!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> ART Module saddles 
module saddles

  implicit none
  
  integer :: natom_displaced    ! # of local atoms displaced
  integer :: KTER_MIN
  integer :: MAXKTER
  integer :: MAXIPERP, MAXKPERP
  integer, dimension(:), allocatable  :: atom_displaced  ! Id of local atoms displaced
  real(kind=8), dimension(:), allocatable, target :: dr

  real(kind=8) :: INITSTEPSIZE 
  real(kind=8) :: LOCAL_CUTOFF
  real(kind=8) :: INCREMENT,BASIN_FACTOR 
  real(kind=8) :: FTHRESHOLD
  real(kind=8) :: EIGEN_THRESH 
  real(kind=8) :: EXITTHRESH

  character(len=20) :: TYPE_EVENTS

  real(kind=8) :: delta_e  ! current_energy - ref_energy
  integer :: m_perp   ! Accepted perpendicular iterations.

  integer :: try      ! Total # of iterationes in a given hyperplane
  real(kind=8) :: ftot     ! Norm of the total force...
  real(kind=8) :: fpar     ! Parallel projection of force
  real(kind=8) :: fperp    ! Norm of the perpendicular force
  real(kind=8) :: delr     ! Distance between two configurations.
  integer :: npart    ! number of particles having moved by more 
                      ! than a THRESHOLD
  logical :: switchDIIS
  logical :: end_activation
  integer :: nsteps_after_eigen_min 
  real(kind=8) :: eigen_min

  !____DEV              
  real(kind=8) :: coord_length
  integer      :: coord_number
  logical      :: cw_try_again ! for clean_wf 

  ! Variables for GUESS_DIRECTION 
  character(len=20)                       :: GUESSFILE
  real(kind=8), dimension(:), allocatable :: g_pos    ! Working positions of the atoms at the presumed direction 
  real(kind=8)                            :: guess_noise

END MODULE saddles

!> ART find_saddle
!!   This subroutine initiates the random displacement at the start
!!   of the ART algorithm. 
!!   After a random escape direction has been selected, the routine call 
!!   saddle_converge which will try to bring the configuration to a saddle point.
!!   If the convergence fails, find_saddle will restart; if it succeeds, it saves 
!!   the event and returns to the main loop.
!!   The displacement can be either LOCAL and involve a single atom
!!   and its neareast neighbours or...
!!
!!   GLOBAL and involve ALL the atoms.
!!
!!   For large cells, it is preferable to use a local initial
!!   displacement to prevent the creation of many trajectories
!!   at the same time in different sections of the cell.
subroutine find_saddle( success, saddle_energy )

  use defs
  use saddles, only : type_events
  implicit none

  !Arguments
  logical, intent(out)      :: success
  real(kind=8), intent(out) :: saddle_energy

  !Local variables
  integer :: ret, ierror, nat

  if ( ( .not. restart ) .and. new_event ) then 

     evalf_number = 0                 ! Initialization of counter.

     nat = 3 * NATOMS                         

     ! We start from the reference configuration.
     pos   = posref  

  ! _______
     ! for dual_search:
     central_atom = 0                 ! default for central_atom 

  !________

     ! These subroutines modify the vector pos and generate a vector
     ! of length 1 indicating the direction of the random
     ! displacement.  There is not case default.

     selectcase( type_events )
         case( 'global' )
             call global_move( )
         case( 'local' )
             call local_move( )
         case( 'list_local' ) 
             call list_and_local( )
         case( 'list' )
             call list_of_atoms( )
         case( 'local_coord' )
             call coord_based_move( )
     end select

  end if

  ! Now, activate per se.
  call saddle_converge( ret, saddle_energy )
  
  ! If the activation did not converge, for whatever reason, we
  ! restart the routine and do not accept the new position.
  call end_report( success, ret, saddle_energy )


END SUBROUTINE find_saddle

!> ART global_move
!!   The initial random direction is taken from the full 3N-dimensional space

!  !Generate a random displacement.

subroutine global_move()

  use defs
  use random
  use saddles

  implicit none 
   
  real(kind=8) :: ran3
  real(kind=8) :: p
  
  selectcase ( search_strategy )

  case ( 0 ) ! default 
        call set_move_random()
  
  case ( 1 ) ! "follow" strategy
        call set_move_follow()

  case ( 2 ) ! "avoid" strategy
        call set_move_avoid()

  case ( 3 ) ! "follow_and_avoid" strategy
        p = ran3()
        if (p > odds_follow_or_avoid) then
                call set_move_follow()
        else
                call set_move_avoid()
        endif

  case ( 4 ) ! "noncovalent"
        call noncovalent()

  case ( 5 ) ! "noncovalent_roll"
        call noncovalent_roll()

  case ( 6 ) ! "noncovalent_attack"
        call noncovalent_attack()

  case ( 7 ) ! "noncovalent_roll_and_attack"
        p = ran3()
        if (p > odds_roll_or_attack) then
                call noncovalent_roll()
        else
                call noncovalent_attack()
        endif

  endselect

  call center_and_norm ( INITSTEPSIZE )

END SUBROUTINE global_move


!> ART local_move
!!   The initial random direction is taken from a restricted space based on 
!!   the local bonding environment. For this, we need to know the list of neighbours
!!   and the cut-off of the potential. Other approaches could be used also.
subroutine local_move( )

  use defs
  use random
  use saddles
  implicit none

  !Local variables
  integer :: i, j, that, ierror
  real(kind=8) :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8) :: dr2
  real(kind=8) :: xi, yi, zi, xij, yij, zij
  real(kind=8) :: ran3
  real(kind=8), dimension(:), pointer  :: dx, dy, dz

                                      ! Select an atom at random.
  if ( preferred_atom < 0 ) then
     ! Only between totally free atoms 
     that = int( NATOMS * ran3() + 1 ) 
  else 
     that = preferred_atom
  end if

  ! for dual_search:
  if ( dual_search ) then
     central_atom = that 
     call neighbours_local( )
  end if 

  call symmetry_break( )              !Breaks the initial symmetry.

                                      ! Write

  open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
  write(FLOG,*) ' '
  write(FLOG,'(1X,A34,I17)') ' - That atom                    : ', that
  close(FLOG)
  write(*,*) 'ARTGAUSS: That atom = ', that


                                      ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF  

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  dr = 0.0d0  
  natom_displaced = 0
  atom_displaced  = 0  

  ! Now we also displace all atoms within a cut-off distance,
  ! LOCAL_CUTOFF.
  xi = x(that)
  yi = y(that)
  zi = z(that)
  do j = 1, NATOMS
        xij = x(j) - xi 
        yij = y(j) - yi
        zij = z(j) - zi
        
        dr2 = xij*xij + yij*yij + zij*zij
        if ( dr2 < lcutoff2 ) then ! Close enough, give a random displacement.
           do
              dx(j) = 0.5d0 - ran3()
              dy(j) = 0.5d0 - ran3()
              dz(j) = 0.5d0 - ran3()

              ! Ensures that the random displacement is isotropic
              dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
              if ( dr2 < 0.25d0 ) exit  
           end do
              
           natom_displaced = natom_displaced + 1
           atom_displaced(j) = 1
        end if
   !  end if
  end do

  call center_and_norm ( INITSTEPSIZE )

END SUBROUTINE local_move


!> ART list_of_atoms
subroutine list_of_atoms ( ) 

  use defs
  use random
  use saddles
  implicit none

  !Local variables  
  logical :: found
  integer :: i, j, i_stat, nlines, sizelist
  integer, dimension(:), allocatable :: list_atoms
  character(len = 128)                            :: filename
  character(len = 150)                            :: line
  character(len = 150), dimension(:), allocatable :: lines
  real(kind=8)                        :: ran3, dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  found    = .false.
  filename = 'list_atoms.dat'
  filename = adjustl( filename )

  ! Is there the file?  
  inquire( file = filename, exist = found )
  if ( .not. found ) then
     write(*,'(a)') 'ERROR: list_atoms.dat not found !! '
     call end_art() 
  end if
  
  open( unit= 99, file= filename, status= 'old' )
   
  ! First pass: to store of the file in a string buffer. 
  ! WARNING: Dimension of lines is 500.
  allocate (lines(500))
  nlines = 1
  do
     read( 99,'(a150)', iostat = i_stat) lines(nlines)
     if (i_stat /= 0) then
        exit
     end if
     nlines = nlines + 1
     if ( nlines > 500 ) then
        write(*,*) 'list_atoms.dat file too long (> 500 lines).'
        write(*,*) 'change this parameter in the source!!'
        call end_art() 
     end if
  end do
  nlines = nlines - 1
  close (99)
 
  if ( nlines < 2 ) then
     write(*,*) 'ERROR: list_atoms, file has less than 2 lines.'
     call end_art() 
  end if
 
  ! Second pass: to determine the correct number of atoms.
  sizelist = 0 
  do i = 1, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (len(trim(line)) /= 0 ) then 
        sizelist = sizelist + 1 
     end if
  end do

  allocate( list_atoms ( sizelist ) )

  ! Last pass. We store the atoms in list_atoms array.
  list_atoms = 0
  do i = 1, sizelist, 1
     write(line, "(a150)") adjustl(lines(i))
     read(line,*, iostat = i_stat)  list_atoms(i) 
  end do
  deallocate(lines)

  call symmetry_break( )              !Breaks the initial symmetry.

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0  
  
  ! Now we apply a random displacement 
  do j = 1, NATOMS
        do i = 1, sizelist 
           if ( j == list_atoms(i) ) then
              do
                 dx(j) = 0.5d0 - ran3()
                 dy(j) = 0.5d0 - ran3()
                 dz(j) = 0.5d0 - ran3()                 
                 dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2

                 ! Ensures that the random displacement is isotropic
                 if ( dr2 < 0.25d0 ) exit
              end do
              natom_displaced = natom_displaced + 1
              atom_displaced(j) = 1
              exit
           end if
        end do
  end do
  
  deallocate(list_atoms)

  call center_and_norm ( INITSTEPSIZE )

END SUBROUTINE list_of_atoms 


!> ART list_and_local
subroutine list_and_local () 

  use defs
  use random
  use saddles
  implicit none

  !Local variables  
  logical :: found
  integer :: i, j, i_stat, nlines, sizelist
  integer :: that, this, ierror
  integer, dimension(:), allocatable :: list_atoms
  character(len = 128)                            :: filename
  character(len = 150)                            :: line
  character(len = 150), dimension(:), allocatable :: lines
  real(kind=8) :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8) :: ran3, dr2
  real(kind=8) :: xi, yi, zi, xij, yij, zij
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  found    = .false.
  filename = 'list_atoms.dat'
  filename = adjustl( filename )

  ! Is there the file?  
  inquire( file = filename, exist = found )
  if ( .not. found ) then
     write(*,'(a)') 'ERROR: list_atoms.dat not found !! '
     call end_art()  
  end if
  
  open( unit= 99, file= filename, status= 'old' )
   
  ! First pass: to store of the file in a string buffer. 
  ! WARNING: Dimension of lines is 500.
  allocate (lines(500))
  nlines = 1
  do
     read( 99,'(a150)', iostat = i_stat) lines(nlines)
     if (i_stat /= 0) then
        exit
     end if
     nlines = nlines + 1
     if ( nlines > 500 ) then
        write(*,*) 'list_atoms.dat file too long (> 500 lines).'
        write(*,*) 'change this parameter in the source!!'
        call end_art() 
     end if
  end do
  nlines = nlines - 1
  close (99)
 
  if ( nlines < 2 ) then
     write(*,*) 'ERROR: list_atoms, file has less than 2 lines.'
     call end_art() 
  end if
 
  ! Second pass: to determine the correct number of atoms.
  sizelist = 0 
  do i = 1, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (len(trim(line)) /= 0 ) then 
        sizelist = sizelist + 1 
     end if
  end do

  allocate( list_atoms ( sizelist ) )

  ! Last pass. We store the atoms in list_atoms array.
  list_atoms = 0
  do i = 1, sizelist, 1
     write(line, "(a150)") adjustl(lines(i))
     read(line,*, iostat = i_stat)  list_atoms(i) 
  end do
  deallocate(lines)


  do 
     this = int( sizelist * ran3() + 1 ) 
     that = list_atoms(this)
  end do


  ! for dual_search:
  if ( dual_search ) then
     central_atom = that 
     call neighbours_local( )
  end if 

  call symmetry_break( )              !Breaks the initial symmetry.

  ! Write
  open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
  write(FLOG,*) ' '
  write(FLOG,'(1X,A34,I17)') ' - That atom                    : ', that
  close(FLOG)
  write(*,*) 'ARTGAUSS: That atom = ', that
  
  
                                      ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF  

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  dr = 0.0d0  
  natom_displaced = 0
  atom_displaced  = 0  

  ! Now we also displace all atoms within a cut-off distance,
  ! LOCAL_CUTOFF.
  xi = x(that)
  yi = y(that)
  zi = z(that)
  do j = 1, NATOMS
        xij = x(j) - xi
        yij = y(j) - yi
        zij = z(j) - zi
        
        dr2 = xij*xij + yij*yij + zij*zij
        if ( dr2 < lcutoff2 ) then ! Close enough, give a random displacement.
           do
              dx(j) = 0.5d0 - ran3()
              dy(j) = 0.5d0 - ran3()
              dz(j) = 0.5d0 - ran3()
              
              ! Ensures that the random displacement is isotropic
              dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
              if ( dr2 < 0.25d0 ) exit  
           end do
           natom_displaced = natom_displaced + 1
           atom_displaced(j) = 1
        end if
  end do

  call center_and_norm ( INITSTEPSIZE )

END SUBROUTINE list_and_local  


!> ART symmetry_break
subroutine symmetry_break( )

  use defs
  use random
  use saddles
  implicit none

  !Local variables
  integer :: i
  real(kind=8) :: dr2
  real(kind=8) :: ran3
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers.
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0                          ! Initialization of dr.

  ! Generate a random displacement. 
  do i = 1, NATOMS
        do
           dx(i) = 0.5d0 - ran3()
           dy(i) = 0.5d0 - ran3()
           dz(i) = 0.5d0 - ran3()

           ! Ensures that the random displacement is isotropic
           dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
           if ( dr2 < 0.25d0 ) exit 
        end do
        natom_displaced = natom_displaced + 1
        atom_displaced(i) = 1
  end do

  call center_and_norm ( sym_break_dist )

END SUBROUTINE symmetry_break


!> ART center_and_norm
subroutine center_and_norm ( step )

  use defs
  use saddles
  implicit none

  !Arguments
  real(kind=8), intent(in) :: step

  !Local variables 
  integer :: nat
  real(kind=8) :: ierror
  real(kind=8) :: xsum, ysum, zsum
  real(kind=8) :: xnorm, ynorm, znorm, norm
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)


  ! We now center only on the atoms that have been randomly selected.
  if ( natom_displaced > 1 ) then
     xsum = sum(dx) / natom_displaced
     ysum = sum(dy) / natom_displaced
     zsum = sum(dz) / natom_displaced

     dx = dx - xsum * atom_displaced
     dy = dy - ysum * atom_displaced
     dz = dz - zsum * atom_displaced
     
  end if
        write(*,*) "This is sum_coords", xsum, ysum, zsum
  
  ! And normalize the total random displacement effected to the value
  ! desired.
  norm = 0.0d0
  norm = dot_product(dr, dr)
  
  ! This renormalizes in angstroems to a displacement INITSTEPSIZE.
  norm = step / sqrt(norm)
  xnorm = norm 
  ynorm = norm
  znorm = norm
  
  dx = dx * xnorm 
  dy = dy * ynorm
  dz = dz * znorm

  ! Redistribute the information.
  nat = 3*NATOMS

  ! Update the position using this random displacement
  pos = pos + dr
  
  ! Now, we normalize dr to get the initial_direction 
  initial_direction = dr 
  norm = dot_product( initial_direction, initial_direction ) 
  norm = 1.0d0 / sqrt(norm)
  initial_direction  = initial_direction * norm

  write(*,*) 'ARTGAUSS: Number of displaced atoms initially: ',natom_displaced

  deallocate(atom_displaced)
  deallocate(dr)

END SUBROUTINE center_and_norm 


! for dual_search:
! Based on the Laurent Karim Beland's subroutine neighbours in neighbour.f9O.
! But here, what we only want is a list of neighbours of a given atom.
subroutine neighbours_local( )

  use defs
  implicit none
 
  !Local variable
  integer i, j
  real(kind=8) :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8) :: dr2
  real(kind=8) :: xi, yi, zi, xij, yij, zij 

  !_______________________

  lcutoff2 = size_system*size_system

  i = central_atom 

  xi = x(i)
  yi = y(i)
  zi = z(i)

  do j = 1, NATOMS
     if ( j == i ) cycle

     xij = x(j) - xi
     yij = y(j) - yi
     zij = z(j) - zi

     dr2 = xij*xij + yij*yij + zij*zij
  end do

END SUBROUTINE neighbours_local


subroutine coord_based_move( ) 

  use defs
  use random
  use saddles
  implicit none
 
  !Local variable
  integer :: i, j, that, ierror
  integer :: numnei 
  real(kind=8) :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8) :: dr2
  real(kind=8) :: xi, yi, zi, xij, yij, zij 
  real(kind=8) :: ran3
  real(kind=8), dimension(:), pointer  :: dx, dy, dz

  integer, dimension(natoms) :: in_list

  !_______________________

  lcutoff2 = coord_length*coord_length

  in_list = 0                       ! Initial vectorial assignment
 
  do i = 1, NATOMS
     
     numnei = 0
     xi = x(i)
     yi = y(i)
     zi = z(i)

     do j = 1, NATOMS
        if ( j == i ) cycle

        xij = x(j) - xi
        yij = y(j) - yi
        zij = z(j) - zi

        dr2 = xij*xij + yij*yij + zij*zij

        if ( dr2 < lcutoff2 ) then
           numnei = numnei + 1
        end if

     end do

     if ( numnei < coord_number ) in_list(i) = 1 

  end do

  ! Only between totally free atoms 
  do 
  that = int( NATOMS * ran3() + 1 ) 
     if ( in_list(that) == 1) then !
        if ( typat(that) == type_sel .and. type_sel/= '' ) then
           exit 
        else if ( type_sel == '' ) then 
           exit
        else 
           cycle
        end if

     else 
        cycle
     end if
  end do
 
  if ( that == 0 ) then
     write(*,*)  'ARTGAUSS ERROR: There is no atom with lower coord than ', coord_number
     stop
  end if

  call symmetry_break( )              !Breaks the initial symmetry.

  open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
  write(FLOG,*) ' '
  write(FLOG,'(1X,A34,I17)') ' - That atom                    : ', that
  close(FLOG)
  write(*,*) 'ARTGAUSS: That atom = ', that

  ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF  

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
  
  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  dr = 0.0d0  
  natom_displaced = 0
  atom_displaced  = 0  

  ! Now we also displace all atoms within a cut-off distance,
  ! LOCAL_CUTOFF.
  xi = x(that)
  yi = y(that)
  zi = z(that)
  do j = 1, NATOMS
        xij = x(j) - xi
        yij = y(j) - yi
        zij = z(j) - zi

        dr2 = xij*xij + yij*yij + zij*zij
        if ( dr2 < lcutoff2 ) then ! Close enough, give a random displacement.
           do
              dx(j) = 0.5d0 - ran3()
              dy(j) = 0.5d0 - ran3()
              dz(j) = 0.5d0 - ran3()
              ! Ensures that the random
              ! displacement is isotropic
              dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
              if ( dr2 < 0.25d0 ) exit  
           end do
           natom_displaced = natom_displaced + 1
           atom_displaced(j) = 1
        end if
  end do

  call center_and_norm ( INITSTEPSIZE )

END SUBROUTINE coord_based_move 


SUBROUTINE set_move_random()

  use defs
  use random
  use saddles

  implicit none

  !Local variables
  integer :: i
  real(kind=8) :: dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8) :: ran3

  allocate(dr(3*natoms)) 
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)
  
  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0

  do i = 1, natoms
      do
         dx(i) = 0.5d0 - ran3()
         dy(i) = 0.5d0 - ran3()
         dz(i) = 0.5d0 - ran3()                             
         ! Ensures that the random displacement is isotropic
         dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
         if ( dr2 < 0.25d0 ) exit 
      end do
      natom_displaced = natom_displaced + 1
      atom_displaced(i) = 1
  end do

END SUBROUTINE set_move_random


SUBROUTINE set_move_follow()

  use defs
  use random
  use saddles

  implicit none

  !Local variables
  integer :: i, success_counter, rand
  real(kind=8) :: dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8) :: ran3
  real(kind=8) :: dr_transformed_list(nsad_read,natoms,3)
  real(kind=8) :: sad_transformed_list(nsad_read,natoms,3)

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)
  
  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0

  call read_and_transform ( dr_transformed_list, sad_transformed_list, success_counter )
  if (success_counter .eq. 0) then
          write(*,*) "None of the read min align well with the current min"
          call end_art()
  endif
  rand = ceiling(success_counter*ran3()) 
  write(*,*) "This is the random pick", rand

  do i = 1, natoms
          dx(i) = dr_transformed_list(rand,i,1)
          dy(i) = dr_transformed_list(rand,i,2)
          dz(i) = dr_transformed_list(rand,i,3)
          natom_displaced = natom_displaced + 1
          atom_displaced(i) = 1
  enddo

END SUBROUTINE set_move_follow


SUBROUTINE set_move_avoid()

  use defs
  use random
  use saddles

  implicit none

  !Local variables
  integer :: i, j, k, success_counter, rand
  real(kind=8) :: dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8) :: cos_theta
  real(kind=8) :: ran3
  real(kind=8) :: dr_transformed_list(nsad_read,natoms,3)
  real(kind=8) :: sad_transformed_list(nsad_read,natoms,3)

  allocate(dr(3*natoms)) 
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)
  
  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0

  call read_and_transform ( dr_transformed_list, sad_transformed_list, success_counter )
  if (success_counter .eq. 0) then
          write(*,*) "None of the read min align well with the current min"
          call end_art()
  endif
  rand = ceiling(success_counter*ran3()) 
  do
      do i = 1, natoms 
          do
                   dx(i) = 0.5d0 - ran3()
                   dy(i) = 0.5d0 - ran3()
                   dz(i) = 0.5d0 - ran3()                             
                   dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
                   if ( dr2 < 0.25d0 ) exit 
          end do
          natom_displaced = natom_displaced + 1
          atom_displaced(i) = 1
      end do
      cos_theta = dot_product(dr,reshape(dr_transformed_list(rand,1:natoms,1:3),(/natoms*3/))) / &
                & sqrt(dot_product(dr,dr)*dot_product(reshape(dr_transformed_list(rand,1:natoms,1:3),(/natoms*3/)),&
                & reshape(dr_transformed_list(rand,1:natoms,1:3),(/natoms*3/))))
      if (cos_theta .LT. 0.8) exit
  enddo

END SUBROUTINE set_move_avoid


SUBROUTINE read_and_transform ( dr_transformed_list, sad_transformed_list, success_count )

        use defs
        
        implicit none

        real(kind=8), intent(out) :: dr_transformed_list(nsad_read,natoms,3)
        real(kind=8), intent(out) :: sad_transformed_list(nsad_read,natoms,3)
        integer, intent(out) :: success_count
        
        integer :: i, j, min_count, sad_count

        real(kind=8) :: read_min(nmin_read,natoms_read,3)
        real(kind=8) :: read_sad(nsad_read,natoms_read,3)
        real(kind=8) :: read_dr(nsad_read,natoms_read,3)
        real(kind=8) :: current_min(natoms,3)
        real(kind=8) :: each_read_min(natoms_read,3)
        real(kind=8) :: each_read_dr(natoms_read,3)
        real(kind=8) :: each_read_sad(natoms_read,3)
        real(kind=8) :: each_dr_transformed(natoms,3)
        real(kind=8) :: each_sad_transformed(natoms,3)
        logical :: align_well
        character(len=20) :: keyword 

        open(VREAD, file=VECREAD, status = 'old', action = 'read')

        min_count = 0
        sad_count = 0

        do
                read(VREAD,*,end = 100) keyword
                if (keyword .eq. 'min') then
                         min_count = min_count + 1
                         do i = 1, natoms_read
                                   read(VREAD,*) typat_read(i), read_min(min_count,i,1), read_min(min_count,i,2), read_min(min_count,i,3)
                         enddo
                endif
                if (keyword .eq. 'sad') then
                         sad_count = sad_count + 1
                         do i = 1, natoms_read
                                   read(VREAD,*) typat_read(i), read_sad(sad_count,i,1), read_sad(sad_count,i,2), read_sad(sad_count,i,3)
                         enddo
                endif
        enddo
        100 rewind(VREAD)

        close(VREAD)

        do j = 1, nsad_read
                do i = 1, natoms_read
                         read_dr(j,i,1:3) = read_sad(j,i,1:3) - read_min(j,i,1:3)
                enddo
        enddo

        do i = 1, natoms
                current_min(i,1) = x(i)
                current_min(i,2) = y(i)
                current_min(i,3) = z(i)
        enddo
        
        success_count = 0
        do j = 1, nsad_read
                do i = 1, natoms_read
                         each_read_min(i,1:3) = read_min(j,i,1:3)
                         each_read_sad(i,1:3) = read_sad(j,i,1:3)
                         each_read_dr(i,1:3) = read_dr(j,i,1:3)
                enddo

                align_well = .false.
                call align ( each_read_min, current_min, align_well, each_read_sad, each_read_dr, each_dr_transformed, each_sad_transformed )
                if ( align_well .eq. .true. ) then                        
                         success_count = success_count + 1
                         do i = 1, natoms
                                 dr_transformed_list(success_count,i,1:3) = each_dr_transformed(i,1:3)
                                 sad_transformed_list(success_count,i,1:3) = each_sad_transformed(i,1:3)
                         enddo
                endif
        enddo

        write(*,*) "Success_count inside subroutine read_and_transform", success_count

END SUBROUTINE read_and_transform


SUBROUTINE detect_fragments ( number_of_fragments, fragment_list )
        
        use defs
        
        implicit none

        !Arguments
        integer, intent(out) :: number_of_fragments 
        integer, dimension(natoms, natoms), intent(out) :: fragment_list

        !Local variables
        integer :: i, j, k
        real(kind=8), dimension(natoms, natoms) :: dist_matrix 
        logical, dimension(natoms, natoms) :: adj_matrix 
        character(len=1), dimension(4) :: atomic_kind
        real(kind=8), dimension(4) :: cov_rad
        real(kind=8), dimension(natoms) :: cov_radius_current
        logical, dimension(natoms) :: visited
        integer, dimension(natoms) :: queue
        integer :: node, q_index, unvisited

        fragment_list = 0      !initializing fragment list

        do i = 1, natoms
                do j = 1, natoms
                      dist_matrix(i,j) = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
                enddo
                write(*,'(<natoms>(f14.8))') dist_matrix(i,1:natoms)
        enddo


        atomic_kind = [ 'C',  'H',  'N',  'O' ]            !Theoretical covalent radii of the
        cov_rad     = [ 0.76, 0.31, 0.71, 0.66 ]           !four most common organic elements 

        !Every atom in the current structure is assigned a covalent radius

        do i = 1, natoms
                do j = 1, 4
                        if (typat(i) == atomic_kind(j)) then
                                cov_radius_current(i) = cov_rad(j)
                        endif
                enddo
        enddo

        do i = 1, natoms
                do j = 1, natoms
                        if ( i .ne. j ) then
                                if (abs(dist_matrix(i,j) - (cov_radius_current(i) + cov_radius_current(j)) < 0.30)) then
                                        adj_matrix(i,j) = .true.
                                else
                                        adj_matrix(i,j) = .false.
                                endif
                        else
                                adj_matrix(i,j) = .false.
                        endif
                enddo
        enddo

        ! Breadth-first search (BFS) algorithm for finding fragments

        do i = 1, natoms                        !Mark all nodes as unvisited
                visited(i) = .false.
        enddo
        
        queue = 0
        visited(1) = .true.                     !Visit the source node (could be any node, I chose number 1)
        queue(1) = 1                            !Enqueue the visited node 
        node = queue(1)                         !The variable "node" holds the first element of the queue     
        q_index = 1
        number_of_fragments = 0
        do
                do j = 1, natoms               
                        if (adj_matrix(node,j) .eq. .true. .and. visited(j) .eq. .false.) then  
                                visited(j) = .true.                         
                                do k = 1, size(queue)
                                        if (queue(k) .eq. 0) then           
                                                queue(k) = j                
                                                exit
                                        endif
                                enddo
                        endif
                enddo
                q_index = q_index + 1
                node = queue(q_index)            
                if (q_index .eq. natoms) then
                        number_of_fragments = number_of_fragments + 1
                        fragment_list(number_of_fragments, 1:natoms) = queue(1:natoms)
                        exit
                endif
                if (queue(q_index) .eq. 0) then  
                        number_of_fragments = number_of_fragments + 1
                        fragment_list(number_of_fragments, 1:natoms) = queue(1:natoms)
                        queue(1:natoms) = 0
                        
                        unvisited = 0
                        do i = 1, natoms
                                if (visited(i) .eq. .false.) then
                                        unvisited = unvisited + 1
                                endif
                        enddo
                        if (unvisited .eq. 0) exit
                        do i = 1, natoms
                                if (visited(i) .eq. .false.) then
                                        visited(i) = .true.
                                        queue(1) = i 
                                        node = queue(1)
                                        q_index = 1
                                        exit
                                endif
                        enddo

                endif

                write(*,'(<natoms>(1X, I2))') queue
        enddo
        
        write(*,*) "fragment"
        do i = 1, natoms
                write(*,'(<natoms>(2X, I2))') fragment_list(i, 1:natoms)
        enddo

END SUBROUTINE detect_fragments


SUBROUTINE fragment_utility ( fragment_to_move, shortest_vec )
        
        use defs

        implicit none

        !Arguments

        real(kind=8), dimension(1,3), intent(out) :: shortest_vec
        integer, dimension(natoms), intent(out) :: fragment_to_move
        
        !Local variables
        integer :: i, j        
        integer :: number_of_fragments
        integer, dimension(natoms, natoms) :: fragment_list
        integer, allocatable, dimension(:) :: frag_size_list
        real(kind=8), allocatable, dimension(:,:) :: dist_between_frag
        integer :: each_frag_size, smallest_fragment, largest_fragment

        fragment_list = 0

        call detect_fragments ( number_of_fragments, fragment_list )
        
        allocate(frag_size_list(number_of_fragments))

        do i = 1, number_of_fragments
                each_frag_size = 0
                do j = 1, natoms
                        if (fragment_list(i,j) .ne. 0) then
                                each_frag_size = each_frag_size + 1
                        endif
                enddo
                frag_size_list(i) = each_frag_size
        enddo

        do i = 1, number_of_fragments
                if (minval(frag_size_list) .eq. frag_size_list(i)) then
                        smallest_fragment = i 
                elseif (maxval(frag_size_list) .eq. frag_size_list(i)) then
                        largest_fragment = i
                endif
        enddo
        
        fragment_to_move(1:natoms) = fragment_list(smallest_fragment, 1:natoms)
        
        allocate(dist_between_frag(frag_size_list(smallest_fragment), frag_size_list(largest_fragment)))        

        do i = 1, frag_size_list(smallest_fragment)
            do j = 1, frag_size_list(largest_fragment)
                dist_between_frag(i, j) = sqrt((x(fragment_list(smallest_fragment, i)) - x(fragment_list(largest_fragment, j)))**2 &
                                           & + (y(fragment_list(smallest_fragment, i)) - y(fragment_list(largest_fragment, j)))**2 &
                                           & + (z(fragment_list(smallest_fragment, i)) - z(fragment_list(largest_fragment, j)))**2) 
            enddo
        enddo

        do i = 1, frag_size_list(smallest_fragment)
                do j = 1, frag_size_list(largest_fragment)
                        if (minval(dist_between_frag) .eq. dist_between_frag(i,j)) then
                                shortest_vec(1,1) = x(fragment_list(smallest_fragment, i)) - x(fragment_list(largest_fragment, j))
                                shortest_vec(1,2) = y(fragment_list(smallest_fragment, i)) - y(fragment_list(largest_fragment, j))
                                shortest_vec(1,3) = z(fragment_list(smallest_fragment, i)) - z(fragment_list(largest_fragment, j))
                        endif
                enddo
        enddo
        
        do i = 1, frag_size_list(smallest_fragment)
                write(*,'(<frag_size_list(largest_fragment)>(f14.3))') dist_between_frag(i, 1:frag_size_list(largest_fragment)) 
        enddo

        deallocate(frag_size_list)
        deallocate(dist_between_frag)

END SUBROUTINE fragment_utility


SUBROUTINE noncovalent()
                
  use defs
  use random
  use saddles

  implicit none

  !Local variables
  integer :: i, j
  real(kind=8) :: dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8) :: ran3
  integer, dimension(natoms) :: fragment_to_move
  real(kind=8), dimension(1,3) :: shortest_vec 

  allocate(dr(3*natoms)) 
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)
  
  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0

  call fragment_utility ( fragment_to_move, shortest_vec )
  
  write(*,*) "This is fragment to move", fragment_to_move

  do i = 1, natoms
      do j = 1, natoms
              if (i .eq. fragment_to_move(j)) then
                do
                        dx(i) = 0.5d0 - ran3()
                        dy(i) = 0.5d0 - ran3()
                        dz(i) = 0.5d0 - ran3()                             
                        ! Ensures that the random displacement is isotropic
                        dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
                        if ( dr2 < 0.25d0 ) exit 
                end do
                natom_displaced = natom_displaced + 1
                atom_displaced(i) = 1
              endif
      enddo
  enddo

  write(*,*) "This is dr_noncovalent"
  do i = 1, natoms
          write(*,'(3(f14.6))') dx(i), dy(i), dz(i)
  enddo

  write(*,*) "Atoms displaced", natom_displaced
        
END SUBROUTINE noncovalent


SUBROUTINE noncovalent_roll()

  use defs
  use random
  use saddles

  implicit none

  !Local variables
  integer :: i, j
  real(kind=8) :: dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8) :: ran3
  integer, dimension(natoms) :: fragment_to_move
  real(kind=8), dimension(1,3) :: shortest_vec 

  allocate(dr(3*natoms)) 
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)
  
  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0

  call fragment_utility ( fragment_to_move, shortest_vec )
  
  do i = 1, natoms
      do j = 1, natoms
              if (i .eq. fragment_to_move(j)) then
                dx(i) = 0.5d0 - ran3()
                dy(i) = 0.5d0 - ran3()
                dz(i) = -1.0*(shortest_vec(1,1) * dx(i) + shortest_vec(1,2) * dy(i)) / shortest_vec(1,3)                             
                natom_displaced = natom_displaced + 1
                atom_displaced(i) = 1
              endif
      enddo
  enddo

  write(*,*) "This is dr_noncovalent_roll"
  do i = 1, natoms
          write(*,'(3(f14.6))') dx(i), dy(i), dz(i)
  enddo

END SUBROUTINE noncovalent_roll


SUBROUTINE noncovalent_attack()

  use defs
  use random
  use saddles

  implicit none

  !Local variables
  integer :: i, j
  real(kind=8) :: dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8) :: ran3
  integer, dimension(natoms) :: fragment_to_move
  real(kind=8), dimension(1,3) :: shortest_vec 

  allocate(dr(3*natoms)) 
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)
  
  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0

  call fragment_utility ( fragment_to_move, shortest_vec )

  do i = 1, natoms
      do j = 1, natoms
              if (i .eq. fragment_to_move(j)) then
                dx(i) = shortest_vec(1,1)
                dy(i) = shortest_vec(1,2)
                dz(i) = shortest_vec(1,3)
              endif
      enddo
      natom_displaced = natom_displaced + 1
      atom_displaced(i) = 1
  enddo

  write(*,*) "This is dr_noncovalent_attack"
  do i = 1, natoms
          write(*,'(3(f14.6))') dx(i), dy(i), dz(i)
  enddo

END SUBROUTINE noncovalent_attack


SUBROUTINE align ( read_min, current_min, align_well, read_sad, read_dr, dr_transformed, sad_transformed )

        ! This subroutine calculates a transformation matrix that minimizes the RMSD between two sets of coordinates.
        ! The first step of this transformation involves a translation. This is simply done by re-centering the molecules so that their
        ! centroids coincide with each other.
        ! Next, a covariance matrix is calculated using the re-centered molecules.
        ! Then, the singular value decomposition (SVD) of this covariance matrix is calculated, which yields the factors of this
        ! covariance, viz U, S and V (Here U and V are the left- and right- singular vectors of the covariance matrix.)
        ! According to the Kabsch algorithm, the product of V and U(transpose) is the optimal rotation matrix 
        ! that minimizes the RMSD between the two sets of vectors.
        ! Here cov is the covariance matrix, U and VT are the left- and
        ! right- singular vectors of the cov, (VT is actually the transpose of V),
        ! S is one of the factors of the cov. R (V*U_transpose) is the rotation matrix
        ! Work is one of the arguments required by the LAPACK library subroutine that calculates the SVD

        use defs

        implicit none

        !Arguments
        real(kind=8), intent(in) :: read_min(natoms_read,3)
        real(kind=8), intent(in) :: current_min(natoms,3)
        real(kind=8), intent(in) :: read_dr(natoms_read,3)
        real(kind=8), intent(in) :: read_sad(natoms_read,3)

        logical, intent(inout) :: align_well

        real(kind=8), intent(out) :: dr_transformed(natoms,3)
        real(kind=8), intent(out) :: sad_transformed(natoms,3)

        !Local variables
        real(kind=8) :: read_min_sliced(natoms_correspond,3)
        real(kind=8) :: current_min_sliced(natoms_correspond,3)
        real(kind=8) :: read_dr_sliced(natoms_correspond,3)
        real(kind=8) :: read_sad_sliced(natoms_correspond,3)
        real(kind=8) :: read_dr_translated(natoms_correspond,3)
        real(kind=8) :: read_sad_translated(natoms_correspond,3)
        real(kind=8) :: read_dr_translated_rotated(natoms_correspond,3)
        real(kind=8) :: read_sad_translated_rotated(natoms_correspond,3)


        integer :: read_atoms_correspond(natoms_correspond)
        integer :: current_atoms_correspond(natoms_correspond)

        real(kind=8) :: read_min_sliced_centered(natoms_correspond,3)
        real(kind=8) :: current_min_sliced_centered(natoms_correspond,3)
        real(kind=8) :: read_min_sliced_transformed(natoms_correspond,3)

        real(kind=8) :: deviation_atom_before(natoms_correspond) ! deviation of each individual atom from current min before alignment
        real(kind=8) :: deviation_atom_after(natoms_correspond) ! deviation of each individual atom from current min after alignment
        real(kind=8) :: sum_deviations_before, rmsd_before, sum_deviations_after, rmsd_after
        real(kind=8) :: sum_read_x, sum_read_y, sum_read_z, cx_read, cy_read, cz_read
        real(kind=8) :: sum_current_x, sum_current_y, sum_current_z, cx_current, cy_current, cz_current
         
        real(kind=8) :: cov(3,3), U(3,3), S(3), VT(3,3), R(3,3), work(20)
        real(kind=8) :: U_transpose(3,3), V(3,3), det_U_transpose, det_V
       
        integer :: i, j, dev_count, info
        character :: dummy
        
        open(ALI, file = 'alignment.txt', status = 'unknown', action = 'read')
        read(ALI,*) dummy, dummy
        do i = 1, natoms_correspond
                read(ALI,*) current_atoms_correspond(i), read_atoms_correspond(i) 
                write(*,*) current_atoms_correspond(i), read_atoms_correspond(i) 

        enddo
        close(ALI)
        
        do i = 1, natoms_correspond
                read_min_sliced(i,1:3) = read_min(read_atoms_correspond(i),1:3)
                read_sad_sliced(i,1:3) = read_sad(read_atoms_correspond(i),1:3)
                read_dr_sliced(i,1:3) = read_dr(read_atoms_correspond(i),1:3)
                current_min_sliced(i,1:3) = current_min(current_atoms_correspond(i),1:3)
        enddo

        
        sum_deviations_before = 0.0

        do i = 1, natoms_correspond     
                deviation_atom_before(i) = sqrt((read_min_sliced(i,1) - current_min_sliced(i,1))**2 &
                                      & + (read_min_sliced(i,2) - current_min_sliced(i,2))**2 &
                                      & + (read_min_sliced(i,3) - current_min_sliced(i,3))**2) 

                write(*,*) deviation_atom_before(i)
                sum_deviations_before = sum_deviations_before + &
                                    & + (read_min_sliced(i,1) - current_min_sliced(i,1))**2 & 
                                    & + (read_min_sliced(i,2) - current_min_sliced(i,2))**2 & 
                                    & + (read_min_sliced(i,3) - current_min_sliced(i,3))**2       
        enddo
        
        rmsd_before = sqrt(sum_deviations_before/natoms_correspond) !Calculating RMSD before structural alignment
        write(*,*) "This is rmsd before alignment: ", rmsd_before

        sum_current_x = 0.0
        sum_current_y = 0.0
        sum_current_z = 0.0

        sum_read_x = 0.0
        sum_read_y = 0.0
        sum_read_z = 0.0

        do i = 1, natoms_correspond
                sum_read_x = sum_read_x + read_min_sliced(i,1)
                sum_read_y = sum_read_y + read_min_sliced(i,2)
                sum_read_z = sum_read_z + read_min_sliced(i,3)

                sum_current_x = sum_current_x + current_min_sliced(i,1)
                sum_current_y = sum_current_y + current_min_sliced(i,2)
                sum_current_z = sum_current_z + current_min_sliced(i,3)
        enddo

        !Calculating the centroids for the two structures

        cx_read = sum_read_x/natoms_correspond
        cy_read = sum_read_y/natoms_correspond
        cz_read = sum_read_z/natoms_correspond

        cx_current = sum_current_x/natoms_correspond
        cy_current = sum_current_y/natoms_correspond
        cz_current = sum_current_z/natoms_correspond

        !Recentering the structures so that their centroids coincide with each other.

        do i = 1, natoms_correspond
                read_min_sliced_centered(i,1) = read_min_sliced(i,1) - cx_read 
                read_min_sliced_centered(i,2) = read_min_sliced(i,2) - cy_read
                read_min_sliced_centered(i,3) = read_min_sliced(i,3) - cz_read 
        enddo

        do i = 1, natoms_correspond
                current_min_sliced_centered(i,1) = current_min_sliced(i,1) - cx_current
                current_min_sliced_centered(i,2) = current_min_sliced(i,2) - cy_current
                current_min_sliced_centered(i,3) = current_min_sliced(i,3) - cz_current 
        enddo
               
        !Calculating the covariance matrix         

        cov = matmul(transpose(read_min_sliced_centered), current_min_sliced_centered)

        !Invoking the LAPACK library subroutine "dgesvd" that calculates the SVD of the covariance matrix

        call dgesvd('S','S',3,3,cov,3,S,U,3,VT,3,work,20,info)
        U_transpose = transpose(U)
        V = transpose(VT)

        det_U_transpose = U_transpose(1,1)*(U_transpose(2,2)*U_transpose(3,3) - U_transpose(2,3)*U_transpose(3,2)) &
                      & - U_transpose(1,2)*(U_transpose(2,1)*U_transpose(3,3) - U_transpose(2,3)*U_transpose(3,1)) & 
                      & + U_transpose(1,3)*(U_transpose(2,1)*U_transpose(3,2) - U_transpose(2,2)*U_transpose(3,1))

        det_V = V(1,1)*(V(2,2)*V(3,3) - V(2,3)*V(3,2)) &
            & - V(1,2)*(V(2,1)*V(3,3) - V(2,3)*V(3,1)) &
            & + V(1,3)*(V(2,1)*V(3,2) - V(2,2)*V(3,1)) 
        
        if(det_V * det_U_transpose .LT. 0) then    ! If det(V*UT) is less than zero, then multiply the last column of V by -1 (reflection)
                V(:,3) = -V(:,3)              
        endif

        R = matmul(V,U_transpose)  !This is the optimal rotation matrix 


        read_min_sliced_transformed = transpose(matmul(R,transpose(read_min_sliced_centered)))

        do i = 1,natoms_correspond
                read_min_sliced_transformed(i,1) = read_min_sliced_transformed(i,1) + cx_current
                read_min_sliced_transformed(i,2) = read_min_sliced_transformed(i,2) + cy_current
                read_min_sliced_transformed(i,3) = read_min_sliced_transformed(i,3) + cz_current
        enddo
        
        sum_deviations_after = 0.0
        do i = 1, natoms_correspond
                deviation_atom_after(i) = sqrt((read_min_sliced_transformed(i,1) - current_min_sliced(i,1))**2 &
                                      & + (read_min_sliced_transformed(i,2) - current_min_sliced(i,2))**2 &
                                      & + (read_min_sliced_transformed(i,3) - current_min_sliced(i,3))**2) 

                sum_deviations_after = sum_deviations_after &
                                   & + (read_min_sliced_transformed(i,1) - current_min_sliced(i,1))**2 & 
                                   & + (read_min_sliced_transformed(i,2) - current_min_sliced(i,2))**2 & 
                                   & + (read_min_sliced_transformed(i,3) - current_min_sliced(i,3))**2 
        enddo

        rmsd_after = sqrt(sum_deviations_after/natoms_correspond) !Calculating RMSD after structural alignment
        write(*,*) "This is rmsd after alignment: ", rmsd_after
        
        dev_count = 0
        do i = 1, natoms_correspond
                write(*,*) deviation_atom_after(i)
                if ( deviation_atom_after(i) .LT. 0.1 ) then
                        dev_count = dev_count + 1
                endif
        enddo
        
        if (rmsd_after .LT. 0.1 .and. dev_count == natoms_correspond) then

                align_well = .true.

                do i = 1,natoms_correspond
                        read_dr_translated(i,1) = read_dr_sliced(i,1) - cx_read
                        read_dr_translated(i,2) = read_dr_sliced(i,2) - cy_read
                        read_dr_translated(i,3) = read_dr_sliced(i,3) - cz_read

                        read_sad_translated(i,1) = read_sad_sliced(i,1) - cx_read
                        read_sad_translated(i,2) = read_sad_sliced(i,2) - cy_read
                        read_sad_translated(i,3) = read_sad_sliced(i,3) - cz_read
                enddo

                read_dr_translated_rotated = transpose(matmul(R,transpose(read_dr_translated)))
                read_sad_translated_rotated = transpose(matmul(R,transpose(read_sad_translated)))
                
                do i = 1,natoms_correspond
                        dr_transformed(i,1) = read_dr_translated_rotated(i,1) + cx_current
                        dr_transformed(i,2) = read_dr_translated_rotated(i,2) + cy_current
                        dr_transformed(i,3) = read_dr_translated_rotated(i,3) + cz_current

                        sad_transformed(i,1) = read_sad_translated_rotated(i,1) + cx_current
                        sad_transformed(i,2) = read_sad_translated_rotated(i,2) + cy_current
                        sad_transformed(i,3) = read_sad_translated_rotated(i,3) + cz_current
                enddo
                
        endif
        
END SUBROUTINE align
