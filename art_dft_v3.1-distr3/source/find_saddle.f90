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
  real(kind=8), dimension(:), allocatable, target :: dr_tried
  real(kind=8), dimension(:), allocatable :: norm_dr, norm_dr_tried

  real(kind=8) :: INITSTEPSIZE 
  real(kind=8) :: LOCAL_CUTOFF
  real(kind=8) :: INCREMENT,BASIN_FACTOR 
  real(kind=8) :: FTHRESHOLD
  real(kind=8) :: EIGEN_THRESH 
  real(kind=8) :: EXITTHRESH

  character(len=20) :: TYPE_EVENTS

  character(len=40) :: SEARCH_STRATEGY 


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

     evalf_number = 0                 ! Initialization of line_counter.

     nat = 3 * NATOMS                         

     ! We copy the reference half box and scaling to the working ones.
     ! We start from the reference configuration.
     scala = scalaref
     box   = boxref
     pos   = posref  

  ! _______
     ! for dual_search:
     central_atom = 0                 ! default for central_atom 

  !________

     ! These subroutines modify the vector pos and generate a vector
     ! of length 1 indicating the direction of the random
     ! displacement.  There is not case default.

     if ( eventtype == "GUESS_DIRECTION" ) then
        call guess_direction ( )  
     else
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

  end if

  ! Now, activate per se.
  call saddle_converge( ret, saddle_energy )
  
  ! If the activation did not converge, for whatever reason, we
  ! restart the routine and do not accept the new position.
  call end_report( success, ret, saddle_energy )


END SUBROUTINE find_saddle


!> ART global_move
!!   The initial random direction is taken from the full 3N-dimensional space
subroutine global_move( )

  use defs
  use random
  use saddles
  implicit none

  !Local variables
  integer :: i, j, number_of_min, s, number_of_dr, k, line_counter, success_count, rand
  real(kind=8) :: dr2
  real(kind=8) :: cos_theta
  real(kind=8) :: ran3
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8), dimension(:), pointer :: dx_tried, dy_tried, dz_tried
  real(kind=8) :: current_min(natoms,3)
  real(kind=8), allocatable :: read_min(:,:,:)
  real(kind=8) :: each_read_min(natoms,3)
  real(kind=8), allocatable :: sad_list(:)
  real(kind=8), allocatable :: dr_list(:)
  real(kind=8), allocatable :: dr_transformed_list(:,:,:)
  real(kind=8), allocatable :: read_dr(:,:,:)
  real(kind=8) :: each_read_dr(natoms,3)
  real(kind=8) :: each_dr_transformed(natoms,3)
  real(kind=8) :: r
  logical :: align_well
  character(len=20) :: keyword 

  allocate(dr(3*natoms)) 
  allocate(dr_tried(3*natoms))
  allocate(norm_dr(3*natoms)) 
  allocate(norm_dr_tried(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  

  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)
  
  dx_tried => dr_tried(1:NATOMS)
  dy_tried => dr_tried(NATOMS+1:2*NATOMS)
  dz_tried => dr_tried(2*NATOMS+1:3*NATOMS)

  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0

!  !Generate a random displacement.

selectcase ( search_strategy )

   case ('0')
        open(VLOG, file=VECLOG, action = 'write', position = 'append')
        write(VLOG,*) "Displacement vector"
        do i = 1, natoms, 1
        ! if ( constr(i) == 0 ) then
            do
               dx(i) = 0.5d0 - ran3()
               dy(i) = 0.5d0 - ran3()
               dz(i) = 0.5d0 - ran3()                             
               ! Ensures that the random displacement is isotropic
               dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
               if ( dr2 < 0.25d0 ) exit 
            end do
            write(VLOG,'((1X,a),3(2x,f16.8))') typat(i), dx(i), dy(i), dz(i)
            natom_displaced = natom_displaced + 1
            atom_displaced(i) = 1
        ! end if
        end do
        close(VLOG)

   case ('1') ! "follow" strategy


      open(VREAD, file = VECREAD, status = 'old', action = 'read')

        number_of_min =0
        s = 0

        do 
                read(VREAD,*,end=50) keyword
 

                if (keyword .eq. 'min') then

                        number_of_min = number_of_min + 1

                endif

                if (keyword .eq. 'sad') then

                        s = s + 1

                endif

        enddo


        50 rewind(VREAD)


        allocate(read_dr(s,natoms,3))
        allocate(read_min(number_of_min,natoms,3))
        allocate(sad_list(s))
        allocate(dr_list(s))
        

        number_of_min = 0
        line_counter = 0
        s= 0

        do
                read(VREAD,*,end=100) keyword

                line_counter = line_counter + 1

                if (keyword .eq. 'min') then

                        number_of_min = number_of_min + 1

                        do i = 1,natoms

                                read(VREAD,*) typat(i), read_min(number_of_min,i,1), read_min(number_of_min,i,2), read_min(number_of_min,i,3)
                                line_counter = line_counter + 1

                        enddo

                endif


                if (keyword .eq. 'sad') then

                        s = s + 1

                        sad_list(s) = line_counter

                endif


        enddo

        100 rewind(VREAD)


        do i=1,s
        dr_list(i) = sad_list(i)-(natoms+1)
        enddo
        
        number_of_dr = 0
        line_counter = 0
        do

                read(VREAD,*,end=200)

                line_counter = line_counter + 1

                do k =1,s
                        if (line_counter .eq. dr_list(k)) then

                        number_of_dr = number_of_dr + 1

                        do i = 1,natoms

                                read(VREAD,*) typat(i), read_dr(number_of_dr,i,1), read_dr(number_of_dr,i,2), read_dr(number_of_dr,i,3)
                                line_counter = line_counter + 1

                        enddo
                
                        endif

               enddo


        enddo

        200 rewind(VREAD)

        allocate(dr_transformed_list(number_of_dr,natoms,3))

        do i = 1,natoms

                current_min(i,1) = x(i)
                current_min(i,2) = y(i)
                current_min(i,3) = z(i)

        enddo
        

        success_count = 0

        align_well = .false.


        do j =1,number_of_dr

                do i = 1,natoms

                        each_read_min(i,1:3) = read_min(j,i,1:3)

                        each_read_dr(i,1:3) = read_dr(j,i,1:3)

                enddo

                call align(each_read_min, current_min, align_well, each_read_dr, each_dr_transformed)

                if(align_well) then
                        
                        success_count = success_count + 1

                        do i =1,natoms

                                dr_transformed_list(j,i,1:3) = each_dr_transformed(i,1:3)

                        enddo

                endif

        enddo


        call random_seed()
        call random_number(r)
        
        rand = ceiling(success_count*r) 

        do i = 1,natoms

                dx(i) =  dr_transformed_list(rand,i,1)
                dy(i) =  dr_transformed_list(rand,i,2)
                dz(i) =  dr_transformed_list(rand,i,3)

        enddo

        open(VLOG, file = VECLOG, status = 'unknown', action = 'write', position = 'append')
                
                write(VLOG,*) "Displacement vector"

                do i =1,NATOMS

                        write(VLOG,'(1X,a,3(2x,f16.8))') typat(i), dx(i), dy(i), dz(i)
        
                enddo

        close(VLOG)

        

   case ('2') ! "avoid" strategy

        open(VREAD, file=VECREAD, action='read', status = 'old')
        do
           do i = 1, natoms, 1
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
         
           read(VREAD,*) keyword

           if (keyword .eq. 'Displacement') then

           do i = 1, NATOMS
           read(VREAD,*) typat(i), dx_tried(i), dy_tried(i), dz_tried(i)
           enddo

           norm_dr = dr/sqrt(dot_product(dr,dr))
           norm_dr_tried = dr_tried/sqrt(dot_product(dr_tried, dr_tried))
           cos_theta = dot_product(norm_dr, norm_dr_tried)

           if (cos_theta .gt. 0.8) continue 

           write(*,*) "This is the angle between the two vectors", cos_theta

           endif
           
           if (keyword .eq. "SADDLE=") exit

           enddo
        
        close(VLOG)

       
        open(VLOG,file=VECLOG, action='write', position = 'append', status = 'unknown')
        
        write(VLOG,*) "Displacement vector"

        do i =1,NATOMS

        write(VLOG,'(1X,a,3(2x,f16.8))') typat(i), dx(i), dy(i), dz(i)
        
        enddo

        close(VLOG)

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
  real(kind=8), dimension(3) :: boxl, invbox
  real(kind=8), dimension(:), pointer  :: dx, dy, dz

                                      ! Select an atom at random.
  if ( preferred_atom < 0 ) then
     ! Only between totally free atoms 
     do 
     that = int( NATOMS * ran3() + 1 ) 
       ! if ( constr(that) == 0 ) exit 
     end do
  else 
     that = preferred_atom
     !do i = 1,NATOMS
      !  write(*,*) 'CONSTR: ', i, constr(i)
     !end do
     !if ( constr(that) .ne. 0 ) then  ! die !
      !  write(*,*) ' ERROR: choosen atom is blocked in the geometry file '
       ! call end_art()                          
     !end if 
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
  write(*,*) 'BART: That atom = ', that


  boxl = box*scala                    ! without periodic boundary conditions box
                                      ! is zero so we set invbox to 1 
  do i = 1, 3
     if ( boxl(i) .lt. 1.d-08 ) then 
        invbox(i) = 1.0d0
     else
        invbox(i) = 1.0d0 / boxl(i)
     end if
  end do
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
    ! if ( constr(j) == 0 ) then
        xij = x(j) - xi - boxl(1) * nint((x(j)-xi) * invbox(1))
        yij = y(j) - yi - boxl(2) * nint((y(j)-yi) * invbox(2))
        zij = z(j) - zi - boxl(3) * nint((z(j)-zi) * invbox(3))
        
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
   !  if ( constr(j) == 0 ) then
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
   !  end if
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
  real(kind=8), dimension(3) :: boxl, invbox
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
    ! if ( constr(that) == 0 ) exit 
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
  write(*,*) 'BART: That atom = ', that
  
  
  boxl = box*scala                    ! without periodic boundary conditions box
                                      ! is zero so we set invbox to 1 
  do i = 1, 3
     if ( boxl(i) .lt. 1.d-08 ) then 
        invbox(i) = 1.0d0
     else
        invbox(i) = 1.0d0 / boxl(i)
     end if
  end do
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
   !  if ( constr(j) == 0 ) then
        xij = x(j) - xi - boxl(1) * nint((x(j)-xi) * invbox(1))
        yij = y(j) - yi - boxl(2) * nint((y(j)-yi) * invbox(2))
        zij = z(j) - zi - boxl(3) * nint((z(j)-zi) * invbox(3))
        
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
    ! end if
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
    ! if ( constr(i) == 0 ) then 
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
   !  end if
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
  
  ! And normalize the total random displacement effected to the value
  ! desired.
  norm = 0.0d0
  norm = dot_product(dr, dr)
  
  ! This renormalizes in angstroems to a displacement INITSTEPSIZE.
  norm = step / sqrt(norm)
  xnorm = norm 
  ynorm = norm
  znorm = norm
  
  ! The displacement is now in box units. 
  dx = dx * xnorm 
  dy = dy * ynorm
  dz = dz * znorm

  ! Redistribute the information.
  nat = 3*NATOMS

  ! Update the position using this random displacement
  pos = pos + dr
  ! DEBUG Bhupinder
  write(*,*) 'Updated position using random displacement' , pos  
  
  ! Now, we normalize dr to get the initial_direction (note that this
  ! had to be done after the transfer into box units.
  initial_direction = dr 
  write(*,*) 'This is dr:'
  write(*,*) dr
  norm = dot_product( initial_direction, initial_direction ) 
  norm = 1.0d0 / sqrt(norm)
  initial_direction  = initial_direction * norm
  write(*,*)'Initial direction is:'
  write(*,*) initial_direction

  write(*,*) 'BART: Number of displaced atoms initially: ',natom_displaced

  deallocate(atom_displaced)
  deallocate(dr)
  deallocate(dr_tried)
  deallocate(norm_dr)
  deallocate(norm_dr_tried)

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
  real(kind=8), dimension(3) :: boxl, invbox

  !_______________________

  lcutoff2 = size_system*size_system
  boxl = box*scala                    ! without periodic boundary conditions box
                                      ! is zero so we set invbox to 1 
  do i = 1, 3
     if ( boxl(i) .lt. 1.d-08 ) then 
        invbox(i) = 1.0d0
     else
        invbox(i) = 1.0d0 / boxl(i)
     end if
  end do

  i = central_atom 

  xi = x(i)
  yi = y(i)
  zi = z(i)

  do j = 1, NATOMS
     if ( j == i ) cycle

     xij = x(j) - xi - boxl(1) * nint((x(j)-xi) * invbox(1))
     yij = y(j) - yi - boxl(2) * nint((y(j)-yi) * invbox(2))
     zij = z(j) - zi - boxl(3) * nint((z(j)-zi) * invbox(3))

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
  real(kind=8), dimension(3) :: boxl, invbox
  real(kind=8), dimension(:), pointer  :: dx, dy, dz

  integer, dimension(natoms) :: in_list

  !_______________________

  lcutoff2 = coord_length*coord_length
  boxl = box*scala                    ! without periodic boundary conditions box
                                      ! is zero so we set invbox to 1 
  do i = 1, 3
     if ( boxl(i) .lt. 1.d-08 ) then 
        invbox(i) = 1.0d0
     else
        invbox(i) = 1.0d0 / boxl(i)
     end if
  end do

  in_list = 0                       ! Initial vectorial assignment
 
  do i = 1, NATOMS
     
     numnei = 0
     xi = x(i)
     yi = y(i)
     zi = z(i)

     do j = 1, NATOMS
        if ( j == i ) cycle

        xij = x(j) - xi - boxl(1) * nint((x(j)-xi) * invbox(1))
        yij = y(j) - yi - boxl(2) * nint((y(j)-yi) * invbox(2))
        zij = z(j) - zi - boxl(3) * nint((z(j)-zi) * invbox(3))

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
     if ( in_list(that) ==1) then ! .and. constr(that) == 0 ) then
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
     write(*,*)  'BART ERROR: There is no atom with lower coord than ', coord_number
     stop
  end if

  call symmetry_break( )              !Breaks the initial symmetry.

  open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
  write(FLOG,*) ' '
  write(FLOG,'(1X,A34,I17)') ' - That atom                    : ', that
  close(FLOG)
  write(*,*) 'BART: That atom = ', that

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
   !  if ( constr(j) == 0 ) then
        xij = x(j) - xi - boxl(1) * nint((x(j)-xi) * invbox(1))
        yij = y(j) - yi - boxl(2) * nint((y(j)-yi) * invbox(2))
        zij = z(j) - zi - boxl(3) * nint((z(j)-zi) * invbox(3))

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
    ! end if
  end do

  call center_and_norm ( INITSTEPSIZE )


END SUBROUTINE coord_based_move 


subroutine guess_direction ( )

  use defs
  use saddles

  !Local variables
  integer :: i, ierror
  real(kind=8) :: ran3
  real(kind=8), dimension(3) :: boxl, invbox
  real(kind=8), dimension(:), pointer :: xa, ya, za
  real(kind=8), dimension(:), pointer :: xb, yb, zb
  real(kind=8), dimension(:), pointer :: dx, dy, dz
  real(kind=8), dimension(:), allocatable, target :: posa
  real(kind=8), dimension(:), allocatable, target :: posb
  real(kind=8)  :: norm

  open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
  write(FLOG,*) ' '
  write(FLOG,'(1X,A47,A20)') ' -following a given initial direction in file: ', GUESSFILE
  write(FLOG,'(1X,A20,F12.6)') '  Noise amplitude : ', guess_noise 
  close(FLOG)

  boxl = box*scala                    ! without periodic boundary conditions box
                                      ! is zero so we set invbox to 1 
  do i = 1, 3
     if ( boxl(i) .lt. 1.d-08 ) then 
        invbox(i) = 1.0d0
     else
        invbox(i) = 1.0d0 / boxl(i)
     end if
  end do

  allocate(posa(3*natoms))
  allocate(posb(3*natoms))
  allocate(dr(3*natoms))
  posa = 0.0d0
  posb = 0.0d0
  dr   = 0.0d0

  posa(:) = pos(:)
  posb(:) = g_pos(:)
                                      ! maybe this is no neccesary.
  call center( posa, 3*natoms )
  call center( posb, 3*natoms )
                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  xa => posa(1:NATOMS)
  ya => posa(NATOMS+1:2*NATOMS)
  za => posa(2*NATOMS+1:3*NATOMS)

  xb => posb(1:NATOMS)
  yb => posb(NATOMS+1:2*NATOMS)
  zb => posb(2*NATOMS+1:3*NATOMS)

  ! a small noise helps if we fail in the attempt
  do i = 1, NATOMS
    ! if ( constr(i) == 0 ) then
        dx(i) = xb(i) - xa(i) - boxl(1) * nint((xb(i)-xa(i)) * invbox(1)) + guess_noise*(0.5d0-ran3()) 
        dy(i) = yb(i) - ya(i) - boxl(2) * nint((yb(i)-ya(i)) * invbox(2)) + guess_noise*(0.5d0-ran3())
        dz(i) = zb(i) - za(i) - boxl(3) * nint((zb(i)-za(i)) * invbox(3)) + guess_noise*(0.5d0-ran3())
   !  end if
  end do
 
  norm = dot_product( dr, dr ) 
  norm = 1.0d0 / sqrt(norm) 
  dr = dr*norm

  call center( dr, 3*natoms )
  norm = dot_product( dr, dr )     ! is it really necesary normalize again?
  norm = 1.0d0 / sqrt(norm) 
  dr = dr*norm
  initial_direction = dr     

  pos = pos + INITSTEPSIZE*initial_direction

  deallocate(dr)
  deallocate(dr_tried)
  deallocate(norm_dr)
  deallocate(norm_dr_tried)
  deallocate(posa)
  deallocate(posb)
  write(*,*) 'BART: Number of displaced atoms initially: ', natoms 

END SUBROUTINE guess_direction 

SUBROUTINE align (each_read_min, current_min, align_well, each_read_dr, each_dr_transformed)


        !! This subroutine calculates a transformation matrix that minimizes the RMSD between two sets of coordinates.
        !The first step of this transformation involves a translation. This is simply done by re-centering the molecules so that their
        !centroids coincide with the origin of the coordinate system.
        !Next, a covariance matrix is calculated using the re-centered molecules.
        !Then, the singular value decomposition (SVD) of this covariance matrix is calculated, which yields the factors of this
        !covariance, viz U, S and V (Here U and V are the left- and right- singular vectors of the covariance matrix.)
        ! According to the Kabsch algorithm, the product of V and U(transpose) is the optimal rotation matrix that minimizes the RMSD between
        ! the two sets of vectors.
        ! Here cov is the covariance matrix, U and VT are the left- and
        !right- singular vectors of the cov, (VT is actually the transpose of V),
        !S is one of the factors of the cov. R (V*U_transpose) is the rotation matrix
        ! Work is one of the arguments required by the LAPACK library subroutine that calculates the SVD
        !Documentation at: 
        !http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html#ga84fdf22a62b12ff364621e4713ce02f2

        use defs

        implicit none
        
        real(kind=8), intent(in) :: each_read_min(natoms,3)
        real(kind=8), intent(in) :: current_min(natoms,3)
        real(kind=8), intent(in) :: each_read_dr(natoms,3)

        logical, intent(inout) :: align_well

        real(kind=8), intent(out) :: each_dr_transformed(natoms,3)

        real(kind=8) :: each_read_min_transformed(natoms,3)

        real(kind=8) :: each_read_min_moved(natoms,3)
        real(kind=8) :: current_min_moved(natoms,3)

        real(kind=8) :: each_read_min_nat4(natoms,4)
        real(kind=8) :: each_read_dr_nat4(natoms,4)

        real(kind=8) :: each_read_min_transformed_nat4(natoms,4)
        real(kind=8) :: each_read_dr_transformed_nat4(natoms,4)
        real(kind=8) :: transformation_vec(4,4)

        real(kind=8) :: deviation_each_atom_before(natoms), deviation_each_atom_after(natoms)
        real(kind=8) :: sum_deviations_before, rmsd_before, sum_deviations_after, rmsd_after
        real(kind=8) :: sum_read_x, sum_read_y, sum_read_z, cx_read, cy_read, cz_read
        real(kind=8) :: sum_current_x, sum_current_y, sum_current_z, cx_current, cy_current, cz_current
         
        real(kind=8) :: cov(3,3), U(3,3), S(3), VT(3,3), R(3,3), work(20)
        real(kind=8) :: U_transpose(3,3), V(3,3), det_U_transpose, det_V
       
        integer :: i, j, info
        

        sum_current_x = 0.0
        sum_current_y = 0.0
        sum_current_z = 0.0

        sum_read_x = 0.0
        sum_read_y = 0.0
        sum_read_z = 0.0
        
        sum_deviations_before = 0.0
        sum_deviations_after = 0.0


        do i=1,natoms

                deviation_each_atom_before(i) = (each_read_min(i,1)-current_min(i,1))**2 + (each_read_min(i,2)-current_min(i,2))**2 + (each_read_min(i,3)-current_min(i,3))**2
       
        enddo

        do i=1,natoms

                sum_deviations_before = sum_deviations_before + (each_read_min(i,1)-current_min(i,1))**2 + (each_read_min(i,2)-current_min(i,2))**2 + (each_read_min(i,3)-current_min(i,3))**2
       
        enddo


        rmsd_before = sqrt(sum_deviations_before/natoms) !Calculating RMSD before structural alignment


        write(*,*) "This is rmsd before alignment: ", rmsd_before


        do i=1,natoms

                sum_read_x = sum_read_x + each_read_min(i,1)
                sum_read_y = sum_read_y + each_read_min(i,2)
                sum_read_z = sum_read_z + each_read_min(i,3)

                sum_current_x = sum_current_x + current_min(i,1)
                sum_current_y = sum_current_y + current_min(i,2)
                sum_current_z = sum_current_z + current_min(i,3)

        enddo

       
        !Calculating the centroids for all the points


        cx_read = sum_read_x/natoms
        cy_read = sum_read_y/natoms
        cz_read = sum_read_z/natoms


        cx_current = sum_current_x/natoms
        cy_current = sum_current_y/natoms
        cz_current = sum_current_z/natoms


        !Recentering the structures so that their centroids coincide with the origin of the coordinate system.


        do i = 1,natoms

                each_read_min_moved(i,1) = each_read_min(i,1) - cx_read
                each_read_min_moved(i,2) = each_read_min(i,2) - cy_read
                each_read_min_moved(i,3) = each_read_min(i,3) - cz_read 

                current_min_moved(i,1) = current_min(i,1) - cx_current
                current_min_moved(i,2) = current_min(i,2) - cy_current
                current_min_moved(i,3) = current_min(i,3) - cz_current

        enddo
       
        
        !Calculating the covariance matrix         


        cov = matmul(transpose(each_read_min_moved),current_min_moved)

        !Invoking the LAPACK library subroutine "dgesvd" that calculates the SVD of the covariance matrix

        call dgesvd('S','S',3,3,cov,3,S,U,3,VT,3,work,20,info)

        U_transpose = transpose(U)

        V = transpose(VT)

        det_U_transpose = U_transpose(1,1)*(U_transpose(2,2)*U_transpose(3,3) - U_transpose(2,3)*U_transpose(3,2)) - U_transpose(1,2)*(U_transpose(2,1)*U_transpose(3,3) - U_transpose(2,3)*U_transpose(3,1)) + U_transpose(3,1)*(U_transpose(2,1)*U_transpose(3,2) - U_transpose(2,2)*U_transpose(3,1))

        det_V = V(1,1)*(V(2,2)*V(3,3) - V(2,3)*V(3,2)) - V(1,2)*(V(2,1)*V(3,3) - V(2,3)*V(3,1)) + V(3,1)*(V(2,1)*V(3,2) - V(2,2)*V(3,1))
        
        if(det_V * det_U_transpose .LT. 0) then       ! A small correction to ensure a right-handed coordinate system

                V(:,3) = -V(:,3)              

        endif

        R = matmul(V,U_transpose)  !This is the optimal rotation matrix 

        do i =1,3

                transformation_vec(i,1:3) = R(i,1:3)

        enddo

        transformation_vec(4,1:4) = 0.0d0

        transformation_vec(1,4) = cx_read
        transformation_vec(2,4) = cy_read
        transformation_vec(3,4) = cz_read


        do i =1,natoms

                each_read_min_nat4(i,1:3) = each_read_min(i,1:3)
                each_read_min_nat4(i,4) = 0.0d0

                each_read_dr_nat4(i,1:3) = each_read_dr(i,1:3)
                each_read_dr_nat4(i,4) = 0.0d0

        enddo

        each_read_min_transformed_nat4 = transpose(matmul(transformation_vec,transpose(each_read_min_nat4)))

        do i = 1,natoms

                each_read_min_transformed(i,1:3) = each_read_min_transformed_nat4(i,1:3)

        enddo

        do i=1,natoms

                deviation_each_atom_after=(each_read_min_transformed(i,1)-current_min(i,1))**2 + (each_read_min_transformed(i,2)-current_min(i,2))**2 + (each_read_min_transformed(i,3)-current_min(i,3))**2 
        enddo
        
        do i=1,natoms

                sum_deviations_after = sum_deviations_after + (each_read_min_transformed(i,1)-current_min(i,1))**2 + (each_read_min_transformed(i,2)-current_min(i,2))**2 + (each_read_min_transformed(i,3)-current_min(i,3))**2 

        enddo


        rmsd_after = sqrt(sum_deviations_after/natoms) !Calculating RMSD after structural alignment

        write(*,*) "This is rmsd after alignment: ", rmsd_after

        do j = 1,natoms

                if (rmsd_after .LT. 0.1 .and. deviation_each_atom_after(j) .LT. 0.1) then

                        align_well = .true.

                        each_read_dr_transformed_nat4 = transpose(matmul(transformation_vec,transpose(each_read_dr_nat4)))

                        do i = 1,natoms

                                each_dr_transformed(i,1:3) = each_read_dr_transformed_nat4(i,1:3)

                        enddo

                endif

        enddo
        
END SUBROUTINE align
