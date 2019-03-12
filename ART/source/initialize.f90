!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! Modified by:
!! -EM 2010, see ~/AUTHORS
!! -Laurent Karim Beland, UdeM, 2011. For working with QM/MM !!
!! -EM 2011, see ~/AUTHORS
!!
!> ART initialize
!!   Initialization of art method
!!   It relaxes it into a local minimum without volume optimization
!!   The initial configuration is a "reference configuration", it will be the reference 
!!   configuration until a new event is accepted. 
!! 
subroutine initialize()

  use defs
  use lanczos_defs,  only : LANCZOS_MIN 
  use saddles,       only : g_pos, GUESSFILE
  implicit none

  !Local variables
  integer :: i, ierror
  character(len=20) :: dummy, fname
  character(len=4)  :: scounter
  logical           :: flag, success

  integer                    :: nat_test
  character(len=3), pointer  :: typ_a(:)    ! Atomic type
  real(kind=8), pointer      :: pos_a(:)    ! Working positions of the atoms
  real(kind=8), dimension(:), pointer :: xa, ya, za

  !_______________________
  
  ! Read the counter in order to continue the run where it stopped or for refine
  ! Format:
  ! Counter:     1000

  if (.not. restart) then 
     inquire( file = COUNTER, exist = flag )
     if ( flag ) then 
        open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
        read(FCOUNTER,'(A12,I6)') dummy, mincounter 
        close(FCOUNTER)
     else
        mincounter = 1000
     end if
  end if 

  ! we read the initial/reference configuration
  allocate(pos_a(vecsize))
  allocate(typ_a(NATOMS))
  xa => pos_a(1:NATOMS)
  ya => pos_a(NATOMS+1:2*NATOMS)
  za => pos_a(2*NATOMS+1:3*NATOMS)  
 
  open(unit=FREFCONFIG,file=REFCONFIG,status='old',action='read',iostat=ierror)
  read(FREFCONFIG,*) dummy, refcounter
  read(FREFCONFIG,*) dummy, ref_energy
  do i = 1, NATOMS
    read(FREFCONFIG,*) typ_a(i), xa(i), ya(i), za(i)
  end do
  close(FREFCONFIG)

  !assign the data from the atomic file
  if ( .not. restart ) then
     typat(:)   = typ_a(:)
     pos(:)     = pos_a(:)
     refcounter = mincounter
  else 
     posref(:) = pos_a(:)
  end if

     deallocate(typ_a)
     deallocate(pos_a)
  

     !atomic type copy in atom
     do i = 1, NATOMS
        Atom(i) = typat(i)
     end do

  ! write
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,'(1X,A34,I17)') ' - Mincounter                   : ', mincounter
  close(FLOG)

  ! We rescale the coordinates. For what ?? 

  call initialize_potential()         ! Initialize Potential (CORE) 


  !DEBUG starts Bhupinder
  write(*,*) 'DEBUG!!! At this point, the coordinates are rescaled by the initialize_potential subroutine'
  !DEBUG ends Bhupinder

  ! for output files
  call convert_to_chain( refcounter, 4, scounter )
  fname = FINAL // scounter
  conf_initial = fname
                                      ! If this is a new event we relax 
  If_ne: if ( new_event .and. (.not. restart) ) then 
     call min_converge( success )     ! Converge the configuration to a local minimum

     posref = pos                     ! New reference configuration.
     ref_energy = total_energy
 

     call write_refconfig( )       ! Write reference in REFCONFIG. 
     call store( fname )           ! Store the configuration into fname.
     open(unit=VLOG, file = VECLOG, status = 'unknown', action='write', position='append')
     write(VLOG,*) FINAL // ' ' // scounter
     do i=1,NATOMS
        write(VLOG,'(1x,a,3(2x,f16.8))') typat(i), x(i), y(i), z(i)
     end do
     close(VLOG)

     open( unit = FLOG, file = LOGFILE, status = 'unknown',&
          & action = 'write', position = 'append', iostat = ierror )
     write(*,*) 'ART: Configuration stored in file ',fname
     write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
     if ( .not. success ) then
        write(FLOG,'(1X,A)') "ERROR: Initial configurations is not a minimum"
        call end_art()  
     end if
     close(FLOG)

     mincounter = mincounter + 1
     
     ! if dual_search we dont do this check at the beginning. It is not
     ! well defined.
     if ( success .and. .not. dual_search ) then
        if ( LANCZOS_MIN .or. setup_initial ) call check_min( 'M' ) 
     end if

  end if If_ne

END SUBROUTINE initialize
