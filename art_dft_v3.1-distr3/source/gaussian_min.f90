! In this routine, we interface between ART and GAUSSIAN
subroutine min_converge_gau(success)
  use defs

  integer :: i, ret

  logical :: success
  integer, parameter :: FGAUSS = 21
  real*8, parameter :: ZERO = 0.0d0
  real*8, dimension(natoms) :: xx, yy, zz
  character(len=20) :: GAUSS   = 'art2gaussian'
  character(len=20) :: GAUSSFORCE = 'gaussian2art'
  character(len=70) :: line
  character(len=10) :: string_natoms
  logical :: read_coordinates_done,read_energy_done
  real(kind=8),dimension(3) :: boxl

  boxl = box * scala

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'start minimization'
  close(flog)

  ! We first write to a file the format requested by GAUSS
  open(unit=FGAUSS,file=GAUSS,status='replace',action='write',iostat=ierror)

  ! Prepares gaussian input file coordinates
  write(FGAUSS,"(i4, f14.8, f14.8, f14.8)")  (typat(i),  x(i), y(i), z(i), i=1, NATOMS )
  write(FGAUSS,*)
  close(FGAUSS)

  !converting the number of atoms to a string value
  write(string_natoms, '(i10)' )  NATOMS

  ! We now call Gaussian do to the minimization
  ! Bash parameters: natoms=$1, optimization=$2
  call system('sh execute_gaussian.sh ' // string_natoms // ' ' // 'opt')

  ! ! We first read the total energy
  ! read_done = .false.
  ! do 
  !   read(FGAUSS,"(A40)") line
  !   if ( line  == "gaussi: Final energy (eV):" ) then
  !     ! Debug
  !     write (*,*) 'test1 after energy'
  !     do i = 1, 7
  !       read(FGAUSS,"(A40)") line
  !     end do
  !     read(FGAUSS,"(A24,f14.6)") line, total_energy
      
  !     read_done = .true.
  !   endif
  !   if(read_done) exit
  ! end do
  ! close(FGAUSS)

  ! write(*,*) 'the last line is : ', line
  ! write(*,*) 'total energy is : ', total_energy  
  
  !NEW SECTION
  
  ! We must now read the energy and positions from gaussian's output file
  open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)

  read_energy_done = .false.
  read_coordinates_done= .false.

  do 
    !Gets the Gaussian output coordinates
    read(FGAUSS,"(A70)") line

    ! Checks for IO errors or EOF marker
    ! if (ierror .eq. end_of_file) then
    !   if (.not.success) then
    !     write(*,*) 'IO error in gaussian_min.f90'
    !     stop
    !   endif
    !   success = .false.
    !   exit
    ! endif

    if ( line  == "outcoor:" ) then
      do i = 1, NATOMS
        ! read(FGAUSS,"(f15.6,f13.6,f11.6)")x(i),y(i),z(i)
        read (FGAUSS,*) x(i), y(i), z(i)
        write (*,*) "pos_min", x(i),y(i),z(i)
      end do
      read_coordinates_done = .true.
    endif
    if(read_coordinates_done ) exit
  end do

  close(FGAUSS)

  open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)
  do 
    read(FGAUSS,"(A40)") line
    ! Checks for IO errors or EOF marker
    ! if (ierror .eq. end_of_file) then
    !   if (.not.success) then
    !     write(*,*) 'IO error in gaussian_min.f90'
    !     stop
    !   endif
    !   success = .false.
    !   exit
    ! endif

    !Gets the final energy
    if ( line  == "energy:" ) then
      read(FGAUSS,*) energy
      read_energy_done = .true.    
    endif
    if(read_energy_done) exit
  end do
  close(FGAUSS)

  ! ! We now do the force
  ! open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)
  ! read_done = .false.
  ! read_energy_done= .false.
  ! do 
  !   read(FGAUSS,"(A70)") line
  !   !TODO change from hardcoded
  !   if ( line(1:8)  == "outcoor:" ) then
  !     do i = 1, NATOMS
  !       read(FGAUSS,"(i7,i12,i12,f15.6,f13.6,f11.6)")x(i),y(i),z(i)
  !       write(*,*) "pos_min", idum,idum1,idum2,x(i),y(i),z(i)
  !     end do
  !     read_done = .true.
  !   endif
  !   if ( line  == "gaussi: Final energy (eV):" ) read_energy_done = .true.
  !   if(read_done .and. read_energy_done ) exit
  ! end do

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'End minimization'
  close(flog)

  success = .true.
  return
 write (*,*) 'Running gaussian_min'
end subroutine
