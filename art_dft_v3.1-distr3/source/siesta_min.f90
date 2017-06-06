! In this routine, we interface between ART and SIESTA
subroutine min_converge_dft(success)
  use defs

  integer :: i, idum, ret 
  logical :: success
  integer, parameter :: FSIESTA = 21
  real*8, parameter :: ZERO = 0.0d0
  real*8, dimension(natoms) :: xx, yy, zz
  character(len=20) :: SIESTA   = 'art2siesta'
  character(len=20) :: SIESTAFORCE = 'siesta2art'
  character(len=70) :: line
  logical :: read_done,read_final
  real(kind=8),dimension(3) :: boxl
  real(8) :: toto

  boxl = box * scala

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'start minimization'
  close(flog)

  ! We first write to a file the format requested by SIESTA
  open(unit=FSIESTA,file=SIESTA,status='replace',action='write',iostat=ierror)

  write(FSIESTA,"('MD.TypeOfRun         cg')")
  write(FSIESTA,"('MD.NumCGsteps       150')")
  write(FSIESTA,"('MD.MaxCGDispl        0.1  Ang')")
  write(FSIESTA,"('MD.MaxForceTol      0.002 eV/Ang')")
  write(FSIESTA,*)
  write(FSIESTA,"('%block Atomic_Coordinates_and_Atomic_Species')")
  write(FSIESTA,"(f14.8, f14.8, f14.8, i4)")  ( x(i), y(i), z(i), typat(i), i=1, NATOMS )
  write(FSIESTA,"('%endblock Atomic_Coordinates_and_Atomic_Species')")
  write(FSIESTA,"('%include basicinfo.fdf')")
  close(FSIESTA)

  ! We now call siesta do to the minimization
  call  system('./execute.sh')
  
  do i=1, 10000
    toto = dexp ( i * 0.001d0)
  end do
  ! We must now read the energy and positions from siesta's output file
  open(unit=FSIESTA,file=SIESTAFORCE,status='old',action='read',iostat=ierror)

  ! We first read the total energy
  read_done = .false.
  do 
    read(FSIESTA,"(A40)") line
    if ( line  == "siesta: Final energy (eV):" ) then
!      do i = 1, 7
      do i = 1, 8

!        read(FSIESTA,"(A40)") line
        read(FSIESTA,"(A38)") line
      end do
      read(FSIESTA,"(A24,f14.6)") line, total_energy
      
      read_done = .true.
    endif
    if(read_done) exit
  end do
  close(FSIESTA)

  write(*,*) 'the last line is : ', line
  write(*,*) 'total energy is : ', total_energy  

  ! We now do the force
  open(unit=FSIESTA,file=SIESTAFORCE,status='old',action='read',iostat=ierror)
  read_done = .false.
  read_final= .false.
  do 
    read(FSIESTA,"(A70)") line
    if ( line(1:8)  == "outcoor:" ) then
      do i = 1, NATOMS
        read(FSIESTA,*) x(i),y(i),z(i)
      end do
      read_done = .true.
    endif
    if ( line  == "siesta: Final energy (eV):" ) read_final = .true.
    if(read_done .and. read_final ) exit
  end do

  close(FSIESTA)
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'End minimization'
  close(flog)

  success = .true.
  return
end subroutine
