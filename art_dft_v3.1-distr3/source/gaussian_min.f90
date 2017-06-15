! In this routine, we interface between ART and GAUSSIAN
subroutine min_converge_gau(success)
  use defs

  integer :: i, idum, ret

  integer :: idum1, idum2
  logical :: success
  integer, parameter :: FGAUSS = 21
  real*8, parameter :: ZERO = 0.0d0
  real*8, dimension(natoms) :: xx, yy, zz
  character(len=20) :: GAUSS   = 'art2gaussian.inp'
  character(len=20) :: GAUSSFORCE = 'log'
  character(len=70) :: line
  character(len=10) :: string_natoms
  logical :: read_done,read_final
  real(kind=8),dimension(3) :: boxl
  real(8) :: toto

  boxl = box * scala

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'start minimization'
  close(flog)

  ! We first write to a file the format requested by GAUSS
  open(unit=FGAUSS,file=GAUSS,status='replace',action='write',iostat=ierror)

 ! character(len=15) :: GAU_mem
 ! character(len=15) :: GAU_nproc
 ! character(len=50) :: GAU_desc
 ! character(len=30) :: GAU_title
 ! integer           :: GAU_charge
 ! integer           :: GAU_multip



 !converting multiplicity to a string value
  write(string_natoms, '(i10)' )  GAU


! TODO get this hardcoded data from parameters
     !write(FGAUSS,"('%chk=temp.chk')")
     !write(FGAUSS,"('#rhf/3-21g nosymm opt')")
     !write(FGAUSS,*)
     !write(FGAUSS,"('message')")
     !write(FGAUSS,*)
     !write(FGAUSS,"('0 1')")
     write(FGAUSS,"(i4, f14.8, f14.8, f14.8)")  (typat(i),  x(i), y(i), z(i), i=1, NATOMS )
     write(FGAUSS,*)
     close(FGAUSS)

  !converting the number of atoms to a string value
  write(string_natoms, '(i10)' )  NATOMS

  ! We now call Gaussian do to the minimization
  !call system('sh execute_gaussian.sh ' // string_natoms // ' ' // 'opt ' // '%mem=8000MB')
  call system('sh execute_gaussian.sh ' // string_natoms // ' ' // 'opt')
  
  do i=1, 10000
    toto = dexp ( i * 0.001d0)
  end do
  ! We must now read the energy and positions from gaussian's output file
  open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)

  ! We first read the total energy
  read_done = .false.
  do 
    read(FGAUSS,"(A40)") line
    if ( line  == "gaussi: Final energy (eV):" ) then
      ! Debug
      write (*,*) 'test1 after energy'
      do i = 1, 7
        read(FGAUSS,"(A40)") line
      end do
      read(FGAUSS,"(A24,f14.6)") line, total_energy
      
      read_done = .true.
    endif
    if(read_done) exit
  end do
  close(FGAUSS)

  write(*,*) 'the last line is : ', line
  write(*,*) 'total energy is : ', total_energy  

  ! We now do the force
  open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)
  read_done = .false.
  read_final= .false.
  do 
    read(FGAUSS,"(A70)") line
    if ( line(1:8)  == "outcoor:" ) then
      do i = 1, NATOMS
!        read(FGAUSS,*) x(i),y(i),z(i)
! bharat starts
        read(FGAUSS,"(i7,i12,i12,f15.6,f13.6,f11.6)")idum,idum1,idum2,x(i),y(i),z(i)
    write(*,*) "pos_min", idum,idum1,idum2,x(i),y(i),z(i)
! bharat ends 
      end do
      read_done = .true.
    endif
    if ( line  == "gaussi: Final energy (eV):" ) read_final = .true.
    if(read_done .and. read_final ) exit
  end do

  close(FGAUSS)
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'End minimization'
  close(flog)

  success = .true.
  return
 write (*,*) 'Running gaussian_min'
end subroutine
