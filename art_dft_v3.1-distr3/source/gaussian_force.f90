! In this routine, is the  interface between ART and Gaussian

subroutine calcforce_gau(nat,typa,posa,boxl,forca,energy)
  use defs
  integer, intent(in) :: nat
  integer, intent(in), dimension(NATOMS) :: typa
  real*8, intent(in), dimension(VECSIZE),target :: posa
  real*8, dimension(VECSIZE), target:: posS
  real*8, intent(in) :: boxl
  real*8, intent(out), dimension(VECSIZE),target :: forca
  real*8, intent(out) :: energy
  real*8 ::  toto

  integer :: i, idum, ierror
  integer :: idum1, idum2
  integer :: end_of_file = -1
  integer, parameter :: FGAUSS = 21
  real*8, parameter :: ZERO = 0.0d0
  character(len=20) :: GAUSS   = 'art2gaussian.inp'
  character(len=20) :: GAUSSFORCE = 'log'
  character(len=40) :: line
  character(len=10) :: string_natoms
  logical :: read_done,read_final,read_doneF,read_doneC, success
  real(8), dimension(:), pointer :: xS, yS, zS
  real(8), dimension(:), pointer :: xa, ya, za
  real(8), dimension(:), pointer :: fax, fay, faz   ! Pointers for working force

  success = .true.

  ! We first set-up pointers for the x, y, z components in the position and forces
  xa => posa(1:NATOMS)
  ya => posa(NATOMS+1:2*NATOMS)
  za => posa(2*NATOMS+1:3*NATOMS)

  posS = posa


  xS => posS(1:NATOMS)
  yS => posS(NATOMS+1:2*NATOMS)
  zS => posS(2*NATOMS+1:3*NATOMS)

  fax => forca(1:NATOMS)
  fay => forca(NATOMS+1:2*NATOMS)
  faz => forca(2*NATOMS+1:3*NATOMS)
  
  ! Transform the coordinates system from box to angstroems

  do

     ! We first write to a file the format requested by Gaussian
     open(unit=FGAUSS,file=GAUSS,status='replace',action='write',iostat=ierror)

! Prepares gaussian input file
! TODO cleanup to take config parameters
!     write(FGAUSS,"('%chk=temp.chk')")
!     write(FGAUSS,"('#rhf/3-21g nosymm force')")
!     write(FGAUSS,*)
!     write(FGAUSS,"('name')")
!     write(FGAUSS,*)
!     write(FGAUSS,"('0 1')")
     write(FGAUSS,"(i4, f14.8, f14.8, f14.8)")  ( typa(i), xa(i), ya(i),za(i), i=1, NATOMS)
     write(FGAUSS,*)
     close(FGAUSS)


     !converting the number of atoms to a string value
     write(string_natoms, '(i10)' )  NATOMS

     !We now call Gaussian do to the minimization
     ! Bash parameters: natoms=$1, nproc=$2, mem=$3
     !call system('sh execute_gaussian.sh ' // string_natoms // ' ' // '%nproc=12 ' // '%mem=8000MB')
     call system('sh execute_gaussian.sh ' // string_natoms // ' ' // 'force')


     do i=1, 10000
        toto = dexp ( i * 0.001d0)
     end do


     ! We must now read the forces from Gaussian's output file
     open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)

     read_done = .false.
     do 
        read(FGAUSS,"(A40)") line
        if ( line  == "gaussi: Final energy (eV):" ) then
           do i = 1, 7
              read(FGAUSS,"(A40)") line
           end do
           read(FGAUSS,"(A24,f14.6)") line, energy
           read_done = .true.
        endif
        if(read_done) exit
     end do
     close(FGAUSS)
     
     ! We now do the force
     open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)
     read_doneF = .false.
     read_final = .false.
     read_doneC = .false.
     do 
        read(FGAUSS,"(A40)") line
        if (ierror .eq. end_of_file) then
           if (.not.success) then
              write(*,*) 'Error with Gaussian force'
              stop
           endif
           success = .false.
           exit
        endif
        if ( line  == "gaussi: Atomic forces (eV/Ang):" ) then
           do i = 1, NATOMS
! bharat              read(FGAUSS,"(i7,3f12.6)") idum,fax(i),fay(i),faz(i)
        read(FGAUSS,"(i7,i10,f21.9,f15.9,f15.9)") idum,idum1, fax(i),fay(i),faz(i)
!              read(FGAUSS,"(i7,i7,3f13.9)") idum, idum,fax(i)*51.4220629786602,fay(i)*51.4220629786602,faz(i)*51.4220629786602
!              read(FGAUSS,"(i8,i10,3f12.9)") idum,idum,fax(i),fay(i),faz(i)

! unit conversion Gaussian format Forces (Hartrees/Bohr) to Gaussian format forces forces (eV/Ang):
! hartree_to_ev=27.2113838668 and bohr_to_angstrom=0.52917721092
! changed sign for forces as suggested by Guillaume on Nov 15,2016 NOT WORKING but did not work so backed to same.
        fax(i) = fax(i)*51.4220629786602
        fay(i) = fay(i)*51.4220629786602
        faz(i) = faz(i)*51.4220629786602
!write (*,*) "Forces", idum,idum1,fax(i)*51.4220629786602,fay(i)*51.4220629786602,faz(i)*51.4220629786602
!write (*,*) "Forces",idum,idum1,fax(i),fay(i),faz(i)
!stop

           end do
           read_doneF = .true.
        endif
        
        if ( line(1:8)  == "outcoor:" ) then
           do i = 1, NATOMS
!              read(FGAUSS,*) xS(i),yS(i),zS(i)
! debug starts
! write (*,*) "position formatting enters"
        read(FGAUSS,"(i7,i12,i12,f15.6,f13.6,f11.6)") idum,idum1,idum2,xS(i),yS(i),zS(i)
! write (*,*) idum,idum1,idum2,xS(i),yS(i),zS(i)
! debug ends


           end do
           read_doneC = .true.
        endif
        
        if ( line == "gaussi: Final energy (eV):" ) read_final = .true.
        if(read_doneF .and. read_doneC .and. read_final ) exit
     end do
     close(FGAUSS)

     if (success) exit
  end do

! debug
 write (*,*)
 write (*,*)
 write (*,*) 'Running with good gaussian_force: '

  call center(forca,VECSIZE)
end subroutine

! This subroutine is not used for siesta
subroutine init_potential_gau()
! placeholder to define default parameters if desired
end subroutine

