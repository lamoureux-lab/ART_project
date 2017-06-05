! In this routine, is the  interface between ART and SIESTA

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
! bharat starts
  integer :: idum1, idum2
! bharat ends
  integer :: end_of_file = -1
  integer, parameter :: FSIESTA = 21
  real*8, parameter :: ZERO = 0.0d0
  character(len=20) :: SIESTA   = 'art2gaussian.inp'
  character(len=20) :: SIESTAFORCE = 'log'
  character(len=40) :: line
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

     ! We first write to a file the format requested by SIESTA
     open(unit=FSIESTA,file=SIESTA,status='replace',action='write',iostat=ierror)

!     write(FSIESTA,"('MD.TypeOfRun         cg')")
!     write(FSIESTA,"('MD.NumCGsteps         0')")
!     write(FSIESTA,"('MD.MaxCGDispl         0.1  Ang')")
!     write(FSIESTA,"('MD.MaxForceTol        0.002  eV/Ang')")
!     write(FSIESTA,*)
!     write(FSIESTA,"('%block Atomic_Coordinates_and_Atomic_Species')")
!     write(FSIESTA,"(f14.8, f14.8, f14.8, i4)")  ( xa(i), ya(i), za(i), typa(i), i=1, NATOMS )
!     write(FSIESTA,"('%endblock Atomic_Coordinates_and_Atomic_Species')")
!     write(FSIESTA,"('%include basicinfo.fdf')")
!     close(FSIESTA)

! bharat starts for gaussian input file preparation
     write(FSIESTA,"('%chk=temp.chk')")
     write(FSIESTA,"('#rhf/3-21g nosymm force')")
     write(FSIESTA,*)
     write(FSIESTA,"('name')")
     write(FSIESTA,*)
     write(FSIESTA,"('0 1')")
     write(FSIESTA,"(i4, f14.8, f14.8, f14.8)")  ( typa(i), xa(i), ya(i),za(i), i=1, NATOMS)
     write(FSIESTA,*)
     close(FSIESTA) 

! bharat ends





     ! We now call siesta do to the minimization
     call  system('sh execute_gaussian.sh')
     do i=1, 10000
        toto = dexp ( i * 0.001d0)
     end do
!OK  write (*,*) 'bharat test1'
     ! We must now read the forces from siesta's output file
     open(unit=FSIESTA,file=SIESTAFORCE,status='old',action='read',iostat=ierror)
!OK  write (*,*) 'bharat test2'
     read_done = .false.
     do 
        read(FSIESTA,"(A40)") line
! OK write (*,*) 'bharat test3', line
        if ( line  == "siesta: Final energy (eV):" ) then
           do i = 1, 7
              read(FSIESTA,"(A40)") line
! bharat debug starts
!OK       write(*,*) 'Line 1: ', line
! bharat debug ends
           end do
           read(FSIESTA,"(A24,f14.6)") line, energy
! bharat debug starts
!       write(*,*) 'Line 2: energy ', energy
!     stop
           
! bharat debug ends
           read_done = .true.
        endif
        if(read_done) exit
     end do
     close(FSIESTA)
     
     ! We now do the force
     open(unit=FSIESTA,file=SIESTAFORCE,status='old',action='read',iostat=ierror)
     read_doneF = .false.
     read_final = .false.
     read_doneC = .false.
     do 
        read(FSIESTA,"(A40)") line
        if (ierror .eq. end_of_file) then
           if (.not.success) then
              write(*,*) 'Error with siesta force'
              stop
           endif
           success = .false.
           exit
        endif
        if ( line  == "siesta: Atomic forces (eV/Ang):" ) then
           do i = 1, NATOMS
! bharat              read(FSIESTA,"(i7,3f12.6)") idum,fax(i),fay(i),faz(i)
        read(FSIESTA,"(i7,i10,f21.9,f15.9,f15.9)") idum,idum1, fax(i),fay(i),faz(i)
!              read(FSIESTA,"(i7,i7,3f13.9)") idum, idum,fax(i)*51.4220629786602,fay(i)*51.4220629786602,faz(i)*51.4220629786602
!              read(FSIESTA,"(i8,i10,3f12.9)") idum,idum,fax(i),fay(i),faz(i)

! unit conversion Gaussian format Forces (Hartrees/Bohr) to Siesta format forces forces (eV/Ang): 
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
!              read(FSIESTA,*) xS(i),yS(i),zS(i)
! bharat debug starts
! write (*,*) "position formatting enters"
        read(FSIESTA,"(i7,i12,i12,f15.6,f13.6,f11.6)") idum,idum1,idum2,xS(i),yS(i),zS(i)
! write (*,*) idum,idum1,idum2,xS(i),yS(i),zS(i)
!stop
! bharat debug ends


           end do
           read_doneC = .true.
        endif
        
        if ( line == "siesta: Final energy (eV):" ) read_final = .true.
        if(read_doneF .and. read_doneC .and. read_final ) exit
     end do
     close(FSIESTA)

     if (success) exit
  end do

! bharat starts
 write (*,*) 'bharat running good gaussian_force'
!stop
! bharat ends
  
  call center(forca,VECSIZE)
end subroutine

! This subroutine is not used for siesta
subroutine init_potential_gau()

end subroutine

