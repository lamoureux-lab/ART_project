! In this routine, is the  interface between ART and Gaussian

subroutine calcforce_gau(nat,typa,posa,boxl,forca,energy)
  use defs
  integer, intent(in) :: nat
  character(len=3) :: temp_string
  !character(len=3), dimension(NATOMS) :: typa
  real*8, intent(in), dimension(VECSIZE),target :: posa
  real*8, dimension(VECSIZE), target:: posS
  real*8, intent(in) :: boxl
  real*8, intent(out), dimension(VECSIZE),target :: forca
  real*8, intent(out) :: energy

  integer :: i, ierror
  integer :: end_of_file = -1
  integer, parameter :: FGAUSS = 21
  real*8, parameter :: ZERO = 0.0d0
  character(len=20) :: GAUSS   = 'art2gaussian.inp'
  character(len=20) :: GAUSSFORCE = 'gaussian2art'
  character(len=40) :: line
  logical :: read_energy_done,read_force_done,read_coordinates_done, success
  real(8), dimension(:), pointer :: xS, yS, zS
  real(8), dimension(:), pointer :: xa, ya, za
  real(8), dimension(:), pointer :: fax, fay, faz   ! Pointers for working force

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
     call system('python create_header.py')
     open(unit=FGAUSS,file=GAUSS,status='unknown',action='write', position = 'append', iostat=ierror)

     ! Prepares gaussian input file coordinates
     do i = 1, NATOMS
        write(FGAUSS,"(a, f14.8, f14.8, f14.8)") typat(i), xa(i), ya(i),za(i)
     end do
     write(FGAUSS,*)    
     
     close(FGAUSS)

     call system('python execute_gaussian.py')
     open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)

     read_force_done = .false.
     read_energy_done = .false.
     read_coordinates_done = .false.

     do 
        read(FGAUSS,"(A40)") line
        ! Checks for IO errors or EOF marker
        if (ierror .eq. end_of_file) then
           if (.not.success) then
              write(*,*) 'IO error in gaussian_force.f90'
              stop
           endif
           success = .false.
           exit
        endif

        !Gets the Gaussian output coordinates
        if ( line  == "outcoor:" ) then
           do i = 1, NATOMS
              read(FGAUSS,*) xS(i),yS(i),zS(i)
              write(*,*) xS(i),yS(i),zS(i)
           end do
           read_coordinates_done = .true.
        endif

        !Gets the final energy
        if ( line  == "energy:" ) then
           read(FGAUSS,*) energy
           read_energy_done = .true.
        endif
        
        !Gets the forces
        if ( line  == "forces:" ) then
           do i = 1, NATOMS
              read(FGAUSS,*) fax(i),fay(i),faz(i)
              ! ! unit conversion Gaussian format Forces (Hartrees/Bohr) to Gaussian format forces forces (eV/Ang):
              ! ! hartree_to_ev=27.2113838668 and bohr_to_angstrom=0.52917721092
              fax(i) = fax(i)*51.4220629786602
              fay(i) = fay(i)*51.4220629786602
              faz(i) = faz(i)*51.4220629786602
           end do
           read_force_done = .true.
        endif
       
        !Checks that 
        if(read_force_done .and. read_coordinates_done .and. read_energy_done ) exit
     end do
     close(FGAUSS)

  if(success) exit
  end do

  call center(forca,VECSIZE)
end subroutine

! This subroutine is not used for siesta
subroutine init_potential_gau()
! placeholder to define default parameters if desired
end subroutine

