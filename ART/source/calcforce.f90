!> @file
!!   Subroutine which calls the right for type inside art
!! @author
!!   Written by Laurent Karim Beland, UdeM 2011!!
!!   Copyright (C) 2010-2011 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Subroutine which calls the right for type inside art
subroutine calcforce(nat, posa, forca, energy, evalf_number, conv )
   use defs, only : energy_type,typat
   implicit none

   !Arguments
   integer,      intent(in)                            :: nat
   integer                                             :: i, j, k
   real(kind=8), intent(in),  dimension(3*nat)         :: posa
   real(kind=8), intent(out), dimension(3*nat)         :: forca
   real(kind=8), intent(out)                           :: energy
   integer,      intent(inout)                         :: evalf_number
   logical,      intent(in)                            :: conv

  if (energy_type == 'GAU') then
      write(*,*) "Calling Gaussian energy forces"
      call gaussian_energy_forces(nat,typat,posa,forca,energy)
  elseif (energy_type == 'CPK') then
      write(*,*) "Calling CP2K energy forces"
      call cp2k_energy_forces(nat,typat,posa,forca,energy)
  endif

  evalf_number = evalf_number + 1

END SUBROUTINE calcforce

! In this routine, is the interface between ART and Gaussian

subroutine gaussian_energy_forces ( nat, typa, posa, forca, energy )
  use defs
  integer, intent(in) :: nat
  character(len=3) :: temp_string
  !character(len=3), dimension(NATOMS) :: typa
  real*8, intent(in), dimension(VECSIZE),target :: posa
  real*8, dimension(VECSIZE), target:: posG
  real*8, intent(out), dimension(VECSIZE),target :: forca
  real*8, intent(out) :: energy

  integer :: i, ierror
  integer, parameter :: FGAUSS = 21
  real*8, parameter :: ZERO = 0.0d0
  character(len=20) :: GAUSS   = 'art2gaussian.inp'
  character(len=20) :: GAUSSFORCE = 'gaussian2art'
  character(len=40) :: line
  logical :: read_energy_done,read_force_done,read_coordinates_done
  real(8), dimension(:), pointer :: xG, yG, zG
  real(8), dimension(:), pointer :: xa, ya, za
  real(8), dimension(:), pointer :: fax, fay, faz   ! Pointers for working force

  ! We first set-up pointers for the x, y, z components in the position and forces
  xa => posa(1:NATOMS)
  ya => posa(NATOMS+1:2*NATOMS)
  za => posa(2*NATOMS+1:3*NATOMS)

  posG = posa


  xG => posG(1:NATOMS)
  yG => posG(NATOMS+1:2*NATOMS)
  zG => posG(2*NATOMS+1:3*NATOMS)

  fax => forca(1:NATOMS)
  fay => forca(NATOMS+1:2*NATOMS)
  faz => forca(2*NATOMS+1:3*NATOMS)
  
  open(PATH, file = 'pathway.xyz', action = 'write', position = 'append')
  write(PATH,*) NATOMS
  write(PATH,'(a10)') "ART_COORDS"
  do i=1,NATOMS
       write(PATH,'((1X,a),3(2X,f15.8))') typat(i), xa(i), ya(i), za(i)
  enddo
  close(PATH)


  ! We first write to a file the format requested by Gaussian
  call system('python update_gaussian_header.py')
  open(unit=FGAUSS,file=GAUSS,status='unknown',action='write', position = 'append', iostat=ierror)

  ! Prepares gaussian input file coordinates
  do i = 1, NATOMS
      write(FGAUSS,"(a, f14.8, f14.8, f14.8)") typat(i), xa(i), ya(i), za(i)
  end do
  write(FGAUSS,*)    
     
  close(FGAUSS)

  call system('python execute_gaussian.py')
  open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)

  read_force_done = .false.
  read_energy_done = .false.
  read_coordinates_done = .false.

  do 
      read(FGAUSS,"(A40)", end=300) line

      !Gets the Gaussian output coordinates
      if ( line  == "outcoor:" ) then
          do i = 1, NATOMS
              read(FGAUSS,*) xG(i),yG(i),zG(i)
              write(*,*) xG(i),yG(i),zG(i)
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

  300 rewind(FGAUSS)

  close(FGAUSS)

  call center(forca,VECSIZE)

end subroutine gaussian_energy_forces

subroutine cp2k_energy_forces(nat,typa,posa,forca,energy)
  use defs
  integer, intent(in) :: nat
  character(len=3) :: temp_string
  !character(len=3), dimension(NATOMS) :: typa
  real*8, intent(in), dimension(VECSIZE),target :: posa
  real*8, dimension(VECSIZE), target:: posG
  real*8, intent(out), dimension(VECSIZE),target :: forca
  real*8, intent(out) :: energy

  integer :: i, ierror
  integer, parameter :: FCPK = 21
  real*8, parameter :: ZERO = 0.0d0
  character(len=20) :: CPK_COORD   = 'cp2k_coords.xyz'
  character(len=20) :: CPK_ENE_FOR = 'cp2k_output'
  character(len=40) :: line
  logical :: read_energy_done,read_forces_done
  real(8), dimension(:), pointer :: xa, ya, za
  real(8), dimension(:), pointer :: fax, fay, faz   ! Pointers for working force

  ! We first set-up pointers for the x, y, z components in the position and forces
  xa => posa(1:NATOMS)
  ya => posa(NATOMS+1:2*NATOMS)
  za => posa(2*NATOMS+1:3*NATOMS)

  fax => forca(1:NATOMS)
  fay => forca(NATOMS+1:2*NATOMS)
  faz => forca(2*NATOMS+1:3*NATOMS)
  
  open(PATH, file = 'pathway.xyz', action = 'write', position = 'append')
  write(PATH,*) NATOMS
  write(PATH,'(a10)') "ART_COORDS"
  do i=1,NATOMS
       write(PATH,'((1X,a),3(2X,f15.8))') typat(i), xa(i), ya(i), za(i)
  enddo
  close(PATH)

  ! We first write to a file the format requested by CP2K
  open(unit=FCPK,file=CPK_COORD,status='unknown', action='write', iostat=ierror)

  write(FCPK,*) NATOMS
  write(FCPK,"(a40)") "coords for cp2k"

  do i = 1, NATOMS
      write(FCPK,"(a, f14.8, f14.8, f14.8)") typat(i), xa(i), ya(i), za(i)
  end do

  !Blank line at the end of coordinates
  write(FCPK,*)    
     
  close(FCPK)

  write(*,*) "single point energy calculation"

  call system('python update_cp2k_header.py -k ef')

  call system('python execute_cp2k.py > cp2k_energy_forces')

  open(unit=FCPK,file=CPK_ENE_FOR,status='old',action='read',iostat=ierror)

  read_energy_done = .false.
  read_forces_done = .false.

  do 
      read(FCPK,"(A40)") line
        
      !Gets the final energy
      if ( line  == "energy:" ) then
          read(FCPK,*) energy
          ! ! hartree_to_ev=27.2113838668 
          energy = energy*27.2113838668
          read_energy_done = .true.
      endif
        
      !Gets the forces
      if ( line  == "forces:" ) then
         do i = 1, NATOMS
            read(FCPK,*) fax(i),fay(i),faz(i)
            ! ! unit conversion CP2K format Forces (Hartrees/Bohr) to CP2K format forces forces (eV/Ang):
            ! ! hartree_to_ev=27.2113838668 and bohr_to_angstrom=0.52917721092
            fax(i) = fax(i)*51.4220629786602
            fay(i) = fay(i)*51.4220629786602
            faz(i) = faz(i)*51.4220629786602
         end do
         read_forces_done = .true.
      endif
       
      !Checks that 
      if(read_energy_done .and. read_forces_done ) exit

  end do

  close(FCPK)

  call center(forca,VECSIZE)

end subroutine cp2k_energy_forces

! This subroutine is not used for siesta
subroutine init_potential_gau()
! placeholder to define default parameters if desired
end subroutine

