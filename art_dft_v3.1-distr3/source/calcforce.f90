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
subroutine calcforce(nat, posa, boxl, forca, energy, evalf_number, conv )
   use defs, only : energy_type,typat
   implicit none

   !Arguments
   integer,      intent(in)                            :: nat
   real(kind=8), intent(in),  dimension(3*nat)         :: posa
   real(kind=8), dimension(3), intent(inout)           :: boxl
   real(kind=8), intent(out), dimension(3*nat)         :: forca
   real(kind=8), intent(out)                           :: energy
   integer,      intent(inout)                         :: evalf_number
   logical,      intent(in)                            :: conv


   if(energy_type == "SWP")  then
      call SWcalcforce(nat,posa,boxl,forca, energy)
      evalf_number = evalf_number +1
   else if (energy_type == "DFT")   then 
      call calcforce_dft(nat,typat,posa,boxl,forca,energy)
      evalf_number = evalf_number +1
   else if (energy_type == "GAU")   then
      call calcforce_gau(nat,typat,posa,boxl,forca,energy)

! debug starts
write(*,*)
write(*,*) 'nat: ',  nat
write(*,*) 'typat: ', typat
write(*,*)'posa:'
do i=1,nat
        write(*,*) posa(i), posa(i+nat), posa(i+2*nat)
end do
!write(*,*) 'boxl: ', boxl
write(*,*) 'forces: ', NEW_LINE('A'), forca
write(*,*) 'energy: ', energy

      evalf_number = evalf_number +1
   endif

END SUBROUTINE calcforce
