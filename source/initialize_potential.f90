!> @file
!!   Initialize the potential
!!
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!  Modified by Laurent Karim Beland, UdeM, 2011. For working with QM/MM !!
subroutine initialize_potential( )
   use defs
   implicit None

   !Local variables
   integer :: ierror

   call init_potential_gau()
   call calcforce( NATOMS, pos, force, total_energy, evalf_number, .false. )

END SUBROUTINE initialize_potential




!> Finalize the potential
subroutine finalise_potential( )
   use defs, only : energy_type

   implicit none  

END SUBROUTINE finalise_potential
