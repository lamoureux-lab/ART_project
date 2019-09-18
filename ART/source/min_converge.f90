!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! Modified by Laurent in 2011 for working with QM/MM and FIRE


!> ART min_converge
!!   Minimizes the energy at constant volume.
!!   This minimization is done with only a minimal knowledge of the 
!!   physics of the problem so that it is portable
subroutine min_converge ( success )

   use defs
   implicit none

   !Arguments
   logical, intent(out) :: success

   !Local variables
   integer :: i, ierror
   real(kind=8), dimension(3*natoms)         :: pos_temp
   integer :: nat
   integer, dimension(natoms) :: numnei
   integer, dimension(natoms,maxnei) :: nei
   integer :: j,k 
   logical, dimension(natoms) :: is_at_quantum
   real(kind=8) :: xij,yij,zij,rij2
   !_______________________
   ! Report
   open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
        & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,'(1X,A)') ' RELAXATION'
   close(FLOG) 

   ! It is possible, here, to use a different minimization routine for each potential.

   if (energy_type == 'GAU') then
       write(*,*) "Calling min converge Gaussian"
       write(*,*) "energy_type ", energy_type
       call min_converge_gau(success)
   elseif (energy_type == 'CPK') then
       write(*,*) "Calling min converge CP2K"
       write(*,*) "energy_type ", energy_type
       call min_converge_cp2k(success)
   endif

   if ( .not. success ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
           & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(1X,A)') "Minimization exited before the geometry optimization converged,"
      write(FLOG,'(1X,A)') "this minimum will be rejected."
      close(FLOG)
   end if

   write(*,"('',' ART: Relaxed energy : ',(1p,e17.10,0p))") total_energy

   return
END SUBROUTINE min_converge


subroutine check_min( stage )
   use defs
   use lanczos_defs
   implicit none

   !Arguments
   character(len=1), intent(in) :: stage

   !Local variables
   integer :: i, ierror, repetition
   logical :: new_projection
   real(kind=8) :: a1
   real(kind=8) :: min_energy ! First reference energy in lanczos


   ! We check how it changes the energy of the system by applying the projection.
   IN_MINIMUN = .True.
   ! Report


   open( unit = FLOG, file = LOGFILE, status = 'unknown',&
        & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,*) ' Starting Lanczos'
   close(FLOG)

   new_projection = .true.            ! We do not use any previously computed direction.

   if ( .not. setup_initial ) then
      ! We call lanczos twice.
      repetition = 2
   else
      ! if not, four times.
      repetition = 3
   end if

   write(*,*) "ART: INIT LANCZOS"  !debug
   do i = 1, repetition
      call lanczos( NVECTOR_LANCZOS_H, new_projection , a1 )
      ! Report
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
           & action = 'write', position = 'append', iostat = ierror )
      if ( i == 1 ) then           ! Our reference energy for the report.
         min_energy = lanc_energy
         write(FLOG,'(2x,A,f12.8)') ' Em= ', min_energy
         write(FLOG,'(A39)') '   Iter     Ep-Em (eV)   Eigenvalue  a1'
      end if
      write(FLOG,'(I6,3X,(1p,e10.2,0p),4X,F12.6,1X,F7.4)') i, proj_energy-min_energy, eigenvalue, a1
      close(FLOG)
      write(*,*) 'ART: Iter ', i, ' : ', lanc_energy, proj_energy,  eigenvalue, a1

      ! Now we start from the previous direction.
      new_projection= .false.
      ! let's see the projection
      if ( setup_initial ) call print_proj ( i, stage, projection, eigenvalue, DEL_LANCZOS )
   end do

   ! Report
   write(*,*) "ART: END  LANCZOS"  !debug
   open( unit = FLOG, file = LOGFILE, status = 'unknown',&
        & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,*) ' Done Lanczos'
   close(FLOG)

   ! Default value in the activation part is false.
   IN_MINIMUN = .False.

END SUBROUTINE check_min


!>  This module defines a number of parameters used during the minimization
MODULE minimization_sd
   implicit none
   save
   integer, parameter :: MAX_ITER = 1000
   real(kind=8), parameter :: FTHRESHOLD = 3.5D-1
   real(kind=8), parameter :: STEPSIZE = 1.0D-4  ! Size in angstroems

   real(kind=8), parameter :: FTHRESH2 = FTHRESHOLD * FTHRESHOLD
END MODULE minimization_sd


subroutine min_converge_sd(minimized)
   use minimization_sd
   use defs
   implicit none

   logical, intent(inout)  :: minimized

   real(kind=8),  dimension(VECSIZE):: forceb, posb, tmp_pos
   integer :: iter, i, npart
   real(kind=8) :: current_ftot2, ftot,ftot2, step, delr, current_energy
   logical :: conv

   conv = .false. !not sure about this variable

   ! We compute at constant volume
   call calcforce(NATOMS, pos, force, total_energy, evalf_number, conv)
   ftot2 = 0.0d0
   do i=1, VECSIZE
      ftot2 = ftot2 + force(i) * force(i)
   end do
   current_ftot2 = ftot2
   current_energy = total_energy ! if quantum, this energy is only quantum
   write(*,*)  'initial energy : ',  total_energy

   step = STEPSIZE
   do iter = 1, MAX_ITER

      ! Apply PBC
      posb = pos + step * force
      tmp_pos = pos
      pos = posb
      call calcforce(NATOMS, pos, forceb, total_energy, evalf_number, conv)
      pos = tmp_pos
      ftot2 = 0.0d0
      do i=1, VECSIZE
         ftot2 = ftot2 + forceb(i) * forceb(i)
      end do

      ! if(ftot2 < current_ftot2 .or. total_energy < current_energy ) then
      if(ftot2 < current_ftot2  ) then

         pos = posb
         force = forceb
         step = 1.2 * step
         current_ftot2 = ftot2
         current_energy = total_energy

      else
         step = 0.6 * step
      endif
      if(ftot2 < FTHRESH2) exit
      call displacement( posref, pos, delr, npart )
      ! Write

      if (modulo(iter,5) == 0 ) then
         call write_step ( 'M', iter, 0.0d0, current_energy )
         if ( SAVE_CONF_INT ) call save_intermediate( 'M' )

      endif

   end do

   ftot = sqrt(ftot2)
   if (ftot < FTHRESHOLD ) then
      write(*,*) 'Minimization successful   ftot : ', ftot
      write(*,*)  'final energy :', total_energy
      minimized = .true.
   else
      write(*,*) 'Minimization failed   ftot : ', ftot
      minimized = .false.
   endif
END SUBROUTINE min_converge_sd


!> Implementation of the damped MD based geometry optimizer FIRE, PRL 97, 170201 (2006)
!! The MD-Integrator is the common velocity verlet, all masses are equal to 1.d0
!! Implemented in August 2010, Maximilian Amsler, Basel University
!! Suggestion for maximal timestep as tmax=2*pi*sqrt(alphaVSSD)*1.2d-1
!! Choose the initial timestep as tinit=tmax*0.5d0
subroutine min_converge_fire(success)
   use defs
   use saddles
   use lanczos_defs

   implicit none

   !Arguments
   logical :: success
   !Local variables
   real(kind=8) :: fnrm
   real(kind=8) :: fmax,vmax
   real(kind=8), parameter :: pi = 3.14159265d0
   integer :: check
   integer :: iat

   logical :: conv

   !Fire parameters:
   real(kind=8) :: alpha,P,finc,fdec,falpha,alphastart,dt,dtmax,vnrm
   real(kind=8) :: velcur(3*natoms), velpred(3*natoms),poscur(3*natoms),pospred(3*natoms)
   real(kind=8) :: fcur(3*natoms),fpred(3*natoms),mass(3*natoms)
   real(kind=8) :: ecur,epred,anoise
   integer:: Nmin,nstep,it
   integer,parameter :: max_iter=100
   real(kind=8),parameter :: norm_criterium = 0.008d0
   real(kind=8),parameter :: fmax_criterium = 0.020d0
   integer :: miter
   real(8) :: initial_energy
   real(8), dimension(VECSIZE) :: initial_pos

   initial_energy = total_energy
   initial_pos = pos

   miter = 0
   conv = .false.


   check=0
   !Set FIRE parameters
   Nmin=5
   finc=1.1d0
   fdec=0.5d0
   !  alphastart=0.25d0
   alphastart =0.1d0

   anoise=1.0d-8


   alpha=alphastart
   falpha=0.99d0
   nstep=1

   ! Optimization for Gaussian or CP2K
   if (energy_type == "GAU" .or. energy_type == 'CPK') then
      dtmax=0.3d0
   else
      dtmax=0.15d0
   end if
   dt = dtmax*0.2d0

   success=.false.
   fnrm=1.d10
   velcur=0.0d0
   poscur=pos

   call calcforce(natoms,pos,fpred,total_energy,evalf_number,conv)

   fcur=force
   mass=1.0d0
   ecur=total_energy
   epred=total_energy


   do it=1,max_iter
      miter = miter + 1
      pas = pas + 1

      do
         do iat=1,3*natoms
            pospred(iat)=poscur(iat)+dt*velcur(iat)+dt*dt*0.5d0*fcur(iat)/mass(iat)
         enddo
         call displacement(pospred, pos, delr, npart)
         if (delr .lt. 0.25d0) exit
         write(*,*) "FIRE: problem with explosion, trying a smaller dt -delr", delr
         dt = 0.5 * dt
         if (dt .lt. 0.005d0*dtmax) then
            return
         end if
      enddo

      pos = pospred
      call calcforce(natoms,pospred,fpred,total_energy,evalf_number,conv)
      force=fpred

      !call fnrmmax_fire(fpred,fnrm,fmax,natoms)
      !  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)
      fmax = dsqrt(MAXVAL(fpred(1:NATOMS)**2 + fpred(NATOMS+1:2*NATOMS)**2 + fpred(2*NATOMS+1:3*NATOMS)**2))
      fnrm = SUM(fpred(1:NATOMS)**2 + fpred(NATOMS+1:2*NATOMS)**2 + fpred(2*NATOMS+1:3*NATOMS)**2)

      do iat=1,3*natoms
         velpred(iat)=velcur(iat)+0.5d0*dt*(fpred(iat))/mass(iat)+0.5d0*dt*fcur(iat)/mass(iat)
      enddo

      P=dot_product(fpred,velpred)
      !call fnrmmax_fire(velpred,vnrm,vmax,natoms)
      vnrm = SUM(velpred(1:NATOMS)**2 + velpred(NATOMS+1:2*NATOMS)**2 + velpred(2*NATOMS+1:3*NATOMS)**2)

      if (modulo(miter,5) == 0 ) then
         call write_step ( 'M', miter, 0.0d0, total_energy )
         write(*,*) "fnrm",fnrm,"fmax",fmax
         pos = pospred
         if ( SAVE_CONF_INT ) call save_intermediate( 'M' )
      endif

      !force = fpred
      !ftot = dsqrt(fnrm)
      !delta_e = total_energy - ref_energy
      !fpar = fmax
      !fperp = 0.0d0
      !eigenvalue = 0.0d0
      ! Magnitude of the displacement (utils.f90).
      !call displacement( posref, poscur, delr, npart )
      ! Write

      if( fnrm < norm_criterium .or. fmax < fmax_criterium) then
         if (total_energy < initial_energy) then
             !pos = pospred
             success = .true.
             exit
         else
             write(*,*) "FIRE: final energy higher than initial energy, try steepest descent instead"
             pos = initial_pos
             call min_converge_sd(success)
             return
         end if
      endif

      !Update variables
      fcur=fpred
      poscur=pospred
      !Normal verlet velocity update
      !  velcur=velpred

      !!FIRE Update
      !call fnrmmax_fire(fpred,fnrm,fmax,natoms)
      fnrm=sqrt(fnrm)
      !call fnrmmax_fire(velpred,vnrm,vmax,natoms)
      vnrm=sqrt(vnrm)
      !Modified velocity update, suggested by Alireza
      !velcur(:)=(1.0d0-alpha)*velpred(:)+fpred(:)*min(alpha*vnrm/fnrm,2.0d0*in%betax)!alpha*fpred(:)/fnrm*vnrm
      !Original FIRE velocitiy update
      velcur(:)=(1.0d0-alpha)*velpred(:)+alpha*fpred(:)/fnrm*vnrm
      if(P.gt.-anoise*vnrm .and. nstep.gt.Nmin) then
         dt=min(dt*finc,dtmax)
         !         alpha=max(alpha*falpha,0.1d0) !Limit the decrease of alpha
         alpha=alpha*falpha
         elseif(P.le.-anoise*vnrm) then
         nstep=0
         dt=dt*fdec
         velcur=0.d0
         alpha=alphastart
      endif
      nstep=nstep+1
   enddo

   if (success .eqv. .false.) then
     write(*,*) "FIRE: failed to minimize forces, trying steepest descent"
     pos = initial_pos
     call min_converge_sd(success)
   else
     write(*,*) "FIRE: minimization successful"
     write(*,*) "fnrm: ", fnrm, "  , fmax: ", fmax," , iter: ", it
     write(*,*) "Final energy: ", total_energy
   endif

   return
END SUBROUTINE min_converge_fire


subroutine fnrmmax_fire(force,fnrm,fmax,natoms)
   implicit none
   integer,intent(in) :: natoms
   real(kind=8),dimension(3*natoms),intent(in) :: force
   real(kind=8),intent(out) :: fnrm,fmax

   integer :: i

   real(kind=8) :: nombre

   fnrm = 0.0d0
   fmax = 0.0d0
   !  force2 = dot_product(force)
   do i = 1,natoms
      nombre = force(i)**2 + force(i+natoms)**2 + force(i+natoms+natoms)**2
      if (nombre > fmax) fmax = nombre
      fnrm = fnrm+nombre
   enddo
   fmax = dsqrt(fmax)

END SUBROUTINE fnrmmax_fire


!> Implementation of the damped MD based geometry optimizer FIRE, PRL 97, 170201
!(2006)
!! The MD-Integrator is the common velocity verlet, all masses are equal to 1.d0
!! Implemented in August 2010, Maximilian Amsler, Basel University
!! Suggestion for maximal timestep as tmax=2*pi*sqrt(alphaVSSD)*1.2d-1
!! Choose the initial timestep as tinit=tmax*0.5d0
subroutine perp_fire(success,max_iter)
   use defs
   use saddles
   use lanczos_defs

   implicit none

   !Arguments
   logical :: success
   integer, intent(in) :: max_iter
   !Local variables
   real(kind=8) :: fnrm
   real(kind=8) :: fmax,vmax
   real(kind=8), parameter :: pi = 3.14159265d0
   integer :: check
   integer :: iat

   logical :: conv

   !Fire parameters:
   real(kind=8) :: alpha,P,finc,fdec,falpha,alphastart,dt,dtmax,vnrm
   real(kind=8) :: velcur(3*natoms),velpred(3*natoms),poscur(3*natoms),pospred(3*natoms)
   real(kind=8) :: fcur(3*natoms),fpred(3*natoms),mass(3*natoms)
   real(kind=8) :: ecur,epred,anoise
   real(kind=8), dimension(VECSIZE) :: perp_force   ! Perpendicular force...
   integer:: Nmin,nstep,it
   real(kind=8),parameter :: norm_criterium = 0.008d0
   real(kind=8),parameter :: fmax_criterium = 0.020d0
   integer :: miter
   real(8), dimension(VECSIZE) :: initial_pos
   real(8) :: initial_energy

   miter = 0
   conv = .false.

   initial_pos = pos
   initial_energy = total_energy

   check=0
   !Set FIRE parameters
   Nmin=5
   finc=1.1d0
   fdec=0.5d0
   !  alphastart=0.25d0
   alphastart =0.1d0
   anoise=1.0d-8
   alpha=alphastart
   falpha=0.99d0
   nstep=1

   ! Optimization for Gaussian or CP2K
   if (energy_type == "GAU" .or. energy_type == 'CPK') then
      write(*,*) "Test1:", energy_type,  "is properly set" 
      dtmax=0.3d0
   else
      dtmax=0.15d0
   end if

   dt = dtmax*0.2d0
   success=.false.
   fnrm=1.d10
   velcur=0.0d0
   poscur=pos

   call calcforce(natoms,pos,fpred,total_energy,evalf_number,conv)

   fpar = dot_product(fpred, projection)
   perp_force = fpred - fpar * projection
   fpred = perp_force

   fcur=fpred
   mass=1.0d0
   ecur=total_energy
   epred=total_energy


   do it=1,max_iter
    
      do

        miter = miter + 1
        pas = pas + 1
        do iat=1,3*natoms
           pospred(iat)=poscur(iat)+dt*velcur(iat)+dt*dt*0.5d0*fcur(iat)/mass(iat)
        enddo
        call displacement(pospred, pos, delr, npart)
        if (delr .lt. 0.25d0) exit
           write(*,*) "PERP_FIRE: problem with explosion, trying a smaller dt-delr: ", delr
           dt = 0.5d0*dt
           if (dt .lt. 0.0050d0*dtmax) then
              return
           endif
      enddo

      pos = pospred
      call calcforce(natoms,pos,fpred,total_energy,evalf_number,conv)
      force=fpred
      fpar = dot_product(fpred, projection)
      perp_force = fpred - fpar * projection
      fpred = perp_force
      fmax = dsqrt(MAXVAL(fpred(1:NATOMS)**2 + fpred(NATOMS+1:2*NATOMS)**2 + fpred(2*NATOMS+1:3*NATOMS)**2))
      fnrm = SUM(fpred(1:NATOMS)**2 + fpred(NATOMS+1:2*NATOMS)**2 + fpred(2*NATOMS+1:3*NATOMS)**2)
      !  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)
      do iat=1,3*natoms
         velpred(iat)=velcur(iat)+0.5d0*dt*(fpred(iat))/mass(iat)+0.5d0*dt*fcur(iat)/mass(iat)
      enddo
      P=dot_product(fpred,velpred)
      !call fnrmmax_fire(velpred,vnrm,vmax,natoms)
      vnrm = SUM(velpred(1:NATOMS)**2 + velpred(NATOMS+1:2*NATOMS)**2 + velpred(2*NATOMS+1:3*NATOMS)**2)

      if( (it > 1) .and.  (fnrm < fpar )) then
        if (total_energy < initial_energy ) then
!           pos = pospred
           success = .true.
           exit
        else
           !FIRE has failed, we will try SD in saddle_converge
           !write(unite,*) "FIRE gave an energy greater than initial. we will
           !try steepest descent"
           pos = initial_pos
           !call min_converge_sd_local(success)
           success = .false.
           exit
        endif
     endif

      !if( fnrm < norm_criterium .or. fmax < fmax_criterium) then
      !   pos = pospred
      !   success = .true.
      !   exit
      !endif

      !Update variables
      fcur=fpred
      poscur=pospred
      !Normal verlet velocity update
      !  velcur=velpred

      !!FIRE Update
!      call fnrmmax_fire(fpred,fnrm,fmax,natoms)
      fnrm=sqrt(fnrm)
!      call fnrmmax_fire(velpred,vnrm,vmax,natoms)
      vnrm=sqrt(vnrm)
      !Modified velocity update, suggested by Alireza
      !velcur(:)=(1.0d0-alpha)*velpred(:)+fpred(:)*min(alpha*vnrm/fnrm,2.0d0*in%betax)!alpha*fpred(:)/fnrm*vnrm
      !Original FIRE velocitiy update
      velcur(:)=(1.0d0-alpha)*velpred(:)+alpha*fpred(:)/fnrm*vnrm

      if(P.gt.-anoise*vnrm .and. nstep.gt.Nmin) then
         dt=min(dt*finc,dtmax)
         !         alpha=max(alpha*falpha,0.1d0) !Limit the decrease of alpha
         alpha=alpha*falpha
      else if(P.le.-anoise*vnrm) then
         nstep=0
         dt=dt*fdec
         velcur=0.d0
         alpha=alphastart
      endif
      nstep=nstep+1
   enddo

   if (it == max_iter+1) success = .true.

   return
END SUBROUTINE perp_fire

! In this routine, we interface between ART and GAUSSIAN
subroutine min_converge_gau(success)
  use defs

  integer :: i, ret

  logical :: success
  integer, parameter :: FGAUSS = 21
  real*8, parameter :: ZERO = 0.0d0
  real*8, dimension(natoms) :: xx, yy, zz
  character(len=20) :: GAUSS   = 'art2gaussian.inp'
  character(len=20) :: GAUSSFORCE = 'gaussian2art'
  character(len=70) :: line
  character(len=10) :: string_natoms
  real*8, dimension(VECSIZE), target :: forcea_opt
  real(8), dimension(:), pointer :: fax_opt, fay_opt, faz_opt   ! Pointers for working force
  logical :: opt_coords_read, opt_energy_read, opt_forces_read

  fax_opt => forcea_opt(1:NATOMS)
  fay_opt => forcea_opt(NATOMS+1:2*NATOMS)
  faz_opt => forcea_opt(2*NATOMS+1:3*NATOMS)

  success = .false.
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'start minimization'
  close(flog)

  ! We first write to a file the format requested by GAUSS

  write(*,*) "Optimizing now"
  call system('python update_gaussian_header.py -k opt')

  open(unit=FGAUSS,file=GAUSS,status='unknown',action='write', position = 'append', iostat=ierror)

  ! Prepares gaussian input file coordinates
  do i = 1,NATOMS
  write(FGAUSS,"(a, f14.8, f14.8, f14.8)") typat(i),  x(i), y(i), z(i)
  enddo
  write(FGAUSS,*) !Blank line... Gaussian expects one blank line at the end of the input file.
  close(FGAUSS)

  ! We now call Gaussian do to the minimization

  call system('python execute_gaussian.py')
  
  ! We must now read the energy and positions from gaussian's output file
  open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)
  opt_coords_read = .false.
  opt_energy_read = .false.
  opt_forces_read = .false.
  do 
    !Gets the Gaussian output coordinates
    read(FGAUSS,'(A40)') line
    if ( line  == "outcoor:" ) then
      do i = 1, NATOMS
        read (FGAUSS,*) x(i), y(i), z(i)
      end do
        opt_coords_read = .true.
    endif

    if (line == 'energy:') then
        read(FGAUSS,*) total_energy
        opt_energy_read = .true.
    endif

      !Gets the forces
    if ( line  == "forces:" ) then
       do i = 1, NATOMS
          read(FGAUSS,*) fax_opt(i),fay_opt(i),faz_opt(i)
            ! ! unit conversion Gaussian format Forces (Hartrees/Bohr) to Gaussian format forces forces (eV/Ang):
            ! ! hartree_to_ev=27.2113838668 and bohr_to_Angstrom=0.52917721092
          fax_opt(i) = fax_opt(i)*51.4220629786602
          fay_opt(i) = fay_opt(i)*51.4220629786602
          faz_opt(i) = faz_opt(i)*51.4220629786602
       end do
       opt_forces_read = .true.
    endif
    if (opt_coords_read .and. opt_energy_read .and. opt_forces_read) exit    
  end do

  close(FGAUSS)

  force_opt = forcea_opt

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'End minimization'
  close(flog)
  success = .true.
  return

end subroutine min_converge_gau

! In this routine, we interface between ART and CP2K
subroutine min_converge_cp2k(success)
  use defs

  integer :: i, ret

  logical :: success
  integer, parameter :: FCPK = 21
  real*8, parameter :: ZERO = 0.0d0
  real*8, dimension(natoms) :: xx, yy, zz
  character(len=20) :: CPK_COORD   = 'cp2k_coords.xyz'
  character(len=20) :: CPK_ENE_FOR = 'cp2k_output'
  character(len=20) :: CPK_OPT_COORD = 'cp2k-pos-1.xyz'
  character(len=70) :: line, dummy
  character(len=10) :: string_natoms

  
  success = .false.
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'start minimization'
  close(flog)

  ! We first write to a file the format requested by CP2KSS

  open(unit=FCPK,file=CPK_COORD,status='unknown',action='write',iostat=ierror)

  write(FCPK,*) NATOMS
  write(FCPK,*) "cp2k coords"
  do i = 1,NATOMS
  write(FCPK,"(a, f14.8, f14.8, f14.8)") typat(i),  x(i), y(i), z(i)
  enddo
  write(FCPK,*)
  close(FCPK)

  ! We now call CP2K do to the minimization

  write(*,*) "Optimizing now"
  call system('python update_cp2k_header.py -k opt')
  call system('python execute_cp2k.py > cp2k_energy_forces')
  
  open(unit=FCPK,file=CPK_ENE_FOR,status='old',action='read',iostat=ierror)
  
  do 
      read(FCPK,"(A40)", end=300) line

      !Gets the final energy
      if ( line  == "energy:" ) then
          read(FCPK,*) energy
      endif
        
  end do

  300 rewind(FCPK)

  close(FCPK)

  ! We must now read the position from cp2k's output file
  open(unit=FCPK,file=CPK_OPT_COORD,status='old',action='read',iostat=ierror)
  do 
    read(FCPK,'(A40)',end=400) dummy
    read(FCPK,'(A40)') dummy
    do i = 1, NATOMS
        read (FCPK,*) typat(i), x(i), y(i), z(i)
        x(i) = x(i) - cell_a/2
        y(i) = y(i) - cell_b/2
        z(i) = z(i) - cell_c/2
    end do
  end do
  400 rewind(FCPK)
  close(FCPK)

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(flog,*) 'End minimization'
  close(flog)
  success = .true.
  return

end subroutine min_converge_cp2k

