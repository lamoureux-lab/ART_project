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
   real(kind=8),dimension(3) :: boxl

   real(kind=8), dimension(3*natoms)         :: pos_temp
   integer :: nat
   integer, dimension(natoms) :: numnei
   integer, dimension(natoms,maxnei) :: nei
   real(kind=8), dimension(3) :: invbox
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
   if (energy_type == "DFT" ) then
      call min_converge_dft(success)
   elseif (energy_type == "SWP") then
      call min_converge_fire(success)
   endif

   if ( .not. success ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
           & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(1X,A)') "Minimization exited before the geometry optimization converged,"
      write(FLOG,'(1X,A)') "this minimum will be rejected."
      close(FLOG) 
   end if

   write(*,"('',' BART: Relaxed energy : ',(1p,e17.10,0p))") total_energy

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

   write(*,*) "BART: INIT LANCZOS"  !debug
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
      write(*,*) 'BART: Iter ', i, ' : ', lanc_energy, proj_energy,  eigenvalue, a1  

      ! Now we start from the previous direction. 
      new_projection= .false.   
      ! let's see the projection
      if ( setup_initial ) call print_proj ( i, stage, projection, eigenvalue, DEL_LANCZOS )
   end do

   if (energy_type == "SWP") call reset_SW_potential()

   ! Report 
   write(*,*) "BART: END  LANCZOS"  !debug
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
   call calcforce(NATOMS, pos, box, force, total_energy, evalf_number, conv)
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
      call calcforce(NATOMS, pos, box, forceb, total_energy, evalf_number, conv)
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

   dtmax=0.15d0
   dt = dtmax*0.2d0

   success=.false.
   fnrm=1.d10
   velcur=0.0d0
   poscur=pos

   call calcforce(natoms,pos,box,fpred,total_energy,evalf_number,conv)

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
      call calcforce(natoms,pospred,box,fpred,total_energy,evalf_number,conv)
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
   dtmax=0.15d0
   dt = dtmax*0.2d0
   success=.false.
   fnrm=1.d10
   velcur=0.0d0
   poscur=pos

   call calcforce(natoms,pos,box,fpred,total_energy,evalf_number,conv)

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
      call calcforce(natoms,pos,box,fpred,total_energy,evalf_number,conv)
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



