!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!******************************************************************************
program main
!******************************************************************************
!
! Main file for lesgo solver
! Contains main time loop
!

use types, only : rprec
use clock_m
use param
use sim_param
use grid_m
use io, only : energy, output_loop, output_final, jt_total
use io, only : write_tau_wall_bot, write_tau_wall_top, kx_energy, kx_energy_fourier
!use io, only : ky_energy
#ifdef PPOUTPUT_CLOCK
use io, only : write_clocks
#endif
use fft
use derivatives, only : filt_da, ddz_uv, ddz_w
use derivatives, only : wave2physF
use derivatives, only : wave2phys, phys2wave
use test_filtermodule
use cfl_util
use sgs_param, only : nu
use sgs_stag_util, only : sgs_stag
use forcing
use functions, only: get_tau_wall_bot, get_tau_wall_top

#ifdef PPMPI
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif

#ifdef PPLVLSET
use level_set, only : level_set_global_CA, level_set_vel_err
use level_set_base, only : global_CA_calc
#endif

#ifdef PPTURBINES
use turbines, only : turbines_forcing, turbine_vel_init
#endif

#ifdef PPSTREAKS
use sgs_param, only: F_LM, F_MM, F_QN, F_NN
#endif

#ifdef PPRNL
use functions, only : x_avg
#endif

#ifdef PPGQL
use functions, only : gql_filter
#endif

#ifdef PPTLWMLES
use tlwmles, only : tlwm_wallstress
#endif

use messages

implicit none

#ifdef PPSTREAKS
real(rprec), dimension(:,:,:), allocatable :: dummyu, dummyv, dummyw
real(rprec), dimension(:,:,:), allocatable :: dummyRHSx, dummyRHSy, dummyRHSz
#endif

character (*), parameter :: prog_name = 'main'

integer :: jt_step, nstart
real(rprec) :: rmsdivvel, maxcfl, tt, maxvisc

type(clock_t) :: clock, clock_total
!type(clock_t) :: clock_forcing
#ifdef PPOUTPUT_CLOCK
type(clock_t) :: clock_deriv, clock_wm, clock_stress, clock_divstress
type(clock_t) :: clock_convec, clock_inter, clock_pres
#endif

! Measure total time in forcing function
!real(rprec) :: clock_total_f = 0.0

#ifdef PPMPI
! Buffers used for MPI communication
real(rprec) :: rbuffer
real(rprec) :: maxdummy ! Used to calculate maximum with mpi_allreduce
real(rprec) :: tau_top   ! Used to write top wall stress at first proc
#endif

real(rprec), dimension(:,:,:), allocatable :: tempRHS
allocate( tempRHS (ld, ny, lbz:nz) )

! Initialize MPI
#ifdef PPMPI
call mpi_init (ierr)
#endif

! Start the clocks, both local and total
call clock%start

! Initialize time variable
tt = 0
jt = 0
jt_total = 0

! Initialize all data
call initialize()

if(coord == 0) then
    call clock%stop
#ifdef PPMPI
    write(*,'(1a,E15.7)') 'Initialization wall time: ', clock % time
#else
    write(*,'(1a,E15.7)') 'Initialization cpu time: ', clock % time
#endif

endif

call clock_total%start

! Initialize starting loop index
! If new simulation jt_total=0 by definition, if restarting jt_total
! provided by total_time.dat
nstart = jt_total+1

! Declare variables for shifting the domain
! This gets rid of streaks in the domain
#ifdef PPSTREAKS
allocate( dummyu     (ld    ,ny, lbz:nz) )
allocate( dummyv     (ld    ,ny, lbz:nz) )
allocate( dummyw     (ld    ,ny, lbz:nz) )
allocate( dummyRHSx  (ld    ,ny, lbz:nz) )
allocate( dummyRHSy  (ld    ,ny, lbz:nz) )
allocate( dummyRHSz  (ld    ,ny, lbz:nz) )
#endif

! BEGIN TIME LOOP
time_loop: do jt_step = nstart, nsteps

    ! Get the starting time for the iteration
    call clock%start

    if (use_cfl_dt) then

        ! Go back to physical space to calculate CFL
        if (fourier) then
            call wave2physF( u, uF )
            call wave2physF( v, vF )
            call wave2physF( w, wF )
        endif

        dt_f = dt
        dt = get_cfl_dt()
        dt_dim = dt * z_i / u_star

        tadv1 = 1._rprec + 0.5_rprec * dt / dt_f
        tadv2 = 1._rprec - tadv1

    end if

    ! Advance time
    jt_total = jt_step
    jt = jt + 1
    total_time = total_time + dt
    total_time_dim = total_time_dim + dt_dim
    tt = tt+dt

    ! Save previous time's right-hand-sides for Adams-Bashforth Integration
    ! NOTE: RHS does not contain the pressure gradient
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz

    ! Attempt to trigger turbulence by lowering the viscosity
    if (trigger) then
        if (jt_total == trig_on) then
            nu_molec = nu_molec / trig_factor
            nu = nu / trig_factor
        endif
        if (jt_total == trig_off) then
            nu_molec = nu_molec*trig_factor
            nu = nu*trig_factor
        endif       
    end if

    ! Fourier check
    if (fourier .and. fourier_check) then
        if ( (mod(jt_total,fourier_nskip)==0) .and.             &
            (jt_total < tavg_nstart) ) then
            ! Transform to physical space
            call wave2phys( u, lbz )
            call wave2phys( v, lbz )
            call wave2phys( w, lbz )
            ! Transform back
            call phys2wave( u, lbz )
            call phys2wave( v, lbz )
            call phys2wave( w, lbz )
        endif
    endif

    ! Calculate velocity derivatives
#ifdef PPOUTPUT_CLOCK
    call clock_deriv%start
#endif
    ! Calculate dudx, dudy, dvdx, dvdy, dwdx, dwdy (in Fourier space)
    call filt_da(u, dudx, dudy, lbz)
    call filt_da(v, dvdx, dvdy, lbz)
    call filt_da(w, dwdx, dwdy, lbz)

    ! Calculate dudz, dvdz using finite differences (for 1:nz on uv-nodes)
    !  except bottom coord, only 2:nz
    call ddz_uv(u, dudz, lbz)
    call ddz_uv(v, dvdz, lbz)

    ! Calculate dwdz using finite differences (for 0:nz-1 on w-nodes)
    !  except bottom coord, only 1:nz-1
    call ddz_w(w, dwdz, lbz)
#ifdef PPOUTPUT_CLOCK
    call clock_deriv%stop
#endif
    ! Calculate wall stress and derivatives at the wall
    ! (txz, tyz, dudz, dvdz at jz=1)
    ! MPI: bottom and top processes only
#ifdef PPTLWMLES
    call tlwm_wallstress()
#else
    if (coord == 0) then
#ifdef PPOUTPUT_CLOCK
        call clock_wm%start
#endif
        call wallstress()
#ifdef PPOUTPUT_CLOCK
        call clock_wm%stop
#endif

#ifdef PPCNDIFF
        ! Add boundary condition to explicit portion
        txz_half2(1:nx,:,1) = txz(1:nx,:,1)
        tyz_half2(1:nx,:,1) = tyz(1:nx,:,1)
#endif
    endif
#endif
    if (coord == nproc-1) then
        call wallstress()
#ifdef PPCNDIFF
        ! Add boundary condition to explicit portion
        txz_half2(1:nx,:,nz) = txz(1:nx,:,nz)
        tyz_half2(1:nx,:,nz) = tyz(1:nx,:,nz)
#endif
    endif

    ! Calculate turbulent (subgrid) stress for entire domain
    !   using the model specified in param.f90 (Smag, LASD, etc)
    !   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
#ifdef PPOUTPUT_CLOCK
    call clock_stress%start
#endif
    call sgs_stag()
#ifdef PPOUTPUT_CLOCK
    call clock_stress%stop
#endif

    ! Exchange ghost node information (since coords overlap) for tau_zz
    !   send info up (from nz-1 below to 0 above)
#ifdef PPMPI
    call mpi_sendrecv (tzz(:,:,nz-1), ld*ny, MPI_RPREC, up, 6,             &
                       tzz(:,:,0), ld*ny, MPI_RPREC, down, 6,              &
                       comm, status, ierr)
#endif

    ! Compute divergence of SGS shear stresses
    ! the divt's and the diagonal elements of t are not equivalenced
    ! in this version. Provides divtz 1:nz-1, except 1:nz at top process
#ifdef PPOUTPUT_CLOCK
    call clock_divstress%start
#endif
#ifdef PPCNDIFF
    call divstress_uv(divtx, divty, txx, txy, txz_half1, tyy, tyz_half1)
    call divstress_w_cndiff(divtz, txz, tyz)
#else
    call divstress_uv(divtx, divty, txx, txy, txz, tyy, tyz)
    call divstress_w(divtz, txz, tyz, tzz)
#endif
#ifdef PPOUTPUT_CLOCK
    call clock_divstress%stop
#endif

    ! Calculates u x (omega) term in physical space. Uses 3/2 rule for
    ! dealiasing. Stores this term in RHS (right hand side) variable
#ifdef PPOUTPUT_CLOCK
    call clock_convec%start
#endif
#ifdef PPRNL
    ! Use RNL equations
    call convec(u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy,RHSx,RHSy,RHSz)

    if (.not. fourier) then
        u_pert = u - x_avg(u)
        v_pert = v - x_avg(v)
        w_pert = w - x_avg(w)
        dudy_pert = dudy - x_avg(dudy)
        dudz_pert = dudz - x_avg(dudz)
        dvdx_pert = dvdx - x_avg(dvdx)
        dvdz_pert = dvdz - x_avg(dvdz)
        dwdx_pert = dwdx - x_avg(dwdx)
        dwdy_pert = dwdy - x_avg(dwdy)

        call convec(u_pert, v_pert, w_pert,                                   &
            dudy_pert, dudz_pert, dvdx_pert, dvdz_pert, dwdx_pert, dwdy_pert, &
            RHSx_pert, RHSy_pert, RHSz_pert)

        RHSx = RHSx - RHSx_pert + x_avg(RHSx_pert)
        RHSy = RHSy - RHSy_pert + x_avg(RHSy_pert)
        RHSz = RHSz - RHSz_pert + x_avg(RHSz_pert)
    endif
    ! if fourier and RNL, just run convec since convolution is used
#else
    ! Use full NS equations --> OR running RNL/GQL in Fourier space
    ! OR if use GQL (not Fourier) equations as below
    call convec(u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy,RHSx,RHSy,RHSz) 

#ifdef PPGQL

    if (.not. fourier) then
        ! Use GQL equations
        RHSx = RHSx - gql_filter( RHSx )
        RHSy = RHSy - gql_filter( RHSy )
        RHSz = RHSz - gql_filter( RHSz )

        u_low = gql_filter( u )
        v_low = gql_filter( v )
        w_low = gql_filter( w )
        dudy_low = gql_filter( dudy )
        dudz_low = gql_filter( dudz )
        dvdx_low = gql_filter( dvdx )
        dvdz_low = gql_filter( dvdz )
        dwdx_low = gql_filter( dwdx )
        dwdy_low = gql_filter( dwdy )

        call convec(u_low, v_low, w_low,                                    &
            dudy_low, dudz_low, dvdx_low, dvdz_low, dwdx_low, dwdy_low,     &
            RHSx_low, RHSy_low, RHSz_low)

        RHSx = RHSx - RHSx_low + 2.0_rprec*gql_filter( RHSx_low )
        RHSy = RHSy - RHSy_low + 2.0_rprec*gql_filter( RHSy_low )
        RHSz = RHSz - RHSz_low + 2.0_rprec*gql_filter( RHSz_low )

        u_high = u - u_low
        v_high = v - v_low
        w_high = w - w_low
        dudy_high = dudy - dudy_low
        dudz_high = dudz - dudz_low
        dvdx_high = dvdx - dvdx_low
        dvdz_high = dvdz - dvdz_low
        dwdx_high = dwdx - dwdx_low
        dwdy_high = dwdy - dwdy_low

        call convec(u_high, v_high, w_high,                                   &
            dudy_high, dudz_high, dvdx_high, dvdz_high, dwdx_high, dwdy_high, &
            RHSx_high, RHSy_high, RHSz_high)

        RHSx = RHSx - RHSx_high + 2.0_rprec*gql_filter( RHSx_high )
        RHSy = RHSy - RHSy_high + 2.0_rprec*gql_filter( RHSy_high )
        RHSz = RHSz - RHSz_high + 2.0_rprec*gql_filter( RHSz_high )
    endif

#endif

#endif
#ifdef PPOUTPUT_CLOCK
call clock_convec%stop
#endif

    ! Add div-tau term to RHS variable
    !   this will be used for pressure calculation
    RHSx(:,:,1:nz-1) = -RHSx(:,:,1:nz-1) - divtx(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = -RHSy(:,:,1:nz-1) - divty(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = -RHSz(:,:,1:nz-1) - divtz(:,:,1:nz-1)
    if (coord == nproc-1) RHSz(:,:,nz) = -RHSz(:,:,nz)-divtz(:,:,nz)

    ! Coriolis: add forcing to RHS
    if (coriolis_forcing) then
        ! This is to put in the coriolis forcing using coriol,ug and vg as
        ! precribed in param.f90. (ug,vg) specfies the geostrophic wind vector
        ! Note that ug and vg are non-dimensional (using u_star in param.f90)
        RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) +                           &
            coriol * v(:,:,1:nz-1) - coriol * vg
        RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) -                           &
            coriol * u(:,:,1:nz-1) + coriol * ug
    end if

    !--calculate u^(*) (intermediate vel field)
    !  at this stage, p, dpdx_i are from previous time step
    !  (assumes old dpdx has NOT been added to RHSx_f, etc)
    !  we add force (mean press forcing) here so that u^(*) is as close
    !  to the final velocity as possible
    if (use_mean_p_force) then
        if (fourier) then
            ! only add to the mean (kx=0) mode
            ! no need to transform mean_p_force to kx space
            RHSx(1,1,1:nz-1) = RHSx(1,1,1:nz-1) + mean_p_force_x
            RHSy(1,1,1:nz-1) = RHSy(1,1,1:nz-1) + mean_p_force_y
        else
            RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + mean_p_force_x
            RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + mean_p_force_y
        endif
    end if

    ! Optional random forcing, i.e. to help prevent relaminarization
    if (use_random_force .and. jt_total < stop_random_force) then
        call forcing_random()
    end if

    !//////////////////////////////////////////////////////
    !/// APPLIED FORCES                                 ///
    !//////////////////////////////////////////////////////
    !  In order to save memory the arrays fxa, fya, and fza are now only defined when needed.
    !  For Levelset RNS all three arrays are assigned.
    !  For turbines at the moment only fxa is assigned.
    !  Look in forcing_applied for calculation of forces.
    !  Look in sim_param.f90 for the assignment of the arrays.

    !  Applied forcing (forces are added to RHS{x,y,z})

    ! Calculate forcing time
    !call clock_forcing%start

    ! Apply forcing. These forces will later go into RHS
    call forcing_applied()

    ! Calculate forcing time
    !call clock_forcing%stop

    ! Calculate the total time of the forcing
    !clock_total_f = clock_total_f + clock_forcing % time

    !  Update RHS with applied forcing
#if defined(PPTURBINES) || defined(PPATM)
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)
#endif

    !//////////////////////////////////////////////////////
    !/// EULER INTEGRATION CHECK                        ///
    !//////////////////////////////////////////////////////
    ! Set RHS*_f if necessary (first timestep)
    if ((jt_total == 1) .and. (.not. initu)) then
        ! if initu, then this is read from the initialization file
        ! else for the first step put RHS_f=RHS
        !--i.e. at first step, take an Euler step
        RHSx_f=RHSx
        RHSy_f=RHSy
        RHSz_f=RHSz
    end if

    !//////////////////////////////////////////////////////
    !/// INTERMEDIATE VELOCITY                          ///
    !//////////////////////////////////////////////////////
    ! Calculate intermediate velocity field
    !   only 1:nz-1 are valid
#ifdef PPOUTPUT_CLOCK
    call clock_inter%start
#endif
#ifdef PPCNDIFF
    call diff_stag_array_uv()
    call diff_stag_array_w()
#else
    u(:,:,1:nz-1) = u(:,:,1:nz-1) +                                     &
        dt * ( tadv1 * RHSx(:,:,1:nz-1) + tadv2 * RHSx_f(:,:,1:nz-1) )
    v(:,:,1:nz-1) = v(:,:,1:nz-1) +                                     &
        dt * ( tadv1 * RHSy(:,:,1:nz-1) + tadv2 * RHSy_f(:,:,1:nz-1) )
    w(:,:,1:nz-1) = w(:,:,1:nz-1) +                                     &
        dt * ( tadv1 * RHSz(:,:,1:nz-1) + tadv2 * RHSz_f(:,:,1:nz-1) )
    if (coord == nproc-1) then
        w(:,:,nz) = w(:,:,nz) +                                         &
            dt * ( tadv1 * RHSz(:,:,nz) + tadv2 * RHSz_f(:,:,nz) )
    end if
#endif
#ifdef PPOUTPUT_CLOCK
    call clock_inter%stop
#endif

    ! Set unused values to BOGUS so unintended uses will be noticable
#ifdef PPSAFETYMODE
#ifdef PPMPI
    u(:,:,0) = BOGUS
    v(:,:,0) = BOGUS
    w(:,:,0) = BOGUS
#endif
    u(:,:,nz) = BOGUS
    v(:,:,nz) = BOGUS
    if(coord < nproc-1) w(:,:,nz) = BOGUS
#endif

    !//////////////////////////////////////////////////////
    !/// PRESSURE SOLUTION                              ///
    !//////////////////////////////////////////////////////
    ! Solve Poisson equation for pressure
    !   div of momentum eqn + continuity (div-vel=0) yields Poisson eqn
    !   do not need to store p --> only need gradient
    !   provides p, dpdx, dpdy at 0:nz-1 and dpdz at 1:nz-1
#ifdef PPOUTPUT_CLOCK
    call clock_pres%start
#endif
    call press_stag_array()
#ifdef PPOUTPUT_CLOCK
    call clock_pres%stop
#endif

    ! SIMPLIFIED LEVEL SET METHOD WITH MAPPING CAPABILITY
#ifdef PPLVLSET_STRETCH
    call level_set2(IBFx,IBFy,IBFz)
#endif

    ! Add pressure gradients to RHS variables (for next time step)
    !   could avoid storing pressure gradients - add directly to RHS
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) - dpdx(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) - dpdy(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) - dpdz(:,:,1:nz-1)
    if(coord == nproc-1) then
        RHSz(:,:,nz) = RHSz(:,:,nz) - dpdz(:,:,nz)
    end if

    !//////////////////////////////////////////////////////
    !/// INDUCED FORCES                                 ///
    !//////////////////////////////////////////////////////
    ! Calculate external forces induced forces. These are
    ! stored in fx,fy,fz arrays. We are calling induced
    ! forces before applied forces as some of the applied
    ! forces (RNS) depend on the induced forces and the
    ! two are assumed independent
    call forcing_induced()

    !//////////////////////////////////////////////////////
    !/// PROJECTION STEP                                ///
    !//////////////////////////////////////////////////////
    ! Projection method provides u,v,w for jz=1:nz
    !   uses fx,fy,fz calculated above
    !   for MPI: syncs 1 -> Nz and Nz-1 -> 0 nodes info for u,v,w
    call project ()

#ifdef PPLVLSET
    if (global_CA_calc) call level_set_global_CA()
#endif

    ! Write output files
    call output_loop()

    ! Check the total time of the simulation up to this point on the master
    ! node and send this to all

    if (modulo (jt_total, wbase) == 0) then

        ! Get the ending time for the iteration
        call clock%stop
        call clock_total%stop

        if (fourier) then
            call wave2physF( u, uF )
            call wave2physF( v, vF )
            call wave2physF( w, wF )
            call wave2physF( txz, txzF )
            call wave2physF( tyz, tyzF )
        endif

        ! Calculate rms divergence of velocity
        ! only written to screen, not used otherwise
        call rmsdiv(rmsdivvel)
        maxcfl = get_max_cfl()
        maxvisc = get_max_visc()

        ! This takes care of the clock times, to obtain the quantities based
        ! on all the processors, not just processor 0
#ifdef PPMPI
        call mpi_allreduce(clock % time, maxdummy, 1, mpi_rprec,            &
            MPI_MAX, comm, ierr)
        clock % time = maxdummy
        call mpi_allreduce(clock_total % time, maxdummy, 1, mpi_rprec,      &
            MPI_MAX, comm, ierr)
        clock_total % time = maxdummy
        !call mpi_allreduce(clock_forcing % time, maxdummy, 1, mpi_rprec,   &
        !    MPI_MAX, comm, ierr)
        !clock_forcing % time = maxdummy
        !call mpi_allreduce(clock_total_f , maxdummy, 1, mpi_rprec,         &
        !    MPI_MAX, comm, ierr)
        !clock_total_f = maxdummy
#endif

        ! Send top wall stress to bottom process
#ifdef PPMPI
        if (ubc_mom >0) then
            if (coord == nproc-1) then
                tau_top = get_tau_wall_top()
            else
                tau_top = 0._rprec
            endif

            call mpi_allreduce(tau_top, maxdummy, 1, mpi_rprec,             &
                MPI_SUM, comm, ierr)
            tau_top = maxdummy
        endif
#endif

        if (coord == 0) then
            write(*,*)
            write(*,'(a)') '========================================================'
            write(*,'(a)') 'Time step information:'
            write(*,'(a,i9)') '  Iteration: ', jt_total
            write(*,'(a,E15.7)') '  Time step: ', dt
            write(*,'(a,E15.7)') '  CFL: ', maxcfl
            write(*,'(a,E15.7)') '  VISC: ', maxvisc
            ! write(*,'(a,2E15.7)') '  AB2 TADV1, TADV2: ', tadv1, tadv2
            write(*,*)
            write(*,'(a)') 'Flow field information:'
            write(*,'(a,E15.7)') '  Velocity divergence metric: ', rmsdivvel
            write(*,'(a,E15.7)') '  Bot wall stress: ', get_tau_wall_bot()
            write(*,'(a,E15.7)') '  Turnovers: ', total_time_dim / ( L_x * z_i / u_star )
#ifdef PPMPI
            if (ubc_mom > 0) then
                write(*,'(a,E15.7)') '  Top wall stress: ', tau_top
            end if
#else
            if (ubc_mom > 0) then
                write(*,'(a,E15.7)') '  Top wall stress: ', get_tau_wall_top()
            end if
#endif
            write(*,*)
            write(*,'(1a)') 'Simulation wall times (s): '
            write(*,'(1a,E15.7)') '  Iteration: ', clock % time
            write(*,'(1a,E15.7)') '  Cumulative: ', clock_total % time
            ! write(*,'(1a,E15.7)') '  Forcing: ', clock_forcing % time
            ! write(*,'(1a,E15.7)') '  Cumulative Forcing: ', clock_total_f
            ! write(*,'(1a,E15.7)') '  Forcing %: ',                        &
            !     clock_total_f /clock_total % time
#ifdef PPOUTPUT_WMLES
            write(*,'(1a,E15.7)') '  WM Iteration: ', clock_wm % time
#endif
            write(*,'(a)') '========================================================'
            call write_tau_wall_bot()
        end if
        if ((coord == nproc-1) .AND. (ubc_mom /= 0)) then
            call write_tau_wall_top()
        end if

        if (fourier) then
            if(coord == 0) then
                write(*,'(a)') '======================= BOTTOM ========================='
                write(*,*) 'u: ', uF(nxp/2,ny/2,1:2)
                write(*,*) 'v: ', vF(nxp/2,ny/2,1:2)
                write(*,*) 'w: ', wF(nxp/2,ny/2,1:2)
                write(*,'(a)') '========================================================'
            end if
            call mpi_barrier(comm, ierr)
            if(coord == nproc-1) then
                write(*,'(a)') '======================== TOP ==========================='
                write(*,*) 'u: ', uF(nxp/2,ny/2,nz-2:nz-1)
                write(*,*) 'v: ', vF(nxp/2,ny/2,nz-2:nz-1)
                write(*,*) 'w: ', wF(nxp/2,ny/2,nz-1:nz)
                write(*,'(a)') '========================================================'
            end if
            call mpi_barrier(comm, ierr)
        else
            if(coord == 0) then
                write(*,'(a)') '======================= BOTTOM ========================='
                write(*,*) 'u: ', u(nx/2,ny/2,1:2)
                write(*,*) 'v: ', v(nx/2,ny/2,1:2)
                write(*,*) 'w: ', w(nx/2,ny/2,1:2)
                write(*,'(a)') '========================================================'
            end if
            call mpi_barrier(comm, ierr)
            if(coord == nproc-1) then
                write(*,'(a)') '======================== TOP ==========================='
                write(*,*) 'u: ', u(nx/2,ny/2,nz-2:nz-1)
                write(*,*) 'v: ', v(nx/2,ny/2,nz-2:nz-1)
                write(*,*) 'w: ', w(nx/2,ny/2,nz-1:nz)
                write(*,'(a)') '========================================================'
            end if
            call mpi_barrier(comm, ierr)
        endif

        if (.not. fourier) then
            call kx_energy()
            !call ky_energy()
        else
            call kx_energy_fourier()
        endif

#ifdef PPOUTPUT_CLOCK
#ifdef PPMPI
        call mpi_allreduce(clock_deriv % time, maxdummy, 1, mpi_rprec,     &
            MPI_MAX, comm, ierr)
        clock_deriv % time = maxdummy
        call mpi_allreduce(clock_stress % time, maxdummy, 1, mpi_rprec,    &
            MPI_MAX, comm, ierr)
        clock_stress % time = maxdummy
        call mpi_allreduce(clock_divstress % time, maxdummy, 1, mpi_rprec, &
            MPI_MAX, comm, ierr)
        clock_divstress % time = maxdummy
        call mpi_allreduce(clock_convec % time, maxdummy, 1, mpi_rprec,    &
            MPI_MAX, comm, ierr)
        clock_convec % time = maxdummy
        call mpi_allreduce(clock_inter % time, maxdummy, 1, mpi_rprec,     &
            MPI_MAX, comm, ierr)
        clock_inter % time = maxdummy
        call mpi_allreduce(clock_pres % time, maxdummy, 1, mpi_rprec,      &
            MPI_MAX, comm, ierr)
        clock_pres % time = maxdummy
        ! Remember clock_wm is only on coord == 0
#endif
if (coord == 0) call write_clocks( clock%time, clock_deriv%time,           &
    clock_stress%time, clock_wm%time, clock_divstress%time,                &
    clock_convec%time, clock_inter%time,  clock_pres%time)
#endif

        ! Check if we are to check the allowable runtime
        if (runtime > 0) then

#ifdef PPMPI
            ! Determine the processor that has used most time and communicate
            ! this. Needed to make sure that all processors stop at the same
            ! time and not just some of them
            call mpi_allreduce(clock_total % time, rbuffer, 1, MPI_RPREC,   &
                MPI_MAX, MPI_COMM_WORLD, ierr)
            clock_total % time = rbuffer
#endif

            ! If maximum time is surpassed go to the end of the program
            if ( clock_total % time >= real(runtime,rprec) ) then
                call mesg( prog_name,                                       &
                    'Specified runtime exceeded. Exiting simulation.')
                exit time_loop
            endif

       endif

    ! Shift the domain in the y (spanwise) direction
#ifdef PPSTREAKS
        if (modulo (jt_total, 1000) == 0) then
            if (coord == 0) then
                write(*,*) '--------------------red shift-------------------'
            endif
            dummyu(:,:,:) = u(:,:,:)
            dummyv(:,:,:) = v(:,:,:)
            dummyw(:,:,:) = w(:,:,:)

            dummyRHSx(:,:,:) = RHSx(:,:,:)
            dummyRHSy(:,:,:) = RHSy(:,:,:)
            dummyRHSz(:,:,:) = RHSz(:,:,:)

            u(:,2:ny,:) = dummyu(:,1:ny-1,:)
            v(:,2:ny,:) = dummyv(:,1:ny-1,:)
            w(:,2:ny,:) = dummyw(:,1:ny-1,:)

            RHSx(:,2:ny,:) = dummyRHSx(:,1:ny-1,:)
            RHSy(:,2:ny,:) = dummyRHSy(:,1:ny-1,:)
            RHSz(:,2:ny,:) = dummyRHSz(:,1:ny-1,:)

            u(:,1,:) = dummyu(:,ny,:)
            v(:,1,:) = dummyv(:,ny,:)
            w(:,1,:) = dummyw(:,ny,:)

            RHSx(:,1,:) = dummyRHSx(:,ny,:)
            RHSy(:,1,:) = dummyRHSy(:,ny,:)
            RHSz(:,1,:) = dummyRHSz(:,ny,:)

            dummyu(:,:,:) = F_LM(:,:,:)
            F_LM(:,2:ny,:) = dummyu(:,1:ny-1,:)
            F_LM(:,1,:) = dummyu(:,ny   ,:)

            dummyu(:,:,:) = F_MM(:,:,:)
            F_MM(:,2:ny,:) = dummyu(:,1:ny-1,:)
            F_MM(:,1,:) = dummyu(:,ny   ,:)

            dummyu(:,:,:) = F_QN(:,:,:)
            F_QN(:,2:ny,:) = dummyu(:,1:ny-1,:)
            F_QN(:,1,:) = dummyu(:,ny   ,:)

            dummyu(:,:,:) = F_NN(:,:,:)
            F_NN(:,2:ny,:) = dummyu(:,1:ny-1,:)
            F_NN(:,1,:) = dummyu(:,ny   ,:)
        endif
#endif

    end if

end do time_loop
! END TIME LOOP

! Finalize
close(2)

! Write total_time.dat and tavg files
call output_final()

! Stop wall clock
call clock_total%stop
#ifdef PPMPI
if (coord == 0) write(*,"(a,e15.7)") 'Simulation wall time (s) : ',      &
    clock_total % time
#else
if (coord == 0) write(*,"(a,e15.7)") 'Simulation cpu time (s) : ',       &
    clock_total % time
#endif

call finalize()

if(coord == 0 ) write(*,'(a)') 'Simulation complete'

end program main
