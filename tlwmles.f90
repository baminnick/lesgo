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

!*******************************************************************************
module tlwmles
!*******************************************************************************
! 
! This module contains all functions necessary for solving a set of equations in
! the inner layer as a wall model for the LES.
! 

use types, only : rprec

implicit none 

save
private
#ifdef PPOUTPUT_WMLES
public tlwm, tlwm_init, tlwm_finalize, nzr, hwm
#else
public tlwm, tlwm_init, tlwm_finalize
#endif

! Wall-model settings to be defined by the user:
! Inputs to specify grid resolution and spacing
real(rprec) :: dzrp = 0.4_rprec !! size of first cell from wall, in wall units
! dzrp <= 0.5 to be wall-resolved
!integer :: ihwm = 4 !! integer value for wall-model height
! if ihwm = 1, then using first LES uv grid point so hwm = dz/2
real(rprec) :: ss = 1.2_rprec !! grid stretching parameter
! if ss = 1, then the grid is uniform, s can only be s >= 1

! Inputs for wall-model eddy viscosity
integer :: ktype = 1 !! how the wall-model coefficient kwm should be specified
! ktype = 1 uses kwm = vonk, ktype = 2 adjusts kwm by nu_sgs
! ktype = 3 uses a combination of 1 and 2 and varies in wall-normal

! Inputs for nonequilibrium wall-model
integer :: ntwm = 1 !! integer value for number of wall-model timesteps from LES
! if ntwm = 1, then wall-model and LES advance together
! if ntwm = 2, then wall-model changes in 2 LES timesteps

! Inputs for 2d3c/rnl equilibrium wall-model
integer, dimension(2) :: kxi = (/1, 6/) !! set of kx mode indices
! Note kxi = 1 refers to the zero mode and should always be retained

! Variables to be defined in tlwm_init and used throughout
real(rprec), dimension(:,:,:), allocatable :: ur, vr, nu_wm
real(rprec), dimension(:,:,:), allocatable :: wr, durdx, durdy, dvrdx, dvrdy
real(rprec), dimension(:), allocatable :: zr
real(rprec) :: dzr, hwm
integer :: nzr, nvel

! Variables to be defined in tlwm_init and used for tlwm_fourier
logical :: tlwm_fourier
real(rprec), dimension(:,:), allocatable :: txzf, tyzf, dudzf, dvdzf
integer :: nxf, nyf

contains

!*******************************************************************************
subroutine tlwm
!*******************************************************************************
!
! This subroutine acts as the main subroutine for all two-layer wall models. It
! calls the wall model that is to be used to give the wall stress, this 
! subroutine also calls functions that output data and statistics.
!
use param, only : lbc_mom
use param, only : checkpoint_data, checkpoint_nskip, jt_total
#ifdef PPOUTPUT_WMLES
use param, only : dt, tavg_nstart, tavg_nend, tavg_nskip, tavg_calc
use param, only : domain_calc, domain_nstart, domain_nskip, domain_nend
use stat_defs, only : tavg_initialized, tavg_wmles_dt
#endif
implicit none

! Calculate wall stress
if (lbc_mom == 5) then
    call ws_tlwm_eq_lbc()
elseif (lbc_mom == 6) then
    call ws_tlwm_noneq_lbc()
elseif (lbc_mom == 7) then
    call ws_rnl_eq_lbc()
    ! call ws_2d3cF_eq_lbc()
elseif (lbc_mom == 8) then
    ! call ws_rnl_noneq_lbc()
end if

! Determine if we are to checkpoint intermediate times
if (checkpoint_data) then
    ! Now check if data should be checkpointed this time step
    if ( modulo (jt_total, checkpoint_nskip ) == 0) call tlwm_checkpoint()
end if

#ifdef PPOUTPUT_WMLES
! OUTPUT LOOP (for TLWMLES)

! Instantaneous Domain Velocities
if (domain_calc) then
    if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and.  &
        (mod(jt_total-domain_nstart,domain_nskip)==0) ) then
        call wmles_inst_write()
    endif
endif

! Time-averaging statistics
! Determine if time summations are to be calculated
if (tavg_calc) then
    ! Are we between the start and stop timesteps?
    if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then
        ! Every timestep (between nstart and nend), add to tavg_wmles_dt
        tavg_wmles_dt = tavg_wmles_dt + dt

        ! Are we at the beginning or a multiple of nstart?
        if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then
            ! Check if we have initialized tavg
            if (.not.tavg_initialized) then
                call tavg_wmles_init()
            else
                call tavg_wmles_compute ()
            end if
        end if
    end if
end if
#endif

end subroutine tlwm

!*******************************************************************************
subroutine tlwm_init
!*******************************************************************************
!
! This subroutine allocates memory used for tlwm and sets up the grid used. 
!
 
use types, only : rprec
use param, only : ld, nx, ny, lbc_mom, path, read_endian
use param, only : tlwm_kxnum
use sim_param, only : txz, tyz
#ifdef PPOUTPUT_WMLES
use stat_defs, only : tavg_wmles
use param, only : tavg_calc
#endif

implicit none

character(64) :: fname
logical :: file_flag
integer :: k
#ifdef PPOUTPUT_WMLES
integer :: i, j
#endif

! Determine if using Fourier mode or not
! Also specify RNL TLWMLES parameters here
tlwm_fourier = .false.
if (lbc_mom == 7) then 
    tlwm_fourier = .true.
    nxf = 2*(tlwm_kxnum - 1)
endif

! Check if checkpoint file exists
fname = path // 'wmles.out'
inquire( file=fname, exist=file_flag )

if (file_flag) then

    ! Read in checkpoint data
    open(12, file=fname, form='unformatted', convert=read_endian)
    read(12) hwm, dzr, nzr, ss, txz(1:nx,1:ny,1), tyz(1:nx,1:ny,1)
    close(12)

    ! Allocate grid variables
    allocate(zr(nzr))

    ! Find number of velocity points on inner-layer grid
    nvel = (nzr+1)/2

    ! Setup inner layer grid as before
    zr(nzr) = 0.0_rprec
    do k = (nzr-1), 1, -1
        zr(k) = zr(k+1) + (ss**(real(nzr-1-k,rprec)))*dzr
    end do

else 
    ! Create inner layer grid
    call tlwm_grid()
end if

! Allocate variables
if (.not. tlwm_fourier) then
    allocate(ur(ld,ny,nvel))
    allocate(vr(ld,ny,nvel))
    allocate(nu_wm(ld,ny,nzr))
else !! fourier mode
    allocate(txzf(nxf+2,ny))
    allocate(tyzf(nxf+2,ny))
    allocate(dudzf(nxf+2,ny))
    allocate(dvdzf(nxf+2,ny))
    txzf(:,:) = 0.0_rprec
    tyzf(:,:) = 0.0_rprec
    dudzf(:,:) = 0.0_rprec
    dvdzf(:,:) = 0.0_rprec
#ifdef PPOUTPUT_WMLES
    allocate(ur(ld,ny,nvel))
    allocate(vr(ld,ny,nvel))
    allocate(nu_wm(ld,ny,nzr))
#endif
endif

! Variables only used in Non-equilibrium model
if (lbc_mom == 6) then 
    allocate(wr(ld,ny,nvel))
    allocate(durdx(ld,ny,nvel))
    allocate(durdy(ld,ny,nvel))
    allocate(dvrdx(ld,ny,nvel))
    allocate(dvrdy(ld,ny,nvel))
end if

! Initialize time-averaged quantites
#ifdef PPOUTPUT_WMLES
if (tavg_calc) then
    allocate(tavg_wmles(1:nx,1:ny,nvel))

    do j = 1, ny
    do i = 1, nx
        do k = 1, nvel
        tavg_wmles(i,j,k) % u = 0._rprec
        tavg_wmles(i,j,k) % v = 0._rprec
        tavg_wmles(i,j,k) % nu = 0._rprec

        tavg_wmles(i,j,k) % uu = 0._rprec
        tavg_wmles(i,j,k) % vv = 0._rprec
        tavg_wmles(i,j,k) % uv = 0._rprec
        end do
    end do
    end do
endif

write(*,*) '--> tlwmles initialized.'
#endif

end subroutine tlwm_init

!*******************************************************************************
subroutine tlwm_grid()
!*******************************************************************************
!
! This subroutine creates the inner layer grid for the tlwm. Currently this is
! done based on the user supplied u_star and nu_molec values. The grid is 
! proportionally stretched by the parameter ss which is defined in the beginning 
! of module tlwmles.
! 
! The code below determines how many points are in the inner layer, then
! corrects ss based on this value. To correct ss, the bisection method is 
! used since a proportionally stretched grid is given by a geometric series.
!

use types, only : rprec
use param, only : u_star, nu_molec, dz, nz, ihwm
#ifdef PPOUTPUT_WMLES
use param, only : dx, dy
#endif

implicit none

integer :: k
real(rprec) :: grid_tol, bia, bib, fs !! parameters to create the grid

! Correct ihwm if needed
if (ihwm > nz) then !! User asked for grid point outside first processor
    ihwm = nz
#ifdef PPOUTPUT_WMLES
    write(*,*) 'TLWM: Changing ihwm to: ', ihwm
#endif
end if

! Determine height of wall-model from ihwm
hwm = real(ihwm-1,rprec)*dz + dz/(2.0_rprec)

! Find physical height of first wall-model cell from wall using dzrp
dzr = dzrp*nu_molec/u_star

! Determine number of inner layer grid points, using stretching parameter
if (ss == 1.0_rprec) then !! uniform grid
    nzr = int( (hwm/dzr) + 1.0_rprec ) !! rounding down
else !! stretched grid
    nzr = int( 1.0_rprec + (1.0_rprec/log(ss))*                          &
        (log(1.0_rprec + ((hwm/dzr)*(ss-1.0_rprec)))) )
end if

! Correct nzr if needed
! It has to be odd so the U-grid lies on the interface and on the wall
! Also must be at least 7 points for the TDMA to work
if (mod(nzr,2) == 0) then !! nzr is even
    nzr = nzr - 1
end if
if (nzr < 7) then
    nzr = 7
end if

! Recompute stretching parameter based on nzr
grid_tol = 10.0_rprec**(-8)
bia = ss !! initial interval bound 
bib = ss + 0.15_rprec !! initial interval bound
fs = 1.0_rprec !! initializing function evaluation, should be == 0
if ( (hwm - real(nzr-1,rprec)*dzr) >= (grid_tol) ) then 
    ! grid is not uniform so ss ~= 1
    ! Using bisection method to find correct ss
    do while (abs(fs) > grid_tol)
        ss = (bia + bib)*0.5_rprec !! find midpoint
        fs = ((1.0_rprec-(ss**real(nzr-1,rprec)))/                     &
            (1.0_rprec-ss))*dzr - hwm !! == f(ss) = 0
        if (fs > 0) then
            bib = ss
        else !! fs < 0
            bia = ss
        end if
    end do
end if

! Check if ktype = 3 can be used
#ifdef PPOUTPUT_WMLES
if ((hwm < max(dx,dy)) .and. (ktype == 3)) then
    write(*,*) 'hwm<(dx or dy): May change wm coef from dynamic to constant!'
end if
#endif

! Allocate grid variables
allocate(zr(nzr))

! Setup inner layer grid
! point 1 on grid is ihwm'th point on LES grid away from the wall
! point nzr on grid is on the wall
zr(nzr) = 0.0_rprec
do k = (nzr-1), 1, -1
    zr(k) = zr(k+1) + (ss**(real(nzr-1-k,rprec)))*dzr
end do
! zr(1) should equal hwm

! Number of points which wall-model velocity lies on, nvel
nvel = (nzr+1)/2

end subroutine tlwm_grid

!*******************************************************************************
subroutine ws_tlwm_eq_lbc()
!*******************************************************************************
! 
! This subroutine solves a 2nd order ODE which assumes equilibrium and a 
! eddy viscosity model that satisfies the log-law and the viscous sublayer if
! the LES grid is sufficiently resolved.
! 
! Currently the model assumes smooth walls, finite Reynolds number, and takes 
! the inner layer starting at the ihwm'th LES grid point from the wall (on the 
! uv grid).
!  
! The ODE is of the form:
!     d/dz[(nu+nu_t)*dUdz] = 0,  U(z=0) = 0,  U(z=hwm) = Ules
! Centered differencing is used. As a result, the quantity (nu+nu_t) is
! evaluated between grid points which U is solved on. The TDMA algorithm is 
! used to invert the resulting matrix equation.
! 

use param, only : ld, nx, ny, nu_molec, z_i, ubot, u_star, jt_total, ihwm
use sim_param, only : u, v, dudz, dvdz, txz, tyz
use test_filtermodule, only : test_filter
implicit none
integer :: i, j, k, cnt
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec), dimension(nx, ny) :: ustar
real(rprec) :: nu1, nu2
real(rprec), dimension(nvel-2) :: Aa, Ab, Ac, Ad, Bd, ur1, vr1 !! for TDMA

! Treatment of LES velocity into wall-model BC
! Average so free stream velocities fluctuate less and obey log-law
u1 = u(:,:,ihwm)
v1 = v(:,:,ihwm)

! Use test filter, see appendix of E Bou-Zeid and Meneveau
call test_filter(u1)
call test_filter(v1)

! Find constant ustar for wall model eddy viscosity
! 1. Take the average of previous timesteps tau value
if (jt_total > 1) then
    ustar = sqrt(abs(txz(1:nx,1:ny,1)) + abs(tyz(1:nx,1:ny,1)))
else !! first time-step, should use nonzero ustar
    ustar = u_star
end if
! 2. Use constant value specified by user - only for debugging
!ustar = u_star 

! Compute eddy viscosity values on inner layer grid
call wm_eddyvisc(ustar)

do j = 1, ny
do i = 1, nx

    ! 1st equation in Matrix, k = 3
    nu1 = (nu_wm(i,j,2) + nu_molec)/((zr(4)-zr(2))*(zr(3)-zr(1)))
    nu2 = (nu_wm(i,j,4) + nu_molec)/((zr(4)-zr(2))*(zr(5)-zr(3)))
    Aa(1) = 0.0_rprec
    Ab(1) = -(nu1 + nu2)
    Ac(1) = nu2

    ! BC from LES
    Ad(1) = -nu1*u1(i,j)
    Bd(1) = -nu1*v1(i,j)

    ! Middle equations in Matrix, k = 5, 7, ..., nzr-4
    cnt = 2
    do k = 5, nzr-4, 2
        nu1 = (nu_wm(i,j,k-1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k)-zr(k-2)))
        nu2 = (nu_wm(i,j,k+1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k+2)-zr(k)))
        Aa(cnt) = nu1
        Ab(cnt) = -(nu1 + nu2)
        Ac(cnt) = nu2
        Ad(cnt) = 0.0_rprec
        Bd(cnt) = 0.0_rprec
        cnt = cnt + 1
    end do

    ! Last equation in Matrix, k = nzr-2
    nu1 = (nu_wm(i,j,nzr-3) + nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr-2)-zr(nzr-4)))
    nu2 = (nu_wm(i,j,nzr-1) + nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr)-zr(nzr-2)))
    Aa((nzr-3)/2) = nu1
    Ab((nzr-3)/2) = -(nu1 + nu2)
    Ac((nzr-3)/2) = 0.0_rprec

    ! no-slip BC 
    Ad((nzr-3)/2) = -nu2*ubot
    Bd((nzr-3)/2) = 0.0_rprec

    ! Solve matrix equation
    call tdma(Aa,Ab,Ac,Ad,ur1)
    call tdma(Aa,Ab,Ac,Bd,vr1)

    ur(i,j,1) = u1(i,j)
    ur(i,j,2:(nvel-1)) = ur1
    ur(i,j,nvel) = ubot
    vr(i,j,1) = v1(i,j)
    vr(i,j,2:(nvel-1)) = vr1
    vr(i,j,nvel) = 0.0_rprec

    ! Output Boundary Conditions for LES
    ! Assuming inner layer is sufficiently resolved down to the wall so that
    ! the closest point within the inner layer is in the viscous sublayer
    ! Using one-sided finite difference
    dudz(i,j,1) = ( ur(i,j,nvel-1) - ubot ) / (zr(nzr-2) - zr(nzr))
    dvdz(i,j,1) = vr(i,j,nvel-1) / (zr(nzr-2) - zr(nzr))
    txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
    tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)

end do
end do

end subroutine ws_tlwm_eq_lbc

!*******************************************************************************
subroutine ws_tlwm_noneq_lbc()
!*******************************************************************************
! 
! 
! 
! TO DO:
! - Should the time-derivative have a different dt?
! - Derivatives not in conservative form
! - Using first order BE in time, improve to CN
! - Output new time-averaging statistics
! - Consider treating advection terms similar to convec.f90
! 

use param, only : ld, nx, ny, nu_molec, z_i, ubot, u_star, dt, jt_total, ihwm
use sim_param, only : u, v, w, dudz, dvdz, txz, tyz, dpdx, dpdy
use sim_param, only : dudx, dvdx, dwdx, dudy, dvdy, dwdy
use test_filtermodule, only : test_filter
implicit none
integer :: i, j, k, cnt
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec), dimension(ld, ny, nvel) :: durdx, durdy, dvrdx, dvrdy, dpdxr, dpdyr
real(rprec), dimension(nx, ny, nvel) :: advux, advuy, advvx, advvy,            &
    advuz, advvz, RHSur, RHSvr
real(rprec), dimension(nx, ny) :: ustar
real(rprec) :: nu1, nu2
real(rprec), dimension(nvel-2) :: Aa, Ab, Ac, Ad, Bd, ur1, vr1 !! for TDMA

! Treatment of LES velocity into wall-model BC
! Average so free stream velocities fluctuate less and obey log-law
u1 = u(:,:,ihwm)
v1 = v(:,:,ihwm)

! Use test filter, see appendix of E Bou-Zeid and Meneveau
call test_filter(u1)
call test_filter(v1)

! Find constant ustar for wall model eddy viscosity
! 1. Take the average of previous timesteps tau value
if (jt_total > 1) then
    ustar = sqrt(abs(txz(1:nx,1:ny,1)) + abs(tyz(1:nx,1:ny,1)))
else !! first time-step, should use nonzero ustar
    ustar = u_star
end if
! 2. Use constant value specified by user - only for debugging
!ustar = u_star 

! Compute eddy viscosity values on inner layer grid
call wm_eddyvisc(ustar)

! Compute velocity derivatives from previous time-step
call tlwm_ddxy(ur,durdx,durdy)
call tlwm_ddxy(vr,dvrdx,dvrdy)

! Find wall-normal velocity using continuity
! Starting from bottom to impose the boundary conditions
wr(1:nx,1:ny,nvel-1:nvel) = 0.0_rprec !! no penetration and no-slip
cnt = nzr !! momentary index for zr since it is indexed differently than wr
do k = (nvel-1), 2, -1
    wr(1:nx,1:ny,k-1) = wr(1:nx,1:ny,k+1) -                     &
        (zr(cnt-4) - zr(cnt))*(durdx(1:nx,1:ny,k) + dvrdy(1:nx,1:ny,k))
    cnt = cnt - 2
end do

! Advection terms that don't involve wall-normal derivative
! NOTE: not in conservative form
advux(1:nx,1:ny,1:nvel) = ur(1:nx,1:ny,1:nvel)*durdx(1:nx,1:ny,1:nvel)
advuy(1:nx,1:ny,1:nvel) = vr(1:nx,1:ny,1:nvel)*durdy(1:nx,1:ny,1:nvel)
advvx(1:nx,1:ny,1:nvel) = ur(1:nx,1:ny,1:nvel)*dvrdx(1:nx,1:ny,1:nvel)
advvy(1:nx,1:ny,1:nvel) = vr(1:nx,1:ny,1:nvel)*dvrdy(1:nx,1:ny,1:nvel)

! Advection term involving wall-normal derivative
! NOTE: in conservative form
! Using first-order forward finite-difference for first point
advuz(1:nx,1:ny,1) = (1.0_rprec/(zr(1) - zr(3)))*                       &
    ( wr(1:nx,1:ny,1)*ur(1:nx,1:ny,1) - wr(1:nx,1:ny,2)*ur(1:nx,1:ny,2) )
advvz(1:nx,1:ny,1) = (1.0_rprec/(zr(1) - zr(3)))*                       &
    ( wr(1:nx,1:ny,1)*vr(1:nx,1:ny,1) - wr(1:nx,1:ny,2)*vr(1:nx,1:ny,2) )
! Using second-order central finite-difference 
cnt = 1
do k = 2, (nvel-1) !! loop through interior points 
    advuz(1:nx,1:ny,k) = (1.0_rprec/(zr(cnt)-zr(cnt+4)))*               &
        ( wr(1:nx,1:ny,k-1)*ur(1:nx,1:ny,k-1) -                         &
        wr(1:nx,1:ny,k+1)*ur(1:nx,1:ny,k+1) )
    advvz(1:nx,1:ny,k) = (1.0_rprec/(zr(cnt)-zr(cnt+4)))*               &
        ( wr(1:nx,1:ny,k-1)*vr(1:nx,1:ny,k-1) -                         &
        wr(1:nx,1:ny,k+1)*vr(1:nx,1:ny,k+1) )
    cnt = cnt + 2
end do
! no-slip and no-penetration
advuz(1:nx,1:ny,nvel) = 0.0_rprec
advvz(1:nx,1:ny,nvel) = 0.0_rprec

! Pressure Gradient
! Removing energy from pressure
! Interpolating w and (dwdx/dwdy) to the uv grid only at point ihwm
dpdxr(1:ld,1:ny,1) = dpdx(1:ld,1:ny,ihwm) -                         &
    ( u(1:ld,1:ny,ihwm)*dudx(1:ld,1:ny,ihwm)                        &
    + v(1:ld,1:ny,ihwm)*dvdx(1:ld,1:ny,ihwm)                        &
    + 0.5_rprec*( w(1:ld,1:ny,ihwm+1)*dwdx(1:ld,1:ny,ihwm+1) +      &
    w(1:ld,1:ny,ihwm-1)*dwdx(1:ld,1:ny,ihwm-1) ) )
dpdyr(1:ld,1:ny,1) = dpdy(1:ld,1:ny,ihwm) -                         &
    ( u(1:ld,1:ny,ihwm)*dudy(1:ld,1:ny,ihwm)                        &
    + v(1:ld,1:ny,ihwm)*dvdy(1:ld,1:ny,ihwm)                        &
    + 0.5_rprec*( w(1:ld,1:ny,ihwm+1)*dwdy(1:ld,1:ny,ihwm+1) +      &
    w(1:ld,1:ny,ihwm-1)*dwdy(1:ld,1:ny,ihwm-1) ) )
! Test-filter for smoothing
call test_filter(dpdxr(1:ld,1:ny,1))
call test_filter(dpdyr(1:ld,1:ny,1))
! Apply constant pressure gradient across wall-normal direction
do k = 2, nvel
    dpdxr(1:nx,1:ny,k) = dpdxr(1:nx,1:ny,1)
    dpdyr(1:nx,1:ny,k) = dpdyr(1:nx,1:ny,1)
end do

! Move nonequilibrium terms (excluding time derivative) to the RHS
RHSur(1:nx,1:ny,1:nvel) = -( advux(1:nx,1:ny,1:nvel) +           &
    advuy(1:nx,1:ny,1:nvel) + advuz(1:nx,1:ny,1:nvel) +          &
    dpdxr(1:nx,1:ny,1:nvel) )
RHSvr(1:nx,1:ny,1:nvel) = -( advvx(1:nx,1:ny,1:nvel) +           &
    advvy(1:nx,1:ny,1:nvel) + advvz(1:nx,1:ny,1:nvel) +          &
    dpdyr(1:nx,1:ny,1:nvel) )

!RHSur(1:nx,1:ny,1:nvel) = -( advux(1:nx,1:ny,1:nvel) +           &
!    advuy(1:nx,1:ny,1:nvel) + advuz(1:nx,1:ny,1:nvel) )
!RHSvr(1:nx,1:ny,1:nvel) = -( advvx(1:nx,1:ny,1:nvel) +           &
!    advvy(1:nx,1:ny,1:nvel) + advvz(1:nx,1:ny,1:nvel) )

!RHSur(1:nx,1:ny,1:nvel) = -( dpdxr(1:nx,1:ny,1:nvel) )
!RHSvr(1:nx,1:ny,1:nvel) = -( dpdyr(1:nx,1:ny,1:nvel) )

!RHSur(1:nx,1:ny,1:nvel) = 0.0_rprec
!RHSvr(1:nx,1:ny,1:nvel) = 0.0_rprec

!DEBUG
!write(*,*) 'ux', advux(15,15,:)
!write(*,*) '------------------------------'
!write(*,*) 'uy', advuy(15,15,:)
!write(*,*) '------------------------------'
!write(*,*) 'uz', advuz(15,15,:)
!write(*,*) '------------------------------'
write(*,*) 'RHS', RHSur(15,15,:)
write(*,*) '------------------------------'

! Diffusive term is moved to the LHS and treated implicitly
do j = 1, ny
do i = 1, nx

    ! 1st equation in Matrix, k = 3
    nu1 = (nu_wm(i,j,2) + nu_molec)/((zr(4)-zr(2))*(zr(3)-zr(1)))
    nu2 = (nu_wm(i,j,4) + nu_molec)/((zr(4)-zr(2))*(zr(5)-zr(3)))
    Aa(1) = 0.0_rprec
    Ab(1) = (nu1 + nu2) + (1.0_rprec/dt)
    Ac(1) = -nu2

    ! BC from LES
    Ad(1) = nu1*u1(i,j) + ((1.0_rprec/dt)*ur(i,j,2)) + RHSur(i,j,2)
    Bd(1) = nu1*v1(i,j) + ((1.0_rprec/dt)*vr(i,j,2)) + RHSvr(i,j,2)

    ! Middle equations in Matrix, k = 5, 7, ..., nzr-4
    cnt = 2
    do k = 5, nzr-4, 2
        nu1 = (nu_wm(i,j,k-1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k)-zr(k-2)))
        nu2 = (nu_wm(i,j,k+1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k+2)-zr(k)))
        Aa(cnt) = -nu1
        Ab(cnt) = (nu1 + nu2) + (1.0_rprec/dt)
        Ac(cnt) = -nu2
        Ad(cnt) = (1.0_rprec/dt)*ur(i,j,cnt+1) + RHSur(i,j,cnt+1)
        Bd(cnt) = (1.0_rprec/dt)*vr(i,j,cnt+1) + RHSvr(i,j,cnt+1)
        cnt = cnt + 1
    end do

    ! Last equation in Matrix, k = nzr-2
    nu1 = (nu_wm(i,j,nzr-3) + nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr-2)-zr(nzr-4)))
    nu2 = (nu_wm(i,j,nzr-1) + nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr)-zr(nzr-2)))
    Aa((nzr-3)/2) = -nu1
    Ab((nzr-3)/2) = (nu1 + nu2) + (1.0_rprec/dt)
    Ac((nzr-3)/2) = 0.0_rprec

    ! no-slip BC 
    Ad((nzr-3)/2) = nu2*ubot + ((1.0_rprec/dt)*ur(i,j,nvel-1)) + RHSur(i,j,nvel-1)
    Bd((nzr-3)/2) = ((1.0_rprec/dt)*vr(i,j,nvel-1)) + RHSvr(i,j,nvel-1)

    ! Solve matrix equation
    call tdma(Aa,Ab,Ac,Ad,ur1)
    call tdma(Aa,Ab,Ac,Bd,vr1)

    ur(i,j,1) = u1(i,j)
    ur(i,j,2:(nvel-1)) = ur1
    ur(i,j,nvel) = ubot
    vr(i,j,1) = v1(i,j)
    vr(i,j,2:(nvel-1)) = vr1
    vr(i,j,nvel) = 0.0_rprec

    ! Output Boundary Conditions for LES
    ! Assuming inner layer is sufficiently resolved down to the wall so that
    ! the closest point within the inner layer is in the viscous sublayer
    ! Using one-sided finite difference
    dudz(i,j,1) = ( ur(i,j,nvel-1) - ubot ) / (zr(nzr-2) - zr(nzr))
    dvdz(i,j,1) = vr(i,j,nvel-1) / (zr(nzr-2) - zr(nzr))
    txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
    tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)

end do
end do

end subroutine ws_tlwm_noneq_lbc

!*******************************************************************************
subroutine ws_2d3cF_eq_lbc()
!*******************************************************************************
! 
! This subroutine solves a 2nd order ODE which assumes equilibrium and a 
! eddy viscosity model that satisfies the log-law and the viscous sublayer if
! the LES grid is sufficiently resolved.
! 
! Currently the model assumes smooth walls, finite Reynolds number, and takes 
! the inner layer starting at the ihwm'th LES grid point from the wall (on the 
! uv grid).
!  
! The ODE is of the form:
!     d/dz[(nu+nu_t)*dUdz] = 0,  U(z=0) = 0,  U(z=hwm) = Ules
! Centered differencing is used. As a result, the quantity (nu+nu_t) is
! evaluated between grid points which U is solved on. The TDMA algorithm is 
! used to invert the resulting matrix equation.
! 

use param, only : ld, nx, ny, nu_molec, z_i, ubot, u_star, jt_total, ihwm
use sim_param, only : u, v, dudz, dvdz, txz, tyz
use test_filtermodule, only : test_filter
use fft
implicit none
integer :: i, j, k, cnt, kcnt
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec), dimension(nx, ny) :: ustar
complex(rprec) :: nu1, nu2
complex(rprec), dimension(nvel-2) :: Aa, Ab, Ac, Ad, Bd, ur1, vr1 !! for TDMA
complex(rprec), dimension(nx/2+1) :: u1f, v1f, nu_temp
complex(rprec), dimension(nx/2+1,nvel) ::urf, vrf
complex(rprec), dimension(nx/2+1,nzr) :: nuf

! Treatment of LES velocity into wall-model BC
u1 = u(:,:,ihwm)
v1 = v(:,:,ihwm)

! Find constant ustar for wall model eddy viscosity
! 1. Take the average of previous timesteps tau value
if (jt_total > 1) then
    ustar = sqrt(abs(txz(1:nx,1:ny,1)) + abs(tyz(1:nx,1:ny,1)))
else !! first time-step, should use nonzero ustar
    ustar = u_star
end if
! 2. Use constant value specified by user - only for debugging
!ustar = u_star 

! Compute eddy viscosity values on inner layer grid
call wm_eddyvisc(ustar)

do j = 1, ny
    ! Take 1D Fourier transform in x
    call dfftw_execute_dft_r2c( forw_x, u1(:,j), u1f )
    call dfftw_execute_dft_r2c( forw_x, v1(:,j), v1f )

    ! Repeat for the wall model eddy viscosity
    do k = 1, nzr
        call dfftw_execute_dft_r2c( forw_x, nu_wm(:,j,k), nu_temp )
        nuf(:,k) = nu_temp / real(nx,rprec)
    end do

    kcnt = 1 !! re-initialize with each spanwise position
    do i = 1, nx/2 + 1

        ! Integrate along wall-normal for each specified kx mode
        if (kcnt <= size(kxi)) then
            if (i == kxi(kcnt)) then

                ! Normalize
                u1f(i) = u1f(i) / real(nx,rprec)
                v1f(i) = v1f(i) / real(nx,rprec)

                ! 1st equation in Matrix, k = 3
                nu1 = (nuf(i,2) + nu_molec)/((zr(4)-zr(2))*(zr(3)-zr(1)))
                nu2 = (nuf(i,4) + nu_molec)/((zr(4)-zr(2))*(zr(5)-zr(3)))
                Aa(1) = 0
                Ab(1) = -(nu1 + nu2)
                Ac(1) = nu2

                ! BC from LES
                Ad(1) = -nu1*u1f(i)
                Bd(1) = -nu1*v1f(i)

                ! Middle equations in Matrix, k = 5, 7, ..., nzr-4
                cnt = 2
                do k = 5, nzr-4, 2
                    nu1 = (nuf(i,k-1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k)-zr(k-2)))
                    nu2 = (nuf(i,k+1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k+2)-zr(k)))
                    Aa(cnt) = nu1
                    Ab(cnt) = -(nu1 + nu2)
                    Ac(cnt) = nu2
                    Ad(cnt) = 0.0_rprec
                    Bd(cnt) = 0.0_rprec
                    cnt = cnt + 1
                end do

                ! Last equation in Matrix, k = nzr-2
                nu1=(nuf(i,nzr-3)+nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr-2)-zr(nzr-4)))
                nu2=(nuf(i,nzr-1)+nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr)-zr(nzr-2)))
                Aa((nzr-3)/2) = nu1
                Ab((nzr-3)/2) = -(nu1 + nu2)
                Ac((nzr-3)/2) = 0.0_rprec

                ! no-slip BC 
                ! Assuming ubot is zero for now, should use its value for only i = 1
                ! and zero out for all nonzero kx modes
                Ad((nzr-3)/2) = -nu2*ubot
                Bd((nzr-3)/2) = 0.0_rprec

                ! Solve matrix equation
                call tdma_complex(Aa,Ab,Ac,Ad,ur1)
                call tdma_complex(Aa,Ab,Ac,Bd,vr1)

                ! Format inner-layer velocity field
                urf(i,1) = u1f(i)
                urf(i,2:(nvel-1)) = ur1
                urf(i,nvel) = ubot !! assuming ubot = 0
                vrf(i,1) = v1f(i)
                vrf(i,2:(nvel-1)) = vr1
                vrf(i,nvel) = 0.0_rprec

                kcnt = kcnt + 1

            ! zero out unselected kx modes
            else
                urf(i,:) = 0.0_rprec
                vrf(i,:) = 0.0_rprec
            endif

        ! all selected kx modes considered, zero out the rest
        else
            urf(i,:) = 0.0_rprec
            vrf(i,:) = 0.0_rprec
        endif

    end do !! end of kx loop

    do k = 1, nvel
        ! Transform back, kx --> x
        call dfftw_execute_dft_c2r( back_x, urf(:,k), ur(:,j,k) )
        call dfftw_execute_dft_c2r( back_x, vrf(:,k), vr(:,j,k) )
    end do

    ! Output Boundary Conditions for LES
    do i = 1, nx
        dudz(i,j,1) = ( ur(i,j,nvel-1) - ubot ) / (zr(nzr-2) - zr(nzr))
        dvdz(i,j,1) = vr(i,j,nvel-1) / (zr(nzr-2) - zr(nzr))
        txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
        tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)
    end do

end do !! end of y-pos loop

end subroutine ws_2d3cF_eq_lbc

!*******************************************************************************
subroutine ws_rnl_eq_lbc()
!*******************************************************************************
! 
! This subroutine builds on the equilibrium wall-model by incorporating the 
! Restricted Nonlinear (RNL) model philosophy in the x-direction. The ODE 
! consists of only wall-normal derivatives and is used to extrapolate the LES
! velocity from the ihwm'th grid point to the wall. The RNL model is used to 
! reduce number of computations in the streamwise direction. Instead of 
! solving the ODE for each x, we solve for particular kx.
! 
! This subroutine solves a 2nd order ODE which assumes equilibrium and a 
! eddy viscosity model that satisfies the log-law and the viscous sublayer if
! the LES grid is sufficiently resolved.
! 
! Currently the model assumes smooth walls, finite Reynolds number, and takes 
! the inner layer starting at the ihwm'th LES grid point from the wall (on the
! uv grid).
!  
! The ODE is of the form:
!     d/dz[(nu+nu_t)*dUdz] = 0,  U(z=0) = 0,  U(z=hwm) = Ules
! Centered differencing is used. As a result, the quantity (nu+nu_t) is
! evaluated between grid points which U is solved on. The TDMA algorithm is 
! used to invert the resulting matrix equation.
! 

use param, only : ld, ny, jt_total, nu_molec, z_i, ubot, u_star, ihwm
use sim_param, only : u, v, dudz, dvdz, txz, tyz
#ifdef PPOUTPUT_WMLES
use param, only : tavg_calc, tavg_nstart, tavg_nend
#endif
implicit none
integer :: i, j, k, cnt
real(rprec), dimension(nxf+2,ny) :: u1f, v1f
real(rprec) :: ustar, nu1, nu2
real(rprec), dimension(nzr) :: nuf
real(rprec), dimension(nvel-2) :: LHSa, LHSb, LHSc, RHSu, RHSv, ur1, vr1 !! for TDMA
real(rprec), dimension(ld,ny) :: dudzp, dvdzp, txzp, tyzp
#ifdef PPOUTPUT_WMLES
real(rprec), dimension(nxf+2,ny,nvel) :: urf, vrf
#endif

#ifdef PPOUTPUT_WMLES
! Initialize variables 
urf(:,:,:) = 0.0_rprec
vrf(:,:,:) = 0.0_rprec
#endif

! Treatment of LES velocity into wall-model BC
call tlwm_phys2waveF( u(1:ld,1:ny,ihwm), u1f )
call tlwm_phys2waveF( v(1:ld,1:ny,ihwm), v1f )

! Find constant ustar for wall model eddy viscosity
! 1. Take the average of previous timesteps tau value
if (jt_total > 1) then
    ustar = sqrt( abs(txzf(1,1)) + abs(tyzf(1,1)) )
else !! first time-step, should use nonzero ustar
    ustar = u_star
end if
! 2. Use constant value specified by user - only for debugging
!ustar = u_star 

! Compute eddy viscosity values on inner layer grid
call wm_eddyviscF(ustar,nuf)

! Integrate along wall-normal for each kx and ky mode
! Since eddy viscosity is independent of x and y, compute LHS outside of loop
! 1st equation in Matrix, k = 3
nu1 = (nuf(2) + nu_molec)/((zr(4)-zr(2))*(zr(3)-zr(1)))
nu2 = (nuf(4) + nu_molec)/((zr(4)-zr(2))*(zr(5)-zr(3)))
LHSa(1) = 0
LHSb(1) = -(nu1 + nu2)
LHSc(1) = nu2
! Only initializing RHS vector, will be overwritten by LES BC
RHSu(1) = 0.0_rprec
RHSv(1) = 0.0_rprec

! Middle equations in Matrix, k = 5, 7, ..., nzr-4
cnt = 2
do k = 5, nzr-4, 2
    nu1 = (nuf(k-1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k)-zr(k-2)))
    nu2 = (nuf(k+1) + nu_molec)/((zr(k+1)-zr(k-1))*(zr(k+2)-zr(k)))
    LHSa(cnt) = nu1
    LHSb(cnt) = -(nu1 + nu2)
    LHSc(cnt) = nu2
    RHSu(cnt) = 0.0_rprec
    RHSv(cnt) = 0.0_rprec
    cnt = cnt + 1
end do

! Last equation in Matrix, k = nzr-2
nu1=(nuf(nzr-3)+nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr-2)-zr(nzr-4)))
nu2=(nuf(nzr-1)+nu_molec)/((zr(nzr-1)-zr(nzr-3))*(zr(nzr)-zr(nzr-2)))
LHSa((nzr-3)/2) = nu1
LHSb((nzr-3)/2) = -(nu1 + nu2)
LHSc((nzr-3)/2) = 0.0_rprec

! no-slip BC 
! Assuming ubot is zero for now, should use its value for only i = 1
! and zero out for all nonzero kx modes
RHSu((nzr-3)/2) = -nu2*ubot
RHSv((nzr-3)/2) = 0.0_rprec

! Implement LES BC on RHS vector for each kx,ky and solve
nu1 = (nuf(2) + nu_molec)/((zr(4)-zr(2))*(zr(3)-zr(1))) !! overwrite for LES BC
do j = 1, ny
do i = 1, nxf

! BC from LES
RHSu(1) = -nu1*u1f(i,j)
RHSv(1) = -nu1*v1f(i,j)

! Solve matrix equation
call tdma(LHSa,LHSb,LHSc,RHSu,ur1)
call tdma(LHSa,LHSb,LHSc,RHSv,vr1)

#ifdef PPOUTPUT_WMLES
! Only if time-averaging
if ((tavg_calc) .and. (jt_total>=tavg_nstart) .and. (jt_total<=tavg_nend)) then
! Format inner-layer velocity field
urf(i,j,1) = u1f(i,j)
urf(i,j,2:(nvel-1)) = ur1
urf(i,j,nvel) = ubot
vrf(i,j,1) = v1f(i,j)
vrf(i,j,2:(nvel-1)) = vr1
vrf(i,j,nvel) = 0.0_rprec
endif
#endif

! Compute derivatives for stress 
dudzf(i,j) = ( ur1(nvel-2) - ubot ) / (zr(nzr-2) - zr(nzr))
dvdzf(i,j) = vr1(nvel-2) / (zr(nzr-2) - zr(nzr))
txzf(i,j) = -nu_molec/(z_i*u_star)*dudzf(i,j)
tyzf(i,j) = -nu_molec/(z_i*u_star)*dvdzf(i,j)

end do !! end of kx loop
end do !! end of ky loop

! Transform back to output BC for LES
! dudzp, txzp, etc. are temp variables
call tlwm_wave2physF(dudzf,dudzp)
call tlwm_wave2physF(dvdzf,dvdzp)
call tlwm_wave2physF(txzf,txzp)
call tlwm_wave2physF(tyzf,tyzp)
dudz(:,:,1) = dudzp(:,:)
dvdz(:,:,1) = dvdzp(:,:)
txz(:,:,1) = txzp(:,:)
tyz(:,:,1) = tyzp(:,:)

#ifdef PPOUTPUT_WMLES
! Only transform back if time-averaging
if ((tavg_calc) .and. (jt_total>=tavg_nstart) .and. (jt_total<=tavg_nend)) then
! Transform each zr-level to output inner-layer velocity field
do k = 1, nvel
    call tlwm_wave2physF(urf(:,:,k),ur(:,:,k))
    call tlwm_wave2physF(vrf(:,:,k),vr(:,:,k))
enddo
! Eddy viscosity is streamwise and spanwise averaged 
! therefore just distribute the same value
do k = 1, nzr
    nu_wm(:,:,k) = nuf(k)
enddo

endif
#endif

end subroutine ws_rnl_eq_lbc

!*******************************************************************************
subroutine wm_eddyvisc(ustar)
!*******************************************************************************
! 
! This subroutine represents the eddy viscosity function used in the two layer
! wall model. It returns the eddy viscosity value for the entire inner layer 
! domain.
! 
! The wall-model coefficient kappa can be either:
!     1. constant, kappa = vonk
!     2. wall-normal constant, sgs dynamic
!     3. wall-normal dynamic, sgs dynamic
!
! The second method was used in the study:
! M. Wang and P. Moin, "Dynamic wall modeling for large-eddy simulation of 
! complex turbulent flows," Phys. Fluids 14, 2043-2051 (2012).
! 
! The third method was used in the study:
! S. Kawai and J. Larsson, "Dynamic non-equilibrium wall-modeling for large
! eddy simulation at high Reynolds numbers," Phys. Fluids 25, (2013).
!  

use types, only : rprec
use param, only : nu_molec, nx, ny, vonk, dx, dy, ihwm
use sgs_param, only : Nu_t
use functions, only : y_avg
implicit none
real(rprec), dimension(nx, ny), intent(in) :: ustar
real(rprec) :: a_plus, alpha, max_dx, zcr, bigK
real(rprec), dimension(nx) :: kwm_temp
real(rprec), dimension(nx, nzr) :: kwm
real(rprec), dimension(nx, ny, nzr) :: decay
integer :: i, j, k

! Model parameters
a_plus = 17.0_rprec !! used in decay portion
alpha = 0.48_rprec !! proportionality constant for kwm in ktype = 3
max_dx = max(dx,dy)
zcr = alpha*max_dx

! Compute decay portion of eddy viscosity
do i = 1, nx
do j = 1, ny
do k = 1, nzr
    decay(i,j,k) = (1.0_rprec - exp(-zr(k)*ustar(i,j)/(a_plus*nu_molec)))**2
end do
end do
end do

! Compute wall-model coefficient kwm
if ((ktype == 1) .or. ((ktype == 3).and.(hwm <= zcr)) )then
    kwm = vonk
    ! the second condition is to ensure that bigK > 0 for ktype == 3
elseif (ktype == 2) then

    kwm_temp = y_avg(Nu_t(1:nx,1:ny,ihwm))/y_avg(ustar(:,:)*hwm*decay(:,:,1))
    do k = 1, nzr
        kwm(:,k) = kwm_temp
    end do

elseif (ktype == 3) then

    kwm_temp = y_avg(Nu_t(1:nx,1:ny,ihwm))/y_avg(ustar(:,:)*hwm*decay(:,:,1))
    do k = 1, nzr
        bigK = min( 1.0_rprec, (hwm-zr(k))/(hwm-zcr) )
        kwm(:,k) = vonk*bigK + kwm_temp*(1.0_rprec - bigK)
    end do

end if

! Compute eddy viscosity
do i = 1, nx
do j = 1, ny
do k = 1, nzr
    nu_wm(i,j,k) = kwm(i,k)*ustar(i,j)*zr(k)*decay(i,j,k)
end do
end do
end do

end subroutine wm_eddyvisc

!*******************************************************************************
subroutine wm_eddyviscF(ustar, nu_wmf)
!*******************************************************************************
! 
! This subroutine outputs the eddy viscosity for the RNL TLWMLES. It is assumed
! that input ustar is streamwise and spanwise averaged and therefore so is the
! eddy viscosity. Hence this subroutine differs from wm_eddyvisc, since the 
! eddy viscosity outputted here is not a parameter in module tlwmles.
!  

use types, only : rprec
use param, only : nu_molec, vonk
implicit none
real(rprec), intent(in) :: ustar
real(rprec), dimension(nzr), intent(out) :: nu_wmf
real(rprec) :: a_plus, decay
integer :: k

! Model parameters
a_plus = 17.0_rprec !! used in decay portion

! Compute eddy viscosity
do k = 1, nzr
    ! Using streamwise and spanwise average
    decay = (1.0_rprec - exp(-zr(k)*ustar/(a_plus*nu_molec)))**2.0_rprec
    nu_wmf(k) = vonk*ustar*zr(k)*decay
end do   

end subroutine wm_eddyviscF

!*******************************************************************************
subroutine tdma(a,b,c,d,x)
!*******************************************************************************
!
! Size of b, d, and x is size of square coefficient matrix, n
! Size of a and c are off-diagonal entries and are therefore of size n-1
! However, a(1) = 0 and c(n) = 0 when inputted into this subroutine
!
! The tridiagonal system is represented as:
!     a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i
! 
use types, only : rprec
implicit none
integer :: n, i
real(rprec), dimension(:), intent(in) :: a, b, c, d
real(rprec), dimension(:), intent(out) :: x
real(rprec), dimension(:), allocatable :: cn, dn

! Find size of tridiagonal matrix equation
n = size(d)

! Forward sweep
allocate(cn(n))
allocate(dn(n))
cn(1) = c(1) / b(1)
dn(1) = d(1) / b(1)

do i = 2, (n-1)
    cn(i) = c(i) / (b(i) - a(i)*cn(i-1))
    dn(i) = ( d(i) - a(i)*dn(i-1) ) / ( b(i) - a(i)*cn(i-1) )
end do

dn(n) = ( d(n) - a(n)*dn(n-1) ) / ( b(n) - a(n)*cn(n-1) )

! Backward substitution
x(n) = dn(n)

do i = (n-1), 1, -1
    x(i) = dn(i) - cn(i)*x(i+1)
end do

end subroutine tdma

!*******************************************************************************
subroutine tdma_complex(a,b,c,d,x)
!*******************************************************************************
!
! Size of b, d, and x is size of square coefficient matrix, n
! Size of a and c are off-diagonal entries and are therefore of size n-1
! However, a(1) = 0 and c(n) = 0 when inputted into this subroutine
!
! The tridiagonal system is represented as:
!     a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i
! 
use types, only : rprec
implicit none
integer :: n, i
complex(rprec), dimension(:), intent(in) :: a, b, c, d
complex(rprec), dimension(:), intent(out) :: x
complex(rprec), dimension(:), allocatable :: cn, dn

! Find size of tridiagonal matrix equation
n = size(d)

! Forward sweep
allocate(cn(n))
allocate(dn(n))
cn(1) = c(1) / b(1)
dn(1) = d(1) / b(1)

do i = 2, (n-1)
    cn(i) = c(i) / (b(i) - a(i)*cn(i-1))
    dn(i) = ( d(i) - a(i)*dn(i-1) ) / ( b(i) - a(i)*cn(i-1) )
end do

dn(n) = ( d(n) - a(n)*dn(n-1) ) / ( b(n) - a(n)*cn(n-1) )

! Backward substitution
x(n) = dn(n)

do i = (n-1), 1, -1
    x(i) = dn(i) - cn(i)*x(i+1)
end do

end subroutine tdma_complex

!*******************************************************************************
subroutine tlwm_ddxy(f,dfdx,dfdy)
!*******************************************************************************
! 
! This subroutine mimics ddxy in derivatives.f90, however instead of looping 
! through horizontal slices of the outer layer domain, this routine loops
! through horizontal slices of the inner layer domain.
! 
use types, only : rprec
use param, only : ld, nx, ny
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

real(rprec), dimension(:,:,:), intent(in) :: f
real(rprec), dimension(:,:,:), intent(inout) :: dfdx, dfdy
real(rprec) :: const
integer :: jz

const = 1._rprec / ( nx * ny )

! Loop through horizontal slices, don't include the wall
do jz = 1, (nvel-1)
    ! Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = const*f(:,:,jz)
    call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz), dfdx(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdx(ld-1:ld,:,jz) = 0._rprec
    dfdx(:,ny/2+1,jz) = 0._rprec

    ! Derivatives: must to y's first here, because we're using dfdx as storage
    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdx(:,:,jz) .MULI. ky
    dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
    call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine tlwm_ddxy

!*******************************************************************************
subroutine tlwm_phys2waveF( u, uhat )
!*******************************************************************************
! 
! Transform inner layer field from physical space to wavenumber space.
! (x,y) --> (kx,ky)
! This function is intended for tlwm_fourier.
! 
! This routine is similar to phys2waveF in derivatives.f90, however phys2waveF
! is to be used on the entire domain. tlwm_phys2waveF only transforms on a
! given z-slice specified by the input. Another difference is that phy2waveF
! also outputs the input variable, this routine makes no changes to the input.
! 

use types, only : rprec
use param, only : tlwm_kxnum, tlwm_kxin, nx, ny
use fft
implicit none

integer :: jx, ii, ir, iih, irh
real(rprec) :: const
real(rprec), dimension(ld, ny), intent(inout) :: u
real(rprec), dimension(ld, ny) :: utemp
real(rprec), dimension(nxf+2, ny), intent(out) :: uhat

const = 1._rprec / (nx*ny)

! Initialize uhat
uhat = 0.0_rprec

! in-place Transform
utemp(:,:) = const * u(:,:) !! normalization
call dfftw_execute_dft_r2c(forw, utemp(:,:), utemp(:,:) )

! Only interested in u values at tlwm_kxin wavenumbers
do jx = 1, tlwm_kxnum
    ! uhat indices
    iih = 2*jx
    irh = iih - 1

    ! u indices
    ii = 2 * int( tlwm_kxin(jx) ) + 2
    ir = ii - 1

    uhat(irh:iih, :) = utemp(ir:ii, :)
enddo

end subroutine tlwm_phys2waveF

!*******************************************************************************
subroutine tlwm_wave2physF( uhat, u )
!*******************************************************************************
! 
! Transform inner layer field from physical space to wavenumber space.
! (kx,ky) --> (x,y)
! This function is intended for tlwm_fourier.
! 
! This routine is similar to wave2physF in derivatives.f90, however wave2physF
! is to be used on the entire domain. tlwm_wave2physF only transforms on a
! given z-slice specified by the input.
! 

use types, only : rprec
use param, only : ld, ny, tlwm_kxnum, tlwm_kxin
use fft
implicit none

integer :: jx, iih, irh, ii, ir
real(rprec), dimension(nxf+2,ny), intent(in) :: uhat
real(rprec), dimension(ld,ny), intent(out) :: u

! Fill physical u with zeros everywhere
u(:,:) = 0.0_rprec

! Place uhat values in appropriate u index before inverse fourier transform
do jx = 1, tlwm_kxnum
    ! uhat indices
    iih = 2*jx
    irh = iih - 1

    ! u indices 
    ii = 2 * int( tlwm_kxin(jx) ) + 2
    ir = ii - 1

    u(ir:ii, :) = uhat( irh:iih, :)
enddo

! in-place Transform
call dfftw_execute_dft_r2c(back, u(1:ld,1:ny), u(1:ld,1:ny) )

end subroutine tlwm_wave2physF

!*******************************************************************************
subroutine tlwm_finalize
!*******************************************************************************
!
! This subroutine outputs wmles results and deallocates memory used for tlwm
!
use param, only : lbc_mom
#ifdef PPOUTPUT_WMLES
use param, only : tavg_calc, jt_total, tavg_nstart
#endif
implicit none

! Perform final checkpoint
call tlwm_checkpoint()

! Write time-averaging files
#ifdef PPOUTPUT_WMLES
if ((tavg_calc) .and. (jt_total > tavg_nstart)) call tavg_wmles_finalize()
#endif

deallocate(zr)
deallocate(ur)
deallocate(vr)
deallocate(nu_wm)

if (lbc_mom == 6) then
    deallocate(wr)
    deallocate(durdx)
    deallocate(durdy)
    deallocate(dvrdx)
    deallocate(dvrdy)
end if

end subroutine tlwm_finalize

!*******************************************************************************
subroutine tlwm_checkpoint
!*******************************************************************************
!
! Checkpoints variables required for tlwm and runs any checkpoint routines
!

use param, only : nx, ny, path, write_endian
use sim_param, only : txz, tyz
#ifdef PPOUTPUT_WMLES
use stat_defs, only : tavg_initialized
use param, only : tavg_calc
#endif
implicit none
character(64) :: fname

! Checkpoint necessary variables for the wall-model
fname = path // 'wmles.out'
open(11, file=fname, form='unformatted', convert=write_endian,      &
    status='unknown', position='rewind')
write(11) hwm, dzr, nzr, ss, txz(1:nx,1:ny,1), tyz(1:nx,1:ny,1)
close(11)

! Checkpoint time averagning restart data
#ifdef PPOUTPUT_WMLES
if (tavg_calc .and. tavg_initialized ) call tavg_wmles_checkpoint()
#endif

end subroutine tlwm_checkpoint

#ifdef PPOUTPUT_WMLES
!*******************************************************************************
subroutine wmles_inst_write()
!*******************************************************************************
! 
! This subroutine writes an instantaneous snapshot of the inner layer velocity
! field. This is similar to io.f90/inst_write, but only for tlwmles and does 
! not write for plane snapshots, only for the full domain.
! 
! Also, only writes binary output, not CGNS.
! 
use param, only : nx, ny, write_endian, jt_total, path
use string_util, only : string_splice, string_concat
implicit none
character(64) :: fname

! write binary output
call string_splice(fname, path // 'output/wmles_vel.', jt_total)
call string_concat(fname, '.bin')
open(unit=13, file=fname, form='unformatted', convert=write_endian,          &
    access='direct', recl=nx*ny*nvel*rprec)
write(13,rec=1) ur(:nx,:ny,1:nvel)
write(13,rec=2) vr(:nx,:ny,1:nvel)
close(13)

end subroutine wmles_inst_write

!*******************************************************************************
subroutine tavg_wmles_init()
!*******************************************************************************
!
! This subroutine loads the tavg_wmles.out file or creates it
!
use types, only : rprec
use param, only : read_endian, path
use stat_defs, only : tavg_wmles_total_time, tavg_wmles
implicit none
character (*), parameter :: ftavg_wmles_in = path // 'tavg_wmles.out'
logical :: exst

! Does a tavg file already exist?
inquire (file=ftavg_wmles_in, exist=exst)
if (.not. exst) then !! No, so initialize tavg_total_time
    tavg_wmles_total_time = 0._rprec
else !! Yes, extract tavg_total_time and tavg data
    open(1, file=ftavg_wmles_in, action='read', position='rewind',         &
        form='unformatted', convert=read_endian)
    read(1) tavg_wmles_total_time
    read(1) tavg_wmles
    close(1)
end if

end subroutine tavg_wmles_init

!*******************************************************************************
subroutine tavg_wmles_compute()
!*******************************************************************************
!
! Running time-average computations
!

use param, only : nx, ny
use stat_defs, only : tavg_wmles, tavg_wmles_total_time, tavg_wmles_dt
implicit none
integer :: i, j, k
real(rprec), dimension(nx,ny,nvel) :: nu_wm_temp

! Extract information at every point onto only the velocity locations
nu_wm_temp(1:nx,:,1:nvel) = nu_wm(1:nx,:,1:nzr:2)

! time-averaging sum
do j = 1, ny
do i = 1, nx
do k = 1, nvel
    tavg_wmles(i,j,k) % u = tavg_wmles(i,j,k) % u + ur(i,j,k) * tavg_wmles_dt
    tavg_wmles(i,j,k) % v = tavg_wmles(i,j,k) % v + vr(i,j,k) * tavg_wmles_dt
    tavg_wmles(i,j,k) % nu = tavg_wmles(i,j,k) % nu + nu_wm_temp(i,j,k) * tavg_wmles_dt

    tavg_wmles(i,j,k) % uu = tavg_wmles(i,j,k) % uu + ur(i,j,k) * ur(i,j,k) * tavg_wmles_dt
    tavg_wmles(i,j,k) % vv = tavg_wmles(i,j,k) % vv + vr(i,j,k) * vr(i,j,k) * tavg_wmles_dt
    tavg_wmles(i,j,k) % uv = tavg_wmles(i,j,k) % uv + ur(i,j,k) * vr(i,j,k) * tavg_wmles_dt

end do
end do
end do

! Update tavg_wmles_total_time for variable time stepping
tavg_wmles_total_time = tavg_wmles_total_time + tavg_wmles_dt

! Set tavg_wmles_dt back to zero for next increment
tavg_wmles_dt = 0._rprec

end subroutine tavg_wmles_compute

!*******************************************************************************
subroutine tavg_wmles_checkpoint
!*******************************************************************************
!
! This subroutine writes the restart data for the wall model statistics. It is 
! called by tavg_wmles_finalize and tlwm_checkpoint.
!
use param, only : checkpoint_tavg_wmles_file, write_endian
use stat_defs, only : tavg_wmles, tavg_wmles_total_time
implicit none

! Write data to tavg_wmles.out
open(1, file=checkpoint_tavg_wmles_file, action='write', position='rewind',   &
    form='unformatted', convert=write_endian)
write(1) tavg_wmles_total_time
write(1) tavg_wmles
close(1)

end subroutine tavg_wmles_checkpoint

!*******************************************************************************
subroutine tavg_wmles_finalize
!*******************************************************************************
!
! Writes time-averaged data to be outputted.
! Currently only uses a binary file
! 

use stat_defs, only : tavg_wmles, tavg_wmles_t, tavg_wmles_total_time
use param, only : nx, ny, write_endian, path
implicit none
character(64) :: fname_grid, fname_wmvel, fname_wmrs
integer :: i, j, k
real(rprec), dimension(nx,ny,nvel) :: upup, vpvp, upvp

fname_grid = path // 'output/wmles_grid.bin'
fname_wmvel = path // 'output/wmvel_avg.bin'
fname_wmrs = path // 'output/wmrs_avg.bin'

! Final checkpoint all restart data
call tavg_wmles_checkpoint()

! Normalize sums by total time
do j = 1, ny
do i = 1, nx
    do k = 1, nvel
    tavg_wmles(i,j,k) % u = tavg_wmles(i,j,k) % u / tavg_wmles_total_time
    tavg_wmles(i,j,k) % v = tavg_wmles(i,j,k) % v / tavg_wmles_total_time
    tavg_wmles(i,j,k) % nu = tavg_wmles(i,j,k) % nu / tavg_wmles_total_time

    tavg_wmles(i,j,k) % uu = tavg_wmles(i,j,k) % uu / tavg_wmles_total_time
    tavg_wmles(i,j,k) % vv = tavg_wmles(i,j,k) % vv / tavg_wmles_total_time
    tavg_wmles(i,j,k) % uv = tavg_wmles(i,j,k) % uv / tavg_wmles_total_time
    end do
end do 
end do

! Write binary data
! write wall-model grid to file
open(unit=13, file=fname_grid, form='unformatted', convert=write_endian, &
    access='direct',recl=nzr*rprec) 
write(13,rec=1) zr(1:nzr)
close(13)

! write variables on velocity grid
open(unit=13, file=fname_wmvel, form='unformatted', convert=write_endian, &
    access='direct',recl=nx*ny*nvel*rprec) 
write(13,rec=1) tavg_wmles(:nx,:ny,1:nvel)%u
write(13,rec=2) tavg_wmles(:nx,:ny,1:nvel)%v
write(13,rec=3) tavg_wmles(:nx,:ny,1:nvel)%nu
close(13)

! Compute Reynolds stresses
do j = 1, ny
do i = 1, nx
do k = 1, nvel
upup(i,j,k) = tavg_wmles(i,j,k) % uu - tavg_wmles(i,j,k) % u * tavg_wmles(i,j,k) % u
vpvp(i,j,k) = tavg_wmles(i,j,k) % vv - tavg_wmles(i,j,k) % v * tavg_wmles(i,j,k) % v
upvp(i,j,k) = tavg_wmles(i,j,k) % uv - tavg_wmles(i,j,k) % u * tavg_wmles(i,j,k) % v
enddo
enddo
enddo


! Output Reynolds stresses
open(unit=13, file=fname_wmrs, form='unformatted', convert=write_endian, &
    access='direct',recl=nx*ny*nvel*rprec) 
write(13,rec=1) upup(:nx,:ny,1:nvel)
write(13,rec=2) vpvp(:nx,:ny,1:nvel)
write(13,rec=3) upvp(:nx,:ny,1:nvel)
close(13)

end subroutine tavg_wmles_finalize
#endif

end module tlwmles
