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
subroutine diff_stag_array_uv(u,v,RHSx,RHSy,RHSx_f,RHSy_f,txz_half2,tyz_half2,txz,tyz)
!*******************************************************************************
!
! Calculate the intermediate xy-velocity from implicit CN scheme on exit.
!
use types, only : rprec
use param
use messages
use sgs_param, only : nu, Nu_t
use derivatives, only : ddz_w
use fft
#ifdef PPMAPPING
use sim_param, only : jaco_w, jaco_uv, mesh_stretch
#endif
use derivatives, only : dft_direct_forw_2d_n_yonlyC, &
    dft_direct_back_2d_n_yonlyC

implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: u, v
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: RHSx, RHSy, RHSx_f, RHSy_f
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: txz_half2, tyz_half2, txz, tyz

real(rprec), dimension(ld,ny,0:nz) :: Rx, usol, Ry, vsol
real(rprec), dimension(nx,ny,0:nz) :: a, b, c
real(rprec), dimension(ld,ny,lbz:nz) :: dtxzdz_rhs, dtyzdz_rhs
real(rprec) :: nu_a, nu_b, nu_c, nu_r, const1, const2, const3
integer :: jx, jy, jz, jz_min, jz_max

! Set constants
const1 = (dt/(2._rprec*dz))
const2 = (1._rprec/dz)
const3 = (1._rprec/(0.5_rprec*dz))

! Get the RHS ready
! Initialize with the explicit terms
Rx(:,:,1:nz-1) = u(:,:,1:nz-1) +                                     &
    dt * ( tadv1 * RHSx(:,:,1:nz-1) + tadv2 * RHSx_f(:,:,1:nz-1) )
Ry(:,:,1:nz-1) = v(:,:,1:nz-1) +                                     &
    dt * ( tadv1 * RHSy(:,:,1:nz-1) + tadv2 * RHSy_f(:,:,1:nz-1) )

! Add explicit portion of Crank-Nicolson
call ddz_w(txz_half2, dtxzdz_rhs, lbz)
call ddz_w(tyz_half2, dtyzdz_rhs, lbz)
dtxzdz_rhs(ld-1:ld, :, 1:nz-1) = 0._rprec
dtyzdz_rhs(ld-1:ld, :, 1:nz-1) = 0._rprec
#ifdef PPSAFETYMODE
#ifdef PPMPI
dtxzdz_rhs(:,:,0) = BOGUS
dtyzdz_rhs(:,:,0) = BOGUS
#endif
dtxzdz_rhs(:,:,nz) = BOGUS
dtyzdz_rhs(:,:,nz) = BOGUS
#endif
! Assuming fixed dt here!
Rx(:,:,1:nz-1) = Rx(:,:,1:nz-1) - dt * 0.5_rprec * dtxzdz_rhs(:,:,1:nz-1)
Ry(:,:,1:nz-1) = Ry(:,:,1:nz-1) - dt * 0.5_rprec * dtyzdz_rhs(:,:,1:nz-1)

! Transform eddy viscosity kx, ky, z --> kx, y, z
! however only going to use the kx=0 mode, assuming Nu_t(y,z) only
if ((fourier) .and. (sgs)) then
    do jz = 1, nz
        call dft_direct_back_2d_n_yonlyC( Nu_t(:,:,jz) )
    end do
endif

! Get bottom row 
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    select case (lbc_mom)
        ! Stress free
        case (0)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on w, not needed here
                    if (fourier) then
                        nu_c = Nu_t(1,jy,2) + nu
                    else
                        nu_c = Nu_t(jx,jy,2) + nu
                    endif
                else
                    nu_c = nu
                end if
                ! txz(jx,jy,1) = 0, so nothing added to RHS
#ifdef PPMAPPING
                b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_uv(1))*        &
                    const2*(1._rprec/jaco_w(2))*nu_c
                c(jx,jy,1) = -const1*(1._rprec/jaco_uv(1))*const2*(1._rprec/jaco_w(2))*nu_c
#else
                b(jx,jy,1) = 1._rprec + const1*const2*nu_c
                c(jx,jy,1) = -const1*const2*nu_c
#endif
            end do
            end do

        ! DNS BC or Wall-resolved
        case (1)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on uvp, not needed here
                    if (fourier) then
                        nu_c = Nu_t(1,jy,2) + nu
                    else
                        nu_c = Nu_t(jx,jy,2) + nu
                    endif
                else
                    nu_c = nu
                end if
                ! Discretized txz(jx,jy,1) as in wallstress,
                ! Therefore BC treated implicitly
                ! vbot = 0 so no changes to Ry
#ifdef PPMAPPING
                b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_uv(1))*        &
                    (const2*(1._rprec/jaco_w(2))*nu_c + (nu/mesh_stretch(1)))
                c(jx,jy,1) = -const1*(1._rprec/jaco_uv(1))*const2*(1._rprec/jaco_w(2))*nu_c
                if (fourier) then
                    Rx(1,1,1) = Rx(1,1,1) + const1*(1._rprec/jaco_uv(1))* &
                        (nu/mesh_stretch(1))*ubot
                else
                    Rx(jx,jy,1) = Rx(jx,jy,1) + const1*(1._rprec/jaco_uv(1))* &
                        (nu/mesh_stretch(1))*ubot
                endif
#else
                b(jx,jy,1) = 1._rprec + const1*(const2*nu_c + const3*nu)
                c(jx,jy,1) = -const1*const2*nu_c
                if (fourier) then
                    Rx(1,1,1) = Rx(1,1,1) + const1*const3*nu*ubot
                else
                    Rx(jx,jy,1) = Rx(jx,jy,1) + const1*const3*nu*ubot
                endif
#endif
            end do
            end do

        ! Wall-model 
        case (2:)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on uvp, not needed here
                    if (fourier) then
                        nu_c = Nu_t(1,jy,2) + nu
                    else
                        nu_c = Nu_t(jx,jy,2) + nu
                    endif
                else
                    nu_c = nu
                end if
                ! Treating txz(jx,jy,1) from wallstress explicitly
#ifdef PPMAPPING
                b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_uv(1))*            &
                    const2*(1._rprec/jaco_w(2))*nu_c
                c(jx,jy,1) = -const1*(1._rprec/jaco_uv(1))*const2*(1._rprec/jaco_w(2))*nu_c
                Rx(jx,jy,1) = Rx(jx,jy,1) + const1*(1._rprec/jaco_uv(1))*txz(jx,jy,1)
                Ry(jx,jy,1) = Ry(jx,jy,1) + const1*(1._rprec/jaco_uv(1))*tyz(jx,jy,1)
#else
                b(jx,jy,1) = 1._rprec + const1*const2*nu_c
                c(jx,jy,1) = -const1*const2*nu_c
                Rx(jx,jy,1) = Rx(jx,jy,1) + const1*txz(jx,jy,1)
                Ry(jx,jy,1) = Ry(jx,jy,1) + const1*tyz(jx,jy,1)
#endif
            end do
            end do

    end select
    jz_min = 2
else
    jz_min = 1
endif

! Get top row
#ifdef PPMPI
if (coord == nproc-1) then
#endif
#ifdef PPSAFETYMODE
    c(:,:,nz-1) = BOGUS
#endif
    select case (ubc_mom)
        ! Stress free
        case (0)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,nz) on w, not needed here
                    if (fourier) then
                        nu_a = Nu_t(1,jy,nz-1) + nu
                    else
                        nu_a = Nu_t(jx,jy,nz-1) + nu
                    endif
                else
                    nu_a = nu
                end if
                ! txz(jx,jy,nz) = 0, so nothing added to RHS
#ifdef PPMAPPING
                a(jx,jy,nz-1) = -const1*(1._rprec/jaco_uv(nz-1))*const2*(1._rprec/jaco_w(nz-1))*nu_a
                b(jx,jy,nz-1) = 1._rprec +                               &
                    const1*(1._rprec/jaco_uv(nz-1))*const2*(1._rprec/jaco_w(nz-1))*nu_a
#else
                a(jx,jy,nz-1) = -const1*const2*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*const2*nu_a
#endif
            end do
            end do

        ! DNS BC or Wall-resolved
        case (1)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,nz) on uvp nz-1, not needed here
                    if (fourier) then
                        nu_a = Nu_t(1,jy,nz-1) + nu
                    else
                        nu_a = Nu_t(jx,jy,nz-1) + nu
                    endif
                else
                    nu_a = nu
                endif
                ! Discretized txz(jx,jy,nz) as in wallstress,
                ! Therefore BC treated implicitly
                ! vtop = 0 so no changes to Ry
#ifdef PPMAPPING
                a(jx,jy,nz-1) = -const1*(1._rprec/jaco_uv(nz-1))*const2*(1._rprec/jaco_w(nz-1))*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/jaco_uv(nz-1))*          &
                    (const2*(1._rprec/jaco_w(nz-1))*nu_a + (nu/(L_z-mesh_stretch(nz-1))))
                if (fourier) then
                    Rx(1,1,nz-1) = Rx(1,1,nz-1) + const1*(1._rprec/jaco_uv(nz-1))* &
                        (nu/(L_z-mesh_stretch(nz-1)))*utop
                else
                    Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) + const1*(1._rprec/jaco_uv(nz-1))* &
                        (nu/(L_z-mesh_stretch(nz-1)))*utop
                endif
#else
                a(jx,jy,nz-1) = -const1*const2*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*(const2*nu_a + const3*nu)
                if (fourier) then
                    Rx(1,1,nz-1) = Rx(1,1,nz-1) + const1*const3*nu*utop
                else
                    Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) + const1*const3*nu*utop
                endif
#endif
            end do
            end do

        ! Wall-model
        case (2:)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,nz) on uvp, not needed here
                    if (fourier) then
                        nu_a = Nu_t(1,jy,nz-1) + nu
                    else
                        nu_a = Nu_t(jx,jy,nz-1) + nu
                    endif
                else
                    nu_a = nu
                end if
                ! Treating txz(jx,jy,nz) from wallstress explicitly
#ifdef PPMAPPING
                a(jx,jy,nz-1) = -const1*(1._rprec/jaco_uv(nz-1))*const2*(1._rprec/jaco_w(nz-1))*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/jaco_uv(nz-1))*            &
                    const2*(1._rprec/jaco_w(nz-1))*nu_a
                Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) - const1*(1._rprec/jaco_uv(nz-1))*txz(jx,jy,nz)
                Ry(jx,jy,nz-1) = Ry(jx,jy,nz-1) - const1*(1._rprec/jaco_uv(nz-1))*tyz(jx,jy,nz)
#else
                a(jx,jy,nz-1) = -const1*const2*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*const2*nu_a
                Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) - const1*txz(jx,jy,nz)
                Ry(jx,jy,nz-1) = Ry(jx,jy,nz-1) - const1*tyz(jx,jy,nz)
#endif
            end do
            end do

    end select
    jz_max = nz-2
#ifdef PPMPI
else
    jz_max = nz-1
endif
#endif

! Compute coefficients in domain
! Treating SGS viscosity explicitly!
do jz = jz_min, jz_max
do jy = 1, ny
do jx = 1, nx
    if (sgs) then
        if (fourier) then
            nu_a = Nu_t(1,jy,jz) + nu
#ifdef PPMAPPING
            nu_b = ((Nu_t(1,jy,jz+1)+nu)/jaco_w(jz+1)) + ((Nu_t(1,jy,jz)+nu)/jaco_w(jz))
#else
            nu_b = Nu_t(1,jy,jz+1) + Nu_t(1,jy,jz) + 2._rprec*nu
#endif
            nu_c = Nu_t(1,jy,jz+1) + nu
        else
            nu_a = Nu_t(jx,jy,jz) + nu
#ifdef PPMAPPING
            nu_b = ((Nu_t(jx,jy,jz+1)+nu)/jaco_w(jz+1)) + ((Nu_t(jx,jy,jz)+nu)/jaco_w(jz))
#else
            nu_b = Nu_t(jx,jy,jz+1) + Nu_t(jx,jy,jz) + 2._rprec*nu
#endif
            nu_c = Nu_t(jx,jy,jz+1) + nu
        endif
    else
        nu_a = nu
#ifdef PPMAPPING
        nu_b = (nu/jaco_w(jz+1)) + (nu/jaco_w(jz))
#else
        nu_b = 2._rprec*nu
#endif
        nu_c = nu
    endif

#ifdef PPMAPPING
    a(jx, jy, jz) = -const1*(1._rprec/jaco_uv(jz))*const2*(1._rprec/jaco_w(jz))*nu_a
    b(jx, jy, jz) = 1._rprec + const1*(1._rprec/jaco_uv(jz))*const2*nu_b
    c(jx, jy, jz) = -const1*(1._rprec/jaco_uv(jz))*const2*(1._rprec/jaco_w(jz+1))*nu_c
#else
    a(jx, jy, jz) = -const1*const2*nu_a
    b(jx, jy, jz) = 1._rprec + const1*const2*nu_b
    c(jx, jy, jz) = -const1*const2*nu_c
#endif
end do
end do
end do

! Done using eddy viscosity, transform kx, y, z --> kx, ky, z
! Now transform Rx and Ry, kx, ky, z --> kx, y, z
if ((fourier) .and. (sgs)) then
    do jz = 1, nz
        call dft_direct_forw_2d_n_yonlyC( Nu_t(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( Rx(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( Ry(:,:,jz) )
    end do
endif

! Find intermediate velocity in TDMA
call tridag_array_diff_uv (a, b, c, Rx, usol)
call tridag_array_diff_uv (a, b, c, Ry, vsol)

! Since a,b,c(y,z) and Rx,Ry(kx,y,z) we have usol,vsol(kx,y,z)
! No need to change Rx, Ry back, will be overwritten in next time-step
if ((fourier) .and. (sgs)) then
    do jz = 1, nz-1
        call dft_direct_forw_2d_n_yonlyC( usol(:,:,jz) )
        call dft_direct_forw_2d_n_yonlyC( vsol(:,:,jz) )
    enddo
endif

! Fill velocity solution
u(:nx,:ny,1:nz-1) = usol(:nx,:ny,1:nz-1)
v(:nx,:ny,1:nz-1) = vsol(:nx,:ny,1:nz-1)

end subroutine diff_stag_array_uv
