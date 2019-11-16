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
subroutine diff_stag_array_w()
!*******************************************************************************
!
! Calculate the intermediate z-velocity from implicit CN scheme on exit.
!
use types, only : rprec
use param
use messages
use sim_param, only : w, RHSz, RHSz_f, tzz
use sgs_param, only : nu, Nu_t
use derivatives, only : ddz_uv
use fft
#ifdef PPMAPPING
use sim_param, only : JACO1, JACO2, mesh_stretch
#endif

implicit none

real(rprec), dimension(ld,ny,0:nz) :: Rz, wsol
real(rprec), dimension(nx,ny,0:nz) :: a, b, c
real(rprec), dimension(ld,ny,lbz:nz) :: dtzzdz_rhs
real(rprec) :: nu_a, nu_b, nu_c, nu_r, const1, const2
real(rprec), dimension(nx,ny) :: Nu_up
integer :: jx, jy, jz, jz_min, jz_max

! Set constants
const1 = (dt/dz)
const2 = (1._rprec/dz)

! Get the RHS ready
! Initialize with the explicit terms
Rz(:,:,1:nz-1) = w(:,:,1:nz-1) +                                     &
    dt * ( tadv1 * RHSz(:,:,1:nz-1) + tadv2 * RHSz_f(:,:,1:nz-1) )
if (coord == nproc-1) then
    Rz(:,:,nz) = w(:,:,nz) +                                         &
        dt * ( tadv1 * RHSz(:,:,nz) + tadv2 * RHSz_f(:,:,nz) )
endif

! Add explicit portion of Crank-Nicolson
call ddz_uv(tzz, dtzzdz_rhs, lbz)
dtzzdz_rhs(ld-1:ld, :, 1:nz-1) = 0._rprec
#ifdef PPSAFETYMODE
#ifdef PPMPI
dtzzdz_rhs(:,:,0) = BOGUS
#endif
#endif
! Assuming fixed dt here!
Rz(:,:,1:nz-1) = Rz(:,:,1:nz-1) - dt * 0.5_rprec * dtzzdz_rhs(:,:,1:nz-1)
! Not adding dtzzdz to Rz at top wall, assumed zero as done in divstress_w.f90

! Unlike uv, w has one less equation on bottom coord
#ifdef PPSAFETYMODE
if (coord == 0) then
    a(:,:,1) = BOGUS
    b(:,:,1) = BOGUS
    c(:,:,1) = BOGUS
endif
#endif

! Get bottom row, at jz = 2 instead of jz = 1
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,2) = BOGUS
#endif
    select case (lbc_mom)
        ! Stress free
        case (0)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    nu_c = 0.5_rprec*(Nu_t(jx,jy,3)+Nu_t(jx,jy,2)) + nu
                    ! Nu_t(jx,jy,1) on w, need to interpolate
#ifdef PPMAPPING
                    nu_b = (nu_c/JACO2(2)) + ((0.5_rprec*(Nu_t(jx,jy,2)+Nu_t(jx,jy,1)) + nu)/JACO2(1))
#else
                    nu_b = nu_c + 0.5_rprec*(Nu_t(jx,jy,2)+Nu_t(jx,jy,1)) + nu
#endif
                else
#ifdef PPMAPPING
                    nu_b = (nu/JACO2(2)) + (nu/JACO2(1))
#else
                    nu_b = 2.0_rprec*nu
#endif
                    nu_c = nu
                end if
#ifdef PPMAPPING
                b(jx,jy,2) = 1._rprec + const1*(1._rprec/JACO1(2))*const2*nu_b
                c(jx,jy,2) = -const1*(1._rprec/JACO1(2))*const2*(1._rprec/JACO2(2))*nu_c
#else
                b(jx,jy,2) = 1._rprec + const1*const2*nu_b
                c(jx,jy,2) = -const1*const2*nu_c
#endif
            end do
            end do

        ! DNS BC of Wall-model
        case (1:)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on uvp, no need to interpolate
#ifdef PPMAPPING
                    nu_b = (nu_c/JACO2(2)) + ((Nu_t(jx,jy,1)+nu)/JACO2(1))
#else
                    nu_b = nu_c + Nu_t(jx,jy,1) + nu
#endif
                    nu_c = 0.5_rprec*(Nu_t(jx,jy,3)+Nu_t(jx,jy,2)) + nu
                else
#ifdef PPMAPPING
                    nu_b = (nu/JACO2(2)) + (nu/JACO2(1))
#else
                    nu_b = 2.0_rprec*nu
#endif
                    nu_c = nu
                end if
#ifdef PPMAPPING
                b(jx,jy,2) = 1._rprec + const1*(1._rprec/JACO1(2))*const2*nu_b
                c(jx,jy,2) = -const1*(1._rprec/JACO1(2))*const2*(1._rprec/JACO2(2))*nu_c
#else
                b(jx,jy,2) = 1._rprec + const1*const2*nu_b
                c(jx,jy,2) = -const1*const2*nu_c
#endif
            end do
            end do
    end select
    jz_min = 3
else
    jz_min = 2
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
                    ! Nu_t(jx,jy,nz) on w, need to interpolate
                    nu_a = Nu_t(jx,jy,nz-1) + nu
                else
                    nu_a = nu
                end if
                ! txz(jx,jy,nz) = 0, so nothing added to RHS
#ifdef PPMAPPING
                a(jx,jy,nz-1) = -const1*(1._rprec/JACO2(nz-1))*const2*(1._rprec/JACO1(nz-1))*nu_a
                b(jx,jy,nz-1) = 1._rprec +                               &
                    const1*(1._rprec/JACO2(nz-1))*const2*(1._rprec/JACO1(nz-1))*nu_a
#else
                a(jx,jy,nz-1) = -const1*const2*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*const2*nu_a
#endif
            end do
            end do

        ! DNS BC or Wall-model
        case (1:)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on uvp, not needed here
                    nu_a = Nu_t(jx,jy,nz-1) + nu
                else
                    nu_a = nu
                endif
                ! Discretized txz(jx,jy,nz) as in wallstress,
                ! Therefore BC treated implicitly
                ! vtop = 0 so no changes to Ry
#ifdef PPMAPPING
                a(jx,jy,nz-1) = -const1*(1._rprec/JACO2(nz-1))*const2*(1._rprec/JACO1(nz-1))*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/JACO2(nz-1))*          &
                    (const2*(1._rprec/JACO1(2))*nu_a + (nu/(L_z-mesh_stretch(nz-1))))
                Rz(jx,jy,nz-1) = Rz(jx,jy,nz-1) + const1*(1._rprec/JACO2(1))*      &
                    (nu/(L_z-mesh_stretch(nz-1)))*utop
#else
                a(jx,jy,nz-1) = -const1*const2*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*(const2*nu_a + const3*nu)
                Rz(jx,jy,nz-1) = Rz(jx,jy,nz-1) + const1*const3*nu*utop
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

! Compute coefficients in domain near bottom interface
! This involves passing SGS viscosity from down coord
#ifdef PPMPI
! Pass SGS info up
if (sgs) then
    call mpi_sendrecv ( Nu_t(1,1,nz-1), nx*ny, MPI_RPREC, up, 9,            &
        Nu_up(1, 1), nx*ny, MPI_RPREC, down, 9, comm, status, ierr)
endif

! Now compute coefficients
if (coord > 0) then
    do jy = 1, ny
    do jx = 1, nx
        if (sgs) then
            ! Interpolate eddy viscosity onto uv-grid
            nu_a = 0.5_rprec*(Nu_t(jx,jy,1) + Nu_up(jx,jy)) + nu
            nu_c = 0.5_rprec*(Nu_t(jx,jy,2) + Nu_t(jx,jy,1)) + nu
#ifdef PPMAPPING
            nu_b = (nu_c/JACO2(1)) + (nu_a/JACO2(0))
#else
            nu_b = nu_a + nu_c
#endif
        else
            nu_a = nu
#ifdef PPMAPPING
            nu_b = (nu/JACO2(1)) + (nu/JACO2(0))
#else
            nu_b = 2._rprec*nu
#endif
            nu_c = nu
        endif

#ifdef PPMAPPING
        a(jx,jy,1) = -const1*(1._rprec/JACO1(1))*const2*(1._rprec/JACO2(0))*nu_a
        b(jx,jy,1) = 1._rprec + const1*(1._rprec/JACO1(1))*const2*nu_b
        c(jx,jy,1) = -const1*(1._rprec/JACO1(1))*const2*(1._rprec/JACO2(1))*nu_a
#else
        a(jx,jy,1) = -const1*const2*nu_a
        b(jx,jy,1) = 1._rprec + const1*const2*nu_b
        c(jx,jy,1) = -const1*const2*nu_c
#endif
    end do
    end do
endif
#endif

! Compute coefficients in domain
! Treating SGS viscosity explicitly!
do jz = jz_min, jz_max
do jy = 1, ny
do jx = 1, nx
    if (sgs) then
        ! Interpolate eddy viscosity onto uv-grid
        nu_a = 0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz-1)) + nu
        nu_c = 0.5_rprec*(Nu_t(jx,jy,jz+1) + Nu_t(jx,jy,jz)) + nu
#ifdef PPMAPPING
        nu_b = (nu_c/JACO2(jz)) + (nu_a/JACO2(jz-1))
#else
        nu_b = nu_a + nu_c
#endif
    else
        nu_a = nu
#ifdef PPMAPPING
        nu_b = (nu/JACO2(jz)) + (nu/JACO2(jz-1))
#else
        nu_b = 2._rprec*nu
#endif
        nu_c = nu
    endif

#ifdef PPMAPPING
    a(jx,jy,jz) = -const1*(1._rprec/JACO1(jz))*const2*(1._rprec/JACO2(jz-1))*nu_a
    b(jx,jy,jz) = 1._rprec + const1*(1._rprec/JACO1(jz))*const2*nu_b
    c(jx,jy,jz) = -const1*(1._rprec/JACO1(jz))*const2*(1._rprec/JACO2(jz))*nu_a
#else
    a(jx,jy,jz) = -const1*const2*nu_a
    b(jx,jy,jz) = 1._rprec + const1*const2*nu_b
    c(jx,jy,jz) = -const1*const2*nu_c
#endif
end do
end do
end do

! Find intermediate velocity in TDMA
call tridag_array_diff_w (a, b, c, Rz, wsol)

! Fill velocity solution
if (coord == 0) then
    w(:,:,1) = w(:,:,1) + dt * ( tadv1 * RHSz(:,:,1) + tadv2 * RHSz_f(:,:,1) )
    w(:nx,:ny,2:nz-1) = wsol(:nx,:ny,2:nz-1)
else
    w(:nx,:ny,1:nz-1) = wsol(:nx,:ny,1:nz-1)
endif

if (coord == nproc-1) then
    w(:,:,nz) = w(:,:,nz) + dt * ( tadv1 * RHSz(:,:,nz) + tadv2 * RHSz_f(:,:,nz) )
endif

end subroutine diff_stag_array_w
