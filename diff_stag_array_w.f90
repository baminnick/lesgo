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
use sim_param, only : jaco_w, jaco_uv, mesh_stretch
#endif
use derivatives, only : dft_direct_forw_2d_n_yonlyC, &
    dft_direct_back_2d_n_yonlyC
#ifdef PPHYBRID
use param, only : nproc_rnl
use derivatives, only : wave2physF, phys2waveF
#endif

implicit none

real(rprec), dimension(ld,ny,0:nz) :: Rz
real(rprec), dimension(ld,ny,lbz:nz) :: dtzzdz_rhs
real(rprec) :: nu_a, nu_b, nu_c, nu_r, const1, const2
#ifdef PPHYBRID
real(rprec), dimension(nxp+2,ny,0:nz) :: RzF, wsol
real(rprec), dimension(nxp,ny,0:nz) :: a, b, c
real(rprec), dimension(nxp,ny) :: Nu_up
#else
real(rprec), dimension(ld,ny,0:nz) :: wsol
real(rprec), dimension(nx,ny,0:nz) :: a, b, c
real(rprec), dimension(nx,ny) :: Nu_up
#endif
integer :: jx, jy, jz, jz_min, jz_max, end_nx

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

#ifdef PPHYBRID
! Move RHS to physical space (kx,ky,z) --> (x,y,z) on RNL coord
if (fourier) then
    call wave2physF( Rz, RzF)
else
    ! non-RNL coords have Rz already in physical space, so no transform needed
    ! Still have to move over values, note: nxp = nx on these coords
    RzF = Rz
endif
end_nx = nxp
#else
end_nx = nx
#endif

! Unlike uv, w has one less equation on bottom coord
#ifdef PPSAFETYMODE
if (coord == 0) then
    a(:,:,1) = BOGUS
    b(:,:,1) = BOGUS
    c(:,:,1) = BOGUS
endif
#endif

! Transform eddy viscosity, kx, ky, z --> kx, y, z
! however only going to use the kx=0 mode, assuming Nu_t(y,z) only
if ((fourier) .and. (sgs)) then
    do jz = 1, nz
        call dft_direct_back_2d_n_yonlyC( Nu_t(:,:,jz) )
    enddo
endif

! Imposing no-penetration BC regardless of whether or not
! it is stress-free or a wall. However still need to select
! case lbc_mom/ubc_mom because the eddy viscosity is in 
! different locations otherwise.

! Get bottom row, at jz = 2 instead of jz = 1
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,2) = BOGUS
#endif
    select case (lbc_mom)
        ! Stress free
        case (0)
            do jy = 1, ny
            do jx = 1, end_nx
                if (sgs) then
                    if (fourier) then
                        nu_c = 0.5_rprec*(Nu_t(1,jy,3)+Nu_t(1,jy,2)) + nu
                        ! Nu_t(jx,jy,1) on w, need to interpolate
#ifdef PPMAPPING
                        nu_b = (nu_c/jaco_uv(2)) + ((0.5_rprec*(Nu_t(1,jy,2)+Nu_t(jx,jy,1)) + nu)/jaco_uv(1))
#else
                        nu_b = nu_c + 0.5_rprec*(Nu_t(1,jy,2)+Nu_t(1,jy,1)) + nu
#endif
                    else
                        nu_c = 0.5_rprec*(Nu_t(jx,jy,3)+Nu_t(jx,jy,2)) + nu
                        ! Nu_t(jx,jy,1) on w, need to interpolate
#ifdef PPMAPPING
                        nu_b = (nu_c/jaco_uv(2)) + ((0.5_rprec*(Nu_t(jx,jy,2)+Nu_t(jx,jy,1)) + nu)/jaco_uv(1))
#else
                        nu_b = nu_c + 0.5_rprec*(Nu_t(jx,jy,2)+Nu_t(jx,jy,1)) + nu
#endif
                    endif
                else
#ifdef PPMAPPING
                    nu_b = (nu/jaco_uv(2)) + (nu/jaco_uv(1))
#else
                    nu_b = 2.0_rprec*nu
#endif
                    nu_c = nu
                end if
#ifdef PPMAPPING
                b(jx,jy,2) = 1._rprec + const1*(1._rprec/jaco_w(2))*const2*nu_b
                c(jx,jy,2) = -const1*(1._rprec/jaco_w(2))*const2*(1._rprec/jaco_uv(2))*nu_c
#else
                b(jx,jy,2) = 1._rprec + const1*const2*nu_b
                c(jx,jy,2) = -const1*const2*nu_c
#endif
            end do
            end do

        ! DNS BC or Wall-model
        case (1:)
            do jy = 1, ny
            do jx = 1, end_nx
                if (sgs) then
                    if (fourier) then
                        nu_c = 0.5_rprec*(Nu_t(1,jy,3)+Nu_t(1,jy,2)) + nu
                        ! Nu_t(jx,jy,1) on uvp, no need to interpolate
#ifdef PPMAPPING
                        nu_b = (nu_c/jaco_uv(2)) + ((Nu_t(1,jy,1)+nu)/jaco_uv(1))
#else
                        nu_b = nu_c + Nu_t(1,jy,1) + nu
#endif
                    else
                        nu_c = 0.5_rprec*(Nu_t(jx,jy,3)+Nu_t(jx,jy,2)) + nu
                        ! Nu_t(jx,jy,1) on uvp, no need to interpolate
#ifdef PPMAPPING
                        nu_b = (nu_c/jaco_uv(2)) + ((Nu_t(jx,jy,1)+nu)/jaco_uv(1))
#else
                        nu_b = nu_c + Nu_t(jx,jy,1) + nu
#endif
                    endif
                else
#ifdef PPMAPPING
                    nu_b = (nu/jaco_uv(2)) + (nu/jaco_uv(1))
#else
                    nu_b = 2.0_rprec*nu
#endif
                    nu_c = nu
                end if
#ifdef PPMAPPING
                b(jx,jy,2) = 1._rprec + const1*(1._rprec/jaco_w(2))*const2*nu_b
                c(jx,jy,2) = -const1*(1._rprec/jaco_w(2))*const2*(1._rprec/jaco_uv(2))*nu_c
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
            do jx = 1, end_nx
                if (sgs) then
                    ! Nu_t(jx,jy,nz) on w, need to interpolate
                    if (fourier) then
                        nu_a = 0.5_rprec*(Nu_t(1,jy,nz-2) + Nu_t(1,jy,nz-1)) + nu
#ifdef PPMAPPING
                        nu_b = (nu_a/jaco_uv(nz-2)) + (( 0.5_rprec*(Nu_t(1,jy,nz)+Nu_t(1,jy,nz-1)) +nu)/jaco_uv(nz-1))
#else
                        nu_b = nu_a + 0.5_rprec*(Nu_t(1,jy,nz)+Nu_t(1,jy,nz-1)) + nu
#endif
                    else
                        nu_a = 0.5_rprec*(Nu_t(jx,jy,nz-2) + Nu_t(jx,jy,nz-1)) + nu
#ifdef PPMAPPING
                        nu_b = (nu_a/jaco_uv(nz-2)) + (( 0.5_rprec*(Nu_t(jx,jy,nz)+Nu_t(jx,jy,nz-1)) +nu)/jaco_uv(nz-1))
#else
                        nu_b = nu_a + 0.5_rprec*(Nu_t(jx,jy,nz)+Nu_t(jx,jy,nz-1)) + nu
#endif
                    endif
                else
                    nu_a = nu
#ifdef PPMAPPING
                    nu_b = (nu/jaco_uv(nz-2)) + (nu/jaco_uv(nz-1))
#else
                    nu_b = 2.0_rprec*nu
#endif
                end if
#ifdef PPMAPPING
                a(jx,jy,nz-1) = -const1*(1._rprec/jaco_w(nz-1))*const2*(1._rprec/jaco_uv(nz-2))*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/jaco_w(nz-1))*const2*nu_b
#else
                a(jx,jy,nz-1) = -const1*const2*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*const2*nu_b
#endif
            end do
            end do

        ! DNS BC or Wall-model
        case (1:)
            do jy = 1, ny
            do jx = 1, end_nx
                if (sgs) then
                    ! Nu_t(jx,jy,nz) on uvp point nz-1, no need to interpolate
                    if (fourier) then
                        nu_a = 0.5_rprec*(Nu_t(1,jy,nz-2) + Nu_t(1,jy,nz-1)) + nu
#ifdef PPMAPPING
                        nu_b = (nu_a/jaco_uv(nz-2)) + ((Nu_t(1,jy,nz)+nu)/jaco_uv(nz-1))
#else
                        nu_b = nu_a + Nu_t(1,jy,nz) + nu
#endif
                    else
                        nu_a = 0.5_rprec*(Nu_t(jx,jy,nz-2) + Nu_t(jx,jy,nz-1)) + nu
#ifdef PPMAPPING
                        nu_b = (nu_a/jaco_uv(nz-2)) + ((Nu_t(jx,jy,nz)+nu)/jaco_uv(nz-1))
#else
                        nu_b = nu_a + Nu_t(jx,jy,nz) + nu
#endif
                    endif
                else
                    nu_a = nu
#ifdef PPMAPPING
                    nu_b = (nu/jaco_uv(nz-2)) + (nu/jaco_uv(nz-1))
#else
                    nu_b = 2.0_rprec*nu
#endif
                endif
#ifdef PPMAPPING
                a(jx,jy,nz-1) = -const1*(1._rprec/jaco_w(nz-1))*const2*(1._rprec/jaco_uv(nz-2))*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/jaco_w(nz-1))*const2*nu_b
#else
                a(jx,jy,nz-1) = -const1*const2*nu_a
                b(jx,jy,nz-1) = 1._rprec + const1*const2*nu_b
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
#ifdef PPHYBRID
if (coord .ne. (nproc_rnl - 1)) then
    call mpi_sendrecv( Nu_t(1,1,nz-1), nxp*ny, MPI_RPREC, up, 9,             &
        Nu_up(1,1), nxp*ny, MPI_RPREC, down, 9, comm, status, ierr)
else !! RNL coord near interface
    Nu_up(:,:) = 0.0_rprec !! fill with zeros initially
    Nu_up(1,:) = Nu_t(1,:,nz-1) !! fill only kx=0
    call mpi_sendrecv(Nu_up, nxp*ny, MPI_RPREC, up, 9,                      &
        Nu_up, nxp*ny, MPI_RPREC, down, 9, comm, status, ierr)
endif
#else
    call mpi_sendrecv ( Nu_t(1,1,nz-1), nx*ny, MPI_RPREC, up, 9,            &
        Nu_up(1, 1), nx*ny, MPI_RPREC, down, 9, comm, status, ierr)
#endif
endif

! Now compute coefficients
if (coord > 0) then
    do jy = 1, ny
    do jx = 1, end_nx
        if (sgs) then
            ! Interpolate eddy viscosity onto uv-grid
            if (fourier) then
                nu_a = 0.5_rprec*(Nu_t(1,jy,1) + Nu_up(1,jy)) + nu
                nu_c = 0.5_rprec*(Nu_t(1,jy,2) + Nu_t(1,jy,1)) + nu
            else
                nu_a = 0.5_rprec*(Nu_t(jx,jy,1) + Nu_up(jx,jy)) + nu
                nu_c = 0.5_rprec*(Nu_t(jx,jy,2) + Nu_t(jx,jy,1)) + nu
            endif
#ifdef PPMAPPING
            nu_b = (nu_c/jaco_uv(1)) + (nu_a/jaco_uv(0))
#else
            nu_b = nu_a + nu_c
#endif
        else
            nu_a = nu
#ifdef PPMAPPING
            nu_b = (nu/jaco_uv(1)) + (nu/jaco_uv(0))
#else
            nu_b = 2._rprec*nu
#endif
            nu_c = nu
        endif

#ifdef PPMAPPING
        a(jx,jy,1) = -const1*(1._rprec/jaco_w(1))*const2*(1._rprec/jaco_uv(0))*nu_a
        b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_w(1))*const2*nu_b
        c(jx,jy,1) = -const1*(1._rprec/jaco_w(1))*const2*(1._rprec/jaco_uv(1))*nu_a
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
do jx = 1, end_nx
    if (sgs) then
        ! Interpolate eddy viscosity onto uv-grid
        if (fourier) then
            nu_a = 0.5_rprec*(Nu_t(1,jy,jz) + Nu_t(1,jy,jz-1)) + nu
            nu_c = 0.5_rprec*(Nu_t(1,jy,jz+1) + Nu_t(1,jy,jz)) + nu
        else
            nu_a = 0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz-1)) + nu
            nu_c = 0.5_rprec*(Nu_t(jx,jy,jz+1) + Nu_t(jx,jy,jz)) + nu
        endif
#ifdef PPMAPPING
        nu_b = (nu_c/jaco_uv(jz)) + (nu_a/jaco_uv(jz-1))
#else
        nu_b = nu_a + nu_c
#endif
    else
        nu_a = nu
#ifdef PPMAPPING
        nu_b = (nu/jaco_uv(jz)) + (nu/jaco_uv(jz-1))
#else
        nu_b = 2._rprec*nu
#endif
        nu_c = nu
    endif

#ifdef PPMAPPING
    a(jx,jy,jz) = -const1*(1._rprec/jaco_w(jz))*const2*(1._rprec/jaco_uv(jz-1))*nu_a
    b(jx,jy,jz) = 1._rprec + const1*(1._rprec/jaco_w(jz))*const2*nu_b
    c(jx,jy,jz) = -const1*(1._rprec/jaco_w(jz))*const2*(1._rprec/jaco_uv(jz))*nu_a
#else
    a(jx,jy,jz) = -const1*const2*nu_a
    b(jx,jy,jz) = 1._rprec + const1*const2*nu_b
    c(jx,jy,jz) = -const1*const2*nu_c
#endif
end do
end do
end do

! Done using eddy viscosity, transform kx, y, z --> kx, ky, z
! Now transform Rz, kx, ky, z --> kx, y, z
if ((fourier) .and. (sgs)) then
    do jz = 1, nz
#ifdef PPHYBRID
        call dft_direct_forw_2d_n_yonlyC( Nu_t(:,:,jz) )
        ! No need to transform Rz, now using RzF
#else
        call dft_direct_forw_2d_n_yonlyC( Nu_t(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( Rz(:,:,jz) )
#endif
    enddo
endif

! Find intermediate velocity in TDMA
#ifdef PPHYBRID
call tridag_array_diff_w (a, b, c, RzF, wsol)
#else
call tridag_array_diff_w (a, b, c, Rz, wsol)
#endif

#ifdef PPHYBRID
if (fourier) then
    ! Transform back and omit non-RNL wavenumbers in fourier coord
    call phys2waveF( wsol, w )

    ! Overwrite level at bottom wall to enforce no-penetration
    if (coord == 0) then
        w(:,:,1) = w(:,:,1) + dt * ( tadv1 * RHSz(:,:,1) + tadv2 * RHSz_f(:,:,1) )
    endif
else
    ! No need to transfrom for physical coords, also nx = nxp here
    w(:nx,:ny,1:nz-1) = wsol(:nx,:ny,1:nz-1)   

    ! Overwrite level at top wall to enforce no-penetration
    if (coord == nproc-1) then
        w(:,:,nz) = w(:,:,nz) + dt * ( tadv1 * RHSz(:,:,nz) + tadv2 * RHSz_f(:,:,nz) )
    endif
endif
#else
! Since a,b,c(y,z) and Rz(kx,y,z) we have wsol(kx,y,z) therefore transform
! No need to change Rz back, will be overwritten in next time-step
if ((fourier) .and. (sgs)) then
    do jz = 1, nz-1
        call dft_direct_forw_2d_n_yonlyC( wsol(:,:,jz) )
    enddo
endif

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
#endif

end subroutine diff_stag_array_w
