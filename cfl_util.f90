!!
!!  Copyright (C) 2010-2017  Johns Hopkins University
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
module cfl_util
!*******************************************************************************
!
! This module provides the subroutines/functions for getting CFL related
! quantities
!
save
private

public get_max_cfl, get_cfl_dt, get_max_visc

contains

!*******************************************************************************
function get_max_cfl() result(cfl)
!*******************************************************************************
!
! This function provides the value of the maximum CFL in the entire
! domain
!
use types, only : rprec
use param, only : dt, dx, dy, dz, nx, ny, nz, fourier, nxp
use sim_param, only : u, v, w
use sim_param, only : uF, vF, wF
#ifdef PPMAPPING
use sim_param, only : JACO1
#endif

#ifdef PPMPI
use mpi
use param, only : ierr, MPI_RPREC
#endif

implicit none
real(rprec) :: cfl
real(rprec) :: cfl_u, cfl_v, cfl_w
#ifdef PPMAPPING
real(rprec), dimension(1:nz-1) :: cfl_w_temp
integer :: jz
#endif

#ifdef PPMPI
real(rprec) :: cfl_buf
#endif

if (fourier) then !! remember dx = L_x / nxp (if fourier=true)
    cfl_u = maxval( abs(uF(1:nxp,1:ny,1:nz-1)) ) / dx
    cfl_v = maxval( abs(vF(1:nxp,1:ny,1:nz-1)) ) / dy
#ifdef PPMAPPING
    do jz = 1, (nz-1)
        cfl_w_temp(jz) = maxval( abs(wF(1:nxp,1:ny,1:nz-1)) ) / (JACO1(jz)*dz)
    enddo
#else
    cfl_w = maxval( abs(wF(1:nxp,1:ny,1:nz-1)) ) / dz
#endif
else
    cfl_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
    cfl_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy
#ifdef PPMAPPING
    do jz = 1, (nz-1)
        cfl_w_temp(jz) = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / (JACO1(jz)*dz)
    enddo
#else
    cfl_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz
#endif
endif

#ifdef PPMAPPING
cfl_w = maxval( cfl_w_temp(1:nz-1) )
#endif
cfl = dt * maxval( (/ cfl_u, cfl_v, cfl_w /) )

#ifdef PPMPI
call mpi_allreduce(cfl, cfl_buf, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)
cfl = cfl_buf
#endif

end function get_max_cfl

!*******************************************************************************
function get_cfl_dt() result(dt)
!*******************************************************************************
!
! This functions determines the maximum allowable time step based on the CFL
! value specified in the param module
!
use types, only : rprec
use param, only : cfl, dx, dy, dz, nx, ny, nz, fourier, nxp
use sim_param, only : u,v,w
use sim_param, only : uF, vF, wF
#ifdef PPMAPPING
use sim_param, only : JACO1
#endif

#ifdef PPMPI
use mpi
use param, only : ierr, MPI_RPREC
#endif

implicit none

real(rprec) :: dt

! dt inverse
real(rprec) :: dt_inv_u, dt_inv_v, dt_inv_w
#ifdef PPMAPPING
real(rprec), dimension(1:nz-1) :: dt_inv_w_temp
integer :: jz
#endif

#ifdef PPMPI
real(rprec) :: dt_buf
#endif

! Avoid division by computing max dt^-1
if (fourier) then
    dt_inv_u = maxval( abs(uF(1:nxp,1:ny,1:nz-1)) ) / dx
    dt_inv_v = maxval( abs(vF(1:nxp,1:ny,1:nz-1)) ) / dy
#ifdef PPMAPPING
    do jz = 1, (nz-1)
        dt_inv_w_temp = maxval( abs(wF(1:nxp,1:ny,1:nz-1)) ) / (JACO1(jz)*dz)
    enddo
#else
    dt_inv_w = maxval( abs(wF(1:nxp,1:ny,1:nz-1)) ) / dz
#endif
else
    dt_inv_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
    dt_inv_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy
#ifdef PPMAPPING
    do jz = 1, (nz-1)
        dt_inv_w_temp = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / (JACO1(jz)*dz)
    enddo
#else
    dt_inv_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz
#endif
endif

#ifdef PPMAPPING
dt_inv_w = maxval( dt_inv_w_temp(1:nz-1) )
#endif
dt = cfl / maxval( (/ dt_inv_u, dt_inv_v, dt_inv_w /) )

#ifdef PPMPI
call mpi_allreduce(dt, dt_buf, 1, MPI_RPREC, MPI_MIN, MPI_COMM_WORLD, ierr)
dt = dt_buf
#endif

end function get_cfl_dt

!*******************************************************************************
function get_max_visc() result(visc)
!*******************************************************************************
!
! This function provides the value of the maximum viscous stability term,
! (nu + nu_t)*dt/(dx^2)
!
use types, only : rprec
use param, only : dt, dx, dy, dz, nx, ny, nz, nu_molec
! use param, only : fourier, nxp
use param, only : coord, lbc_mom, sgs, molec
use sgs_param, only : Nu_t
! use sgs_param, only : Nu_tF
#ifdef PPMAPPING
use sim_param, only : JACO1, JACO2
#endif

#ifdef PPMPI
use mpi
use param, only : ierr, MPI_RPREC
#endif

implicit none
real(rprec) :: visc, visc_x, visc_y, visc_z, nu_eff
#ifdef PPMAPPING
real(rprec), dimension(1:nz-1) :: visc_z_temp
integer :: jz
#endif

#ifdef PPMPI
real(rprec) :: visc_buf
#endif

nu_eff = 0.0_rprec
if (sgs) then
! if (fourier) then
!    nu_eff = nu_eff + maxval( abs(Nu_tF(1:nxp,1:ny,1:nz-1)) )
! else
    nu_eff = nu_eff + maxval( abs(Nu_t(1:nx,1:ny,1:nz-1)) )
! endif
endif
if (molec) then
    nu_eff = nu_eff + nu_molec
endif

visc_x = nu_eff / (dx**2)
visc_y = nu_eff / (dy**2)
#ifdef PPMAPPING
do jz = 1, (nz-1)
    visc_z_temp(jz) = nu_eff / ((JACO1(jz)*dz)**2)
enddo
if ((coord == 0) .and. (lbc_mom > 0)) then !! Nu_t on uv-grid at wall
    visc_z_temp(1) = nu_eff / ((JACO2(jz)*dz)**2)
endif
!! not considering upper wall since at nz
#else
visc_z = nu_eff / (dz**2)
#endif

#ifdef PPMAPPING
visc_z = maxval( visc_z_temp(1:nz-1) )
#endif
visc = dt * maxval( (/ visc_x, visc_y, visc_z /) )

#ifdef PPMPI
call mpi_allreduce(visc, visc_buf, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)
visc = visc_buf
#endif

end function get_max_visc

end module cfl_util
