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

public get_max_cfl, get_cfl_dt

contains

!*******************************************************************************
function get_max_cfl() result(cfl)
!*******************************************************************************
!
! This function provides the value of the maximum CFL in the entire
! domain
!
use types, only : rprec
use param, only : dt, dx, dy, nx, ny, nz, fourier, nxp
use sim_param, only : u,v,w
use sim_param, only : uF, vF, wF
use grid_m

#ifdef PPMPI
use mpi
use param, only : ierr, MPI_RPREC
#endif

implicit none
real(rprec) :: cfl
real(rprec) :: cfl_u, cfl_v, cfl_w
real(rprec), dimension(1:nz-1) :: cfl_w_temp
real(rprec), pointer, dimension(:) :: zw
integer :: jz

#ifdef PPMPI
real(rprec) :: cfl_buf
#endif

! Nullify pointers
nullify(zw)

zw => grid % zw

if (fourier) then !! remember dx = L_x / nxp (if fourier=true)
    cfl_u = maxval( abs(uF(1:nxp,1:ny,1:nz-1)) ) / dx
    cfl_v = maxval( abs(vF(1:nxp,1:ny,1:nz-1)) ) / dy

    do jz = 1, (nz-1)
        cfl_w_temp(jz) = maxval( abs(wF(1:nxp,1:ny,jz)) ) / (zw(jz+1) - zw(jz))
    end do
else
    cfl_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
    cfl_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy

    do jz = 1, (nz-1)
        cfl_w_temp(jz) = maxval( abs(w(1:nx,1:ny,jz)) ) / (zw(jz+1) - zw(jz))
    end do
endif

nullify(zw)
cfl_w = maxval( cfl_w_temp(1:nz-1) )

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
use param, only : cfl, dx, dy, nx, ny, nz, fourier, nxp
use sim_param, only : u,v,w
use sim_param, only : uF, vF, wF
use grid_m

#ifdef PPMPI
use mpi
use param, only : ierr, MPI_RPREC
#endif

implicit none

real(rprec) :: dt
real(rprec), pointer, dimension(:) :: zw
real(rprec), dimension(1:nz-1) :: dt_inv_w_temp
integer :: jz

! dt inverse
real(rprec) :: dt_inv_u, dt_inv_v, dt_inv_w

#ifdef PPMPI
real(rprec) :: dt_buf
#endif

! Nullify pointer
nullify(zw)

zw => grid % zw

! Avoid division by computing max dt^-1
if (fourier) then
    dt_inv_u = maxval( abs(uF(1:nxp,1:ny,1:nz-1)) ) / dx
    dt_inv_v = maxval( abs(vF(1:nxp,1:ny,1:nz-1)) ) / dy

    do jz = 1, (nz-1)
        dt_inv_w_temp(jz) = maxval( abs(wF(1:nxp,1:ny,jz)) ) / (zw(jz+1) - zw(jz))
    end do
else
    dt_inv_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
    dt_inv_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy

    do jz = 1, (nz-1)
        dt_inv_w_temp(jz) = maxval( abs(w(1:nx,1:ny,jz)) ) / (zw(jz+1) - zw(jz))
    end do
endif

nullify(zw)

dt_inv_w = maxval( dt_inv_w_temp(1:nz-1) )

dt = cfl / maxval( (/ dt_inv_u, dt_inv_v, dt_inv_w /) )

#ifdef PPMPI
call mpi_allreduce(dt, dt_buf, 1, MPI_RPREC, MPI_MIN, MPI_COMM_WORLD, ierr)
dt = dt_buf
#endif

end function get_cfl_dt

end module cfl_util
