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
module mfm
!******************************************************************************
!
! This module contains all of the subroutines associated with the macroscopic
! forcing method, including solving the generalized momentum transport (GMT)
!

use types, only : rprec
use param, only : path

implicit none

save
private

public :: gmt

! Main simulation variables for the GMT
real(rprec), dimension(:,:,:), allocatable :: gmx, gmy, gmz
real(rprec), dimension(:,:,:), allocatable :: gmtxx, gmtxy, gmtyy, &
    gmtxz, gmtyz, gmtzz, div_gmtx, div_gmty, div_gmtz
#ifdef PPCNDIFF
real(rprec), dimension(:,:,:), allocatable :: gmtxz_half1, gmtxz_half2
real(rprec), dimension(:,:,:), allocatable :: gmtyz_half1, gmtyz_half2
#endif
real(rprec), dimension(:,:,:), allocatable :: u_big, v_big, w_big

contains 

!******************************************************************************
subroutine mfm_init
!******************************************************************************
!
! This subroutine initializes the variables for the mfm module
!
use param, only : lbz, ld, nx, ny, nz, nxp

! Allocate simulations variables
allocate ( gmx(nxp+2, ny, lbz:nz) ); gmx = 0._rprec
allocate ( gmy(nxp+2, ny, lbz:nz) ); gmy = 0._rprec
allocate ( gmz(nxp+2, ny, lbz:nz) ); gmz = 0._rprec
allocate ( gmtxx(nxp+2, ny, lbz:nz) ); gmtxx = 0._rprec
allocate ( gmtxy(nxp+2, ny, lbz:nz) ); gmtxy = 0._rprec
allocate ( gmtyy(nxp+2, ny, lbz:nz) ); gmtyy = 0._rprec
allocate ( gmtxz(nxp+2, ny, lbz:nz) ); gmtxz = 0._rprec
allocate ( gmtyz(nxp+2, ny, lbz:nz) ); gmtyz = 0._rprec
allocate ( gmtzz(nxp+2, ny, lbz:nz) ); gmtzz = 0._rprec
allocate ( div_gmtx(nxp+2, ny, lbz:nz) ); div_gmtx = 0._rprec
allocate ( div_gmty(nxp+2, ny, lbz:nz) ); div_gmty = 0._rprec
allocate ( div_gmtz(nxp+2, ny, lbz:nz) ); div_gmtz = 0._rprec
#ifdef PPCNDIFF
allocate ( gmtxz_half1(nxp+2, ny, lbz:nz) ); gmtxz_half1 = 0.0_rprec
allocate ( gmtxz_half2(nxp+2, ny, lbz:nz) ); gmtxz_half2 = 0.0_rprec
allocate ( gmtyz_half1(nxp+2, ny, lbz:nz) ); gmtyz_half1 = 0.0_rprec
allocate ( gmtyz_half2(nxp+2, ny, lbz:nz) ); gmtyz_half2 = 0.0_rprec
#endif
allocate ( u_big(ld_big, ny2, lbz:nz) ); u_big = 0._rprec
allocate ( v_big(ld_big, ny2, lbz:nz) ); v_big = 0._rprec
allocate ( w_big(ld_big, ny2, lbz:nz) ); w_big = 0._rprec

end subroutine mfm_init

!******************************************************************************
subroutine ic_gmt
!******************************************************************************
!
! This subroutine initializes the initial velocity profile for the GMT equation
!

end subroutine ic_gmt

!******************************************************************************
subroutine mfm_checkpoint
!******************************************************************************
!
! This subroutine saves checkpoint variables for MFM analysis
!

end subroutine mfm_checkpoint

!******************************************************************************
subroutine gmt_wallstress
!******************************************************************************
! 
! This subroutine computed txz, tyz, dudz, dvdz at walls for the GMT
! 
use types, only : rprec
use param, only : lbc_mom, ubc_mom, coord, nproc, nz, ny, nx, dz, Lz
#ifdef PPMAPPING 
use sim_param, only : mesh_stretch
#endif
implicit none
real(rprec) :: denom
integer :: i, j

! Lower boundary condition
if (coord == 0) then

#ifdef PPMAPPING
    denom = mesh_stretch(1)
#else
    denom = 0.5_rprec*dz
#endif

    select case (lbc_mom)
        ! Stress free
        case (0)
            dgmxdz(:,:,1) = 0.0_rprec
            dgmydz(:,:,1) = 0.0_rprec
            gmtxz(:,:,1) = 0.0_rprec
            gmtyz(:,:,1) = 0.0_rprec

        ! DNS wall
        case (1)
            do j = 1, ny
            do i = 1, nx
                dgmxdz(i,j,1) = ( gmx(i,j,1) - ubot ) / denom
                dgmydz(i,j,1) = gmy(i,j,1) / denom
                gmtxz(i,j,1) = -nu_molec/(z_i*u_star)*dgmxdz(i,j,1)
                gmtyz(i,j,1) = -nu_molec/(z_i*u_star)*dgmydz(i,j,1)
            enddo
            enddo

    end select
endif

! Upper boundary condition
if (coord == nproc-1) then

#ifdef PPMAPPING
    denom = L_z - mesh_stretch(nz-1)
#else
    denom = 0.5_rprec*dz
#endif

    select case (ubc_mom)
        ! Stress free
        case (0)
            dgmxdz(:,:,nz) = 0.0_rprec
            dgmydz(:,:,nz) = 0.0_rprec
            gmtxz(:,:,nz) = 0.0_rprec
            gmtyz(:,:,nz) = 0.0_rprec

        ! DNS wall
        case (1)
            do j = 1, ny
            do i = 1, nx
                dgmxdz(i,j,nz) = ( utop - gmx(i,j,nz-1) ) / denom
                dgmydz(i,j,nz) = -gmy(i,j,nz-1) / denom
                gmtxz(i,j,nz) = -nu_molec/(z_i*u_star)*dgmxdz(i,j,nz)
                gmtyz(i,j,nz) = -nu_molec/(z_i*u_star)*dgmydz(i,j,nz)
            enddo
            enddo

    end select
endif

end subroutine gmt_wallstress

!******************************************************************************
subroutine to_big(a, a_big)
!******************************************************************************
! 
! Add padding to variable a for multiplication in physical space and dealiasing
! 
! Note: This subroutine is not yet ready for fourier mode
! 

use fft
use param, only : lbz, ny, nz, nxp

real(rprec), dimension(ld, ny, lbz:nz), intent(inout) :: a
real(rprec), dimension(ld_big, ny2, lbz:nz), intent(inout) :: a_big
real(rprec), dimension(ld, ny) :: temp

integer :: jz
real(rprec) :: const

const = 1._rprec/(nx*ny)

do jz = lbz, nz
    temp(:,:) = const*a(:,:,jz)
    call dfftw_execute_dft_r2c(forw, temp, temp)
    call padd(a_big(:,:,jz), temp)
    call dfftw_execute_dft_c2r(back_big, a_big(:,:,jz), a_big(:,:,jz))
enddo

end subroutine

!******************************************************************************
subroutine gmt
!******************************************************************************
!
! This subroutine solves the generalized momentum transport (GMT) equation
!
use param, only : nu_molec
use derivatives, only : filt_da, ddx, ddy, ddz_uv, ddz_w
use mpi

real(rprec) :: diff_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate derivatives of generalized momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dudx, dudy, dvdx, dvdy, dwdx, dwdy derivatives (in Fourier space)
call filt_da(gmx, dgmxdx, dgmxdy, lbz)
call filt_da(gmy, dgmydx, dgmydy, lbz)
call filt_da(gmz, dgmzdx, dgmzdy, lbz)

! dudz, dvdz using finite differences (for 1:nz on uv-nodes)
! except bottom coord, only 2:nz
call ddz_uv(gmx, dgmxdz, lbz)
call ddz_uv(gmy, dgmydz, lbz)

! dwdz using finite differences (for 0:nz-1 on w-nodes)
! except bottom coord, only 1:nz-1
call ddz_w(gmz, dgmzdz, lbz)

! Wall stress and derivatives at the wall
! (txz, tyz, dgmxdz, dgmydz at jz=1)
if (coord == 0) then
    call gmt_wallstress()
endif
if (coord == nproc-1) then
    call gmt_wallstress()
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Stress of generalized momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This part is written out explicitly instead of using sgs_stag from the
! main loop since sgs_stag is messy, however set up in a similar manner
diff_coef = -2.0_rprec*(nu_molec/(u_star*z_i))

! Calculate stress for bottom level
if (coord == 0) then
    gmtxx(:,:,1) = diff_coef*dgmxdx(:,:,1)
    gmtxy(:,:,1) = diff_coef*0.5_rprec*( dgmxdy(:,:,1) + dgmydx(:,:,1) )
    gmtyy(:,:,1) = diff_coef*dgmydy(:,:,1)
    gmtzz(:,:,1) = diff_coef*dgmzdz(:,:,1)
    ! Remember txz & tyz already calculated in gmt_wallstress

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
endif

! Unlike sgs_stag, calculate all tij for nz-1 on top coord with 
! all other levels isntead of separately

! Calculate stress for entire domain
do jz = jz_min, nz-1
    gmtxx(:,:,jz) = diff_coef*dgmxdx(:,:,jz)
    gmtxy(:,:,jz) = diff_coef*0.5_rprec*( dgmxdy(:,:,jz) + dgmydx(:,:,jz) )
    gmtyy(:,:,jz) = diff_coef*dgmydy(:,:,jz)
    gmtzz(:,:,jz) = diff_coef*dgmzdz(:,:,jz)
#ifdef PPCNDIFF
    gmtxz(:,:,jz) = diff_coef*0.5_rprec*( dgmxdz(:,:,jz) + dgmzdx(:,:,jz) )
    gmtxz_half1(:,:,jz) = diff_coef*0.5_rprec*( dgmzdx(:,:,jz) )
    gmtxz_half2(:,:,jz) = diff_coef*0.5_rprec*( dgmxdz(:,:,jz) )
    gmtyz(:,:,jz) = diff_coef*0.5_rprec*( dgmydz(:,:,jz) + dgmzdy(:,:,jz) )
    gmtyz_half1(:,:,jz) = diff_coef*0.5_rprec*( dgmzdy(:,:,jz) )
    gmtyz_half2(:,:,jz) = diff_coef*0.5_rprec*( dgmydz(:,:,jz) )
#else
    gmtxz(:,:,jz) = diff_coef*0.5_rprec*( dgmxdz(:,:,jz) + dgmzdx(:,:,jz) )
    gmtyz(:,:,jz) = diff_coef*0.5_rprec*( dgmydz(:,:,jz) + dgmzdy(:,:,jz) )
#endif
enddo

! Exchange information between processors to set values at nz 
! from jz=1 above to jz=nz below
#ifdef PPCNDIFF
call mpi_sync_real_array( gmtxz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( gmtxz_half1, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( gmtxz_half2, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( gmtyz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( gmtyz_half1, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( gmtyz_half2, 0, MPI_SYNC_DOWN )
#else
call mpi_sync_real_array( gmtxz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( gmtyz, 0, MPI_SYNC_DOWN )
#endif

! Exchange ghost node information for tzz
! send info up (from nz-1 below to 0 above)
#ifdef PPMPI
    call mpi_sendrecv (gmtzz(:,:,nz-1), (nxp+2)*ny, MPI_RPREC, up, 6,       &
        gmtzz(:,:,0), (nxp+2)*ny, MPI_RPREC, down, 6, comm, status, ierr)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diffusive term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate divergence of Stress for GMT
! Using the same divstress routines used in the main loop
#ifdef PPCNDIFF
call divstress_uv(div_gmtx, div_gmty, gmtxx, gmtxy, gmtxz_half1, gmtyy, gmtyz_half1)
call divstress_w_cndiff(div_gmtz, gmtxz, gmtyz, gmtzz)
#else
call divstress_uv(div_gmtx, div_gmty, gmtxx, gmtxy, gmtxz, gmtyy, gmtyz)
call divstress_w(div_gmtz, gmtxz, gmtyz, gmtzz)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advective term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Using simulation velocities (u,v,w) here to advect GMT

! Move variables to big domain
call to_big(u, u_big)
call to_big(v, v_big)
call to_big(w, w_big)

end subroutine gmt

end module mfm
