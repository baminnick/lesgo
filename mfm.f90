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

!*****************************************************************************
module mfm
!*****************************************************************************
!
! This module contains all of the subroutines associated with the macroscopic
! forcing method, including solving the generalized momentum transport (GMT)
!

use types, only : rprec
use param, only : path

implicit none

save
private

public :: mfm_init, ic_gmt, gm_transport, gm_2d3c_transport, mfm_checkpoint

! Main simulation variables for the GMT
real(rprec), public, dimension(:,:,:), allocatable :: gmu, gmv, gmw
real(rprec), dimension(:,:,:), allocatable :: dgmudx, dgmudy, dgmudz,   &
    dgmvdx, dgmvdy, dgmvdz, dgmwdx, dgmwdy, dgmwdz,                     &
    gmtxx, gmtxy, gmtyy, gmtxz, gmtyz, gmtzz,                           &
    div_gmtx, div_gmty, div_gmtz,                                       &
    rhs_gmx, rhs_gmy, rhs_gmz, rhs_gmx_f, rhs_gmy_f, rhs_gmz_f,         &
    gmp, dgmpdx, dgmpdy, dgmpdz
#ifdef PPCNDIFF
real(rprec), dimension(:,:,:), allocatable :: gmtxz_half1, gmtxz_half2
real(rprec), dimension(:,:,:), allocatable :: gmtyz_half1, gmtyz_half2
#endif
real(rprec), dimension(:,:,:), allocatable :: u_big, v_big, w_big,      &
    dgmudx_big, dgmudy_big, dgmudz_big,                                 & 
    dgmvdx_big, dgmvdy_big, dgmvdz_big,                                 & 
    dgmwdx_big, dgmwdy_big, dgmwdz_big,                                 & 
    temp_big

! Whether to initialize gmt field
logical, public :: init_gmt = .true.
! Name of file for restarting
character(64) :: fname

! Initial Condition / Forcing Method
integer, public :: ic_mfm = 1

! Boundary conditions for ic_mfm = 1
real(rprec), public :: gmu_bot = 0._rprec
real(rprec), public :: gmu_top = 0._rprec
real(rprec), public :: gmv_bot = 0._rprec
real(rprec), public :: gmv_top = 0._rprec

! Step location for ic_mfm = 2 brute force method
real(rprec), public :: bf_loc = 1._rprec

! Initial noise
real(rprec), public :: initial_noise_gmt = 0._rprec

! Whether to advect GMT field with full velocity field or streamwise mean
logical, public :: total_advec = .true.

! Target Macroscopic field
real(rprec), dimension(:), allocatable :: gmu_bar, gmv_bar, gmw_bar

contains 

!*****************************************************************************
subroutine mfm_init
!*****************************************************************************
!
! This subroutine initializes the variables for the mfm module
!
use param

! Allocate simulations variables
allocate ( gmu(nxp+2, ny, lbz:nz) ); gmu = 0._rprec
allocate ( gmv(nxp+2, ny, lbz:nz) ); gmv = 0._rprec
allocate ( gmw(nxp+2, ny, lbz:nz) ); gmw = 0._rprec
allocate ( dgmudx(nxp+2, ny, lbz:nz) ); dgmudx = 0._rprec
allocate ( dgmudy(nxp+2, ny, lbz:nz) ); dgmudy = 0._rprec
allocate ( dgmudz(nxp+2, ny, lbz:nz) ); dgmudz = 0._rprec
allocate ( dgmvdx(nxp+2, ny, lbz:nz) ); dgmvdx = 0._rprec
allocate ( dgmvdy(nxp+2, ny, lbz:nz) ); dgmvdy = 0._rprec
allocate ( dgmvdz(nxp+2, ny, lbz:nz) ); dgmvdz = 0._rprec
allocate ( dgmwdx(nxp+2, ny, lbz:nz) ); dgmwdx = 0._rprec
allocate ( dgmwdy(nxp+2, ny, lbz:nz) ); dgmwdy = 0._rprec
allocate ( dgmwdz(nxp+2, ny, lbz:nz) ); dgmwdz = 0._rprec
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
allocate ( rhs_gmx(nxp+2, ny, lbz:nz) ); rhs_gmx = 0.0_rprec
allocate ( rhs_gmy(nxp+2, ny, lbz:nz) ); rhs_gmy = 0.0_rprec
allocate ( rhs_gmz(nxp+2, ny, lbz:nz) ); rhs_gmz = 0.0_rprec
allocate ( rhs_gmx_f(nxp+2, ny, lbz:nz) ); rhs_gmx_f = 0.0_rprec
allocate ( rhs_gmy_f(nxp+2, ny, lbz:nz) ); rhs_gmy_f = 0.0_rprec
allocate ( rhs_gmz_f(nxp+2, ny, lbz:nz) ); rhs_gmz_f = 0.0_rprec
allocate ( gmp(nxp+2, ny, lbz:nz) ); gmp = 0._rprec
allocate ( dgmpdx(nxp+2, ny, nz) ); dgmpdx = 0._rprec
allocate ( dgmpdy(nxp+2, ny, nz) ); dgmpdy = 0._rprec
allocate ( dgmpdz(nxp+2, ny, nz) ); dgmpdz = 0._rprec

! Big variables
! Remember nx2 = 3*nx/2, lh_big = nx2/2 + 1, ld_big= 2*lh_big
! Therefore big variables using nxp are size 3*nxp/2 + 2
allocate ( u_big(3*nxp/2 + 2, ny2, lbz:nz) ); u_big = 0._rprec
allocate ( v_big(3*nxp/2 + 2, ny2, lbz:nz) ); v_big = 0._rprec
allocate ( w_big(3*nxp/2 + 2, ny2, lbz:nz) ); w_big = 0._rprec
allocate ( dgmudx_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmudx_big = 0._rprec
allocate ( dgmudy_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmudy_big = 0._rprec
allocate ( dgmudz_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmudz_big = 0._rprec
allocate ( dgmvdx_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmvdx_big = 0._rprec
allocate ( dgmvdy_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmvdy_big = 0._rprec
allocate ( dgmvdz_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmvdz_big = 0._rprec
allocate ( dgmwdx_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmwdx_big = 0._rprec
allocate ( dgmwdy_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmwdy_big = 0._rprec
allocate ( dgmwdz_big(3*nxp/2 + 2, ny2, lbz:nz) ); dgmwdz_big = 0._rprec
allocate ( temp_big(3*nxp/2 + 2, ny2, lbz:nz) ); temp_big = 0._rprec

! MFM variables
allocate ( gmu_bar(nz) ); gmu_bar = 0._rprec
allocate ( gmv_bar(nz) ); gmv_bar = 0._rprec
allocate ( gmw_bar(nz) ); gmw_bar = 0._rprec

end subroutine mfm_init

!*****************************************************************************
subroutine ic_gmt
!*****************************************************************************
! Set initial profile for GMT equation, determine if there is a file to read
! in or a new profile needs to be generated
use param, only : coord, lbc_mom, BOGUS, lbz, nz
use string_util
use grid_m
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
#ifdef PPMAPPING
use sim_param, only : mesh_stretch
#endif

integer :: jz
real(rprec) :: z

fname = path // 'gmt.out'
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif
inquire (file=fname, exist=init_gmt)

! Initialize target macroscopic field
do jz = 1, nz
#ifdef PPMPI
#ifdef PPMAPPING
    z = mesh_stretch(jz)
#else
    z = (coord*(nz-1) + jz - 0.5_rprec) * dz
#endif
#else
    z = (jz - 0.5_rprec) * dz
#endif

    select case (ic_mfm)
        ! Linear profile in u, zero for v and w
        case (1)
        gmu_bar(jz) = z
        gmv_bar(jz) = 0.0_rprec
        gmw_bar(jz) = 0.0_rprec

        ! Heaviside profile in u, zero for v and w
        case (2)
        if (z < bf_loc) then
            gmu_bar(jz) = 0.0_rprec
        else
            gmu_bar(jz) = 1.0_rprec
        endif
        gmv_bar(jz) = 0.0_rprec
        gmw_bar(jz) = 0.0_rprec

    end select
end do

! Initialize GMT field
if (init_gmt) then
    if (coord == 0) write(*,*) "--> Reading initial GMT field from file"
    call ic_gmt_file
else
    if (coord == 0) write(*,*) "--> Creating initial GMT field by prescribed macroscopic forcing"
    call ic_gmt_mfm
! Uncomment this for debugging purposes... to initialize the same as initial.f90
!    if (coord == 0) write(*,*) "--> Creating initial GMT field as blended profile"
!    call ic_gmt_blend
endif

#ifdef PPMPI
! Exchange ghost node information for u, v, and w
call mpi_sync_real_array( gmu, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( gmv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( gmw, 0, MPI_SYNC_DOWNUP )

!--set 0-level velocities to BOGUS
if (coord == 0) then
    gmu(:, :, lbz) = BOGUS
    gmv(:, :, lbz) = BOGUS
    gmw(:, :, lbz) = BOGUS
end if
#endif

end subroutine ic_gmt

!*****************************************************************************
subroutine ic_gmt_file
!*****************************************************************************
! Read initial profile for GMT equation from a file
use param, only : nz, read_endian
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use grid_m

open(12, file=fname, form='unformatted', convert=read_endian)
read(12) gmu(:,:,1:nz), gmv(:,:,1:nz), gmw(:,:,1:nz),                 &
    rhs_gmx(:,:,1:nz), rhs_gmy(:,:,1:nz), rhs_gmz(:,:,1:nz)
close(12)

#ifdef PPMPI
call mpi_sync_real_array(gmu, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(gmv, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(gmw, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(rhs_gmx, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(rhs_gmy, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(rhs_gmz, 0, MPI_SYNC_DOWNUP)
#endif

end subroutine ic_gmt_file

!*****************************************************************************
subroutine ic_gmt_mfm
!*****************************************************************************
!
! This subroutine initializes the initial generalized momentum profile to 
! achieve a specified macroscopic forcing
!
use param
#ifdef PPMAPPING
use sim_param, only : mesh_stretch
#endif
implicit none
integer :: jz, jz_abs
real(rprec) :: z
real(rprec) :: wall_noise, decay, rms, sigma_rv

sigma_rv = 0.289_rprec
wall_noise = 10 !! dictates how strong the noise is at the wall
! The higher wall_noise is, the stronger the noise is

! Fill u, v, and w with uniformly distributed random numbers between 0 and 1
call init_random_seed
call random_number(gmu)
call random_number(gmv)
call random_number(gmw)

! Center random number about 0 and rescale
gmu = initial_noise_gmt / sigma_rv * (gmu - 0.5_rprec)
gmv = initial_noise_gmt / sigma_rv * (gmv - 0.5_rprec)
gmw = initial_noise_gmt / sigma_rv * (gmw - 0.5_rprec)

! Rescale noise depending on distance from wall and mean log profile
! z is in meters
do jz = 1, nz
#ifdef PPMPI
    jz_abs = coord * (nz-1) + jz
#ifdef PPMAPPING
    z = mesh_stretch(jz)
#else
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
#endif
#else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
#endif
    ! Change z-value for different boundary conditions
    ! For channel flow, choose closest wall
    if(lbc_mom > 0 .and. ubc_mom > 0) z = min(z, L_z - z)
    ! For half-channel (lbc_mom > 0 .and. ubc_mom .eq 0) leave as is

    decay = (1-exp(-wall_noise*z))/(1-exp(-wall_noise))
    gmu(:,:,jz) = gmu(:,:,jz) * decay + gmu_bar(jz)
    gmv(:,:,jz) = gmv(:,:,jz) * decay + gmv_bar(jz)
    gmw(:,:,jz) = gmw(:,:,jz) * decay + gmw_bar(jz)
end do

! Bottom boundary conditions
if (coord == 0) then
    gmw(:, :, 1) = 0._rprec !! Always no-penetration
#ifdef PPMPI
    gmu(:, :, 0) = gmu_bot
    gmv(:, :, 0) = gmv_bot
    gmw(:, :, 0) = 0._rprec !! Always no-penetration
#endif
end if

! Set upper boundary condition as zero for u, v, and w
#ifdef PPMPI
if (coord == nproc-1) then
#endif
    gmw(1:nxp, 1:ny, nz) = 0._rprec !! Always no-penetration
    gmu(1:nxp, 1:ny, nz) = gmu_top
    gmv(1:nxp, 1:ny, nz) = gmv_top
#ifdef PPMPI
end if
#endif

end subroutine ic_gmt_mfm

!*****************************************************************************
subroutine ic_gmt_blend
!*****************************************************************************
!
! This subroutine initializes the initial velocity profile for the GMT 
! equation as a blended turbulent profile
!
use param
#ifdef PPMAPPING
use sim_param, only : mesh_stretch
#endif
implicit none
integer :: jz, jz_abs
real(rprec), dimension(nz) :: ubar
real(rprec) :: sigma_rv, z_plus, uturb, z, nu_eff
real(rprec) :: wall_noise, decay
real(rprec) :: kappa1, kappa2, beta

! Set fitted constants for initial velocity profile
kappa2 = 8.0_rprec
kappa1 = (1.0_rprec/vonk)*log(kappa2) + 5.2_rprec
beta = 2.0_rprec

do jz = 1, nz

#ifdef PPMPI
#ifdef PPMAPPING
    z = mesh_stretch(jz)
#else
    z = (coord*(nz-1) + jz - 0.5_rprec) * dz
#endif
#else
    z = (jz - 0.5_rprec) * dz
#endif

    ! Change z-value for different boundary conditions
    ! For channel flow, choose closest wall
    if(lbc_mom > 0 .and. ubc_mom > 0) z = min(z, L_z - z)
    ! For half-channel (lbc_mom > 0 .and. ubc_mom .eq 0) leave as is

    if (trigger) then
        nu_eff = nu_molec / trig_factor
    else
        nu_eff = nu_molec
    endif

    ! Old Blended profile 
    !z_plus = z * z_i * u_star / nu_eff !! ~ z^+
    !uturb = (1.0_rprec/vonk)*log(z_plus) + 5.0_rprec !! log-law
    !gam = - 0.01_rprec*(z_plus**4)/(1.0_rprec + 5.0_rprec*(z_plus))
    !ubar(jz) = exp(gam)*z_plus + exp(1.0_rprec/gam)*uturb

    z_plus = z * z_i * u_star / nu_eff !! effective z^+
    uturb = (1.0_rprec/vonk)*log(kappa2 + z_plus) + 5.2_rprec !! log-law
    ubar(jz) = uturb*((1.0_rprec + ((z_plus/kappa1)**(-beta)))**(-1.0_rprec/beta))

end do

sigma_rv = 0.289_rprec
wall_noise = 10 !! dictates how strong the noise is at the wall
! The higher wall_noise is, the stronger the noise is

! Fill u, v, and w with uniformly distributed random numbers between 0 and 1
call init_random_seed
call random_number(gmu)
call random_number(gmv)
call random_number(gmw)

! Center random number about 0 and rescale
gmu = initial_noise / sigma_rv * (gmu - 0.5_rprec)
gmv = initial_noise / sigma_rv * (gmv - 0.5_rprec)
gmw = initial_noise / sigma_rv * (gmw - 0.5_rprec)

! Rescale noise depending on distance from wall and mean log profile
! z is in meters
do jz = 1, nz
#ifdef PPMPI
    jz_abs = coord * (nz-1) + jz
#ifdef PPMAPPING
    z = mesh_stretch(jz)
#else
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
#endif
#else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
#endif
    ! Change z-value for different boundary conditions
    ! For channel flow, choose closest wall
    if(lbc_mom > 0 .and. ubc_mom > 0) z = min(z, L_z - z)
    ! For half-channel (lbc_mom > 0 .and. ubc_mom .eq 0) leave as is

    decay = (1-exp(-wall_noise*z))/(1-exp(-wall_noise))
    gmu(:,:,jz) = gmu(:,:,jz) * decay + ubar(jz)
    gmv(:,:,jz) = gmv(:,:,jz) * decay
    gmw(:,:,jz) = gmw(:,:,jz) * decay
end do

! Bottom boundary conditions
if (coord == 0) then
    gmw(:, :, 1) = 0._rprec
#ifdef PPMPI
    gmu(:, :, 0) = 0._rprec
    gmv(:, :, 0) = 0._rprec
    gmw(:, :, 0) = 0._rprec
#endif
end if

! Set upper boundary condition as zero for u, v, and w
#ifdef PPMPI
if (coord == nproc-1) then
#endif
    gmw(1:nxp, 1:ny, nz) = 0._rprec
    gmu(1:nxp, 1:ny, nz) = 0._rprec
    gmv(1:nxp, 1:ny, nz) = 0._rprec
#ifdef PPMPI
end if
#endif

end subroutine ic_gmt_blend

!*****************************************************************************
subroutine mfm_checkpoint
!*****************************************************************************
!
! This subroutine saves checkpoint variables for MFM analysis
!
use param, only : nz, write_endian

open(11, file=fname, form='unformatted', convert=write_endian,              &
    status='unknown', position='rewind')
write (11) gmu(:,:,1:nz), gmv(:,:,1:nz), gmw(:,:,1:nz),                     &
    rhs_gmx(:,:,1:nz), rhs_gmy(:,:,1:nz), rhs_gmz(:,:,1:nz)
close(11)

end subroutine mfm_checkpoint

!*****************************************************************************
subroutine gmt_wallstress
!*****************************************************************************
! 
! This subroutine computed txz, tyz, dudz, dvdz at walls for the GMT
! 
use types, only : rprec
use param, only : lbc_mom, ubc_mom, coord, nproc, nz, ny, nxp, dz, L_z
use param, only : nu_molec, u_star, z_i
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
    ! Dirichlet BCs
    do j = 1, ny
    do i = 1, nxp
        dgmudz(i,j,1) = ( gmu(i,j,1) - gmu_bot ) / denom
        dgmvdz(i,j,1) = ( gmv(i,j,1) - gmv_bot ) / denom
        gmtxz(i,j,1) = -nu_molec/(z_i*u_star)*dgmudz(i,j,1)
        gmtyz(i,j,1) = -nu_molec/(z_i*u_star)*dgmvdz(i,j,1)
    enddo
    enddo
endif

! Upper boundary condition
if (coord == nproc-1) then
#ifdef PPMAPPING
    denom = L_z - mesh_stretch(nz-1)
#else
    denom = 0.5_rprec*dz
#endif
    ! Dirichlet BCs
    do j = 1, ny
    do i = 1, nxp
        dgmudz(i,j,nz) = ( gmu_top - gmu(i,j,nz-1) ) / denom
        dgmvdz(i,j,nz) = ( gmv_top - gmv(i,j,nz-1) ) / denom
        gmtxz(i,j,nz) = -nu_molec/(z_i*u_star)*dgmudz(i,j,nz)
        gmtyz(i,j,nz) = -nu_molec/(z_i*u_star)*dgmvdz(i,j,nz)
    enddo
    enddo
endif

end subroutine gmt_wallstress

!*****************************************************************************
subroutine to_big_fourier(a, a_big)
!*****************************************************************************
! 
! Add padding to variable a for multiplication in physical space and 
! dealiasing
! 
! This is similar to the to_big subroutine, however input a should be in 
! kx, ky space, output a_big should be in kx, y space
! 

use fft
use param
use derivatives, only : dft_direct_back_2d_n_yonlyC_big

real(rprec), dimension(ld, ny, lbz:nz), intent(in) :: a
real(rprec), dimension(ld_big, ny2, lbz:nz), intent(inout) :: a_big

integer :: jz

do jz = lbz, nz
    ! padd a while in kx, ky space
    call padd(a_big(:,:,jz), a(:,:,jz))
    ! transform ky --> y, keep in kx space
    call dft_direct_back_2d_n_yonlyC_big( a_big(:,:,jz) )
enddo

end subroutine to_big_fourier

!*****************************************************************************
subroutine to_big(a, a_big)
!*****************************************************************************
! 
! Add padding to variable a for multiplication in physical space and 
! dealiasing
! 

use fft
use param

real(rprec), dimension(nxp+2, ny, lbz:nz), intent(in) :: a
real(rprec), dimension(3*nxp/2+2, ny2, lbz:nz), intent(inout) :: a_big
real(rprec), dimension(nxp+2, ny) :: temp

integer :: jz
real(rprec) :: const

const = 1._rprec/(nxp*ny)

do jz = lbz, nz
    temp(:,:) = const*a(:,:,jz)
    call dfftw_execute_dft_r2c(gmt_forw, temp, temp)
    call gmt_padd(a_big(:,:,jz), temp)
    call dfftw_execute_dft_c2r(gmt_back_big, a_big(:,:,jz), a_big(:,:,jz))
enddo

end subroutine to_big

!*****************************************************************************
subroutine to_small(a_big, a)
!*****************************************************************************
! 
! Undo padding to variable after multiplication in physical space
! 

use fft
use param, only : lbz, ny, nz, nxp

real(rprec), dimension(nxp+2, ny, lbz:nz), intent(inout) :: a
real(rprec), dimension(3*nxp/2+2, ny2, lbz:nz), intent(inout) :: a_big

integer :: jz

! Normalization constant for transforms is accounted for outside
! this routine

do jz = 1, nz-1
    call dfftw_execute_dft_r2c(gmt_forw_big, a_big(:,:,jz), a_big(:,:,jz))
    call gmt_unpadd(a(:,:,jz), a_big(:,:,jz))
    call dfftw_execute_dft_c2r(gmt_back, a(:,:,jz), a(:,:,jz))
enddo

end subroutine to_small

!*****************************************************************************
subroutine gm_transport
!*****************************************************************************
!
! This subroutine marches the generalized momentum transport (GMT) equation
! one time-step
!
use param
use sim_param
use derivatives, only : ddz_uv, ddz_w
use derivatives, only : wave2physF
#ifdef PPMPI
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
use forcing, only : project

real(rprec) :: diff_coef, rmsdivvel, const
integer :: jx, jy, jz, jz_min, jz_max, nxp2

real(rprec), dimension(nz) :: gmu_avg, gmv_avg, gmw_avg, du, dv, dw

real(rprec), dimension(nxp+2,ny,lbz:nz) :: dtxdx, dtydy, dtzdz,         &
    dtxdx2, dtydy2, dtzdz2, dtxdx3, dtydy3, dtzdz3

real(rprec), dimension(nxp+2,ny,lbz:nz) :: uF_avg, vF_avg, wF_avg

! Save previous timestep's RHS
rhs_gmx_f = rhs_gmx
rhs_gmy_f = rhs_gmy
rhs_gmz_f = rhs_gmz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate derivatives of generalized momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dudx, dudy, dvdx, dvdy, dwdx, dwdy derivatives (in Fourier space)
call gmt_filt_da(gmu, dgmudx, dgmudy, lbz)
call gmt_filt_da(gmv, dgmvdx, dgmvdy, lbz)
call gmt_filt_da(gmw, dgmwdx, dgmwdy, lbz)

! dudz, dvdz using finite differences (for 1:nz on uv-nodes)
! except bottom coord, only 2:nz
call ddz_uv(gmu, dgmudz, lbz)
call ddz_uv(gmv, dgmvdz, lbz)

! dwdz using finite differences (for 0:nz-1 on w-nodes)
! except bottom coord, only 1:nz-1
call ddz_w(gmw, dgmwdz, lbz)

! Wall stress and derivatives at the wall
! (txz, tyz, dgmudz, dgmvdz at jz=1)
if (coord == 0) then
    call gmt_wallstress()
#ifdef PPCNDIFF
    ! Add boundary condition to explicit portion
    gmtxz_half2(1:nxp,:,1) = gmtxz(1:nxp,:,1)
    gmtyz_half2(1:nxp,:,1) = gmtyz(1:nxp,:,1)
#endif
endif
if (coord == nproc-1) then
    call gmt_wallstress()
#ifdef PPCNDIFF
    ! Add boundary condition to explicit portion
    gmtxz_half2(1:nxp,:,nz) = gmtxz(1:nxp,:,nz)
    gmtyz_half2(1:nxp,:,nz) = gmtyz(1:nxp,:,nz)
#endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Stress of generalized momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This part is written out explicitly instead of using sgs_stag from the
! main loop since sgs_stag is messy, however set up in a similar manner
diff_coef = -2.0_rprec*(nu_molec/(u_star*z_i))

! Calculate stress for bottom level
if (coord == 0) then
    gmtxx(:,:,1) = diff_coef*dgmudx(:,:,1)
    gmtxy(:,:,1) = diff_coef*0.5_rprec*( dgmudy(:,:,1) + dgmvdx(:,:,1) )
    gmtyy(:,:,1) = diff_coef*dgmvdy(:,:,1)
    gmtzz(:,:,1) = diff_coef*dgmwdz(:,:,1)
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
    gmtxx(:,:,jz) = diff_coef*dgmudx(:,:,jz)
    gmtxy(:,:,jz) = diff_coef*0.5_rprec*( dgmudy(:,:,jz) + dgmvdx(:,:,jz) )
    gmtyy(:,:,jz) = diff_coef*dgmvdy(:,:,jz)
    gmtzz(:,:,jz) = diff_coef*dgmwdz(:,:,jz)
#ifdef PPCNDIFF
    gmtxz(:,:,jz) = diff_coef*0.5_rprec*( dgmudz(:,:,jz) + dgmwdx(:,:,jz) )
    gmtxz_half1(:,:,jz) = diff_coef*0.5_rprec*( dgmwdx(:,:,jz) )
    gmtxz_half2(:,:,jz) = diff_coef*0.5_rprec*( dgmudz(:,:,jz) )
    gmtyz(:,:,jz) = diff_coef*0.5_rprec*( dgmvdz(:,:,jz) + dgmwdy(:,:,jz) )
    gmtyz_half1(:,:,jz) = diff_coef*0.5_rprec*( dgmwdy(:,:,jz) )
    gmtyz_half2(:,:,jz) = diff_coef*0.5_rprec*( dgmvdz(:,:,jz) )
#else
    gmtxz(:,:,jz) = diff_coef*0.5_rprec*( dgmudz(:,:,jz) + dgmwdx(:,:,jz) )
    gmtyz(:,:,jz) = diff_coef*0.5_rprec*( dgmvdz(:,:,jz) + dgmwdy(:,:,jz) )
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diffusive term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate divergence of Stress for GMT
! Compute stress gradients
call gmt_ddx(gmtxx, dtxdx, lbz)
call gmt_ddy(gmtyy, dtydy2, lbz)
call gmt_ddxy(gmtxy, dtxdx2, dtydy, lbz)
call gmt_ddx(gmtxz, dtxdx3, lbz)
call gmt_ddy(gmtyz, dtydy3, lbz)
#ifdef PPCNDIFF
call ddz_w(gmtxz_half1, dtzdz, lbz)
call ddz_w(gmtyz_half1, dtzdz2, lbz)
#else
call ddz_w(gmtxz, dtzdz, lbz)
call ddz_w(gmtyz, dtzdz2, lbz)
call ddz_uv(gmtzz, dtzdz3, lbz)
#endif

! Take sum, remember only 1:nz-1 are valid
div_gmtx(:,:,1:nz-1) = dtxdx(:,:,1:nz-1) + dtydy(:,:,1:nz-1) + dtzdz(:,:,1:nz-1)
div_gmty(:,:,1:nz-1) = dtxdx2(:,:,1:nz-1) + dtydy2(:,:,1:nz-1) + dtzdz2(:,:,1:nz-1)
#ifdef PPCNDIFF
div_gmtz(:,:,1:nz) = dtxdx3(:,:,1:nz) + dtydy3(:,:,1:nz)
#else
! As in divstress_w, assume that dz(tzz)=0.0 at walls
if (coord == 0) dtzdz3(:,:,1) = 0._rprec
if (coord == nproc-1) dtzdz3(:,:,nz) = 0._rprec
div_gmtz(:,:,1:nz) = dtxdx3(:,:,1:nz) + dtydy3(:,:,1:nz) + dtzdz3(:,:,1:nz)
#endif

! Set ld-1 and ld oddballs to 0
div_gmtx(nxp+1:nxp+2,:,1:nz-1) = 0._rprec
div_gmty(nxp+1:nxp+2,:,1:nz-1) = 0._rprec
div_gmtz(nxp+1:nxp+2,:,1:nz-1) = 0._rprec

! Set boundary points to BOGUS
#ifdef PPSAFETYMODE
#ifdef PPMPI
div_gmtx(:,:,0) = BOGUS
div_gmty(:,:,0) = BOGUS
div_gmtz(:,:,0) = BOGUS
#endif
div_gmtx(:,:,nz) = BOGUS
div_gmty(:,:,nz) = BOGUS
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advective term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This differs from the convec subroutine in the main loop in two main ways:
! 1. Uses coupling between simulation velocities (u,v,w) and GMT velocities
! 2. Is not u x omega, instead computed as u[sim]*grad(u[gmt])

! Move variables to big domain
if (fourier) then
    ! Transform RNL velocities to physical domain, same grid size as GMT
    call wave2physF( u, uF )
    call wave2physF( v, vF )
    call wave2physF( w, wF )
    if (.not. total_advec) then !! using streamwise mean to advec
        ! Overwrite uF,vF,wF with streamwise mean, will be overwritten later
        call x_avg_mfm( uF )
        call x_avg_mfm( vF )
        call x_avg_mfm( wF )
    endif
    call to_big(uF, u_big)
    call to_big(vF, v_big)
    call to_big(wF, w_big)
else !! not fourier
    call to_big(u, u_big)
    call to_big(v, v_big)
    call to_big(w, w_big)
endif
call to_big(dgmudx, dgmudx_big)
call to_big(dgmudy, dgmudy_big)
call to_big(dgmudz, dgmudz_big)
call to_big(dgmvdx, dgmvdx_big)
call to_big(dgmvdy, dgmvdy_big)
call to_big(dgmvdz, dgmvdz_big)
call to_big(dgmwdx, dgmwdx_big)
call to_big(dgmwdy, dgmwdy_big)
call to_big(dgmwdz, dgmwdz_big)

! Normalization for FFTs
nxp2 = 3*nxp/2
const = 1._rprec/(nxp2*ny2)

! Compute advective term in x-GMT
! Interpolate w and dudz onto uv-grid
if (coord == 0) then
    ! Bottom wall take w(jz=1) = 0
    temp_big(:,:,1) = const*(u_big(:,:,1)*dgmudx_big(:,:,1) +      &
        v_big(:,:,1)*dgmudy_big(:,:,1) +                           &
        0.5_rprec*w_big(:,:,2)*dgmudz_big(:,:,2))
    jz_min = 2
else
    jz_min = 1
endif

if (coord == nproc-1) then
    ! Top wall take w(jz=nz) = 0
    temp_big(:,:,nz-1) = const*(u_big(:,:,nz-1)*dgmudx_big(:,:,nz-1) +   &
        v_big(:,:,nz-1)*dgmudy_big(:,:,nz-1) +                           &
        0.5_rprec*w_big(:,:,nz-1)*dgmudz_big(:,:,nz-1))
    jz_max = nz-2
else
    jz_max = nz-1
endif

! For entire domain
do jz = jz_min, jz_max
    temp_big(:,:,jz) = const*(u_big(:,:,jz)*dgmudx_big(:,:,jz) +      &
        v_big(:,:,jz)*dgmudy_big(:,:,jz) +                            &
        0.5_rprec*(w_big(:,:,jz+1)*dgmudz_big(:,:,jz+1) +             &
        w_big(:,:,jz)*dgmudz_big(:,:,jz)))
enddo

! Move temp_big into RHSx for GMT and make small
call to_small(temp_big, rhs_gmx)

! Compute advective term in y-GMT
! Interpolate w and dvdz onto uv-grid
if (coord == 0) then
    ! Bottom wall take w(jz=1) = 0
    temp_big(:,:,1) = const*(u_big(:,:,1)*dgmvdx_big(:,:,1) +      &
        v_big(:,:,1)*dgmvdy_big(:,:,1) +                           &
        0.5_rprec*w_big(:,:,2)*dgmvdz_big(:,:,2))
    jz_min = 2
else
    jz_min = 1
endif

if (coord == nproc-1) then
    ! Top wall take w(jz=nz) = 0
    temp_big(:,:,nz-1) = const*(u_big(:,:,nz-1)*dgmvdx_big(:,:,nz-1) +   &
        v_big(:,:,nz-1)*dgmvdy_big(:,:,nz-1) +                           &
        0.5_rprec*w_big(:,:,nz-1)*dgmvdz_big(:,:,nz-1))
    jz_max = nz-2
else
    jz_max = nz-1
endif

! For entire domain
do  jz = jz_min, jz_max
    temp_big(:,:,jz) = const*(u_big(:,:,jz)*dgmvdx_big(:,:,jz) +      &
        v_big(:,:,jz)*dgmvdy_big(:,:,jz) +                            &
        0.5_rprec*(w_big(:,:,jz+1)*dgmvdz_big(:,:,jz+1) +             &
        w_big(:,:,jz)*dgmvdz_big(:,:,jz)))
enddo

! Move temp_big into RHSy for GMT and make small
call to_small(temp_big, rhs_gmy)

! Compute advective term in z-GMT
! Interpolate u, v, and dwdz onto w-grid
if (coord == 0) then
    ! Bottom wall take w(jz=1) = dwdx(jz=1) = dwdy(jz=1) = 0
    temp_big(:,:,1) = 0._rprec
    jz_min = 2
else
    jz_min = 1
endif

if (coord == nproc-1) then
    ! Top wall take w(jz=nz) = dwdx(jz=1) = dwdy(jz=1) = 0
    temp_big(:,:,nz) = 0._rprec
    jz_max = nz-1
else
    jz_max = nz-1
endif

! For entire domain
do jz = jz_min, jz_max
    temp_big(:,:,jz) = const*(                                           &
        0.5_rprec*(u_big(:,:,jz)+u_big(:,:,jz-1))*dgmwdx_big(:,:,jz) +   &
        0.5_rprec*(v_big(:,:,jz)+v_big(:,:,jz-1))*dgmwdy_big(:,:,jz) +   &
        w_big(:,:,jz)*0.5_rprec*(dgmwdz_big(:,:,jz)+dgmwdz_big(:,:,jz-1)))
enddo

! Move temp_big into RHSz for GMT and make small
call to_small(temp_big, rhs_gmz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add terms to the RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add div-tau term 
rhs_gmx(:,:,1:nz-1) = -rhs_gmx(:,:,1:nz-1) - div_gmtx(:,:,1:nz-1)
rhs_gmy(:,:,1:nz-1) = -rhs_gmy(:,:,1:nz-1) - div_gmty(:,:,1:nz-1)
rhs_gmz(:,:,1:nz-1) = -rhs_gmz(:,:,1:nz-1) - div_gmtz(:,:,1:nz-1)
if (coord == nproc-1) rhs_gmz(:,:,nz) = -rhs_gmz(:,:,nz) - div_gmtz(:,:,nz)

! Add pressure forcing -- only use for debugging purposes
!if (use_mean_p_force) rhs_gmx(:,:,1:nz-1)=rhs_gmx(:,:,1:nz-1)+mean_p_force_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Intermediate Velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PPCNDIFF
call gmt_diff_stag_array_uv(gmu,gmv,rhs_gmx,rhs_gmy,rhs_gmx_f,rhs_gmy_f,    &
    gmtxz_half2,gmtyz_half2,gmtxz,gmtyz)
call gmt_diff_stag_array_w(gmw,rhs_gmz,rhs_gmz_f,gmtzz)
#else
gmu(:,:,1:nz-1) = gmu(:,:,1:nz-1) +                                         &
    dt * ( tadv1 * rhs_gmx(:,:,1:nz-1) + tadv2 * rhs_gmx_f(:,:,1:nz-1) )
gmv(:,:,1:nz-1) = gmv(:,:,1:nz-1) +                                         &
    dt * ( tadv1 * rhs_gmy(:,:,1:nz-1) + tadv2 * rhs_gmy_f(:,:,1:nz-1) )
gmw(:,:,1:nz-1) = gmw(:,:,1:nz-1) +                                         &
    dt * ( tadv1 * rhs_gmz(:,:,1:nz-1) + tadv2 * rhs_gmz_f(:,:,1:nz-1) )
if (coord == nproc-1) then
    gmw(:,:,nz) = gmw(:,:,nz) +                                             &
        dt * ( tadv1 * rhs_gmz(:,:,nz) + tadv2 * rhs_gmz_f(:,:,nz) )
endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inverse Macroscopic Forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute difference in current macroscopic field and target
! and update the intermediate velocity field with this difference
! 
! Adding this difference acts as the source term for the GMT and 
! ensures the macroscopic field is unchanging with time-step
do jz = 1, nz-1
    ! Initialize averages
    gmu_avg(jz) = 0.0_rprec
    gmv_avg(jz) = 0.0_rprec
    gmw_avg(jz) = 0.0_rprec

    ! Average in x and y
    do jx = 1, nxp
    do jy = 1, ny
        gmu_avg(jz) = gmu_avg(jz) + gmu(jx,jy,jz)
        gmv_avg(jz) = gmv_avg(jz) + gmv(jx,jy,jz)
        gmw_avg(jz) = gmw_avg(jz) + gmw(jx,jy,jz)
    end do
    end do
    gmu_avg(jz) = gmu_avg(jz) / (nxp*ny)
    gmv_avg(jz) = gmv_avg(jz) / (nxp*ny)
    gmw_avg(jz) = gmw_avg(jz) / (nxp*ny)

    ! Find difference
    du(jz) = gmu_bar(jz) - gmu_avg(jz)
    dv(jz) = gmv_bar(jz) - gmv_avg(jz)
    dw(jz) = gmw_bar(jz) - gmw_avg(jz)

    ! Update intermediate velocity
    do jx = 1, nxp
    do jy = 1, ny
        gmu(jx,jy,jz) = gmu(jx,jy,jz) + du(jz)
        gmv(jx,jy,jz) = gmv(jx,jy,jz) + dv(jz)
        gmw(jx,jy,jz) = gmw(jx,jy,jz) + dw(jz)
    end do
    end do
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pressure Solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call gmt_press_stag_array(gmu,gmv,gmw,div_gmtz,gmp,dgmpdx,dgmpdy,dgmpdz)

! Add pressure gradients to RHS variables (for next time step)
rhs_gmx(:,:,1:nz-1) = rhs_gmx(:,:,1:nz-1) - dgmpdx(:,:,1:nz-1)
rhs_gmy(:,:,1:nz-1) = rhs_gmy(:,:,1:nz-1) - dgmpdy(:,:,1:nz-1)
rhs_gmz(:,:,1:nz-1) = rhs_gmz(:,:,1:nz-1) - dgmpdz(:,:,1:nz-1)
if (coord == nproc-1) then
    rhs_gmz(:,:,nz) = rhs_gmz(:,:,nz) - dgmpdz(:,:,nz)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Projection step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call project(gmu,gmv,gmw,dgmpdx,dgmpdy,dgmpdz)

! Output updated information on GMT to screen
if (modulo (jt_total, wbase) == 0) then
    call check_gmt(gmu,gmv,gmw,du,dv,dw)
    call rmsdiv(dgmudx,dgmvdy,dgmwdz,rmsdivvel,nxp)

    if(coord == 0) then
        call write_gmt_tau_wall_bot()
        write(*,*)
        write(*,'(a,E15.7)') ' GM  Velocity divergence metric: ', rmsdivvel
        write(*,'(a)') '===================== GMT BOTTOM ======================='
        write(*,*) 'u: ', gmu(nxp/2,ny/2,1:2)
        write(*,*) 'v: ', gmv(nxp/2,ny/2,1:2)
        write(*,*) 'w: ', gmw(nxp/2,ny/2,1:2)
        write(*,'(a)') '========================================================'
    end if
    call mpi_barrier(comm, ierr)
    if(coord == nproc-1) then
        call write_gmt_tau_wall_top()
        write(*,'(a)') '====================== GMT TOP ========================='
        write(*,*) 'u: ', gmu(nxp/2,ny/2,nz-2:nz-1)
        write(*,*) 'v: ', gmv(nxp/2,ny/2,nz-2:nz-1)
        write(*,*) 'w: ', gmw(nxp/2,ny/2,nz-1:nz)
        write(*,'(a)') '========================================================'
    end if
    call mpi_barrier(comm, ierr)
end if

end subroutine gm_transport

!*****************************************************************************
subroutine gm_2d3c_transport
!*****************************************************************************
!
! This subroutine marches the generalized momentum transport (GMT) equation
! one time-step using the streamwise mean to advect. This simplifies the GMT
! equation to have only one velocity component, gmu. Therefore, gmv and gmw
! are unchanging and no poisson equation is solved.
! 
! This subroutine is used when total_advec = false
!
use param
use sim_param
use derivatives, only : ddz_uv, ddz_w
use derivatives, only : wave2physF
#ifdef PPMPI
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
use forcing, only : project

real(rprec) :: diff_coef, rmsdivvel, const
integer :: jx, jy, jz, jz_min, jz_max, nxp2

real(rprec), dimension(nz) :: gmu_avg, du

real(rprec), dimension(nxp+2,ny,lbz:nz) :: dtydy, dtzdz

real(rprec), dimension(nxp+2,ny,lbz:nz) :: uF_avg, vF_avg, wF_avg

! Save previous timestep's RHS
rhs_gmx_f = rhs_gmx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate derivatives of generalized momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dudy derivatives (in Fourier space)
call gmt_ddy(gmu, dgmudy, lbz)

! dudz, dvdz using finite differences (for 1:nz on uv-nodes)
! except bottom coord, only 2:nz
call ddz_uv(gmu, dgmudz, lbz)

! Wall stress and derivatives at the wall
! (txz, tyz, dgmudz, dgmvdz at jz=1)
if (coord == 0) then
    call gmt_wallstress()
#ifdef PPCNDIFF
    ! Add boundary condition to explicit portion
    gmtxz_half2(1:nxp,:,1) = gmtxz(1:nxp,:,1)
#endif
endif
if (coord == nproc-1) then
    call gmt_wallstress()
#ifdef PPCNDIFF
    ! Add boundary condition to explicit portion
    gmtxz_half2(1:nxp,:,nz) = gmtxz(1:nxp,:,nz)
#endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Stress of generalized momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This part is written out explicitly instead of using sgs_stag from the
! main loop since sgs_stag is messy, however set up in a similar manner
diff_coef = -2.0_rprec*(nu_molec/(u_star*z_i))

! Calculate stress for bottom level
if (coord == 0) then
    gmtxy(:,:,1) = diff_coef*0.5_rprec*( dgmudy(:,:,1) )
    ! Remember txz already calculated in gmt_wallstress

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
endif

! Unlike sgs_stag, calculate all tij for nz-1 on top coord with 
! all other levels isntead of separately

! Calculate stress for entire domain
do jz = jz_min, nz-1
    gmtxy(:,:,jz) = diff_coef*0.5_rprec*( dgmudy(:,:,jz) )
#ifdef PPCNDIFF
    gmtxz(:,:,jz) = diff_coef*0.5_rprec*( dgmudz(:,:,jz) )
    gmtxz_half2(:,:,jz) = diff_coef*0.5_rprec*( dgmudz(:,:,jz) )
#else
    gmtxz(:,:,jz) = diff_coef*0.5_rprec*( dgmudz(:,:,jz) )
#endif
enddo

! Exchange information between processors to set values at nz 
! from jz=1 above to jz=nz below
#ifdef PPCNDIFF
call mpi_sync_real_array( gmtxz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( gmtxz_half2, 0, MPI_SYNC_DOWN )
#else
call mpi_sync_real_array( gmtxz, 0, MPI_SYNC_DOWN )
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diffusive term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate divergence of Stress for GMT
! Compute stress gradients
call gmt_ddy(gmtxy, dtydy, lbz)
#ifdef PPCNDIFF
! Do nothing, gmtxz_half was the gmwdx derivative, which is now zero
! call ddz_w(gmtxz_half1, dtzdz, lbz)
dtzdz = 0.0_rprec
#else
call ddz_w(gmtxz, dtzdz, lbz)
#endif

! Take sum, remember only 1:nz-1 are valid
div_gmtx(:,:,1:nz-1) = dtydy(:,:,1:nz-1) + dtzdz(:,:,1:nz-1)

! Set ld-1 and ld oddballs to 0
div_gmtx(nxp+1:nxp+2,:,1:nz-1) = 0._rprec

! Set boundary points to BOGUS
#ifdef PPSAFETYMODE
#ifdef PPMPI
div_gmtx(:,:,0) = BOGUS
#endif
div_gmtx(:,:,nz) = BOGUS
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advective term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This differs from the convec subroutine in the main loop in two main ways:
! 1. Uses coupling between simulation velocities (u,v,w) and GMT velocities
! 2. Is not u x omega, instead computed as u[sim]*grad(u[gmt])

! Move variables to big domain
if (fourier) then 
    ! Transform RNL velocities to physical domain, same grid size as GMT
    call wave2physF( v, vF )
    call wave2physF( w, wF )
else
    ! v and vf should be the same size, just carry values over
    vF = v
    wF = w
endif
! Overwrite vF,wF with streamwise mean, will be overwritten later
call x_avg_mfm( vF )
call x_avg_mfm( wF )
call to_big(vF, v_big)
call to_big(wF, w_big)
call to_big(dgmudy, dgmudy_big)
call to_big(dgmudz, dgmudz_big)

! Normalization for FFTs
nxp2 = 3*nxp/2
const = 1._rprec/(nxp2*ny2)

! Compute advective term in x-GMT
! Interpolate w and dudz onto uv-grid
if (coord == 0) then
    ! Bottom wall take w(jz=1) = 0
    temp_big(:,:,1) = const*( v_big(:,:,1)*dgmudy_big(:,:,1) +          &
        0.5_rprec*w_big(:,:,2)*dgmudz_big(:,:,2))
    jz_min = 2
else
    jz_min = 1
endif

if (coord == nproc-1) then
    ! Top wall take w(jz=nz) = 0
    temp_big(:,:,nz-1) = const*( v_big(:,:,nz-1)*dgmudy_big(:,:,nz-1) + &
        0.5_rprec*w_big(:,:,nz-1)*dgmudz_big(:,:,nz-1))
    jz_max = nz-2
else
    jz_max = nz-1
endif

! For entire domain
do jz = jz_min, jz_max
    temp_big(:,:,jz) = const*( v_big(:,:,jz)*dgmudy_big(:,:,jz) +     &
        0.5_rprec*(w_big(:,:,jz+1)*dgmudz_big(:,:,jz+1) +             &
        w_big(:,:,jz)*dgmudz_big(:,:,jz)))
enddo

! Move temp_big into RHSx for GMT and make small
call to_small(temp_big, rhs_gmx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add terms to the RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add div-tau term 
rhs_gmx(:,:,1:nz-1) = -rhs_gmx(:,:,1:nz-1) - div_gmtx(:,:,1:nz-1)

! Add pressure forcing -- only use for debugging purposes
!if (use_mean_p_force) rhs_gmx(:,:,1:nz-1)=rhs_gmx(:,:,1:nz-1)+mean_p_force_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Intermediate Velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PPCNDIFF
call gmt_2d3c_diff_stag_array_uv(gmu,rhs_gmx,rhs_gmx_f,gmtxz_half2,gmtxz)
#else
gmu(:,:,1:nz-1) = gmu(:,:,1:nz-1) +                                         &
    dt * ( tadv1 * rhs_gmx(:,:,1:nz-1) + tadv2 * rhs_gmx_f(:,:,1:nz-1) )
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inverse Macroscopic Forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute difference in current macroscopic field and target
! and update the intermediate velocity field with this difference
! 
! Adding this difference acts as the source term for the GMT and 
! ensures the macroscopic field is unchanging with time-step
do jz = 1, nz-1
    ! Initialize averages
    gmu_avg(jz) = 0.0_rprec

    ! Average in x and y
    do jx = 1, nxp
    do jy = 1, ny
        gmu_avg(jz) = gmu_avg(jz) + gmu(jx,jy,jz)
    end do
    end do
    gmu_avg(jz) = gmu_avg(jz) / (nxp*ny)

    ! Find difference
    du(jz) = gmu_bar(jz) - gmu_avg(jz)

    ! Update intermediate velocity
    do jx = 1, nxp
    do jy = 1, ny
        gmu(jx,jy,jz) = gmu(jx,jy,jz) + du(jz)
    end do
    end do
enddo

! Output updated information on GMT to screen
if (modulo (jt_total, wbase) == 0) then
    call check_gmt_2d3c(gmu,du)

    if(coord == 0) then
        call write_gmt_tau_wall_bot()
        write(*,*)
        write(*,'(a,E15.7)') ' GM  Velocity divergence metric: ', rmsdivvel
        write(*,'(a)') '===================== GMT BOTTOM ======================='
        write(*,*) 'u: ', gmu(nxp/2,ny/2,1:2)
        write(*,'(a)') '========================================================'
    end if
    call mpi_barrier(comm, ierr)
    if(coord == nproc-1) then
        call write_gmt_tau_wall_top()
        write(*,'(a)') '====================== GMT TOP ========================='
        write(*,*) 'u: ', gmu(nxp/2,ny/2,nz-2:nz-1)
        write(*,'(a)') '========================================================'
    end if
    call mpi_barrier(comm, ierr)
end if

end subroutine gm_2d3c_transport

!*****************************************************************************
subroutine write_gmt_tau_wall_bot()
!*****************************************************************************
!
! Measure and write spatially average bottom wall stress
! 
use types, only : rprec
use param, only : jt_total, wbase, nxp, ny
implicit none
real(rprec) :: twall, txavg, tyavg
integer :: jx, jy

! -------------------------- Wall Stress Statistics --------------------------
txavg = 0._rprec
tyavg = 0._rprec
do jx = 1, nxp
do jy = 1, ny
    txavg = txavg + gmtxz(jx,jy,1)
    tyavg = tyavg + gmtyz(jx,jy,1)
end do
end do
txavg = txavg / (nxp*ny)
tyavg = tyavg / (nxp*ny)
twall = sqrt( (txavg)**2 + (tyavg)**2 )

! ------------------------------ Write to file -------------------------------
open(2,file=path // 'output/gmt_tau_wall_bot.dat', status='unknown',        &
    form='formatted', position='append')

! one time header output
if (jt_total==wbase) write(2,*) 'jt_total, txavg, tyavg, twall'

! continual time-related output
write(2,*) jt_total, txavg, tyavg, twall
close(2)

end subroutine write_gmt_tau_wall_bot

!*****************************************************************************
subroutine write_gmt_tau_wall_top()
!*****************************************************************************
!
! Measure and write spatially average top wall stress
! 
use types, only : rprec
use param, only : jt_total, wbase, nxp, ny, nz
implicit none
real(rprec) :: twall, txavg, tyavg
integer :: jx, jy

! -------------------------- Wall Stress Statistics --------------------------
txavg = 0._rprec
tyavg = 0._rprec
do jx = 1, nxp
do jy = 1, ny
    txavg = txavg + gmtxz(jx,jy,nz)
    tyavg = tyavg + gmtyz(jx,jy,nz)
end do
end do
txavg = txavg / (nxp*ny)
tyavg = tyavg / (nxp*ny)
twall = sqrt( (txavg)**2 + (tyavg)**2 )

! ------------------------------ Write to file -------------------------------
open(2,file=path // 'output/gmt_tau_wall_top.dat', status='unknown',        &
    form='formatted', position='append')

! one time header output
if (jt_total==wbase) write(2,*) 'jt_total, txavg, tyavg, twall'

! continual time-related output
write(2,*) jt_total, txavg, tyavg, twall
close(2)

end subroutine write_gmt_tau_wall_top

!*****************************************************************************
subroutine check_gmt(u,v,w,du,dv,dw)
!*****************************************************************************
! 
! Measure and write the average difference from the target macroscopic value.
! Remember, these differences were measured in the main gm_transport routine
! and were found by averaging in the x & y directions. Therefore, only 
! averaging in the z direction here to get a global estimate.
! 
! Added computation of the energy of the u, v, and w components
! 
use types, only : rprec
use param
use messages
implicit none
integer :: jx, jy, jz
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: u, v, w
real(rprec), dimension(nz), intent(in) :: du, dv, dw
real(rprec) :: uu, vv, ww, du_avg, dv_avg, dw_avg
#ifdef PPMPI
real(rprec) :: uu_global, vv_global, ww_global
real(rprec) :: du_global, dv_global, dw_global
#endif

! Initialize variables
uu = 0._rprec
vv = 0._rprec
ww = 0._rprec
du_avg = 0._rprec
dv_avg = 0._rprec
dw_avg = 0._rprec

! Compute spatial average of squared velocities
do jz = 1, nz-1
do jy = 1, ny
do jx = 1, nxp
    uu = uu + u(jx,jy,jz)**2
    vv = vv + v(jx,jy,jz)**2
    ww = ww + w(jx,jy,jz)**2
enddo
enddo
enddo
uu = uu / (nxp*ny*(nz-1))
vv = vv / (nxp*ny*(nz-1))
ww = ww / (nxp*ny*(nz-1))

! Perform spatial averaging for differences
! Remember these are already averaged in the x & y directions
do jz = 1, nz-1
    du_avg = du_avg + du(jz)
    dv_avg = dv_avg + dv(jz)
    dw_avg = dw_avg + dw(jz)
enddo
du_avg = du_avg / (nz-1)
dv_avg = dv_avg / (nz-1)
dw_avg = dw_avg / (nz-1)

#ifdef PPMPI
call mpi_reduce (uu, uu_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (vv, vv_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (ww, ww_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (du_avg, du_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (dv_avg, dv_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (dw_avg, dw_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then ! note that it's rank here, not coord
    uu = uu_global/nproc
    vv = vv_global/nproc
    ww = ww_global/nproc
    du_avg = du_global/nproc
    dv_avg = dv_global/nproc
    dw_avg = dw_global/nproc
#endif
    open(2,file=path // 'output/check_gmt.dat', status='unknown',      &
        form='formatted', position='append')

    ! one time header output
    if (jt_total==wbase) write(2,*) 'jt_total, uu, vv, ww, du_avg, dv_avg, dw_avg'

    ! continual time-related output
    write(2,*) jt_total, uu, vv, ww, du_avg, dv_avg, dw_avg
    close(2)
#ifdef PPMPI
endif
#endif

end subroutine check_gmt

!*****************************************************************************
subroutine check_gmt_2d3c(u,du)
!*****************************************************************************
! 
! Measure and write the average difference from the target macroscopic value.
! Remember, these differences were measured in the main gm_2d3c routine
! and were found by averaging in the x & y directions. Therefore, only 
! averaging in the z direction here to get a global estimate.
! 
! Added computation of the energy of the u component
! 
use types, only : rprec
use param
use messages
implicit none
integer :: jx, jy, jz
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: u
real(rprec), dimension(nz), intent(in) :: du
real(rprec) :: uu, du_avg
#ifdef PPMPI
real(rprec) :: uu_global, du_global
#endif

! Initialize variables
uu = 0._rprec
du_avg = 0._rprec

! Compute spatial average of squared velocities
do jz = 1, nz-1
do jy = 1, ny
do jx = 1, nxp
    uu = uu + u(jx,jy,jz)**2
enddo
enddo
enddo
uu = uu / (nxp*ny*(nz-1))

! Perform spatial averaging for differences
! Remember these are already averaged in the x & y directions
do jz = 1, nz-1
    du_avg = du_avg + du(jz)
enddo
du_avg = du_avg / (nz-1)

#ifdef PPMPI
call mpi_reduce (uu, uu_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (du_avg, du_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then ! note that it's rank here, not coord
    uu = uu_global/nproc
    du_avg = du_global/nproc
#endif
    open(2,file=path // 'output/check_gmt.dat', status='unknown',      &
        form='formatted', position='append')

    ! one time header output
    if (jt_total==wbase) write(2,*) 'jt_total, uu, du_avg'

    ! continual time-related output
    write(2,*) jt_total, uu, du_avg
    close(2)
#ifdef PPMPI
endif
#endif

end subroutine check_gmt_2d3c

!*****************************************************************************
subroutine gmt_ddx(f, dfdx, lbz)
!*****************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x using spectral decomposition.
!
use types, only : rprec
use param, only : nxp, ny, nz
use fft
use emul_complex, only : OPERATOR(.GMTMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx
integer :: jz

! Loop through horizontal slices
do jz = lbz, nz
    !  Use dfdx to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = f(:,:,jz) / (nxp*ny)
    call dfftw_execute_dft_r2c(gmt_forw, dfdx(:,:,jz), dfdx(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdx(nxp+1:nxp+2,:,jz) = 0._rprec
    dfdx(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdx to perform complex multiplication
    ! Optimized version for real(eye*kx)=0
    ! only passing imaginary part of eye*kx
    dfdx(:,:,jz) = dfdx(:,:,jz) .GMTMULI. gmt_kx

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(gmt_back, dfdx(:,:,jz), dfdx(:,:,jz))
enddo

end subroutine gmt_ddx

!*****************************************************************************
subroutine gmt_ddy(f, dfdy, lbz)
!*****************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! y using spectral decomposition.
!
use types, only : rprec
use param, only : nxp, ny, nz
use fft
use emul_complex, only : OPERATOR(.GMTMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
integer :: jz

! Loop through horizontal slices
do jz = lbz, nz
    !  Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdy(:,:,jz) = f(:,:,jz) / (nxp*ny)
    call dfftw_execute_dft_r2c(gmt_forw, dfdy(:,:,jz), dfdy(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdy(nxp+1:nxp+2,:,jz) = 0._rprec
    dfdy(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdy(:,:,jz) .GMTMULI. gmt_ky

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(gmt_back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine gmt_ddy

!*****************************************************************************
subroutine gmt_ddxy (f, dfdx, dfdy, lbz)
!*****************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x and y using spectral decomposition.
!
use types, only : rprec
use param, only : nxp, ny, nz
use fft
use emul_complex, only : OPERATOR(.GMTMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
integer :: jz

! Loop through horizontal slices
do jz = lbz, nz
    ! Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = f(:,:,jz) / (nxp * ny)
    call dfftw_execute_dft_r2c(gmt_forw, dfdx(:,:,jz), dfdx(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdx(nxp+1:nxp+2,:,jz) = 0._rprec
    dfdx(:,ny/2+1,jz) = 0._rprec

    ! Derivatives: must to y's first here, because we're using dfdx as storage
    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdx(:,:,jz) .GMTMULI. gmt_ky
    dfdx(:,:,jz) = dfdx(:,:,jz) .GMTMULI. gmt_kx

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(gmt_back, dfdx(:,:,jz), dfdx(:,:,jz))
    call dfftw_execute_dft_c2r(gmt_back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine gmt_ddxy

!*****************************************************************************
subroutine gmt_filt_da(f, dfdx, dfdy, lbz)
!*****************************************************************************
!
! This subroutine kills the oddball components in f and computes the partial
! derivative of f with respect to x and y using spectral decomposition.
!
use types, only : rprec
use param, only : nxp, ny, nz
use fft
use emul_complex, only : OPERATOR(.GMTMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(inout) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
integer :: jz

! loop through horizontal slices
do jz = lbz, nz
    ! Calculate FFT in place
    f(:,:,jz) = f(:,:,jz) / (nxp*ny)
    call dfftw_execute_dft_r2c(gmt_forw, f(:,:,jz), f(:,:,jz))

    ! Kill oddballs in zero padded region and Nyquist frequency
    f(nxp+1:nxp+2,:,jz) = 0._rprec
    f(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdx(:,:,jz) = f(:,:,jz) .GMTMULI. gmt_kx
    dfdy(:,:,jz) = f(:,:,jz) .GMTMULI. gmt_ky

    ! Perform inverse transform to get pseudospectral derivative
    ! The oddballs for derivatives should already be dead, since they are for f
    ! inverse transform
    call dfftw_execute_dft_c2r(gmt_back, f(:,:,jz), f(:,:,jz))
    call dfftw_execute_dft_c2r(gmt_back, dfdx(:,:,jz), dfdx(:,:,jz))
    call dfftw_execute_dft_c2r(gmt_back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine gmt_filt_da

!*****************************************************************************
subroutine gmt_diff_stag_array_uv(u,v,RHSx,RHSy,RHSx_f,RHSy_f,txz_half2,tyz_half2,txz,tyz)
!*****************************************************************************
!
! Calculate the intermediate xy-velocity from implicit CN scheme on exit.
!
use types, only : rprec
use param
use messages
use sgs_param, only : nu
use derivatives, only : ddz_w
use fft
#ifdef PPMAPPING
use sim_param, only : jaco_w, jaco_uv, mesh_stretch
#endif

implicit none

real(rprec), dimension(nxp+2,ny,lbz:nz), intent(inout) :: u, v
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: RHSx, RHSy, RHSx_f, RHSy_f
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: txz_half2, tyz_half2, txz, tyz

real(rprec), dimension(nxp+2,ny,0:nz) :: Rx, usol, Ry, vsol
real(rprec), dimension(nxp,ny,0:nz) :: a, b, c
real(rprec), dimension(nxp+2,ny,lbz:nz) :: dtxzdz_rhs, dtyzdz_rhs
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
dtxzdz_rhs(nxp+1:nxp+2, :, 1:nz-1) = 0._rprec
dtyzdz_rhs(nxp+1:nxp+2, :, 1:nz-1) = 0._rprec
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

! Get bottom row 
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    ! Dirichlet BCs
    do jy = 1, ny
    do jx = 1, nxp
        nu_c = nu
        ! Discretized txz(jx,jy,1) as in wallstress,
        ! Therefore BC treated implicitly
#ifdef PPMAPPING
        b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_uv(1))*        &
            (const2*(1._rprec/jaco_w(2))*nu_c + (nu/mesh_stretch(1)))
        c(jx,jy,1) = -const1*(1._rprec/jaco_uv(1))*const2*(1._rprec/jaco_w(2))*nu_c
        Rx(jx,jy,1) = Rx(jx,jy,1) + const1*(1._rprec/jaco_uv(1))*   &
            (nu/mesh_stretch(1))*gmu_bot
        Ry(jx,jy,1) = Ry(jx,jy,1) + const1*(1._rprec/jaco_uv(1))*   &
            (nu/mesh_stretch(1))*gmv_bot
#else
        b(jx,jy,1) = 1._rprec + const1*(const2*nu_c + const3*nu)
        c(jx,jy,1) = -const1*const2*nu_c
        Rx(jx,jy,1) = Rx(jx,jy,1) + const1*const3*nu*gmu_bot
        Ry(jx,jy,1) = Ry(jx,jy,1) + const1*const3*nu*gmv_bot
#endif
    end do
    end do

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
    ! Dirichlet BCs
    do jy = 1, ny
    do jx = 1, nxp
        nu_a = nu
        ! Discretized txz(jx,jy,nz) as in wallstress,
        ! Therefore BC treated implicitly
#ifdef PPMAPPING
        a(jx,jy,nz-1) = -const1*(1._rprec/jaco_uv(nz-1))*const2*(1._rprec/jaco_w(nz-1))*nu_a
        b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/jaco_uv(nz-1))*        &
            (const2*(1._rprec/jaco_w(nz-1))*nu_a + (nu/(L_z-mesh_stretch(nz-1))))
        Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) + const1*(1._rprec/jaco_uv(nz-1))* &
            (nu/(L_z-mesh_stretch(nz-1)))*gmu_top
        Ry(jx,jy,nz-1) = Ry(jx,jy,nz-1) + const1*(1._rprec/jaco_uv(nz-1))* &
            (nu/(L_z-mesh_stretch(nz-1)))*gmv_top
#else
        a(jx,jy,nz-1) = -const1*const2*nu_a
        b(jx,jy,nz-1) = 1._rprec + const1*(const2*nu_a + const3*nu)
        Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) + const1*const3*nu*gmu_top
        Ry(jx,jy,nz-1) = Ry(jx,jy,nz-1) + const1*const3*nu*gmv_top
#endif
    end do
    end do

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
do jx = 1, nxp
    nu_a = nu
#ifdef PPMAPPING
    nu_b = (nu/jaco_w(jz+1)) + (nu/jaco_w(jz))
#else
    nu_b = 2._rprec*nu
#endif
    nu_c = nu

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

! Find intermediate velocity in TDMA
call tridag_array_diff_uv (a, b, c, Rx, usol, nxp)
call tridag_array_diff_uv (a, b, c, Ry, vsol, nxp)

! Fill velocity solution
u(:nxp,:ny,1:nz-1) = usol(:nxp,:ny,1:nz-1)
v(:nxp,:ny,1:nz-1) = vsol(:nxp,:ny,1:nz-1)

end subroutine gmt_diff_stag_array_uv

!*****************************************************************************
subroutine gmt_2d3c_diff_stag_array_uv(u,RHSx,RHSx_f,txz_half2,txz)
!*****************************************************************************
!
! Calculate the intermediate xy-velocity from implicit CN scheme on exit.
!
use types, only : rprec
use param
use messages
use sgs_param, only : nu
use derivatives, only : ddz_w
use fft
#ifdef PPMAPPING
use sim_param, only : jaco_w, jaco_uv, mesh_stretch
#endif

implicit none

real(rprec), dimension(nxp+2,ny,lbz:nz), intent(inout) :: u
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: RHSx, RHSx_f
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: txz_half2, txz

real(rprec), dimension(nxp+2,ny,0:nz) :: Rx, usol
real(rprec), dimension(nxp,ny,0:nz) :: a, b, c
real(rprec), dimension(nxp+2,ny,lbz:nz) :: dtxzdz_rhs
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

! Add explicit portion of Crank-Nicolson
call ddz_w(txz_half2, dtxzdz_rhs, lbz)
dtxzdz_rhs(nxp+1:nxp+2, :, 1:nz-1) = 0._rprec
#ifdef PPSAFETYMODE
#ifdef PPMPI
dtxzdz_rhs(:,:,0) = BOGUS
#endif
dtxzdz_rhs(:,:,nz) = BOGUS
#endif
! Assuming fixed dt here!
Rx(:,:,1:nz-1) = Rx(:,:,1:nz-1) - dt * 0.5_rprec * dtxzdz_rhs(:,:,1:nz-1)

! Get bottom row 
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    ! Dirichlet BCs
    do jy = 1, ny
    do jx = 1, nxp
        nu_c = nu
        ! Discretized txz(jx,jy,1) as in wallstress,
        ! Therefore BC treated implicitly
#ifdef PPMAPPING
        b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_uv(1))*        &
            (const2*(1._rprec/jaco_w(2))*nu_c + (nu/mesh_stretch(1)))
        c(jx,jy,1) = -const1*(1._rprec/jaco_uv(1))*const2*(1._rprec/jaco_w(2))*nu_c
        Rx(jx,jy,1) = Rx(jx,jy,1) + const1*(1._rprec/jaco_uv(1))*   &
            (nu/mesh_stretch(1))*gmu_bot
#else
        b(jx,jy,1) = 1._rprec + const1*(const2*nu_c + const3*nu)
        c(jx,jy,1) = -const1*const2*nu_c
        Rx(jx,jy,1) = Rx(jx,jy,1) + const1*const3*nu*gmu_bot
#endif
    end do
    end do

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
    ! Dirichlet BCs
    do jy = 1, ny
    do jx = 1, nxp
        nu_a = nu
        ! Discretized txz(jx,jy,nz) as in wallstress,
        ! Therefore BC treated implicitly
#ifdef PPMAPPING
        a(jx,jy,nz-1) = -const1*(1._rprec/jaco_uv(nz-1))*const2*(1._rprec/jaco_w(nz-1))*nu_a
        b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/jaco_uv(nz-1))*        &
            (const2*(1._rprec/jaco_w(nz-1))*nu_a + (nu/(L_z-mesh_stretch(nz-1))))
        Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) + const1*(1._rprec/jaco_uv(nz-1))* &
            (nu/(L_z-mesh_stretch(nz-1)))*gmu_top
#else
        a(jx,jy,nz-1) = -const1*const2*nu_a
        b(jx,jy,nz-1) = 1._rprec + const1*(const2*nu_a + const3*nu)
        Rx(jx,jy,nz-1) = Rx(jx,jy,nz-1) + const1*const3*nu*gmu_top
#endif
    end do
    end do

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
do jx = 1, nxp
    nu_a = nu
#ifdef PPMAPPING
    nu_b = (nu/jaco_w(jz+1)) + (nu/jaco_w(jz))
#else
    nu_b = 2._rprec*nu
#endif
    nu_c = nu

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

! Find intermediate velocity in TDMA
call tridag_array_diff_uv (a, b, c, Rx, usol, nxp)

! Fill velocity solution
u(:nxp,:ny,1:nz-1) = usol(:nxp,:ny,1:nz-1)

end subroutine gmt_2d3c_diff_stag_array_uv

!*****************************************************************************
subroutine gmt_diff_stag_array_w(w,RHSz,RHSz_f,tzz)
!*****************************************************************************
!
! Calculate the intermediate z-velocity from implicit CN scheme on exit.
!
use types, only : rprec
use param
use messages
use sgs_param, only : nu
use derivatives, only : ddz_uv
use fft
#ifdef PPMAPPING
use sim_param, only : jaco_w, jaco_uv, mesh_stretch
#endif

implicit none

real(rprec), dimension(nxp+2,ny,lbz:nz), intent(inout) :: w
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: RHSz, RHSz_f, tzz

real(rprec), dimension(nxp+2,ny,0:nz) :: Rz, wsol
real(rprec), dimension(nxp,ny,0:nz) :: a, b, c
real(rprec), dimension(nxp+2,ny,lbz:nz) :: dtzzdz_rhs
real(rprec) :: nu_a, nu_b, nu_c, nu_r, const1, const2
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
dtzzdz_rhs(nxp+1:nxp+2, :, 1:nz-1) = 0._rprec
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

! Imposing no-penetration BC regardless of whether or not
! it is stress-free or a wall. However still need to select
! case lbc_mom/ubc_mom because the eddy viscosity is in 
! different locations otherwise.

! Get bottom row, at jz = 2 instead of jz = 1
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,2) = BOGUS
#endif
    ! Dirichlet BC, wbot = 0
    do jy = 1, ny
    do jx = 1, nxp
#ifdef PPMAPPING
        nu_b = (nu/jaco_uv(2)) + (nu/jaco_uv(1))
#else
        nu_b = 2.0_rprec*nu
#endif
        nu_c = nu
#ifdef PPMAPPING
        b(jx,jy,2) = 1._rprec + const1*(1._rprec/jaco_w(2))*const2*nu_b
        c(jx,jy,2) = -const1*(1._rprec/jaco_w(2))*const2*(1._rprec/jaco_uv(2))*nu_c
#else
        b(jx,jy,2) = 1._rprec + const1*const2*nu_b
        c(jx,jy,2) = -const1*const2*nu_c
#endif
    end do
    end do
    jz_min = 3
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
    ! Dirichlet BC, wtop = 0
    do jy = 1, ny
    do jx = 1, nxp
        nu_a = nu
#ifdef PPMAPPING
        nu_b = (nu/jaco_uv(nz-2)) + (nu/jaco_uv(nz-1))
#else
        nu_b = 2.0_rprec*nu
#endif
#ifdef PPMAPPING
        a(jx,jy,nz-1) = -const1*(1._rprec/jaco_w(nz-1))*const2*(1._rprec/jaco_uv(nz-2))*nu_a
        b(jx,jy,nz-1) = 1._rprec + const1*(1._rprec/jaco_w(nz-1))*const2*nu_b
#else
        a(jx,jy,nz-1) = -const1*const2*nu_a
        b(jx,jy,nz-1) = 1._rprec + const1*const2*nu_b
#endif
    end do
    end do
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
do jx = 1, nxp
    nu_a = nu
#ifdef PPMAPPING
    nu_b = (nu/jaco_uv(jz)) + (nu/jaco_uv(jz-1))
#else
    nu_b = 2._rprec*nu
#endif
    nu_c = nu

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

! Find intermediate velocity in TDMA
call tridag_array_diff_w (a, b, c, Rz, wsol, nxp)

! Fill velocity solution
if (coord == 0) then
    w(:,:,1) = w(:,:,1) + dt * ( tadv1 * RHSz(:,:,1) + tadv2 * RHSz_f(:,:,1) )
    w(:nxp,:ny,2:nz-1) = wsol(:nxp,:ny,2:nz-1)
else
    w(:nxp,:ny,1:nz-1) = wsol(:nxp,:ny,1:nz-1)
endif

if (coord == nproc-1) then
    w(:,:,nz) = w(:,:,nz) + dt * ( tadv1 * RHSz(:,:,nz) + tadv2 * RHSz_f(:,:,nz) )
endif

end subroutine gmt_diff_stag_array_w

!*****************************************************************************
subroutine gmt_press_stag_array(u,v,w,divtz,p,dpdx,dpdy,dpdz)
!*****************************************************************************
!
! Calculate the pressure and its derivatives on exit. Everything is in 
! physical space on exit.
!
use types, only : rprec
use param
use messages
use fft
#ifdef PPMAPPING
use sim_param, only : jaco_w, jaco_uv
#endif

implicit none

real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz
integer :: ir, ii
integer :: jz_min
integer :: end_kx

real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: u, v, w, divtz
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(out) :: p
real(rprec), dimension(nxp+2,ny,nz), intent(out) :: dpdx, dpdy, dpdz

real(rprec), save, dimension(:,:,:), allocatable :: rH_x, rH_y, rH_z
real(rprec), save, dimension(:,:), allocatable :: rtopw, rbottomw
real(rprec), save, dimension(:,:,:), allocatable :: RHS_col
real(rprec), save, dimension(:,:,:), allocatable :: a, b, c

logical, save :: arrays_allocated = .false.

real(rprec), dimension(2) :: aH_x, aH_y

! Specifiy cached constants
const = 1._rprec/(nxp*ny)
const2 = const/tadv1/dt
const3 = 1._rprec/(dz**2)
const4 = 1._rprec/(dz)

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( rH_x(nxp+2,ny,lbz:nz), rH_y(nxp+2,ny,lbz:nz), rH_z(nxp+2,ny,lbz:nz) )
    allocate ( rtopw(nxp+2,ny), rbottomw(nxp+2,ny) )
    allocate ( RHS_col(nxp+2,ny,nz+1) )
    allocate ( a(nxp/2+1,ny,nz+1), b(nxp/2+1,ny,nz+1), c(nxp/2+1,ny,nz+1) )

    arrays_allocated = .true.
endif

if (coord == 0) then
    p(:,:,0) = 0._rprec
#ifdef PPSAFETYMODE
else
    p(:,:,0) = BOGUS
#endif
end if

! Get the right hand side ready
! Loop over levels
! Recall that the old timestep guys already contain the pressure
do jz = 1, nz-1
    rH_x(:,:,jz) = const2 * u(:,:,jz)
    rH_y(:,:,jz) = const2 * v(:,:,jz)
    rH_z(:,:,jz) = const2 * w(:,:,jz)

    call dfftw_execute_dft_r2c(gmt_forw, rH_x(:,:,jz), rH_x(:,:,jz))
    call dfftw_execute_dft_r2c(gmt_forw, rH_y(:,:,jz), rH_y(:,:,jz))
    call dfftw_execute_dft_r2c(gmt_forw, rH_z(:,:,jz), rH_z(:,:,jz))
end do

#if defined(PPMPI) && defined(PPSAFETYMODE)
  !Careful - only update real values (odd indicies)
  rH_x(1:(nxp+2):2,:,0) = BOGUS
  rH_y(1:(nxp+2):2,:,0) = BOGUS
  rH_z(1:(nxp+2):2,:,0) = BOGUS
#endif

#ifdef PPSAFETYMODE
!Careful - only update real values (odd indicies)
rH_x(1:(nxp+2):2,:,nz) = BOGUS
rH_y(1:(nxp+2):2,:,nz) = BOGUS
#endif

#ifdef PPMPI
if (coord == nproc-1) then
    rH_z(:,:,nz) = const2 * w(:,:,nz)
    call dfftw_execute_dft_r2c(gmt_forw, rH_z(:,:,nz), rH_z(:,:,nz))
#ifdef PPSAFETYMODE
else
    rH_z(1:(nxp+2):2,:,nz) = BOGUS !--perhaps this should be 0 on top process?
#endif
endif
#else
rH_z(:,:,nz) = const2 * w(:,:,nz)
call dfftw_execute_dft_r2c(gmt_forw, rH_z(:,:,nz), rH_z(:,:,nz))
#endif

if (coord == 0) then
    rbottomw(:,:) = const * divtz(:,:,1)
    call dfftw_execute_dft_r2c(gmt_forw, rbottomw, rbottomw )
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    rtopw(:,:) = const * divtz(:,:,nz)
    call dfftw_execute_dft_r2c(gmt_forw, rtopw, rtopw)
#ifdef PPMPI
endif
#endif

! set oddballs to 0
rH_x(nxp+1:nxp+2,:,1:nz-1) = 0._rprec
rH_y(nxp+1:nxp+2,:,1:nz-1) = 0._rprec
rH_z(nxp+1:nxp+2,:,1:nz-1) = 0._rprec
rH_x(:,ny/2+1,1:nz-1) = 0._rprec
rH_y(:,ny/2+1,1:nz-1) = 0._rprec
rH_z(:,ny/2+1,1:nz-1) = 0._rprec
! should also set to zero for rH_z (nz) on coord == nproc-1
if (coord == nproc-1) then
    rH_z(nxp+1:nxp+2,:,nz) = 0._rprec
    rH_z(:,ny/2+1,nz) = 0._rprec
end if

! with MPI; topw and bottomw are only on top & bottom processes
rtopw(nxp+1:nxp+2, :) = 0._rprec
rtopw(:, ny/2+1) = 0._rprec
rbottomw(nxp+1:nxp+2, :) = 0._rprec
rbottomw(:, ny/2+1) = 0._rprec

! Loop over (Kx,Ky) to solve for Pressure amplitudes
if (coord == 0) then
    !  a, b, and c are treated as the real part of a complex array
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    b(:,:,1) = -1._rprec
    c(:,:,1) = 1._rprec
#ifdef PPMAPPING
    RHS_col(:,:,1) = - jaco_w(1)*dz* rbottomw(:,:)
#else
    RHS_col(:,:,1) = -dz * rbottomw(:,:)
#endif
    jz_min = 2
else
  jz_min = 1
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    !--top nodes
    a(:,:,nz+1) = -1._rprec
    b(:,:,nz+1) = 1._rprec
#ifdef PPSAFETYMODE
    c(:,:,nz+1) = BOGUS
#endif
#ifdef PPMAPPING
    RHS_col(:,:,nz+1) = -jaco_w(nz)*dz * rtopw(:,:)
#else
    RHS_col(:,:,nz+1) = -dz * rtopw(:,:)
#endif
#ifdef PPMPI
endif
#endif

#ifdef PPMPI
    call mpi_sendrecv (rH_x(1, 1, nz-1), (nxp+2)*ny, MPI_RPREC, up, 1,      &
        rH_x(1, 1, 0), (nxp+2)*ny, MPI_RPREC, down, 1, comm, status, ierr)
    call mpi_sendrecv (rH_y(1, 1, nz-1), (nxp+2)*ny, MPI_RPREC, up, 2,      &
        rH_y(1, 1, 0), (nxp+2)*ny, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, nz-1), (nxp+2)*ny, MPI_RPREC, up, 3,      &
        rH_z(1, 1, 0), (nxp+2)*ny, MPI_RPREC, down, 3, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, 1), (nxp+2)*ny, MPI_RPREC, down, 6,       &
        rH_z(1, 1, nz), (nxp+2)*ny, MPI_RPREC, up, 6, comm, status, ierr)
#endif

! end_kx = lh-1
end_kx = (nxp/2 + 1) - 1

do jz = jz_min, nz
do jy = 1, ny
    if (jy == ny/2 + 1) cycle

    do jx = 1, end_kx

        if (jx*jy == 1) cycle

        ii = 2*jx   ! imaginary index
        ir = ii - 1 ! real index

        ! JDA dissertation, eqn(2.85) a,b,c=coefficients and RHS_col=r_m
#ifdef PPMAPPING
        a(jx, jy, jz) = const3*(1._rprec/(jaco_uv(jz-1)))*(1._rprec/(jaco_w(jz-1)))
        b(jx, jy, jz) = -(gmt_kx(jx,jy)**2 + gmt_ky(jx,jy)**2               &
            + const3*(1._rprec/(jaco_uv(jz-1)))*                            &
            (1._rprec/(jaco_w(jz-1))+1._rprec/(jaco_w(jz))))
        c(jx, jy, jz) = const3*(1._rprec/(jaco_uv(jz-1)))*(1._rprec/(jaco_w(jz)))
#else
        a(jx, jy, jz) = const3
        b(jx, jy, jz) = -(gmt_kx(jx,jy)**2 + gmt_ky(jx,jy)**2 + 2._rprec*const3)
        c(jx, jy, jz) = const3
#endif

        !  Compute eye * kx * H_x
        aH_x(1) = -rH_x(ii,jy,jz-1) * gmt_kx(jx,jy)
        aH_x(2) =  rH_x(ir,jy,jz-1) * gmt_kx(jx,jy)
        aH_y(1) = -rH_y(ii,jy,jz-1) * gmt_ky(jx,jy)
        aH_y(2) =  rH_y(ir,jy,jz-1) * gmt_ky(jx,jy)

#ifdef PPMAPPING
        RHS_col(ir:ii,jy,jz) = aH_x + aH_y + (rH_z(ir:ii, jy, jz) -         &
            rH_z(ir:ii, jy, jz-1))*const4/jaco_uv(jz-1)
#else
        RHS_col(ir:ii,jy,jz) =  aH_x + aH_y + (rH_z(ir:ii, jy, jz) -        &
            rH_z(ir:ii, jy, jz-1))*const4
#endif

    end do
end do
end do

! this skips zero wavenumber solution, nyquist freqs
call tridag_array (a, b, c, RHS_col, p, nxp)

! zero-wavenumber solution
#ifdef PPMPI
! wait for p(1, 1, 1) from "down"
call mpi_recv (p(1:2, 1, 1), 2, MPI_RPREC, down, 8, comm, status, ierr)
#endif

if (coord == 0) then
    p(1:2, 1, 0) = 0._rprec !! BC, arbitrary pressure
#ifdef PPMAPPING
    p(1:2, 1, 1) = p(1:2,1,0) - jaco_w(1)*dz * rbottomw(1:2,1)
#else
    p(1:2, 1, 1) = p(1:2,1,0) - dz * rbottomw(1:2,1)
#endif
end if

do jz = 2, nz
    ! JDA dissertation, eqn(2.88)
#ifdef PPMAPPING
    p(1:2, 1, jz) = p(1:2, 1, jz-1) + rH_z(1:2, 1, jz) * dz * jaco_w(jz)
#else
    p(1:2, 1, jz) = p(1:2, 1, jz-1) + rH_z(1:2, 1, jz) * dz
#endif
end do

#ifdef PPMPI
! send p(1, 1, nz) to "up"
call mpi_send (p(1:2, 1, nz), 2, MPI_RPREC, up, 8, comm, ierr)
#endif

#ifdef PPMPI
! make sure 0 <-> nz-1 are syncronized
! 1 <-> nz should be in sync already
call mpi_sendrecv (p(1, 1, nz-1), (nxp+2)*ny, MPI_RPREC, up, 2,             &
    p(1, 1, 0), (nxp+2)*ny, MPI_RPREC, down, 2, comm, status, ierr)
#endif

! zero the nyquist freqs
p(nxp+1:nxp+2,:,:) = 0._rprec
p(:,ny/2+1,:) = 0._rprec

! Now need to get p(wave,level) to physical p(jx,jy,jz)
! Loop over height levels
call dfftw_execute_dft_c2r(gmt_back,p(:,:,0),p(:,:,0))

do jz = 1, nz-1
    do jy = 1, ny
    do jx = 1, end_kx
        ii = 2*jx   ! imaginary index
        ir = ii - 1 ! real index
        dpdx(ir,jy,jz) = -p(ii,jy,jz) * gmt_kx(jx,jy)
        dpdx(ii,jy,jz) =  p(ir,jy,jz) * gmt_kx(jx,jy)
        dpdy(ir,jy,jz) = -p(ii,jy,jz) * gmt_ky(jx,jy)
        dpdy(ii,jy,jz) =  p(ir,jy,jz) * gmt_ky(jx,jy)
    end do
    end do

    ! note the oddballs of p are already 0, so we should be OK here
    call dfftw_execute_dft_c2r(gmt_back,dpdx(:,:,jz), dpdx(:,:,jz))
    call dfftw_execute_dft_c2r(gmt_back,dpdy(:,:,jz), dpdy(:,:,jz))
    call dfftw_execute_dft_c2r(gmt_back,p(:,:,jz), p(:,:,jz))
end do

if (coord==nproc-1) then
    call dfftw_execute_dft_c2r(gmt_back,p(:,:,nz),p(:,:,nz))
endif

! nz level is not needed elsewhere (although its valid)
#ifdef PPSAFETYMODE
dpdx(:,:,nz) = BOGUS
dpdy(:,:,nz) = BOGUS
if(coord<nproc-1) p(:,:,nz) = BOGUS
#endif

! Final step compute the z-derivative of p
! note: p has additional level at z=-dz/2 for this derivative
#ifdef PPMAPPING
do jz = 1, nz-1
    dpdz(1:nxp, 1:ny, jz) = (p(1:nxp, 1:ny, jz) - p(1:nxp, 1:ny, jz-1))     &
        / dz / jaco_w(jz)
end do
#else
dpdz(1:nxp, 1:ny, 1:nz-1) = (p(1:nxp, 1:ny, 1:nz-1) - p(1:nxp, 1:ny, 0:nz-2)) / dz
#endif
#ifdef PPSAFETYMODE
if(coord<nproc-1)  dpdz(:,:,nz) = BOGUS
#endif
#ifdef PPMAPPING
if(coord==nproc-1) dpdz(1:nxp,1:ny,nz) = (p(1:nxp,1:ny,nz)-p(1:nxp,1:ny,nz-1))  &
    / dz / jaco_w(nz)
#else
if(coord==nproc-1) dpdz(1:nxp,1:ny,nz) = (p(1:nxp,1:ny,nz)-p(1:nxp,1:ny,nz-1)) / dz
#endif

end subroutine gmt_press_stag_array

!*****************************************************************************
subroutine gmt_padd(u_big,u)
!*****************************************************************************
! puts arrays into larger, zero-padded arrays
! automatically zeroes the oddballs
use types, only : rprec
use param, only : nxp,ny,ny2
implicit none

!  u and u_big are interleaved as complex arrays
real(rprec), dimension(nxp+2,ny), intent(in) :: u
real(rprec), dimension(3*nxp/2 + 2,ny2), intent(out) :: u_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

! make sure the big array is zeroed!
u_big(:,:) = 0._rprec

! note: split access in an attempt to maintain locality
u_big(:nxp,:ny_h) = u(:nxp,:ny_h)

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

u_big(:nxp,j_big_s:ny2) = u(:nxp,j_s:ny)

end subroutine gmt_padd

!*****************************************************************************
subroutine gmt_unpadd(cc,cc_big)
!*****************************************************************************
use types, only : rprec
use param, only : nxp,ny,ny2
implicit none

!  cc and cc_big are interleaved as complex arrays
real(rprec), dimension( nxp+2, ny ) :: cc
real(rprec), dimension( 3*nxp/2 + 2, ny2 ) :: cc_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

cc(:nxp,:ny_h) = cc_big(:nxp,:ny_h)

! oddballs
cc(nxp+1:nxp+2,:) = 0._rprec
cc(:,ny_h+1) = 0._rprec

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2
cc(:nxp,j_s:ny) = cc_big(:nxp,j_big_s:ny2)

end subroutine gmt_unpadd

!*****************************************************************************
subroutine x_avg_mfm( vel )
!*****************************************************************************
!
! This function provides the streamwise averaged value of the velocity field.
! 
! This function is copied from x_avg in functions.f90 with modifications to
! input/outputs.
! 
! The input should be in physical space with size nxp+2, also the input is
! over-written and outputted.
!

use types, only: rprec
use param, only: nxp, ny, lbz, nz

implicit none
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(inout)  :: vel
real(rprec), dimension(ny,lbz:nz) :: vel_sum
integer :: i

vel_sum = 0._rprec

do i = 1, nxp
    vel_sum = vel_sum + vel(i,:,:)
enddo

vel_sum = vel_sum / real(nxp,rprec)

! Overwrite input, then output streamwise mean field
do i = 1, nxp
    vel(i,:,:) = vel_sum
enddo


end subroutine x_avg_mfm

end module mfm
