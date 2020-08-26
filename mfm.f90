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

public :: mfm_init, ic_gmt, gm_transport, mfm_checkpoint

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

! Boundary conditions
real(rprec), public :: gmu_bot = 0._rprec
real(rprec), public :: gmu_top = 0._rprec
real(rprec), public :: gmv_bot = 0._rprec
real(rprec), public :: gmv_top = 0._rprec

! Initial condition
integer, public :: ic_mfm = 1
real(rprec), public :: initial_noise_gmt = 0._rprec

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
! If .not. fourier, then nx = nxp and uvw are the same size as gmu
! If fourier, then nx < nxp and uvw and gmu are different sizes
allocate ( u_big(ld_big, ny2, lbz:nz) ); u_big = 0._rprec
allocate ( v_big(ld_big, ny2, lbz:nz) ); v_big = 0._rprec
allocate ( w_big(ld_big, ny2, lbz:nz) ); w_big = 0._rprec
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

    select case (lbc_mom)
        ! Stress free
        case (0)
            dgmudz(:,:,1) = 0.0_rprec
            dgmvdz(:,:,1) = 0.0_rprec
            gmtxz(:,:,1) = 0.0_rprec
            gmtyz(:,:,1) = 0.0_rprec

        ! DNS wall
        case (1)
            do j = 1, ny
            do i = 1, nxp
                dgmudz(i,j,1) = ( gmu(i,j,1) - gmu_bot ) / denom
                dgmvdz(i,j,1) = ( gmv(i,j,1) - gmv_bot ) / denom
                gmtxz(i,j,1) = -nu_molec/(z_i*u_star)*dgmudz(i,j,1)
                gmtyz(i,j,1) = -nu_molec/(z_i*u_star)*dgmvdz(i,j,1)
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
            dgmudz(:,:,nz) = 0.0_rprec
            dgmvdz(:,:,nz) = 0.0_rprec
            gmtxz(:,:,nz) = 0.0_rprec
            gmtyz(:,:,nz) = 0.0_rprec

        ! DNS wall
        case (1)
            do j = 1, ny
            do i = 1, nxp
                dgmudz(i,j,nz) = ( gmu_top - gmu(i,j,nz-1) ) / denom
                dgmvdz(i,j,nz) = ( gmv_top - gmv(i,j,nz-1) ) / denom
                gmtxz(i,j,nz) = -nu_molec/(z_i*u_star)*dgmudz(i,j,nz)
                gmtyz(i,j,nz) = -nu_molec/(z_i*u_star)*dgmvdz(i,j,nz)
            enddo
            enddo

    end select
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
    call dfftw_execute_dft_r2c(forw, temp, temp)
    call padd(a_big(:,:,jz), temp)
    call dfftw_execute_dft_c2r(back_big, a_big(:,:,jz), a_big(:,:,jz))
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
    call dfftw_execute_dft_r2c(forw_big, a_big(:,:,jz), a_big(:,:,jz))
    call unpadd(a(:,:,jz), a_big(:,:,jz))
    call dfftw_execute_dft_c2r(back, a(:,:,jz), a(:,:,jz))
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
use derivatives, only : filt_da, ddx, ddy, ddz_uv, ddz_w
#ifdef PPMPI
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
use forcing, only : project

real(rprec) :: diff_coef, rmsdivvel, const
integer :: jx, jy, jz, jz_min, jz_max, nxp2

real(rprec), dimension(nz) :: gmu_avg, gmv_avg, gmw_avg, du, dv, dw

! Save previous timestep's RHS
rhs_gmx_f = rhs_gmx
rhs_gmy_f = rhs_gmy
rhs_gmz_f = rhs_gmz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate derivatives of generalized momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dudx, dudy, dvdx, dvdy, dwdx, dwdy derivatives (in Fourier space)
call filt_da(gmu, dgmudx, dgmudy, lbz)
call filt_da(gmv, dgmvdx, dgmvdy, lbz)
call filt_da(gmw, dgmwdx, dgmwdy, lbz)

! dudz, dvdz using finite differences (for 1:nz on uv-nodes)
! except bottom coord, only 2:nz
call ddz_uv(gmu, dgmudz, lbz)
call ddz_uv(gmv, dgmvdz, lbz)

! dwdz using finite differences (for 0:nz-1 on w-nodes)
! except bottom coord, only 1:nz-1
call ddz_w(gmw, dgmwdz, lbz)

!debug
write(*,*) coord, dgmudx(1,1,:)

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
! Using the same divstress routines used in the main loop
#ifdef PPCNDIFF
call divstress_uv(div_gmtx, div_gmty, gmtxx, gmtxy, gmtxz_half1, gmtyy, gmtyz_half1)
call divstress_w_cndiff(div_gmtz, gmtxz, gmtyz, gmtzz)
#else
call divstress_uv(div_gmtx, div_gmty, gmtxx, gmtxy, gmtxz, gmtyy, gmtyz)
call divstress_w(div_gmtz, gmtxz, gmtyz, gmtzz)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advective term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This differs from the convec subroutine in the main loop in two main ways:
! 1. Uses coupling between simulation velocities (u,v,w) and GMT velocities
! 2. Is not u x omega, instead computed as u[sim]*grad(u[gmt])

! Move variables to big domain
if (fourier) then
    call to_big_fourier(u, u_big)
    call to_big_fourier(v, v_big)
    call to_big_fourier(w, w_big)
else
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
if (fourier) then
    ! Remember u,v,w(kx,y,z) after to_big_fourier routine
    ! Using streamwise constant mean (kx=0) to advect generalized momentum
    do jy = 1, ny2
        if (coord == 0) then
            ! Bottom wall take w(jz=1) = 0
            temp_big(:,jy,1) = const*(u_big(1,jy,1)*dgmudx_big(:,jy,1) +     &
                v_big(1,jy,1)*dgmudy_big(:,jy,1) +                           &
                0.5_rprec*w_big(1,jy,2)*dgmudz_big(:,jy,2))
            jz_min = 2
        else
            jz_min = 1
        endif

        if (coord == nproc-1) then
            ! Top wall take w(jz=nz) = 0
            temp_big(:,jy,nz-1) = const*(u_big(1,jy,nz-1)*dgmudx_big(:,jy,nz-1) +   &
                v_big(1,jy,nz-1)*dgmudy_big(:,jy,nz-1) +                            &
                0.5_rprec*w_big(1,jy,nz-1)*dgmudz_big(:,jy,nz-1))
            jz_max = nz-2
        else
            jz_max = nz-1
        endif

        ! For entire domain
        do jz = jz_min, jz_max
            temp_big(:,jy,jz) = const*(u_big(1,jy,jz)*dgmudx_big(:,jy,jz) +  &
                v_big(1,jy,jz)*dgmudy_big(:,jy,jz) +                         &
                0.5_rprec*(w_big(1,jy,jz+1)*dgmudz_big(:,jy,jz+1) +          &
                w_big(1,jy,jz)*dgmudz_big(:,jy,jz)))
!debug - shows NaN for all jx positions except for the last two
!if ((coord==0) .and. (jy == 1)) write(*,*) jz, temp_big(:,jy,jz)
!debug - filled with numbers and zeros in the wrong places
!if ((coord==0) .and. (jy == 1)) write(*,*) jz, dgmudx_big(:,jy,jz)
        enddo
    enddo
else
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
endif

! Move temp_big into RHSx for GMT and make small
call to_small(temp_big, rhs_gmx)

!debug -- shows zeroes on some z-levels, but others are filled
!if ((coord == 0)) write(*,*) rhs_gmx(1,1,:)

! Compute advective term in y-GMT
! Interpolate w and dvdz onto uv-grid
if (fourier) then
    ! Remember u,v,w(kx,y,z) after to_big_fourier routine
    ! Using streamwise constant mean (kx=0) to advect generalized momentum
    do jy = 1, ny2
        if (coord == 0) then
            ! Bottom wall take w(jz=1) = 0
            temp_big(:,jy,1) = const*(u_big(1,jy,1)*dgmvdx_big(:,jy,1) +    &
                v_big(1,jy,1)*dgmvdy_big(:,jy,1) +                          &
                0.5_rprec*w_big(1,jy,2)*dgmvdz_big(:,jy,2))
            jz_min = 2
        else
            jz_min = 1
        endif

        if (coord == nproc-1) then
            ! Top wall take w(jz=nz) = 0
            temp_big(:,jy,nz-1) = const*(u_big(1,jy,nz-1)*dgmvdx_big(:,jy,nz-1) +   &
                v_big(1,jy,nz-1)*dgmvdy_big(:,jy,nz-1) +                            &
                0.5_rprec*w_big(1,jy,nz-1)*dgmvdz_big(:,jy,nz-1))
            jz_max = nz-2
        else
            jz_max = nz-1
        endif

        ! For entire domain
        do jz = jz_min, jz_max
            temp_big(:,jy,jz) = const*(u_big(1,jy,jz)*dgmvdx_big(:,jy,jz) +      &
                v_big(1,jy,jz)*dgmvdy_big(:,jy,jz) +                             &
                0.5_rprec*(w_big(1,jy,jz+1)*dgmvdz_big(:,jy,jz+1) +              &
                w_big(1,jy,jz)*dgmvdz_big(:,jy,jz)))
        enddo
    enddo
else
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
    do jz = jz_min, jz_max
        temp_big(:,:,jz) = const*(u_big(:,:,jz)*dgmvdx_big(:,:,jz) +      &
            v_big(:,:,jz)*dgmvdy_big(:,:,jz) +                            &
            0.5_rprec*(w_big(:,:,jz+1)*dgmvdz_big(:,:,jz+1) +             &
            w_big(:,:,jz)*dgmvdz_big(:,:,jz)))
    enddo
endif

! Move temp_big into RHSy for GMT and make small
call to_small(temp_big, rhs_gmy)

! Compute advective term in z-GMT
! Interpolate u, v, and dwdz onto w-grid
if (fourier) then
    ! Remember u,v,w(kx,y,z) after to_big_fourier routine
    ! Using streamwise constant mean (kx=0) to advect generalized momentum
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
    do jy = 1, ny2
        do jz = jz_min, jz_max
            temp_big(:,jy,jz) = const*(                                                &
                0.5_rprec*(u_big(1,jy,jz)+u_big(1,jy,jz-1))*dgmwdx_big(:,jy,jz) +      &
                0.5_rprec*(v_big(1,jy,jz)+v_big(1,jy,jz-1))*dgmwdy_big(:,jy,jz) +      &
                w_big(1,jy,jz)*0.5_rprec*(dgmwdz_big(:,jy,jz)+dgmwdz_big(:,jy,jz-1)))
        enddo
    enddo
else
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
endif

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
call diff_stag_array_uv(gmu,gmv,rhs_gmx,rhs_gmy,rhs_gmx_f,rhs_gmy_f,        &
    gmtxz_half2,gmtyz_half2,gmtxz,gmtyz)
call diff_stag_array_w(gmw,rhs_gmz,rhs_gmz_f,gmtzz)
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
call press_stag_array(gmu,gmv,gmw,div_gmtz,gmp,dgmpdx,dgmpdy,dgmpdz)

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
    call rmsdiv(dgmudx,dgmvdy,dgmwdz,rmsdivvel)

    if(coord == 0) then
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
        write(*,'(a)') '====================== GMT TOP ========================='
        write(*,*) 'u: ', gmu(nxp/2,ny/2,nz-2:nz-1)
        write(*,*) 'v: ', gmv(nxp/2,ny/2,nz-2:nz-1)
        write(*,*) 'w: ', gmw(nxp/2,ny/2,nz-1:nz)
        write(*,'(a)') '========================================================'
    end if
    call mpi_barrier(comm, ierr)
end if

end subroutine gm_transport

end module mfm
