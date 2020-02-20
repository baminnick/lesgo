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
subroutine vreman
!*******************************************************************************
!
! This subroutine calculates the eddy viscosity predicted by
! Vreman, PF, 16: 3670-3681 (2004) 
!
use types, only : rprec
use param
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
use sgs_param, only : Nu_t
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
#ifdef PPMAPPING
use sim_param, only : JACO1, JACO2
#endif

implicit none

real(rprec), dimension(ld,ny) :: dudx_p, dudy_p, dudz_p, &
    dvdx_p, dvdy_p, dvdz_p, dwdx_p, dwdy_p, dwdz_p
real(rprec), dimension(ld,ny) :: alpha2, beta
real(rprec), dimension(ld,ny) :: b11, b12, b13, b22, b23, b33
real(rprec) :: dz_p, cvre, eps
integer :: jz, jz_min, jz_max

! Set constants
cvre = 0.07_rprec !! to be tuned
eps = 1e-12 !! some small number in case denominator is zero
dz_p = dz !! to be overwritten in mapping to stretched grid

! Calculate Nu_t for jz = 1 (coord==0 only)
!   stored on uvp-nodes (this level only) for 'wall'
!   stored on w-nodes (all) for 'stress free'
! See sgs_stag_util/calc_Sij for why derivatives are interpolated this way
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free
        case(0)
            ! Interpolate derivatives onto w-node
            dudx_p(:,:) = dudx(:,:,1)
            dudy_p(:,:) = dudy(:,:,1)
            dudz_p(:,:) = dudz(:,:,1)
            dvdx_p(:,:) = dvdx(:,:,1)
            dvdy_p(:,:) = dvdy(:,:,1)
            dvdz_p(:,:) = dvdz(:,:,1)
            dwdx_p(:,:) = dwdx(:,:,1)
            dwdy_p(:,:) = dwdy(:,:,1)
            dwdz_p(:,:) = 0.5_rprec*(dwdz(:,:,1) + 0.0_rprec)

#ifdef PPMAPPING
            dz_p = JACO1(1)*dz
#endif

        ! Wall
        case(1:)
            ! Interpolate derivatives onto uvp-node
            dudx_p(:,:) = dudx(:,:,1)
            dudy_p(:,:) = dudy(:,:,1)
            dudz_p(:,:) = 0.5_rprec*(dudz(:,:,1) + dudz(:,:,2))
            dvdx_p(:,:) = dvdx(:,:,1)
            dvdy_p(:,:) = dvdy(:,:,1)
            dvdz_p(:,:) = 0.5_rprec*(dvdz(:,:,1) + dvdz(:,:,2))
            dwdx_p(:,:) = 0.5_rprec*(dwdx(:,:,1) + dwdx(:,:,2))
            dwdy_p(:,:) = 0.5_rprec*(dwdy(:,:,1) + dwdy(:,:,2))
            dwdz_p(:,:) = dwdz(:,:,1)

#ifdef PPMAPPING
            dz_p = JACO2(1)*dz
#endif

        ! Wall-Model
        case(2:)
            ! Interpolate derivatives onto uvp-node
            dudx_p(:,:) = dudx(:,:,1)
            dudy_p(:,:) = dudy(:,:,1)
            ! dudz from wallstress, therefore no interpolation
            dudz_p(:,:) = dudz(:,:,1)
            dvdx_p(:,:) = dvdx(:,:,1)
            dvdy_p(:,:) = dvdy(:,:,1)
            ! dvdz from wallstress, therefore no interpolation
            dvdz_p(:,:) = dvdz(:,:,1)
            dwdx_p(:,:) = 0.5_rprec*(dwdx(:,:,1) + dwdx(:,:,2))
            dwdy_p(:,:) = 0.5_rprec*(dwdy(:,:,1) + dwdy(:,:,2))
            dwdz_p(:,:) = dwdz(:,:,1)

#ifdef PPMAPPING
            dz_p = JACO2(1)*dz
#endif

    end select

    ! Compute inner velocity gradient product
    alpha2 = dudx_p**2 + dudy_p**2 + dudz_p**2 +          &
        dvdx_p**2 + dvdy_p**2 + dvdz_p**2 +               &
        dwdx_p**2 + dwdy_p**2 + dwdz_p**2

    ! Compute outer velocity gradient product with grid size
    b11 = (dx**2)*(dudx_p**2) + (dy**2)*(dudy_p**2) +     &
        (dz_p**2)*(dudz_p**2)
    b12 = (dx**2)*dudx_p*dvdx_p + (dy**2)*dudy_p*dvdy_p + &
        (dz_p**2)*dudz_p*dvdz_p
    b13 = (dx**2)*dudx_p*dwdx_p + (dy**2)*dudy_p*dwdy_p + &
        (dz_p**2)*dudz_p*dwdz_p
    b22 = (dx**2)*(dvdx_p**2) + (dy**2)*(dvdy_p**2) +     &
        (dz_p**2)*(dvdz_p**2)
    b23 = (dx**2)*dvdx_p*dwdx_p + (dy**2)*dvdy_p*dwdy_p + &
        (dz_p**2)*dvdz_p*dwdz_p
    b33 = (dx**2)*(dwdx_p**2) + (dy**2)*(dwdy_p**2) +     &
        (dz_p**2)*(dwdz_p**2)

    ! Compute beta product
    beta = b11*b22 - (b12**2) +                           &
        b11*b33 - (b13**2) + b22*b33 - (b23**2)

    ! Compute eddy viscosity
    ! abs here is a safety precaution, only negative if beta ~ 0
    Nu_t(:,:,1) = cvre*sqrt( abs( beta(:,:) / (alpha2(:,:)+eps) ) )

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if

! Calculate Nu_t for jz = nz (coord==nproc-1 only)
!   stored on nz-1 uvp-node (this level only) for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        case(0)
            ! Interpolate derivatives onto w-node
            dudx_p(:,:) = dudx(:,:,nz-1)
            dudy_p(:,:) = dudy(:,:,nz-1)
            dudz_p(:,:) = dudz(:,:,nz)
            dvdx_p(:,:) = dvdx(:,:,nz-1)
            dvdy_p(:,:) = dvdy(:,:,nz-1)
            dvdz_p(:,:) = dvdz(:,:,nz)
            dwdx_p(:,:) = dwdx(:,:,nz)
            dwdy_p(:,:) = dwdy(:,:,nz)
            dwdz_p(:,:) = 0.5_rprec*(dwdz(:,:,nz-1) + 0.0_rprec)

#ifdef PPMAPPING
            dz_p = JACO1(nz)*dz
#endif

        ! Wall
        case(1:)
            ! Interpolate derivatives onto uvp-node
            dudx_p(:,:) = dudx(:,:,nz-1)
            dudy_p(:,:) = dudy(:,:,nz-1)
            dudz_p(:,:) = 0.5_rprec*(dudz(:,:,nz-1) + dudz(:,:,nz))
            dvdx_p(:,:) = dvdx(:,:,nz-1)
            dvdy_p(:,:) = dvdy(:,:,nz-1)
            dvdz_p(:,:) = 0.5_rprec*(dvdz(:,:,nz-1) + dvdz(:,:,nz))
            dwdx_p(:,:) = 0.5_rprec*(dwdx(:,:,nz-1) + dwdx(:,:,nz))
            dwdy_p(:,:) = 0.5_rprec*(dwdy(:,:,nz-1) + dwdy(:,:,nz))
            dwdz_p(:,:) = dwdz(:,:,nz-1)

#ifdef PPMAPPING
            dz_p = JACO2(nz-1)*dz
#endif

        ! Wall-Model
        case(2:)
            ! Interpolate derivatives onto uvp-node
            dudx_p(:,:) = dudx(:,:,nz-1)
            dudy_p(:,:) = dudy(:,:,nz-1)
            ! dudz from wallstress, therefore no interpolation
            dudz_p(:,:) = dudz(:,:,nz)
            dvdx_p(:,:) = dvdx(:,:,nz-1)
            dvdy_p(:,:) = dvdy(:,:,nz-1)
            ! dvdz from wallstress, therefore no interpolation
            dvdz_p(:,:) = dvdz(:,:,nz)
            dwdx_p(:,:) = 0.5_rprec*(dwdx(:,:,nz-1) + dwdx(:,:,nz))
            dwdy_p(:,:) = 0.5_rprec*(dwdy(:,:,nz-1) + dwdy(:,:,nz))
            dwdz_p(:,:) = dwdz(:,:,nz-1)

#ifdef PPMAPPING
            dz_p = JACO2(nz-1)*dz
#endif

    end select

    ! Compute inner velocity gradient product
    alpha2 = dudx_p**2 + dudy_p**2 + dudz_p**2 +          &
        dvdx_p**2 + dvdy_p**2 + dvdz_p**2 +               &
        dwdx_p**2 + dwdy_p**2 + dwdz_p**2

    ! Compute outer velocity gradient product with grid size
    b11 = (dx**2)*(dudx_p**2) + (dy**2)*(dudy_p**2) +     &
        (dz_p**2)*(dudz_p**2)
    b12 = (dx**2)*dudx_p*dvdx_p + (dy**2)*dudy_p*dvdy_p + &
        (dz_p**2)*dudz_p*dvdz_p
    b13 = (dx**2)*dudx_p*dwdx_p + (dy**2)*dudy_p*dwdy_p + &
        (dz_p**2)*dudz_p*dwdz_p
    b22 = (dx**2)*(dvdx_p**2) + (dy**2)*(dvdy_p**2) +     &
        (dz_p**2)*(dvdz_p**2)
    b23 = (dx**2)*dvdx_p*dwdx_p + (dy**2)*dvdy_p*dwdy_p + &
        (dz_p**2)*dvdz_p*dwdz_p
    b33 = (dx**2)*(dwdx_p**2) + (dy**2)*(dwdy_p**2) +     &
        (dz_p**2)*(dwdz_p**2)

    ! Compute beta product
    beta = b11*b22 - (b12**2) +                           &
        b11*b33 - (b13**2) + b22*b33 - (b23**2)

    ! Compute eddy viscosity
    ! abs here is a safety precaution, only negative if beta ~ 0
    Nu_t(:,:,nz) = cvre*sqrt( abs( beta(:,:) / (alpha2(:,:)+eps) ) )

    ! since last level already calculated
    jz_max = nz-1
else
    jz_max = nz
end if

#ifdef PPMPI
! dudz calculated for 0:nz-1 (on w-nodes) except bottom process
! (only 1:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
#endif

! Calculate Nu_t for the rest of the domain
!   values are stored on w-nodes
do jz = jz_min, jz_max

    ! Interpolate derivatives onto w-nodes
    dudx_p(:,:) = 0.5_rprec*(dudx(:,:,jz) + dudx(:,:,jz-1))
    dudy_p(:,:) = 0.5_rprec*(dudy(:,:,jz) + dudy(:,:,jz-1))
    dudz_p(:,:) = dudz(:,:,jz)
    dvdx_p(:,:) = 0.5_rprec*(dvdx(:,:,jz) + dvdx(:,:,jz-1))
    dvdy_p(:,:) = 0.5_rprec*(dvdy(:,:,jz) + dvdy(:,:,jz-1))
    dvdz_p(:,:) = dvdz(:,:,jz)
    dwdx_p(:,:) = dwdx(:,:,jz)
    dwdy_p(:,:) = dwdy(:,:,jz)
    dwdz_p(:,:) = 0.5_rprec*(dwdz(:,:,jz) + dwdz(:,:,jz-1))

    ! Compute inner velocity gradient product
    alpha2 = dudx_p**2 + dudy_p**2 + dudz_p**2 +          &
        dvdx_p**2 + dvdy_p**2 + dvdz_p**2 +               &
        dwdx_p**2 + dwdy_p**2 + dwdz_p**2

    ! Compute outer velocity gradient product with grid size
#ifdef PPMAPPING
    dz_p = JACO1(jz)*dz
#endif
    b11 = (dx**2)*(dudx_p**2) + (dy**2)*(dudy_p**2) +     &
        (dz_p**2)*(dudz_p**2)
    b12 = (dx**2)*dudx_p*dvdx_p + (dy**2)*dudy_p*dvdy_p + &
        (dz_p**2)*dudz_p*dvdz_p
    b13 = (dx**2)*dudx_p*dwdx_p + (dy**2)*dudy_p*dwdy_p + &
        (dz_p**2)*dudz_p*dwdz_p
    b22 = (dx**2)*(dvdx_p**2) + (dy**2)*(dvdy_p**2) +     &
        (dz_p**2)*(dvdz_p**2)
    b23 = (dx**2)*dvdx_p*dwdx_p + (dy**2)*dvdy_p*dwdy_p + &
        (dz_p**2)*dvdz_p*dwdz_p
    b33 = (dx**2)*(dwdx_p**2) + (dy**2)*(dwdy_p**2) +     &
        (dz_p**2)*(dwdz_p**2)

    ! Compute beta product
    beta = b11*b22 - (b12**2) +                           &
        b11*b33 - (b13**2) + b22*b33 - (b23**2)

    ! Compute eddy viscosity
    ! abs here is a safety precaution, only negative if beta ~ 0
    Nu_t(:,:,jz) = cvre*sqrt( abs( beta(:,:) / (alpha2(:,:)+eps) ) )

end do

end subroutine vreman
