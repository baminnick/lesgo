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
subroutine wale
!*******************************************************************************
!
! This subroutine calculates the Wall-Adapting Local Eddy-viscosity (WALE) by
! Nicoud & Ducros, Flow, Turbulence and Combustion, 1999, 62 (2), pp.183-200.
!
use types, only : rprec
use param
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
use sgs_param, only : Nu_t
use derivatives, only : dft_direct_back_2d_n_yonlyC, dft_direct_forw_2d_n_yonlyC
#ifdef PPHYBRID
use derivatives, only : mpi_sync_hybrid
use mpi_defs, only : MPI_SYNC_DOWN
#else
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
#endif
#ifdef PPMAPPING
use sim_param, only : delta_stretch
#endif

implicit none

real(rprec), dimension(ld,ny) :: dudx_p, dudy_p, dudz_p, &
    dvdx_p, dvdy_p, dvdz_p, dwdx_p, dwdy_p, dwdz_p
real(rprec), dimension(ld,ny) :: alpha2, beta
real(rprec), dimension(ld,ny) :: S11, S12, S13, S22, S23, S33
real(rprec), dimension(ld,ny) :: O12, O13, O23
real(rprec), dimension(ld,ny) :: SS, OO, SdSd
real(rprec), dimension(ld,ny) :: SS11, SS12, SS13, SS22, SS23, SS33
real(rprec), dimension(ld,ny) :: OO11, OO12, OO13, OO22, OO23, OO33
real(rprec), dimension(ld,ny) :: Sd11, Sd12, Sd13, Sd22, Sd23, Sd33
real(rprec) :: cwale, eps
integer :: jz, jz_min, jz_max

! Set constants
cwale = 0.325_rprec*Co/0.16_rprec !! wale constant
eps = 1e-12 !! some small number in case denominator is zero

#ifdef PPHYBRID
call mpi_sync_hybrid( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
#else
#ifdef PPMPI
! dwdz calculated for 0:nz-1 (on w-nodes) except bottom process
! (only 1:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
#endif
#endif

if (fourier) then

do jz = 1, nz

    ! Calculate Nu_t for jz = 1 (coord==0 only)
    !   stored on uvp-nodes (this level only) for 'wall'
    !   stored on w-nodes (all) for 'stress free'
    ! See sgs_stag_util/calc_Sij for why derivatives are interpolated this way
    if ( (coord == 0) .and. (jz == 1) ) then
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

            ! Wall
            case(1)
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

        end select

    ! Calculate Nu_t for jz = nz (coord==nproc-1 only)
    !   stored on nz-1 uvp-node (this level only) for 'wall'
    !   stored on w-nodes (all) for 'stress free'
    elseif ((coord == nproc-1) .and. (jz == nz)) then
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

            ! Wall
            case(1)
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

        end select

    ! Calculate Nu_t for the rest of the domain
    !   values are stored on w-nodes
    else
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

    end if

    ! Currently duidxj(kx,ky,z), need duidxj(kx,y,z)
    call dft_direct_back_2d_n_yonlyC( dudx_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dudy_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dudz_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dvdx_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dvdy_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dvdz_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dwdx_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dwdy_p(:,:) )
    call dft_direct_back_2d_n_yonlyC( dwdz_p(:,:) )
    ! Only use kx=0 mode in calculations to avoid convolution

    ! Compute strain-rate tensor
    ! Remember calc_Sij is not called for this SGS model
    S11(1,:) = dudx_p(1,:)
    S12(1,:) = 0.5_rprec*(dudy_p(1,:) + dvdx_p(1,:))
    S13(1,:) = 0.5_rprec*(dudz_p(1,:) + dwdx_p(1,:))
    S22(1,:) = dvdy_p(1,:)
    S23(1,:) = 0.5_rprec*(dvdz_p(1,:) + dwdy_p(1,:))
    S33(1,:) = dwdz_p(1,:)

    ! Compute anti-symmetric part of gradient tensor
    O12(1,:) = 0.5_rprec*(dudy_p(1,:) - dvdx_p(1,:))
    O13(1,:) = 0.5_rprec*(dudz_p(1,:) - dwdx_p(1,:))
    O23(1,:) = 0.5_rprec*(dvdz_p(1,:) - dwdy_p(1,:))

    ! Compute inner product of Sij
    SS(1,:) = S11(1,:)**2 + S22(1,:)**2 + S33(1,:)**2 +       &
        2.0_rprec*(S12(1,:)**2 + S13(1,:)**2 + S23(1,:)**2)

    ! Compute inner product of Oij
    OO(1,:) = 2.0_rprec*(O12(1,:)**2 + O13(1,:)**2 + O23(1,:)**2)

    ! Compute outer product of Sij: Sik*Skj
    SS11(1,:) = S11(1,:)**2 + S12(1,:)**2 + S13(1,:)**2
    SS12(1,:) = S11(1,:)*S12(1,:) + S12(1,:)*S22(1,:) + S13(1,:)*S23(1,:)
    SS13(1,:) = S11(1,:)*S13(1,:) + S12(1,:)*S23(1,:) + S13(1,:)*S33(1,:)
    SS22(1,:) = S12(1,:)**2 + S22(1,:)**2 + S23(1,:)**2
    SS23(1,:) = S12(1,:)*S13(1,:) + S22(1,:)*S23(1,:) + S23(1,:)*S33(1,:)
    SS33(1,:) = S13(1,:)**2 + S23(1,:)**2 + S33(1,:)**2

    ! Compute outer product of Oij: Oik*Okj
    OO11(1,:) = -O12(1,:)**2 - O13(1,:)**2
    OO12(1,:) = -O13(1,:)*O23(1,:)
    OO13(1,:) = O12(1,:)*O23(1,:)
    OO22(1,:) = -O12(1,:)**2 - O23(1,:)**2
    OO23(1,:) = -O12(1,:)*O13(1,:)
    OO33(1,:) = -O13(1,:)**2 - O23(1,:)**2

    ! Compute fancy tensor: Sdij = SSij+OOij-(1/3)*kron_ij*(SS-OO)
    Sd11(1,:) = SS11(1,:) + OO11(1,:) - ((SS(1,:) - OO(1,:))/3._rprec)
    Sd12(1,:) = SS12(1,:) + OO12(1,:)
    Sd13(1,:) = SS13(1,:) + OO13(1,:)
    Sd22(1,:) = SS22(1,:) + OO22(1,:) - ((SS(1,:) - OO(1,:))/3._rprec)
    Sd23(1,:) = SS23(1,:) + OO23(1,:)
    Sd33(1,:) = SS33(1,:) + OO33(1,:) - ((SS(1,:) - OO(1,:))/3._rprec)

    ! Compute inner product of Sdij
    SdSd(1,:) = Sd11(1,:)**2 + Sd22(1,:)**2 + Sd33(1,:)**2 +           &
        2.0_rprec*(Sd12(1,:)**2 + Sd13(1,:)**2 + Sd23(1,:)**2)

    ! Compute eddy viscosity
    Nu_t(1,:,jz) = ((cwale*delta_stretch(jz))**2)*                    &
        ( (SdSd(1,:)**1.5_rprec) / ((SS(1,:)**2.5_rprec) + (SdSd(1,:)**1.25_rprec) + eps) )

    ! Zero-out non-zero kx modes
    Nu_t(2:ld,:,jz) = 0.0_rprec

    ! Transform eddy viscosity, (kx,y,z) --> (kx,ky,z)
    call dft_direct_forw_2d_n_yonlyC( Nu_t(:,:,jz) )

end do

else !! .not. fourier

do jz = 1, nz

    ! Calculate Nu_t for jz = 1 (coord==0 only)
    !   stored on uvp-nodes (this level only) for 'wall'
    !   stored on w-nodes (all) for 'stress free'
    ! See sgs_stag_util/calc_Sij for why derivatives are interpolated this way
    if ( (coord == 0) .and. (jz == 1) ) then
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

            ! Wall
            case(1)
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

        end select

    ! Calculate Nu_t for jz = nz (coord==nproc-1 only)
    !   stored on nz-1 uvp-node (this level only) for 'wall'
    !   stored on w-nodes (all) for 'stress free'
    elseif ((coord == nproc-1) .and. (jz == nz)) then
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

            ! Wall
            case(1)
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

        end select

    ! Calculate Nu_t for the rest of the domain
    !   values are stored on w-nodes
    else
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

    end if

    ! Compute strain-rate tensor
    ! Remember calc_Sij is not called for this SGS model
    S11 = dudx_p
    S12 = 0.5_rprec*(dudy_p + dvdx_p)
    S13 = 0.5_rprec*(dudz_p + dwdx_p)
    S22 = dvdy_p
    S23 = 0.5_rprec*(dvdz_p + dwdy_p)
    S33 = dwdz_p

    ! Compute anti-symmetric part of gradient tensor
    O12 = 0.5_rprec*(dudy_p - dvdx_p)
    O13 = 0.5_rprec*(dudz_p - dwdx_p)
    O23 = 0.5_rprec*(dvdz_p - dwdy_p)

    ! Compute inner product of Sij
    SS = S11**2 + S22**2 + S33**2 +               &
        2.0_rprec*(S12**2 + S13**2 + S23**2)

    ! Compute inner product of Oij
    OO = 2.0_rprec*(O12**2 + O13**2 + O23**2)

    ! Compute outer product of Sij: Sik*Skj
    SS11 = S11**2 + S12**2 + S13**2
    SS12 = S11*S12 + S12*S22 + S13*S23
    SS13 = S11*S13 + S12*S23 + S13*S33
    SS22 = S12**2 + S22**2 + S23**2
    SS23 = S12*S13 + S22*S23 + S23*S33
    SS33 = S13**2 + S23**2 + S33**2

    ! Compute outer product of Oij: Oik*Okj
    OO11 = -O12**2 - O13**2
    OO12 = -O13*O23
    OO13 = O12*O23
    OO22 = -O12**2 - O23**2
    OO23 = -O12*O13
    OO33 = -O13**2 - O23**2

    ! Compute fancy tensor: Sdij = SSij+OOij-(1/3)*kron_ij*(SS-OO)
    Sd11 = SS11 + OO11 - ((SS - OO)/3._rprec)
    Sd12 = SS12 + OO12
    Sd13 = SS13 + OO13
    Sd22 = SS22 + OO22 - ((SS - OO)/3._rprec)
    Sd23 = SS23 + OO23
    Sd33 = SS33 + OO33 - ((SS - OO)/3._rprec)

    ! Compute inner product of Sdij
    SdSd = Sd11**2 + Sd22**2 + Sd33**2 +           &
        2.0_rprec*(Sd12**2 + Sd13**2 + Sd23**2)

    ! Compute eddy viscosity
    Nu_t(:,:,jz) = ((cwale*delta_stretch(jz))**2)*                    &
        ( (SdSd**1.5_rprec) / ((SS**2.5_rprec) + (SdSd**1.25_rprec) + eps) )

end do

endif

end subroutine wale
