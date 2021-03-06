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
module sgs_stag_util
!*******************************************************************************
implicit none

save
private

public sgs_stag, rtnewt

contains

!*******************************************************************************
subroutine sgs_stag ()
!*******************************************************************************
!
! Calculates turbulent (subgrid) stress for entire domain
!   using the model specified in param.f90 (Smag, LASD, etc)
!   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
!   txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)
!   Sij values are stored on w-nodes (1:nz)
!
!   module is used to share Sij values b/w subroutines
!     (avoids memory error when arrays are very large)
!
! put everything onto w-nodes, follow original version

use types, only : rprec
use param
use sim_param, only : txx, txy, txz, tyy, tyz, tzz
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
#ifdef PPMAPPING
use sim_param, only : jaco_uv, mesh_stretch, delta_stretch
#endif
use sgs_param
use messages
#ifdef PPCNDIFF
use sim_param, only : txz_half1, txz_half2, tyz_half1, tyz_half2
#endif
use derivatives, only : convolve_rnl, dft_direct_back_2d_n_yonlyC_big
use derivatives, only : dft_direct_forw_2d_n_yonlyC_big
use derivatives, only : dft_direct_back_2d_n_yonlyC
use derivatives, only : dft_direct_forw_2d_n_yonlyC
use fft, only : padd, unpadd

#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif

#ifdef PPLVLSET
use level_set, only : level_set_BC, level_set_Cs
#endif

implicit none

character (*), parameter :: sub_name = 'sgs_stag'

real(rprec), dimension(nz) :: l, ziko, zz
integer :: jz, jz_min, jz_max

! sgs models of the form nu_t = cs_opt2*(l**2)*S
if ((sgs) .and. (sgs_model < 6)) then
    ! Cs is Smagorinsky's constant. l is a filter size (non-dim.)
    call calc_Sij ()

    ! This approximates the sum displacement during cs_count timesteps
    ! This is used with the lagrangian model only
#ifdef PPCFL_DT
    if (sgs_model == 4 .OR. sgs_model==5) then
        if ( ( jt .GE. DYN_init-cs_count + 1 ) .OR.  initu ) then
            lagran_dt = lagran_dt + dt
        endif
    endif
#else
    lagran_dt = cs_count*dt
#endif

    ! Traditional Smagorinsky model
    if (sgs_model == 1) then

#ifdef PPLVLSET
        l = delta
        call level_set_Cs (delta)
#else
        ! Parameters (Co and nn) for wallfunction defined in param.f90
        Cs_opt2 = Co**2  ! constant coefficient

        ! both Stress free
        if (lbc_mom == 0 .and. ubc_mom == 0) then
            l = delta

        ! top Stress free, bottom wall
        else if (lbc_mom > 0 .and. ubc_mom == 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == 0) then
                ! z's nondimensional, l here is on uv-nodes
#ifdef PPMAPPING
                if (damp_model == 1) then
                    zz(1) = mesh_stretch(1)
                    l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)    &
                        + (delta_stretch(1))**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                else
                    zz(1) = mesh_stretch(1)*u_star/nu_molec !! plus units
                    l(1) = delta_stretch(1)*(1.0_rprec - exp(-zz(1)/25.0_rprec)) !! A+ = 25.0
                end if
#else
                if (damp_model == 1) then
                    zz(1) = 0.5_rprec * dz
                    l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)    &
                        + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                else
                    zz(1) = 0.5_rprec * dz * u_star / nu_molec !! plus units
                    l(1) = delta*(1.0_rprec - exp(-zz(1)/25.0_rprec)) !! A+ = 25.0
                end if
#endif
                jz_min = 2
            else
                jz_min = 1
            end if

            do jz = jz_min, nz
                ! z's nondimensional, l here is on w-nodes
#ifdef PPMAPPING
                if (damp_model == 1) then
                    zz(jz) = mesh_stretch(jz)-0.5*jaco_uv(jz)*dz
                    l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)  &
                        + (delta_stretch(jz))**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                else
                    zz(jz) = (mesh_stretch(jz)-0.5*jaco_uv(jz)*dz) * u_star / nu_molec !! plus units
                    l(jz) = delta_stretch(jz)*(1.0_rprec - exp(-zz(jz)/25.0_rprec)) !! A= = 25.0
                endif
#else
                if (damp_model == 1) then
                    zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz
                    l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)  &
                        + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                else
                    zz(jz) = (((jz - 1) + coord * (nz - 1)) * dz) * u_star / nu_molec !! plus units
                    l(jz) = delta*(1.0_rprec - exp(-zz(jz)/25.0_rprec)) !! A= = 25.0
                end if
#endif
            end do

        ! both top and bottom walls, zz is distance to nearest wall
        else if (lbc_mom > 0 .and. ubc_mom > 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == 0) then
                ! z's nondimensional, l here is on uv-nodes
                zz(1) = 0.5_rprec * dz
                l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)&
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_min = 2
            else
                jz_min = 1
            end if

            if (coord == nproc-1) then
                ! z's nondimensional, l here is on uv-nodes
                zz(nz) = 0.5_rprec * dz
                l(nz) = (Co**(wall_damp_exp)*(vonk*zz(nz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_max = nz-1
            else
                jz_max = nz
            end if

            do jz = jz_min, jz_max
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz
                zz(jz) = min( zz(jz), (nz-1)*nproc*dz - zz(jz) )
                l(jz) = (Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! top wall, bottom stress free, zz is distance to top
        else if (lbc_mom == 0 .and. ubc_mom > 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == nproc-1) then
                ! z's nondimensional, l here is on uv-nodes
                zz(nz) = 0.5_rprec * dz
                l(nz) = (Co**(wall_damp_exp)*(vonk*zz(nz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_max = nz-1
            else
                jz_max = nz
            end if

            do jz = 1, jz_max
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((nproc - coord)*(nz - 1) - (jz - 1)) * dz
                l(jz) = (Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! Invalid combination
        else
            call error (sub_name, 'invalid b.c. combination')
        end if
#endif

    ! Dynamic procedures: modify/set Sij and Cs_opt2 (specific to sgs_model)
    else
        ! recall: l is the filter size
#ifdef PPMAPPING
        l(:) = delta_stretch(1:nz) ! delta_stretch from lbz:nz, l from 1:nz
#else
        l = delta
#endif

        ! Use the Smagorinsky model until DYN_init timestep
        if ((jt == 1) .and. (inilag)) then
#ifdef PPOUTPUT_SGS
            if (coord == 0) write(*,*) 'CS_opt2 initialiazed'
#endif
            Cs_opt2 = 0.03_rprec

        ! Update Sij, Cs every cs_count timesteps (specified in param)
        elseif ( ((jt.GE.DYN_init).OR.(initu)) .AND.                           &
            (mod(jt_total,cs_count)==0) ) then

#ifdef PPOUTPUT_SGS
            if ((jt == DYN_init) .and. (coord == 0)) then
                write(*,*) 'running dynamic sgs_model = ', sgs_model
            end if
#endif

            ! Standard dynamic model
            if (sgs_model == 2) then
                call std_dynamic(ziko)
                forall (jz = 1:nz) Cs_opt2(:, :, jz) = ziko(jz)

            ! Plane average dynamic model
            else if (sgs_model==3) then
                call scaledep_dynamic(ziko)
                do jz = 1, nz
                    Cs_opt2(:, :, jz) = ziko(jz)
                end do
            ! Lagrangian scale similarity model
            else if (sgs_model==4.) then
                call lagrange_Ssim()

            ! Lagrangian scale dependent model
            elseif (sgs_model==5) then
                call lagrange_Sdep()
            end if

        end if
    end if !! for Smagorinsky/Dynamic SGS

    ! Define |S| and eddy viscosity (nu_t= c_s^2 l^2 |S|) for entire domain
    !   stored on w-nodes (on uvp node for jz=1 or nz for 'wall' BC only)
    if (fourier) then
        do jz = 1, nz

            ! Remember Sij(kx,y,z) at the end of calc_Sij
            ! Use only kx = 0 when computing strain-rate magnitude (SRM)
            ! Assuming streamwise average
            S(1,:) = sqrt( 2.0_rprec*(S11(1,:,jz)**2 +          &
                S22(1,:,jz)**2 + S33(1,:,jz)**2 +               &
                2.0_rprec*(S12(1,:,jz)**2 + S13(1,:,jz)**2 +    &
                S23(1,:,jz)**2 )))
            S(2:ld,:) = 0.0_rprec

            ! Transform SRM y --> ky
            ! No need to transform Sij, gets overwritten in calc_Sij
            call dft_direct_forw_2d_n_yonlyC( S(:,:) )

            ! Commented code below used Bretheim et al. 2018, assumes Sij(kx,ky,z) at this point
            !! Padd to prepare for convolution
            !call padd(S11_big(:,:,jz), S11(:,:,jz))
            !call padd(S22_big(:,:,jz), S22(:,:,jz))
            !call padd(S33_big(:,:,jz), S33(:,:,jz))
            !call padd(S12_big(:,:,jz), S12(:,:,jz))
            !call padd(S13_big(:,:,jz), S13(:,:,jz))
            !call padd(S23_big(:,:,jz), S23(:,:,jz))

            !! ky --> y
            !call dft_direct_back_2d_n_yonlyC_big(S11_big(:,:,jz))
            !call dft_direct_back_2d_n_yonlyC_big(S22_big(:,:,jz))
            !call dft_direct_back_2d_n_yonlyC_big(S33_big(:,:,jz))
            !call dft_direct_back_2d_n_yonlyC_big(S12_big(:,:,jz))
            !call dft_direct_back_2d_n_yonlyC_big(S13_big(:,:,jz))
            !call dft_direct_back_2d_n_yonlyC_big(S23_big(:,:,jz))

            !! Convolve to find strain-rate magnitude
            !S_big(:,:) = 2._rprec*(convolve_rnl(S11_big(:,:,jz),S11_big(:,:,jz)) + &
            !                       convolve_rnl(S22_big(:,:,jz),S22_big(:,:,jz)) + &
            !                       convolve_rnl(S33_big(:,:,jz),S33_big(:,:,jz)) + &
            !             2._rprec*(convolve_rnl(S12_big(:,:,jz),S12_big(:,:,jz)) + &
            !                       convolve_rnl(S13_big(:,:,jz),S13_big(:,:,jz)) + &
            !                       convolve_rnl(S23_big(:,:,jz),S23_big(:,:,jz))))

            !! Code below uses a streamwise average instead of streamwise and 
            !! spanwise average as used in Bretheim et al. 2018
            !S_big(1,1:ny2) = sqrt( abs( S_big(1,1:ny2) ) )
            !! S_big(1,1:ny2) = sqrt( S_big(1,1:ny2) ) !! sometimes gives sqrt(negative number)!
            !S_big(2:ld,1:ny2) = 0._rprec 

            !! y --> ky
            !call dft_direct_forw_2d_n_yonlyC_big( S_big(:,:) )
            !call unpadd( S(:,:), S_big(:,:) )

            ! ! Commented code below is the SGS model used in Bretheim et al. 2018
            ! ! It uses a streamwise and spanwise average of the strain-rate magnitude
            ! ! Only consider streamwise and spanwise average
            ! S(1,1) = sqrt( S(1,1) )
            ! S(2:ld, 1:ny) = 0._rprec
            ! S(1,2:ny) = 0._rprec

            ! Compute eddy viscosity
            Nu_t(:,:,jz) = S(:,:)*Cs_opt2(:,:,jz)*l(jz)**2

        end do
    else !! not fourier
        do jz = 1, nz
            S(1:nx,:) = sqrt( 2.0_rprec*(S11(1:nx,:,jz)**2 +          &
                S22(1:nx,:,jz)**2 + S33(1:nx,:,jz)**2 +               &
                2.0_rprec*(S12(1:nx,:,jz)**2 + S13(1:nx,:,jz)**2 +    &
                S23(1:nx,:,jz)**2 )))
            Nu_t(1:nx,:,jz) = S(1:nx,:)*Cs_opt2(1:nx,:,jz)*l(jz)**2
        end do
    endif

! sgs models of the form, nu_t = c*OP(duidxj) for some operator OP
! Vreman (2004) sgs model
elseif ((sgs) .and. (sgs_model == 6)) then

    call vreman()

elseif ((sgs) .and. (sgs_model == 7)) then

    if ( ((jt.GE.DYN_init).OR.(initu)) .AND. (mod(jt_total,cs_count)==0) ) then
        call calc_Sij ()
        call dyn_vreman(cvre)
    else
        call vreman()
    endif

elseif ((sgs) .and. (sgs_model == 8)) then

    call wale()

else !! not sgs, molec only

    ! define nu_coefs here since it does not change case-by-case for DNS
    ! if sgs, nu_coefs defined case-by-case since it is dependent on Nu_t
    if (.not. sgs) then
        nu_coef(1:nx,:) = 2.0_rprec*nu
        nu_coef2(1:nx,:) = 2.0_rprec*nu
    endif

end if !! if (sgs)

! Padd velocity gradients before convolution with Nu_t
if ((fourier) .and. (sgs)) then !! RNL-LES
    do jz = 1, nz
        call padd( dudx_big(:,:,jz), dudx(:,:,jz) )
        call padd( dudy_big(:,:,jz), dudy(:,:,jz) )
        call padd( dudz_big(:,:,jz), dudz(:,:,jz) )
        call padd( dvdx_big(:,:,jz), dvdx(:,:,jz) )
        call padd( dvdy_big(:,:,jz), dvdy(:,:,jz) )
        call padd( dvdz_big(:,:,jz), dvdz(:,:,jz) )
        call padd( dwdx_big(:,:,jz), dwdx(:,:,jz) )
        call padd( dwdy_big(:,:,jz), dwdy(:,:,jz) )
        call padd( dwdz_big(:,:,jz), dwdz(:,:,jz) )

        ! ky --> y
        call dft_direct_back_2d_n_yonlyC_big(dudx_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dudy_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dudz_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dvdx_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dvdy_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dvdz_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dwdx_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dwdy_big(:,:,jz))
        call dft_direct_back_2d_n_yonlyC_big(dwdz_big(:,:,jz))
    end do
end if

! Calculate txx, txy, tyy, tzz for bottom level: jz=1 node (coord==0 only)
if (coord == 0) then
    if (sgs) then
        select case (lbc_mom)
            ! Stress free
            case (0)
                nu_coef(1:nx,:) = 2.0_rprec*0.5_rprec*(Nu_t(1:nx,:,1) + Nu_t(1:nx,:,2)) + nu

            ! Wall
            case (1:)
                if (fourier) then
                    nu_coef(:,:) = 2._rprec*Nu_t(:,:,1) !! Initialize
                    ! overwrite and add nu only to the mean --> (kx,ky) = (0,0)
                    nu_coef(1,1) = 2._rprec*(Nu_t(1,1,1)+nu)
                else
                    nu_coef(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,1)+nu) !! Nu_t on uvp-node(1) here
                end if

        end select
    endif
    if ((fourier) .and. (sgs)) then !! only convolve if sgs since nu_coef is a function of kx
        call padd( nu_coef_big(:,:), nu_coef(:,:) )
        call dft_direct_back_2d_n_yonlyC_big( nu_coef_big(:,:) )

        txx_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dudx_big(:,:,1) )
        txy_big(:,:) = convolve_rnl( -nu_coef_big(:,:), 0.5_rprec*(dudy_big(:,:,1)+dvdx_big(:,:,1)) )
        tyy_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dvdy_big(:,:,1) )
        tzz_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dwdz_big(:,:,1) )

        call dft_direct_forw_2d_n_yonlyC_big( txx_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( txy_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyy_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tzz_big(:,:) )

        call unpadd( txx(:,:,1), txx_big(:,:) )
        call unpadd( txy(:,:,1), txy_big(:,:) )
        call unpadd( tyy(:,:,1), tyy_big(:,:) )
        call unpadd( tzz(:,:,1), tzz_big(:,:) )

    else
        txx(1:nx,:,1) = -nu_coef(1:nx,:)*dudx(1:nx,:,1) !! uvp-node(1)
        txy(1:nx,:,1) = -nu_coef(1:nx,:)*(0.5_rprec*(dudy(1:nx,:,1)+dvdx(1:nx,:,1))) !! uvp-node(1)
        tyy(1:nx,:,1) = -nu_coef(1:nx,:)*dvdy(1:nx,:,1) !! uvp-node(1)
        tzz(1:nx,:,1) = -nu_coef(1:nx,:)*dwdz(1:nx,:,1) !! uvp-node(1)
    endif

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if

! Calculate txx, txy, tyy, tzz for top level: jz=nz node (coord==nproc-1)
if (coord == nproc-1) then
    if (sgs) then
        select case (ubc_mom)
            ! Stress free
            case (0)
                if (fourier) then
                    nu_coef(:,:) = 2.0_rprec*(0.5_rprec*(Nu_t(:,:,nz-1) + Nu_t(:,:,nz))) !! initialize
                    nu_coef2(:,:) = 2.0_rprec*Nu_t(:,:,nz-1) !! initialize
                    ! overwrite and add nu only to the mean --> (kx,ky) = (0,0)
                    nu_coef(1,1) = 2.0_rprec*(0.5_rprec*(Nu_t(1,1,nz-1) + Nu_t(1,1,nz)) + nu) !! uvp-node(nz-1)
                    nu_coef2(1,1) = 2.0_rprec*(Nu_t(1,1,nz-1) + nu) !! w-node(nz-1)
                else
                    nu_coef(1:nx,:) = 2.0_rprec*(0.5_rprec*(Nu_t(1:nx,:,nz-1) + Nu_t(1:nx,:,nz)) + nu) !! uvp-node(nz-1)
                    nu_coef2(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,nz-1) + nu) !! w-node(nz-1)
                end if

            ! Wall
            case (1:)
                nu_coef(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,nz) + nu) !! uvp-node
                nu_coef2(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,nz-1) + nu) !! w-node

        end select
    endif

    if ((fourier) .and. (sgs)) then !! only convolve if sgs since nu_coef is a function of kx
        call padd( nu_coef_big(:,:), nu_coef(:,:) )
        call padd( nu_coef2_big(:,:), nu_coef2(:,:) )
        call dft_direct_back_2d_n_yonlyC_big( nu_coef_big(:,:) )
        call dft_direct_back_2d_n_yonlyC_big( nu_coef2_big(:,:) )

        txx_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dudx_big(:,:,nz-1) )
        txy_big(:,:) = convolve_rnl( -nu_coef_big(:,:), 0.5_rprec*(dudy_big(:,:,nz-1)+dvdx_big(:,:,nz-1)) )
        tyy_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dvdy_big(:,:,nz-1) )
        tzz_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dwdz_big(:,:,nz-1) )
        ! for top wall, include w-grid stress since we touched nz-1
#ifdef PPCNDIFF
        txz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dudz_big(:,:,nz-1)+dwdx_big(:,:,nz-1)) )
        txz_half1_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dwdx_big(:,:,nz-1)) )
        txz_half2_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dudz_big(:,:,nz-1)) )
        tyz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dvdz_big(:,:,nz-1)+dwdy_big(:,:,nz-1)) )
        tyz_half1_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dwdy_big(:,:,nz-1)) )
        tyz_half2_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dvdz_big(:,:,nz-1)) )
#else
        txz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dudz_big(:,:,nz-1)+dwdx_big(:,:,nz-1)) )
        tyz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dvdz_big(:,:,nz-1)+dwdy_big(:,:,nz-1)) )
#endif

        call dft_direct_forw_2d_n_yonlyC_big( txx_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( txy_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyy_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tzz_big(:,:) )
#ifdef PPCNDIFF
        call dft_direct_forw_2d_n_yonlyC_big( txz_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( txz_half1_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( txz_half2_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_half1_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_half2_big(:,:) )
#else
        call dft_direct_forw_2d_n_yonlyC_big( txz_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_big(:,:) )
#endif

        call unpadd( txx(:,:,nz-1), txx_big(:,:) )
        call unpadd( txy(:,:,nz-1), txy_big(:,:) )
        call unpadd( tyy(:,:,nz-1), tyy_big(:,:) )
        call unpadd( tzz(:,:,nz-1), tzz_big(:,:) )
#ifdef PPCNDIFF
        call unpadd( txz(:,:,nz-1), txz_big(:,:) )
        call unpadd( txz_half1(:,:,nz-1), txz_half1_big(:,:) )
        call unpadd( txz_half2(:,:,nz-1), txz_half2_big(:,:) )
        call unpadd( tyz(:,:,nz-1), tyz_big(:,:) )
        call unpadd( tyz_half1(:,:,nz-1), tyz_half1_big(:,:) )
        call unpadd( tyz_half2(:,:,nz-1), tyz_half2_big(:,:) )
#else
        call unpadd( txz(:,:,nz-1), txz_big(:,:) )
        call unpadd( tyz(:,:,nz-1), tyz_big(:,:) )
#endif

    else
        txx(1:nx,:,nz-1) = -nu_coef(1:nx,:)*dudx(1:nx,:,nz-1) !! uvp-node(nz-1)
        txy(1:nx,:,nz-1) = -nu_coef(1:nx,:)*(0.5_rprec*(dudy(1:nx,:,nz-1)+dvdx(1:nx,:,nz-1))) !! uvp-node(nz-1)
        tyy(1:nx,:,nz-1) = -nu_coef(1:nx,:)*dvdy(1:nx,:,nz-1) !! uvp-node(nz-1)
        tzz(1:nx,:,nz-1) = -nu_coef(1:nx,:)*dwdz(1:nx,:,nz-1) !! uvp-node(nz-1)
        ! for top wall, include w-grid stress since we touched nz-1
#ifdef PPCNDIFF
        txz(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,nz-1)+dwdx(1:nx,:,nz-1))) !! w-node(nz-1)
        txz_half1(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dwdx(1:nx,:,nz-1))) !! w-node(nz-1)
        txz_half2(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,nz-1))) !! w-node(nz-1)
        tyz(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,nz-1)+dwdy(1:nx,:,nz-1))) !! w-node(nz-1)
        tyz_half1(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dwdy(1:nx,:,nz-1))) !! w-node(nz-1)
        tyz_half2(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,nz-1))) !! w-node(nz-1)
#else
        txz(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,nz-1)+dwdx(1:nx,:,nz-1))) !! w-node(nz-1)
        tyz(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,nz-1)+dwdy(1:nx,:,nz-1))) !! w-node(nz-1)
#endif
    end if

    ! since last level already calculated
    jz_max = nz-2
else
    jz_max = nz-1
end if

! Calculate all tau for the rest of the domain
!   txx, txy, tyy, tzz not needed at nz (so they aren't calculated)
!     txz, tyz at nz will be done later
!   txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)
do jz = jz_min, jz_max
    if (sgs) then
        if (fourier) then
            nu_coef(:,:) = 2.0_rprec*(0.5_rprec*(Nu_t(:,:,jz) + Nu_t(:,:,jz+1))) !! initialize
            nu_coef2(:,:) = 2.0_rprec*Nu_t(:,:,jz) !! initialize
            ! add nu only to the mean --> (kx,ky) = (0,0)
            nu_coef(1,1) = 2.0_rprec*(0.5_rprec*(Nu_t(1,1,jz) + Nu_t(1,1,jz+1)) + nu) !! uvp-node(jz)
            nu_coef2(1,1) = 2.0_rprec*(Nu_t(1,1,jz) + nu) !! w-node(jz)
        else
            nu_coef(1:nx,:) = 2.0_rprec*(0.5_rprec*(Nu_t(1:nx,:,jz) + Nu_t(1:nx,:,jz+1)) + nu) !! uvp-node(jz)
            nu_coef2(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,jz) + nu) !! w-node(jz)
        endif
    endif
    if ((fourier) .and. (sgs)) then !! only convolve if sgs since nu_coef is a function of kx
        call padd( nu_coef_big(:,:), nu_coef(:,:) )
        call padd( nu_coef2_big(:,:), nu_coef2(:,:) )
        call dft_direct_back_2d_n_yonlyC_big( nu_coef_big(:,:) )
        call dft_direct_back_2d_n_yonlyC_big( nu_coef2_big(:,:) )

        txx_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dudx_big(:,:,jz) )
        txy_big(:,:) = convolve_rnl( -nu_coef_big(:,:), 0.5_rprec*(dudy_big(:,:,jz)+dvdx_big(:,:,jz)) )
        tyy_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dvdy_big(:,:,jz) )
        tzz_big(:,:) = convolve_rnl( -nu_coef_big(:,:), dwdz_big(:,:,jz) )
        ! for top wall, include w-grid stress since we touched nz-1
#ifdef PPCNDIFF
        txz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dudz_big(:,:,jz)+dwdx_big(:,:,jz)) )
        txz_half1_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dwdx_big(:,:,jz)) )
        txz_half2_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dudz_big(:,:,jz)) )
        tyz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dvdz_big(:,:,jz)+dwdy_big(:,:,jz)) )
        tyz_half1_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dwdy_big(:,:,jz)) )
        tyz_half2_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dvdz_big(:,:,jz)) )
#else
        txz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dudz_big(:,:,jz)+dwdx_big(:,:,jz)) )
        tyz_big(:,:) = convolve_rnl( -nu_coef2_big(:,:), 0.5_rprec*(dvdz_big(:,:,jz)+dwdy_big(:,:,jz)) )
#endif

        call dft_direct_forw_2d_n_yonlyC_big( txx_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( txy_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyy_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tzz_big(:,:) )
#ifdef PPCNDIFF
        call dft_direct_forw_2d_n_yonlyC_big( txz_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( txz_half1_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( txz_half2_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_half1_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_half2_big(:,:) )
#else
        call dft_direct_forw_2d_n_yonlyC_big( txz_big(:,:) )
        call dft_direct_forw_2d_n_yonlyC_big( tyz_big(:,:) )
#endif

        call unpadd( txx(:,:,jz), txx_big(:,:) )
        call unpadd( txy(:,:,jz), txy_big(:,:) )
        call unpadd( tyy(:,:,jz), tyy_big(:,:) )
        call unpadd( tzz(:,:,jz), tzz_big(:,:) )
#ifdef PPCNDIFF
        call unpadd( txz(:,:,jz), txz_big(:,:) )
        call unpadd( txz_half1(:,:,jz), txz_half1_big(:,:) )
        call unpadd( txz_half2(:,:,jz), txz_half2_big(:,:) )
        call unpadd( tyz(:,:,jz), tyz_big(:,:) )
        call unpadd( tyz_half1(:,:,jz), tyz_half1_big(:,:) )
        call unpadd( tyz_half2(:,:,jz), tyz_half2_big(:,:) )
#else
        call unpadd( txz(:,:,jz), txz_big(:,:) )
        call unpadd( tyz(:,:,jz), tyz_big(:,:) )
#endif

    else
        txx(1:nx,:,jz)=-nu_coef(1:nx,:)*dudx(1:nx,:,jz) !! uvp-node(jz)
        txy(1:nx,:,jz)=-nu_coef(1:nx,:)*(0.5_rprec*(dudy(1:nx,:,jz)+dvdx(1:nx,:,jz))) !! uvp-node(jz)
        tyy(1:nx,:,jz)=-nu_coef(1:nx,:)*dvdy(1:nx,:,jz) !! uvp-node(jz)
        tzz(1:nx,:,jz)=-nu_coef(1:nx,:)*dwdz(1:nx,:,jz) !! uvp-node(jz)
#ifdef PPCNDIFF
        txz(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,jz)+dwdx(1:nx,:,jz))) !! w-node(jz)
        txz_half1(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dwdx(1:nx,:,jz))) !! w-node(jz)
        txz_half2(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,jz))) !! w-node(jz)
        tyz(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,jz)+dwdy(1:nx,:,jz))) !! w-node(jz)
        tyz_half1(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dwdy(1:nx,:,jz))) !! w-node(jz)
        tyz_half2(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,jz))) !! w-node(jz)
#else
        txz(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,jz)+dwdx(1:nx,:,jz))) !! w-node(jz)
        tyz(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,jz)+dwdy(1:nx,:,jz))) !! w-node(jz)
#endif
    endif
enddo

#ifdef PPLVLSET
!--at this point tij are only set for 1:nz-1
!--at this point u, v, w are set for 0:nz, except bottom process is 1:nz
!--some MPI synchronizing may be done in here, but this will be kept
!  separate from the rest of the code (at the risk of some redundancy)
call level_set_BC ()
#endif

#ifdef PPMPI
! txz,tyz calculated for 1:nz-1 (on w-nodes) except bottom process
! (only 2:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
#ifdef PPCNDIFF
call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( txz_half1, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( txz_half2, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( tyz_half1, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( tyz_half2, 0, MPI_SYNC_DOWN )
#else
call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
#endif
#ifdef PPSAFETYMODE
! Set bogus values (easier to catch if there's an error)
txx(:, :, 0) = BOGUS
txy(:, :, 0) = BOGUS
#ifdef PPCNDIFF
txz(:, :, 0) = BOGUS
txz_half1(:, :, 0) = BOGUS
txz_half2(:, :, 0) = BOGUS
tyz(:, :, 0) = BOGUS
tyz_half1(:, :, 0) = BOGUS
tyz_half2(:, :, 0) = BOGUS
#else
txz(:, :, 0) = BOGUS
tyz(:, :, 0) = BOGUS
#endif
tyy(:, :, 0) = BOGUS
tzz(:, :, 0) = BOGUS
#endif
#endif

! Set bogus values (easier to catch if there's an error)
#ifdef PPSAFETYMODE
txx(:, :, nz) = BOGUS
txy(:, :, nz) = BOGUS
tyy(:, :, nz) = BOGUS
tzz(:, :, nz) = BOGUS
#endif

end subroutine sgs_stag

!*******************************************************************************
subroutine calc_Sij()
!*******************************************************************************
! Calculate the resolved strain rate tensor, Sij = 0.5(djui + diuj)
!   values are stored on w-nodes (1:nz) except for near the wall, then stored
!   on uv1 node
! 
! Sij is only used to calculated nu_t, not tau_ij
! 

use types, only : rprec
use param
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
use sgs_param
use derivatives, only : dft_direct_back_2d_n_yonlyC
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
implicit none

integer :: jz, jz_min, jz_max

! Calculate Sij for jz=1 (coord==0 only)
!   stored on uvp-nodes (this level only) for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free
        case (0)
            ! Sij values are supposed to be on w-nodes for this case
            ! dudx, dudy, dvdx, dvdy at boundary (w) are the same as
            !   on closest uvp because of slip
            ! dwdz should be 0 beyond boundary due to no-penetration
            S11(1:nx,:,1) = dudx(1:nx,:,1)
            S12(1:nx,:,1) = 0.5_rprec*(dudy(1:nx,:,1)+dvdx(1:nx,:,1))
            S13(1:nx,:,1) = 0.5_rprec*(dudz(1:nx,:,1)+dwdx(1:nx,:,1))
            S22(1:nx,:,1) = dvdy(1:nx,:,1)
            S23(1:nx,:,1) = 0.5_rprec*(dvdy(1:nx,:,1)+dwdy(1:nx,:,1))
            S33(1:nx,:,1) = 0.5_rprec*(dwdz(1:nx,:,1) + 0.0_rprec) ! interpolating here

        ! Wall
        ! dudz, dvdz, dwdx, and dwdy on w-nodes... need to be interpolated
        case (1)
            ! these values stored on uvp-nodes
            S11(1:nx,:,1) = dudx(1:nx,:,1) !! uvp_node(1)
            S12(1:nx,:,1) = 0.5_rprec*(dudy(1:nx,:,1)+dvdx(1:nx,:,1)) !! uvp_node(1)
            S13(1:nx,:,1) = 0.5_rprec*( (0.5_rprec*(dudz(1:nx,:,1)+dudz(1:nx,:,2))) + & 
                (0.5_rprec*(dwdx(1:nx,:,1)+dwdx(1:nx,:,2))) ) !! uvp_node(1)
            S22(1:nx,:,1) = dvdy(1:nx,:,1) !! uvp_node(1)
            S23(1:nx,:,1) = 0.5_rprec*( (0.5_rprec*(dvdz(1:nx,:,1)+dvdz(1:nx,:,2))) + &
                (0.5_rprec*(dwdy(1:nx,:,1)+dwdy(1:nx,:,2))) ) !! uvp_node(1)
            S33(1:nx,:,1) = dwdz(1:nx,:,1) !! uvp_node(1)

        ! Wall-Model
        ! Only differs from case (1) by definition of dudz and dvdz
        ! before these had to be interpolated, however if wall-modeling,
        ! dudz and dvdz are defined on the first grid point (uvp) from the wall
        case (2:)
            ! these values stored on uvp-nodes
            S11(1:nx,:,1) = dudx(1:nx,:,1) !! uvp_node(1)
            S12(1:nx,:,1) = 0.5_rprec*(dudy(1:nx,:,1)+dvdx(1:nx,:,1)) !! uvp_node(1)
            S13(1:nx,:,1) = 0.5_rprec*(dudz(1:nx,:,1) +                            & 
                (0.5_rprec*(dwdx(1:nx,:,1)+dwdx(1:nx,:,2))) ) !! uvp_node(1)
            S22(1:nx,:,1) = dvdy(1:nx,:,1) !! uvp_node(1)
            S23(1:nx,:,1) = 0.5_rprec*(dvdz(1:nx,:,1) +                            &
                (0.5_rprec*(dwdy(1:nx,:,1)+dwdy(1:nx,:,2))) ) !! uvp_node(1)
            S33(1:nx,:,1) = dwdz(1:nx,:,1) !! uvp_node(1)

    end select

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if

! Calculate Sij for jz=nz (coord==nproc-1 only)
!   stored on uvp-nodes (this level only nz on w-grid --> nz-1 on uvp-grid)
!       for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        case (0)
            ! Sij values are supposed to be on w-nodes for this case
            ! dudx, dudy, dvdx, dvdy at boundary (w) are the same as
            !   on closest uvp because of slip
            ! dwdz should be 0 beyond boundary due to no-penetration
            S11(1:nx,:,nz) = dudx(1:nx,:,nz-1)
            S12(1:nx,:,nz) = 0.5_rprec*(dudy(1:nx,:,nz-1)+dvdx(1:nx,:,nz-1))
            S13(1:nx,:,nz) = 0.5_rprec*(dudz(1:nx,:,nz)+dwdx(1:nx,:,nz))
            S22(1:nx,:,nz) = dvdy(1:nx,:,nz-1)
            S23(1:nx,:,nz) = 0.5_rprec*(dvdz(1:nx,:,nz)+dwdy(1:nx,:,nz))
            S33(1:nx,:,nz) = 0.5_rprec*(dwdz(1:nx,:,nz-1) + 0.0_rprec)

        ! Wall
        ! dudz, dvdz, dwdx, and dwdy on w-nodes... need to be interpolated
        case (1)
            ! these values stored on uvp-nodes
            S11(1:nx,:,nz) = dudx(1:nx,:,nz-1) !! uvp_node(nz-1)
            S12(1:nx,:,nz) = 0.5_rprec*(dudy(1:nx,:,nz-1)+dvdx(1:nx,:,nz-1)) !! uvp_node(nz-1)
            S13(1:nx,:,nz) = 0.5_rprec*( (0.5_rprec*(dudz(1:nx,:,nz-1)+dudz(1:nx,:,nz))) +   &
                (0.5_rprec*(dwdx(1:nx,:,nz-1)+dwdx(1:nx,:,nz))) ) !! uvp_node(nz-1)
            S22(1:nx,:,nz) = dvdy(1:nx,:,nz-1) !! uvp_node(nz-1)
            S23(1:nx,:,nz) = 0.5_rprec*( (0.5_rprec*(dvdz(1:nx,:,nz-1)+dvdz(1:nx,:,nz))) +   &
                (0.5_rprec*(dwdy(1:nx,:,nz-1)+dwdy(1:nx,:,nz))) ) !! uvp_node(nz-1)
            S33(1:nx,:,nz) = dwdz(1:nx,:,nz-1) !! uvp_node(nz-1)

        ! Wall-Model
        ! Only differs from case (1) by definition of dudz and dvdz
        ! before these had to be interpolated, however if wall-modeling,
        ! dudz and dvdz are defined on the first grid point (uvp) from the wall
        case (2:)
            ! these values stored on uvp-nodes
            S11(1:nx,:,nz) = dudx(1:nx,:,nz-1) !! uvp_node(nz-1)
            S12(1:nx,:,nz) = 0.5_rprec*(dudy(1:nx,:,nz-1)+dvdx(1:nx,:,nz-1)) !! uvp_node(nz-1)
            ! dudz from wallstress()
            S13(1:nx,:,nz) = 0.5_rprec*(dudz(1:nx,:,nz) +                                    &
                (0.5_rprec*(dwdx(1:nx,:,nz-1)+dwdx(1:nx,:,nz))) ) !! uvp_node(nz-1)
            S22(1:nx,:,nz) = dvdy(1:nx,:,nz-1) !! uvp_node(nz-1)
            ! dvdz from wallstress()
            S23(1:nx,:,nz) = 0.5_rprec*(dvdz(1:nx,:,nz) +                                    &
                (0.5_rprec*(dwdy(1:nx,:,nz-1)+dwdy(1:nx,:,nz))) ) !! uvp_node(nz-1)
            S33(1:nx,:,nz) = dwdz(1:nx,:,nz-1) !! uvp_node(nz-1)

    end select

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

! Calculate Sij for the rest of the domain
!   values are stored on w-nodes
!   dudz, dvdz, dwdx, dwdy are already stored on w-nodes
do jz = jz_min, jz_max
    S11(1:nx,:,jz) = 0.5_rprec*(dudx(1:nx,:,jz) + dudx(1:nx,:,jz-1)) !! w-node(jz)
    S12(1:nx,:,jz) = 0.25_rprec*(dudy(1:nx,:,jz) + dudy(1:nx,:,jz-1) +           &
        dvdx(1:nx,:,jz) + dvdx(1:nx,:,jz-1))                         !! w-node(jz)
    S13(1:nx,:,jz) = 0.5_rprec*(dudz(1:nx,:,jz) + dwdx(1:nx,:,jz))   !! w-node(jz)
    S22(1:nx,:,jz) = 0.5_rprec*(dvdy(1:nx,:,jz) + dvdy(1:nx,:,jz-1)) !! w-node(jz)
    S23(1:nx,:,jz) = 0.5_rprec*(dvdz(1:nx,:,jz) + dwdy(1:nx,:,jz))   !! w-node(jz)
    S33(1:nx,:,jz) = 0.5_rprec*(dwdz(1:nx,:,jz) + dwdz(1:nx,:,jz-1)) !! w-node(jz)
end do

if (fourier) then
    do jz = 1, nz
        ! Transform ky --> y
        call dft_direct_back_2d_n_yonlyC( S11(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( S22(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( S33(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( S12(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( S13(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC( S23(:,:,jz) )
    end do
endif

end subroutine calc_Sij

!*******************************************************************************
real(rprec) function rtnewt(A, jz)
!*******************************************************************************
use types, only : rprec
integer, parameter :: jmax=100
real(rprec) :: x1, x2, xacc
real(rprec) :: df,dx,f
integer :: j, jz
real(rprec), dimension(0:5) :: A
x1 = 0._rprec
x2 = 15._rprec                  ! try to find the largest root first
xacc = 0.001_rprec              ! doesn't need to be that accurate
rtnewt = 0.5_rprec*(x1+x2)
do j = 1, jmax
    f = A(0)+rtnewt*(A(1)+rtnewt*(A(2)+rtnewt*(A(3)+rtnewt*(A(4)+rtnewt*A(5)))))
    df = A(1) + rtnewt*(2._rprec*A(2) + rtnewt*(3._rprec*A(3) +                &
        rtnewt*(4._rprec*A(4) + rtnewt*(5._rprec*A(5)))))
    dx = f/df
    rtnewt = rtnewt - dx
    if (abs(dx) < xacc) return
end do
rtnewt = 1._rprec  ! if don't converge fast enough
write(6,*) 'using beta=1 at jz= ', jz

end function rtnewt

end module sgs_stag_util
