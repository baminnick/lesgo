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
use sim_param, only : JACO2, mesh_stretch, delta_stretch
#endif
use sgs_param
use messages
#ifdef PPCNDIFF
use sim_param, only : txz_half1, txz_half2
#endif

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

if (sgs) then
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
                zz(1) = mesh_stretch(1)
                l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)    &
                    + (delta_stretch(1))**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
#else
                zz(1) = 0.5_rprec * dz
                l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)    &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
#endif
                jz_min = 2
            else
                jz_min = 1
            end if

            do jz = jz_min, nz
                ! z's nondimensional, l here is on w-nodes
#ifdef PPMAPPING
                zz(jz) = mesh_stretch(jz)-0.5*JACO2(jz)*dz
                l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)  &
                    + (delta_stretch(jz))**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
#else
                zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz
                l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)  &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
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
        l = delta

        ! Use the Smagorinsky model until DYN_init timestep
        if ((jt == 1) .and. (inilag)) then
#ifdef PPOUTPUT_SGS
            write(*,*) 'CS_opt2 initialiazed'
#endif
            Cs_opt2 = 0.03_rprec

        ! Update Sij, Cs every cs_count timesteps (specified in param)
        elseif ( ((jt.GE.DYN_init).OR.(initu)) .AND.                           &
            (mod(jt_total,cs_count)==0) ) then

#ifdef PPOUTPUT_SGS
            if (jt == DYN_init) then
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
    do jz = 1, nz
        S(1:nx,:) = sqrt( 2.0_rprec*(S11(1:nx,:,jz)**2 + S22(1:nx,:,jz)**2 +       &
            S33(1:nx,:,jz)**2 + 2.0_rprec*(S12(1:nx,:,jz)**2 +                     &
            S13(1:nx,:,jz)**2 + S23(1:nx,:,jz)**2 )))
        Nu_t(1:nx,:,jz) = S(1:nx,:)*Cs_opt2(1:nx,:,jz)*l(jz)**2
    end do

else

    ! define nu_coefs here since it does not change case-by-case for DNS
    ! if sgs, nu_coefs defined case-by-case since it is dependent on Nu_t
    if (.not. sgs) then
        nu_coef(1:nx,:) = 2.0_rprec*nu
        nu_coef2(1:nx,:) = 2.0_rprec*nu
    endif

end if !! if (sgs)

! Calculate txx, txy, tyy, tzz for bottom level: jz=1 node (coord==0 only)
if (coord == 0) then
    if (sgs) then
        select case (lbc_mom)
            ! Stress free
            case (0)
                nu_coef(1:nx,:) = 2.0_rprec*0.5_rprec*(Nu_t(1:nx,:,1) + Nu_t(1:nx,:,2)) + nu

            ! Wall
            case (1:)
                nu_coef(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,1)+nu) !! Nu_t on uvp-node(1) here

        end select
    endif
    txx(1:nx,:,1) = -nu_coef(1:nx,:)*dudx(1:nx,:,1) !! uvp-node(1)
    txy(1:nx,:,1) = -nu_coef(1:nx,:)*(0.5_rprec*(dudy(1:nx,:,1)+dvdx(1:nx,:,1))) !! uvp-node(1)
    tyy(1:nx,:,1) = -nu_coef(1:nx,:)*dvdy(1:nx,:,1) !! uvp-node(1)
    tzz(1:nx,:,1) = -nu_coef(1:nx,:)*dwdz(1:nx,:,1) !! uvp-node(1)

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if

! Calculate txx, txy, tyy, tzz for bottom level: jz=nz node (coord==nproc-1)
if (coord == nproc-1) then
    if (sgs) then
        select case (lbc_mom)
            ! Stress free
            case (0)
                nu_coef(1:nx,:) = 2.0_rprec*(0.5_rprec*(Nu_t(1:nx,:,nz-1) + Nu_t(1:nx,:,nz)) + nu) !! uvp-node(nz-1)
                nu_coef2(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,nz-1) + nu) !! w-node(nz-1)

            ! Wall
            case (1:)
                nu_coef(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,nz) + nu) !! uvp-node
                nu_coef2(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,nz-1) + nu) !! w-node

        end select
    endif
    txx(1:nx,:,nz-1) = -nu_coef(1:nx,:)*dudx(1:nx,:,nz-1) !! uvp-node(nz-1)
    txy(1:nx,:,nz-1) = -nu_coef(1:nx,:)*(0.5_rprec*(dudy(1:nx,:,nz-1)+dvdx(1:nx,:,nz-1))) !! uvp-node(nz-1)
    tyy(1:nx,:,nz-1) = -nu_coef(1:nx,:)*dvdy(1:nx,:,nz-1) !! uvp-node(nz-1)
    tzz(1:nx,:,nz-1) = -nu_coef(1:nx,:)*dwdz(1:nx,:,nz-1) !! uvp-node(nz-1)
    ! for top wall, include w-grid stress since we touched nz-1
#ifdef PPCNDIFF
    txz(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,nz-1)+dwdx(1:nx,:,nz-1))) !! w-node(nz-1)
    txz_half1(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dwdx(1:nx,:,nz-1))) !! w-node(nz-1)
    txz_half2(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,nz-1))) !! w-node(nz-1)
#else
    txz(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,nz-1)+dwdx(1:nx,:,nz-1))) !! w-node(nz-1)
#endif
    tyz(1:nx,:,nz-1) = -nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,nz-1)+dwdy(1:nx,:,nz-1))) !! w-node(nz-1)

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
        nu_coef(1:nx,:) = 2.0_rprec*(0.5_rprec*(Nu_t(1:nx,:,jz) + Nu_t(1:nx,:,jz+1)) + nu) !! uvp-node(jz)
        nu_coef2(1:nx,:) = 2.0_rprec*(Nu_t(1:nx,:,jz) + nu) !! w-node(jz)
    endif
    txx(1:nx,:,jz)=-nu_coef(1:nx,:)*dudx(1:nx,:,jz) !! uvp-node(jz)
    txy(1:nx,:,jz)=-nu_coef(1:nx,:)*(0.5_rprec*(dudy(1:nx,:,jz)+dvdx(1:nx,:,jz))) !! uvp-node(jz)
    tyy(1:nx,:,jz)=-nu_coef(1:nx,:)*dvdy(1:nx,:,jz) !! uvp-node(jz)
    tzz(1:nx,:,jz)=-nu_coef(1:nx,:)*dwdz(1:nx,:,jz) !! uvp-node(jz)
#ifdef PPCNDIFF
    txz(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,jz)+dwdx(1:nx,:,jz))) !! w-node(jz)
    txz_half1(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dwdx(1:nx,:,jz))) !! w-node(jz)
    txz_half2(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,jz))) !! w-node(jz)
#else
    txz(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dudz(1:nx,:,jz)+dwdx(1:nx,:,jz))) !! w-node(jz)
#endif
    tyz(1:nx,:,jz)=-nu_coef2(1:nx,:)*(0.5_rprec*(dvdz(1:nx,:,jz)+dwdy(1:nx,:,jz))) !! w-node(jz)
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
#else
call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
#endif
call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
#ifdef PPSAFETYMODE
! Set bogus values (easier to catch if there's an error)
txx(:, :, 0) = BOGUS
txy(:, :, 0) = BOGUS
#ifdef PPCNDIFF
txz(:, :, 0) = BOGUS
txz_half1(:, :, 0) = BOGUS
txz_half2(:, :, 0) = BOGUS
#else
txz(:, :, 0) = BOGUS
#endif
tyy(:, :, 0) = BOGUS
tyz(:, :, 0) = BOGUS
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
            !   does that mean they (Sij) should all be zero?
            ! these values are stored on w-nodes
            S11(1:nx,:,1) = dudx(1:nx,:,1) ! was 0.5_rprec*(dudx(1:nx,:,1) + dudx(1:nx,:,1))
            S12(1:nx,:,1) = 0.5_rprec*(dudy(1:nx,:,1)+dvdx(1:nx,:,1))
            S13(1:nx,:,1) = 0.5_rprec*(dudz(1:nx,:,1)+dwdx(1:nx,:,1))
            S22(1:nx,:,1) = dvdy(1:nx,:,1)
            S23(1:nx,:,1) = 0.5_rprec*(dvdy(1:nx,:,1)+dwdy(1:nx,:,1))
            S33(1:nx,:,1) = 0.5_rprec*(dwdz(1:nx,:,1) + 0.0_rprec)

        ! Wall
        ! recall dudz and dvdz are stored on uvp-nodes for first level only,
        !   'wall' only
        ! recall dwdx and dwdy are stored on w-nodes (always)
        ! All Sij near wall are put on uv1 node
        case (1:)
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
            !   does that mean they (Sij) should all be zero?
            ! these values are stored on w-nodes

            ! dudx(1:nx,:,nz-1) on uvp_node(nz-1)
            S11(1:nx,:,nz) = dudx(1:nx,:,nz-1)                               ! w_node(nz)
            ! dudy(1:nx,:,nz-1) on uvp_node(nz-1)
            ! dvdx(1:nx,:,nz-1) on uvp_node(nz-1)
            S12(1:nx,:,nz) = 0.5_rprec*(dudy(1:nx,:,nz-1)+dvdx(1:nx,:,nz-1)) ! w_node(nz)
            ! dudz(1:nx,:,nz) comes from wallstress() i.e. zero, on w_node(nz)
            ! dwdx(1:nx,:,nz) on uvp_node(nz-1)
            S13(1:nx,:,nz) = 0.5_rprec*(dudz(1:nx,:,nz)+dwdx(1:nx,:,nz))     ! w_node(nz)
            ! dvdy(1:nx,:,nz-1) on uvp_node(nz-1)
            S22(1:nx,:,nz) = dvdy(1:nx,:,nz-1)                               ! w_node(nz)
            ! dvdz(1:nx,:,nz) comes from wallstress() i.e. zero, on w_node(nz)
            ! dwdy(1:nx,:,nz) on uvp_node(nz-1)
            S23(1:nx,:,nz) = 0.5_rprec*(dvdz(1:nx,:,nz)+dwdy(1:nx,:,nz))     ! w_node(nz)
            ! dwdz(1:nx,:,nz-1) + 0._rprec ! w_node(nz)
            S33(1:nx,:,nz) = 0.5_rprec*(dwdz(1:nx,:,nz-1) + 0.0_rprec)       ! w_node(nz)           

        ! Wall
        ! recall dudz and dvdz are stored on uvp-nodes for first level only,
        !   'wall' only
        ! recall dwdx and dwdy are stored on w-nodes (always)
        case (1:)
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
