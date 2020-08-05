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
subroutine dyn_vreman(dyn_cvre)
!*******************************************************************************
!
! This subroutine calculates the eddy viscosity predicted by
! Vreman dynamically according to the Germano identity.
!
! This subroutine only outputs the dynamic constant (dyn_cvre)
! which is then used in the vreman subroutine.
!
use types, only : rprec
use param
use sim_param, only : u, v, w
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
#ifdef PPMAPPING
use sim_param, only : jaco_w, jaco_uv
#endif
use test_filtermodule
use sgs_param, only : S11, S12, S13, S22, S23, S33,           &
    u_bar, v_bar, w_bar, L11, L12, L13, L22, L23, L33,        &
    M11, M12, M13, M22, M23, M33,                             &
    S11_bar, S12_bar, S13_bar, S22_bar, S23_bar, S33_bar

implicit none

real(rprec), intent(out) :: dyn_cvre
real(rprec), dimension(ld,ny) :: dudx_p, dudy_p, dudz_p,         &
    dvdx_p, dvdy_p, dvdz_p, dwdx_p, dwdy_p, dwdz_p,              &
    dudx_bar, dudy_bar, dudz_bar, dvdx_bar, dvdy_bar, dvdz_bar,  &
    dwdx_bar, dwdy_bar, dwdz_bar
real(rprec), dimension(ld,ny) :: alpha2_bar, beta_bar, alpha2, beta, Mtemp
real(rprec), dimension(ld,ny) :: b11, b12, b13, b22, b23, b33,   &
    b11_bar, b12_bar, b13_bar, b22_bar, b23_bar, b33_bar
real(rprec) :: dz_p, eps, ctop, cbot, ctop_global, cbot_global
real(rprec), dimension(nz) :: cbot_z, ctop_z
integer :: jz

! Set constants
eps = 1e-12 !! some small number in case denominator is zero
dz_p = dz !! to be overwritten in mapping to stretched grid

#ifdef PPMPI
! dwdz calculated for 0:nz-1 (on w-nodes) except bottom process
! (only 1:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
#endif

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

#ifdef PPMAPPING
                dz_p = jaco_w(1)*dz
#endif

                u_bar(:,:) = 0.5_rprec*(u(:,:,1) + u(:,:,0))
                v_bar(:,:) = 0.5_rprec*(v(:,:,1) + v(:,:,0))
                w_bar(:,:) = w(:,:,1)
                L11(:,:) = u_bar(:,:)*u_bar(:,:)
                L12(:,:) = u_bar(:,:)*v_bar(:,:)
                L13(:,:) = u_bar(:,:)*w_bar(:,:)
                L22(:,:) = v_bar(:,:)*v_bar(:,:)
                L23(:,:) = v_bar(:,:)*w_bar(:,:)
                L33(:,:) = w_bar(:,:)*w_bar(:,:)

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

#ifdef PPMAPPING
                dz_p = jaco_uv(1)*dz
#endif

                L11(:,:) = u(:,:,1)*u(:,:,1)
                L12(:,:) = u(:,:,1)*v(:,:,1)
                L13(:,:) = u(:,:,1)*0.25_rprec*w(:,:,2) ! parabolic interp.
                L22(:,:) = v(:,:,1)*v(:,:,1)
                L23(:,:) = v(:,:,1)*0.25_rprec*w(:,:,2) ! parabolic interp.
                L33(:,:) = (0.25_rprec*w(:,:,2))**2     ! parabolic interp.
                u_bar(:,:) = u(:,:,1)
                v_bar(:,:) = v(:,:,1)
                w_bar(:,:) = 0.25_rprec*w(:,:,2) ! parabolic interp.

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
                dz_p = jaco_uv(1)*dz
#endif

                L11(:,:) = u(:,:,1)*u(:,:,1)
                L12(:,:) = u(:,:,1)*v(:,:,1)
                L13(:,:) = u(:,:,1)*0.25_rprec*w(:,:,2) ! parabolic interp.
                L22(:,:) = v(:,:,1)*v(:,:,1)
                L23(:,:) = v(:,:,1)*0.25_rprec*w(:,:,2) ! parabolic interp.
                L33(:,:) = (0.25_rprec*w(:,:,2))**2     ! parabolic interp.
                u_bar(:,:) = u(:,:,1)
                v_bar(:,:) = v(:,:,1)
                w_bar(:,:) = 0.25_rprec*w(:,:,2) ! parabolic interp.

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

#ifdef PPMAPPING
                dz_p = jaco_w(nz)*dz
#endif

                u_bar(:,:) = 0.5_rprec*(u(:,:,nz) + u(:,:,nz-1))
                v_bar(:,:) = 0.5_rprec*(v(:,:,nz) + v(:,:,nz-1))
                w_bar(:,:) = w(:,:,nz)
                L11(:,:) = u_bar(:,:)*u_bar(:,:)
                L12(:,:) = u_bar(:,:)*v_bar(:,:)
                L13(:,:) = u_bar(:,:)*w_bar(:,:)
                L22(:,:) = v_bar(:,:)*v_bar(:,:)
                L23(:,:) = v_bar(:,:)*w_bar(:,:)
                L33(:,:) = w_bar(:,:)*w_bar(:,:)

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

#ifdef PPMAPPING
                dz_p = jaco_uv(nz-1)*dz
#endif

                L11(:,:) = u(:,:,nz)*u(:,:,nz)
                L12(:,:) = u(:,:,nz)*v(:,:,nz)
                L13(:,:) = u(:,:,nz)*0.25_rprec*w(:,:,nz-1) ! parabolic interp.
                L22(:,:) = v(:,:,nz)*v(:,:,nz)
                L23(:,:) = v(:,:,nz)*0.25_rprec*w(:,:,nz-1) ! parabolic interp.
                L33(:,:) = (0.25_rprec*w(:,:,nz-1))**2     ! parabolic interp.
                u_bar(:,:) = u(:,:,nz)
                v_bar(:,:) = v(:,:,nz)
                w_bar(:,:) = 0.25_rprec*w(:,:,nz-1) ! parabolic interp.

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
                dz_p = jaco_uv(nz-1)*dz
#endif

                L11(:,:) = u(:,:,nz)*u(:,:,nz)
                L12(:,:) = u(:,:,nz)*v(:,:,nz)
                L13(:,:) = u(:,:,nz)*0.25_rprec*w(:,:,nz-1) ! parabolic interp.
                L22(:,:) = v(:,:,nz)*v(:,:,nz)
                L23(:,:) = v(:,:,nz)*0.25_rprec*w(:,:,nz-1) ! parabolic interp.
                L33(:,:) = (0.25_rprec*w(:,:,nz-1))**2     ! parabolic interp.
                u_bar(:,:) = u(:,:,nz)
                v_bar(:,:) = v(:,:,nz)
                w_bar(:,:) = 0.25_rprec*w(:,:,nz-1) ! parabolic interp.

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

#ifdef PPMAPPING
        dz_p = jaco_w(jz)*dz
#endif

        u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        w_bar(:,:) = w(:,:,jz)
        L11(:,:) = u_bar(:,:)*u_bar(:,:)
        L12(:,:) = u_bar(:,:)*v_bar(:,:)
        L13(:,:) = u_bar(:,:)*w_bar(:,:)
        L22(:,:) = v_bar(:,:)*v_bar(:,:)
        L23(:,:) = v_bar(:,:)*w_bar(:,:)
        L33(:,:) = w_bar(:,:)*w_bar(:,:)

    end if

    ! in-place filtering
    call test_filter ( u_bar )
    call test_filter ( v_bar )
    call test_filter ( w_bar )

    call test_filter ( L11 )
    L11 = L11 - u_bar*u_bar
    call test_filter ( L12 )
    L12 = L12 - u_bar*v_bar
    call test_filter ( L13 )
    L13 = L13 - u_bar*w_bar
    call test_filter ( L22 )
    L22 = L22 - v_bar*v_bar
    call test_filter ( L23 )
    L23 = L23 - v_bar*w_bar
    call test_filter ( L33 )
    L33 = L33 - w_bar*w_bar

    ! S_ij already on w-nodes
    S11_bar(:,:) = S11(:,:,jz)
    S12_bar(:,:) = S12(:,:,jz)
    S13_bar(:,:) = S13(:,:,jz)
    S22_bar(:,:) = S22(:,:,jz)
    S23_bar(:,:) = S23(:,:,jz)
    S33_bar(:,:) = S33(:,:,jz)

    call test_filter ( S11_bar )
    call test_filter ( S12_bar )
    call test_filter ( S13_bar )
    call test_filter ( S22_bar )
    call test_filter ( S23_bar )
    call test_filter ( S33_bar )

    dudx_bar(:,:) = dudx_p(:,:)
    dudy_bar(:,:) = dudy_p(:,:)
    dudz_bar(:,:) = dudz_p(:,:)
    dvdx_bar(:,:) = dvdx_p(:,:)
    dvdy_bar(:,:) = dvdy_p(:,:)
    dvdz_bar(:,:) = dvdz_p(:,:)
    dwdx_bar(:,:) = dwdx_p(:,:)
    dwdy_bar(:,:) = dwdy_p(:,:)
    dwdz_bar(:,:) = dwdz_p(:,:)

    call test_filter ( dudx_bar )
    call test_filter ( dudy_bar )
    call test_filter ( dudz_bar )
    call test_filter ( dvdx_bar )
    call test_filter ( dvdy_bar )
    call test_filter ( dvdz_bar )
    call test_filter ( dwdx_bar )
    call test_filter ( dwdy_bar )
    call test_filter ( dwdz_bar )

    ! Compute inner velocity gradient product
    alpha2(:,:) = dudx_p**2 + dudy_p**2 + dudz_p**2 +       &
        dvdx_p**2 + dvdy_p**2 + dvdz_p**2 +                 &
        dwdx_p**2 + dwdy_p**2 + dwdz_p**2
    alpha2_bar = dudx_bar**2 + dudy_bar**2 + dudz_bar**2 +  &
        dvdx_bar**2 + dvdy_bar**2 + dvdz_bar**2 +           &
        dwdx_bar**2 + dwdy_bar**2 + dwdz_bar**2

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
    b11_bar = 4._rprec*(dx**2)*(dudx_bar**2) +                            &
        4._rprec*(dy**2)*(dudy_bar**2) + (dz_p**2)*(dudz_bar**2)
    b12_bar = 4._rprec*(dx**2)*dudx_bar*dvdx_bar +                        &
        4._rprec*(dy**2)*dudy_bar*dvdy_bar + (dz_p**2)*dudz_bar*dvdz_bar
    b13_bar = 4._rprec*(dx**2)*dudx_bar*dwdx_bar +                        &
        4._rprec*(dy**2)*dudy_bar*dwdy_bar + (dz_p**2)*dudz_bar*dwdz_bar
    b22_bar = 4._rprec*(dx**2)*(dvdx_bar**2) +                            &
        4._rprec*(dy**2)*(dvdy_bar**2) + (dz_p**2)*(dvdz_bar**2)
    b23_bar = 4._rprec*(dx**2)*dvdx_bar*dwdx_bar +                        &
        4._rprec*(dy**2)*dvdy_bar*dwdy_bar + (dz_p**2)*dvdz_bar*dwdz_bar
    b33_bar = 4._rprec*(dx**2)*(dwdx_bar**2) +                            &
        4._rprec*(dy**2)*(dwdy_bar**2) + (dz_p**2)*(dwdz_bar**2)

    ! Compute beta product
    beta(:,:) = b11*b22 - (b12**2) +                                     &
        b11*b33 - (b13**2) + b22*b33 - (b23**2)
    beta_bar = b11_bar*b22_bar - (b12_bar**2) +                          &
        b11_bar*b33_bar - (b13_bar**2) + b22_bar*b33_bar - (b23_bar**2)

    ! Use M_ij terms as temp variable to perform test_filtered part
    Mtemp(:,:) = sqrt( abs( beta(:,:) / (alpha2(:,:) + eps) ) )
    M11 = Mtemp(:,:)*S11(:,:,jz)
    M12 = Mtemp(:,:)*S12(:,:,jz)
    M13 = Mtemp(:,:)*S13(:,:,jz)
    M22 = Mtemp(:,:)*S22(:,:,jz)
    M23 = Mtemp(:,:)*S23(:,:,jz)
    M33 = Mtemp(:,:)*S33(:,:,jz)

    call test_filter ( M11 )
    call test_filter ( M12 )
    call test_filter ( M13 )
    call test_filter ( M22 )
    call test_filter ( M23 )
    call test_filter ( M33 )

    ! Compute actual M_ij terms
    Mtemp(:,:) = sqrt( abs( beta_bar(:,:) / (alpha2_bar(:,:) + eps) ) )
    M11(:,:) = M11(:,:) - Mtemp(:,:)*S11_bar(:,:)
    M12(:,:) = M12(:,:) - Mtemp(:,:)*S12_bar(:,:)
    M13(:,:) = M13(:,:) - Mtemp(:,:)*S13_bar(:,:)
    M22(:,:) = M22(:,:) - Mtemp(:,:)*S22_bar(:,:)
    M23(:,:) = M23(:,:) - Mtemp(:,:)*S23_bar(:,:)
    M33(:,:) = M33(:,:) - Mtemp(:,:)*S33_bar(:,:)

    ! Perform planar average for points on this coord
    ctop_z(jz) = sum(L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23))
    cbot_z(jz) = sum(M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2))

end do

! Sum across z for points on this coord
ctop = sum( ctop_z )
cbot = sum( cbot_z )

! Sum across coord
#ifdef PPMPI
call mpi_allreduce(ctop,ctop_global,1,MPI_RPREC,MPI_SUM,MPI_COMM_WORLD,ierr)
call mpi_allreduce(cbot,cbot_global,1,MPI_RPREC,MPI_SUM,MPI_COMM_WORLD,ierr)
ctop = ctop_global
cbot = cbot_global
#endif

! Compute model constant
dyn_cvre = 0.5_rprec*(ctop/cbot)

! Clip
dyn_cvre = max(0._rprec,dyn_cvre)

#ifdef PPOUTPUT_SGS
if (coord == 0) write(*,*) 'Dynamic Vreman Model Constant: ', dyn_cvre
#endif

end subroutine dyn_vreman
