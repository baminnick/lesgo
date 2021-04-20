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
module tlwmles
!*****************************************************************************
!
! This module contains all of the subroutines associated with the two-layer
! wall-model equations for LES
!

use types, only : rprec
use param, only : path

implicit none

save
private

public :: tlwm_init, ic_tlwm, tlwm_checkpoint, tlwm_wallstress, tavg_tlwm_finalize

! Main simulation variables for the TLWMLES
real(rprec), public, dimension(:,:,:), allocatable :: ur, vr, wr
real(rprec), dimension(:,:,:), allocatable :: dudxr, dudyr, dudzr,      &
    dvdxr, dvdyr, dvdzr, dwdxr, dwdyr, dwdzr,                           &
    txxr, txyr, tyyr, txzr, tyzr,                                       &
    txzr_half1, txzr_half2, tyzr_half1, tyzr_half2,                     &
    div_txr, div_tyr,                                                   &
    rhs_xr, rhs_yr, rhs_xr_f, rhs_yr_f
real(rprec), dimension(:,:,:), allocatable :: nu_uvr, nu_wr
real(rprec), dimension(:,:,:), allocatable :: ur_big, vr_big, wr_big,   &
    dudxr_big, dudyr_big, dudzr_big,                                    & 
    dvdxr_big, dvdyr_big, dvdzr_big,                                    & 
    temp_big
real(rprec), dimension(:,:), allocatable :: dpdxr, dpdyr

! Grid variables
real(rprec), dimension(:), allocatable :: zuvr, zwr, jaco_uvr, jaco_wr

#ifdef PPTLWM_LVLSET
! Level-set parameters
real(rprec), dimension(:,:,:), allocatable :: phi_uvr, phi_wr
real(rprec), dimension(:,:), allocatable :: dudzr_wall, dvdzr_wall,     &
    txzr_wall, tyzr_wall, phi_wall
integer, dimension(:,:), allocatable :: coord_wall, k_wall
#endif

! Whether to initialize tlwm field
logical, public :: init_tlwm = .true.
! Name of file for restarting
character(64) :: fname

! Initial noise
real(rprec), public :: initial_noise_wm = 0._rprec

contains 

!*****************************************************************************
subroutine tlwm_init
!*****************************************************************************
!
! This subroutine initializes the variables for the tlwmles module
!
use param

! Allocate simulations variables
allocate ( ur(nxr+2, nyr, lbz:nzr) ); ur = 0._rprec
allocate ( vr(nxr+2, nyr, lbz:nzr) ); vr = 0._rprec
allocate ( wr(nxr+2, nyr, lbz:nzr) ); wr = 0._rprec
allocate ( dudxr(nxr+2, nyr, lbz:nzr) ); dudxr = 0._rprec
allocate ( dudyr(nxr+2, nyr, lbz:nzr) ); dudyr = 0._rprec
allocate ( dudzr(nxr+2, nyr, lbz:nzr) ); dudzr = 0._rprec
allocate ( dvdxr(nxr+2, nyr, lbz:nzr) ); dvdxr = 0._rprec
allocate ( dvdyr(nxr+2, nyr, lbz:nzr) ); dvdyr = 0._rprec
allocate ( dvdzr(nxr+2, nyr, lbz:nzr) ); dvdzr = 0._rprec
allocate ( dwdxr(nxr+2, nyr, lbz:nzr) ); dwdxr = 0._rprec
allocate ( dwdyr(nxr+2, nyr, lbz:nzr) ); dwdyr = 0._rprec
allocate ( dwdzr(nxr+2, nyr, lbz:nzr) ); dwdzr = 0._rprec
allocate ( txxr(nxr+2, nyr, lbz:nzr) ); txxr = 0._rprec
allocate ( txyr(nxr+2, nyr, lbz:nzr) ); txyr = 0._rprec
allocate ( tyyr(nxr+2, nyr, lbz:nzr) ); tyyr = 0._rprec
allocate ( txzr(nxr+2, nyr, lbz:nzr) ); txzr = 0._rprec
allocate ( txzr_half1(nxr+2, nyr, lbz:nzr) ); txzr_half1 = 0._rprec
allocate ( txzr_half2(nxr+2, nyr, lbz:nzr) ); txzr_half2 = 0._rprec
allocate ( tyzr(nxr+2, nyr, lbz:nzr) ); tyzr = 0._rprec
allocate ( tyzr_half1(nxr+2, nyr, lbz:nzr) ); tyzr_half1 = 0._rprec
allocate ( tyzr_half2(nxr+2, nyr, lbz:nzr) ); tyzr_half2 = 0._rprec
allocate ( div_txr(nxr+2, nyr, lbz:nzr) ); div_txr = 0._rprec
allocate ( div_tyr(nxr+2, nyr, lbz:nzr) ); div_tyr = 0._rprec
allocate ( nu_wr(nxr,nyr,lbz:nzr) ); nu_wr = 0._rprec
allocate ( nu_uvr(nxr,nyr,lbz:nzr) ); nu_uvr = 0._rprec
allocate ( rhs_xr(nxr+2, nyr, lbz:nzr) ); rhs_xr = 0.0_rprec
allocate ( rhs_yr(nxr+2, nyr, lbz:nzr) ); rhs_yr = 0.0_rprec
allocate ( rhs_xr_f(nxr+2, nyr, lbz:nzr) ); rhs_xr_f = 0.0_rprec
allocate ( rhs_yr_f(nxr+2, nyr, lbz:nzr) ); rhs_yr_f = 0.0_rprec
allocate ( dpdxr(nxr, nyr) ); dpdxr = 0._rprec
allocate ( dpdyr(nxr, nyr) ); dpdyr = 0._rprec

! Big variables
! Remember nx2 = 3*nx/2, lh_big = nx2/2 + 1, ld_big= 2*lh_big
! Therefore big variables using nxr are size 3*nxr/2 + 2
allocate ( ur_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); ur_big = 0._rprec
allocate ( vr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); vr_big = 0._rprec
allocate ( wr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); wr_big = 0._rprec
allocate ( dudxr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dudxr_big = 0._rprec
allocate ( dudyr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dudyr_big = 0._rprec
allocate ( dudzr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dudzr_big = 0._rprec
allocate ( dvdxr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dvdxr_big = 0._rprec
allocate ( dvdyr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dvdyr_big = 0._rprec
allocate ( dvdzr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dvdzr_big = 0._rprec
allocate ( temp_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); temp_big = 0._rprec

! Create TLWM inner layer wall-normal grid
if (str_on) then
    call tlwm_grid_stretched()
else
    call tlwm_grid_uniform()
endif

#ifdef PPTLWM_LVLSET
! Simulation variables for level-set
allocate ( dudzr_wall(nxr+2,nyr) ); dudzr_wall = 0._rprec
allocate ( dvdzr_wall(nxr+2,nyr) ); dvdzr_wall = 0._rprec
allocate ( txzr_wall(nxr+2,nyr) ); txzr_wall = 0._rprec
allocate ( tyzr_wall(nxr+2,nyr) ); tyzr_wall = 0._rprec

! Load topography and allocate parameters for level-set
call load_tlwm_topography()
#endif

end subroutine tlwm_init

!*****************************************************************************
subroutine tlwm_wallstress
!*****************************************************************************
!
! This subroutine acts as the main subroutine for all two-layer wall-models.
! It calls the wall-model that is to be used to give the wall stress, this
! subroutine also calls functions that output data and statistics.
!
use param, only : lbc_mom, coord
use param, only : checkpoint_data, checkpoint_nskip, jt_total
use messages, only : error
#ifdef PPOUTPUT_WMLES
use param, only : dt, tavg_nstart, tavg_nend, tavg_nskip, tavg_calc
use param, only : domain_calc, domain_nstart, domain_nskip, domain_nend
use stat_defs, only : tavg_initialized, tavg_tlwm_dt
#endif
implicit none
character(*), parameter :: sub_name = 'tlwm_wallstress'

! Calculate wall stress
select case (lbc_mom)
    ! Standard equilibrium wall-model
    case (5)
        call tlwm_eq_ubc() !! Compute u and v upper BC for tlwm
        call tlwm_eq_solve() !! Solve for wall-model velocities

    ! Standard non-equilibrium wall-model
    case (6)
        call tlwm_noneq_ubc() !! compute u, v, w, dudz, dvdz, dwdz upper BC
        call tlwm_noneq_solve() !! Solve for wall-model velocities

    ! Otherwise, invalid
    case default
        call error (sub_name, 'invalid lbc_mom')
end select

! Compute and filter wall-stress answer for LES using TLWM sol
#ifdef PPTLWM_LVLSET
call tlwm_lvlset_wallstress()
call lvlset_les_lbc()
#else
if (coord == 0) then
    call inner_layer_wallstress()
    call les_lbc()
endif
#endif

#ifdef PPOUTPUT_WMLES
! OUTPUT LOOP (for TLWMLES)
! Instantaneous Domain Velocities
if (domain_calc) then
    if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and.  &
        (mod(jt_total-domain_nstart,domain_nskip)==0) ) then
        call tlwm_inst_write()
    endif
endif

! Time-averaging statistics
! Determine if time summations are to be calculated
if (tavg_calc) then
    ! Are we between the start and stop timesteps?
    if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then
        ! Every timestep (between nstart and nend), add to tavg_tlwm_dt
        tavg_tlwm_dt = tavg_tlwm_dt + dt

        ! Are we at the beginning or a multiple of nstart?
        if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then
            ! Check if we have initialized tavg

            if (.not.tavg_initialized) then
                call tavg_tlwm_init()
            else
                call tavg_tlwm_compute ()
            end if
        end if
    end if
end if
#endif

end subroutine tlwm_wallstress

!*****************************************************************************
subroutine tlwm_grid_uniform
!*****************************************************************************
!
! This subroutine creates a inner-layer wall-normal uniform grid for the 
! uv and w variables for the TLWM.
!
use param, only : nzr, dzr, lbz, coord
implicit none
integer :: jz

allocate(zuvr(lbz:nzr), zwr(lbz:nzr), jaco_uvr(lbz:nzr), jaco_wr(lbz:nzr))


do jz = lbz, nzr
    zuvr(jz) = (coord*(nzr-1) + jz - 0.5_rprec) * dzr
enddo

zwr = zuvr - dzr/2._rprec

! Define Jacobian values as 1, no stretching
jaco_uvr(:) = 1._rprec
jaco_wr(:) = 1._rprec

end subroutine tlwm_grid_uniform

!*****************************************************************************
subroutine tlwm_grid_stretched
!*****************************************************************************
!
! This subroutine creates a inner-layer wall-normal stretched grid for the 
! uv and w variables for the TLWM. This is done in a similar fashion to
! what is done in load_jacobian.f90, where a uniform grid is mapped to the
! stretched grid.
!
use param, only : nzr, nzr_tot, dzr, lbz, coord, str_r, L_zr
implicit none
integer :: jz
real(rprec), dimension(nzr_tot) :: z_w, z_uv, f1, f2, f3, f4 !! temp grid variables

allocate(zuvr(lbz:nzr), zwr(lbz:nzr), jaco_uvr(lbz:nzr), jaco_wr(lbz:nzr))

! Create unstretched grid to be mapped to new grid
! NOTE: using nzr_tot here
do jz = 1, nzr_tot
    z_uv(jz) = (jz-0.5_rprec)*dzr !! z-locations on uv-grid
enddo
z_w(:) = z_uv(:) - dzr/2._rprec !! z-locations on w-grid

! Mapped to stretched grid
! For the uv-grid
f3(:) = L_zr*(1.0_rprec+(tanh(str_r*(z_uv(:)/L_zr-1.0_rprec)) &
    /tanh(str_r)))
! For the w-grid
f4(:) = L_zr*(1.0_rprec+(tanh(str_r*(z_w(:)/L_zr-1.0_rprec))  &
    /tanh(str_r)))

! Compute Jacobian values for both w- and uv-grids
! Using analytical derivative expression
f1(:) = L_zr*(str_r/L_zr)*                                 &
    (1-(tanh(str_r*(z_w(:)/L_zr-1)))**2)/tanh(str_r)
f2(:) = L_zr*(str_r/L_zr)*                                 &
    (1-(tanh(str_r*(z_uv(:)/L_zr-1)))**2)/tanh(str_r)

! Store variables into what TLWM will use
do jz = 1, nzr
    jaco_wr(jz) = f1(coord*(nzr-1)+jz)
    jaco_uvr(jz) = f2(coord*(nzr-1)+jz)
    zuvr(jz) = f3(coord*(nzr-1)+jz)
    zwr(jz) = f4(coord*(nzr-1)+jz)
enddo

if (coord == 0) then
    jaco_wr(lbz) = jaco_wr(1)
    jaco_uvr(lbz) = jaco_uvr(1)
    zuvr(lbz) = -zuvr(1)
    zwr(lbz) = -zwr(lbz)
    write(*,*) '--> Inner-layer grid stretched using mapping function'
else
    jaco_wr(lbz) = f1((coord-1)*(nzr-1)+nzr-1)
    jaco_uvr(lbz) = f2((coord-1)*(nzr-1)+nzr-1)
    zuvr(lbz) = f3((coord-1)*(nzr-1)+nzr-1)
    zwr(lbz) = f4((coord-1)*(nzr-1)+nzr-1)
endif

end subroutine tlwm_grid_stretched

!*****************************************************************************
subroutine ic_tlwm
!*****************************************************************************
! Set initial profile for wall-model equation, determine if there is a file
! to read in or a new profile needs to be generated
!
use param, only : coord, lbc_mom, BOGUS, lbz, nz
use string_util

fname = path // 'tlwm.out'
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif
inquire (file=fname, exist=init_tlwm)

! Initialize any additional TLWM simulation parameters here

! Initialize TLWM field
if (init_tlwm) then
    if (coord == 0) write(*,*) "--> Reading initial TLWM field from file"
    call ic_tlwm_file
else
    if (coord == 0) write(*,*) "--> Creating initial TLWM field from prescribed turbulent profile"
    call ic_tlwm_vel
endif

!--set 0-level velocities to BOGUS
if (coord == 0) then
    ur(:, :, lbz) = BOGUS
    vr(:, :, lbz) = BOGUS
end if

end subroutine ic_tlwm

!*****************************************************************************
subroutine ic_tlwm_file
!*****************************************************************************
! Read initial profile for TLWM equation from a file
use param, only : nzr, read_endian

open(12, file=fname, form='unformatted', convert=read_endian)
read(12) ur(:,:,1:nzr), vr(:,:,1:nzr), rhs_xr(:,:,1:nzr), rhs_yr(:,:,1:nzr)
close(12)

end subroutine ic_tlwm_file

!*****************************************************************************
subroutine ic_tlwm_vel
!*****************************************************************************
!
! This subroutine initializes the initial tlwm profile 
!
use param
implicit none

! Enforce interpolated LES velocity field as upper BC for TLWM
call tlwm_eq_ubc()

! Solve equilibrium wall-model ODE to initialize the TLWM field
call tlwm_eq_solve()

#ifdef PPTLWM_LVLSET
! Compute and filter wall-stress for LES
call tlwm_lvlset_wallstress()
call lvlset_les_lbc()
#else
! Filter wall-stress answer for LES
if (coord == 0) then
    call inner_layer_wallstress()
    call les_lbc()
endif
#endif

end subroutine ic_tlwm_vel

!*****************************************************************************
subroutine tlwm_eq_ubc
!*****************************************************************************
!
! This subroutine creates the upper boundary condition for the equilibrium
! inner-layer wall-model. This requires grabbing the LES velocity at z = hwm
! on a lower coord, then moving the data to coord = nproc-1 and 
! interpolating the data to the TLWM grid.
! 
! The equilibrium wall-model only requires the u, v velocity from the LES.
! 
use param, only : nx, ny, jz_coord, jz_r, nzr
use param, only : coord, nproc, comm, ierr, MPI_RPREC
use param, only : nxr, nyr
use sim_param, only : u, v
use mpi
implicit none
real(rprec), dimension(nx,ny) :: dummy_in, u_les, v_les
real(rprec), dimension(nxr,nyr) :: dummy_out

if (coord == jz_coord) then
    u_les = u(1:nx,1:ny,jz_r)
    v_les = v(1:nx,1:ny,jz_r)
else
    u_les = 0.0_rprec
    v_les = 0.0_rprec
endif

! Gather data and interpolate on top coord
call mpi_allreduce(u_les, dummy_in, nx*ny, mpi_rprec,                      &
    MPI_SUM, comm, ierr)
if (coord == nproc-1) then
    call interp_les_to_tlwm(dummy_in,dummy_out)
    ur(1:nxr,1:nyr,nzr) = dummy_out
endif

call mpi_allreduce(v_les, dummy_in, nx*ny, mpi_rprec,                      &
    MPI_SUM, comm, ierr)
if (coord == nproc-1) then
    call interp_les_to_tlwm(dummy_in,dummy_out)
    vr(1:nxr,1:nyr,nzr) = dummy_out
endif

!if (coord == jz_coord) write(*,*) 'LES velocity: ', u(1,1,jz_r), v(1,1,jz_r)
!if (coord == nproc-1) write(*,*) 'Interpolated: ', ur(2,2,nzr), vr(2,2,nzr)

end subroutine tlwm_eq_ubc

!*****************************************************************************
subroutine tlwm_noneq_ubc
!*****************************************************************************
!
! This subroutine follows a similar procedure to tlwm_eq_ubc, however moves
! additional flow variables around because more information is needed for the
! non-equilibrium inner-layer wall-model.
! 
! u, v, w, dudz, dvdz, dwdz, dpdx, dpdy at z=hwm (a uvp-node) are gathered
! for the TLWM non-equilibrium upper-boundary condition. Necessary 
! interpolations are taken.
!
use param, only : ld, nx, ny, jz_coord, jz_r, nzr, ld_big, nx2, ny2
use param, only : coord, nproc, comm, ierr, MPI_RPREC
use param, only : nxr, nyr
use sim_param, only : u, v, w, dudx, dvdx, dwdx, dudy, dvdy, dwdy, dpdx, dpdy
use mpi
use fft
implicit none
real(rprec), dimension(nx,ny) :: dummy_in, u_les, v_les, dpdx_les, dpdy_les
real(rprec), dimension(ld,ny) :: ut, vt, wt, dudxt, dvdxt, dwdxt,            &
    dudyt, dvdyt, dwdyt, dEdx, dEdy
real(rprec), dimension(ld_big,ny2) :: ut_big, vt_big, wt_big,                &
    dudxt_big, dvdxt_big, dwdxt_big, dudyt_big, dvdyt_big, dwdyt_big,        &
    dEdx_big, dEdy_big
real(rprec), dimension(nxr,nyr) :: dummy_out
real(rprec) :: const

if (coord == jz_coord) then
    u_les = u(1:nx,1:ny,jz_r)
    v_les = v(1:nx,1:ny,jz_r)

    ! Compute real pressure gradients by extracting energy
    ! This requires taking products, so take to big domain
    ! Also interpolate w-grid variables: w, dwdx, dwdy
    const = 1._rprec/(nx*ny)

    ut = const*u(:,:,jz_r)
    vt = const*v(:,:,jz_r)
    wt = const*0.5_rprec*( w(:,:,jz_r+1) + w(:,:,jz_r) )
    dudxt = const*dudx(:,:,jz_r)
    dvdxt = const*dvdx(:,:,jz_r)
    dwdxt = const*0.5_rprec*( dwdx(:,:,jz_r+1) + dwdx(:,:,jz_r) )
    dudyt = const*dudy(:,:,jz_r)
    dvdyt = const*dvdy(:,:,jz_r)
    dwdyt = const*0.5_rprec*( dwdy(:,:,jz_r+1) + dwdy(:,:,jz_r) )

    call dfftw_execute_dft_r2c(forw, ut(:,:), ut(:,:))
    call dfftw_execute_dft_r2c(forw, vt(:,:), vt(:,:))
    call dfftw_execute_dft_r2c(forw, wt(:,:), wt(:,:))
    call dfftw_execute_dft_r2c(forw, dudxt(:,:), dudxt(:,:))
    call dfftw_execute_dft_r2c(forw, dvdxt(:,:), dvdxt(:,:))
    call dfftw_execute_dft_r2c(forw, dwdxt(:,:), dwdxt(:,:))
    call dfftw_execute_dft_r2c(forw, dudyt(:,:), dudyt(:,:))
    call dfftw_execute_dft_r2c(forw, dvdyt(:,:), dvdyt(:,:))
    call dfftw_execute_dft_r2c(forw, dwdyt(:,:), dwdyt(:,:))

    call padd( ut_big(:,:), ut(:,:) )
    call padd( vt_big(:,:), vt(:,:) )
    call padd( wt_big(:,:), wt(:,:) )
    call padd( dudxt_big(:,:), dudxt(:,:) )
    call padd( dvdxt_big(:,:), dvdxt(:,:) )
    call padd( dwdxt_big(:,:), dwdxt(:,:) )
    call padd( dudyt_big(:,:), dudyt(:,:) )
    call padd( dvdyt_big(:,:), dvdyt(:,:) )
    call padd( dwdyt_big(:,:), dwdyt(:,:) )

    call dfftw_execute_dft_c2r(back_big, ut_big(:,:), ut_big(:,:))
    call dfftw_execute_dft_c2r(back_big, vt_big(:,:), vt_big(:,:))
    call dfftw_execute_dft_c2r(back_big, wt_big(:,:), wt_big(:,:))
    call dfftw_execute_dft_c2r(back_big, dudxt_big(:,:), dudxt_big(:,:))
    call dfftw_execute_dft_c2r(back_big, dvdxt_big(:,:), dvdxt_big(:,:))
    call dfftw_execute_dft_c2r(back_big, dwdxt_big(:,:), dwdxt_big(:,:))
    call dfftw_execute_dft_c2r(back_big, dudyt_big(:,:), dudyt_big(:,:))
    call dfftw_execute_dft_c2r(back_big, dvdyt_big(:,:), dvdyt_big(:,:))
    call dfftw_execute_dft_c2r(back_big, dwdyt_big(:,:), dwdyt_big(:,:))

    ! Compute horizontal derivatives of energy    
    ! Redefinition of const
    const = 1._rprec/(nx2*ny2)
    dEdx_big = const*( ut_big(:,:)*dudxt_big(:,:) +                    &
        vt_big(:,:)*dvdxt_big(:,:) + wt_big(:,:)*dwdxt_big(:,:) )
    dEdy_big = const*( ut_big(:,:)*dudyt_big(:,:) +                    &
        vt_big(:,:)*dvdyt_big(:,:) + wt_big(:,:)*dwdyt_big(:,:) )

    ! Take back to small
    call dfftw_execute_dft_r2c(forw_big, dEdx_big(:,:), dEdx_big(:,:))
    call dfftw_execute_dft_r2c(forw_big, dEdy_big(:,:), dEdy_big(:,:))

    call unpadd(dEdx(:,:), dEdx_big(:,:))
    call unpadd(dEdy(:,:), dEdy_big(:,:))

    call dfftw_execute_dft_c2r(back, dEdx(:,:), dEdx(:,:))
    call dfftw_execute_dft_c2r(back, dEdy(:,:), dEdy(:,:))

    ! Finally compute real pressure gradients dpdx and dpdy
    dpdx_les(:,:) = dpdx(1:nx,:,jz_r) - dEdx(1:nx,:)
    dpdy_les(:,:) = dpdy(1:nx,:,jz_r) - dEdy(1:nx,:)

else
    u_les = 0.0_rprec
    v_les = 0.0_rprec
endif

! Gather u & v velocity data and interpolate on top coord
call mpi_allreduce(u_les, dummy_in, nx*ny, mpi_rprec, MPI_SUM, comm, ierr)
if (coord == nproc-1) then
    call interp_les_to_tlwm(dummy_in,dummy_out)
    ur(1:nxr,1:nyr,nzr) = dummy_out
endif

call mpi_allreduce(v_les, dummy_in, nx*ny, mpi_rprec, MPI_SUM, comm, ierr)
if (coord == nproc-1) then
    call interp_les_to_tlwm(dummy_in,dummy_out)
    vr(1:nxr,1:nyr,nzr) = dummy_out
endif

! Since dpdx and dpdy will be used on all processors, interpolate 
! on coord where it was computed, then send out to all processors
if (coord == jz_coord) then
    dummy_in = dpdx_les
    call interp_les_to_tlwm(dummy_in,dummy_out)
else
    dummy_out = 0._rprec
endif
call mpi_allreduce(dummy_out, dpdxr, nxr*nyr, mpi_rprec, MPI_SUM, comm, ierr)

if (coord == jz_coord) then
    dummy_in = dpdy_les
    call interp_les_to_tlwm(dummy_in,dummy_out)
else
    dummy_out = 0._rprec
endif
call mpi_allreduce(dummy_out, dpdyr, nxr*nyr, mpi_rprec, MPI_SUM, comm, ierr)

end subroutine tlwm_noneq_ubc

!*****************************************************************************
subroutine les_lbc
!*****************************************************************************
!
! This subroutine creates the lower boundary condition for the LES using the
! solution of the inner-layer wall-model. The derivatives and wall-stress
! at z=0 are assumed to have already been computed, then in this subroutine,
! these data are filtered onto the LES grid.
!
! This routine should only be accessed by coord=0
!
use param, only : nx, ny, ubot, nu_molec, z_i, u_star
use param, only : nxr, nyr, dzr
use sim_param, only : dudz, dvdz, txz, tyz
use mpi
use test_filtermodule, only : tlwm_filter
implicit none
integer :: jx, jy
real(rprec), dimension(nxr+2,nyr) :: dummy_in
real(rprec), dimension(nx,ny) :: dummy_out

! Filter data to LES
dummy_in = dudzr(:,:,1)
call tlwm_filter(dummy_in)
call interp_tlwm_to_les(dummy_in,dummy_out)
dudz(1:nx,1:ny,1) = dummy_out

dummy_in = dvdzr(:,:,1)
call tlwm_filter(dummy_in)
call interp_tlwm_to_les(dummy_in,dummy_out)
dvdz(1:nx,1:ny,1) = dummy_out

dummy_in = txzr(:,:,1)
call tlwm_filter(dummy_in)
call interp_tlwm_to_les(dummy_in,dummy_out)
txz(1:nx,1:ny,1) = dummy_out

dummy_in = tyzr(:,:,1)
call tlwm_filter(dummy_in)
call interp_tlwm_to_les(dummy_in,dummy_out)
tyz(1:nx,1:ny,1) = dummy_out

end subroutine les_lbc

!*****************************************************************************
subroutine tlwm_eq_solve
!*****************************************************************************
! 
! This subroutine solves a 2nd order ODE which assumes equilibrium and a 
! eddy viscosity model that satisfies the log-law and the viscous sublayer if
! the LES grid is sufficiently resolved.
! 
! Currently the model assumes smooth walls, finite Reynolds number, and takes 
! the inner layer starting at the ihwm'th LES grid point from the wall (on the 
! uv grid).
!  
! The ODE is of the form:
!     d/dz[(nu+nu_t)*dUdz] = 0,  U(z=0) = 0,  U(z=hwm) = Ules
! Centered differencing is used. As a result, the quantity (nu+nu_t) is
! evaluated between grid points which U is solved on. The TDMA algorithm is 
! used to invert the resulting matrix equation.
! 
use messages
use param, only : nxr, nyr, nzr, nu_molec, z_i, ubot, u_star, jt_total, dzr, L_zr
use param, only : coord, nproc, comm, ierr, MPI_RPREC
#ifdef PPSAFETYMODE
use param, only : BOGUS
#endif
use mpi
implicit none
integer :: jx, jy, jz, jz_min, jz_max
real(rprec), dimension(nxr, nyr) :: ustar, dummy
real(rprec) :: const1
real(rprec), dimension(nxr,nyr,0:nzr) :: a, b, c
real(rprec), dimension(nxr+2,nyr,0:nzr) :: Rx, usol, Ry, vsol

const1 = 1._rprec/dzr
! Commented out old code for uniform grid which used const1=1/(dzr**2)
! Now using stretched grid with Jacobians

! Find ustar for wall model eddy viscosity
! 1. Use wall-stress value
if (jt_total > 1) then
#ifdef PPTLWM_LVLSET
    call tlwm_lvlset_wallstress
    ustar = sqrt(abs(txzr_wall(1:nxr,1:nyr)) + abs(tyzr_wall(1:nxr,1:nyr)))
    ! from tlwm_lvlset_wallstress, txzr_wall and tyzr_wall are 
    ! already zero on coords where the wall interface is not located
#else
    ! Compute wall-stress on only bottom coord
    if (coord == 0) then
        ustar = sqrt(abs(txzr(1:nxr,1:nyr,1)) + abs(tyzr(1:nxr,1:nyr,1)))
    else
        ustar = 0._rprec
    endif
#endif

    ! Send wall-stress to all coords
    call mpi_allreduce(ustar, dummy, nxr*nyr, mpi_rprec,                  &
        MPI_SUM, comm, ierr)
    ustar = dummy
else !! first time-step, use assumed nonzero ustar
    ustar = u_star
endif
! 2. Use constant value specified by user
! ustar = u_star

! Compute eddy viscosity values on inner layer grid
call wm_eddyvisc(ustar)

! Prepare the RHS
Rx(:,:,:) = 0._rprec
Ry(:,:,:) = 0._rprec

! Get bottom row
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    do jy = 1, nyr
    do jx = 1, nxr
#ifdef PPTLWM_LVLSET
        if (phi_uvr(jx,jy,1) > 0) then !! in fluid
             b(jx,jy,1) = -const1*(1._rprec/jaco_uvr(1))*                        &
                (const1*(nu_wr(jx,jy,2)+nu_molec)/jaco_wr(2) + nu_molec/phi_uvr(jx,jy,1))
             c(jx,jy,1) = const1*(1._rprec/jaco_uvr(1))*                         &
                const1*(nu_wr(jx,jy,2)+nu_molec)/jaco_wr(2)

            Rx(jx,jy,1) = Rx(jx,jy,1) - const1*(1._rprec/jaco_uvr(1))*          &
                (nu_molec/zuvr(1))*ubot
            ! Ry(jx,jy,1) = Ry(jx,jy,1) !! vbot = 0, so do nothing           
        else !! in solid
            ! Change coefficient to manually specify what the velocity is
            b(jx,jy,1) = 1._rprec
            c(jx,jy,1) = 0._rprec
            ! Force velocity values to be known values
            Rx(jx,jy,1) = ubot
            ! Leave Ry as zero
        endif
#else
        ! Coefficients are different on bottom row because using one-sided
        ! finite difference at wall, also nu_t(z=0)=0
        b(jx,jy,1) = -const1*(1._rprec/jaco_uvr(1))*                        &
            (const1*(nu_wr(jx,jy,2)+nu_molec)/jaco_wr(2) + nu_molec/zuvr(1))
        c(jx,jy,1) = const1*(1._rprec/jaco_uvr(1))*                         &
            const1*(nu_wr(jx,jy,2)+nu_molec)/jaco_wr(2)

        Rx(jx,jy,1) = Rx(jx,jy,1) - const1*(1._rprec/jaco_uvr(1))*          &
            (nu_molec/zuvr(1))*ubot
        ! Ry(jx,jy,1) = Ry(jx,jy,1) !! vbot = 0, so do nothing
#endif
    enddo
    enddo
    jz_min = 2
else
    jz_min = 1
endif

! Get top row
if (coord == nproc-1) then
#ifdef PPSAFETYMODE
    c(:,:,nzr-1) = BOGUS
#endif
    do jy = 1, nyr
    do jx = 1, nxr
        a(jx,jy,nzr-1) = const1*(1._rprec/jaco_uvr(nzr-1))*                 &
            const1*(nu_wr(jx,jy,nzr-1)+nu_molec)/jaco_wr(nzr-1)
        b(jx,jy,nzr-1) = -const1*(1._rprec/jaco_uvr(nzr-1))*                &
            (const1*(nu_wr(jx,jy,nzr-1)+nu_molec)/jaco_wr(nzr-1) +          &
            (nu_wr(jx,jy,nzr)+nu_molec)/(L_zr-zuvr(nzr-1)))

        !! interpolated LES velocity as upper BC
        Rx(jx,jy,nzr-1) = Rx(jx,jy,nzr-1) - const1*(1._rprec/jaco_uvr(nzr-1))*   &
            ((nu_wr(jx,jy,nzr)+nu_molec)/(L_zr-zuvr(nzr-1)))*ur(jx,jy,nzr)
        Ry(jx,jy,nzr-1) = Ry(jx,jy,nzr-1) - const1*(1._rprec/jaco_uvr(nzr-1))*   &
            ((nu_wr(jx,jy,nzr)+nu_molec)/(L_zr-zuvr(nzr-1)))*vr(jx,jy,nzr)
    enddo
    enddo
    jz_max = nzr-2
else
    jz_max = nzr-1
endif

! Coefficients in domain
do jz = jz_min, jz_max
do jy = 1, nyr
do jx = 1, nxr
#ifdef PPTLWM_LVLSET
    if (k_wall(jx,jy) == jz) then !! fluid point near wall
        a(jx,jy,jz) = 0._rprec !! moved to rhs
        b(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*                        &
            (const1*(nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1) +              &
            (nu_wr(jx,jy,jz)+nu_molec)/phi_wall(jx,jy))
        c(jx,jy,jz) = const1*(1._rprec/jaco_uvr(jz))*                         &
            const1*(nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1)
        Rx(jx,jy,jz) = Rx(jx,jy,jz) - const1*(1._rprec/jaco_uvr(jz))*         &
            ((nu_wr(jx,jy,jz)+nu_molec)/phi_wall(jx,jy))*ubot
        ! Leave Ry as is
    elseif (phi_uvr(jx,jy,jz) > 0) then !! in fluid away from wall
        a(jx,jy,jz) = const1*(1._rprec/jaco_uvr(jz))*                         &
            const1*(nu_wr(jx,jy,jz)+nu_molec)/jaco_wr(jz)
        b(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*                        &
            (const1*(nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1) +              &
            const1*(nu_wr(jx,jy,jz)+nu_molec)/jaco_wr(jz))
        c(jx,jy,jz) = const1*(1._rprec/jaco_uvr(jz))*                         &
            const1*(nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1)
    else !! in solid
        a(jx,jy,jz) = 0._rprec
        b(jx,jy,jz) = 1._rprec
        c(jx,jy,jz) = 0._rprec
        Rx(jx,jy,jz) = ubot
        ! Leave Ry as is
    endif
#else
    a(jx,jy,jz) = const1*(1._rprec/jaco_uvr(jz))*                         &
        const1*(nu_wr(jx,jy,jz)+nu_molec)/jaco_wr(jz)
    b(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*                        &
        (const1*(nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1) +              &
        const1*(nu_wr(jx,jy,jz)+nu_molec)/jaco_wr(jz))
    c(jx,jy,jz) = const1*(1._rprec/jaco_uvr(jz))*                         &
        const1*(nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1)
#endif
enddo
enddo
enddo

! Solve ur, vr velocities using TDMA
call tridag_array_tlwm_diff(a,b,c,Rx,usol)
call tridag_array_tlwm_diff(a,b,c,Ry,vsol)

! Fill velocity solution
ur(:nxr,:nyr,1:nzr-1) = usol(:nxr,:nyr,1:nzr-1)
vr(:nxr,:nyr,1:nzr-1) = vsol(:nxr,:nyr,1:nzr-1)

end subroutine tlwm_eq_solve

!*****************************************************************************
subroutine tlwm_noneq_solve
!*****************************************************************************
! 
! This subroutine solves the streamwise, spanwise momentum equations, and the
! equation of continuity for a non-equilibrium wall-model boundary condition
! for the outer-layer LES.
! 
use messages
use param, only : nxr, nyr, nzr, lbz, dzr, L_zr
use param, only : nu_molec, z_i, ubot, u_star, jt_total, wbase
use param, only : coord, nproc, comm, ierr, MPI_RPREC, up, down, status
use param, only : dt, tadv1, tadv2, dtr, tadvr1, tadvr2, wm_count, use_tlwm_visc_dt
#ifdef PPSAFETYMODE
use param, only : BOGUS
#endif
use mpi
implicit none
integer :: jx, jy, jz, jz_min, jz_max, jtr
real(rprec), dimension(nxr, nyr) :: ustar, dummy
real(rprec) :: const, const1, maxvisc, rmsdivvel
real(rprec), dimension(nxr,nyr,lbz:nzr) :: diff_coef_uv, diff_coef_w
real(rprec), dimension(nxr+2,nyr,lbz:nzr) :: dtxdx, dtydy, dtzdz,         &
    dtxdx2, dtydy2, dtzdz2
real(rprec), dimension(nxr,nyr,0:nzr) :: a, b, c
real(rprec), dimension(nxr+2,nyr,0:nzr) :: Rx, usol, Ry, vsol

! Determine time-step size, wall-model time-steps are fixed regardless of use_cfl_dt
if (use_tlwm_visc_dt) call get_tlwm_dt(wm_count)
dtr = dt / wm_count
tadvr1 = 1.5_rprec
tadvr2 = -0.5_rprec

do jtr = 1, wm_count

    ! Find ustar for wall model eddy viscosity
    ! 1. Use wall-stress value
    if (jt_total > 1) then
#ifdef PPTLWM_LVLSET
        call tlwm_lvlset_wallstress
        ustar = sqrt(abs(txzr_wall(1:nxr,1:nyr)) + abs(tyzr_wall(1:nxr,1:nyr)))
        ! from tlwm_lvlset_wallstress, txzr_wall and tyzr_wall are
        ! already zero on coords where the wall interface is not located
#else
        ! Compute wall-stress on only bottom coord
        if (coord == 0) then
            ustar = sqrt(abs(txzr(1:nxr,1:nyr,1)) + abs(tyzr(1:nxr,1:nyr,1)))
        else
            ustar = 0._rprec
        endif
#endif
        ! Send wall-stress to all coords
        call mpi_allreduce(ustar, dummy, nxr*nyr, mpi_rprec,                &
            MPI_SUM, comm, ierr)
        ustar = dummy
    else !! first time-step, use assumed nonzero ustar
        ustar = u_star
    endif
    ! 2. Use constant value specified by user
    ! ustar = u_star

    ! Compute eddy viscosity values on inner layer grid
    !call wm_eddyvisc_blend(ustar)
    call wm_eddyvisc(ustar)

    ! Exchange ghost node information
    ! send info down from jz=1 on coord to jz=nzr on coord-1
    call mpi_sendrecv( ur(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,          &
        ur(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)
    call mpi_sendrecv( vr(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,          &
        vr(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)
    ! send info up from jz=nzr-1 on coord to jz=0 on coord+1
    call mpi_sendrecv( ur(:,:,nzr-1), (nxr+2)*nyr, MPI_RPREC, up, 2,        &
        ur(:,:,0), (nxr+2)*nyr, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_sendrecv( vr(:,:,nzr-1), (nxr+2)*nyr, MPI_RPREC, up, 2,        &
        vr(:,:,0), (nxr+2)*nyr, MPI_RPREC, down, 2, comm, status, ierr)

    ! Save previous timestep's RHS
    rhs_xr_f = rhs_xr
    rhs_yr_f = rhs_yr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate derivatives
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute horizontal derivatives of u and v
    call tlwm_filt_da(ur, dudxr, dudyr, lbz)
    call tlwm_filt_da(vr, dvdxr, dvdyr, lbz)

    ! Compute the velocity w that satisfies continuity for given u and v
    call div_calc_w()

    ! Report divergence calculation
    if (( modulo (jt_total, wbase) == 0) .and. (jtr == wm_count) ) then
        call tlwm_rmsdiv(rmsdivvel)
    endif

    ! Now compute horizontal derivatives of w
    call tlwm_filt_da(wr, dwdxr, dwdyr, lbz)

    ! Compute vertical derivatives of u and v
    call tlwm_ddz_uv(ur, dudzr, lbz)
    call tlwm_ddz_uv(vr, dvdzr, lbz)

    ! Compute dudz, dvdz, txz, and tyz at wall on bottom coord
    if (coord == 0) then
        call inner_layer_wallstress()
        ! Add boundary condition to explicit portion
        txzr_half2(1:nxr,:,1) = txzr(1:nxr,:,1)
        tyzr_half2(1:nxr,:,1) = tyzr(1:nxr,:,1)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate stress
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Diffusion coefficient
    diff_coef_uv = -2.0_rprec*(nu_molec + nu_uvr(:,:,:))
    diff_coef_w = -2.0_rprec*(nu_molec + nu_wr(:,:,:))

    ! Calculate stress for bottom level
    if (coord == 0) then
        txxr(1:nxr,:,1) = diff_coef_uv(:,:,1)*dudxr(1:nxr,:,1)
        txyr(1:nxr,:,1) = diff_coef_uv(:,:,1)*0.5_rprec*(dudyr(1:nxr,:,1)+dvdxr(1:nxr,:,1))
        tyyr(1:nxr,:,1) = diff_coef_uv(:,:,1)*dvdyr(1:nxr,:,1)
        ! Remember txz & tyz already calculated in routine inner_layer_wallstress

        ! since first level already calculated
        jz_min = 2
    else
        jz_min = 1
    endif

    ! Calculate stress at top level
    ! in main.f90 this is done in wallstress, however here the only upper BC
    ! is to use the LES velocity, so use derivatives as normal
    !     dudzr and dvdzr are defined as normal on top coord at jz=nzr
    ! txx, txy, and tyy are not needed to step in time, but calculated anyway
    if (coord == nproc-1) then
        ! Stresses on uv-grid
        txxr(1:nxr,:,nzr) = diff_coef_uv(:,:,nzr)*dudxr(1:nxr,:,nzr)
        txyr(1:nxr,:,nzr) = diff_coef_uv(:,:,nzr)*0.5_rprec*(dudyr(1:nxr,:,nzr)+dvdxr(1:nxr,:,nzr))
        tyyr(1:nxr,:,nzr) = diff_coef_uv(:,:,nzr)*dvdyr(1:nxr,:,nzr)

        ! Stresses on w-grid
        txzr(1:nxr,:,nzr) = diff_coef_w(:,:,nzr)*0.5_rprec*(dudzr(1:nxr,:,nzr)+dwdxr(1:nxr,:,nzr))
        txzr_half1(1:nxr,:,nzr) = diff_coef_w(:,:,nzr)*0.5_rprec*( dwdxr(1:nxr,:,nzr) )
        txzr_half2(1:nxr,:,nzr) = diff_coef_w(:,:,nzr)*0.5_rprec*( dudzr(1:nxr,:,nzr) )
        tyzr(1:nxr,:,nzr) = diff_coef_w(:,:,nzr)*0.5_rprec*(dvdzr(1:nxr,:,nzr)+dwdyr(1:nxr,:,nzr))
        tyzr_half1(1:nxr,:,nzr) = diff_coef_w(:,:,nzr)*0.5_rprec*( dwdyr(1:nxr,:,nzr) )
        tyzr_half2(1:nxr,:,nzr) = diff_coef_w(:,:,nzr)*0.5_rprec*( dvdzr(1:nxr,:,nzr) )
    endif

    ! Calculate stress for entire domain
    do jz = jz_min, nzr-1
        ! Stresses on uv-grid
        txxr(1:nxr,:,jz) = diff_coef_uv(:,:,jz)*dudxr(1:nxr,:,jz)
        txyr(1:nxr,:,jz) = diff_coef_uv(:,:,jz)*0.5_rprec*(dudyr(1:nxr,:,jz)+dvdxr(1:nxr,:,jz))
        tyyr(1:nxr,:,jz) = diff_coef_uv(:,:,jz)*dvdyr(1:nxr,:,jz)

        ! Stresses on w-grid
        txzr(1:nxr,:,jz) = diff_coef_w(:,:,jz)*0.5_rprec*(dudzr(1:nxr,:,jz)+dwdxr(1:nxr,:,jz))
        txzr_half1(1:nxr,:,jz) = diff_coef_w(:,:,jz)*0.5_rprec*( dwdxr(1:nxr,:,jz) )
        txzr_half2(1:nxr,:,jz) = diff_coef_w(:,:,jz)*0.5_rprec*( dudzr(1:nxr,:,jz) )
        tyzr(1:nxr,:,jz) = diff_coef_w(:,:,jz)*0.5_rprec*(dvdzr(1:nxr,:,jz)+dwdyr(1:nxr,:,jz))
        tyzr_half1(1:nxr,:,jz) = diff_coef_w(:,:,jz)*0.5_rprec*( dwdyr(1:nxr,:,jz) )
        tyzr_half2(1:nxr,:,jz) = diff_coef_w(:,:,jz)*0.5_rprec*( dvdzr(1:nxr,:,jz) )
    enddo

    ! Exchange information between processors to set values at nz
    ! from jz=1 above to jz=nz below
    call mpi_sendrecv( txzr(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,        &
        txzr(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)
    call mpi_sendrecv( txzr_half1(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,  &
        txzr_half1(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)
    call mpi_sendrecv( txzr_half2(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,  &
        txzr_half2(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)
    call mpi_sendrecv( tyzr(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,        &
        tyzr(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)
    call mpi_sendrecv( tyzr_half1(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,  &
        tyzr_half1(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)
    call mpi_sendrecv( tyzr_half2(:,:,1), (nxr+2)*nyr, MPI_RPREC, down, 1,  &
        tyzr_half2(:,:,nzr), (nxr+2)*nyr, MPI_RPREC, up, 1, comm, status, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Diffusive term
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate divergence of stress
    ! Compute stress gradients
    call tlwm_ddx(txxr, dtxdx, lbz)
    call tlwm_ddy(tyyr, dtydy2, lbz)
    call tlwm_ddxy(txyr, dtxdx2, dtydy, lbz)
    call tlwm_ddz_w(txzr_half1, dtzdz, lbz)
    call tlwm_ddz_w(tyzr_half1, dtzdz2, lbz)
    ! For Adams-Bashforth time-advancement uncomment:
    ! call tlwm_ddz_w(txzr, dtzdz, lbz)
    ! call tlwm_ddz_w(tyzr, dtzdz2, lbz)

    ! Take sum, remember only 1:nzr-1 are valid
    div_txr(:,:,1:nzr-1) = dtxdx(:,:,1:nzr-1) + dtydy(:,:,1:nzr-1)+ dtzdz(:,:,1:nzr-1)
    div_tyr(:,:,1:nzr-1) = dtxdx2(:,:,1:nzr-1) + dtydy2(:,:,1:nzr-1)+ dtzdz2(:,:,1:nzr-1)

    ! Set ld-1 and ld oddballs to 0
    div_txr(nxr+1:nxr+2,:,1:nzr-1) = 0._rprec
    div_tyr(nxr+1:nxr+2,:,1:nzr-1) = 0._rprec

    ! Set boundary points to BOGUS
#ifdef PPSAFETYMODE
    div_txr(:,:,0) = BOGUS
    div_tyr(:,:,0) = BOGUS
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advective term
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Take variables to big domain
    call to_big(ur, ur_big)
    call to_big(vr, vr_big)
    call to_big(dudxr, dudxr_big)
    call to_big(dudyr, dudyr_big)
    call to_big(dudzr, dudzr_big)
    call to_big(dvdxr, dvdxr_big)
    call to_big(dvdyr, dvdyr_big)
    call to_big(dvdzr, dvdzr_big)

    ! Normalization for FFTs
    const = 1._rprec/( (3*nxr/2) * (3*nyr/2) )

    ! Compute advective term in x-momentum
    ! Interpolate w and dudz onto uv-grid 
    if (coord == 0) then
        ! Bottom wall take w(jz=1) = 0
        temp_big(:,:,1) = const*(ur_big(:,:,1)*dudxr_big(:,:,1) +         &
            vr_big(:,:,1)*dudyr_big(:,:,1) +                              &
            0.5_rprec*wr_big(:,:,2)*dudzr_big(:,:,2))
        jz_min = 2
    else
        jz_min = 1
    endif

    ! In other advec routine, would take w(jz=nz)=0 on top coord
    ! however here that is not necessarily true

    ! For entire domain
    do jz = jz_min, nzr-1
        temp_big(:,:,jz) = const*(ur_big(:,:,jz)*dudxr_big(:,:,jz) +    &
            vr_big(:,:,jz)*dudyr_big(:,:,jz) +                            &
            0.5_rprec*(wr_big(:,:,jz+1)*dudzr_big(:,:,jz+1) +             &
            wr_big(:,:,jz)*dudzr_big(:,:,jz)))
    enddo

    ! Move temp_big into RHSx and make small
    call to_small(temp_big, rhs_xr)

    ! Compute advective term in y-momentum
    ! Interpolate w and dvdz onto uv-grid 
    if (coord == 0) then
        ! Bottom wall take w(jz=1) = 0
        temp_big(:,:,1) = const*(ur_big(:,:,1)*dvdxr_big(:,:,1) +         &
            vr_big(:,:,1)*dvdyr_big(:,:,1) +                              &
            0.5_rprec*wr_big(:,:,2)*dvdzr_big(:,:,2))
        jz_min = 2
    else
        jz_min = 1
    endif

    ! In other advec routine, would take w(jz=nz)=0 on top coord
    ! however here that is not necessarily true

    ! For entire domain
    do jz = jz_min, nzr-1
        temp_big(:,:,jz) = const*(ur_big(:,:,jz)*dvdxr_big(:,:,jz) +      &
            vr_big(:,:,jz)*dvdyr_big(:,:,jz) +                            &
            0.5_rprec*(wr_big(:,:,jz+1)*dvdzr_big(:,:,jz+1) +             &
            wr_big(:,:,jz)*dvdzr_big(:,:,jz)))
    enddo

    ! Move temp_big into RHSx and make small
    call to_small(temp_big, rhs_yr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Add terms to the RHS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Add div-tau term
    rhs_xr(:,:,1:nzr-1) = -rhs_xr(:,:,1:nzr-1) - div_txr(:,:,1:nzr-1)
    rhs_yr(:,:,1:nzr-1) = -rhs_yr(:,:,1:nzr-1) - div_tyr(:,:,1:nzr-1)

    ! Add pressure forcing
    ! Use pressure gradients from LES to every wall-normal point
    do jz = 1, nzr-1
        rhs_xr(1:nxr,:,jz) = rhs_xr(1:nxr,:,jz) + dpdxr(1:nxr,1:nyr)
        rhs_yr(1:nxr,:,jz) = rhs_yr(1:nxr,:,jz) + dpdyr(1:nxr,1:nyr)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Euler Integration check
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (jt_total == 1) then
        ! Make rhs_f = rhs for Euler integration
        rhs_xr_f = rhs_xr
        rhs_yr_f = rhs_yr
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate velocity at next time-step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare matrix equation and solve using TDMA
    call tlwm_diff_stag_array_uv()

enddo !! end of main wall-model time-loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write time-iteration info to the screen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (modulo (jt_total, wbase) == 0) then
    call max_tlwm_visc(maxvisc)

    if (coord == 0) then
        write(*,*) 
        write(*,'(a,E15.7)') ' TLWM DIV: ', rmsdivvel
        write(*,'(a,E15.7)') ' TLWM VISC: ', maxvisc
        write(*,'(a,i9)')    ' TLWM Nsteps: ', wm_count
    endif
endif

end subroutine tlwm_noneq_solve

!*******************************************************************************
subroutine tlwm_diff_stag_array_uv
!*******************************************************************************
! 
! Calculate xy-velocity from CN scheme using TDMA. This routine is only for 
! the non-equilibrium wall-model, where a time-derivative is involved.
! 
use messages
use param, only : nxr, nyr, nzr, nu_molec, z_i, ubot, u_star, dzr, L_zr, lbz
use param, only : dtr, tadvr1, tadvr2
use param, only : coord, nproc, comm, ierr, MPI_RPREC
#ifdef PPSAFETYMODE
use param, only : BOGUS
#endif
use mpi
implicit none

integer :: jx, jy, jz, jz_min, jz_max
real(rprec) :: const1, const2, nu_a, nu_b, nu_c
real(rprec), dimension(nxr,nyr,0:nzr) :: a, b, c
real(rprec), dimension(nxr+2,nyr,0:nzr) :: usol, vsol, Rx, Ry
real(rprec), dimension(nxr+2,nyr,lbz:nzr) :: dtxzdz_rhs, dtyzdz_rhs

const1 = (dtr/(2._rprec*dzr))
const2 = (1._rprec/dzr)

! Get the RHS ready for CN
! Initialize with the explicit terms
Rx(:,:,1:nzr-1) = ur(:,:,1:nzr-1) +                                      &
    dtr*(tadvr1*rhs_xr(:,:,1:nzr-1) + tadvr2*rhs_xr_f(:,:,1:nzr-1))
Ry(:,:,1:nzr-1) = vr(:,:,1:nzr-1) +                                      &
    dtr*(tadvr1*rhs_yr(:,:,1:nzr-1) + tadvr2*rhs_yr_f(:,:,1:nzr-1))

! Add explicit portion of CN
call tlwm_ddz_w(txzr_half2, dtxzdz_rhs, lbz)
call tlwm_ddz_w(tyzr_half2, dtyzdz_rhs, lbz)
#ifdef PPSAFETYMODE
dtxzdz_rhs(:,:,0) = BOGUS
dtyzdz_rhs(:,:,0) = BOGUS
#endif
Rx(:,:,1:nzr-1) = Rx(:,:,1:nzr-1) - dtr*0.5_rprec*dtxzdz_rhs(:,:,1:nzr-1)
Ry(:,:,1:nzr-1) = Ry(:,:,1:nzr-1) - dtr*0.5_rprec*dtyzdz_rhs(:,:,1:nzr-1)

! Get bottom row
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    do jy = 1, nyr
    do jx = 1, nxr
#ifdef PPTLWM_LVLSET
        if (phi_uvr(jx,jy,1) > 0) then !! in fluid
            nu_c = nu_wr(jx,jy,2) + nu_molec
            b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_uvr(1))*              &
                (const2*(1._rprec/jaco_wr(2))*nu_c + (nu_molec/phi_uvr(jx,jy,1)))
            c(jx,jy,1) = -const1*(1._rprec/jaco_uvr(1))*                        &
                const2*(1._rprec/jaco_wr(2))*nu_c

            Rx(jx,jy,1) = Rx(jx,jy,1) + const1*(1._rprec/jaco_uvr(1))*          &
                (nu_molec/zuvr(1))*ubot
            ! Ry(jx,jy,1) = Ry(jx,jy,1) !! vbot = 0, so do nothing
        else !! in solid
            ! Change coefficient to manually specify what the velocity is
            b(jx,jy,1) = 1._rprec
            c(jx,jy,1) = 0._rprec
            ! Force velocity values to be known values
            Rx(jx,jy,1) = ubot
            ! Leave Ry as zero
        endif
#else
        nu_c = nu_wr(jx,jy,2) + nu_molec
        b(jx,jy,1) = 1._rprec + const1*(1._rprec/jaco_uvr(1))*              &
            (const2*(1._rprec/jaco_wr(2))*nu_c + (nu_molec/zuvr(1)))
        c(jx,jy,1) = -const1*(1._rprec/jaco_uvr(1))*                        &
            const2*(1._rprec/jaco_wr(2))*nu_c

        Rx(jx,jy,1) = Rx(jx,jy,1) + const1*(1._rprec/jaco_uvr(1))*          &
            (nu_molec/zuvr(1))*ubot
        ! Ry(jx,jy,1) = Ry(jx,jy,1) !! vbot = 0, so do nothing
#endif
    enddo
    enddo
    jz_min = 2
else
    jz_min = 1
endif

! Get top row
if (coord == nproc-1) then
#ifdef PPSAFETYMODE
    c(:,:,nzr-1) = BOGUS
#endif
    do jy = 1, nyr
    do jx = 1, nxr
        nu_a = nu_wr(jx,jy,nzr-1) + nu_molec
        nu_c = nu_wr(jx,jy,nzr) + nu_molec

        a(jx,jy,nzr-1) = -const1*(1._rprec/jaco_uvr(nzr-1))*                    &
            const2*(1._rprec/jaco_wr(nzr-1))*nu_a
        b(jx,jy,nzr-1) = 1._rprec + const1*(1._rprec/jaco_uvr(nzr-1))*          &
            (const2*(1._rprec/jaco_wr(nzr-1))*nu_a + nu_c/(L_zr-zuvr(nzr-1)))

        !! interpolated LES velocity as upper BC
        Rx(jx,jy,nzr-1) = Rx(jx,jy,nzr-1) + const1*(1._rprec/jaco_uvr(nzr-1))*  &
            (nu_c/(L_zr-zuvr(nzr-1)))*ur(jx,jy,nzr)
        Ry(jx,jy,nzr-1) = Ry(jx,jy,nzr-1) + const1*(1._rprec/jaco_uvr(nzr-1))*  &
            (nu_c/(L_zr-zuvr(nzr-1)))*vr(jx,jy,nzr)
    enddo
    enddo
    jz_max = nzr-2
else
    jz_max = nzr-1
endif

! Coefficients in domain
do jz = jz_min, jz_max
do jy = 1, nyr
do jx = 1, nxr
#ifdef PPTLWM_LVLSET
    if (k_wall(jx,jy) == jz) then !! fluid point near wall
        nu_a = (nu_wr(jx,jy,jz) + nu_molec)/phi_wall(jx,jy)
        nu_b = (const2*(nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1)) +          &
            ((nu_wr(jx,jy,jz)+nu_molec)/phi_wall(jx,jy))
        nu_c = nu_wr(jx,jy,jz+1) + nu_molec

        a(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*const2*(1._rprec/jaco_wr(jz))*nu_a
        a(jx,jy,jz) = 0._rprec !! moved to rhs
        b(jx,jy,jz) = 1._rprec + const1*(1._rprec/jaco_uvr(jz))*nu_b
        c(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*const2*(1._rprec/jaco_wr(jz+1))*nu_c
        Rx(jx,jy,jz) = Rx(jx,jy,jz) + const1*(1._rprec/jaco_uvr(jz))*nu_a*ubot
        ! Leave Ry as is
    elseif (phi_uvr(jx,jy,jz) > 0) then !! in fluid, away from wall
        nu_a = nu_wr(jx,jy,jz) + nu_molec
        nu_b = ((nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1)) +                 &
            ((nu_wr(jx,jy,jz)+nu_molec)/jaco_wr(jz))
        nu_c = nu_wr(jx,jy,jz+1) + nu_molec

        a(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*const2*(1._rprec/jaco_wr(jz))*nu_a
        b(jx,jy,jz) = 1._rprec + const1*(1._rprec/jaco_uvr(jz))*const2*nu_b
        c(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*const2*(1._rprec/jaco_wr(jz+1))*nu_c
    else !! in solid
        a(jx,jy,jz) = 0._rprec
        b(jx,jy,jz) = 1._rprec
        c(jx,jy,jz) = 0._rprec
        Rx(jx,jy,jz) = ubot
        ! Leave Ry as is
    endif
#else
    nu_a = nu_wr(jx,jy,jz) + nu_molec
    nu_b = ((nu_wr(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1)) +                 &
        ((nu_wr(jx,jy,jz)+nu_molec)/jaco_wr(jz))
    nu_c = nu_wr(jx,jy,jz+1) + nu_molec

    a(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*const2*(1._rprec/jaco_wr(jz))*nu_a
    b(jx,jy,jz) = 1._rprec + const1*(1._rprec/jaco_uvr(jz))*const2*nu_b
    c(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*const2*(1._rprec/jaco_wr(jz+1))*nu_c
#endif
enddo
enddo
enddo

! Solve ur, vr velocities using TDMA
call tridag_array_tlwm_diff(a,b,c,Rx,usol)
call tridag_array_tlwm_diff(a,b,c,Ry,vsol)

! Fill velocity solution
ur(:nxr,:nyr,1:nzr-1) = usol(:nxr,:nyr,1:nzr-1)
vr(:nxr,:nyr,1:nzr-1) = vsol(:nxr,:nyr,1:nzr-1)

end subroutine tlwm_diff_stag_array_uv

!*****************************************************************************
subroutine tlwm_ddx(f, dfdx, lbz)
!*****************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x using spectral decomposition.
!
use types, only : rprec
use param, only : nxr, nyr, nzr
use fft
use emul_complex, only : OPERATOR(.TLWMMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx
integer :: jz

! Loop through horizontal slices
do jz = lbz, nzr
    !  Use dfdx to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = f(:,:,jz) / (nxr*nyr)
    call dfftw_execute_dft_r2c(tlwm_forw, dfdx(:,:,jz),dfdx(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdx(nxr+1:nxr+2,:,jz) = 0._rprec
    dfdx(:,nyr/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdx to perform complex multiplication
    ! Optimized version for real(eye*kx)=0
    ! only passing imaginary part of eye*k
    dfdx(:,:,jz) = dfdx(:,:,jz) .TLWMMULI. kxr

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(tlwm_back, dfdx(:,:,jz), dfdx(:,:,jz))
enddo

end subroutine tlwm_ddx

!*****************************************************************************
subroutine tlwm_ddy(f, dfdy, lbz)
!*****************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! y using spectral decomposition.
!
use types, only : rprec
use param, only : nxr, nyr, nzr
use fft
use emul_complex, only : OPERATOR(.TLWMMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
integer :: jz

! Loop through horizontal slices
do jz = lbz, nzr
    !  Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdy(:,:,jz) = f(:,:,jz) / (nxr*nyr)
    call dfftw_execute_dft_r2c(tlwm_forw, dfdy(:,:,jz), dfdy(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdy(nxr+1:nxr+2,:,jz) = 0._rprec
    dfdy(:,nyr/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdy(:,:,jz) .TLWMMULI. kyr

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(tlwm_back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine tlwm_ddy

!*****************************************************************************
subroutine tlwm_ddxy (f, dfdx, dfdy, lbz)
!*****************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x and y using spectral decomposition.
!
use types, only : rprec
use param, only : nxr, nyr, nzr
use fft
use emul_complex, only : OPERATOR(.TLWMMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
integer :: jz

! Loop through horizontal slices
do jz = lbz, nzr
    ! Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = f(:,:,jz) / (nxr*nyr)
    call dfftw_execute_dft_r2c(tlwm_forw, dfdx(:,:,jz), dfdx(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdx(nxr+1:nxr+2,:,jz) = 0._rprec
    dfdx(:,nyr/2+1,jz) = 0._rprec

    ! Derivatives: must to y's first here, because we're using dfdx as storage
    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdx(:,:,jz) .TLWMMULI. kyr
    dfdx(:,:,jz) = dfdx(:,:,jz) .TLWMMULI. kxr

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(tlwm_back, dfdx(:,:,jz), dfdx(:,:,jz))
    call dfftw_execute_dft_c2r(tlwm_back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine tlwm_ddxy

!*******************************************************************************
subroutine tlwm_filt_da(f, dfdx, dfdy, lbz)
!*******************************************************************************
!
! This subroutine kills the oddball components in f and computes the partial
! derivative of f with respect to x and y using spectral decomposition.
!
use types, only : rprec
use param, only : nxr, nyr, nzr
use fft
use emul_complex, only : OPERATOR(.TLWMMULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(inout) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
integer :: jz

! loop through horizontal slices
do jz = lbz, nzr
    ! Calculate FFT in place
    f(:,:,jz) = f(:,:,jz) / (nxr*nyr)
    call dfftw_execute_dft_r2c(tlwm_forw, f(:,:,jz), f(:,:,jz))

    ! Kill oddballs in zero padded region and Nyquist frequency
    f(nxr+1:nxr+2,:,jz) = 0._rprec
    f(:,nyr/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdx(:,:,jz) = f(:,:,jz) .TLWMMULI. kxr
    dfdy(:,:,jz) = f(:,:,jz) .TLWMMULI. kyr

    ! Perform inverse transform to get pseudospectral derivative
    ! The oddballs for derivatives should already be dead, since they are for
    ! inverse transform
    call dfftw_execute_dft_c2r(tlwm_back, f(:,:,jz), f(:,:,jz))
    call dfftw_execute_dft_c2r(tlwm_back, dfdx(:,:,jz), dfdx(:,:,jz))
    call dfftw_execute_dft_c2r(tlwm_back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine tlwm_filt_da

!*******************************************************************************
subroutine tlwm_ddz_uv(f, dfdz, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to z using
! 2nd order finite differencing. f is on the uv grid and dfdz is on the w grid.
! The serial version provides dfdz(:,:,2:nz), and the value at jz=1 is not
! touched. The MPI version provides dfdz(:,:,1:nz), except at the bottom
! process it only supplies 2:nz
!
use types, only : rprec
use param, only : nxr, nyr, nzr, dzr, BOGUS
#ifdef PPSAFETYMODE
use param, only : nproc, coord
#endif
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
integer :: jx, jy, jz
real(rprec) :: const

#if defined(PPMPI) && defined(PPSAFETYMODE)
dfdz(:,:,0) = BOGUS
#endif

! Calculate derivative
! The ghost node information is available here
! if coord == 0, dudz(1) will be set elsewhere
do jz = lbz+1, nzr
do jy = 1, nyr
do jx = 1, nxr
    dfdz(jx,jy,jz) = (1._rprec/(dzr*jaco_wr(jz)))*(f(jx,jy,jz)-f(jx,jy,jz-1))
end do
end do
end do

! Not necessarily accurate at bottom boundary
! Set to BOGUS just to be safe
#ifdef PPSAFETYMODE
if (coord == 0) then
    dfdz(1:nxr,:,1) = BOGUS
end if
#endif
! Unlike ddz_uv in the main routine, here dfdz(:,:,nzr) is accurate
! on top coord

end subroutine tlwm_ddz_uv

!*******************************************************************************
subroutine tlwm_ddz_w(f, dfdz, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to z using
! 2nd order finite differencing. f is on the w grid and dfdz is on the uv grid.
! The serial version provides dfdz(:,:,1:nz-1), and the value at jz=1 is not
! touched. The MPI version provides dfdz(:,:,0:nz-1), except at the top and
! bottom processes, which each has has 0:nz, and 1:nz-1, respectively.
!
use types, only : rprec
use param, only : nxr, nyr, nzr, dzr, BOGUS
#ifdef PPSAFETYMODE
#ifdef PPMPI
use param, only : coord
#endif
#endif
implicit none

real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
integer, intent(in) :: lbz
integer :: jx, jy, jz

do jz = lbz, nzr-1
do jy = 1, nyr
do jx = 1, nxr
    dfdz(jx,jy,jz) = (1/(jaco_uvr(jz)*dzr))*(f(jx,jy,jz+1)-f(jx,jy,jz))
end do
end do
end do

#ifdef PPSAFETYMODE
#ifdef PPMPI
! bottom process cannot calculate dfdz(jz=0)
if (coord == 0) then
    dfdz(1:nxr,:,lbz) = BOGUS
endif
#endif
! All processes cannot calculate dfdz(jz=nz)
dfdz(1:nxr,:,nzr) = BOGUS
#endif

end subroutine tlwm_ddz_w

!*****************************************************************************
subroutine inner_layer_wallstress
!*****************************************************************************
!
! This subroutine computes the wall derivatives and stress values for the
! inner-layer velocity field. No-slip is assumed.
! 
! This routine should only be accessed on coord=0.
! 
use param, only : nxr, nyr, nu_molec, z_i, u_star, ubot
integer :: jx, jy

do jy = 1, nyr
do jx = 1, nxr
    ! This is at it is in wallstress.f90/ws_dns_lbc
    dudzr(jx,jy,1) = ( ur(jx,jy,1) - ubot ) / zuvr(1)
    dvdzr(jx,jy,1) = vr(jx,jy,1) / zuvr(1) !! vbot = 0
    txzr(jx,jy,1) = -nu_molec/(z_i*u_star)*dudzr(jx,jy,1)
    tyzr(jx,jy,1) = -nu_molec/(z_i*u_star)*dvdzr(jx,jy,1)
enddo
enddo

end subroutine inner_layer_wallstress

!*****************************************************************************
subroutine to_big(a, a_big)
!*****************************************************************************
! 
! Add padding to variable a for multiplication in physical space and 
! dealiasing
! 

use fft
use param, only : nxr, nyr, lbz, nzr

real(rprec), dimension(nxr+2, nyr, lbz:nzr), intent(in) :: a
real(rprec), dimension(3*nxr/2+2, 3*nyr/2, lbz:nzr), intent(inout) :: a_big
real(rprec), dimension(nxr+2, nyr) :: temp

integer :: jz
real(rprec) :: const

const = 1._rprec/(nxr*nyr)

do jz = lbz, nzr
    temp(:,:) = const*a(:,:,jz)
    call dfftw_execute_dft_r2c(tlwm_forw, temp, temp)
    call tlwm_padd(a_big(:,:,jz), temp)
    call dfftw_execute_dft_c2r(tlwm_back_big, a_big(:,:,jz), a_big(:,:,jz))
enddo

end subroutine to_big

!*****************************************************************************
subroutine to_small(a_big, a)
!*****************************************************************************
! 
! Undo padding to variable after multiplication in physical space
! 

use fft
use param, only : nxr, nyr, lbz, nzr

real(rprec), dimension(nxr+2, nyr, lbz:nzr), intent(inout) :: a
real(rprec), dimension(3*nxr/2+2, 3*nyr/2, lbz:nzr), intent(inout) :: a_big

integer :: jz

! Normalization constant for transforms is accounted for outside
! this routine
do jz = 1, nzr-1
    call dfftw_execute_dft_r2c(tlwm_forw_big, a_big(:,:,jz), a_big(:,:,jz))
    call tlwm_unpadd(a(:,:,jz), a_big(:,:,jz))
    call dfftw_execute_dft_c2r(tlwm_back, a(:,:,jz), a(:,:,jz))
enddo

end subroutine to_small

!*****************************************************************************
subroutine tlwm_padd(u_big,u)
!*****************************************************************************
! puts arrays into larger, zero-padded arrays
! automatically zeroes the oddballs
use types, only : rprec
use param, only : nxr, nyr
implicit none

!  u and u_big are interleaved as complex arrays
real(rprec), dimension(nxr+2,nyr), intent(in) :: u
real(rprec), dimension(3*nxr/2+2,3*nyr/2), intent(out) :: u_big

integer :: ny_h, j_s, j_big_s

ny_h = nyr/2

! make sure the big array is zeroed!
u_big(:,:) = 0._rprec

! note: split access in an attempt to maintain locality
u_big(:nxr,:ny_h) = u(:nxr,:ny_h)

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = (3*nyr/2) - ny_h + 2

u_big(:nxr,j_big_s:(3*nyr/2)) = u(:nxr,j_s:nyr)

end subroutine tlwm_padd

!*****************************************************************************
subroutine tlwm_unpadd(cc,cc_big)
!*****************************************************************************
use types, only : rprec
use param, only : nxr, nyr
implicit none

!  cc and cc_big are interleaved as complex arrays
real(rprec), dimension( nxr+2, nyr ) :: cc
real(rprec), dimension( 3*nxr/2+2, 3*nyr/2 ) :: cc_big

integer :: ny_h, j_s, j_big_s

ny_h = nyr/2

cc(:nxr,:ny_h) = cc_big(:nxr,:ny_h)

! oddballs
cc(nxr+1:nxr+2,:) = 0._rprec
cc(:,ny_h+1) = 0._rprec

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = (3*nyr/2) - ny_h + 2
cc(:nxr,j_s:nyr) = cc_big(:nxr,j_big_s:(3*nyr/2))

end subroutine tlwm_unpadd

!*****************************************************************************
subroutine div_calc_w
!*****************************************************************************
!
! Compute w that satisfies continuity for a given u and v. Only one boundary
! condition is imposed, no-penetration at wall.
! 
! This subroutine assumes dudx and dvdy have already been calculated.
! w is on the w-grid upon exit.
! 
use messages
use param, only : nxr, nyr, nzr, dzr, coord, nproc
use param, only : status, comm, ierr, MPI_RPREC, up, down
use mpi
implicit none
integer :: jx, jy, jz
real(rprec), dimension(nxr+2,nyr,0:nzr) :: rhs

! Prepare the RHS
rhs(:,:,:) = -( dudxr(:,:,:) + dvdyr(:,:,:) )

! No-penetration boundary condition
if (coord == 0) then
    wr(1:nxr,1:nyr,1) = 0._rprec
endif

! This ODE is solved starting at the bottom wall, then moving up
! wait for wr(:,:,1) from "down"
if (coord /= 0) then
    call mpi_recv(wr(1,1,1), (nxr+2)*nyr, MPI_RPREC, down, 8, comm, status, ierr)
endif

! Solve ODE on given coord with known lower boundary (at jz=1)
do jz = 2, nzr
    wr(1:nxr,1:nyr,jz) = (jaco_uvr(jz-1)*dzr)*rhs(1:nxr,1:nyr,jz-1) +     &
        wr(1:nxr,1:nyr,jz-1)
enddo

! Send wr(:,:,nzr) to "up" 
if (coord /= nproc-1) then
    call mpi_send(wr(1,1,nzr), (nxr+2)*nyr, MPI_RPREC, up, 8, comm, ierr)
endif

end subroutine div_calc_w

!*****************************************************************************
subroutine interp_les_to_tlwm(u1,u2)
!*****************************************************************************
!
! This subroutine takes a flow variable from the LES field (u1) at plane 
! z=hwm, then interpolates onto the TLWM field (u2) as the upper BC
!
use param, only : ld, nx, ny, L_x, L_y, dx, dy
use param, only : nxr, nyr, dxr, dyr
use functions, only : binary_search
implicit none
integer :: i, j, i1, i2, j1, j2
real(rprec), allocatable, dimension(:) :: x1, y1, x2, y2
real(rprec) :: ax, ay, bx, by
real(rprec), dimension(nx,ny), intent(in) :: u1
real(rprec), dimension(nxr,nyr), intent(out) :: u2

! Recreate LES grid
allocate(x1(nx), y1(ny))
do i = 1, nx
    x1(i) = (i-1) * L_x/nx
enddo
do j = 1, ny
    y1(j) = (j-1) * L_y/ny
enddo

! Create TLWM grid
allocate(x2(nxr), y2(nyr))
do i = 1, nxr
    x2(i) = (i-1) * L_x/nxr
enddo
do j = 1, nyr
    y2(j) = (j-1) * L_y/nyr
enddo

! Initialize output
u2(:,:) = 0._rprec

! Perform bilinear interpolation
do i = 1, nxr
    i1 = binary_search(x1,x2(i))
    i2 = mod(i1,nx) + 1
    bx = (x2(i) - x1(i1))/dx
    ax = 1._rprec - bx
    do j = 1, nyr
        j1 = binary_search(y1,y2(j))
        j2 = mod(j1,ny) + 1
        by = (y2(j) - y1(j1))/dy
        ay = 1._rprec - by

        u2(i,j) = ax*ay*u1(i1,j1) + bx*ay*u1(i2,j1) +                    &
            ax*by*u1(i1,j2) + bx*by*u1(i2,j2)
    enddo
enddo

end subroutine interp_les_to_tlwm

!*****************************************************************************
subroutine interp_tlwm_to_les(u1,u2)
!*****************************************************************************
!
! This subroutine takes a flow variable from the TLWM field (u1) at plane 
! z=0, then interpolates onto the LES field (u2) as the lower BC.
! 
! It should be noted that this routine is simpler than interp_les_to_tlwm
! because both TLWM and the LES share wall information on coord = 0
!
use param, only : ld, nx, ny, L_x, L_y, dx, dy
use param, only : nxr, nyr, dxr, dyr
use functions, only : binary_search
implicit none
integer :: i, j, i1, i2, j1, j2
real(rprec), allocatable, dimension(:) :: x1, y1, x2, y2
real(rprec) :: ax, ay, bx, by
real(rprec), dimension(nxr+2,nyr), intent(in) :: u1
real(rprec), dimension(nx,ny), intent(out) :: u2

! Recreate TLWM grid
allocate(x1(nxr), y1(nyr))
do i = 1, nxr
    x1(i) = (i-1) * L_x/nxr
enddo
do j = 1, nyr
    y1(j) = (j-1) * L_y/nyr
enddo

! Create LES grid
allocate(x2(nx), y2(ny))
do i = 1, nx
    x2(i) = (i-1) * L_x/nx
enddo
do j = 1, ny
    y2(j) = (j-1) * L_y/ny
enddo

! Initialize output
u2(:,:) = 0.0_rprec

! Perform bilinear interpolation
do i = 1, nx
    i1 = binary_search(x1,x2(i))
    i2 = mod(i1,nxr) + 1
    bx = (x2(i) - x1(i1))/dxr
    ax = 1._rprec - bx
    do j = 1, ny
        j1 = binary_search(y1,y2(j))
        j2 = mod(j1,nyr) + 1
        by = (y2(j) - y1(j1))/dyr
        ay = 1._rprec - by

        u2(i,j) = ax*ay*u1(i1,j1) + bx*ay*u1(i2,j1) +                    &
            ax*by*u1(i1,j2) + bx*by*u1(i2,j2)
    enddo
enddo

end subroutine interp_tlwm_to_les

!******************************************************************************
subroutine wm_eddyvisc(ustar)
!******************************************************************************
! 
! This subroutine represents the eddy viscosity function used in the two layer
! wall model. It returns the eddy viscosity value for the entire inner layer 
! domain on the uv-grid: nu_urv and on the w-grid: nu_wr
!

use types, only : rprec
use param, only : nu_molec, nxr, nyr, nzr, vonk, lbz
implicit none
real(rprec), dimension(nxr, nyr), intent(in) :: ustar
real(rprec) :: a_plus
integer :: i, j, k

! Model parameters
a_plus = 17.0_rprec !! used in decay portion

! Compute eddy viscosity
do i = 1, nxr
do j = 1, nyr
do k = lbz, nzr
#ifdef PPTLWM_LVLSET
! Using phi values to specify distance from wall
    if (phi_uvr(i,j,k) > 0._rprec) then !! in fluid
        nu_uvr(i,j,k) = vonk*ustar(i,j)*phi_uvr(i,j,k)*                          &
            ((1.0_rprec - exp(-phi_uvr(i,j,k)*ustar(i,j)/(a_plus*nu_molec)))**2)
    else !! phi_uvr <= 0 !! in solid
        nu_uvr(i,j,k) = 0._rprec
    endif
    if (phi_wr(i,j,k) > 0._rprec) then !! in fluid
        nu_wr(i,j,k) = vonk*ustar(i,j)*phi_wr(i,j,k)*                            &
            ((1.0_rprec - exp(-phi_wr(i,j,k)*ustar(i,j)/(a_plus*nu_molec)))**2)
    else !! phi_wr <= 0 !! in solid
        nu_wr(i,j,k) = 0._rprec
    endif
#else
    nu_uvr(i,j,k) = vonk*ustar(i,j)*zuvr(k)*                             &
        ((1.0_rprec - exp(-zuvr(k)*ustar(i,j)/(a_plus*nu_molec)))**2)
    nu_wr(i,j,k) = vonk*ustar(i,j)*zwr(k)*                               &
        ((1.0_rprec - exp(-zwr(k)*ustar(i,j)/(a_plus*nu_molec)))**2)
#endif
end do
end do
end do

end subroutine wm_eddyvisc

!******************************************************************************
subroutine wm_eddyvisc_blend(ustar)
!******************************************************************************
! 
! This subroutine outputs the eddy viscosity for the inner-layer wall-model
! on both the uv- and w-grid. A blended approach is used, if the wall-model 
! grid is coarse, a RANS mixing length is used (vonk*z), if not a LES mixing 
! length is used (C0*delta) with delta = (dx*dy*dz)^(1/3). A van Driest
! decay portion is included to ensure the eddy viscosity approaches zero
! at the wall.
!

use types, only : rprec
use param, only : nu_molec, nxr, nyr, nzr, vonk, lbz, Co, dxr, dyr, dzr
implicit none
real(rprec), dimension(nxr, nyr), intent(in) :: ustar
real(rprec) :: a_plus, decay, Lm
integer :: i, j, k

! Model parameters
a_plus = 17.0_rprec !! used in decay portion

! Compute eddy viscosity
do i = 1, nxr
do j = 1, nyr
do k = lbz, nzr
#ifdef PPTLWM_LVLSET
    if (phi_uvr(i,j,k) > 0._rprec) then !! in fluid
        decay = (1.0_rprec - exp(-phi_uvr(i,j,k)*ustar(i,j)/(a_plus*nu_molec)))**2
        Lm = min(vonk*phi_uvr(i,j,k),Co*(dxr*dyr*jaco_uvr(k)*dzr)**(1._rprec/3._rprec))
        nu_uvr(i,j,k) = ustar(i,j)*Lm*decay
    else !! phi_uvr <= 0 !! in solid
        nu_uvr(i,j,k) = 0._rprec
    endif
    if (phi_wr(i,j,k) > 0._rprec) then !! in fluid
        decay = (1.0_rprec - exp(-phi_wr(i,j,k)*ustar(i,j)/(a_plus*nu_molec)))**2
        Lm = min(vonk*phi_wr(i,j,k),Co*(dxr*dyr*jaco_wr(k)*dzr)**(1._rprec/3._rprec))
        nu_wr(i,j,k) = ustar(i,j)*Lm*decay
    else !! phi_wr <= 0 !! in solid
        nu_wr(i,j,k) = 0._rprec
    endif
#else
    decay = (1.0_rprec - exp(-zuvr(k)*ustar(i,j)/(a_plus*nu_molec)))**2
    Lm = min(vonk*zuvr(k),Co*(dxr*dyr*jaco_uvr(k)*dzr)**(1._rprec/3._rprec))
    nu_uvr(i,j,k) = ustar(i,j)*Lm*decay

    decay = (1.0_rprec - exp(-zwr(k)*ustar(i,j)/(a_plus*nu_molec)))**2
    Lm = min(vonk*zwr(k),Co*(dxr*dyr*jaco_wr(k)*dzr)**(1._rprec/3._rprec))
    nu_wr(i,j,k) = ustar(i,j)*Lm*decay
#endif
end do
end do
end do

end subroutine wm_eddyvisc_blend

!******************************************************************************
subroutine tlwm_rmsdiv(rms_global)
!******************************************************************************
! 
! This subroutine calculates the velocity divergence metric for the inner-layer
! wall-model. Uses the L_inf norm of the velocity divergence.
! 
use types, only : rprec
use param, only : nxr, nyr, nzr, lbz
use mpi
use param, only : ierr, MPI_RPREC
implicit none
real(rprec) :: rms
real(rprec), intent(out) :: rms_global

! Horizontal derivatives of u and v already calculated
! Calculate wall-normal derivative of w
call tlwm_ddz_w(wr, dwdzr, lbz)

rms = maxval( abs( dudxr(1:nxr,1:nyr,1:nzr-1) +             &
    dvdyr(1:nxr,1:nyr,1:nzr-1) + dwdzr(1:nxr,1:nyr,1:nzr-1) ) )

! Transfer between processors
call mpi_allreduce(rms, rms_global, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)

end subroutine tlwm_rmsdiv

!******************************************************************************
subroutine get_tlwm_dt(wm_count)
!******************************************************************************
!
! This subroutine determines the smallest number of wall-model time-steps
! are needed to satisfy the specified viscous stability ratio
!
use param, only : dxr, dyr, dzr, nxr, nyr, nzr, nu_molec, visc_ratio, dt
use mpi
use param, only : ierr, MPI_RPREC
implicit none
integer, intent(out) :: wm_count
real(rprec) :: visc_x, visc_y, visc_z, nu_eff, visc, visc_global
real(rprec), dimension(1:nzr-1) :: visc_z_temp
integer :: jz

! In wall-model, nu_t is defined on uv and w grids
! Only using nu_t on the uv grid for simplicity
nu_eff = maxval( abs(nu_uvr(1:nxr,1:nyr,1:nzr-1)) ) + nu_molec
visc_x = nu_eff / (dxr**2)
visc_y = nu_eff / (dyr**2)
do jz = 1, nzr-1
    visc_z_temp(jz) = nu_eff / ((jaco_uvr(jz)*dzr)**2)
enddo
visc_z = maxval( visc_z_temp(1:nzr-1) )
visc = maxval( (/ visc_x, visc_y, visc_z /) )
call mpi_allreduce(visc, visc_global, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)
wm_count = max( nint( visc_global*dt/visc_ratio ), 1 )

end subroutine get_tlwm_dt

!******************************************************************************
subroutine max_tlwm_visc(visc_global)
!******************************************************************************
! 
! This subroutine provides the value of the maximum viscous stability term,
! (nu + nu_t)*dt / (dx^2) for the non-equilibrium wall-model
! 
use param, only : dtr, dxr, dyr, dzr, nxr, nyr, nzr, nu_molec
use mpi
use param, only : ierr, MPI_RPREC
implicit none
real(rprec), intent(out) :: visc_global
real(rprec) :: visc_x, visc_y, visc_z, nu_eff, visc
real(rprec), dimension(1:nzr-1) :: visc_z_temp
integer :: jz

! In wall-model, nu_t is defined on uv and w grids
! Only using nu_t on the uv grid for simplicity
nu_eff = maxval( abs(nu_uvr(1:nxr,1:nyr,1:nzr-1)) ) + nu_molec
visc_x = nu_eff / (dxr**2)
visc_y = nu_eff / (dyr**2)
do jz = 1, nzr-1
    visc_z_temp(jz) = nu_eff / ((jaco_uvr(jz)*dzr)**2)
enddo
visc_z = maxval( visc_z_temp(1:nzr-1) )
visc = dtr * maxval( (/ visc_x, visc_y, visc_z /) )
call mpi_allreduce(visc, visc_global, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)

end subroutine max_tlwm_visc

!******************************************************************************
subroutine tridag_array_tlwm_diff(a, b, c, r, u)
!******************************************************************************
use types, only : rprec
use param
implicit none

real(rprec), dimension(nxr,nyr,0:nzr), intent(in) :: a, b, c
real(rprec), dimension(nxr+2,nyr,0:nzr), intent(in) :: r
real(rprec), dimension(nxr+2,nyr,0:nzr), intent(out) :: u

integer :: nchunks, chunksize
integer :: cstart, cend
integer :: jx, jy, j
integer :: tag0, q
integer :: j_minf, j_minb, j_maxf, j_maxb

real(rprec) :: bet(nxr, nyr)
real(rprec), dimension(nxr,nyr,0:nzr) :: gam

! Initialize variables
nchunks = nyr

! make sure nchunks divides ny evenly
chunksize = nyr / nchunks

if (coord == 0) then
    do jy = 1, nyr
        do jx = 1, nxr
#ifdef PPSAFETYMODE
            if (b(jx, jy, 1) == 0._rprec) then
                write (*, *) 'tridag_array_tlwm_diff: rewrite eqs, jx, jy= ', jx, jy
                stop
            end if
#endif
            u(jx,jy,1) = r(jx,jy,1) / b(jx,jy,1)
        end do
    end do
    bet = b(:, :, 1)
    j_minf = 2 ! this is only for forward pass
else
    j_minf = 1 ! this is only for forward pass
end if

j_maxf = nzr-1 ! this is only for forward pass
j_minb = 1 ! this is only for backward pass

if (coord == nproc-1) then
    j_maxb = nzr-2 ! this is only for backward pass
else
    j_maxb = nzr-1 ! this is only for backward pass
end if

do q = 1, nchunks
    cstart = 1 + (q - 1) * chunksize
    cend = cstart + chunksize - 1
    tag0 = 0 + 10 * (q - 1)

    if (coord /= 0) then
        ! wait for c(:,:,0), bet(:,:), u(:,:,0) from "down"
        call mpi_recv(c(1, cstart, 0), nxr*chunksize, MPI_RPREC, down, tag0+1,  &
                       comm, status, ierr)
        call mpi_recv(bet(1, cstart), nxr*chunksize, MPI_RPREC, down, tag0+2,   &
                       comm, status, ierr)
        call mpi_recv(u(1, cstart, 0), (nxr+2)*chunksize, MPI_RPREC, down, tag0+3,  &
                       comm, status, ierr)
    end if

    do j = j_minf, j_maxf
        do jy = cstart, cend
            do jx = 1, nxr
                gam(jx, jy, j) = c(jx, jy, j-1) / bet(jx, jy)
                bet(jx, jy) = b(jx, jy, j) - a(jx, jy, j)*gam(jx, jy, j)

#ifdef PPSAFETYMODE
                if (bet(jx, jy) == 0._rprec) then
                    write (*, *) 'tridag_array_tlwm_diff failed at jx,jy,j=', jx, jy, j
                    write (*, *) 'a,b,c,gam,bet=', a(jx, jy, j), b(jx, jy, j), &
                        c(jx, jy, j), gam(jx, jy, j), bet(jx, jy)
                    stop
                end if
#endif
                u(jx, jy, j) = (r(jx, jy, j) - a(jx, jy, j) *            &
                u(jx, jy, j-1)) / bet(jx, jy)
            end do
        end do
    end do

    if (coord /= nproc - 1) then
        ! send c(nzr-1), bet, u(nzr-1) to "up"
        call mpi_send (c(1, cstart, nzr-1), nxr*chunksize, MPI_RPREC, up, tag0+1, &
                       comm, ierr)
        call mpi_send (bet(1, cstart), nxr*chunksize, MPI_RPREC, up, tag0+2,    &
                       comm, ierr)
        call mpi_send (u(1, cstart, nzr-1), (nxr+2)*chunksize, MPI_RPREC, up, tag0+3, &
                       comm, ierr)                   
    end if
end do

do q = 1, nchunks
    cstart = 1 + (q - 1) * chunksize
    cend = cstart + chunksize - 1
    tag0 = 0 + 10 * (q - 1)

    if (coord /= nproc - 1) then  
        ! wait for u(n), gam(n) from "up"
        call mpi_recv (u(1, cstart, nzr), (nxr+2)*chunksize, MPI_RPREC, up, tag0+4,   &
                       comm, status, ierr)
        call mpi_recv (gam(1, cstart, nzr), nxr*chunksize, MPI_RPREC, up, tag0+5, &
                       comm, status, ierr)
    end if

    do j = j_maxb, j_minb, -1
        do jy = cstart, cend
            do jx = 1, nxr
                u(jx, jy, j) = u(jx, jy, j) - gam(jx, jy, j+1) *         &
                    u(jx, jy, j+1)
            end do
        end do
    end do

    ! send u(2), gam(2) to "down"
    call mpi_send (u(1, cstart, 1), (nxr+2)*chunksize, MPI_RPREC, down, tag0+4,     &
                   comm, ierr)
    call mpi_send (gam(1, cstart, 1), nxr*chunksize, MPI_RPREC, down, tag0+5,   &
                   comm, ierr)

end do

end subroutine tridag_array_tlwm_diff

!*****************************************************************************
subroutine tlwm_checkpoint
!*****************************************************************************
!
! This subroutine saves checkpoint variables for MFM analysis
!
use param, only : nzr, write_endian

open(11, file=fname, form='unformatted', convert=write_endian,              &
    status='unknown', position='rewind')
write (11) ur(:,:,1:nzr), vr(:,:,1:nzr), rhs_xr(:,:,1:nzr), rhs_yr(:,:,1:nzr)
close(11)

end subroutine tlwm_checkpoint

#ifdef PPOUTPUT_WMLES
!*****************************************************************************
subroutine tlwm_inst_write()
!*****************************************************************************
! 
! This subroutine writes an instantaneous snapshot of the inner layer velocity
! field. This is similar to io.f90/inst_write, but only for tlwmles and does 
! not write for plane snapshots, only for the full domain.
! 
! Also, only writes binary output, not CGNS.
! 
use param, only : nxr, nyr, nzr, lbz, write_endian, jt_total, path
use param, only : coord, nproc, comm, ierr, MPI_RPREC, up, down, status
use string_util, only : string_splice, string_concat
use mpi
implicit none
integer :: jz
character(64) :: fname, bin_ext
real(rprec), dimension(1:nxr,1:nyr,1:nzr) :: wr_uv, vortxr, vortyr, vortzr, vortzr_w

call string_splice(bin_ext, '.c', coord, '.bin')

! Solve for wr that satisfies continuity for given u and v before writing
call div_calc_w()

! Also take horizontal derivatives of this w to compute vorticity
call tlwm_filt_da(wr, dwdxr, dwdyr, lbz)

! Interpolate wr onto uv-grid
wr_uv(1:nxr,:,1:nzr-1) = 0.5_rprec*(wr(1:nxr,:,2:nzr) + wr(1:nxr,:,1:nzr-1))
! Take care of top boundary
if (coord == nproc-1) wr_uv(1:nxr,:,nzr) = wr_uv(1:nxr,:,nzr-1)
! Sync overlapping data
call mpi_sendrecv( wr_uv(:,:,1), nxr*nyr, MPI_RPREC, down, 1,              &
        wr_uv(:,:,nzr), nxr*nyr, MPI_RPREC, up, 1, comm, status, ierr)

! Compute vorticity
vortxr(1:nxr,:,1:nzr) = dwdyr(1:nxr,:,1:nzr) - dvdzr(1:nxr,:,1:nzr)
vortyr(1:nxr,:,1:nzr) = dudzr(1:nxr,:,1:nzr) - dwdxr(1:nxr,:,1:nzr)
vortzr(1:nxr,:,1:nzr) = dvdxr(1:nxr,:,1:nzr) - dudyr(1:nxr,:,1:nzr)

! Interpolate vortzr onto w-grid
vortzr_w(1:nxr,:,2:nzr) = 0.5_rprec*(vortzr(1:nxr,:,1:nzr-1) + vortzr(1:nxr,:,2:nzr))
! Take care of bottom boundary
if (coord == 0) vortzr_w(1:nxr,:,1) = 0.0_rprec !! no-slip
! Sync overlapping data
call mpi_sendrecv( vortzr_w(:,:,nzr), nxr*nyr, MPI_RPREC, up, 2,           &
    vortzr_w(:,:,1), nxr*nyr, MPI_RPREC, down, 2, comm, status, ierr)

! write binary output
! Instantaneous velocity output (all components on uv-grid)
call string_splice(fname, path // 'output/tlwm_vel.', jt_total)
call string_concat(fname, bin_ext)
open(unit=13, file=fname, form='unformatted', convert=write_endian,         &
    access='direct', recl=nxr*nyr*nzr*rprec)
write(13,rec=1) ur(:nxr,:nyr,1:nzr)
write(13,rec=2) vr(:nxr,:nyr,1:nzr)
write(13,rec=3) wr_uv(:nxr,:nyr,1:nzr)
close(13)

! Instantaneous vorticity output (all components on w-grid)
call string_splice(fname, path // 'output/tlwm_vort.', jt_total)
call string_concat(fname, bin_ext)
open(unit=13, file=fname, form='unformatted', convert=write_endian,         &
    access='direct', recl=nxr*nyr*nzr*rprec)
write(13,rec=1) vortxr(:nxr,:nyr,1:nzr)
write(13,rec=2) vortyr(:nxr,:nyr,1:nzr)
write(13,rec=3) vortzr_w(:nxr,:nyr,1:nzr)
close(13)

! Instantaneous wall-stress output
if (coord == 0) then
    call string_splice(fname, path // 'output/tlwm_wallstress.', jt_total)
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,         &
        access='direct', recl=nxr*nyr*rprec)
#ifdef PPTLWM_LVLSET
    write(13,rec=1) txzr_wall(:nxr,:nyr)
    write(13,rec=2) tyzr_wall(:nxr,:nyr)
#else
    write(13,rec=1) txzr(:nxr,:nyr,1)
    write(13,rec=2) tyzr(:nxr,:nyr,1)
#endif
    close(13)
endif

end subroutine tlwm_inst_write

!*****************************************************************************
subroutine tavg_tlwm_init()
!*****************************************************************************
!
! This subroutine loads the tavg_tlwm.out file or creates it
!
use types, only : rprec
use param, only : read_endian, path, coord
use stat_defs, only : tavg_tlwm_total_time, tavg_tlwm
use string_util
implicit none
character (*), parameter :: ftavg_tlwm_in = path // 'tavg_tlwm.out'
character (128) :: fname
logical :: exst

! Does a tavg file already exist?
inquire (file=ftavg_tlwm_in, exist=exst)
if (.not. exst) then !! No, so initialize tavg_tlwm_total_time
    tavg_tlwm_total_time = 0._rprec
else !! Yes, extract tavg_total_time and tavg data
    fname = ftavg_tlwm_in
    call string_concat( fname, '.c', coord )
    open(1, file=fname, action='read', position='rewind',         &
        form='unformatted', convert=read_endian)
    read(1) tavg_tlwm_total_time
    read(1) tavg_tlwm
    close(1)
end if

end subroutine tavg_tlwm_init

!*****************************************************************************
subroutine tavg_tlwm_compute()
!*****************************************************************************
!
! Running time-average computations for inner-layer wall-model
!

use param, only : nxr, nyr, nzr, lbz
use param, only : coord, nproc, comm, ierr, MPI_RPREC, up, down, status
use stat_defs, only : tavg_tlwm, tavg_tlwm_total_time, tavg_tlwm_dt
implicit none
integer :: i, j, k
real(rprec), dimension(1:nxr,1:nyr,1:nzr) :: wr_uv, vortxr, vortyr, vortzr, vortzr_w

! Solve for wr that satisfies continuity for given u and v before writing
call div_calc_w()

! Also take horizontal derivatives of this w to compute vorticity
call tlwm_filt_da(wr, dwdxr, dwdyr, lbz)

! Interpolate wr onto uv-grid
wr_uv(1:nxr,:,1:nzr-1) = 0.5_rprec*(wr(1:nxr,:,2:nzr) + wr(1:nxr,:,1:nzr-1))
! Take care of top boundary
if (coord == nproc-1) wr_uv(1:nxr,:,nzr) = wr_uv(1:nxr,:,nzr-1)
! Sync overlapping data
call mpi_sendrecv( wr_uv(:,:,1), nxr*nyr, MPI_RPREC, down, 1,              &
        wr_uv(:,:,nzr), nxr*nyr, MPI_RPREC, up, 1, comm, status, ierr)

! Compute vorticity
vortxr(1:nxr,:,1:nzr) = dwdyr(1:nxr,:,1:nzr) - dvdzr(1:nxr,:,1:nzr)
vortyr(1:nxr,:,1:nzr) = dudzr(1:nxr,:,1:nzr) - dwdxr(1:nxr,:,1:nzr)
vortzr(1:nxr,:,1:nzr) = dvdxr(1:nxr,:,1:nzr) - dudyr(1:nxr,:,1:nzr)

! Interpolate vortzr onto w-grid
vortzr_w(1:nxr,:,2:nzr) = 0.5_rprec*(vortzr(1:nxr,:,1:nzr-1) + vortzr(1:nxr,:,2:nzr))
! Take care of bottom boundary
if (coord == 0) vortzr_w(1:nxr,:,1) = 0.0_rprec !! no-slip
! Sync overlapping data
call mpi_sendrecv( vortzr_w(:,:,nzr), nxr*nyr, MPI_RPREC, up, 2,           &
    vortzr_w(:,:,1), nxr*nyr, MPI_RPREC, down, 2, comm, status, ierr)

! time-averaging sum
do j = 1, nyr
do i = 1, nxr
do k = 1, nzr
    tavg_tlwm(i,j,k) % u  = tavg_tlwm(i,j,k) % u + ur(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % v  = tavg_tlwm(i,j,k) % v + vr(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % w  = tavg_tlwm(i,j,k) % w + wr_uv(i,j,k) * tavg_tlwm_dt

    tavg_tlwm(i,j,k) % nu_t = tavg_tlwm(i,j,k) % nu_t + nu_uvr(i,j,k) * tavg_tlwm_dt

    tavg_tlwm(i,j,k) % uu = tavg_tlwm(i,j,k) % uu + ur(i,j,k) * ur(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % vv = tavg_tlwm(i,j,k) % vv + vr(i,j,k) * vr(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % ww = tavg_tlwm(i,j,k) % ww + wr_uv(i,j,k) * wr_uv(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % uv = tavg_tlwm(i,j,k) % uv + ur(i,j,k) * vr(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % uw = tavg_tlwm(i,j,k) % uw + ur(i,j,k) * wr_uv(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % vw = tavg_tlwm(i,j,k) % vw + vr(i,j,k) * wr_uv(i,j,k) * tavg_tlwm_dt

    tavg_tlwm(i,j,k) % vortx = tavg_tlwm(i,j,k) % vortx + vortxr(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % vorty = tavg_tlwm(i,j,k) % vorty + vortyr(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % vortz = tavg_tlwm(i,j,k) % vortz + vortzr_w(i,j,k) * tavg_tlwm_dt

    tavg_tlwm(i,j,k) % vortx2 = tavg_tlwm(i,j,k) % vortx2 + vortxr(i,j,k) * vortxr(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % vorty2 = tavg_tlwm(i,j,k) % vorty2 + vortyr(i,j,k) * vortyr(i,j,k) * tavg_tlwm_dt
    tavg_tlwm(i,j,k) % vortz2 = tavg_tlwm(i,j,k) % vortz2 + vortzr_w(i,j,k) * vortzr_w(i,j,k) * tavg_tlwm_dt

end do
end do
end do

! Update tavg_tlwm_total_time for variable time stepping
tavg_tlwm_total_time = tavg_tlwm_total_time + tavg_tlwm_dt

! Set tavg_wmles_dt back to zero for next increment
tavg_tlwm_dt = 0._rprec

end subroutine tavg_tlwm_compute

!*****************************************************************************
subroutine tavg_tlwm_checkpoint
!*****************************************************************************
!
! This subroutine writes the restart data for the wall model statistics. It
! is called by tavg_tlwm_finalize and tlwm_checkpoint.
!
use param, only : checkpoint_tavg_tlwm_file, write_endian, coord
use stat_defs, only : tavg_tlwm, tavg_tlwm_total_time
use string_util
implicit none
character(64) :: fname

! Write data to tavg_tlwm.out
fname = checkpoint_tavg_tlwm_file
call string_concat( fname, '.c', coord )
open(1, file=fname, action='write', position='rewind', &
    form='unformatted', convert=write_endian)
write(1) tavg_tlwm_total_time
write(1) tavg_tlwm
close(1)

end subroutine tavg_tlwm_checkpoint

!*****************************************************************************
subroutine tavg_tlwm_finalize
!*****************************************************************************
!
! Writes time-averaged data to be outputted.
! Currently only uses a binary file
! 

use stat_defs, only : tavg_tlwm, tavg_tlwm_t, tavg_tlwm_total_time
use param, only : nxr, nyr, nzr, write_endian, path, coord
use string_util
implicit none
character(64) :: fname_tlwm_vel, fname_tlwm_rs, fname_tlwm_nu,              &
    fname_tlwm_vort, fname_tlwm_vortrms
integer :: i, j, k
real(rprec), dimension(nxr,nyr,nzr) :: upup, vpvp, wpwp, upvp, upwp, vpwp,  &
    vortxrms, vortyrms, vortzrms
character(64) :: bin_ext

fname_tlwm_vel = path // 'output/tlwm_veluv_avg'
fname_tlwm_nu = path // 'output/tlwm_nu_t_avg'
fname_tlwm_rs = path // 'output/tlwm_rs'
fname_tlwm_vort = path // 'output/tlwm_vort_avg'
fname_tlwm_vortrms = path // 'output/tlwm_vortrms'

call string_splice(bin_ext, '.c', coord, '.bin')

! Final checkpoint all restart data
call tavg_tlwm_checkpoint()

! Normalize sums by total time
do j = 1, nyr
do i = 1, nxr
do k = 1, nzr
    tavg_tlwm(i,j,k) % u = tavg_tlwm(i,j,k) % u / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % v = tavg_tlwm(i,j,k) % v / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % w = tavg_tlwm(i,j,k) % w / tavg_tlwm_total_time

    tavg_tlwm(i,j,k) % nu_t = tavg_tlwm(i,j,k) % nu_t / tavg_tlwm_total_time

    tavg_tlwm(i,j,k) % uu = tavg_tlwm(i,j,k) % uu / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % vv = tavg_tlwm(i,j,k) % vv / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % ww = tavg_tlwm(i,j,k) % ww / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % uv = tavg_tlwm(i,j,k) % uv / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % uw = tavg_tlwm(i,j,k) % uw / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % vw = tavg_tlwm(i,j,k) % vw / tavg_tlwm_total_time

    tavg_tlwm(i,j,k) % vortx = tavg_tlwm(i,j,k) % vortx / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % vorty = tavg_tlwm(i,j,k) % vorty / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % vortz = tavg_tlwm(i,j,k) % vortz / tavg_tlwm_total_time

    tavg_tlwm(i,j,k) % vortx2 = tavg_tlwm(i,j,k) % vortx2 / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % vorty2 = tavg_tlwm(i,j,k) % vorty2 / tavg_tlwm_total_time
    tavg_tlwm(i,j,k) % vortz2 = tavg_tlwm(i,j,k) % vortz2 / tavg_tlwm_total_time
end do
end do 
end do

! Write binary data
! write time-averaged velocity
call string_concat(fname_tlwm_vel, bin_ext)
open(unit=13, file=fname_tlwm_vel, form='unformatted', convert=write_endian, &
    access='direct',recl=nxr*nyr*nzr*rprec) 
write(13,rec=1) tavg_tlwm(:nxr,:nyr,1:nzr)%u
write(13,rec=2) tavg_tlwm(:nxr,:nyr,1:nzr)%v
write(13,rec=3) tavg_tlwm(:nxr,:nyr,1:nzr)%w
close(13)

! Write time-averaged wall-model eddy viscosity
call string_concat(fname_tlwm_nu, bin_ext)
open(unit=13, file=fname_tlwm_nu, form='unformatted', convert=write_endian,  &
    access='direct',recl=nxr*nyr*nzr*rprec) 
write(13,rec=1) tavg_tlwm(:nxr,:nyr,1:nzr)%nu_t
close(13)

! Compute Reynolds stresses
do j = 1, nyr
do i = 1, nxr
do k = 1, nzr
upup(i,j,k) = tavg_tlwm(i,j,k) % uu - tavg_tlwm(i,j,k) % u * tavg_tlwm(i,j,k) % u
vpvp(i,j,k) = tavg_tlwm(i,j,k) % vv - tavg_tlwm(i,j,k) % v * tavg_tlwm(i,j,k) % v
wpwp(i,j,k) = tavg_tlwm(i,j,k) % ww - tavg_tlwm(i,j,k) % w * tavg_tlwm(i,j,k) % w
upvp(i,j,k) = tavg_tlwm(i,j,k) % uv - tavg_tlwm(i,j,k) % u * tavg_tlwm(i,j,k) % v
upwp(i,j,k) = tavg_tlwm(i,j,k) % uw - tavg_tlwm(i,j,k) % u * tavg_tlwm(i,j,k) % w
vpwp(i,j,k) = tavg_tlwm(i,j,k) % vw - tavg_tlwm(i,j,k) % v * tavg_tlwm(i,j,k) % w
enddo
enddo
enddo

! Output Reynolds stresses
call string_concat(fname_tlwm_rs, bin_ext)
open(unit=13, file=fname_tlwm_rs, form='unformatted', convert=write_endian,  &
    access='direct',recl=nxr*nyr*nzr*rprec)
write(13,rec=1) upup(:nxr,:nyr,1:nzr)
write(13,rec=2) vpvp(:nxr,:nyr,1:nzr)
write(13,rec=3) wpwp(:nxr,:nyr,1:nzr)
write(13,rec=4) upvp(:nxr,:nyr,1:nzr)
write(13,rec=5) upwp(:nxr,:nyr,1:nzr)
write(13,rec=6) vpwp(:nxr,:nyr,1:nzr)
close(13)

! Write time-averaged vorticity
call string_concat(fname_tlwm_vort, bin_ext)
open(unit=13, file=fname_tlwm_vort, form='unformatted', convert=write_endian, &
    access='direct',recl=nxr*nyr*nzr*rprec) 
write(13,rec=1) tavg_tlwm(:nxr,:nyr,1:nzr)%vortx
write(13,rec=2) tavg_tlwm(:nxr,:nyr,1:nzr)%vorty
write(13,rec=3) tavg_tlwm(:nxr,:nyr,1:nzr)%vortz
close(13)

! Compute Vorticity RMS
do j = 1, nyr
do i = 1, nxr
do k = 1, nzr
vortxrms(i,j,k) = tavg_tlwm(i,j,k) % vortx2 - tavg_tlwm(i,j,k) % vortx * tavg_tlwm(i,j,k) % vortx
vortyrms(i,j,k) = tavg_tlwm(i,j,k) % vorty2 - tavg_tlwm(i,j,k) % vorty * tavg_tlwm(i,j,k) % vorty
vortzrms(i,j,k) = tavg_tlwm(i,j,k) % vortz2 - tavg_tlwm(i,j,k) % vortz * tavg_tlwm(i,j,k) % vortz
enddo
enddo
enddo

! Output Vorticity RMS 
call string_concat(fname_tlwm_vortrms, bin_ext)
open(unit=13, file=fname_tlwm_vortrms, form='unformatted', convert=write_endian, &
    access='direct',recl=nxr*nyr*nzr*rprec) 
write(13,rec=1) vortxrms(:nxr,:nyr,1:nzr)
write(13,rec=2) vortyrms(:nxr,:nyr,1:nzr)
write(13,rec=3) vortzrms(:nxr,:nyr,1:nzr)
close(13)

end subroutine tavg_tlwm_finalize
#endif

#ifdef PPTLWM_LVLSET
!*****************************************************************************
subroutine load_tlwm_topography
!*****************************************************************************
!
! Loads topography specified from "hij.dat" and outputs grid locations
! relative to topography, phi>0 in fluid, phi<0 in solid
!
use param, only : nxr, nyr, lbz, nzr, coord
implicit none
real(rprec), dimension(nxr*nyr) :: FIELD1
real(rprec), dimension(nxr,nyr) :: tlwm_hij
integer :: i, j, k

allocate(phi_uvr(nxr,nyr,lbz:nzr),phi_wr(nxr,nyr,lbz:nzr))
allocate(coord_wall(nxr,nyr),k_wall(nxr,nyr),phi_wall(nxr,nyr))

! Read hij.dat file
open(1,file=path//'hij.dat')
do i = 1, nxr*nyr
read(1,*) FIELD1(i)
enddo
close(1)

! Re-sort input data
do i = 1, nxr
do j = 1, nyr
    tlwm_hij(i,j) = FIELD1((i-1)*nyr+j)
enddo
enddo

! Compute output phi values
do i = 1, nxr
do j = 1, nyr
do k = lbz, nzr
    phi_uvr(i,j,k) = zuvr(k) - tlwm_hij(i,j)
    phi_wr(i,j,k) = zwr(k) - tlwm_hij(i,j)
enddo
enddo
enddo

! Output dz distance from topography and first grid point
! Also output coord and wall-normal location of this distance
coord_wall = -1
k_wall = -1
phi_wall = -1

do i = 1, nxr
do j = 1, nyr
    do k = 1, (nzr-1)
        if ((phi_uvr(i,j,k-1)<=0) .and. (phi_uvr(i,j,k)>0)) then
            coord_wall(i,j) = coord
            k_wall(i,j) = k
            phi_wall(i,j) = phi_uvr(i,j,k)
        endif
    enddo
enddo
enddo

end subroutine load_tlwm_topography

!*****************************************************************************
subroutine tlwm_lvlset
!*****************************************************************************
!
! Enforce topography in tlwmles using level-set method. Velocity are zeroed
! depending on the sign of phi
!
use param, only : nxr, nyr, lbz, nzr
implicit none
integer :: i, j, k

do k = lbz, nzr
do j = 1, nyr
do i = 1, nxr
    if (phi_uvr(i,j,k) <= 0) then
        ur(i,j,k) = 0._rprec
        vr(i,j,k) = 0._rprec
    endif
enddo
enddo
enddo
! Don't need to zero out wr-velocity b/c its found by solving continuity

end subroutine tlwm_lvlset

!*****************************************************************************
subroutine tlwm_lvlset_wallstress
!*****************************************************************************
!
! This subroutine computes the wallstress given a roughness topography
! enforced by the level-set method
!
use param, only : nxr, nyr, lbz, nzr, coord
use param, only : ubot, z_i, u_star, nu_molec
implicit none
integer :: i, j, k

do i = 1, nxr
do j = 1, nyr
    if (coord == coord_wall(i,j)) then
        k = k_wall(i,j)
        dudzr_wall(i,j) = ( ur(i,j,k) - ubot ) / phi_wall(i,j)
        dvdzr_wall(i,j) = vr(i,j,k) / phi_wall(i,j)
        txzr_wall(i,j) = -nu_molec/(z_i*u_star)*dudzr_wall(i,j)
        tyzr_wall(i,j) = -nu_molec/(z_i*u_star)*dvdzr_wall(i,j)
    else
        ! Zero these values on the other coords so it can summed 
        ! later to output for LES
        dudzr_wall(i,j) = 0._rprec
        dvdzr_wall(i,j) = 0._rprec
        txzr_wall(i,j) = 0._rprec
        tyzr_wall(i,j) = 0._rprec
    endif
enddo
enddo

end subroutine tlwm_lvlset_wallstress

!*****************************************************************************
subroutine lvlset_les_lbc
!*****************************************************************************
!
! Similar to les_lbc, this subroutine creates the lower boundary condition
! for the LES using the solution of the inner-layer wall-model. Unlike 
! les_lbc, this subroutine needs to be called on all coords since the 
! topography may extend beyond coord=0.
! 
use param, only : nxr, nyr, lbz, nzr, nx, ny
use param, only : coord, comm, ierr, MPI_RPREC
use sim_param, only : dudz, dvdz, txz, tyz
use mpi
use test_filtermodule, only : tlwm_filter
implicit none
integer :: i, j, k
real(rprec), dimension(nxr+2,nyr) :: dummy_in
real(rprec), dimension(nx,ny) :: dummy_out

! Gather data and interpolate on bottom coord
call mpi_allreduce(dudzr_wall, dummy_in, (nxr+2)*nyr, mpi_rprec, MPI_SUM, comm, ierr)
if (coord == 0) then
    call tlwm_filter(dummy_in)
    call interp_tlwm_to_les(dummy_in,dummy_out)
    dudz(1:nx,1:ny,1) = dummy_out
endif

call mpi_allreduce(dvdzr_wall, dummy_in, (nxr+2)*nyr, mpi_rprec, MPI_SUM, comm, ierr)
if (coord == 0) then
    call tlwm_filter(dummy_in)
    call interp_tlwm_to_les(dummy_in,dummy_out)
    dvdz(1:nx,1:ny,1) = dummy_out
endif

call mpi_allreduce(txzr_wall, dummy_in, (nxr+2)*nyr, mpi_rprec, MPI_SUM, comm, ierr)
if (coord == 0) then
    call tlwm_filter(dummy_in)
    call interp_tlwm_to_les(dummy_in,dummy_out)
    txz(1:nx,1:ny,1) = dummy_out
endif

call mpi_allreduce(tyzr_wall, dummy_in, (nxr+2)*nyr, mpi_rprec, MPI_SUM, comm, ierr)
if (coord == 0) then
    call tlwm_filter(dummy_in)
    call interp_tlwm_to_les(dummy_in,dummy_out)
    tyz(1:nx,1:ny,1) = dummy_out
endif

end subroutine lvlset_les_lbc
#endif

end module tlwmles
