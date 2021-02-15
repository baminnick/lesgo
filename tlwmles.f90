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

public :: tlwm_init, ic_tlwm, tlwm_checkpoint, tlwm_wallstress

! Main simulation variables for the TLWMLES
real(rprec), public, dimension(:,:,:), allocatable :: ur, vr, wr
real(rprec), dimension(:,:,:), allocatable :: dudxr, dudyr, dudzr,      &
    dvdxr, dvdyr, dvdzr, dwdxr, dwdyr, dwdzr,                           &
    txxr, txyr, tyyr, txzr, tyzr, tzzr,                                 &
    div_txr, div_tyr, div_tzr,                                          &
    rhs_xr, rhs_yr, rhs_zr, rhs_xr_f, rhs_yr_f, rhs_zr_f,               &
    pr, dprdx, dprdy, dprdz
real(rprec), dimension(:,:,:), allocatable :: nu_r
real(rprec), dimension(:,:,:), allocatable :: ur_big, vr_big, wr_big,   &
    dudxr_big, dudyr_big, dudzr_big,                                    & 
    dvdxr_big, dvdyr_big, dvdzr_big,                                    & 
    dwdxr_big, dwdyr_big, dwdzr_big,                                    & 
    temp_big

! Grid variables
real(rprec), dimension(:), allocatable :: zuvr, zwr, jaco_uvr, jaco_wr

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
allocate ( tyzr(nxr+2, nyr, lbz:nzr) ); tyzr = 0._rprec
allocate ( tzzr(nxr+2, nyr, lbz:nzr) ); tzzr = 0._rprec
allocate ( div_txr(nxr+2, nyr, lbz:nzr) ); div_txr = 0._rprec
allocate ( div_tyr(nxr+2, nyr, lbz:nzr) ); div_tyr = 0._rprec
allocate ( div_tzr(nxr+2, nyr, lbz:nzr) ); div_tzr = 0._rprec
allocate ( nu_r(nxr,nyr,lbz:nzr) ); nu_r = 0._rprec
allocate ( rhs_xr(nxr+2, nyr, lbz:nzr) ); rhs_xr = 0.0_rprec
allocate ( rhs_yr(nxr+2, nyr, lbz:nzr) ); rhs_yr = 0.0_rprec
allocate ( rhs_zr(nxr+2, nyr, lbz:nzr) ); rhs_zr = 0.0_rprec
allocate ( rhs_xr_f(nxr+2, nyr, lbz:nzr) ); rhs_xr_f = 0.0_rprec
allocate ( rhs_yr_f(nxr+2, nyr, lbz:nzr) ); rhs_yr_f = 0.0_rprec
allocate ( rhs_zr_f(nxr+2, nyr, lbz:nzr) ); rhs_zr_f = 0.0_rprec
allocate ( pr(nxr+2, nyr, lbz:nzr) ); pr = 0._rprec
allocate ( dprdx(nxr+2, nyr, nzr) ); dprdx = 0._rprec
allocate ( dprdy(nxr+2, nyr, nzr) ); dprdy = 0._rprec
allocate ( dprdz(nxr+2, nyr, nzr) ); dprdz = 0._rprec

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
allocate ( dwdxr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dwdxr_big = 0._rprec
allocate ( dwdyr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dwdyr_big = 0._rprec
allocate ( dwdzr_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); dwdzr_big = 0._rprec
allocate ( temp_big(3*nxr/2 + 2, 3*nyr/2, lbz:nzr) ); temp_big = 0._rprec

! Create TLWM inner layer wall-normal grid
if (str_on) then
    call tlwm_grid_stretched()
else
    call tlwm_grid_uniform()
endif

! Additional TLWMLES Variables

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
use stat_defs, only : tavg_initialized, tavg_wmles_dt
#endif
implicit none
character(*), parameter :: sub_name = 'tlwm_wallstress'

! Calculate wall stress
select case (lbc_mom)
    ! Standard equilibrium wall-model
    case (5)
        call tlwm_eq_ubc() !! Compute u and v upper BC for tlwm
        call tlwm_eq_solve() !! Solve for wall-model velocities

    ! Otherwise, invalid
    case default
        call error (sub_name, 'invalid lbc_mom')
end select

! Compute and interpolate/filter wall-stress answer for LES using TLWM sol
if (coord == 0) call les_lbc()

! Determine if we are to checkpoint intermediate times
if (checkpoint_data) then
    ! Now check if data should be checkpointed this time step
    if ( modulo (jt_total, checkpoint_nskip ) == 0) call tlwm_checkpoint()
end if

#ifdef PPOUTPUT_WMLES
! OUTPUT LOOP (for TLWMLES)

! Instantaneous Domain Velocities
if (domain_calc) then
    if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and.  &
        (mod(jt_total-domain_nstart,domain_nskip)==0) ) then
        call wmles_inst_write()
    endif
endif

! Time-averaging statistics
! Determine if time summations are to be calculated
if (tavg_calc) then
    ! Are we between the start and stop timesteps?
    if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then
        ! Every timestep (between nstart and nend), add to tavg_wmles_dt
        tavg_wmles_dt = tavg_wmles_dt + dt

        ! Are we at the beginning or a multiple of nstart?
        if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then
            ! Check if we have initialized tavg
            if (.not.tavg_initialized) then
                call tavg_wmles_init()
            else
                call tavg_wmles_compute ()
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
    zwr(lbz) = - zwr(lbz)
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
use grid_m
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
#ifdef PPMAPPING
use sim_param, only : mesh_stretch
#endif

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

#ifdef PPMPI
! Exchange ghost node information for u, v, and w
!call mpi_sync_real_array( ur, 0, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( vr, 0, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( wr, 0, MPI_SYNC_DOWNUP )

!--set 0-level velocities to BOGUS
if (coord == 0) then
    ur(:, :, lbz) = BOGUS
    vr(:, :, lbz) = BOGUS
    wr(:, :, lbz) = BOGUS
end if
#endif

end subroutine ic_tlwm

!*****************************************************************************
subroutine ic_tlwm_file
!*****************************************************************************
! Read initial profile for TLWM equation from a file
use param, only : nzr, read_endian
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use grid_m

open(12, file=fname, form='unformatted', convert=read_endian)
read(12) ur(:,:,1:nzr), vr(:,:,1:nzr), wr(:,:,1:nzr),           &
    rhs_xr(:,:,1:nzr), rhs_yr(:,:,1:nzr), rhs_zr(:,:,1:nzr)
close(12)

#ifdef PPMPI
!call mpi_sync_real_array(ur, 0, MPI_SYNC_DOWNUP)
!call mpi_sync_real_array(vr, 0, MPI_SYNC_DOWNUP)
!call mpi_sync_real_array(wr, 0, MPI_SYNC_DOWNUP)
!call mpi_sync_real_array(rhs_xr, 0, MPI_SYNC_DOWNUP)
!call mpi_sync_real_array(rhs_yr, 0, MPI_SYNC_DOWNUP)
!call mpi_sync_real_array(rhs_zr, 0, MPI_SYNC_DOWNUP)
#endif

end subroutine ic_tlwm_file

!*****************************************************************************
subroutine ic_tlwm_vel
!*****************************************************************************
!
! This subroutine initializes the initial tlwm profile 
!
use param
implicit none

! Enforce interpolated/filtered LES velocity field as upper BC for TLWM
call tlwm_eq_ubc()

! Solve equilibrium wall-model ODE to initialize the TLWM field
call tlwm_eq_solve()

! Interpolate/filter wall-stress answer for LES
if (coord == 0) call les_lbc()

end subroutine ic_tlwm_vel

!*****************************************************************************
subroutine tlwm_checkpoint
!*****************************************************************************
!
! This subroutine saves checkpoint variables for MFM analysis
!
use param, only : nzr, write_endian

open(11, file=fname, form='unformatted', convert=write_endian,              &
    status='unknown', position='rewind')
write (11) ur(:,:,1:nzr), vr(:,:,1:nzr), wr(:,:,1:nzr),                     &
    rhs_xr(:,:,1:nzr), rhs_yr(:,:,1:nzr), rhs_zr(:,:,1:nzr)
close(11)

end subroutine tlwm_checkpoint

!*****************************************************************************
subroutine tlwm_eq_ubc
!*****************************************************************************
!
! This subroutine creates the upper boundary condition for the equilibrium
! inner-layer wall-model. This requires grabbing the LES velocity at z = hwm
! on a lower coord, then moving the data to coord = nproc-1 and 
! interpolating/filtering the data to the TLWM grid.
! 
! The equilibrium wall-model only requires the u, v velocity from the LES.
! 
use param, only : nx, ny, jz_r, nzr
use param, only : coord, nproc, comm, ierr, MPI_RPREC
use param, only : nxr, nyr
use sim_param, only : u, v
use mpi
implicit none
real(rprec), dimension(nx,ny) :: dummy_in, u_les, v_les
real(rprec), dimension(nxr,nyr) :: dummy_out

if (coord == 0) then
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

end subroutine tlwm_eq_ubc

!*****************************************************************************
subroutine tlwm_noneq_ubc
!*****************************************************************************
!
! This subroutine follows a similar procedure to tlwm_eq_ubc, however moves
! additional flow variables around because more information is needed for the
! non-equilibrium inner-layer wall-model.
!
use param, only : nx, ny, jz_r, nzr
use param, only : coord, nproc, comm, ierr, MPI_RPREC
use param, only : nxr, nyr
use sim_param, only : u, v, w
use mpi
implicit none
real(rprec), dimension(nx,ny) :: dummy_in, u_les, v_les, w_les
real(rprec), dimension(nxr,nyr) :: dummy_out

if (coord == 0) then
    u_les = u(1:nx,1:ny,jz_r)
    v_les = v(1:nx,1:ny,jz_r)
    ! Interpolate w onto uv-grid at z=hwm
    w_les = 0.5_rprec*(w(1:nx,1:ny,jz_r+1)+w(1:nx,1:ny,jz_r))
else
    u_les = 0.0_rprec
    v_les = 0.0_rprec
    w_les = 0.0_rprec
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

call mpi_allreduce(w_les, dummy_in, nx*ny, mpi_rprec,                      &
    MPI_SUM, comm, ierr)
if (coord == nproc-1) then
    call interp_les_to_tlwm(dummy_in,dummy_out)
    wr(1:nxr,1:nyr,nzr) = dummy_out
endif

end subroutine tlwm_noneq_ubc

!*****************************************************************************
subroutine les_lbc
!*****************************************************************************
!
! This subroutine creates the lower boundary condition for the LES using the
! solution of the inner-layer wall-model. This requires computing the 
! derivatives and wall-stress at z=0, then interpolating/filtering the data
! to the LES grid
!
! This routine should only be accessed by coord=0
!
use param, only : nx, ny, ubot, nu_molec, z_i, u_star
use param, only : nxr, nyr,dzr
use sim_param, only : dudz, dvdz, txz, tyz
use mpi
implicit none
integer :: jx, jy
real(rprec), dimension(nxr,nyr) :: dummy_in
real(rprec), dimension(nx,ny) :: dummy_out

do jy = 1, nyr
do jx = 1, nxr
    ! This is at it is in wallstress.f90/ws_dns_lbc
!    dudzr(jx,jy,1) = ( ur(jx,jy,1) - ubot ) / (0.5_rprec*dzr)
!    dvdzr(jx,jy,1) = vr(jx,jy,1) / (0.5_rprec*dzr) !! vbot = 0
    dudzr(jx,jy,1) = ( ur(jx,jy,1) - ubot ) / zuvr(1)
    dvdzr(jx,jy,1) = vr(jx,jy,1) / zuvr(1) !! vbot = 0
    txzr(jx,jy,1) = -nu_molec/(z_i*u_star)*dudzr(jx,jy,1)
    tyzr(jx,jy,1) = -nu_molec/(z_i*u_star)*dvdzr(jx,jy,1)
enddo
enddo

! Interpolate/filter the data to LES
dummy_in = dudzr(1:nxr,1:nyr,1)
call interp_tlwm_to_les(dummy_in,dummy_out)
dudz(1:nx,1:ny,1) = dummy_out

dummy_in = dvdzr(1:nxr,1:nyr,1)
call interp_tlwm_to_les(dummy_in,dummy_out)
dvdz(1:nx,1:ny,1) = dummy_out

dummy_in = txzr(1:nxr,1:nyr,1)
call interp_tlwm_to_les(dummy_in,dummy_out)
txz(1:nx,1:ny,1) = dummy_out

dummy_in = tyzr(1:nxr,1:nyr,1)
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
use param, only : coord, nproc
#ifdef PPSAFETYMODE
use param, only : BOGUS
#endif
implicit none
integer :: jx, jy, jz, jz_min, jz_max
real(rprec), dimension(nxr, nyr) :: ustar
real(rprec) :: const1
real(rprec), dimension(nxr,nyr,0:nzr) :: a, b, c
real(rprec), dimension(nxr+2,nyr,0:nzr) :: Rx, usol, Ry, vsol

!const1 = 1._rprec/(dzr**2)
const1 = 1._rprec/dzr
! Commented out old code for uniform grid which used const1=1/(dzr**2)
! Now using stretched grid with Jacobians

! Find constant ustar for wall model eddy viscosity
! 1. Take the average of previous timesteps tau value
if (jt_total > 1) then
    ustar = sqrt(abs(txzr(1:nxr,1:nyr,1)) + abs(tyzr(1:nxr,1:nyr,1)))
else !! first time-step, should use nonzero ustar
    ustar = u_star
end if
! 2. Use constant value specified by user - only for debugging
!ustar = u_star 

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
        !! Coefficients are different on bottom row because using one-sided
        !! finite difference at wall, also nu_t(z=0)=0
        !b(jx,jy,1) = -const1*((nu_r(jx,jy,2)+nu_molec)+(nu_molec/0.5_rprec))
        !c(jx,jy,1) = const1*(nu_r(jx,jy,2)+nu_molec)
        !Rx(jx,jy,1) = Rx(jx,jy,1) - const1*(nu_molec/0.5_rprec)*ubot
        !Ry(jx,jy,1) = Ry(jx,jy,1) !! vbot = 0

        b(jx,jy,1) = -const1*(1._rprec/jaco_uvr(1))*                        &
            (const1*(nu_r(jx,jy,2)+nu_molec)/jaco_wr(2) + nu_molec/zuvr(1))
        c(jx,jy,1) = const1*(1._rprec/jaco_uvr(1))*                         &
            const1*(nu_r(jx,jy,2)+nu_molec)/jaco_wr(2)

        Rx(jx,jy,1) = Rx(jx,jy,1) - const1*(1._rprec/jaco_uvr(1))*          &
            (nu_molec/zuvr(1))*ubot
        Ry(jx,jy,1) = Ry(jx,jy,1) !! vbot = 0
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
        !a(jx,jy,nzr-1) = const1*(nu_r(jx,jy,nzr-1)+nu_molec)
        !b(jx,jy,nzr-1) = -const1*(nu_r(jx,jy,nzr)+nu_r(jx,jy,nzr-1)+2._rprec*nu_molec)
        !! interpolated/filtered LES velocity as upper BC
        !Rx(jx,jy,nzr-1) = Rx(jx,jy,nzr-1) - const1*(nu_r(jx,jy,nzr)+nu_molec)*ur(jx,jy,nzr)
        !Ry(jx,jy,nzr-1) = Ry(jx,jy,nzr-1) - const1*(nu_r(jx,jy,nzr)+nu_molec)*vr(jx,jy,nzr)

        a(jx,jy,nzr-1) = const1*(1._rprec/jaco_uvr(nzr-1))*                 &
            const1*(nu_r(jx,jy,nzr-1)+nu_molec)/jaco_wr(nzr-1)
        b(jx,jy,nzr-1) = -const1*(1._rprec/jaco_uvr(nzr-1))*                &
            (const1*(nu_r(jx,jy,nzr-1)+nu_molec)/jaco_wr(nzr-1) +           &
            (nu_r(jx,jy,nzr)+nu_molec)/(L_zr-zwr(nzr)))

        !! interpolated/filtered LES velocity as upper BC
        Rx(jx,jy,nzr-1) = Rx(jx,jy,nzr-1) - const1*(1._rprec/jaco_uvr(nzr-1))*   &
            ((nu_r(jx,jy,nzr)+nu_molec)/(L_zr-zwr(nzr)))*ur(jx,jy,nzr)
        Rx(jx,jy,nzr-1) = Rx(jx,jy,nzr-1) - const1*(1._rprec/jaco_uvr(nzr-1))*   &
            ((nu_r(jx,jy,nzr)+nu_molec)/(L_zr-zwr(nzr)))*vr(jx,jy,nzr)
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
    !a(jx,jy,jz) = const1*(nu_r(jx,jy,jz)+nu_molec)
    !b(jx,jy,jz) = -const1*(nu_r(jx,jy,jz+1)+nu_r(jx,jy,jz)+2._rprec*nu_molec)
    !c(jx,jy,jz) = const1*(nu_r(jx,jy,jz+1)+nu_molec)

    a(jx,jy,jz) = const1*(1._rprec/jaco_uvr(jz))*                         &
        const1*(nu_r(jx,jy,jz)+nu_molec)/jaco_wr(jz)
    b(jx,jy,jz) = -const1*(1._rprec/jaco_uvr(jz))*                        &
        (const1*(nu_r(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1) +               &
        const1*(nu_r(jx,jy,jz)+nu_molec)/jaco_wr(jz))
    c(jx,jy,jz) = const1*(1._rprec/jaco_uvr(jz))*                         &
        const1*(nu_r(jx,jy,jz+1)+nu_molec)/jaco_wr(jz+1)
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
real(rprec), dimension(nxr,nyr), intent(in) :: u1
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
! domain.
! 
! The wall-model coefficient kappa can be either:
!     1. constant, kappa = vonk
!     2. wall-normal constant, sgs dynamic
!     3. wall-normal dynamic, sgs dynamic
!
! The second method was used in the study:
! M. Wang and P. Moin, "Dynamic wall modeling for large-eddy simulation of 
! complex turbulent flows," Phys. Fluids 14, 2043-2051 (2012).
! 
! The third method was used in the study:
! S. Kawai and J. Larsson, "Dynamic non-equilibrium wall-modeling for large
! eddy simulation at high Reynolds numbers," Phys. Fluids 25, (2013).
!  
! Currently methods 2 and 3 are commented out.
!

use types, only : rprec
use param, only : nu_molec, nxr, nyr, nzr, vonk, lbz
! use param, only : dx, dy, ihwm
!use sgs_param, only : Nu_t
!use functions, only : y_avg
implicit none
real(rprec), dimension(nxr, nyr), intent(in) :: ustar
real(rprec) :: a_plus
!real(rprec) :: alpha, max_dx, zcr, bigK
!real(rprec), dimension(nx) :: kwm_temp
!real(rprec), dimension(nx, nzr) :: kwm
real(rprec), dimension(nxr, nyr, lbz:nzr) :: decay
integer :: i, j, k

! Model parameters
a_plus = 17.0_rprec !! used in decay portion
!alpha = 0.48_rprec !! proportionality constant for kwm in ktype = 3
!max_dx = max(dx,dy)
!zcr = alpha*max_dx

! Compute decay portion of eddy viscosity
do i = 1, nxr
do j = 1, nyr
do k = lbz, nzr
    decay(i,j,k) = (1.0_rprec - exp(-zwr(k)*ustar(i,j)/(a_plus*nu_molec)))**2
end do
end do
end do

! Compute wall-model coefficient kwm
!kwm = vonk

!if ((ktype == 1) .or. ((ktype == 3).and.(hwm <= zcr)) )then
!    kwm = vonk
!    ! the second condition is to ensure that bigK > 0 for ktype == 3
!elseif (ktype == 2) then
!
!    kwm_temp = y_avg(Nu_t(1:nx,1:ny,ihwm))/y_avg(ustar(:,:)*hwm*decay(:,:,1))
!    do k = 1, nzr
!        kwm(:,k) = kwm_temp
!    end do
!
!elseif (ktype == 3) then
!
!    kwm_temp = y_avg(Nu_t(1:nx,1:ny,ihwm))/y_avg(ustar(:,:)*hwm*decay(:,:,1))
!    do k = 1, nzr
!        bigK = min( 1.0_rprec, (hwm-zr(k))/(hwm-zcr) )
!        kwm(:,k) = vonk*bigK + kwm_temp*(1.0_rprec - bigK)
!    end do

!end if

! Compute eddy viscosity
do i = 1, nxr
do j = 1, nyr
do k = lbz, nzr
!    nu_r(i,j,k) = kwm(i,k)*ustar(i,j)*zwr(k)*decay(i,j,k)
    nu_r(i,j,k) = vonk*ustar(i,j)*zwr(k)*decay(i,j,k)
end do
end do
end do

end subroutine wm_eddyvisc

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

end module tlwmles
