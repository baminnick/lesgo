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
subroutine press_stag_array_hybrid()
!*******************************************************************************
!
! Calculate the pressure and its derivatives on exit. Everything is in physical
! space on exit.
!
use types, only : rprec
use param
use messages
use sim_param, only : u, v, w, divtz, p, dpdx, dpdy, dpdz
use fft
use derivatives, only : mpi_sync_hybrid
use mpi_defs, only : MPI_SYNC_UP
#ifdef PPMAPPING
use sim_param, only : jaco_w, jaco_uv
#endif

implicit none

real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz
integer :: ir, ii
integer :: jz_min
integer :: end_kx
integer :: jj, kk
integer :: lhp, ldp

real(rprec), save, dimension(:,:,:), allocatable :: rH_x, rH_y, rH_z
real(rprec), save, dimension(:,:), allocatable :: rtopw, rbottomw
real(rprec), save, dimension(:,:,:), allocatable :: RHS_col
real(rprec), save, dimension(:,:,:), allocatable :: a, b, c
real(rprec), save, dimension(:,:,:), allocatable :: p_rnl, p_phys

logical, save :: arrays_allocated = .false.

real(rprec), dimension(2) :: aH_x, aH_y

real(rprec), dimension(1:(nxp+2),1:ny) :: ptemp

! Specifiy cached constants
if (fourier) then
    const = 1._rprec
else
    const = 1._rprec/(nxp*ny)
endif
const2 = const/tadv1/dt
const3 = 1._rprec/(dz**2)
const4 = 1._rprec/(dz)
! Fourier mode uses const = 1._rprec, Physical uses const = 1._rprec/(nx*ny)
! Both modes use const2 = const/tadv1/dt which is different

lhp = nxp/2 + 1
ldp = 2*lhp

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( rH_x(ldp,ny,lbz:nz), rH_y(ldp,ny,lbz:nz), rH_z(ldp,ny,lbz:nz) )
    allocate ( rtopw(ldp,ny), rbottomw(ldp,ny) )
    allocate ( RHS_col(ldp,ny,nz+1) )
    allocate ( a(lhp,ny,nz+1), b(lhp,ny,nz+1), c(lhp,ny,nz+1) )

    allocate ( p_rnl(nxf,ny,0:nz) )
    allocate ( p_phys(nxp-nxf,ny,0:nz) )

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
    if (fourier) then !! RNL WM coords
        ! Initialize with all zeros
        rH_x(:,:,jz) = 0.0_rprec
        rH_y(:,:,jz) = 0.0_rprec
        rH_z(:,:,jz) = 0.0_rprec

        ! Add only RNL wavenumbers, ignore Nyquist
        rH_x(kxi,:,jz) = const2 * u(1:nx,:,jz)
        rH_y(kxi,:,jz) = const2 * v(1:nx,:,jz)
        rH_z(kxi,:,jz) = const2 * w(1:nx,:,jz)
    else !! non Fourier coords
        ! Transform 
        rH_x(:,:,jz) = const2 * u(:,:,jz)
        rH_y(:,:,jz) = const2 * v(:,:,jz)
        rH_z(:,:,jz) = const2 * w(:,:,jz)
        call dfftw_execute_dft_r2c(forw, rH_x(:,:,jz), rH_x(:,:,jz))
        call dfftw_execute_dft_r2c(forw, rH_y(:,:,jz), rH_y(:,:,jz))
        call dfftw_execute_dft_r2c(forw, rH_z(:,:,jz), rH_z(:,:,jz))
    endif

end do

#ifdef PPSAFETYMODE
  !Careful - only update real values (odd indicies)
  rH_x(1:ldp:2,:,0) = BOGUS
  rH_y(1:ldp:2,:,0) = BOGUS
  rH_z(1:ldp:2,:,0) = BOGUS
#endif

#ifdef PPSAFETYMODE
!Careful - only update real values (odd indicies)
rH_x(1:ldp:2,:,nz) = BOGUS
rH_y(1:ldp:2,:,nz) = BOGUS
#endif

if (coord == nproc-1) then
    ! Top coord in physical space, i.e. fourier = false
    rH_z(:,:,nz) = const2 * w(:,:,nz)
    call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,nz))
#ifdef PPSAFETYMODE
else
    rH_z(1:ld:2,:,nz) = BOGUS !--perhaps this should be 0 on top process?
#endif
endif

if (coord == 0) then
    ! Bottom coord in fourier space, i.e. fourier = true
    rbottomw(kxi,:) = const * divtz(1:nx,:,1)
end if

if (coord == nproc-1) then
    ! Top coord in physical space, i.e. fourier = false
    rtopw(:,:) = const * divtz(:,:,nz)
    call dfftw_execute_dft_r2c(forw, rtopw, rtopw)
endif

! set oddballs to 0
rH_x(ldp-1:ldp,:,1:nz-1) = 0._rprec
rH_y(ldp-1:ldp,:,1:nz-1) = 0._rprec
rH_z(ldp-1:ldp,:,1:nz-1) = 0._rprec
rH_x(:,ny/2+1,1:nz-1) = 0._rprec
rH_y(:,ny/2+1,1:nz-1) = 0._rprec
rH_z(:,ny/2+1,1:nz-1) = 0._rprec
! should also set to zero for rH_z (nz) on coord == nproc-1
if (coord == nproc-1) then
    rH_z(ldp-1:ldp,:,nz) = 0._rprec
    rH_z(:,ny/2+1,nz) = 0._rprec
end if

! with MPI; topw and bottomw are only on top & bottom processes
rtopw(ldp-1:ldp, :) = 0._rprec
rtopw(:, ny/2+1) = 0._rprec
rbottomw(ldp-1:ldp, :) = 0._rprec
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
    call mpi_sendrecv (rH_x(1, 1, nz-1), ldp*ny, MPI_RPREC, up, 1,          &
        rH_x(1, 1, 0), ldp*ny, MPI_RPREC, down, 1, comm, status, ierr)
    call mpi_sendrecv (rH_y(1, 1, nz-1), ldp*ny, MPI_RPREC, up, 2,          &
        rH_y(1, 1, 0), ldp*ny, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, nz-1), ldp*ny, MPI_RPREC, up, 3,          &
        rH_z(1, 1, 0), ldp*ny, MPI_RPREC, down, 3, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, 1), ldp*ny, MPI_RPREC, down, 6,           &
        rH_z(1, 1, nz), ldp*ny, MPI_RPREC, up, 6, comm, status, ierr)
#endif

! If RNL WM coord, loop through RNL wavenumbers only
! If non Fourier coord, loop through all wavenumbers
if (fourier) then
    end_kx = kx_num-1
else
    end_kx = lh-1
endif

do jz = jz_min, nz
do jy = 1, ny
    if (jy == ny/2 + 1) cycle

    kk = 0 !! index for wavenumbers
    do jx = 1, end_kx

        ! indices for a, b, c
        if (fourier) then
            jj = int(kxs_in(jx)) + 1
        else
            jj = jx
        endif

        kk = kk + 1
        if (jj*jy == 1) cycle !! skip kx = ky = 0 

        ! indices for rH_x, rH_y, rH_z, and RHS_col
        ii = 2*jj   ! imaginary index
        ir = ii - 1 ! real index

        ! JDA dissertation, eqn(2.85) a,b,c=coefficients and RHS_col=r_m
#ifdef PPMAPPING
        a(jj, jy, jz) = const3*(1._rprec/(jaco_uv(jz-1)))*(1._rprec/(jaco_w(jz-1)))
        b(jj, jy, jz) = -(kx(kk,jy)**2 + ky(kk,jy)**2                      &
            + const3*(1._rprec/(jaco_uv(jz-1)))*                           &
            (1._rprec/(jaco_w(jz-1))+1._rprec/(jaco_w(jz))))
        c(jj, jy, jz) = const3*(1._rprec/(jaco_uv(jz-1)))*(1._rprec/(jaco_w(jz)))
#else
        a(jj, jy, jz) = const3
        b(jj, jy, jz) = -(kx(kk, jy)**2 + ky(kk, jy)**2 + 2._rprec*const3)
        c(jj, jy, jz) = const3
#endif

        !  Compute eye * kx * H_x
        aH_x(1) = -rH_x(ii,jy,jz-1) * kx(kk,jy)
        aH_x(2) =  rH_x(ir,jy,jz-1) * kx(kk,jy)
        aH_y(1) = -rH_y(ii,jy,jz-1) * ky(kk,jy)
        aH_y(2) =  rH_y(ir,jy,jz-1) * ky(kk,jy)

#ifdef PPMAPPING
        RHS_col(ir:ii,jy,jz) = aH_x + aH_y + (rH_z(ir:ii, jy, jz) -        &
            rH_z(ir:ii, jy, jz-1))*const4/jaco_uv(jz-1)
#else
        RHS_col(ir:ii,jy,jz) =  aH_x + aH_y + (rH_z(ir:ii, jy, jz) -       &
            rH_z(ir:ii, jy, jz-1))*const4
#endif

    end do
end do
end do

! Correct coefficients for wavenumbers NOT included in RNL near the interface.
! This is done after values are placed for every wavenumber. So only need to 
! give zeroes to enforce interface BC. This acts as the pressure boundary 
! condition for the LES.
if (coord == nproc_rnl) then
do jy = 1, ny
    if (jy == ny/2 + 1) cycle
    jj = 1
    do jx = 1, end_kx
        if (jj .le. (kx_num-1)) then
        if (jx .eq. (int(kxs_in(jj))+1)) then !! skip RNL wavenumbers
            jj = jj + 1
            cycle
        endif
        endif

        if (jx*jy == 1) cycle !! skip kx = ky = 0
        a(jx,jy,1) = 0.0_rprec
    end do
end do
endif

! tridag_array routines skips zero wavenumber solution, nyquist freqs
! --> Only accounted for in tridag_array_hybrid_rnl, 
!     assuming kx=ky=0 already included in kxs_in
! Find pressure modes for rnl wavenumbers by solving system across all z-levels
call tridag_array_hybrid_rnl (a(int(kxs_in)+1,:,:), b(int(kxs_in)+1,:,:), c(int(kxs_in)+1,:,:), RHS_col(kxi,:,:), p_rnl)

if (.not. fourier) then
! Find pressure modes for non-rnl wavenumbers by solving system above wall-model
call tridag_array_hybrid_phys (a(int(kxs_phys)+1,:,:), b(int(kxs_phys)+1,:,:), c(int(kxs_phys)+1,:,:), RHS_col(kxpi,:,:), p_phys)
endif

! Put pressure values into sim_param pressure
if (fourier) then
    ! Only store RNL related wavenumbers
    p(1:nx,:,:) = p_rnl(:,:,:)
else
    ! Store all wavenumbers
    p(kxi,:,:) = p_rnl(:,:,:)
    p(kxpi,:,:) = p_phys(:,:,:)
endif

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
if (coord .ne. (nproc_rnl - 1)) then
    call mpi_sendrecv( p(1,1,nz-1), ld*ny, MPI_RPREC, up, 2,      &
        p(1,1,0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
else !! RNL coord near interface
    ptemp(:,:) = 0.0_rprec !! fill with zeros initially
    ptemp(kxi,:) = p(1:nx,:,nz-1) !! fill RNL wavenumber spots only
    call mpi_sendrecv( ptemp, (nxp+2)*ny, MPI_RPREC, up, 2,       &
        p(:,:,0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
endif
#endif

! zero the nyquist freqs
p(ld-1:ld,:,:) = 0._rprec
p(:,ny/2+1,:) = 0._rprec

! Now need to get p(wave,level) to physical p(jx,jy,jz)
! Loop over height levels
if (.not. fourier) then
    call dfftw_execute_dft_c2r(back,p(:,:,0), p(:,:,0))
endif

do jz = 1, nz-1
    do jy = 1, ny
    do jx = 1, end_kx
        ii = 2*jx   ! imaginary index
        ir = ii - 1 ! real index
        ! Take ddx and ddy derivatives of pressure
        dpdx(ir,jy,jz) = -p(ii,jy,jz) * kx(jx,jy)
        dpdx(ii,jy,jz) =  p(ir,jy,jz) * kx(jx,jy)
        dpdy(ir,jy,jz) = -p(ii,jy,jz) * ky(jx,jy)
        dpdy(ii,jy,jz) =  p(ir,jy,jz) * ky(jx,jy)
    end do
    end do

    ! note the oddballs of p are already 0, so we should be OK here
    if (.not. fourier) then
        call dfftw_execute_dft_c2r(back,dpdx(:,:,jz), dpdx(:,:,jz))
        call dfftw_execute_dft_c2r(back,dpdy(:,:,jz), dpdy(:,:,jz))
        call dfftw_execute_dft_c2r(back,p(:,:,jz), p(:,:,jz))
    endif
end do

if ((coord==nproc-1) .and. (.not. fourier)) then !! top coord, fourier = false
    call dfftw_execute_dft_c2r(back,p(:,:,nz),p(:,:,nz))
endif

! nz level is not needed elsewhere (although its valid)
#ifdef PPSAFETYMODE
dpdx(:,:,nz) = BOGUS
dpdy(:,:,nz) = BOGUS
if(coord<nproc-1) p(:,:,nz) = BOGUS
#endif

! Final step compute the z-derivative of p
! note: p has additional level at z=-dz/2 for this derivative
do jz = 1, nz-1
#ifdef PPMAPPING
    dpdz(1:nx,1:ny,jz) = (p(1:nx,1:ny,jz) - p(1:nx,1:ny,jz-1))      &
        / dz / jaco_w(jz)
#else
    dpdz(1:nx,1:ny,jz) = (p(1:nx,1:ny,jz) - p(1:nx,1:ny,jz-1)) / dz
#endif
end do

#ifdef PPSAFETYMODE
if(coord<nproc-1)  dpdz(:,:,nz) = BOGUS
#endif
#ifdef PPMAPPING
if(coord==nproc-1) dpdz(1:nx,1:ny,nz) = (p(1:nx,1:ny,nz)-p(1:nx,1:ny,nz-1))  &
    / dz / jaco_w(nz)
#else
if(coord==nproc-1) dpdz(1:nx,1:ny,nz) = (p(1:nx,1:ny,nz)-p(1:nx,1:ny,nz-1)) / dz
#endif

end subroutine press_stag_array_hybrid
