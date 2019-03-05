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
use sim_param, only : zhyb, coord_int
use fft
#ifdef PPMAPPING
use sim_param, only : JACO1, JACO2
#endif

implicit none

real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz
integer :: ir, ii
integer :: jz_min
integer :: jj, kcnt

real(rprec), save, dimension(:,:,:), allocatable :: rH_x, rH_y, rH_z
real(rprec), save, dimension(:,:), allocatable :: rtopw, rbottomw
real(rprec), save, dimension(:,:,:), allocatable :: RHS_col
real(rprec), save, dimension(:,:,:), allocatable :: a, b, c
real(rprec), save, dimension(:,:,:), allocatable :: p_rnl, p_phys

logical, save :: arrays_allocated = .false.

real(rprec), dimension(2) :: aH_x, aH_y

real(rprec), dimension(1:ld,1:ny) :: ptemp !! for computing dpdz

! Specifiy cached constants
const = 1._rprec/(nx*ny)
const2 = 1._rprec/tadv1/dt
const3 = 1._rprec/(dz**2)
const4 = 1._rprec/(dz)
! Fourier mode uses const = 1._rprec, Physical uses const = 1._rprec/(nx*ny)
! Both modes use const2 = const/tadv1/dt which is different

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( rH_x(ld,ny,lbz:nz), rH_y(ld,ny,lbz:nz), rH_z(ld,ny,lbz:nz) )
    allocate ( rtopw(ld,ny), rbottomw(ld,ny) )
    allocate ( RHS_col(ld,ny,nz+1) )
    allocate ( a(lh,ny,nz+1), b(lh,ny,nz+1), c(lh,ny,nz+1) )

    allocate ( p_rnl(2*kx_num,ny,0:nz) )
    allocate ( p_phys(ld-2*kx_num,ny,0:nz) )

    ! zeroing for debugging purposes
    !a(:,:,:) = 0.0_rprec
    !b(:,:,:) = 0.0_rprec
    !c(:,:,:) = 0.0_rprec
    !RHS_col(:,:,:) = 0.0_rprec

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

    if (.not. zhyb(jz)) then
        rH_x(:,:,jz) = const * rH_x(:,:,jz)
        rH_y(:,:,jz) = const * rH_y(:,:,jz)
        rH_z(:,:,jz) = const * rH_z(:,:,jz)

        call dfftw_execute_dft_r2c(forw, rH_x(:,:,jz), rH_x(:,:,jz))
        call dfftw_execute_dft_r2c(forw, rH_y(:,:,jz), rH_y(:,:,jz))
        call dfftw_execute_dft_r2c(forw, rH_z(:,:,jz), rH_z(:,:,jz))
    endif
end do

#if defined(PPMPI) && defined(PPSAFETYMODE)
  !Careful - only update real values (odd indicies)
  rH_x(1:ld:2,:,0) = BOGUS
  rH_y(1:ld:2,:,0) = BOGUS
  rH_z(1:ld:2,:,0) = BOGUS
#endif

#ifdef PPSAFETYMODE
!Careful - only update real values (odd indicies)
rH_x(1:ld:2,:,nz) = BOGUS
rH_y(1:ld:2,:,nz) = BOGUS
#endif

#ifdef PPMPI
if (coord == nproc-1) then
    rH_z(:,:,nz) = const * const2 * w(:,:,nz)
    call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,nz))
#ifdef PPSAFETYMODE
else
    rH_z(1:ld:2,:,nz) = BOGUS !--perhaps this should be 0 on top process?
#endif
endif
#else
rH_z(:,:,nz) = const * const2 * w(:,:,nz)
call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,nz))
endif
#endif

if (coord == 0) then
    rbottomw(:,:) = divtz(:,:,1)
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    rtopw(:,:) = const * divtz(:,:,nz)
    call dfftw_execute_dft_r2c(forw, rtopw, rtopw)
#ifdef PPMPI
endif
#endif

! set oddballs to 0
rH_x(ld-1:ld,:,1:nz-1) = 0._rprec
rH_y(ld-1:ld,:,1:nz-1) = 0._rprec
rH_z(ld-1:ld,:,1:nz-1) = 0._rprec
rH_x(:,ny/2+1,1:nz-1) = 0._rprec
rH_y(:,ny/2+1,1:nz-1) = 0._rprec
rH_z(:,ny/2+1,1:nz-1) = 0._rprec
! should also set to zero for rH_z (nz) on coord == nproc-1
if (coord == nproc-1) then
    rH_z(ld-1:ld,:,nz) = 0._rprec
    rH_z(:,ny/2+1,nz) = 0._rprec
end if

! with MPI; topw and bottomw are only on top & bottom processes
rtopw(ld-1:ld, :) = 0._rprec
rtopw(:, ny/2+1) = 0._rprec
rbottomw(ld-1:ld, :) = 0._rprec
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
    RHS_col(:,:,1) = - JACO1(1)*dz* rbottomw(:,:)
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
    RHS_col(:,:,nz+1) = -JACO1(nz)*dz * rtopw(:,:)
#else
    RHS_col(:,:,nz+1) = -dz * rtopw(:,:)
#endif
#ifdef PPMPI
endif
#endif

#ifdef PPMPI
    call mpi_sendrecv (rH_x(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,              &
        rH_x(1, 1, 0), ld*ny, MPI_RPREC, down, 1, comm, status, ierr)
    call mpi_sendrecv (rH_y(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,              &
        rH_y(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,              &
        rH_z(1, 1, 0), ld*ny, MPI_RPREC, down, 3, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, 1), ld*ny, MPI_RPREC, down, 6,               &
        rH_z(1, 1, nz), ld*ny, MPI_RPREC, up, 6, comm, status, ierr)
#endif

do jz = jz_min, nz

if (zhyb(jz-1)) then !! Fourier levels, only include RNL wavenumbers
do jy = 1, ny
    if (jy == ny/2 + 1) cycle

    do jx = 1, kx_num

        ! indices for a, b, c, kx, and ky
        jj = int(kxs_in(jx)) + 1

        if (jj*jy == 1) cycle !! skip kx = ky = 0 

        ! indices for rH_x, rH_y, rH_z, and RHS_col
        ii = 2*jj   ! imaginary index
        ir = ii - 1 ! real index

        ! JDA dissertation, eqn(2.85) a,b,c=coefficients and RHS_col=r_m
#ifdef PPMAPPING
        a(jj, jy, jz) = const3*(1._rprec/(JACO2(jz-1)))*(1._rprec/(JACO1(jz-1)))
        b(jj, jy, jz) = -(kx(jj,jy)**2 + ky(jj,jy)**2                          &
            + const3*(1._rprec/(JACO2(jz-1)))*                                 &
            (1._rprec/(JACO1(jz-1))+1._rprec/(JACO1(jz))))
        c(jj, jy, jz) = const3*(1._rprec/(JACO2(jz-1)))*(1._rprec/(JACO1(jz)))
#else
        a(jj, jy, jz) = const3
        b(jj, jy, jz) = -(kx(jj, jy)**2 + ky(jj, jy)**2 + 2._rprec*const3)
        c(jj, jy, jz) = const3
#endif

        !  Compute eye * kx * H_x
        aH_x(1) = -rH_x(ii,jy,jz-1) * kx(jj,jy)
        aH_x(2) =  rH_x(ir,jy,jz-1) * kx(jj,jy)
        aH_y(1) = -rH_y(ii,jy,jz-1) * ky(jj,jy)
        aH_y(2) =  rH_y(ir,jy,jz-1) * ky(jj,jy)

#ifdef PPMAPPING
        RHS_col(ir:ii,jy,jz) = aH_x + aH_y + (rH_z(ir:ii, jy, jz) -            &
            rH_z(ir:ii, jy, jz-1))*const4/JACO2(jz-1)
#else
        RHS_col(ir:ii,jy,jz) =  aH_x + aH_y + (rH_z(ir:ii, jy, jz) -           &
            rH_z(ir:ii, jy, jz-1))*const4
#endif

    end do
end do

else !! Physical levels, include all wavenumbers

do jy = 1, ny
    if (jy == ny/2 + 1) cycle

    do jx = 1, lh-1

        if (jx*jy == 1) cycle !! skip kx = ky = 0 

        ! indices for rH_x, rH_y, rH_z, and RHS_col
        ii = 2*jx   ! imaginary index
        ir = ii - 1 ! real index

        ! JDA dissertation, eqn(2.85) a,b,c=coefficients and RHS_col=r_m
#ifdef PPMAPPING
        a(jx, jy, jz) = const3*(1._rprec/(JACO2(jz-1)))*(1._rprec/(JACO1(jz-1)))
        b(jx, jy, jz) = -(kx(jx,jy)**2 + ky(jx,jy)**2                          &
            + const3*(1._rprec/(JACO2(jz-1)))*                                 &
            (1._rprec/(JACO1(jz-1))+1._rprec/(JACO1(jz))))
        c(jx, jy, jz) = const3*(1._rprec/(JACO2(jz-1)))*(1._rprec/(JACO1(jz)))
#else
        a(jx, jy, jz) = const3
        b(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2 + 2._rprec*const3)
        c(jx, jy, jz) = const3
#endif

        !  Compute eye * kx * H_x
        aH_x(1) = -rH_x(ii,jy,jz-1) * kx(jx,jy)
        aH_x(2) =  rH_x(ir,jy,jz-1) * kx(jx,jy)
        aH_y(1) = -rH_y(ii,jy,jz-1) * ky(jx,jy)
        aH_y(2) =  rH_y(ir,jy,jz-1) * ky(jx,jy)

#ifdef PPMAPPING
        RHS_col(ir:ii,jy,jz) = aH_x + aH_y + (rH_z(ir:ii, jy, jz) -            &
            rH_z(ir:ii, jy, jz-1))*const4/JACO2(jz-1)
#else
        RHS_col(ir:ii,jy,jz) =  aH_x + aH_y + (rH_z(ir:ii, jy, jz) -           &
            rH_z(ir:ii, jy, jz-1))*const4
#endif

    end do
end do

endif

end do

! Correct coefficients for wavenumbers NOT included in RNL at first uv point
! above the interface. This is done after values are placed for every wavenumber.
! So only need to give zeros to enforce interface BC. This acts as the pressure 
! boundary condition for the LES. Currently using Dirichlet BC, may use Neumann.
do jz = jz_min, nz
if ((zhyb(jz-1)) .and. (.not. zhyb(jz))) then !! interface at jz-1, want jz

do jy = 1, ny
    if (jy == ny/2 + 1) cycle

    kcnt = 1
    do jx = 1, lh-1

        if (kcnt .le. kx_num) then
        if (jx .eq. (int(kxs_in(kcnt))+1)) then !! skip RNL wavenumbers
            kcnt = kcnt + 1
            cycle
        endif
        endif

        if (jx*jy == 1) cycle !! skip kx = ky = 0 

        a(jx, jy, jz+1) = 0.0_rprec

    end do
end do

endif
end do

! DEBUG
!if (coord == 5) then
!do jz = 1, nz+1
!do jy = 1, ny
!write(*,*) jy, jz, RHS_col(kxi,jy,jz)
!enddo
!enddo
!endif

! tridag_array routines skips zero wavenumber solution, nyquist freqs
! --> Only accounted for in tridag_array_hybrid_rnl, 
!     assuming kx=ky=0 only included in kxs_in

! call mpi_barrier( comm, ierr )

! Find pressure modes for rnl wavenumbers by solving system across all z-levels
call tridag_array_hybrid_rnl (a(int(kxs_in)+1,:,:), b(int(kxs_in)+1,:,:), c(int(kxs_in)+1,:,:), RHS_col(kxi,:,:), p_rnl)

if (coord_int /= -1) then
! Find pressure modes for non-rnl wavenumbers by solving system above wall-model
call tridag_array_hybrid_phys (a(int(kxs_phys)+1,:,:), b(int(kxs_phys)+1,:,:), c(int(kxs_phys)+1,:,:), RHS_col(kxpi,:,:), p_phys)
endif

! zeroing for debugging purposes
p(:,:,:) = 0.0_rprec
dpdx(:,:,:) = 0.0_rprec
dpdy(:,:,:) = 0.0_rprec
dpdz(:,:,:) = 0.0_rprec

! Put pressure values into sim_param pressure
p(kxi,:,:) = p_rnl(:,:,:)
p(kxpi,:,:) = p_phys(:,:,:)

! DEBUG
!if (coord == 6) then
!write(*,*) 'B', p(kxi,1,4)
!endif

! zero-wavenumber solution
#ifdef PPMPI
! wait for p(1, 1, 1) from "down"
call mpi_recv (p(1:2, 1, 1), 2, MPI_RPREC, down, 8, comm, status, ierr)
#endif

if (coord == 0) then
    p(1:2, 1, 0) = 0._rprec !! BC, arbitrary pressure
#ifdef PPMAPPING
    p(1:2, 1, 1) = p(1:2,1,0) - JACO1(1)*dz * rbottomw(1:2,1)
#else
    p(1:2, 1, 1) = p(1:2,1,0) - dz * rbottomw(1:2,1)
#endif
end if

do jz = 2, nz
    ! JDA dissertation, eqn(2.88)
#ifdef PPMAPPING
    p(1:2, 1, jz) = p(1:2, 1, jz-1) + rH_z(1:2, 1, jz) * dz * JACO1(jz)
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
call mpi_sendrecv (p(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,                     &
    p(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
#endif

! zero the nyquist freqs
p(ld-1:ld,:,:) = 0._rprec
p(:,ny/2+1,:) = 0._rprec

! Now need to get p(wave,level) to physical p(jx,jy,jz)
! Loop over height levels
if (.not. zhyb(0)) then
    call dfftw_execute_dft_c2r(back,p(:,:,0), p(:,:,0))
endif

do jz = 1, nz-1
    if (zhyb(jz)) then
        do jy = 1, ny
        do jx = 1, kx_num
            ii = 2*int( kxs_in(jx) ) + 2   ! imaginary index
            ir = ii - 1 ! real index
            dpdx(ir,jy,jz) = -p(ii,jy,jz) * kxrnl(jx,jy)
            dpdx(ii,jy,jz) =  p(ir,jy,jz) * kxrnl(jx,jy)
            dpdy(ir,jy,jz) = -p(ii,jy,jz) * kyrnl(jx,jy)
            dpdy(ii,jy,jz) =  p(ir,jy,jz) * kyrnl(jx,jy)
        end do
        end do
    else
        do jy = 1, ny
        do jx = 1, lh-1
            ii = 2*jx   ! imaginary index
            ir = ii - 1 ! real index
            dpdx(ir,jy,jz) = -p(ii,jy,jz) * kx(jx,jy)
            dpdx(ii,jy,jz) =  p(ir,jy,jz) * kx(jx,jy)
            dpdy(ir,jy,jz) = -p(ii,jy,jz) * ky(jx,jy)
            dpdy(ii,jy,jz) =  p(ir,jy,jz) * ky(jx,jy)
        end do
        end do
    endif

    ! note the oddballs of p are already 0, so we should be OK here
    if (.not. zhyb(jz)) then
        call dfftw_execute_dft_c2r(back,dpdx(:,:,jz), dpdx(:,:,jz))
        call dfftw_execute_dft_c2r(back,dpdy(:,:,jz), dpdy(:,:,jz))
        call dfftw_execute_dft_c2r(back,p(:,:,jz), p(:,:,jz))
    endif
end do

if (coord==nproc-1) then
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
    if (zhyb(jz) .and. zhyb(jz-1)) then !! both Fourier levels

#ifdef PPMAPPING
        dpdz(kxi, 1:ny, jz) = (p(kxi, 1:ny, jz) - p(kxi, 1:ny, jz-1))        &
            / dz / JACO1(jz)
#else
        dpdz(kxi, 1:ny, jz) = (p(kxi, 1:ny, jz) - p(kxi, 1:ny, jz-1)) / dz
#endif

    elseif ( (.not. zhyb(jz)) .and. (.not. zhyb(jz-1)) ) then !! both Physical

#ifdef PPMAPPING
        dpdz(1:nx, 1:ny, jz) = (p(1:nx, 1:ny, jz) - p(1:nx, 1:ny, jz-1))         &
            / dz / JACO1(jz)
#else
        dpdz(1:nx, 1:ny, jz) = (p(1:nx, 1:ny, jz) - p(1:nx, 1:ny, jz-1)) / dz
#endif

    else !! jz-1 is at interface in Fourier, jz in Physical, dpdz in Physical

        ptemp(:,:) = p(:,:,jz-1)
        call dfftw_execute_dft_c2r(back, ptemp(:,:), ptemp(:,:))

#ifdef PPMAPPING
        dpdz(1:nx, 1:ny, jz) = (p(1:nx, 1:ny, jz) - ptemp(1:nx, 1:ny))         &
            / dz / JACO1(jz)
#else
        dpdz(1:nx, 1:ny, jz) = (p(1:nx, 1:ny, jz) - ptemp(1:nx, 1:ny)) / dz
#endif

    end if
end do


#ifdef PPSAFETYMODE
if(coord<nproc-1)  dpdz(:,:,nz) = BOGUS
#endif
#ifdef PPMAPPING
if(coord==nproc-1) dpdz(1:nx,1:ny,nz) = (p(1:nx,1:ny,nz)-p(1:nx,1:ny,nz-1))  &
    / dz / JACO1(nz)
#else
if(coord==nproc-1) dpdz(1:nx,1:ny,nz) = (p(1:nx,1:ny,nz)-p(1:nx,1:ny,nz-1)) / dz
#endif

! DEBUG
!if (coord == 0) then
!write(*,*) 'C', dpdy(kxi,2,4)
!write(*,*) '---------------------------------'
!endif

end subroutine press_stag_array_hybrid
