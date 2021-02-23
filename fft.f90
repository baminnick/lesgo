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

!******************************************************************************
! fftw 3.X version
!******************************************************************************
module fft
use types, only : rprec
use param, only : ld, lh, ny, ld_big, ny2
use iso_c_binding
implicit none
include 'fftw3.f'
save

public :: padd, unpadd, init_fft

public :: kx, ky, k2
public :: forw, back, forw_big, back_big
public :: forw_x, back_x, forw_y, forw_x_fourier
public :: forw_fourier, back_fourier
public :: ycomp_forw_big, ycomp_back_big
public :: ycomp_forw, ycomp_back
#ifdef PPTLWMLES
public :: tlwm_forw, tlwm_back, tlwm_forw_big, tlwm_back_big
#endif

real(rprec), allocatable, dimension(:,:) :: kx, ky, k2
#ifdef PPTLWMLES
real(rprec), allocatable, dimension(:,:) :: kxr, kyr, k2r
#endif

integer*8 :: forw, back, forw_big, back_big
integer*8 :: forw_x, back_x, forw_y, forw_x_fourier
integer*8 :: forw_fourier, back_fourier
integer*8 :: ycomp_forw_big, ycomp_back_big
integer*8 :: ycomp_forw, ycomp_back
#ifdef PPTLWMLES
integer*8 :: tlwm_forw, tlwm_back, tlwm_forw_big, tlwm_back_big
#endif

real (rprec), dimension (:, :), allocatable :: data, data_big

real (rprec), dimension (:), allocatable :: data_x_in, data_y_in, data_x_fourier_in
complex (rprec), dimension (:), allocatable :: data_x_out, data_y_out, data_x_fourier_out

real (rprec), dimension (:, :), allocatable :: data_fourier

complex (rprec), dimension(:), allocatable :: ycomp_data_big

complex (rprec), dimension(:), allocatable :: ycomp_data

#ifdef PPTLWMLES
real(rprec), dimension(:,:), allocatable :: tlwm_data, tlwm_data_big
#endif

contains

!******************************************************************************
subroutine padd (u_big,u)
!******************************************************************************
! puts arrays into larger, zero-padded arrays
! automatically zeroes the oddballs
use types, only : rprec
use param, only : ld,ld_big,nx,ny,ny2
implicit none

!  u and u_big are interleaved as complex arrays
real(rprec), dimension(ld,ny), intent(in) :: u
real(rprec), dimension(ld_big,ny2), intent(out) :: u_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

! make sure the big array is zeroed!
u_big(:,:) = 0._rprec

! note: split access in an attempt to maintain locality
u_big(:nx,:ny_h) = u(:nx,:ny_h)

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

u_big(:nx,j_big_s:ny2) = u(:nx,j_s:ny)

end subroutine padd

!******************************************************************************
subroutine unpadd(cc,cc_big)
!******************************************************************************
use types, only : rprec
use param, only : ld,nx,ny,ny2,ld_big
implicit none

!  cc and cc_big are interleaved as complex arrays
real(rprec), dimension( ld, ny ) :: cc
real(rprec), dimension( ld_big, ny2 ) :: cc_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

cc(:nx,:ny_h) = cc_big(:nx,:ny_h)

! oddballs
cc(ld-1:ld,:) = 0._rprec
cc(:,ny_h+1) = 0._rprec

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2
cc(:nx,j_s:ny) = cc_big(:nx,j_big_s:ny2)

end subroutine unpadd

!******************************************************************************
subroutine init_fft()
!******************************************************************************
use param, only : nx, ny, nx2, ny2, nxp
#ifdef PPTLWMLES
use param, only : nxr, nyr
#endif
implicit none

! Allocate temporary arrays for creating the FFTW plans
allocate( data(ld, ny) )
allocate( data_big(ld_big, ny2) )

allocate( data_x_in(nx) )
allocate( data_x_out(nx/2+1) )
allocate( data_y_in(ny) )
allocate( data_y_out(ny/2+1) )
allocate( data_x_fourier_in(nxp) )
allocate( data_x_fourier_out(nxp/2+1) )

allocate( data_fourier(nxp+2, ny) )

allocate( ycomp_data_big(ny2) ) !! complex stored as real

allocate( ycomp_data(ny) ) !! complex stored as real

#ifdef PPTLWMLES
allocate( tlwm_data(nxr+2,nyr) )
allocate( tlwm_data_big(3*nxr/2+2, 3*nyr/2) )
#endif

! Create the forward and backward plans for the unpadded and padded
! domains. Notice we are using FFTW_UNALIGNED since the arrays used will not be
! guaranteed to be memory aligned.
call dfftw_plan_dft_r2c_2d(forw, nx, ny, data,                               &
    data, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back, nx, ny, data,                               &
    data, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_2d(forw_big, nx2, ny2, data_big,                     &
    data_big, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back_big, nx2, ny2, data_big,                     &
    data_big, FFTW_PATIENT, FFTW_UNALIGNED)

call dfftw_plan_dft_r2c_1d(forw_x, nx, data_x_in,                            &
    data_x_out, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(back_x, nx, data_x_out,                           &
    data_x_in, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_y, ny, data_y_in,                            &
    data_y_out, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_x_fourier, nxp, data_x_fourier_in,           &
    data_x_fourier_out, FFTW_PATIENT, FFTW_UNALIGNED)

call dfftw_plan_dft_r2c_2d(forw_fourier, nxp, ny, data_fourier,              &
    data_fourier, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back_fourier, nxp, ny, data_fourier,              &
    data_fourier, FFTW_PATIENT, FFTW_UNALIGNED)

call dfftw_plan_dft_1d(ycomp_forw_big, ny2, ycomp_data_big,                  &
    ycomp_data_big, FFTW_FORWARD, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_1d(ycomp_back_big, ny2, ycomp_data_big,                  &
    ycomp_data_big, FFTW_BACKWARD, FFTW_PATIENT, FFTW_UNALIGNED)

call dfftw_plan_dft_1d(ycomp_forw, ny, ycomp_data,                           &
    ycomp_data, FFTW_FORWARD, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_1d(ycomp_back, ny, ycomp_data,                           &
    ycomp_data, FFTW_BACKWARD, FFTW_PATIENT, FFTW_UNALIGNED)

#ifdef PPTLWMLES
call dfftw_plan_dft_r2c_2d(tlwm_forw, nxr, nyr, tlwm_data,                   &
    tlwm_data, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(tlwm_back, nxr, nyr, tlwm_data,                   &
    tlwm_data, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_2d(tlwm_forw_big, 3*nxr/2, 3*nyr/2, tlwm_data_big,   &
    tlwm_data_big, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(tlwm_back_big, 3*nxr/2, 3*nyr/2, tlwm_data_big,   &
    tlwm_data_big, FFTW_PATIENT, FFTW_UNALIGNED)
#endif

deallocate(data)
deallocate(data_big)

deallocate(data_x_in)
deallocate(data_x_out)
deallocate(data_y_in)
deallocate(data_y_out)
deallocate(data_x_fourier_in)
deallocate(data_x_fourier_out)


deallocate(data_fourier)

deallocate(ycomp_data_big)

deallocate(ycomp_data)

#ifdef PPTLWMLES
deallocate(tlwm_data)
deallocate(tlwm_data_big)
#endif

call init_wavenumber()
end subroutine init_fft

!******************************************************************************
subroutine init_wavenumber()
!******************************************************************************
use param, only : lh, ny, L_x, L_y, pi
use param, only : kxs_in, fourier, kx_num, coord
#ifdef PPTLWMLES
use param, only : nxr, nyr
#endif
implicit none
integer :: jx, jy

! Allocate wavenumbers
allocate( kx(lh,ny), ky(lh,ny), k2(lh,ny) )
#ifdef PPTLWMLES
allocate( kxr(nxr/2+1,nyr), kyr(nxr/2+1,nyr), k2r(nxr/2+1,nyr) )
#endif

do jx = 1, lh-1
    kx(jx,:) = real(jx-1,kind=rprec)
end do

if (fourier) then
    do jx = 1, kx_num
        kx(jx,:) = kxs_in(jx)
    end do
endif

do jy = 1, ny
    ky(:,jy) = real(modulo(jy - 1 + ny/2,ny) - ny/2,kind=rprec)
end do

! Nyquist: makes doing derivatives easier
kx(lh,:) = 0._rprec
ky(lh,:) = 0._rprec
kx(:,ny/2+1) = 0._rprec
ky(:,ny/2+1) = 0._rprec

! for the aspect ratio change
kx = 2._rprec*pi/L_x*kx
ky = 2._rprec*pi/L_y*ky

! magnitude squared: will have 0's around the edge
k2 = kx*kx + ky*ky

#ifdef PPTLWMLES
do jx = 1, nxr/2
    kxr(jx,:) = real(jx-1,kind=rprec)
end do

do jy = 1, nyr
    kyr(:,jy) = real(modulo(jy - 1 + nyr/2,nyr) - nyr/2,kind=rprec)
end do

! Nyquist: makes doing derivatives easier
kxr(nxr/2+1,:) = 0._rprec
kyr(nxr/2+1,:) = 0._rprec
kxr(:,nyr/2+1) = 0._rprec
kyr(:,nyr/2+1) = 0._rprec

! for the aspect ratio change
kxr = 2._rprec*pi/L_x*kxr
kyr = 2._rprec*pi/L_y*kyr

! magnitude squared: will have 0's around the edge
k2r = kxr*kxr + kyr*kyr
#endif

if ((coord == 0) .and. fourier) then
    write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    write(*,*) 'SIMULATING IN FOURIER SPACE'
    WRITE(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    write(*,*) 'SOLVING ON KX MODES:'
    do jx = 1, lh-1
        write(*,*) kx(jx,1)
    enddo
endif

end subroutine init_wavenumber

end module fft
