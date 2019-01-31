!!
!!  Copyright (C) 2010-2017  Johns Hopkins University
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
module derivatives
!*******************************************************************************
!
! This module contains all of the major subroutines used for computing
! derivatives and taking Fourier transforms.
!
implicit none

save
private

public ddx, ddy, ddxy, filt_da, ddz_uv, ddz_w,                          &
    phys2wave, wave2phys, phys2waveF, wave2physF,                       &
    dft_direct_forw_2d_n_yonlyC_big, dft_direct_back_2d_n_yonlyC_big,   &
    convolve_rnl

contains

!*******************************************************************************
subroutine ddx(f,dfdx,lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz, fourier
use fft
use emul_complex, only : OPERATOR(.MULI.)
!use param, only : coord
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx
real(rprec) :: const
integer :: jz

!real(rprec) :: xx, yy
!real(rprec), dimension(ld,ny) :: FF
!integer :: jx

!if (coord == 0) then
!    write(*,*) '----------------------------------'
!    ! Create FF
!    do jx = 1, nx
!    do jz = 1, ny
!        xx = (2.0_rprec*3.1415926535_rprec)*real(jx-1,rprec)/real(nx-1,rprec)
!        yy = (2.0_rprec*3.1415926535_rprec)*real(jz-1,rprec)/real(ny-1,rprec)
!        ! FF(jx,jz) = cos(5.0_rprec*xx) + cos(2.0_rprec*yy) + 5.0_rprec
!        FF(jx,jz) = cos(5.0_rprec*xx) + cos(2.0_rprec*yy)
!        ! FF(jx,jz) = cos(5.0_rprec*xx) 
!        ! FF(jx,jz) = cos(5.0_rprec*yy)
!        ! FF(jx,jz) = 5.0_rprec
!    enddo
!    enddo
!
!    call dfftw_execute_dft_r2c(forw,FF(:,:),FF(:,:))
!    FF(:,:) = FF(:,:) / ( nx * ny )
!    write(*,*) FF(:,1)
!    write(*,*) '***'
!    write(*,*) FF(:,2)
!    write(*,*) '***'
!    write(*,*) FF(:,3)
!    write(*,*) '***'
!    write(*,*) FF(:,4)
!    write(*,*) '***'
!    write(*,*) FF(:,5)
!    write(*,*) '***'
!    write(*,*) FF(:,6)
!    write(*,*) '***'
!    write(*,*) FF(:,7)
!    write(*,*) '***'
!    write(*,*) FF(:,8)
!    write(*,*) '***'
!    write(*,*) FF(:,9)
!    write(*,*) '***'
!    write(*,*) FF(:,10)
!
!    write(*,*) '----------------------------------'
!endif

const = 1._rprec !! if fourier
if (.not. fourier) const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz = lbz, nz
    !  Use dfdx to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = const*f(:,:,jz)
    if (.not. fourier) then 
        call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz),dfdx(:,:,jz))
    endif

    ! Zero padded region and Nyquist frequency
    dfdx(ld-1:ld,:,jz) = 0._rprec
    dfdx(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdx to perform complex multiplication
    ! Optimized version for real(eye*kx)=0
    ! only passing imaginary part of eye*kx
    dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

    ! Perform inverse transform to get pseudospectral derivative
    if (.not. fourier) then
        call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
    endif
enddo

end subroutine ddx

!*******************************************************************************
subroutine ddy(f,dfdy, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! y using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz, fourier
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
real(rprec) :: const
integer :: jz

const = 1._rprec !! if fourier
if (.not. fourier) const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz = lbz, nz
    !  Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdy(:,:,jz) = const * f(:,:,jz)
    if (.not. fourier) then
        call dfftw_execute_dft_r2c(forw, dfdy(:,:,jz), dfdy(:,:,jz))
    endif

    ! Zero padded region and Nyquist frequency
    dfdy(ld-1:ld,:,jz) = 0._rprec
    dfdy(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdy(:,:,jz) .MULI. ky

    ! Perform inverse transform to get pseudospectral derivative
    if (.not. fourier) then
        call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
    endif
end do

end subroutine ddy

!*******************************************************************************
subroutine ddxy (f, dfdx, dfdy, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x and y using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz, fourier
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
real(rprec) :: const
integer :: jz

const = 1._rprec !! if fourier
if (.not. fourier) const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz = lbz, nz
    ! Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = const*f(:,:,jz)
    if (.not. fourier) then
        call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz), dfdx(:,:,jz))
    endif

    ! Zero padded region and Nyquist frequency
    dfdx(ld-1:ld,:,jz) = 0._rprec
    dfdx(:,ny/2+1,jz) = 0._rprec

    ! Derivatives: must to y's first here, because we're using dfdx as storage
    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdx(:,:,jz) .MULI. ky
    dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

    ! Perform inverse transform to get pseudospectral derivative
    if (.not. fourier) then
        call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
        call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
    endif
end do

end subroutine ddxy

!*******************************************************************************
subroutine filt_da(f,dfdx,dfdy, lbz)
!*******************************************************************************
!
! This subroutine kills the oddball components in f and computes the partial
! derivative of f with respect to x and y using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz, fourier
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(inout) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
real(rprec) :: const
integer :: jz

const = 1._rprec !! if fourier
if (.not. fourier) const = 1._rprec/(nx*ny)

! loop through horizontal slices
do jz = lbz, nz
    ! Calculate FFT in place
    f(:,:,jz) = const*f(:,:,jz)
    if (.not. fourier) then
        call dfftw_execute_dft_r2c(forw, f(:,:,jz), f(:,:,jz))
    endif

    ! Kill oddballs in zero padded region and Nyquist frequency
    f(ld-1:ld,:,jz) = 0._rprec
    f(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
    dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

    ! Perform inverse transform to get pseudospectral derivative
    ! The oddballs for derivatives should already be dead, since they are for f
    ! inverse transform
    if (.not. fourier) then
        call dfftw_execute_dft_c2r(back, f(:,:,jz), f(:,:,jz))
        call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
        call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
    endif
end do

end subroutine filt_da

!*******************************************************************************
subroutine ddz_uv(f, dfdz, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to z using
! 2nd order finite differencing. f is on the uv grid and dfdz is on the w grid.
! The serial version provides dfdz(:,:,2:nz), and the value at jz=1 is not
! touched. The MPI version provides dfdz(:,:,1:nz), except at the bottom
! process it only supplies 2:nz
!
use types, only : rprec
use param, only : nx, ny, nz, dz, BOGUS
#ifdef PPSAFETYMODE
use param, only : nproc, coord
#endif
use sim_param, only : JACO1
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
integer :: jx, jy, jz
real(rprec) :: const

const = 1._rprec/dz

#if defined(PPMPI) && defined(PPSAFETYMODE)
dfdz(:,:,0) = BOGUS
#endif

! Calculate derivative.
! The ghost node information is available here
! if coord == 0, dudz(1) will be set in wallstress
do jz = lbz+1, nz
do jy = 1, ny
do jx = 1, nx !! jb has ld instead of nx fourier
    dfdz(jx,jy,jz) = (1/JACO1(jz))*const*(f(jx,jy,jz)-f(jx,jy,jz-1))
end do
end do
end do

! Not necessarily accurate at top and bottom boundary
! Set to BOGUS just to be safe
#ifdef PPSAFETYMODE
if (coord == 0) then
    dfdz(:,:,1) = BOGUS
end if
if (coord == nproc-1) then
    dfdz(:,:,nz) = BOGUS
end if
#endif

end subroutine ddz_uv

!*******************************************************************************
subroutine ddz_w(f, dfdz, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to z using
! 2nd order finite differencing. f is on the w grid and dfdz is on the uv grid.
! The serial version provides dfdz(:,:,1:nz-1), and the value at jz=1 is not
! touched. The MPI version provides dfdz(:,:,0:nz-1), except at the top and
! bottom processes, which each has has 0:nz, and 1:nz-1, respectively.
!
use types, only : rprec
use param, only : ld, ny, nz, dz, BOGUS
#ifdef PPSAFETYMODE
#ifdef PPMPI
use param, only : coord
#endif
#endif
use sim_param, only : JACO2
implicit none

real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
integer, intent(in) :: lbz
real(rprec)::const
integer :: jx, jy, jz

const = 1._rprec/dz
do jz = lbz, nz-1
do jy = 1, ny
do jx = 1, ld !! nx fourier
    dfdz(jx,jy,jz) = (1/JACO2(jz))*const*(f(jx,jy,jz+1)-f(jx,jy,jz))
end do
end do
end do

#ifdef PPSAFETYMODE
#ifdef PPMPI
! bottom process cannot calculate dfdz(jz=0)
if (coord == 0) then
    dfdz(:,:,lbz) = BOGUS
endif
#endif
! All processes cannot calculate dfdz(jz=nz)
dfdz(:,:,nz) = BOGUS
#endif

end subroutine ddz_w

!*******************************************************************************
subroutine phys2wave( f, lb )
!*******************************************************************************
!
! Transform field from physical space to wavenumber space. 
! (x,y) --> (kx,ky)
!
! lb - lower bound index for z-dimension
!     lb = lbz for everything except pressure gradients (which should be lb = 1)
! 
use types, only : rprec
use param, only : nx, ny, nz
use fft
implicit none

real(rprec), dimension(:,:,lb:), intent(inout) :: f
real(rprec) :: const
integer :: jz
integer, intent(in) :: lb

const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz = lb, nz
    f(:,:,jz) = const*f(:,:,jz) !! normalization
    call dfftw_execute_dft_r2c(forw, f(:,:,jz), f(:,:,jz))
enddo

end subroutine phys2wave

!*******************************************************************************
subroutine wave2phys( f, lb )
!*******************************************************************************
!
! Transform field from wavenumber space to physical space. 
! (kx,ky) --> (x,y)
!
! lb - lower bound index for z-dimension
!     lb = lbz for everything except pressure gradients (which should be lb = 1)
! 
use types, only : rprec
use param, only : nz
use fft
implicit none

real(rprec), dimension(:,:,lb:), intent(inout) :: f
integer :: jz
integer, intent(in) :: lb

! Loop through horizontal slices
do jz = lb, nz
    call dfftw_execute_dft_c2r(back, f(:,:,jz), f(:,:,jz))
enddo

end subroutine wave2phys

!*******************************************************************************
subroutine phys2waveF( u, uhat )
!*******************************************************************************
!
! Transform field from physical space to wavenumber space.
! (x,y) --> (kx,ky)
! This function is intended for fourier
!
use types, only : rprec
use param, only : nxp, ny, nz, lbz, kx_num, kxs_in
use fft
implicit none

integer :: jx, jz, ii, ir, iih, irh
real(rprec) :: const
real(rprec), dimension(nxp+2, ny, lbz:nz), intent(inout) :: u
real(rprec), dimension(ld, ny, lbz:nz), intent(out) :: uhat
! u is inout since it is being modified

const = 1._rprec / (nxp*ny)

! Initialize uhat
uhat = 0.0_rprec

! Loop through horizontal slices
do jz = lbz, nz
    u(:,:,jz) = const * u(:,:,jz) !! normalization
    call dfftw_execute_dft_r2c(forw_fourier, u(:,:,jz), u(:,:,jz) )
enddo

! Only interested in u values at kxs_in wavenumbers
do jx = 1, kx_num
    ! uhat indices
    iih = 2*jx
    irh = iih - 1

    ! u indices
    ii = 2 * int( kxs_in(jx) ) + 2
    ir = ii - 1

    uhat(irh:iih, :, :) = u(ir:ii, :, :)
end do

end subroutine phys2waveF

!*******************************************************************************
subroutine wave2physF( uhat, u )
!*******************************************************************************
!
! Transform field from wavenumber space to physical space.
! (kx,ky) --> (x,y)
! This function is intended to be used for fourier
!
use types, only : rprec
use param, only : nxp, ny, nz, lbz, kx_num, kxs_in
use fft
implicit none

integer :: jx, jz, ii, ir, iih, irh
real(rprec), dimension(ld, ny, lbz:nz), intent(in) :: uhat
real(rprec), dimension(nxp+2, ny, lbz:nz), intent(out) :: u

! Fill physical u with zeroes everywhere
u(:,:,:) = 0.0_rprec

! Place uhat values in appropriate u index before inverse fourier transform
do jx = 1, kx_num
    ! uhat indices
    iih = 2*jx
    irh = iih - 1

    ! u indices
    ii = 2 * int( kxs_in(jx) ) + 2
    ir = ii - 1

    u(ir:ii, :, :) = uhat( irh:iih, :, :)
end do

! Loop through horizontal slices
do jz = lbz, nz
    call dfftw_execute_dft_c2r(back_fourier, u(1:nxp+2,1:ny,jz), u(1:nxp+2,1:ny,jz))
enddo

end subroutine wave2physF

!*******************************************************************************
subroutine dft_direct_forw_2d_n_yonlyC_big(f)
!*******************************************************************************
! 
! Computes 2D DFT directly, no FFT in x-direction.
! 
! This function is inteded to be used for fourier
! 
use types, only: rprec
use param, only: ny2, kx_num
use fft
implicit none

integer :: jx, jy, ii, ir, jx_s
real(rprec), dimension(:,:), intent(inout) :: f
real(rprec) :: const
complex(rprec), dimension(ny2) :: fhat

const = 1._rprec / real(ny2,rprec)
f(:,:) = const * f(:,:)
fhat(:) = ( 0._rprec, 0._rprec )

do jx = 1, kx_num
    ! un-interleave the real array into a complex array
    jx_s = kx_veci( jx )
    ii = 2*jx_s ! imag index
    ir = ii-1 ! real index

    do jy = 1, ny2
        fhat(jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
    enddo

    ! remains in complex space
    call dfftw_execute_dft(ycomp_forw_big, fhat(1:ny2), fhat(1:ny2) )

    ! interleave the complex array into a real array
    do jy = 1, ny2
        f(ir,jy) = real( fhat(jy), rprec )
        f(ii,jy) = aimag( fhat(jy) ) !! or dimag?
    enddo
enddo

end subroutine dft_direct_forw_2d_n_yonlyC_big

!*******************************************************************************
subroutine dft_direct_back_2d_n_yonlyC_big(f)
!*******************************************************************************
! 
! Computes 2D inverse DFT directly, no FFT in x-direction.
! 
! This function is inteded to be used for fourier
! 
use types, only: rprec
use param, only: ny2, kx_num
use fft
implicit none

integer :: jx, jy, ii, ir, jx_s
real(rprec), dimension(:,:), intent(inout) :: f
complex(rprec), dimension(ny2) :: fhat

fhat(:) = ( 0._rprec, 0._rprec )

do jx = 1, kx_num
    ! un-interleave the real array into a complex array
    jx_s = kx_veci( jx )
    ii = 2*jx_s ! imag index
    ir = ii-1 ! real index

    do jy = 1, ny2
        fhat(jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
    enddo

    ! remains in complex space
    call dfftw_execute_dft(ycomp_back_big, fhat(1:ny2), fhat(1:ny2) )

    ! interleave the complex array into a real array
    do jy = 1, ny2
        f(ir,jy) = real( fhat(jy), rprec )
        f(ii,jy) = aimag( fhat(jy) ) !! or dimag?
    enddo
enddo

end subroutine dft_direct_back_2d_n_yonlyC_big

!*******************************************************************************
function convolve_rnl(f,g) result(out)
!*******************************************************************************
! 
! This function computes the convolution of f and g, which are both in Fourier
! space. The f and g arrays contain complex numbers stored as interleaved real
! arrays.
! 
! This function is intended to be used for fourier.
! 
! This function mimics convolve_rnl_old as was done in Joel's work, however
! reduces amount of space allocated for arrays. This function should eventually
! replace convolve_rnl_old.
! 
! Now edited to perform band-limited GQL in Fourier space.
! 
use types, only: rprec
use param, only: nx
use fft
use functions, only: interleave_r2c, interleave_c2r
#ifdef PPGQL
use param, only: thrx
#endif
implicit none

real(rprec), dimension(:,:), intent(in) :: f, g
real(rprec), allocatable, dimension(:,:) :: out

complex(rprec), allocatable, dimension(:,:) :: fc, gc, outc
integer :: jx
#ifdef PPGQL
integer :: ix, jx_end
#endif

integer :: ldh, nxh, nyh

ldh = size(f,1)    !! either ld or ld_big
nxh = ldh - 2      !! either nx or nx2
nyh = size(f,2)    !! either ny or ny2

allocate( out(ldh, nyh) )
allocate( fc(nxh, nyh) )
allocate( gc(nxh, nyh) )
allocate( outc(nxh, nyh) )

fc(:,:) = (0._rprec, 0._rprec )
gc(:,:) = (0._rprec, 0._rprec )
outc(:,:) = (0._rprec, 0._rprec )
out(:,:) = 0._rprec

! f, g are in kx space, structured as real arrays
! now re-structure as complex arrays
fc(1:nx,:) = interleave_r2c( f(:,:) )
gc(1:nx,:) = interleave_r2c( g(:,:) )

! UV part
outc(1,:) = outc(1,:) + fc(1,:) * gc(1,:)

! Uv, vU parts
do jx = 2, nxh/2
    outc(jx,:) = outc(jx,:) + fc(1,:) * gc(jx,:)
    outc(jx,:) = outc(jx,:) + fc(jx,:) * gc(1,:)
enddo

! < uv > part
do jx = 2, nxh/2
    outc(1,:) = outc(1,:) + fc(jx,:) * conjg( gc(jx,:) )
    outc(1,:) = outc(1,:) + conjg( fc(jx,: )) * gc(jx,:)
enddo

#ifdef PPGQL
! GQL has the same interactions as RNL above, only adding more interactions
! These interactions involve only nonzero wavenumbers

! Kxi+Kxi -> Kx, Add large scales to get large scales (Kxi,Kxi) only
do ix = 2, floor((thrx+2.0_rprec)/2.0_rprec)
    outc(ix+ix-1,:) = outc(ix+ix-1,:) + fc(ix,:) * gc(ix,:)
enddo

! Kxi+Kxj -> Kx, Add large scales to get large scales (Kxi,Kxj) Kxi ~= Kxj only
do ix = 2, thrx
    do jx = (ix+1), (thrx-ix+2)
        outc(ix+jx-1,:) = outc(ix+jx-1,:) + fc(ix,:) * gc(jx,:)
        outc(ix+jx-1,:) = outc(ix+jx-1,:) + fc(jx,:) * gc(ix,:)
    enddo
enddo

! Kxj-Kxi -> Kx, Subtract large scales to get large scales
do ix = 2, (thrx)
    do jx = (ix+1), thrx+1
        outc(jx-ix+1,:) = outc(jx-ix+1,:) + conjg( fc(ix,: )) * gc(jx,:)
        outc(jx-ix+1,:) = outc(jx-ix+1,:) + fc(jx,:) * conjg( gc(ix,:) )
    enddo
enddo

! kxj-kxi -> Kx, Subtract small scales to get large scales
do ix = (thrx+2), ((nxh/2)-1)
    jx_end = ix + thrx
    if ( jx_end > (nxh/2) ) then
        jx_end = (nxh/2)
    endif
    do jx = (ix+1), jx_end
        outc(jx-ix+1,:) = outc(jx-ix+1,:) + conjg( fc(ix,:) ) * gc(jx,:)
        outc(jx-ix+1,:) = outc(jx-ix+1,:) + fc(jx,:) * conjg( gc(ix,:) )
    enddo
enddo

! kxi+Kxj -> kx, Add large and small scales to get small scales
do ix = 2, (thrx+1)
    do jx = (thrx+2), ( (nxh/2)-ix+1 )
        outc(ix+jx-1,:) = outc(ix+jx-1,:) + fc(ix,:) * gc(jx,:)
        outc(ix+jx-1,:) = outc(ix+jx-1,:) + fc(jx,:) * gc(ix,:)
    enddo
enddo

! kxj-Kxi -> kx, Subtract large and small scales to get small scales
do ix = 2, (thrx+1)
    do jx = (thrx+ix+1), (nxh/2)
        outc(jx-ix+1,:) = outc(jx-ix+1,:) + conjg( fc(ix,:) ) * gc(jx,:)
        outc(jx-ix+1,:) = outc(jx-ix+1,:) + fc(jx,:) * conjg( gc(ix,:) )
    enddo
enddo

#endif

! re-structure the complex array into a real array
out(:,:) = interleave_c2r( outc(1:nx,:) )

return

end function convolve_rnl








!*******************************************************************************
function convolve_rnl_old(f,g) result(out)
!*******************************************************************************
! 
! This function computes the convolution of f and g, which are both in Fourier
! space. The f and g arrays contain complex numbers stored as interleaved real
! arrays.
! 
! This function is intended to be used for fourier
! 
use types, only: rprec
use param, only: nx
use fft
use functions, only: interleave_r2c, interleave_c2r
implicit none

real(rprec), dimension(:,:), intent(in) :: f, g
real(rprec), allocatable, dimension(:,:) :: out

complex(rprec), allocatable, dimension(:,:) :: fc, gc, outc
integer :: jx

integer :: ldh, nxh, nyh

ldh = size(f,1)    !! either ld or ld_big
nxh = ldh - 2      !! either nx or nx2
nyh = size(f,2)    !! either ny or ny2

allocate( out(ldh, nyh) )
allocate( fc(2*nxh - 1, nyh) )
allocate( gc(2*nxh - 1, nyh) )
allocate( outc(2*nxh - 1, nyh) )

fc(:,:) = (0._rprec, 0._rprec )
gc(:,:) = (0._rprec, 0._rprec )
outc(:,:) = (0._rprec, 0._rprec )
out(:,:) = 0._rprec

! f, g are in kx space, structured as real arrays
! now re-structure as complex arrays
fc(1:nx,:) = interleave_r2c( f(:,:) )
gc(1:nx,:) = interleave_r2c( g(:,:) )

! get the -kx modes from the +kx modes
do jx = 2, nxh/2
    fc(nxh-jx+2,:) = conjg(fc(jx,:))
    gc(nxh-jx+2,:) = conjg(gc(jx,:))
enddo

! UV part
outc(1,:) = outc(1,:) + fc(1,:) * gc(1,:)

! Uv, vU parts
do jx = 2, (2*nxh - 1)
    outc(jx,:) = outc(jx,:) + fc(1,:) * gc(jx,:)
    outc(jx,:) = outc(jx,:) + fc(jx,:) * gc(1,:)
enddo

! < uv > part
do jx = 2, nxh/2
    outc(1,:) = outc(1,:) + fc(jx,:) * conjg( gc(jx,:) )
    outc(1,:) = outc(1,:) + conjg( fc(jx,: )) * gc(jx,:)
enddo

outc(1:nxh-1,:) = outc(1:nxh-1,:) + outc(nxh+1:2*nxh-1,:)

! re-structure the complex array into a real array
out(:,:) = interleave_c2r( outc(1:nx,:) )

return

end function convolve_rnl_old

end module derivatives



