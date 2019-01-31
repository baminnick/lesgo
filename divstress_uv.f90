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
subroutine divstress_uv (divtx, divty, txx, txy, txz, tyy, tyz)
!*******************************************************************************
!
! This subroutine provides divt for 1:nz. MPI provides 1:nz-1,
! except at top, where 1:nz is provided
!
use types, only : rprec
use param, only : ld, ny, nz, BOGUS, lbz
use derivatives, only : ddx, ddy, ddz_w, ddxy
implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: divtx, divty
real(rprec), dimension(ld, ny, lbz:nz), intent (in) :: txx, txy, txz, tyy, tyz
real(rprec), dimension(ld,ny,lbz:nz) :: dtxdx, dtydy, dtzdz
real(rprec), dimension(ld,ny,lbz:nz) :: dtxdx2, dtydy2, dtzdz2

! compute stress gradients
! MPI: tx 1:nz-1 => dtxdx 1:nz-1
call ddx(txx, dtxdx, lbz)  ! really should replace with ddxy (save an fft)
call ddy(tyy, dtydy2, lbz)
call ddxy(txy , dtxdx2, dtydy, lbz)

! MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(txz, dtzdz, lbz)

! MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(tyz, dtzdz2, lbz)

! MPI following comment only true at bottom process
! the following gives bad results...but it seems like i the
! issue should be taken care of somewhere
! need to correct wall level, since tz will be on uv-node there
!      dtzdz(:,:,1) = (tz(:,:,2)-tz(:,:,1))/(0.5*dz)

! only 1:nz-1 are valid
divtx(:,:,1:nz-1) = dtxdx(:,:,1:nz-1) + dtydy(:,:,1:nz-1) + dtzdz(:,:,1:nz-1)

! Set ld-1, ld to 0 (or could do BOGUS)
divtx(ld-1:ld, :, 1:nz-1) = 0._rprec

#ifdef PPSAFETYMODE
#ifdef PPMPI
divtx(:,:,0) = BOGUS
#endif
divtx(:,:,nz) = BOGUS
#endif

! only 1:nz-1 are valid
divty(:,:,1:nz-1) = dtxdx2(:,:,1:nz-1) + dtydy2(:,:,1:nz-1) + dtzdz2(:,:,1:nz-1)

! Set ld-1, ld to 0 (or could do BOGUS)
divty(ld-1:ld,:,1:nz-1) = 0._rprec

#ifdef PPSAFETYMODE
#ifdef PPMPI
divty(:,:,0) = BOGUS
#endif
divty(:,:,nz) = BOGUS
#endif

end subroutine divstress_uv

