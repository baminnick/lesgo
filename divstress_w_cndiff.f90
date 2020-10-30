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
subroutine divstress_w_cndiff(divt, tx, ty)
!******************************************************************************
!
! This subroutine provides divt for 1:nz. MPI provides 1:nz-1,
! except at top, where 1:nz is provided
! 
! This subroutine is copied from divstress_w, however does not compute
! dtzdz, as this is computed in diff_stag_array.f90 since it is to be 
! treated implicitly
!
use types, only : rprec
use param, only : ld, nx, ny, nz, BOGUS, lbz
use derivatives, only : ddx, ddy, ddz_uv
implicit none

real(rprec),dimension(ld,ny,lbz:nz),intent(out)::divt
real(rprec),dimension (ld, ny, lbz:nz), intent (in) :: tx, ty
real(rprec),dimension(ld,ny,lbz:nz) :: dtxdx, dtydy
integer :: jx, jy, jz

! compute stress gradients
! tx 1:nz => dtxdx 1:nz
! ty 1:nz => dtydy 1:nz
call ddx(tx, dtxdx, lbz)
call ddy(ty, dtydy, lbz)

#ifdef PPSAFETYMODE
#ifdef PPMPI
dtxdx(:, :, 0) = BOGUS
dtydy(:, :, 0) = BOGUS
divt(:, :, 0) = BOGUS
#endif
#endif

do jz = 1, nz
do jy = 1, ny
do jx = 1, nx
    divt(jx,jy,jz) = dtxdx(jx,jy,jz) + dtydy(jx,jy,jz)
end do 
end do
end do

! set ld-1, ld to 0 (could maybe do BOGUS)
divt(ld-1:ld, :, 1:nz-1) = 0._rprec

end subroutine divstress_w_cndiff
