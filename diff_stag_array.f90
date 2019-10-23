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
subroutine diff_stag_array()
!*******************************************************************************
!
! Calculate the intermediate velocity from implicit CN scheme on exit.
!
use types, only : rprec
use param
use messages
use sim_param, only : u, txz_half2, RHSx, RHSx_f
use sgs_param, only : nu, Nu_t
use fft
#ifdef PPMAPPING
use sim_param, only : JACO1, JACO2
#endif

implicit none

real(rprec), dimension(ld,ny,1:nz-1) :: Rx
real(rprec), dimension(ld,ny,lbz:nz) :: dtxzdz_rhs
real(rprec), dimension(nx,ny,1:nz-1) :: ax, bx, cx
real(rprec) :: nu_a, nu_b, nu_c, nu_r
integer :: jz_min, jz_max

! Get the RHS ready
! Initialize with the explicit terms
Rx(:,:,1:nz-1) = u(:,:,1:nz-1) +                                     &
    dt * ( tadv1 * RHSx(:,:,nz-1) + tadv2 * RHSx_f(:,:,1:nz-1) )

! Add explicit portion of Crank-Nicolson
call ddz_w(txz_half2, dtxzdz_rhs, lbz)
dtxzdz_rhs(ld-1,ld, :, 1:nz-1) = 0._rprec
#ifdef PPSAFETYMODE
#ifdef PPMPI
dtxzdz_rhs(:,:,0) = BOGUS
#endif
dtxzdz_rhs(:,:,nz) = BOGUS
#endif
! Assuming fixed dt here!
Rx(:,:,1:nz-1) = Rx(:,:,1:nz-1) - dt * 0.5_rprec * dtxzdz_rhs(:,:,1:nz-1)

! Get bottom row 
if (coord == 0) then
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    select case (lbc_mom)
        ! Stress free
        case (0)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on w, not needed here
                    nu_c = Nu_t(jx,jy,2) + nu
                else
                    nu_c = nu
                end if
                ! txz(jx,jy,1) = 0, so nothing added to RHS
                b(jx,jy,1) = 1._rprec + (dt/(2._rprec*dz))*(nu_c/dz)
                c(jx,jy,1) = -(dt/(2._rprec*dz))*(1._rprec/dz)*nu_c
            end do
            end do

        ! DNS BC or Wall-resolved
        case (1)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on uvp, not needed here
                    nu_c = Nu_t(jx,jy,2) + nu
                else
                    nu_c = nu
                end if
                ! Discretized txz(jx,jy,1) as in wallstress,
                ! Therefore BC treated implicitly
                b(jx,jy,1) = 1._rprec + (dt/(2._rprec*dz))*        &
                    ((nu_c/dz) + (nu/(0.5_rprec*dz)))
                c(jx,jy,1) = -(dt/(2._rprec*dz))*(1._rprec/dz)*nu_c
                Rx(jx,jy,1) = Rx(jx,jy,1) + (dt/(2._rprec*dz))*    &
                    (nu/(0.5_rprec*dz))*ubot
            end do
            end do

        ! Wall-model 
        case (2:)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                    ! Nu_t(jx,jy,1) on uvp, not needed here
                    nu_c = Nu_t(jx,jy,2) + nu
                else
                    nu_c = nu
                end if
                ! Treating txz(jx,jy,1) from wallstress explicitly
                b(jx,jy,1) = 1._rprec + (dt/(2._rprec*dz))*(nu_c/dz)
                c(jx,jy,1) = -(dt/(2._rprec*dz))*(1._rprec/dz)*nu_c
                Rx(jx,jy,1) = Rx(jx,jy,1) - (dt/(2._rprec*dz))*txz(jx,jy,1)
            end do
            end do

    end select
    jz_min = 2
else
    jz_min = 1
endif

! Get top row
#ifdef PPMPI
if (coord == nproc-1) then
#endif
#ifdef PPSAFETYMODE
    c(:,:,nz-1) = BOGUS
#endif
    select case (ubc_mom)
        ! Stress free
        case (0)
            do jy = 1, ny
            do jx = 1, nx
                if (sgs) then
                ! Nu_t(jx,jy,nz) on w, not needed here

                else

                end if
                ! txz(jx,jy,1) = 0, so nothing added to RHS


            end do
            end do

        ! DNS BC or Wall-resolved
        case (1)


        ! Wall-model
        case (2)


    end select
    jz_max = nz-2
#ifdef PPMPI
else
    jz_max = nz-1
endif
#endif

! Compute coefficients
! Treating SGS viscosity explicitly!
do jz = jz_min, jz_max
do jy = 1, ny
do jx = 1, nx
    if (sgs) then
        a_nu = Nu_t(jx,jy,jz) + nu
        b_nu = Nu_t(jx,jy,jz+1) + Nu_t(jx,jy,jz) + 2._rprec*nu
        c_nu = Nu_t(jx,jy,jz+1) + nu
    else
        a_nu = nu
        b_nu = 2._rprec*nu
        c_nu = nu
    endif

    ax(jx, jy, jz) = -(dt/(2._rprec*dz))*(1._rprec/dz)*a_nu
    bx(jx, jy, jz) = 1._rprec + (dt/(2._rprec*dz))*(1._rprec/dz)*b_nu
    cx(jx, jy, jz) = -(dt/(2._rprec*dz))*(1._rprec/dz)*c_nu
end do
end do
end do

! Find intermediate velocity in TDMA
call tridag_array_diff (a, b, c, Rx, u)

end subroutine diff_stag_array
