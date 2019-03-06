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

!*******************************************************************************
module sim_param
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz, nxp, BOGUS
implicit none

save
public

logical :: sim_param_initialized = .false.

real(rprec), dimension(:,:,:), allocatable :: u, v, w,                         &
    dudx, dudy, dudz, dvdx, dvdy, dvdz,  dwdx, dwdy, dwdz,                     &
    RHSx, RHSy, RHSz, RHSx_f, RHSy_f, RHSz_f,                                  &
    dpdx, dpdy, dpdz, txx, txy, tyy,                                           &
    txz, tyz, tzz, divtx, divty, divtz,                                        &
    fx, fy, fz, fxa, fya, fza

#ifdef PPRNL
real(rprec), dimension(:,:,:), allocatable :: u_pert, v_pert, w_pert,          &
    dudy_pert, dudz_pert, dvdx_pert, dvdz_pert, dwdx_pert, dwdy_pert,          &
    RHSx_pert, RHSy_pert, RHSz_pert
#endif

#ifdef PPGQL
real(rprec), dimension(:,:,:), allocatable :: u_low, v_low, w_low,             &
    dudy_low, dudz_low, dvdx_low, dvdz_low, dwdx_low, dwdy_low,                &
    RHSx_low, RHSy_low, RHSz_low,                                              &
    u_high, v_high, w_high,                                                    &
    dudy_high, dudz_high, dvdx_high, dvdz_high, dwdx_high, dwdy_high,          &
    RHSx_high, RHSy_high, RHSz_high
#endif

real(rprec), target, dimension(:,:,:), allocatable :: p

#ifdef PPMAPPING
real(rprec), dimension(:), allocatable :: JACO1, JACO2
real(rprec), dimension(:), allocatable :: mesh_stretch, delta_stretch
#endif

real(rprec), dimension(:,:,:), allocatable :: uF, vF, wF
real(rprec), target, dimension(:,:,:), allocatable :: pF
real(rprec), dimension(:,:,:), allocatable :: dudxF, dudyF, dudzF,             &
    dvdxF, dvdyF, dvdzF, dwdxF, dwdyF, dwdzF
real(rprec), dimension(:,:,:), allocatable :: txzF, tyzF

logical, dimension(:), allocatable :: zhyb
integer :: jz_int = -1
integer :: coord_int = -1
integer, dimension(:), allocatable :: coord_intv

contains

!*******************************************************************************
subroutine sim_param_init ()
!*******************************************************************************
!
! This subroutine initilizes all global arrays defined in the sim_param
! module. Here they are allocated and initialized to zero.
!
implicit none

allocate ( u(ld, ny, lbz:nz) ); u = 0.0_rprec
allocate ( v(ld, ny, lbz:nz) ); v = 0.0_rprec
allocate ( w(ld, ny, lbz:nz) ); w = 0.0_rprec
allocate ( dudx(ld, ny, lbz:nz) ); dudx = 0.0_rprec
allocate ( dudy(ld, ny, lbz:nz) ); dudy = 0.0_rprec
allocate ( dudz(ld, ny, lbz:nz) ); dudz = 0.0_rprec
allocate ( dvdx(ld, ny, lbz:nz) ); dvdx = 0.0_rprec
allocate ( dvdy(ld, ny, lbz:nz) ); dvdy = 0.0_rprec
allocate ( dvdz(ld, ny, lbz:nz) ); dvdz = 0.0_rprec
allocate ( dwdx(ld, ny, lbz:nz) ); dwdx = 0.0_rprec
allocate ( dwdy(ld, ny, lbz:nz) ); dwdy = 0.0_rprec
allocate ( dwdz(ld, ny, lbz:nz) ); dwdz = 0.0_rprec
allocate ( RHSx(ld, ny, lbz:nz) ); RHSx = 0.0_rprec
allocate ( RHSy(ld, ny, lbz:nz) ); RHSy = 0.0_rprec
allocate ( RHSz(ld, ny, lbz:nz) ); RHSz = 0.0_rprec
allocate ( RHSx_f(ld, ny, lbz:nz) ); RHSx_f = 0.0_rprec
allocate ( RHSy_f(ld, ny, lbz:nz) ); RHSy_f = 0.0_rprec
allocate ( RHSz_f(ld, ny, lbz:nz) ); RHSz_f = 0.0_rprec
allocate ( dpdx(ld, ny, nz) ); dpdx = 0.0_rprec
allocate ( dpdy(ld, ny, nz) ); dpdy = 0.0_rprec
allocate ( dpdz(ld, ny, nz) ); dpdz = 0.0_rprec
allocate ( txx(ld, ny, lbz:nz) ); txx = 0.0_rprec
allocate ( txy(ld, ny, lbz:nz) ); txy = 0.0_rprec
allocate ( tyy(ld, ny, lbz:nz) ); tyy = 0.0_rprec
allocate ( txz(ld, ny, lbz:nz) ); txz = 0.0_rprec
allocate ( tyz(ld, ny, lbz:nz) ); tyz = 0.0_rprec
allocate ( tzz(ld, ny, lbz:nz) ); tzz = 0.0_rprec
allocate ( p(ld, ny, 0:nz) ); p = 0.0_rprec
allocate ( divtx(ld, ny, lbz:nz) ); divtx = 0.0_rprec
allocate ( divty(ld, ny, lbz:nz) ); divty = 0.0_rprec
allocate ( divtz(ld, ny, lbz:nz) ); divtz = 0.0_rprec

#ifdef PPRNL
allocate ( u_pert(ld, ny, lbz:nz) ); u_pert = 0.0_rprec
allocate ( v_pert(ld, ny, lbz:nz) ); v_pert = 0.0_rprec
allocate ( w_pert(ld, ny, lbz:nz) ); w_pert = 0.0_rprec
! note dudx_pert, dvdy_pert, and dwdz_pert are not here 
! since they are not needed for convec subroutine
allocate ( dudy_pert(ld, ny, lbz:nz) ); dudy_pert = 0.0_rprec
allocate ( dudz_pert(ld, ny, lbz:nz) ); dudz_pert = 0.0_rprec
allocate ( dvdx_pert(ld, ny, lbz:nz) ); dvdx_pert = 0.0_rprec
allocate ( dvdz_pert(ld, ny, lbz:nz) ); dvdz_pert = 0.0_rprec
allocate ( dwdx_pert(ld, ny, lbz:nz) ); dwdx_pert = 0.0_rprec
allocate ( dwdy_pert(ld, ny, lbz:nz) ); dwdy_pert = 0.0_rprec
allocate ( RHSx_pert(ld, ny, lbz:nz) ); RHSx_pert = 0.0_rprec
allocate ( RHSy_pert(ld, ny, lbz:nz) ); RHSy_pert = 0.0_rprec
allocate ( RHSz_pert(ld, ny, lbz:nz) ); RHSz_pert = 0.0_rprec
#endif

#ifdef PPGQL
allocate ( u_low(ld, ny, lbz:nz) ); u_low = 0.0_rprec
allocate ( v_low(ld, ny, lbz:nz) ); v_low = 0.0_rprec
allocate ( w_low(ld, ny, lbz:nz) ); w_low = 0.0_rprec
allocate ( dudy_low(ld, ny, lbz:nz) ); dudy_low = 0.0_rprec
allocate ( dudz_low(ld, ny, lbz:nz) ); dudz_low = 0.0_rprec
allocate ( dvdx_low(ld, ny, lbz:nz) ); dvdx_low = 0.0_rprec
allocate ( dvdz_low(ld, ny, lbz:nz) ); dvdz_low = 0.0_rprec
allocate ( dwdx_low(ld, ny, lbz:nz) ); dwdx_low = 0.0_rprec
allocate ( dwdy_low(ld, ny, lbz:nz) ); dwdy_low = 0.0_rprec
allocate ( RHSx_low(ld, ny, lbz:nz) ); RHSx_low = 0.0_rprec
allocate ( RHSy_low(ld, ny, lbz:nz) ); RHSy_low = 0.0_rprec
allocate ( RHSz_low(ld, ny, lbz:nz) ); RHSz_low = 0.0_rprec
allocate ( u_high(ld, ny, lbz:nz) ); u_high = 0.0_rprec
allocate ( v_high(ld, ny, lbz:nz) ); v_high = 0.0_rprec
allocate ( w_high(ld, ny, lbz:nz) ); w_high = 0.0_rprec
allocate ( dudy_high(ld, ny, lbz:nz) ); dudy_high = 0.0_rprec
allocate ( dudz_high(ld, ny, lbz:nz) ); dudz_high = 0.0_rprec
allocate ( dvdx_high(ld, ny, lbz:nz) ); dvdx_high = 0.0_rprec
allocate ( dvdz_high(ld, ny, lbz:nz) ); dvdz_high = 0.0_rprec
allocate ( dwdx_high(ld, ny, lbz:nz) ); dwdx_high = 0.0_rprec
allocate ( dwdy_high(ld, ny, lbz:nz) ); dwdy_high = 0.0_rprec
allocate ( RHSx_high(ld, ny, lbz:nz) ); RHSx_high = 0.0_rprec
allocate ( RHSy_high(ld, ny, lbz:nz) ); RHSy_high = 0.0_rprec
allocate ( RHSz_high(ld, ny, lbz:nz) ); RHSz_high = 0.0_rprec
#endif

allocate ( uF(nxp+2, ny, lbz:nz) ); uF = 0.0_rprec
allocate ( vF(nxp+2, ny, lbz:nz) ); vF = 0.0_rprec
allocate ( wF(nxp+2, ny, lbz:nz) ); wF = 0.0_rprec
allocate ( pF(nxp+2, ny, 0:nz) ); pF = 0.0_rprec
allocate ( dudxF(nxp+2, ny, lbz:nz) ); dudxF = 0.0_rprec
allocate ( dudyF(nxp+2, ny, lbz:nz) ); dudyF = 0.0_rprec
allocate ( dudzF(nxp+2, ny, lbz:nz) ); dudzF = 0.0_rprec
allocate ( dvdxF(nxp+2, ny, lbz:nz) ); dvdxF = 0.0_rprec
allocate ( dvdyF(nxp+2, ny, lbz:nz) ); dvdyF = 0.0_rprec
allocate ( dvdzF(nxp+2, ny, lbz:nz) ); dvdzF = 0.0_rprec
allocate ( dwdxF(nxp+2, ny, lbz:nz) ); dwdxF = 0.0_rprec
allocate ( dwdyF(nxp+2, ny, lbz:nz) ); dwdyF = 0.0_rprec
allocate ( dwdzF(nxp+2, ny, lbz:nz) ); dwdzF = 0.0_rprec
allocate ( txzF(nxp+2, ny, lbz:nz) ); txzF = 0.0_rprec
allocate ( tyzF(nxp+2, ny, lbz:nz) ); tyzF = 0.0_rprec

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
allocate ( fxa(ld, ny, lbz:nz) ); fxa = 0.0_rprec
allocate ( fya(ld, ny, lbz:nz) ); fya = 0.0_rprec
allocate ( fza(ld, ny, lbz:nz) ); fza = 0.0_rprec
#endif

#if defined(PPLVLSET) || defined(PPATM)
allocate ( fx(ld, ny, nz) ); fx = 0.0_rprec
allocate ( fy(ld, ny, nz) ); fy = 0.0_rprec
allocate ( fz(ld, ny, nz) ); fz = 0.0_rprec
#endif

#ifdef PPMAPPING
allocate ( JACO1(lbz:nz) ); JACO1 = 1/BOGUS
allocate ( JACO2(lbz:nz) ); JACO2 = 1/BOGUS
allocate ( mesh_stretch(lbz:nz)); mesh_stretch = BOGUS
allocate ( delta_stretch(lbz:nz)); delta_stretch = BOGUS
#endif

allocate( zhyb(lbz:nz)); zhyb = .false.

sim_param_initialized = .true.

end subroutine sim_param_init

end module sim_param
