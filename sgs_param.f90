!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
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
module sgs_param
!*******************************************************************************
use types, only : rprec
use param, only : fourier

save
private rprec
public

! For all sgs models
integer ::jt_count
real(rprec) :: delta, nu
real(rprec), dimension(:,:,:),allocatable :: S11, S12, S22, S33, S13, S23
real(rprec), dimension(:,:,:),allocatable :: Nu_t !! eddy viscosity
real(rprec), dimension(:,:,:),allocatable :: Cs_opt2 !! (C_s)^2, Dynamic Smag Coef
real(rprec), dimension(:,:),  allocatable :: S
real(rprec), dimension(:,:),  allocatable :: nu_coef, nu_coef2

! RNL-LES only for sgs = 1
real(rprec), dimension(:,:,:), allocatable :: S11_big, S22_big, S33_big
real(rprec), dimension(:,:,:), allocatable :: S12_big, S13_big, S23_big
real(rprec), dimension(:,:),   allocatable :: S_big
real(rprec), dimension(:,:),   allocatable :: txx_big, tyy_big, tzz_big
real(rprec), dimension(:,:),   allocatable :: txy_big, txz_big, tyz_big
real(rprec), dimension(:,:),   allocatable :: nu_coef_big, nu_coef2_big
real(rprec), dimension(:,:,:), allocatable :: dudx_big, dudy_big, dudz_big
real(rprec), dimension(:,:,:), allocatable :: dvdx_big, dvdy_big, dvdz_big
real(rprec), dimension(:,:,:), allocatable :: dwdx_big, dwdy_big, dwdz_big
#ifdef PPCNDIFF
real(rprec), dimension(:,:),   allocatable :: txz_half1_big, txz_half2_big
real(rprec), dimension(:,:),   allocatable :: tyz_half1_big, tyz_half2_big
#endif

! For all dynamic models (2-5)
real(rprec), dimension(:,:,:),allocatable :: ee_now
real(rprec), dimension(:,:),  allocatable :: L11, L12, L13, L22, L23, L33
real(rprec), dimension(:,:),  allocatable :: M11, M12, M13, M22, M23, M33

real(rprec), dimension(:,:),  allocatable :: S11_bar, S12_bar, S13_bar
real(rprec), dimension(:,:),  allocatable :: S22_bar, S23_bar, S33_bar
real(rprec), dimension(:,:),  allocatable :: S_S11_bar, S_S12_bar, S_S13_bar
real(rprec), dimension(:,:),  allocatable :: S_S22_bar, S_S23_bar, S_S33_bar
real(rprec), dimension(:,:),  allocatable :: u_bar, v_bar, w_bar, S_bar

! For Lagrangian models (4,5)
real(rprec) :: lagran_dt = 0._rprec
real(rprec), parameter :: opftime = 1.5_rprec ! (Meneveau, Lund, Cabot; JFM 1996)
real(rprec), dimension(:,:,:), allocatable :: F_LM, F_MM, F_QN, F_NN
real(rprec), dimension(:,:,:), allocatable :: Beta, Tn_all

! For scale dependent model 5
real(rprec), dimension(:,:), allocatable :: Q11, Q12, Q13, Q22, Q23, Q33
real(rprec), dimension(:,:), allocatable :: N11, N12, N13, N22, N23, N33

! For scale dependent models (3,5)
real(rprec), dimension(:,:), allocatable :: S11_hat, S12_hat, S13_hat
real(rprec), dimension(:,:), allocatable :: S22_hat, S23_hat, S33_hat
real(rprec), dimension(:,:), allocatable :: S_S11_hat, S_S12_hat, S_S13_hat
real(rprec), dimension(:,:), allocatable :: S_S22_hat, S_S23_hat, S_S33_hat
real(rprec), dimension(:,:), allocatable :: u_hat, v_hat, w_hat, S_hat

! For Vreman SGS model
real(rprec) :: cvre

! The following are for dynamically updating T, the timescale for Lagrangian averaging
!   F_ee2 is the running average of (eij*eij)^2
!   F_deedt2 is the running average of [d(eij*eij)/dt]^2
!   ee_past is the array (eij*eij) for the past timestep
#ifdef PPDYN_TN
real(rprec), dimension(:,:,:), allocatable :: F_ee2, F_deedt2, ee_past
#endif

contains

!*******************************************************************************
subroutine sgs_param_init ()
!*******************************************************************************
use param, only : ld, ny, nz, lbz, molec, nu_molec, u_star,          &
    z_i, dx, dy, dz, sgs_model
use param, only : fourier, ld_big, ny2
use test_filtermodule, only : filter_size

implicit none

! For all sgs models:
allocate ( S11(ld,ny,nz) ); S11 = 0._rprec
allocate ( S12(ld,ny,nz) ); S12 = 0._rprec
allocate ( S13(ld,ny,nz) ); S13 = 0._rprec
allocate ( S22(ld,ny,nz) ); S22 = 0._rprec
allocate ( S23(ld,ny,nz) ); S23 = 0._rprec
allocate ( S33(ld,ny,nz) ); S33 = 0._rprec
allocate ( Nu_t(ld,ny,nz) ); Nu_t = 0._rprec
allocate ( Cs_opt2(ld,ny,nz) ); Cs_opt2 = 0._rprec
allocate ( S(ld,ny) ); S = 0._rprec
allocate ( nu_coef(ld,ny) ); nu_coef = 0._rprec
allocate ( nu_coef2(ld,ny) ); nu_coef2 = 0._rprec

! RNL-LES only for sgs = 1
if (fourier) then
    allocate ( S11_big(ld_big,ny2,nz) ); S11_big = 0.0_rprec
    allocate ( S22_big(ld_big,ny2,nz) ); S22_big = 0.0_rprec
    allocate ( S33_big(ld_big,ny2,nz) ); S33_big = 0.0_rprec
    allocate ( S12_big(ld_big,ny2,nz) ); S12_big = 0.0_rprec
    allocate ( S13_big(ld_big,ny2,nz) ); S13_big = 0.0_rprec
    allocate ( S23_big(ld_big,ny2,nz) ); S23_big = 0.0_rprec
    allocate ( S_big(ld_big,ny2) ); S_big = 0.0_rprec
    allocate ( txx_big(ld_big,ny2) ); txx_big = 0.0_rprec
    allocate ( tyy_big(ld_big,ny2) ); tyy_big = 0.0_rprec
    allocate ( tzz_big(ld_big,ny2) ); tzz_big = 0.0_rprec
    allocate ( txy_big(ld_big,ny2) ); txy_big = 0.0_rprec
    allocate ( txz_big(ld_big,ny2) ); txz_big = 0.0_rprec
    allocate ( tyz_big(ld_big,ny2) ); tyz_big = 0.0_rprec
    allocate ( nu_coef_big(ld_big,ny2) ); nu_coef_big = 0.0_rprec
    allocate ( nu_coef2_big(ld_big,ny2) ); nu_coef2_big = 0.0_rprec
    allocate ( dudx_big(ld_big,ny2,nz) ); dudx_big = 0.0_rprec
    allocate ( dudy_big(ld_big,ny2,nz) ); dudy_big = 0.0_rprec
    allocate ( dudz_big(ld_big,ny2,nz) ); dudz_big = 0.0_rprec
    allocate ( dvdx_big(ld_big,ny2,nz) ); dvdx_big = 0.0_rprec
    allocate ( dvdy_big(ld_big,ny2,nz) ); dvdy_big = 0.0_rprec
    allocate ( dvdz_big(ld_big,ny2,nz) ); dvdz_big = 0.0_rprec
    allocate ( dwdx_big(ld_big,ny2,nz) ); dwdx_big = 0.0_rprec
    allocate ( dwdy_big(ld_big,ny2,nz) ); dwdy_big = 0.0_rprec
    allocate ( dwdz_big(ld_big,ny2,nz) ); dwdz_big = 0.0_rprec
#ifdef PPCNDIFF
    allocate ( txz_half1_big(ld_big,ny2) ); txz_half1_big = 0.0_rprec
    allocate ( txz_half2_big(ld_big,ny2) ); txz_half2_big = 0.0_rprec
    allocate ( tyz_half1_big(ld_big,ny2) ); tyz_half1_big = 0.0_rprec
    allocate ( tyz_half2_big(ld_big,ny2) ); tyz_half2_big = 0.0_rprec
#endif
endif

! For dynamic models:
if (sgs_model .ne. 1) then
    allocate ( ee_now(ld,ny,lbz:nz) ); ee_now = 0.0_rprec
    allocate ( L11(ld,ny) ); L11 = 0._rprec
    allocate ( L12(ld,ny) ); L12 = 0._rprec
    allocate ( L13(ld,ny) ); L13 = 0._rprec
    allocate ( L22(ld,ny) ); L22 = 0._rprec
    allocate ( L23(ld,ny) ); L23 = 0._rprec
    allocate ( L33(ld,ny) ); L33 = 0._rprec
    allocate ( M11(ld,ny) ); M11 = 0._rprec
    allocate ( M12(ld,ny) ); M12 = 0._rprec
    allocate ( M13(ld,ny) ); M13 = 0._rprec
    allocate ( M22(ld,ny) ); M22 = 0._rprec
    allocate ( M23(ld,ny) ); M23 = 0._rprec
    allocate ( M33(ld,ny) ); M33 = 0._rprec

    allocate ( S11_bar(ld,ny) ); S11_bar = 0._rprec
    allocate ( S12_bar(ld,ny) ); S12_bar = 0._rprec
    allocate ( S13_bar(ld,ny) ); S13_bar = 0._rprec
    allocate ( S22_bar(ld,ny) ); S22_bar = 0._rprec
    allocate ( S23_bar(ld,ny) ); S23_bar = 0._rprec
    allocate ( S33_bar(ld,ny) ); S33_bar = 0._rprec
    allocate ( S_S11_bar(ld,ny) ); S_S11_bar = 0._rprec
    allocate ( S_S12_bar(ld,ny) ); S_S12_bar = 0._rprec
    allocate ( S_S13_bar(ld,ny) ); S_S13_bar = 0._rprec
    allocate ( S_S22_bar(ld,ny) ); S_S22_bar = 0._rprec
    allocate ( S_S23_bar(ld,ny) ); S_S23_bar = 0._rprec
    allocate ( S_S33_bar(ld,ny) ); S_S33_bar = 0._rprec

    allocate ( S_bar(ld,ny) ); S_bar(ld,ny) = 0.0_rprec
    allocate ( u_bar(ld,ny) ); u_bar = 0.0_rprec
    allocate ( v_bar(ld,ny) ); v_bar = 0.0_rprec
    allocate ( w_bar(ld,ny) ); w_bar = 0.0_rprec

endif

! For Lagrangian models:
allocate ( F_LM(ld,ny,lbz:nz) ); F_LM = 0.0_rprec
allocate ( F_MM(ld,ny,lbz:nz) ); F_MM = 0.0_rprec
allocate ( F_QN(ld,ny,lbz:nz) ); F_QN = 0.0_rprec
allocate ( F_NN(ld,ny,lbz:nz) ); F_NN = 0.0_rprec
if ((sgs_model .eq. 4).or.(sgs_model .eq. 5)) then
    allocate ( Beta(ld,ny,lbz:nz) ); Beta = 0.0_rprec
    allocate ( Tn_all(ld,ny,lbz:nz) ); Tn_all  = 0._rprec

#ifdef PPDYN_TN
    ! Lagrangian zero-crossing time scale variables
    allocate ( F_ee2(ld,ny,lbz:nz) ) ; F_ee2 = 0._rprec
    allocate ( F_deedt2(ld,ny,lbz:nz) ); F_deedt2 = 0._rprec
    allocate ( ee_past(ld,ny,lbz:nz) );  ee_past = 0.0_rprec
#endif

endif

! For scale dependent models:
if ((sgs_model .eq. 3).or.(sgs_model .eq. 5)) then
    allocate ( S11_hat(ld,ny) ); S11_hat = 0._rprec
    allocate ( S12_hat(ld,ny) ); S12_hat = 0._rprec
    allocate ( S13_hat(ld,ny) ); S13_hat = 0._rprec
    allocate ( S22_hat(ld,ny) ); S22_hat = 0._rprec
    allocate ( S23_hat(ld,ny) ); S23_hat = 0._rprec
    allocate ( S33_hat(ld,ny) ); S33_hat = 0._rprec
    allocate ( S_S11_hat(ld,ny) ); S_S11_hat = 0._rprec
    allocate ( S_S12_hat(ld,ny) ); S_S12_hat = 0._rprec
    allocate ( S_S13_hat(ld,ny) ); S_S13_hat = 0._rprec
    allocate ( S_S22_hat(ld,ny) ); S_S22_hat = 0._rprec
    allocate ( S_S23_hat(ld,ny) ); S_S23_hat = 0._rprec
    allocate ( S_S33_hat(ld,ny) ); S_S33_hat = 0._rprec
    allocate ( S_hat(ld,ny) ); S_hat = 0._rprec
    allocate ( u_hat(ld,ny) ); u_hat = 0._rprec
    allocate ( v_hat(ld,ny) ); v_hat = 0._rprec
    allocate ( w_hat(ld,ny) ); w_hat = 0._rprec
endif

if (sgs_model .eq. 5) then
    allocate ( Q11(ld,ny) ); Q11 = 0._rprec
    allocate ( Q12(ld,ny) ); Q12 = 0._rprec
    allocate ( Q13(ld,ny) ); Q13 = 0._rprec
    allocate ( Q22(ld,ny) ); Q22 = 0._rprec
    allocate ( Q23(ld,ny) ); Q23 = 0._rprec
    allocate ( Q33(ld,ny) ); Q33 = 0._rprec
    allocate ( N11(ld,ny) ); N11 = 0._rprec
    allocate ( N12(ld,ny) ); N12 = 0._rprec
    allocate ( N13(ld,ny) ); N13 = 0._rprec
    allocate ( N22(ld,ny) ); N22 = 0._rprec
    allocate ( N23(ld,ny) ); N23 = 0._rprec
    allocate ( N33(ld,ny) ); N33 = 0._rprec
endif

if (sgs_model >= 6) then
    cvre = 0.07_rprec
endif

! Set dimensionless constants
delta = filter_size*(dy*dz)**(1._rprec/2._rprec)
!delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)
if (molec) then
    nu = (nu_molec/(u_star*z_i))
else
    nu = 0._rprec
end if

end subroutine sgs_param_init

end module sgs_param
