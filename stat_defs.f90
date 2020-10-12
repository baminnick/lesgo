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
module stat_defs
!*****************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lh
#ifdef PPTURBINES
use turbine_indicator
#endif

save
public

type point_t
    integer :: istart, jstart, kstart, coord
    real(rprec) :: xdiff, ydiff, zdiff
    integer :: fid
end type point_t

type plane_t
    integer :: istart
    real(rprec) :: ldiff
end type plane_t

type zplane_t
    integer :: istart, coord
    real(rprec) :: ldiff
end type zplane_t

type spectra_t
    real(rprec), dimension(:), allocatable :: power
    integer :: istart, coord
    real(rprec) :: ldiff
end type spectra_t

real(rprec) :: spectra_total_time
real(rprec) :: tavg_total_time

! Time between calls of tavg_compute, built by summing dt
real(rprec) :: tavg_dt
! Switch for determining if time averaging has been initialized
logical :: tavg_initialized = .false.

!  Sums performed over time
type tavg_t
    real(rprec) :: u, v, w, u_w, v_w, w_uv, p
    real(rprec) :: txx, tyy, tzz, txy, txz, tyz
    real(rprec) :: u2, v2, w2, uv, uw, vw, fx, fy, fz
    ! real(rprec) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
end type tavg_t

type rs_t
    real(rprec) :: up2, vp2, wp2, upvp, upwp, vpwp
end type rs_t

#ifdef PPSCALARS
type tavg_scal_t
    real(rprec) :: theta, theta2, utheta, vtheta, wtheta, pi_z
end type tavg_scal_t

type rs_scal_t
    real(rprec) :: thetap2, upthetap, vpthetap, wpthetap
end type rs_scal_t

#ifdef PPOUTPUT_TURBSPEC
type tavg_scal_turbspec_t
    complex(rprec) :: thetaf
    real(rprec) :: theta2, utheta, vtheta, wtheta
end type tavg_scal_turbspec_t

type scal_turbspec_t
    real(rprec) :: thetap2, upthetap, vpthetap, wpthetap
end type scal_turbspec_t
#endif
#endif

type tavg_vort_t
    real(rprec) :: vortx, vorty, vortz
    real(rprec) :: vortx2, vorty2, vortz2
end type tavg_vort_t

type vortrms_t
    real(rprec) :: vortxrms, vortyrms, vortzrms
end type vortrms_t

#ifdef PPOUTPUT_SGS
type tavg_sgs_t
    real(rprec) :: cs_opt2, Nu_t
end type tavg_sgs_t
#endif

#ifdef PPOUTPUT_BUDGET
type tavg_budget_t
    real(rprec) :: p
    real(rprec) :: uu, vv, ww, uv, uw, vw
    real(rprec) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(rprec) :: dpdx, dpdy, dpdz
    real(rprec) :: ududx, ududy, ududz, udvdx, udvdy, udvdz, udwdx, udwdy, udwdz
    real(rprec) :: vdudx, vdudy, vdudz, vdvdx, vdvdy, vdvdz, vdwdx, vdwdy, vdwdz
    real(rprec) :: wdudx, wdudy, wdudz, wdvdx, wdvdy, wdvdz, wdwdx, wdwdy, wdwdz
    real(rprec) :: uududx, uvdudy, uwdudz, uudvdx, uvdvdy, uwdvdz, uudwdx, uvdwdy, uwdwdz
    real(rprec) :: vududx, vvdudy, vwdudz, vudvdx, vvdvdy, vwdvdz, vudwdx, vvdwdy, vwdwdz
    real(rprec) :: wududx, wvdudy, wwdudz, wudvdx, wvdvdy, wwdvdz, wudwdx, wvdwdy, wwdwdz
    real(rprec) :: uxux, uyuy, uzuz, vxvx, vyvy, vzvz, wxwx, wywy, wzwz
    real(rprec) :: uxvx, uyvy, uzvz, uxwx, uywy, uzwz, vxwx, vywy, vzwz
    ! real(rprec) :: uyvx, uzwx, vzwy
    real(rprec) :: udpdx, udpdy, udpdz, vdpdx, vdpdy, vdpdz, wdpdx, wdpdy, wdpdz
    real(rprec) :: pdudx, pdudy, pdudz, pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz
    real(rprec) :: lapu, lapv, lapw
    real(rprec) :: ulapu, ulapv, ulapw, vlapu, vlapv, vlapw, wlapu, wlapv, wlapw
end type tavg_budget_t

type budget_t
    real(rprec) :: advxx, advyy, advzz, advxy, advxz, advyz
    real(rprec) :: tflucxx, tflucyy, tfluczz, tflucxy, tflucxz, tflucyz
    real(rprec) :: tpresxx, tpresyy, tpreszz, tpresxy, tpresxz, tpresyz
    real(rprec) :: pstrainxx, pstrainyy, pstrainzz, pstrainxy, pstrainxz, pstrainyz
    real(rprec) :: tviscxx, tviscyy, tvisczz, tviscxy, tviscxz, tviscyz
    real(rprec) :: prodxx, prodyy, prodzz, prodxy, prodxz, prodyz
    real(rprec) :: pdissxx, pdissyy, pdisszz, pdissxy, pdissxz, pdissyz
end type budget_t
#endif

#ifdef PPOUTPUT_TURBSPEC
type tavg_turbspec_t
    complex(rprec) :: uf, vf, wf, vortxf, vortyf, vortzf
    complex(rprec) :: uu, vv, ww, uv, uw, vw 
    complex(rprec) :: vortx2, vorty2, vortz2
    !real(rprec) :: vel2, vort2
end type tavg_turbspec_t

type turbspec_t
    real(rprec) :: upup, vpvp, wpwp, upvp, upwp, vpwp
    real(rprec) :: vortxp2, vortyp2, vortzp2
end type turbspec_t
#endif

#ifdef PPOUTPUT_SPECBUDG
type tavg_specbudg_t
    complex(rprec) :: p, dpdx, dpdy, dpdz
    complex(rprec) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    complex(rprec) :: uh_dudxh, uh_dudyh, uh_dudzh, uh_dvdxh, uh_dvdyh, uh_dvdzh, uh_dwdxh, uh_dwdyh, uh_dwdzh
    complex(rprec) :: vh_dudxh, vh_dudyh, vh_dudzh, vh_dvdxh, vh_dvdyh, vh_dvdzh, vh_dwdxh, vh_dwdyh, vh_dwdzh
    complex(rprec) :: wh_dudxh, wh_dudyh, wh_dudzh, wh_dvdxh, wh_dvdyh, wh_dvdzh, wh_dwdxh, wh_dwdyh, wh_dwdzh
    complex(rprec) :: ududxh, vdudyh, wdudzh, udvdxh, vdvdyh, wdvdzh, udwdxh, vdwdyh, wdwdzh 
    complex(rprec) :: uh_ududxh, uh_vdudyh, uh_wdudzh, vh_ududxh, vh_vdudyh, vh_wdudzh, wh_ududxh, wh_vdudyh, wh_wdudzh
    complex(rprec) :: uh_udvdxh, uh_vdvdyh, uh_wdvdzh, vh_udvdxh, vh_vdvdyh, vh_wdvdzh, wh_udvdxh, wh_vdvdyh, wh_wdvdzh
    complex(rprec) :: uh_udwdxh, uh_vdwdyh, uh_wdwdzh, vh_udwdxh, vh_vdwdyh, vh_wdwdzh, wh_udwdxh, wh_vdwdyh, wh_wdwdzh
    complex(rprec) :: pdudx, pdudy, pdudz, pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz
    complex(rprec) :: udpdx, udpdy, udpdz, vdpdx, vdpdy, vdpdz, wdpdx, wdpdy, wdpdz
    complex(rprec) :: uxux, uyuy, uzuz, vxvx, vyvy, vzvz, wxwx, wywy, wzwz
    complex(rprec) :: uxvx, uyvy, uzvz, uxwx, uywy, uzwz, vxwx, vywy, vzwz
    complex(rprec) :: lapu, lapv, lapw
    complex(rprec) :: ulapu, ulapv, ulapw, vlapu, vlapv, vlapw, wlapu, wlapv, wlapw
end type tavg_specbudg_t

type specbudg_t
    real(rprec) :: advxx, advyy, advzz, advxy, advxz, advyz
    real(rprec) :: tflucxx, tflucyy, tfluczz, tflucxy, tflucxz, tflucyz
    real(rprec) :: tpresxx, tpresyy, tpreszz, tpresxy, tpresxz, tpresyz
    real(rprec) :: pstrainxx, pstrainyy, pstrainzz, pstrainxy, pstrainxz, pstrainyz
    real(rprec) :: tviscxx, tviscyy, tvisczz, tviscxy, tviscxz, tviscyz
    real(rprec) :: prodxx, prodyy, prodzz, prodxy, prodxz, prodyz
    real(rprec) :: pdissxx, pdissyy, pdisszz, pdissxy, pdissxz, pdissyz
end type specbudg_t
#endif

#ifdef PPOUTPUT_WMLES
real(rprec) :: tavg_wmles_total_time, tavg_wmles_dt
type tavg_wmles_t
    real(rprec) :: u, v, nu, uu, vv, uv
end type tavg_wmles_t
#endif

! Types for including wind turbines as drag disks
#ifdef PPTURBINES
! Single turbines
type turbine_t
    real(rprec) :: xloc, yloc, height, dia, thk
    ! term used for volume correction
    ! real(rprec) :: vol_c
    ! angle CCW(from above) from -x direction [degrees]
    real(rprec) :: theta1
    ! angle above the horizontal, from -x dir [degrees]
    real(rprec) :: theta2
    ! number of nodes associated with each turbine
    integer :: num_nodes
    ! location of turbine center (local k)
    integer :: icp, jcp, kcp
    ! true if the center is in the processor
    logical :: center_in_proc
    ! thrust coefficient
    real(rprec) :: Ct_prime
    ! running time-average of mean disk velocity
    real(rprec) :: u_d, u_d_T
    ! normal force on turbine disk
    real(rprec) :: f_n
    ! (nx,ny,nz) of unit normal for each turbine
    real(rprec), dimension(3) :: nhat
    ! indicator function - weighting of each node
    real(rprec), dimension(50000) :: ind
    ! object to calculate indicator function
    type(turb_ind_func_t) :: turb_ind_func
    ! (i,j,k) of each included node
    integer, dimension(50000,3) :: nodes
    ! search area for nearby nodes
    integer, dimension(6) :: nodes_max
end type turbine_t

! A collection of wind turbines
type wind_farm_t
    type(turbine_t), pointer, dimension(:) :: turbine
end type wind_farm_t

! The wind farm
type(wind_farm_t) :: wind_farm
#endif

! Create types for outputting data (instantaneous or averaged)
type(point_t), allocatable, dimension(:) :: point
type(plane_t), allocatable, dimension(:) :: xplane, yplane
type(zplane_t), allocatable, dimension(:) :: zplane

type(tavg_t), allocatable, dimension(:,:,:) :: tavg
type(tavg_t), allocatable, dimension(:) :: tavg_zplane

type(rs_t), allocatable, dimension(:,:,:) :: rs
type(rs_t), allocatable, dimension(:) :: rs_zplane, cnpy_zplane

#ifdef PPSCALARS
type(tavg_scal_t), allocatable, dimension(:,:,:) :: tavg_scal
type(rs_scal_t), allocatable, dimension(:,:,:) :: rs_scal
#ifdef PPOUTPUT_TURBSPEC
type(tavg_scal_turbspec_t), allocatable, dimension(:,:,:) :: tavg_scal_turbspecx
type(tavg_scal_turbspec_t), allocatable, dimension(:,:,:) :: tavg_scal_turbspecy
type(scal_turbspec_t), allocatable, dimension(:,:,:) :: scal_turbspecx
type(scal_turbspec_t), allocatable, dimension(:,:,:) :: scal_turbspecy
#endif
#endif

type(tavg_vort_t), allocatable, dimension(:,:,:) :: tavg_vort
type(vortrms_t), allocatable, dimension(:,:,:) :: vortrms

#ifdef PPOUTPUT_SGS
type(tavg_sgs_t), allocatable, dimension(:,:,:) :: tavg_sgs
#endif

#ifdef PPOUTPUT_BUDGET
type(tavg_budget_t), allocatable, dimension(:,:,:) :: tavg_budget
type(budget_t), allocatable, dimension(:,:,:) :: budget
#endif

#ifdef PPOUTPUT_TURBSPEC
type(tavg_turbspec_t), allocatable, dimension(:,:,:) :: tavg_turbspecx
type(tavg_turbspec_t), allocatable, dimension(:,:,:) :: tavg_turbspecy
type(turbspec_t), allocatable, dimension(:,:,:) :: turbspecx
type(turbspec_t), allocatable, dimension(:,:,:) :: turbspecy
#endif

#ifdef PPOUTPUT_SPECBUDG
type(tavg_specbudg_t), allocatable, dimension(:,:,:) :: tavg_specbudgx
type(specbudg_t), allocatable, dimension(:,:,:) :: specbudgx
type(tavg_specbudg_t), allocatable, dimension(:,:,:) :: tavg_specbudgy
type(specbudg_t), allocatable, dimension(:,:,:) :: specbudgy
#endif

#ifdef PPOUTPUT_WMLES
type(tavg_wmles_t), allocatable, dimension(:,:,:) :: tavg_wmles
#endif

contains

#ifdef PPTURBINES
!*****************************************************************************
function val(this, r, x) result(Rval)
!*****************************************************************************
use functions, only : linear_interp
implicit none
class(turb_ind_func_t), intent(in) :: this
real(rprec), intent(in) :: r, x
real(rprec) :: R1, R23, Rval

R23 = linear_interp(this%r, this%R23, r)
R1 = erf(this%sqrt6overdelta*(this%t_half + x)) +                         &
    erf(this%sqrt6overdelta*(this%t_half - x))
Rval = 0.5 * R1 * R23

end function val

!*****************************************************************************
subroutine init(this, delta2, thk, dia, N)
!*****************************************************************************
use param, only : write_endian, path, pi
use functions, only : bilinear_interp
implicit none
include'fftw3.f'

class(turb_ind_func_t), intent(inout) :: this
real(rprec), intent(in) :: delta2, thk, dia
integer, intent(in) :: N

real(rprec) :: L, d, R
integer, dimension(:), allocatable :: ind
real(rprec), dimension(:), allocatable :: yz
real(rprec), dimension(:,:), allocatable :: g, f, h
real(rprec), dimension(:), allocatable :: xi
real(rprec) :: dr, Lr
integer :: i, j

integer*8 plan
complex(rprec), dimension(:,:), allocatable :: ghat, fhat, hhat

L = 4 * dia
d = L / N
R = 0.5 * dia;

allocate(yz(N))
allocate(ind(N))
allocate(g(N, N))
allocate(h(N, N))
allocate(f(N, N))
allocate(ghat(N/2+1, N))
allocate(hhat(N/2+1, N))
allocate(fhat(N/2+1, N))

! Calculate constants
this%t_half = 0.5 * thk
this%sqrt6overdelta = sqrt(6._rprec) / sqrt(delta2)

! Calculate yz and indices to sort the result
do i = 1, N/2
    yz(i) = d*(i-0.5)
    ind(i) = N/2+i
end do
do i = N/2+1, N
    yz(i) = -L + d*(i-0.5)
    ind(i) = i-N/2
end do

! Calculate g and f
do j = 1, N
    do i = 1, N
        g(i,j) = exp(-6*(yz(i)**2+yz(j)**2)/delta2)
        if (sqrt(yz(i)**2 + yz(j)**2) < R) then
            h(i,j) = 1.0
        else
            h(i,j) = 0.0
        end if
    end do
end do

! Do the convolution f = g*h in fourier space
call dfftw_plan_dft_r2c_2d(plan, N, N, g, ghat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, g, ghat)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan, N, N, h, hhat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, h, hhat)
call dfftw_destroy_plan(plan)

fhat = ghat*hhat

! Compute the inverse fft of fhat
call dfftw_plan_dft_c2r_2d(plan, N, N, fhat, f, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fhat, f)
call dfftw_destroy_plan(plan)

! Normalize
f = f / N**2 * d**2

! Sort the results
f = f(ind,ind)
yz = yz(ind);

! Interpolate onto the lookup table
allocate(xi(N))
if (allocated(this%r) ) then
    deallocate(this%r)
end if
allocate( this%r(N) )
allocate( this%R23(N) )

Lr = R + 2 * sqrt(delta2)
dr = Lr / (N - 1)
do i = 1,N
    this%r(i) = (i-1)*dr
    xi(i) = 0
end do
this%R23 = bilinear_interp(yz, yz, f, xi, this%r)
this%R23 = this%R23 / this%R23(1)

end subroutine init
#endif

!/////////////////////////////////////////////////////////////////////////////
!/// Spectral RS operators
!/////////////////////////////////////////////////////////////////////////////

!*****************************************************************************
function rs_compute( a , lbz2) result(c)
!*****************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: a
type(rs_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % up2 = a % u2 - a % u * a % u
c % vp2 = a % v2 - a % v * a % v
c % wp2 = a % w2 - a % w * a % w
c % upvp = a % uv - a % u * a % v
!! using u_w and v_w below instead of u and v ensures that the Reynolds
!! stresses are on the same grid as the squared velocities (i.e., w-grid)
c % upwp = a % uw - a % u_w * a % w   !!pj
c % vpwp = a % vw - a % v_w * a % w   !!pj

end function rs_compute

#ifdef PPSCALARS
!*******************************************************************************
function rs_scal_compute( a, b, lbz2) result(c)
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_scal_t), dimension(:,:,lbz2:), intent(in) :: a
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: b
type(rs_scal_t), allocatable, dimension(:,:,:) :: c
integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % thetap2 = a % theta2 - a % theta * a % theta
c % upthetap = a % utheta - b % u * a % theta
c % vpthetap = a % vtheta - b % v * a % theta
c % wpthetap = a % wtheta - b % w * a % theta

end function rs_scal_compute

#ifdef PPOUTPUT_TURBSPEC
!*******************************************************************************
function scal_turbspec_compute( a , b, lbz2 ) result( c )
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_scal_turbspec_t), dimension(:,:,lbz2:), intent(in) :: a
type(tavg_turbspec_t), dimension(:,:,lbz2:), intent(in) :: b
type(scal_turbspec_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % thetap2  = a % theta2 - real( a % thetaf * conjg( a % thetaf ) )
c % upthetap = a % utheta - real( b % uf * conjg( a % thetaf ) )
c % vpthetap = a % vtheta - real( b % vf * conjg( a % thetaf ) )
c % wpthetap = a % wtheta - real( b % wf * conjg( a % thetaf ) )

end function scal_turbspec_compute
#endif
#endif

!*******************************************************************************
function vortrms_compute( a , lbz2) result(c)
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_vort_t), dimension(:,:,lbz2:), intent(in) :: a
type(vortrms_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % vortxrms = a % vortx2 - a % vortx * a % vortx
c % vortyrms = a % vorty2 - a % vorty * a % vorty
c % vortzrms = a % vortz2 - a % vortz * a % vortz

end function vortrms_compute

#ifdef PPOUTPUT_BUDGET
!*****************************************************************************
function budget_compute( a, b, lbz2) result(c)
!*****************************************************************************
use param, only: nu_molec
implicit none
integer, intent(in) :: lbz2
type(tavg_budget_t), dimension(:,:,lbz2:), intent(in) :: a
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: b
type(budget_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz, i, j, k
real(rprec) :: Rxx, Ryy, Rzz, Rxy, Rxz, Ryz
real(rprec) :: Cududx, Cududy, Cududz, Cudvdx, Cudvdy, Cudvdz, Cudwdx, Cudwdy, Cudwdz
real(rprec) :: Cvdudx, Cvdudy, Cvdudz, Cvdvdx, Cvdvdy, Cvdvdz, Cvdwdx, Cvdwdy, Cvdwdz 
real(rprec) :: Cwdudx, Cwdudy, Cwdudz, Cwdvdx, Cwdvdy, Cwdvdz, Cwdwdx, Cwdwdy, Cwdwdz 
real(rprec) :: Cuududx,Cuvdudy,Cuwdudz,Cuudvdx,Cuvdvdy,Cuwdvdz
real(rprec) :: Cuudwdx, Cuvdwdy, Cuwdwdz
real(rprec) :: Cvududx, Cvvdudy, Cvwdudz, Cvudvdx, Cvvdvdy, Cvwdvdz
real(rprec) :: Cvudwdx, Cvvdwdy, Cvwdwdz
real(rprec) :: Cwududx, Cwvdudy, Cwwdudz, Cwudvdx, Cwvdvdy, Cwwdvdz
real(rprec) :: Cwudwdx, Cwvdwdy, Cwwdwdz
real(rprec) :: Cudpdx,Cudpdy, Cudpdz, Cvdpdx, Cvdpdy, Cvdpdz, Cwdpdx, Cwdpdy, Cwdpdz
real(rprec) :: Cpdudx,Cpdudy, Cpdudz, Cpdvdx, Cpdvdy, Cpdvdz, Cpdwdx, Cpdwdy, Cpdwdz
real(rprec) :: Cuxux, Cuyuy, Cuzuz, Cvxvx, Cvyvy, Cvzvz, Cwxwx, Cwywy, Cwzwz
real(rprec) :: Cuxvx, Cuyvy, Cuzvz, Cuxwx, Cuywy, Cuzwz, Cvxwx, Cvywy, Cvzwz
real(rprec) :: Culapu,Culapv, Culapw, Cvlapu, Cvlapv, Cvlapw, Cwlapu, Cwlapv, Cwlapw

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

! -------------------------- Reynolds Stress Budget --------------------------

do i = 1, ubx
do j = 1, uby
do k = 1, ubz

! Compute intermediate terms used in the budget
! Velocity-Velocity correlation, Reynolds stress
! Similar to rs_compute, except all are on w-grid
Rxx = a(i,j,k) % uu - b(i,j,k) % u_w * b(i,j,k) % u_w
Ryy = a(i,j,k) % vv - b(i,j,k) % v_w * b(i,j,k) % v_w
Rzz = a(i,j,k) % ww - b(i,j,k) % w   * b(i,j,k) % w
Rxy = a(i,j,k) % uv - b(i,j,k) % u_w * b(i,j,k) % v_w
Rxz = a(i,j,k) % uw - b(i,j,k) % u_w * b(i,j,k) % w
Ryz = a(i,j,k) % vw - b(i,j,k) % v_w * b(i,j,k) % w

! Velocity-Velocity Gradient correlation, ui*dujdxk
Cududx = a(i,j,k)%ududx - b(i,j,k)%u_w * a(i,j,k)%dudx
Cududy = a(i,j,k)%ududy - b(i,j,k)%u_w * a(i,j,k)%dudy
Cududz = a(i,j,k)%ududz - b(i,j,k)%u_w * a(i,j,k)%dudz
Cudvdx = a(i,j,k)%udvdx - b(i,j,k)%u_w * a(i,j,k)%dvdx
Cudvdy = a(i,j,k)%udvdy - b(i,j,k)%u_w * a(i,j,k)%dvdy
Cudvdz = a(i,j,k)%udvdz - b(i,j,k)%u_w * a(i,j,k)%dvdz
Cudwdx = a(i,j,k)%udwdx - b(i,j,k)%u_w * a(i,j,k)%dwdx
Cudwdy = a(i,j,k)%udwdy - b(i,j,k)%u_w * a(i,j,k)%dwdy
Cudwdz = a(i,j,k)%udwdz - b(i,j,k)%u_w * a(i,j,k)%dwdz

Cvdudx = a(i,j,k)%vdudx - b(i,j,k)%v_w * a(i,j,k)%dudx
Cvdudy = a(i,j,k)%vdudy - b(i,j,k)%v_w * a(i,j,k)%dudy
Cvdudz = a(i,j,k)%vdudz - b(i,j,k)%v_w * a(i,j,k)%dudz
Cvdvdx = a(i,j,k)%vdvdx - b(i,j,k)%v_w * a(i,j,k)%dvdx
Cvdvdy = a(i,j,k)%vdvdy - b(i,j,k)%v_w * a(i,j,k)%dvdy
Cvdvdz = a(i,j,k)%vdvdz - b(i,j,k)%v_w * a(i,j,k)%dvdz
Cvdwdx = a(i,j,k)%vdwdx - b(i,j,k)%v_w * a(i,j,k)%dwdx
Cvdwdy = a(i,j,k)%vdwdy - b(i,j,k)%v_w * a(i,j,k)%dwdy
Cvdwdz = a(i,j,k)%vdwdz - b(i,j,k)%v_w * a(i,j,k)%dwdz

Cwdudx = a(i,j,k)%wdudx - b(i,j,k)%w * a(i,j,k)%dudx
Cwdudy = a(i,j,k)%wdudy - b(i,j,k)%w * a(i,j,k)%dudy
Cwdudz = a(i,j,k)%wdudz - b(i,j,k)%w * a(i,j,k)%dudz
Cwdvdx = a(i,j,k)%wdvdx - b(i,j,k)%w * a(i,j,k)%dvdx
Cwdvdy = a(i,j,k)%wdvdy - b(i,j,k)%w * a(i,j,k)%dvdy
Cwdvdz = a(i,j,k)%wdvdz - b(i,j,k)%w * a(i,j,k)%dvdz
Cwdwdx = a(i,j,k)%wdwdx - b(i,j,k)%w * a(i,j,k)%dwdx
Cwdwdy = a(i,j,k)%wdwdy - b(i,j,k)%w * a(i,j,k)%dwdy
Cwdwdz = a(i,j,k)%wdwdz - b(i,j,k)%w * a(i,j,k)%dwdz

! vel-vel-velGrad triple correlation, ui*uk*dujdxk
Cuududx = a(i,j,k)%uududx - b(i,j,k)%u_w*b(i,j,k)%u_w*a(i,j,k)%dudx -        &
    Rxx*a(i,j,k)%dudx - Cududx*b(i,j,k)%u_w - Cududx*b(i,j,k)%u_w
Cuvdudy = a(i,j,k)%uvdudy - b(i,j,k)%u_w*b(i,j,k)%v_w*a(i,j,k)%dudy -        &
    Rxy*a(i,j,k)%dudy - Cududy*b(i,j,k)%v_w - Cvdudy*b(i,j,k)%u_w
Cuwdudz = a(i,j,k)%uwdudz - b(i,j,k)%u_w*b(i,j,k)%w*a(i,j,k)%dudz   -        &
    Rxz*a(i,j,k)%dudz - Cududz*b(i,j,k)%w   - Cwdudz*b(i,j,k)%u_w
Cuudvdx = a(i,j,k)%uudvdx - b(i,j,k)%u_w*b(i,j,k)%u_w*a(i,j,k)%dvdx -        &
    Rxx*a(i,j,k)%dvdx - Cudvdx*b(i,j,k)%u_w - Cudvdx*b(i,j,k)%u_w
Cuvdvdy = a(i,j,k)%uvdvdy - b(i,j,k)%u_w*b(i,j,k)%v_w*a(i,j,k)%dvdy -        &
    Rxy*a(i,j,k)%dvdy - Cudvdy*b(i,j,k)%v_w - Cvdvdy*b(i,j,k)%u_w
Cuwdvdz = a(i,j,k)%uwdvdz - b(i,j,k)%u_w*b(i,j,k)%w*a(i,j,k)%dvdz   -        &
    Rxz*a(i,j,k)%dvdz - Cudvdz*b(i,j,k)%w   - Cwdvdz*b(i,j,k)%u_w
Cuudwdx = a(i,j,k)%uudwdx - b(i,j,k)%u_w*b(i,j,k)%u_w*a(i,j,k)%dwdx -        &
    Rxx*a(i,j,k)%dwdx - Cudwdx*b(i,j,k)%u_w - Cudwdx*b(i,j,k)%u_w
Cuvdwdy = a(i,j,k)%uvdwdy - b(i,j,k)%u_w*b(i,j,k)%v_w*a(i,j,k)%dwdy -        &
    Rxy*a(i,j,k)%dwdy - Cudwdy*b(i,j,k)%v_w - Cvdwdy*b(i,j,k)%u_w
Cuwdwdz = a(i,j,k)%uwdwdz - b(i,j,k)%u_w*b(i,j,k)%w*a(i,j,k)%dwdz   -        &
    Rxz*a(i,j,k)%dwdz - Cudwdz*b(i,j,k)%w   - Cwdwdz*b(i,j,k)%u_w

Cvududx = a(i,j,k)%vududx - b(i,j,k)%v_w*b(i,j,k)%u_w*a(i,j,k)%dudx -        &
    Rxy*a(i,j,k)%dudx - Cvdudx*b(i,j,k)%u_w - Cududx*b(i,j,k)%v_w
Cvvdudy = a(i,j,k)%vvdudy - b(i,j,k)%v_w*b(i,j,k)%v_w*a(i,j,k)%dudy -        &
    Ryy*a(i,j,k)%dudy - Cvdudy*b(i,j,k)%v_w - Cvdudy*b(i,j,k)%v_w
Cvwdudz = a(i,j,k)%vwdudz - b(i,j,k)%v_w*b(i,j,k)%w*a(i,j,k)%dudz   -        &
    Ryz*a(i,j,k)%dudz - Cvdudz*b(i,j,k)%w   - Cwdudz*b(i,j,k)%v_w
Cvudvdx = a(i,j,k)%vudvdx - b(i,j,k)%v_w*b(i,j,k)%u_w*a(i,j,k)%dvdx -        &
    Rxy*a(i,j,k)%dvdx - Cvdvdx*b(i,j,k)%u_w - Cudvdx*b(i,j,k)%v_w
Cvvdvdy = a(i,j,k)%vvdvdy - b(i,j,k)%v_w*b(i,j,k)%v_w*a(i,j,k)%dvdy -        &
    Ryy*a(i,j,k)%dvdy - Cvdvdy*b(i,j,k)%v_w - Cvdvdy*b(i,j,k)%v_w
Cvwdvdz = a(i,j,k)%vwdvdz - b(i,j,k)%v_w*b(i,j,k)%w*a(i,j,k)%dvdz   -        &
    Ryz*a(i,j,k)%dvdz - Cvdvdz*b(i,j,k)%w   - Cwdvdz*b(i,j,k)%v_w
Cvudwdx = a(i,j,k)%vudwdx - b(i,j,k)%v_w*b(i,j,k)%u_w*a(i,j,k)%dwdx -        &
    Rxy*a(i,j,k)%dwdx - Cvdwdx*b(i,j,k)%u_w - Cudwdx*b(i,j,k)%v_w
Cvvdwdy = a(i,j,k)%vvdwdy - b(i,j,k)%v_w*b(i,j,k)%v_w*a(i,j,k)%dwdy -        &
    Ryy*a(i,j,k)%dwdy - Cvdwdy*b(i,j,k)%v_w - Cvdwdy*b(i,j,k)%v_w
Cvwdwdz = a(i,j,k)%vwdwdz - b(i,j,k)%v_w*b(i,j,k)%w*a(i,j,k)%dwdz   -        &
    Ryz*a(i,j,k)%dwdz - Cvdwdz*b(i,j,k)%w   - Cwdwdz*b(i,j,k)%v_w

Cwududx = a(i,j,k)%wududx - b(i,j,k)%w*b(i,j,k)%u_w*a(i,j,k)%dudx -          &
    Rxz*a(i,j,k)%dudx - Cwdudx*b(i,j,k)%u_w - Cududx*b(i,j,k)%w
Cwvdudy = a(i,j,k)%wvdudy - b(i,j,k)%w*b(i,j,k)%v_w*a(i,j,k)%dudy -          &
    Ryz*a(i,j,k)%dudy - Cwdudy*b(i,j,k)%v_w - Cvdudy*b(i,j,k)%w
Cwwdudz = a(i,j,k)%wwdudz - b(i,j,k)%w*b(i,j,k)%w*a(i,j,k)%dudz   -          &
    Rzz*a(i,j,k)%dudz - Cwdudz*b(i,j,k)%w   - Cwdudz*b(i,j,k)%w
Cwudvdx = a(i,j,k)%wudvdx - b(i,j,k)%w*b(i,j,k)%u_w*a(i,j,k)%dvdx -          &
    Rxz*a(i,j,k)%dvdx - Cwdvdx*b(i,j,k)%u_w - Cudvdx*b(i,j,k)%w
Cwvdvdy = a(i,j,k)%wvdvdy - b(i,j,k)%w*b(i,j,k)%v_w*a(i,j,k)%dvdy -          &
    Ryz*a(i,j,k)%dvdy - Cwdvdy*b(i,j,k)%v_w - Cvdvdy*b(i,j,k)%w
Cwwdvdz = a(i,j,k)%wwdvdz - b(i,j,k)%w*b(i,j,k)%w*a(i,j,k)%dvdz   -          &
    Rzz*a(i,j,k)%dvdz - Cwdvdz*b(i,j,k)%w   - Cwdvdz*b(i,j,k)%w
Cwudwdx = a(i,j,k)%wudwdx - b(i,j,k)%w*b(i,j,k)%u_w*a(i,j,k)%dwdx -          &
    Rxz*a(i,j,k)%dwdx - Cwdwdx*b(i,j,k)%u_w - Cudwdx*b(i,j,k)%w
Cwvdwdy = a(i,j,k)%wvdwdy - b(i,j,k)%w*b(i,j,k)%v_w*a(i,j,k)%dwdy -          &
    Ryz*a(i,j,k)%dwdy - Cwdwdy*b(i,j,k)%v_w - Cvdwdy*b(i,j,k)%w
Cwwdwdz = a(i,j,k)%wwdwdz - b(i,j,k)%w*b(i,j,k)%w*a(i,j,k)%dwdz   -          &
    Rzz*a(i,j,k)%dwdz - Cwdwdz*b(i,j,k)%w   - Cwdwdz*b(i,j,k)%w

! vel-presGrad correlation, ui*dpdxj
Cudpdx = a(i,j,k)%udpdx - b(i,j,k)%u_w * a(i,j,k)%dpdx
Cudpdy = a(i,j,k)%udpdy - b(i,j,k)%u_w * a(i,j,k)%dpdy
Cudpdz = a(i,j,k)%udpdz - b(i,j,k)%u_w * a(i,j,k)%dpdz
Cvdpdx = a(i,j,k)%vdpdx - b(i,j,k)%v_w * a(i,j,k)%dpdx
Cvdpdy = a(i,j,k)%vdpdy - b(i,j,k)%v_w * a(i,j,k)%dpdy
Cvdpdz = a(i,j,k)%vdpdz - b(i,j,k)%v_w * a(i,j,k)%dpdz
Cwdpdx = a(i,j,k)%wdpdx - b(i,j,k)%w * a(i,j,k)%dpdx
Cwdpdy = a(i,j,k)%wdpdy - b(i,j,k)%w * a(i,j,k)%dpdy
Cwdpdz = a(i,j,k)%wdpdz - b(i,j,k)%w * a(i,j,k)%dpdz

! pres-velGrad correlation, p*duidxj
Cpdudx = a(i,j,k)%pdudx - a(i,j,k)%p * a(i,j,k)%dudx
Cpdudy = a(i,j,k)%pdudy - a(i,j,k)%p * a(i,j,k)%dudy
Cpdudz = a(i,j,k)%pdudz - a(i,j,k)%p * a(i,j,k)%dudz
Cpdvdx = a(i,j,k)%pdvdx - a(i,j,k)%p * a(i,j,k)%dvdx
Cpdvdy = a(i,j,k)%pdvdy - a(i,j,k)%p * a(i,j,k)%dvdy
Cpdvdz = a(i,j,k)%pdvdz - a(i,j,k)%p * a(i,j,k)%dvdz
Cpdwdx = a(i,j,k)%pdwdx - a(i,j,k)%p * a(i,j,k)%dwdx
Cpdwdy = a(i,j,k)%pdwdy - a(i,j,k)%p * a(i,j,k)%dwdy
Cpdwdz = a(i,j,k)%pdwdz - a(i,j,k)%p * a(i,j,k)%dwdz

! velGrad-velGrad correlation, duidxk*dujdxk, i=j
Cuxux = a(i,j,k)%uxux - (a(i,j,k)%dudx**2)
Cuyuy = a(i,j,k)%uyuy - (a(i,j,k)%dudy**2)
Cuzuz = a(i,j,k)%uzuz - (a(i,j,k)%dudz**2)
Cvxvx = a(i,j,k)%vxvx - (a(i,j,k)%dvdx**2)
Cvyvy = a(i,j,k)%vyvy - (a(i,j,k)%dvdy**2)
Cvzvz = a(i,j,k)%vzvz - (a(i,j,k)%dvdz**2)
Cwxwx = a(i,j,k)%wxwx - (a(i,j,k)%dwdx**2)
Cwywy = a(i,j,k)%wywy - (a(i,j,k)%dwdy**2)
Cwzwz = a(i,j,k)%wzwz - (a(i,j,k)%dwdz**2)

! velGrad-velGrad correlation, duidxk*dujdxk, i/=j
Cuxvx = a(i,j,k)%uxvx - a(i,j,k)%dudx * a(i,j,k)%dvdx
Cuyvy = a(i,j,k)%uyvy - a(i,j,k)%dudy * a(i,j,k)%dvdy
Cuzvz = a(i,j,k)%uzvz - a(i,j,k)%dudz * a(i,j,k)%dvdz
Cuxwx = a(i,j,k)%uxwx - a(i,j,k)%dudx * a(i,j,k)%dwdx
Cuywy = a(i,j,k)%uywy - a(i,j,k)%dudy * a(i,j,k)%dwdy
Cuzwz = a(i,j,k)%uzwz - a(i,j,k)%dudz * a(i,j,k)%dwdz
Cvxwx = a(i,j,k)%vxwx - a(i,j,k)%dvdx * a(i,j,k)%dwdx
Cvywy = a(i,j,k)%vywy - a(i,j,k)%dvdy * a(i,j,k)%dwdy
Cvzwz = a(i,j,k)%vzwz - a(i,j,k)%dvdz * a(i,j,k)%dwdz

! vel-Laplacian correlation, ui*lap(uj)
Culapu = a(i,j,k)%ulapu - b(i,j,k)%u_w * a(i,j,k)%lapu
Culapv = a(i,j,k)%ulapv - b(i,j,k)%u_w * a(i,j,k)%lapv
Culapw = a(i,j,k)%ulapw - b(i,j,k)%u_w * a(i,j,k)%lapw
Cvlapu = a(i,j,k)%vlapu - b(i,j,k)%v_w * a(i,j,k)%lapu
Cvlapv = a(i,j,k)%vlapv - b(i,j,k)%v_w * a(i,j,k)%lapv
Cvlapw = a(i,j,k)%vlapw - b(i,j,k)%v_w * a(i,j,k)%lapw
Cwlapu = a(i,j,k)%wlapu - b(i,j,k)%w * a(i,j,k)%lapu
Cwlapv = a(i,j,k)%wlapv - b(i,j,k)%w * a(i,j,k)%lapv
Cwlapw = a(i,j,k)%wlapw - b(i,j,k)%w * a(i,j,k)%lapw

! Now compute terms to be outputted
! Advection, Uk*dRijdxk
c(i,j,k) % advxx = 2.0_rprec*( b(i,j,k)%u_w*Cududx + b(i,j,k)%v_w*Cududy + b(i,j,k)%w*Cududz)
c(i,j,k) % advyy = 2.0_rprec*( b(i,j,k)%u_w*Cvdvdx + b(i,j,k)%v_w*Cvdvdy + b(i,j,k)%w*Cvdvdz)
c(i,j,k) % advzz = 2.0_rprec*( b(i,j,k)%u_w*Cwdwdx + b(i,j,k)%v_w*Cwdwdy + b(i,j,k)%w*Cwdwdz)
c(i,j,k) % advxy = (b(i,j,k)%u_w*Cudvdx + b(i,j,k)%v_w*Cudvdy + b(i,j,k)%w*Cudvdz)            &
    + (b(i,j,k)%u_w*Cvdudx + b(i,j,k)%v_w*Cvdudy + b(i,j,k)%w*Cvdudz)
c(i,j,k) % advxz = (b(i,j,k)%u_w*Cudwdx + b(i,j,k)%v_w*Cudwdy + b(i,j,k)%w*Cudwdz)            &
    + (b(i,j,k)%u_w*Cwdudx + b(i,j,k)%v_w*Cwdudy + b(i,j,k)%w*Cwdudz)
c(i,j,k) % advyz = (b(i,j,k)%u_w*Cvdwdx + b(i,j,k)%v_w*Cvdwdy + b(i,j,k)%w*Cvdwdz)            &
    + (b(i,j,k)%u_w*Cwdvdx + b(i,j,k)%v_w*Cwdvdy + b(i,j,k)%w*Cwdvdz)

! Transport by fluctuations, d(ui*uj*uk)dxk = ui*uk*dujdxk + uj*uk*duidxk
c(i,j,k) % tflucxx = 2.0_rprec*(Cuududx + Cuvdudy + Cuwdudz)
c(i,j,k) % tflucyy = 2.0_rprec*(Cvudvdx + Cvvdvdy + Cvwdvdz)
c(i,j,k) % tfluczz = 2.0_rprec*(Cwudwdx + Cwvdwdy + Cwwdwdz)
c(i,j,k) % tflucxy = (Cuudvdx + Cuvdvdy + Cuwdvdz) + (Cvududx + Cvvdudy + Cvwdudz)
c(i,j,k) % tflucxz = (Cuudwdx + Cuvdwdy + Cuwdwdz) + (Cwududx + Cwvdudy + Cwwdudz)
c(i,j,k) % tflucyz = (Cvudwdx + Cvvdwdy + Cvwdwdz) + (Cwudvdx + Cwvdvdy + Cwwdvdz)

! Transport by pressure-velocity, d(p*(uj*Iik+ui*Ijk))dxk = d(p*uj)dxi + d(p*ui)dxj
c(i,j,k) % tpresxx = 2.0_rprec*(Cpdudx + Cudpdx)
c(i,j,k) % tpresyy = 2.0_rprec*(Cpdvdy + Cvdpdy)
c(i,j,k) % tpreszz = 2.0_rprec*(Cpdwdz + Cwdpdz)
c(i,j,k) % tpresxy = (Cpdvdx + Cvdpdx) + (Cpdudy + Cudpdy)
c(i,j,k) % tpresxz = (Cpdwdx + Cwdpdx) + (Cpdudz + Cudpdz)
c(i,j,k) % tpresyz = (Cpdwdy + Cwdpdy) + (Cpdvdz + Cvdpdz)

! Pressure-Strain, 2*p*sij = p*(dujdxi + duidxj)
c(i,j,k) % pstrainxx = 2.0_rprec*Cpdudx
c(i,j,k) % pstrainyy = 2.0_rprec*Cpdvdy
c(i,j,k) % pstrainzz = 2.0_rprec*Cpdwdz
c(i,j,k) % pstrainxy = Cpdudy + Cpdvdx
c(i,j,k) % pstrainxz = Cpdudz + Cpdwdx
c(i,j,k) % pstrainyz = Cpdvdz + Cpdwdy

! Production Rate, -(Rik*dujdxk + Rjk*duidxk)
c(i,j,k) % prodxx = -2.0_rprec*(Rxx*a(i,j,k)%dudx + Rxy*a(i,j,k)%dudy + Rxz*a(i,j,k)%dudz)
c(i,j,k) % prodyy = -2.0_rprec*(Rxy*a(i,j,k)%dvdx + Ryy*a(i,j,k)%dvdy + Ryz*a(i,j,k)%dvdz)
c(i,j,k) % prodzz = -2.0_rprec*(Rxz*a(i,j,k)%dwdx + Ryz*a(i,j,k)%dwdy + Rzz*a(i,j,k)%dwdz)
c(i,j,k) % prodxy = -(Rxx*a(i,j,k)%dvdx + Rxy*a(i,j,k)%dvdy + Rxz*a(i,j,k)%dvdz +                 &
    Rxy*a(i,j,k)%dudx + Ryy*a(i,j,k)%dudy + Ryz*a(i,j,k)%dudz)
c(i,j,k) % prodxz = -(Rxx*a(i,j,k)%dwdx + Rxy*a(i,j,k)%dwdy + Rxz*a(i,j,k)%dwdz +                 &
    Rxz*a(i,j,k)%dudx + Ryz*a(i,j,k)%dudy + Rzz*a(i,j,k)%dudz)
c(i,j,k) % prodyz = -(Rxy*a(i,j,k)%dwdx + Ryy*a(i,j,k)%dwdy + Ryz*a(i,j,k)%dwdz +                 &
    Rxz*a(i,j,k)%dvdx + Ryz*a(i,j,k)%dvdy + Rzz*a(i,j,k)%dvdz)

! Pseudo-dissipation, 2*nu*duidxk*dujdxk
c(i,j,k) % pdissxx = 2.0_rprec*nu_molec*(Cuxux + Cuyuy + Cuzuz)
c(i,j,k) % pdissyy = 2.0_rprec*nu_molec*(Cvxvx + Cvyvy + Cvzvz)
c(i,j,k) % pdisszz = 2.0_rprec*nu_molec*(Cwxwx + Cwywy + Cwzwz)
c(i,j,k) % pdissxy = 2.0_rprec*nu_molec*(Cuxvx + Cuyvy + Cuzvz)
c(i,j,k) % pdissxz = 2.0_rprec*nu_molec*(Cuxwx + Cuywy + Cuzwz)
c(i,j,k) % pdissyz = 2.0_rprec*nu_molec*(Cvxwx + Cvywy + Cvzwz)

! Transport by viscous diffusion
! Cuilapuj already multiplied by viscosity because divtxj was used
! Change sign of Cuilapuj because of how divtxj is set up
c(i,j,k) % tviscxx = c(i,j,k) % pdissxx - 2.0_rprec*Culapu
c(i,j,k) % tviscyy = c(i,j,k) % pdissyy - 2.0_rprec*Cvlapv
c(i,j,k) % tvisczz = c(i,j,k) % pdisszz - 2.0_rprec*Cwlapw
c(i,j,k) % tviscxy = c(i,j,k) % pdissxy - (Culapv + Cvlapu)
c(i,j,k) % tviscxz = c(i,j,k) % pdissxz - (Culapw + Cwlapu)
c(i,j,k) % tviscyz = c(i,j,k) % pdissyz - (Cvlapw + Cwlapv)

! ---------------------- Mean kinetic energy balance -------------------------

end do 
end do
end do

! ---------------------------- Prepare Output --------------------------------
! Move all terms to the RHS
c % advxx = - c % advxx
c % advyy = - c % advyy
c % advzz = - c % advzz
c % advxy = - c % advxy
c % advxz = - c % advxz
c % advyz = - c % advyz

c % tflucxx = - c % tflucxx
c % tflucyy = - c % tflucyy
c % tfluczz = - c % tfluczz
c % tflucxy = - c % tflucxy
c % tflucxz = - c % tflucxz
c % tflucyz = - c % tflucyz

c % tpresxx = - c % tpresxx
c % tpresyy = - c % tpresyy
c % tpreszz = - c % tpreszz
c % tpresxy = - c % tpresxy
c % tpresxz = - c % tpresxz
c % tpresyz = - c % tpresyz

c % pdissxx = - c % pdissxx
c % pdissyy = - c % pdissyy
c % pdisszz = - c % pdisszz
c % pdissxy = - c % pdissxy
c % pdissxz = - c % pdissxz
c % pdissyz = - c % pdissyz

end function budget_compute
#endif

#ifdef PPOUTPUT_TURBSPEC
!*****************************************************************************
function turbspec_compute( a, lbz2 ) result( c )
!*****************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_turbspec_t), dimension(:,:,lbz2:), intent(in) :: a
type(turbspec_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % upup = real( a % uu - a % uf * conjg( a % uf ) )
c % vpvp = real( a % vv - a % vf * conjg( a % vf ) )
c % wpwp = real( a % ww - a % wf * conjg( a % wf ) )

c % upvp = real( a % uv - a % uf * conjg( a % vf ) )
c % upwp = real( a % uw - a % uf * conjg( a % wf ) )
c % vpwp = real( a % vw - a % vf * conjg( a % wf ) )

c % vortxp2 = real( a % vortx2 - a % vortxf * conjg( a % vortxf ) )
c % vortyp2 = real( a % vorty2 - a % vortyf * conjg( a % vortyf ) )
c % vortzp2 = real( a % vortz2 - a % vortzf * conjg( a % vortzf ) )

end function turbspec_compute
#endif

#ifdef PPOUTPUT_SPECBUDG
!*****************************************************************************
function specbudgx_compute( a, b, bb, c, d, lbz2 ) result( e )
!*****************************************************************************
! 
! Compute the spectral budget. All quantities should be on the w-grid.
! 
use param, only: nu_molec
implicit none
integer, intent(in) :: lbz2
type(tavg_specbudg_t), dimension(:,:,lbz2:), intent(in) :: a
type(turbspec_t), dimension(:,:,lbz2:), intent(in) :: b
type(tavg_turbspec_t), dimension(:,:,lbz2:), intent(in) :: bb
type(tavg_budget_t), dimension(:,:,lbz2:), intent(in) :: c
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: d
type(specbudg_t), allocatable, dimension(:,:,:) :: e

integer :: ubx, uby, ubz, vbx, vby, vbz, i, j, k

real(rprec), allocatable, dimension(:,:,:) :: u_avg, v_avg, w_avg,          &
    dudx_avg, dvdx_avg, dwdx_avg, dudy_avg, dvdy_avg, dwdy_avg,             &
    dudz_avg, dvdz_avg, dwdz_avg

complex(rprec) :: Cududx, Cududy, Cududz, Cudvdx, Cudvdy, Cudvdz, Cudwdx, Cudwdy, Cudwdz,  &
    Cvdudx, Cvdudy, Cvdudz, Cvdvdx, Cvdvdy, Cvdvdz, Cvdwdx, Cvdwdy, Cvdwdz,                &
    Cwdudx, Cwdudy, Cwdudz, Cwdvdx, Cwdvdy, Cwdvdz, Cwdwdx, Cwdwdy, Cwdwdz 

complex(rprec) :: Cuududx, Cuvdudy, Cuwdudz, Cvududx, Cvvdudy, Cvwdudz, Cwududx, Cwvdudy, Cwwdudz, &
    Cuudvdx, Cuvdvdy, Cuwdvdz, Cvudvdx, Cvvdvdy, Cvwdvdz, Cwudvdx, Cwvdvdy, Cwwdvdz,               &
    Cuudwdx, Cuvdwdy, Cuwdwdz, Cvudwdx, Cvvdwdy, Cvwdwdz, Cwudwdx, Cwvdwdy, Cwwdwdz

complex(rprec) :: Cpdudx, Cpdudy, Cpdudz, Cpdvdx, Cpdvdy, Cpdvdz, &
    Cpdwdx, Cpdwdy, Cpdwdz

complex(rprec) :: Cudpdx, Cudpdy, Cudpdz, Cvdpdx, Cvdpdy, Cvdpdz, &
    Cwdpdx, Cwdpdy, Cwdpdz

complex(rprec) :: Cuxux, Cuyuy, Cuzuz, Cvxvx, Cvyvy, Cvzvz, Cwxwx, Cwywy, Cwzwz
complex(rprec) :: Cuxvx, Cuyvy, Cuzvz, Cuxwx, Cuywy, Cuzwz, Cvxwx, Cvywy, Cvzwz
complex(rprec) :: Culapu, Culapv, Culapw, Cvlapu, Cvlapv, Cvlapw, Cwlapu, Cwlapv, Cwlapw

! Size of a, b, and e
ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

! Size of c and d
vbx = ubound(c,1)
vby = ubound(c,2)
vbz = ubound(c,3)

allocate(e(ubx,uby,lbz2:ubz))

allocate(u_avg(ubx,vby,lbz2:vbz))
allocate(v_avg(ubx,vby,lbz2:vbz))
allocate(w_avg(ubx,vby,lbz2:vbz))
allocate(dudx_avg(ubx,vby,lbz2:vbz))
allocate(dudy_avg(ubx,vby,lbz2:vbz))
allocate(dudz_avg(ubx,vby,lbz2:vbz))
allocate(dvdx_avg(ubx,vby,lbz2:vbz))
allocate(dvdy_avg(ubx,vby,lbz2:vbz))
allocate(dvdz_avg(ubx,vby,lbz2:vbz))
allocate(dwdx_avg(ubx,vby,lbz2:vbz))
allocate(dwdy_avg(ubx,vby,lbz2:vbz))
allocate(dwdz_avg(ubx,vby,lbz2:vbz))

! -------------------------- Streamwise Average ------------------------------
! When deriving the spectral budget equation, there are transforms of the 
! time-averaged mean multiplied with other quantities. To avoid convolutions
! in the budget equation, it is assumed that the time-averaged mean is not a 
! function of the direction in which the transform is being applied.
! 
! Because of this, we average the time-averaged quantities in that direction
! before removing the time-averaged mean component from the instantaneous 
! running sums.
!
! Quantities used in average are in x,y,z space, however only kx=0 mode
! is filled, kx~=0 modes are left as zeros
u_avg(:,:,:) = 0._rprec
v_avg(:,:,:) = 0._rprec
w_avg(:,:,:) = 0._rprec
dudx_avg(:,:,:) = 0._rprec
dudy_avg(:,:,:) = 0._rprec
dudz_avg(:,:,:) = 0._rprec
dvdx_avg(:,:,:) = 0._rprec
dvdy_avg(:,:,:) = 0._rprec
dvdz_avg(:,:,:) = 0._rprec
dwdx_avg(:,:,:) = 0._rprec
dwdy_avg(:,:,:) = 0._rprec
dwdz_avg(:,:,:) = 0._rprec
do i = 1, vbx
    u_avg(1,:,:) = u_avg(1,:,:) + d(i,:,:) % u_w
    v_avg(1,:,:) = v_avg(1,:,:) + d(i,:,:) % v_w
    w_avg(1,:,:) = w_avg(1,:,:) + d(i,:,:) % w
    dudx_avg(1,:,:) = dudx_avg(1,:,:) + c(i,:,:) % dudx
    dudy_avg(1,:,:) = dudy_avg(1,:,:) + c(i,:,:) % dudy
    dudz_avg(1,:,:) = dudz_avg(1,:,:) + c(i,:,:) % dudz
    dvdx_avg(1,:,:) = dvdx_avg(1,:,:) + c(i,:,:) % dvdx
    dvdy_avg(1,:,:) = dvdy_avg(1,:,:) + c(i,:,:) % dvdy
    dvdz_avg(1,:,:) = dvdz_avg(1,:,:) + c(i,:,:) % dvdz
    dwdx_avg(1,:,:) = dwdx_avg(1,:,:) + c(i,:,:) % dwdx
    dwdy_avg(1,:,:) = dwdy_avg(1,:,:) + c(i,:,:) % dwdy
    dwdz_avg(1,:,:) = dwdz_avg(1,:,:) + c(i,:,:) % dwdz
enddo
u_avg(1,:,:) = u_avg(1,:,:) / vbx
v_avg(1,:,:) = v_avg(1,:,:) / vbx
w_avg(1,:,:) = w_avg(1,:,:) / vbx
dudx_avg(1,:,:) = dudx_avg(1,:,:) / vbx
dudy_avg(1,:,:) = dudy_avg(1,:,:) / vbx
dudz_avg(1,:,:) = dudz_avg(1,:,:) / vbx
dvdx_avg(1,:,:) = dvdx_avg(1,:,:) / vbx
dvdy_avg(1,:,:) = dvdy_avg(1,:,:) / vbx
dvdz_avg(1,:,:) = dvdz_avg(1,:,:) / vbx
dwdx_avg(1,:,:) = dwdx_avg(1,:,:) / vbx
dwdy_avg(1,:,:) = dwdy_avg(1,:,:) / vbx
dwdz_avg(1,:,:) = dwdz_avg(1,:,:) / vbx

do i = 1, ubx
do j = 1, uby
do k = 1, ubz

! ---------------------- Compute Intermediate Terms --------------------------
! Compute intermediate terms used in the budget
! Velocity-Velocity Gradient correlation, ui_hat*dujdxk_hat
Cududx = a(i,j,k)%uh_dudxh - bb(i,j,k)%uf * conjg( a(i,j,k)%dudx )
Cududy = a(i,j,k)%uh_dudyh - bb(i,j,k)%uf * conjg( a(i,j,k)%dudy )
Cududz = a(i,j,k)%uh_dudzh - bb(i,j,k)%uf * conjg( a(i,j,k)%dudz )
Cudvdx = a(i,j,k)%uh_dvdxh - bb(i,j,k)%uf * conjg( a(i,j,k)%dvdx )
Cudvdy = a(i,j,k)%uh_dvdyh - bb(i,j,k)%uf * conjg( a(i,j,k)%dvdy )
Cudvdz = a(i,j,k)%uh_dvdzh - bb(i,j,k)%uf * conjg( a(i,j,k)%dvdz )
Cudwdx = a(i,j,k)%uh_dwdxh - bb(i,j,k)%uf * conjg( a(i,j,k)%dwdx )
Cudwdy = a(i,j,k)%uh_dwdyh - bb(i,j,k)%uf * conjg( a(i,j,k)%dwdy )
Cudwdz = a(i,j,k)%uh_dwdzh - bb(i,j,k)%uf * conjg( a(i,j,k)%dwdz )

Cvdudx = a(i,j,k)%vh_dudxh - bb(i,j,k)%vf * conjg( a(i,j,k)%dudx )
Cvdudy = a(i,j,k)%vh_dudyh - bb(i,j,k)%vf * conjg( a(i,j,k)%dudy )
Cvdudz = a(i,j,k)%vh_dudzh - bb(i,j,k)%vf * conjg( a(i,j,k)%dudz )
Cvdvdx = a(i,j,k)%vh_dvdxh - bb(i,j,k)%vf * conjg( a(i,j,k)%dvdx )
Cvdvdy = a(i,j,k)%vh_dvdyh - bb(i,j,k)%vf * conjg( a(i,j,k)%dvdy )
Cvdvdz = a(i,j,k)%vh_dvdzh - bb(i,j,k)%vf * conjg( a(i,j,k)%dvdz )
Cvdwdx = a(i,j,k)%vh_dwdxh - bb(i,j,k)%vf * conjg( a(i,j,k)%dwdx )
Cvdwdy = a(i,j,k)%vh_dwdyh - bb(i,j,k)%vf * conjg( a(i,j,k)%dwdy )
Cvdwdz = a(i,j,k)%vh_dwdzh - bb(i,j,k)%vf * conjg( a(i,j,k)%dwdz )

Cwdudx = a(i,j,k)%wh_dudxh - bb(i,j,k)%wf * conjg( a(i,j,k)%dudx )
Cwdudy = a(i,j,k)%wh_dudyh - bb(i,j,k)%wf * conjg( a(i,j,k)%dudy )
Cwdudz = a(i,j,k)%wh_dudzh - bb(i,j,k)%wf * conjg( a(i,j,k)%dudz )
Cwdvdx = a(i,j,k)%wh_dvdxh - bb(i,j,k)%wf * conjg( a(i,j,k)%dvdx )
Cwdvdy = a(i,j,k)%wh_dvdyh - bb(i,j,k)%wf * conjg( a(i,j,k)%dvdy )
Cwdvdz = a(i,j,k)%wh_dvdzh - bb(i,j,k)%wf * conjg( a(i,j,k)%dvdz )
Cwdwdx = a(i,j,k)%wh_dwdxh - bb(i,j,k)%wf * conjg( a(i,j,k)%dwdx )
Cwdwdy = a(i,j,k)%wh_dwdyh - bb(i,j,k)%wf * conjg( a(i,j,k)%dwdy )
Cwdwdz = a(i,j,k)%wh_dwdzh - bb(i,j,k)%wf * conjg( a(i,j,k)%dwdz )

! vel-vel-vel Grad correlation, uj_hat*(uk*duidxk)_hat
! Here the time-average mean is non-zero for i=1 (kx=0) and zero otherwise (kx~=0)
Cuududx = conjg(a(i,j,k)%uh_ududxh)           - bb(i,j,k)%uu * dudx_avg(1,j,k)             &
    - conjg(a(i,j,k)%uh_dudxh) * u_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%ududxh             &
    + 2.0_rprec*u_avg(i,j,k)*u_avg(i,j,k)*dudx_avg(i,j,k)
Cuvdudy = conjg(a(i,j,k)%uh_vdudyh)           - conjg(bb(i,j,k)%uv) * dudy_avg(1,j,k)      &
    - conjg(a(i,j,k)%uh_dudyh) * v_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%vdudyh             &
    + 2.0_rprec*u_avg(i,j,k)*v_avg(i,j,k)*dudy_avg(i,j,k)
Cuwdudz = conjg(a(i,j,k)%uh_wdudzh)           - conjg(bb(i,j,k)%uw) * dudz_avg(1,j,k)      &
    - conjg(a(i,j,k)%uh_dudzh) * w_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%wdudzh             &
    + 2.0_rprec*u_avg(i,j,k)*w_avg(i,j,k)*dudz_avg(i,j,k)
Cvududx = conjg(a(i,j,k)%vh_ududxh)           - bb(i,j,k)%uv * dudx_avg(1,j,k)             &
    - conjg(a(i,j,k)%vh_dudxh) * u_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%ududxh             &
    + 2.0_rprec*v_avg(i,j,k)*u_avg(i,j,k)*dudx_avg(i,j,k)
Cvvdudy = conjg(a(i,j,k)%vh_vdudyh)           - bb(i,j,k)%vv * dudy_avg(1,j,k)             &
    - conjg(a(i,j,k)%vh_dudyh) * v_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%vdudyh             &
    + 2.0_rprec*v_avg(i,j,k)*v_avg(i,j,k)*dudy_avg(i,j,k)
Cvwdudz = conjg(a(i,j,k)%vh_wdudzh)           - conjg(bb(i,j,k)%vw) * dudz_avg(1,j,k)      &
    - conjg(a(i,j,k)%vh_dudzh) * w_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%wdudzh             &
    + 2.0_rprec*v_avg(i,j,k)*w_avg(i,j,k)*dudz_avg(i,j,k)
Cwududx = conjg(a(i,j,k)%wh_ududxh)           - bb(i,j,k)%uw * dudx_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dudxh) * u_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%ududxh             &
    + 2.0_rprec*w_avg(i,j,k)*u_avg(i,j,k)*dudx_avg(i,j,k)
Cwvdudy = conjg(a(i,j,k)%wh_vdudyh)           - bb(i,j,k)%vw * dudy_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dudyh) * v_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%vdudyh             &
    + 2.0_rprec*w_avg(i,j,k)*v_avg(i,j,k)*dudy_avg(i,j,k)
Cwwdudz = conjg(a(i,j,k)%wh_wdudzh)           - bb(i,j,k)%ww * dudz_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dudzh) * w_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%wdudzh             &
    + 2.0_rprec*w_avg(i,j,k)*w_avg(i,j,k)*dudz_avg(i,j,k)

Cuudvdx = conjg(a(i,j,k)%uh_udvdxh)           - bb(i,j,k)%uu * dvdx_avg(1,j,k)             &
    - conjg(a(i,j,k)%uh_dvdxh) * u_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%udvdxh             &
    + 2.0_rprec*u_avg(i,j,k)*u_avg(i,j,k)*dvdx_avg(i,j,k)
Cuvdvdy = conjg(a(i,j,k)%uh_vdvdyh)           - conjg(bb(i,j,k)%uv) * dvdy_avg(1,j,k)      &
    - conjg(a(i,j,k)%uh_dvdyh) * v_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%vdvdyh             &
    + 2.0_rprec*u_avg(i,j,k)*v_avg(i,j,k)*dvdy_avg(i,j,k)
Cuwdvdz = conjg(a(i,j,k)%uh_wdvdzh)           - conjg(bb(i,j,k)%uw) * dvdz_avg(1,j,k)      &
    - conjg(a(i,j,k)%uh_dvdzh) * w_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%wdvdzh             &
    + 2.0_rprec*u_avg(i,j,k)*w_avg(i,j,k)*dvdz_avg(i,j,k)
Cvudvdx = conjg(a(i,j,k)%vh_udvdxh)           - bb(i,j,k)%uv * dvdx_avg(1,j,k)             &
    - conjg(a(i,j,k)%vh_dvdxh) * u_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%udvdxh             &
    + 2.0_rprec*v_avg(i,j,k)*u_avg(i,j,k)*dvdx_avg(i,j,k)
Cvvdvdy = conjg(a(i,j,k)%vh_vdvdyh)           - bb(i,j,k)%vv * dvdy_avg(1,j,k)             &
    - conjg(a(i,j,k)%vh_dvdyh) * v_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%vdvdyh             &
    + 2.0_rprec*v_avg(i,j,k)*v_avg(i,j,k)*dvdy_avg(i,j,k)
Cvwdvdz = conjg(a(i,j,k)%vh_wdvdzh)           - conjg(bb(i,j,k)%vw) * dvdz_avg(1,j,k)      &
    - conjg(a(i,j,k)%vh_dvdzh) * w_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%wdvdzh             &
    + 2.0_rprec*v_avg(i,j,k)*w_avg(i,j,k)*dvdz_avg(i,j,k)
Cwudvdx = conjg(a(i,j,k)%wh_udvdxh)           - bb(i,j,k)%uw * dvdx_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dvdxh) * u_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%udvdxh             &
    + 2.0_rprec*w_avg(i,j,k)*u_avg(i,j,k)*dvdx_avg(i,j,k)
Cwvdvdy = conjg(a(i,j,k)%wh_vdvdyh)           - bb(i,j,k)%vw * dvdy_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dvdyh) * v_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%vdvdyh             &
    + 2.0_rprec*w_avg(i,j,k)*v_avg(i,j,k)*dvdy_avg(i,j,k)
Cwwdvdz = conjg(a(i,j,k)%wh_wdvdzh)           - bb(i,j,k)%ww * dvdz_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dvdzh) * w_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%wdvdzh             &
    + 2.0_rprec*w_avg(i,j,k)*w_avg(i,j,k)*dvdz_avg(i,j,k)

Cuudwdx = conjg(a(i,j,k)%uh_udwdxh)           - bb(i,j,k)%uu * dwdx_avg(1,j,k)             &
    - conjg(a(i,j,k)%uh_dwdxh) * u_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%udwdxh             &
    + 2.0_rprec*u_avg(i,j,k)*u_avg(i,j,k)*dwdx_avg(i,j,k)
Cuvdwdy = conjg(a(i,j,k)%uh_vdwdyh)           - conjg(bb(i,j,k)%uv) * dwdy_avg(1,j,k)      &
    - conjg(a(i,j,k)%uh_dwdyh) * v_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%vdwdyh             &
    + 2.0_rprec*u_avg(i,j,k)*v_avg(i,j,k)*dwdy_avg(i,j,k)
Cuwdwdz = conjg(a(i,j,k)%uh_wdwdzh)           - conjg(bb(i,j,k)%uw) * dwdz_avg(1,j,k)      &
    - conjg(a(i,j,k)%uh_dwdzh) * w_avg(1,j,k) - u_avg(1,j,k) * a(i,j,k)%wdwdzh             &
    + 2.0_rprec*u_avg(i,j,k)*w_avg(i,j,k)*dwdz_avg(i,j,k)
Cvudwdx = conjg(a(i,j,k)%vh_udwdxh)           - bb(i,j,k)%uv * dwdx_avg(1,j,k)             &
    - conjg(a(i,j,k)%vh_dwdxh) * u_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%udwdxh             &
    + 2.0_rprec*v_avg(i,j,k)*u_avg(i,j,k)*dwdx_avg(i,j,k)
Cvvdwdy = conjg(a(i,j,k)%vh_vdwdyh)           - bb(i,j,k)%vv * dwdy_avg(1,j,k)             &
    - conjg(a(i,j,k)%vh_dwdyh) * v_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%vdwdyh             &
    + 2.0_rprec*v_avg(i,j,k)*v_avg(i,j,k)*dwdy_avg(i,j,k)
Cvwdwdz = conjg(a(i,j,k)%vh_wdwdzh)           - conjg(bb(i,j,k)%vw) * dwdz_avg(1,j,k)      &
    - conjg(a(i,j,k)%vh_dwdzh) * w_avg(1,j,k) - v_avg(1,j,k) * a(i,j,k)%wdwdzh             &
    + 2.0_rprec*v_avg(i,j,k)*w_avg(i,j,k)*dwdz_avg(i,j,k)
Cwudwdx = conjg(a(i,j,k)%wh_udwdxh)           - bb(i,j,k)%uw * dwdx_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dwdxh) * u_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%udwdxh             &
    + 2.0_rprec*w_avg(i,j,k)*u_avg(i,j,k)*dwdx_avg(i,j,k)
Cwvdwdy = conjg(a(i,j,k)%wh_vdwdyh)           - bb(i,j,k)%vw * dwdy_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dwdyh) * v_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%vdwdyh             &
    + 2.0_rprec*w_avg(i,j,k)*v_avg(i,j,k)*dwdy_avg(i,j,k)
Cwwdwdz = conjg(a(i,j,k)%wh_wdwdzh)           - bb(i,j,k)%ww * dwdz_avg(1,j,k)             &
    - conjg(a(i,j,k)%wh_dwdzh) * w_avg(1,j,k) - w_avg(1,j,k) * a(i,j,k)%wdwdzh             &
    + 2.0_rprec*w_avg(i,j,k)*w_avg(i,j,k)*dwdz_avg(i,j,k)

! Pressure-Velocity Gradient correlation, p*dujdxk
Cpdudx = a(i,j,k)%pdudx - a(i,j,k)%p * conjg( a(i,j,k)%dudx )
Cpdudy = a(i,j,k)%pdudy - a(i,j,k)%p * conjg( a(i,j,k)%dudy )
Cpdudz = a(i,j,k)%pdudz - a(i,j,k)%p * conjg( a(i,j,k)%dudz )
Cpdvdx = a(i,j,k)%pdvdx - a(i,j,k)%p * conjg( a(i,j,k)%dvdx )
Cpdvdy = a(i,j,k)%pdvdy - a(i,j,k)%p * conjg( a(i,j,k)%dvdy )
Cpdvdz = a(i,j,k)%pdvdz - a(i,j,k)%p * conjg( a(i,j,k)%dvdz )
Cpdwdx = a(i,j,k)%pdwdx - a(i,j,k)%p * conjg( a(i,j,k)%dwdx )
Cpdwdy = a(i,j,k)%pdwdy - a(i,j,k)%p * conjg( a(i,j,k)%dwdy )
Cpdwdz = a(i,j,k)%pdwdz - a(i,j,k)%p * conjg( a(i,j,k)%dwdz )

! Velocity-Pressure Gradient correlation, uj*dpdxk
Cudpdx = a(i,j,k)%udpdx - bb(i,j,k)%uf * conjg( a(i,j,k)%dpdx )
Cudpdy = a(i,j,k)%udpdy - bb(i,j,k)%uf * conjg( a(i,j,k)%dpdy )
Cudpdz = a(i,j,k)%udpdz - bb(i,j,k)%uf * conjg( a(i,j,k)%dpdz )
Cvdpdx = a(i,j,k)%vdpdx - bb(i,j,k)%vf * conjg( a(i,j,k)%dpdx )
Cvdpdy = a(i,j,k)%vdpdy - bb(i,j,k)%vf * conjg( a(i,j,k)%dpdy )
Cvdpdz = a(i,j,k)%vdpdz - bb(i,j,k)%vf * conjg( a(i,j,k)%dpdz )
Cwdpdx = a(i,j,k)%wdpdx - bb(i,j,k)%wf * conjg( a(i,j,k)%dpdx )
Cwdpdy = a(i,j,k)%wdpdy - bb(i,j,k)%wf * conjg( a(i,j,k)%dpdy )
Cwdpdz = a(i,j,k)%wdpdz - bb(i,j,k)%wf * conjg( a(i,j,k)%dpdz )

! velGrad-velGrad correlation, duidxk*dujdxk, i=j
Cuxux = a(i,j,k)%uxux - a(i,j,k)%dudx * conjg(a(i,j,k)%dudx)
Cuyuy = a(i,j,k)%uyuy - a(i,j,k)%dudy * conjg(a(i,j,k)%dudy)
Cuzuz = a(i,j,k)%uzuz - a(i,j,k)%dudz * conjg(a(i,j,k)%dudz)
Cvxvx = a(i,j,k)%vxvx - a(i,j,k)%dvdx * conjg(a(i,j,k)%dvdx)
Cvyvy = a(i,j,k)%vyvy - a(i,j,k)%dvdy * conjg(a(i,j,k)%dvdy)
Cvzvz = a(i,j,k)%vzvz - a(i,j,k)%dvdz * conjg(a(i,j,k)%dvdz)
Cwxwx = a(i,j,k)%wxwx - a(i,j,k)%dwdx * conjg(a(i,j,k)%dwdx)
Cwywy = a(i,j,k)%wywy - a(i,j,k)%dwdy * conjg(a(i,j,k)%dwdy)
Cwzwz = a(i,j,k)%wzwz - a(i,j,k)%dwdz * conjg(a(i,j,k)%dwdz)

! velGrad-velGrad correlation, duidxk*dujdxk, i\=j
Cuxvx = a(i,j,k)%uxvx - a(i,j,k)%dudx * conjg(a(i,j,k)%dvdx)
Cuyvy = a(i,j,k)%uyvy - a(i,j,k)%dudy * conjg(a(i,j,k)%dvdy)
Cuzvz = a(i,j,k)%uzvz - a(i,j,k)%dudz * conjg(a(i,j,k)%dvdz)
Cuxwx = a(i,j,k)%uxwx - a(i,j,k)%dudx * conjg(a(i,j,k)%dwdx)
Cuywy = a(i,j,k)%uywy - a(i,j,k)%dudy * conjg(a(i,j,k)%dwdy)
Cuzwz = a(i,j,k)%uzwz - a(i,j,k)%dudz * conjg(a(i,j,k)%dwdz)
Cvxwx = a(i,j,k)%vxwx - a(i,j,k)%dvdx * conjg(a(i,j,k)%dwdx)
Cvywy = a(i,j,k)%vywy - a(i,j,k)%dvdy * conjg(a(i,j,k)%dwdy)
Cvzwz = a(i,j,k)%vzwz - a(i,j,k)%dvdz * conjg(a(i,j,k)%dwdz)

! vel-Laplacian correlation, ui*lap(uj)
Culapu = a(i,j,k)%ulapu - bb(i,j,k)%uf * conjg(a(i,j,k)%lapu)
Culapv = a(i,j,k)%ulapv - bb(i,j,k)%uf * conjg(a(i,j,k)%lapv)
Culapw = a(i,j,k)%ulapw - bb(i,j,k)%uf * conjg(a(i,j,k)%lapw)
Cvlapu = a(i,j,k)%vlapu - bb(i,j,k)%vf * conjg(a(i,j,k)%lapu)
Cvlapv = a(i,j,k)%vlapv - bb(i,j,k)%vf * conjg(a(i,j,k)%lapv)
Cvlapw = a(i,j,k)%vlapw - bb(i,j,k)%vf * conjg(a(i,j,k)%lapw)
Cwlapu = a(i,j,k)%wlapu - bb(i,j,k)%wf * conjg(a(i,j,k)%lapu)
Cwlapv = a(i,j,k)%wlapv - bb(i,j,k)%wf * conjg(a(i,j,k)%lapv)
Cwlapw = a(i,j,k)%wlapw - bb(i,j,k)%wf * conjg(a(i,j,k)%lapw)

! ------------------------- Compute Budget Terms -----------------------------
! Mean Advection, Uk*dRijdxk
! Multiplying by kx=0 mode time-average here
e(i,j,k) % advxx = u_avg(1,j,k)*real(conjg(Cududx) + Cududx)       &
    + v_avg(1,j,k)*real(conjg(Cududy) + Cududy)                    &
    + w_avg(1,j,k)*real(conjg(Cududz) + Cududz)
e(i,j,k) % advyy = u_avg(1,j,k)*real(conjg(Cvdvdx) + Cvdvdx)       &
    + v_avg(1,j,k)*real(conjg(Cvdvdy) + Cvdvdy)                    &
    + w_avg(1,j,k)*real(conjg(Cvdvdz) + Cvdvdz)
e(i,j,k) % advzz = u_avg(1,j,k)*real(conjg(Cwdwdx) + Cwdwdx)       &
    + v_avg(1,j,k)*real(conjg(Cwdwdy) + Cwdwdy)                    &
    + w_avg(1,j,k)*real(conjg(Cwdwdz) + Cwdwdz)
e(i,j,k) % advxy = u_avg(1,j,k)*real(conjg(Cvdudx) + Cudvdx)       &
    + v_avg(1,j,k)*real(conjg(Cvdudy) + Cudvdy)                    &
    + w_avg(1,j,k)*real(conjg(Cvdudz) + Cudvdz)
e(i,j,k) % advxz = u_avg(1,j,k)*real(conjg(Cwdudx) + Cudwdx)       &
    + v_avg(1,j,k)*real(conjg(Cwdudy) + Cudwdy)                    & 
    + w_avg(1,j,k)*real(conjg(Cwdudz) + Cudwdz)
e(i,j,k) % advyz = u_avg(1,j,k)*real(conjg(Cwdvdx) + Cvdwdx)       &
    + v_avg(1,j,k)*real(conjg(Cwdvdy) + Cvdwdy)                    &
    + w_avg(1,j,k)*real(conjg(Cwdvdz) + Cvdwdz)

! Transport by fluctuations, d(ui*uj*uk)dxk = ui*uk*dujdxk + uj*uk*duidxk
e(i,j,k) % tflucxx = real( Cuududx + conjg(Cuududx) +               &
    Cuvdudy + conjg(Cuvdudy) + Cuwdudz + conjg(Cuwdudz) )
e(i,j,k) % tflucyy = real( Cvudvdx + conjg(Cvudvdx) +               &
    Cvvdvdy + conjg(Cvvdvdy) + Cvwdvdz + conjg(Cvwdvdz) )
e(i,j,k) % tfluczz = real( Cwudwdx + conjg(Cwudwdx) +               &
    Cwvdwdy + conjg(Cwvdwdy) + Cwwdwdz + conjg(Cwwdwdz) )
e(i,j,k) % tflucxy = real( Cvududx + conjg(Cuudvdx) +               &
    Cvvdudy + conjg(Cuvdvdy) + Cvwdudz + conjg(Cuwdvdz) )
e(i,j,k) % tflucxz = real( Cwududx + conjg(Cuudwdx) +               &
    Cwvdudy + conjg(Cuvdwdy) + Cwwdudz + conjg(Cuwdwdz) )
e(i,j,k) % tflucyz = real( Cwudvdx + conjg(Cvudwdx) +               &
    Cwvdvdy + conjg(Cvvdwdy) + Cwwdvdz + conjg(Cvwdwdz) )

! Pressure-Strain, 2*p*sij = p*(dujdxi + duidxj)
! x-derivative terms have no transport, only pressure-strain
! yz-derivative terms have both transport and strain components
! after having applied the product rule
e(i,j,k) % pstrainxx = real( conjg(Cudpdx) + Cudpdx )
e(i,j,k) % pstrainyy = -real( Cpdvdy + conjg(Cpdvdy) )
e(i,j,k) % pstrainzz = -real( Cpdwdz + conjg(Cpdwdz) )
e(i,j,k) % pstrainxy = real( conjg(Cvdpdx) - conjg(Cpdudy) )
e(i,j,k) % pstrainxz = real( conjg(Cwdpdx) - conjg(Cpdudz) )
e(i,j,k) % pstrainyz = -real( Cpdwdy + conjg(Cpdvdz) )

! Transport by pressure-velocity, d(p*(Iik+ui*Ijk))dxk = d(p*uj)dxi + d(p*ui)dxj
! Terms that involve ddx have been transformed and contribute as a source 
! term not a transport term and are therefore included in the pressure-strain
e(i,j,k) % tpresxx = 0.0_rprec 
e(i,j,k) % tpresyy = real( Cpdvdy + conjg(Cpdvdy) + conjg(Cvdpdy) + Cvdpdy )
e(i,j,k) % tpreszz = real( Cpdwdz + conjg(Cpdwdz) + conjg(Cwdpdz) + Cwdpdz )
e(i,j,k) % tpresxy = real( conjg(Cpdudy) + Cudpdy )
e(i,j,k) % tpresxz = real( conjg(Cpdudz) + Cudpdz )
e(i,j,k) % tpresyz = real( Cpdwdy + conjg(Cpdvdz) + conjg(Cwdpdy) + Cvdpdz )

! Production Rate, -(Rik*dUjdxk + Rjk*dUidxk)
! Multiplying by kx=0 mode time-average here
e(i,j,k) % prodxx = 2.0_rprec*( b(i,j,k)%upup * dudx_avg(1,j,k) +                        &
    b(i,j,k)%upvp * dudy_avg(1,j,k) + b(i,j,k)%upwp * dudz_avg(1,j,k) )
e(i,j,k) % prodyy = 2.0_rprec*( b(i,j,k)%upvp * dvdx_avg(1,j,k) +                        &
    b(i,j,k)%vpvp * dvdy_avg(1,j,k) + b(i,j,k)%vpwp * dvdz_avg(1,j,k) )
e(i,j,k) % prodzz = 2.0_rprec*( b(i,j,k)%upwp * dwdx_avg(1,j,k) +                        &
    b(i,j,k)%vpwp * dwdy_avg(1,j,k) + b(i,j,k)%wpwp * dwdz_avg(1,j,k) )
e(i,j,k) % prodxy = (b(i,j,k)%upup * dvdx_avg(1,j,k) + b(i,j,k)%upvp*dudx_avg(1,j,k) +   &
    b(i,j,k)%upvp * dvdy_avg(1,j,k) + b(i,j,k)%vpvp * dudy_avg(1,j,k) +                  &
    b(i,j,k)%upwp * dvdz_avg(1,j,k) + b(i,j,k)%vpwp * dudz_avg(1,j,k) )
e(i,j,k) % prodxz = (b(i,j,k)%upup * dwdx_avg(1,j,k) + b(i,j,k)%upwp*dudx_avg(1,j,k) +   &
    b(i,j,k)%upvp * dwdy_avg(1,j,k) + b(i,j,k)%vpwp * dudy_avg(1,j,k) +                  &
    b(i,j,k)%upwp * dwdz_avg(1,j,k) + b(i,j,k)%wpwp * dudz_avg(1,j,k) )
e(i,j,k) % prodyz = (b(i,j,k)%upvp * dwdx_avg(1,j,k) + b(i,j,k)%upwp*dvdx_avg(1,j,k) +   &
    b(i,j,k)%vpvp * dwdy_avg(1,j,k) + b(i,j,k)%vpwp * dvdy_avg(1,j,k) +                  &
    b(i,j,k)%vpwp * dwdz_avg(1,j,k) + b(i,j,k)%wpwp * dvdz_avg(1,j,k) )

! Pseudo-dissipation, 2*nu*duidxk*dujdxk
e(i,j,k) % pdissxx = 2.0_rprec*nu_molec*real(Cuxux + Cuyuy + Cuzuz)
e(i,j,k) % pdissyy = 2.0_rprec*nu_molec*real(Cvxvx + Cvyvy + Cvzvz)
e(i,j,k) % pdisszz = 2.0_rprec*nu_molec*real(Cwxwx + Cwywy + Cwzwz)
e(i,j,k) % pdissxy = 2.0_rprec*nu_molec*real(Cuxvx + Cuyvy + Cuzvz)
e(i,j,k) % pdissxz = 2.0_rprec*nu_molec*real(Cuxwx + Cuywy + Cuzwz)
e(i,j,k) % pdissyz = 2.0_rprec*nu_molec*real(Cvxwx + Cvywy + Cvzwz)

! Transport by viscous diffusion
! Cuilapuj already multiplied by viscosity because divtxj was used
! Change sign of Cuilapuj because of how divtxj is set up
! Note: lapui terms include x-deriveration terms from pseudo-dissipation
!       which need to be canceled, which is why Cuxux terms have different sign
e(i,j,k) % tviscxx = 2.0_rprec*nu_molec*real(-Cuxux + Cuyuy + Cuzuz) - real( conjg(Culapu) + Culapu )
e(i,j,k) % tviscyy = 2.0_rprec*nu_molec*real(-Cvxvx + Cvyvy + Cvzvz) - real( conjg(Cvlapv) + Cvlapv )
e(i,j,k) % tvisczz = 2.0_rprec*nu_molec*real(-Cwxwx + Cwywy + Cwzwz) - real( conjg(Cwlapw) + Cwlapw )
e(i,j,k) % tviscxy = 2.0_rprec*nu_molec*real(-Cuxvx + Cuyvy + Cuzvz) - real( conjg(Cvlapu) + Culapv )
e(i,j,k) % tviscxz = 2.0_rprec*nu_molec*real(-Cuxwx + Cuywy + Cuzwz) - real( conjg(Cwlapu) + Culapw )
e(i,j,k) % tviscyz = 2.0_rprec*nu_molec*real(-Cvxwx + Cvywy + Cvzwz) - real( conjg(Cwlapv) + Cvlapw )

end do
end do
end do

! ---------------------------- Prepare Output --------------------------------
! Move all terms to the RHS
e % advxx = - e % advxx
e % advyy = - e % advyy
e % advzz = - e % advzz
e % advxy = - e % advxy
e % advxz = - e % advxz
e % advyz = - e % advyz

e % tflucxx = - e % tflucxx
e % tflucyy = - e % tflucyy
e % tfluczz = - e % tfluczz
e % tflucxy = - e % tflucxy
e % tflucxz = - e % tflucxz
e % tflucyz = - e % tflucyz

e % prodxx = - e % prodxx
e % prodyy = - e % prodyy
e % prodzz = - e % prodzz
e % prodxy = - e % prodxy
e % prodxz = - e % prodxz
e % prodyz = - e % prodyz

e % pstrainxx = - e % pstrainxx
e % pstrainyy = - e % pstrainyy
e % pstrainzz = - e % pstrainzz
e % pstrainxy = - e % pstrainxy
e % pstrainxz = - e % pstrainxz
e % pstrainyz = - e % pstrainyz

e % tpresxx = - e % tpresxx
e % tpresyy = - e % tpresyy
e % tpreszz = - e % tpreszz
e % tpresxy = - e % tpresxy
e % tpresxz = - e % tpresxz
e % tpresyz = - e % tpresyz

e % pdissxx = - e % pdissxx
e % pdissyy = - e % pdissyy
e % pdisszz = - e % pdisszz
e % pdissxy = - e % pdissxy
e % pdissxz = - e % pdissxz
e % pdissyz = - e % pdissyz

end function specbudgx_compute

!*****************************************************************************
function specbudgy_compute( a, b, bb, c, d, lbz2 ) result( e )
!*****************************************************************************
! 
! Compute the spectral budget. All quantities should be on the w-grid.
! 
use param, only: nu_molec
implicit none
integer, intent(in) :: lbz2
type(tavg_specbudg_t), dimension(:,:,lbz2:), intent(in) :: a
type(turbspec_t), dimension(:,:,lbz2:), intent(in) :: b
type(tavg_turbspec_t), dimension(:,:,lbz2:), intent(in) :: bb
type(tavg_budget_t), dimension(:,:,lbz2:), intent(in) :: c
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: d
type(specbudg_t), allocatable, dimension(:,:,:) :: e

integer :: ubx, uby, ubz, vbx, vby, vbz, i, j, k

real(rprec), allocatable, dimension(:,:,:) :: u_avg, v_avg, w_avg,          &
    dudx_avg, dvdx_avg, dwdx_avg, dudy_avg, dvdy_avg, dwdy_avg,             &
    dudz_avg, dvdz_avg, dwdz_avg

complex(rprec) :: Cududx, Cududy, Cududz, Cudvdx, Cudvdy, Cudvdz, Cudwdx, Cudwdy, Cudwdz,  &
    Cvdudx, Cvdudy, Cvdudz, Cvdvdx, Cvdvdy, Cvdvdz, Cvdwdx, Cvdwdy, Cvdwdz,                &
    Cwdudx, Cwdudy, Cwdudz, Cwdvdx, Cwdvdy, Cwdvdz, Cwdwdx, Cwdwdy, Cwdwdz 

complex(rprec) :: Cuududx, Cuvdudy, Cuwdudz, Cvududx, Cvvdudy, Cvwdudz, Cwududx, Cwvdudy, Cwwdudz, &
    Cuudvdx, Cuvdvdy, Cuwdvdz, Cvudvdx, Cvvdvdy, Cvwdvdz, Cwudvdx, Cwvdvdy, Cwwdvdz,               &
    Cuudwdx, Cuvdwdy, Cuwdwdz, Cvudwdx, Cvvdwdy, Cvwdwdz, Cwudwdx, Cwvdwdy, Cwwdwdz

complex(rprec) :: Cpdudx, Cpdudy, Cpdudz, Cpdvdx, Cpdvdy, Cpdvdz, &
    Cpdwdx, Cpdwdy, Cpdwdz

complex(rprec) :: Cudpdx, Cudpdy, Cudpdz, Cvdpdx, Cvdpdy, Cvdpdz, &
    Cwdpdx, Cwdpdy, Cwdpdz

complex(rprec) :: Cuxux, Cuyuy, Cuzuz, Cvxvx, Cvyvy, Cvzvz, Cwxwx, Cwywy, Cwzwz
complex(rprec) :: Cuxvx, Cuyvy, Cuzvz, Cuxwx, Cuywy, Cuzwz, Cvxwx, Cvywy, Cvzwz
complex(rprec) :: Culapu, Culapv, Culapw, Cvlapu, Cvlapv, Cvlapw, Cwlapu, Cwlapv, Cwlapw

! Size of a, b, and e
ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

! Size of c and d
vbx = ubound(c,1)
vby = ubound(c,2)
vbz = ubound(c,3)

allocate(e(ubx,uby,lbz2:ubz))

allocate(u_avg(vbx,uby,lbz2:vbz))
allocate(v_avg(vbx,uby,lbz2:vbz))
allocate(w_avg(vbx,uby,lbz2:vbz))
allocate(dudx_avg(vbx,uby,lbz2:vbz))
allocate(dudy_avg(vbx,uby,lbz2:vbz))
allocate(dudz_avg(vbx,uby,lbz2:vbz))
allocate(dvdx_avg(vbx,uby,lbz2:vbz))
allocate(dvdy_avg(vbx,uby,lbz2:vbz))
allocate(dvdz_avg(vbx,uby,lbz2:vbz))
allocate(dwdx_avg(vbx,uby,lbz2:vbz))
allocate(dwdy_avg(vbx,uby,lbz2:vbz))
allocate(dwdz_avg(vbx,uby,lbz2:vbz))

! -------------------------- Spanwise Average ------------------------------
! When deriving the spectral budget equation, there are transforms of the 
! time-averaged mean multiplied with other quantities. To avoid convolutions
! in the budget equation, it is assumed that the time-averaged mean is not a 
! function of the direction in which the transform is being applied.
! 
! Because of this, we average the time-averaged quantities in that direction
! before removing the time-averaged mean component from the instantaneous 
! running sums.
!
! Quantities used in average are in x,y,z space, however only ky=0 mode
! is filled, ky~=0 modes are left as zeros
u_avg(:,:,:) = 0._rprec
v_avg(:,:,:) = 0._rprec
w_avg(:,:,:) = 0._rprec
dudx_avg(:,:,:) = 0._rprec
dudy_avg(:,:,:) = 0._rprec
dudz_avg(:,:,:) = 0._rprec
dvdx_avg(:,:,:) = 0._rprec
dvdy_avg(:,:,:) = 0._rprec
dvdz_avg(:,:,:) = 0._rprec
dwdx_avg(:,:,:) = 0._rprec
dwdy_avg(:,:,:) = 0._rprec
dwdz_avg(:,:,:) = 0._rprec
do i = 1, vby
    u_avg(:,1,:) = u_avg(:,1,:) + d(:,i,:) % u_w
    v_avg(:,1,:) = v_avg(:,1,:) + d(:,i,:) % v_w
    w_avg(:,1,:) = w_avg(:,1,:) + d(:,i,:) % w
    dudx_avg(:,1,:) = dudx_avg(:,1,:) + c(:,i,:) % dudx
    dudy_avg(:,1,:) = dudy_avg(:,1,:) + c(:,i,:) % dudy
    dudz_avg(:,1,:) = dudz_avg(:,1,:) + c(:,i,:) % dudz
    dvdx_avg(:,1,:) = dvdx_avg(:,1,:) + c(:,i,:) % dvdx
    dvdy_avg(:,1,:) = dvdy_avg(:,1,:) + c(:,i,:) % dvdy
    dvdz_avg(:,1,:) = dvdz_avg(:,1,:) + c(:,i,:) % dvdz
    dwdx_avg(:,1,:) = dwdx_avg(:,1,:) + c(:,i,:) % dwdx
    dwdy_avg(:,1,:) = dwdy_avg(:,1,:) + c(:,i,:) % dwdy
    dwdz_avg(:,1,:) = dwdz_avg(:,1,:) + c(:,i,:) % dwdz
enddo
u_avg(:,1,:) = u_avg(:,1,:) / vby
v_avg(:,1,:) = v_avg(:,1,:) / vby
w_avg(:,1,:) = w_avg(:,1,:) / vby
dudx_avg(:,1,:) = dudx_avg(:,1,:) / vby
dudy_avg(:,1,:) = dudy_avg(:,1,:) / vby
dudz_avg(:,1,:) = dudz_avg(:,1,:) / vby
dvdx_avg(:,1,:) = dvdx_avg(:,1,:) / vby
dvdy_avg(:,1,:) = dvdy_avg(:,1,:) / vby
dvdz_avg(:,1,:) = dvdz_avg(:,1,:) / vby
dwdx_avg(:,1,:) = dwdx_avg(:,1,:) / vby
dwdy_avg(:,1,:) = dwdy_avg(:,1,:) / vby
dwdz_avg(:,1,:) = dwdz_avg(:,1,:) / vby

do i = 1, ubx
do j = 1, uby
do k = 1, ubz

! ---------------------- Compute Intermediate Terms --------------------------
! Compute intermediate terms used in the budget
! Velocity-Velocity Gradient correlation, ui_hat*dujdxk_hat
Cududx = a(i,j,k)%uh_dudxh - bb(i,j,k)%uf * conjg( a(i,j,k)%dudx )
Cududy = a(i,j,k)%uh_dudyh - bb(i,j,k)%uf * conjg( a(i,j,k)%dudy )
Cududz = a(i,j,k)%uh_dudzh - bb(i,j,k)%uf * conjg( a(i,j,k)%dudz )
Cudvdx = a(i,j,k)%uh_dvdxh - bb(i,j,k)%uf * conjg( a(i,j,k)%dvdx )
Cudvdy = a(i,j,k)%uh_dvdyh - bb(i,j,k)%uf * conjg( a(i,j,k)%dvdy )
Cudvdz = a(i,j,k)%uh_dvdzh - bb(i,j,k)%uf * conjg( a(i,j,k)%dvdz )
Cudwdx = a(i,j,k)%uh_dwdxh - bb(i,j,k)%uf * conjg( a(i,j,k)%dwdx )
Cudwdy = a(i,j,k)%uh_dwdyh - bb(i,j,k)%uf * conjg( a(i,j,k)%dwdy )
Cudwdz = a(i,j,k)%uh_dwdzh - bb(i,j,k)%uf * conjg( a(i,j,k)%dwdz )

Cvdudx = a(i,j,k)%vh_dudxh - bb(i,j,k)%vf * conjg( a(i,j,k)%dudx )
Cvdudy = a(i,j,k)%vh_dudyh - bb(i,j,k)%vf * conjg( a(i,j,k)%dudy )
Cvdudz = a(i,j,k)%vh_dudzh - bb(i,j,k)%vf * conjg( a(i,j,k)%dudz )
Cvdvdx = a(i,j,k)%vh_dvdxh - bb(i,j,k)%vf * conjg( a(i,j,k)%dvdx )
Cvdvdy = a(i,j,k)%vh_dvdyh - bb(i,j,k)%vf * conjg( a(i,j,k)%dvdy )
Cvdvdz = a(i,j,k)%vh_dvdzh - bb(i,j,k)%vf * conjg( a(i,j,k)%dvdz )
Cvdwdx = a(i,j,k)%vh_dwdxh - bb(i,j,k)%vf * conjg( a(i,j,k)%dwdx )
Cvdwdy = a(i,j,k)%vh_dwdyh - bb(i,j,k)%vf * conjg( a(i,j,k)%dwdy )
Cvdwdz = a(i,j,k)%vh_dwdzh - bb(i,j,k)%vf * conjg( a(i,j,k)%dwdz )

Cwdudx = a(i,j,k)%wh_dudxh - bb(i,j,k)%wf * conjg( a(i,j,k)%dudx )
Cwdudy = a(i,j,k)%wh_dudyh - bb(i,j,k)%wf * conjg( a(i,j,k)%dudy )
Cwdudz = a(i,j,k)%wh_dudzh - bb(i,j,k)%wf * conjg( a(i,j,k)%dudz )
Cwdvdx = a(i,j,k)%wh_dvdxh - bb(i,j,k)%wf * conjg( a(i,j,k)%dvdx )
Cwdvdy = a(i,j,k)%wh_dvdyh - bb(i,j,k)%wf * conjg( a(i,j,k)%dvdy )
Cwdvdz = a(i,j,k)%wh_dvdzh - bb(i,j,k)%wf * conjg( a(i,j,k)%dvdz )
Cwdwdx = a(i,j,k)%wh_dwdxh - bb(i,j,k)%wf * conjg( a(i,j,k)%dwdx )
Cwdwdy = a(i,j,k)%wh_dwdyh - bb(i,j,k)%wf * conjg( a(i,j,k)%dwdy )
Cwdwdz = a(i,j,k)%wh_dwdzh - bb(i,j,k)%wf * conjg( a(i,j,k)%dwdz )

! vel-vel-vel Grad correlation, uj_hat*(uk*duidxk)_hat
! Here the time-average mean is non-zero for i=1 (kx=0) and zero otherwise (kx~=0)
Cuududx = conjg(a(i,j,k)%uh_ududxh)           - bb(i,j,k)%uu * dudx_avg(i,1,k)             &
    - conjg(a(i,j,k)%uh_dudxh) * u_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%ududxh             &
    + 2.0_rprec*u_avg(i,j,k)*u_avg(i,j,k)*dudx_avg(i,j,k)
Cuvdudy = conjg(a(i,j,k)%uh_vdudyh)           - conjg(bb(i,j,k)%uv) * dudy_avg(i,1,k)      &
    - conjg(a(i,j,k)%uh_dudyh) * v_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%vdudyh             &
    + 2.0_rprec*u_avg(i,j,k)*v_avg(i,j,k)*dudy_avg(i,j,k)
Cuwdudz = conjg(a(i,j,k)%uh_wdudzh)           - conjg(bb(i,j,k)%uw) * dudz_avg(i,1,k)      &
    - conjg(a(i,j,k)%uh_dudzh) * w_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%wdudzh             &
    + 2.0_rprec*u_avg(i,j,k)*w_avg(i,j,k)*dudz_avg(i,j,k)
Cvududx = conjg(a(i,j,k)%vh_ududxh)           - bb(i,j,k)%uv * dudx_avg(i,1,k)             &
    - conjg(a(i,j,k)%vh_dudxh) * u_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%ududxh             &
    + 2.0_rprec*v_avg(i,j,k)*u_avg(i,j,k)*dudx_avg(i,j,k)
Cvvdudy = conjg(a(i,j,k)%vh_vdudyh)           - bb(i,j,k)%vv * dudy_avg(i,1,k)             &
    - conjg(a(i,j,k)%vh_dudyh) * v_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%vdudyh             &
    + 2.0_rprec*v_avg(i,j,k)*v_avg(i,j,k)*dudy_avg(i,j,k)
Cvwdudz = conjg(a(i,j,k)%vh_wdudzh)           - conjg(bb(i,j,k)%vw) * dudz_avg(i,1,k)      &
    - conjg(a(i,j,k)%vh_dudzh) * w_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%wdudzh             &
    + 2.0_rprec*v_avg(i,j,k)*w_avg(i,j,k)*dudz_avg(i,j,k)
Cwududx = conjg(a(i,j,k)%wh_ududxh)           - bb(i,j,k)%uw * dudx_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dudxh) * u_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%ududxh             &
    + 2.0_rprec*w_avg(i,j,k)*u_avg(i,j,k)*dudx_avg(i,j,k)
Cwvdudy = conjg(a(i,j,k)%wh_vdudyh)           - bb(i,j,k)%vw * dudy_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dudyh) * v_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%vdudyh             &
    + 2.0_rprec*w_avg(i,j,k)*v_avg(i,j,k)*dudy_avg(i,j,k)
Cwwdudz = conjg(a(i,j,k)%wh_wdudzh)           - bb(i,j,k)%ww * dudz_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dudzh) * w_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%wdudzh             &
    + 2.0_rprec*w_avg(i,j,k)*w_avg(i,j,k)*dudz_avg(i,j,k)

Cuudvdx = conjg(a(i,j,k)%uh_udvdxh)           - bb(i,j,k)%uu * dvdx_avg(i,1,k)             &
    - conjg(a(i,j,k)%uh_dvdxh) * u_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%udvdxh             &
    + 2.0_rprec*u_avg(i,j,k)*u_avg(i,j,k)*dvdx_avg(i,j,k)
Cuvdvdy = conjg(a(i,j,k)%uh_vdvdyh)           - conjg(bb(i,j,k)%uv) * dvdy_avg(i,1,k)      &
    - conjg(a(i,j,k)%uh_dvdyh) * v_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%vdvdyh             &
    + 2.0_rprec*u_avg(i,j,k)*v_avg(i,j,k)*dvdy_avg(i,j,k)
Cuwdvdz = conjg(a(i,j,k)%uh_wdvdzh)           - conjg(bb(i,j,k)%uw) * dvdz_avg(i,1,k)      &
    - conjg(a(i,j,k)%uh_dvdzh) * w_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%wdvdzh             &
    + 2.0_rprec*u_avg(i,j,k)*w_avg(i,j,k)*dvdz_avg(i,j,k)
Cvudvdx = conjg(a(i,j,k)%vh_udvdxh)           - bb(i,j,k)%uv * dvdx_avg(i,1,k)             &
    - conjg(a(i,j,k)%vh_dvdxh) * u_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%udvdxh             &
    + 2.0_rprec*v_avg(i,j,k)*u_avg(i,j,k)*dvdx_avg(i,j,k)
Cvvdvdy = conjg(a(i,j,k)%vh_vdvdyh)           - bb(i,j,k)%vv * dvdy_avg(i,1,k)             &
    - conjg(a(i,j,k)%vh_dvdyh) * v_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%vdvdyh             &
    + 2.0_rprec*v_avg(i,j,k)*v_avg(i,j,k)*dvdy_avg(i,j,k)
Cvwdvdz = conjg(a(i,j,k)%vh_wdvdzh)           - conjg(bb(i,j,k)%vw) * dvdz_avg(i,1,k)      &
    - conjg(a(i,j,k)%vh_dvdzh) * w_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%wdvdzh             &
    + 2.0_rprec*v_avg(i,j,k)*w_avg(i,j,k)*dvdz_avg(i,j,k)
Cwudvdx = conjg(a(i,j,k)%wh_udvdxh)           - bb(i,j,k)%uw * dvdx_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dvdxh) * u_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%udvdxh             &
    + 2.0_rprec*w_avg(i,j,k)*u_avg(i,j,k)*dvdx_avg(i,j,k)
Cwvdvdy = conjg(a(i,j,k)%wh_vdvdyh)           - bb(i,j,k)%vw * dvdy_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dvdyh) * v_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%vdvdyh             &
    + 2.0_rprec*w_avg(i,j,k)*v_avg(i,j,k)*dvdy_avg(i,j,k)
Cwwdvdz = conjg(a(i,j,k)%wh_wdvdzh)           - bb(i,j,k)%ww * dvdz_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dvdzh) * w_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%wdvdzh             &
    + 2.0_rprec*w_avg(i,j,k)*w_avg(i,j,k)*dvdz_avg(i,j,k)

Cuudwdx = conjg(a(i,j,k)%uh_udwdxh)           - bb(i,j,k)%uu * dwdx_avg(i,1,k)             &
    - conjg(a(i,j,k)%uh_dwdxh) * u_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%udwdxh             &
    + 2.0_rprec*u_avg(i,j,k)*u_avg(i,j,k)*dwdx_avg(i,j,k)
Cuvdwdy = conjg(a(i,j,k)%uh_vdwdyh)           - conjg(bb(i,j,k)%uv) * dwdy_avg(i,1,k)      &
    - conjg(a(i,j,k)%uh_dwdyh) * v_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%vdwdyh             &
    + 2.0_rprec*u_avg(i,j,k)*v_avg(i,j,k)*dwdy_avg(i,j,k)
Cuwdwdz = conjg(a(i,j,k)%uh_wdwdzh)           - conjg(bb(i,j,k)%uw) * dwdz_avg(i,1,k)      &
    - conjg(a(i,j,k)%uh_dwdzh) * w_avg(i,1,k) - u_avg(i,1,k) * a(i,j,k)%wdwdzh             &
    + 2.0_rprec*u_avg(i,j,k)*w_avg(i,j,k)*dwdz_avg(i,j,k)
Cvudwdx = conjg(a(i,j,k)%vh_udwdxh)           - bb(i,j,k)%uv * dwdx_avg(i,1,k)             &
    - conjg(a(i,j,k)%vh_dwdxh) * u_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%udwdxh             &
    + 2.0_rprec*v_avg(i,j,k)*u_avg(i,j,k)*dwdx_avg(i,j,k)
Cvvdwdy = conjg(a(i,j,k)%vh_vdwdyh)           - bb(i,j,k)%vv * dwdy_avg(i,1,k)             &
    - conjg(a(i,j,k)%vh_dwdyh) * v_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%vdwdyh             &
    + 2.0_rprec*v_avg(i,j,k)*v_avg(i,j,k)*dwdy_avg(i,j,k)
Cvwdwdz = conjg(a(i,j,k)%vh_wdwdzh)           - conjg(bb(i,j,k)%vw) * dwdz_avg(i,1,k)      &
    - conjg(a(i,j,k)%vh_dwdzh) * w_avg(i,1,k) - v_avg(i,1,k) * a(i,j,k)%wdwdzh             &
    + 2.0_rprec*v_avg(i,j,k)*w_avg(i,j,k)*dwdz_avg(i,j,k)
Cwudwdx = conjg(a(i,j,k)%wh_udwdxh)           - bb(i,j,k)%uw * dwdx_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dwdxh) * u_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%udwdxh             &
    + 2.0_rprec*w_avg(i,j,k)*u_avg(i,j,k)*dwdx_avg(i,j,k)
Cwvdwdy = conjg(a(i,j,k)%wh_vdwdyh)           - bb(i,j,k)%vw * dwdy_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dwdyh) * v_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%vdwdyh             &
    + 2.0_rprec*w_avg(i,j,k)*v_avg(i,j,k)*dwdy_avg(i,j,k)
Cwwdwdz = conjg(a(i,j,k)%wh_wdwdzh)           - bb(i,j,k)%ww * dwdz_avg(i,1,k)             &
    - conjg(a(i,j,k)%wh_dwdzh) * w_avg(i,1,k) - w_avg(i,1,k) * a(i,j,k)%wdwdzh             &
    + 2.0_rprec*w_avg(i,j,k)*w_avg(i,j,k)*dwdz_avg(i,j,k)

! Pressure-Velocity Gradient correlation, p*dujdxk
Cpdudx = a(i,j,k)%pdudx - a(i,j,k)%p * conjg( a(i,j,k)%dudx )
Cpdudy = a(i,j,k)%pdudy - a(i,j,k)%p * conjg( a(i,j,k)%dudy )
Cpdudz = a(i,j,k)%pdudz - a(i,j,k)%p * conjg( a(i,j,k)%dudz )
Cpdvdx = a(i,j,k)%pdvdx - a(i,j,k)%p * conjg( a(i,j,k)%dvdx )
Cpdvdy = a(i,j,k)%pdvdy - a(i,j,k)%p * conjg( a(i,j,k)%dvdy )
Cpdvdz = a(i,j,k)%pdvdz - a(i,j,k)%p * conjg( a(i,j,k)%dvdz )
Cpdwdx = a(i,j,k)%pdwdx - a(i,j,k)%p * conjg( a(i,j,k)%dwdx )
Cpdwdy = a(i,j,k)%pdwdy - a(i,j,k)%p * conjg( a(i,j,k)%dwdy )
Cpdwdz = a(i,j,k)%pdwdz - a(i,j,k)%p * conjg( a(i,j,k)%dwdz )

! Velocity-Pressure Gradient correlation, uj*dpdxk
Cudpdx = a(i,j,k)%udpdx - bb(i,j,k)%uf * conjg( a(i,j,k)%dpdx )
Cudpdy = a(i,j,k)%udpdy - bb(i,j,k)%uf * conjg( a(i,j,k)%dpdy )
Cudpdz = a(i,j,k)%udpdz - bb(i,j,k)%uf * conjg( a(i,j,k)%dpdz )
Cvdpdx = a(i,j,k)%vdpdx - bb(i,j,k)%vf * conjg( a(i,j,k)%dpdx )
Cvdpdy = a(i,j,k)%vdpdy - bb(i,j,k)%vf * conjg( a(i,j,k)%dpdy )
Cvdpdz = a(i,j,k)%vdpdz - bb(i,j,k)%vf * conjg( a(i,j,k)%dpdz )
Cwdpdx = a(i,j,k)%wdpdx - bb(i,j,k)%wf * conjg( a(i,j,k)%dpdx )
Cwdpdy = a(i,j,k)%wdpdy - bb(i,j,k)%wf * conjg( a(i,j,k)%dpdy )
Cwdpdz = a(i,j,k)%wdpdz - bb(i,j,k)%wf * conjg( a(i,j,k)%dpdz )

! velGrad-velGrad correlation, duidxk*dujdxk, i=j
Cuxux = a(i,j,k)%uxux - a(i,j,k)%dudx * conjg(a(i,j,k)%dudx)
Cuyuy = a(i,j,k)%uyuy - a(i,j,k)%dudy * conjg(a(i,j,k)%dudy)
Cuzuz = a(i,j,k)%uzuz - a(i,j,k)%dudz * conjg(a(i,j,k)%dudz)
Cvxvx = a(i,j,k)%vxvx - a(i,j,k)%dvdx * conjg(a(i,j,k)%dvdx)
Cvyvy = a(i,j,k)%vyvy - a(i,j,k)%dvdy * conjg(a(i,j,k)%dvdy)
Cvzvz = a(i,j,k)%vzvz - a(i,j,k)%dvdz * conjg(a(i,j,k)%dvdz)
Cwxwx = a(i,j,k)%wxwx - a(i,j,k)%dwdx * conjg(a(i,j,k)%dwdx)
Cwywy = a(i,j,k)%wywy - a(i,j,k)%dwdy * conjg(a(i,j,k)%dwdy)
Cwzwz = a(i,j,k)%wzwz - a(i,j,k)%dwdz * conjg(a(i,j,k)%dwdz)

! velGrad-velGrad correlation, duidxk*dujdxk, i\=j
Cuxvx = a(i,j,k)%uxvx - a(i,j,k)%dudx * conjg(a(i,j,k)%dvdx)
Cuyvy = a(i,j,k)%uyvy - a(i,j,k)%dudy * conjg(a(i,j,k)%dvdy)
Cuzvz = a(i,j,k)%uzvz - a(i,j,k)%dudz * conjg(a(i,j,k)%dvdz)
Cuxwx = a(i,j,k)%uxwx - a(i,j,k)%dudx * conjg(a(i,j,k)%dwdx)
Cuywy = a(i,j,k)%uywy - a(i,j,k)%dudy * conjg(a(i,j,k)%dwdy)
Cuzwz = a(i,j,k)%uzwz - a(i,j,k)%dudz * conjg(a(i,j,k)%dwdz)
Cvxwx = a(i,j,k)%vxwx - a(i,j,k)%dvdx * conjg(a(i,j,k)%dwdx)
Cvywy = a(i,j,k)%vywy - a(i,j,k)%dvdy * conjg(a(i,j,k)%dwdy)
Cvzwz = a(i,j,k)%vzwz - a(i,j,k)%dvdz * conjg(a(i,j,k)%dwdz)

! vel-Laplacian correlation, ui*lap(uj)
Culapu = a(i,j,k)%ulapu - bb(i,j,k)%uf * conjg(a(i,j,k)%lapu)
Culapv = a(i,j,k)%ulapv - bb(i,j,k)%uf * conjg(a(i,j,k)%lapv)
Culapw = a(i,j,k)%ulapw - bb(i,j,k)%uf * conjg(a(i,j,k)%lapw)
Cvlapu = a(i,j,k)%vlapu - bb(i,j,k)%vf * conjg(a(i,j,k)%lapu)
Cvlapv = a(i,j,k)%vlapv - bb(i,j,k)%vf * conjg(a(i,j,k)%lapv)
Cvlapw = a(i,j,k)%vlapw - bb(i,j,k)%vf * conjg(a(i,j,k)%lapw)
Cwlapu = a(i,j,k)%wlapu - bb(i,j,k)%wf * conjg(a(i,j,k)%lapu)
Cwlapv = a(i,j,k)%wlapv - bb(i,j,k)%wf * conjg(a(i,j,k)%lapv)
Cwlapw = a(i,j,k)%wlapw - bb(i,j,k)%wf * conjg(a(i,j,k)%lapw)

! ------------------------- Compute Budget Terms -----------------------------
! Mean Advection, Uk*dRijdxk
! Multiplying by kx=0 mode time-average here
e(i,j,k) % advxx = u_avg(i,1,k)*real(conjg(Cududx) + Cududx)       &
    + v_avg(i,1,k)*real(conjg(Cududy) + Cududy)                    &
    + w_avg(i,1,k)*real(conjg(Cududz) + Cududz)
e(i,j,k) % advyy = u_avg(i,1,k)*real(conjg(Cvdvdx) + Cvdvdx)       &
    + v_avg(i,1,k)*real(conjg(Cvdvdy) + Cvdvdy)                    &
    + w_avg(i,1,k)*real(conjg(Cvdvdz) + Cvdvdz)
e(i,j,k) % advzz = u_avg(i,1,k)*real(conjg(Cwdwdx) + Cwdwdx)       &
    + v_avg(i,1,k)*real(conjg(Cwdwdy) + Cwdwdy)                    &
    + w_avg(i,1,k)*real(conjg(Cwdwdz) + Cwdwdz)
e(i,j,k) % advxy = u_avg(i,1,k)*real(conjg(Cvdudx) + Cudvdx)       &
    + v_avg(i,1,k)*real(conjg(Cvdudy) + Cudvdy)                    &
    + w_avg(i,1,k)*real(conjg(Cvdudz) + Cudvdz)
e(i,j,k) % advxz = u_avg(i,1,k)*real(conjg(Cwdudx) + Cudwdx)       &
    + v_avg(i,1,k)*real(conjg(Cwdudy) + Cudwdy)                    & 
    + w_avg(i,1,k)*real(conjg(Cwdudz) + Cudwdz)
e(i,j,k) % advyz = u_avg(i,1,k)*real(conjg(Cwdvdx) + Cvdwdx)       &
    + v_avg(i,1,k)*real(conjg(Cwdvdy) + Cvdwdy)                    &
    + w_avg(i,1,k)*real(conjg(Cwdvdz) + Cvdwdz)

! Transport by fluctuations, d(ui*uj*uk)dxk = ui*uk*dujdxk + uj*uk*duidxk
e(i,j,k) % tflucxx = real( Cuududx + conjg(Cuududx) +               &
    Cuvdudy + conjg(Cuvdudy) + Cuwdudz + conjg(Cuwdudz) )
e(i,j,k) % tflucyy = real( Cvudvdx + conjg(Cvudvdx) +               &
    Cvvdvdy + conjg(Cvvdvdy) + Cvwdvdz + conjg(Cvwdvdz) )
e(i,j,k) % tfluczz = real( Cwudwdx + conjg(Cwudwdx) +               &
    Cwvdwdy + conjg(Cwvdwdy) + Cwwdwdz + conjg(Cwwdwdz) )
e(i,j,k) % tflucxy = real( Cvududx + conjg(Cuudvdx) +               &
    Cvvdudy + conjg(Cuvdvdy) + Cvwdudz + conjg(Cuwdvdz) )
e(i,j,k) % tflucxz = real( Cwududx + conjg(Cuudwdx) +               &
    Cwvdudy + conjg(Cuvdwdy) + Cwwdudz + conjg(Cuwdwdz) )
e(i,j,k) % tflucyz = real( Cwudvdx + conjg(Cvudwdx) +               &
    Cwvdvdy + conjg(Cvvdwdy) + Cwwdvdz + conjg(Cvwdwdz) )

! Pressure-Strain, 2*p*sij = p*(dujdxi + duidxj)
! y-derivative terms have no transport, only pressure-strain
! xz-derivative terms have both transport and strain components
! after having applied the product rule
e(i,j,k) % pstrainxx = -real( Cpdudx + conjg(Cpdudx) )
e(i,j,k) % pstrainyy = real( conjg(Cvdpdy) + Cvdpdy )
e(i,j,k) % pstrainzz = -real( Cpdwdz + conjg(Cpdwdz) )
e(i,j,k) % pstrainxy = real( Cudpdy - Cpdvdx )
e(i,j,k) % pstrainxz = -real(Cpdwdx + conjg(Cpdudz) )
e(i,j,k) % pstrainyz = real( conjg(Cwdpdy) - conjg(Cpdvdz) )

! Transport by pressure-velocity, d(p*(Iik+ui*Ijk))dxk = d(p*uj)dxi + d(p*ui)dxj
! Terms that involve ddy have been transformed and contribute as a source 
! term not a transport term and are therefore included in the pressure-strain
e(i,j,k) % tpresxx = real( Cpdudx + conjg(Cpdudx) + conjg(Cudpdx) + Cudpdx )
e(i,j,k) % tpresyy = 0.0_rprec
e(i,j,k) % tpreszz = real( Cpdwdz + conjg(Cpdwdz) + conjg(Cwdpdz) + Cwdpdz )
e(i,j,k) % tpresxy = real( Cpdvdx + conjg(Cvdpdx) )
e(i,j,k) % tpresxz = real( Cpdwdx + conjg(Cpdudz) + conjg(Cwdpdx) + Cudpdz )
e(i,j,k) % tpresyz = real( conjg(Cpdvdz) + Cvdpdz )

! Production Rate, -(Rik*dUjdxk + Rjk*dUidxk)
! Multiplying by kx=0 mode time-average here
e(i,j,k) % prodxx = 2.0_rprec*( b(i,j,k)%upup * dudx_avg(i,1,k) +                        &
    b(i,j,k)%upvp * dudy_avg(i,1,k) + b(i,j,k)%upwp * dudz_avg(i,1,k) )
e(i,j,k) % prodyy = 2.0_rprec*( b(i,j,k)%upvp * dvdx_avg(i,1,k) +                        &
    b(i,j,k)%vpvp * dvdy_avg(i,1,k) + b(i,j,k)%vpwp * dvdz_avg(i,1,k) )
e(i,j,k) % prodzz = 2.0_rprec*( b(i,j,k)%upwp * dwdx_avg(i,1,k) +                        &
    b(i,j,k)%vpwp * dwdy_avg(i,1,k) + b(i,j,k)%wpwp * dwdz_avg(i,1,k) )
e(i,j,k) % prodxy = (b(i,j,k)%upup * dvdx_avg(i,1,k) + b(i,j,k)%upvp*dudx_avg(i,1,k) +   &
    b(i,j,k)%upvp * dvdy_avg(i,1,k) + b(i,j,k)%vpvp * dudy_avg(i,1,k) +                  &
    b(i,j,k)%upwp * dvdz_avg(i,1,k) + b(i,j,k)%vpwp * dudz_avg(i,1,k) )
e(i,j,k) % prodxz = (b(i,j,k)%upup * dwdx_avg(i,1,k) + b(i,j,k)%upwp*dudx_avg(i,1,k) +   &
    b(i,j,k)%upvp * dwdy_avg(i,1,k) + b(i,j,k)%vpwp * dudy_avg(i,1,k) +                  &
    b(i,j,k)%upwp * dwdz_avg(i,1,k) + b(i,j,k)%wpwp * dudz_avg(i,1,k) )
e(i,j,k) % prodyz = (b(i,j,k)%upvp * dwdx_avg(i,1,k) + b(i,j,k)%upwp*dvdx_avg(i,1,k) +   &
    b(i,j,k)%vpvp * dwdy_avg(i,1,k) + b(i,j,k)%vpwp * dvdy_avg(i,1,k) +                  &
    b(i,j,k)%vpwp * dwdz_avg(i,1,k) + b(i,j,k)%wpwp * dvdz_avg(i,1,k) )

! Pseudo-dissipation, 2*nu*duidxk*dujdxk
e(i,j,k) % pdissxx = 2.0_rprec*nu_molec*real(Cuxux + Cuyuy + Cuzuz)
e(i,j,k) % pdissyy = 2.0_rprec*nu_molec*real(Cvxvx + Cvyvy + Cvzvz)
e(i,j,k) % pdisszz = 2.0_rprec*nu_molec*real(Cwxwx + Cwywy + Cwzwz)
e(i,j,k) % pdissxy = 2.0_rprec*nu_molec*real(Cuxvx + Cuyvy + Cuzvz)
e(i,j,k) % pdissxz = 2.0_rprec*nu_molec*real(Cuxwx + Cuywy + Cuzwz)
e(i,j,k) % pdissyz = 2.0_rprec*nu_molec*real(Cvxwx + Cvywy + Cvzwz)

! Transport by viscous diffusion
! Cuilapuj already multiplied by viscosity because divtxj was used
! Change sign of Cuilapuj because of how divtxj is set up
! Note: lapui terms include x-deriveration terms from pseudo-dissipation
!       which need to be canceled, which is why Cuxux terms have different sign
e(i,j,k) % tviscxx = 2.0_rprec*nu_molec*real(-Cuxux + Cuyuy + Cuzuz) - real( conjg(Culapu) + Culapu )
e(i,j,k) % tviscyy = 2.0_rprec*nu_molec*real(-Cvxvx + Cvyvy + Cvzvz) - real( conjg(Cvlapv) + Cvlapv )
e(i,j,k) % tvisczz = 2.0_rprec*nu_molec*real(-Cwxwx + Cwywy + Cwzwz) - real( conjg(Cwlapw) + Cwlapw )
e(i,j,k) % tviscxy = 2.0_rprec*nu_molec*real(-Cuxvx + Cuyvy + Cuzvz) - real( conjg(Cvlapu) + Culapv )
e(i,j,k) % tviscxz = 2.0_rprec*nu_molec*real(-Cuxwx + Cuywy + Cuzwz) - real( conjg(Cwlapu) + Culapw )
e(i,j,k) % tviscyz = 2.0_rprec*nu_molec*real(-Cvxwx + Cvywy + Cvzwz) - real( conjg(Cwlapv) + Cvlapw )

end do
end do
end do

! ---------------------------- Prepare Output --------------------------------
! Move all terms to the RHS
e % advxx = - e % advxx
e % advyy = - e % advyy
e % advzz = - e % advzz
e % advxy = - e % advxy
e % advxz = - e % advxz
e % advyz = - e % advyz

e % tflucxx = - e % tflucxx
e % tflucyy = - e % tflucyy
e % tfluczz = - e % tfluczz
e % tflucxy = - e % tflucxy
e % tflucxz = - e % tflucxz
e % tflucyz = - e % tflucyz

e % prodxx = - e % prodxx
e % prodyy = - e % prodyy
e % prodzz = - e % prodzz
e % prodxy = - e % prodxy
e % prodxz = - e % prodxz
e % prodyz = - e % prodyz

e % pstrainxx = - e % pstrainxx
e % pstrainyy = - e % pstrainyy
e % pstrainzz = - e % pstrainzz
e % pstrainxy = - e % pstrainxy
e % pstrainxz = - e % pstrainxz
e % pstrainyz = - e % pstrainyz

e % tpresxx = - e % tpresxx
e % tpresyy = - e % tpresyy
e % tpreszz = - e % tpreszz
e % tpresxy = - e % tpresxy
e % tpresxz = - e % tpresxz
e % tpresyz = - e % tpresyz

e % pdissxx = - e % pdissxx
e % pdissyy = - e % pdissyy
e % pdisszz = - e % pdisszz
e % pdissxy = - e % pdissxy
e % pdissxz = - e % pdissxz
e % pdissyz = - e % pdissyz

end function specbudgy_compute

#endif

end module stat_defs
