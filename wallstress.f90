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
subroutine wallstress
!*******************************************************************************
! 
! This subroutine calculates the wall stress txz, tyz (w-nodes) and dudz,
! dvdz (w-nodes) at the first z-location k = 1. The wall stress is calculated
! depending on lower boundary condition lbc_mom. This subroutine should only
! be called after ensuring coord==0
!
! Options for lbc_mom:
!   0 - stress free
!       txz, tyz, dudz, and dvdz are all 0
!
!   1 - DNS wall boundary conditions 
!       calculates wall stress values from the first grid point
!
!   2 - Equilibirum wall model
!       See John D. Albertson's dissertation, eqns (2.46)-(2.52)
!       Also see E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent 
!           Lagrangian dynamic model for large eddy simulation of complex 
!           turbulent flows" (2005) -- Appendix
!
!   3 - Integral wall model
!       See X.I.A. Yang, J. Sadique, R. Mittal & C. Meneveau, "Integral wall 
!           model for large eddy simulations of wall-bounded turbulent flows." (2015)
! 
!   4 - Smooth Equilibrium wall model with blending function
! 

use types, only : rprec
use param, only : lbc_mom
use param, only : ubc_mom, coord, nproc, nz ! these necessary only for upper bc
use messages, only : error
use iwmles, only : iwm_wallstress
use tlwmles, only : tlwm
use sim_param, only : txz, tyz, dudz, dvdz
implicit none
character(*), parameter :: sub_name = 'wallstress'

! Lower boundary condition
if (coord == 0) then
    select case (lbc_mom)
        ! Stress free
        case (0)
            call ws_free_lbc

        ! DNS wall
        case (1)
            call ws_dns_lbc

        ! Equilibrium wall model
        case (2)
            call ws_equilibrium_lbc

        ! Integral wall model (not implemented for top wall)
        case (3)
            call iwm_wallstress()

        ! Smooth equilibrium wall model
        case (4)
            call ws_smoothequil_lbc

        ! two-layer equilibrium wall model
        case (5)
            call tlwm()

        ! two-layer nonequilibrium wall model
        case (6)
            call tlwm()

        ! two-layer equilibrium wall model in (kx,ky) space
        case (7)
            call tlwm()

        ! two-layer nonequilibrium wall model in (kx,ky) space
        case (8)
            call tlwm()

        ! Equilibrium wall model for RNL-LES, assumes horizontal average
        case (9)
            call ws_equilibrium_lbc_fourier_xyavg

        ! Equilibrium wall model for RNL-LES, assumes streamwise average
        case (10)
            call ws_equilibrium_lbc_fourier_xavg

        ! Otherwise, invalid
        case default
            call error (sub_name, 'invalid lbc_mom')
    end select
end if

if (coord == nproc-1) then
    select case (ubc_mom)
        ! Stress free
        case (0)                        
            call ws_free_ubc

        ! DNS wall
        case (1)                        
            call ws_dns_ubc

        ! Equilibrium wall model
        case (2)                       
            call ws_equilibrium_ubc

        ! Integral wall model (not implemented for top wall)
        case (3)                        
            call error(sub_name, 'invalid ubc_mom')

        ! Otherwise, invalid
        case default       
            call error(sub_name, 'invalid ubc_mom')
    end select
end if

contains

!*******************************************************************************
subroutine ws_free_lbc
!*******************************************************************************
implicit none

txz(:, :, 1) = 0._rprec
tyz(:, :, 1) = 0._rprec
dudz(:, :, 1) = 0._rprec
dvdz(:, :, 1) = 0._rprec

end subroutine ws_free_lbc

!*******************************************************************************
subroutine ws_free_ubc
!*******************************************************************************
implicit none

txz(:, :,nz) = 0._rprec 
tyz(:, :,nz) = 0._rprec
dudz(:,:,nz) = 0._rprec
dvdz(:,:,nz) = 0._rprec

end subroutine ws_free_ubc

!*******************************************************************************
subroutine ws_dns_lbc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, ubot 
#ifdef PPMAPPING
use sim_param, only : mesh_stretch
#else
use param, only : dz
#endif
use sim_param , only : u, v
implicit none
integer :: i, j

do j = 1, ny
    do i = 1, nx
#ifdef PPMAPPING
        dudz(i,j,1) = ( u(i,j,1) - ubot ) / (mesh_stretch(1))
        dvdz(i,j,1) = v(i,j,1) / (mesh_stretch(1))
#else
        dudz(i,j,1) = ( u(i,j,1) - ubot ) / ( 0.5_rprec*dz )
        dvdz(i,j,1) = v(i,j,1) / ( 0.5_rprec*dz )
#endif
        txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
        tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)
    end do
end do

end subroutine ws_dns_lbc

!*******************************************************************************
subroutine ws_dns_ubc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, utop, L_z
#ifdef PPMAPPING
use sim_param, only : mesh_stretch
#else
use param, only : dz
#endif
use sim_param , only : u, v
implicit none
integer :: i, j

do j = 1, ny
    do i = 1, nx
#ifdef PPMAPPING
        dudz(i,j,nz) = ( utop - u(i,j,nz-1) ) / (L_z - mesh_stretch(nz-1))
        dvdz(i,j,nz) = -v(i,j,nz-1) / (L_z - mesh_stretch(nz-1))
#else
        dudz(i,j,nz) = ( utop - u(i,j,nz-1) ) / (0.5_rprec*dz)
        dvdz(i,j,nz) = -v(i,j,nz-1) / (0.5_rprec*dz)
#endif
        txz(i,j,nz) = -nu_molec/(z_i*u_star)*dudz(i,j,nz)
        tyz(i,j,nz) = -nu_molec/(z_i*u_star)*dvdz(i,j,nz)
    end do
end do

end subroutine ws_dns_ubc

!*******************************************************************************
subroutine ws_equilibrium_lbc
!*******************************************************************************
use param, only : dz, ld, nx, ny, vonk, zo
use sim_param, only : u, v
#ifdef PPMAPPING
use sim_param, only : jaco_uv
#endif
use test_filtermodule
implicit none
integer :: i, j
real(rprec), dimension(nx, ny) :: denom, u_avg, ustar
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const, const2

u1 = u(:,:,1)
v1 = v(:,:,1)
call test_filter(u1)
call test_filter(v1)

#ifdef PPMAPPING
denom = log(0.5_rprec*jaco_uv(1)*dz/zo)
const2 = jaco_uv(1)
#else
denom = log(0.5_rprec*dz/zo)
const2 = 1._rprec
#endif
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
ustar = u_avg*vonk/denom

do j = 1, ny
    do i = 1, nx
        const = -(ustar(i,j)**2)/u_avg(i,j)
        txz(i,j,1) = const*u1(i,j)
        tyz(i,j,1) = const*v1(i,j)
        !this is as in Moeng 84
        dudz(i,j,1) = ustar(i,j)/(0.5_rprec*const2*dz*vonK)*u(i,j,1)/u_avg(i,j)
        dvdz(i,j,1) = ustar(i,j)/(0.5_rprec*const2*dz*vonK)*v(i,j,1)/u_avg(i,j)

        dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
        dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_lbc

!*******************************************************************************
subroutine ws_equilibrium_lbc_fourier_xyavg
!*******************************************************************************
! 
! This subroutine is only to be used if sgs = true, molec = false, and
! when fourier = true. This routine mimics ws_equilibrium_lbc, but with
! some additional assumptions since working in Fourier space.
! 
! This is the wall-model used in Bretheim et al. 2018, and assumes averaging
! in the streamwise and spanwise direction for the prefactor, the outputted
! wallstress is in (kx,y) space because of u1 and v1.
! 
! Note the routine assumes a uniform grid, i.e. mapping not used.
! 
use param, only : dz, ld, nx, ny, vonk, zo, ld_big, ny2
use sim_param, only : u, v
use derivatives, only : dft_direct_back_2d_n_yonlyC_big, dft_direct_forw_2d_n_yonlyC_big
use derivatives, only : convolve_rnl
use fft, only : padd, unpadd
!use test_filtermodule
implicit none
integer :: i, j
real(rprec), dimension(ld, ny) :: u1, v1, u_avg
!real(rprec), dimension(ld, ny) :: ustar
real(rprec), dimension(ld_big,ny2) :: u1_big, v1_big, s1_big, txz_big, tyz_big
real(rprec) :: const, denom

u1 = u(:,:,1)
v1 = v(:,:,1)
! Not test filtering for fourier
!call test_filter(u1)
!call test_filter(v1)
denom = log(0.5_rprec*dz/zo)

! Compute u_avg = sqrt( u**2 + v**2 ) and ustar
u1_big = 0._rprec
v1_big = 0._rprec
s1_big = 0._rprec
call padd(u1_big(:,:), u1(:,:))
call padd(v1_big(:,:), v1(:,:))
call dft_direct_back_2d_n_yonlyC_big(u1_big(:,:))
call dft_direct_back_2d_n_yonlyC_big(v1_big(:,:))
s1_big = convolve_rnl(u1_big(:,:),u1_big(:,:)) + convolve_rnl(v1_big(:,:),v1_big(:,:))
call dft_direct_forw_2d_n_yonlyC_big(s1_big(:,:))
call unpadd(u_avg(:,:), s1_big(:,:))
! Use horizontal average of u_avg, zero out kx & ky varying modes
u_avg(1,1) = sqrt( u_avg(1,1) )
u_avg(2:ld,1:ny) = 0._rprec
u_avg(1,2:ny) = 0._rprec
!ustar = u_avg*vonk/denom ! not used, already made simplifications below

! More simplified here than original ws_equilibrium_lbc because of above assumptions
const = -((vonk/denom)**2) * u_avg(1,1)

! Rewrite s1_big variable to convolve const with u1/v1 to find txz/tyz
s1_big = 0._rprec
s1_big(1,1) = const
call dft_direct_back_2d_n_yonlyC_big(s1_big(:,:))

! Compute txz & tyz
txz_big(:,:) = convolve_rnl(s1_big(:,:), u1_big(:,:))
tyz_big(:,:) = convolve_rnl(s1_big(:,:), v1_big(:,:))
call dft_direct_forw_2d_n_yonlyC_big(txz_big(:,:))
call dft_direct_forw_2d_n_yonlyC_big(tyz_big(:,:))
call unpadd(txz(:,:,1), txz_big(:,:))
call unpadd(tyz(:,:,1), tyz_big(:,:))

! Compute dudz & dvdz, again more simplified than original ws_equilibrium_lbc
dudz(:,:,1) = 1._rprec/(0.5_rprec*dz*denom) * u(:,:,1)
dvdz(:,:,1) = 1._rprec/(0.5_rprec*dz*denom) * v(:,:,1)

do j = 1, ny
    do i = 1, nx
        dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
        dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_lbc_fourier_xyavg

!*******************************************************************************
subroutine ws_equilibrium_lbc_fourier_xavg
!*******************************************************************************
! 
! Similar to ws_equilibrium_lbc_fourier_xyavg, this wall-model is only for
! sgs = true, molec = false, and fourier = true. Instead of assuming a 
! streamwise and spanwise average, only a streamwise average is used.
! 
! Note the routine assumes a uniform grid, i.e. mapping not used.
! 
use param, only : dz, ld, nx, ny, vonk, zo
use sim_param, only : u, v
use derivatives, only : dft_direct_back_2d_n_yonlyC
use derivatives, only : dft_direct_forw_2d_n_yonlyC
!use test_filtermodule, only : test_filter_fourier
implicit none
integer :: i, j
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec), dimension(ny) :: u1_avg, v1_avg, u_avg, ustar
real(rprec) :: denom, const, const2

u1 = u(:,:,1)
v1 = v(:,:,1)

! Instead of test-filtering, use streamwise average
! Transform (kx,ky) --> (kx,y)
call dft_direct_back_2d_n_yonlyC(u1(:,:))
call dft_direct_back_2d_n_yonlyC(v1(:,:))

! Test-filter in only the spanwise direction
!call test_filter_fourier( u1 )
!call test_filter_fourier( v1 )

! Take only kx=0 mode, effectively streamwise averaging
u1_avg(1:ny) = u1(1,:)
v1_avg(1:ny) = v1(1,:)

! Compute prefactor
denom = log(0.5_rprec*dz/zo)
u_avg = sqrt(u1_avg(1:ny)**2+v1_avg(1:ny)**2)
ustar = u_avg*vonk/denom

! remember ustar(y) and u_avg(y), but u1(kx,y) and v1(kx,y)
do j = 1, ny
    const = -(ustar(j)**2)/u_avg(j)
    const2 = ustar(j)/(0.5_rprec*dz*vonK)/u_avg(j)
    do i = 1, nx
        ! Using u1 & v1, still has non-zero kx modes
        txz(i,j,1) = const*u1(i,j)
        tyz(i,j,1) = const*v1(i,j)
        !this is as in Moeng 84
        dudz(i,j,1) = const2*u1(i,j)
        dvdz(i,j,1) = const2*v1(i,j)

        dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
        dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
    end do
end do

! Transform (kx,y) --> (kx,ky)
call dft_direct_forw_2d_n_yonlyC(txz(:,:,1))
call dft_direct_forw_2d_n_yonlyC(tyz(:,:,1))
call dft_direct_forw_2d_n_yonlyC(dudz(:,:,1))
call dft_direct_forw_2d_n_yonlyC(dvdz(:,:,1))

end subroutine ws_equilibrium_lbc_fourier_xavg

!*******************************************************************************
subroutine ws_equilibrium_ubc
!*******************************************************************************
use param, only : dz, ld, nx, ny, vonk, zo
use sim_param, only : u, v
use test_filtermodule
implicit none
integer :: i, j
real(rprec), dimension(nx, ny) :: denom, u_avg, ustar
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const


u1 = u(:,:,nz-1)
v1 = v(:,:,nz-1)
call test_filter(u1)
call test_filter(v1)
denom = log(0.5_rprec*dz/zo)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
ustar = u_avg*vonk/denom

do j = 1, ny
    do i = 1, nx
        const = (ustar(i,j)**2)/u_avg(i,j) ! diff sign for upper b.c.
        txz(i,j,nz) = const*u1(i,j)
        tyz(i,j,nz) = const*v1(i,j)
        !this is as in Moeng 84
        dudz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*u(i,j,nz-1)/u_avg(i,j)
        dvdz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*v(i,j,nz-1)/u_avg(i,j)
        dudz(i,j,nz) = merge(0._rprec,dudz(i,j,nz),u(i,j,nz-1).eq.0._rprec)
        dvdz(i,j,nz) = merge(0._rprec,dvdz(i,j,nz),v(i,j,nz-1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_ubc

!*******************************************************************************
subroutine ws_smoothequil_lbc
!*******************************************************************************
! This wall model computes the stress at the wall given that the wall is smooth.
! Unlike ws_equilibrium_lbc, this wall model does not require a roughness 
! length. In its place, the viscous length is used. As a result an iterative 
! method (Newton's method) is required to find the wall stress.
! 
! A blending function is also used for low Reynolds number LES cases. The 
! blending function used was proposed in
!       B. Kader. "Temperature and Concentration Profiles in Fully 
!           Turbulent Boundary Layers". Int. J. Heat Mass Transfer. 
!           24(9). 1541-1544. 1981.
! and is the wall model used in commercial software ANSYS FLUENT.
!
! TO DO:
! - Currently using u_star as initial guess for Newton's method, should
!     consider using previous wall stress value 
! - Call this wall model if user specifies zo = 0? This could then be 
!     incorporated with original equilibrium wall model. If molec is 
!     false however (i.e. infinite Reynolds number) then would have to
!     call original equilibrium wall model with given roughness length
! - Using 10 iterations with each calculation. Is this sufficient?
!     This was tested in MATLAB before implementation. Could vary number of 
!     iterations depending on dy_plus
! 

use param, only: dz, ld, nx, ny, vonk, nu_molec, u_star
use sim_param, only : u, v, txz, tyz
use test_filtermodule
implicit none
integer :: i, j, iter, Niter
real(rprec), dimension(nx, ny) :: u_avg, ustar
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const, B0, aa, bb, lam, turb, bigG, numer, y_plus
real(rprec) :: Dlam, Dturb, DbigG1, DbigG2, DbigG, denom1, denom2, denom

u1 = u(:,:,1)
v1 = v(:,:,1)
call test_filter(u1)
call test_filter(v1)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)

! Using previous txz and tyz value as an initial guess
! Consider using u_avg and (u1,v1) as well as (txz,tyz)

! ustar = sqrt( abs(txz(1:nx,1:ny,1) + tyz(1:nx,1:ny,1)) )
ustar = u_star

! Model parameters
B0 = 5.3_rprec
aa = 0.01_rprec
bb = 5.0_rprec

! Iterative solver settings
Niter = 10

! Find stress at wall as before
do j = 1, ny
    do i = 1, nx

        ! Perform Newton's method
        do iter = 1, Niter
            y_plus = (0.5_rprec*dz)*ustar(i,j)/nu_molec

            ! Functions
            lam = ustar(i,j)*y_plus
            turb = ustar(i,j)*( (1/vonk)*log(y_plus) + B0 )
            bigG = - aa*((y_plus)**4)/(1 + bb*y_plus)
            numer = exp(bigG)*lam + exp(1/bigG)*turb - u_avg(i,j)

            ! Analytical derivatives of root function wrt ustar
            Dlam = 2*y_plus
            Dturb = (1/vonk)*log(y_plus) + B0 + 1.0_rprec
            DbigG1 = ((aa*bb*(y_plus**4)*(0.5_rprec*dz/nu_molec))/((bb*y_plus+1)**2))
            DbigG2 = - (4.0_rprec*aa*(y_plus**3)*(0.5_rprec*dz/nu_molec))/(bb*y_plus+1)
            DbigG = DbigG1 + DbigG2
            denom1 = exp(bigG)*Dlam + exp(1/bigG)*Dturb
            denom2 = (exp(bigG)*lam - (exp(1/bigG)/(bigG**2))*turb)*DbigG
            denom = denom1 + denom2

            ! Applying Newton's method
            ustar(i,j) = ustar(i,j) - numer/denom

        end do

        const = -(ustar(i,j)**2)/u_avg(i,j)
        txz(i,j,1) = const*u1(i,j)
        tyz(i,j,1) = const*v1(i,j)

        ! Functions
        lam = ustar(i,j)*y_plus
        turb = ustar(i,j)*( (1/vonk)*log(y_plus) + B0 )
        bigG = - aa*((y_plus)**4)/(1 + bb*y_plus)

        ! Using derivatives of average velocity wrt z
        Dlam = (ustar(i,j)**2)/nu_molec
        Dturb = (ustar(i,j))/(vonk*0.5_rprec*dz)
        DbigG = -aa*(ustar(i,j)**4)*((0.5_rprec*dz)**3)*                       &
            (4.0_rprec*nu_molec+3.0_rprec*bb*ustar(i,j)*0.5_rprec*dz)/         &
            ( (nu_molec**3)*((nu_molec + bb*ustar(i,j)*0.5_rprec*dz)**3) )
        const = exp(bigG)*Dlam + exp(1/bigG)*Dturb +                           &
            lam*exp(bigG)*DbigG - (1/(bigG**2))*exp(1/bigG)*DbigG
        dudz(i,j,1) = -const*u1(i,j)/u_avg(i,j)
        dvdz(i,j,1) = -const*v1(i,j)/u_avg(i,j)
        dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
        dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)

    end do
end do

end subroutine ws_smoothequil_lbc

end subroutine wallstress
