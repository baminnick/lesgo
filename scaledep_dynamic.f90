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
subroutine scaledep_dynamic(Cs_1D)
!*******************************************************************************
!
! Subroutine uses the scale dependent dynamic model to calculate the Smagorinsky
! coefficient Cs_1D and |S|. This is done layer-by-layer to save memory.
!
use types, only : rprec
use param, only : path, ld, nx, ny, nz, coord
use sim_param, only : u, v, w
use test_filtermodule
use sgs_stag_util, only : rtnewt
use sgs_param, only : S11, S12, S13, S22, S23, S33, delta, S, ee_now,          &
    u_bar, v_bar, w_bar, L11, L12, L13, L22, L23, L33, S_bar,                  &
    S11_bar, S12_bar, S13_bar, S22_bar, S23_bar, S33_bar,                      &
    S_S11_bar, S_S12_bar, S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar,          &
    u_hat, v_hat, w_hat,                                                       &
    S_hat, S11_hat, S12_hat, S13_hat, S22_hat, S23_hat, S33_hat,               &
    S_S11_hat, S_S12_hat, S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat
use param, only : fourier
use derivatives, only : dft_direct_back_2d_n_yonlyC
implicit none

integer :: jz
real(rprec), dimension(nz), intent (inout) :: Cs_1D
real(rprec), save, allocatable, target, dimension(:,:) :: Q11, Q12, Q13, Q22,  &
    Q23, Q33
real(rprec), pointer, dimension(:,:) :: M11, M12, M13, M22, M23, M33
real(rprec), save, dimension(:), allocatable :: beta
logical, save :: arrays_allocated = .false.
real(rprec) :: const
real(rprec), dimension(0:5) :: A
real(rprec) :: a1, b1, c1, d1, e1, a2, b2, c2, d2, e2
real(rprec), dimension(ld,ny) :: LM,MM
real(rprec), dimension(ld,ny) :: um, vm, wm

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( Q11(ld,ny), Q12(ld,ny), Q13(ld,ny), &
        Q22(ld,ny), Q23(ld,ny), Q33(ld,ny) )
    allocate ( beta(nz) )
    arrays_allocated = .true.
endif

! Associate pointers
M11 => Q11
M12 => Q12
M13 => Q13
M22 => Q22
M23 => Q23
M33 => Q33

if (fourier) then
do jz = 1, nz
    ! Write interpolated variables into temp variables
    if ( (coord == 0) .and. (jz == 1) ) then
        ! put on uvp-nodes
        ! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
        um(:,:) = u(:,:,1)
        vm(:,:) = v(:,:,1)
        wm(:,:) = 0.25*w(:,:,2) ! parabolic interp.
    else
        ! put on w-nodes
        um(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        vm(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        wm(:,:) = w(:,:,jz)
    endif

    ! Transform variables, ky --> y, to take products
    call dft_direct_back_2d_n_yonlyC( um(:,:) )
    call dft_direct_back_2d_n_yonlyC( vm(:,:) )
    call dft_direct_back_2d_n_yonlyC( wm(:,:) )

    ! Using L_ij as temp storage here
    ! Use streamwise average when taking products in physical space
    L11(1,:) = um(1,:)*um(1,:)
    L12(1,:) = um(1,:)*vm(1,:)
    L13(1,:) = um(1,:)*wm(1,:)
    L22(1,:) = vm(1,:)*vm(1,:)
    L23(1,:) = vm(1,:)*wm(1,:)
    L33(1,:) = wm(1,:)*wm(1,:)
    L11(2:ld,:) = 0.0_rprec
    L12(2:ld,:) = 0.0_rprec
    L13(2:ld,:) = 0.0_rprec
    L22(2:ld,:) = 0.0_rprec
    L23(2:ld,:) = 0.0_rprec
    L33(2:ld,:) = 0.0_rprec
    u_bar(1,:) = um(1,:)
    v_bar(1,:) = vm(1,:)
    w_bar(1,:) = wm(1,:)
    u_bar(2:ld,:) = 0.0_rprec
    v_bar(2:ld,:) = 0.0_rprec
    w_bar(2:ld,:) = 0.0_rprec
    u_hat = u_bar
    v_hat = v_bar
    w_hat = w_bar

    ! Filter first term and add the second term to get the final value
    ! in-place filtering, only on y not x
    call test_filter_fourier ( u_bar )
    call test_filter_fourier ( v_bar )
    call test_filter_fourier ( w_bar )
    call test_filter_fourier ( L11 )
    L11 = L11 - u_bar*u_bar
    call test_filter_fourier ( L12 )
    L12 = L12 - u_bar*v_bar
    call test_filter_fourier ( L13 )
    L13 = L13 - u_bar*w_bar
    call test_filter_fourier ( L22 )
    L22 = L22 - v_bar*v_bar
    call test_filter_fourier ( L23 )
    L23 = L23 - v_bar*w_bar
    call test_filter_fourier ( L33 )
    L33 = L33 - w_bar*w_bar

    Q11 = u_bar*u_bar
    Q12 = u_bar*v_bar
    Q13 = u_bar*w_bar
    Q22 = v_bar*v_bar
    Q23 = v_bar*w_bar
    Q33 = w_bar*w_bar

    call test_test_filter_fourier ( u_hat )
    call test_test_filter_fourier ( v_hat )
    call test_test_filter_fourier ( w_hat )
    call test_test_filter_fourier ( Q11 )
    Q11 = Q11 - u_hat*u_hat
    call test_test_filter_fourier ( Q12 )
    Q12 = Q12 - u_hat*v_hat
    call test_test_filter_fourier ( Q13 )
    Q13 = Q13 - u_hat*w_hat
    call test_test_filter_fourier ( Q22 )
    Q22 = Q22 - v_hat*v_hat
    call test_test_filter_fourier ( Q23 )
    Q23 = Q23 - v_hat*w_hat
    call test_test_filter_fourier ( Q33 )
    Q33 = Q33 - w_hat*w_hat

    ! Remember Sij(kx,y,z) after calc_Sij
    ! calculate |S| using streamwise average
    S(1,:) = sqrt(2._rprec*(S11(1,:,jz)**2 + S22(1,:,jz)**2 +     &
        S33(1,:,jz)**2 + 2._rprec*(S12(1,:,jz)**2 +               &
        S13(1,:,jz)**2 + S23(1,:,jz)**2)))
    S(2:ld,:) = 0.0_rprec

    ! S_ij already on w-nodes
    ! Take streamwise average, zero-out other modes
    S11_bar(1,:) = S11(1,:,jz)
    S12_bar(1,:) = S12(1,:,jz)
    S13_bar(1,:) = S13(1,:,jz)
    S22_bar(1,:) = S22(1,:,jz)
    S23_bar(1,:) = S23(1,:,jz)
    S33_bar(1,:) = S33(1,:,jz)
    S11_bar(2:ld,:) = 0.0_rprec
    S12_bar(2:ld,:) = 0.0_rprec
    S13_bar(2:ld,:) = 0.0_rprec
    S22_bar(2:ld,:) = 0.0_rprec
    S23_bar(2:ld,:) = 0.0_rprec
    S33_bar(2:ld,:) = 0.0_rprec

    S11_hat = S11_bar
    S12_hat = S12_bar
    S13_hat = S13_bar
    S22_hat = S22_bar
    S23_hat = S23_bar
    S33_hat = S33_bar

    call test_filter_fourier ( S11_bar )
    call test_filter_fourier ( S12_bar )
    call test_filter_fourier ( S13_bar )
    call test_filter_fourier ( S22_bar )
    call test_filter_fourier ( S23_bar )
    call test_filter_fourier ( S33_bar )

    call test_test_filter_fourier ( S11_hat )
    call test_test_filter_fourier ( S12_hat )
    call test_test_filter_fourier ( S13_hat )
    call test_test_filter_fourier ( S22_hat )
    call test_test_filter_fourier ( S23_hat )
    call test_test_filter_fourier ( S33_hat )

    S_bar(1,:) = sqrt(2._rprec*(S11_bar(1,:)**2 +               &
        S22_bar(1,:)**2 + S33_bar(1,:)**2 +                     &
        2._rprec*(S12_bar(1,:)**2 + S13_bar(1,:)**2 +           &
        S23_bar(1,:)**2)))
    S_bar(2:ld,:) = 0.0_rprec

    S_hat(1,:) = sqrt(2._rprec*(S11_hat(1,:)**2 +               &
        S22_hat(1,:)**2 + S33_hat(1,:)**2 +                     &
        2._rprec*(S12_hat(1,:)**2 + S13_hat(1,:)**2 +           &
        S23_hat(1,:)**2)))
    S_hat(2:ld,:) = 0.0_rprec

    S_S11_bar(1,:) = S(1,:)*S11(1,:,jz)
    S_S12_bar(1,:) = S(1,:)*S12(1,:,jz)
    S_S13_bar(1,:) = S(1,:)*S13(1,:,jz)
    S_S22_bar(1,:) = S(1,:)*S22(1,:,jz)
    S_S23_bar(1,:) = S(1,:)*S23(1,:,jz)
    S_S33_bar(1,:) = S(1,:)*S33(1,:,jz)
    S_S11_bar(2:ld,:) = 0.0_rprec
    S_S12_bar(2:ld,:) = 0.0_rprec
    S_S13_bar(2:ld,:) = 0.0_rprec
    S_S22_bar(2:ld,:) = 0.0_rprec
    S_S23_bar(2:ld,:) = 0.0_rprec
    S_S33_bar(2:ld,:) = 0.0_rprec

    S_S11_hat = S_S11_bar
    S_S12_hat = S_S12_bar
    S_S13_hat = S_S13_bar
    S_S22_hat = S_S22_bar
    S_S23_hat = S_S23_bar
    S_S33_hat = S_S33_bar

    call test_filter_fourier ( S_S11_bar )
    call test_filter_fourier ( S_S12_bar )
    call test_filter_fourier ( S_S13_bar )
    call test_filter_fourier ( S_S22_bar )
    call test_filter_fourier ( S_S23_bar )
    call test_filter_fourier ( S_S33_bar )

    call test_test_filter_fourier ( S_S11_hat )
    call test_test_filter_fourier ( S_S12_hat )
    call test_test_filter_fourier ( S_S13_hat )
    call test_test_filter_fourier ( S_S22_hat )
    call test_test_filter_fourier ( S_S23_hat )
    call test_test_filter_fourier ( S_S33_hat )

    ! note: check that the Nyquist guys are zero!
    ! Averaging only over spanwise direction here, using the streamwise mean
    a1 = -2._rprec*(delta**2)*4._rprec*sum( S_bar(1,:)*(S11_bar(1,:)*L11(1,:) + S22_bar(1,:)*L22(1,:) + &
        S33_bar(1,:)*L33(1,:) + 2._rprec*( S12_bar(1,:)*L12(1,:) + S13_bar(1,:)*L13(1,:) + S23_bar(1,:)*L23(1,:))))/   &
        (ny)
    b1 = -2._rprec*(delta**2)*sum( S_S11_bar(1,:)*L11(1,:) + S_S22_bar(1,:)*L22(1,:) +             &
        S_S33_bar(1,:)*L33(1,:) + 2._rprec*( S_S12_bar(1,:)*L12(1,:) + S_S13_bar(1,:)*L13(1,:) +             &
        S_S23_bar(1,:)*L23(1,:)))/(ny)
    c1 = (2._rprec*delta**2)**2 * sum( S_S11_bar(1,:)**2 + S_S22_bar(1,:)**2 +           &
        S_S33_bar(1,:)**2 + &
        2._rprec*( S_S12_bar(1,:)**2 + S_S13_bar(1,:)**2 + S_S23_bar(1,:)**2))/(ny)
    d1 = (2._rprec*delta**2)**2 * 16._rprec*sum( 0.5_rprec*S_bar(1,:)**4 )/(ny)
    e1 = 2._rprec*(2._rprec*delta**2)**2*4._rprec*sum( S_bar(1,:)*( S11_bar(1,:)*        &
        S_S11_bar(1,:) +S22_bar(1,:)*S_S22_bar(1,:) + S33_bar(1,:)*S_S33_bar(1,:) + 2._rprec*(          &
        S12_bar(1,:)*S_S12_bar(1,:) + S13_bar(1,:)*S_S13_bar(1,:) + S23_bar(1,:)*S_S23_bar(1,:))))/          &
        (ny)

    a2 = -2._rprec*(delta**2)*16._rprec*sum( S_hat(1,:)*(S11_hat(1,:)*Q11(1,:) + S22_hat(1,:)*Q22(1,:) +&
        S33_hat(1,:)*Q33(1,:) + 2._rprec*( S12_hat(1,:)*Q12(1,:) + S13_hat(1,:)*Q13(1,:) + S23_hat(1,:)*Q23(1,:))))/&
        (ny)
    b2 = -2._rprec*(delta**2)*sum( S_S11_hat(1,:)*Q11(1,:) + S_S22_hat(1,:)*Q22(1,:) +&
        S_S33_hat(1,:)*Q33(1,:) + 2._rprec*( S_S12_hat(1,:)*Q12(1,:) + S_S13_hat(1,:)*Q13(1,:) + &
        S_S23_hat(1,:)*Q23(1,:)))/(ny)
    c2 = (2._rprec*delta**2)**2 * sum( S_S11_hat(1,:)**2 + S_S22_hat(1,:)**2 +&
        S_S33_hat(1,:)**2 + &
        2._rprec*( S_S12_hat(1,:)**2 + S_S13_hat(1,:)**2 + S_S23_hat(1,:)**2))/(ny)
    d2 = (2._rprec*delta**2)**2 *256._rprec*sum( 0.5_rprec*S_hat(1,:)**4 )/(ny)
    e2 = 2._rprec*(2._rprec*delta**2)**2*16._rprec*sum( S_hat(1,:)*( S11_hat(1,:)*       &
        S_S11_hat(1,:) + S22_hat(1,:)*S_S22_hat(1,:) + S33_hat(1,:)*S_S33_hat(1,:) + 2.*(               &
        S12_hat(1,:)*S_S12_hat(1,:) + S13_hat(1,:)*S_S13_hat(1,:) + S23_hat(1,:)*S_S23_hat(1,:))))/          &
        (ny)

    A(0) = b2*c1 - b1*c2
    A(1) = a1*c2 - b2*e1
    A(2) = b2*d1 + b1*e2 - a2*c1
    A(3) = a2*e1 - a1*e2
    A(4) = -a2*d1 - b1*d2
    A(5) = a1*d2

    beta(jz) = rtnewt(A,jz)

    ! now put beta back into M_ij: using Q_ij as storage
    const = 2._rprec*delta**2
    M11(1,:) = const*(S_S11_bar(1,:) - 4._rprec*beta(jz)*S_bar(1,:)*S11_bar(1,:))
    M12(1,:) = const*(S_S12_bar(1,:) - 4._rprec*beta(jz)*S_bar(1,:)*S12_bar(1,:))
    M13(1,:) = const*(S_S13_bar(1,:) - 4._rprec*beta(jz)*S_bar(1,:)*S13_bar(1,:))
    M22(1,:) = const*(S_S22_bar(1,:) - 4._rprec*beta(jz)*S_bar(1,:)*S22_bar(1,:))
    M23(1,:) = const*(S_S23_bar(1,:) - 4._rprec*beta(jz)*S_bar(1,:)*S23_bar(1,:))
    M33(1,:) = const*(S_S33_bar(1,:) - 4._rprec*beta(jz)*S_bar(1,:)*S33_bar(1,:))

    Cs_1D(jz) =                                                                   &
        sum(L11(1,:)*M11(1,:) + L22(1,:)*M22(1,:) + L33(1,:)*M33(1,:) +           &
        2._rprec*(L12(1,:)*M12(1,:) + L13(1,:)*M13(1,:) + L23(1,:)*M23(1,:))) /   &
        sum(M11(1,:)**2 + M22(1,:)**2 + M33(1,:)**2 +                             &
        2._rprec*(M12(1,:)**2 + M13(1,:)**2 + M23(1,:)**2))
    Cs_1D(jz) = max(0._rprec, real(Cs_1D(jz),rprec))

    ! Calculate ee_now (the current value of eij*eij)
    LM(1,:) = L11(1,:)*M11(1,:) + L22(1,:)*M22(1,:) + L33(1,:)*M33(1,:) +         &
        2._rprec*(L12(1,:)*M12(1,:) + L13(1,:)*M13(1,:) + L23(1,:)*M23(1,:))
    LM(2:ld,:) = 0.0_rprec
    MM(1,:) = M11(1,:)**2 + M22(1,:)**2 + M33(1,:)**2 +                           &
        2._rprec*(M12(1,:)**2 + M13(1,:)**2 + M23(1,:)**2)
    MM(2:ld,:) = 0.0_rprec
    ee_now(1,:,jz) = L11(1,:)**2 + L22(1,:)**2 + L33(1,:)**2 +                    &
        2._rprec*(L12(1,:)**2 + L13(1,:)**2 + L23(1,:)**2) -                      &
        2._rprec*LM(1,:)*Cs_1D(jz) + MM(1,:)*Cs_1D(jz)**2
    ee_now(2:ld,:,jz) = 0.0_rprec

end do
else !! .not. fourier
do jz = 1, nz
    ! using L_ij as temp storage here
    if ( (coord == 0) .and. (jz == 1) ) then
        !!! watch the 0.25's: w = c*z**z near wall, so get 0.25
        ! put on uvp-nodes
        L11(:,:) = u(:,:,1)*u(:,:,1)  ! uv-node
        L12(:,:) = u(:,:,1)*v(:,:,1)  ! uv-node
        L13(:,:) = u(:,:,1)*0.25_rprec*w(:,:,2)  ! assume parabolic near wall
        L22(:,:) = v(:,:,1)*v(:,:,1)  ! uv-node
        L23(:,:) = v(:,:,jz)*0.25_rprec*w(:,:,2)  ! uv-node
        L33(:,:) = (0.25_rprec*w(:,:,2))**2  ! uv-node
        u_bar(:,:) = u(:,:,1)
        v_bar(:,:) = v(:,:,1)
        w_bar(:,:) = 0.25_rprec*w(:,:,2)
    else  ! w-nodes
        L11(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*                        &
                0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        L12(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*                        &
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        L13(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*w(:,:,jz)
        L22(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*                        &
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        L23(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*w(:,:,jz)
        L33(:,:) = w(:,:,jz)*w(:,:,jz)
        u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        w_bar(:,:) = w(:,:,jz)
    end if
    u_hat = u_bar
    v_hat = v_bar
    w_hat = w_bar

    ! Filter first term and add the second term to get the final value
    ! in-place filtering
    call test_filter ( u_bar )
    call test_filter ( v_bar )
    call test_filter ( w_bar )
    call test_filter ( L11 )
    L11 = L11 - u_bar*u_bar
    call test_filter ( L12 )
    L12 = L12 - u_bar*v_bar
    call test_filter ( L13 )
    L13 = L13 - u_bar*w_bar
    call test_filter ( L22 )
    L22 = L22 - v_bar*v_bar
    call test_filter ( L23 )
    L23 = L23 - v_bar*w_bar
    call test_filter ( L33 )
    L33 = L33 - w_bar*w_bar

    Q11 = u_bar*u_bar
    Q12 = u_bar*v_bar
    Q13 = u_bar*w_bar
    Q22 = v_bar*v_bar
    Q23 = v_bar*w_bar
    Q33 = w_bar*w_bar

    call test_test_filter ( u_hat )
    call test_test_filter ( v_hat )
    call test_test_filter ( w_hat )
    call test_test_filter ( Q11 )
    Q11 = Q11 - u_hat*u_hat
    call test_test_filter ( Q12 )
    Q12 = Q12 - u_hat*v_hat
    call test_test_filter ( Q13 )
    Q13 = Q13 - u_hat*w_hat
    call test_test_filter ( Q22 )
    Q22 = Q22 - v_hat*v_hat
    call test_test_filter ( Q23 )
    Q23 = Q23 - v_hat*w_hat
    call test_test_filter ( Q33 )
    Q33 = Q33 - w_hat*w_hat

    ! calculate |S|
    S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 +&
        S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 + &
        S13(:,:,jz)**2 + S23(:,:,jz)**2)))

    ! already on w-nodes
    S11_bar(:,:) = S11(:,:,jz)
    S12_bar(:,:) = S12(:,:,jz)
    S13_bar(:,:) = S13(:,:,jz)
    S22_bar(:,:) = S22(:,:,jz)
    S23_bar(:,:) = S23(:,:,jz)
    S33_bar(:,:) = S33(:,:,jz)

    S11_hat = S11_bar
    S12_hat = S12_bar
    S13_hat = S13_bar
    S22_hat = S22_bar
    S23_hat = S23_bar
    S33_hat = S33_bar

    call test_filter ( S11_bar )
    call test_filter ( S12_bar )
    call test_filter ( S13_bar )
    call test_filter ( S22_bar )
    call test_filter ( S23_bar )
    call test_filter ( S33_bar )

    call test_test_filter ( S11_hat )
    call test_test_filter ( S12_hat )
    call test_test_filter ( S13_hat )
    call test_test_filter ( S22_hat )
    call test_test_filter ( S23_hat )
    call test_test_filter ( S33_hat )

    S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

    S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 + &
        2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

    S_S11_bar = S*S11(:,:,jz)
    S_S12_bar = S*S12(:,:,jz)
    S_S13_bar = S*S13(:,:,jz)
    S_S22_bar = S*S22(:,:,jz)
    S_S23_bar = S*S23(:,:,jz)
    S_S33_bar = S*S33(:,:,jz)

    S_S11_hat = S_S11_bar
    S_S12_hat = S_S12_bar
    S_S13_hat = S_S13_bar
    S_S22_hat = S_S22_bar
    S_S23_hat = S_S23_bar
    S_S33_hat = S_S33_bar

    call test_filter ( S_S11_bar )
    call test_filter ( S_S12_bar )
    call test_filter ( S_S13_bar )
    call test_filter ( S_S22_bar )
    call test_filter ( S_S23_bar )
    call test_filter ( S_S33_bar )

    call test_test_filter ( S_S11_hat )
    call test_test_filter ( S_S12_hat )
    call test_test_filter ( S_S13_hat )
    call test_test_filter ( S_S22_hat )
    call test_test_filter ( S_S23_hat )
    call test_test_filter ( S_S33_hat )

    ! note: check that the Nyquist guys are zero!
    ! the 1./(nx*ny) is not really neccessary, but in practice it does
    ! affect the results.
    a1 = -2._rprec*(delta**2)*4._rprec*sum( S_bar*(S11_bar*L11 + S22_bar*L22 + &
        S33_bar*L33 + 2._rprec*( S12_bar*L12 + S13_bar*L13 + S23_bar*L23)))/   &
        (nx*ny)
    b1 = -2._rprec*(delta**2)*sum( S_S11_bar*L11 + S_S22_bar*L22 +             &
        S_S33_bar*L33 + 2._rprec*( S_S12_bar*L12 + S_S13_bar*L13 +             &
        S_S23_bar*L23))/(nx*ny)
    c1 = (2._rprec*delta**2)**2 * sum( S_S11_bar**2 + S_S22_bar**2 +           &
        S_S33_bar**2 + &
        2._rprec*( S_S12_bar**2 + S_S13_bar**2 + S_S23_bar**2))/(nx*ny)
    d1 = (2._rprec*delta**2)**2 * 16._rprec*sum( 0.5_rprec*S_bar**4 )/(nx*ny)
    e1 = 2._rprec*(2._rprec*delta**2)**2*4._rprec*sum( S_bar*( S11_bar*        &
        S_S11_bar +S22_bar*S_S22_bar + S33_bar*S_S33_bar + 2._rprec*(          &
        S12_bar*S_S12_bar + S13_bar*S_S13_bar + S23_bar*S_S23_bar)))/          &
        (nx*ny)

    a2 = -2._rprec*(delta**2)*16._rprec*sum( S_hat*(S11_hat*Q11 + S22_hat*Q22 +&
        S33_hat*Q33 + 2._rprec*( S12_hat*Q12 + S13_hat*Q13 + S23_hat*Q23)))/&
        (nx*ny)
    b2 = -2._rprec*(delta**2)*sum( S_S11_hat*Q11 + S_S22_hat*Q22 +&
        S_S33_hat*Q33 + 2._rprec*( S_S12_hat*Q12 + S_S13_hat*Q13 + &
        S_S23_hat*Q23))/(nx*ny)
    c2 = (2._rprec*delta**2)**2 * sum( S_S11_hat**2 + S_S22_hat**2 +&
        S_S33_hat**2 + &
        2._rprec*( S_S12_hat**2 + S_S13_hat**2 + S_S23_hat**2))/(nx*ny)
    d2 = (2._rprec*delta**2)**2 *256._rprec*sum( 0.5_rprec*S_hat**4 )/(nx*ny)
    e2 = 2._rprec*(2._rprec*delta**2)**2*16._rprec*sum( S_hat*( S11_hat*       &
        S_S11_hat + S22_hat*S_S22_hat + S33_hat*S_S33_hat + 2.*(               &
        S12_hat*S_S12_hat + S13_hat*S_S13_hat + S23_hat*S_S23_hat)))/          &
        (nx*ny)

    A(0) = b2*c1 - b1*c2
    A(1) = a1*c2 - b2*e1
    A(2) = b2*d1 + b1*e2 - a2*c1
    A(3) = a2*e1 - a1*e2
    A(4) = -a2*d1 - b1*d2
    A(5) = a1*d2

    beta(jz) = rtnewt(A,jz)

    ! now put beta back into M_ij: using Q_ij as storage
    const = 2._rprec*delta**2
    M11 = const*(S_S11_bar - 4._rprec*beta(jz)*S_bar*S11_bar)
    M12 = const*(S_S12_bar - 4._rprec*beta(jz)*S_bar*S12_bar)
    M13 = const*(S_S13_bar - 4._rprec*beta(jz)*S_bar*S13_bar)
    M22 = const*(S_S22_bar - 4._rprec*beta(jz)*S_bar*S22_bar)
    M23 = const*(S_S23_bar - 4._rprec*beta(jz)*S_bar*S23_bar)
    M33 = const*(S_S33_bar - 4._rprec*beta(jz)*S_bar*S33_bar)

    Cs_1D(jz) = sum(L11*M11 + L22*M22 + L33*M33 + 2._rprec*(L12*M12 +          &
        L13*M13 + L23*M23)) /                                                  &
        sum(M11**2 + M22**2 + M33**2 + 2._rprec*(M12**2 + M13**2 + M23**2))
    Cs_1D(jz) = max(0._rprec, real(Cs_1D(jz),rprec))

    ! Calculate ee_now (the current value of eij*eij)
    LM = L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
    MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)
    ee_now(:,:,jz) = L11**2+L22**2+L33**2+2._rprec*(L12**2+L13**2+L23**2)      &
                    -2._rprec*LM*Cs_1D(jz) + MM*Cs_1D(jz)**2
end do
endif

! Nullify pointers
nullify( M11, M12, M13, M22, M23, M33 )

end subroutine scaledep_dynamic
