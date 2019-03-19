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
subroutine convec(u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy,RHSx,RHSy,RHSz)
!*******************************************************************************
!
! Computes the rotation convective term in physical space
!       c = - (u X vort)
! Uses 3/2-rule for dealiasing for more info see Canuto 1991 Spectral Methods,
! chapter 7
!
use types, only : rprec
use param
use fft
use derivatives, only : convolve_rnl, dft_direct_forw_2d_n_yonlyC_big,      &
    dft_direct_back_2d_n_yonlyC_big
use sim_param, only : zhyb

implicit none

integer :: jz, jz_min
integer :: jzLo, jzHi, jz_max  ! added for full channel capabilities

real(rprec), save, allocatable, dimension(:,:,:) :: cc_big,                    &
    u_big, v_big, w_big, vort1_big, vort2_big, vort3_big
logical, save :: arrays_allocated = .false.

real(rprec) :: const

!! big temp variables for hybrid_fourier only
real(rprec), dimension(ld_big,ny2) :: u_temp, v_temp, w_temp, vort1_temp, vort2_temp

real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: u, v, w,                &
    dudy, dudz, dvdx, dvdz, dwdx, dwdy

real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: RHSx, RHSy, RHSz

if (sgs) then
    jzLo = 2        !! necessary for LES or else blows up ....?
    jzHi = nz-1     !! can remove after testing
else
    jzLo = 1        !! for DNS
    jzHi = nz-1     !! can remove after testing
endif

if( .not. arrays_allocated ) then
   allocate( cc_big( ld_big,ny2,nz ) )
   allocate( u_big(ld_big, ny2, lbz:nz) )
   allocate( v_big(ld_big, ny2, lbz:nz) )
   allocate( w_big(ld_big, ny2, lbz:nz) )
   allocate( vort1_big( ld_big,ny2,nz ) )
   allocate( vort2_big( ld_big,ny2,nz ) )
   allocate( vort3_big( ld_big,ny2,nz ) )
   arrays_allocated = .true.
endif

! Recall dudz, and dvdz are on UVP node for k=1 only
! So du2 does not vary from arg2a to arg2b in 1st plane (k=1)

! Loop through horizontal slices
! MPI: u_big, v_big needed at jz = 0, w_big not needed though
! MPI: could get u{1,2}_big
if (fourier) then
    const = 1._rprec
else !! keep const .neq. 1 for hybrid_fourier
    const = 1._rprec/(nx*ny)
endif

do jz = lbz, nz

    if (hybrid_fourier) then

    if (.not. zhyb(jz)) then
        ! use RHSx,RHSy,RHSz for temp storage
        RHSx(:,:,jz)=const*u(:,:,jz)
        RHSy(:,:,jz)=const*v(:,:,jz)
        RHSz(:,:,jz)=const*w(:,:,jz)

        ! do forward fft on normal-size arrays
        call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
    else
        ! use RHSx,RHSy,RHSz for temp storage
        RHSx(kxi,:,jz)=u(kxi,:,jz)
        RHSy(kxi,:,jz)=v(kxi,:,jz)
        RHSz(kxi,:,jz)=w(kxi,:,jz)

        ! no need to transform, already in fourier space
    endif

    ! zero pad: padd takes care of the oddballs
    ! no changes needed for fourier here
    call padd(u_big(:,:,jz), RHSx(:,:,jz))
    call padd(v_big(:,:,jz), RHSy(:,:,jz))
    call padd(w_big(:,:,jz), RHSz(:,:,jz))

    ! Keep w_big used for interface calculation of RHSx and RHSy in
    ! fourier space, store as a temp variable
    if (jz .ne. lbz) then
    if ( (zhyb(jz-1)) .and. (.not. zhyb(jz)) ) then
        w_temp(:,:) = w_big(:,:,jz)
    endif
    endif

    if (.not. zhyb(jz)) then
        ! Back to physical space
        call dfftw_execute_dft_c2r(back_big, u_big(:,:,jz), u_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, v_big(:,:,jz), v_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, w_big(:,:,jz), w_big(:,:,jz))
    endif

    ! Transform u_big and v_big to physical space at interface for 
    ! calculation of RHSz, store as a temp variable
    if (jz .ne. nz) then
    if ( (zhyb(jz)) .and. (.not. zhyb(jz+1)) ) then
        u_temp(:,:) = u_big(:,:,jz)
        v_temp(:,:) = v_big(:,:,jz)
        call dfftw_execute_dft_c2r(back_big, u_temp(:,:), u_temp(:,:))
        call dfftw_execute_dft_c2r(back_big, v_temp(:,:), v_temp(:,:))
    endif
    endif

    else !! not hybrid_fourier

    ! use RHSx,RHSy,RHSz for temp storage
    RHSx(:,:,jz)=const*u(:,:,jz)
    RHSy(:,:,jz)=const*v(:,:,jz)
    RHSz(:,:,jz)=const*w(:,:,jz)

    if (.not. fourier) then
        ! do forward fft on normal-size arrays
        call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
    endif
    ! else (fourier) no need to transform, already in fourier space

    ! zero pad: padd takes care of the oddballs
    ! no changes needed for fourier here
    call padd(u_big(:,:,jz), RHSx(:,:,jz))
    call padd(v_big(:,:,jz), RHSy(:,:,jz))
    call padd(w_big(:,:,jz), RHSz(:,:,jz))

    if (.not. fourier) then
        ! Back to physical space
        call dfftw_execute_dft_c2r(back_big, u_big(:,:,jz), u_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, v_big(:,:,jz), v_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, w_big(:,:,jz), w_big(:,:,jz))
    endif
    ! else (fourier) no need to transform, want to stay in fourier space

    endif

end do

! Do the same thing with the vorticity
do jz = 1, nz
    ! if dudz, dvdz are on u-nodes for jz=1, then we need a special
    ! definition of the vorticity in that case which also interpolates
    ! dwdx, dwdy to the u-node at jz=1
    if ( (coord == 0) .and. (jz == 1) ) then

        select case (lbc_mom)
        ! Stress free
        case (0)
            RHSx(:, :, 1) = 0._rprec
            RHSy(:, :, 1) = 0._rprec

        ! Wall (all cases >= 1)
        case (1:)
            if (hybrid_fourier) then
            !! remove const since it should be 1.0_rprec for fourier
            if (sgs) then
                ! dwdy(jz=1) should be 0, so we can use this
                RHSx(kxi, :, 1) = ( 0.5_rprec * (dwdy(kxi, :, 1) +             &
                    dwdy(kxi, :, 2))  - dvdz(kxi, :, 1) )
                ! dwdx(jz=1) should be 0, so we can use this
                RHSy(kxi, :, 1) = ( dudz(kxi, :, 1) -                          &
                    0.5_rprec * (dwdx(kxi, :, 1) + dwdx(kxi, :, 2)) )
            else ! for DNS, dudz(1) and dvdz(1) located at the w-node(1)
                RHSx(kxi, :, 1) = ( 0.5_rprec *(dwdy(kxi, :, 1)+ dwdy(kxi, :,2)) &
                    - 0.5_rprec *(dvdz(kxi, :, 1)+dvdz(kxi, :, 2)))
                ! RHSx(:, :, 1) located at uv-node(1)
                RHSy(kxi, :, 1) = (0.5_rprec *(dudz(kxi, :, 1)+ dudz(kxi, :, 2)) &
                    -0.5_rprec *(dwdx(kxi, :, 1)+ dwdx(kxi, :, 2)))
                ! RHSy(:, :, 1) located at uv-node(1)
            endif
            else !! fourier or not fourier
            if (sgs) then
                ! dwdy(jz=1) should be 0, so we can use this
                RHSx(:, :, 1) = const * ( 0.5_rprec * (dwdy(:, :, 1) +             &
                    dwdy(:, :, 2))  - dvdz(:, :, 1) )
                ! dwdx(jz=1) should be 0, so we can use this
                RHSy(:, :, 1) = const * ( dudz(:, :, 1) -                          &
                    0.5_rprec * (dwdx(:, :, 1) + dwdx(:, :, 2)) )
            else ! for DNS, dudz(1) and dvdz(1) located at the w-node(1)
                RHSx(:, :, 1) = const * ( 0.5_rprec *(dwdy(:, :, 1)+ dwdy(:, :,2)) &
                    - 0.5_rprec *(dvdz(:, :, 1)+dvdz(:, :, 2)))
                ! RHSx(:, :, 1) located at uv-node(1)
                RHSy(:, :, 1) = const * (0.5_rprec *(dudz(:, :, 1)+ dudz(:, :, 2)) &
                    -0.5_rprec *(dwdx(:, :, 1)+ dwdx(:, :, 2)))
                ! RHSy(:, :, 1) located at uv-node(1)
            endif
            endif
        end select
    endif

    if ( (coord == nproc-1) .and. (jz == nz) ) then

        select case (ubc_mom)

        ! Stress free
        case (0)
            RHSx(:, :, nz) = 0._rprec
            RHSy(:, :, nz) = 0._rprec

        ! No-slip and wall model
        case (1:)
            ! dwdy(jz=1) should be 0, so we could use this
            ! this RHSx = vort1 is actually uvp nz-1 but stored as w nz
            RHSx(:, :, nz) = const * ( 0.5_rprec * (dwdy(:, :, nz-1) +         &
                dwdy(:, :, nz)) - dvdz(:, :, nz-1) )
            ! dwdx(jz=1) should be 0, so we could use this
            ! this RHSy = vort2 is actually uvp nz-1 but stored as w nz
            RHSy(:, :, nz) = const * ( dudz(:, :, nz-1) -                      &
                0.5_rprec * (dwdx(:, :, nz-1) + dwdx(:, :, nz)) )

        end select
    endif

    ! very kludgy -- fix later      !! channel
    if (.not.(coord==0 .and. jz==1) .and. .not. (ubc_mom>0 .and.               &
        coord==nproc-1 .and. jz==nz)  ) then
        if (hybrid_fourier) then !! be careful of the const here!
            if (zhyb(jz)) then !! const = 1
                RHSx(kxi,:,jz)=(dwdy(kxi,:,jz)-dvdz(kxi,:,jz))
                RHSy(kxi,:,jz)=(dudz(kxi,:,jz)-dwdx(kxi,:,jz))
            else !! const = 1 / (nx*ny)
                RHSx(:,:,jz)=const*(dwdy(:,:,jz)-dvdz(:,:,jz))
                RHSy(:,:,jz)=const*(dudz(:,:,jz)-dwdx(:,:,jz))
            endif
        else !! not fourier and fourier
            RHSx(:,:,jz)=const*(dwdy(:,:,jz)-dvdz(:,:,jz))
            RHSy(:,:,jz)=const*(dudz(:,:,jz)-dwdx(:,:,jz))
        endif
    end if

    if (hybrid_fourier) then
        if (zhyb(jz)) then !! const = 1
            RHSz(kxi,:,jz)=(dvdx(kxi,:,jz)-dudy(kxi,:,jz))
        else !! const = 1 / (nx*ny)
            RHSz(:,:,jz)=const*(dvdx(:,:,jz)-dudy(:,:,jz))
        endif
    else !! fourier and not fourier
        RHSz(:,:,jz)=const*(dvdx(:,:,jz)-dudy(:,:,jz))
    endif

    if (hybrid_fourier) then

    if (.not. zhyb(jz)) then
        ! do forward fft on normal-size arrays
        call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
    endif 
    ! else (fourier) no need to transform, already in fourier space

    ! no changes needed for fourier here
    call padd(vort1_big(:,:,jz), RHSx(:,:,jz))
    call padd(vort2_big(:,:,jz), RHSy(:,:,jz))
    call padd(vort3_big(:,:,jz), RHSz(:,:,jz))

    ! Keep variables used for interface calculation of RHSx and RHSy in
    ! fourier space, store as a temp variable before transform
    if (jz .ne. 1) then
    if ( (zhyb(jz-1)) .and. (.not. zhyb(jz)) ) then
        vort1_temp(:,:) = vort1_big(:,:,jz)
        vort2_temp(:,:) = vort2_big(:,:,jz)
    endif
    endif

    if (.not. zhyb(jz)) then
        ! Back to physical space
        call dfftw_execute_dft_c2r(back_big, vort1_big(:,:,jz), vort1_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, vort2_big(:,:,jz), vort2_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, vort3_big(:,:,jz), vort3_big(:,:,jz))
    endif
    ! else (fourier) no need to transform, want to stay in fourier space

    else !! not hybrid_fourier

    if (.not. fourier) then
        ! do forward fft on normal-size arrays
        call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
        call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
    endif 
    ! else (fourier) no need to transform, already in fourier space

    ! no changes needed for fourier here
    call padd(vort1_big(:,:,jz), RHSx(:,:,jz))
    call padd(vort2_big(:,:,jz), RHSy(:,:,jz))
    call padd(vort3_big(:,:,jz), RHSz(:,:,jz))

    if (.not. fourier) then
        ! Back to physical space
        call dfftw_execute_dft_c2r(back_big, vort1_big(:,:,jz), vort1_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, vort2_big(:,:,jz), vort2_big(:,:,jz))
        call dfftw_execute_dft_c2r(back_big, vort3_big(:,:,jz), vort3_big(:,:,jz))
    endif
    ! else (fourier) no need to transform, want to stay in fourier space

    endif
end do

! For Fourier modes transform from (kx,ky) to (kx,y)
if (fourier) then
    ! Transform ky --> y
    do jz = lbz, nz
        call dft_direct_back_2d_n_yonlyC_big( u_big(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC_big( v_big(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC_big( w_big(:,:,jz) )
    enddo
    do jz = 1, nz
        call dft_direct_back_2d_n_yonlyC_big( vort1_big(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC_big( vort2_big(:,:,jz) )
        call dft_direct_back_2d_n_yonlyC_big( vort3_big(:,:,jz) )
    enddo
elseif (hybrid_fourier) then
    ! Transform ky --> y
    do jz = lbz, nz
        if (zhyb(jz)) then
            call dft_direct_back_2d_n_yonlyC_big( u_big(:,:,jz) )
            call dft_direct_back_2d_n_yonlyC_big( v_big(:,:,jz) )
            call dft_direct_back_2d_n_yonlyC_big( w_big(:,:,jz) )
        endif
    enddo
    do jz = 1, nz
        if (zhyb(jz)) then
            call dft_direct_back_2d_n_yonlyC_big( vort1_big(:,:,jz) )
            call dft_direct_back_2d_n_yonlyC_big( vort2_big(:,:,jz) )
            call dft_direct_back_2d_n_yonlyC_big( vort3_big(:,:,jz) )
        endif
        ! Variables used in convolution at interface
        if ( (zhyb(jz-1)) .and. (.not. zhyb(jz)) ) then
            call dft_direct_back_2d_n_yonlyC_big( w_temp(:,:) )
            call dft_direct_back_2d_n_yonlyC_big( vort1_temp(:,:) )
            call dft_direct_back_2d_n_yonlyC_big( vort2_temp(:,:) )
        endif
    enddo
endif

! RHSx
! redefinition of const
if (.not. fourier) then !! for not fourier and hybrid_fourier
    const=1._rprec/(nx2*ny2)
endif
! else (fourier) keep const = 1._rprec

! (fourier) u_big and vort_big in (kx,y,z) space, need to convolve because of kx
if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    if ((fourier) .or. (hybrid_fourier)) then
        cc_big(:,:,1)=(convolve_rnl(v_big(:,:,1),-vort3_big(:,:,1))    &
            +0.5_rprec*convolve_rnl(w_big(:,:,2),vort2_big(:,:,jzLo)))
    else !! not fourier or hybrid_fourier
        cc_big(:,:,1)=const*(v_big(:,:,1)*(-vort3_big(:,:,1))                &
            +0.5_rprec*w_big(:,:,2)*(vort2_big(:,:,jzLo)))   ! (default index was 2)
    endif
    !--vort2(jz=1) is located on uvp-node        ^  try with 1 (experimental)
    !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    !  above the wall (could arguably be 0.25 * w(:,:,2))
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1 ) then  ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    if (fourier) then
        cc_big(:,:,nz-1)=const*(convolve_rnl(v_big(:,:,nz-1),-vort3_big(:,:,nz-1))       &
            +0.5_rprec*convolve_rnl(w_big(:,:,nz-1),vort2_big(:,:,jzHi)))
    else !! hybrid_fourier or not fourier
        cc_big(:,:,nz-1)=const*(v_big(:,:,nz-1)*(-vort3_big(:,:,nz-1))       &
            +0.5_rprec*w_big(:,:,nz-1)*(vort2_big(:,:,jzHi)))
    endif
    !--vort2(jz=1) is located on uvp-node        ^  try with nz-1 (experimental)
    !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    !  below the wall (could arguably be 0.25 * w(:,:,2))

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max    !nz-1   ! channel
    if (fourier) then
        cc_big(:,:,jz)=const*(convolve_rnl(v_big(:,:,jz),-vort3_big(:,:,jz))       &
            +0.5_rprec*(convolve_rnl(w_big(:,:,jz+1),vort2_big(:,:,jz+1))          &
            +convolve_rnl(w_big(:,:,jz),vort2_big(:,:,jz))))
    elseif (hybrid_fourier) then
        if ( (zhyb(jz+1)) .and. (zhyb(jz)) ) then !! all points in Fourier
            cc_big(:,:,jz)=(convolve_rnl(v_big(:,:,jz),-vort3_big(:,:,jz))       &
                +0.5_rprec*(convolve_rnl(w_big(:,:,jz+1),vort2_big(:,:,jz+1))    &
                +convolve_rnl(w_big(:,:,jz),vort2_big(:,:,jz))))
        elseif ( (.not. zhyb(jz+1)) .and. (.not. zhyb(jz)) ) then !! all Physical
            cc_big(:,:,jz)=const*(v_big(:,:,jz)*(-vort3_big(:,:,jz))       &
                +0.5_rprec*(w_big(:,:,jz+1)*(vort2_big(:,:,jz+1))          &
                +w_big(:,:,jz)*(vort2_big(:,:,jz))))
        else !! jz+1 w node in Physical, jz w node in Fourier, Convolve at interface in Fourier
            cc_big(:,:,jz)=(convolve_rnl(v_big(:,:,jz),-vort3_big(:,:,jz))       &
                +0.5_rprec*(convolve_rnl(w_temp(:,:),vort2_temp(:,:))            &
                +convolve_rnl(w_big(:,:,jz),vort2_big(:,:,jz))))
        endif
    else !! not hybrid_fourier or fourier
        cc_big(:,:,jz)=const*(v_big(:,:,jz)*(-vort3_big(:,:,jz))       &
            +0.5_rprec*(w_big(:,:,jz+1)*(vort2_big(:,:,jz+1))          &
            +w_big(:,:,jz)*(vort2_big(:,:,jz))))
    endif
end do

! Loop through horizontal slices
do jz=1,nz-1
    if (fourier) then  !! transform y --> ky
         call dft_direct_forw_2d_n_yonlyC_big( cc_big(:,:,jz) )
    elseif (hybrid_fourier) then
        if (zhyb(jz)) then
            call dft_direct_forw_2d_n_yonlyC_big( cc_big(:,:,jz) )
        else 
            call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz),cc_big(:,:,jz))
        endif
    else !! not fourier or hybrid_fourier
        call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz),cc_big(:,:,jz))
    endif

    ! un-zero pad
    ! note: cc_big is going into RHSx
    ! (fourier) no changes needed here
    call unpadd(RHSx(:,:,jz),cc_big(:,:,jz))

    ! Back to physical space
    if (hybrid_fourier) then
        if (.not. zhyb(jz)) then
            call dfftw_execute_dft_c2r(back, RHSx(:,:,jz), RHSx(:,:,jz))
        endif
    elseif (.not. fourier) then !! not fourier or hybrid_fourier
        call dfftw_execute_dft_c2r(back, RHSx(:,:,jz), RHSx(:,:,jz))
    endif
    ! else (fourier) do nothing, exit with RHSx in fourier space
end do

! RHSy
! const should be 1./(nx2*ny2) here -> for fourier const = 1.0_rprec
if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    if ((fourier) .or. (hybrid_fourier)) then
        cc_big(:,:,1)=(convolve_rnl(u_big(:,:,1),vort3_big(:,:,1))   &
            +0.5_rprec*convolve_rnl(w_big(:,:,2),-vort1_big(:,:,jzLo)))
    else !! not fourier or hybrid_fourier
        cc_big(:,:,1)=const*(u_big(:,:,1)*(vort3_big(:,:,1))               &
            +0.5_rprec*w_big(:,:,2)*(-vort1_big(:,:,jzLo)))
    endif
    !--vort1(jz=1) is uvp-node                    ^ try with 1 (experimental)
    !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    !  above the wall (could arguably be 0.25 * w(:,:,2))
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then   ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    if (fourier) then
        cc_big(:,:,nz-1)=const*(convolve_rnl(u_big(:,:,nz-1),vort3_big(:,:,nz-1)) &
            +0.5_rprec*convolve_rnl(w_big(:,:,nz-1),-vort1_big(:,:,jzHi)))
    else !! hybrid_fourier or not fourier
        cc_big(:,:,nz-1)=const*(u_big(:,:,nz-1)*(vort3_big(:,:,nz-1))             &
            +0.5_rprec*w_big(:,:,nz-1)*(-vort1_big(:,:,jzHi)))
    endif
    !--vort1(jz=1) is uvp-node                    ^ try with nz-1 (experimental)
    !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    !  below the wall

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max  !nz - 1   ! channel
    if (fourier) then
        cc_big(:,:,jz)=(convolve_rnl(u_big(:,:,jz),vort3_big(:,:,jz))  &
            +0.5_rprec*(convolve_rnl(w_big(:,:,jz+1),-vort1_big(:,:,jz+1))   &
            +convolve_rnl(w_big(:,:,jz),-vort1_big(:,:,jz))))
    elseif (hybrid_fourier) then
        if ( (zhyb(jz+1)) .and. (zhyb(jz)) ) then !! all points in Fourier
            cc_big(:,:,jz)=(convolve_rnl(u_big(:,:,jz),vort3_big(:,:,jz))  &
                +0.5_rprec*(convolve_rnl(w_big(:,:,jz+1),-vort1_big(:,:,jz+1))   &
                +convolve_rnl(w_big(:,:,jz),-vort1_big(:,:,jz))))
        elseif ( (.not. zhyb(jz+1)) .and. (.not. zhyb(jz)) ) then !! all Physical
            cc_big(:,:,jz)=const*(u_big(:,:,jz)*(vort3_big(:,:,jz))              &
                +0.5_rprec*(w_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))               &
                +w_big(:,:,jz)*(-vort1_big(:,:,jz))))
        else !! jz+1 w node in Physical, jz w node in Fourier, Convolve at interface in Fourier
            cc_big(:,:,jz)=(convolve_rnl(u_big(:,:,jz),vort3_big(:,:,jz))  &
                +0.5_rprec*(convolve_rnl(w_temp(:,:),-vort1_temp(:,:))     &
                +convolve_rnl(w_big(:,:,jz),-vort1_big(:,:,jz))))
        endif
    else !! not hybrid_fourier or fourier
        cc_big(:,:,jz)=const*(u_big(:,:,jz)*(vort3_big(:,:,jz))              &
            +0.5_rprec*(w_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))               &
            +w_big(:,:,jz)*(-vort1_big(:,:,jz))))
    endif
end do

do jz=1,nz-1
    if (fourier) then !! transform y --> ky
        call dft_direct_forw_2d_n_yonlyC_big( cc_big(:,:,jz) )
    elseif (hybrid_fourier) then
        if (zhyb(jz)) then
            call dft_direct_forw_2d_n_yonlyC_big( cc_big(:,:,jz) )
        else
            call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz), cc_big(:,:,jz))
        endif
    else !! no fourier or hybrid_fourier
        call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz), cc_big(:,:,jz))
    endif

    ! un-zero pad
    ! note: cc_big is going into RHSy
    ! (fourier) no changes needed here
    call unpadd(RHSy(:,:,jz), cc_big(:,:,jz))

    ! Back to physical space
    if (hybrid_fourier) then
        if (.not. zhyb(jz)) then
            call dfftw_execute_dft_c2r(back, RHSy(:,:,jz), RHSy(:,:,jz))
        endif
    elseif (.not. fourier) then !! not fourier or hybrid_fourier
        call dfftw_execute_dft_c2r(back, RHSy(:,:,jz), RHSy(:,:,jz))
    endif
    ! else (fourier) do nothing, exit with RHSy in fourier space
end do

! RHSz

if (coord == 0) then
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--wrong, we do solve z-eqn (w-eqn) at bottom wall --pj
    !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
    ! slip BC
    cc_big(:,:,1)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then     ! channel
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--but now we do solve z-eqn (w-eqn) at top wall --pj
    !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
    ! slip BC
    cc_big(:,:,nz)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_max = nz-1
else
    jz_max = nz-1   !! or nz ?       ! channel
end if

!#ifdef PPMPI
!  if (coord == nproc-1) then
!    cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!    jz_max = nz - 1
!  else
!    jz_max = nz
!  endif
!#else
!  cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!  jz_max = nz - 1
!#endif

! channel
do jz = jz_min, jz_max    !nz - 1
    if (fourier) then
        cc_big(:,:,jz) = const*0.5_rprec*(                                         &
            convolve_rnl(u_big(:,:,jz)+u_big(:,:,jz-1),-vort2_big(:,:,jz))         &
            +convolve_rnl(v_big(:,:,jz)+v_big(:,:,jz-1),vort1_big(:,:,jz)))
    elseif (hybrid_fourier) then
        if ( (zhyb(jz)) .and. (zhyb(jz-1)) ) then !! all points in Fourier
            cc_big(:,:,jz) = 0.5_rprec*(                                         &
                convolve_rnl(u_big(:,:,jz)+u_big(:,:,jz-1),-vort2_big(:,:,jz))   &
                +convolve_rnl(v_big(:,:,jz)+v_big(:,:,jz-1),vort1_big(:,:,jz)))
        elseif ( (.not. zhyb(jz)) .and. (.not. zhyb(jz-1)) ) then !! all Physical
            cc_big(:,:,jz) = const*0.5_rprec*(                                   &
                (u_big(:,:,jz)+u_big(:,:,jz-1))*(-vort2_big(:,:,jz))             &
                +(v_big(:,:,jz)+v_big(:,:,jz-1))*(vort1_big(:,:,jz)))
        else !! jz uv node in Physical, jz-1 uv node at interface in Fourier, Product in Physical
            cc_big(:,:,jz) = const*0.5_rprec*(                                   &
                (u_big(:,:,jz)+u_temp(:,:))*(-vort2_big(:,:,jz))             &
                +(v_big(:,:,jz)+v_temp(:,:))*(vort1_big(:,:,jz)))
        endif
    else !! not fourier or hybrid_fourier
        cc_big(:,:,jz) = const*0.5_rprec*(                                         &
            (u_big(:,:,jz)+u_big(:,:,jz-1))*(-vort2_big(:,:,jz))                   &
            +(v_big(:,:,jz)+v_big(:,:,jz-1))*(vort1_big(:,:,jz)))
    endif
end do

! Loop through horizontal slices
do jz=1,nz !nz - 1
    if (fourier) then
        call dft_direct_forw_2d_n_yonlyC_big( cc_big(:,:,jz) )
    elseif (hybrid_fourier) then
        if (zhyb(jz)) then
            call dft_direct_forw_2d_n_yonlyC_big( cc_big(:,:,jz) )
        else
            call dfftw_execute_dft_r2c(forw_big,cc_big(:,:,jz),cc_big(:,:,jz))
        endif
    else !! not fourier or hybrid_fourier
        call dfftw_execute_dft_r2c(forw_big,cc_big(:,:,jz),cc_big(:,:,jz))
    endif

    ! un-zero pad
    ! note: cc_big is going into RHSz!!!!
    ! (fourier) no changes needed here
    call unpadd(RHSz(:,:,jz),cc_big(:,:,jz))

    ! Back to physical space
    if (hybrid_fourier) then
        if (.not. zhyb(jz)) then
            call dfftw_execute_dft_c2r(back,RHSz(:,:,jz),   RHSz(:,:,jz))
        endif
    elseif (.not. fourier) then !! not fourier or hybrid_fourier
        call dfftw_execute_dft_c2r(back,RHSz(:,:,jz),   RHSz(:,:,jz))
    endif
    ! else (fourier) do nothing, exit with RHSz in fourier space
end do

#ifdef PPMPI
#ifdef PPSAFETYMODE
RHSx(:, :, 0) = BOGUS
RHSy(:, :, 0) = BOGUS
RHSz(: ,:, 0) = BOGUS
#endif
#endif

!--top level is not valid
#ifdef PPSAFETYMODE
RHSx(:, :, nz) = BOGUS
RHSy(:, :, nz) = BOGUS
if(coord<nproc-1) RHSz(:, :, nz) = BOGUS
#endif

end subroutine convec
