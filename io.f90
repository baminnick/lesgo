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

!*****************************************************************************
module io
!*****************************************************************************
use types, only : rprec
use param, only : ld, nx, ny, nz, nz_tot, path, coord, rank, nproc, jt_total
use param, only : total_time, total_time_dim, lbz, jzmin, jzmax
use param, only : cumulative_time, fcumulative_time
use sim_param, only : w, dudz, dvdz
use sgs_param, only : cs_opt2
use string_util
use messages

#ifdef PPMPI
use mpi
#endif

#ifdef PPCGNS
use cgns
#ifdef PPMPI
use param, only: ierr
#endif
#endif

implicit none
save
private

public jt_total, openfiles, energy, output_loop, output_final, output_init, &
    write_tau_wall_bot, write_tau_wall_top, kx_energy, kx_energy_fourier,   &
    ky_energy
#ifdef PPOUTPUT_CLOCK
public write_clocks
#endif

! Where to end with nz index.
integer :: nz_end

contains

!*****************************************************************************
subroutine openfiles()
!*****************************************************************************
use param, only : use_cfl_dt, dt, cfl_f
implicit none
logical :: exst

! Temporary values used to read time step and CFL from file
real(rprec) :: dt_r, cfl_r

if (cumulative_time) then
    inquire (file=fcumulative_time, exist=exst)
    if (exst) then
        open (1, file=fcumulative_time)
        read(1, *) jt_total, total_time, total_time_dim, dt_r, cfl_r
        close (1)
    else
        ! assume this is the first run on cumulative time
        if ( coord == 0 ) then
            write (*, *) '--> Assuming jt_total = 0, total_time = 0.0'
        end if
        jt_total = 0
        total_time = 0._rprec
        total_time_dim = 0._rprec
        ! Initialize temp values to be overwritten
        dt_r = 0.0_rprec
        cfl_r = 0.0_rprec
    end if
end if

! Update dynamic time stepping info if required; otherwise discard.
if ( use_cfl_dt ) then
    dt = dt_r
    cfl_f = cfl_r
end if

end subroutine openfiles

!*****************************************************************************
subroutine energy (ke)
!*****************************************************************************
use types, only : rprec
use param
use sim_param, only : u, v, w
use messages
implicit none
integer :: jx, jy, jz, nan_count
real(rprec)::KE,temp_w
#ifdef PPMPI
real(rprec) :: ke_global
#endif

! Initialize variables
nan_count = 0
ke = 0._rprec

do jz = 1, nz-1
do jy = 1, ny
do jx = 1, nx
    temp_w = 0.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
    ke = ke + (u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)
end do
end do
end do

! Perform spatial averaging
ke = ke*0.5_rprec/(nx*ny*(nz-1))

#ifdef PPMPI
call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then  ! note that it's rank here, not coord
    ke = ke_global/nproc
#endif
    open(2,file=path // 'output/check_ke.dat', status='unknown',               &
        form='formatted', position='append')
    write(2,*) total_time,ke
    close(2)
#ifdef PPMPI
end if
#endif

end subroutine energy

!*****************************************************************************
subroutine xpert (u,upert)
!*****************************************************************************
! 
! This routine takes in field variable u and outputs the streamwise fluctuating 
! component upert
! 
! This function is intended for .not. fourier
! 
use types, only: rprec
use param, only: ld, nx, ny, lbz, nz

implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: u
real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: upert
real(rprec), dimension(ny,lbz:nz) :: uavg
integer :: jx

! Initialize the temporary average variable
uavg = 0.0_rprec

! Compute streamwise average
do jx = 1, nx
    uavg(:,:) = uavg(:,:) + u(jx,:,:)
enddo
uavg = uavg / nx

! Compute streamwise fluctuating component
do jx = 1, nx
    upert(jx,:,:) = u(jx,:,:) - uavg(:,:)
enddo

return
end subroutine xpert

!*****************************************************************************
subroutine ypert (u,upert)
!*****************************************************************************
! 
! This routine takes in field variable u and outputs the spanwise fluctuating 
! component upert
! 
! This function is intended for .not. fourier
! 
use types, only: rprec
use param, only: ld, ny, lbz, nz

implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: u
real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: upert
real(rprec), dimension(ld,lbz:nz) :: uavg
integer :: jy

! Initialize the temporary average variable
uavg = 0.0_rprec

! Compute spanwise average
do jy = 1, ny
    uavg(:,:) = uavg(:,:) + u(:,jy,:)
enddo
uavg = uavg / ny

! Compute spanwise fluctuating component
do jy = 1, ny
    upert(:,jy,:) = u(:,jy,:) - uavg(:,:)
enddo

return
end subroutine ypert

!*****************************************************************************
subroutine ypert_fourier (u,upert)
!*****************************************************************************
! 
! This routine takes in field variable u and outputs the spanwise fluctuating 
! component upert
! 
! This function is intended for fourier
! 
use types, only: rprec
use param, only: ny, lbz, nz, nxp

implicit none

real(rprec), dimension(nxp+2,ny,lbz:nz), intent(in) :: u
real(rprec), dimension(nxp+2,ny,lbz:nz), intent(out) :: upert
real(rprec), dimension(nxp+2,lbz:nz) :: uavg
integer :: jy

! Initialize the temporary average variable
uavg = 0.0_rprec

! Compute spanwise average
do jy = 1, ny
    uavg(:,:) = uavg(:,:) + u(:,jy,:)
enddo
uavg = uavg / ny

! Compute spanwise fluctuating component
do jy = 1, ny
    upert(:,jy,:) = u(:,jy,:) - uavg(:,:)
enddo

return
end subroutine ypert_fourier

!*****************************************************************************
subroutine zpert (u,upert)
!*****************************************************************************
! 
! This routine takes in field variable u and outputs the wall-normal fluctuating 
! component upert
! 
use types, only: rprec
use param, only: ld, ny, lbz, nz

implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: u
real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: upert
real(rprec), dimension(ld,ny) :: uavg
integer :: jz

! Initialize the temporary average variable
uavg = 0.0_rprec

! Compute wall-normal average
do jz = lbz, nz
    uavg(:,:) = uavg(:,:) + u(:,:,jz)
enddo
uavg = uavg / (nz - lbz)

! Compute wall-normal fluctuating component
do jz = lbz, nz
    upert(:,:,jz) = u(:,:,jz) - uavg(:,:)
enddo

return
end subroutine zpert

!*****************************************************************************
subroutine ypert_by_z (u,upert)
!*****************************************************************************
! 
! This routine takes in field variable u and outputs the spanwise fluctuating 
! component upert. This subroutine differs from ypert by assuming the input
! variable u is not the entire flow domain, but only a particular z-level
! 
! This function is intended for .not. fourier
! 
use types, only: rprec
use param, only: ld, ny

implicit none

real(rprec), dimension(ld,ny), intent(in) :: u
real(rprec), dimension(ld,ny), intent(out) :: upert
real(rprec), dimension(ld) :: uavg
integer :: jy

! Initialize the temporary average variable
uavg = 0.0_rprec

! Compute spanwise average
do jy = 1, ny
    uavg(:) = uavg(:) + u(:,jy)
enddo
uavg = uavg / ny

! Compute spanwise fluctuating component
do jy = 1, ny
    upert(:,jy) = u(:,jy) - uavg(:)
enddo

return
end subroutine ypert_by_z

!*****************************************************************************
subroutine kx_energy ()
!*****************************************************************************
! 
! Computes and writes the energy of the streamwise modes
! 
! This function is intended for .not. fourier
! 
use types, only: rprec
use param
use sim_param, only: u, v, w
use messages
use fft
use functions, only: int2str
#ifdef PPSCALARS
use scalars, only: theta
#endif

implicit none

integer :: jx, jy, jz
real(rprec), dimension(ld,ny,lbz:nz) :: upert, vpert, wpert
complex(rprec), dimension(nx/2+1) :: uhat, vhat, what
real(rprec), dimension(nx/2+1) :: ke, uu, vv, ww
#ifdef PPOUTPUT_SSP
complex(rprec), dimension(nx/2+1) :: bigvhat, bigwhat
real(rprec), dimension(nx/2+1) :: roll, bigvv, bigww
#endif
#ifdef PPSCALARS
real(rprec), dimension(ld,ny,lbz:nz) :: theta_pert
complex(rprec), dimension(nx/2+1) :: theta_hat
real(rprec), dimension(nx/2+1) :: theta2
#endif
character(len = 10), dimension(nx/2+1) :: fhead !! file header variable

#ifdef PPMPI
real(rprec), dimension(nx/2+1) :: ke_global
#ifdef PPOUTPUT_SSP
real(rprec), dimension(nx/2+1) :: uu_global, roll_global
#endif
#ifdef PPSCALARS
real(rprec), dimension(nx/2+1) :: theta2_global
#endif
#endif

! Initialize variables
ke = 0.0_rprec
uu = 0.0_rprec
vv = 0.0_rprec
ww = 0.0_rprec
#ifdef PPOUTPUT_SSP
bigvv = 0.0_rprec
bigww = 0.0_rprec
#endif
#ifdef PPSCALARS
theta2 = 0.0_rprec
#endif

#ifdef PPMPI
ke_global = 0.0_rprec
#ifdef PPOUTPUT_SSP
uu_global = 0.0_rprec
roll_global = 0.0_rprec
#endif
#ifdef PPSCALARS
theta2_global = 0.0_rprec
#endif
#endif

! Initialize file, write header
if (jt_total == wbase) then
    do jx = 1, nx/2 + 1
        fhead(jx) = 'kx='//trim(int2str(jx-1))//' '
    enddo
endif

! Take out spanwise mean, ky = 0
call ypert(u,upert)
call ypert(v,vpert)
call ypert(w,wpert)
#ifdef PPSCALARS
call ypert(theta,theta_pert)
#endif

! Consider each y and z location
do jy = 1, ny
do jz = 1, nz-1

    ! Take 1D Fourier transform
    call dfftw_execute_dft_r2c( forw_x, upert(:,jy,jz), uhat)
    call dfftw_execute_dft_r2c( forw_x, vpert(:,jy,jz), vhat)
    call dfftw_execute_dft_r2c( forw_x, wpert(:,jy,jz), what)
#ifdef PPOUTPUT_SSP
    call dfftw_execute_dft_r2c( forw_x, v(:,jy,jz), bigvhat)
    call dfftw_execute_dft_r2c( forw_x, w(:,jy,jz), bigwhat)
#endif
#ifdef PPSCALARS
    call dfftw_execute_dft_r2c( forw_x, theta_pert(:,jy,jz), theta_hat)
#endif

    ! Normalize transformed variables
    uhat = uhat / nx
    vhat = vhat / nx
    what = what / nx
#ifdef PPOUTPUT_SSP
    bigvhat = bigvhat / nx
    bigwhat = bigwhat / nx
#endif
#ifdef PPSCALARS
    theta_hat = theta_hat / nx
#endif

    ! Sum over boundary points
    uu(1) = uu(1) + 0.5d0*real(uhat(1)*conjg(uhat(1)))
    vv(1) = vv(1) + 0.5d0*real(vhat(1)*conjg(vhat(1)))
    ww(1) = ww(1) + 0.5d0*real(what(1)*conjg(what(1)))
#ifdef PPOUTPUT_SSP
    bigvv(1) = bigvv(1) + 0.5d0*real(bigvhat(1)*conjg(bigvhat(1)))
    bigww(1) = bigww(1) + 0.5d0*real(bigwhat(1)*conjg(bigwhat(1)))
#endif
#ifdef PPSCALARS
    theta2(1) = theta2(1) + 0.5d0*real(theta_hat(1)*conjg(theta_hat(1)))
#endif

    uu(lh) = uu(lh) + 0.5d0*real(uhat(lh)*conjg(uhat(lh)))
    vv(lh) = vv(lh) + 0.5d0*real(vhat(lh)*conjg(vhat(lh)))
    ww(lh) = ww(lh) + 0.5d0*real(what(lh)*conjg(what(lh)))
#ifdef PPOUTPUT_SSP
    bigvv(lh) = bigvv(lh) + 0.5d0*real(bigvhat(lh)*conjg(bigvhat(lh)))
    bigww(lh) = bigww(lh) + 0.5d0*real(bigwhat(lh)*conjg(bigwhat(lh)))
#endif
#ifdef PPSCALARS
    theta2(lh) = theta2(lh) + 0.5d0*real(theta_hat(lh)*conjg(theta_hat(lh)))
#endif

    ! Sum over interior points
    do jx = 2, lh-1
        uu(jx) = uu(jx) + real(uhat(jx)*conjg(uhat(jx)))
        vv(jx) = vv(jx) + real(vhat(jx)*conjg(vhat(jx)))
        ww(jx) = ww(jx) + real(what(jx)*conjg(what(jx)))
#ifdef PPOUTPUT_SSP
        bigvv(jx) = bigvv(jx) + 0.5d0*real(bigvhat(jx)*conjg(bigvhat(jx)))
        bigww(jx) = bigww(jx) + 0.5d0*real(bigwhat(jx)*conjg(bigwhat(jx)))
#endif
#ifdef PPSCALARS
        theta2(jx) = theta2(jx) + 0.5d0*real(theta_hat(jx)*conjg(theta_hat(jx)))
#endif
    enddo

enddo
enddo

! Normalize by spanwise length - treating as a spanwise average
uu = uu / L_y
vv = vv / L_y
ww = ww / L_y
#ifdef PPOUTPUT_SSP
bigvv = bigvv / L_y
bigww = bigww / L_y
#endif
#ifdef PPSCALARS
theta2 = theta2 / L_y
#endif

! Make summations
ke = uu + vv + ww
#ifdef PPOUTPUT_SSP
roll = bigvv + bigww
#endif

! Write data to file
#ifdef PPMPI
call mpi_reduce (ke, ke_global, lh, MPI_RPREC, MPI_SUM, 0, comm, ierr)
#ifdef PPOUTPUT_SSP
call mpi_reduce (uu, uu_global, lh, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (roll, roll_global, lh, MPI_RPREC, MPI_SUM, 0, comm, ierr)
#endif
#ifdef PPSCALARS
call mpi_reduce ( theta2, theta2_global, lh, MPI_RPREC, MPI_SUM, 0, comm, ierr)
#endif
if (rank == 0) then  ! note that it's rank here, not coord
    ke = ke_global
#ifdef PPOUTPUT_SSP
    uu = uu_global
    roll = roll_global
#endif
#ifdef PPSCALARS
    theta2 = theta2_global
#endif
#endif
    ! kinetic energy of spanwise perturbations
    open(2,file=path // 'output/ke_kx.dat', status='unknown',               &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, ke
    close(2)

#ifdef PPOUTPUT_SSP
    ! Streamwise energy of spanwise perturbations - streak energy
    open(2,file=path // 'output/streak_kx.dat', status='unknown',           &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, uu
    close(2)

    ! Spanwise and Wall-normal energy - roll energy
    open(2,file=path // 'output/roll_kx.dat', status='unknown',             &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, roll
    close(2)
#endif

#ifdef PPSCALARS
    ! temperature "energy" of spanwise perturbations
    open(2,file=path // 'output/theta_kx.dat', status='unknown',            &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, theta2
    close(2)
#endif

#ifdef PPMPI
end if
#endif

end subroutine kx_energy

!*****************************************************************************
subroutine ky_energy ()
!*****************************************************************************
! 
! Computes and writes the energy of the spanwise modes
! 
! This function is intended for .not. fourier
! 
use types, only: rprec
use param
use sim_param, only: u, v, w
use messages
use fft
use functions, only: int2str

implicit none

integer :: jx, jy, jz
real(rprec), dimension(ld,ny,lbz:nz) :: upert, vpert, wpert
real(rprec), dimension(ny) :: utemp, vtemp, wtemp
complex(rprec), dimension(ny/2+1) :: uhat, vhat, what
real(rprec), dimension(ny/2+1) :: ke, uu, vv, ww
character(len = 10), dimension(ny/2+1) :: fhead !! file header variable

#ifdef PPMPI
    real(rprec), dimension(ny/2+1) :: ke_global
#endif

! Initialize variables
ke = 0.0_rprec
uu = 0.0_rprec
vv = 0.0_rprec
ww = 0.0_rprec

#ifdef PPMPI
    ke_global = 0.0_rprec
#endif

! Initialize file, write header
if (jt_total == wbase) then
    do jy = 1, ny/2 + 1
        fhead(jy) = 'ky='//trim(int2str(jy-1))//' '
    enddo
endif

! Take out streamwise mean, kx = 0
call xpert(u,upert)
call xpert(v,vpert)
call xpert(w,wpert)

! Consider each y and z location
do jx = 1, nx
do jz = 1, nz-1

    ! Store perturbations as temp variable to take fft
    ! This is only so a temp array is not allocated
    utemp = upert(jx,:,jz)
    vtemp = vpert(jx,:,jz)
    wtemp = wpert(jx,:,jz)

    ! Take 1D Fourier transform
    call dfftw_execute_dft_r2c( forw_y, utemp, uhat)
    call dfftw_execute_dft_r2c( forw_y, vtemp, vhat)
    call dfftw_execute_dft_r2c( forw_y, wtemp, what)

    ! Normalize transformed variables
    uhat = uhat / ny
    vhat = vhat / ny
    what = what / ny

    ! Sum over boundary points
    uu(1) = uu(1) + 0.5d0*real(uhat(1)*conjg(uhat(1)))
    vv(1) = vv(1) + 0.5d0*real(vhat(1)*conjg(vhat(1)))
    ww(1) = ww(1) + 0.5d0*real(what(1)*conjg(what(1)))

    uu(ny/2+1) = uu(ny/2+1) + 0.5d0*real(uhat(ny/2+1)*conjg(uhat(ny/2+1)))
    vv(ny/2+1) = vv(ny/2+1) + 0.5d0*real(vhat(ny/2+1)*conjg(vhat(ny/2+1)))
    ww(ny/2+1) = ww(ny/2+1) + 0.5d0*real(what(ny/2+1)*conjg(what(ny/2+1)))

    ! Sum over interior points
    do jy = 2, (ny/2+1-1)
        uu(jy) = uu(jy) + real(uhat(jy)*conjg(uhat(jy)))
        vv(jy) = vv(jy) + real(vhat(jy)*conjg(vhat(jy)))
        ww(jy) = ww(jy) + real(what(jy)*conjg(what(jy)))
    enddo

enddo
enddo

! Normalize by streamwise length - treating as a streamwise average
uu = uu / L_x
vv = vv / L_x
ww = ww / L_x

! Make summations
ke = uu + vv + ww

! Write data to file
#ifdef PPMPI
call mpi_reduce (ke, ke_global, lh, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then  ! note that it's rank here, not coord
    ke = ke_global
#endif
    ! kinetic energy in ky space
    open(2,file=path // 'output/ke_ky.dat', status='unknown',               &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, ke
    close(2)
#ifdef PPMPI
end if
#endif

end subroutine ky_energy

!*****************************************************************************
subroutine kx_energy_fourier ()
!*****************************************************************************
! 
! Computes and writes the energy of the streamwise modes
! 
! This function is intended for fourier
! 
use types, only: rprec
use param
use sim_param, only: uF, vF, wF
use messages
use fft
use functions, only: int2str
#ifdef PPSCALARS
use scalars, only: thetaF
#endif

implicit none

integer :: jx, jy, jz
real(rprec), dimension(nxp+2,ny,lbz:nz) :: upert, vpert, wpert
complex(rprec), dimension(nxp/2+1) :: uhat, vhat, what
real(rprec), dimension(nxp/2+1) :: ke, uu, vv, ww
character(len = 10), dimension(nxp/2+1) :: fhead !! file header variable
#ifdef PPOUTPUT_SSP
complex(rprec), dimension(nxp/2+1) :: bigvhat, bigwhat
real(rprec), dimension(nxp/2+1) :: roll, bigvv, bigww
#endif
#ifdef PPSCALARS
real(rprec), dimension(nxp+2,ny,lbz:nz) :: theta_pert
complex(rprec), dimension(nxp/2+1) :: theta_hat
real(rprec), dimension(nxp/2+1) :: theta2
#endif

#ifdef PPMPI
    real(rprec), dimension(nxp/2+1) :: ke_global
#ifdef PPOUTPUT_SSP
    real(rprec), dimension(nxp/2+1) :: uu_global, roll_global
#endif
#ifdef PPSCALARS
    real(rprec), dimension(nxp/2+1) :: theta2_global
#endif
#endif

! Initialize variables
ke = 0.0_rprec
uu = 0.0_rprec
vv = 0.0_rprec
ww = 0.0_rprec
#ifdef PPOUTPUT_SSP
bigvv = 0.0_rprec
bigww = 0.0_rprec
#endif

#ifdef PPMPI
    ke_global = 0.0_rprec
#ifdef PPOUTPUT_SSP
    uu_global = 0.0_rprec
    roll_global = 0.0_rprec
#endif
#ifdef PPSCALARS
theta2_global = 0.0_rprec
#endif
#endif

! Initialize file, write header
if (jt_total == wbase) then
    do jx = 1, nxp/2 + 1
        fhead(jx) = 'kx='//trim(int2str(jx-1))//' '
    enddo
endif

! Take out spanwise mean, ky = 0
call ypert_fourier(uF,upert)
call ypert_fourier(vF,vpert)
call ypert_fourier(wF,wpert)
#ifdef PPSCALARS
call ypert_fourier(thetaF,theta_pert)
#endif

! Consider each y and z location
do jy = 1, ny
do jz = 1, nz-1

    ! Take 1D Fourier transform
    call dfftw_execute_dft_r2c( forw_x_fourier, upert(:,jy,jz), uhat)
    call dfftw_execute_dft_r2c( forw_x_fourier, vpert(:,jy,jz), vhat)
    call dfftw_execute_dft_r2c( forw_x_fourier, wpert(:,jy,jz), what)
#ifdef PPOUTPUT_SSP
    call dfftw_execute_dft_r2c( forw_x_fourier, vF(:,jy,jz), bigvhat)
    call dfftw_execute_dft_r2c( forw_x_fourier, wF(:,jy,jz), bigwhat)
#endif
#ifdef PPSCALARS
    call dfftw_execute_dft_r2c(forw_x_fourier, theta_pert(:,jy,jz), theta_hat)
#endif

    ! Normalize transformed variables
    uhat = uhat / nxp
    vhat = vhat / nxp
    what = what / nxp
#ifdef PPOUTPUT_SSP
    bigvhat = bigvhat / nxp
    bigwhat = bigwhat / nxp
#endif
#ifdef PPSCALARS
    theta_hat = theta_hat / nxp
#endif

    ! Sum over boundary points
    uu(1) = uu(1) + 0.5d0*real(uhat(1)*conjg(uhat(1)))
    vv(1) = vv(1) + 0.5d0*real(vhat(1)*conjg(vhat(1)))
    ww(1) = ww(1) + 0.5d0*real(what(1)*conjg(what(1)))
#ifdef PPOUTPUT_SSP
    bigvv(1) = bigvv(1) + 0.5d0*real(bigvhat(1)*conjg(bigvhat(1)))
    bigww(1) = bigww(1) + 0.5d0*real(bigwhat(1)*conjg(bigwhat(1)))
#endif
#ifdef PPSCALARS
    theta2(1) = theta2(1) + 0.50d0*real(theta_hat(1)*conjg(theta_hat(1)))
#endif

    uu(nxp/2+1) = uu(nxp/2+1) + 0.5d0*real(uhat(nxp/2+1)*conjg(uhat(nxp/2+1)))
    vv(nxp/2+1) = vv(nxp/2+1) + 0.5d0*real(vhat(nxp/2+1)*conjg(vhat(nxp/2+1)))
    ww(nxp/2+1) = ww(nxp/2+1) + 0.5d0*real(what(nxp/2+1)*conjg(what(nxp/2+1)))
#ifdef PPOUTPUT_SSP
    bigvv(nxp/2+1) = bigvv(nxp/2+1) + 0.5d0*real(bigvhat(nxp/2+1)*conjg(bigvhat(nxp/2+1)))
    bigww(nxp/2+1) = bigww(nxp/2+1) + 0.5d0*real(bigwhat(nxp/2+1)*conjg(bigwhat(nxp/2+1)))
#endif
#ifdef PPSCALARS
    theta2(nxp/2+1) = theta2(nxp/2+1) + 0.5d0*real(theta_hat(nxp/2+1)*conjg(theta_hat(nxp/2+1)))
#endif

    ! Sum over interior points
    do jx = 2, ( nxp/2 + 1 - 1 )
        uu(jx) = uu(jx) + real(uhat(jx)*conjg(uhat(jx)))
        vv(jx) = vv(jx) + real(vhat(jx)*conjg(vhat(jx)))
        ww(jx) = ww(jx) + real(what(jx)*conjg(what(jx)))
#ifdef PPOUTPUT_SSP
        bigvv(jx) = bigvv(jx) + real(bigvhat(jx)*conjg(bigvhat(jx)))
        bigww(jx) = bigww(jx) + real(bigwhat(jx)*conjg(bigwhat(jx)))
#endif
#ifdef PPSCALARS
        theta2(jx) = theta2(jx) + 0.5d0*real(theta_hat(jx)*conjg(theta_hat(jx)))
#endif
    enddo

enddo
enddo

! Normalize by spanwise length - treating as a spanwise average
uu = uu / L_y
vv = vv / L_y
ww = ww / L_y
#ifdef PPOUTPUT_SSP
bigvv = bigvv / L_y
bigww = bigww / L_y
#endif
#ifdef PPSCALARS
theta2 = theta2 / L_y
#endif

! Compute total kinetic energy
ke = uu + vv + ww
#ifdef PPOUTPUT_SSP
roll = bigvv + bigww
#endif

! Write data to file
#ifdef PPMPI
call mpi_reduce (ke, ke_global, nxp/2+1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
#ifdef PPOUTPUT_SSP
call mpi_reduce (uu, uu_global, nxp/2+1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
call mpi_reduce (roll, roll_global, nxp/2+1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
#endif
#ifdef PPSCALARS
call mpi_reduce (theta2, theta2_global, nxp/2+1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
#endif
if (rank == 0) then  ! note that it's rank here, not coord
    ke = ke_global
#ifdef PPOUTPUT_SSP
    uu = uu_global
    roll = roll_global
#endif
#ifdef PPSCALARS
    theta2 = theta2_global
#endif
#endif
    ! kinetic energy of spanwise perturbations
    open(2,file=path // 'output/ke_kx.dat', status='unknown',               &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, ke
    close(2)

#ifdef PPOUTPUT_SSP
    ! Streamwise energy of spanwise perturbations - streak energy
    open(2,file=path // 'output/streak_kx.dat', status='unknown',       &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, uu
    close(2)

    ! Spanwise and Wall-normal energy - roll energy
    open(2,file=path // 'output/roll_kx.dat', status='unknown',         &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, roll
    close(2)
#endif

#ifdef PPSCALARS
    ! temperature "energy" of spanwise perturbations
     open(2,file=path // 'output/theta_kx.dat', status='unknown',       &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, theta2
    close(2)   
#endif

#ifdef PPMPI
end if
#endif

end subroutine kx_energy_fourier

!*****************************************************************************
subroutine kx_energy_by_z_fourier ()
!*****************************************************************************
! 
! Computes and writes the energy of the streamwise modes for each z-level.
! 
! This function is intended for fourier
! 
use types, only: rprec
use param
use sim_param, only: uF, vF, wF
use messages
use fft
use functions, only: int2str

implicit none

integer :: jx, jy, jz
real(rprec), dimension(nxp+2,ny,lbz:nz) :: upert, vpert, wpert
complex(rprec), dimension(nxp/2+1) :: uhat, vhat, what
real(rprec), dimension(nxp/2+1) :: ke, uu, vv, ww
character(len = 10), dimension(nxp/2+1) :: fhead !! file header variable
character(len = 20) :: fname

! Initialize variables
ke = 0.0_rprec
uu = 0.0_rprec
vv = 0.0_rprec
ww = 0.0_rprec

! Initialize file, write header
if (jt_total == wbase) then
    do jx = 1, nxp/2 + 1
        fhead(jx) = 'kx='//trim(int2str(jx-1))//' '
    enddo
endif

! Take out spanwise mean, ky = 0
call ypert_fourier(uF,upert)
call ypert_fourier(vF,vpert)
call ypert_fourier(wF,wpert)

! Sum across y, record at each z
do jz = 1, nz-1
    do jy = 1, ny

        ! Take 1D Fourier transform
        call dfftw_execute_dft_r2c( forw_x_fourier, upert(:,jy,jz), uhat)
        call dfftw_execute_dft_r2c( forw_x_fourier, vpert(:,jy,jz), vhat)
        call dfftw_execute_dft_r2c( forw_x_fourier, wpert(:,jy,jz), what)

        ! Normalize transformed variables
        uhat = uhat / nxp
        vhat = vhat / nxp
        what = what / nxp

        ! Sum over boundary points
        uu(1) = uu(1) + 0.5d0*real(uhat(1)*conjg(uhat(1)))
        vv(1) = vv(1) + 0.5d0*real(vhat(1)*conjg(vhat(1)))
        ww(1) = ww(1) + 0.5d0*real(what(1)*conjg(what(1)))

        uu(nxp/2+1) = uu(nxp/2+1) + 0.5d0*real(uhat(nxp/2+1)*conjg(uhat(nxp/2+1)))
        vv(nxp/2+1) = vv(nxp/2+1) + 0.5d0*real(vhat(nxp/2+1)*conjg(vhat(nxp/2+1)))
        ww(nxp/2+1) = ww(nxp/2+1) + 0.5d0*real(what(nxp/2+1)*conjg(what(nxp/2+1)))

        ! Sum over interior points
        do jx = 2, ( nxp/2 + 1 - 1 )
            uu(jx) = uu(jx) + real(uhat(jx)*conjg(uhat(jx)))
            vv(jx) = vv(jx) + real(vhat(jx)*conjg(vhat(jx)))
            ww(jx) = ww(jx) + real(what(jx)*conjg(what(jx)))
        enddo

    enddo

    ! Normalize by spanwise length - treating as a spanwise average
    uu = uu / L_y
    vv = vv / L_y
    ww = ww / L_y

    ! Compute total kinetic energy
    ke = uu + vv + ww

    ! Create file name
    fname = 'kx_y' !! refresh with each z
    call string_concat( fname,'_c',coord )
    call string_concat( fname,'_zi',jz )
    call string_concat( fname,'.dat' )

    ! Write data to file
    open(2,file=path // fname, status='unknown',               &
        form='formatted', position='append')
    if (jt_total==wbase) write(2,*) 'jt_total ', fhead
    write(2,*) jt_total, ke
    close(2)

enddo

end subroutine kx_energy_by_z_fourier

!*****************************************************************************
subroutine write_tau_wall_bot()
!*****************************************************************************
! 
! Write spatially-averaged statistics on wall stress
! 
use types, only : rprec
use param, only : wbase, jt_total, total_time_dim, L_x, z_i, u_star, nx, ny
use param, only : nxp, fourier
use sim_param, only : txz, tyz, txzF, tyzF
use functions, only : int2str
use param, only : lh, L_y
use fft
implicit none

real(rprec) :: turnovers, txzavg, tyzavg, twall, txzsamp, tyzsamp
integer :: jx, jy
real(rprec), dimension(ld,ny) :: txzpert, tyzpert
complex(rprec), dimension(nx/2+1) :: txzhat, tyzhat
real(rprec), dimension(nx/2+1) :: txz2, tyz2, twall2
character(len = 10), dimension(nx/2+1) :: fhead !! file header variable for kx

! ------------------------- Wall Stress Statistics  -------------------------
! Compute output to write to file
turnovers = total_time_dim / (L_x * z_i / u_star)
txzavg = 0._rprec
tyzavg = 0._rprec

if (fourier) then
    ! compute spatial average
    do jx = 1, nxp
    do jy = 1, ny
        txzavg = txzavg + txzF(jx,jy,1)
        tyzavg = tyzavg + tyzF(jx,jy,1)
    end do
    end do
    txzavg = txzavg/(nxp*ny)
    tyzavg = tyzavg/(nxp*ny)
    ! Consider fluctuations about txzavg and tyzavg at a given sample point
    txzsamp = txzF(nxp/2,ny/2,1) - txzavg
    tyzsamp = tyzF(nxp/2,ny/2,1) - tyzavg
else
    ! compute spatial average
    do jx = 1, nx
    do jy = 1, ny
        txzavg = txzavg + txz(jx,jy,1)
        tyzavg = tyzavg + tyz(jx,jy,1)
    end do
    end do
    txzavg = txzavg/(nx*ny)
    tyzavg = tyzavg/(nx*ny)
    ! Consider fluctuations about txzavg and tyzavg at a given sample point
    txzsamp = txz(nx/2,ny/2,1) - txzavg
    tyzsamp = tyz(nx/2,ny/2,1) - tyzavg
endif

twall = sqrt( (txzavg**2) + (tyzavg**2) )
! Make txzavg and tyzavg positive to compare to twall
txzavg = abs(txzavg)
tyzavg = abs(tyzavg)
! Have the sample points fluctuate about positive txzavg and tyzavg
txzsamp = txzavg + txzsamp
tyzsamp = tyzavg + tyzsamp

! ------------------------------ kx Wall Stress ------------------------------
if (.not. fourier) then

! Initialize variables
txz2 = 0.0_rprec
tyz2 = 0.0_rprec

! Initialize file, write header 
if (jt_total == wbase) then
    do jx = 1, nx/2 + 1
        fhead(jx) = 'kx='//trim(int2str(jx-1))//' '
    enddo
endif

! Take out spanwise mean, ky = 0
call ypert_by_z(txz(:,:,1),txzpert)
call ypert_by_z(tyz(:,:,1),tyzpert)

! Consider each y location
do jy = 1, ny
    ! Take 1D Fourier Transform
    call dfftw_execute_dft_r2c( forw_x, txzpert(:,jy), txzhat)
    call dfftw_execute_dft_r2c( forw_x, tyzpert(:,jy), tyzhat)

    ! Normalize transformed variables
    txzhat = txzhat / nx
    tyzhat = tyzhat / nx

    ! Take product and integrate
    do jx = 1, lh
        txz2(jx) = txz2(jx) + real(txzhat(jx)*conjg(txzhat(jx)))
        tyz2(jx) = tyz2(jx) + real(tyzhat(jx)*conjg(tyzhat(jx)))
    enddo

enddo

! Normalize by spanwise length - treating as a spanwise average
txz2 = txz2 / L_y
tyz2 = tyz2 / L_y

! Make summations
twall2 = txz2 + tyz2

endif

! ------------------------------ Write to file ------------------------------
! Wall Stress Statistics File
open(2,file=path // 'output/tau_wall_bot.dat', status='unknown',               &
    form='formatted', position='append')
!! one time header output
if (jt_total==wbase) write(2,*)                                                &
    'jt_total, turnovers, 1.0, tau_wall, txzavg, tyzavg, txzsamp, tyzsamp'
!! continual time-related output
write(2,*) jt_total, turnovers, 1.0, twall, txzavg, tyzavg, txzsamp, tyzsamp
close(2)
!! Note: sample points are at nx/2 and ny/2

! kx Wall Stress
open(2,file=path // 'output/twall_kx.dat', status='unknown',                     &
    form='formatted', position='append')
if (jt_total==wbase) write(2,*) 'jt_total ', fhead
write(2,*) jt_total, twall2
close(2)

end subroutine write_tau_wall_bot

!*****************************************************************************
subroutine write_tau_wall_top()
!*****************************************************************************
use types, only : rprec
use param, only : jt_total, total_time, total_time_dim, dt, dt_dim, wbase
use param, only : L_x, z_i, u_star
use functions, only : get_tau_wall_top
implicit none

real(rprec) :: turnovers

turnovers = total_time_dim / (L_x * z_i / u_star)

open(2,file=path // 'output/tau_wall_top.dat', status='unknown',               &
    form='formatted', position='append')

! one time header output
if (jt_total==wbase) write(2,*)                                                &
    'jt_total, total_time, total_time_dim, turnovers, dt, dt_dim, 1.0, tau_wall'

! continual time-related output
write(2,*) jt_total, total_time, total_time_dim, turnovers, dt, dt_dim,        &
    1.0, get_tau_wall_top()
close(2)

end subroutine write_tau_wall_top

#ifdef PPOUTPUT_CLOCK
!*****************************************************************************
subroutine write_clocks(t1,t2,t3,t4,t5,t6,t7,t8)
!*****************************************************************************
use param, only : jt_total, wbase
implicit none

real(rprec), intent(in) :: t1,t2,t3,t4,t5,t6,t7,t8

open(2,file=path // 'output/clocks.dat', status='unknown',                &
    form='formatted', position='append')

! one time header output
if (jt_total==wbase) write(2,*)                                           &
    'jt_total, total, deriv, stress, wall_model, divstress, convec, inter, pres'

! continual time-related output
write(2,*) jt_total, t1,t2,t3,t4,t5,t6,t7,t8

end subroutine write_clocks
#endif

#ifdef PPCGNS
#ifdef PPMPI
!*****************************************************************************
subroutine write_parallel_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,  &
    end_n_in, xin, yin, zin, num_fields, fieldNames, input )
!*****************************************************************************
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Data to be written
real(rprec), intent(in), dimension(:) :: input
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1,                                 &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2,   &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 3,   &
                            start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i=1,num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
        field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n, &
        input((i-1)*nnodes+1:(i)*nnodes), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

end subroutine write_parallel_cgns

!*****************************************************************************
subroutine write_null_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,     &
    end_n_in, xin, yin, zin, num_fields, fieldNames )
!*****************************************************************************
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 3, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i = 1, num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),  &
                           field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n, &
                                %VAL(0), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

write(*,*) "end of write_null_cgns"

end subroutine write_null_cgns
#endif
#endif

!*****************************************************************************
subroutine output_loop()
!*****************************************************************************
!
!  This subroutine is called every time step and acts as a driver for
!  computing statistics and outputing instantaneous data. No actual
!  calculations are performed here.
!
use param, only : jt_total, dt
use param, only : checkpoint_data, checkpoint_nskip
use param, only : tavg_calc, tavg_nstart, tavg_nend, tavg_nskip
use param, only : point_calc, point_nstart, point_nend, point_nskip
use param, only : domain_calc, domain_nstart, domain_nend, domain_nskip
use param, only : xplane_calc, xplane_nstart, xplane_nend, xplane_nskip
use param, only : yplane_calc, yplane_nstart, yplane_nend, yplane_nskip
use param, only : zplane_calc, zplane_nstart, zplane_nend, zplane_nskip
use stat_defs, only : tavg_initialized,tavg_dt
use derivatives, only : wave2phys, phys2wave
use param, only : fourier
use sim_param, only : u, v, w, txz, tyz
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
#ifdef PPOUTPUT_BUDGET
use sim_param, only : p, dpdx, dpdy, dpdz, divtx, divty, divtz
#ifdef PPSCALARS
use scalars, only : dTdx, dTdy, dTdz
#endif
#ifdef PPCNDIFF
use sim_param, only : txx, txy, tyy, tzz
#ifdef PPSCALARS
use derivatives, only : ddz_w
use scalars, only : temp_var, div_pi
#endif
#endif
#endif
#ifdef PPSCALARS
use scalars, only : theta, pi_z
#endif
implicit none

! Determine if we are to checkpoint intermediate times
if( checkpoint_data ) then
    ! Now check if data should be checkpointed this time step
    if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint()
end if

!  Determine if time summations are to be calculated
if (tavg_calc) then
    ! Are we between the start and stop timesteps?
    if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then
        ! Every timestep (between nstart and nend), add to tavg_dt
        tavg_dt = tavg_dt + dt

        ! Are we at the beginning or a multiple of nstart?
        if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then

#ifdef PPOUTPUT_BUDGET 
#ifdef PPCNDIFF
            ! Re-calculate stress-divergence terms since this was only partially done
            ! when using the Crank-Nicolson scheme
            call divstress_uv(divtx, divty, txx, txy, txz, tyy, tyz)
            call divstress_w(divtz, txz, tyz, tzz)

#ifdef PPSCALARS
            call ddz_w(pi_z, temp_var, lbz)
            div_pi = div_pi + temp_var
#endif

#endif
#endif

            if (fourier) then
                call wave2phys( u, lbz )
                call wave2phys( v, lbz )
                call wave2phys( w, lbz )
                call wave2phys( txz, lbz )
                call wave2phys( tyz, lbz )
                call wave2phys( dudx, lbz )
                call wave2phys( dudy, lbz )
                call wave2phys( dudz, lbz )
                call wave2phys( dvdx, lbz )
                call wave2phys( dvdy, lbz )
                call wave2phys( dvdz, lbz )
                call wave2phys( dwdx, lbz )
                call wave2phys( dwdy, lbz )
                call wave2phys( dwdz, lbz )
#ifdef PPOUTPUT_BUDGET
                call wave2phys( p, lbz )
                call wave2phys( dpdx, 1 )
                call wave2phys( dpdy, 1 )
                call wave2phys( dpdz, 1 )
                call wave2phys( divtx, lbz )
                call wave2phys( divty, lbz )
                call wave2phys( divtz, lbz )
#endif
#ifdef PPSCALARS
                call wave2phys( theta, lbz )
                call wave2phys( pi_z, lbz )
#ifdef PPOUTPUT_SPECBUDG
                call wave2phys( dTdx, lbz )
                call wave2phys( dTdy, lbz )
                call wave2phys( dTdz, lbz )
                call wave2phys( div_pi, lbz )
#endif
#endif
            endif

            ! Check if we have initialized tavg
            if (.not.tavg_initialized) then
                if (coord == 0) then
                    write(*,*) '-------------------------------'
                    write(*,"(1a,i9,1a,i9)")                               &
                        'Starting running time summation from ',           &
                        tavg_nstart, ' to ', tavg_nend
                    write(*,*) '-------------------------------'
                end if

                call tavg_init()
            else
                call tavg_compute()
            end if

            if (fourier) then
                call phys2wave( u, lbz )
                call phys2wave( v, lbz )
                call phys2wave( w, lbz )
                call phys2wave( txz, lbz )
                call phys2wave( tyz, lbz )
                call phys2wave( dudx, lbz )
                call phys2wave( dudy, lbz )
                call phys2wave( dudz, lbz )
                call phys2wave( dvdx, lbz )
                call phys2wave( dvdy, lbz )
                call phys2wave( dvdz, lbz )
                call phys2wave( dwdx, lbz )
                call phys2wave( dwdy, lbz )
                call phys2wave( dwdz, lbz )
#ifdef PPOUTPUT_BUDGET
                call phys2wave( p, lbz )
                call phys2wave( dpdx, 1 )
                call phys2wave( dpdy, 1 )
                call phys2wave( dpdz, 1 )
                call phys2wave( divtx, lbz )
                call phys2wave( divty, lbz )
                call phys2wave( divtz, lbz )
#endif
#ifdef PPSCALARS
                call phys2wave( theta, lbz )
                call phys2wave( pi_z, lbz )
#ifdef PPOUTPUT_SPECBUDG
                call phys2wave( dTdx, lbz )
                call phys2wave( dTdy, lbz )
                call phys2wave( dTdz, lbz )
                call phys2wave( div_pi, lbz )
#endif
#endif
            endif

        end if ! Are we at the beginning or a multiple of nstart?

    end if
end if

!  Determine if instantaneous point velocities are to be recorded
if(point_calc) then
    if (jt_total >= point_nstart .and. jt_total <= point_nend .and.        &
        ( mod(jt_total-point_nstart,point_nskip)==0) ) then
        if (jt_total == point_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                   &
                    'Writing instantaneous point velocities from ',        &
                    point_nstart, ' to ', point_nend
                write(*,"(1a,i9)") 'Iteration skip:', point_nskip
                write(*,*) '-------------------------------'
            end if
        end if
        call inst_write(1)
    end if
end if

!  Determine if instantaneous domain velocities are to be recorded
if(domain_calc) then
    if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and.       &
        ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
        if (jt_total == domain_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                    &
                    'Writing instantaneous domain velocities from ',        &
                    domain_nstart, ' to ', domain_nend
                write(*,"(1a,i9)") 'Iteration skip:', domain_nskip
                write(*,*) '-------------------------------'
            end if

        end if
        call inst_write(2)
    end if
end if

!  Determine if instantaneous x-plane velocities are to be recorded
if(xplane_calc) then
    if (jt_total >= xplane_nstart .and. jt_total <= xplane_nend .and.       &
        ( mod(jt_total-xplane_nstart,xplane_nskip)==0) ) then
    if (jt_total == xplane_nstart) then
        if (coord == 0) then
            write(*,*) '-------------------------------'
            write(*,"(1a,i9,1a,i9)")                                        &
                'Writing instantaneous x-plane velocities from ',           &
                xplane_nstart, ' to ', xplane_nend
            write(*,"(1a,i9)") 'Iteration skip:', xplane_nskip
            write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(3)
    end if
end if

!  Determine if instantaneous y-plane velocities are to be recorded
if(yplane_calc) then
    if (jt_total >= yplane_nstart .and. jt_total <= yplane_nend .and.       &
        ( mod(jt_total-yplane_nstart,yplane_nskip)==0) ) then
        if (jt_total == yplane_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                    &
                    'Writing instantaneous y-plane velocities from ',       &
                    yplane_nstart, ' to ', yplane_nend
                write(*,"(1a,i9)") 'Iteration skip:', yplane_nskip
                write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(4)
    end if
end if

!  Determine if instantaneous z-plane velocities are to be recorded
if(zplane_calc) then
    if (jt_total >= zplane_nstart .and. jt_total <= zplane_nend .and.       &
        ( mod(jt_total-zplane_nstart,zplane_nskip)==0) ) then
        if (jt_total == zplane_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                    &
                    'Writing instantaneous z-plane velocities from ',       &
                    zplane_nstart, ' to ', zplane_nend
                write(*,"(1a,i9)") 'Iteration skip:', zplane_nskip
                write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(5)
    end if
end if

end subroutine output_loop

!*****************************************************************************
subroutine inst_write(itype)
!*****************************************************************************
!
! This subroutine is used to write all of the instantaneous data from
! lesgo to file. The types of data written are:
!
!   points   : itype=1
!   domain   : itype=2
!   x-planes : itype=3
!   y-planes : itype=4
!   z-planes : itype=5
!
! For the points and planar data, this subroutine writes using the
! locations specfied from the param module.
! If additional instantenous values are
! desired to be written, they should be done so using this subroutine.
!
use functions, only : linear_interp, trilinear_interp, interp_to_uv_grid
use param, only : point_nloc, point_loc
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use param, only : dx, dy
use param, only : write_endian
use param, only : fourier
use grid_m
use sim_param, only : u, v, w, p
use sim_param, only : dwdy, dwdx, dvdx, dudy
use functions, only : interp_to_w_grid
use derivatives, only : wave2physF
use sim_param, only : uF, vF, wF, pF
use sim_param, only : dudyF, dudzF, dvdxF, dvdzF, dwdxF, dwdyF
use param, only : nxp
#ifdef PPSCALARS
use scalars, only : theta, thetaF
#endif

use stat_defs, only : xplane, yplane
#ifdef PPMPI
use stat_defs, only : zplane, point
use param, only : ny, nz, dz
#endif
#ifdef PPLVLSET
use level_set_base, only : phi
use sim_param, only : fx, fy, fz, fxa, fya, fza
#endif

!debug
use derivatives, only : wave2phys, phys2wave

implicit none

!debug
real(rprec), dimension(ny,lbz:nz) :: u_sum, v_sum, w_sum, u_avg, v_avg, w_avg

integer, intent(in) :: itype
character (64) :: fname
integer :: n, i, j, k
real(rprec), allocatable, dimension(:,:,:) :: ui, vi, wi, w_uv, wF_uv
real(rprec), pointer, dimension(:) :: x, y, z, zw
#ifdef PPSCALARS
real(rprec), allocatable, dimension(:,:,:) :: thetai
#endif

#ifndef PPCGNS
character(64) :: bin_ext

#ifdef PPLVLSET
real(rprec), allocatable, dimension(:,:,:) :: fx_tot, fy_tot, fz_tot
#endif

! Vorticity
real(rprec), dimension (:,:,:), allocatable :: vortx, vorty, vortz
real(rprec), dimension (:,:,:), allocatable :: vortxF, vortyF, vortzF

! Pressure
real(rprec), dimension(:,:,:), allocatable :: pres_real
real(rprec), dimension(:,:,:), allocatable :: presF_real

#ifdef PPMPI
call string_splice(bin_ext, '.c', coord, '.bin')
#else
bin_ext = '.bin'
#endif
#endif

! Nullify pointers
nullify(x,y,z,zw)

! Set grid pointers
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

if (fourier) then
    call wave2physF( u, uF )
    call wave2physF( v, vF )
    call wave2physF( w, wF )
    call wave2physF( p, pF )
    call wave2physF( dudy, dudyF )
    call wave2physF( dudz, dudzF )
    call wave2physF( dvdx, dvdxF )
    call wave2physF( dvdz, dvdzF )
    call wave2physF( dwdx, dwdxF )
    call wave2physF( dwdy, dwdyF )
#ifdef PPSCALARS
    call wave2physF( theta, thetaF )
#endif
endif

!  Allocate space for the interpolated w values
allocate(w_uv(nx,ny,lbz:nz))
allocate(wF_uv(nxp,ny,lbz:nz))

!  Make sure w has been interpolated to uv-grid
if (fourier) then
    wF_uv = interp_to_uv_grid(wF(1:nxp,1:ny,lbz:nz), lbz)
else
    w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)
endif

!  Instantaneous velocity sampled at point
if(itype==1) then
    do n = 1, point_nloc
        ! Common file name for all output types
        call string_splice(fname, path // 'output/vel.x-', point_loc(n)%xyz(1),&
            '.y-', point_loc(n)%xyz(2), '.z-', point_loc(n)%xyz(3), '.dat')

#ifdef PPMPI
        if(point(n) % coord == coord) then
#endif
            open(unit=13, position="append", file=fname)
            write(13,*) total_time,                                          &
            trilinear_interp(u(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz),    &
            trilinear_interp(v(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz),    &
            trilinear_interp(w_uv(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz)
            close(13)
#ifdef PPMPI
        end if
#endif
    end do

!  Instantaneous write for entire domain
elseif(itype==2) then
    ! Common file name for all output types
    call string_splice(fname, path //'output/vel.', jt_total)

if (fourier) then
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nxp, ny, nz - nz_end, nz_tot,           &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                     &
        (/ nxp, ny, (nz-1)*(coord+1) + 1 - nz_end /),                       &
        x(1:nxp) , y(1:ny) , z(1:(nz-nz_end) ),                             &
        3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                     &
        (/ uF(1:nxp,1:ny,1:(nz-nz_end)), vF(1:nxp,1:ny,1:(nz-nz_end)),      &
         wF_uv(1:nxp,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,     &
        access='direct', recl=nxp*ny*nz*rprec)
    write(13,rec=1) uF(:nxp,:ny,1:nz)
    write(13,rec=2) vF(:nxp,:ny,1:nz)
    write(13,rec=3) wF_uv(:nxp,:ny,1:nz)
    close(13)
#endif
else !! .not. fourier
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,            &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                     &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                        &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                              &
        3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                     &
        (/ u(1:nx,1:ny,1:(nz-nz_end)), v(1:nx,1:ny,1:(nz-nz_end)),          &
         w_uv(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,     &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) u(:nx,:ny,1:nz)
    write(13,rec=2) v(:nx,:ny,1:nz)
    write(13,rec=3) w_uv(:nx,:ny,1:nz)
    close(13)
#endif
endif

    ! Compute vorticity
    allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))
    allocate(vortxF(nxp,ny,lbz:nz), vortyF(nxp,ny,lbz:nz), vortzF(nxp,ny,lbz:nz))

    ! Use vorticityx as an intermediate step for performing uv-w interpolation
    ! Vorticity is written in w grid
    if (fourier) then
        vortxF(1:nxp,1:ny,lbz:nz) = 0._rprec
        vortyF(1:nxp,1:ny,lbz:nz) = 0._rprec
        vortzF(1:nxp,1:ny,lbz:nz) = 0._rprec

        vortxF(1:nxp,1:ny,lbz:nz) = dvdxF(1:nxp,1:ny,lbz:nz) -            &
            dudyF(1:nxp,1:ny,lbz:nz)
        vortzF(1:nxp,1:ny,lbz:nz) =                                       &
            interp_to_w_grid( vortxF(1:nxp,1:ny,lbz:nz), lbz)
        vortxF(1:nxp,1:ny,lbz:nz) = dwdyF(1:nxp,1:ny,lbz:nz) -            &
            dvdzF(1:nxp,1:ny,lbz:nz)
        vortyF(1:nxp,1:ny,lbz:nz) = dudzF(1:nxp,1:ny,lbz:nz) -            &
            dwdxF(1:nxp,1:ny,lbz:nz)

        if (coord == 0) then
            vortzF(1:nxp,1:ny, 1) = 0._rprec
        end if
    else
        vortx(1:nx,1:ny,lbz:nz) = 0._rprec
        vorty(1:nx,1:ny,lbz:nz) = 0._rprec
        vortz(1:nx,1:ny,lbz:nz) = 0._rprec

        vortx(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) - dudy(1:nx,1:ny,lbz:nz)
        vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid( vortx(1:nx,1:ny,lbz:nz), lbz)
        vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) - dvdz(1:nx,1:ny,lbz:nz)
        vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) - dwdx(1:nx,1:ny,lbz:nz)

        if (coord == 0) then
            vortz(1:nx,1:ny, 1) = 0._rprec
        end if
    endif

    ! Common file name for all output types
    call string_splice(fname, path //'output/vort.', jt_total)

if (fourier) then
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname,nxp,ny, nz - nz_end, nz_tot,             &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                     &
        (/ nxp, ny, (nz-1)*(coord+1) + 1 - nz_end /),                       &
        x(1:nxp) , y(1:ny) , zw(1:(nz-nz_end) ),                            &
        3, (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                  &
        (/ vortxF(1:nxp,1:ny,1:(nz-nz_end)), vortyF(1:nxp,1:ny,1:(nz-nz_end)),  &
        vortzF(1:nxp,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,     &
        access='direct', recl=nxp*ny*nz*rprec)
    write(13,rec=1) vortxF(:nxp,:ny,1:nz)
    write(13,rec=2) vortyF(:nxp,:ny,1:nz)
    write(13,rec=3) vortzF(:nxp,:ny,1:nz)
    close(13)
#endif
else !! .not. fourier
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname,nx,ny, nz - nz_end, nz_tot,              &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                     &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                        &
        x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                             &
        3, (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                  &
        (/ vortx(1:nx,1:ny,1:(nz-nz_end)), vorty(1:nx,1:ny,1:(nz-nz_end)),  &
        vortz(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,     &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) vortx(:nx,:ny,1:nz)
    write(13,rec=2) vorty(:nx,:ny,1:nz)
    write(13,rec=3) vortz(:nx,:ny,1:nz)
    close(13)
#endif
endif !! fourier

    deallocate(vortx, vorty, vortz)
    deallocate(vortxF, vortyF, vortzF)

    ! Compute pressure
    allocate(pres_real(nx,ny,lbz:nz))
    allocate(presF_real(nxp,ny,lbz:nz))
    if (fourier) then
        presF_real(1:nxp,1:ny,lbz:nz) = 0._rprec

        presF_real(1:nxp,1:ny,lbz:nz) = pF(1:nxp,1:ny,lbz:nz)               &
            - 0.5 * ( uF(1:nxp,1:ny,lbz:nz)**2                              &
            + wF_uv(1:nxp,1:ny,lbz:nz)**2                                   &
            + vF(1:nxp,1:ny,lbz:nz)**2 )
    else !! .not. fourier
        pres_real(1:nx,1:ny,lbz:nz) = 0._rprec

        pres_real(1:nx,1:ny,lbz:nz) = p(1:nx,1:ny,lbz:nz)                   &
            - 0.5 * ( u(1:nx,1:ny,lbz:nz)**2                                &
            + w_uv(1:nx,1:ny,lbz:nz)**2                                     &
            + v(1:nx,1:ny,lbz:nz)**2 )
    endif

    ! Common file name for all output types
    call string_splice(fname, path //'output/pres.', jt_total)

if (fourier) then
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nxp, ny, nz - nz_end, nz_tot,           &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                     &
        (/ nxp, ny, (nz-1)*(coord+1) + 1 - nz_end /),                       &
        x(1:nxp) , y(1:ny) , z(1:(nz-nz_end) ),                             &
        1, (/ 'Pressure' /), (/ presF_real(1:nxp,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,     &
        access='direct', recl=nxp*ny*nz*rprec)
    write(13,rec=1) presF_real(:nxp,:ny,1:nz)
    close(13)
#endif
else !! .not. fourier
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,            &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                     &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                        &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                              &
        1, (/ 'Pressure' /), (/ pres_real(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,     &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) pres_real(:nx,:ny,1:nz)
    close(13)
#endif
endif

     deallocate(pres_real)
     deallocate(presF_real)

#ifdef PPSCALARS
    ! Common file name for all output types
    call string_splice(fname, path //'output/theta.', jt_total)
if (fourier) then
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nxp, ny, nz - nz_end, nz_tot,              &
     (/ 1, 1,   (nz-1)*coord + 1 /),                                           &
     (/ nxp, ny, (nz-1)*(coord+1) + 1 - nz_end /),                             &
     x(1:nxp) , y(1:ny) , z(1:(nz-nz_end) ),                                   &
     1, (/ 'Theta' /), (/ thetaF(1:nxp,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
     access='direct', recl=nxp*ny*nz*rprec)
    write(13,rec=1) thetaF(:nxp,:ny,1:nz)
    close(13)
#endif
else !! .not. fourier
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
     (/ 1, 1,   (nz-1)*coord + 1 /),                                           &
     (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                              &
     x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                    &
     1, (/ 'Theta' /), (/ theta(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
     access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) theta(:nx,:ny,1:nz)
    close(13)
#endif
endif
#endif

!  Write instantaneous x-plane values
elseif(itype==3) then

if (fourier) then

    ! Not going to interpolate since physical grid is probably coarse
    ! Instead going to use arbitrary xplane location, ignoring user-input
    ! Assuming user only wants a single xplane location

    ! Common file name portion for all output types
    call string_splice(fname, path // 'output/vel.x-', xplane_loc(1), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        ! Write CGNS Output
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,1,ny, nz - nz_end, nz_tot,     &
                        (/ 1, 1,   (nz-1)*coord + 1 /),                &
                        (/ 1, ny, (nz-1)*(coord+1) + 1 - nz_end /),    &
                    xplane_loc(i:i) , y(1:ny) , z(1:(nz-nz_end) ),     &
              3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),          &
              (/ uF(1,1:ny,1:(nz-nz_end)), vF(1,1:ny,1:(nz-nz_end)),   &
                 wF(1,1:ny,1:(nz-nz_end)) /) )

#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=ny*nz*rprec)
        write(13,rec=1) uF(1,:ny,1:nz)
        write(13,rec=2) vF(1,:ny,1:nz)
        write(13,rec=3) wF(1,:ny,1:nz)
        close(13)
#endif

    ! Now to compute vorticity
    allocate(vortxF(nxp,ny,lbz:nz), vortyF(nxp,ny,lbz:nz), vortzF(nxp,ny,lbz:nz))

    ! Use vorticityx as an intermediate step for performing uv-w interpolation
    ! Vorticity is written in w grid
    vortxF(1:nxp,1:ny,lbz:nz) = 0._rprec
    vortyF(1:nxp,1:ny,lbz:nz) = 0._rprec
    vortzF(1:nxp,1:ny,lbz:nz) = 0._rprec

    vortxF(1:nxp,1:ny,lbz:nz) = dvdxF(1:nxp,1:ny,lbz:nz) -            &
        dudyF(1:nxp,1:ny,lbz:nz)
    vortzF(1:nxp,1:ny,lbz:nz) =                                       &
        interp_to_w_grid( vortxF(1:nxp,1:ny,lbz:nz), lbz)
    vortxF(1:nxp,1:ny,lbz:nz) = dwdyF(1:nxp,1:ny,lbz:nz) -            &
        dvdzF(1:nxp,1:ny,lbz:nz)
    vortyF(1:nxp,1:ny,lbz:nz) = dudzF(1:nxp,1:ny,lbz:nz) -            &
        dwdxF(1:nxp,1:ny,lbz:nz)

    if (coord == 0) then
        vortzF(1:nxp,1:ny, 1) = 0._rprec
    end if
 
    ! Common file name for all output types
    call string_splice(fname, path // 'output/vort.x-', xplane_loc(1), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname,1,ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                     &
        (/ nxp, ny, (nz-1)*(coord+1) + 1 - nz_end /),                       &
        x(1:nxp) , y(1:ny) , zw(1:(nz-nz_end) ),                            &
        3, (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                  &
        (/ vortxF(1,1:ny,1:(nz-nz_end)), vortyF(1,1:ny,1:(nz-nz_end)),      &
        vortzF(1,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,     &
        access='direct', recl=ny*nz*rprec)
    write(13,rec=1) vortxF(1,:ny,1:nz)
    write(13,rec=2) vortyF(1,:ny,1:nz)
    write(13,rec=3) vortzF(1,:ny,1:nz)
    close(13)
#endif

    deallocate(vortxF, vortyF, vortzF)

#ifdef PPSCALARS
    ! Common file name portion for all output types
    call string_splice(fname, path // 'output/theta.x-', xplane_loc(1), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        ! Write CGNS Output
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,1,ny, nz - nz_end, nz_tot,     &
                        (/ 1, 1,   (nz-1)*coord + 1 /),                &
                        (/ 1, ny, (nz-1)*(coord+1) + 1 - nz_end /),    &
                    xplane_loc(i:i) , y(1:ny) , z(1:(nz-nz_end) ),     &
              1, (/ 'Theta' /),                                        &
              (/ thetaF(1,1:ny,1:(nz-nz_end)) /) )

#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=ny*nz*rprec)
        write(13,rec=1) thetaF(1,:ny,1:nz)
        close(13)
#endif
#endif

else !! not fourier
    allocate(ui(1,ny,nz), vi(1,ny,nz), wi(1,ny,nz))

    !  Loop over all xplane locations
    do i = 1, xplane_nloc
        do k = 1, nz
            do j = 1, ny
                ui(1,j,k) = linear_interp(u(xplane(i) % istart,j,k),    &
                     u(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
                vi(1,j,k) = linear_interp(v(xplane(i) % istart,j,k),    &
                     v(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
                wi(1,j,k) = linear_interp(w_uv(xplane(i) % istart,j,k), &
                     w_uv(xplane(i) % istart+1,j,k), dx, &
                     xplane(i) % ldiff)
            end do
        end do

        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.x-', xplane_loc(i), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        ! Write CGNS Output
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,1,ny, nz - nz_end, nz_tot,     &
                        (/ 1, 1,   (nz-1)*coord + 1 /),                &
                        (/ 1, ny, (nz-1)*(coord+1) + 1 - nz_end /),    &
                    xplane_loc(i:i) , y(1:ny) , z(1:(nz-nz_end) ),     &
              3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),          &
              (/ ui(1,1:ny,1:(nz-nz_end)), vi(1,1:ny,1:(nz-nz_end)),   &
                 wi(1,1:ny,1:(nz-nz_end)) /) )

#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=ny*nz*rprec)
        write(13,rec=1) ui
        write(13,rec=2) vi
        write(13,rec=3) wi
        close(13)
#endif
    end do

    deallocate(ui,vi,wi)

    ! Output vorticity on arbitrary x=0 plane same as fourier, not interpolating
    ! Now to compute vorticity
    allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))

    ! Use vorticityx as an intermediate step for performing uv-w interpolation
    ! Vorticity is written in w grid
    vortx(1:nx,1:ny,lbz:nz) = 0._rprec
    vorty(1:nx,1:ny,lbz:nz) = 0._rprec
    vortz(1:nx,1:ny,lbz:nz) = 0._rprec

    vortx(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) -            &
        dudy(1:nx,1:ny,lbz:nz)
    vortz(1:nx,1:ny,lbz:nz) =                                     &
        interp_to_w_grid( vortx(1:nx,1:ny,lbz:nz), lbz)
    vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) -            &
        dvdz(1:nx,1:ny,lbz:nz)
    vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) -            &
        dwdx(1:nx,1:ny,lbz:nz)

    if (coord == 0) then
        vortz(1:nx,1:ny, 1) = 0._rprec
    end if
 
    ! Common file name for all output types
    call string_splice(fname, path // 'output/vort.x-', xplane_loc(1), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname,1,ny, nz - nz_end, nz_tot,                   &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
        (/ nxp, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nxp) , y(1:ny) , zw(1:(nz-nz_end) ),                                &
        3, (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                      &
        (/ vortx(1,1:ny,1:(nz-nz_end)), vorty(1,1:ny,1:(nz-nz_end)),  &
        vortz(1,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=ny*nz*rprec)
    write(13,rec=1) vortx(1,:ny,1:nz)
    write(13,rec=2) vorty(1,:ny,1:nz)
    write(13,rec=3) vortz(1,:ny,1:nz)
    close(13)
#endif

    deallocate(vortx, vorty, vortz)

#ifdef PPSCALARS
    ! Common file name portion for all output types
    call string_splice(fname, path // 'output/theta.x-', xplane_loc(1), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        ! Write CGNS Output
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,1,ny, nz - nz_end, nz_tot,     &
                        (/ 1, 1,   (nz-1)*coord + 1 /),                &
                        (/ 1, ny, (nz-1)*(coord+1) + 1 - nz_end /),    &
                    xplane_loc(i:i) , y(1:ny) , z(1:(nz-nz_end) ),     &
              1, (/ 'Theta' /),                                        &
              (/ theta(1,1:ny,1:(nz-nz_end)) /) )

#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=ny*nz*rprec)
        write(13,rec=1) theta(1,:ny,1:nz)
        close(13)
#endif
#endif

endif !! to fourier or not to fourier

!  Write instantaneous y-plane values
elseif(itype==4) then

    allocate(ui(nx,1,nz), vi(nx,1,nz), wi(nx,1,nz))

    !  Loop over all yplane locations
    do j = 1, yplane_nloc
        do k = 1, nz
            do i = 1, nx

                ui(i,1,k) = linear_interp(u(i,yplane(j) % istart,k),        &
                     u(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
                vi(i,1,k) = linear_interp(v(i,yplane(j) % istart,k),        &
                     v(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
                wi(i,1,k) = linear_interp(w_uv(i,yplane(j) % istart,k),     &
                     w_uv(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
            end do
        end do

        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.y-', yplane_loc(j), '.', &
             jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,nx,1, nz - nz_end, nz_tot,          &
            (/ 1, 1,   (nz-1)*coord + 1 /),                                 &
            (/ nx, 1, (nz-1)*(coord+1) + 1 - nz_end /),                     &
            x(1:nx) , yplane_loc(j:j) , z(1:(nz-nz_end) ),                  &
            3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                 &
            (/ ui(1:nx,1,1:(nz-nz_end)), vi(1:nx,1,1:(nz-nz_end)),          &
            wi(1:nx,1,1:(nz-nz_end)) /) )
#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=nx*nz*rprec)
        write(13,rec=1) ui
        write(13,rec=2) vi
        write(13,rec=3) wi
        close(13)
#endif

    end do

    deallocate(ui,vi,wi)

!  Write instantaneous z-plane values
elseif (itype==5) then

    allocate(ui(nx,ny,1), vi(nx,ny,1), wi(nx,ny,1))

    !  Loop over all zplane locations
    do k = 1, zplane_nloc
        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.z-',                  &
                zplane_loc(k), '.', jt_total)

#ifdef PPCGNS
        call string_concat(fname, '.cgns')
#endif

#ifdef PPMPI
        if(zplane(k) % coord == coord) then
            do j = 1, Ny
                do i = 1, Nx
                    ui(i,j,1) = linear_interp(u(i,j,zplane(k) % istart),     &
                         u(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                    vi(i,j,1) = linear_interp(v(i,j,zplane(k) % istart),     &
                         v(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                    wi(i,j,1) = linear_interp(w_uv(i,j,zplane(k) % istart),  &
                         w_uv(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                end do
            end do

#ifdef PPCGNS
            call warn("inst_write","Z plane writting is currently disabled.")
!            ! Write CGNS Data
!            ! Only the processor with data writes, the other one is written
!            ! using null arguments with 'write_null_cgns'
!            call write_parallel_cgns (fname ,nx, ny, 1, 1,                  &
!                (/ 1, 1,   1 /),                                            &
!                (/ nx, ny, 1 /),                                            &
!                x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                     &
!                (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                &
!                (/ ui(1:nx,1:ny,1), vi(1:nx,1:ny,1), wi(1:nx,1:ny,1) /) )
#else
            call string_concat(fname, bin_ext)
            open(unit=13,file=fname,form='unformatted',convert=write_endian, &
                            access='direct',recl=nx*ny*1*rprec)
            write(13,rec=1) ui(1:nx,1:ny,1)
            write(13,rec=2) vi(1:nx,1:ny,1)
            write(13,rec=3) wi(1:nx,1:ny,1)
            close(13)
#endif
!
! #ifdef PPMPI
!         else
! #ifdef PPCGNS
!            write(*,*) "At write_null_cgns"
!            call write_null_cgns (fname ,nx, ny, 1, 1,                      &
!            (/ 1, 1,   1 /),                                                &
!            (/ nx, ny, 1 /),                                                &
!            x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                         &
!            (/ 'VelocityX', 'VelocityY', 'VelocityZ' /) )
!#endif
        end if
#endif
    end do
    deallocate(ui,vi,wi)

#ifdef PPSCALARS
    allocate(thetai(nx,ny,1))

    !  Loop over all zplane locations
    do k = 1, zplane_nloc
        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/theta.z-',                     &
                zplane_loc(k), '.', jt_total)

#ifdef PPCGNS
        call string_concat(fname, '.cgns')
#endif

#ifdef PPMPI
        if(zplane(k) % coord == coord) then
            do j = 1, Ny
                do i = 1, Nx
                    thetai(i,j,1) = linear_interp(theta(i,j,zplane(k) % istart),       &
                         theta(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                end do
            end do

#ifdef PPCGNS
            call warn("inst_write","Z plane writting is currently disabled.")
!            ! Write CGNS Data
!            ! Only the processor with data writes, the other one is written
!            ! using null arguments with 'write_null_cgns'
!            call write_parallel_cgns (fname ,nx, ny, 1, 1,                     &
!                (/ 1, 1,   1 /),                                               &
!                (/ nx, ny, 1 /),                                               &
!                x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                        &
!                (/ 'Theta' /),                   &
!                (/ thetai(1:nx,1:ny,1) /) )
#else
            call string_concat(fname, bin_ext)
            open(unit=13,file=fname,form='unformatted',convert=write_endian,   &
                            access='direct',recl=nx*ny*1*rprec)
            write(13,rec=1) thetai(1:nx,1:ny,1)
            close(13)
#endif
!
! #ifdef PPMPI
!         else
! #ifdef PPCGNS
!            write(*,*) "At write_null_cgns"
!            call write_null_cgns (fname ,nx, ny, 1, 1,                         &
!            (/ 1, 1,   1 /),                                                   &
!            (/ nx, ny, 1 /),                                                   &
!            x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                            &
!            (/ 'Theta' /) )
!#endif
        end if
#endif
    end do
    deallocate(thetai)
#endif

else
    write(*,*) 'Error: itype not specified properly to inst_write!'
    stop
end if

deallocate(w_uv)
nullify(x,y,z,zw)

#ifdef PPLVLSET
contains
!*****************************************************************************
subroutine force_tot()
!*****************************************************************************
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
implicit none

! Zero bogus values
fx(:,:,nz) = 0._rprec
fy(:,:,nz) = 0._rprec
fz(:,:,nz) = 0._rprec

!  Sum both the induced and applied forces
allocate(fx_tot(nx,ny,nz), fy_tot(nx,ny,nz), fz_tot(nx,ny,nz))

#ifdef PPTURBINES
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = fya(1:nx,1:ny,1:nz)
fz_tot = fza(1:nx,1:ny,1:nz)

#elif PPATM
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = fya(1:nx,1:ny,1:nz)
fz_tot = fza(1:nx,1:ny,1:nz)

#elif PPLVLSET
fx_tot = fx(1:nx,1:ny,1:nz)+fxa(1:nx,1:ny,1:nz)
fy_tot = fy(1:nx,1:ny,1:nz)+fya(1:nx,1:ny,1:nz)
fz_tot = fz(1:nx,1:ny,1:nz)+fza(1:nx,1:ny,1:nz)
#else
fx_tot = 0._rprec
fy_tot = 0._rprec
fz_tot = 0._rprec
#endif

#ifdef PPMPI
!  Sync forces
call mpi_sync_real_array( fx_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fy_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fz_tot, 1, MPI_SYNC_DOWN )
#endif

! Put fz_tot on uv-grid
fz_tot(1:nx,1:ny,1:nz) = interp_to_uv_grid( fz_tot(1:nx,1:ny,1:nz), 1 )

return
end subroutine force_tot
#endif

!!****************************************************************************
!subroutine pressure_sync()
!!****************************************************************************
!use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
!use param, only : ld
!implicit none
!
!! Reset bogus values
!p(:,:,nz) = p(:,:,nz-1)
!dpdx(:,:,nz) = dpdx(:,:,nz-1)
!dpdy(:,:,nz) = dpdy(:,:,nz-1)
!dpdz(:,:,nz) = dpdz(:,:,nz-1)
!
!#ifdef PPMPI
!!  Sync pressure
!call mpi_sync_real_array( p, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdx, 1 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdy, 1 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdz, 1 , MPI_SYNC_DOWN )
!#endif
!
!return
!end subroutine pressure_sync
!
!!****************************************************************************
!subroutine RHS_sync()
!!****************************************************************************
!use param, only : ld
!use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
!implicit none
!
!! Reset bogus values
!RHSx(:,:,nz) = RHSx(:,:,nz-1)
!RHSy(:,:,nz) = RHSy(:,:,nz-1)
!RHSz(:,:,nz) = RHSz(:,:,nz-1)
!
!#ifdef PPMPI
!!  Sync RHS
!call mpi_sync_real_array( RHSx, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( RHSy, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( RHSz, 0 , MPI_SYNC_DOWN )
!#endif
!
!return
!end subroutine RHS_sync

end subroutine inst_write

!*****************************************************************************
subroutine checkpoint ()
!*****************************************************************************
use iwmles
use param, only : nz, checkpoint_file, tavg_calc, lbc_mom, L_x, L_y, L_z
#ifdef PPMPI
use param, only : comm, ierr
#endif
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
use param, only : jt_total, total_time, total_time_dim, dt,                 &
    use_cfl_dt, cfl, write_endian
use cfl_util, only : get_max_cfl
use stat_defs, only : tavg_initialized
use string_util, only : string_concat
#if PPUSE_TURBINES
use turbines, only : turbines_checkpoint
#endif
#ifdef PPSCALARS
use scalars, only : scalars_checkpoint
#endif

! HIT Inflow
#ifdef PPHIT
use hit_inflow, only : hit_write_restart
#endif

implicit none
character(64) :: fname
real(rprec) :: cfl_w

fname = checkpoint_file
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif

!  Open vel.out (lun_default in io) for final output
open(11, file=fname, form='unformatted', convert=write_endian,              &
    status='unknown', position='rewind')
write (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),                     &
    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),                   &
    Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),                      &
    F_QN(:,:,1:nz), F_NN(:,:,1:nz)
close(11)

! Open grid.out for final output
if (coord == 0) then
    open(11, file='grid.out', form='unformatted', convert=write_endian)
    write(11) nproc, Nx, Ny, Nz, L_x, L_y, L_z
    close(11)
end if

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

! Checkpoint time averaging restart data
if ( tavg_calc .and. tavg_initialized ) call tavg_checkpoint()

! Write time and current simulation state
! Set the current cfl to a temporary (write) value based whether CFL is
! specified or must be computed
if( use_cfl_dt ) then
    cfl_w = cfl
else
    cfl_w = get_max_cfl()
end if

!xiang check point for iwm
if(lbc_mom==3)then
    if (coord == 0) call iwm_checkPoint()
end if

#ifdef PPHIT
    if (coord == 0) call hit_write_restart()
#endif

#if PPUSE_TURBINES
call turbines_checkpoint
#endif

#ifdef PPSCALARS
call scalars_checkpoint
#endif

!  Update total_time.dat after simulation
if (coord == 0) then
    !--only do this for true final output, not intermediate recording
    open (1, file=fcumulative_time)
    write(1, *) jt_total, total_time, total_time_dim, dt, cfl_w
    close(1)
end if

end subroutine checkpoint

!*****************************************************************************
subroutine output_final()
!*****************************************************************************
use stat_defs, only : tavg_initialized
use param, only : tavg_calc
implicit none

! Perform final checkpoint
call checkpoint()

!  Check if average quantities are to be recorded
if (tavg_calc .and. tavg_initialized ) call tavg_finalize()

end subroutine output_final

!*****************************************************************************
subroutine output_init ()
!*****************************************************************************
!
!  This subroutine allocates the memory for arrays used for statistical
!  calculations
!
use param, only : dx, dy, dz, nz, lbz
use param, only : point_calc, point_nloc, point_loc
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
use param, only : tavg_calc
use grid_m
use functions, only : cell_indx
use stat_defs, only : point, xplane, yplane, zplane
use stat_defs, only : tavg, tavg_zplane
#ifdef PPSCALARS
use stat_defs, only : tavg_scal
#ifdef PPOUTPUT_TURBSPEC
use stat_defs, only : tavg_scal_turbspecx, tavg_scal_turbspecy
#endif
#ifdef PPOUTPUT_SPECBUDG
use stat_defs, only : tavg_scal_specbudgx
#endif
#endif
use stat_defs, only : tavg_vort
#ifdef PPOUTPUT_SGS
use stat_defs, only : tavg_sgs
#endif
#ifdef PPOUTPUT_BUDGET
use stat_defs, only : tavg_budget
#endif
#ifdef PPOUTPUT_TURBSPEC
use stat_defs, only : tavg_turbspecx, tavg_turbspecy
#endif
#ifdef PPOUTPUT_SPECBUDG
use stat_defs, only : tavg_specbudgx, tavg_specbudgy
#endif

implicit none

integer :: i,j,k
real(rprec), pointer, dimension(:) :: x,y,z


#ifdef PPMPI
! This adds one more element to the last processor (which contains an extra one)
! Processor nproc-1 has data from 1:nz
! Rest of processors have data from 1:nz-1
if ( coord == nproc-1 ) then
    nz_end = 0
else
    nz_end = 1
end if
#else
nz_end = 0
#endif

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

if( tavg_calc ) then

    allocate(tavg(nx,ny,lbz:nz))
    ! allocate(tavg_zplane(nz))
#ifdef PPSCALARS
    allocate(tavg_scal(nx,ny,lbz:nz))
#ifdef PPOUTPUT_TURBSPEC
    allocate(tavg_scal_turbspecx(nx/2+1,ny,lbz:nz))
    allocate(tavg_scal_turbspecy(nx,ny/2+1,lbz:nz))
#endif
#ifdef PPOUTPUT_SPECBUDG
    allocate(tavg_scal_specbudgx(nx/2+1,ny,lbz:nz))
#endif
#endif
    allocate(tavg_vort(nx,ny,lbz:nz))
#ifdef PPOUTPUT_SGS
    allocate(tavg_sgs(nx,ny,lbz:nz))
#endif
#ifdef PPOUTPUT_BUDGET
    allocate(tavg_budget(nx,ny,lbz:nz))
#endif
#ifdef PPOUTPUT_TURBSPEC
    allocate(tavg_turbspecx(nx/2+1,ny,lbz:nz))
    allocate(tavg_turbspecy(nx,ny/2+1,lbz:nz))
#endif
#ifdef PPOUTPUT_SPECBUDG
    allocate(tavg_specbudgx(nx/2+1,ny,lbz:nz))
    allocate(tavg_specbudgy(nx,ny/2+1,lbz:nz))
#endif

  ! Initialize the derived types tavg and tavg_zplane
    do k = 1, Nz
        do j = 1, Ny
        do i = 1, Nx
            tavg(i,j,k) % u    = 0._rprec
            tavg(i,j,k) % v    = 0._rprec
            tavg(i,j,k) % w_uv = 0._rprec
            tavg(i,j,k) % u_w  = 0._rprec
            tavg(i,j,k) % v_w  = 0._rprec
            tavg(i,j,k) % w    = 0._rprec
            tavg(i,j,k) % u2   = 0._rprec
            tavg(i,j,k) % v2   = 0._rprec
            tavg(i,j,k) % w2   = 0._rprec
            tavg(i,j,k) % uv   = 0._rprec
            tavg(i,j,k) % uw   = 0._rprec
            tavg(i,j,k) % vw   = 0._rprec
            tavg(i,j,k) % txx  = 0._rprec
            tavg(i,j,k) % tyy  = 0._rprec
            tavg(i,j,k) % tzz  = 0._rprec
            tavg(i,j,k) % txy  = 0._rprec
            tavg(i,j,k) % txz  = 0._rprec
            tavg(i,j,k) % tyz  = 0._rprec
            tavg(i,j,k) % fx   = 0._rprec
            tavg(i,j,k) % fy   = 0._rprec
            tavg(i,j,k) % fz   = 0._rprec

            tavg_vort(i,j,k) % vortx = 0._rprec
            tavg_vort(i,j,k) % vorty = 0._rprec
            tavg_vort(i,j,k) % vortz = 0._rprec
            tavg_vort(i,j,k) % vortx2 = 0._rprec
            tavg_vort(i,j,k) % vorty2 = 0._rprec
            tavg_vort(i,j,k) % vortz2 = 0._rprec

        end do
        end do
        ! type_set removed, is tavg_zplane even used?
        ! to replace commented line below, set all types of tavg_zplane = 0._rprec
        ! call type_set( tavg_zplane(k), 0._rprec )
    end do

#ifdef PPSCALARS
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx
        tavg_scal(i,j,k) % theta  = 0._rprec
        tavg_scal(i,j,k) % theta2 = 0._rprec
        tavg_scal(i,j,k) % utheta = 0._rprec
        tavg_scal(i,j,k) % vtheta = 0._rprec
        tavg_scal(i,j,k) % wtheta = 0._rprec
        tavg_scal(i,j,k) % pi_z   = 0._rprec
    end do
    end do
    end do

#ifdef PPOUTPUT_TURBSPEC
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx/2 + 1

        tavg_scal_turbspecx(i,j,k) % thetaf = 0.d0

        tavg_scal_turbspecx(i,j,k) % theta2 = 0.d0
        tavg_scal_turbspecx(i,j,k) % utheta = 0.d0
        tavg_scal_turbspecx(i,j,k) % vtheta = 0.d0
        tavg_scal_turbspecx(i,j,k) % wtheta = 0.d0

    end do
    end do
    end do

    do k = 1, Nz
    do j = 1, Ny/2 + 1
    do i = 1, Nx

        tavg_scal_turbspecy(i,j,k) % thetaf = 0.d0

        tavg_scal_turbspecy(i,j,k) % theta2 = 0.d0
        tavg_scal_turbspecy(i,j,k) % utheta = 0.d0
        tavg_scal_turbspecy(i,j,k) % vtheta = 0.d0
        tavg_scal_turbspecy(i,j,k) % wtheta = 0.d0

    end do
    end do
    end do
#endif

#ifdef PPOUTPUT_SPECBUDG
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx/2 + 1

        ! Scalar gradients
        tavg_scal_specbudgx(i,j,k) % dTdx = 0.d0
        tavg_scal_specbudgx(i,j,k) % dTdy = 0.d0
        tavg_scal_specbudgx(i,j,k) % dTdz = 0.d0

        ! Scalar-Scalar Gradient, Th_dTdxjh
        tavg_scal_specbudgx(i,j,k) % Th_dTdxh = 0.d0
        tavg_scal_specbudgx(i,j,k) % Th_dTdyh = 0.d0
        tavg_scal_specbudgx(i,j,k) % Th_dTdzh = 0.d0

        ! Velocity-Scalar Gradient, (uj*dTdxj)h
        tavg_scal_specbudgx(i,j,k) % udTdxh = 0.d0
        tavg_scal_specbudgx(i,j,k) % vdTdyh = 0.d0
        tavg_scal_specbudgx(i,j,k) % wdTdzh = 0.d0

        ! Scalar-Velocity-Scalar Gradient, Th_(uj*dTdxj)h
        tavg_scal_specbudgx(i,j,k) % Th_udTdxh = 0.d0
        tavg_scal_specbudgx(i,j,k) % Th_vdTdyh = 0.d0
        tavg_scal_specbudgx(i,j,k) % Th_wdTdzh = 0.d0

        ! Scalar Gradient-Scalar Gradient
        tavg_scal_specbudgx(i,j,k) % TxTx = 0.d0
        tavg_scal_specbudgx(i,j,k) % TyTy = 0.d0
        tavg_scal_specbudgx(i,j,k) % TzTz = 0.d0

        ! Laplacian, (nu/Pr)*lap(theta)
        tavg_scal_specbudgx(i,j,k) % lapT = 0.d0

        ! Scalar-Laplacian, (theta*(nu/Pr)*lap(theta)
        tavg_scal_specbudgx(i,j,k) % TlapT = 0.d0

    end do
    end do
    end do
#endif

#endif

#ifdef PPOUTPUT_SGS
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx
        tavg_sgs(i,j,k) % cs_opt2  = 0._rprec
        ! tavg_sgs(i,j,k) % Tn       = 0._rprec
        tavg_sgs(i,j,k) % Nu_t     = 0._rprec
        ! tavg_sgs(i,j,k) % F_LM     = 0._rprec
        ! tavg_sgs(i,j,k) % F_MM     = 0._rprec
        ! tavg_sgs(i,j,k) % F_QN     = 0._rprec
        ! tavg_sgs(i,j,k) % F_NN     = 0._rprec
        ! tavg_sgs(i,j,k) % ee_now   = 0._rprec
! #ifdef PPDYN_TN
        ! tavg_sgs(i,j,k) % F_ee2    = 0._rprec
        ! tavg_sgs(i,j,k) % F_deedt2 = 0._rprec
! #endif
    end do
    end do
    end do
#endif

#ifdef PPOUTPUT_BUDGET
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx
        ! Mean pressure on w-grid
        tavg_budget(i,j,k) % p = 0._rprec

        ! Mean velocity product, ui*uj, on w-grid
        tavg_budget(i,j,k) % uu = 0._rprec
        tavg_budget(i,j,k) % vv = 0._rprec
        tavg_budget(i,j,k) % ww = 0._rprec
        tavg_budget(i,j,k) % uv = 0._rprec
        tavg_budget(i,j,k) % uw = 0._rprec
        tavg_budget(i,j,k) % vw = 0._rprec

        ! Mean Velocity Gradients, duidxj
        tavg_budget(i,j,k) % dudx = 0._rprec
        tavg_budget(i,j,k) % dudy = 0._rprec
        tavg_budget(i,j,k) % dudz = 0._rprec
        tavg_budget(i,j,k) % dvdx = 0._rprec
        tavg_budget(i,j,k) % dvdy = 0._rprec
        tavg_budget(i,j,k) % dvdz = 0._rprec
        tavg_budget(i,j,k) % dwdx = 0._rprec
        tavg_budget(i,j,k) % dwdy = 0._rprec
        tavg_budget(i,j,k) % dwdz = 0._rprec

        ! Mean Pressure Gradients, dpdxi
        tavg_budget(i,j,k) % dpdx = 0._rprec
        tavg_budget(i,j,k) % dpdy = 0._rprec
        tavg_budget(i,j,k) % dpdz = 0._rprec

        ! Mean vel-velGrad product, ui*dujdxk
        tavg_budget(i,j,k) % ududx = 0._rprec
        tavg_budget(i,j,k) % ududy = 0._rprec
        tavg_budget(i,j,k) % ududz = 0._rprec
        tavg_budget(i,j,k) % udvdx = 0._rprec
        tavg_budget(i,j,k) % udvdy = 0._rprec
        tavg_budget(i,j,k) % udvdz = 0._rprec
        tavg_budget(i,j,k) % udwdx = 0._rprec
        tavg_budget(i,j,k) % udwdy = 0._rprec
        tavg_budget(i,j,k) % udwdz = 0._rprec

        tavg_budget(i,j,k) % vdudx = 0._rprec
        tavg_budget(i,j,k) % vdudy = 0._rprec
        tavg_budget(i,j,k) % vdudz = 0._rprec
        tavg_budget(i,j,k) % vdvdx = 0._rprec
        tavg_budget(i,j,k) % vdvdy = 0._rprec
        tavg_budget(i,j,k) % vdvdz = 0._rprec
        tavg_budget(i,j,k) % vdwdx = 0._rprec
        tavg_budget(i,j,k) % vdwdy = 0._rprec
        tavg_budget(i,j,k) % vdwdz = 0._rprec

        tavg_budget(i,j,k) % wdudx = 0._rprec
        tavg_budget(i,j,k) % wdudy = 0._rprec
        tavg_budget(i,j,k) % wdudz = 0._rprec
        tavg_budget(i,j,k) % wdvdx = 0._rprec
        tavg_budget(i,j,k) % wdvdy = 0._rprec
        tavg_budget(i,j,k) % wdvdz = 0._rprec
        tavg_budget(i,j,k) % wdwdx = 0._rprec
        tavg_budget(i,j,k) % wdwdy = 0._rprec
        tavg_budget(i,j,k) % wdwdz = 0._rprec

        ! Mean vel-vel-velGrad product, ui*uk*dujdxk
        tavg_budget(i,j,k) % uududx = 0._rprec
        tavg_budget(i,j,k) % uvdudy = 0._rprec
        tavg_budget(i,j,k) % uwdudz = 0._rprec
        tavg_budget(i,j,k) % uudvdx = 0._rprec
        tavg_budget(i,j,k) % uvdvdy = 0._rprec
        tavg_budget(i,j,k) % uwdvdz = 0._rprec
        tavg_budget(i,j,k) % uudwdx = 0._rprec
        tavg_budget(i,j,k) % uvdwdy = 0._rprec
        tavg_budget(i,j,k) % uwdwdz = 0._rprec

        tavg_budget(i,j,k) % vududx = 0._rprec
        tavg_budget(i,j,k) % vvdudy = 0._rprec
        tavg_budget(i,j,k) % vwdudz = 0._rprec
        tavg_budget(i,j,k) % vudvdx = 0._rprec
        tavg_budget(i,j,k) % vvdvdy = 0._rprec
        tavg_budget(i,j,k) % vwdvdz = 0._rprec
        tavg_budget(i,j,k) % vudwdx = 0._rprec
        tavg_budget(i,j,k) % vvdwdy = 0._rprec
        tavg_budget(i,j,k) % vwdwdz = 0._rprec

        tavg_budget(i,j,k) % wududx = 0._rprec
        tavg_budget(i,j,k) % wvdudy = 0._rprec
        tavg_budget(i,j,k) % wwdudz = 0._rprec
        tavg_budget(i,j,k) % wudvdx = 0._rprec
        tavg_budget(i,j,k) % wvdvdy = 0._rprec
        tavg_budget(i,j,k) % wwdvdz = 0._rprec
        tavg_budget(i,j,k) % wudwdx = 0._rprec
        tavg_budget(i,j,k) % wvdwdy = 0._rprec
        tavg_budget(i,j,k) % wwdwdz = 0._rprec

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
        tavg_budget(i,j,k) % uxux = 0._rprec
        tavg_budget(i,j,k) % uyuy = 0._rprec
        tavg_budget(i,j,k) % uzuz = 0._rprec
        tavg_budget(i,j,k) % vxvx = 0._rprec
        tavg_budget(i,j,k) % vyvy = 0._rprec
        tavg_budget(i,j,k) % vzvz = 0._rprec
        tavg_budget(i,j,k) % wxwx = 0._rprec
        tavg_budget(i,j,k) % wywy = 0._rprec
        tavg_budget(i,j,k) % wzwz = 0._rprec

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
        tavg_budget(i,j,k) % uxvx = 0._rprec
        tavg_budget(i,j,k) % uyvy = 0._rprec
        tavg_budget(i,j,k) % uzvz = 0._rprec
        tavg_budget(i,j,k) % uxwx = 0._rprec
        tavg_budget(i,j,k) % uywy = 0._rprec
        tavg_budget(i,j,k) % uzwz = 0._rprec
        tavg_budget(i,j,k) % vxwx = 0._rprec
        tavg_budget(i,j,k) % vywy = 0._rprec
        tavg_budget(i,j,k) % vzwz = 0._rprec

        ! duidxj*dujdxi, i /= j
        ! tavg_budget(i,j,k) % uyvx = 0._rprec
        ! tavg_budget(i,j,k) % uzwx = 0._rprec
        ! tavg_budget(i,j,k) % vzwy = 0._rprec

        ! Mean vel-presGrad product, ui*dpdxj
        tavg_budget(i,j,k) % udpdx = 0._rprec
        tavg_budget(i,j,k) % udpdy = 0._rprec
        tavg_budget(i,j,k) % udpdz = 0._rprec
        tavg_budget(i,j,k) % vdpdx = 0._rprec
        tavg_budget(i,j,k) % vdpdy = 0._rprec
        tavg_budget(i,j,k) % vdpdz = 0._rprec
        tavg_budget(i,j,k) % wdpdx = 0._rprec
        tavg_budget(i,j,k) % wdpdy = 0._rprec
        tavg_budget(i,j,k) % wdpdz = 0._rprec

        ! Mean pres-velGrad product, p*duidxj
        tavg_budget(i,j,k) % pdudx = 0._rprec
        tavg_budget(i,j,k) % pdudy = 0._rprec
        tavg_budget(i,j,k) % pdudz = 0._rprec
        tavg_budget(i,j,k) % pdvdx = 0._rprec
        tavg_budget(i,j,k) % pdvdy = 0._rprec
        tavg_budget(i,j,k) % pdvdz = 0._rprec
        tavg_budget(i,j,k) % pdwdx = 0._rprec
        tavg_budget(i,j,k) % pdwdy = 0._rprec
        tavg_budget(i,j,k) % pdwdz = 0._rprec

        ! Mean Laplacian, nu*lap(uj)
        tavg_budget(i,j,k) % lapu = 0._rprec
        tavg_budget(i,j,k) % lapv = 0._rprec
        tavg_budget(i,j,k) % lapw = 0._rprec

        ! Mean Vel-Laplacian, nu*ui*lap(uj)
        tavg_budget(i,j,k) % ulapu = 0._rprec
        tavg_budget(i,j,k) % ulapv = 0._rprec
        tavg_budget(i,j,k) % ulapw = 0._rprec
        tavg_budget(i,j,k) % vlapu = 0._rprec
        tavg_budget(i,j,k) % vlapv = 0._rprec
        tavg_budget(i,j,k) % vlapw = 0._rprec
        tavg_budget(i,j,k) % wlapu = 0._rprec
        tavg_budget(i,j,k) % wlapv = 0._rprec
        tavg_budget(i,j,k) % wlapw = 0._rprec
    end do
    end do
    end do
#endif

#ifdef PPOUTPUT_TURBSPEC
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx/2 + 1

        tavg_turbspecx(i,j,k) % uf = 0.d0
        tavg_turbspecx(i,j,k) % vf = 0.d0
        tavg_turbspecx(i,j,k) % wf = 0.d0

        tavg_turbspecx(i,j,k) % uu = 0.d0
        tavg_turbspecx(i,j,k) % vv = 0.d0
        tavg_turbspecx(i,j,k) % ww = 0.d0
        !tavg_turbspecx(i,j,k) % vel2 = 0.d0

        tavg_turbspecx(i,j,k) % uv = 0.d0
        tavg_turbspecx(i,j,k) % uw = 0.d0
        tavg_turbspecx(i,j,k) % vw = 0.d0

        tavg_turbspecx(i,j,k) % vortxf = 0.d0
        tavg_turbspecx(i,j,k) % vortyf = 0.d0
        tavg_turbspecx(i,j,k) % vortzf = 0.d0

        tavg_turbspecx(i,j,k) % vortx2 = 0.d0
        tavg_turbspecx(i,j,k) % vorty2 = 0.d0
        tavg_turbspecx(i,j,k) % vortz2 = 0.d0
        !tavg_turbspecx(i,j,k) % vort2 = 0.d0

    end do
    end do
    end do

    do k = 1, Nz
    do j = 1, Ny/2 + 1
    do i = 1, Nx

        tavg_turbspecy(i,j,k) % uf = 0.d0
        tavg_turbspecy(i,j,k) % vf = 0.d0
        tavg_turbspecy(i,j,k) % wf = 0.d0

        tavg_turbspecy(i,j,k) % uu = 0.d0
        tavg_turbspecy(i,j,k) % vv = 0.d0
        tavg_turbspecy(i,j,k) % ww = 0.d0
        !tavg_turbspecy(i,j,k) % vel2 = 0.d0

        tavg_turbspecy(i,j,k) % uv = 0.d0
        tavg_turbspecy(i,j,k) % uw = 0.d0
        tavg_turbspecy(i,j,k) % vw = 0.d0

        tavg_turbspecy(i,j,k) % vortxf = 0.d0
        tavg_turbspecy(i,j,k) % vortyf = 0.d0
        tavg_turbspecy(i,j,k) % vortzf = 0.d0

        tavg_turbspecy(i,j,k) % vortx2 = 0.d0
        tavg_turbspecy(i,j,k) % vorty2 = 0.d0
        tavg_turbspecy(i,j,k) % vortz2 = 0.d0
        !tavg_turbspecy(i,j,k) % vort2 = 0.d0
    end do
    end do
    end do
#endif

#ifdef PPOUTPUT_SPECBUDG
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx/2 + 1
        ! Mean pressure
        tavg_specbudgx(i,j,k) % p = 0.d0

        ! Mean pressure gradients
        tavg_specbudgx(i,j,k) % dpdx = 0.d0
        tavg_specbudgx(i,j,k) % dpdy = 0.d0
        tavg_specbudgx(i,j,k) % dpdz = 0.d0

        ! Mean velocity gradients
        tavg_specbudgx(i,j,k) % dudx = 0.d0
        tavg_specbudgx(i,j,k) % dudy = 0.d0
        tavg_specbudgx(i,j,k) % dudz = 0.d0
        tavg_specbudgx(i,j,k) % dvdx = 0.d0
        tavg_specbudgx(i,j,k) % dvdy = 0.d0
        tavg_specbudgx(i,j,k) % dvdz = 0.d0
        tavg_specbudgx(i,j,k) % dwdx = 0.d0
        tavg_specbudgx(i,j,k) % dwdy = 0.d0
        tavg_specbudgx(i,j,k) % dwdz = 0.d0

        ! Mean vel-velGrad product, ui_hat*dujdxk_hat
        tavg_specbudgx(i,j,k) % uh_dudxh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dudyh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dudzh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dvdxh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dvdyh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dvdzh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dwdxh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dwdyh = 0.d0
        tavg_specbudgx(i,j,k) % uh_dwdzh = 0.d0

        tavg_specbudgx(i,j,k) % vh_dudxh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dudyh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dudzh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dvdxh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dvdyh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dvdzh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dwdxh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dwdyh = 0.d0
        tavg_specbudgx(i,j,k) % vh_dwdzh = 0.d0

        tavg_specbudgx(i,j,k) % wh_dudxh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dudyh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dudzh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dvdxh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dvdyh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dvdzh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dwdxh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dwdyh = 0.d0
        tavg_specbudgx(i,j,k) % wh_dwdzh = 0.d0

        ! Mean vel-velGrad product, uk*duidxk_hat
        tavg_specbudgx(i,j,k) % ududxh = 0.d0
        tavg_specbudgx(i,j,k) % vdudyh = 0.d0
        tavg_specbudgx(i,j,k) % wdudzh = 0.d0
        tavg_specbudgx(i,j,k) % udvdxh = 0.d0
        tavg_specbudgx(i,j,k) % vdvdyh = 0.d0
        tavg_specbudgx(i,j,k) % wdvdzh = 0.d0
        tavg_specbudgx(i,j,k) % udwdxh = 0.d0
        tavg_specbudgx(i,j,k) % vdwdyh = 0.d0
        tavg_specbudgx(i,j,k) % wdwdzh = 0.d0

        ! Mean vel-vel-velGrad product, uj_hat*(uk*duidxk)_hat
        tavg_specbudgx(i,j,k) % uh_ududxh = 0.d0
        tavg_specbudgx(i,j,k) % uh_vdudyh = 0.d0
        tavg_specbudgx(i,j,k) % uh_wdudzh = 0.d0
        tavg_specbudgx(i,j,k) % vh_ududxh = 0.d0
        tavg_specbudgx(i,j,k) % vh_vdudyh = 0.d0
        tavg_specbudgx(i,j,k) % vh_wdudzh = 0.d0
        tavg_specbudgx(i,j,k) % wh_ududxh = 0.d0
        tavg_specbudgx(i,j,k) % wh_vdudyh = 0.d0
        tavg_specbudgx(i,j,k) % wh_wdudzh = 0.d0

        tavg_specbudgx(i,j,k) % uh_udvdxh = 0.d0
        tavg_specbudgx(i,j,k) % uh_vdvdyh = 0.d0
        tavg_specbudgx(i,j,k) % uh_wdvdzh = 0.d0
        tavg_specbudgx(i,j,k) % vh_udvdxh = 0.d0
        tavg_specbudgx(i,j,k) % vh_vdvdyh = 0.d0
        tavg_specbudgx(i,j,k) % vh_wdvdzh = 0.d0
        tavg_specbudgx(i,j,k) % wh_udvdxh = 0.d0
        tavg_specbudgx(i,j,k) % wh_vdvdyh = 0.d0
        tavg_specbudgx(i,j,k) % wh_wdvdzh = 0.d0

        tavg_specbudgx(i,j,k) % uh_udwdxh = 0.d0
        tavg_specbudgx(i,j,k) % uh_vdwdyh = 0.d0
        tavg_specbudgx(i,j,k) % uh_wdwdzh = 0.d0
        tavg_specbudgx(i,j,k) % vh_udwdxh = 0.d0
        tavg_specbudgx(i,j,k) % vh_vdwdyh = 0.d0
        tavg_specbudgx(i,j,k) % vh_wdwdzh = 0.d0
        tavg_specbudgx(i,j,k) % wh_udwdxh = 0.d0
        tavg_specbudgx(i,j,k) % wh_vdwdyh = 0.d0
        tavg_specbudgx(i,j,k) % wh_wdwdzh = 0.d0

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
        tavg_specbudgx(i,j,k) % uxux = 0.d0
        tavg_specbudgx(i,j,k) % uyuy = 0.d0
        tavg_specbudgx(i,j,k) % uzuz = 0.d0
        tavg_specbudgx(i,j,k) % vxvx = 0.d0
        tavg_specbudgx(i,j,k) % vyvy = 0.d0
        tavg_specbudgx(i,j,k) % vzvz = 0.d0
        tavg_specbudgx(i,j,k) % wxwx = 0.d0
        tavg_specbudgx(i,j,k) % wywy = 0.d0
        tavg_specbudgx(i,j,k) % wzwz = 0.d0

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
        tavg_specbudgx(i,j,k) % uxvx = 0.d0
        tavg_specbudgx(i,j,k) % uyvy = 0.d0
        tavg_specbudgx(i,j,k) % uzvz = 0.d0
        tavg_specbudgx(i,j,k) % uxwx = 0.d0
        tavg_specbudgx(i,j,k) % uywy = 0.d0
        tavg_specbudgx(i,j,k) % uzwz = 0.d0
        tavg_specbudgx(i,j,k) % vxwx = 0.d0
        tavg_specbudgx(i,j,k) % vywy = 0.d0
        tavg_specbudgx(i,j,k) % vzwz = 0.d0

        ! Mean pressure * velocity gradient product, p*dujdxk
        tavg_specbudgx(i,j,k) % pdudx = 0.d0
        tavg_specbudgx(i,j,k) % pdudy = 0.d0
        tavg_specbudgx(i,j,k) % pdudz = 0.d0
        tavg_specbudgx(i,j,k) % pdvdx = 0.d0
        tavg_specbudgx(i,j,k) % pdvdy = 0.d0
        tavg_specbudgx(i,j,k) % pdvdz = 0.d0
        tavg_specbudgx(i,j,k) % pdwdx = 0.d0
        tavg_specbudgx(i,j,k) % pdwdy = 0.d0
        tavg_specbudgx(i,j,k) % pdwdz = 0.d0

        ! Mean velocity * pressure gradient product, uj*dpdxk
        tavg_specbudgx(i,j,k) % udpdx = 0.d0
        tavg_specbudgx(i,j,k) % udpdy = 0.d0
        tavg_specbudgx(i,j,k) % udpdz = 0.d0
        tavg_specbudgx(i,j,k) % vdpdx = 0.d0
        tavg_specbudgx(i,j,k) % vdpdy = 0.d0
        tavg_specbudgx(i,j,k) % vdpdz = 0.d0
        tavg_specbudgx(i,j,k) % wdpdx = 0.d0
        tavg_specbudgx(i,j,k) % wdpdy = 0.d0
        tavg_specbudgx(i,j,k) % wdpdz = 0.d0

        ! Mean Laplacian, nu*lap(uj)
        tavg_specbudgx(i,j,k) % lapu = 0._rprec
        tavg_specbudgx(i,j,k) % lapv = 0._rprec
        tavg_specbudgx(i,j,k) % lapw = 0._rprec

        ! Mean Vel-Laplacian, nu*ui*lap(uj)
        tavg_specbudgx(i,j,k) % ulapu = 0._rprec
        tavg_specbudgx(i,j,k) % ulapv = 0._rprec
        tavg_specbudgx(i,j,k) % ulapw = 0._rprec
        tavg_specbudgx(i,j,k) % vlapu = 0._rprec
        tavg_specbudgx(i,j,k) % vlapv = 0._rprec
        tavg_specbudgx(i,j,k) % vlapw = 0._rprec
        tavg_specbudgx(i,j,k) % wlapu = 0._rprec
        tavg_specbudgx(i,j,k) % wlapv = 0._rprec
        tavg_specbudgx(i,j,k) % wlapw = 0._rprec

    end do
    end do
    end do

    do k = 1, Nz
    do j = 1, Ny/2+1
    do i = 1, Nx
        ! Mean pressure
        tavg_specbudgy(i,j,k) % p = 0.d0

        ! Mean pressure gradients
        tavg_specbudgy(i,j,k) % dpdx = 0.d0
        tavg_specbudgy(i,j,k) % dpdy = 0.d0
        tavg_specbudgy(i,j,k) % dpdz = 0.d0

        ! Mean velocity gradients
        tavg_specbudgy(i,j,k) % dudx = 0.d0
        tavg_specbudgy(i,j,k) % dudy = 0.d0
        tavg_specbudgy(i,j,k) % dudz = 0.d0
        tavg_specbudgy(i,j,k) % dvdx = 0.d0
        tavg_specbudgy(i,j,k) % dvdy = 0.d0
        tavg_specbudgy(i,j,k) % dvdz = 0.d0
        tavg_specbudgy(i,j,k) % dwdx = 0.d0
        tavg_specbudgy(i,j,k) % dwdy = 0.d0
        tavg_specbudgy(i,j,k) % dwdz = 0.d0

        ! Mean vel-velGrad product, ui_hat*dujdxk_hat
        tavg_specbudgy(i,j,k) % uh_dudxh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dudyh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dudzh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dvdxh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dvdyh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dvdzh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dwdxh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dwdyh = 0.d0
        tavg_specbudgy(i,j,k) % uh_dwdzh = 0.d0

        tavg_specbudgy(i,j,k) % vh_dudxh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dudyh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dudzh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dvdxh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dvdyh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dvdzh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dwdxh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dwdyh = 0.d0
        tavg_specbudgy(i,j,k) % vh_dwdzh = 0.d0

        tavg_specbudgy(i,j,k) % wh_dudxh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dudyh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dudzh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dvdxh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dvdyh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dvdzh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dwdxh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dwdyh = 0.d0
        tavg_specbudgy(i,j,k) % wh_dwdzh = 0.d0

        ! Mean vel-velGrad product, uk*duidxk_hat
        tavg_specbudgy(i,j,k) % ududxh = 0.d0
        tavg_specbudgy(i,j,k) % vdudyh = 0.d0
        tavg_specbudgy(i,j,k) % wdudzh = 0.d0
        tavg_specbudgy(i,j,k) % udvdxh = 0.d0
        tavg_specbudgy(i,j,k) % vdvdyh = 0.d0
        tavg_specbudgy(i,j,k) % wdvdzh = 0.d0
        tavg_specbudgy(i,j,k) % udwdxh = 0.d0
        tavg_specbudgy(i,j,k) % vdwdyh = 0.d0
        tavg_specbudgy(i,j,k) % wdwdzh = 0.d0

        ! Mean vel-vel-velGrad product, uj_hat*(uk*duidxk)_hat
        tavg_specbudgy(i,j,k) % uh_ududxh = 0.d0
        tavg_specbudgy(i,j,k) % uh_vdudyh = 0.d0
        tavg_specbudgy(i,j,k) % uh_wdudzh = 0.d0
        tavg_specbudgy(i,j,k) % vh_ududxh = 0.d0
        tavg_specbudgy(i,j,k) % vh_vdudyh = 0.d0
        tavg_specbudgy(i,j,k) % vh_wdudzh = 0.d0
        tavg_specbudgy(i,j,k) % wh_ududxh = 0.d0
        tavg_specbudgy(i,j,k) % wh_vdudyh = 0.d0
        tavg_specbudgy(i,j,k) % wh_wdudzh = 0.d0

        tavg_specbudgy(i,j,k) % uh_udvdxh = 0.d0
        tavg_specbudgy(i,j,k) % uh_vdvdyh = 0.d0
        tavg_specbudgy(i,j,k) % uh_wdvdzh = 0.d0
        tavg_specbudgy(i,j,k) % vh_udvdxh = 0.d0
        tavg_specbudgy(i,j,k) % vh_vdvdyh = 0.d0
        tavg_specbudgy(i,j,k) % vh_wdvdzh = 0.d0
        tavg_specbudgy(i,j,k) % wh_udvdxh = 0.d0
        tavg_specbudgy(i,j,k) % wh_vdvdyh = 0.d0
        tavg_specbudgy(i,j,k) % wh_wdvdzh = 0.d0

        tavg_specbudgy(i,j,k) % uh_udwdxh = 0.d0
        tavg_specbudgy(i,j,k) % uh_vdwdyh = 0.d0
        tavg_specbudgy(i,j,k) % uh_wdwdzh = 0.d0
        tavg_specbudgy(i,j,k) % vh_udwdxh = 0.d0
        tavg_specbudgy(i,j,k) % vh_vdwdyh = 0.d0
        tavg_specbudgy(i,j,k) % vh_wdwdzh = 0.d0
        tavg_specbudgy(i,j,k) % wh_udwdxh = 0.d0
        tavg_specbudgy(i,j,k) % wh_vdwdyh = 0.d0
        tavg_specbudgy(i,j,k) % wh_wdwdzh = 0.d0

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
        tavg_specbudgy(i,j,k) % uxux = 0.d0
        tavg_specbudgy(i,j,k) % uyuy = 0.d0
        tavg_specbudgy(i,j,k) % uzuz = 0.d0
        tavg_specbudgy(i,j,k) % vxvx = 0.d0
        tavg_specbudgy(i,j,k) % vyvy = 0.d0
        tavg_specbudgy(i,j,k) % vzvz = 0.d0
        tavg_specbudgy(i,j,k) % wxwx = 0.d0
        tavg_specbudgy(i,j,k) % wywy = 0.d0
        tavg_specbudgy(i,j,k) % wzwz = 0.d0

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
        tavg_specbudgy(i,j,k) % uxvx = 0.d0
        tavg_specbudgy(i,j,k) % uyvy = 0.d0
        tavg_specbudgy(i,j,k) % uzvz = 0.d0
        tavg_specbudgy(i,j,k) % uxwx = 0.d0
        tavg_specbudgy(i,j,k) % uywy = 0.d0
        tavg_specbudgy(i,j,k) % uzwz = 0.d0
        tavg_specbudgy(i,j,k) % vxwx = 0.d0
        tavg_specbudgy(i,j,k) % vywy = 0.d0
        tavg_specbudgy(i,j,k) % vzwz = 0.d0

        ! Mean pressure * velocity gradient product, p*dujdxk
        tavg_specbudgy(i,j,k) % pdudx = 0.d0
        tavg_specbudgy(i,j,k) % pdudy = 0.d0
        tavg_specbudgy(i,j,k) % pdudz = 0.d0
        tavg_specbudgy(i,j,k) % pdvdx = 0.d0
        tavg_specbudgy(i,j,k) % pdvdy = 0.d0
        tavg_specbudgy(i,j,k) % pdvdz = 0.d0
        tavg_specbudgy(i,j,k) % pdwdx = 0.d0
        tavg_specbudgy(i,j,k) % pdwdy = 0.d0
        tavg_specbudgy(i,j,k) % pdwdz = 0.d0

        ! Mean velocity * pressure gradient product, uj*dpdxk
        tavg_specbudgy(i,j,k) % udpdx = 0.d0
        tavg_specbudgy(i,j,k) % udpdy = 0.d0
        tavg_specbudgy(i,j,k) % udpdz = 0.d0
        tavg_specbudgy(i,j,k) % vdpdx = 0.d0
        tavg_specbudgy(i,j,k) % vdpdy = 0.d0
        tavg_specbudgy(i,j,k) % vdpdz = 0.d0
        tavg_specbudgy(i,j,k) % wdpdx = 0.d0
        tavg_specbudgy(i,j,k) % wdpdy = 0.d0
        tavg_specbudgy(i,j,k) % wdpdz = 0.d0

        ! Mean Laplacian, nu*lap(uj)
        tavg_specbudgy(i,j,k) % lapu = 0._rprec
        tavg_specbudgy(i,j,k) % lapv = 0._rprec
        tavg_specbudgy(i,j,k) % lapw = 0._rprec

        ! Mean Vel-Laplacian, nu*ui*lap(uj)
        tavg_specbudgy(i,j,k) % ulapu = 0._rprec
        tavg_specbudgy(i,j,k) % ulapv = 0._rprec
        tavg_specbudgy(i,j,k) % ulapw = 0._rprec
        tavg_specbudgy(i,j,k) % vlapu = 0._rprec
        tavg_specbudgy(i,j,k) % vlapv = 0._rprec
        tavg_specbudgy(i,j,k) % vlapw = 0._rprec
        tavg_specbudgy(i,j,k) % wlapu = 0._rprec
        tavg_specbudgy(i,j,k) % wlapv = 0._rprec
        tavg_specbudgy(i,j,k) % wlapw = 0._rprec

    end do
    end do
    end do
#endif

end if

! Initialize information for x-planar stats/data
if (xplane_calc) then
    allocate(xplane(xplane_nloc))
    xplane(:) % istart = -1
    xplane(:) % ldiff = 0.

    !  Compute istart and ldiff
    do i = 1, xplane_nloc
        xplane(i) % istart = cell_indx('i', dx, xplane_loc(i))
        xplane(i) % ldiff = xplane_loc(i) - x(xplane(i) % istart)
    end do
end if

! Initialize information for y-planar stats/data
if (yplane_calc) then
    allocate(yplane(yplane_nloc))
    yplane(:) % istart = -1
    yplane(:) % ldiff = 0.

    !  Compute istart and ldiff
    do j = 1, yplane_nloc
        yplane(j) % istart = cell_indx('j', dy, yplane_loc(j))
        yplane(j) % ldiff = yplane_loc(j) - y(yplane(j) % istart)
    end do
end if

! Initialize information for z-planar stats/data
if(zplane_calc) then
    allocate(zplane(zplane_nloc))

    !  Initialize
    zplane(:) % istart = -1
    zplane(:) % ldiff = 0.
    zplane(:) % coord = -1

    !  Compute istart and ldiff
    do k = 1, zplane_nloc

#ifdef PPMPI
        if (zplane_loc(k) >= z(1) .and. zplane_loc(k) < z(nz)) then
            zplane(k) % coord = coord
            zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
            zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
        end if
#else
        zplane(k) % coord = 0
        zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
        zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
#endif
    end do
end if

!  Open files for instantaneous writing
if (point_calc) then
    allocate(point(point_nloc))

    !  Intialize the coord values
    ! (-1 shouldn't be used as coord so initialize to this)
    point % coord=-1
    point % fid = -1

    do i = 1, point_nloc
        !  Find the processor in which this point lives
#ifdef PPMPI
        if (point_loc(i)%xyz(3) >= z(1) .and. point_loc(i)%xyz(3) < z(nz)) then
#endif

            point(i) % coord = coord

            point(i) % istart = cell_indx('i',dx,point_loc(i)%xyz(1))
            point(i) % jstart = cell_indx('j',dy,point_loc(i)%xyz(2))
            point(i) % kstart = cell_indx('k',dz,point_loc(i)%xyz(3))

            point(i) % xdiff = point_loc(i)%xyz(1) - x(point(i) % istart)
            point(i) % ydiff = point_loc(i)%xyz(2) - y(point(i) % jstart)
            point(i) % zdiff = point_loc(i)%xyz(3) - z(point(i) % kstart)

#ifdef PPMPI
        end if
#endif
    end do
end if

nullify(x,y,z)

end subroutine output_init

!*******************************************************************************
subroutine tavg_init()
!*******************************************************************************
!
!  This subroutine loads the tavg.out files
!
use messages
use param, only : read_endian
use stat_defs, only : tavg, tavg_total_time, tavg_dt, tavg_initialized
use stat_defs, only : tavg_vort
#ifdef PPOUTPUT_SGS
use stat_defs, only : tavg_sgs
#endif
#ifdef PPOUTPUT_BUDGET
use stat_defs, only : tavg_budget
#endif
#ifdef PPOUTPUT_TURBSPEC
use stat_defs, only : tavg_turbspecx, tavg_turbspecy
#endif
#ifdef PPSCALARS
use stat_defs, only : tavg_scal
#ifdef PPOUTPUT_TURBSPEC
use stat_defs, only : tavg_scal_turbspecx, tavg_scal_turbspecy
#endif
#ifdef PPOUTPUT_SPECBUDG
use stat_defs, only : tavg_scal_specbudgx
#endif
#endif
#ifdef PPOUTPUT_SPECBUDG
use stat_defs, only : tavg_specbudgx, tavg_specbudgy
#endif

implicit none

character (*), parameter :: ftavg_in = path // 'tavg.out'
character (*), parameter :: ftavg_vort_in = path // 'tavg_vort.out'
#ifdef PPOUTPUT_SGS
character (*), parameter :: ftavg_sgs_in = path // 'tavg_sgs.out'
#endif
#ifdef PPOUTPUT_BUDGET
character (*), parameter :: ftavg_budget_in = path // 'tavg_budget.out'
#endif
#ifdef PPOUTPUT_TURBSPEC
character (*), parameter :: ftavg_turbspec_in = path // 'tavg_turbspec.out'
#endif
#ifdef PPSCALARS
character (*), parameter :: ftavg_scal_in = path // 'tavg_scal.out'
#ifdef PPOUTPUT_TURBSPEC
character (*), parameter :: ftavg_scal_turbspec_in = path // 'tavg_scal_turbspec.out'
#endif
#ifdef PPOUTPUT_SPECBUDG
character (*), parameter :: ftavg_scal_specbudg_in = path // 'tavg_scal_specbudg.out'
#endif
#endif
#ifdef PPOUTPUT_SPECBUDG
character (*), parameter :: ftavg_specbudg_in = path // 'tavg_specbudg.out'
#endif
#ifdef PPMPI
character (*), parameter :: MPI_suffix = '.c'
#endif
character (128) :: fname

logical :: exst

fname = ftavg_in
#ifdef PPMPI
call string_concat( fname, MPI_suffix, coord )
#endif

inquire (file=fname, exist=exst)
if (.not. exst) then
    !  Nothing to read in
    if (coord == 0) then
        write(*,*) ' '
        write(*,*)'No previous time averaged data - starting from scratch.'
    end if
    ! note: tavg was already initialized to zero in output_init routine
    tavg_total_time = 0._rprec
else
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg
    close(1)

    fname = ftavg_vort_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_vort
    close(1)

#ifdef PPOUTPUT_SGS
    fname = ftavg_sgs_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_sgs
    close(1)
#endif

#ifdef PPOUTPUT_BUDGET
    fname = ftavg_budget_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_budget
    close(1)
#endif

#ifdef PPOUTPUT_TURBSPEC
    fname = ftavg_turbspec_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_turbspecx
    read(1) tavg_turbspecy
    close(1)
#endif

#ifdef PPSCALARS
    fname = ftavg_scal_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_scal
    close(1)

#ifdef PPOUTPUT_TURBSPEC
    fname = ftavg_scal_turbspec_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_scal_turbspecx
    read(1) tavg_scal_turbspecy
    close(1)
#endif

#ifdef PPOUTPUT_SPECBUDG
    fname = ftavg_scal_specbudg_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_scal_specbudgx
    close(1)
#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
    fname = ftavg_specbudg_in
#ifdef PPMPI
    call string_concat( fname, MPI_suffix, coord )
#endif
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) tavg_total_time
    read(1) tavg_specbudgx
    read(1) tavg_specbudgy
    close(1)
#endif

end if

! Initialize tavg_dt
tavg_dt = 0._rprec

! Set global switch that tavg as been initialized
tavg_initialized = .true.

end subroutine tavg_init

!*******************************************************************************
subroutine tavg_compute()
!*******************************************************************************
!
!  This subroutine collects the stats for each flow
!  variable quantity
!
use stat_defs, only : tavg, tavg_total_time, tavg_dt
use stat_defs, only : tavg_vort
use param, only : nx, ny, nz, lbz, jzmax, ubc_mom, lbc_mom
use sim_param, only : u, v, w, p
use sim_param, only : txx, txy, tyy, txz, tyz, tzz
use sim_param, only : dudy, dudz, dvdx, dvdz, dwdx, dwdy
!use sim_param, only : dudx, dvdy, dwdz

#ifdef PPOUTPUT_SGS
use stat_defs, only : tavg_sgs
use sgs_param
#endif

#ifdef PPTURBINES
use sim_param, only : fxa, fya, fza
#endif

#ifdef PPSCALARS
use stat_defs, only : tavg_scal
use scalars, only : theta, pi_z
#endif

use functions, only : interp_to_uv_grid, interp_to_w_grid

implicit none

integer :: i, j, k
real(rprec) :: u_p, u_p2, v_p, v_p2, w_p, w_p2
real(rprec), allocatable, dimension(:,:,:) :: w_uv, u_w, v_w
real(rprec), allocatable, dimension(:,:,:) :: pres_real
real(rprec), allocatable, dimension(:,:,:) :: dwdx_uv, dwdy_uv, dudz_uv, dvdz_uv
#ifdef PPSCALARS
real(rprec) :: theta_p, pi_zp
#endif
real(rprec), allocatable, dimension(:,:,:) :: vortx, vorty, vortz
real(rprec) :: vortx_p, vorty_p, vortz_p

allocate(w_uv(nx,ny,lbz:nz), u_w(nx,ny,lbz:nz), v_w(nx,ny,lbz:nz))
allocate(pres_real(nx,ny,lbz:nz))
allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))

! Compute vorticity on w-grid from definition
!vortx(1:nx,1:ny,lbz:nz) = 0._rprec
!vorty(1:nx,1:ny,lbz:nz) = 0._rprec
!vortz(1:nx,1:ny,lbz:nz) = 0._rprec
vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) - dvdz(1:nx,1:ny,lbz:nz)
vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) - dwdx(1:nx,1:ny,lbz:nz)
vortz(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) - dudy(1:nx,1:ny,lbz:nz) !! on uv-grid, to be interpoalted soon

! Prepare variables for time-averaging
w_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz )
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(v(1:nx,1:ny,lbz:nz), lbz )

pres_real(1:nx,1:ny,lbz:nz) = 0._rprec
pres_real(1:nx,1:ny,lbz:nz) = p(1:nx,1:ny,lbz:nz)                              &
    - 0.5 * ( u(1:nx,1:ny,lbz:nz)**2 + w_uv(1:nx,1:ny,lbz:nz)**2               &
    + v(1:nx,1:ny,lbz:nz)**2 )

vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid( vortz(1:nx,1:ny,lbz:nz), lbz )

! note: u_w not necessarily zero on walls, but only mult by w=0 vu u'w', so OK
! can zero u_w at BC anyway:
if(coord==0       .and. lbc_mom>0) u_w(:,:,1)  = 0._rprec
if(coord==nproc-1 .and. ubc_mom>0) u_w(:,:,nz) = 0._rprec
if(coord==0       .and. lbc_mom>0) v_w(:,:,1)  = 0._rprec
if(coord==nproc-1 .and. ubc_mom>0) v_w(:,:,nz) = 0._rprec

! Begin time-averaging
do k = lbz, jzmax     ! lbz = 0 for mpi runs, otherwise lbz = 1
do j = 1, ny
do i = 1, nx
    u_p = u(i,j,k)       !! uv grid
    u_p2= u_w(i,j,k)     !! w grid
    v_p = v(i,j,k)       !! uv grid
    v_p2= v_w(i,j,k)     !! w grid
    w_p = w(i,j,k)       !! w grid
    w_p2= w_uv(i,j,k)    !! uv grid

    tavg(i,j,k) % u = tavg(i,j,k) % u + u_p * tavg_dt !! uv grid
    tavg(i,j,k) % v = tavg(i,j,k) % v + v_p * tavg_dt !! uv grid
    tavg(i,j,k) % w_uv = tavg(i,j,k) % w_uv + w_p2 * tavg_dt !! uv grid
    tavg(i,j,k) % u_w = tavg(i,j,k) % u_w + u_p2 * tavg_dt !! w grid
    tavg(i,j,k) % v_w = tavg(i,j,k) % v_w + v_p2 * tavg_dt !! w grid
    tavg(i,j,k) % w = tavg(i,j,k) % w + w_p * tavg_dt !! w grid

    ! Note: compute u'w' on w-grid because stresses on w-grid --pj
    tavg(i,j,k) % u2 = tavg(i,j,k) % u2 + u_p * u_p * tavg_dt !! uv grid
    tavg(i,j,k) % v2 = tavg(i,j,k) % v2 + v_p * v_p * tavg_dt !! uv grid
    tavg(i,j,k) % w2 = tavg(i,j,k) % w2 + w_p * w_p * tavg_dt !! w grid
    tavg(i,j,k) % uv = tavg(i,j,k) % uv + u_p * v_p * tavg_dt !! uv grid
    tavg(i,j,k) % uw = tavg(i,j,k) % uw + u_p2 * w_p * tavg_dt !! w grid
    tavg(i,j,k) % vw = tavg(i,j,k) % vw + v_p2 * w_p * tavg_dt !! w grid

    tavg(i,j,k) % txx = tavg(i,j,k) % txx + txx(i,j,k) * tavg_dt !! uv grid
    tavg(i,j,k) % tyy = tavg(i,j,k) % tyy + tyy(i,j,k) * tavg_dt !! uv grid
    tavg(i,j,k) % tzz = tavg(i,j,k) % tzz + tzz(i,j,k) * tavg_dt !! uv grid
    tavg(i,j,k) % txy = tavg(i,j,k) % txy + txy(i,j,k) * tavg_dt !! uv grid
    tavg(i,j,k) % txz = tavg(i,j,k) % txz + txz(i,j,k) * tavg_dt !! w grid
    tavg(i,j,k) % tyz = tavg(i,j,k) % tyz + tyz(i,j,k) * tavg_dt !! w grid

    tavg(i,j,k) % p = tavg(i,j,k) % p + pres_real(i,j,k) * tavg_dt !! uv grid

    ! All vorticity components are on the w grid
    tavg_vort(i,j,k)%vortx = tavg_vort(i,j,k)%vortx + vortx(i,j,k)*tavg_dt
    tavg_vort(i,j,k)%vorty = tavg_vort(i,j,k)%vorty + vorty(i,j,k)*tavg_dt
    tavg_vort(i,j,k)%vortz = tavg_vort(i,j,k)%vortz + vortz(i,j,k)*tavg_dt

    tavg_vort(i,j,k)%vortx2 = tavg_vort(i,j,k)%vortx2 +            &
        vortx(i,j,k)*vortx(i,j,k)*tavg_dt
    tavg_vort(i,j,k)%vorty2 = tavg_vort(i,j,k)%vorty2 +            &
        vorty(i,j,k)*vorty(i,j,k)*tavg_dt
    tavg_vort(i,j,k)%vortz2 = tavg_vort(i,j,k)%vortz2 +            &
        vortz(i,j,k)*vortz(i,j,k)*tavg_dt

#ifdef PPSCALARS
    theta_p = theta(i,j,k) !! uv grid
    pi_zp = pi_z(i,j,k) !! w grid

    tavg_scal(i,j,k) % theta = tavg_scal(i,j,k) % theta + theta_p * tavg_dt

    tavg_scal(i,j,k) % theta2 = tavg_scal(i,j,k) % theta2 + theta_p * theta_p * tavg_dt
    tavg_scal(i,j,k) % utheta = tavg_scal(i,j,k) % utheta + u_p * theta_p * tavg_dt
    tavg_scal(i,j,k) % vtheta = tavg_scal(i,j,k) % vtheta + v_p * theta_p * tavg_dt
    tavg_scal(i,j,k) % wtheta = tavg_scal(i,j,k) % wtheta + w_p2 * theta_p * tavg_dt
    tavg_scal(i,j,k) % pi_z = tavg_scal(i,j,k) % pi_z + pi_zp * tavg_dt
#endif
end do
end do
end do

#ifdef PPOUTPUT_SGS
do k = 1, jzmax
do j = 1, ny
do i = 1, nx
    tavg_sgs(i,j,k)%cs_opt2 = tavg_sgs(i,j,k)%cs_opt2 + Cs_opt2(i,j,k) * tavg_dt

    ! w-grid variables
    tavg_sgs(i,j,k)%Nu_t = tavg_sgs(i,j,k)%Nu_t + Nu_t(i,j,k) * tavg_dt
end do
end do
end do
#endif

#ifdef PPOUTPUT_BUDGET
call tavg_budget_compute
#endif

#ifdef PPOUTPUT_TURBSPEC
call tavg_turbspec_compute
#endif

#ifdef PPTURBINES
do k = 1, jzmax     ! lbz = 0 for mpi runs, otherwise lbz = 1
do j = 1, ny
do i = 1, nx
    tavg(i,j,k)%fx = tavg(i,j,k)%fx + fxa(i,j,k) * tavg_dt
    tavg(i,j,k)%fy = tavg(i,j,k)%fy + fya(i,j,k) * tavg_dt
    tavg(i,j,k)%fz = tavg(i,j,k)%fz + fza(i,j,k) * tavg_dt
end do
end do
end do
#endif

! Update tavg_total_time for variable time stepping
tavg_total_time = tavg_total_time + tavg_dt

! Set tavg_dt back to zero for next increment
tavg_dt = 0._rprec

end subroutine tavg_compute

#ifdef PPOUTPUT_BUDGET
!*****************************************************************************
subroutine tavg_budget_compute
!*****************************************************************************
!
! Computes terms that are used for budgets and time-averages accordingly.
! All quantities used in the averaging are moved to the w-grid.
!
! Has not been tested on MARCC, without MPI, or for full channel configuration.
!
! Assumes no-slip applies even for LES.
! 

use stat_defs, only : tavg_dt, tavg_budget
use param, only : nx, ny, nz, lbz, jzmax, lbc_mom, ubc_mom
use sim_param, only : u, v, w, p
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
use sim_param, only : dpdx, dpdy, dpdz
use sim_param, only : divtx, divty, divtz
use functions, only : interp_to_w_grid
use mpi_defs, only: mpi_sync_real_array, MPI_SYNC_DOWNUP

implicit none

integer :: i, j, k
real(rprec) :: u_p, v_p, w_p, p_p
real(rprec) :: dpdx_p, dpdy_p, dpdz_p
real(rprec) :: dudx_p, dudy_p, dudz_p, dvdx_p, dvdy_p, dvdz_p
real(rprec) :: dwdx_p, dwdy_p, dwdz_p, divtx_p, divty_p, divtz_p
real(rprec), allocatable, dimension(:,:,:) :: u_w, v_w, p_w
real(rprec), allocatable, dimension(:,:,:) :: dudx_w, dudy_w, dvdx_w, dvdy_w
real(rprec), allocatable, dimension(:,:,:) :: dwdz_w, divtx_w, divty_w
real(rprec), allocatable, dimension(:,:,:) :: pres_real
real(rprec), allocatable, dimension(:,:,:) :: dpdx_real, dpdy_real, dpdz_real

allocate(u_w(nx,ny,lbz:nz), v_w(nx,ny,lbz:nz), p_w(nx,ny,lbz:nz))
allocate(dudx_w(nx,ny,lbz:nz), dudy_w(nx,ny,lbz:nz), dvdx_w(nx,ny,lbz:nz))
allocate(dvdy_w(nx,ny,lbz:nz), dwdz_w(nx,ny,lbz:nz))
allocate(divtx_w(nx,ny,lbz:nz), divty_w(nx,ny,lbz:nz))
allocate(pres_real(nx,ny,lbz:nz))
allocate(dpdx_real(nx,ny,lbz:nz), dpdy_real(nx,ny,lbz:nz), dpdz_real(nx,ny,lbz:nz))

! Prepare variables that need to be interpolated onto the w-grid
#ifdef PPMPI
! Remove BOGUS values at processor interfaces
call mpi_sync_real_array( p, lbz, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( divtx, lbz, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( divty, lbz, MPI_SYNC_DOWNUP )

! Remove BOGUS value within boundary conditions as well
if (coord == 0) then 
    u(:,:,lbz) = 0._rprec
    v(:,:,lbz) = 0._rprec
    dwdz(:,:,lbz) = 0._rprec
    divtx(:,:,lbz) = 0._rprec
    divty(:,:,lbz) = 0._rprec
end if
if (coord == nproc-1) then
    dwdz(:,:,nz) = 0._rprec
    divtx(:,:,nz) = 0._rprec
    divty(:,:,nz) = 0._rprec
end if
#endif

! Now interpolate variables to the w-grid
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(v(1:nx,1:ny,lbz:nz), lbz )
p_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(p(1:nx,1:ny,lbz:nz), lbz)

dudx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dudx(1:nx,1:ny,lbz:nz), lbz )
dudy_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dudy(1:nx,1:ny,lbz:nz), lbz )
dvdx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dvdx(1:nx,1:ny,lbz:nz), lbz )
dvdy_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dvdy(1:nx,1:ny,lbz:nz), lbz )
dwdz_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dwdz(1:nx,1:ny,lbz:nz), lbz )

divtx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(divtx(1:nx,1:ny,lbz:nz), lbz)
divty_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(divty(1:nx,1:ny,lbz:nz), lbz)

! Remove energy from dynamic simulation pressure for static pressure
! This is different than in tavg_compute since on w-grid
pres_real(1:nx,1:ny,lbz:nz) = 0._rprec
pres_real(1:nx,1:ny,lbz:nz) = p_w(1:nx,1:ny,lbz:nz)                         &
    - 0.5 * ( u_w(1:nx,1:ny,lbz:nz)**2 + v_w(1:nx,1:ny,lbz:nz)**2           &
    + w(1:nx,1:ny,lbz:nz)**2 )

! dpdx and dpdy are treated differently since the 0 index is empty 
! and the nz index is BOGUS
! Initialize dpdx_real and dpdy_real from dpdx and dpdy
! BOGUS values of dpdx_real and dpdy_real are removed
! on the uv grid then brought over
! to the w grid
dpdx_real(1:nx,1:ny,1:nz) = dpdx(1:nx,1:ny,1:nz)
dpdy_real(1:nx,1:ny,1:nz) = dpdy(1:nx,1:ny,1:nz)

! Remove BOGUS value above ubc_mom
if (coord == nproc - 1) then
    dpdx_real(:,:,nz) = 0._rprec
    dpdy_real(:,:,nz) = 0._rprec
end if

#ifdef PPMPI
! Fill empty 0 index to be overwritten
dpdx_real(1:nx,1:ny,0) = 0._rprec
dpdy_real(1:nx,1:ny,0) = 0._rprec

! Remove intermediate BOGUS values (at nz) and zeros (at 0 index)
call mpi_sync_real_array( dpdx_real, lbz, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( dpdy_real, lbz, MPI_SYNC_DOWNUP )
#endif

! Now bring dpdx_real and dpdy_real actually onto the w grid
dpdx_real(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dpdx_real(1:nx,1:ny,lbz:nz), lbz)
dpdy_real(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dpdy_real(1:nx,1:ny,lbz:nz), lbz)

! Extract energy from pressure
dpdx_real(1:nx,1:ny,lbz:nz) = dpdx_real(1:nx,1:ny,lbz:nz)                   &
    - ( u_w(1:nx,1:ny,lbz:nz)*dudx_w(1:nx,1:ny,lbz:nz)                      &
    + v_w(1:nx,1:ny,lbz:nz)*dvdx_w(1:nx,1:ny,lbz:nz)                        &
    + w(1:nx,1:ny,lbz:nz)*dwdx(1:nx,1:ny,lbz:nz) )

dpdy_real(1:nx,1:ny,lbz:nz) = dpdy_real(1:nx,1:ny,lbz:nz)                   &
    - ( u_w(1:nx,1:ny,lbz:nz)*dudy_w(1:nx,1:ny,lbz:nz)                      &
    + v_w(1:nx,1:ny,lbz:nz)*dvdy_w(1:nx,1:ny,lbz:nz)                        &
    + w(1:nx,1:ny,lbz:nz)*dwdy(1:nx,1:ny,lbz:nz) )

! dpdx and dpdy were already interpolated, still need to consider dpdz_real
dpdz_real(1:nx,1:ny,1:nz) = dpdz(1:nx,1:ny,1:nz)

#ifdef PPMPI
dpdz_real(1:nx,1:ny,0) = dpdz(1:nx,1:ny,1)

! Remove intermediate BOGUS values (at nz) and zeros (at 0 index)
call mpi_sync_real_array( dpdz_real, lbz, MPI_SYNC_DOWNUP )
#endif

! Extract energy from pressure
dpdz_real(1:nx,1:ny,lbz:nz) = dpdz_real(1:nx,1:ny,lbz:nz)                   &
    - ( u_w(1:nx,1:ny,lbz:nz)*dudz(1:nx,1:ny,lbz:nz)                        &
    + v_w(1:nx,1:ny,lbz:nz)*dvdz(1:nx,1:ny,lbz:nz)                          &
    + w(1:nx,1:ny,lbz:nz)*dwdz_w(1:nx,1:ny,lbz:nz) )

! Enforce no penetration and no-slip
if (coord == 0 .and. lbc_mom > 0) then
    ! No-slip
    u_w(:,:,1) = 0._rprec
    v_w(:,:,1) = 0._rprec
    dudx_w(:,:,1) = 0._rprec
    dudy_w(:,:,1) = 0._rprec
    dvdx_w(:,:,1) = 0._rprec
    dvdy_w(:,:,1) = 0._rprec
    dwdz_w(:,:,1) = 0._rprec

    ! No penetration
    w(:,:,1) = 0._rprec
    dwdx(:,:,1) = 0._rprec
    dwdy(:,:,1) = 0._rprec
end if
if (coord == nproc-1 .and. ubc_mom > 0) then
    ! No-slip
    u_w(:,:,nz) = 0._rprec
    v_w(:,:,nz) = 0._rprec
    dudx_w(:,:,nz) = 0._rprec
    dudy_w(:,:,nz) = 0._rprec
    dvdx_w(:,:,nz) = 0._rprec
    dvdy_w(:,:,nz) = 0._rprec
    dwdz_w(:,:,nz) = 0._rprec

    ! No penetration
    w(:,:,nz) = 0._rprec
    dwdx(:,:,nz) = 0._rprec
    dwdy(:,:,nz) = 0._rprec
 end if

! Begin time-averaging
do k = lbz, jzmax     ! lbz = 0 for mpi runs, otherwise lbz = 1
do j = 1, ny
do i = 1, nx

    ! Rewrite temp variables
    u_p = u_w(i,j,k)
    v_p = v_w(i,j,k)
    w_p = w(i,j,k)
    p_p = pres_real(i,j,k)

    dudx_p = dudx_w(i,j,k)
    dudy_p = dudy_w(i,j,k)
    dudz_p = dudz(i,j,k)
    dvdx_p = dvdx_w(i,j,k)
    dvdy_p = dvdy_w(i,j,k)
    dvdz_p = dvdz(i,j,k)
    dwdx_p = dwdx(i,j,k)
    dwdy_p = dwdy(i,j,k)
    dwdz_p = dwdz_w(i,j,k)

    dpdx_p = dpdx_real(i,j,k)
    dpdy_p = dpdy_real(i,j,k)
    dpdz_p = dpdz_real(i,j,k)

    divtx_p = divtx_w(i,j,k)
    divty_p = divty_w(i,j,k)
    divtz_p = divtz(i,j,k)

    ! mean pressure on w-grid
    tavg_budget(i,j,k) % p = tavg_budget(i,j,k) % p + p_p * tavg_dt

    ! mean velocity-velocity product, ui*uj, on w-grid
    tavg_budget(i,j,k) % uu = tavg_budget(i,j,k) % uu + u_p * u_p * tavg_dt
    tavg_budget(i,j,k) % vv = tavg_budget(i,j,k) % vv + v_p * v_p * tavg_dt
    tavg_budget(i,j,k) % ww = tavg_budget(i,j,k) % ww + w_p * w_p * tavg_dt
    tavg_budget(i,j,k) % uv = tavg_budget(i,j,k) % uv + u_p * v_p * tavg_dt
    tavg_budget(i,j,k) % uw = tavg_budget(i,j,k) % uw + u_p * w_p * tavg_dt
    tavg_budget(i,j,k) % vw = tavg_budget(i,j,k) % vw + v_p * w_p * tavg_dt

    ! mean velocity gradients, duidxj
    tavg_budget(i,j,k) % dudx = tavg_budget(i,j,k) % dudx + dudx_p * tavg_dt
    tavg_budget(i,j,k) % dudy = tavg_budget(i,j,k) % dudy + dudy_p * tavg_dt 
    tavg_budget(i,j,k) % dudz = tavg_budget(i,j,k) % dudz + dudz_p * tavg_dt 
    tavg_budget(i,j,k) % dvdx = tavg_budget(i,j,k) % dvdx + dvdx_p * tavg_dt
    tavg_budget(i,j,k) % dvdy = tavg_budget(i,j,k) % dvdy + dvdy_p * tavg_dt
    tavg_budget(i,j,k) % dvdz = tavg_budget(i,j,k) % dvdz + dvdz_p * tavg_dt 
    tavg_budget(i,j,k) % dwdx = tavg_budget(i,j,k) % dwdx + dwdx_p * tavg_dt 
    tavg_budget(i,j,k) % dwdy = tavg_budget(i,j,k) % dwdy + dwdy_p * tavg_dt 
    tavg_budget(i,j,k) % dwdz = tavg_budget(i,j,k) % dwdz + dwdz_p * tavg_dt 

    ! mean pressure gradient, dpdxi
    tavg_budget(i,j,k) % dpdx = tavg_budget(i,j,k) % dpdx + dpdx_p * tavg_dt
    tavg_budget(i,j,k) % dpdy = tavg_budget(i,j,k) % dpdy + dpdy_p * tavg_dt
    tavg_budget(i,j,k) % dpdz = tavg_budget(i,j,k) % dpdz + dpdz_p * tavg_dt

    ! mean vel-velGrad product, ui*dujdxk
    tavg_budget(i,j,k)% ududx = tavg_budget(i,j,k) % ududx + u_p * dudx_p * tavg_dt
    tavg_budget(i,j,k)% ududy = tavg_budget(i,j,k) % ududy + u_p * dudy_p * tavg_dt
    tavg_budget(i,j,k)% ududz = tavg_budget(i,j,k) % ududz + u_p * dudz_p * tavg_dt
    tavg_budget(i,j,k)% udvdx = tavg_budget(i,j,k) % udvdx + u_p * dvdx_p * tavg_dt
    tavg_budget(i,j,k)% udvdy = tavg_budget(i,j,k) % udvdy + u_p * dvdy_p * tavg_dt
    tavg_budget(i,j,k)% udvdz = tavg_budget(i,j,k) % udvdz + u_p * dvdz_p * tavg_dt
    tavg_budget(i,j,k)% udwdx = tavg_budget(i,j,k) % udwdx + u_p * dwdx_p * tavg_dt
    tavg_budget(i,j,k)% udwdy = tavg_budget(i,j,k) % udwdy + u_p * dwdy_p * tavg_dt
    tavg_budget(i,j,k)% udwdz = tavg_budget(i,j,k) % udwdz + u_p * dwdz_p * tavg_dt

    tavg_budget(i,j,k)% vdudx = tavg_budget(i,j,k) % vdudx + v_p * dudx_p * tavg_dt
    tavg_budget(i,j,k)% vdudy = tavg_budget(i,j,k) % vdudy + v_p * dudy_p * tavg_dt
    tavg_budget(i,j,k)% vdudz = tavg_budget(i,j,k) % vdudz + v_p * dudz_p * tavg_dt
    tavg_budget(i,j,k)% vdvdx = tavg_budget(i,j,k) % vdvdx + v_p * dvdx_p * tavg_dt
    tavg_budget(i,j,k)% vdvdy = tavg_budget(i,j,k) % vdvdy + v_p * dvdy_p * tavg_dt
    tavg_budget(i,j,k)% vdvdz = tavg_budget(i,j,k) % vdvdz + v_p * dvdz_p * tavg_dt
    tavg_budget(i,j,k)% vdwdx = tavg_budget(i,j,k) % vdwdx + v_p * dwdx_p * tavg_dt
    tavg_budget(i,j,k)% vdwdy = tavg_budget(i,j,k) % vdwdy + v_p * dwdy_p * tavg_dt
    tavg_budget(i,j,k)% vdwdz = tavg_budget(i,j,k) % vdwdz + v_p * dwdz_p * tavg_dt

    tavg_budget(i,j,k)% wdudx = tavg_budget(i,j,k) % wdudx + w_p * dudx_p * tavg_dt
    tavg_budget(i,j,k)% wdudy = tavg_budget(i,j,k) % wdudy + w_p * dudy_p * tavg_dt
    tavg_budget(i,j,k)% wdudz = tavg_budget(i,j,k) % wdudz + w_p * dudz_p * tavg_dt
    tavg_budget(i,j,k)% wdvdx = tavg_budget(i,j,k) % wdvdx + w_p * dvdx_p * tavg_dt
    tavg_budget(i,j,k)% wdvdy = tavg_budget(i,j,k) % wdvdy + w_p * dvdy_p * tavg_dt
    tavg_budget(i,j,k)% wdvdz = tavg_budget(i,j,k) % wdvdz + w_p * dvdz_p * tavg_dt
    tavg_budget(i,j,k)% wdwdx = tavg_budget(i,j,k) % wdwdx + w_p * dwdx_p * tavg_dt
    tavg_budget(i,j,k)% wdwdy = tavg_budget(i,j,k) % wdwdy + w_p * dwdy_p * tavg_dt
    tavg_budget(i,j,k)% wdwdz = tavg_budget(i,j,k) % wdwdz + w_p * dwdz_p * tavg_dt

    ! Mean vel-vel-velGrad product, ui*uk*dujdxk
    tavg_budget(i,j,k)%uududx = tavg_budget(i,j,k)%uududx+ u_p*u_p*dudx_p * tavg_dt
    tavg_budget(i,j,k)%uvdudy = tavg_budget(i,j,k)%uvdudy+ u_p*v_p*dudy_p * tavg_dt
    tavg_budget(i,j,k)%uwdudz = tavg_budget(i,j,k)%uwdudz+ u_p*w_p*dudz_p * tavg_dt
    tavg_budget(i,j,k)%uudvdx = tavg_budget(i,j,k)%uudvdx+ u_p*u_p*dvdx_p * tavg_dt
    tavg_budget(i,j,k)%uvdvdy = tavg_budget(i,j,k)%uvdvdy+ u_p*v_p*dvdy_p * tavg_dt
    tavg_budget(i,j,k)%uwdvdz = tavg_budget(i,j,k)%uwdvdz+ u_p*w_p*dvdz_p * tavg_dt
    tavg_budget(i,j,k)%uudwdx = tavg_budget(i,j,k)%uudwdx+ u_p*u_p*dwdx_p * tavg_dt
    tavg_budget(i,j,k)%uvdwdy = tavg_budget(i,j,k)%uvdwdy+ u_p*v_p*dwdy_p * tavg_dt
    tavg_budget(i,j,k)%uwdwdz = tavg_budget(i,j,k)%uwdwdz+ u_p*w_p*dwdz_p * tavg_dt

    tavg_budget(i,j,k)%vududx=tavg_budget(i,j,k)%vududx + v_p*u_p*dudx_p * tavg_dt
    tavg_budget(i,j,k)%vvdudy=tavg_budget(i,j,k)%vvdudy + v_p*v_p*dudy_p * tavg_dt
    tavg_budget(i,j,k)%vwdudz=tavg_budget(i,j,k)%vwdudz + v_p*w_p*dudz_p * tavg_dt
    tavg_budget(i,j,k)%vudvdx=tavg_budget(i,j,k)%vudvdx + v_p*u_p*dvdx_p * tavg_dt
    tavg_budget(i,j,k)%vvdvdy=tavg_budget(i,j,k)%vvdvdy + v_p*v_p*dvdy_p * tavg_dt
    tavg_budget(i,j,k)%vwdvdz=tavg_budget(i,j,k)%vwdvdz + v_p*w_p*dvdz_p * tavg_dt
    tavg_budget(i,j,k)%vudwdx=tavg_budget(i,j,k)%vudwdx + v_p*u_p*dwdx_p * tavg_dt
    tavg_budget(i,j,k)%vvdwdy=tavg_budget(i,j,k)%vvdwdy + v_p*v_p*dwdy_p * tavg_dt
    tavg_budget(i,j,k)%vwdwdz=tavg_budget(i,j,k)%vwdwdz + v_p*w_p*dwdz_p * tavg_dt

    tavg_budget(i,j,k)%wududx=tavg_budget(i,j,k)%wududx + w_p*u_p*dudx_p * tavg_dt
    tavg_budget(i,j,k)%wvdudy=tavg_budget(i,j,k)%wvdudy + w_p*v_p*dudy_p * tavg_dt
    tavg_budget(i,j,k)%wwdudz=tavg_budget(i,j,k)%wwdudz + w_p*w_p*dudz_p * tavg_dt
    tavg_budget(i,j,k)%wudvdx=tavg_budget(i,j,k)%wudvdx + w_p*u_p*dvdx_p * tavg_dt
    tavg_budget(i,j,k)%wvdvdy=tavg_budget(i,j,k)%wvdvdy + w_p*v_p*dvdy_p * tavg_dt
    tavg_budget(i,j,k)%wwdvdz=tavg_budget(i,j,k)%wwdvdz + w_p*w_p*dvdz_p * tavg_dt
    tavg_budget(i,j,k)%wudwdx=tavg_budget(i,j,k)%wudwdx + w_p*u_p*dwdx_p * tavg_dt
    tavg_budget(i,j,k)%wvdwdy=tavg_budget(i,j,k)%wvdwdy + w_p*v_p*dwdy_p * tavg_dt
    tavg_budget(i,j,k)%wwdwdz=tavg_budget(i,j,k)%wwdwdz + w_p*w_p*dwdz_p * tavg_dt

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
    tavg_budget(i,j,k)%uxux = tavg_budget(i,j,k)%uxux + dudx_p * dudx_p * tavg_dt
    tavg_budget(i,j,k)%uyuy = tavg_budget(i,j,k)%uyuy + dudy_p * dudy_p * tavg_dt
    tavg_budget(i,j,k)%uzuz = tavg_budget(i,j,k)%uzuz + dudz_p * dudz_p * tavg_dt
    tavg_budget(i,j,k)%vxvx = tavg_budget(i,j,k)%vxvx + dvdx_p * dvdx_p * tavg_dt
    tavg_budget(i,j,k)%vyvy = tavg_budget(i,j,k)%vyvy + dvdy_p * dvdy_p * tavg_dt
    tavg_budget(i,j,k)%vzvz = tavg_budget(i,j,k)%vzvz + dvdz_p * dvdz_p * tavg_dt
    tavg_budget(i,j,k)%wxwx = tavg_budget(i,j,k)%wxwx + dwdx_p * dwdx_p * tavg_dt
    tavg_budget(i,j,k)%wywy = tavg_budget(i,j,k)%wywy + dwdy_p * dwdy_p * tavg_dt
    tavg_budget(i,j,k)%wzwz = tavg_budget(i,j,k)%wzwz + dwdz_p * dwdz_p * tavg_dt

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
    tavg_budget(i,j,k)%uxvx = tavg_budget(i,j,k)%uxvx + dudx_p * dvdx_p * tavg_dt
    tavg_budget(i,j,k)%uyvy = tavg_budget(i,j,k)%uyvy + dudy_p * dvdy_p * tavg_dt
    tavg_budget(i,j,k)%uzvz = tavg_budget(i,j,k)%uzvz + dudz_p * dvdz_p * tavg_dt
    tavg_budget(i,j,k)%uxwx = tavg_budget(i,j,k)%uxwx + dudx_p * dwdx_p * tavg_dt
    tavg_budget(i,j,k)%uywy = tavg_budget(i,j,k)%uywy + dudy_p * dwdy_p * tavg_dt
    tavg_budget(i,j,k)%uzwz = tavg_budget(i,j,k)%uzwz + dudz_p * dwdz_p * tavg_dt
    tavg_budget(i,j,k)%vxwx = tavg_budget(i,j,k)%vxwx + dvdx_p * dwdx_p * tavg_dt
    tavg_budget(i,j,k)%vywy = tavg_budget(i,j,k)%vywy + dvdy_p * dwdy_p * tavg_dt
    tavg_budget(i,j,k)%vzwz = tavg_budget(i,j,k)%vzwz + dvdz_p * dwdz_p * tavg_dt

    ! another velocity gradient-velocity gradient correlation, duixj*dUjdxi
    !tavg_budget(i,j,k)%uyvx = tavg_budget(i,j,k)%uyvx + dudy_p * dvdx_p * tavg_dt
    !tavg_budget(i,j,k)%uzwx = tavg_budget(i,j,k)%uzwx + dudz_p * dwdx_p * tavg_dt
    !tavg_budget(i,j,k)%vzwy = tavg_budget(i,j,k)%vzwy + dvdz_p * dwdy_p * tavg_dt

    ! Mean vel-presGrad product, ui*dpdxj
    tavg_budget(i,j,k)%udpdx = tavg_budget(i,j,k) % udpdx + u_p * dpdx_p * tavg_dt
    tavg_budget(i,j,k)%udpdy = tavg_budget(i,j,k) % udpdy + u_p * dpdy_p * tavg_dt
    tavg_budget(i,j,k)%udpdz = tavg_budget(i,j,k) % udpdz + u_p * dpdz_p * tavg_dt
    tavg_budget(i,j,k)%vdpdx = tavg_budget(i,j,k) % vdpdx + v_p * dpdx_p * tavg_dt
    tavg_budget(i,j,k)%vdpdy = tavg_budget(i,j,k) % vdpdy + v_p * dpdy_p * tavg_dt
    tavg_budget(i,j,k)%vdpdz = tavg_budget(i,j,k) % vdpdz + v_p * dpdz_p * tavg_dt
    tavg_budget(i,j,k)%wdpdx = tavg_budget(i,j,k) % wdpdx + w_p * dpdx_p * tavg_dt
    tavg_budget(i,j,k)%wdpdy = tavg_budget(i,j,k) % wdpdy + w_p * dpdy_p * tavg_dt
    tavg_budget(i,j,k)%wdpdz = tavg_budget(i,j,k) % wdpdz + w_p * dpdz_p * tavg_dt

    ! Mean pres-velGrad product, p*duidxj
    tavg_budget(i,j,k)%pdudx = tavg_budget(i,j,k) % pdudx + p_p * dudx_p * tavg_dt
    tavg_budget(i,j,k)%pdudy = tavg_budget(i,j,k) % pdudy + p_p * dudy_p * tavg_dt
    tavg_budget(i,j,k)%pdudz = tavg_budget(i,j,k) % pdudz + p_p * dudz_p * tavg_dt
    tavg_budget(i,j,k)%pdvdx = tavg_budget(i,j,k) % pdvdx + p_p * dvdx_p * tavg_dt
    tavg_budget(i,j,k)%pdvdy = tavg_budget(i,j,k) % pdvdy + p_p * dvdy_p * tavg_dt
    tavg_budget(i,j,k)%pdvdz = tavg_budget(i,j,k) % pdvdz + p_p * dvdz_p * tavg_dt
    tavg_budget(i,j,k)%pdwdx = tavg_budget(i,j,k) % pdwdx + p_p * dwdx_p * tavg_dt
    tavg_budget(i,j,k)%pdwdy = tavg_budget(i,j,k) % pdwdy + p_p * dwdy_p * tavg_dt
    tavg_budget(i,j,k)%pdwdz = tavg_budget(i,j,k) % pdwdz + p_p * dwdz_p * tavg_dt

    ! Mean Laplacian, nu*lap(uj)
    tavg_budget(i,j,k) % lapu = tavg_budget(i,j,k) % lapu + divtx_p * tavg_dt
    tavg_budget(i,j,k) % lapv = tavg_budget(i,j,k) % lapv + divty_p * tavg_dt
    tavg_budget(i,j,k) % lapw = tavg_budget(i,j,k) % lapw + divtz_p * tavg_dt

    ! Mean Vel-Laplacian, nu*ui*lap(uj)
    tavg_budget(i,j,k)%ulapu = tavg_budget(i,j,k)%ulapu + u_p * divtx_p * tavg_dt
    tavg_budget(i,j,k)%ulapv = tavg_budget(i,j,k)%ulapv + u_p * divty_p * tavg_dt
    tavg_budget(i,j,k)%ulapw = tavg_budget(i,j,k)%ulapw + u_p * divtz_p * tavg_dt
    tavg_budget(i,j,k)%vlapu = tavg_budget(i,j,k)%vlapu + v_p * divtx_p * tavg_dt
    tavg_budget(i,j,k)%vlapv = tavg_budget(i,j,k)%vlapv + v_p * divty_p * tavg_dt
    tavg_budget(i,j,k)%vlapw = tavg_budget(i,j,k)%vlapw + v_p * divtz_p * tavg_dt
    tavg_budget(i,j,k)%wlapu = tavg_budget(i,j,k)%wlapu + w_p * divtx_p * tavg_dt
    tavg_budget(i,j,k)%wlapv = tavg_budget(i,j,k)%wlapv + w_p * divty_p * tavg_dt
    tavg_budget(i,j,k)%wlapw = tavg_budget(i,j,k)%wlapw + w_p * divtz_p * tavg_dt

end do
end do
end do

end subroutine tavg_budget_compute
#endif

#ifdef PPOUTPUT_TURBSPEC
!*****************************************************************************
subroutine tavg_turbspec_compute
!*****************************************************************************
!
! Computes various spectral densities, including energy (velocity) and 
! dissipation (vorticity), then time-averages these quantities.
! 
! Additional quantities for the spectral budget are also computed here.
! 
use types, only: rprec 
use param, only: nx, ny, nz, lbz, coord, ubc_mom, lbc_mom
use sim_param, only: u, v, w, dwdy, dvdz, dudz, dwdx, dvdx, dudy
use stat_defs, only: tavg_dt, tavg_turbspecx, tavg_turbspecy
use fft
use functions, only : interp_to_w_grid
#ifdef PPSCALARS
use scalars, only: theta
use stat_defs, only: tavg_scal_turbspecx, tavg_scal_turbspecy
#ifdef PPOUTPUT_SPECBUDG
use scalars, only: dTdx, dTdy, dTdz, div_pi
use stat_defs, only: tavg_scal_specbudgx
#endif
#endif
#ifdef PPOUTPUT_SPECBUDG
use stat_defs, only: tavg_specbudgx, tavg_specbudgy
use sim_param, only: p, dpdx, dpdy, dpdz, dudx, dvdy, dwdz
use sim_param, only: divtx, divty, divtz
use mpi_defs, only: mpi_sync_real_array, MPI_SYNC_DOWNUP
#endif

implicit none

integer :: jx, jy, jz
real(rprec), dimension(ld,ny,lbz:nz) :: u_w, v_w
real(rprec), dimension(ld,ny,lbz:nz) :: vortx, vorty, vortz
real(rprec), dimension(nx) :: upx, vpx, wpx
real(rprec), dimension(nx) :: vortxpx, vortypx, vortzpx
complex(rprec), dimension(nx/2+1) :: ufx, vfx, wfx
complex(rprec), dimension(nx/2+1) :: vortxfx, vortyfx, vortzfx
real(rprec), dimension(ny) :: upy, vpy, wpy
real(rprec), dimension(ny) :: vortxpy, vortypy, vortzpy
complex(rprec), dimension(ny/2+1) :: ufy, vfy, wfy
complex(rprec), dimension(ny/2+1) :: vortxfy, vortyfy, vortzfy
real(rprec) :: const1, const2
#ifdef PPOUTPUT_SPECBUDG
real(rprec), dimension(ld,ny,lbz:nz) :: p_w, pres_real
real(rprec), dimension(ld,ny,lbz:nz) :: dpdx_real, dpdy_real, dpdz_real
real(rprec), dimension(ld,ny,lbz:nz) :: dudx_w, dudy_w, dvdx_w, dvdy_w, dwdz_w
real(rprec), dimension(ld,ny,lbz:nz) :: divtx_w, divty_w
real(rprec), dimension(nx) :: dudxpx, dudypx, dudzpx,                     &
    dvdxpx, dvdypx, dvdzpx, dwdxpx, dwdypx, dwdzpx,                       &
    ppx, dpdxpx, dpdypx, dpdzpx, divtxpx, divtypx, divtzpx,               &
    ududxpx, vdudypx, wdudzpx, udvdxpx, vdvdypx, wdvdzpx,                 &
    udwdxpx, vdwdypx, wdwdzpx
complex(rprec), dimension(nx/2+1) :: dudxfx, dudyfx, dudzfx,              &
    dvdxfx, dvdyfx, dvdzfx, dwdxfx, dwdyfx, dwdzfx,                       &
    pfx, dpdxfx, dpdyfx, dpdzfx, divtxfx, divtyfx, divtzfx,               &
    ududxfx, vdudyfx, wdudzfx, udvdxfx, vdvdyfx, wdvdzfx,                 &
    udwdxfx, vdwdyfx, wdwdzfx
real(rprec), dimension(ny) :: dudxpy, dudypy, dudzpy,                     &
    dvdxpy, dvdypy, dvdzpy, dwdxpy, dwdypy, dwdzpy,                       &
    ppy, dpdxpy, dpdypy, dpdzpy, divtxpy, divtypy, divtzpy,               &
    ududxpy, vdudypy, wdudzpy, udvdxpy, vdvdypy, wdvdzpy,                 &
    udwdxpy, vdwdypy, wdwdzpy
complex(rprec), dimension(ny/2+1) :: dudxfy, dudyfy, dudzfy,              &
    dvdxfy, dvdyfy, dvdzfy, dwdxfy, dwdyfy, dwdzfy,                       &
    pfy, dpdxfy, dpdyfy, dpdzfy, divtxfy, divtyfy, divtzfy,               &
    ududxfy, vdudyfy, wdudzfy, udvdxfy, vdvdyfy, wdvdzfy,                 &
    udwdxfy, vdwdyfy, wdwdzfy
#endif

!real(rprec), dimension(ld,ny,lbz:nz) :: vel, vort
!real(rprec), dimension(nx) :: velpx, vortpx
!complex(rprec), dimension(nx/2+1) :: velfx, vortfx
!real(rprec), dimension(nx) :: velpy, vortpy
!complex(rprec), dimension(ny/2+1) :: velfy, vortfy

#ifdef PPSCALARS
real(rprec), dimension(ld,ny,lbz:nz) :: theta_w
real(rprec), dimension(nx) :: thetapx
complex(rprec), dimension(nx/2+1) :: thetafx
real(rprec), dimension(ny) :: thetapy
complex(rprec), dimension(ny/2+1) :: thetafy
#ifdef PPOUTPUT_SPECBUDG
real(rprec), dimension(ld,ny,lbz:nz) :: dTdx_w, dTdy_w
real(rprec), dimension(nx) :: dTdxpx, dTdypx, dTdzpx
complex(rprec), dimension(nx/2+1) :: dTdxfx, dTdyfx, dTdzfx
real(rprec), dimension(nx) :: udTdxpx, vdTdypx, wdTdzpx
complex(rprec), dimension(nx/2+1) :: udTdxfx, vdTdyfx, wdTdzfx
real(rprec), dimension(ld,ny,lbz:nz) :: lapT_w
real(rprec), dimension(nx) :: lapTpx
complex(rprec), dimension(nx/2+1) :: lapTfx
#endif
#endif

! Set normalization constants
const1 = 1._rprec / nx
const2 = 1._rprec / ny

#ifdef PPOUTPUT_SPECBUDG
! Prepare variables that need to be interpolated onto the w-grid
! Remove BOGUS values at processor interfaces
call mpi_sync_real_array( p, lbz, MPI_SYNC_DOWNUP )

! Remove BOGUS value within boundary conditions as well
if (coord == 0) then
    u(:,:,lbz) = 0._rprec
    v(:,:,lbz) = 0._rprec
    dwdz(:,:,lbz) = 0._rprec
endif
if (coord == nproc-1) then
    dwdz(:,:,nz) = 0._rprec
endif
#endif

! Compute vorticity on w-grid from definition
vortx(:,:,:) = dwdy(:,:,:) - dvdz(:,:,:)
vorty(:,:,:) = dudz(:,:,:) - dwdx(:,:,:)
vortz(:,:,:) = dvdx(:,:,:) - dudy(:,:,:)

! Interpolate quantities to w-grid
u_w(:,:,:) = u(:,:,:) !! carry over ld-1:ld values
v_w(:,:,:) = v(:,:,:)
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(u_w(1:nx,1:ny,lbz:nz), lbz)
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(v_w(1:nx,1:ny,lbz:nz), lbz)
vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid(vortz(1:nx,1:ny,lbz:nz), lbz)
! dvdx and dudy both on uv-grid, just using vortz, so only interpolate once
! rewriting vortz variable

#ifdef PPOUTPUT_SPECBUDG
dudx_w(:,:,:) = dudx(:,:,:)
dudy_w(:,:,:) = dudy(:,:,:)
dvdx_w(:,:,:) = dvdx(:,:,:)
dvdy_w(:,:,:) = dvdy(:,:,:)
dwdz_w(:,:,:) = dwdz(:,:,:)
dudx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dudx_w(1:nx,1:ny,lbz:nz), lbz)
dudy_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dudy_w(1:nx,1:ny,lbz:nz), lbz)
dvdx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dvdx_w(1:nx,1:ny,lbz:nz), lbz)
dvdy_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dvdy_w(1:nx,1:ny,lbz:nz), lbz)
dwdz_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dwdz_w(1:nx,1:ny,lbz:nz), lbz)
p_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(p(1:nx,1:ny,lbz:nz), lbz)
divtx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(divtx(1:nx,1:ny,lbz:nz), lbz)
divty_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(divty(1:nx,1:ny,lbz:nz), lbz)

! Enforce no penetration and no-slip
if (coord == 0 .and. lbc_mom > 0) then
    ! No-slip
    u_w(:,:,1) = 0._rprec
    v_w(:,:,1) = 0._rprec
    dudx_w(:,:,1) = 0._rprec
    dudy_w(:,:,1) = 0._rprec
    dvdx_w(:,:,1) = 0._rprec
    dvdy_w(:,:,1) = 0._rprec
    dwdz_w(:,:,1) = 0._rprec

    ! No penetration
    w(:,:,1) = 0._rprec
    dwdx(:,:,1) = 0._rprec
    dwdy(:,:,1) = 0._rprec
end if
if (coord == nproc-1 .and. ubc_mom > 0) then
    ! No-slip
    u_w(:,:,nz) = 0._rprec
    v_w(:,:,nz) = 0._rprec
    dudx_w(:,:,nz) = 0._rprec
    dudy_w(:,:,nz) = 0._rprec
    dvdx_w(:,:,nz) = 0._rprec
    dvdy_w(:,:,nz) = 0._rprec
    dwdz_w(:,:,nz) = 0._rprec

    ! No penetration
    w(:,:,nz) = 0._rprec
    dwdx(:,:,nz) = 0._rprec
    dwdy(:,:,nz) = 0._rprec
 end if

! Remove energy from dynamic simulation pressure for static pressure
pres_real(:,:,:) = 0._rprec
pres_real(1:nx,1:ny,lbz:nz) = p_w(1:nx,1:ny,lbz:nz)                   &
    - 0.5 * ( u_w(1:nx,1:ny,lbz:nz)**2 + v_w(1:nx,1:ny,lbz:nz)**2     &
    + w(1:nx,1:ny,lbz:nz)**2 )

! dpdx and dpdy are treated differently since the 0 index is empty
! and the nz index is BOGUS
! Initialize dpdx_real and dpdy_real from dpdx and dpdy
! BOGUS values of dpdx_real and dpy_real are removed
! on the uv grid then brought over to the w grid
dpdx_real(:,:,1:nz) = dpdx(:,:,1:nz)
dpdy_real(:,:,1:nz) = dpdy(:,:,1:nz)

! Remove BOGUS value above ubc_mom
if (coord == nproc-1) then
    dpdx_real(:,:,nz) = 0._rprec
    dpdy_real(:,:,nz) = 0._rprec
endif

! Fill empty 0 index to be overwritten
dpdx_real(1:nx,1:ny,0) = 0._rprec
dpdy_real(1:nx,1:ny,0) = 0._rprec

! Remove intermediate BOGUS values (at nz) and zeros (at 0 index)
call mpi_sync_real_array( dpdx_real, lbz, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( dpdy_real, lbz, MPI_SYNC_DOWNUP )

! Now bring dpdx_real and dpdy_real onto the w grid
dpdx_real(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dpdx_real(1:nx,1:ny,lbz:nz), lbz)
dpdy_real(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dpdy_real(1:nx,1:ny,lbz:nz), lbz)

! Extract energy from pressure
dpdx_real(1:nx,1:ny,lbz:nz) = dpdx_real(1:nx,1:ny,lbz:nz)                   &
    - ( u_w(1:nx,1:ny,lbz:nz)*dudx_w(1:nx,1:ny,lbz:nz)                      &
    + v_w(1:nx,1:ny,lbz:nz)*dvdx_w(1:nx,1:ny,lbz:nz)                        &
    + w(1:nx,1:ny,lbz:nz)*dwdx(1:nx,1:ny,lbz:nz) )

dpdy_real(1:nx,1:ny,lbz:nz) = dpdy_real(1:nx,1:ny,lbz:nz)                   &
    - ( u_w(1:nx,1:ny,lbz:nz)*dudy_w(1:nx,1:ny,lbz:nz)                      &
    + v_w(1:nx,1:ny,lbz:nz)*dvdy_w(1:nx,1:ny,lbz:nz)                        &
    + w(1:nx,1:ny,lbz:nz)*dwdy(1:nx,1:ny,lbz:nz) )

! dpdx and dpdy were already interpolated, still need to consider dpdz_real
dpdz_real(1:nx,1:ny,1:nz) = dpdz(1:nx,1:ny,1:nz)

#ifdef PPMPI
dpdz_real(1:nx,1:ny,0) = dpdz(1:nx,1:ny,1)

! Remove intermediate BOGUS values (at nz) and zeros (at 0 index)
call mpi_sync_real_array( dpdz_real, lbz, MPI_SYNC_DOWNUP )
#endif

! Extract energy from pressure
dpdz_real(1:nx,1:ny,lbz:nz) = dpdz_real(1:nx,1:ny,lbz:nz)                   &
    - ( u_w(1:nx,1:ny,lbz:nz)*dudz(1:nx,1:ny,lbz:nz)                        &
    + v_w(1:nx,1:ny,lbz:nz)*dvdz(1:nx,1:ny,lbz:nz)                          &
    + w(1:nx,1:ny,lbz:nz)*dwdz_w(1:nx,1:ny,lbz:nz) )
#endif

! Compute magnitudes
! vel(:,:,:) = sqrt( u_w(:,:,:)**2 + v_w(:,:,:)**2 + w(:,:,:)**2 )
! vort(:,:,:) = sqrt( vortx(:,:,:)**2 + vorty(:,:,:)**2 + vortz(:,:,:)**2 )

#ifdef PPSCALARS
theta_w(:,:,:) = theta(:,:,:)
theta_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(theta_w(1:nx,1:ny,lbz:nz), lbz)
dTdx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dTdx(1:nx,1:ny,lbz:nz), lbz)
dTdy_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(dTdy(1:nx,1:ny,lbz:nz), lbz)
lapT_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(div_pi(1:nx,1:ny,lbz:nz), lbz)
#endif

do jy = 1, ny
do jz = 1, nz

    ! Normalize
    upx = const1 * u_w(1:nx,jy,jz)
    vpx = const1 * v_w(1:nx,jy,jz)
    wpx = const1 * w(1:nx,jy,jz)
    !velpx = const1 * vel(1:nx,jy,jz)
    vortxpx = const1 * vortx(1:nx,jy,jz)
    vortypx = const1 * vorty(1:nx,jy,jz)
    vortzpx = const1 * vortz(1:nx,jy,jz)
    !vortpx = const1 * vort(1:nx,jy,jz)
#ifdef PPSCALARS
    thetapx = const1 * theta_w(1:nx,jy,jz)
#ifdef PPOUTPUT_SPECBUDG
    dTdxpx  = const1 * dTdx_w(1:nx,jy,jz)
    dTdypx  = const1 * dTdy_w(1:nx,jy,jz)
    dTdzpx  = const1 * dTdz(1:nx,jy,jz)
    lapTpx  = const1 * lapT_w(1:nx,jy,jz)
    udTdxpx = const1 * u_w(1:nx,jy,jz)*dTdx_w(1:nx,jy,jz)
    vdTdypx = const1 * v_w(1:nx,jy,jz)*dTdy_w(1:nx,jy,jz)
    wdTdzpx = const1 * w(1:nx,jy,jz)*dTdz(1:nx,jy,jz)
#endif
#endif

    ! Take 1d transform, x --> kx
    call dfftw_execute_dft_r2c( forw_x, upx,     ufx )
    call dfftw_execute_dft_r2c( forw_x, vpx,     vfx )
    call dfftw_execute_dft_r2c( forw_x, wpx,     wfx )
    !call dfftw_execute_dft_r2c( forw_x, velpx,   velfx )
    call dfftw_execute_dft_r2c( forw_x, vortxpx, vortxfx )
    call dfftw_execute_dft_r2c( forw_x, vortypx, vortyfx )
    call dfftw_execute_dft_r2c( forw_x, vortzpx, vortzfx )
    !call dfftw_execute_dft_r2c( forw_x, vortpx,  vortfx )
#ifdef PPSCALARS
    call dfftw_execute_dft_r2c( forw_x, thetapx, thetafx )
#ifdef PPOUTPUT_SPECBUDG
    call dfftw_execute_dft_r2c( forw_x, dTdxpx, dTdxfx )
    call dfftw_execute_dft_r2c( forw_x, dTdypx, dTdyfx )
    call dfftw_execute_dft_r2c( forw_x, dTdzpx, dTdzfx )
    call dfftw_execute_dft_r2c( forw_x, lapTpx, lapTfx )
    call dfftw_execute_dft_r2c( forw_x, udTdxpx, udTdxfx )
    call dfftw_execute_dft_r2c( forw_x, vdTdypx, vdTdyfx )
    call dfftw_execute_dft_r2c( forw_x, wdTdzpx, wdTdzfx )
#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
    ! Normalize
    ppx = const1 * pres_real(1:nx,jy,jz)
    dpdxpx = const1 * dpdx_real(1:nx,jy,jz)
    dpdypx = const1 * dpdx_real(1:nx,jy,jz)
    dpdzpx = const1 * dpdx_real(1:nx,jy,jz)
    dudxpx = const1 * dudx_w(1:nx,jy,jz)
    dudypx = const1 * dudy_w(1:nx,jy,jz)
    dudzpx = const1 * dudz(1:nx,jy,jz)
    dvdxpx = const1 * dvdx_w(1:nx,jy,jz)
    dvdypx = const1 * dvdy_w(1:nx,jy,jz)
    dvdzpx = const1 * dvdz(1:nx,jy,jz)
    dwdxpx = const1 * dwdx(1:nx,jy,jz)
    dwdypx = const1 * dwdy(1:nx,jy,jz)
    dwdzpx = const1 * dwdz_w(1:nx,jy,jz)
    divtxpx = const1 * divtx_w(1:nx,jy,jz)
    divtypx = const1 * divty_w(1:nx,jy,jz)
    divtzpx = const1 * divtz(1:nx,jy,jz)
    ududxpx = const1 * u_w(1:nx,jy,jz)*dudx_w(1:nx,jy,jz)
    vdudypx = const1 * v_w(1:nx,jy,jz)*dudy_w(1:nx,jy,jz)
    wdudzpx = const1 * w(1:nx,jy,jz)*dudz(1:nx,jy,jz)
    udvdxpx = const1 * u_w(1:nx,jy,jz)*dvdx_w(1:nx,jy,jz)
    vdvdypx = const1 * v_w(1:nx,jy,jz)*dvdy_w(1:nx,jy,jz)
    wdvdzpx = const1 * w(1:nx,jy,jz)*dvdz(1:nx,jy,jz)
    udwdxpx = const1 * u_w(1:nx,jy,jz)*dwdx(1:nx,jy,jz)
    vdwdypx = const1 * v_w(1:nx,jy,jz)*dwdy(1:nx,jy,jz)
    wdwdzpx = const1 * w(1:nx,jy,jz)*dwdz_w(1:nx,jy,jz)

    ! Take 1d transform, x --> kx
    call dfftw_execute_dft_r2c( forw_x, ppx, pfx )
    call dfftw_execute_dft_r2c( forw_x, dpdxpx, dpdxfx )
    call dfftw_execute_dft_r2c( forw_x, dpdypx, dpdyfx )
    call dfftw_execute_dft_r2c( forw_x, dpdzpx, dpdzfx )
    call dfftw_execute_dft_r2c( forw_x, dudxpx, dudxfx )
    call dfftw_execute_dft_r2c( forw_x, dudypx, dudyfx )
    call dfftw_execute_dft_r2c( forw_x, dudzpx, dudzfx )
    call dfftw_execute_dft_r2c( forw_x, dvdxpx, dvdxfx )
    call dfftw_execute_dft_r2c( forw_x, dvdypx, dvdyfx )
    call dfftw_execute_dft_r2c( forw_x, dvdzpx, dvdzfx )
    call dfftw_execute_dft_r2c( forw_x, dwdxpx, dwdxfx )
    call dfftw_execute_dft_r2c( forw_x, dwdypx, dwdyfx )
    call dfftw_execute_dft_r2c( forw_x, dwdzpx, dwdzfx )
    call dfftw_execute_dft_r2c( forw_x, divtxpx, divtxfx )
    call dfftw_execute_dft_r2c( forw_x, divtypx, divtyfx )
    call dfftw_execute_dft_r2c( forw_x, divtzpx, divtzfx )
    call dfftw_execute_dft_r2c( forw_x, ududxpx, ududxfx )
    call dfftw_execute_dft_r2c( forw_x, vdudypx, vdudyfx )
    call dfftw_execute_dft_r2c( forw_x, wdudzpx, wdudzfx )
    call dfftw_execute_dft_r2c( forw_x, udvdxpx, udvdxfx )
    call dfftw_execute_dft_r2c( forw_x, vdvdypx, vdvdyfx )
    call dfftw_execute_dft_r2c( forw_x, wdvdzpx, wdvdzfx )
    call dfftw_execute_dft_r2c( forw_x, udwdxpx, udwdxfx )
    call dfftw_execute_dft_r2c( forw_x, vdwdypx, vdwdyfx )
    call dfftw_execute_dft_r2c( forw_x, wdwdzpx, wdwdzfx )
#endif

    ! Multiply together and time-average
    do jx = 1, nx/2 + 1
        tavg_turbspecx(jx,jy,jz)%uf = tavg_turbspecx(jx,jy,jz)%uf +       &
            ufx(jx)*tavg_dt
        tavg_turbspecx(jx,jy,jz)%vf = tavg_turbspecx(jx,jy,jz)%vf +       &
            vfx(jx)*tavg_dt
        tavg_turbspecx(jx,jy,jz)%wf = tavg_turbspecx(jx,jy,jz)%wf +       &
            wfx(jx)*tavg_dt

        tavg_turbspecx(jx,jy,jz)%uu = tavg_turbspecx(jx,jy,jz)%uu +       &
            ufx(jx)*conjg(ufx(jx))*tavg_dt
        tavg_turbspecx(jx,jy,jz)%vv = tavg_turbspecx(jx,jy,jz)%vv +       &
            vfx(jx)*conjg(vfx(jx))*tavg_dt
        tavg_turbspecx(jx,jy,jz)%ww = tavg_turbspecx(jx,jy,jz)%ww +       &
            wfx(jx)*conjg(wfx(jx))*tavg_dt
        !tavg_turbspecx(jx,jy,jz)%vel2 = tavg_turbspecx(jx,jy,jz)%vel2 +   &
        !    velfx(jx)*conjg(velfx(jx))*tavg_dt

        tavg_turbspecx(jx,jy,jz)%uv = tavg_turbspecx(jx,jy,jz)%uv +       &
            ufx(jx)*conjg(vfx(jx))*tavg_dt
        tavg_turbspecx(jx,jy,jz)%uw = tavg_turbspecx(jx,jy,jz)%uw +       &
            ufx(jx)*conjg(wfx(jx))*tavg_dt
        tavg_turbspecx(jx,jy,jz)%vw = tavg_turbspecx(jx,jy,jz)%vw +       &
            vfx(jx)*conjg(wfx(jx))*tavg_dt

        tavg_turbspecx(jx,jy,jz)%vortxf = tavg_turbspecx(jx,jy,jz)%vortxf + &
            vortxfx(jx)*tavg_dt
        tavg_turbspecx(jx,jy,jz)%vortyf = tavg_turbspecx(jx,jy,jz)%vortyf + &
            vortyfx(jx)*tavg_dt
        tavg_turbspecx(jx,jy,jz)%vortzf = tavg_turbspecx(jx,jy,jz)%vortzf + &
            vortzfx(jx)*tavg_dt

        tavg_turbspecx(jx,jy,jz)%vortx2 = tavg_turbspecx(jx,jy,jz)%vortx2 + &
            vortxfx(jx)*conjg(vortxfx(jx))*tavg_dt
        tavg_turbspecx(jx,jy,jz)%vorty2 = tavg_turbspecx(jx,jy,jz)%vorty2 + &
            vortyfx(jx)*conjg(vortyfx(jx))*tavg_dt
        tavg_turbspecx(jx,jy,jz)%vortz2 = tavg_turbspecx(jx,jy,jz)%vortz2 + &
            vortzfx(jx)*conjg(vortzfx(jx))*tavg_dt
        !tavg_turbspecx(jx,jy,jz)%vort2 = tavg_turbspecx(jx,jy,jz)%vort2 +   &
        !    real(vortfx(jx)*conjg(vortfx(jx)))*tavg_dt

#ifdef PPSCALARS
        tavg_scal_turbspecx(jx,jy,jz)%thetaf = tavg_scal_turbspecx(jx,jy,jz)%thetaf + &
            thetafx(jx)*tavg_dt

        tavg_scal_turbspecx(jx,jy,jz)%theta2 = tavg_scal_turbspecx(jx,jy,jz)%theta2 + &
            thetafx(jx)*conjg(thetafx(jx))*tavg_dt
        tavg_scal_turbspecx(jx,jy,jz)%utheta = tavg_scal_turbspecx(jx,jy,jz)%utheta + &
            ufx(jx)*conjg(thetafx(jx))*tavg_dt
        tavg_scal_turbspecx(jx,jy,jz)%vtheta = tavg_scal_turbspecx(jx,jy,jz)%vtheta + &
            vfx(jx)*conjg(thetafx(jx))*tavg_dt
        tavg_scal_turbspecx(jx,jy,jz)%wtheta = tavg_scal_turbspecx(jx,jy,jz)%wtheta + &
            wfx(jx)*conjg(thetafx(jx))*tavg_dt

#ifdef PPOUTPUT_SPECBUDG
        ! Scalar gradient, dTdxj_hat
        tavg_scal_specbudgx(jx,jy,jz)%dTdx = tavg_scal_specbudgx(jx,jy,jz)%dTdx + &
            dTdxfx(jx)*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%dTdy = tavg_scal_specbudgx(jx,jy,jz)%dTdy + &
            dTdyfx(jx)*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%dTdz = tavg_scal_specbudgx(jx,jy,jz)%dTdz + &
            dTdzfx(jx)*tavg_dt

        ! Scalar-Scalar gradient, T_hat*dTdxj_hat
        tavg_scal_specbudgx(jx,jy,jz)%Th_dTdxh = tavg_scal_specbudgx(jx,jy,jz)%Th_dTdxh + &
            thetafx(jx)*conjg(dTdxfx(jx))*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%Th_dTdyh = tavg_scal_specbudgx(jx,jy,jz)%Th_dTdyh + &
            thetafx(jx)*conjg(dTdyfx(jx))*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%Th_dTdzh = tavg_scal_specbudgx(jx,jy,jz)%Th_dTdzh + &
            thetafx(jx)*conjg(dTdzfx(jx))*tavg_dt

        ! Velocity-Scalar Gradient, (uj*dTdxj)_hat
        tavg_scal_specbudgx(jx,jy,jz)%udTdxh = tavg_scal_specbudgx(jx,jy,jz)%udTdxh + &
            udTdxfx(jx)*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%vdTdyh = tavg_scal_specbudgx(jx,jy,jz)%vdTdyh + &
            vdTdyfx(jx)*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%wdTdzh = tavg_scal_specbudgx(jx,jy,jz)%wdTdzh + &
            wdTdzfx(jx)*tavg_dt

        ! Scalar-Velocity-Scalar Gradient, T_hat*((uj*dTdxj)_hat)
        tavg_scal_specbudgx(jx,jy,jz)%Th_udTdxh = tavg_scal_specbudgx(jx,jy,jz)%Th_udTdxh + &
            thetafx(jx)*conjg(udTdxfx(jx))*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%Th_vdTdyh = tavg_scal_specbudgx(jx,jy,jz)%Th_vdTdyh + &
            thetafx(jx)*conjg(vdTdyfx(jx))*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%Th_wdTdzh = tavg_scal_specbudgx(jx,jy,jz)%Th_wdTdzh + &
            thetafx(jx)*conjg(wdTdzfx(jx))*tavg_dt

        ! Scalar gradient-scalar gradient, dTdxj_hat*dTdxj_hat
        tavg_scal_specbudgx(jx,jy,jz)%TxTx = tavg_scal_specbudgx(jx,jy,jz)%TxTx + &
            dTdxfx(jx)*conjg(dTdxfx(jx))*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%TyTy = tavg_scal_specbudgx(jx,jy,jz)%TyTy + &
            dTdyfx(jx)*conjg(dTdyfx(jx))*tavg_dt
        tavg_scal_specbudgx(jx,jy,jz)%TzTz = tavg_scal_specbudgx(jx,jy,jz)%TzTz + &
            dTdzfx(jx)*conjg(dTdzfx(jx))*tavg_dt

        ! Laplacian, (nu/Pr)*lap(theta)
        tavg_scal_specbudgx(jx,jy,jz)%lapT = tavg_scal_specbudgx(jx,jy,jz)%lapT + &
            lapTfx(jx)*tavg_dt

        ! Scalar-Laplacian, theta*(nu/Pr)*lap(theta)
        tavg_scal_specbudgx(jx,jy,jz)%TlapT = tavg_scal_specbudgx(jx,jy,jz)%TlapT + &
            thetafx(jx)*conjg(lapTfx(jx))*tavg_dt

#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
        ! Mean pressure
        tavg_specbudgx(jx,jy,jz)%p = tavg_specbudgx(jx,jy,jz)%p + &
            pfx(jx)*tavg_dt 

        ! Mean velocity gradients, duidxj
        tavg_specbudgx(jx,jy,jz)%dudx = tavg_specbudgx(jx,jy,jz)%dudx + &
            dudxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dudy = tavg_specbudgx(jx,jy,jz)%dudy + &
            dudyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dudz = tavg_specbudgx(jx,jy,jz)%dudz + &
            dudzfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dvdx = tavg_specbudgx(jx,jy,jz)%dvdx + &
            dvdxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dvdy = tavg_specbudgx(jx,jy,jz)%dvdy + &
            dvdyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dvdz = tavg_specbudgx(jx,jy,jz)%dvdz + &
            dvdzfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dwdx = tavg_specbudgx(jx,jy,jz)%dwdx + &
            dwdxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dwdy = tavg_specbudgx(jx,jy,jz)%dwdy + &
            dwdyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dwdz = tavg_specbudgx(jx,jy,jz)%dwdz + &
            dwdzfx(jx)*tavg_dt

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
        tavg_specbudgx(jx,jy,jz)%uxux = tavg_specbudgx(jx,jy,jz)%uxux + &
            dudxfx(jx)*conjg(dudxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uyuy = tavg_specbudgx(jx,jy,jz)%uyuy + &
            dudyfx(jx)*conjg(dudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uzuz = tavg_specbudgx(jx,jy,jz)%uzuz + &
            dudzfx(jx)*conjg(dudzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vxvx = tavg_specbudgx(jx,jy,jz)%vxvx + &
            dvdxfx(jx)*conjg(dvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vyvy = tavg_specbudgx(jx,jy,jz)%vyvy + &
            dvdyfx(jx)*conjg(dvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vzvz = tavg_specbudgx(jx,jy,jz)%vzvz + &
            dvdzfx(jx)*conjg(dvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wxwx = tavg_specbudgx(jx,jy,jz)%wxwx + &
            dwdxfx(jx)*conjg(dwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wywy = tavg_specbudgx(jx,jy,jz)%wywy + &
            dwdyfx(jx)*conjg(dwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wzwz = tavg_specbudgx(jx,jy,jz)%wzwz + &
            dwdzfx(jx)*conjg(dwdzfx(jx))*tavg_dt

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
        tavg_specbudgx(jx,jy,jz)%uxvx = tavg_specbudgx(jx,jy,jz)%uxvx + &
            dudxfx(jx)*conjg(dvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uyvy = tavg_specbudgx(jx,jy,jz)%uyvy + &
            dudyfx(jx)*conjg(dvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uzvz = tavg_specbudgx(jx,jy,jz)%uzvz + &
            dudzfx(jx)*conjg(dvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uxwx = tavg_specbudgx(jx,jy,jz)%uxwx + &
            dudxfx(jx)*conjg(dwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uywy = tavg_specbudgx(jx,jy,jz)%uywy + &
            dudyfx(jx)*conjg(dwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uzwz = tavg_specbudgx(jx,jy,jz)%uzwz + &
            dudzfx(jx)*conjg(dwdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vxwx = tavg_specbudgx(jx,jy,jz)%vxwx + &
            dvdxfx(jx)*conjg(dwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vywy = tavg_specbudgx(jx,jy,jz)%vywy + &
            dvdyfx(jx)*conjg(dwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vzwz = tavg_specbudgx(jx,jy,jz)%vzwz + &
            dvdzfx(jx)*conjg(dwdzfx(jx))*tavg_dt

        ! Mean pressure gradients, dpdxi
        tavg_specbudgx(jx,jy,jz)%dpdx = tavg_specbudgx(jx,jy,jz)%dpdx + &
            dpdxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dpdy = tavg_specbudgx(jx,jy,jz)%dpdy + &
            dpdyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%dpdz = tavg_specbudgx(jx,jy,jz)%dpdz + &
            dpdzfx(jx)*tavg_dt

        ! Mean velocity * velocity gradient product, ui_hat*dujdxk_hat
        tavg_specbudgx(jx,jy,jz)%uh_dudxh = tavg_specbudgx(jx,jy,jz)%uh_dudxh + &
            ufx(jx)*conjg(dudxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dudyh = tavg_specbudgx(jx,jy,jz)%uh_dudyh + &
            ufx(jx)*conjg(dudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dudzh = tavg_specbudgx(jx,jy,jz)%uh_dudzh + &
            ufx(jx)*conjg(dudzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dvdxh = tavg_specbudgx(jx,jy,jz)%uh_dvdxh + &
            ufx(jx)*conjg(dvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dvdyh = tavg_specbudgx(jx,jy,jz)%uh_dvdyh + &
            ufx(jx)*conjg(dvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dvdzh = tavg_specbudgx(jx,jy,jz)%uh_dvdzh + &
            ufx(jx)*conjg(dvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dwdxh = tavg_specbudgx(jx,jy,jz)%uh_dwdxh + &
            ufx(jx)*conjg(dwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dwdyh = tavg_specbudgx(jx,jy,jz)%uh_dwdyh + &
            ufx(jx)*conjg(dwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_dwdzh = tavg_specbudgx(jx,jy,jz)%uh_dwdzh + &
            ufx(jx)*conjg(dwdzfx(jx))*tavg_dt

        tavg_specbudgx(jx,jy,jz)%vh_dudxh = tavg_specbudgx(jx,jy,jz)%vh_dudxh + &
            vfx(jx)*conjg(dudxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dudyh = tavg_specbudgx(jx,jy,jz)%vh_dudyh + &
            vfx(jx)*conjg(dudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dudzh = tavg_specbudgx(jx,jy,jz)%vh_dudzh + &
            vfx(jx)*conjg(dudzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dvdxh = tavg_specbudgx(jx,jy,jz)%vh_dvdxh + &
            vfx(jx)*conjg(dvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dvdyh = tavg_specbudgx(jx,jy,jz)%vh_dvdyh + &
            vfx(jx)*conjg(dvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dvdzh = tavg_specbudgx(jx,jy,jz)%vh_dvdzh + &
            vfx(jx)*conjg(dvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dwdxh = tavg_specbudgx(jx,jy,jz)%vh_dwdxh + &
            vfx(jx)*conjg(dwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dwdyh = tavg_specbudgx(jx,jy,jz)%vh_dwdyh + &
            vfx(jx)*conjg(dwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_dwdzh = tavg_specbudgx(jx,jy,jz)%vh_dwdzh + &
            vfx(jx)*conjg(dwdzfx(jx))*tavg_dt

        tavg_specbudgx(jx,jy,jz)%wh_dudxh = tavg_specbudgx(jx,jy,jz)%wh_dudxh + &
            wfx(jx)*conjg(dudxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dudyh = tavg_specbudgx(jx,jy,jz)%wh_dudyh + &
            wfx(jx)*conjg(dudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dudzh = tavg_specbudgx(jx,jy,jz)%wh_dudzh + &
            wfx(jx)*conjg(dudzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dvdxh = tavg_specbudgx(jx,jy,jz)%wh_dvdxh + &
            wfx(jx)*conjg(dvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dvdyh = tavg_specbudgx(jx,jy,jz)%wh_dvdyh + &
            wfx(jx)*conjg(dvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dvdzh = tavg_specbudgx(jx,jy,jz)%wh_dvdzh + &
            wfx(jx)*conjg(dvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dwdxh = tavg_specbudgx(jx,jy,jz)%wh_dwdxh + &
            wfx(jx)*conjg(dwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dwdyh = tavg_specbudgx(jx,jy,jz)%wh_dwdyh + &
            wfx(jx)*conjg(dwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_dwdzh = tavg_specbudgx(jx,jy,jz)%wh_dwdzh + &
            wfx(jx)*conjg(dwdzfx(jx))*tavg_dt

        ! Mean velocity * velocity gradient product, uk*duidxk_hat
        tavg_specbudgx(jx,jy,jz)%ududxh = tavg_specbudgx(jx,jy,jz)%ududxh + &
            ududxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vdudyh = tavg_specbudgx(jx,jy,jz)%vdudyh + &
            vdudyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wdudzh = tavg_specbudgx(jx,jy,jz)%wdudzh + &
            wdudzfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%udvdxh = tavg_specbudgx(jx,jy,jz)%udvdxh + &
            udvdxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vdvdyh = tavg_specbudgx(jx,jy,jz)%vdvdyh + &
            vdvdyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wdvdzh = tavg_specbudgx(jx,jy,jz)%wdvdzh + &
            wdvdzfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%udwdxh = tavg_specbudgx(jx,jy,jz)%udwdxh + &
            udwdxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vdwdyh = tavg_specbudgx(jx,jy,jz)%vdwdyh + &
            vdwdyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wdwdzh = tavg_specbudgx(jx,jy,jz)%wdwdzh + &
            wdwdzfx(jx)*tavg_dt

        ! Mean vel * vel * vel gradient product, uj_hat*(uk*duidxk)_hat
        tavg_specbudgx(jx,jy,jz)%uh_ududxh = tavg_specbudgx(jx,jy,jz)%uh_ududxh + &
            ufx(jx)*conjg(ududxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_vdudyh = tavg_specbudgx(jx,jy,jz)%uh_vdudyh + &
            ufx(jx)*conjg(vdudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_wdudzh = tavg_specbudgx(jx,jy,jz)%uh_wdudzh + &
            ufx(jx)*conjg(wdudzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_ududxh = tavg_specbudgx(jx,jy,jz)%vh_ududxh + &
            vfx(jx)*conjg(ududxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_vdudyh = tavg_specbudgx(jx,jy,jz)%vh_vdudyh + &
            vfx(jx)*conjg(vdudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_wdudzh = tavg_specbudgx(jx,jy,jz)%vh_wdudzh + &
            vfx(jx)*conjg(wdudzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_ududxh = tavg_specbudgx(jx,jy,jz)%wh_ududxh + &
            wfx(jx)*conjg(ududxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_vdudyh = tavg_specbudgx(jx,jy,jz)%wh_vdudyh + &
            wfx(jx)*conjg(vdudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_wdudzh = tavg_specbudgx(jx,jy,jz)%wh_wdudzh + &
            wfx(jx)*conjg(wdudzfx(jx))*tavg_dt

        tavg_specbudgx(jx,jy,jz)%uh_udvdxh = tavg_specbudgx(jx,jy,jz)%uh_udvdxh + &
            ufx(jx)*conjg(udvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_vdvdyh = tavg_specbudgx(jx,jy,jz)%uh_vdvdyh + &
            ufx(jx)*conjg(vdvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_wdvdzh = tavg_specbudgx(jx,jy,jz)%uh_wdvdzh + &
            ufx(jx)*conjg(wdvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_udvdxh = tavg_specbudgx(jx,jy,jz)%vh_udvdxh + &
            vfx(jx)*conjg(udvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_vdvdyh = tavg_specbudgx(jx,jy,jz)%vh_vdvdyh + &
            vfx(jx)*conjg(vdvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_wdvdzh = tavg_specbudgx(jx,jy,jz)%vh_wdvdzh + &
            vfx(jx)*conjg(wdvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_udvdxh = tavg_specbudgx(jx,jy,jz)%wh_udvdxh + &
            wfx(jx)*conjg(udvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_vdvdyh = tavg_specbudgx(jx,jy,jz)%wh_vdvdyh + &
            wfx(jx)*conjg(vdvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_wdvdzh = tavg_specbudgx(jx,jy,jz)%wh_wdvdzh + &
            wfx(jx)*conjg(wdvdzfx(jx))*tavg_dt

        tavg_specbudgx(jx,jy,jz)%uh_udwdxh = tavg_specbudgx(jx,jy,jz)%uh_udwdxh + &
            ufx(jx)*conjg(udwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_vdwdyh = tavg_specbudgx(jx,jy,jz)%uh_vdwdyh + &
            ufx(jx)*conjg(vdwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%uh_wdwdzh = tavg_specbudgx(jx,jy,jz)%uh_wdwdzh + &
            ufx(jx)*conjg(wdwdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_udwdxh = tavg_specbudgx(jx,jy,jz)%vh_udwdxh + &
            vfx(jx)*conjg(udwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_vdwdyh = tavg_specbudgx(jx,jy,jz)%vh_vdwdyh + &
            vfx(jx)*conjg(vdwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vh_wdwdzh = tavg_specbudgx(jx,jy,jz)%vh_wdwdzh + &
            vfx(jx)*conjg(wdwdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_udwdxh = tavg_specbudgx(jx,jy,jz)%wh_udwdxh + &
            wfx(jx)*conjg(udwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_vdwdyh = tavg_specbudgx(jx,jy,jz)%wh_vdwdyh + &
            wfx(jx)*conjg(vdwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wh_wdwdzh = tavg_specbudgx(jx,jy,jz)%wh_wdwdzh + &
            wfx(jx)*conjg(wdwdzfx(jx))*tavg_dt

        ! Mean velocity * pressure gradient product, ui*dpdxj
        tavg_specbudgx(jx,jy,jz)%udpdx = tavg_specbudgx(jx,jy,jz)%udpdx + &
            ufx(jx)*conjg(dpdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%udpdy = tavg_specbudgx(jx,jy,jz)%udpdy + &
            ufx(jx)*conjg(dpdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%udpdz = tavg_specbudgx(jx,jy,jz)%udpdz + &
            ufx(jx)*conjg(dpdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vdpdx = tavg_specbudgx(jx,jy,jz)%vdpdx + &
            vfx(jx)*conjg(dpdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vdpdy = tavg_specbudgx(jx,jy,jz)%vdpdy + &
            vfx(jx)*conjg(dpdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vdpdz = tavg_specbudgx(jx,jy,jz)%vdpdz + &
            vfx(jx)*conjg(dpdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wdpdx = tavg_specbudgx(jx,jy,jz)%wdpdx + &
            wfx(jx)*conjg(dpdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wdpdy = tavg_specbudgx(jx,jy,jz)%wdpdy + &
            wfx(jx)*conjg(dpdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wdpdz = tavg_specbudgx(jx,jy,jz)%wdpdz + &
            wfx(jx)*conjg(dpdzfx(jx))*tavg_dt

        ! Mean pressure * velocity gradient product, p*duidxj
        tavg_specbudgx(jx,jy,jz)%pdudx = tavg_specbudgx(jx,jy,jz)%pdudx + &
            pfx(jx)*conjg(dudxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdudy = tavg_specbudgx(jx,jy,jz)%pdudy + &
            pfx(jx)*conjg(dudyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdudz = tavg_specbudgx(jx,jy,jz)%pdudz + &
            pfx(jx)*conjg(dudzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdvdx = tavg_specbudgx(jx,jy,jz)%pdvdx + &
            pfx(jx)*conjg(dvdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdvdy = tavg_specbudgx(jx,jy,jz)%pdvdy + &
            pfx(jx)*conjg(dvdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdvdz = tavg_specbudgx(jx,jy,jz)%pdvdz + &
            pfx(jx)*conjg(dvdzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdwdx = tavg_specbudgx(jx,jy,jz)%pdwdx + &
            pfx(jx)*conjg(dwdxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdwdy = tavg_specbudgx(jx,jy,jz)%pdwdy + &
            pfx(jx)*conjg(dwdyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%pdwdz = tavg_specbudgx(jx,jy,jz)%pdwdz + &
            pfx(jx)*conjg(dwdzfx(jx))*tavg_dt

        ! Mean Laplacian, nu*lap(uj)
        tavg_specbudgx(jx,jy,jz)%lapu = tavg_specbudgx(jx,jy,jz)%lapu + &
            divtxfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%lapv = tavg_specbudgx(jx,jy,jz)%lapv + &
            divtyfx(jx)*tavg_dt
        tavg_specbudgx(jx,jy,jz)%lapw = tavg_specbudgx(jx,jy,jz)%lapw + &
            divtzfx(jx)*tavg_dt

        ! Mean Vel-Laplacian, nu*ui*lap(uj)
        tavg_specbudgx(jx,jy,jz)%ulapu = tavg_specbudgx(jx,jy,jz)%ulapu + &
            ufx(jx)*conjg(divtxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%ulapv = tavg_specbudgx(jx,jy,jz)%ulapv + &
            ufx(jx)*conjg(divtyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%ulapw = tavg_specbudgx(jx,jy,jz)%ulapw + &
            ufx(jx)*conjg(divtzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vlapu = tavg_specbudgx(jx,jy,jz)%vlapu + &
            vfx(jx)*conjg(divtxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vlapv = tavg_specbudgx(jx,jy,jz)%vlapv + &
            vfx(jx)*conjg(divtyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%vlapw = tavg_specbudgx(jx,jy,jz)%vlapw + &
            vfx(jx)*conjg(divtzfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wlapu = tavg_specbudgx(jx,jy,jz)%wlapu + &
            wfx(jx)*conjg(divtxfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wlapv = tavg_specbudgx(jx,jy,jz)%wlapv + &
            wfx(jx)*conjg(divtyfx(jx))*tavg_dt
        tavg_specbudgx(jx,jy,jz)%wlapw = tavg_specbudgx(jx,jy,jz)%wlapw + &
            wfx(jx)*conjg(divtzfx(jx))*tavg_dt
#endif
    end do

end do
end do

do jx = 1, nx
do jz = 1, nz

    ! Normalize
    upy = const2 * u_w(jx,1:ny,jz)
    vpy = const2 * v_w(jx,1:ny,jz)
    wpy = const2 * w(jx,1:ny,jz)
    !velpy = const2 * vel(jx,1:ny,jz)
    vortxpy = const2 * vortx(jx,1:ny,jz)
    vortypy = const2 * vorty(jx,1:ny,jz)
    vortzpy = const2 * vortz(jx,1:ny,jz)
    !vortpy = const2 * vort(jx,1:ny,jz)
#ifdef PPSCALARS
    thetapy = const2 * theta_w(jx,1:ny,jz)
#endif

    ! Take 1d transform, y --> ky
    call dfftw_execute_dft_r2c( forw_y, upy,     ufy )
    call dfftw_execute_dft_r2c( forw_y, vpy,     vfy )
    call dfftw_execute_dft_r2c( forw_y, wpy,     wfy )
    !call dfftw_execute_dft_r2c( forw_y, velpy,   velfy )
    call dfftw_execute_dft_r2c( forw_y, vortxpy, vortxfy )
    call dfftw_execute_dft_r2c( forw_y, vortypy, vortyfy )
    call dfftw_execute_dft_r2c( forw_y, vortzpy, vortzfy )
    !call dfftw_execute_dft_r2c( forw_y, vortpy,  vortfy )
#ifdef PPSCALARS
    call dfftw_execute_dft_r2c( forw_y, thetapy, thetafy )
#endif

#ifdef PPOUTPUT_SPECBUDG
    ! Normalize
    ppy = const2 * pres_real(jx,1:ny,jz)
    dpdxpy = const2 * dpdx_real(jx,1:ny,jz)
    dpdypy = const2 * dpdx_real(jx,1:ny,jz)
    dpdzpy = const2 * dpdx_real(jx,1:ny,jz)
    dudxpy = const2 * dudx_w(jx,1:ny,jz)
    dudypy = const2 * dudy_w(jx,1:ny,jz)
    dudzpy = const2 * dudz(jx,1:ny,jz)
    dvdxpy = const2 * dvdx_w(jx,1:ny,jz)
    dvdypy = const2 * dvdy_w(jx,1:ny,jz)
    dvdzpy = const2 * dvdz(jx,1:ny,jz)
    dwdxpy = const2 * dwdx(jx,1:ny,jz)
    dwdypy = const2 * dwdy(jx,1:ny,jz)
    dwdzpy = const2 * dwdz_w(jx,1:ny,jz)
    divtxpy = const2 * divtx_w(jx,1:ny,jz)
    divtypy = const2 * divty_w(jx,1:ny,jz)
    divtzpy = const2 * divtz(jx,1:ny,jz)
    ududxpy = const2 * u_w(jx,1:ny,jz)*dudx_w(jx,1:ny,jz)
    vdudypy = const2 * v_w(jx,1:ny,jz)*dudy_w(jx,1:ny,jz)
    wdudzpy = const2 * w(jx,1:ny,jz)*dudz(jx,1:ny,jz)
    udvdxpy = const2 * u_w(jx,1:ny,jz)*dvdx_w(jx,1:ny,jz)
    vdvdypy = const2 * v_w(jx,1:ny,jz)*dvdy_w(jx,1:ny,jz)
    wdvdzpy = const2 * w(jx,1:ny,jz)*dvdz(jx,1:ny,jz)
    udwdxpy = const2 * u_w(jx,1:ny,jz)*dwdx(jx,1:ny,jz)
    vdwdypy = const2 * v_w(jx,1:ny,jz)*dwdy(jx,1:ny,jz)
    wdwdzpy = const2 * w(jx,1:ny,jz)*dwdz_w(jx,1:ny,jz)

    ! Take 1d transform, y --> ky
    call dfftw_execute_dft_r2c( forw_y, ppy, pfy )
    call dfftw_execute_dft_r2c( forw_y, dpdxpy, dpdxfy )
    call dfftw_execute_dft_r2c( forw_y, dpdypy, dpdyfy )
    call dfftw_execute_dft_r2c( forw_y, dpdzpy, dpdzfy )
    call dfftw_execute_dft_r2c( forw_y, dudxpy, dudxfy )
    call dfftw_execute_dft_r2c( forw_y, dudypy, dudyfy )
    call dfftw_execute_dft_r2c( forw_y, dudzpy, dudzfy )
    call dfftw_execute_dft_r2c( forw_y, dvdxpy, dvdxfy )
    call dfftw_execute_dft_r2c( forw_y, dvdypy, dvdyfy )
    call dfftw_execute_dft_r2c( forw_y, dvdzpy, dvdzfy )
    call dfftw_execute_dft_r2c( forw_y, dwdxpy, dwdxfy )
    call dfftw_execute_dft_r2c( forw_y, dwdypy, dwdyfy )
    call dfftw_execute_dft_r2c( forw_y, dwdzpy, dwdzfy )
    call dfftw_execute_dft_r2c( forw_y, divtxpy, divtxfy )
    call dfftw_execute_dft_r2c( forw_y, divtypy, divtyfy )
    call dfftw_execute_dft_r2c( forw_y, divtzpy, divtzfy )
    call dfftw_execute_dft_r2c( forw_y, ududxpy, ududxfy )
    call dfftw_execute_dft_r2c( forw_y, vdudypy, vdudyfy )
    call dfftw_execute_dft_r2c( forw_y, wdudzpy, wdudzfy )
    call dfftw_execute_dft_r2c( forw_y, udvdxpy, udvdxfy )
    call dfftw_execute_dft_r2c( forw_y, vdvdypy, vdvdyfy )
    call dfftw_execute_dft_r2c( forw_y, wdvdzpy, wdvdzfy )
    call dfftw_execute_dft_r2c( forw_y, udwdxpy, udwdxfy )
    call dfftw_execute_dft_r2c( forw_y, vdwdypy, vdwdyfy )
    call dfftw_execute_dft_r2c( forw_y, wdwdzpy, wdwdzfy )
#endif

    ! Multiply together and time-average
    do jy = 1, ny/2 + 1
        tavg_turbspecy(jx,jy,jz)%uf = tavg_turbspecy(jx,jy,jz)%uf +       &
            ufy(jy)*tavg_dt
        tavg_turbspecy(jx,jy,jz)%vf = tavg_turbspecy(jx,jy,jz)%vf +       &
            vfy(jy)*tavg_dt
        tavg_turbspecy(jx,jy,jz)%wf = tavg_turbspecy(jx,jy,jz)%wf +       &
            wfy(jy)*tavg_dt

        tavg_turbspecy(jx,jy,jz)%uu = tavg_turbspecy(jx,jy,jz)%uu +       &
            ufy(jy)*conjg(ufy(jy))*tavg_dt
        tavg_turbspecy(jx,jy,jz)%vv = tavg_turbspecy(jx,jy,jz)%vv +       &
            vfy(jy)*conjg(vfy(jy))*tavg_dt
        tavg_turbspecy(jx,jy,jz)%ww = tavg_turbspecy(jx,jy,jz)%ww +       &
            wfy(jy)*conjg(wfy(jy))*tavg_dt
        !tavg_turbspecy(jx,jy,jz)%vel2 = tavg_turbspec(jx,jy,jz)%vel2 +   &
        !    velfy(jy)*conjg(velfy(jy))*tavg_dt

        tavg_turbspecy(jx,jy,jz)%uv = tavg_turbspecy(jx,jy,jz)%uv +        &
            ufy(jy)*conjg(vfy(jy))*tavg_dt
        tavg_turbspecy(jx,jy,jz)%uw = tavg_turbspecy(jx,jy,jz)%uw +        &
            ufy(jy)*conjg(wfy(jy))*tavg_dt
        tavg_turbspecy(jx,jy,jz)%vw = tavg_turbspecy(jx,jy,jz)%vw +        &
            vfy(jy)*conjg(wfy(jy))*tavg_dt

        tavg_turbspecy(jx,jy,jz)%vortxf = tavg_turbspecy(jx,jy,jz)%vortxf + &
            vortxfy(jy)*tavg_dt
        tavg_turbspecy(jx,jy,jz)%vortyf = tavg_turbspecy(jx,jy,jz)%vortyf + &
            vortyfy(jy)*tavg_dt
        tavg_turbspecy(jx,jy,jz)%vortzf = tavg_turbspecy(jx,jy,jz)%vortzf + &
            vortzfy(jy)*tavg_dt

        tavg_turbspecy(jx,jy,jz)%vortx2 = tavg_turbspecy(jx,jy,jz)%vortx2 +  &
            vortxfy(jy)*conjg(vortxfy(jy))*tavg_dt
        tavg_turbspecy(jx,jy,jz)%vorty2 = tavg_turbspecy(jx,jy,jz)%vorty2 +  &
            vortyfy(jy)*conjg(vortyfy(jy))*tavg_dt
        tavg_turbspecy(jx,jy,jz)%vortz2 = tavg_turbspecy(jx,jy,jz)%vortz2 +  &
            real(vortzfy(jy)*conjg(vortzfy(jy)))*tavg_dt
        !tavg_turbspecy(jx,jy,jz)%vort2 = tavg_turbspecy(jx,jy,jz)%vort2 +    &
        !    real(vortfy(jy)*conjg(vortfy(jy)))*tavg_dt
#ifdef PPSCALARS
        tavg_scal_turbspecy(jx,jy,jz)%thetaf = tavg_scal_turbspecy(jx,jy,jz)%thetaf + &
            thetafy(jy)*tavg_dt

        tavg_scal_turbspecy(jx,jy,jz)%theta2 = tavg_scal_turbspecy(jx,jy,jz)%theta2 + &
            thetafy(jy)*conjg(thetafy(jy))*tavg_dt
        tavg_scal_turbspecy(jx,jy,jz)%utheta = tavg_scal_turbspecy(jx,jy,jz)%utheta + &
            ufy(jy)*conjg(thetafy(jy))*tavg_dt
        tavg_scal_turbspecy(jx,jy,jz)%vtheta = tavg_scal_turbspecy(jx,jy,jz)%vtheta + &
            vfy(jy)*conjg(thetafy(jy))*tavg_dt
        tavg_scal_turbspecy(jx,jy,jz)%wtheta = tavg_scal_turbspecy(jx,jy,jz)%wtheta + &
            wfy(jy)*conjg(thetafy(jy))*tavg_dt
#endif

#ifdef PPOUTPUT_SPECBUDG
        ! Mean pressure
        tavg_specbudgy(jx,jy,jz)%p = tavg_specbudgy(jx,jy,jz)%p + &
            pfy(jy)*tavg_dt 

        ! Mean velocity gradients, duidxj
        tavg_specbudgy(jx,jy,jz)%dudx = tavg_specbudgy(jx,jy,jz)%dudx + &
            dudxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dudy = tavg_specbudgy(jx,jy,jz)%dudy + &
            dudyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dudz = tavg_specbudgy(jx,jy,jz)%dudz + &
            dudzfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dvdx = tavg_specbudgy(jx,jy,jz)%dvdx + &
            dvdxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dvdy = tavg_specbudgy(jx,jy,jz)%dvdy + &
            dvdyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dvdz = tavg_specbudgy(jx,jy,jz)%dvdz + &
            dvdzfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dwdx = tavg_specbudgy(jx,jy,jz)%dwdx + &
            dwdxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dwdy = tavg_specbudgy(jx,jy,jz)%dwdy + &
            dwdyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dwdz = tavg_specbudgy(jx,jy,jz)%dwdz + &
            dwdzfy(jy)*tavg_dt

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
        tavg_specbudgy(jx,jy,jz)%uxux = tavg_specbudgy(jx,jy,jz)%uxux + &
            dudxfy(jy)*conjg(dudxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uyuy = tavg_specbudgy(jx,jy,jz)%uyuy + &
            dudyfy(jy)*conjg(dudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uzuz = tavg_specbudgy(jx,jy,jz)%uzuz + &
            dudzfy(jy)*conjg(dudzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vxvx = tavg_specbudgy(jx,jy,jz)%vxvx + &
            dvdxfy(jy)*conjg(dvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vyvy = tavg_specbudgy(jx,jy,jz)%vyvy + &
            dvdyfy(jy)*conjg(dvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vzvz = tavg_specbudgy(jx,jy,jz)%vzvz + &
            dvdzfy(jy)*conjg(dvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wxwx = tavg_specbudgy(jx,jy,jz)%wxwx + &
            dwdxfy(jy)*conjg(dwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wywy = tavg_specbudgy(jx,jy,jz)%wywy + &
            dwdyfy(jy)*conjg(dwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wzwz = tavg_specbudgy(jx,jy,jz)%wzwz + &
            dwdzfy(jy)*conjg(dwdzfy(jy))*tavg_dt

        ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
        tavg_specbudgy(jx,jy,jz)%uxvx = tavg_specbudgy(jx,jy,jz)%uxvx + &
            dudxfy(jy)*conjg(dvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uyvy = tavg_specbudgy(jx,jy,jz)%uyvy + &
            dudyfy(jy)*conjg(dvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uzvz = tavg_specbudgy(jx,jy,jz)%uzvz + &
            dudzfy(jy)*conjg(dvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uxwx = tavg_specbudgy(jx,jy,jz)%uxwx + &
            dudxfy(jy)*conjg(dwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uywy = tavg_specbudgy(jx,jy,jz)%uywy + &
            dudyfy(jy)*conjg(dwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uzwz = tavg_specbudgy(jx,jy,jz)%uzwz + &
            dudzfy(jy)*conjg(dwdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vxwx = tavg_specbudgy(jx,jy,jz)%vxwx + &
            dvdxfy(jy)*conjg(dwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vywy = tavg_specbudgy(jx,jy,jz)%vywy + &
            dvdyfy(jy)*conjg(dwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vzwz = tavg_specbudgy(jx,jy,jz)%vzwz + &
            dvdzfy(jy)*conjg(dwdzfy(jy))*tavg_dt

        ! Mean pressure gradients, dpdxi
        tavg_specbudgy(jx,jy,jz)%dpdx = tavg_specbudgy(jx,jy,jz)%dpdx + &
            dpdxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dpdy = tavg_specbudgy(jx,jy,jz)%dpdy + &
            dpdyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%dpdz = tavg_specbudgy(jx,jy,jz)%dpdz + &
            dpdzfy(jy)*tavg_dt

        ! Mean velocity * velocity gradient product, ui_hat*dujdxk_hat
        tavg_specbudgy(jx,jy,jz)%uh_dudxh = tavg_specbudgy(jx,jy,jz)%uh_dudxh + &
            ufy(jy)*conjg(dudxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dudyh = tavg_specbudgy(jx,jy,jz)%uh_dudyh + &
            ufy(jy)*conjg(dudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dudzh = tavg_specbudgy(jx,jy,jz)%uh_dudzh + &
            ufy(jy)*conjg(dudzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dvdxh = tavg_specbudgy(jx,jy,jz)%uh_dvdxh + &
            ufy(jy)*conjg(dvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dvdyh = tavg_specbudgy(jx,jy,jz)%uh_dvdyh + &
            ufy(jy)*conjg(dvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dvdzh = tavg_specbudgy(jx,jy,jz)%uh_dvdzh + &
            ufy(jy)*conjg(dvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dwdxh = tavg_specbudgy(jx,jy,jz)%uh_dwdxh + &
            ufy(jy)*conjg(dwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dwdyh = tavg_specbudgy(jx,jy,jz)%uh_dwdyh + &
            ufy(jy)*conjg(dwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_dwdzh = tavg_specbudgy(jx,jy,jz)%uh_dwdzh + &
            ufy(jy)*conjg(dwdzfy(jy))*tavg_dt

        tavg_specbudgy(jx,jy,jz)%vh_dudxh = tavg_specbudgy(jx,jy,jz)%vh_dudxh + &
            vfy(jy)*conjg(dudxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dudyh = tavg_specbudgy(jx,jy,jz)%vh_dudyh + &
            vfy(jy)*conjg(dudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dudzh = tavg_specbudgy(jx,jy,jz)%vh_dudzh + &
            vfy(jy)*conjg(dudzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dvdxh = tavg_specbudgy(jx,jy,jz)%vh_dvdxh + &
            vfy(jy)*conjg(dvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dvdyh = tavg_specbudgy(jx,jy,jz)%vh_dvdyh + &
            vfy(jy)*conjg(dvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dvdzh = tavg_specbudgy(jx,jy,jz)%vh_dvdzh + &
            vfy(jy)*conjg(dvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dwdxh = tavg_specbudgy(jx,jy,jz)%vh_dwdxh + &
            vfy(jy)*conjg(dwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dwdyh = tavg_specbudgy(jx,jy,jz)%vh_dwdyh + &
            vfy(jy)*conjg(dwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_dwdzh = tavg_specbudgy(jx,jy,jz)%vh_dwdzh + &
            vfy(jy)*conjg(dwdzfy(jy))*tavg_dt

        tavg_specbudgy(jx,jy,jz)%wh_dudxh = tavg_specbudgy(jx,jy,jz)%wh_dudxh + &
            wfy(jy)*conjg(dudxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dudyh = tavg_specbudgy(jx,jy,jz)%wh_dudyh + &
            wfy(jy)*conjg(dudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dudzh = tavg_specbudgy(jx,jy,jz)%wh_dudzh + &
            wfy(jy)*conjg(dudzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dvdxh = tavg_specbudgy(jx,jy,jz)%wh_dvdxh + &
            wfy(jy)*conjg(dvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dvdyh = tavg_specbudgy(jx,jy,jz)%wh_dvdyh + &
            wfy(jy)*conjg(dvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dvdzh = tavg_specbudgy(jx,jy,jz)%wh_dvdzh + &
            wfy(jy)*conjg(dvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dwdxh = tavg_specbudgy(jx,jy,jz)%wh_dwdxh + &
            wfy(jy)*conjg(dwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dwdyh = tavg_specbudgy(jx,jy,jz)%wh_dwdyh + &
            wfy(jy)*conjg(dwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_dwdzh = tavg_specbudgy(jx,jy,jz)%wh_dwdzh + &
            wfy(jy)*conjg(dwdzfy(jy))*tavg_dt

        ! Mean velocity * velocity gradient product, uk*duidxk_hat
        tavg_specbudgy(jx,jy,jz)%ududxh = tavg_specbudgy(jx,jy,jz)%ududxh + &
            ududxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vdudyh = tavg_specbudgy(jx,jy,jz)%vdudyh + &
            vdudyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wdudzh = tavg_specbudgy(jx,jy,jz)%wdudzh + &
            wdudzfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%udvdxh = tavg_specbudgy(jx,jy,jz)%udvdxh + &
            udvdxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vdvdyh = tavg_specbudgy(jx,jy,jz)%vdvdyh + &
            vdvdyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wdvdzh = tavg_specbudgy(jx,jy,jz)%wdvdzh + &
            wdvdzfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%udwdxh = tavg_specbudgy(jx,jy,jz)%udwdxh + &
            udwdxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vdwdyh = tavg_specbudgy(jx,jy,jz)%vdwdyh + &
            vdwdyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wdwdzh = tavg_specbudgy(jx,jy,jz)%wdwdzh + &
            wdwdzfy(jy)*tavg_dt

        ! Mean vel * vel * vel gradient product, uj_hat*(uk*duidxk)_hat
        tavg_specbudgy(jx,jy,jz)%uh_ududxh = tavg_specbudgy(jx,jy,jz)%uh_ududxh + &
            ufy(jy)*conjg(ududxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_vdudyh = tavg_specbudgy(jx,jy,jz)%uh_vdudyh + &
            ufy(jy)*conjg(vdudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_wdudzh = tavg_specbudgy(jx,jy,jz)%uh_wdudzh + &
            ufy(jy)*conjg(wdudzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_ududxh = tavg_specbudgy(jx,jy,jz)%vh_ududxh + &
            vfy(jy)*conjg(ududxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_vdudyh = tavg_specbudgy(jx,jy,jz)%vh_vdudyh + &
            vfy(jy)*conjg(vdudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_wdudzh = tavg_specbudgy(jx,jy,jz)%vh_wdudzh + &
            vfy(jy)*conjg(wdudzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_ududxh = tavg_specbudgy(jx,jy,jz)%wh_ududxh + &
            wfy(jy)*conjg(ududxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_vdudyh = tavg_specbudgy(jx,jy,jz)%wh_vdudyh + &
            wfy(jy)*conjg(vdudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_wdudzh = tavg_specbudgy(jx,jy,jz)%wh_wdudzh + &
            wfy(jy)*conjg(wdudzfy(jy))*tavg_dt

        tavg_specbudgy(jx,jy,jz)%uh_udvdxh = tavg_specbudgy(jx,jy,jz)%uh_udvdxh + &
            ufy(jy)*conjg(udvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_vdvdyh = tavg_specbudgy(jx,jy,jz)%uh_vdvdyh + &
            ufy(jy)*conjg(vdvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_wdvdzh = tavg_specbudgy(jx,jy,jz)%uh_wdvdzh + &
            ufy(jy)*conjg(wdvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_udvdxh = tavg_specbudgy(jx,jy,jz)%vh_udvdxh + &
            vfy(jy)*conjg(udvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_vdvdyh = tavg_specbudgy(jx,jy,jz)%vh_vdvdyh + &
            vfy(jy)*conjg(vdvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_wdvdzh = tavg_specbudgy(jx,jy,jz)%vh_wdvdzh + &
            vfy(jy)*conjg(wdvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_udvdxh = tavg_specbudgy(jx,jy,jz)%wh_udvdxh + &
            wfy(jy)*conjg(udvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_vdvdyh = tavg_specbudgy(jx,jy,jz)%wh_vdvdyh + &
            wfy(jy)*conjg(vdvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_wdvdzh = tavg_specbudgy(jx,jy,jz)%wh_wdvdzh + &
            wfy(jy)*conjg(wdvdzfy(jy))*tavg_dt

        tavg_specbudgy(jx,jy,jz)%uh_udwdxh = tavg_specbudgy(jx,jy,jz)%uh_udwdxh + &
            ufy(jy)*conjg(udwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_vdwdyh = tavg_specbudgy(jx,jy,jz)%uh_vdwdyh + &
            ufy(jy)*conjg(vdwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%uh_wdwdzh = tavg_specbudgy(jx,jy,jz)%uh_wdwdzh + &
            ufy(jy)*conjg(wdwdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_udwdxh = tavg_specbudgy(jx,jy,jz)%vh_udwdxh + &
            vfy(jy)*conjg(udwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_vdwdyh = tavg_specbudgy(jx,jy,jz)%vh_vdwdyh + &
            vfy(jy)*conjg(vdwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vh_wdwdzh = tavg_specbudgy(jx,jy,jz)%vh_wdwdzh + &
            vfy(jy)*conjg(wdwdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_udwdxh = tavg_specbudgy(jx,jy,jz)%wh_udwdxh + &
            wfy(jy)*conjg(udwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_vdwdyh = tavg_specbudgy(jx,jy,jz)%wh_vdwdyh + &
            wfy(jy)*conjg(vdwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wh_wdwdzh = tavg_specbudgy(jx,jy,jz)%wh_wdwdzh + &
            wfy(jy)*conjg(wdwdzfy(jy))*tavg_dt

        ! Mean velocity * pressure gradient product, ui*dpdxj
        tavg_specbudgy(jx,jy,jz)%udpdx = tavg_specbudgy(jx,jy,jz)%udpdx + &
            ufy(jy)*conjg(dpdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%udpdy = tavg_specbudgy(jx,jy,jz)%udpdy + &
            ufy(jy)*conjg(dpdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%udpdz = tavg_specbudgy(jx,jy,jz)%udpdz + &
            ufy(jy)*conjg(dpdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vdpdx = tavg_specbudgy(jx,jy,jz)%vdpdx + &
            vfy(jy)*conjg(dpdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vdpdy = tavg_specbudgy(jx,jy,jz)%vdpdy + &
            vfy(jy)*conjg(dpdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vdpdz = tavg_specbudgy(jx,jy,jz)%vdpdz + &
            vfy(jy)*conjg(dpdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wdpdx = tavg_specbudgy(jx,jy,jz)%wdpdx + &
            wfy(jy)*conjg(dpdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wdpdy = tavg_specbudgy(jx,jy,jz)%wdpdy + &
            wfy(jy)*conjg(dpdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wdpdz = tavg_specbudgy(jx,jy,jz)%wdpdz + &
            wfy(jy)*conjg(dpdzfy(jy))*tavg_dt

        ! Mean pressure * velocity gradient product, p*duidxj
        tavg_specbudgy(jx,jy,jz)%pdudx = tavg_specbudgy(jx,jy,jz)%pdudx + &
            pfy(jy)*conjg(dudxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdudy = tavg_specbudgy(jx,jy,jz)%pdudy + &
            pfy(jy)*conjg(dudyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdudz = tavg_specbudgy(jx,jy,jz)%pdudz + &
            pfy(jy)*conjg(dudzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdvdx = tavg_specbudgy(jx,jy,jz)%pdvdx + &
            pfy(jy)*conjg(dvdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdvdy = tavg_specbudgy(jx,jy,jz)%pdvdy + &
            pfy(jy)*conjg(dvdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdvdz = tavg_specbudgy(jx,jy,jz)%pdvdz + &
            pfy(jy)*conjg(dvdzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdwdx = tavg_specbudgy(jx,jy,jz)%pdwdx + &
            pfy(jy)*conjg(dwdxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdwdy = tavg_specbudgy(jx,jy,jz)%pdwdy + &
            pfy(jy)*conjg(dwdyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%pdwdz = tavg_specbudgy(jx,jy,jz)%pdwdz + &
            pfy(jy)*conjg(dwdzfy(jy))*tavg_dt

        ! Mean Laplacian, nu*lap(uj)
        tavg_specbudgy(jx,jy,jz)%lapu = tavg_specbudgy(jx,jy,jz)%lapu + &
            divtxfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%lapv = tavg_specbudgy(jx,jy,jz)%lapv + &
            divtyfy(jy)*tavg_dt
        tavg_specbudgy(jx,jy,jz)%lapw = tavg_specbudgy(jx,jy,jz)%lapw + &
            divtzfy(jy)*tavg_dt

        ! Mean Vel-Laplacian, nu*ui*lap(uj)
        tavg_specbudgy(jx,jy,jz)%ulapu = tavg_specbudgy(jx,jy,jz)%ulapu + &
            ufy(jy)*conjg(divtxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%ulapv = tavg_specbudgy(jx,jy,jz)%ulapv + &
            ufy(jy)*conjg(divtyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%ulapw = tavg_specbudgy(jx,jy,jz)%ulapw + &
            ufy(jy)*conjg(divtzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vlapu = tavg_specbudgy(jx,jy,jz)%vlapu + &
            vfy(jy)*conjg(divtxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vlapv = tavg_specbudgy(jx,jy,jz)%vlapv + &
            vfy(jy)*conjg(divtyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%vlapw = tavg_specbudgy(jx,jy,jz)%vlapw + &
            vfy(jy)*conjg(divtzfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wlapu = tavg_specbudgy(jx,jy,jz)%wlapu + &
            wfy(jy)*conjg(divtxfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wlapv = tavg_specbudgy(jx,jy,jz)%wlapv + &
            wfy(jy)*conjg(divtyfy(jy))*tavg_dt
        tavg_specbudgy(jx,jy,jz)%wlapw = tavg_specbudgy(jx,jy,jz)%wlapw + &
            wfy(jy)*conjg(divtzfy(jy))*tavg_dt
#endif
    end do

end do 
end do

end subroutine tavg_turbspec_compute

#endif

!*****************************************************************************
subroutine tavg_finalize()
!*****************************************************************************
use grid_m
use stat_defs, only : tavg_t, tavg_total_time, tavg
use stat_defs, only : rs_compute, rs
use stat_defs, only : tavg_vort, vortrms_compute, vortrms
use param, only : write_endian
use param, only : ny,nz
use param, only : coord

#ifdef PPOUTPUT_SGS
use stat_defs, only : tavg_sgs
#endif

#ifdef PPOUTPUT_BUDGET
use stat_defs, only : tavg_budget, budget_compute, budget
#endif

#ifdef PPOUTPUT_TURBSPEC
use stat_defs, only : tavg_turbspecx, tavg_turbspecy, turbspecx, turbspecy
use stat_defs, only : turbspec_compute
#endif

#ifdef PPSCALARS
use stat_defs, only : tavg_scal_t, tavg_scal, rs_scal_compute, rs_scal
#ifdef PPOUTPUT_TURBSPEC
use stat_defs, only : tavg_scal_turbspecx, tavg_scal_turbspecy
use stat_defs, only : scal_turbspecx, scal_turbspecy
use stat_defs, only : scal_turbspec_compute
#endif
#ifdef PPOUTPUT_SPECBUDG
use stat_defs, only : tavg_scal_specbudgx, scal_specbudgx, scal_specbudgx_compute
#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
use stat_defs, only : tavg_specbudgx, specbudgx, specbudgx_compute
use stat_defs, only : tavg_specbudgy, specbudgy, specbudgy_compute
#endif

#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array,MPI_SYNC_DOWNUP
use param, only : ierr,comm
#endif

implicit none

#ifndef PPCGNS
character(64) :: bin_ext
#endif

character(64) :: fname_vel, fname_velw, fname_tau, fname_pres, fname_rs
character(64) :: fname_vort, fname_vortrms
! character(64) :: fname_f, fname_vel2

#ifdef PPOUTPUT_SGS
character(64) :: fname_sgs, fname_cs
#endif

#ifdef PPOUTPUT_BUDGET
character(64) :: fname_rxx, fname_ryy, fname_rzz, fname_rxy, fname_rxz, fname_ryz
#endif

#ifdef PPOUTPUT_TURBSPEC
character(64) :: fname_sxvel, fname_syvel, fname_sxvort, fname_syvort
#endif

#ifdef PPSCALARS
character(64) :: fname_scal
#ifdef PPOUTPUT_TURBSPEC
character(64) :: fname_sxscal, fname_syscal
#endif
#ifdef PPOUTPUT_SPECBUDG
character(64) :: fname_sxscal_budg
#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
character(64) :: fname_sxuu, fname_sxvv, fname_sxww, fname_sxuv, fname_sxuw, fname_sxvw
character(64) :: fname_syuu, fname_syvv, fname_syww, fname_syuv, fname_syuw, fname_syvw
#endif

integer :: i,j,k

real(rprec), pointer, dimension(:) :: x,y,z,zw

nullify(x,y,z,zw)

x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! Common file name
fname_vel = path     // 'output/veluv_avg'
fname_velw = path    // 'output/velw_avg'
! fname_vel2 = path    // 'output/vel2_avg'
fname_tau = path     // 'output/tau_avg'
! fname_f = path       // 'output/force_avg'
fname_pres = path    // 'output/pres_avg'
fname_rs = path      // 'output/rs'
fname_vort = path // 'output/vort_avg'
fname_vortrms = path // 'output/vortrms'
#ifdef PPOUTPUT_SGS
fname_cs = path // 'output/cs_opt2'
fname_sgs = path // 'output/sgs'
#endif
#ifdef PPOUTPUT_BUDGET
fname_rxx = path // 'output/rxx'
fname_ryy = path // 'output/ryy'
fname_rzz = path // 'output/rzz'
fname_rxy = path // 'output/rxy'
fname_rxz = path // 'output/rxz'
fname_ryz = path // 'output/ryz'
#endif
#ifdef PPOUTPUT_TURBSPEC
fname_sxvel = path // 'output/sxvel'
fname_syvel = path // 'output/syvel'
fname_sxvort = path // 'output/sxvort'
fname_syvort = path // 'output/syvort'
#endif

#ifdef PPSCALARS
fname_scal = path // 'output/scal'
#ifdef PPOUTPUT_TURBSPEC
fname_sxscal = path // 'output/sxscal'
fname_syscal = path // 'output/syscal'
#endif
#ifdef PPOUTPUT_SPECBUDG
fname_sxscal_budg = path // 'output/sxscal_budg'
#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
fname_sxuu = path // 'output/sxuu_budg'
fname_sxvv = path // 'output/sxvv_budg'
fname_sxww = path // 'output/sxww_budg'
fname_sxuv = path // 'output/sxuv_budg'
fname_sxuw = path // 'output/sxuw_budg'
fname_sxvw = path // 'output/sxvw_budg'

fname_syuu = path // 'output/syuu_budg'
fname_syvv = path // 'output/syvv_budg'
fname_syww = path // 'output/syww_budg'
fname_syuv = path // 'output/syuv_budg'
fname_syuw = path // 'output/syuw_budg'
fname_syvw = path // 'output/syvw_budg'
#endif

! CGNS
#ifdef PPCGNS
call string_concat(fname_vel, '.cgns')
call string_concat(fname_velw, '.cgns')
! call string_concat(fname_vel2, '.cgns')
call string_concat(fname_tau, '.cgns')
call string_concat(fname_pres, '.cgns')
! call string_concat(fname_f, '.cgns')
call string_concat(fname_rs, '.cgns')
call string_concat(fname_vort, '.cgns')
call string_concat(fname_vortrms, '.cgns')
#ifdef PPOUTPUT_SGS
call string_concat(fname_cs, '.cgns')
call string_concat(fname_sgs, '.cgns')
#endif
#ifdef PPOUTPUT_BUDGET
call string_concat(fname_rxx, '.cgns')
call string_concat(fname_ryy, '.cgns')
call string_concat(fname_rzz, '.cgns')
call string_concat(fname_rxy, '.cgns')
call string_concat(fname_rxz, '.cgns')
call string_concat(fname_ryz, '.cgns')
#endif
#ifdef PPOUTPUT_TURBSPEC
call string_concat(fname_sxvel, '.cgns')
call string_concat(fname_syvel, '.cgns')
call string_concat(fname_sxvort, '.cgns')
call string_concat(fname_syvort, '.cgns')
#endif

#ifdef PPSCALARS
call string_concat(fname_scal, '.cgns')
#ifdef PPOUTPUT_TURBSPEC
call string_concat(fname_sxscal, '.cgns')
call string_concat(fname_syscal, '.cgns')
#endif
#ifdef PPOUTPUT_SPECBUDG
call string_concat(fname_sxscal_budg, '.cgns')
#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
call string_concat(fname_sxuu, '.cgns')
call string_concat(fname_sxvv, '.cgns')
call string_concat(fname_sxww, '.cgns')
call string_concat(fname_sxuv, '.cgns')
call string_concat(fname_sxuw, '.cgns')
call string_concat(fname_sxvw, '.cgns')

call string_concat(fname_syuu, '.cgns')
call string_concat(fname_syvv, '.cgns')
call string_concat(fname_syww, '.cgns')
call string_concat(fname_syuv, '.cgns')
call string_concat(fname_syuw, '.cgns')
call string_concat(fname_syvw, '.cgns')
#endif

! Binary
#else
#ifdef PPMPI
call string_splice(bin_ext, '.c', coord, '.bin')
#else
bin_ext = '.bin'
#endif
call string_concat(fname_vel, bin_ext)
call string_concat(fname_velw, bin_ext)
! call string_concat(fname_vel2, bin_ext)
call string_concat(fname_tau, bin_ext)
call string_concat(fname_pres, bin_ext)
! call string_concat(fname_f, bin_ext)
call string_concat(fname_rs, bin_ext)
call string_concat(fname_vort, bin_ext)
call string_concat(fname_vortrms, bin_ext)
#ifdef PPOUTPUT_SGS
call string_concat(fname_cs, bin_ext)
call string_concat(fname_sgs, bin_ext)
#endif
#ifdef PPOUTPUT_BUDGET
call string_concat(fname_rxx, bin_ext)
call string_concat(fname_ryy, bin_ext)
call string_concat(fname_rzz, bin_ext)
call string_concat(fname_rxy, bin_ext)
call string_concat(fname_rxz, bin_ext)
call string_concat(fname_ryz, bin_ext)
#endif
#ifdef PPOUTPUT_TURBSPEC
call string_concat(fname_sxvel, bin_ext)
call string_concat(fname_syvel, bin_ext)
call string_concat(fname_sxvort, bin_ext)
call string_concat(fname_syvort, bin_ext)
#endif
#ifdef PPSCALARS
call string_concat(fname_scal, bin_ext)
#ifdef PPOUTPUT_TURBSPEC
call string_concat(fname_sxscal, bin_ext)
call string_concat(fname_syscal, bin_ext)
#endif
#ifdef PPOUTPUT_SPECBUDG
call string_concat(fname_sxscal_budg, bin_ext)
#endif
#endif
#ifdef PPOUTPUT_SPECBUDG
call string_concat(fname_sxuu, bin_ext)
call string_concat(fname_sxvv, bin_ext)
call string_concat(fname_sxww, bin_ext)
call string_concat(fname_sxuv, bin_ext)
call string_concat(fname_sxuw, bin_ext)
call string_concat(fname_sxvw, bin_ext)

call string_concat(fname_syuu, bin_ext)
call string_concat(fname_syvv, bin_ext)
call string_concat(fname_syww, bin_ext)
call string_concat(fname_syuv, bin_ext)
call string_concat(fname_syuw, bin_ext)
call string_concat(fname_syvw, bin_ext)
#endif
#endif

! Final checkpoint all restart data
call tavg_checkpoint()

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

!  Perform time averaging operation
do k = jzmin, jzmax
do j = 1, Ny
do i = 1, Nx
    tavg(i,j,k) % u    = tavg(i,j,k) % u    / tavg_total_time
    tavg(i,j,k) % v    = tavg(i,j,k) % v    / tavg_total_time
    tavg(i,j,k) % w_uv = tavg(i,j,k) % w_uv / tavg_total_time
    tavg(i,j,k) % u_w  = tavg(i,j,k) % u_w  / tavg_total_time
    tavg(i,j,k) % v_w  = tavg(i,j,k) % v_w  / tavg_total_time
    tavg(i,j,k) % w    = tavg(i,j,k) % w    / tavg_total_time
    tavg(i,j,k) % u2   = tavg(i,j,k) % u2   / tavg_total_time
    tavg(i,j,k) % v2   = tavg(i,j,k) % v2   / tavg_total_time
    tavg(i,j,k) % w2   = tavg(i,j,k) % w2   / tavg_total_time
    tavg(i,j,k) % uv   = tavg(i,j,k) % uv   / tavg_total_time
    tavg(i,j,k) % uw   = tavg(i,j,k) % uw   / tavg_total_time
    tavg(i,j,k) % vw   = tavg(i,j,k) % vw   / tavg_total_time
    tavg(i,j,k) % txx  = tavg(i,j,k) % txx  / tavg_total_time
    tavg(i,j,k) % tyy  = tavg(i,j,k) % tyy  / tavg_total_time
    tavg(i,j,k) % tzz  = tavg(i,j,k) % tzz  / tavg_total_time
    tavg(i,j,k) % txy  = tavg(i,j,k) % txy  / tavg_total_time
    tavg(i,j,k) % txz  = tavg(i,j,k) % txz  / tavg_total_time
    tavg(i,j,k) % tyz  = tavg(i,j,k) % tyz  / tavg_total_time
    tavg(i,j,k) % fx   = tavg(i,j,k) % fx   / tavg_total_time
    tavg(i,j,k) % fy   = tavg(i,j,k) % fy   / tavg_total_time
    tavg(i,j,k) % fz   = tavg(i,j,k) % fz   / tavg_total_time

#ifdef PPSCALARS
    tavg_scal(i,j,k) % theta  = tavg_scal(i,j,k) % theta / tavg_total_time
    tavg_scal(i,j,k) % theta2 = tavg_scal(i,j,k) % theta2 / tavg_total_time
    tavg_scal(i,j,k) % utheta = tavg_scal(i,j,k) % utheta / tavg_total_time
    tavg_scal(i,j,k) % vtheta = tavg_scal(i,j,k) % vtheta / tavg_total_time
    tavg_scal(i,j,k) % wtheta = tavg_scal(i,j,k) % wtheta / tavg_total_time
    tavg_scal(i,j,k) % pi_z = tavg_scal(i,j,k) % pi_z / tavg_total_time
#endif

    tavg_vort(i,j,k)%vortx = tavg_vort(i,j,k)%vortx / tavg_total_time
    tavg_vort(i,j,k)%vorty = tavg_vort(i,j,k)%vorty / tavg_total_time
    tavg_vort(i,j,k)%vortz = tavg_vort(i,j,k)%vortz / tavg_total_time
    tavg_vort(i,j,k)%vortx2 = tavg_vort(i,j,k)%vortx2 / tavg_total_time
    tavg_vort(i,j,k)%vorty2 = tavg_vort(i,j,k)%vorty2 / tavg_total_time
    tavg_vort(i,j,k)%vortz2 = tavg_vort(i,j,k)%vortz2 / tavg_total_time

end do
end do
end do

#ifdef PPOUTPUT_SGS
do k = jzmin, jzmax
do j = 1, Ny
do i = 1, Nx
    tavg_sgs(i,j,k) % cs_opt2 = tavg_sgs(i,j,k) % cs_opt2 / tavg_total_time
    ! tavg_sgs(i,j,k) % Tn       = tavg_sgs(i,j,k) % Tn / tavg_total_time
    tavg_sgs(i,j,k) % Nu_t     = tavg_sgs(i,j,k) % Nu_t / tavg_total_time
    ! tavg_sgs(i,j,k) % F_LM     = tavg_sgs(i,j,k) % F_LM / tavg_total_time
    ! tavg_sgs(i,j,k) % F_MM     = tavg_sgs(i,j,k) % F_MM / tavg_total_time
    ! tavg_sgs(i,j,k) % F_QN     = tavg_sgs(i,j,k) % F_QN / tavg_total_time
    ! tavg_sgs(i,j,k) % F_NN     = tavg_sgs(i,j,k) % F_NN / tavg_total_time
    ! tavg_sgs(i,j,k) % ee_now   = tavg_sgs(i,j,k) % ee_now / tavg_total_time
! #ifdef PPDYN_TN
    ! tavg_sgs(i,j,k) % F_ee2    = tavg_sgs(i,j,k) % F_ee2 / tavg_total_time
    ! tavg_sgs(i,j,k) % F_deedt2 = tavg_sgs(i,j,k) % F_deedt2 / tavg_total_time
! #endif
end do
end do
end do
#endif

#ifdef PPOUTPUT_BUDGET
do k = jzmin, jzmax
do j = 1, Ny
do i = 1, Nx
    ! Mean pressure on w-grid
    tavg_budget(i,j,k) % p = tavg_budget(i,j,k) % p / tavg_total_time

    ! Mean velocity-velocity product, ui*uj
    tavg_budget(i,j,k) % uu = tavg_budget(i,j,k) % uu / tavg_total_time
    tavg_budget(i,j,k) % vv = tavg_budget(i,j,k) % vv / tavg_total_time
    tavg_budget(i,j,k) % ww = tavg_budget(i,j,k) % ww / tavg_total_time
    tavg_budget(i,j,k) % uv = tavg_budget(i,j,k) % uv / tavg_total_time
    tavg_budget(i,j,k) % uw = tavg_budget(i,j,k) % uw / tavg_total_time
    tavg_budget(i,j,k) % vw = tavg_budget(i,j,k) % vw / tavg_total_time

    ! Mean velocity gradients, duidxj
    tavg_budget(i,j,k) % dudx = tavg_budget(i,j,k) % dudx / tavg_total_time
    tavg_budget(i,j,k) % dudy = tavg_budget(i,j,k) % dudy / tavg_total_time
    tavg_budget(i,j,k) % dudz = tavg_budget(i,j,k) % dudz / tavg_total_time
    tavg_budget(i,j,k) % dvdx = tavg_budget(i,j,k) % dvdx / tavg_total_time
    tavg_budget(i,j,k) % dvdy = tavg_budget(i,j,k) % dvdy / tavg_total_time
    tavg_budget(i,j,k) % dvdz = tavg_budget(i,j,k) % dvdz / tavg_total_time
    tavg_budget(i,j,k) % dwdx = tavg_budget(i,j,k) % dwdx / tavg_total_time
    tavg_budget(i,j,k) % dwdy = tavg_budget(i,j,k) % dwdy / tavg_total_time
    tavg_budget(i,j,k) % dwdz = tavg_budget(i,j,k) % dwdz / tavg_total_time

    ! Mean pressure gradients, dpdxi
    tavg_budget(i,j,k) % dpdx = tavg_budget(i,j,k) % dpdx / tavg_total_time
    tavg_budget(i,j,k) % dpdy = tavg_budget(i,j,k) % dpdy / tavg_total_time
    tavg_budget(i,j,k) % dpdz = tavg_budget(i,j,k) % dpdz / tavg_total_time

    ! Mean vel-velGrad product, ui*dujdxk
    tavg_budget(i,j,k) % ududx = tavg_budget(i,j,k) % ududx / tavg_total_time
    tavg_budget(i,j,k) % ududy = tavg_budget(i,j,k) % ududy / tavg_total_time
    tavg_budget(i,j,k) % ududz = tavg_budget(i,j,k) % ududz / tavg_total_time
    tavg_budget(i,j,k) % udvdx = tavg_budget(i,j,k) % udvdx / tavg_total_time
    tavg_budget(i,j,k) % udvdy = tavg_budget(i,j,k) % udvdy / tavg_total_time
    tavg_budget(i,j,k) % udvdz = tavg_budget(i,j,k) % udvdz / tavg_total_time
    tavg_budget(i,j,k) % udwdx = tavg_budget(i,j,k) % udwdx / tavg_total_time
    tavg_budget(i,j,k) % udwdy = tavg_budget(i,j,k) % udwdy / tavg_total_time
    tavg_budget(i,j,k) % udwdz = tavg_budget(i,j,k) % udwdz / tavg_total_time

    tavg_budget(i,j,k) % vdudx = tavg_budget(i,j,k) % vdudx / tavg_total_time
    tavg_budget(i,j,k) % vdudy = tavg_budget(i,j,k) % vdudy / tavg_total_time
    tavg_budget(i,j,k) % vdudz = tavg_budget(i,j,k) % vdudz / tavg_total_time
    tavg_budget(i,j,k) % vdvdx = tavg_budget(i,j,k) % vdvdx / tavg_total_time
    tavg_budget(i,j,k) % vdvdy = tavg_budget(i,j,k) % vdvdy / tavg_total_time
    tavg_budget(i,j,k) % vdvdz = tavg_budget(i,j,k) % vdvdz / tavg_total_time
    tavg_budget(i,j,k) % vdwdx = tavg_budget(i,j,k) % vdwdx / tavg_total_time
    tavg_budget(i,j,k) % vdwdy = tavg_budget(i,j,k) % vdwdy / tavg_total_time
    tavg_budget(i,j,k) % vdwdz = tavg_budget(i,j,k) % vdwdz / tavg_total_time

    tavg_budget(i,j,k) % wdudx = tavg_budget(i,j,k) % wdudx / tavg_total_time
    tavg_budget(i,j,k) % wdudy = tavg_budget(i,j,k) % wdudy / tavg_total_time
    tavg_budget(i,j,k) % wdudz = tavg_budget(i,j,k) % wdudz / tavg_total_time
    tavg_budget(i,j,k) % wdvdx = tavg_budget(i,j,k) % wdvdx / tavg_total_time
    tavg_budget(i,j,k) % wdvdy = tavg_budget(i,j,k) % wdvdy / tavg_total_time
    tavg_budget(i,j,k) % wdvdz = tavg_budget(i,j,k) % wdvdz / tavg_total_time
    tavg_budget(i,j,k) % wdwdx = tavg_budget(i,j,k) % wdwdx / tavg_total_time
    tavg_budget(i,j,k) % wdwdy = tavg_budget(i,j,k) % wdwdy / tavg_total_time
    tavg_budget(i,j,k) % wdwdz = tavg_budget(i,j,k) % wdwdz / tavg_total_time

    ! Mean vel-vel-velGrad product, ui*uk*dujdxk
    tavg_budget(i,j,k) % uududx = tavg_budget(i,j,k) % uududx / tavg_total_time
    tavg_budget(i,j,k) % uvdudy = tavg_budget(i,j,k) % uvdudy / tavg_total_time
    tavg_budget(i,j,k) % uwdudz = tavg_budget(i,j,k) % uwdudz / tavg_total_time
    tavg_budget(i,j,k) % uudvdx = tavg_budget(i,j,k) % uudvdx / tavg_total_time
    tavg_budget(i,j,k) % uvdvdy = tavg_budget(i,j,k) % uvdvdy / tavg_total_time
    tavg_budget(i,j,k) % uwdvdz = tavg_budget(i,j,k) % uwdvdz / tavg_total_time
    tavg_budget(i,j,k) % uudwdx = tavg_budget(i,j,k) % uudwdx / tavg_total_time
    tavg_budget(i,j,k) % uvdwdy = tavg_budget(i,j,k) % uvdwdy / tavg_total_time
    tavg_budget(i,j,k) % uwdwdz = tavg_budget(i,j,k) % uwdwdz / tavg_total_time

    tavg_budget(i,j,k) % vududx = tavg_budget(i,j,k) % vududx / tavg_total_time
    tavg_budget(i,j,k) % vvdudy = tavg_budget(i,j,k) % vvdudy / tavg_total_time
    tavg_budget(i,j,k) % vwdudz = tavg_budget(i,j,k) % vwdudz / tavg_total_time
    tavg_budget(i,j,k) % vudvdx = tavg_budget(i,j,k) % vudvdx / tavg_total_time
    tavg_budget(i,j,k) % vvdvdy = tavg_budget(i,j,k) % vvdvdy / tavg_total_time
    tavg_budget(i,j,k) % vwdvdz = tavg_budget(i,j,k) % vwdvdz / tavg_total_time
    tavg_budget(i,j,k) % vudwdx = tavg_budget(i,j,k) % vudwdx / tavg_total_time
    tavg_budget(i,j,k) % vvdwdy = tavg_budget(i,j,k) % vvdwdy / tavg_total_time
    tavg_budget(i,j,k) % vwdwdz = tavg_budget(i,j,k) % vwdwdz / tavg_total_time

    tavg_budget(i,j,k) % wududx = tavg_budget(i,j,k) % wududx / tavg_total_time
    tavg_budget(i,j,k) % wvdudy = tavg_budget(i,j,k) % wvdudy / tavg_total_time
    tavg_budget(i,j,k) % wwdudz = tavg_budget(i,j,k) % wwdudz / tavg_total_time
    tavg_budget(i,j,k) % wudvdx = tavg_budget(i,j,k) % wudvdx / tavg_total_time
    tavg_budget(i,j,k) % wvdvdy = tavg_budget(i,j,k) % wvdvdy / tavg_total_time
    tavg_budget(i,j,k) % wwdvdz = tavg_budget(i,j,k) % wwdvdz / tavg_total_time
    tavg_budget(i,j,k) % wudwdx = tavg_budget(i,j,k) % wudwdx / tavg_total_time
    tavg_budget(i,j,k) % wvdwdy = tavg_budget(i,j,k) % wvdwdy / tavg_total_time
    tavg_budget(i,j,k) % wwdwdz = tavg_budget(i,j,k) % wwdwdz / tavg_total_time

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
    tavg_budget(i,j,k) % uxux = tavg_budget(i,j,k) % uxux / tavg_total_time
    tavg_budget(i,j,k) % uyuy = tavg_budget(i,j,k) % uyuy / tavg_total_time
    tavg_budget(i,j,k) % uzuz = tavg_budget(i,j,k) % uzuz / tavg_total_time
    tavg_budget(i,j,k) % vxvx = tavg_budget(i,j,k) % vxvx / tavg_total_time
    tavg_budget(i,j,k) % vyvy = tavg_budget(i,j,k) % vyvy / tavg_total_time
    tavg_budget(i,j,k) % vzvz = tavg_budget(i,j,k) % vzvz / tavg_total_time
    tavg_budget(i,j,k) % wxwx = tavg_budget(i,j,k) % wxwx / tavg_total_time
    tavg_budget(i,j,k) % wywy = tavg_budget(i,j,k) % wywy / tavg_total_time
    tavg_budget(i,j,k) % wzwz = tavg_budget(i,j,k) % wzwz / tavg_total_time

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
    tavg_budget(i,j,k) % uxvx = tavg_budget(i,j,k) % uxvx / tavg_total_time
    tavg_budget(i,j,k) % uyvy = tavg_budget(i,j,k) % uyvy / tavg_total_time
    tavg_budget(i,j,k) % uzvz = tavg_budget(i,j,k) % uzvz / tavg_total_time
    tavg_budget(i,j,k) % uxwx = tavg_budget(i,j,k) % uxwx / tavg_total_time
    tavg_budget(i,j,k) % uywy = tavg_budget(i,j,k) % uywy / tavg_total_time
    tavg_budget(i,j,k) % uzwz = tavg_budget(i,j,k) % uzwz / tavg_total_time
    tavg_budget(i,j,k) % vxwx = tavg_budget(i,j,k) % vxwx / tavg_total_time
    tavg_budget(i,j,k) % vywy = tavg_budget(i,j,k) % vywy / tavg_total_time
    tavg_budget(i,j,k) % vzwz = tavg_budget(i,j,k) % vzwz / tavg_total_time

    ! duidxj*dujdxk, i /= j
    ! tavg_budget(i,j,k) % uyvx = tavg_budget(i,j,k) % uyvx / tavg_total_time
    ! tavg_budget(i,j,k) % uzwx = tavg_budget(i,j,k) % uzwx / tavg_total_time
    ! tavg_budget(i,j,k) % vzwy = tavg_budget(i,j,k) % vzwy / tavg_total_time

    ! Mean vel-presGrad product, ui*dpdxj
    tavg_budget(i,j,k) % udpdx = tavg_budget(i,j,k) % udpdx / tavg_total_time
    tavg_budget(i,j,k) % udpdy = tavg_budget(i,j,k) % udpdy / tavg_total_time
    tavg_budget(i,j,k) % udpdz = tavg_budget(i,j,k) % udpdz / tavg_total_time
    tavg_budget(i,j,k) % vdpdx = tavg_budget(i,j,k) % vdpdx / tavg_total_time
    tavg_budget(i,j,k) % vdpdy = tavg_budget(i,j,k) % vdpdy / tavg_total_time
    tavg_budget(i,j,k) % vdpdz = tavg_budget(i,j,k) % vdpdz / tavg_total_time
    tavg_budget(i,j,k) % wdpdx = tavg_budget(i,j,k) % wdpdx / tavg_total_time
    tavg_budget(i,j,k) % wdpdy = tavg_budget(i,j,k) % wdpdy / tavg_total_time
    tavg_budget(i,j,k) % wdpdz = tavg_budget(i,j,k) % wdpdz / tavg_total_time

    ! Mean pres-velGrad product, p*duidxj
    tavg_budget(i,j,k) % pdudx = tavg_budget(i,j,k) % pdudx / tavg_total_time
    tavg_budget(i,j,k) % pdudy = tavg_budget(i,j,k) % pdudy / tavg_total_time
    tavg_budget(i,j,k) % pdudz = tavg_budget(i,j,k) % pdudz / tavg_total_time
    tavg_budget(i,j,k) % pdvdx = tavg_budget(i,j,k) % pdvdx / tavg_total_time
    tavg_budget(i,j,k) % pdvdy = tavg_budget(i,j,k) % pdvdy / tavg_total_time
    tavg_budget(i,j,k) % pdvdz = tavg_budget(i,j,k) % pdvdz / tavg_total_time
    tavg_budget(i,j,k) % pdwdx = tavg_budget(i,j,k) % pdwdx / tavg_total_time
    tavg_budget(i,j,k) % pdwdy = tavg_budget(i,j,k) % pdwdy / tavg_total_time
    tavg_budget(i,j,k) % pdwdz = tavg_budget(i,j,k) % pdwdz / tavg_total_time

    ! Mean Laplacian, nu*lap(uj)
    tavg_budget(i,j,k) % lapu = tavg_budget(i,j,k) % lapu / tavg_total_time
    tavg_budget(i,j,k) % lapv = tavg_budget(i,j,k) % lapv / tavg_total_time
    tavg_budget(i,j,k) % lapw = tavg_budget(i,j,k) % lapw / tavg_total_time

    ! Mean Vel-Laplacian, nu*ui*lap(uj)
    tavg_budget(i,j,k) % ulapu = tavg_budget(i,j,k) % ulapu / tavg_total_time
    tavg_budget(i,j,k) % ulapv = tavg_budget(i,j,k) % ulapv / tavg_total_time
    tavg_budget(i,j,k) % ulapw = tavg_budget(i,j,k) % ulapw / tavg_total_time
    tavg_budget(i,j,k) % vlapu = tavg_budget(i,j,k) % vlapu / tavg_total_time
    tavg_budget(i,j,k) % vlapv = tavg_budget(i,j,k) % vlapv / tavg_total_time
    tavg_budget(i,j,k) % vlapw = tavg_budget(i,j,k) % vlapw / tavg_total_time
    tavg_budget(i,j,k) % wlapu = tavg_budget(i,j,k) % wlapu / tavg_total_time
    tavg_budget(i,j,k) % wlapv = tavg_budget(i,j,k) % wlapv / tavg_total_time
    tavg_budget(i,j,k) % wlapw = tavg_budget(i,j,k) % wlapw / tavg_total_time

end do
end do
end do
#endif

#ifdef PPOUTPUT_TURBSPEC
do k = jzmin, jzmax

do j = 1, Ny
do i = 1, Nx/2 + 1
    tavg_turbspecx(i,j,k) % uf = tavg_turbspecx(i,j,k) % uf / tavg_total_time
    tavg_turbspecx(i,j,k) % vf = tavg_turbspecx(i,j,k) % vf / tavg_total_time
    tavg_turbspecx(i,j,k) % wf = tavg_turbspecx(i,j,k) % wf / tavg_total_time

    tavg_turbspecx(i,j,k) % uu = tavg_turbspecx(i,j,k) % uu / tavg_total_time
    tavg_turbspecx(i,j,k) % vv = tavg_turbspecx(i,j,k) % vv / tavg_total_time
    tavg_turbspecx(i,j,k) % ww = tavg_turbspecx(i,j,k) % ww / tavg_total_time
    !tavg_turbspecx(i,j,k) % vel2 = tavg_turbspecx(i,j,k) % vel2 / tavg_total_time

    tavg_turbspecx(i,j,k) % uv = tavg_turbspecx(i,j,k) % uv / tavg_total_time
    tavg_turbspecx(i,j,k) % uw = tavg_turbspecx(i,j,k) % uw / tavg_total_time
    tavg_turbspecx(i,j,k) % vw = tavg_turbspecx(i,j,k) % vw / tavg_total_time

    tavg_turbspecx(i,j,k) % vortxf = tavg_turbspecx(i,j,k) % vortxf / tavg_total_time
    tavg_turbspecx(i,j,k) % vortyf = tavg_turbspecx(i,j,k) % vortyf / tavg_total_time
    tavg_turbspecx(i,j,k) % vortzf = tavg_turbspecx(i,j,k) % vortzf / tavg_total_time

    tavg_turbspecx(i,j,k) % vortx2 = tavg_turbspecx(i,j,k) % vortx2 / tavg_total_time
    tavg_turbspecx(i,j,k) % vorty2 = tavg_turbspecx(i,j,k) % vorty2 / tavg_total_time
    tavg_turbspecx(i,j,k) % vortz2 = tavg_turbspecx(i,j,k) % vortz2 / tavg_total_time
    !tavg_turbspecx(i,j,k) % vort2 = tavg_turbspecx(i,j,k) % vort2 / tavg_total_time

#ifdef PPSCALARS
    tavg_scal_turbspecx(i,j,k) % thetaf = tavg_scal_turbspecx(i,j,k) % thetaf / tavg_total_time

    tavg_scal_turbspecx(i,j,k) % theta2 = tavg_scal_turbspecx(i,j,k) % theta2 / tavg_total_time
    tavg_scal_turbspecx(i,j,k) % utheta = tavg_scal_turbspecx(i,j,k) % utheta / tavg_total_time
    tavg_scal_turbspecx(i,j,k) % vtheta = tavg_scal_turbspecx(i,j,k) % vtheta / tavg_total_time
    tavg_scal_turbspecx(i,j,k) % wtheta = tavg_scal_turbspecx(i,j,k) % wtheta / tavg_total_time

#ifdef PPOUTPUT_SPECBUDG
    ! Scalar gradient, dTdxj_hat
    tavg_scal_specbudgx(i,j,k) % dTdx = tavg_scal_specbudgx(i,j,k) % dTdx / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % dTdy = tavg_scal_specbudgx(i,j,k) % dTdy / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % dTdz = tavg_scal_specbudgx(i,j,k) % dTdz / tavg_total_time

    ! Scalar-Scalar gradient, T_hat*dTdxj_hat
    tavg_scal_specbudgx(i,j,k) % Th_dTdxh = tavg_scal_specbudgx(i,j,k) % Th_dTdxh / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % Th_dTdyh = tavg_scal_specbudgx(i,j,k) % Th_dTdyh / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % Th_dTdzh = tavg_scal_specbudgx(i,j,k) % Th_dTdzh / tavg_total_time

    ! Velocity-Scalar Gradient, (uj*dTdxj)_hat
    tavg_scal_specbudgx(i,j,k) % udTdxh = tavg_scal_specbudgx(i,j,k) % udTdxh / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % vdTdyh = tavg_scal_specbudgx(i,j,k) % vdTdyh / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % wdTdzh = tavg_scal_specbudgx(i,j,k) % wdTdzh / tavg_total_time

    ! Scalar-Velocity-Scalar Gradient, T_hat*((uj*dTdxj)_hat)
    tavg_scal_specbudgx(i,j,k) % Th_udTdxh = tavg_scal_specbudgx(i,j,k) % Th_udTdxh / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % Th_vdTdyh = tavg_scal_specbudgx(i,j,k) % Th_vdTdyh / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % Th_wdTdzh = tavg_scal_specbudgx(i,j,k) % Th_wdTdzh / tavg_total_time

    ! Scalar gradient-scalar gradient, dTdxj_hat*dTdxj_hat
    tavg_scal_specbudgx(i,j,k) % TxTx = tavg_scal_specbudgx(i,j,k) % TxTx / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % TyTy = tavg_scal_specbudgx(i,j,k) % TyTy / tavg_total_time
    tavg_scal_specbudgx(i,j,k) % TzTz = tavg_scal_specbudgx(i,j,k) % TzTz / tavg_total_time

    ! Laplacian, (nu/Pr)*lap(theta)
    tavg_scal_specbudgx(i,j,k) % lapT = tavg_scal_specbudgx(i,j,k) % lapT / tavg_total_time

    ! Scalar-Laplacian, theta*(nu/Pr)*lap(theta)
    tavg_scal_specbudgx(i,j,k) % TlapT = tavg_scal_specbudgx(i,j,k) % TlapT / tavg_total_time

#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
    ! Mean pressure
    tavg_specbudgx(i,j,k) % p = tavg_specbudgx(i,j,k) % p / tavg_total_time

    ! Mean velocity gradient
    tavg_specbudgx(i,j,k) % dudx = tavg_specbudgx(i,j,k) % dudx / tavg_total_time
    tavg_specbudgx(i,j,k) % dudy = tavg_specbudgx(i,j,k) % dudy / tavg_total_time
    tavg_specbudgx(i,j,k) % dudz = tavg_specbudgx(i,j,k) % dudz / tavg_total_time
    tavg_specbudgx(i,j,k) % dvdx = tavg_specbudgx(i,j,k) % dvdx / tavg_total_time
    tavg_specbudgx(i,j,k) % dvdy = tavg_specbudgx(i,j,k) % dvdy / tavg_total_time
    tavg_specbudgx(i,j,k) % dvdz = tavg_specbudgx(i,j,k) % dvdz / tavg_total_time
    tavg_specbudgx(i,j,k) % dwdx = tavg_specbudgx(i,j,k) % dwdx / tavg_total_time
    tavg_specbudgx(i,j,k) % dwdy = tavg_specbudgx(i,j,k) % dwdy / tavg_total_time
    tavg_specbudgx(i,j,k) % dwdz = tavg_specbudgx(i,j,k) % dwdz / tavg_total_time

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
    tavg_specbudgx(i,j,k) % uxux = tavg_specbudgx(i,j,k) % uxux / tavg_total_time
    tavg_specbudgx(i,j,k) % uyuy = tavg_specbudgx(i,j,k) % uyuy / tavg_total_time
    tavg_specbudgx(i,j,k) % uzuz = tavg_specbudgx(i,j,k) % uzuz / tavg_total_time
    tavg_specbudgx(i,j,k) % vxvx = tavg_specbudgx(i,j,k) % vxvx / tavg_total_time
    tavg_specbudgx(i,j,k) % vyvy = tavg_specbudgx(i,j,k) % vyvy / tavg_total_time
    tavg_specbudgx(i,j,k) % vzvz = tavg_specbudgx(i,j,k) % vzvz / tavg_total_time
    tavg_specbudgx(i,j,k) % wxwx = tavg_specbudgx(i,j,k) % wxwx / tavg_total_time
    tavg_specbudgx(i,j,k) % wywy = tavg_specbudgx(i,j,k) % wywy / tavg_total_time
    tavg_specbudgx(i,j,k) % wzwz = tavg_specbudgx(i,j,k) % wzwz / tavg_total_time

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
    tavg_specbudgx(i,j,k) % uxvx = tavg_specbudgx(i,j,k) % uxvx / tavg_total_time
    tavg_specbudgx(i,j,k) % uyvy = tavg_specbudgx(i,j,k) % uyvy / tavg_total_time
    tavg_specbudgx(i,j,k) % uzvz = tavg_specbudgx(i,j,k) % uzvz / tavg_total_time
    tavg_specbudgx(i,j,k) % uxwx = tavg_specbudgx(i,j,k) % uxwx / tavg_total_time
    tavg_specbudgx(i,j,k) % uywy = tavg_specbudgx(i,j,k) % uywy / tavg_total_time
    tavg_specbudgx(i,j,k) % uzwz = tavg_specbudgx(i,j,k) % uzwz / tavg_total_time
    tavg_specbudgx(i,j,k) % vxwx = tavg_specbudgx(i,j,k) % vxwx / tavg_total_time
    tavg_specbudgx(i,j,k) % vywy = tavg_specbudgx(i,j,k) % vywy / tavg_total_time
    tavg_specbudgx(i,j,k) % vzwz = tavg_specbudgx(i,j,k) % vzwz / tavg_total_time

    ! Mean pressure gradients, dpdxi
    tavg_specbudgx(i,j,k) % dpdx = tavg_specbudgx(i,j,k) % dpdx / tavg_total_time
    tavg_specbudgx(i,j,k) % dpdy = tavg_specbudgx(i,j,k) % dpdy / tavg_total_time
    tavg_specbudgx(i,j,k) % dpdz = tavg_specbudgx(i,j,k) % dpdz / tavg_total_time

    ! Mean velocity * velocity gradient product, ui_hat*dujdxk_hat
    tavg_specbudgx(i,j,k) % uh_dudxh = tavg_specbudgx(i,j,k) % uh_dudxh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dudyh = tavg_specbudgx(i,j,k) % uh_dudyh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dudzh = tavg_specbudgx(i,j,k) % uh_dudzh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dvdxh = tavg_specbudgx(i,j,k) % uh_dvdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dvdyh = tavg_specbudgx(i,j,k) % uh_dvdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dvdzh = tavg_specbudgx(i,j,k) % uh_dvdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dwdxh = tavg_specbudgx(i,j,k) % uh_dwdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dwdyh = tavg_specbudgx(i,j,k) % uh_dwdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_dwdzh = tavg_specbudgx(i,j,k) % uh_dwdzh / tavg_total_time

    tavg_specbudgx(i,j,k) % vh_dudxh = tavg_specbudgx(i,j,k) % vh_dudxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dudyh = tavg_specbudgx(i,j,k) % vh_dudyh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dudzh = tavg_specbudgx(i,j,k) % vh_dudzh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dvdxh = tavg_specbudgx(i,j,k) % vh_dvdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dvdyh = tavg_specbudgx(i,j,k) % vh_dvdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dvdzh = tavg_specbudgx(i,j,k) % vh_dvdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dwdxh = tavg_specbudgx(i,j,k) % vh_dwdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dwdyh = tavg_specbudgx(i,j,k) % vh_dwdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_dwdzh = tavg_specbudgx(i,j,k) % vh_dwdzh / tavg_total_time

    tavg_specbudgx(i,j,k) % wh_dudxh = tavg_specbudgx(i,j,k) % wh_dudxh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dudyh = tavg_specbudgx(i,j,k) % wh_dudyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dudzh = tavg_specbudgx(i,j,k) % wh_dudzh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dvdxh = tavg_specbudgx(i,j,k) % wh_dvdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dvdyh = tavg_specbudgx(i,j,k) % wh_dvdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dvdzh = tavg_specbudgx(i,j,k) % wh_dvdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dwdxh = tavg_specbudgx(i,j,k) % wh_dwdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dwdyh = tavg_specbudgx(i,j,k) % wh_dwdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_dwdzh = tavg_specbudgx(i,j,k) % wh_dwdzh / tavg_total_time

    ! Mean velocity * velocity gradient product, uk*duidxk_hat
    tavg_specbudgx(i,j,k) % ududxh = tavg_specbudgx(i,j,k) % ududxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vdudyh = tavg_specbudgx(i,j,k) % vdudyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wdudzh = tavg_specbudgx(i,j,k) % wdudzh / tavg_total_time
    tavg_specbudgx(i,j,k) % udvdxh = tavg_specbudgx(i,j,k) % udvdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vdvdyh = tavg_specbudgx(i,j,k) % vdvdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wdvdzh = tavg_specbudgx(i,j,k) % wdvdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % udwdxh = tavg_specbudgx(i,j,k) % udwdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vdwdyh = tavg_specbudgx(i,j,k) % vdwdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wdwdzh = tavg_specbudgx(i,j,k) % wdwdzh / tavg_total_time

    ! Mean vel * vel * vel gradient product, uj_hat*(uk*duidxk)_hat
    tavg_specbudgx(i,j,k) % uh_ududxh = tavg_specbudgx(i,j,k)%uh_ududxh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_vdudyh = tavg_specbudgx(i,j,k)%uh_vdudyh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_wdudzh = tavg_specbudgx(i,j,k)%uh_wdudzh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_ududxh = tavg_specbudgx(i,j,k)%vh_ududxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_vdudyh = tavg_specbudgx(i,j,k)%vh_vdudyh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_wdudzh = tavg_specbudgx(i,j,k)%vh_wdudzh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_ududxh = tavg_specbudgx(i,j,k)%wh_ududxh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_vdudyh = tavg_specbudgx(i,j,k)%wh_vdudyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_wdudzh = tavg_specbudgx(i,j,k)%wh_wdudzh / tavg_total_time

    tavg_specbudgx(i,j,k) % uh_udvdxh = tavg_specbudgx(i,j,k)%uh_udvdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_vdvdyh = tavg_specbudgx(i,j,k)%uh_vdvdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_wdvdzh = tavg_specbudgx(i,j,k)%uh_wdvdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_udvdxh = tavg_specbudgx(i,j,k)%vh_udvdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_vdvdyh = tavg_specbudgx(i,j,k)%vh_vdvdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_wdvdzh = tavg_specbudgx(i,j,k)%vh_wdvdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_udvdxh = tavg_specbudgx(i,j,k)%wh_udvdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_vdvdyh = tavg_specbudgx(i,j,k)%wh_vdvdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_wdvdzh = tavg_specbudgx(i,j,k)%wh_wdvdzh / tavg_total_time

    tavg_specbudgx(i,j,k) % uh_udwdxh = tavg_specbudgx(i,j,k)%uh_udwdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_vdwdyh = tavg_specbudgx(i,j,k)%uh_vdwdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % uh_wdwdzh = tavg_specbudgx(i,j,k)%uh_wdwdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_udwdxh = tavg_specbudgx(i,j,k)%vh_udwdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_vdwdyh = tavg_specbudgx(i,j,k)%vh_vdwdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % vh_wdwdzh = tavg_specbudgx(i,j,k)%vh_wdwdzh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_udwdxh = tavg_specbudgx(i,j,k)%wh_udwdxh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_vdwdyh = tavg_specbudgx(i,j,k)%wh_vdwdyh / tavg_total_time
    tavg_specbudgx(i,j,k) % wh_wdwdzh = tavg_specbudgx(i,j,k)%wh_wdwdzh / tavg_total_time

    ! Mean velocity * pressure gradient product, ui*dpdxj
    tavg_specbudgx(i,j,k) % udpdx = tavg_specbudgx(i,j,k) % udpdx / tavg_total_time
    tavg_specbudgx(i,j,k) % udpdy = tavg_specbudgx(i,j,k) % udpdy / tavg_total_time
    tavg_specbudgx(i,j,k) % udpdz = tavg_specbudgx(i,j,k) % udpdz / tavg_total_time
    tavg_specbudgx(i,j,k) % vdpdx = tavg_specbudgx(i,j,k) % vdpdx / tavg_total_time
    tavg_specbudgx(i,j,k) % vdpdy = tavg_specbudgx(i,j,k) % vdpdy / tavg_total_time
    tavg_specbudgx(i,j,k) % vdpdz = tavg_specbudgx(i,j,k) % vdpdz / tavg_total_time
    tavg_specbudgx(i,j,k) % wdpdx = tavg_specbudgx(i,j,k) % wdpdx / tavg_total_time
    tavg_specbudgx(i,j,k) % wdpdy = tavg_specbudgx(i,j,k) % wdpdy / tavg_total_time
    tavg_specbudgx(i,j,k) % wdpdz = tavg_specbudgx(i,j,k) % wdpdz / tavg_total_time

    ! Mean pressure * velocity gradient product, p*duidxj
    tavg_specbudgx(i,j,k) % pdudx = tavg_specbudgx(i,j,k) % pdudx / tavg_total_time
    tavg_specbudgx(i,j,k) % pdudy = tavg_specbudgx(i,j,k) % pdudy / tavg_total_time
    tavg_specbudgx(i,j,k) % pdudz = tavg_specbudgx(i,j,k) % pdudz / tavg_total_time
    tavg_specbudgx(i,j,k) % pdvdx = tavg_specbudgx(i,j,k) % pdvdx / tavg_total_time
    tavg_specbudgx(i,j,k) % pdvdy = tavg_specbudgx(i,j,k) % pdvdy / tavg_total_time
    tavg_specbudgx(i,j,k) % pdvdz = tavg_specbudgx(i,j,k) % pdvdz / tavg_total_time
    tavg_specbudgx(i,j,k) % pdwdx = tavg_specbudgx(i,j,k) % pdwdx / tavg_total_time
    tavg_specbudgx(i,j,k) % pdwdy = tavg_specbudgx(i,j,k) % pdwdy / tavg_total_time
    tavg_specbudgx(i,j,k) % pdwdz = tavg_specbudgx(i,j,k) % pdwdz / tavg_total_time

    ! Mean Laplacian, nu*lap(uj)
    tavg_specbudgx(i,j,k) % lapu = tavg_specbudgx(i,j,k) % lapu / tavg_total_time
    tavg_specbudgx(i,j,k) % lapv = tavg_specbudgx(i,j,k) % lapv / tavg_total_time
    tavg_specbudgx(i,j,k) % lapw = tavg_specbudgx(i,j,k) % lapw / tavg_total_time

    ! Mean Vel-Laplacian, nu*ui*lap(uj)
    tavg_specbudgx(i,j,k) % ulapu = tavg_specbudgx(i,j,k) % ulapu / tavg_total_time
    tavg_specbudgx(i,j,k) % ulapv = tavg_specbudgx(i,j,k) % ulapv / tavg_total_time
    tavg_specbudgx(i,j,k) % ulapw = tavg_specbudgx(i,j,k) % ulapw / tavg_total_time
    tavg_specbudgx(i,j,k) % vlapu = tavg_specbudgx(i,j,k) % vlapu / tavg_total_time
    tavg_specbudgx(i,j,k) % vlapv = tavg_specbudgx(i,j,k) % vlapv / tavg_total_time
    tavg_specbudgx(i,j,k) % vlapw = tavg_specbudgx(i,j,k) % vlapw / tavg_total_time
    tavg_specbudgx(i,j,k) % wlapu = tavg_specbudgx(i,j,k) % wlapu / tavg_total_time
    tavg_specbudgx(i,j,k) % wlapv = tavg_specbudgx(i,j,k) % wlapv / tavg_total_time
    tavg_specbudgx(i,j,k) % wlapw = tavg_specbudgx(i,j,k) % wlapw / tavg_total_time
#endif

end do 
end do

do j = 1, Ny/2 + 1
do i = 1, Nx
    tavg_turbspecy(i,j,k) % uf = tavg_turbspecy(i,j,k) % uf / tavg_total_time
    tavg_turbspecy(i,j,k) % vf = tavg_turbspecy(i,j,k) % vf / tavg_total_time
    tavg_turbspecy(i,j,k) % wf = tavg_turbspecy(i,j,k) % wf / tavg_total_time

    tavg_turbspecy(i,j,k) % uu = tavg_turbspecy(i,j,k) % uu / tavg_total_time
    tavg_turbspecy(i,j,k) % vv = tavg_turbspecy(i,j,k) % vv / tavg_total_time
    tavg_turbspecy(i,j,k) % ww = tavg_turbspecy(i,j,k) % ww / tavg_total_time
    !tavg_turbspecy(i,j,k) % vel2 = tavg_turbspecy(i,j,k) % vel2 / tavg_total_time

    tavg_turbspecy(i,j,k) % uv = tavg_turbspecy(i,j,k) % uv / tavg_total_time
    tavg_turbspecy(i,j,k) % uw = tavg_turbspecy(i,j,k) % uw / tavg_total_time
    tavg_turbspecy(i,j,k) % vw = tavg_turbspecy(i,j,k) % vw / tavg_total_time

    tavg_turbspecy(i,j,k) % vortxf = tavg_turbspecy(i,j,k) % vortxf / tavg_total_time
    tavg_turbspecy(i,j,k) % vortyf = tavg_turbspecy(i,j,k) % vortyf / tavg_total_time
    tavg_turbspecy(i,j,k) % vortzf = tavg_turbspecy(i,j,k) % vortzf / tavg_total_time

    tavg_turbspecy(i,j,k) % vortx2 = tavg_turbspecy(i,j,k) % vortx2 / tavg_total_time
    tavg_turbspecy(i,j,k) % vorty2 = tavg_turbspecy(i,j,k) % vorty2 / tavg_total_time
    tavg_turbspecy(i,j,k) % vortz2 = tavg_turbspecy(i,j,k) % vortz2 / tavg_total_time
    !tavg_turbspecy(i,j,k) % vort2 = tavg_turbspecy(i,j,k) % vort2 / tavg_total_time

#ifdef PPSCALARS
    tavg_scal_turbspecy(i,j,k) % thetaf = tavg_scal_turbspecy(i,j,k) % thetaf / tavg_total_time

    tavg_scal_turbspecy(i,j,k) % theta2 = tavg_scal_turbspecy(i,j,k) % theta2 / tavg_total_time
    tavg_scal_turbspecy(i,j,k) % utheta = tavg_scal_turbspecy(i,j,k) % utheta / tavg_total_time
    tavg_scal_turbspecy(i,j,k) % vtheta = tavg_scal_turbspecy(i,j,k) % vtheta / tavg_total_time
    tavg_scal_turbspecy(i,j,k) % wtheta = tavg_scal_turbspecy(i,j,k) % wtheta / tavg_total_time
#endif

#ifdef PPOUTPUT_SPECBUDG
    ! Mean pressure
    tavg_specbudgy(i,j,k) % p = tavg_specbudgy(i,j,k) % p / tavg_total_time

    ! Mean velocity gradient
    tavg_specbudgy(i,j,k) % dudx = tavg_specbudgy(i,j,k) % dudx / tavg_total_time
    tavg_specbudgy(i,j,k) % dudy = tavg_specbudgy(i,j,k) % dudy / tavg_total_time
    tavg_specbudgy(i,j,k) % dudz = tavg_specbudgy(i,j,k) % dudz / tavg_total_time
    tavg_specbudgy(i,j,k) % dvdx = tavg_specbudgy(i,j,k) % dvdx / tavg_total_time
    tavg_specbudgy(i,j,k) % dvdy = tavg_specbudgy(i,j,k) % dvdy / tavg_total_time
    tavg_specbudgy(i,j,k) % dvdz = tavg_specbudgy(i,j,k) % dvdz / tavg_total_time
    tavg_specbudgy(i,j,k) % dwdx = tavg_specbudgy(i,j,k) % dwdx / tavg_total_time
    tavg_specbudgy(i,j,k) % dwdy = tavg_specbudgy(i,j,k) % dwdy / tavg_total_time
    tavg_specbudgy(i,j,k) % dwdz = tavg_specbudgy(i,j,k) % dwdz / tavg_total_time

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
    tavg_specbudgy(i,j,k) % uxux = tavg_specbudgy(i,j,k) % uxux / tavg_total_time
    tavg_specbudgy(i,j,k) % uyuy = tavg_specbudgy(i,j,k) % uyuy / tavg_total_time
    tavg_specbudgy(i,j,k) % uzuz = tavg_specbudgy(i,j,k) % uzuz / tavg_total_time
    tavg_specbudgy(i,j,k) % vxvx = tavg_specbudgy(i,j,k) % vxvx / tavg_total_time
    tavg_specbudgy(i,j,k) % vyvy = tavg_specbudgy(i,j,k) % vyvy / tavg_total_time
    tavg_specbudgy(i,j,k) % vzvz = tavg_specbudgy(i,j,k) % vzvz / tavg_total_time
    tavg_specbudgy(i,j,k) % wxwx = tavg_specbudgy(i,j,k) % wxwx / tavg_total_time
    tavg_specbudgy(i,j,k) % wywy = tavg_specbudgy(i,j,k) % wywy / tavg_total_time
    tavg_specbudgy(i,j,k) % wzwz = tavg_specbudgy(i,j,k) % wzwz / tavg_total_time

    ! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
    tavg_specbudgy(i,j,k) % uxvx = tavg_specbudgy(i,j,k) % uxvx / tavg_total_time
    tavg_specbudgy(i,j,k) % uyvy = tavg_specbudgy(i,j,k) % uyvy / tavg_total_time
    tavg_specbudgy(i,j,k) % uzvz = tavg_specbudgy(i,j,k) % uzvz / tavg_total_time
    tavg_specbudgy(i,j,k) % uxwx = tavg_specbudgy(i,j,k) % uxwx / tavg_total_time
    tavg_specbudgy(i,j,k) % uywy = tavg_specbudgy(i,j,k) % uywy / tavg_total_time
    tavg_specbudgy(i,j,k) % uzwz = tavg_specbudgy(i,j,k) % uzwz / tavg_total_time
    tavg_specbudgy(i,j,k) % vxwx = tavg_specbudgy(i,j,k) % vxwx / tavg_total_time
    tavg_specbudgy(i,j,k) % vywy = tavg_specbudgy(i,j,k) % vywy / tavg_total_time
    tavg_specbudgy(i,j,k) % vzwz = tavg_specbudgy(i,j,k) % vzwz / tavg_total_time

    ! Mean pressure gradients, dpdxi
    tavg_specbudgy(i,j,k) % dpdx = tavg_specbudgy(i,j,k) % dpdx / tavg_total_time
    tavg_specbudgy(i,j,k) % dpdy = tavg_specbudgy(i,j,k) % dpdy / tavg_total_time
    tavg_specbudgy(i,j,k) % dpdz = tavg_specbudgy(i,j,k) % dpdz / tavg_total_time

    ! Mean velocity * velocity gradient product, ui_hat*dujdxk_hat
    tavg_specbudgy(i,j,k) % uh_dudxh = tavg_specbudgy(i,j,k) % uh_dudxh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dudyh = tavg_specbudgy(i,j,k) % uh_dudyh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dudzh = tavg_specbudgy(i,j,k) % uh_dudzh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dvdxh = tavg_specbudgy(i,j,k) % uh_dvdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dvdyh = tavg_specbudgy(i,j,k) % uh_dvdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dvdzh = tavg_specbudgy(i,j,k) % uh_dvdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dwdxh = tavg_specbudgy(i,j,k) % uh_dwdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dwdyh = tavg_specbudgy(i,j,k) % uh_dwdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_dwdzh = tavg_specbudgy(i,j,k) % uh_dwdzh / tavg_total_time

    tavg_specbudgy(i,j,k) % vh_dudxh = tavg_specbudgy(i,j,k) % vh_dudxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dudyh = tavg_specbudgy(i,j,k) % vh_dudyh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dudzh = tavg_specbudgy(i,j,k) % vh_dudzh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dvdxh = tavg_specbudgy(i,j,k) % vh_dvdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dvdyh = tavg_specbudgy(i,j,k) % vh_dvdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dvdzh = tavg_specbudgy(i,j,k) % vh_dvdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dwdxh = tavg_specbudgy(i,j,k) % vh_dwdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dwdyh = tavg_specbudgy(i,j,k) % vh_dwdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_dwdzh = tavg_specbudgy(i,j,k) % vh_dwdzh / tavg_total_time

    tavg_specbudgy(i,j,k) % wh_dudxh = tavg_specbudgy(i,j,k) % wh_dudxh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dudyh = tavg_specbudgy(i,j,k) % wh_dudyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dudzh = tavg_specbudgy(i,j,k) % wh_dudzh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dvdxh = tavg_specbudgy(i,j,k) % wh_dvdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dvdyh = tavg_specbudgy(i,j,k) % wh_dvdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dvdzh = tavg_specbudgy(i,j,k) % wh_dvdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dwdxh = tavg_specbudgy(i,j,k) % wh_dwdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dwdyh = tavg_specbudgy(i,j,k) % wh_dwdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_dwdzh = tavg_specbudgy(i,j,k) % wh_dwdzh / tavg_total_time

    ! Mean velocity * velocity gradient product, uk*duidxk_hat
    tavg_specbudgy(i,j,k) % ududxh = tavg_specbudgy(i,j,k) % ududxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vdudyh = tavg_specbudgy(i,j,k) % vdudyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wdudzh = tavg_specbudgy(i,j,k) % wdudzh / tavg_total_time
    tavg_specbudgy(i,j,k) % udvdxh = tavg_specbudgy(i,j,k) % udvdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vdvdyh = tavg_specbudgy(i,j,k) % vdvdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wdvdzh = tavg_specbudgy(i,j,k) % wdvdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % udwdxh = tavg_specbudgy(i,j,k) % udwdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vdwdyh = tavg_specbudgy(i,j,k) % vdwdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wdwdzh = tavg_specbudgy(i,j,k) % wdwdzh / tavg_total_time

    ! Mean vel * vel * vel gradient product, uj_hat*(uk*duidxk)_hat
    tavg_specbudgy(i,j,k) % uh_ududxh = tavg_specbudgy(i,j,k)%uh_ududxh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_vdudyh = tavg_specbudgy(i,j,k)%uh_vdudyh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_wdudzh = tavg_specbudgy(i,j,k)%uh_wdudzh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_ududxh = tavg_specbudgy(i,j,k)%vh_ududxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_vdudyh = tavg_specbudgy(i,j,k)%vh_vdudyh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_wdudzh = tavg_specbudgy(i,j,k)%vh_wdudzh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_ududxh = tavg_specbudgy(i,j,k)%wh_ududxh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_vdudyh = tavg_specbudgy(i,j,k)%wh_vdudyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_wdudzh = tavg_specbudgy(i,j,k)%wh_wdudzh / tavg_total_time

    tavg_specbudgy(i,j,k) % uh_udvdxh = tavg_specbudgy(i,j,k)%uh_udvdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_vdvdyh = tavg_specbudgy(i,j,k)%uh_vdvdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_wdvdzh = tavg_specbudgy(i,j,k)%uh_wdvdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_udvdxh = tavg_specbudgy(i,j,k)%vh_udvdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_vdvdyh = tavg_specbudgy(i,j,k)%vh_vdvdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_wdvdzh = tavg_specbudgy(i,j,k)%vh_wdvdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_udvdxh = tavg_specbudgy(i,j,k)%wh_udvdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_vdvdyh = tavg_specbudgy(i,j,k)%wh_vdvdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_wdvdzh = tavg_specbudgy(i,j,k)%wh_wdvdzh / tavg_total_time

    tavg_specbudgy(i,j,k) % uh_udwdxh = tavg_specbudgy(i,j,k)%uh_udwdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_vdwdyh = tavg_specbudgy(i,j,k)%uh_vdwdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % uh_wdwdzh = tavg_specbudgy(i,j,k)%uh_wdwdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_udwdxh = tavg_specbudgy(i,j,k)%vh_udwdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_vdwdyh = tavg_specbudgy(i,j,k)%vh_vdwdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % vh_wdwdzh = tavg_specbudgy(i,j,k)%vh_wdwdzh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_udwdxh = tavg_specbudgy(i,j,k)%wh_udwdxh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_vdwdyh = tavg_specbudgy(i,j,k)%wh_vdwdyh / tavg_total_time
    tavg_specbudgy(i,j,k) % wh_wdwdzh = tavg_specbudgy(i,j,k)%wh_wdwdzh / tavg_total_time

    ! Mean velocity * pressure gradient product, ui*dpdxj
    tavg_specbudgy(i,j,k) % udpdx = tavg_specbudgy(i,j,k) % udpdx / tavg_total_time
    tavg_specbudgy(i,j,k) % udpdy = tavg_specbudgy(i,j,k) % udpdy / tavg_total_time
    tavg_specbudgy(i,j,k) % udpdz = tavg_specbudgy(i,j,k) % udpdz / tavg_total_time
    tavg_specbudgy(i,j,k) % vdpdx = tavg_specbudgy(i,j,k) % vdpdx / tavg_total_time
    tavg_specbudgy(i,j,k) % vdpdy = tavg_specbudgy(i,j,k) % vdpdy / tavg_total_time
    tavg_specbudgy(i,j,k) % vdpdz = tavg_specbudgy(i,j,k) % vdpdz / tavg_total_time
    tavg_specbudgy(i,j,k) % wdpdx = tavg_specbudgy(i,j,k) % wdpdx / tavg_total_time
    tavg_specbudgy(i,j,k) % wdpdy = tavg_specbudgy(i,j,k) % wdpdy / tavg_total_time
    tavg_specbudgy(i,j,k) % wdpdz = tavg_specbudgy(i,j,k) % wdpdz / tavg_total_time

    ! Mean pressure * velocity gradient product, p*duidxj
    tavg_specbudgy(i,j,k) % pdudx = tavg_specbudgy(i,j,k) % pdudx / tavg_total_time
    tavg_specbudgy(i,j,k) % pdudy = tavg_specbudgy(i,j,k) % pdudy / tavg_total_time
    tavg_specbudgy(i,j,k) % pdudz = tavg_specbudgy(i,j,k) % pdudz / tavg_total_time
    tavg_specbudgy(i,j,k) % pdvdx = tavg_specbudgy(i,j,k) % pdvdx / tavg_total_time
    tavg_specbudgy(i,j,k) % pdvdy = tavg_specbudgy(i,j,k) % pdvdy / tavg_total_time
    tavg_specbudgy(i,j,k) % pdvdz = tavg_specbudgy(i,j,k) % pdvdz / tavg_total_time
    tavg_specbudgy(i,j,k) % pdwdx = tavg_specbudgy(i,j,k) % pdwdx / tavg_total_time
    tavg_specbudgy(i,j,k) % pdwdy = tavg_specbudgy(i,j,k) % pdwdy / tavg_total_time
    tavg_specbudgy(i,j,k) % pdwdz = tavg_specbudgy(i,j,k) % pdwdz / tavg_total_time

    ! Mean Laplacian, nu*lap(uj)
    tavg_specbudgy(i,j,k) % lapu = tavg_specbudgy(i,j,k) % lapu / tavg_total_time
    tavg_specbudgy(i,j,k) % lapv = tavg_specbudgy(i,j,k) % lapv / tavg_total_time
    tavg_specbudgy(i,j,k) % lapw = tavg_specbudgy(i,j,k) % lapw / tavg_total_time

    ! Mean Vel-Laplacian, nu*ui*lap(uj)
    tavg_specbudgy(i,j,k) % ulapu = tavg_specbudgy(i,j,k) % ulapu / tavg_total_time
    tavg_specbudgy(i,j,k) % ulapv = tavg_specbudgy(i,j,k) % ulapv / tavg_total_time
    tavg_specbudgy(i,j,k) % ulapw = tavg_specbudgy(i,j,k) % ulapw / tavg_total_time
    tavg_specbudgy(i,j,k) % vlapu = tavg_specbudgy(i,j,k) % vlapu / tavg_total_time
    tavg_specbudgy(i,j,k) % vlapv = tavg_specbudgy(i,j,k) % vlapv / tavg_total_time
    tavg_specbudgy(i,j,k) % vlapw = tavg_specbudgy(i,j,k) % vlapw / tavg_total_time
    tavg_specbudgy(i,j,k) % wlapu = tavg_specbudgy(i,j,k) % wlapu / tavg_total_time
    tavg_specbudgy(i,j,k) % wlapv = tavg_specbudgy(i,j,k) % wlapv / tavg_total_time
    tavg_specbudgy(i,j,k) % wlapw = tavg_specbudgy(i,j,k) % wlapw / tavg_total_time
#endif

end do
end do

end do
#endif

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

!  Sync entire tavg structure
#ifdef PPMPI
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%u, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%v, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%w_uv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%u_w, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%v_w, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%w, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%u2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%v2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%w2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%uw, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%vw, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%uv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%p, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%fx, 0, MPI_SYNC_DOWNUP )
#ifdef PPSCALARS
call mpi_sync_real_array( tavg_scal(1:nx,1:ny,lbz:nz)%theta, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_scal(1:nx,1:ny,lbz:nz)%theta2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_scal(1:nx,1:ny,lbz:nz)%utheta, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_scal(1:nx,1:ny,lbz:nz)%vtheta, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_scal(1:nx,1:ny,lbz:nz)%wtheta, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_scal(1:nx,1:ny,lbz:nz)%pi_z, 0, MPI_SYNC_DOWNUP )
#endif
#ifdef PPOUTPUT_SGS
call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%cs_opt2, 0, MPI_SYNC_DOWNUP )
! call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%Tn, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%Nu_t, 0, MPI_SYNC_DOWNUP )
! call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%F_LM, 0, MPI_SYNC_DOWNUP )
! call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%F_MM, 0, MPI_SYNC_DOWNUP )
! call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%F_QN, 0, MPI_SYNC_DOWNUP )
! call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%F_NN, 0, MPI_SYNC_DOWNUP )
! call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%ee_now, 0, MPI_SYNC_DOWNUP )
! #ifdef PPDYN_TN
! call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%F_ee2, 0, MPI_SYNC_DOWNUP )
! call mpi_sync_real_array(tavg_sgs(1:nx,1:ny,lbz:nz)%F_deedt2, 0, MPI_SYNC_DOWNUP)
! #endif
#endif
#ifdef PPOUTPUT_BUDGET
! Mean pressure on w-grid
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%p, 0, MPI_SYNC_DOWNUP )

! Mean velocity-velocity product, ui*uj, on w-grid
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uu, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%ww, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uw, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vw, 0, MPI_SYNC_DOWNUP )

! Mean velocity gradients, duidxj
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dudx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dudy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dudz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dwdz, 0, MPI_SYNC_DOWNUP )

! Mean pressure gradients, dpdxi
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dpdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dpdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%dpdz, 0, MPI_SYNC_DOWNUP )

! Mean vel-velGrad product, ui*dujdxk
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%ududx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%ududy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%ududz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udwdz, 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdudx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdudy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdudz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdwdz, 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdudx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdudy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdudz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdwdz, 0, MPI_SYNC_DOWNUP )

! Mean vel-vel-velGrad product, ui*uk*dujdxk
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uududx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uvdudy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uwdudz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uudvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uvdvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uwdvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uudwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uvdwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uwdwdz, 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vududx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vvdudy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vwdudz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vudvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vvdvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vwdvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vudwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vvdwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vwdwdz, 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wududx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wvdudy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wwdudz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wudvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wvdvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wwdvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wudwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wvdwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%wwdwdz, 0, MPI_SYNC_DOWNUP )

! Mean velGrad-velGrad product, duidxk*dujdxk, i=j
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uxux, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uyuy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uzuz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vxvx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vyvy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vzvz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wxwx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wywy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wzwz, 0, MPI_SYNC_DOWNUP )

! Mean velGrad-velGrad product, duidxk*dujdxk, i/=j
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uxvx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uyvy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uzvz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uxwx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uywy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%uzwz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vxwx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vywy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vzwz, 0, MPI_SYNC_DOWNUP )

! duidxj*dujdxi, i /= j
! call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uyvx, 0, MPI_SYNC_DOWNUP)
! call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%uzwx, 0, MPI_SYNC_DOWNUP)
! call mpi_sync_real_array(tavg_budget(1:nx,1:ny,lbz:nz)%vzwy, 0, MPI_SYNC_DOWNUP)

! Mean vel-presGrad product, ui*dpdxj
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udpdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udpdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%udpdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdpdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdpdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vdpdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdpdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdpdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wdpdz, 0, MPI_SYNC_DOWNUP )

! Mean pres-velGrad product, p*duidxj
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdudx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdudy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdudz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdvdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdvdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdvdz, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdwdx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdwdy, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%pdwdz, 0, MPI_SYNC_DOWNUP )

! Mean Laplacian, nu*lap(uj)
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%lapu, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%lapv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%lapw, 0, MPI_SYNC_DOWNUP )

! Mean Vel-Laplacian, nu*ui*lap(uj)
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%ulapu, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%ulapv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%ulapw, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vlapu, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vlapv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%vlapw, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wlapu, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wlapv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg_budget(1:nx,1:ny,lbz:nz)%wlapw, 0, MPI_SYNC_DOWNUP )

#endif

#endif

! Write all the 3D data
#ifdef PPCGNS
! Write CGNS Data
call write_parallel_cgns (fname_vel ,nx, ny, nz - nz_end, nz_tot,           &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 3,                               &
    (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                            &
    (/ tavg(1:nx,1:ny,1:nz - nz_end) % u,                                   &
       tavg(1:nx,1:ny,1:nz - nz_end) % v,                                   &
       tavg(1:nx,1:ny,1:nz- nz_end) % w_uv /) )

call write_parallel_cgns (fname_velw ,nx, ny, nz - nz_end, nz_tot,          &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                                 &
    1, (/ 'VelocityZ' /), (/ tavg(1:nx,1:ny,1:nz- nz_end) % w /) )

!call write_parallel_cgns(fname_vel2,nx,ny,nz- nz_end,nz_tot,               &
!    (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
!    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
!    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                             &
!    (/ 'Mean--uu', 'Mean--vv', 'Mean--ww','Mean--uw','Mean--vw','Mean--uv'/),  &
!    (/ tavg(1:nx,1:ny,1:nz- nz_end) % u2,                                  &
!       tavg(1:nx,1:ny,1:nz- nz_end) % v2,                                  &
!       tavg(1:nx,1:ny,1:nz- nz_end) % w2,                                  &
!       tavg(1:nx,1:ny,1:nz- nz_end) % uw,                                  &
!       tavg(1:nx,1:ny,1:nz- nz_end) % vw,                                  &
!       tavg(1:nx,1:ny,1:nz- nz_end) % uv /) )

call write_parallel_cgns(fname_tau,nx,ny,nz- nz_end,nz_tot,                 &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                              &
    (/ 'Tau--txx', 'Tau--txy', 'Tau--tyy','Tau--txz','Tau--tyz','Tau--tzz'/),&
    (/ tavg(1:nx,1:ny,1:nz- nz_end) % txx,                                  &
       tavg(1:nx,1:ny,1:nz- nz_end) % txy,                                  &
       tavg(1:nx,1:ny,1:nz- nz_end) % tyy,                                  &
       tavg(1:nx,1:ny,1:nz- nz_end) % txz,                                  &
       tavg(1:nx,1:ny,1:nz- nz_end) % tyz,                                  &
       tavg(1:nx,1:ny,1:nz- nz_end) % tzz /) )

call write_parallel_cgns(fname_pres,nx,ny,nz- nz_end,nz_tot,                &
   (/ 1, 1,   (nz-1)*coord + 1 /),                                          &
   (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                             &
   x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                               &
   (/ 'pressure' /),                                                        &
   (/ tavg(1:nx,1:ny,1:nz- nz_end) % p /) )

!call write_parallel_cgns(fname_f,nx,ny,nz- nz_end,nz_tot,                  &
!    (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
!    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
!    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 3,                             &
!    (/ 'bodyForX', 'bodyForY', 'bodyForZ' /),                              &
!    (/ tavg(1:nx,1:ny,1:nz- nz_end) % fx,                                  &
!       tavg(1:nx,1:ny,1:nz- nz_end) % fy,                                  &
!       tavg(1:nx,1:ny,1:nz- nz_end) % fz /) )

call write_parallel_cgns (fname_vort ,nx, ny, nz - nz_end, nz_tot,          &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 3,                               &
    (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                         &
    (/ tavg(1:nx,1:ny,1:nz - nz_end) % vortx,                               &
       tavg(1:nx,1:ny,1:nz - nz_end) % vorty,                               &
       tavg(1:nx,1:ny,1:nz- nz_end) % vortz /) )

!call write_parallel_cgns(fname_velgrad,nx,ny,nz- nz_end,nz_tot,            &
!    (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
!    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
!    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 9,                              &
!    (/ 'dudx', 'dudy', 'dudz','dvdx','dvdy','dvdz','dwdx','dwdy','dwdz'/), &
!    (/ tavg(1:nx,1:ny,1:nz- nz_end) % dudx,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dudy,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dudz,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dvdx,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dvdy,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dvdz,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dwdx,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dwdy,                                &
!       tavg(1:nx,1:ny,1:nz- nz_end) % dwdz /) )

#else
! Write binary data
open(unit=13, file=fname_vel, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%u
write(13,rec=2) tavg(:nx,:ny,1:nz)%v
write(13,rec=3) tavg(:nx,:ny,1:nz)%w_uv
close(13)

! Write binary data
open(unit=13, file=fname_velw, form='unformatted', convert=write_endian,    &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%w
close(13)

!open(unit=13, file=fname_vel2, form='unformatted', convert=write_endian,   &
!    access='direct', recl=nx*ny*nz*rprec)
!write(13,rec=1) tavg(:nx,:ny,1:nz)%u2
!write(13,rec=2) tavg(:nx,:ny,1:nz)%v2
!write(13,rec=3) tavg(:nx,:ny,1:nz)%w2
!write(13,rec=4) tavg(:nx,:ny,1:nz)%uw
!write(13,rec=5) tavg(:nx,:ny,1:nz)%vw
!write(13,rec=6) tavg(:nx,:ny,1:nz)%uv
!close(13)

open(unit=13, file=fname_tau, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%txx
write(13,rec=2) tavg(:nx,:ny,1:nz)%txy
write(13,rec=3) tavg(:nx,:ny,1:nz)%tyy
write(13,rec=4) tavg(:nx,:ny,1:nz)%txz
write(13,rec=5) tavg(:nx,:ny,1:nz)%tyz
write(13,rec=6) tavg(:nx,:ny,1:nz)%tzz
close(13)

open(unit=13, file=fname_pres, form='unformatted', convert=write_endian,    &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%p
close(13)

!open(unit=13, file=fname_f, form='unformatted', convert=write_endian,       &
!    access='direct', recl=nx*ny*nz*rprec)
!write(13,rec=1) tavg(:nx,:ny,1:nz)%fx
!write(13,rec=2) tavg(:nx,:ny,1:nz)%fy
!write(13,rec=3) tavg(:nx,:ny,1:nz)%fz
!close(13)


open(unit=13, file=fname_vort, form='unformatted', convert=write_endian,    &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) tavg_vort(:nx,:ny,1:nz)%vortx
write(13,rec=2) tavg_vort(:nx,:ny,1:nz)%vorty
write(13,rec=3) tavg_vort(:nx,:ny,1:nz)%vortz
close(13)

!open(unit=13, file=fname_velgrad, form='unformatted', convert=write_endian, &
!    access='direct', recl=nx*ny*nz*rprec)
!write(13,rec=1) tavg(:nx,:ny,1:nz)%dudx
!write(13,rec=2) tavg(:nx,:ny,1:nz)%dudy
!write(13,rec=3) tavg(:nx,:ny,1:nz)%dudz
!write(13,rec=4) tavg(:nx,:ny,1:nz)%dvdx
!write(13,rec=5) tavg(:nx,:ny,1:nz)%dvdy
!write(13,rec=6) tavg(:nx,:ny,1:nz)%dvdz
!write(13,rec=7) tavg(:nx,:ny,1:nz)%dwdx
!write(13,rec=8) tavg(:nx,:ny,1:nz)%dwdy
!write(13,rec=9) tavg(:nx,:ny,1:nz)%dwdz
!close(13)

#endif

#ifdef PPMPI
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
#endif

! Do the Reynolds stress calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate(rs(nx,ny,lbz:nz))
rs = rs_compute(tavg , lbz)

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rs,nx,ny,nz- nz_end,nz_tot,                  &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                               &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),&
    (/ rs(1:nx,1:ny,1:nz- nz_end) % up2,                                    &
    rs(1:nx,1:ny,1:nz- nz_end) % vp2,                                       &
    rs(1:nx,1:ny,1:nz- nz_end) % wp2,                                       &
    rs(1:nx,1:ny,1:nz- nz_end) % upwp,                                      &
    rs(1:nx,1:ny,1:nz- nz_end) % vpwp,                                      &
    rs(1:nx,1:ny,1:nz- nz_end) % upvp /) )
#else
! Write binary data
open(unit=13, file=fname_rs, form='unformatted', convert=write_endian,      &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) rs(:nx,:ny,1:nz)%up2
write(13,rec=2) rs(:nx,:ny,1:nz)%vp2
write(13,rec=3) rs(:nx,:ny,1:nz)%wp2
write(13,rec=4) rs(:nx,:ny,1:nz)%upwp
write(13,rec=5) rs(:nx,:ny,1:nz)%vpwp
write(13,rec=6) rs(:nx,:ny,1:nz)%upvp
close(13)
#endif

deallocate(rs)

#ifdef PPSCALARS
allocate(rs_scal(nx,ny,lbz:nz))
rs_scal = rs_scal_compute(tavg_scal, tavg, lbz)

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_scal,nx,ny,nz- nz_end,nz_tot,              &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                       &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                          &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                             &
    (/ 'Meantheta', 'Meanthetap', 'Meanupthetap','Meanvpthetap','Meanwpthetap'/),  &
    (/ tavg_scal(1:nx,1:ny,1:nz- nz_end) % theta,                         &
    rs_scal(1:nx,1:ny,1:nz- nz_end) % thetap2,                            &
    rs_scal(1:nx,1:ny,1:nz- nz_end) % upthetap,                           &
    rs_scal(1:nx,1:ny,1:nz- nz_end) % vpthetap,                           &
    rs_scal(1:nx,1:ny,1:nz- nz_end) % wpthetap,                           &
    tavg_scal(1:nx,1:ny,1:nz- nz_end) % pi_z /) )
#else
! Write binary data
open(unit=13, file=fname_scal, form='unformatted', convert=write_endian,           &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tavg_scal(:nx,:ny,1:nz)%theta
write(13,rec=2) rs_scal(:nx,:ny,1:nz)%thetap2
write(13,rec=3) rs_scal(:nx,:ny,1:nz)%upthetap
write(13,rec=4) rs_scal(:nx,:ny,1:nz)%vpthetap
write(13,rec=5) rs_scal(:nx,:ny,1:nz)%wpthetap
write(13,rec=6) tavg_scal(:nx,:ny,1:nz)%pi_z
close(13)
#endif

deallocate(rs_scal)
#endif

! Do the same for vorticity rms
allocate(vortrms(nx,ny,lbz:nz))
vortrms = vortrms_compute(tavg_vort, lbz)

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_vortrms,nx,ny,nz- nz_end,nz_tot,                &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 3,                                  &
    (/ 'Vortxrms', 'Vortyrms', 'Vortzrms'/),                                   &
    (/ vortrms(1:nx,1:ny,1:nz- nz_end) % vortxrms,                             &
    vortrms(1:nx,1:ny,1:nz- nz_end) % vortyrms,                                &
    vortrms(1:nx,1:ny,1:nz- nz_end) % vortzrms /) )
#else
! Write binary data
open(unit=13, file=fname_vortrms, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) vortrms(:nx,:ny,1:nz)%vortxrms
write(13,rec=2) vortrms(:nx,:ny,1:nz)%vortyrms
write(13,rec=3) vortrms(:nx,:ny,1:nz)%vortzrms
close(13)
#endif

deallocate(vortrms)

#ifdef PPOUTPUT_SGS
#ifdef PPCGNS
call write_parallel_cgns(fname_cs,nx,ny,nz- nz_end,nz_tot,                &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                       &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                          &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                            &
    (/ 'Cs_Coeff'/),  (/ tavg(1:nx,1:ny,1:nz- nz_end) % cs_opt2 /) )

call write_parallel_cgns(fname_sgs,nx,ny,nz- nz_end,nz_tot,               &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                       &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                          &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                            &
    (/ 'SGS_eddy '/),  (/ tavg_sgs(1:nx,1:ny,1:nz- nz_end) % Nu_t /) )
#else
open(unit=13, file=fname_cs, form='unformatted', convert=write_endian,    &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) tavg_sgs(:nx,:ny,1:nz)%cs_opt2
close(13)

open(unit=13, file=fname_sgs, form='unformatted', convert=write_endian,   &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) tavg_sgs(:nx,:ny,1:nz)%Nu_t
close(13)
#endif
#endif

#ifdef PPOUTPUT_BUDGET
! Do the budget calculations after time-averaging.
allocate(budget(nx,ny,lbz:nz))
budget = budget_compute(tavg_budget, tavg , lbz)

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rxx,nx,ny,nz- nz_end,nz_tot,                 &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 7,                               &
    (/ 'advxx', 'tflucxx', 'tpresxx', 'pstrainxx', 'tviscxx', 'prodxx', 'pdissxx'/),     &
    (/ budget(1:nx,1:ny,1:nz- nz_end) % advxx,                              &
    budget(1:nx,1:ny,1:nz- nz_end) % tflucxx,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % tpresxx,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % pstrainxx,                             &
    budget(1:nx,1:ny,1:nz- nz_end) % tviscxx,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % prodxx,                                &
    budget(1:nx,1:ny,1:nz- nz_end) % pdissxx /) )

call write_parallel_cgns(fname_ryy,nx,ny,nz- nz_end,nz_tot,                 &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 7,                               &
    (/ 'advyy', 'tflucyy', 'tpresyy', 'pstrainyy', 'tviscyy', 'prodyy', 'pdissyy'/),     &
    (/ budget(1:nx,1:ny,1:nz- nz_end) % advyy,                              &
    budget(1:nx,1:ny,1:nz- nz_end) % tflucyy,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % tpresyy,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % pstrainyy,                             &
    budget(1:nx,1:ny,1:nz- nz_end) % tviscyy,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % prodyy,                                &
    budget(1:nx,1:ny,1:nz- nz_end) % pdissyy /) )

call write_parallel_cgns(fname_rzz,nx,ny,nz- nz_end,nz_tot,                 &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 7,                               &
    (/ 'advzz', 'tfluczz', 'tpreszz', 'pstrainzz', 'tvisczz', 'prodzz', 'pdisszz'/),     &
    (/ budget(1:nx,1:ny,1:nz- nz_end) % advzz,                              &
    budget(1:nx,1:ny,1:nz- nz_end) % tfluczz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % tpreszz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % pstrainzz,                             &
    budget(1:nx,1:ny,1:nz- nz_end) % tvisczz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % prodzz,                                &
    budget(1:nx,1:ny,1:nz- nz_end) % pdisszz /) )

call write_parallel_cgns(fname_rxy,nx,ny,nz- nz_end,nz_tot,                 &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 7,                               &
    (/ 'advxy', 'tflucxy', 'tpresxy', 'pstrainxy', 'tviscxy', 'prodxy', 'pdissxy'/),     &
    (/ budget(1:nx,1:ny,1:nz- nz_end) % advxy,                              &
    budget(1:nx,1:ny,1:nz- nz_end) % tflucxy,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % tpresxy,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % pstrainxy,                             &
    budget(1:nx,1:ny,1:nz- nz_end) % tviscxy,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % prodxy,                                &
    budget(1:nx,1:ny,1:nz- nz_end) % pdissxy /) )

call write_parallel_cgns(fname_rxz,nx,ny,nz- nz_end,nz_tot,                 &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 7,                               &
    (/ 'advxz', 'tflucxz', 'tpresxz', 'pstrainxz', 'tviscxz', 'prodxz', 'pdissxz'/),     &
    (/ budget(1:nx,1:ny,1:nz- nz_end) % advxz,                              &
    budget(1:nx,1:ny,1:nz- nz_end) % tflucxz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % tpresxz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % pstrainxz,                             &
    budget(1:nx,1:ny,1:nz- nz_end) % tviscxz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % prodxz,                                &
    budget(1:nx,1:ny,1:nz- nz_end) % pdissxz /) )

call write_parallel_cgns(fname_ryz,nx,ny,nz- nz_end,nz_tot,                 &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                         &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                            &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 7,                               &
    (/ 'advyz', 'tflucyz', 'tpresyz', 'pstrainyz', 'tviscyz', 'prodyz', 'pdissyz'/),     &
    (/ budget(1:nx,1:ny,1:nz- nz_end) % advyz,                              &
    budget(1:nx,1:ny,1:nz- nz_end) % tflucyz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % tpresyz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % pstrainyz,                             &
    budget(1:nx,1:ny,1:nz- nz_end) % tviscyz,                               &
    budget(1:nx,1:ny,1:nz- nz_end) % prodyz,                                &
    budget(1:nx,1:ny,1:nz- nz_end) % pdissyz /) )

!call write_parallel_cgns(fname_mke,nx,ny,nz- nz_end,nz_tot,               &
!    (/ 1, 1,   (nz-1)*coord + 1 /),                                       &
!    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                          &
!    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                             &
!    (/ 'madv', 'mtfluc', 'mtpres', 'mtvisc', 'mpdiss', 'mdiss'/),         &
!    (/ budget(1:nx,1:ny,1:nz- nz_end) % madv,                             &
!    budget(1:nx,1:ny,1:nz- nz_end) % mtfluc,                              &
!    budget(1:nx,1:ny,1:nz- nz_end) % mtpres,                              &
!    budget(1:nx,1:ny,1:nz- nz_end) % mtvisc,                              &
!    budget(1:nx,1:ny,1:nz- nz_end) % mpdiss,                              &
!    budget(1:nx,1:ny,1:nz- nz_end) % mdiss /) )

#else
! Write binary data
open(unit=13, file=fname_rxx, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) budget(:nx,:ny,1:nz)%advxx
write(13,rec=2) budget(:nx,:ny,1:nz)%tflucxx
write(13,rec=3) budget(:nx,:ny,1:nz)%tpresxx
write(13,rec=4) budget(:nx,:ny,1:nz)%pstrainxx
write(13,rec=5) budget(:nx,:ny,1:nz)%tviscxx
write(13,rec=6) budget(:nx,:ny,1:nz)%prodxx
write(13,rec=7) budget(:nx,:ny,1:nz)%pdissxx
close(13)

open(unit=13, file=fname_ryy, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) budget(:nx,:ny,1:nz)%advyy
write(13,rec=2) budget(:nx,:ny,1:nz)%tflucyy
write(13,rec=3) budget(:nx,:ny,1:nz)%tpresyy
write(13,rec=4) budget(:nx,:ny,1:nz)%pstrainyy
write(13,rec=5) budget(:nx,:ny,1:nz)%tviscyy
write(13,rec=6) budget(:nx,:ny,1:nz)%prodyy
write(13,rec=7) budget(:nx,:ny,1:nz)%pdissyy
close(13)

open(unit=13, file=fname_rzz, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) budget(:nx,:ny,1:nz)%advzz
write(13,rec=2) budget(:nx,:ny,1:nz)%tfluczz
write(13,rec=3) budget(:nx,:ny,1:nz)%tpreszz
write(13,rec=4) budget(:nx,:ny,1:nz)%pstrainzz
write(13,rec=5) budget(:nx,:ny,1:nz)%tvisczz
write(13,rec=6) budget(:nx,:ny,1:nz)%prodzz
write(13,rec=7) budget(:nx,:ny,1:nz)%pdisszz
close(13)

open(unit=13, file=fname_rxy, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) budget(:nx,:ny,1:nz)%advxy
write(13,rec=2) budget(:nx,:ny,1:nz)%tflucxy
write(13,rec=3) budget(:nx,:ny,1:nz)%tpresxy
write(13,rec=4) budget(:nx,:ny,1:nz)%pstrainxy
write(13,rec=5) budget(:nx,:ny,1:nz)%tviscxy
write(13,rec=6) budget(:nx,:ny,1:nz)%prodxy
write(13,rec=7) budget(:nx,:ny,1:nz)%pdissxy
close(13)

open(unit=13, file=fname_rxz, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) budget(:nx,:ny,1:nz)%advxz
write(13,rec=2) budget(:nx,:ny,1:nz)%tflucxz
write(13,rec=3) budget(:nx,:ny,1:nz)%tpresxz
write(13,rec=4) budget(:nx,:ny,1:nz)%pstrainxz
write(13,rec=5) budget(:nx,:ny,1:nz)%tviscxz
write(13,rec=6) budget(:nx,:ny,1:nz)%prodxz
write(13,rec=7) budget(:nx,:ny,1:nz)%pdissxz
close(13)

open(unit=13, file=fname_ryz, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) budget(:nx,:ny,1:nz)%advyz
write(13,rec=2) budget(:nx,:ny,1:nz)%tflucyz
write(13,rec=3) budget(:nx,:ny,1:nz)%tpresyz
write(13,rec=4) budget(:nx,:ny,1:nz)%pstrainyz
write(13,rec=5) budget(:nx,:ny,1:nz)%tviscyz
write(13,rec=6) budget(:nx,:ny,1:nz)%prodyz
write(13,rec=7) budget(:nx,:ny,1:nz)%pdissyz
close(13)

!open(unit=13, file=fname_mke, form='unformatted', convert=write_endian,    &
!    access='direct', recl=nx*ny*nz*rprec)
!write(13,rec=1) budget(:nx,:ny,1:nz)%madv
!write(13,rec=2) budget(:nx,:ny,1:nz)%mtfluc
!write(13,rec=3) budget(:nx,:ny,1:nz)%mtpres
!write(13,rec=4) budget(:nx,:ny,1:nz)%mtvisc
!write(13,rec=5) budget(:nx,:ny,1:nz)%mpdiss
!write(13,rec=6) budget(:nx,:ny,1:nz)%mdiss
!close(13)
#endif

deallocate(budget)

#endif

#ifdef PPOUTPUT_TURBSPEC
allocate(turbspecx(nx/2+1,ny,lbz:nz))
allocate(turbspecy(nx,ny/2+1,lbz:nz))
turbspecx = turbspec_compute(tavg_turbspecx, lbz)
turbspecy = turbspec_compute(tavg_turbspecy, lbz)

#ifdef PPSCALARS
allocate(scal_turbspecx(nx/2+1,ny,lbz:nz))
allocate(scal_turbspecy(nx,ny/2+1,lbz:nz))
scal_turbspecx = scal_turbspec_compute(tavg_scal_turbspecx, tavg_turbspecx, lbz)
scal_turbspecy = scal_turbspec_compute(tavg_scal_turbspecy, tavg_turbspecy, lbz)
#endif

#ifdef PPCGNS
! Write CGNS data

! Position variable needs to be corrected!
call write_parallel_cgns(fname_sxvel,nx/2+1,ny,nz- nz_end,nz_tot,        &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                      &
    (/ nx/2+1, ny, (nz-1)*(coord+1) + 1 - nz_end /),                     &
    x(1:nx/2+1) , y(1:ny) , z(1:(nz-nz_end) ), 6,                        &
    (/ 'sxuu', 'sxvv', 'sxww', 'sxuv', 'sxuw', 'sxvw'/),                 &
    (/ turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % upup,                     &
    turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % vpvp,                        &
    turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % wpwp,                        &
    turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % upvp,                        &
    turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % upwp,                        &
    turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % vpwp /) )

call write_parallel_cgns(fname_sxvort,nx/2+1,ny,nz- nz_end,nz_tot,       &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                      &
    (/ nx/2+1, ny, (nz-1)*(coord+1) + 1 - nz_end /),                     &
    x(1:nx/2+1) , y(1:ny) , z(1:(nz-nz_end) ), 3,                        &
    (/ 'sxvortx', 'sxvorty', 'sxvortz'/),                                &
    (/ turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % vortxp2,                  &
    turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % vortyp2,                     &
    turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % vortzp2 /) )

call write_parallel_cgns(fname_syvel,nx,ny/2+1,nz- nz_end,nz_tot,        &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                      &
    (/ nx, ny/2+1, (nz-1)*(coord+1) + 1 - nz_end /),                     &
    x(1:nx) , y(1:ny/2+1) , z(1:(nz-nz_end) ), 6,                        &
    (/ 'syuu', 'syvv', 'syww', 'syuv', 'syuw', 'syvw'/),                 &
    (/ turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % upup,                     &
    turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % vpvp,                        &
    turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % wpwp,                        &
    turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % upvp,                        &
    turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % upwp,                        &
    turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % vpwp /) )

call write_parallel_cgns(fname_syvort,nx,ny/2+1,nz- nz_end,nz_tot,       &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                      &
    (/ nx, ny/2+1, (nz-1)*(coord+1) + 1 - nz_end /),                     &
    x(1:nx) , y(1:ny/2+1) , z(1:(nz-nz_end) ), 3,                        &
    (/ 'syvortx', 'syvorty', 'syvortz'/),                                &
    (/ turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % vortxp2,                  &
    turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % vortyp2,                     &
    turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % vortzp2 /) )

#ifdef PPSCALARS
call write_parallel_cgns(fname_sxscal,nx/2+1,ny,nz- nz_end,nz_tot,       &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                      &
    (/ nx/2+1, ny, (nz-1)*(coord+1) + 1 - nz_end /),                     &
    x(1:nx/2+1) , y(1:ny) , z(1:(nz-nz_end) ), 4,                        &
    (/ 'sxtheta2', 'sxutheta', 'sxvtheta', 'sxwtheta'/),                 &
    (/ scal_turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % thetap2,             &
    scal_turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % upthetap,               &
    scal_turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % vpthetap,               &
    scal_turbspecx(1:nx/2+1,1:ny,1:nz- nz_end) % wpthetap /) )

call write_parallel_cgns(fname_syscal,nx,ny/2+1,nz- nz_end,nz_tot,       &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                      &
    (/ nx, ny/2+1, (nz-1)*(coord+1) + 1 - nz_end /),                     &
    x(1:nx) , y(1:ny/2+1) , z(1:(nz-nz_end) ), 4,                        &
    (/ 'sytheta2', 'syutheta', 'syvtheta', 'sywtheta'/),                 &
    (/ scal_turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % thetap2,             &
    scal_turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % upthetap,               &
    scal_turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % vpthetap,               &
    scal_turbspecy(1:nx,1:ny/2+1,1:nz- nz_end) % wpthetap /) )
#endif

#else
! Write binary data

open(unit=13, file=fname_sxvel, form='unformatted', convert=write_endian,   &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) turbspecx(:nx/2+1,:ny,1:nz)%upup
write(13,rec=2) turbspecx(:nx/2+1,:ny,1:nz)%vpvp
write(13,rec=3) turbspecx(:nx/2+1,:ny,1:nz)%wpwp
write(13,rec=4) turbspecx(:nx/2+1,:ny,1:nz)%upvp
write(13,rec=5) turbspecx(:nx/2+1,:ny,1:nz)%upwp
write(13,rec=6) turbspecx(:nx/2+1,:ny,1:nz)%vpwp
close(13)

open(unit=13, file=fname_sxvort, form='unformatted', convert=write_endian,  &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) turbspecx(:nx/2+1,:ny,1:nz)%vortxp2
write(13,rec=2) turbspecx(:nx/2+1,:ny,1:nz)%vortyp2
write(13,rec=3) turbspecx(:nx/2+1,:ny,1:nz)%vortzp2
close(13)

open(unit=13, file=fname_syvel, form='unformatted', convert=write_endian,   &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) turbspecy(:nx,:ny/2+1,1:nz)%upup
write(13,rec=2) turbspecy(:nx,:ny/2+1,1:nz)%vpvp
write(13,rec=3) turbspecy(:nx,:ny/2+1,1:nz)%wpwp
write(13,rec=4) turbspecy(:nx,:ny/2+1,1:nz)%upvp
write(13,rec=5) turbspecy(:nx,:ny/2+1,1:nz)%upwp
write(13,rec=6) turbspecy(:nx,:ny/2+1,1:nz)%vpwp
close(13)

open(unit=13, file=fname_syvort, form='unformatted', convert=write_endian,  &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) turbspecy(:nx,:ny/2+1,1:nz)%vortxp2
write(13,rec=2) turbspecy(:nx,:ny/2+1,1:nz)%vortyp2
write(13,rec=3) turbspecy(:nx,:ny/2+1,1:nz)%vortzp2
close(13)

#ifdef PPSCALARS
open(unit=13, file=fname_sxscal, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) scal_turbspecx(:nx/2+1,:ny,1:nz)%thetap2
write(13,rec=2) scal_turbspecx(:nx/2+1,:ny,1:nz)%upthetap
write(13,rec=3) scal_turbspecx(:nx/2+1,:ny,1:nz)%vpthetap
write(13,rec=4) scal_turbspecx(:nx/2+1,:ny,1:nz)%wpthetap
close(13)

open(unit=13, file=fname_syscal, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) scal_turbspecy(:nx,:ny/2+1,1:nz)%thetap2
write(13,rec=2) scal_turbspecy(:nx,:ny/2+1,1:nz)%upthetap
write(13,rec=3) scal_turbspecy(:nx,:ny/2+1,1:nz)%vpthetap
write(13,rec=4) scal_turbspecy(:nx,:ny/2+1,1:nz)%wpthetap
close(13)

#ifdef PPOUTPUT_SPECBUDG
! Write spectral budget data
allocate(scal_specbudgx(nx/2+1,ny,lbz:nz))
scal_specbudgx = scal_specbudgx_compute(tavg_scal_specbudgx, scal_turbspecx, tavg_scal_turbspecx, tavg_budget, tavg, lbz)

#ifdef PPCGNS
! Not ready
#else
! Write binary data
open(unit=13, file=fname_sxscal_budg, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) scal_specbudgx(:nx/2+1,:ny,1:nz)%adv
write(13,rec=2) scal_specbudgx(:nx/2+1,:ny,1:nz)%tfluc
write(13,rec=3) scal_specbudgx(:nx/2+1,:ny,1:nz)%tvisc
write(13,rec=4) scal_specbudgx(:nx/2+1,:ny,1:nz)%prod
write(13,rec=5) scal_specbudgx(:nx/2+1,:ny,1:nz)%pdiss
close(13)
#endif
#endif

deallocate(scal_turbspecx)
deallocate(scal_turbspecy)
#endif

#endif

#endif

#ifdef PPOUTPUT_SPECBUDG
! Write spectral budget data
allocate(specbudgx(nx/2+1,ny,lbz:nz))
specbudgx = specbudgx_compute(tavg_specbudgx, turbspecx, tavg_turbspecx, tavg_budget, tavg, lbz)

#ifdef PPCGNS
! Not ready
#else
! Write binary data
open(unit=13, file=fname_sxuu, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) specbudgx(:nx/2+1,:ny,1:nz)%advxx
write(13,rec=2) specbudgx(:nx/2+1,:ny,1:nz)%tflucxx
write(13,rec=3) specbudgx(:nx/2+1,:ny,1:nz)%tpresxx
write(13,rec=4) specbudgx(:nx/2+1,:ny,1:nz)%pstrainxx
write(13,rec=5) specbudgx(:nx/2+1,:ny,1:nz)%tviscxx
write(13,rec=6) specbudgx(:nx/2+1,:ny,1:nz)%prodxx
write(13,rec=7) specbudgx(:nx/2+1,:ny,1:nz)%pdissxx
close(13)

open(unit=13, file=fname_sxvv, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) specbudgx(:nx/2+1,:ny,1:nz)%advyy
write(13,rec=2) specbudgx(:nx/2+1,:ny,1:nz)%tflucyy
write(13,rec=3) specbudgx(:nx/2+1,:ny,1:nz)%tpresyy
write(13,rec=4) specbudgx(:nx/2+1,:ny,1:nz)%pstrainyy
write(13,rec=5) specbudgx(:nx/2+1,:ny,1:nz)%tviscyy
write(13,rec=6) specbudgx(:nx/2+1,:ny,1:nz)%prodyy
write(13,rec=7) specbudgx(:nx/2+1,:ny,1:nz)%pdissyy
close(13)

open(unit=13, file=fname_sxww, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) specbudgx(:nx/2+1,:ny,1:nz)%advzz
write(13,rec=2) specbudgx(:nx/2+1,:ny,1:nz)%tfluczz
write(13,rec=3) specbudgx(:nx/2+1,:ny,1:nz)%tpreszz
write(13,rec=4) specbudgx(:nx/2+1,:ny,1:nz)%pstrainzz
write(13,rec=5) specbudgx(:nx/2+1,:ny,1:nz)%tvisczz
write(13,rec=6) specbudgx(:nx/2+1,:ny,1:nz)%prodzz
write(13,rec=7) specbudgx(:nx/2+1,:ny,1:nz)%pdisszz
close(13)

open(unit=13, file=fname_sxuv, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) specbudgx(:nx/2+1,:ny,1:nz)%advxy
write(13,rec=2) specbudgx(:nx/2+1,:ny,1:nz)%tflucxy
write(13,rec=3) specbudgx(:nx/2+1,:ny,1:nz)%tpresxy
write(13,rec=4) specbudgx(:nx/2+1,:ny,1:nz)%pstrainxy
write(13,rec=5) specbudgx(:nx/2+1,:ny,1:nz)%tviscxy
write(13,rec=6) specbudgx(:nx/2+1,:ny,1:nz)%prodxy
write(13,rec=7) specbudgx(:nx/2+1,:ny,1:nz)%pdissxy
close(13)

open(unit=13, file=fname_sxuw, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) specbudgx(:nx/2+1,:ny,1:nz)%advxz
write(13,rec=2) specbudgx(:nx/2+1,:ny,1:nz)%tflucxz
write(13,rec=3) specbudgx(:nx/2+1,:ny,1:nz)%tpresxz
write(13,rec=4) specbudgx(:nx/2+1,:ny,1:nz)%pstrainxz
write(13,rec=5) specbudgx(:nx/2+1,:ny,1:nz)%tviscxz
write(13,rec=6) specbudgx(:nx/2+1,:ny,1:nz)%prodxz
write(13,rec=7) specbudgx(:nx/2+1,:ny,1:nz)%pdissxz
close(13)

open(unit=13, file=fname_sxvw, form='unformatted', convert=write_endian,     &
    access='direct', recl=(nx/2+1)*ny*nz*rprec)
write(13,rec=1) specbudgx(:nx/2+1,:ny,1:nz)%advyz
write(13,rec=2) specbudgx(:nx/2+1,:ny,1:nz)%tflucyz
write(13,rec=3) specbudgx(:nx/2+1,:ny,1:nz)%tpresyz
write(13,rec=4) specbudgx(:nx/2+1,:ny,1:nz)%pstrainyz
write(13,rec=5) specbudgx(:nx/2+1,:ny,1:nz)%tviscyz
write(13,rec=6) specbudgx(:nx/2+1,:ny,1:nz)%prodyz
write(13,rec=7) specbudgx(:nx/2+1,:ny,1:nz)%pdissyz
close(13)
#endif

deallocate(specbudgx)

! Write spectral budget data
allocate(specbudgy(nx,ny/2+1,lbz:nz))
specbudgy = specbudgy_compute(tavg_specbudgy, turbspecy, tavg_turbspecy, tavg_budget, tavg, lbz)

#ifdef PPCGNS
! Not ready
#else
! Write binary data
open(unit=13, file=fname_syuu, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) specbudgy(:nx,:ny/2+1,1:nz)%advxx
write(13,rec=2) specbudgy(:nx,:ny/2+1,1:nz)%tflucxx
write(13,rec=3) specbudgy(:nx,:ny/2+1,1:nz)%tpresxx
write(13,rec=4) specbudgy(:nx,:ny/2+1,1:nz)%pstrainxx
write(13,rec=5) specbudgy(:nx,:ny/2+1,1:nz)%tviscxx
write(13,rec=6) specbudgy(:nx,:ny/2+1,1:nz)%prodxx
write(13,rec=7) specbudgy(:nx,:ny/2+1,1:nz)%pdissxx
close(13)

open(unit=13, file=fname_syvv, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) specbudgy(:nx,:ny/2+1,1:nz)%advyy
write(13,rec=2) specbudgy(:nx,:ny/2+1,1:nz)%tflucyy
write(13,rec=3) specbudgy(:nx,:ny/2+1,1:nz)%tpresyy
write(13,rec=4) specbudgy(:nx,:ny/2+1,1:nz)%pstrainyy
write(13,rec=5) specbudgy(:nx,:ny/2+1,1:nz)%tviscyy
write(13,rec=6) specbudgy(:nx,:ny/2+1,1:nz)%prodyy
write(13,rec=7) specbudgy(:nx,:ny/2+1,1:nz)%pdissyy
close(13)

open(unit=13, file=fname_syww, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) specbudgy(:nx,:ny/2+1,1:nz)%advzz
write(13,rec=2) specbudgy(:nx,:ny/2+1,1:nz)%tfluczz
write(13,rec=3) specbudgy(:nx,:ny/2+1,1:nz)%tpreszz
write(13,rec=4) specbudgy(:nx,:ny/2+1,1:nz)%pstrainzz
write(13,rec=5) specbudgy(:nx,:ny/2+1,1:nz)%tvisczz
write(13,rec=6) specbudgy(:nx,:ny/2+1,1:nz)%prodzz
write(13,rec=7) specbudgy(:nx,:ny/2+1,1:nz)%pdisszz
close(13)

open(unit=13, file=fname_syuv, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) specbudgy(:nx,:ny/2+1,1:nz)%advxy
write(13,rec=2) specbudgy(:nx,:ny/2+1,1:nz)%tflucxy
write(13,rec=3) specbudgy(:nx,:ny/2+1,1:nz)%tpresxy
write(13,rec=4) specbudgy(:nx,:ny/2+1,1:nz)%pstrainxy
write(13,rec=5) specbudgy(:nx,:ny/2+1,1:nz)%tviscxy
write(13,rec=6) specbudgy(:nx,:ny/2+1,1:nz)%prodxy
write(13,rec=7) specbudgy(:nx,:ny/2+1,1:nz)%pdissxy
close(13)

open(unit=13, file=fname_syuw, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) specbudgy(:nx,:ny/2+1,1:nz)%advxz
write(13,rec=2) specbudgy(:nx,:ny/2+1,1:nz)%tflucxz
write(13,rec=3) specbudgy(:nx,:ny/2+1,1:nz)%tpresxz
write(13,rec=4) specbudgy(:nx,:ny/2+1,1:nz)%pstrainxz
write(13,rec=5) specbudgy(:nx,:ny/2+1,1:nz)%tviscxz
write(13,rec=6) specbudgy(:nx,:ny/2+1,1:nz)%prodxz
write(13,rec=7) specbudgy(:nx,:ny/2+1,1:nz)%pdissxz
close(13)

open(unit=13, file=fname_syvw, form='unformatted', convert=write_endian,     &
    access='direct', recl=nx*(ny/2+1)*nz*rprec)
write(13,rec=1) specbudgy(:nx,:ny/2+1,1:nz)%advyz
write(13,rec=2) specbudgy(:nx,:ny/2+1,1:nz)%tflucyz
write(13,rec=3) specbudgy(:nx,:ny/2+1,1:nz)%tpresyz
write(13,rec=4) specbudgy(:nx,:ny/2+1,1:nz)%pstrainyz
write(13,rec=5) specbudgy(:nx,:ny/2+1,1:nz)%tviscyz
write(13,rec=6) specbudgy(:nx,:ny/2+1,1:nz)%prodyz
write(13,rec=7) specbudgy(:nx,:ny/2+1,1:nz)%pdissyz
close(13)
#endif

deallocate(specbudgy)

#endif

#ifdef PPOUTPUT_TURBSEC
deallocate(turbspecx)
deallocate(turbspecy)
#endif

#ifdef PPMPI
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
#endif

end subroutine tavg_finalize

!*******************************************************************************
subroutine tavg_checkpoint()
!*******************************************************************************
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'tavg_finalize' at the end of the
! simulation.
!
use param, only : checkpoint_tavg_file, write_endian
use param, only : checkpoint_tavg_vort_file
use stat_defs, only : tavg_total_time, tavg
use stat_defs, only : tavg_vort
#ifdef PPOUTPUT_SGS
use param, only : checkpoint_tavg_sgs_file
use stat_defs, only : tavg_sgs
#endif
#ifdef PPOUTPUT_BUDGET
use param, only : checkpoint_tavg_budget_file
use stat_defs, only : tavg_budget
#endif
#ifdef PPOUTPUT_TURBSPEC
use param, only : checkpoint_tavg_turbspec_file
use stat_defs, only : tavg_turbspecx, tavg_turbspecy
#endif
#ifdef PPSCALARS
use param, only : checkpoint_tavg_scal_file
use stat_defs, only : tavg_scal
#ifdef PPOUTPUT_TURBSPEC
use param, only : checkpoint_tavg_scal_turbspec_file
use stat_defs, only : tavg_scal_turbspecx, tavg_scal_turbspecy
#endif
#ifdef PPOUTPUT_SPECBUDG
use param, only : checkpoint_tavg_scal_specbudg_file
use stat_defs, only : tavg_scal_specbudgx
#endif
#endif
#ifdef PPOUTPUT_SPECBUDG
use param, only : checkpoint_tavg_specbudg_file
use stat_defs, only : tavg_specbudgx, tavg_specbudgy
#endif

implicit none

character(64) :: fname

fname = checkpoint_tavg_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif

!  Write data to tavg.out
open(1, file=fname, action='write', position='rewind',form='unformatted',   &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg
close(1)

fname = checkpoint_tavg_vort_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
! Write data to tavg_vort.out
open(1, file=fname, action='write', position='rewind',form='unformatted',      &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_vort
close(1)

#ifdef PPOUTPUT_SGS
fname = checkpoint_tavg_sgs_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
!  Write data to tavg_sgs.out
open(1, file=fname, action='write', position='rewind',form='unformatted',   &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_sgs
close(1)
#endif

#ifdef PPOUTPUT_BUDGET
fname = checkpoint_tavg_budget_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
!  Write data to tavg_budget.out
open(1, file=fname, action='write', position='rewind',form='unformatted',   &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_budget
close(1)
#endif

#ifdef PPOUTPUT_TURBSPEC
fname = checkpoint_tavg_turbspec_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
!  Write data to tavg_turbspec.out
open(1, file=fname, action='write', position='rewind',form='unformatted',   &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_turbspecx
write(1) tavg_turbspecy
close(1)
#endif

#ifdef PPSCALARS
fname = checkpoint_tavg_scal_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
!  Write data to tavg_sgs.out
open(1, file=fname, action='write', position='rewind',form='unformatted',      &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_scal
close(1)

#ifdef PPOUTPUT_TURBSPEC
fname = checkpoint_tavg_scal_turbspec_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
!  Write data to tavg_turbspec.out
open(1, file=fname, action='write', position='rewind',form='unformatted',      &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_scal_turbspecx
write(1) tavg_scal_turbspecy
close(1)
#endif

#ifdef PPOUTPUT_SPECBUDG
fname = checkpoint_tavg_scal_specbudg_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
!  Write data to tavg_turbspec.out
open(1, file=fname, action='write', position='rewind',form='unformatted',      &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_scal_specbudgx
close(1)
#endif
#endif

#ifdef PPOUTPUT_SPECBUDG
fname = checkpoint_tavg_specbudg_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif
!  Write data to tavg_specbudg.out
open(1, file=fname, action='write', position='rewind',form='unformatted',   &
    convert=write_endian)
write(1) tavg_total_time
write(1) tavg_specbudgx
write(1) tavg_specbudgy
close(1)
#endif

end subroutine tavg_checkpoint

end module io
