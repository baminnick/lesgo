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
subroutine init_hybrid()
!*******************************************************************************
! 
! This subroutine initializes variables used for hybrid_fourier mode, 
! specifically the interface is defined.
! 
use param, only : hwm, lbz, nz, coord, nproc
use sim_param, only : zhyb, coord_int, jz_int, coord_intv
#ifdef PPMAPPING
use sim_param, only : mesh_stretch
#else
use grid_m
#endif
use param, only : comm, ierr, mpi_rprec
use mpi
use types, only : rprec

implicit none

integer:: jz
real(rprec), dimension(lbz:nz) :: z_uv
real(rprec), dimension(:), allocatable :: coord_intv_temp

allocate( coord_intv_temp(nproc) )
allocate( coord_intv(nproc) )

! Specify uv-locations to be used
#ifdef PPMAPPING
z_uv(:) = mesh_stretch(:)
#else
! Nullify pointers
nullify(z)

! Set grid pointers
z_uv => grid % z
#endif

! Specify which z-levels should be run in RNL-fourier mode
! NOTE: the interface always exists on uv-node.
! 
! coord_int = -1 if coord below interface, coord_int = 1 if 
! coord above interface, coord_int = 0 if coord owns interface.
! jz_int = -1 for all z-levels, except where the interface is.
! coord_int and jz_int is only used in tridag_array_hybrid_phys.f90
! 
! Currently hybrid_fourier initialized to run where one coord
! owns the interface, however seems to work OK if interface
! is at jz = 1 and jz = nz for two processors. Does not work 
! if interface is at jz = 0 and jz = nz-1.
do jz = lbz, nz
    if (z_uv(jz) <= hwm) then !! z-level in fourier mode
        zhyb(jz) = .true.
        coord_int = -1
    else !! z-level in physical mode
        zhyb(jz) = .false.
        if (jz == lbz) then
            coord_int = 1
        endif
    endif
enddo

! Correct location of interface if being shared between processors
! This forces the interface to be at jz = 2 for the processor above
if (zhyb(nz-1) .and. (.not. zhyb(nz))) then
    zhyb(nz) = .true. !! entire processor now in fourier mode
elseif (zhyb(lbz) .and. (.not. zhyb(1))) then
    zhyb(1) = .true.
    zhyb(2) = .true. !! new location of interface
elseif (zhyb(1) .and. (.not. zhyb(2))) then
    zhyb(2) = .true. !! new location of interface
endif

! Enforce part of the grid to be physical/fourier
! if specified hwm is too small or too large
if ((coord == 0) .and. (.not. zhyb(2))) then !! hwm too small
    zhyb(lbz) = .true.
    zhyb(1) = .true.
    zhyb(2) = .true.
elseif ((coord == nproc-1) .and. (zhyb(nz-1))) then !! hwm too large
    zhyb(nz) = .false.
    zhyb(nz-1) = .false.
endif

! Find the processor which owns the interface and jz-location
! jz_int is initialized as -1 on all processors
do jz = 2, nz-2
    if ( (zhyb(jz)) .and. (.not. zhyb(jz+1)) ) then
        jz_int = jz
        coord_int = 0
        write(*,*) '--> Hybrid Fourier: interface at, hwm = ', z_uv(jz)
    endif
enddo

! Let coord = 0 know which processor owns the interface
! This is only for reporting iteration time of wall-model
call mpi_gather(real(coord_int,rprec), 1, MPI_RPREC,                     &
     coord_intv_temp, 1, MPI_RPREC, 0, comm, ierr)
if (coord == 0) then
    coord_intv = int( coord_intv_temp(:) )
endif

! Output which z-levels are in fourier mode for debugging
!do jz = lbz, nz
!write(*,*) coord, jz, z_uv(jz), zhyb(jz)
!enddo
!write(*,*) '----------------------'

end subroutine init_hybrid

