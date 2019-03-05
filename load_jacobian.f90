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
subroutine load_jacobian()
!*******************************************************************************
! 
! This subroutine gives the variables necessary to run LESGO with 
! a stretched grid in the wall-normal direction. Based on user
! input, the grid may be loaded via external files or created.
! 
! JACO1        - Jacobian on the w-grid
! JACO2        - Jacobian on the uv-grid
! mesh_stretch - z-locations on the uv-grid
! 
use sim_param
use param
use types, only : rprec

implicit none
real(rprec),dimension(nz_tot) ::FIELD1
real(rprec),dimension(nz_tot) ::FIELD2
real(rprec),dimension(nz_tot) ::FIELD3
real(rprec),dimension(nz_tot) ::FIELD4

real(rprec), dimension(nz_tot) :: z_w, z_uv

integer:: i,jz

if (load_stretch) then
    ! Load user-supplied files
    open(1,file=path//'jaco101.dat')
    do i=1,nz_tot
        read(1,*) FIELD1(i)
    end do
    close(1)

    open(2,file=path//'jaco102.dat')
    do i=1,nz_tot
        read(2,*) FIELD2(i)
    end do
    close(2)

    open(3,file=path//'jaco103.dat')
    do i=1,nz_tot
        read(3,*) FIELD3(i)
    end do
    close(3)

    open(4,file=path//'jaco104.dat')
    do i=1,nz_tot
        read(4,*) FIELD4(i)
    end do
    close(4)

else
    ! Use tanh stretched grid
    ! Create unstretched grid to be mapped to new grid
    do i = 1, nz_tot
        z_w(i) = (i-1)*dz !! z-locations on w-grid
    enddo
    z_uv(:) = z_w(:) + 0.5_rprec*dz !! z-locations on uv-grid

    ! Map to stretched grid
    ! Only need to do this for the uv-grid
    FIELD3(:) = L_z*(1.0_rprec+(tanh(str_factor*(z_uv(:)/L_z-1.0_rprec))    &
        /tanh(str_factor)))

    ! Compute Jacobian values for both w- and uv-grids
    ! Using analytical derivative expression
    FIELD1(:) = L_z*(str_factor/L_z)*                                       &
        (1-(tanh(str_factor*(z_w(:)/L_z-1)))**2)/tanh(str_factor)
    FIELD2(:) = L_z*(str_factor/L_z)*                                       &
        (1-(tanh(str_factor*(z_uv(:)/L_z-1)))**2)/tanh(str_factor)

    ! Compute d(1/J)/dz term for PPE, only on uv-grid
    FIELD4(:) = 2.0_rprec*tanh(str_factor)*                                 &
        tanh(str_factor*(z_uv(:)/L_z-1))/                                   &
        ((1-tanh(str_factor*(z_uv(:)/L_z-1))**2))

endif

! Store variables into what LESGO will use
do jz=1,nz
    JACO1(jz) = FIELD1(coord*(nz-1)+jz)
end do

do jz=1,nz
    JACO2(jz) = FIELD2(coord*(nz-1)+jz)
end do

do jz=1,nz
    mesh_stretch(jz) = FIELD3(coord*(nz-1)+jz)
end do

if (coord == 0) then
    JACO1(lbz)=JACO1(1)    
    JACO2(lbz)=JACO2(1)
    mesh_stretch(lbz)=mesh_stretch(1)     
    write(*,*) '--> Grid stretched using mapping function'
else
    JACO1(lbz)=FIELD1((coord-1)*(nz-1)+nz-1)
    JACO2(lbz)=FIELD2((coord-1)*(nz-1)+nz-1)
    mesh_stretch(lbz)=FIELD3((coord-1)*(nz-1)+nz-1)
end if

end subroutine load_jacobian

