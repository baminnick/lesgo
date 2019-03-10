subroutine load_jacobian()
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
use types

implicit none
real(kind=rprec),dimension(nz_tot) ::FIELD1
real(kind=rprec),dimension(nz_tot) ::FIELD2
real(kind=rprec),dimension(nz_tot) ::FIELD3

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

! Specify which z-levels should be run in RNL-fourier mode
! NOTE: since mesh-stretch is on uv-grid, the interface is
! always located at a uv-point
if ((hybrid_baseline) .or. (hybrid_natural)) then
    do jz = lbz, nz
        if (mesh_stretch(jz) .le. hwm) then !! z-level in RNL-Fourier mode
            zhyb(jz) = .true.
        else !! z-level in LES-Physical mode
            zhyb(jz) = .false.
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

    ! Tell user location of interface
    do jz = 2, nz-2
        if ( (zhyb(jz)) .and. (.not. zhyb(jz+1)) ) then
            write(*,*) '--> Hybrid Fourier: interface at, hwm = ', mesh_stretch(jz)
        endif
    enddo

endif

end subroutine load_jacobian

