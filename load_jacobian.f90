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
#ifdef PPLVLSET_STRETCH
real(kind=rprec),dimension(nz_tot) ::FIELD4
#endif

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

#ifdef PPLVLSET_STRETCH
    open(4,file=path//'jaco104.dat')
    do i=1,nz_tot
        read(4,*) FIELD4(i)
    end do
    close(4)
#endif

else
    ! Create unstretched grid to be mapped to new grid
    do i = 1, nz_tot
        z_w(i) = (i-1)*dz !! z-locations on w-grid
    enddo
    z_uv(:) = z_w(:) + 0.5_rprec*dz !! z-locations on uv-grid

    if (ubc_mom == 0) then
        ! Map to stretched grid
        ! For the uv-grid
        FIELD3(:) = L_z*(1.0_rprec+(tanh(str_factor*(z_uv(:)/L_z-1.0_rprec)) &
            /tanh(str_factor)))
        ! For the w-grid
#ifdef PPLVLSET_STRETCH
        FIELD4(:) = L_z*(1.0_rprec+(tanh(str_factor*(z_w(:)/L_z-1.0_rprec))  &
            /tanh(str_factor)))
#endif

        ! Compute Jacobian values for both w- and uv-grids
        ! Using analytical derivative expression
        FIELD1(:) = L_z*(str_factor/L_z)*                                 &
            (1-(tanh(str_factor*(z_w(:)/L_z-1)))**2)/tanh(str_factor)
        FIELD2(:) = L_z*(str_factor/L_z)*                                 &
            (1-(tanh(str_factor*(z_uv(:)/L_z-1)))**2)/tanh(str_factor)
    else !! full channel
        ! Map to stretched grid
        ! For the uv-grid
        FIELD3(:) = L_z*0.5_rprec*(1.0_rprec+(tanh(str_factor*(z_uv(:)/L_z-0.5_rprec))    &
            /tanh(0.5_rprec*str_factor)))
        ! For the w-grid
#ifdef PPLVLSET_STRETCH
        FIELD4(:) = L_z*0.5_rprec*(1.0_rprec+(tanh(str_factor*(z_w(:)/L_z-0.5_rprec))     &
            /tanh(0.5_rprec*str_factor)))
#endif

        ! Compute Jacobian values for both w- and uv-grids
        ! Using analytical derivative expression
        FIELD1(:) = L_z*0.5_rprec*(str_factor/L_z)*                    &
            (1-(tanh(str_factor*(z_w(:)/L_z-0.5_rprec)))**2)/tanh(0.5_rprec*str_factor)
        FIELD2(:) = L_z*0.5_rprec*(str_factor/L_z)*                    &
            (1-(tanh(str_factor*(z_uv(:)/L_z-0.5_rprec)))**2)/tanh(0.5_rprec*str_factor)

    endif !! half/full channel

endif

! Store variables into what LESGO will use
do jz=1,nz
#ifdef PPHYBRID
    if (fourier) then !! RNL WM coords, nz is larger here
        i = coord*(nz-1) + jz
    else !! non RNL coords, need to shift accordingly
        i = nz_tot_rnl + (coord-nproc_rnl)*(nz-1) + jz
    endif
#else
    i = coord*(nz-1) + jz
#endif
    JACO1(jz) = FIELD1(i)
    JACO2(jz) = FIELD2(i)
    mesh_stretch(jz) = FIELD3(i)
#ifdef PPLVLSET_STRETCH
    mesh_stretch_w(jz) = FIELD4(i)
#endif
end do

if (coord == 0) then
    JACO1(lbz)=JACO1(1)
    JACO2(lbz)=JACO2(1)
    mesh_stretch(lbz)=-mesh_stretch(1)
#ifdef PPLVLSET_STRETCH
    mesh_stretch_w(lbz)=-mesh_stretch_w(1)
#endif
    write(*,*) '--> Grid stretched using mapping function'
else
    JACO1(lbz)=FIELD1((coord-1)*(nz-1)+nz-1)
    JACO2(lbz)=FIELD2((coord-1)*(nz-1)+nz-1)
    mesh_stretch(lbz)=FIELD3((coord-1)*(nz-1)+nz-1)
#ifdef PPLVLSET_STRETCH
    mesh_stretch_w(lbz)=FIELD4((coord-1)*(nz-1)+nz-1)
#endif
end if

#ifdef PPHYBRID
! Tell user where the interface is located
if (coord == (nproc_rnl-1)) then
    write(*,*) 'HYBRID: Interface located at z =', FIELD3(coord*(nz-1) + nz)
endif
#endif

end subroutine load_jacobian

