subroutine level_set2()
!subroutine level_set2(IBFx1,IBFy1,IBFz1)
use param, only : tadv1, dt, BOGUS, fourier, lbz, nx, ny, nz  !--in addition to param vars above
use sim_param
use derivatives, only : wave2phys, phys2wave
implicit none
integer :: i, j, k
real(rprec) :: Rx, Ry, Rz
!real(rprec) :: const
!real(rprec), dimension(nx,ny,nz-1), intent(inout) :: IBFx1,IBFy1,IBFz1

! Transform velocity and pressure gradients to impose force in physical
if (fourier) then
    call wave2phys( u, lbz )
    call wave2phys( v, lbz )
    call wave2phys( w, lbz )
    call wave2phys( dpdx, 1 )
    call wave2phys( dpdy, 1 )
    call wave2phys( dpdz, 1 )
endif

do k = 1, nz
do j = 1, ny
do i = 1, nx

    if (phi_uv(i, j, k) <= 0._rprec) then  !--uv-nodes
        Rx = -tadv1 * dpdx(i, j, k)
        Ry = -tadv1 * dpdy(i, j, k)

        IBFx(i, j, k) = (-u(i, j, k)/dt - Rx)
        IBFy(i, j, k) = (-v(i, j, k)/dt - Ry)
     endif

     if (phi_w(i, j, k) <= 0._rprec) then  !--w-nodes
        Rz = -tadv1 * dpdz(i, j, k)
        IBFz(i, j, k) = (-w(i, j, k)/dt - Rz)
     endif

enddo
enddo
enddo

! Transform velocity, pressure gradient, and force to add to RHS in fourier
if (fourier) then
    call phys2wave( u, lbz )
    call phys2wave( v, lbz )
    call phys2wave( w, lbz )
    call phys2wave( dpdx, 1 )
    call phys2wave( dpdy, 1 )
    call phys2wave( dpdz, 1 )
    call phys2wave( IBFx, 1 )
    call phys2wave( IBFy, 1 )
    call phys2wave( IBFz, 1 )
endif

end subroutine level_set2
