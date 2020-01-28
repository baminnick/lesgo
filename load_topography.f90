subroutine load_topography ()
use sim_param
use param
use types

implicit none
real(kind=rprec),dimension(nx*ny)  ::FIELD1
!real(kind=rprec),dimension(nx*ny/4)::FIELD2
!character(120)::ffname1
integer:: i, j, k
real(rprec) :: z_uv, z_w

open(1,file=path//'hij.dat')

do i=1,nx*ny
read(1,*) FIELD1(i)
end do
close(1)
do i=1,nx
do j=1,ny
    hij(i,j) = FIELD1((i-1)*ny+j)
end do
end do

do i=1,nx
do j=1,ny
do k=lbz,nz
#ifdef PPMAPPING
    phi_uv(i,j,k) = mesh_stretch(k)-hij(i,j)
    phi_w(i,j,k) = mesh_stretch_w(k)-hij(i,j) 
#else
    z_uv = (coord*(nz-1) + k - 0.5_rprec) * dz
    z_w = (coord*(nz-1) + k - 1.0_rprec) * dz
    phi_uv(i,j,k) = z_uv - hij(i,j)
    phi_w(i,j,k) = z_w - hij(i,j)
#endif
end do
end do
end do

end subroutine load_topography
