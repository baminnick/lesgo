subroutine load_jacobian()
use sim_param
use param
use types
!xiaowei jan20_2019. check ok.
implicit none
real(kind=rprec),dimension(nz_tot) ::FIELD1
real(kind=rprec),dimension(nz_tot) ::FIELD2
real(kind=rprec),dimension(nz_tot) ::FIELD3
real(kind=rprec),dimension(nz_tot) ::FIELD4

integer:: i,jz

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

do jz=1,nz
 JACO1(jz) = FIELD1(coord*(nz-1)+jz)
end do

do jz=1,nz
 JACO2(jz) = FIELD2(coord*(nz-1)+jz)
end do

do jz=1,nz
  mesh_stretch(jz) = FIELD3(coord*(nz-1)+jz)
end do

do jz=1,nz
  dj_dzeta(jz) = FIELD4(coord*(nz-1)+jz)
end do

if (coord == 0) then
JACO1(lbz)=JACO1(1)    
JACO2(lbz)=JACO2(1)
mesh_stretch(lbz)=mesh_stretch(1)     
dj_dzeta(lbz)=dj_dzeta(1)
else
JACO1(lbz)=FIELD1((coord-1)*(nz-1)+nz-1)
JACO2(lbz)=FIELD2((coord-1)*(nz-1)+nz-1)
mesh_stretch(lbz)=FIELD3((coord-1)*(nz-1)+nz-1)
dj_dzeta(lbz)=FIELD4((coord-1)*(nz-1)+nz-1)
end if
end subroutine load_jacobian
