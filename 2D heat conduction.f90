! Solves the 2D transient heat conduction equation
! case: square plate, subjected to:
! initial conditions: 20ºC (everywhere)
! boundary conditions: 60ºC (top), 30ºC (bottom), 90ºC (left), 5ºC (right)

! it adittionally creates vtk files for post-processing (using ParaView)

program conduction

integer, parameter :: nx=100 ,ny=100,nz=1
integer :: i,j
integer :: tt,ts,totalt
integer :: i_sample, it
real :: dx,dy
real :: k
real :: T0,Ttop,Tbottom,Tright,Tleft
real,allocatable,dimension (:,:,:) :: T
real,allocatable,dimension (:) :: x,y,z

character(16) outputfield
character(3) sample

i_sample=0
it=0

lx=1.0
ly=1.0

k=0.005

T0=20.0

Ttop=60.0
Tbottom=30.0
Tleft=90.0
Tright=5.0

totalt=90.0
ts=5000
dt=0.001 !totalt/ts

allocate(T(nx,ny,ts))
allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

dx=lx/float(nx-1)
dy=ly/float(ny-1)

! position vectors
do i=1,nx
  x(i)=dx*float(i-1)
end do
do i=1,ny
  y(i)=dy*float(i-1)
end do
do i=1,nz
  z(i)=0.0
end do

!applies initial conditions
do j=1,ny
  do i=1,nx
    T(i,j,1)=T0
  end do
end do

! applies boundary conditions
do tt=1,ts
  do j=1,ny-1
    T(1,j,tt)=Tleft
    T(nx,j,tt)=Tright
  end do
  do i=2,nx
    T(i,1,tt)=Tbottom
    T(i,ny,tt)=Ttop
  end do
end do

! applies a finite difference scheme
do tt=1,ts-1
  do j=2,ny-1
    do i=2,nx-1
      T(i,j,tt+1)=T(i,j,tt)+k*dt/(dx**2)*(T(i+1,j,tt)-2*T(i,j,tt)+T(i-1,j,tt)) &
                           +k*dt/(dy**2)*(T(i,j+1,tt)-2*T(i,j,tt)+T(i,j-1,tt))
    end do
  end do
end do

! creates vtk files
do tt=1,ts
  if (it.eq.50) then
    i_sample=i_sample+1
    write(sample,'(i3.3)')(i_sample)
    outputfield='temp_'//sample//'.vtk'

    open(78, file = Outputfield, form = 'formatted')

      write(78,110) '# vtk DataFile Version 2.3'
      write(78,110) '3D Mesh'
      write(78,111) 'ASCII'

      write(78,110) 'DATASET STRUCTURED_GRID'
      write(78,120) 'DIMENSIONS',nx, ny, nz
      write(78,130) 'POINTS', nx*ny, ' float'

      do j=1,ny
        do i=1,nx
          write(78,100) x(i), y(j), z(1)
        enddo
      enddo

      write(78,140) 'POINT_DATA', nx*ny

      write(78,110) 'SCALARS Temperature float'
      write(78,110) 'LOOKUP_TABLE default'

        do j=1,ny
          do i=1,nx
            write(78,100) T(i,j,tt)
          end do
        end do

    close(78)

    it=1
  else
    it=it+1
  end if

end do

100  FORMAT(3(f12.6));
110  FORMAT(A);
111  FORMAT(A,/);
120  FORMAT(A,I4,I4,I4);
130  FORMAT(A,I10,A);
140  FORMAT(A,I10);

end program conduction
