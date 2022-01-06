! Solves the 2-dimensional diffusion equation, using a finite difference method

! with the initial condition: T(x,y,0)=20ºC,

program diffusion
implicit none

! sets up variables used
integer :: nx, ny, i, tt, lx, xl, i_sample,n_serie,it,j
real :: lambdax, dx, nx1, ny1, ly,lambday, yl, dy
real, allocatable, dimension (:) :: x,y
real, allocatable, dimension (:,:) :: T,T0

character(16) outputfield
character(2) serie
character(3) sample

nx=1000.0 ! mesh size, x-direction
ny=1000.0 ! mesh size, y-direction
nx1=500.0 ! mesh size, x-direction (float)
ny1=500.0 ! mesh size, y-direction (float)
lx=1.0; ! length in x
ly=1.0 ! length in y
xl=0 ! left physical coordinate
yl=0
lambdax=0.2
lambday=0.2

dx=lx/(nx1-1)
dy=ly/(ny1-1)

i_sample=0
n_serie=10.0
it=1.0

allocate(T(nx,ny))
allocate(T0(nx,ny))
allocate(x(nx))
allocate(y(ny))

! initial condition
do i=1,nx
  do j=1,ny
    T0(i,j)=20.0 ! todo el plato a 20ºC
  end do
end do

! time
do tt=1,10000

  ! position vector, x and y
  do i=1,nx
    x(i)=xl+i*lx/(nx1)
  end do

  do i=1,ny
    y(i)=yl+i*ly/(ny1)
  end do

  ! writes to disk 100 files
  if (it.eq.100) then
    i_sample=i_sample+1
    write(sample,'(i3.3)')i_sample
    write(serie,'(i2.2)')n_serie
    Outputfield='field_'//serie//'.'//sample//'.txt'
    open(25,file=Outputfield,form='formatted')
    do i=1,nx
      do j=1,ny
        write(25,100)x(i),y(j),T(i,j)
      end do
    end do
    close(25)
    it=1
  else
    it=it+1
  end if

  ! explicit method
  do i=2,nx-1
    do j=2,ny-1
      T(i,j)=T0(i,j)+lambdax*(T0(i-1,j)-2*T0(i,j)+T0(i+1,j))+lambday*(T0(i,j-1)-2*T0(i,j)+T0(i,j+1))
    end do
  end do

  ! applies boundary conditions
  do i=1,nx
    T(i,1)=30
    T(i,ny)=60
  end do
  do i=1,ny
    T(1,j)=20
    T(ny,j)=50
  end do

  do i=1,nx
    do j=1,ny
      T0(i,j)=T(i,j)
    end do
  end do

end do

100 format(f12.8,f12.8,f12.8)

end program diffusion
