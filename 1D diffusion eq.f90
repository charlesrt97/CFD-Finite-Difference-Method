! Solves the 1-dimensional diffusion equation, using a finite difference method

! with the initial condition: T(x,0)=20ºC
! boundary conditions: left boundary T(1,0)=60ºC, right boundary T(nx,0)=20ºC

program diffusion
implicit none

! sets up variables used
integer :: nx, i, tt, lx, xl, i_sample,n_serie,it
real :: lambda, dx, nx1
real, allocatable, dimension (:) :: T,T0,x

character(16) outputfield
character(2) serie
character(3) sample

nx=100.0 ! mesh size
nx1=100.0 ! mesh size (float)
lx=10.0; ! wire's length
xl=0 ! left physical coordinate
lambda=0.4

i_sample=0
n_serie=10.0
it=1.0

allocate(T(nx))
allocate(T0(nx))
allocate(x(nx))

! initial conditions
do i=1,nx
  T0(i)=20.0
end do

! time
do tt=1,1000

  ! position vector, x
  do i=1,nx
    x(i)=xl+i*lx/(nx1)
  end do

  ! writes to disk 100 files
  if (it.eq.100) then
    i_sample=i_sample+1
    write(sample,'(i3.3)')i_sample
    write(serie,'(i2.2)')n_serie
    Outputfield='field_'//serie//'.'//sample//'.txt'
    open(25,file=Outputfield,form='formatted')
    do i=1,nx
      write(25,100)x(i),T(i)
    end do
    close(25)
    it=1
  else
    it=it+1
  end if

  ! explicit method
  do i=2,nx-1
    T(i)=lambda*T0(i-1)+(1.0-2.0*lambda)*T0(i)+lambda*T0(i+1)
  end do

  ! applies boundary conditions
  T(1)=60.0
  T(nx)=20.0
  
  do i=1,nx
    T0(i)=T(i)
  end do

end do

100 format(f12.8,f12.8)

end program diffusion
