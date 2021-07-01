! Solves the 2-dimensional diffusion equation

! sujeta a condicion inicial: T(x,0)=20ºC (en todo el plato),

! definicion de variables
program difusion
implicit none

! definicion de variables
integer :: nx, ny, i, tt, lx, xl, i_sample,n_serie,it,j
real :: lambdax, dx, nx1, ny1, ly,lambday, yl, dy
real, allocatable, dimension (:) :: x,y
real, allocatable, dimension (:,:) :: T,T0

character(16) outputfield
character(2) serie
character(3) sample

nx=1000.0 ! tamaño de la malla direccion x
ny=1000.0 ! tamaño de la malla direccion y
nx1=500.0 ! tamaño de la malla direccion x flotante
ny1=500.0 ! tamaño de la malla direccion y flotante
lx=1.0; ! longitud x
ly=1.0 ! longitud y
xl=0 ! coordenada fisica extremo izquierdo
yl=0 !coordenada fisica extremo
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

! condiciones iniciales
do i=1,nx
  do j=1,ny
    T0(i,j)=20.0 ! todo el plato a 20ºC
  end do
end do

! tiempo
do tt=1,10000

  ! vector x (posicion)
  do i=1,nx
    x(i)=xl+i*lx/(nx1)
  end do

  do i=1,ny
    y(i)=yl+i*ly/(ny1)
  end do

  ! escribir archivos cada 100 unidades de tiempo
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

  ! metodo explicito
  do i=2,nx-1
    do j=2,ny-1
      T(i,j)=T0(i,j)+lambdax*(T0(i-1,j)-2*T0(i,j)+T0(i+1,j))+lambday*(T0(i,j-1)-2*T0(i,j)+T0(i,j+1))
    end do
  end do

  ! aplicamos condiciones de frontera
  do i=1,nx
    T(i,1)=30
    T(i,ny)=60
  end do
  do i=1,ny
    T(1,j)=20
    T(ny,j)=50
  end do

  ! sobreescribe T a T0. Para tener a T0 como T^(n), y a T como T^(n+1)
  do i=1,nx
    do j=1,ny
      T0(i,j)=T(i,j)
    end do
  end do

end do

100 format(f12.8,f12.8,f12.8)

end program difusion
