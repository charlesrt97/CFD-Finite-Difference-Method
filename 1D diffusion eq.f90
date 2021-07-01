! Solves the 1-dimensional diffusion equation, using a finite difference method

! sujeta a condicion inicial: T(x,0)=20ºC (en todo el alambre),
! condiciones de frontera: extremo izquierdo T(1,0)=60ºC, derecho T(nx,0)=20ºC

! definicion de variables
program difusion
implicit none

! definicion de variables
integer :: nx, i, tt, lx, xl, i_sample,n_serie,it
real :: lambda, dx, nx1
real, allocatable, dimension (:) :: T,T0,x

character(16) outputfield
character(2) serie
character(3) sample

nx=100.0 ! tamaño de la malla
nx1=100.0 ! tamaño de la malla flotante
lx=10.0; ! longitud del alambre
xl=0 ! coordenada fisica extremo izquierdo
lambda=0.4

i_sample=0
n_serie=10.0
it=1.0

allocate(T(nx))
allocate(T0(nx))
allocate(x(nx))

! condiciones iniciales
do i=1,nx
  T0(i)=20.0
end do

! tiempo
do tt=1,1000

  ! vector x (posicion)
  do i=1,nx
    x(i)=xl+i*lx/(nx1)
  end do

  ! escribir archivos cada 100 unidades de tiempo
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

  ! metodo explicito
  do i=2,nx-1
    T(i)=lambda*T0(i-1)+(1.0-2.0*lambda)*T0(i)+lambda*T0(i+1)
  end do

  ! aplicamos condiciones de frontera
  T(1)=60.0
  T(nx)=20.0

  ! sobreescribe T a T0. Para tener a T0 como T^(n), y a T como T^(n+1)
  do i=1,nx
    T0(i)=T(i)
  end do

end do

100 format(f12.8,f12.8)

end program difusion
