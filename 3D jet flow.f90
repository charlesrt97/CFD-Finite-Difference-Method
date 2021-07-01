! INCOMPLETE

!-----------------------------------------------------------------------
module dimensiones

  integer nx,ny,nz,nd,ne

end  module dimensiones

!-----------------------------------------------------------------------
module deltas

  real :: deltax,deltay,deltaz

end  module deltas

!----------------------------------------------------------------------

module consderper

  real, allocatable, dimension (:) :: axp,ayp,azp
  real, allocatable, dimension (:) :: bxp,byp,bzp
  real, allocatable, dimension (:) :: cxp,cyp,czp

end  module consderper

!----------------------------------------------------------------------

module consdernper

  real, allocatable, dimension (:) :: axnp,aynp,aznp
  real, allocatable, dimension (:) :: bxnp,bynp,bznp
  real, allocatable, dimension (:) :: cxnp,cynp,cznp

end  module consdernper

!----------------------------------------------------------------------

module consrs

  real ars,brs,ars1,brs1,crs1,ars2,drs1

end module consrs

!----------------------------------------------------------------------

module derivtools

  real, allocatable, dimension (:) :: du,dv,dw
  real, allocatable, dimension (:) :: u,dup,dunp

end module derivtools

!-----------------------------------------------------------------------

subroutine deriv_alloc()

        use dimensiones
        use consderper
        use consdernper
        use derivtools

        IMPLICIT NONE


        allocate(u(nx))
        allocate(dup(nx))
        allocate(dunp(nx))

        allocate(axp(nx))
        allocate(bxp(nx))
        allocate(cxp(nx))
        allocate(ayp(ny))
        allocate(byp(ny))
        allocate(cyp(ny))
        allocate(azp(nz))
        allocate(bzp(nz))
        allocate(czp(nz))
        allocate(axnp(nx))
        allocate(bxnp(nx))
        allocate(cxnp(nx))
        allocate(aynp(ny))
        allocate(bynp(ny))
        allocate(cynp(ny))
        allocate(aznp(nz))
        allocate(bznp(nz))
        allocate(cznp(nz))
        allocate(du(nx))
        allocate(dv(ny))
        allocate(dw(nz))

        return
        end subroutine deriv_alloc

!-----------------------------------------------------------------------
        subroutine deriv_dealloc()

        use dimensiones
        use consderper
        use consdernper
        use derivtools

        IMPLICIT NONE


        deallocate(u)
        deallocate(dup)
        deallocate(dunp)

        deallocate(axp)
        deallocate(bxp)
        deallocate(cxp)
        deallocate(ayp)
        deallocate(byp)
        deallocate(cyp)
        deallocate(azp)
        deallocate(bzp)
        deallocate(czp)
        deallocate(axnp)
        deallocate(bxnp)
        deallocate(cxnp)
        deallocate(aynp)
        deallocate(bynp)
        deallocate(cynp)
        deallocate(aznp)
        deallocate(bznp)
        deallocate(cznp)
        deallocate(du)
        deallocate(dv)
        deallocate(dw)

        return
        end subroutine deriv_dealloc

!-----------------------------------------------------------------------

        SUBROUTINE inideriv()

        use consderper
        use consdernper
        use consrs
        use dimensiones
        IMPLICIT NONE
        integer i,j,k,l,m

        do i=1,nx
         bxnp(i)=1.0
        enddo
        do i=1,ny
         bynp(i)=1.0
        enddo
        do i=1,nz
         bznp(i)=1.0
        enddo

        axnp(1)=0.0
        axnp(2)=0.25
        axnp(nx-1)=0.25
        axnp(nx)=3.0
        do i=3,nx-2
         axnp(i)=1.0/3.0
        enddo
        aynp(1)=0.0
        aynp(2)=0.25
        aynp(ny-1)=0.25
        aynp(ny)=3.0
        do i=3,ny-2
         aynp(i)=1.0/3.0
        enddo
        aznp(1)=0.0
        aznp(2)=0.25
        aznp(nz-1)=0.25
        aznp(nz)=3.0
        do i=3,nz-2
         aznp(i)=1.0/3.0
        enddo

        cxnp(1)=3.0
        cxnp(2)=0.25
        cxnp(nx-1)=0.25
        cxnp(nx)=0.0
        do i=3,nx-2
         cxnp(i)=1.0/3.0
        enddo
        cynp(1)=3.0
        cynp(2)=0.25
        cynp(ny-1)=0.25
        cynp(ny)=0.0
        do i=3,ny-2
         cynp(i)=1.0/3.0
        enddo
        cznp(1)=3.0
        cznp(2)=0.25
        cznp(nz-1)=0.25
        cznp(nz)=0.0
        do i=3,nz-2
         cznp(i)=1.0/3.0
        enddo

        do i=1,nx
         bxp(i)=1.0
         axp(i)=1.0/3.0
         cxp(i)=1.0/3.0
        enddo
        do i=1,ny
         byp(i)=1.0
         ayp(i)=1.0/3.0
         cyp(i)=1.0/3.0
        enddo
        do i=1,nz
         bzp(i)=1.0
         azp(i)=1.0/3.0
         czp(i)=1.0/3.0
        enddo


       ars=7./9.
       brs=1./36.
  !     ars1=-2.5
  !     brs1=2.
  !     crs1=0.5
       ars1=-17./6.
       brs1=3./2.
       crs1=3./2.
       drs1=-1./6.
       ars2=0.75


       return
       end subroutine inideriv

!-----------------------------------------------------------------------

        SUBROUTINE derivper(n,u,a,b,c,delta)

        use consrs
        IMPLICIT NONE
        integer :: n
        integer i,j,k,l,m
        real :: delta,bet,alpha,beta,gamma,fact
        real,dimension (n) ::  u,r,gam,a,b,c,bb
        real, dimension (n) :: x,z


          DO l=3,n-2
           r(l)=delta*(ars*(u(l+1)-u(l-1))+    &
                     brs*(u(l+2)-u(l-2)))
          ENDDO
           r(1)=delta*(ars*(u(2)-u(n))+        &
                     brs*(u(3)-u(n-1)))
           r(n)=delta*(ars*(u(1)-u(n-1))+      &
                     brs*(u(2)-u(n-2)))
           r(2)=delta*(ars*(u(3)-u(1))+        &
                     brs*(u(4)-u(n)))
           r(n-1)=delta*(ars*(u(n)-u(n-2))+    &
                     brs*(u(1)-u(n-3)))



          alpha=a(1)
          beta=c(n)
          gamma=-b(1)
          bb(1)=b(1)-gamma
          bb(n)=b(n)-alpha*beta/gamma
          do l=2,n-1
           bb(l)=b(l)
          enddo


          bet=bb(1)
          x(1)=r(1)/bet


          do l=2,n
           gam(l)=c(l-1)/bet
           bet=bb(l)-a(l)*gam(l)
           x(l)=(r(l)-a(l)*x(l-1))/bet
          enddo

          do l=n-1,1,-1
           x(l)=x(l)-gam(l+1)*x(l+1)
          enddo


          u(1)=gamma
          u(n)=alpha
          do l=2,n-1
           u(l)=0.0
          enddo

          bet=bb(1)
          z(1)=u(1)/bet

          do l=2,n
           gam(l)=c(l-1)/bet
           bet=bb(l)-a(l)*gam(l)
           z(l)=(u(l)-a(l)*z(l-1))/bet
          enddo

          do l=n-1,1,-1
           z(l)=z(l)-gam(l+1)*z(l+1)
          enddo

          fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)

          do l=1,n
           u(l)=x(l)-fact*z(l)
          enddo

          return
          end subroutine derivper
! ------------------------------------------------------------------------------

SUBROUTINE derivnper(n,u,a,b,c,delta)

use consrs
IMPLICIT NONE
integer i,j,k,l,m
integer :: n
real :: delta,bet
real, dimension (n) ::  u,r,gam,a,b,c


  DO l=3,n-2
   r(l)=delta*(ars*(u(l+1)-u(l-1))+     &
             brs*(u(l+2)-u(l-2)))
  ENDDO
   r(1)=delta*(ars1*u(1)+brs1*u(2)+crs1*u(3)+drs1*u(4))
   r(n)=delta*(-ars1*u(n)-brs1*u(n-1)-crs1*u(n-2)-drs1*u(n-3))
   r(2)=delta*(ars2*(u(3)-u(1)))
   r(n-1)=delta*(ars2*(u(n)-u(n-2)))

  bet=b(1)
  u(1)=r(1)/bet

  do l=2,n
   gam(l)=c(l-1)/bet
   bet=b(l)-a(l)*gam(l)
   u(l)=(r(l)-a(l)*u(l-1))/bet
  enddo

  do l=n-1,1,-1
   u(l)=u(l)-gam(l+1)*u(l+1)
  enddo

  return
  end subroutine derivnper

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module constants

  real, parameter :: Re = 5000
  real, parameter :: Pr = 0.7
  real, parameter :: M = 0.3
  real, parameter :: gam = 1.4

end module constants

!-------------------------------------------------------------------------------

module constdim

  integer :: nx = 100
  integer :: ny = 100
  integer :: nz = 100
  integer :: ne = 5     ! numero de ecuaciones

  real, parameter :: lx=3.0
  real, parameter :: ly=3.0
  real, parameter :: lz=3.0

end module constdim

!-------------------------------------------------------------------------------

subroutine u2prim(U,P) ! entra U, regresa P

use constdim   ! nx, ny, nz, ne, lx, ly, lz
use constants  ! Re, Pr, M, gamma

implicit none

real, dimension (ne,nx,ny,nz) :: U,P
integer :: i,j,k

! allocate (U(ne,nx,ny,nz))
! allocate (P(ne,nx,ny,nz))

do k=1,nz
  do j=1,ny
    do i=1,nx
      P(1,i,j,k)=U(1,i,j,k)
      P(2,i,j,k)=U(2,i,j,k)/U(1,i,j,k)
      P(3,i,j,k)=U(3,i,j,k)/U(1,i,j,k)
      P(4,i,j,k)=U(4,i,j,k)/U(1,i,j,k)
      P(5,i,j,k)=(U(5,i,j,k)-gam*M**2/2*(U(2,i,j,k)**2/U(1,i,j,k)+U(3,i,j,k)**2/U(1,i,j,k)+U(4,i,j,k)**2/U(1,i,j,k)))*(gam-1)
    end do
  end do
end do


return

end subroutine u2prim

!-------------------------------------------------------------------------------

subroutine prim2u(P,U) ! entra P, regresa U

use constdim   ! nx, ny, nz, ne, lx, ly, lz
use constants  ! Re, Pr, M, gamma

implicit none

real, dimension (ne,nx,ny,nz) :: U,P
integer :: i,j,k

! allocate (P(ne,nx,ny,nz))
! allocate (U(ne,nx,ny,nz))

do k=1,nz
  do j=1,ny
    do i=1,nx
      U(1,i,j,k)=P(1,i,j,k)
      U(2,i,j,k)=P(1,i,j,k)*P(2,i,j,k)
      U(3,i,j,k)=P(1,i,j,k)*P(3,i,j,k)
      U(4,i,j,k)=P(1,i,j,k)*P(4,i,j,k)
      U(5,i,j,k)=P(5,i,j,k)/(gam-1)+gam*M**2*P(1,i,j,k)/2*(P(2,i,j,k)**2+P(3,i,j,k)**2+P(4,i,j,k)**2)
    end do
  end do
end do

return
end subroutine prim2u

!-------------------------------------------------------------------------------

subroutine derivQx(Q,dQx) ! entra Q (una funcion) y sale su derivada dQx

use constdim   ! nx, ny, nz, ne, lx, ly, lz
use constants  ! Re, Pr, M, gamma

implicit none

real :: dx
real, dimension (nx,ny,nz) :: Q,dQx
integer :: i,j,k

dx=lx/float(nx)

do k=1,nz
  do j=1,ny
    do i=2,nx-1
      dQx(i,j,k) = (Q(i+1,j,k)-Q(i-1,j,k))/(2.0*dx)
    end do
  end do
end do

do k=1,nz
  do j=1,ny
    dQx(1,j,k) = (Q(2,j,k)-Q(nx,j,k))/(2.0*dx)
    dQx(nx,j,k) = (Q(1,j,k)-Q(nx-1,j,k))/(2.0*dx)
  end do
end do

return
end subroutine derivQx

!-------------------------------------------------------------------------------

subroutine derivQy(Q,dQy) ! entra Q (una funcion) y sale su derivada dQy

use constdim   ! nx, ny, nz, ne, lx, ly, lz
use constants  ! Re, Pr, M, gamma

implicit none

real :: dy
real, dimension (nx,ny,nz) :: Q,dQy
integer :: i,j,k

dy=ly/float(ny)

do k=1,nz
  do i=1,nx
    do j=2,ny-1
      dQy(i,j,k) = (Q(i,j+1,k)-Q(i,j-1,k))/(2.0*dy)
    end do
  end do
end do

do k=1,nz
  do i=1,nx
    dQy(i,1,k) = (Q(i,2,k)-Q(i,ny,k))/(2.0*dy)
    dQy(i,ny,k) = (Q(i,1,k)-Q(i,ny-1,k))/(2.0*dy)
  end do
end do

return
end subroutine derivQy

!-------------------------------------------------------------------------------

subroutine derivQz(Q,dQz) ! entra Q (una funcion) y sale su derivada dQz

use constdim   ! nx, ny, nz, ne, lx, ly, lz
use constants  ! Re, Pr, M, gamma

implicit none

real :: dz
real, dimension (nx,ny,nz) :: Q,dQz
integer :: i,j,k

dz=lz/float(nz)

do j=1,ny
  do i=1,nx
    do k=2,nz-1
      dQz(i,j,k) = (Q(i,j,k+1)-Q(i,j,k-1))/(2.0*dz)
    end do
  end do
end do

do j=1,ny
  do i=1,nx
    dQz(i,j,1) = (Q(i,j,2)-Q(i,j,nz))/(2.0*dz)
    dQz(i,j,nz) = (Q(i,j,1)-Q(i,j,nz-1))/(2.0*dz)
  end do
end do

return

end subroutine derivQz

!-------------------------------------------------------------------------------

subroutine Sij()

use constdim   ! nx, ny, nz, ne, lx, ly, lz
use constants  ! Re, Pr, M, gamma

implicit none

!real, dimension (nx,ny,nz) :: S1j
integer :: i,j,k
real, dimension (nx,ny,nz) :: dux,duy,duz
real, dimension (nx,ny,nz) :: dvx,dvy,dvz
real, dimension (nx,ny,nz) :: dwx,dwy,dwz
real, dimension (nx,ny,nz) :: S11,S12,S13,S21,S22,S23,S31,S32,S33

call derivQx(u,dux)
call derivQy(u,duy)
call derivQz(u,duz)
call derivQx(v,dvx)
call derivQy(v,dvy)
call derivQz(v,dvz)
call derivQx(w,dwx)
call derivQy(w,dwy)
call derivQz(w,dwz)

do k=1,nz
  do j=1,ny
    do i=1,nx
      S11(i,j,k)=dux(i,j,k)+dux(i,j,k)-2/3*(dux(i,j,k)+dvy(i,j,k)+dwz(i,j,k))
      S12(i,j,k)=duy(i,j,k)+dvx(i,j,k)
      S13(i,j,k)=duz(i,j,k)+dwx(i,j,k)
      S21(i,j,k)=dvx(i,j,k)+duy(i,j,k)
      S22(i,j,k)=dvy(i,j,k)+dvy(i,j,k)-2/3*(dux(i,j,k)+dvy(i,j,k)+dwz(i,j,k))
      S23(i,j,k)=dvz(i,j,k)+dwy(i,j,k)
      S31(i,j,k)=dwx+duz(i,j,k)
      S32(i,j,k)=dwy(i,j,k)+dvz(i,j,k)
      S33(i,j,k)=dwz(i,j,k)+dwz(i,j,k)-2/3*(dux(i,j,k)+dvy(i,j,k)+dwz(i,j,k))
    end do
  end do
end do

end subroutine Sij

!-------------------------------------------------------------------------------

subroutine fluxF(P,F,T) ! damos P, sale F

use constdim   ! nx, ny, nz, ne, lx, ly, lz
use constants  ! Re, Pr, M, gamma

implicit none

real, allocatable, dimension (nx,ny,nz) :: P,F
integer :: i,j,k
real, allocatable, dimension (nx,ny,nz) :: T,dTx

! allocate (F(ne,nx,ny,nz))

call derivQx(T,dTx)

do k=1,nz
  do j=1,ny
    do i=1,nx
      F(1,i,j,k)=-P(1,i,j,k)*P(2,i,j,k)
      F(2,i,j,k)=-P(1,i,j,k)*P(2,i,j,k)-1/(gam*M**2)*P(5,i,j,k)+1/Re*S11
      F(3,i,j,k)=-P(1,i,j,k)*P(2,i,j,k)*P(3,i,j,k)+1/Re*S12
      F(4,i,j,k)=-P(1,i,j,k)*P(2,i,j,k)*P(4,i,j,k)+1/Re*S13
      F(5,i,j,k)=-P(2,i,j,k)*(P(5,i,j,k)/(gam-1) + &
      gam*M**2*P(5,i,j,k)/2*(P(2,i,j,k)**2+P(3,i,j,k)**2+P(4,i,j,k)**2)+P(5,i,j,k)) + &
      gam*M**2/Re*(S11*P(2,i,j,k)+S12*P(3,i,j,k)+S13*P(4,i,j,k)) + &
      gam/(gam-1)*(1/(Pr*Re))*dTx
    end do
  end do
end do

end subroutine fluxF

!-------------------------------------------------------------------------------

program derivada

use dimensiones ! si se usa, para nx,ny,nz
use deltas      ! si se usa, para deltax,deltay,deltaz
use consderper
use derivtools

implicit none

real :: lx,ly,lz
real :: dx,dy,dz
real :: x,y,z
integer :: i,j,k
real, allocatable, dimension (:) :: dux,duy,duz
real, allocatable, dimension (:) :: dvx,dvy,dvz
real, allocatable, dimension (:) :: dwx,dwy,dwz

nx=100
ny=100
nz=100

lx=1.0
ly=1.0
lz=1.0

dx=lx/nx
dy=ly/ny
dz=lz/nz

deltax=1.0/dx
deltay=1.0/dy
deltaz=1.0/dz



call deriv_alloc() ! aqui se declara u, y le da dimension a u, u(nx)
call inideriv() ! aqui definen parametros que se usan en la deriv periodica


do k=1,nz
  z=dz*float(k-1)
  do j=1,ny
    y=dy*float(j-1)
    do i=1,nx
      x=dx*float(i-1)
      u(i)=x*x*x
      dux(i)=u(i)
    end do
  end do
end do



call derivper(nx,dux,axp,bxp,cxp,deltax)  ! derivada, regresa dux (derivada de u, respecto a x)



open(78,file = 'deriv.txt', form = 'formatted')
do i=1,nx
  write(78,105)dx*float(i),u(i),dux(i)
end do
close(78)

105 FORMAT(5(f12.6))



end program derivada
