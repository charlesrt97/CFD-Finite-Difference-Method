! Simulation of a 3D compressible turbulent jet flow
! solves the 3D compressible Navier-Stokes equations
! using a compact finite difference scheme, a turbulence model and a 3rd order Runge-Kutta scheme

! For a higher convergence, a refined mesh or an unstructured mesh must be used
!------------------------------------------------------------------------------

module variables

integer :: nx,ny,nz,nt
integer :: i,j,k,m
integer :: imp,condimp
real :: u0,v0,w0,p0,t0,rdm
real :: dx,dy,dz,dt
real :: x1,y1,z1,te
real :: re,pr,ma,gm ! parameters
character(20) :: readdata,savedata
character(3) :: count

real,allocatable,dimension(:) :: x,y,z

real,allocatable,dimension(:,:,:) :: u,v,w,T,P
real,allocatable,dimension(:,:,:) :: ro,ru,rv,rw,rn,ec
real,allocatable,dimension(:,:,:) :: turb_visc

real,allocatable,dimension(:,:,:) :: dudx,dudy,dudz
real,allocatable,dimension(:,:,:) :: dvdx,dvdy,dvdz
real,allocatable,dimension(:,:,:) :: dwdx,dwdy,dwdz
real,allocatable,dimension(:,:,:) :: dTdx,dTdy,dTdz

real,allocatable,dimension(:,:,:) :: Fmass,Gmass,Hmass
real,allocatable,dimension(:,:,:) :: Fmomx,Gmomx,Hmomx
real,allocatable,dimension(:,:,:) :: Fmomy,Gmomy,Hmomy
real,allocatable,dimension(:,:,:) :: Fmomz,Gmomz,Hmomz
real,allocatable,dimension(:,:,:) :: Fnrg,Gnrg,Hnrg

real,allocatable,dimension(:,:,:) :: Unmass,Unmomx,Unmomy,Unmomz,Unnrg
real,allocatable,dimension(:,:,:) :: U1mass,U1momx,U1momy,U1momz,U1nrg
real,allocatable,dimension(:,:,:) :: U2mass,U2momx,U2momy,U2momz,U2nrg
real,allocatable,dimension(:,:,:) :: Unnmass,Unnmomx,Unnmomy,Unnmomz,Unnnrg
real,allocatable,dimension(:,:,:) :: divmass,divmomx,divmomy,divmomz,divnrg

end module variables

!------------------------------------------------------------------------------

subroutine vel_bounds()

use variables
implicit none
real :: u_upper,u_lower,v_upper,v_lower,w_upper,w_lower

u_upper = 1.3
u_lower = -u_upper
v_upper = 0.1
v_lower = -v_upper
w_upper = 0.1
w_lower = -w_upper

do k=1,nz
  do j=1,ny
    do i=1,nx
      if (u(i,j,k).ge.u_upper) then
        u(i,j,k)=u_upper
      end if
      if (u(i,j,k).le.u_lower) then
        u(i,j,k)=u_lower
      end if

      if (v(i,j,k).ge.v_upper) then
        v(i,j,k)=v_upper
      end if
      if (v(i,j,k).le.v_lower) then
        v(i,j,k)=v_lower
      end if

      if (w(i,j,k).ge.w_upper) then
        w(i,j,k)=w_upper
      end if
      if (w(i,j,k).le.w_lower) then
        w(i,j,k)=w_lower
      end if

    end do
  end do
end do

return
end subroutine vel_bounds

!------------------------------------------------------------------------------

subroutine press_bounds()

use variables
implicit none
real :: T_upper,T_lower,P_upper,P_lower

T_upper = 1.1
T_lower = 0.9
P_upper = 1.1
P_lower = 0.9

do k=1,nz
  do j=1,ny
    do i=1,nx

      if (T(i,j,k).ge.T_upper) then
        T(i,j,k)=T_upper
      end if
      if (T(i,j,k).le.T_lower) then
        T(i,j,k)=T_lower
      end if

      if (P(i,j,k).ge.P_upper) then
        P(i,j,k)=P_upper
      end if
      if (P(i,j,k).le.P_lower) then
        P(i,j,k)=P_lower
      end if
    end do
  end do
end do

return
end subroutine press_bounds

!------------------------------------------------------------------------------

subroutine conservatives()

use variables
implicit none

ro=P/T
ru=ro*u
rv=ro*v
rw=ro*w
rn=ro*(T/(gm-1)+gm*(ma**2)*(u**2+v**2+w**2)/2.)

return
end subroutine conservatives

!------------------------------------------------------------------------------

subroutine primitives(Umass,Umomx,Umomy,Umomz,Unrg)

use variables
implicit none
real,dimension(nx,ny,nz) :: Umass,Umomx,Umomy,Umomz,Unrg

ro=Umass

u=Umomx/ro
v=Umomy/ro
w=Umomz/ro
T=(Unrg/ro-gm*(ma**2)*(u**2+v**2+w**2)/2.)*(gm-1)
P=ro*T

call press_bounds()
call vel_bounds()

return
end subroutine primitives

!------------------------------------------------------------------------------

subroutine fluxes()

use variables
implicit none
real,dimension(nx,ny,nz) :: duv,duw,dvw,div

duv=dvdx+dudy
duw=dwdx+dudz
div=dudx+dvdy+dwdz

! conserv. of mass
Fmass=-ru
Gmass=-rv
Hmass=-rw

! x-momentum
Fmomx = -(ru*u+(P/(gm*ma**2))-(2./re)*(dudx-(1./3.)*div)-turb_visc*2.*(dudx-(1./3.)*div))
Gmomx = -(rv*u-(1./re)*duv-turb_visc*duv)
Hmomx = -(rw*u-(1./re)*duw-turb_visc*duw)

! y-momentum
Fmomy = -(ru*v-(1./re)*duv-turb_visc*duv)
Gmomy = -(rv*v+(P/(gm*ma**2))-(2./re)*(dvdy-(1./3.)*div)-turb_visc*2.*(dvdy-(1./3.)*div))
Hmomy = -(rw*v-(1./re)*dvw-turb_visc*dvw)

! z-momentum
Fmomz = -(ru*w-(1./re)*duw-turb_visc*duw)
Gmomz = -(rv*w-(1./re)*dvw-turb_visc*dvw)
Hmomz = -(rw*w+(P/(gm*ma**2))-(2./re)*(dwdz-(1./3.)*div)-turb_visc*2.*(dwdz-(1./3.)*div))

! energy eq.
Fnrg = -((rn+P)*u-(gm*ma**2./re)*                  &
       (2.*u*(dudx-(1./3.)*div)+v*duv+w*duw)       &
       -(gm/(pr*re*(gm-1)))*dTdx)
Gnrg = -((rn+P)*v-(gm*ma**2./re)*                  &
       (u*duv+2.*v*(dvdy-(1./3.)*div)+w*dvw)       &
       -(gm/(pr*re*(gm-1)))*dTdy)
Hnrg = -((rn+P)*w-(gm*ma**2./re)*                  &
       (u*duw+v*dvw+2.*w*(dwdz-(1./3.)*div))       &
       -(gm/(pr*re*(gm-1)))*dTdz)

call periodic_diverg(Fmass,Gmass,Hmass)
call periodic_diverg(Fmomx,Gmomx,Hmomx)
call periodic_diverg(Fmomy,Gmomy,Hmomy)
call periodic_diverg(Fmomz,Gmomz,Hmomz)
call periodic_diverg(Fnrg,Gnrg,Hnrg)

divmass=Fmass+Gmass+Hmass
divmomx=Fmomx+Gmomx+Hmomx
divmomy=Fmomy+Gmomy+Hmomy
divmomz=Fmomz+Gmomz+Hmomz
divnrg=Fnrg+Gnrg+Hnrg

return
end subroutine fluxes

!------------------------------------------------------------------------------

subroutine timestep()

use variables
implicit none
real :: u_max,v_max,w_max,c_max, CFL, lambda
real,dimension (nx,ny,nz) :: c

CFL=0.5

c=1./ma*(T**(1./2.))
u_max=maXVAL(abs(u)+c)
v_max=maXVAL(abs(v)+c)
w_max=maXVAL(abs(w)+c)
lambda=MIN(dx/u_max,dy/v_max,dz/w_max)
dt=lambda*CFL

return
end subroutine timestep

!------------------------------------------------------------------------------

subroutine turb_viscosity()

use variables
implicit none
real,dimension(nx,ny,nz,6) :: v1,v2,v3

do k=2,nz-1
  do j=2,ny-1
    do i=2,nx-1
      v1(i,j,k,1) = (u(i,j,k)-u(i+1,j,k))*dx
      v1(i,j,k,2) = (u(i,j,k)-u(i-1,j,k))*dx
      v1(i,j,k,3) = (u(i,j,k)-u(i,j+1,k))*dy
      v1(i,j,k,4) = (u(i,j,k)-u(i,j-1,k))*dy
      v1(i,j,k,5) = (u(i,j,k)-u(i,j,k+1))*dz
      v1(i,j,k,6) = (u(i,j,k)-u(i,j,k-1))*dz

      v2(i,j,k,1) = (v(i,j,k)-v(i+1,j,k))*dx
      v2(i,j,k,2) = (v(i,j,k)-v(i-1,j,k))*dx
      v2(i,j,k,3) = (v(i,j,k)-v(i,j+1,k))*dy
      v2(i,j,k,4) = (v(i,j,k)-v(i,j-1,k))*dy
      v2(i,j,k,5) = (v(i,j,k)-v(i,j,k+1))*dz
      v2(i,j,k,6) = (v(i,j,k)-v(i,j,k-1))*dz

      v3(i,j,k,1) = (w(i,j,k)-w(i+1,j,k))*dx
      v3(i,j,k,2) = (w(i,j,k)-w(i-1,j,k))*dx
      v3(i,j,k,3) = (w(i,j,k)-w(i,j+1,k))*dy
      v3(i,j,k,4) = (w(i,j,k)-w(i,j-1,k))*dy
      v3(i,j,k,5) = (w(i,j,k)-w(i,j,k+1))*dz
      v3(i,j,k,6) = (w(i,j,k)-w(i,j,k-1))*dz
    end do
  end do
end do

v1(1,:,:,:) = v1(2,:,:,:)
v1(:,1,:,:) = v1(:,2,:,:)
v1(:,:,1,:) = v1(:,:,2,:)
v1(nx,:,:,:) = v1(nx-1,:,:,:)
v1(:,ny,:,:) = v1(:,ny-1,:,:)
v1(:,:,nz,:) = v1(:,:,nz-1,:)

v2(1,:,:,:) = v2(2,:,:,:)
v2(:,1,:,:) = v2(:,2,:,:)
v2(:,:,1,:) = v2(:,:,2,:)
v2(nx,:,:,:) = v2(nx-1,:,:,:)
v2(:,ny,:,:) = v2(:,ny-1,:,:)
v2(:,:,nz,:) = v2(:,:,nz-1,:)

v3(1,:,:,:) = v3(2,:,:,:)
v3(:,1,:,:) = v3(:,2,:,:)
v3(:,:,1,:) = v3(:,:,2,:)
v3(nx,:,:,:) = v3(nx-1,:,:,:)
v3(:,ny,:,:) = v3(:,ny-1,:,:)
v3(:,:,nz,:) = v3(:,:,nz-1,:)

do k=1,nz
  do j=1,ny
    do i=1,nx
      turb_visc(i,j,k) = v1(i,j,k,1)**2+v2(i,j,k,1)**2+v3(i,j,k,1)**2  &
                       +v1(i,j,k,2)**2+v2(i,j,k,2)**2+v3(i,j,k,2)**2   &
                       +v1(i,j,k,3)**2+v2(i,j,k,3)**2+v3(i,j,k,3)**2   &
                       +v1(i,j,k,4)**2+v2(i,j,k,4)**2+v3(i,j,k,4)**2   &
                       +v1(i,j,k,5)**2+v2(i,j,k,5)**2+v3(i,j,k,5)**2   &
                       +v1(i,j,k,6)**2+v2(i,j,k,6)**2+v3(i,j,k,6)**2
    end do
  end do
end do

turb_visc = turb_visc/6.0
turb_visc = 0.063*((dx*dy*dz)**(1./3.))*(turb_visc**(1./2.))

return
end subroutine turb_viscosity

!------------------------------------------------------------------------------

subroutine VT_deriv()

use variables
implicit none

dudx=u
dudy=u
dudz=u
dvdx=v
dvdy=v
dvdz=v
dwdz=w
dwdy=w
dwdz=w
call periodic_diverg(dudx,dudy,dudz)
call periodic_diverg(dvdx,dvdy,dvdz)
call periodic_diverg(dwdx,dwdy,dwdz)

dTdx=T
dTdy=T
dTdz=T

call periodic_diverg(dTdx,dTdy,dTdz)
call turb_viscosity()

return
end subroutine VT_deriv

!------------------------------------------------------------------------------

subroutine finite_deriv(n,derivf,delta)

implicit none
integer :: n,i
real :: delta
real, dimension(n) :: derivf,y

do i=2, n-1
  derivf(i) = (y(i+1)-y(i-1))/(2*delta)
end do
derivf(1) = (y(2)-y(n))/(2*delta)
derivf(n) = (y(1)-y(n-1))/(2*delta)

derivf=y

return
end subroutine finite_deriv

!------------------------------------------------------------------------------

subroutine periodic_deriv(n,u,delta)

implicit none
integer :: n,i
real :: alpha,beta,gamma,a,b,delta,fact
real, dimension(n) :: u,y
real, dimension(n) :: f,gam,x,z

alpha=1./3.
beta=1./3.
gamma=-1.
a=14./9.*1./2.
b=1./9.*1./4.

do i=3,n-2
  y(i) = (a*(u(i+1)-u(i-1))+    &
          b*(u(i+2)-u(i-2)))/delta
end do

y(1) = (a*(u(2)-u(n))+          &
        b*(u(3)-u(n-1)))/delta
y(n) = (a*(u(1)-u(n-1))+        &
        b*(u(2)-u(n-2)))/delta
y(2) = (a*(u(3)-u(1))+          &
        b*(u(4)-u(n)))/delta
y(n-1) = (a*(u(n)-u(n-2))+      &
        b*(u(1)-u(n-3)))/delta

f(1) = 2.
f(n) = 1.-alpha*beta/gamma

do i=2,n-1
  f(i)=1.0
end do

x(1) = y(1)/f(1)

do i=2,n
  gam(i)=(alpha/f(1))
  f(1)=(f(i)-alpha*gam(i))
  x(i)=(y(i)-alpha*x(i-1))/f(1)
end do

do i=n-1,1,-1
  x(i)=x(i)-gam(i+1)*x(i+1)
end do

u(1)=gamma
u(n)=alpha
do i=2,n-1
 u(i)=0.0
enddo

f(1)=2
z(1)=u(1)/f(1)

do i=2,n
  gam(i) = (alpha/f(1))
  f(1)=(f(i)-alpha*gam(i))
  z(i)=(u(i)-alpha*z(i-1))/f(1)
end do

do i=n-1,1,-1
  z(i)=z(i)-gam(i+1)*z(i+1)
end do

fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)

do i=1,n
  u(i)=x(i)-z(i)*fact
end do

return
end subroutine periodic_deriv

!------------------------------------------------------------------------------

subroutine periodic_diverg(F,G,H)
use variables
implicit none
real, dimension(nx,ny,nz) :: F,G,H

do k=1,nz
  do j=1,ny
    call periodic_deriv(nx,F(:,j,k),dx)
  end do
end do

do k=1,nz
  do i=1,nx
    call periodic_deriv(ny,G(i,:,k),dy)
  end do
end do

do j=1,ny
  do i=1,nx
    call periodic_deriv(nz,H(i,j,:),dz)
  end do
end do

return
end subroutine periodic_diverg

!------------------------------------------------------------------------------

subroutine finite_diverg(F,G,H)

use variables
implicit none
real, dimension(nx,ny,nz) :: F,G,H

do k=1,nz
  do j=1,ny
    call finite_deriv(nx,F(:,j,k),dx)
  end do
end do

do k=1,nz
  do i=1,nx
    call finite_deriv(ny,G(i,:,k),dy)
  end do
end do

do j=1,ny
  do i=1,nx
    call finite_deriv(nz,H(i,j,:),dz)
  end do
end do

return
end subroutine finite_diverg

!------------------------------------------------------------------------------
! only saves the x-component of the velocity field
subroutine save_data()
use variables
implicit none

write(count,'(i6)')m
savedata = 'jet_'//count//'.vtk'

open(78, file=savedata, form='formatted')

write(78,110) '# vtk DataFile Version 2.3'
write(78,110) '3D Mesh'
write(78,111) 'ASCII'

write(78,110) 'DATASET STRUCTUreD_GRID'
write(78,120) 'DIMENSIONS', nx, ny, nz
write(78,130) 'POINTS', nx*ny*nz,' float'

do k=1,nz
  do j=1,ny
    do i=1,nx
      write(78,100) x(i),y(j),z(k)
    end do
  end do
end do

write(78,140) 'POINT_DATA',nx*ny*nz

write(78,110) 'SCALARS VEL_U float'
write(78,110) 'LOOKUP_TABLE default'
do k=1,nz
  do j=1,ny
    do i=1,nx
      write(78,100) u(i,j,k)
    end do
  end do
end do

close(78)

100 format(5(f12.6))
110 format(A)
111 format(A,/)
120 format(A,I4,I4,I4)
130 format(A,I10,A)
140 format(A,I10)

return
end subroutine save_data

!------------------------------------------------------------------------------

program main

use variables
implicit none

open(25,file='data.input',form='unformatted')
read(25)re
read(25)pr
read(25)ma
read(25)gm

read(25)nx
read(25)ny
read(25)nz

read(25)u0
read(25)v0
read(25)w0
read(25)p0
read(25)t0
read(25)x1
read(25)y1
read(25)z1

close(25)

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

allocate(turb_visc(nx,ny,nz))

allocate(ro(nx,ny,nz))
allocate(ru(nx,ny,nz))
allocate(rv(nx,ny,nz))
allocate(rw(nx,ny,nz))
allocate(rn(nx,ny,nz))
allocate(ec(nx,ny,nz))

allocate(u(nx,ny, nz))
allocate(v(nx,ny, nz))
allocate(w(nx,ny, nz))
allocate(P(nx,ny, nz))
allocate(T(nx,ny, nz))

allocate(dudx(nx,ny,nz))
allocate(dudy(nx,ny,nz))
allocate(dudz(nx,ny,nz))
allocate(dvdx(nx,ny,nz))
allocate(dvdy(nx,ny,nz))
allocate(dvdz(nx,ny,nz))
allocate(dwdx(nx,ny,nz))
allocate(dwdy(nx,ny,nz))
allocate(dwdz(nx,ny,nz))
allocate(dTdx(nx,ny,nz))
allocate(dTdy(nx,ny,nz))
allocate(dTdz(nx,ny,nz))

allocate(Fmass(nx,ny,nz))
allocate(Gmass(nx,ny,nz))
allocate(Hmass(nx,ny,nz))
allocate(Fmomx(nx,ny,nz))
allocate(Gmomx(nx,ny,nz))
allocate(Hmomx(nx,ny,nz))
allocate(Fmomy(nx,ny,nz))
allocate(Gmomy(nx,ny,nz))
allocate(Hmomy(nx,ny,nz))
allocate(Fmomz(nx,ny,nz))
allocate(Gmomz(nx,ny,nz))
allocate(Hmomz(nx,ny,nz))
allocate(Fnrg(nx,ny,nz))
allocate(Gnrg(nx,ny,nz))
allocate(Hnrg(nx,ny,nz))

allocate(Unmass(nx,ny,nz))
allocate(Unmomx(nx,ny,nz))
allocate(Unmomy(nx,ny,nz))
allocate(Unmomz(nx,ny,nz))
allocate(Unnrg(nx,ny,nz))

allocate(U1mass(nx,ny,nz))
allocate(U1momx(nx,ny,nz))
allocate(U1momy(nx,ny,nz))
allocate(U1momz(nx,ny,nz))
allocate(U1nrg(nx,ny,nz))
allocate(U2mass(nx,ny,nz))
allocate(U2momx(nx,ny,nz))
allocate(U2momy(nx,ny,nz))
allocate(U2momz(nx,ny,nz))
allocate(U2nrg(nx,ny,nz))
allocate(Unnmass(nx,ny,nz))
allocate(Unnmomx(nx,ny,nz))
allocate(Unnmomy(nx,ny,nz))
allocate(Unnmomz(nx,ny,nz))
allocate(Unnnrg(nx,ny,nz))
allocate(divmass(nx,ny,nz))
allocate(divmomx(nx,ny,nz))
allocate(divmomy(nx,ny,nz))
allocate(divmomz(nx,ny,nz))
allocate(divnrg(nx,ny,nz))

dx = x1/float(nx-1)
dy = y1/float(ny-1)
dz = z1/float(nz-1)

do i=1,nx
  x(i) = dx*float(i-1)
end do
do j=1,ny
  y(j) = dy*float(j-1)
end do
do k=1,nz
  z(k) = dz*float(k-1)
end do

imp=0
m=1
nt=2000

do j=1,ny/2
  u(1,j,1) = u0*(erf((y(j)-0.2*y1)*7.))
  u(1,ny-(j-1),1) = u(1,j,1)
end do

do i=1,nx
  do k=1,nz
    u(i,:,k) = u(1,:,1)
  end do
end do

do k=1,nz
  do i=1,nx
    do j=1,ny
      rdm = j+k-i
      call random_number(rdm)
      u(i,j,k) = u(i,j,k)+0.05*u0*(rand()-0.5)
      v(i,j,k) = 0.01*v0*(rand()-0.5)
      w(i,j,k) = 0.01*w0*(rand()-0.5)
    end do
  end do
end do

P(:,:,:) = p0
T(:,:,:) = t0

call conservatives()
call VT_deriv()
call fluxes()
call timestep()

Unmass=ro
Unmomx=ru
Unmomy=rv
Unmomz=rw
Unnrg=rn

do j=1,ny
  write(*,*)y(j),u(20,j,50),u(3,j,5),dudy(3,j,5),w(9,j,4),dwdy(9,j,4)
end do

condimp = 100

do m=m,nt
  U1mass = Unmass+dt*divmass
  U1momx = Unmomx+dt*divmomx
  U1momy = Unmomy+dt*divmomy
  U1momz = Unmomz+dt*divmomz
  U1nrg = Unnrg+dt*divnrg

  call primitives(U1mass,U1momx,U1momy,U1momz,U1nrg)
  call conservatives()
  call VT_deriv()
  call fluxes()

  U2mass = 3./4.*Unmass+1./4.*U1mass+1./4.*dt*divmass
  U2momx = 3./4.*Unmomx+1./4.*U1momx+1./4.*dt*divmomx
  U2momy = 3./4.*Unmomy+1./4.*U1momy+1./4.*dt*divmomy
  U2momz = 3./4.*Unmomz+1./4.*U1momz+1./4.*dt*divmomz
  U2nrg = 3./4.*Unnrg+1./4.*U1nrg+1./4.*dt*divnrg

  call primitives(U2mass, U2momx, U2momy, U2momz, U2nrg)
  call conservatives()
  call VT_deriv()
  call fluxes()

  Unnmass = 1./3.*Unmass+2./3.*U2mass+2./3.*dt*divmass
  Unnmomx = 1./3.*Unmomx+2./3.*U2momx+2./3.*dt*divmomx
  Unnmomy = 1./3.*Unmomy+2./3.*U2momy+2./3.*dt*divmomy
  Unnmomz = 1./3.*Unmomz+2./3.*U2momz+2./3.*dt*divmomz
  Unnnrg = 1./3.*Unnrg+2./3.*U2nrg+2./3.*dt*divnrg

  call primitives(Unnmass, Unnmomx, Unnmomy, Unnmomz, Unnnrg)
  call conservatives()

  imp = imp+1

  if (imp.eq.condimp) then
    call save_data
    write(*,*)dt,m,Unmass(nx/2,ny/2,nz/2),Unmomx(nx/2,ny/2,nz/2),          &
              Unmomy(nx/2,ny/2,nz/2),Unmomz(nx/2,ny/2,nz/2),               &
              Unnrg(nx/2,ny/2,nz/2),P(nx/2,ny/2,nz/2),T(nx/2,ny/2,nz/2),   &
              turb_visc(nx/2,ny/2,nz/2)
    imp=0
  end if

  call VT_deriv()
  call fluxes()
  call timestep()

  Unmass = Unnmass
  Unmomx = Unnmomx
  Unmomy = Unnmomy
  Unmomz = Unnmomz
  Unnrg = Unnnrg

end do

200 format(7(e16.8))

end program main

!------------------------------------------------------------------------------
