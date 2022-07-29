! uses an alternating direction implicit (ADI) scheme to solve the 2D transient heat conduction eq.
! subjected to neumann and dirichlet boundary conditions
! uses a subroutine for inverting the tridiagonal matrix
! it as well computes the local and average nusselt numbers at the horizontal center of the domain

program heatcond2d

integer,parameter :: mi=150, nj=120
integer :: i_sample,it

real(kind=8) :: x(mi),y(nj),lx,ly
real(kind=8) :: prom

real(kind=8) :: T(mi,nj),Ai(mi),Bi(mi),Ci(mi),Ri(mi),Ui(mi)
real(kind=8) :: Aj(nj),Bj(nj),Cj(nj),Rj(nj),Uj(nj)

real (kind=8) :: dx,dy,dt
real (kind=8) :: kappa,T_initial,T_left,T_right,q_bot,q_top

real (kind=8) :: Tadim(mi,nj)
real (kind=8) :: Nux(nj)
real (kind=8) :: Nu

integer :: ii,jj,tt,time

character(16) outputfield
character(3) sample
i_sample=0
it=0

lx=1.d0
ly=1.d0

do ii=1,mi
  x(ii)=lx*dfloat(ii-1)/dfloat(mi-1)
end do

do jj=1,nj
  y(jj)=ly*dfloat(jj-1)/dfloat(nj-1)
end do

kappa=0.005d0
tt=2000
T_initial=20.d0
T_left=90.d0
T_right=5.d0
q_bot=0.d0  ! no heat flux (adiabatic walls)
q_top=0.d0  ! no heat flux (adiabatic walls)

T=T_initial

dt=0.005d0      ! time step
dx=lx/dfloat(mi-1)
dy=ly/dfloat(nj-1)

Ui=0
Uj=0

open(25,file='TempSEC.txt',form='formatted')
open(26,file='NuSEC.txt',form='formatted')

do time=1,tt
  d1=kappa*dt/(2.d0*dx*dx)
  d2=kappa*dt/(2.d0*dy*dy)
  ! x-sweep
  do jj=2,nj-1
    Bi(1)=1.d0   ! dirichlet BC
    Ci(1)=0.d0   ! dirichlet BC
    Ri(1)=T_left
    do ii=2,mi-1
      Ai(ii)=-d1
      Bi(ii)=1.d0+2.d0*d1
      Ci(ii)=-d1
      Ri(ii)=d2*T(ii,jj+1)+(1.d0-2.d0*d2)*T(ii,jj)+d2*T(ii,jj-1)
    end do
    Ai(mi)=0.d0
    Bi(mi)=1.d0
    Ri(mi)=T_right
    call tri(Ai,Bi,Ci,Ri,Ui,mi)
    T(:,jj)=Ui
  end do

  ! y-sweep
  do ii=2,mi-1
    Bj(1)=-1.d0 ! neumann BC
    Cj(1)=1.d0  ! neumann BC
    Rj(1)=q_bot
    do jj=2,nj-1
      Aj(jj)=-d2
      Bj(jj)=1.d0+2.d0*d2
      Cj(jj)=-d2
      Rj(jj)=d1*T(ii+1,jj)+(1.d0-2.d0*d1)*T(ii,jj)+d1*T(ii-1,jj)
    end do
    Aj(nj)=-1.d0 ! neumann BC
    Bj(nj)=1.d0  ! neumann BC
    Rj(nj)=q_top
    call tri(Aj,Bj,Cj,Rj,Uj,nj)
    T(ii,:)=Uj
  end do

  ! --- post-processing ---
  ! dimensionless temperature field
  do ii=1,mi
    do jj=2,nj-1
      Tadim(ii,jj)=(T(ii,jj)-T_right)/(T_left-T_right)
    end do
  end do

  ! local nusselt number at the horizontal center of the domain
  do jj=2,nj-1
    Nux(jj)=-(Tadim(mi/2+1,jj)-Tadim(mi/2,jj))/dx*lx
  end do

  ! average nusselt number at the horizontal center of the domain
  Nu=0
  do jj=2,nj-1
    Nu=Nux(jj)*dy+Nu
  end do

  write(25,160) dt*time,Tadim(10,nj/2),Tadim(20,nj/2),Tadim(30,nj/2), &
                 & Tadim(45,nj/2),Tadim(60,nj/2),Tadim(75,nj/2), &
                 & Tadim(90,nj/2),Tadim(100,nj/2),Tadim(115,nj/2), Tadim(135,nj/2), &
                 & Nu

  ! creates vtk files for post-processing
  if (mod(time,50)==0) then
    i_sample=i_sample+1
    write(sample,'(i3.3)')(i_sample)
    outputfield='temp_'//sample//'.vtk'

    open(78, file = Outputfield, form = 'formatted')

      write(78,110) '# vtk DataFile Version 2.3'
      write(78,110) '3D Mesh'
      write(78,111) 'ASCII'

      write(78,110) 'DATASET STRUCTURED_GRID'
      write(78,120) 'DIMENSIONS',mi, nj, 1
      write(78,130) 'POINTS', mi*nj, ' float'

      do jj=1,nj
        do ii=1,mi
          write(78,100) x(ii), y(jj), 0.d0
        enddo
      enddo

      write(78,140) 'POINT_DATA', mi*nj

      write(78,110) 'SCALARS Temperature float'
      write(78,110) 'LOOKUP_TABLE default'

        do jj=1,nj
          do ii=1,mi
            write(78,100) T(ii,jj)
          end do
        end do

    close(78)
  end if

end do
close(25)
close(26)

100  FORMAT(3(f12.6));
110  FORMAT(A);
111  FORMAT(A,/);
120  FORMAT(A,I4,I4,I4);
130  FORMAT(A,I10,A);
140  FORMAT(A,I10);

160 FORMAT(f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8)


end program heatcond2d
