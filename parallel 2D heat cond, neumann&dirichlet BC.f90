! uses an ADI scheme to solve the 2d transient heat conduction eq.
! subjected to neumann and dirichlet boundary conditions

! uses a domain decomposition method (neumann/dirichlet BC at the interface)
! and openMP to parallelize loops and sections

program heatsec

implicit none

! integer, parameter :: mi1=200, mi2=201, nj=200 ! mi=50, nj=20
! integer, parameter :: mi=400

integer, parameter :: mi1=75, mi2=76, nj=100 ! mi=50, nj=20
integer, parameter :: mi=150

integer :: ii, jj, tt, time, i_sample, coupling, iter

real(kind=8) :: x(mi),x1(mi1),x2(mi2),y(nj),lx1,lx2,lx,ly
real(kind=8) :: d1,d2
real(kind=8) :: Nu

real(kind=8) :: Ai1(mi1-1), Bi1(mi1), Ci1(mi1-1), Ri1(mi1)
real(kind=8) :: Ai2(mi2-1), Bi2(mi2), Ci2(mi2-1), Ri2(mi2)

real(kind=8) :: Aj1(nj-1), Bj1(nj), Cj1(nj-1), Rj1(nj)
real(kind=8) :: Aj2(nj-1), Bj2(nj), Cj2(nj-1), Rj2(nj)

real(kind=8) :: Nux(nj)

real(kind=8) :: T1(mi1,nj), T2(mi2,nj)
real(kind=8) :: T1n12(mi1,nj), T2n12(mi2,nj)
real(kind=8) :: Tast1(mi1,nj), Tast2(mi2,nj)
real(kind=8) :: Tbuena1(mi1,nj), Tbuena2(mi2,nj)
real(kind=8) :: Ttotal(mi,nj)
real(kind=8) :: Tadim(mi,nj)

real(kind=8) :: dx, dy, dt

real(kind=8) :: kappa, T_inicial, T_izq, T_der, q_bot, q_top

real(kind=8) :: T_m(nj), q_m(nj), q_mant(nj), T_mant(nj)

character(3):: sample
character(16):: outputfield

lx1 = 1.d0
lx2 = 1.d0
lx=lx1+lx2
ly = 1.d0

do ii = 1, mi
   x(ii) = lx*dfloat(ii-1)/dfloat(mi-1)
end do
do ii = 1, mi1
   x1(ii) = lx1*dfloat(ii-1)/dfloat(mi1-1)
end do
do ii = 1, mi2
   x2(ii) = lx2*dfloat(ii-1)/dfloat(mi2-1)
end do
do jj = 1, nj
   y(jj) = ly*dfloat(jj-1)/dfloat(nj-1)
end do

dx = x1(2)-x1(1)
dy = y(2)-y(1)

kappa = 0.005d0
tt = 30000
dt = 0.005d0
T_inicial = 0.d0
T_izq = 90.d0
T_der = 5.d0
q_bot = 0.d0
q_top = 0.d0

Ttotal=0

T1=T_inicial
T2=T_inicial
T1n12=T1
T2n12=T2

q_m = 0
T_m = T_inicial

q_mant=0

Tadim=0

open(25,file='TempSEC.txt',form='formatted')
open(26,file='NuSEC.txt',form='formatted')

do time=1,tt
  d1 = kappa*dt/(2.d0*dx*dx)
  d2 = kappa*dt/(2.d0*dy*dy)
  coupling=1
  iter = 1

  do while (coupling==1)
    !$OMP PARALLEL NUM_THREADS(2) DEFAULT(NONE) SHARED(T1,T2,T1n12,T2n12,Tast1,Tast2,Tbuena1,Tbuena2, &
    !$OMP q_m,T_m,d1,d2,T_der,T_izq,q_bot,q_top,dx,T_mant,q_mant,coupling, iter)

    !$OMP SECTIONS

    !$OMP SECTION

    !$OMP PARALLEL NUM_THREADS(4) DEFAULT(NONE) SHARED(T1,T1n12,Tast1,Tbuena1, &
    !$OMP q_m,T_m,d1,d2,T_izq,T_der,q_bot,q_top,dx,iter) &
    !$OMP PRIVATE(Ai1,Bi1,Ci1,Ri1,Aj1,Bj1,Cj1,Rj1)

    ! section 1 --------------------------------------------------------------

    ! x-sweep
    !$OMP DO
    do jj=2, nj-1
       Bi1(1) = 1.d0
       Ci1(1) = 0.d0
       Ri1(1) = T_izq

       do ii=2, mi1-1
          Ai1(ii-1) =-d1
          Bi1(ii)   = 1.d0+2.d0*d1
          Ci1(ii)   =-d1
          Ri1(ii)   = d2*T1(ii,jj+1)+(1.d0-2.d0*d2)*T1(ii,jj)+d2*T1(ii,jj-1)
       end do

       Ai1(mi1-1) =-1.d0
       Bi1(mi1)   = 1.d0
       Ri1(mi1)   = q_m(jj)*dx
       ! Ri1(mi1) = 100.d0*dx

       call tri(Ai1,Bi1,Ci1,Ri1,mi1)
       do ii = 1, mi1
          T1n12(ii,jj) = Ri1(ii)
       end do

     end do
     !$OMP END DO

     ! y-sweep
     !$OMP DO
     do ii = 2, mi1-1
        Bj1(1) =-1.d0
        Cj1(1) = 1.d0
        Rj1(1) = q_bot

        do jj = 2, nj-1
           Aj1(jj-1) =-d2
           Bj1(jj)   = 1.d0+2*d2
           Cj1(jj)   =-d2
           Rj1(jj)   = d1*T1n12(ii+1,jj)+(1.d0-2.d0*d1)*T1n12(ii,jj)+d1*T1n12(ii-1,jj)
        end do

        Aj1(nj-1) =-1.d0
        Bj1(nj) = 1.d0
        Rj1(nj) = q_top

        call tri(Aj1,Bj1,Cj1,Rj1,nj)
        do jj = 1, nj
           Tast1(ii,jj) = Rj1(jj)
        end do
     end do
     !$OMP END DO
     !$OMP DO
     do jj=1,nj
       T_m(jj)=T1n12(mi1,jj)
     end do
     !$OMP END DO
     !$OMP END PARALLEL

     ! -----------------------------------------------------------------------
     !$OMP SECTION
     !$OMP PARALLEL NUM_THREADS(4) DEFAULT(NONE) SHARED(T2,T2n12,Tast2,Tbuena2, &
     !$OMP q_m,T_m,d1,d2,T_der,q_bot,q_top,dx,q_mant) &
     !$OMP PRIVATE(Ai2,Bi2,Ci2,Ri2,Aj2,Bj2,Cj2,Rj2)
     ! section 2 -------------------------------------------------------------

     ! x-sweep
     !$OMP DO
     do jj=2, nj-1
        Bi2(1) = 1.d0
        Ci2(1) = 0.d0
        Ri2(1) = T_m(jj)

        do ii=2, mi2-1
           Ai2(ii-1) =-d1
           Bi2(ii)   = 1.d0+2.d0*d1
           Ci2(ii)   =-d1
           Ri2(ii)   = d2*T2(ii,jj+1)+(1.d0-2.d0*d2)*T2(ii,jj)+d2*T2(ii,jj-1)
        end do

        Ai2(mi2-1) = 0.d0
        Bi2(mi2)   = 1.d0
        Ri2(mi2)   = T_der

        call tri(Ai2,Bi2,Ci2,Ri2,mi2)
        do ii = 1, mi2
           T2n12(ii,jj) = Ri2(ii)
        end do
     end do
     !$OMP END DO

     ! y-sweep
     !$OMP DO
     do ii = 2, mi2-1
        Bj2(1) =-1.d0
        Cj2(1) = 1.d0
        Rj2(1) = q_bot

        do jj = 2, nj-1
           Aj2(jj-1) =-d2
           Bj2(jj)   = 1.d0+2*d2
           Cj2(jj)   =-d2
           Rj2(jj)   = d1*T2n12(ii+1,jj)+(1.d0-2.d0*d1)*T2n12(ii,jj)+d1*T2n12(ii-1,jj)
        end do

        Aj2(nj-1) =-1.d0
        Bj2(nj) = 1.d0
        Rj2(nj) = q_top

        call tri(Aj2,Bj2,Cj2,Rj2,nj)
        do jj = 1, nj
           Tast2(ii,jj) = Rj2(jj)
        end do
     end do

     !$OMP END DO
     !$OMP DO
     do jj=1,nj
       q_m(jj)=(Tast2(2,jj)-T2n12(1,jj))/dx
       q_m(jj)=0.4d0*q_mant(jj)+0.6d0*q_m(jj)
     end do
     !$OMP END DO
     !$OMP END PARALLEL

     !$OMP END SECTIONS

     ! -----------------------------------------------------------------------

     ! coupling
     !$OMP SINGLE

     if (maxval(dabs(T_m-T_mant)) > 0.001) then
       coupling=1
       iter = iter + 1
    else
       Tbuena1=Tast1
       Tbuena2=Tast2
       coupling=0
       iter = 1
!       print*, coupling, maxval(dabs(T_m-T_mant))
     end if
     write(101,*) iter, maxval(abs(T_m-T_mant))
     q_mant=q_m
     T_mant=T_m

     !$OMP END SINGLE

     !$OMP END PARALLEL

   end do

   T1=Tbuena1
   T2=Tbuena2

   ! --- post-processing ---

   ! one single temperature field
   do ii=1,mi1-1 ! 1 to 74
     Ttotal(ii,:)=T1(ii,:)       ! 1 to 74
   end do

   do ii=1,mi2 ! 1 to 76
     Ttotal(ii+mi1-1,:)=T2(ii,:) ! 75 a 150
   end do

   Ttotal(1,:)=T1n12(1,:)
   Ttotal(mi1,:)=T2n12(1,:)

   ! dimensionless temperature field
   do ii=1,mi-1
     do jj=2,nj-1
       Tadim(ii,jj)=(Ttotal(ii,jj)-T_der)/(T_izq-T_der)
     end do
   end do

   ! local nusselt number at the interface
   do jj=2,nj-1
     Nux(jj)=-(Tadim(mi/2+1,jj)-Tadim(mi/2,jj))/dx*lx
   end do

   ! if (time==10000) then
   ! print*, time, Tadim(75,:)
   ! end if

   ! average nusselt number at the interface
   Nu=0
   do jj=2,nj-1
     Nu=Nux(jj)*dy+Nu
   end do

   write(25,160) dt*time,Tadim(10,nj/2),Tadim(20,nj/2),Tadim(30,nj/2), &
                  & Tadim(45,nj/2),Tadim(60,nj/2),Tadim(75,nj/2), &
                  & Tadim(90,nj/2),Tadim(100,nj/2),Tadim(115,nj/2), Tadim(135,nj/2), &
                  & Nu

   if (mod(time,30000)==0) then
     do jj=2,nj-1
       write(26,170) y(jj),Nux(jj)
     end do

      i_sample=i_sample+1
      write(sample,'(i3.3)')(i_sample)
      outputfield='TempSEC_'//sample//'.vtk'

      open(78, file = Outputfield, form = 'formatted')

      write(78,110) '# vtk DataFile Version 2.3'
      write(78,110) '3D Mesh'
      write(78,111) 'ASCII'

      write(78,110) 'DATASET STRUCTURED_GRID'
      write(78,120) 'DIMENSIONS',mi, nj, 1
      write(78,130) 'POINTS', mi*nj, ' float'

      do jj=1,nj
         do ii=1,mi
            write(78,100) x(ii), y(jj), 0.0
         enddo
      enddo

      write(78,140) 'POINT_DATA', mi*nj

      write(78,110) 'SCALARS Temperature float'
      write(78,110) 'LOOKUP_TABLE default'

      do jj=1,nj
         do ii=1,mi
            write(78,100) Ttotal(ii,jj)
         end do
      end do

      close(78)
   end if

end do
close(25)
close(26)

100 FORMAT(3(f12.6));
110 FORMAT(A);
111 FORMAT(A,/);
120 FORMAT(A,I4,I4,I4);
130 FORMAT(A,I10,A);
140 FORMAT(A,I10);

150 format(f12.5)
160 FORMAT(f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8)
170 FORMAT(f12.8,f12.8)
180 FORMAT(I6,f12.8)

end program heatsec
