!-----------------------!
 module systemparameters
!-----------------------!
 implicit none

 real(8), parameter :: pi=3.141592653589793d0      ! pi
 real(8), parameter :: gme=3.9869598d14            ! G*M for earth
 real(8), parameter :: gmm=gme*0.0123d0            ! G*m for the moon
 real(8), parameter :: ts=86164d0                  ! siderial day in seconds
 real(8), parameter :: rmoon=3.844d8               ! distance from earth to moon
 real(8), parameter :: wmoon=2.d0*pi/(ts*27.25d0)  ! angular speed of moon

 real(8) :: incl                                   ! inclination (radians)
 real(8) :: cinc                                   ! cos(inclination)
 real(8) :: sinc                                   ! sin(inclination)

 end module systemparameters
!---------------------------!

!---------------------!
 program geostationary
!---------------------!
 use systemparameters; implicit none

 integer    :: i,j,nper,nstp,wstp
 real(8)    :: r,r0,dt,t,x,y,z,vx,vy,vz,tin,phi,phi0,theta,tmax

 write(*,*)'Inclination of moon orbit'
 read(*,*)incl
 incl=2.d0*pi*incl/360.d0
 cinc=cos(incl)
 sinc=sin(incl)

 write(*,*)'Period of initial orbit; T/Ts'
 read(*,*)tin
 tin=tin*ts

 write(*,*)'Number of integration days, steps/day, write-frequency'
 read(*,*)nper,nstp,wstp

 tmax=ts*nper
 dt=ts/dble(nstp)
 x=0.d0
 y=0.d0
 z=0.d0
 vx=0.d0
 vy=0.d0
 vz=0.d0
 call gorbit(tin,x,vy)

 nper=0
 phi0=0.d0
 open(10,file='orb.dat',status='replace')
 do while (t<=tmax)
    call polarcoordinates(x,y,z,r,phi,theta)   
    if (phi<phi0) nper=nper+1
    phi0=phi         
    write(10,10)t/ts,nper*2.d0*pi+phi-2.d0*pi*t/ts,theta,r/1000.d0
    10 format (f15.3,'  ',f15.6,'  ',f15.6,'  ',f15.6)
    do i=1,wstp
       call rkstep(dt,t,x,y,z,vx,vy,vz)
    enddo
 enddo
 close(10)

 end program geostationary
!-------------------------!

!---------------------------------------------!
 subroutine rkstep(dt,t0,x0,y0,z0,vx0,vy0,vz0)
!--------------------------------------------------------------------------------!
! Integrates the equations of motion one time step, using the Runge-Kutta method !
!--------------------------------------------------------------------------------!
 implicit none

 real(8) :: ax,ay,az,x0,y0,z0,x1,y1,z1,vx0,vy0,vz0,vx1,vy1,vz1,dt,t0,th,t1,dt2
 real(8) :: kx1,kx2,kx3,kx4,ky1,ky2,ky3,ky4,kz1,kz2,kz3,kz4
 real(8) :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4,lz1,lz2,lz3,lz4

 dt2=dt/2.d0
 t1=t0+dt; th=t0+dt2
 call accel(x0,y0,z0,t0,ax,ay,az) 
 kx1=dt2*ax 
 ky1=dt2*ay
 kz1=dt2*az
 lx1=dt2*vx0 
 ly1=dt2*vy0
 lz1=dt2*vz0
 call accel(x0+lx1,y0+ly1,z0+lz1,th,ax,ay,az)
 kx2=dt2*ax; ky2=dt2*ay; kz2=dt2*az
 lx2=dt2*(vx0+kx1)
 ly2=dt2*(vy0+ky1)
 lz2=dt2*(vz0+kz1)
 call accel(x0+lx2,y0+ly2,z0+lz2,th,ax,ay,az)
 kx3=dt*ax
 ky3=dt*ay
 kz3=dt*az
 lx3=dt*(vx0+kx2)
 ly3=dt*(vy0+ky2)
 lz3=dt*(vz0+kz2)
 call accel(x0+lx3,y0+ly3,z0+lz3,t1,ax,ay,az)
 kx4=dt2*ax
 ky4=dt2*ay
 kz4=dt2*az
 lx4=dt2*(vx0+kx3)
 ly4=dt2*(vy0+ky3)
 lz4=dt2*(vz0+kz3)
 x1=x0+(lx1+2.d0*lx2+lx3+lx4)/3.d0
 y1=y0+(ly1+2.d0*ly2+ly3+ly4)/3.d0
 z1=z0+(lz1+2.d0*lz2+lz3+lz4)/3.d0
 vx1=vx0+(kx1+2.d0*kx2+kx3+kx4)/3.d0
 vy1=vy0+(ky1+2.d0*ky2+ky3+ky4)/3.d0
 vz1=vz0+(kz1+2.d0*kz2+kz3+kz4)/3.d0

 t0=t1
 x0=x1; y0=y1; z0=z1
 vx0=vx1; vy0=vy1; vz0=vz1

 end subroutine rkstep
!---------------------!

!----------------------------------!
 subroutine accel(x,y,z,t,ax,ay,az)
!----------------------------------!
 use systemparameters; implicit none

 real(8) :: x,y,z,t,ax,ay,az,xm,ym,zm,re3,rm3

 re3=gme*(x**2+y**2+z**2)**(-1.5d0)
 xm=rmoon*cinc*cos(wmoon*t)
 ym=rmoon*sin(wmoon*t)
 zm=rmoon*sinc*cos(wmoon*t)
 xm=x-xm
 ym=y-ym
 zm=z-zm
 rm3=gmm*(xm**2+ym**2+zm**2)**(-1.5)
 
 ax=-x*re3-xm*rm3
 ay=-y*re3-ym*rm3
 az=-z*re3-zm*rm3

 end subroutine accel
!--------------------!

!----------------------------------------------!
 subroutine polarcoordinates(x,y,z,r,phi,theta)
!----------------------------------------------!
 use systemparameters; implicit none
 
 real(8) :: x,y,z,r,phi,phi0,theta

 r=sqrt(x**2+y**2+z**2)
 if (y>=0.d0) then
    phi=acos(x/r)
 else
    phi=2.d0*pi-acos(x/r)
 endif
 theta=asin(z/r)

 end subroutine polarcoordinates
!--------------------------------!

!----------------------------------------------------------------------!
!Calculates the radius rgeo and velocity vgeo of an orbit with period t!
!----------------------------------------------------------------------!
 subroutine gorbit(t,r,v)
!------------------------!
 use systemparameters; implicit none

 real(8) :: t,r,v

 r=((gme*t**2)/(4.d0*pi**2))**(1.d0/3.d0) 
 v=sqrt(gme/r)               

 end subroutine gorbit
!---------------------!



