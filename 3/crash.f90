!------------------!
 module systemparam
 implicit none

 real(8), parameter :: pi=3.141592653589793d0 
 real(8), parameter :: gm=3.987d14  ! Earth's mass times constant of gravity
 real(8), parameter :: arocket=5.d0 ! deceleration due to rocket motor
 real(8), parameter :: dragc=8.d-4  ! Air drag coefficient
 real(8), parameter :: re=6.378d6   ! Earth's radius

 real(8) :: dt,dt2  ! time step, half of the time step
 real(8) :: tbrake  ! run-time of rocket motor

 end module systemparam
!----------------------!


!======================!
 program satellitecrash
!------------------------------------------------------------------!
!calculates the trajectory of a satellite that is initially in a   !
!circular orbit and begins to brake at the start of the calculation!
!------------------------------------------------------------------!
 use systemparam
 implicit none

 integer :: i,nstp,wstp
 real(8) :: a,r,t,x,y,vx,vy,r0,tmax

 print*,'Initial altitude of satellite (km)'; read*,r0
 r0=r0*1.d3+re
 print*,'Rocket motor run-time (seconds)'; read*,tbrake
 print*,'Time step delta-t for RK integration (seconds)';read*,dt
 dt2=dt/2.d0
 print*,'Writing results every Nth step; give N';read*,wstp
 print*,'Maximum integration time (hours)';read*,tmax
 tmax=tmax*3600.d0

 x=r0
 y=0.d0
 vx=0.d0
 vy=sqrt(gm/r0)
 nstp=int(tmax/dt)

 open(1,file='sat.dat',status='replace')
 do i=0,nstp
   call polarposition(x,y,r,a)
   if (r > re) then
     t=dble(i)*dt
     if (mod(i,wstp)==0) write(1,1)t,a,(r-re)/1.d3,sqrt(vx**2+vy**2)
     1 format(f12.3,'  ',f12.8,'  ',f14.6,'  ',f12.4)
     call rkstep(t,x,y,vx,vy)
   else
     print*,'The satellite has successfully crashed!'
     goto 2
   end if
 end do
 print*,'The satellite did not crash within the specified time.'
 2 close(1)

 end program satellitecrash
!==========================!

!-----------------------------------!
 subroutine rkstep(t0,x0,y0,vx0,vy0)
!---------------------------------------!
!Integrates the equations of motion one !
!time step, using the Runge-Kutta method!
!---------------------------------------!
 use systemparam
 implicit none

 real(8) :: x0,y0,x1,y1,vx0,vy0,vx1,vy1,t0,th,t1,ax,ay
 real(8) :: kx1,kx2,kx3,kx4,ky1,ky2,ky3,ky4
 real(8) :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4

 t1=t0+dt; th=t0+dt2   
 call accel(x0,y0,vx0,vy0,t0,ax,ay) 
 kx1=dt2*ax 
 ky1=dt2*ay
 lx1=dt2*vx0 
 ly1=dt2*vy0
 call accel(x0+lx1,y0+ly1,vx0+kx1,vy0+ky1,th,ax,ay)
 kx2=dt2*ax; ky2=dt2*ay
 lx2=dt2*(vx0+kx1)
 ly2=dt2*(vy0+ky1)
 call accel(x0+lx2,y0+ly2,vx0+kx2,vy0+ky2,th,ax,ay)
 kx3=dt*ax
 ky3=dt*ay
 lx3=dt*(vx0+kx2)
 ly3=dt*(vy0+ky2)
 call accel(x0+lx3,y0+ly3,vx0+kx3,vy0+ky3,t1,ax,ay)
 kx4=dt2*ax
 ky4=dt2*ay
 lx4=dt2*(vx0+kx3)
 ly4=dt2*(vy0+ky3)
 x1=x0+(lx1+2.d0*lx2+lx3+lx4)/3.d0
 y1=y0+(ly1+2.d0*ly2+ly3+ly4)/3.d0
 vx1=vx0+(kx1+2.d0*kx2+kx3+kx4)/3.d0
 vy1=vy0+(ky1+2.d0*ky2+ky3+ky4)/3.d0

 x0=x1; y0=y1; vx0=vx1; vy0=vy1

 end subroutine rkstep
!---------------------!

!-----------------------------------!
 subroutine accel(x,y,vx,vy,t,ax,ay)
!------------------------------------------------------!
!calculates the x and y components of the acceleration,!
!due to gravitation, air drag, and brake-motor thrust  !
!------------------------------------------------------!
 use systemparam
 implicit none

 real(8) :: x,y,vx,vy,t,ax,ay,r,v1,v2,r3,ad
 real(8), external :: airdens

 r=sqrt(x**2+y**2)       
 v2=vx**2+vy**2
 v1=sqrt(v2)

   !*** evaluates the acceleration due to gravitation
 r3=1.d0/r**3
 ax=-gm*x*r3           
 ay=-gm*y*r3  

   !*** evaluates the acceleration due to air drag
 if (v1 > 1.d-12) then
   ad=dragc*airdens(r)*v2
   ax=ax-ad*vx/v1
   ay=ay-ad*vy/v1
 end if

   !*** evaluates the acceleration due to rocket motor thrust
 if (t < tbrake .and. v1 > 1.d-12) then
   ax=ax-arocket*vx/v1
   ay=ay-arocket*vy/v1
 end if

 end subroutine accel
!--------------------!

!---------------------------!
 real(8) function airdens(r)
!-------------------------------------------------!
!the assumed air-density as a function of altitude!
!-------------------------------------------------!
 use systemparam
 implicit none

 real(8), parameter :: k1=1.2d4
 real(8), parameter :: k2=2.2d4

 real(8) :: r,rr

 if (r > re) then
    airdens=1.225d0*exp(-((r-re)/k1+((r-re)/k2)**1.5))
 else
    airdens=0.d0
    return
 endif

 end function airdens
!--------------------!

!---------------------------------!
 subroutine polarposition(x,y,r,a)
!------------------------------------------------------!
!converts the x and y coordinates to a distance r from !
!the center of the earth and a normalized angle a (0-1)!
!------------------------------------------------------!
 use systemparam
 implicit none
 
 real(8) :: x,y,r,a

 r=sqrt(x**2+y**2)
 if (y >= 0.d0) then
   a=acos(x/r)/(2.d0*pi)
 else
   a=1.d0-acos(x/r)/(2.d0*pi)
 end if

 end subroutine polarposition
!----------------------------!
