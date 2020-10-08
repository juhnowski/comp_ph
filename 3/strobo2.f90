!---------------------------------------------------------------------------!
!This program produces a stroboscopic sampling plot for a driven, damped    !
!pendulum. The output is given in the form of a postscript file 'strobo.ps'.!
!Input is read from the file 'read.in', containing:                         !
!  ksprng,mass,gamma,dramp,drfrq                                            !
!  ndt,ntmax,ntrans                                                         !
!  x0,v0                                                                    !
!The first line contains the pendulum parameters.                           !
!The second line contains:                                                        !
!  ndt = number of time-steps for each period of the driving force          !
!  ntmax = the total integration time as a multiple of the driving period   !
!  ntrans = number of periods to skip before writing (to avoid transients)  !
!The third line contains the initial conditions.                            !
!---------------------------------------------------------------------------!

!----------------------!
 program strobopendulum
 implicit none

 real(8), parameter :: pi=3.14159265358979328d0
 real(8), parameter :: pi2=2.d0*pi

 real(8) ::  ksprng,mass,gamma,drfrq,dramp
 common/block1/ksprng,mass,gamma,drfrq,dramp

 integer :: i,ip,it,ndt,ntmax,ntrans,ndramp
 real(8) ::  x0,v0,x1,v1,t0,t1,th,k1,k2,k3,k4,l1,l2,l3,l4,dt,dt2,accel

 open(1,file='read.in',status='old')
 read(1,*)ksprng,mass,gamma,dramp,drfrq
 read(1,*)ndt,ntmax,ntrans
 read(1,*)x0,v0
 close(1)

 ksprng=ksprng/mass
 gamma=gamma/mass
 dramp=dramp/mass
 dt=pi2/(drfrq*dble(ndt))
 dt2=dt/2.d0

 call initgraphic(drfrq)
 t0=0.d0
 do ip=0,ntmax-1
   if (ip.ge.ntrans) then
     write(1,1)x0*30.d0,v0*30.d0
     1 format(f7.2,' ',f7.2,' a')
   endif
   do it=1,ndt
     th=(pi2/drfrq)*dble(ip)+dt*(dble(it)-0.5d0)   
     t1=(pi2/drfrq)*dble(ip)+dt*(dble(it))   
     k1=dt2*accel(x0,v0,t0) 
     l1=dt2*v0
     k2=dt2*accel(x0+l1,v0+k1,th)
     l2=dt2*(v0+k1)
     k3=dt*accel(x0+l2,v0+k2,th)
     l3=dt*(v0+k2)
     k4=dt2*accel(x0+l3,v0+k3,t1)
     l4=dt2*(v0+k3)
     x1=x0+(l1+2.d0*l2+l3+l4)/3.d0
     v1=v0+(k1+2.d0*k2+k3+k4)/3.d0
     if (x1<pi) x1=x1+pi2
     if (x1>pi) x1=x1-pi2
     x0=x1
     v0=v1
     t0=t1
   end do
 end do
 call closegraphic
 
 end program strobopendulum
!--------------------------!

!--------------------------!
 real(8) function accel(x,v,t)
 implicit none

 real(8) :: ksprng,mass,gamma,drfrq,dramp
 common/block1/ksprng,mass,gamma,drfrq,dramp

 real(8) :: x,v,t

 accel=dramp*sin(t*drfrq)-ksprng*sin(x)-gamma*v

 end function accel
!------------------!

!------------------------!
 subroutine initgraphic(q)
 implicit none

 real(8) :: q

 open(1,file='strobo.ps',status='replace')
 write(1,*)'%!'
 write(1,*)'1 1 scale 0 0 translate 0.5 0.5 moveto'
 write(1,*)'600 0 rlineto 0 600 rlineto'
 write(1,*)'-600 0 rlineto 0 -600 rlineto'
 write(1,*)'stroke'
 write(1,*)'6 6 scale 0.2 setlinewidth'
 write(1,*)'50 50 translate 0.45 0.45 scale'
 write(1,*)'0 -100 moveto 0 200 rlineto stroke'
 write(1,*)'-100 0 moveto 200  0 rlineto stroke'
 write(1,*)'/TimesRoman findfont 14 scalefont setfont'
 write(1,*)'-3.25 102 moveto (v) show 102 -3 moveto (x) show'
 write(1,*)'/a {newpath 0.5 0 360 arc fill} def'
 write(1,*)'/TimesRoman findfont 14 scalefont setfont'
 write(1,1)' -107 99 moveto (Q=',q,') show'
 1 format(a,f7.5,a)

 end subroutine initgraphic
!--------------------------!


!-----------------------!
 subroutine closegraphic
 implicit none

 write(1,*)'showpage'
 close(1)

 end subroutine closegraphic
!---------------------------!
