!-------------!
 module system
!-------------!

 real(8), parameter :: rcore=0.1d0        !hard-core radius
 real(8), parameter :: beta=0.053667d0    !energy parameter, corr. to 2.226MeV
 real(8), parameter :: eps=1.d-7          !accuracy of r=r0 boundary condition

 integer :: nr
 real(8) :: v0,a
 real(8) :: h,h2,h12,rmax
 real(8) :: q0,q1,p1,f1

 real(8), allocatable :: psi(:)

 end module system
!-----------------!

!============================================================================!
!This program calculates the wave function of the deuteron, using a Yukawa   !
!potential of depth V0 and width a. For a range of a-values, the depth V0    !
!is found which gives a valid wave function when the energy E=2.226MeV, the  !
!binding energy of the deuteron. The deuteron radius is calculated vs a.     !
!The potential depth V0 and the radius r are written to the file 'vr.dat'.
!============================================================================!
 program deuteron
 use system
 implicit none

 integer :: i,na
 real(8) :: a1,a2,va,dv,dev,r

 write(*,*)'Maximum integration radius'
 read(*,*)rmax
 write(*,*)'Number of integration points'
 read(*,*)nr
 write(*,*)'Delta-V in course search'
 read(*,*)dv
 write(*,*)'Potential widths from a1 to a2'
 read(*,*)a1,a2
 write(*,*)'Number of widths to study'
 read(*,*)na

 allocate(psi(0:nr))

 h=(rmax-rcore)/dble(nr)       
 h2=h**2
 h12=h2/12.d0

 dev=10.d0
 open(1,file='vr.dat',status='replace')
 do i=0,na
   write(*,*)'Doing width ',i
   a=a1+dble(i)*(a2-a1)/dble(na)
   call potentialdeapth(dv)
   call radius(r)
   write(1,1)a,v0,r
   1 format(3f12.6)
   if (abs(r-2.d0) < dev) then
      va=a
      dev=abs(r-2.d0)
   endif
 enddo
 close(1)

 write(*,*)'Best parameter a:',va

 call potentialdeapth(dv)
 open(1,file='psi.dat',status='replace')
 do i=0,nr,10
   write(1,1)rcore+dble(i)*h,psi(i)
 end do
 close(1)

 deallocate(psi)

 end program deuteron
!====================!

!-----------------------------------------!
!Bisection search for the potential depth.!
!-----------------------------------------!
 subroutine potentialdeapth(dv)
 use system; implicit none

 integer :: i
 real(8) :: dv,b0,b1,b2,v1,v2

 v0=0.d0; v1=v0
 call integrate; b1=psi(0)
 do 
   v2=v1+dv; v0=v2
   call integrate; b2=psi(0)
   if (b1*b2 < 0.d0) exit
   v1=v2; b1=b2
 end do
 do
   v0=(v1+v2)/2.0d0
   call integrate; b0=psi(0)
   if (abs(b0) <= eps) exit
   if (b0*b1 <= 0.0d0) then
     v2=v0; b2=b0 
   else
     v1=v0; b1=b0
   end if
 end do

 end subroutine potentialdeapth
!------------------------------!

!----------------------------------------------------!
!Integration of the wave function, starting from rmax!
!and integrating toward the hard-core radius r0.     !
!----------------------------------------------------!
 subroutine integrate()
 use system; implicit none

 integer :: i,n
 real(8) :: r,vv

 call setinitcond
 do i=nr-2,0,-1
   r=rcore+dble(i)*h
   call numerovstep(r)
   psi(i)=p1
 end do 
 call normalize()

 end subroutine integrate 
!------------------------!

!--------------------------------------------!
!Numerov's algorith for one integration step.!
!--------------------------------------------!
 subroutine numerovstep(r)
 use system; implicit none
       
 real(8) :: r,q2,potentialminuse

 q2=h2*f1*p1+2.d0*q1-q0
 q0=q1
 q1=q2
 f1=potentialminuse(r)        
 p1=q1/(1.d0-h12*f1)
         
 end subroutine numerovstep
!--------------------------!

!-------------------------------------------------------------!
!Sets the boundary conditions Psi(rmax)=psi0, Psi(rmax-h)=psi1!
!-------------------------------------------------------------!
 subroutine setinitcond
 use system; implicit none

 real(8) :: potentialminuse

! psi(nr)=1.d0
! psi(nr-1)=exp(h*sqrt(potentialminuse(rmax-h)))
 psi(nr)=0.001
 psi(nr-1)=0.002

 p1=psi(nr-1)
 f1=potentialminuse(rmax)        
 q0=psi(nr)*(1.d0-h12*f1)
 f1=potentialminuse(rmax-h)
 q1=p1*(1.d0-h12*f1)
	
 end subroutine setinitcond
!--------------------------!

!----------------------------------------------!
!V-E for Yukawa potential with parameters v0,a.!
!----------------------------------------------!
 real(8) function potentialminuse(r)
 use system; implicit none

 real(8) :: r

 potentialminuse=beta*(1.d0-v0*exp(-r/a)/(r/a))
 
 end function potentialminuse
!----------------------------!

!----------------------------!
!Normalizes the wave function!
!----------------------------!
 subroutine normalize
 use system; implicit none

 integer :: i
 real(8) :: norm

 norm=psi(0)**2
 do i=1,nr-3,2
   norm=norm+4.d0*psi(i)**2+2.d0*psi(i+1)**2
 end do
 norm=norm+4.d0*psi(nr-1)**2+psi(nr)**2
 norm=1.d0/sqrt(norm*h/3.d0)
 do i=0,nr
   psi(i)=psi(i)*norm
 end do

 end subroutine normalize
!------------------------!

!---------------------------------------------!
!Calculates the deuteron radius, using Simpson! 
!integration of the radial wave function.     !
!---------------------------------------------!
 subroutine radius(r)
!--------------------!
 use system; implicit none

 integer :: i
 real(8) :: r,r1,r2,r3,rn

 rn=(psi(0)**2)
 r3=(psi(0)**2)*rcore**2
 do i=1,nr-3,2
   r1=rcore+h*dble(i)
   rn=rn+4.d0*(psi(i)**2)
   r3=r3+4.d0*(psi(i)**2)*r1**2
   r1=rcore+h*dble(i+1)
   rn=rn+2.d0*(psi(i+1)**2)
   r3=r3+2.d0*(psi(i+1)**2)*r1**2
 end do
 r1=rmax-h
 rn=rn+4.d0*(psi(nr-1)**2)
 r3=r3+4.d0*(psi(nr-1)**2)*r1**2
 rn=rn+(psi(nr)**2)
 r3=r3+(psi(nr)**2)*rmax**2
 r=sqrt(r3/rn)/2.d0

 end subroutine radius
!---------------------!
