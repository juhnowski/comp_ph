!--------------------------------------------------------------------------------!
! - Solves the Schrodinger equation in one dimension, for a potential given      !
!   by the function 'potential(x)', for x in the range (-xmax,xmax).             !
! - The  integration of the wave functipn Psi(x) starts from boundary conditions !
!   applied to Psi(-xmax) and Psi(-xmax+h), where h is the integration step      !
!   (this boundary condition is set in the subroutine 'setinitcond').            !
! - The boundary condition at x=xmax is taken to be Psi(xmax)=0.                 !
! - Bisection is used to find the energy (preceded by a course search            !
!   for two energies bracketing a solution).                                     !
!--------------------------------------------------------------------------------!
 program schrodinger1d
!---------------------!
 implicit none

 real(8), parameter :: eps=1.d-6

 integer :: i,ni,nx
 real(8) :: xmax,de,b0,b1,b2,e0,e1,e2

 real(8), allocatable :: psi(:)

 write(*,*)'integrating between -xmax and xmax; give xmax:'
 read(*,*)xmax
 write(*,*)'Number of points on either side of 0'
 read(*,*)nx
 write(*,*)'searching for E>E0; give E0'
 read(*,*)e1
 write(*,*)'Delta-E in initial course search'
 read(*,*)de
 
 allocate(psi(-nx:nx))

 ni=0
 write(*,*)'Starting course search'
 call integrate(xmax,nx,e1,psi)
 b1=psi(nx)
 call writemessage(ni,e1,b1) 
 do 
   ni=ni+1
   e2=e1+de
   call integrate(xmax,nx,e2,psi)
   b2=psi(nx)
   call writemessage(ni,e2,b2) 
   if (b1*b2 < 0.d0) exit
   e1=e2
   b1=b2
 enddo

 ni=0
 write(*,*)'Starting bisection'
 do
   ni=ni+1
   e0=(e1+e2)/2.0d0
   call integrate(xmax,nx,e0,psi)
   b0=psi(nx)
   call writemessage(ni,e0,b0) 
   if (abs(b0) <= eps) exit
   if (b0*b1 <= 0.0d0) then
     e2=e0      
     b2=b0 
   else
     e1=e0
     b1=b0
   endif
 enddo

 open(1,file='psi.dat',status='replace')
 do i=-nx,nx
   write(1,*)dble(i)*xmax/dble(nx),psi(i)
 enddo
 close(1)

 end program schrodinger1d
!-------------------------!

!------------------------------------!
 subroutine integrate(xmax,nx,e1,psi)
!------------------------------------!
 implicit none

 real(8) :: q0,q1,p1,f1,ee,h2,h12
 common/block1/q0,q1,p1,f1,ee,h2,h12

 integer :: i,n,nx
 real(8) :: psi(-nx:nx),xmax,x,h,de,e1

 ee=e1
 h=xmax/dble(nx)      
 h2=h**2
 h12=h2/12.d0
 call setinitcond(xmax,h,psi(-nx),psi(-nx+1))
 do i=-nx+2,nx
   x=dble(i)*h
   call numerovstep(x)
   psi(i)=p1
 end do 
 call normalize(psi,nx,h)

 end subroutine integrate 
!------------------------!

!----------------------------------------------------------------!
!Integrates the wave function one step, using Numerov's method.  !
!On entry, the wave-function at the previous step is p1, on exit !
!it is the wave-function after the integration step. The modified! 
!wave function, Psi(x)*[1-(h**2)*f1], at the two previous steps  !
!are q0,q1 on entry, and f1=2.d0*(potential(x1)-energy).         !
!----------------------------------------------------------------!
 subroutine numerovstep(x)
!---------------------------!
 implicit none

 real(8) :: q0,q1,p1,f1,ee,h2,h12
 common/block1/q0,q1,p1,f1,ee,h2,h12
       
 real(8) :: x,q2
 real(8), external :: potential

 q2=h2*f1*p1+2.d0*q1-q0
 q0=q1; q1=q2
 f1=2.d0*(potential(x)-ee)        
 p1=q1/(1.d0-h12*f1)
         
 end subroutine numerovstep
!--------------------------!

!---------------------------------------------------!
!Potential energy as a function of the coordinate x !
!Three examples are given (two commented out).      !
!---------------------------------------------------!
 real(8) function potential(x)
!-----------------------------!
 implicit none

 real(8) :: x

!----- square well with infinite walls       
! potential=0.d0

!----- square well with soft walls V(|x|>1)=10
!if (abs(x) < 1.d0) then
!  potential=0.d0
!else
!  potential=10.d0
!endif 

!----- square well with hard walls internal and box between -0.5 and 0.5
 if (x > -0.5d0 .and. x < 0.5d0) then
   potential=5.d0
 else
   potential=0.d0
 endif 

 end function potential
!----------------------!

!----------------------------------------------------------------!
!Sets the boundary conditions Psi(-xmax)=0, Psi(-xmax+h)=0.0001  !
!- note that the boundary condition is strictly not correct for  !
!  a "sof" potential, but OK if xmax is large enough             !
!----------------------------------------------------------------!
 subroutine setinitcond(xmax,h,psi0,psi1)
!----------------------------------------!
 implicit none

 real(8) :: q0,q1,p1,f1,ee,h2,h12
 common/block1/q0,q1,p1,f1,ee,h2,h12

 real(8) :: h,vv,xmax,psi0,psi1
 real(8), external :: potential

 psi0=0.d0
 psi1=0.0001d0

 p1=psi0
 f1=2.d0*(potential(-xmax)-ee)        
 q0=psi0*(1.d0-h12*f1)
 f1=2.d0*(potential(-xmax+h)-ee)        
 q1=psi1*(1.d0-h12*f1)
	
 end subroutine setinitcond
!--------------------------!

!------------------------------!
 subroutine normalize(psi,nx,h)
!------------------------------!
 implicit none

 integer :: i,nx
 real(8) :: psi(-nx:nx),norm,h

 norm=psi(-nx)**2+psi(nx)**2
 do i=-nx+1,nx-3,2
   norm=norm+4.d0*psi(i)**2+2.d0*psi(i+1)**2
 end do
 norm=norm+4.d0*psi(nx-1)**2
 norm=1.d0/sqrt(norm*h/3.d0)
 do i=-nx,nx
   psi(i)=psi(i)*norm
 end do

 end subroutine normalize
!------------------------!

!------------------------------!
 subroutine writemessage(n,e,b) 
!------------------------------!
 implicit none

 integer :: n
 real(8) :: e,b

 write(*,1)n,e,b
 1 format(i4,'    E = ',f10.7,'    Boundary deviation ',f10.7)

 end subroutine writemessage
!---------------------------!
