!---------------!
 program romberg
!---------------------------------------------------------!
! Romberg integration of a function implemented as func(x)
!---------------------------------------------------------!
 implicit none

 integer :: i,n,steps
 real(8) :: x0,xn,t

 real(8), external    :: polext0
 real(8), allocatable :: trap(:)

 write(*,'(a)',advance='no')'Integration bounds x0,xn: '
 read*,x0,xn
 write(*,'(a)',advance='no')'Number of Romberg steps: '
 read*,steps
 write(*,'(a)',advance='no')'Initial number of intervals: '
 read*,n

 allocate(trap(0:steps))

! loop over trapezodial integrations using n*2**i terms in the
! sum over function values. After each loop iteration, t contains
! the value after the last trapezoidal integration. The function
! trapezoidal() needs the previous value of t in order to evaluate
! the following one.

 do i=0,steps-1
   call trapezoidal(x0,xn,n,t,i)   
   trap(i)=t
   if (i/=0) write(*,*)i,trap(i),polext0(i,trap)
 enddo

 deallocate(trap)

 end program romberg
!-------------------!

!-----------------------------------------!
 subroutine trapezoidal(x0,xn,n,trap,step)
!----------------------------------------------------------------!
! Carries out an integration of the function func(x) between
! x=x0 and x=xn. The number of integration intervals is n*2**step.
! If step>0, the result of the previous summation (at step-1) is
! used, so that only the new n/2 new points added between the old
! grid points have to be used.
!-----------------------------------------------------------------
 implicit none

 integer :: i,n,step
 real(8) :: h,h2,x0,xn,x0h,trap,func

 h=(xn-x0)/dble(n*2**step)  ! the step size
 if (step == 0) then
   trap=0.5d0*(func(x0)+func(xn))
   do i=1,n-1
     trap=trap+func(x0+h*dble(i))
   enddo
 else   
   h2=2.d0*h                ! the step size at the previous step
   x0h=x0+h
   trap=trap/h2             ! the value of the previous sum (without h)
   do i=0,n*2**(step-1)-1
      trap=trap+func(x0h+h2*dble(i))   ! the new sum is constructed,
   enddo                               ! skipping the points already included
 endif                                         
 trap=trap*h

 end subroutine trapezoidal
!--------------------------!

!---------------------!
 function polext0(n,y)
!----------------------------------------------------------------!
! Evaluates at x=0 the n:th order polynomial going throigh points 
! y(xi), i=0,...n, for xi of the form x0/2**i, using Legendre's
! formula for the interpolating polynomial.
!----------------------------------------------------------------!
 implicit none

 integer :: j,k,n
 real(8) :: y(0:n),p,polext0

 polext0=0.d0
 do k=0,n
    p=1.d0
    do j=0,n
       if (j/=k) p=p/(2.d0**(2*(j-k))-1.d0)
    enddo
   polext0=polext0+p*y(k)*(-1)**n
 enddo

 end function polext0
!--------------------!

!----------------!
 function func(x)
!---------------------------------------------------!
! The function to be integrated, defined by the user
!---------------------------------------------------!
 implicit none

 real(8) :: x,func

 func=log(x+sqrt(x**2+1))*x**4

 end function func
!-----------------!
