!-----------------------------------------------------------------------------------!
! Variational calculation for a one-dimensional "particle in a box" with a square   !
! potential barrier of height V between points x1,x2. The locations of the infinite !
! walls are x=-1 and x=+1 (so therefore |x1|,|x2| < 1). The basis set used is the   !
! one of eigenfunctions in the absence of internal barrier (i.e., the solution of   !
! the problem when V=0). The basis is truncated to the N first eigenfunctions. One  !
! of the states k=1,...,N is written to files; the expansion coefficients in the    !
! unperturbaed basis in "psik.dat" and the real-space wave function in "psir.dat".  !
! The energy eigenvalues are written to the file "eig.dat".                         !
!-----------------------------------------------------------------------------------!

!-------------------!
 program variational
!-------------------!
 implicit none

 integer, parameter :: npsi=1000

 integer :: i,j,pn,nbasis
 real(8) :: x1,x2,vv,energ0,velement,psi(0:npsi)
 real(8), allocatable :: eig(:)
 real(8), allocatable :: ham(:,:)

 print*,'Internal barrier between x1 and x2; give x1,x2'
 read(*,*)x1,x2
 print*, 'Internal potential V'
 read(*,*)vv 

 print*,'Number of basis states'
 read(*,*)nbasis
 print*,'Write eigenfunction #'
 read(*,*)pn

 allocate(eig(nbasis))
 allocate(ham(nbasis,nbasis))
 
 ham=0.d0
 do j=1,nbasis
   ham(j,j)=energ0(j)+velement(j,j,x1,x2,vv)
   do i=1,j-1
      ham(i,j)=velement(i,j,x1,x2,vv)
      ham(j,i)=ham(i,j)
   end do
 end do
 call diasym(ham,eig,nbasis)

 open(1,file='eig.dat',status='replace')
 do i=1,nbasis
    write(1,1)i,eig(i)
 end do
 1 format(i4,' ',f12.7)
 close (1)

 open(1,file='psik.dat',status='replace')
 do i=1,nbasis
    write(1,1)i,ham(i,pn)
 end do
 close (1)

 call fullwavefunction(ham(:,pn),nbasis,psi,npsi)
 open(1,file='psir.dat',status='replace')
 do i=0,npsi
    write(1,2)dble(2*i)/dble(npsi)-1.d0,psi(i)
 end do
 2 format(f8.5,' ',f12.7)
 close (1)

 end program variational
!-----------------------!

!---------------------------------------!
 real(8) function velement(p,k,x1,x2,vv)
 implicit none

 integer,parameter :: romsteps=3

 integer :: i,p,k,n0
 real(8) :: x1,x2,vv,tint,polext0,trap(0:romsteps)

 n0=4*max(p,k)
 do i=0,romsteps
   call trapezoidal(p,k,x1,x2,n0,tint,i)
   trap(i)=tint
 end do
 velement=vv*polext0(romsteps,trap)
 
 end function velement
!---------------------!

!--------------------------------------------!
 subroutine trapezoidal(p,k,x0,xn,n,trap,step)
 implicit none

 integer :: i,p,k,n,step
 real(8) :: h,h2,x,x0,xn,x0h,trap,wavefunc

 h=(xn-x0)/dble(n*2**step)
 if (step == 0) then
   trap=0.5d0*wavefunc(p,x0)*wavefunc(k,x0)
   do i=1,n-1
     x=x0+h*dble(i)
     trap=trap+wavefunc(p,x)*wavefunc(k,x)
   end do
   trap=trap+0.5d0*wavefunc(p,xn)*wavefunc(k,xn)
 else   
   h2=2.d0*h
   x0h=x0+h
   trap=trap/h2
   do i=0,n*2**(step-1)-1
     x=x0h+h2*dble(i)
     trap=trap+wavefunc(p,x)*wavefunc(k,x)
   end do 
 end if
 trap=trap*h

 end subroutine trapezoidal
!--------------------------!

!-----------------------------!
 real(8) function polext0(n,y)
 implicit none

 integer :: j,k,n
 real(8) :: y(0:n),p

 polext0=0.d0
 do k=0,n
    p=1.d0
    do j=0,n
       if (j /= k) p=p/(2.d0**(2*(j-k))-1.d0)
    end do
    polext0=polext0+p*y(k)*(-1)**n
 end do

 end function polext0
!--------------------!

!-----------------------------!
 real(8) function wavefunc(n,x)
 implicit none

 real(8), parameter :: pi2=3.14159265358979328d0/2.d0

 integer :: n
 real(8) :: x

 if (mod(n,2)==0) then
    wavefunc=sin(dble(n)*pi2*x)
 else
    wavefunc=cos(dble(n)*pi2*x)
 end if

 end function wavefunc
!---------------------!

!-----------------------------------------------!
 subroutine fullwavefunction(vec,nbasis,psi,npsi)
 implicit none
 
 integer :: i,k,nbasis,npsi
 real(8) :: vec(nbasis),psi(0:npsi),wavefunc,x

 psi=0.d0
 do i=0,npsi
   x=-1.d0+dble(2*i)/dble(npsi)
   do k=1,nbasis
     psi(i)=psi(i)+vec(k)*wavefunc(k,x)
   end do
 end do

 end subroutine fullwavefunction
!-------------------------------! 

!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!

!--------------------------!
 subroutine diasym(a,eig,n)
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!---------------------!

!-------------------------!
 real(8) function energ0(n)
 implicit none

 real(8), parameter :: pi=3.14159265358979328d0
 integer :: n

 energ0=0.125d0*(pi*dble(n))**2

 end function energ0
!-------------------!

