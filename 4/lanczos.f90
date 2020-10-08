!---------------------------------------------------------------------------!
! Lanczos calculation for a cubic d-dimensional box, discretized using L**d ! 
! elements. Energies are written to a file 'e.dat', and a selected state    !
! (its wave function squared) is written to 'p.dat'. The potential energy   !
! is stored in a vector vpot, which is constructed in the subroutine        !
! potentialenergy. The basis construction starts from a random state, for   !
! which a file seed.in containing four random seed integers is needed.      !
!---------------------------------------------------------------------------!

!-----------------------!
 module systemparameters
!-----------------------!

 
 integer :: d
 integer :: l
 real(8), allocatable :: vpot(:)

 end module systemparameters
!---------------------------!

!---------------------!
 module lanczosvectors
!---------------------!

 real(8), save, allocatable :: f0(:)
 real(8), save, allocatable :: f1(:)
 real(8), save, allocatable :: f2(:)
 real(8), save, allocatable :: aa(:)
 real(8), save, allocatable :: nn(:)

 end module lanczosvectors
!-------------------------!

!===================!
 program lanczoscube
!=======================================================!
 use systemparameters; use lanczosvectors; implicit none

 integer :: i,niter,st
 real(8) :: len,delta,escale

 real(8), allocatable :: psi(:)
 real(8), allocatable :: eig(:)
 real(8), allocatable :: vec(:,:)

 write(*,*)'Dimensionality (1-3)'
 read(*,*)d
 write(*,*)'Box length'
 read(*,*)len
 write(*,*)'Linear number of elements'
 read(*,*)l
 write(*,*)'Lanczos iterations'
 read(*,*)niter
 write(*,*)'State n to write to file'
 read(*,*)st

 allocate(f0(l**d))
 allocate(f1(l**d))
 allocate(f2(l**d))
 allocate(aa(0:niter))
 allocate(nn(0:niter))

 allocate(psi(l**d))
 allocate(eig(0:niter-1))
 allocate(vec(0:niter-1,0:niter-1))

 delta=len/dble(l)
 escale=2.d0*delta**2  ! scale the energy so that the hopping t=1

 call potentialenergy(escale)

 call initstate(psi,l**d) 
 call lanczos1(l**d,niter,psi,eig,vec)
 call lanczos2(l**d,niter,psi,vec(:,st))

 open(10,file='e.dat',status='replace')
 do i=0,niter-1
    write(10,*)i,eig(i)/escale
 enddo
 close(10)

 open(10,file='p.dat',status='replace')
 do i=1,l**d
    write(10,*)psi(i)**2
 enddo
 close(10)

 deallocate(psi)
 deallocate(eig)
 deallocate(vec)
 deallocate(f0)
 deallocate(f1)
 deallocate(f2)
 deallocate(aa)
 deallocate(nn)
 deallocate(vpot)

 end program lanczoscube
!-----------------------!

!---------------------------------------!
 subroutine lanczos1(n,niter,p0,eig,vec)
!--------------------------------------------------------------!
! This subroutine is constructing the normalized basis states
! directly, with the coefficients aa(m) and nn(m) being the
! matrix elements of the tri-diagonal hamiltonian matrix.
!--------------------------------------------------------------!
 use systemparameters; use lanczosvectors; implicit none

 integer :: m,n,niter
 real(8) :: t,eig(0:niter-1),vec(0:niter-1,0:niter-1),p0(n)

 f0=p0
 nn(0)=1.d0
 call hamoperation(l**d,f0,f1)
 aa(0)=dot_product(f0,f1)
 f1=f1-aa(0)*f0
 nn(1)=sqrt(dot_product(f1,f1))
 f1=f1/nn(1)
 do m=2,niter
    write(*,*)'Initial Lanczos iteration: ',m
    call hamoperation(l**d,f1,f2)
    aa(m-1)=dot_product(f1,f2)
    f2=f2-aa(m-1)*f1-nn(m-1)*f0
    nn(m)=sqrt(dot_product(f2,f2))
    f2=f2/nn(m)
    f0=f1
    f1=f2
 enddo
 call diatri(niter,eig,vec)

 end subroutine lanczos1
!-----------------------!

!------------------------------------!
 subroutine lanczos2(n,niter,psi,vec)
!-------------------------------------------------------------------!
! This subroutine re-constructs the normalized Lanczos basis states
! using the coefficients aa(m) and nn(m) previously computed in the 
! subroutine lanczos1(). The vector vec() is an eigenstate in the 
! Lanczos basis and it's tranformed to the real-space basis state 
! stored as the vector psi().
!-------------------------------------------------------------------!
 use systemparameters; use lanczosvectors; implicit none

 integer :: n,m,niter
 real(8) :: psi(n),vec(0:niter-1)
 
 f0=psi
 psi=psi*vec(0)
 call hamoperation(l**d,f0,f1)
 f1=(f1-aa(0)*f0)/nn(1)
 psi=psi+vec(1)*f1
 do m=2,niter-1
    write(*,*)'Second Lanczos iteration: ',m
    call hamoperation(l**d,f1,f2)
    f2=(f2-aa(m-1)*f1-nn(m-1)*f0)/nn(m)
    psi=psi+vec(m)*f2
    f0=f1
    f1=f2
 enddo
 
 end subroutine lanczos2
!-----------------------!

!--------------------------------!
 subroutine hamoperation(n,f1,f2)
!-------------------------------------------------------!
! Acting with H on a state vector f1(), leading to f2().
!-------------------------------------------------------!
 use systemparameters; implicit none

 integer :: j,k,n,x,y,z,ll
 real(8) :: f1(n),f2(n)

 f2=vpot*f1
 if (d==1) then
    do j=1,n
       if (j.ne.1) f2(j)=f2(j)-f1(j-1)
       if (j.ne.n) f2(j)=f2(j)-f1(j+1)
    end do
 else if (d==2) then
    do j=1,l**2
       x=1+mod(j-1,l)
       y=1+(j-1)/l    
       if (x.ne.1) f2(j)=f2(j)-f1(j-1)
       if (x.ne.l) f2(j)=f2(j)-f1(j+1)
       if (y.ne.1) f2(j)=f2(j)-f1(j-l)
       if (y.ne.l) f2(j)=f2(j)-f1(j+l)
    end do
 else if (d==3) then
    ll=l**2
    do j=1,n
       x=1+mod(j-1,l)
       y=1+mod(j-1,ll)/l
       z=1+(j-1)/ll    
       if (x.ne.1) f2(j)=f2(j)-f1(j-1)
       if (x.ne.l) f2(j)=f2(j)-f1(j+1)
       if (y.ne.1) f2(j)=f2(j)-f1(j-l)
       if (y.ne.l) f2(j)=f2(j)-f1(j+l)
       if (z.ne.1) f2(j)=f2(j)-f1(j-ll)
       if (z.ne.l) f2(j)=f2(j)-f1(j+ll)
    end do
 end if

 end subroutine hamoperation
!---------------------------!

!----------------------------!
 subroutine diatri(n,eig,vec)
!----------------------------!
 use lanczosvectors
 implicit none

 integer :: n,inf
 real(8) :: d(n),e(n),eig(n),vec(n,n),work(max(1,2*n-2))

 d=aa(0:n-1)
 e=nn(1:n)
 call dstev('V',n,d,e,vec,n,work,inf)
 eig=d

 end subroutine diatri
!---------------------!

!-------------------------!
 subroutine initstate(f,n)
!-------------------------!
 implicit none

 integer :: i,n
 real(8) :: f(n),norm

 real(8), external :: ran
 
 call initran(1)
 do i=1,n
    f(i)=ran()-0.5d0
 end do
 norm=1.d0/sqrt(dot_product(f,f))
 do i=1,n
    f(i)=f(i)*norm
 end do
 
 end subroutine initstate
!------------------------!

!----------------------------------!
 subroutine potentialenergy(escale)
!-------------------------------------------!
!multiply actual potential energy by escale !
!-------------------------------------------!
 use systemparameters; implicit none

 real(8) :: escale

 allocate(vpot(l**d))
 vpot(:)=0.d0  ! no potential here; this is implemented by user
 vpot(:)=vpot(:)*escale
 vpot(:)=vpot(:)+2.d0*d

 end subroutine potentialenergy
!------------------------------!

!----------------------!
 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!----------------------!
