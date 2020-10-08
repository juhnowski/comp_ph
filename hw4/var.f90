!---------------------!
 program variational2d
!---------------------!
 implicit none

 real(8), parameter :: lx=5.d0
 real(8), parameter :: ly=10.d0
 real(8), parameter :: beta=13.11d0

 integer :: i,j,nx,ny,np
 real(8) :: vv,x1,x2,y1,y2,psi(5),psirealspace
 real(8), allocatable :: ham(:,:)
 real(8), allocatable :: eig(:)

 open(1,file='vparam.in',status='old')
 read(1,*)vv
 read(1,*)x1,x2,y1,y2
 read(1,*)nx,ny
 close(1)

 allocate(ham(nx*ny,nx*ny))
 allocate(eig(nx*ny))

 call makeh(nx,ny,x1,x2,y1,y2,lx,ly,vv,ham)
 call diasym(ham,eig,nx*ny)
 
 open(1,file='eig.dat',status='replace')
 do i=1,nx*ny
    write(1,*)i,eig(i)/beta
 end do
 close(1)

 open(1,file='psi0.dat',status='replace')
 open(2,file='psi1.dat',status='replace')
 open(3,file='psi2.dat',status='replace')
 open(4,file='psi3.dat',status='replace')
 do j=0,200
    y1=ly*dble(j)/200.d0
    do i=0,100
       x1=lx*dble(i)/100.d0
       write(1,*)psirealspace(x1,y1,nx,ny,lx,ly,ham(:,1))**2
       write(2,*)psirealspace(x1,y1,nx,ny,lx,ly,ham(:,2))**2
       write(3,*)psirealspace(x1,y1,nx,ny,lx,ly,ham(:,3))**2
       write(4,*)psirealspace(x1,y1,nx,ny,lx,ly,ham(:,4))**2
    end do
 end do
 close(1)

 end program variational2D
!-------------------------!

!------------------------------------------------!
 subroutine makeh(nx,ny,x1,x2,y1,y2,lx,ly,vv,ham)
!------------------------------------------------!
 implicit none

 real(8), parameter :: beta=13.11d0

 integer :: i,j,nx,ny,ix,iy,jx,jy
 real(8) :: x1,x2,y1,y2,lx,ly,vv,xint,yint1,yint2,integ1
 real(8) :: ham(nx*ny,nx*ny),energ0

 do j=1,nx*ny
    jx=1+mod(j-1,nx)
    jy=1+(j-1)/nx
    xint=integ1(jx,jx,x1,x2,lx)
    yint1=integ1(jy,jy,y1,y2,ly)
    yint2=integ1(jy,jy,ly-y2,ly-y1,ly)
    ham(j,j)=energ0(jx,jy,lx,ly)+vv*beta*xint*(yint1+yint2)
    do i=1,j-1
       ix=1+mod(i-1,nx)
       iy=1+(i-1)/nx
       xint=integ1(ix,jx,x1,x2,lx)
       yint1=integ1(iy,jy,y1,y2,ly)
       yint2=integ1(iy,jy,ly-y2,ly-y1,ly)
       ham(i,j)=vv*beta*xint*(yint1+yint2)
       ham(j,i)=ham(i,j)
    end do
 end do

 end subroutine makeh
!--------------------!

!--------------------------------------------------!
 real(8) function psirealspace(x,y,nx,ny,lx,ly,vec)
!--------------------------------------------------!
 implicit none

 real(8), parameter :: pi=3.14159265358979328d0

 integer :: i,kx,ky,nx,ny
 real(8) :: x,y,lx,ly,psi,phix,phiy,vec(nx*ny)

 psi=0.d0
 do i=1,nx*ny
    kx=1+mod(i-1,nx)
    ky=1+(i-1)/nx
    phix=sin(dble(kx)*pi*x/lx)
    phiy=sin(dble(ky)*pi*y/ly)
    psi=psi+vec(i)*phix*phiy
 end do
 psirealspace=psi/sqrt(lx*ly)
 
 end function psirealspace
!-------------------------!

!------------------------------------!
 real(8) function energ0(kx,ky,lx,ly)
!------------------------------------!
 implicit none

 real(8), parameter :: pi2=(3.14159265358979328d0**2)/2.d0

 integer :: kx,ky
 real(8) :: lx,ly

 energ0=pi2*((dble(kx)/lx)**2+(dble(ky)/ly)**2)

 end function energ0
!-------------------!

!---------------------------------------!
 real(8) function integ1(k1,k2,x1,x2,lx)
!---------------------------------------!
 implicit none

 real(8), parameter :: pi=3.14159265358979328d0

 integer :: k1,k2
 real(8) :: x1,x2,lx,pilx,i1,i2

 pilx=pi/lx
 if (k1==k2) then
    i1=x1/2.d0-sin(dble(2*k1)*pilx*x1)/(dble(4*k1)*pilx)
    i2=x2/2.d0-sin(dble(2*k1)*pilx*x2)/(dble(4*k1)*pilx)
 else
    i1=sin(dble(k1-k2)*pilx*x1)/(2.d0*dble(k1-k2)*pilx)
    i1=i1-sin(dble(k1+k2)*pilx*x1)/(2.d0*dble(k1+k2)*pilx)
    i2=sin(dble(k1-k2)*pilx*x2)/(2.d0*dble(k1-k2)*pilx)
    i2=i2-sin(dble(k1+k2)*pilx*x2)/(2.d0*dble(k1+k2)*pilx)
 endif
 integ1=2.d0*(i2-i1)/lx

 end function integ1
!-------------------!

!--------------------------!
 subroutine diasym(a,eig,n)
!--------------------------!
 implicit none

 integer :: n,l,inf
 real(8) :: a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!----------------------!


