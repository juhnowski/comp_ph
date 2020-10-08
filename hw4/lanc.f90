!-----------------------!
 module systemparameters
!-----------------------!
 
 real(8), parameter :: beta=13.11d0
 integer :: nx
 integer :: ny
 real(8) :: thop
 real(8) :: delta
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

!--------------------------------------------------------------!
!Reads from the file 'lparam.in:                               !
! vv                                                           !
! x1,x2,y1,y2                                                  !
! lx,ly,delta,niter                                            !
! nplot,st                                                     !
!--------------------------------------------------------------!
! vv    : potential in the two patches                         !
! x1,x2 : patches between these x-values                       !
! y1,y2 : patches between these y-values                       !
! nx,ny : number of volume elements in the x and y directions  !
! delta : discretization (same for x and y)                    !
! nplot : number of state to be graphed                        !
! st    : st=0: start from random state                        !
!       : st=1; start from saved state for the same size       !
!       : st=2; start from saved state for half the size       !
!--------------------------------------------------------------!
!Writes the following files:                                   !
!eig.dat   : Eigenvalues (energies)                            !
!plsi0.dat : Linear combination of the nplot first states      !
!            (to be used as initial state in subsequent runs)  !
!psin.ps   : n=0,...,nplot-1; squared wave function graphs     !
!--------------------------------------------------------------!

!-----------------!
 program lanczos2d
!-----------------!
 use systemparameters; use lanczosvectors; implicit none

 integer :: i,niter,step,nplot
 real(8) :: vv,lx,ly,x1,x2,y1,y2,escale
 real(8), allocatable :: psi(:)
 real(8), allocatable :: psi0(:)
 real(8), allocatable :: psisave(:)
 real(8), allocatable :: eig(:)
 real(8), allocatable :: vec(:,:)

 open(1,file='lparam.in',status='old')
 read(1,*)vv
 read(1,*)lx,ly
 read(1,*)x1,x2,y1,y2
 read(1,*)delta,niter
 read(1,*)nplot
 close(1)

 nx=int(lx/delta+0.1d0)
 ny=int(ly/delta+0.1d0)

 call allocatelanczosvectors(nx*ny,niter)
 allocate(psi(nx*ny))
 allocate(psi0(nx*ny))
 allocate(psisave(nx*ny))
 allocate(eig(0:niter-1))
 allocate(vec(0:niter-1,0:niter-1))
 call potentialenergy(vv,x1,x2,y1,y2,nx*ny)
 escale=0.5d0/delta**2
 vpot(:)=vpot(:)*beta/escale
 thop=1.d0

 call initstate(psi,nx*ny)
 call lanczos1(nx*ny,niter,psi,eig,vec)
 psi0(:)=psi(:)

 open(1,file='eig.dat',status='replace')
 do i=0,niter-1
    write(1,*)i,(escale*eig(i)+4.d0*escale)/beta
 end do
 close(1) 

 psisave(:)=0.d0
 do i=1,nplot
    psi(:)=psi0(:)
    call lanczos2(i-1,nx*ny,niter,psi,vec(:,i-1))
    call initgraph(i-1,300.,nx,ny)
    call plotstate(psi,nx,ny)
    call closegraph()
    psisave(:)=psisave(:)+psi(:)
 enddo
 call writestate(psisave,nx*ny)

 call cleanup
 deallocate(psi)
 deallocate(eig)
 deallocate(vec)

 end program lanczos2d
!---------------------!

!---------------------------------------!
 subroutine lanczos1(n,niter,p0,eig,vec)
!---------------------------------------!
 use systemparameters; use lanczosvectors; implicit none

 integer :: m,n,niter
 real(8) :: t,eig(0:niter-1),vec(0:niter-1,0:niter-1),p0(n)

 f0=p0
 nn(0)=1.d0
 call hamoperation(nx*ny,f0,f1)
 aa(0)=dot_product(f0,f1)
 f1=f1-aa(0)*f0
 nn(1)=sqrt(dot_product(f1,f1))
 f1=f1/nn(1)
 do m=2,niter
    write(*,*)'Initial Lanczos iteration: ',m
    call hamoperation(nx*ny,f1,f2)
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

!---------------------------------------!
 subroutine lanczos2(st,n,niter,psi,vec)
!---------------------------------------!
 use systemparameters; use lanczosvectors; implicit none

 integer :: n,m,st,niter
 real(8) :: psi(n),vec(0:niter-1)
 
 f0=psi
 psi=psi*vec(0)
 call hamoperation(nx*ny,f0,f1)
 f1=(f1-aa(0)*f0)/nn(1)
 psi=psi+vec(1)*f1
 do m=2,niter-1
    write(*,*)st,'  Second Lanczos iteration: ',m
    call hamoperation(nx*ny,f1,f2)
    f2=(f2-aa(m-1)*f1-nn(m-1)*f0)/nn(m)
    psi=psi+vec(m)*f2
    f0=f1
    f1=f2
 enddo
 
 end subroutine lanczos2
!-----------------------!
 
!--------------------------------!
 subroutine hamoperation(n,f1,f2)
!--------------------------------!
 use systemparameters
 implicit none

 integer :: j,k,n,x,y
 real(8) :: f1(n),f2(n)

 f2=vpot*f1
 do j=1,n
    x=1+mod(j-1,nx)
    y=1+(j-1)/nx    
    if (x.ne.1) f2(j)=f2(j)-thop*f1(j-1)
    if (x.ne.nx) f2(j)=f2(j)-thop*f1(j+1)
    if (y.ne.1) f2(j)=f2(j)-thop*f1(j-nx)
    if (y.ne.ny) f2(j)=f2(j)-thop*f1(j+nx)
 end do
 
 end subroutine hamoperation
!---------------------------!

!----------------------------!
 subroutine diatri(n,eig,vec)
!----------------------------!
 use lanczosvectors; implicit none

 integer :: n,inf
 real(8) :: d(n),e(n),eig(n),vec(n,n),work(max(1,2*n-2))

 d=aa(0:n-1)
 e=nn(1:n)
 call dstev('V',n,d,e,vec,n,work,inf)
 eig=d

 end subroutine diatri
!---------------------!

!--------------------------------------------!
 subroutine potentialenergy(vv,x1,x2,y1,y2,n)
!--------------------------------------------!
!Constructs the potential-energy array       !
!--------------------------------------------!
 use systemparameters
 implicit none

 integer :: i,j,k,n
 real(8) :: vv,x,y,x1,x2,y1,y2,lly

 allocate(vpot(n))
 lly=dble(ny)*delta
 k=0
 do j=0,ny-1
    y=delta*(dble(j)+0.5d0)
    do i=0,nx-1
       k=k+1
       x=delta*(dble(i)+0.5d0)
       if (x>=x1.and.x<=x2.and.y>=y1.and.y<=y2) then
          vpot(k)=vv
       else if (x>=x1.and.x<=x2.and.y>=lly-y2.and.y<=lly-y1) then
          vpot(k)=vv
       else
          vpot(k)=0.d0
       end if
    end do
 end do
       
 end subroutine potentialenergy
!------------------------------!

!-------------------------!
 subroutine initstate(f,n)
!---------------------------------!
!Constructs a random initial state!
!---------------------------------!
 implicit none

 integer :: i,n
 real(8) :: f(n),rand,norm

 call initrand(0)
 do i=1,n
    f(i)=rand()-0.5d0
 end do
 norm=1.d0/sqrt(dot_product(f,f))
 f(:)=f(:)*norm
 
 end subroutine initstate
!------------------------!

!--------------------------!
 subroutine writestate(f,n)
!------------------------------------------------!
!Writes the linear combination f of all graphed  ! 
!states to the file psi0.dat (for use as starting!
!state in subsequent rruns)                      !
!------------------------------------------------!
 implicit none

 integer :: i,n
 real(8) :: f(n),norm

 norm=1.d0/sqrt(dot_product(f,f))
 f(:)=f(:)*norm
 open(10,file='psi0.dat',status='replace')
 do i=1,n
    write(10,*)f(i)
 enddo
 close(10)

 end subroutine writestate
!-------------------------!

!-----------------------------!
 subroutine plotstate(a,nx,ny)
!----------------------------------------------------!
!Writes the postscript file for the state stored in a!
!----------------------------------------------------!
 implicit none
 
 integer :: i,j,nx,ny
 real(8) :: a(nx*ny),amax,amin,aa,ar,ag,ab

 a=a**2
 amax=maxval(a)
 amin=minval(a)
 do j=1,ny
    do i=1,nx
       aa=(a(i+(j-1)*nx)-amin)/(amax-amin)
       call setcolor(aa,ar,ag,ab)
       write(2,1)i,j,ar,ag,ab,' p'
       1 format(2i5,' ',3f5.2,a)
    end do
 end do

 end subroutine plotstate
!------------------------!

!---------------------------------!
 subroutine initgraph(st,lx,nx,ny)
!-------------------------------------------------------!
!Initializes a postscript file called psin.ps, where    !
!n=0,1,2,..., is the state label (0 = ground state, etc)! 
!-------------------------------------------------------!
 implicit none

 integer :: nx,ny,st
 real    :: lx,p
 character(10) :: fname

 fname='psi0.ps'
 fname(4:4)=achar(48+st)
 open(2,file=fname,status='replace') 
 p=100.d0/real(nx) 
 1 format(a)
 write(2,1)'%!'
 write(2,*)'100 200 translate'
 write(2,*)lx/nx,lx/nx,' scale'
 write(2,*)'0 0 moveto ',nx+1,' 0 rlineto 0 ',ny+1,' rlineto'
 write(2,*)-nx-1,' 0 rlineto closepath fill'
 write(2,*)'1.1 setlinewidth -0.55 -0.05 translate' 
 write(2,*)'/p {setrgbcolor moveto 1.1 0 rlineto stroke} def'
 
 end subroutine initgraph
!------------------------! 

!-----------------------! 
 subroutine closegraph()
!-----------------------! 

 write(2,*)'showpage'
 close(2)

 end subroutine closegraph
!-------------------------! 

!-------------------------------!
 subroutine setcolor(a,rc,gc,bc)
!-------------------------------!
 implicit none

 real(8), parameter :: f1=1.d0/4.d0
 real(8), parameter :: f2=2.d0/4.d0
 real(8), parameter :: f3=3.d0/4.d0
 
 real(8) :: a,r,g,b,rc,bc,gc,aa

 if (a > 1.d0) then
    rc=1.d0
    gc=1.d0
    bc=1.d0
 else if (a < f1) then
    r=0.d0
    g=0.d0
    b=0.75d0*a/f1
 else if (a < f2) then
    r=((a-f1)/f1)**0.5d0
    g=0.d0
    b=0.75d0*(1-(a-f1)/f1)
 else if (a < f3) then
    r=1.d0
    g=(a-f2)/f1
    b=0.d0
 else
    r=1.d0
    g=1.d0
    b=(a-f3)/f1
 end if
 aa=a**0.2d0
 rc=r*aa
 gc=g*aa
 bc=b*aa

 end subroutine setcolor
!-----------------------! 

!--------------------------------------------!
 subroutine allocatelanczosvectors(nst,nlanc)
!--------------------------------------------!
 use lanczosvectors
 implicit none

 integer :: nst,nlanc

 allocate(f0(nst))
 allocate(f1(nst))
 allocate(f2(nst))
 allocate(aa(0:nlanc))
 allocate(nn(0:nlanc))

 end subroutine allocatelanczosvectors
!-------------------------------------!

!------------------!
 subroutine cleanup
 use systemparameters
 use lanczosvectors

 deallocate(f0)
 deallocate(f1)
 deallocate(f2)
 deallocate(aa)
 deallocate(nn)
 deallocate(vpot)

 end subroutine cleanup
!----------------------!

!-----------------------!
 real(8) function rand()
!-----------------------!
 implicit none

 integer :: iir,jjr,kkr,nnr,mzran
 common/brand/iir,jjr,kkr,nnr,mzran

 mzran=iir-kkr
 if (mzran < 0) mzran=mzran+2147483579
 iir=jjr; jjr=kkr; kkr=mzran
 nnr=69069*nnr+1013904243
 mzran=mzran+nnr
 rand=0.5d0+mzran*0.23283064d-9

 end function rand
!-----------------!

!----------------------!
 subroutine initrand(w)
!----------------------!
 implicit none

 integer :: iir,jjr,kkr,nnr,mzran
 common/brand/iir,jjr,kkr,nnr,mzran

 integer :: w
 real(8) :: rand

 open(1,file='seed.in',status='old')
 read(1,*)iir 
 read(1,*)jjr
 read(1,*)kkr
 read(1,*)nnr
 close(1)
 iir=abs(iir)+1; jjr=abs(jjr)+1; kkr=abs(kkr)+1
 if (w /= 0) then
   open(1,file='seed.in',status='replace')
   write(1,*)abs(int((rand()-0.5d0)/0.23283064d-9))
   write(1,*)abs(int((rand()-0.5d0)/0.23283064d-9))
   write(1,*)abs(int((rand()-0.5d0)/0.23283064d-9))
   write(1,*)abs(int((rand()-0.5d0)/0.23283064d-9))
   close(1)
 end if

 end subroutine initrand
!-----------------------!
