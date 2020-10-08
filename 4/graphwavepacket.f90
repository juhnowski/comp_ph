!-------------------!
 module boxpotential
!-------------------!

 real(8), parameter :: len=40.d0   ! length of the box
 real(8), parameter :: vbeg=3.d0   ! start of the potential 
 real(8), parameter :: vend=13.d0   ! start of the potential 
 real(8), parameter :: vmax=10.d0  ! slope of the potential 

 end module boxpotential
!-----------------------!

!--------------------!
 program timevolution
!--------------------!
 use boxpotential
 implicit none

 integer :: i,j,n,nt,gf
 real(8) :: a,dt,dx,k0,potential
 real(8), allocatable :: plot(:)
 complex(8), allocatable :: psi(:),pot(:),kin(:)


 write(*,*)'Number of lattice points (power of 2 recommended):'
 read(*,*)n
 write(*,*)'Initial wave packet width'
 read(*,*)a
 write(*,*)'Average momentum'
 read(*,*)k0
 write(*,*)'Propagation time step'
 read(*,*)dt
 write(*,*)'Number of propagation steps'
 read(*,*)nt
 write(*,*)'graph frequency'
 read(*,*)gf

 allocate(psi(n))
 allocate(pot(n))
 allocate(kin(n))
 allocate(plot(n))
 dx=len/dble(n)

 call initpsi(a,k0,n,len,psi)
 call operators(n,dt,pot,kin)
 call fourier(0,n,psi)
 call fourier(-1,n,psi)

 do i=1,nt
    write(*,*)'step ',i
    psi=psi*pot
    call fourier(1,n,psi)
    psi=psi*kin
    call fourier(-1,n,psi)
    if (mod(i,gf)==1) then
       call initgraph(i/gf,(i-1)*dt)
       plot=abs(dx*psi)**2
       call graphpsi(dx,n,plot)
    endif
 end do

 deallocate(psi)
 deallocate(pot)
 deallocate(kin)
 deallocate(plot)

 end program timevolution
!------------------------!

!----------------------------------!
 subroutine initpsi(a,k0,n,len,psi)
!----------------------------------!
 implicit none

 real(8),parameter :: pi=3.141592653589793d0

 integer :: i,j,n
 real(8) :: a,k,k0,dk,nr,len
 complex(8) :: psi(n)
 
 dk=2.d0*pi/len
 do i=-n/2+1,n/2
    k=dble(i)*dk
    if (i>0) then
       j=i
    else     
       j=i+n
    endif
    psi(j)=exp(-(a*(k-k0)/2.d0)**2)
    nr=nr+abs(psi(j))**2
 end do
 nr=1.d0/sqrt(nr*dk)
 do i=1,n
   psi(i)=psi(i)*nr
 end do

 end subroutine initpsi
!----------------------!

!----------------------------!
 subroutine fourier(dir,n,psi)
!----------------------------!
 implicit none

 real(8),parameter :: pi=3.141592653589793d0

 integer :: i,n,dir
 real(8) :: nr
 complex(8) :: psi(n)
 real(8), allocatable, save :: wsave(:)       


 if (dir==1) then
    call dcfftf(n,psi,wsave)
    nr=1.d0/DFLOAT(n)
    do i=1,n
       psi(i)=psi(i)*nr
    end do
 elseif (dir==-1) then
    call dcfftb(n,psi,wsave)
 elseif (dir==0) then
    allocate(wsave(4*n+20))
    call dcffti(n,wsave)
 endif

 end subroutine fourier
!----------------------!

!--------------------------------------!
 subroutine operators(n,dt,pot,kin)
!--------------------------------------!
 use boxpotential
 implicit none

 real(8),parameter :: pi=3.141592653589793d0

 integer :: i,j,n
 real(8) :: x,k,b,dt,dx,dk,vmid,vslo
 complex(8) :: pot(n),kin(n)

 dx=len/dble(n)
 dk=2.d0*pi/len
 vmid=(vbeg+vend)/2.d0
 vslo=vmax/(vmid-vbeg)
 do i=-n/2+1,n/2
    x=dble(i)*dx
    k=dble(i)*dk
    if (i>0) then
       j=i
    else     
       j=i+n
    endif
    if (x>=vbeg.and.x<=vmid) then
       pot(j)=exp(-dt*(0,1)*(x-vbeg)*vslo*vmax)          
    elseif (x>vmid.and.x<=vend) then
       pot(j)=exp(-dt*(0,1)*(vmax-(x-vmid))*vslo*vmax)          
    else
       pot(j)=1.d0
    end if
    kin(j)=exp(-dt*0.5d0*(0,1)*k**2)
 end do

 end subroutine operators
!------------------------!

!-------------------------------!
 subroutine graphpsi(dx,n,psi)
!-------------------------------!
 implicit none

 integer :: i,j,n
 real(8) :: x,dx,psi(n)

 psi=min(700.d0*psi,300.d0)
 write(1,*)'gsave 1 0 0 setrgbcolor 2 setlinewidth'
 do i=-n/2+1,n/2
    if (i>0) then
       j=i
    else     
       j=i+n
    endif
    x=dble(i)*dx
    if (i/=-n/2+1) then
       write(1,10)10.d0*x,psi(j)
       10 format(f7.2,' ',f7.2,' l')
    else
       write(1,20)10.d0*x,psi(j)
       20 format(f7.2,' ',f7.2,' m')
    endif
 end do
 write(1,*)'stroke grestore'
 close(1)

 end subroutine graphpsi
!-----------------------!

!-----------------------------!
 subroutine initgraph(frame,t)
!-----------------------------!
 implicit none

 integer :: i1,i2,i3,frame
 real(8) :: t
 character(9) :: fname

 fname='psi000.ps'
 i1=mod(frame,10)
 i2=mod(frame,100)/10
 i3=frame/100
 fname(4:4)=achar(48+i3)
 fname(5:5)=achar(48+i2)
 fname(6:6)=achar(48+i1)
 open(1,file=fname)
 write(1,*)'%!'
 write(1,*)'/m {moveto} def'
 write(1,*)'/l {lineto} def'
 write(1,*)'/TimesRoman findfont 24 scalefont setfont'
 write(1,*)'300 200 translate'
 write(1,*)'-200 0 moveto 400 0 rlineto 0 300 rlineto -400 0 rlineto'
 write(1,*)'closepath stroke'
 write(1,*)'gsave 0 0 1 setrgbcolor'
 write(1,*)'30 0 moveto 50 100 rlineto 50 -100 rlineto' 
 write(1,*)'stroke grestore'
 write(1,*)'-180 268 moveto'
 write(1,10)t
 10 format ('(t=',f5.3,') show')
 
 end subroutine initgraph 
!------------------------!

