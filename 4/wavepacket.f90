!-------------------!
 module boxpotential
!-------------------!

 real(8), parameter :: len=40.d0   ! length of the box
 real(8), parameter :: vbeg=2.d0   ! start of the potential 
 real(8), parameter :: vend=6.d0   ! start of the potential 
 real(8), parameter :: vmax=50.d0  ! height of the potential 

 end module boxpotential
!-----------------------!

!--------------------!
 program timevolution
!---------------------------------------------------------------------!
!Split-operator propagation of a Gaussian wave packet in a 1D periodic!
!box with a square potential barrier.                                 !
!---------------------------------------------------------------------!
 use boxpotential
 implicit none

 integer :: i,j,n,nt,fout
 real(8) :: a,dt,dx,k0,potential
 complex(8), allocatable :: psi(:),pot(:),kin(:)

 write(*,*)'Number of lattice points (power of 2 recommended):'
 read(*,*)n
 write(*,*)'Initial wave packet width'
 read(*,*)a
 write(*,*)'Average momentum'
 read(*,*)k0
 write(*,*)'Propagation time step size'
 read(*,*)dt
 write(*,*)'Number of propagation steps'
 read(*,*)nt
 write(*,*)'Psi-output frequency'
 read(*,*)fout

 allocate(psi(n))
 allocate(pot(n))
 allocate(kin(n))

 call initpsi(a,k0,n,len,psi)
 call operators(n,dt,pot,kin)
 call fourier(0,n,psi)
 call fourier(-1,n,psi)

 dx=len/dble(n)
 do i=0,nt
    psi=psi*pot
    call fourier(1,n,psi)
    psi=psi*kin
    call fourier(-1,n,psi)
    if (mod(i,fout)==0) then
       call writepsi(i,n,len,abs(psi))
    endif
 end do

 deallocate(psi)
 deallocate(pot)
 deallocate(kin)

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
!-----------------------------!
 implicit none

 real(8),parameter :: pi=3.141592653589793d0

 integer :: i,n,dir
 real(8) :: nr
 complex(8) :: psi(n)
 real(8), allocatable, save :: wsave(:)       

 if (dir==1) then
    call dcfftf(n,psi,wsave)
    nr=1.d0/dble(n)
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
 real(8) :: x,k,b,dt,dx,dk
 complex(8) :: pot(n),kin(n)

 dx=len/dble(n)
 dk=2.d0*pi/len
 do i=-n/2+1,n/2
    x=dble(i)*dx
    k=dble(i)*dk
    if (i>0) then
       j=i
    else     
       j=i+n
    endif
    if (x>=vbeg.and.x<=vend) then
       pot(j)=exp(-dt*(0,1)*vmax)          
    else
       pot(j)=1.d0
    end if
    kin(j)=exp(-dt*0.5d0*(0,1)*k**2)
 end do

 end subroutine operators
!------------------------!

!------------------------------------!
 subroutine writepsi(step,n,len,psi2)
!------------------------------------!
 implicit none

 integer :: i,j,n,i1,i2,i3,i4,step
 real(8) :: dx,len,psi2(n)
 character(11) :: fname

 dx=len/dble(n)
 psi2=(psi2*dx)**2

 i1=mod(step,10)
 i2=mod(step,100)/10
 i3=mod(step,1000)/100
 i4=step/1000
 fname='psi0000.dat'
 fname(4:4)=achar(48+i4)
 fname(5:5)=achar(48+i3)
 fname(6:6)=achar(48+i2)
 fname(7:7)=achar(48+i1)

 open(1,file=fname)
 do i=-n/2+1,n/2          
    if (i>0) then
       j=i
    else     
       j=i+n
    endif
    write(1,*)dble(i)*dx,psi2(j)
 end do
 close(1)
 write(*,*)'wrote step ',step

 end subroutine writepsi
!-----------------------!
