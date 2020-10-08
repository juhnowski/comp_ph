!----------------------!
 module systemvariables
!----------------------!
 implicit none

 integer :: d
 integer :: n
 integer :: nbors
 real(8) :: bprob

 integer, allocatable :: l(:)
 integer, allocatable :: neighbor(:,:)
 integer, allocatable :: bondspin(:,:)
 integer, allocatable :: spinbond(:,:)
 logical, allocatable :: notvisited(:)
 integer, allocatable :: stack(:)
 integer, allocatable :: spin(:)
 logical, allocatable :: bond(:)

 end module systemvariables
!--------------------------!

!-------------------------------------------------------------------!
!Monte carlo simulation of the d-simensional (d=1,2,3) Ising model  !
!based on the Swendesen-Wang cluster algorithm. An arbitrary number !
!of temperature can be studied. Input is read from 'ising.in';      !
!   d          = dimensionality (1,2,3)                             !
!   lx,..,lz   = lengths in the x,y,z directions                    !
!   nt         = number of temperatures                             !
!   t_i,bins_i,steps_i = T/J, number of bins, MC step per bin       !
!   (one line for each temperarure)                                 !
!Random number seed (4 integers on separate lines) from 'seed.in'   !
!Output is given in 'bindata.dat', which should be processed by     !
!the program 'multiaverage.f90'. The quantities calculated are:     !
!  - the energy                                                     !
!  - the energy squared                                             !
!  - absolute value of the magnetization                            !
!  - squared magnetization (using standard and improved estimators) !
!-------------------------------------------------------------------!
 program isingsw
!---------------!
 use systemvariables
 implicit none

 integer :: i,j,it,nt,bins,steps
 real(8) :: mag2
 real(8), allocatable :: temp(:)

 open(unit=1,file='ising.in',status='old')
read(1,*)d; allocate(l(d))
 read(1,*)l
 read(1,*)bins,steps
 read(1,*)nt
 allocate(temp(nt))
 do i=1,nt
    read(1,*)temp(i)
 end do
 close(1)

 call initran(1)
 call initarrays
 call lattice
 do it=1,nt
    bprob=1.d0-exp(-2.d0/temp(it))
    do i=1,steps/4
       call castbonds
       call flipclusters(mag2)
    end do
    do j=1,bins
       call resetbindata
       do i=1,steps
          call castbonds
          call flipclusters(mag2)
          call measure(mag2)
       end do
       call writebindata(it,n,steps)
    end do
 end do

 deallocate(temp)
 call cleanup

 end program isingsw
!-------------------!

!--------------------!
 subroutine castbonds
!--------------------!
 use systemvariables
 implicit none

 integer :: b
 real(8),external :: ran

 do b=0,d*n-1
    if (spin(bondspin(1,b))==spin(bondspin(2,b))) then
       if (ran()<bprob) then
          bond(b)=.true.    
       else
          bond(b)=.false.    
       end if
    else
       bond(b)=.false.    
    end if
 end do

 end subroutine castbonds
!------------------------!

!-----------------------------!
 subroutine flipclusters(mag2)
!-----------------------------!
 use systemvariables
 implicit none

 integer :: i,s0,s1,cseed,csize,nstack
 logical :: flipclus
 real(8) :: ran,mag2
 external :: ran

 notvisited=.true.
 mag2=0.d0 
 cseed=0

 1 if (ran()<0.5d0) then
    flipclus=.true.
 else
    flipclus=.false.
 end if
 csize=0

 notvisited(cseed)=.false.
 if (flipclus) spin(cseed)=-spin(cseed)
 nstack=1
 stack(1)=cseed

 do
    if (nstack==0) exit
    s0=stack(nstack)
    nstack=nstack-1
    csize=csize+1
    do i=1,nbors
       s1=neighbor(i,s0)
       if (bond(spinbond(i,s0)).and.notvisited(s1)) then          
          notvisited(s1)=.false.
          if (flipclus) spin(s1)=-spin(s1)
          nstack=nstack+1
          stack(nstack)=s1
       end if    
    end do
 end do
 mag2=mag2+dble(csize)**2
  
 do i=cseed+1,n-1
    if (notvisited(i)) then
       cseed=i
       goto 1
    end if
 end do

 end subroutine flipclusters
!---------------------------!

!--------------------------------------------------------!
!Measures the energy (enrg1), the squared energy (enrg2),! 
!the absolute value of the magnetization (magn1), and the!
!squared magnetization (magn2). An improved-estimator    !
!value (m2clus) obtained in the flipclusters routine is  !
!passed as an input argument.                            !
!--------------------------------------------------------!
 subroutine measure(m2clus)
!--------------------------!
 use systemvariables
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,magn2c
 common/measurments/enrg1,enrg2,magn1,magn2,magn2c

 integer :: i,s,e,m,ss
 real(8) :: m2clus

 e=0
 m=0
 do s=0,n-1
    ss=0
    do i=1,d
       ss=ss+spin(neighbor(i,s))
    end do
    e=e-spin(s)*ss
    m=m+spin(s)
 end do
 enrg1=enrg1+dble(e)
 enrg2=enrg2+dble(e)**2
 magn1=magn1+abs(m)
 magn2=magn2+abs(m)**2
 magn2c=magn2c+m2clus

 end subroutine measure
!----------------------!

!-------------------------------------------------!
!Writes accumulated data to the file 'bindata.dat'!
!-------------------------------------------------!
 subroutine writebindata(it,n,steps)
!-----------------------------------!
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,magn2c
 common/measurments/enrg1,enrg2,magn1,magn2,magn2c

 integer :: n,it,steps
 
 open(1,file='bindata.dat',status='unknown',position='append')
 enrg1=enrg1/(dble(steps)*dble(n))
 enrg2=enrg2/(dble(steps)*dble(n)**2)
 magn1=magn1/(dble(steps)*dble(n))
 magn2=magn2/(dble(steps)*dble(n)**2)
 magn2c=magn2c/(dble(steps)*dble(n)**2)
 write(1,1)it,enrg1,enrg2,magn1,magn2,magn2c
 1 format(i3,' ',5f15.10)
 close(1)

 end subroutine writebindata
!---------------------------!

!--------------------------------------------!
!Sets the data accumulation variables to zero!
!--------------------------------------------!
 subroutine resetbindata
!-----------------------!
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,magn2c
 common/measurments/enrg1,enrg2,magn1,magn2,magn2c

 enrg1=0.d0
 enrg2=0.d0
 magn1=0.d0
 magn2=0.d0
 magn2c=0.d0

 end subroutine resetbindata
!---------------------------!

!--------------------------------------------------------------!
!Constructs lattice tables for a simple cubic lattice (d=1,2,3)! 
!neighbor(i,s) = neighbor i of site s (i=1,...,2d)             !
!bondspin(i,b) = spin i connected to bond b (i=1,...,2d)       !
!spinbond(i,s) = bond i connected to site s (i=1,...,2d)       !
!--------------------------------------------------------------!
 subroutine lattice
!------------------!
 use systemvariables
 implicit none

 integer :: s0,s1,s2,s3,s4,s5,s6,x0,x1,x2,y0,y1,y2,z0,z1,z2

 select case(d)
 case (1)
    do x0=0,n-1
       x1=mod(x0+1,l(1))
       x2=mod(x0-1+l(1),l(1))
       neighbor(1,s0)=x1
       neighbor(2,s0)=x2
       bondspin(1,x0)=x0
       bondspin(2,x0)=x1
       spinbond(1,x0)=x0
       spinbond(2,x1)=x0
    end do
 case (2)
    do s0=0,n-1
       x0=mod(s0,l(1))
       y0=s0/l(1)
       x1=mod(x0+1,l(1))
       x2=mod(x0-1+l(1),l(1))
       y1=mod(y0+1,l(2))
       y2=mod(y0-1+l(2),l(2))
       s1=x1+y0*l(1)
       s2=x0+y1*l(1)
       s3=x2+y0*l(1)
       s4=x0+y2*l(1)
       neighbor(1,s0)=s1
       neighbor(2,s0)=s2
       neighbor(3,s0)=s3
       neighbor(4,s0)=s4
       bondspin(1,2*s0)=s0
       bondspin(2,2*s0)=s1
       bondspin(1,2*s0+1)=s0
       bondspin(2,2*s0+1)=s2
       spinbond(1,s0)=2*s0
       spinbond(2,s0)=2*s0+1
       spinbond(3,s1)=2*s0
       spinbond(4,s2)=2*s0+1
    end do
 case (3)
    do s0=0,n-1
       x0=mod(s0,l(1))
       y0=mod(s0,l(1)*l(2))/l(1)
       z0=s0/(l(1)*l(2))
       x1=mod(x0+1,l(1))
       x2=mod(x0-1+l(1),l(1))
       y1=mod(y0+1,l(2))
       y2=mod(y0-1+l(2),l(2))
       z1=mod(z0+1,l(3))
       z2=mod(z0-1+l(3),l(3))
       s1=x1+y0*l(1)+z0*l(1)*l(2)
       s2=x0+y1*l(1)+z0*l(1)*l(2)
       s3=x0+y0*l(1)+z1*l(1)*l(2)
       s4=x2+y0*l(1)+z0*l(1)*l(2)
       s5=x0+y2*l(1)+z0*l(1)*l(2)
       s6=x0+y0*l(1)+z2*l(1)*l(2)
       neighbor(1,s0)=s1
       neighbor(2,s0)=s2
       neighbor(3,s0)=s3
       neighbor(4,s0)=s4
       neighbor(5,s0)=s5
       neighbor(6,s0)=s6
       bondspin(1,3*s0)=s0
       bondspin(2,3*s0)=s1
       bondspin(1,3*s0+1)=s0
       bondspin(2,3*s0+1)=s2
       bondspin(1,3*s0+2)=s0
       bondspin(2,3*s0+2)=s3
       spinbond(1,s0)=3*s0
       spinbond(2,s0)=3*s0+1
       spinbond(3,s0)=3*s0+2
       spinbond(4,s1)=3*s0
       spinbond(5,s2)=3*s0+1
       spinbond(6,s3)=3*s0+2
    end do
 end select

 end subroutine lattice
!----------------------!

!---------------------!
 subroutine initarrays
!---------------------!
 use systemvariables
 implicit none

 integer :: i
 real(8), external :: ran

 nbors=2*d
 n=product(l)
 allocate(neighbor(nbors,0:n-1))
 allocate(bondspin(2,0:d*n-1))
 allocate(spinbond(nbors,0:n-1))
 allocate(notvisited(0:n-1))
 allocate(stack(n))
 allocate(spin(0:n-1))
 allocate(bond(0:d*n-1))
 do i=0,n-1
    spin(i)=2*int(ran()*2.d0)-1
 end do

 end subroutine initarrays
!-------------------------!

!------------------!
 subroutine cleanup
!------------------!
 use systemvariables
 implicit none

 deallocate(l)
 deallocate(neighbor)
 deallocate(bondspin)
 deallocate(spinbond)
 deallocate(notvisited)
 deallocate(stack)
 deallocate(spin)
 deallocate(bond)

 end subroutine cleanup
!----------------------!
 
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
