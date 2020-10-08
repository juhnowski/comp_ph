!--------------!
 module system
!--------------!
 save

 integer :: nn
 integer :: nb
 integer :: zz

 integer, allocatable :: site(:,:)
 integer, allocatable :: hamx(:)
 integer, allocatable :: nxop(:)
 real(8), allocatable :: hamz(:)
 real(8), allocatable :: jbnd(:)

 complex(8), allocatable :: ff(:)
 complex(8), allocatable :: f0(:)
 complex(8), allocatable :: f1(:)
 complex(8), allocatable :: f2(:)

 end module system
!-----------------!

!--------------------!
 program timevolution
!--------------------!
 use system; implicit none

 integer :: t,nt,wf,a0
 real(8) :: hh,jj,dh1,dh2,dj1,dj2,tt,dt

 open(10,file='read.in',status='old')
 read(10,*)nn
 read(10,*)tt,dt,wf
 close(10)

! nn = number of spins
! tt = integration time
! dt = integration time step
! wf = output results every wf time step

 call lattice()
 call makehamiltonian()
 call couplings(a0)

 call initstate()
 write(*,*)'d'
 nt=int(tt/dt+0.1d0)    ! total number of RK integration steps

 open(10,file='r.dat',status='replace')
 close(10)

 do t=0,nt-1
    call ramp(t,nt,hh,jj,dh1,dh2,dj1,dj2)
    call timestep(dt,hh,jj,dh1,dh2,dj1,dj2)               
    if (mod(t+1,wf)==0) call writedata(t+1,nt,a0)
 enddo

 end program timevolution
!------------------------!

!-----------------------------!
 subroutine writedata(t,nt,a0)
!-----------------------------!
 use system; implicit none

 integer :: a,t,nt,a0
 real(8) :: ez,p0

 ez=sum(hamz(:)*abs(ff(:))**2)-hamz(a0)
 p0=2.d0*abs(ff(a0))**2
 open(10,file='r.dat',status='old',position='append')
 write(10,'(3f18.12)')dble(t)/dble(nt),p0,ez
 close(10)

 end subroutine writedata
!------------------------!

!-------------------------------------------!
 subroutine ramp(t,nt,hh,jj,dh1,dh2,dj1,dj2)
!------------------------------------------------------------------!
! Gives the Hamiltonian parameters jj and hh as a function of time
! according to a linear ramp as a function of time from hh=1, jj=0
! at t=0 to hh=0,jj=1 at t=nt (t is the integer time step in the
! RK integration process. The other output numbers are the parameter
! values at other time points (t+1/2 and t=1) needed in the RK
! integration.
!-------------------------------------------------------------------
 use system; implicit none

 integer :: t,nt
 real(8) :: ss,ds,hh,jj,dh1,dh2,dj1,dj2

 ss=dble(t)/dble(nt)
 jj=ss
 hh=1.d0-ss
 ds=1.d0/dble(nt)    
 dj1=0.5d0*ds
 dh1=-dj1
 dj2=ds
 dh2=-dj2

 end subroutine ramp
!-------------------!

!---------------------------------------------!
 subroutine timestep(dt,hh,jj,dh1,dh2,dj1,dj2)
!------------------------------------------------------------------------------!
! Uses the RK method to integrate the Schrodinger equation by one time step dt !
!------------------------------------------------------------------------------!
 use system; implicit none

 integer :: i,j
 real(8) :: hh,jj,dh1,dh2,dj1,dj2,dt

 f2=ff
 f0=ff
 call hoperation(hh,jj)
 f1=-dt*(0,1)*f1
 f2=f2+f1/6.d0
 f0=ff+0.5d0*f1
 call hoperation(hh+dh1,jj+dj1)
 f1=-dt*(0,1)*f1
 f2=f2+f1/3.d0
 f0=ff+0.5d0*f1
 call hoperation(hh+dh1,jj+dj1)
 f1=-dt*(0,1)*f1
 f2=f2+f1/3.d0
 f0=ff+f1
 call hoperation(hh+dh2,jj+dj2)
 f1=-dt*(0,1)*f1
 f2=f2+f1/6.d0
 ff=f2/sqrt(dot_product(f2,f2))

 ! Explicit symmetrization of the state
 do i=0,2**nn-1
    j=ieor(i,zz)
    if (i<j) ff(j)=ff(i)
 enddo

 end subroutine timestep
!-----------------------!

!----------------------------!
 subroutine hoperation(hh,jj)
!----------------------------!
 use system; implicit none

 integer :: i,j,a,b,c
 real(8) :: hh,jj

 f1(:)=jj*hamz(:)*f0(:)
 j=0
 do a=0,2**nn-1
    f1(hamx(j+1:j+nxop(a)))=f1(hamx(j+1:j+nxop(a)))-hh*f0(a)
    j=j+nxop(a)
 enddo

! Take care of spin-reflection related states not included above, 
! using spin-reflection symmetry

 do b=0,2**nn-1           
    c=ieor(b,zz)          
    if (b<c) f1(c)=f1(b)  
 enddo                    

 end subroutine hoperation
!-------------------------!

!----------------------------!
 subroutine makehamiltonian()
!-------------------------------------------------------!
! allocates the storage arrays for the hamiltoniand and !
! makes the list of the non-zero off-diagonal elements  !
! (their locations in the matrix)                       !
!-------------------------------------------------------!
 use system; implicit none

 integer :: i,j,a,b,c,hz

 allocate(nxop(0:2**nn-1))    ! number of off-diagonal elements for each basis state
 allocate(hamx(nn*2**(nn-1))) ! location of off-diagonal H elements compact-stored here
 allocate(hamz(0:2**nn-1))    ! diagonal H elements will be stored here

! zz is the number used to spin-reflect a state with ieor().
 zz=0
 do i=0,nn-1
    zz=ibset(zz,i)
 enddo

 j=0
 do a=0,2**nn-1
    nxop(a)=0
    do i=0,nn-1
       if (btest(a,i)) then
          b=ibclr(a,i)
       else
          b=ibset(a,i)
       endif

      ! Spin-reflect the state b, giving c.
      ! Only include Hamiltonian term if b<c.
       c=ieor(b,zz)
       if (b<c) then
          j=j+1
          hamx(j)=b
          nxop(a)=nxop(a)+1
       endif
    enddo

    ! Sort the Hamiltonian terms, can speed up the execution
    if (nxop(a)>1) call sort(nxop(a),hamx(j-nxop(a)+1:j))

 enddo

 end subroutine makehamiltonian
!------------------------------!

!--------------------!                                                                                           
 subroutine sort(n,a)
!----------------------------!
! Based on Numerical Recipes.
!----------------------------!
 implicit none

 integer :: i,j,l,n,ir,rra,a(n)

 l=n/2+1
 ir=n
 do
    if (l>1) then
       l=l-1
          rra=a(l)
       else
       rra=a(ir)
       a(ir)=a(1)
       ir=ir-1
       if (ir==1) then
          a(1)=rra
          return
       endif
    endif
    i=l
    j=l+l
    do
       if (j<=ir) then
          if (j<ir) then
             if (a(j)<a(j+1)) j=j+1
          endif
          if (rra<a(j)) then
             a(i)=a(j)
             i=j
             j=j+j
          else
             j=ir+1
          endif
       else
          exit
       endif
    enddo
    a(i)=rra
 enddo

 end subroutine sort
!-------------------!                                                                                                            
!----------------------!
 subroutine initstate()
!----------------------!
 use system; implicit none

 allocate(ff(0:2**nn-1))
 allocate(f0(0:2**nn-1))
 allocate(f1(0:2**nn-1)) 
 allocate(f2(0:2**nn-1)) 

 ff(:)=1.d0
 ff(:)=ff(:)/sqrt(dot_product(ff(:),ff(:)))
 
 end subroutine initstate
!------------------------!

!--------------------!
 subroutine lattice()
!------------------------------------!
! generates the lattice connectivity !
! (sites connected by interactions)  !
!------------------------------------!
 use system; implicit none

 integer :: i,j,b

 nb=(nn*(nn-1))/2
 allocate(jbnd(nb))
 allocate(site(2,nb))
 b=0
 do i=0,nn-1
    do j=i+1,nn-1
       b=b+1
       site(1,b)=i
       site(2,b)=j
    enddo
 enddo

 end subroutine lattice
!----------------------!

!------------------------!
 subroutine couplings(a0)
!------------------------!
 use system; implicit none

 integer :: i,a,b,a0
 real(8) :: hz

 jbnd=1.d0
 do a=0,2**nn-1
    hz=0.d0
    do i=1,nb
       if (btest(a,site(1,i)).eqv.btest(a,site(2,i))) then
          hz=hz-jbnd(i)
       else
          hz=hz+jbnd(i)
       endif
    enddo
    hamz(a)=hz/(dble(nn-1)/2.d0)
 enddo

 hz=0.d0
 do a=0,2**nn-1
    b=ieor(a,zz)
    if (a<b.and.hamz(a)<hz) then
       hz=hamz(a)
       a0=a
    endif
 enddo

 end subroutine couplings
!------------------------!

!--------------------!
 subroutine cleanup()
!--------------------!
 use system; implicit none

 deallocate(site)
 deallocate(jbnd)
 deallocate(hamz)
 deallocate(hamx)
 deallocate(nxop)
 deallocate(ff)
 deallocate(f0)
 deallocate(f1)
 deallocate(f2)

 end subroutine cleanup
!----------------------!

