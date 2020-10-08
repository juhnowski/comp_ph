!----------------------!
 module systemvariables
!----------------------!

 real, parameter :: r2mol1=1.0       ! square of r1
 real, parameter :: r2mol2=9.0       ! square of r2
 real, parameter :: potential=-1.0  

 integer :: n            ! number of particles
 integer :: n0           ! number of overlapping potentials
 real    :: l            ! system length
 real    :: temp         ! temperature
 real    :: delta        ! max coordinate change = +/- delta/2

 real, allocatable :: xyz(:,:)    ! particle coordinates

 end module systemvariables
!--------------------------!

!--------------------------!
 program moleculemontecarlo
!--------------------------!
 use systemvariables
 implicit none

 integer :: i,j,it,nt,bins,steps
 real    :: tmax,dt,arate,accepted

 open(1,file='param.in',status='old')
 read(1,*)n,l
 read(1,*)tmax,dt,nt
 read(1,*)bins,steps
 close(1)

 allocate(xyz(3,n)) 
 call initrand(1)
 call initialstate
 delta=sqrt(r2mol2)

 do it=0,nt-1
    temp=tmax-dt*it

    arate=0.
    do i=1,steps
       call mcstep(accepted)
       arate=arate+accepted
       if (mod(i,steps/20)==0) call adjustedelta(steps/20,arate,delta,l)
    enddo

    do j=1,bins

       call resetbindata
       do i=1,steps
          call mcstep(accepted)
          call measure(accepted)
       enddo
       call writebindata(it,steps)

       call openlog
       write(10,*)it,j
       call closelog

    enddo

 enddo

 end program moleculemontecarlo
!------------------------------!

!---------------------------!
 subroutine mcstep(accepted)
!---------------------------!
 use systemvariables
 implicit none

 integer :: i,j,n1 
 real    :: p,r,r2,x0,y0,z0,x1,y1,z1,dx,dy,dz,accepted
 logical :: accept

 accepted=0.
 do j=1,n
    x0=xyz(1,j); y0=xyz(2,j); z0=xyz(3,j)
    call newcoordinates(l,delta,x0,y0,z0,x1,y1,z1)
    n0=0; n1=0
    do i=1,n
       if (i/=j) then
          dx=abs(x1-xyz(1,i)); dx=min(dx,l-dx)
          dy=abs(y1-xyz(2,i)); dy=min(dy,l-dy)
          dz=abs(z1-xyz(3,i)); dz=min(dz,l-dz)
          r2=dx**2+dy**2+dz**2
          if (r2 < r2mol1) goto 1
          if (r2 < r2mol2) n1=n1+1
       endif
    enddo
    accept=.true.
    if (n1 < n0) then 
       call random_number(r)
       if (r > exp(-potential*real(n1-n0)/temp)) accept=.false.
    endif
    if (accept) then
       n0=n1
       xyz(1,j)=x1; xyz(2,j)=y1; xyz(3,j)=z1
       accepted=accepted+1.
    endif
    1 continue
 enddo
 accepted=accepted/real(n)

 end subroutine mcstep
!---------------------!

!----------------------------!
 subroutine measure(accepted)
!----------------------------!
 use systemvariables
 implicit none

 real(8) :: energy,arate
 common/measurements/energy,arate

 real :: accepted

 energy=energy+potential*dble(n0)/dble(n)
 arate=arate+dble(accepted)

 end subroutine measure
!----------------------!

!---------------------------------!
 subroutine writebindata(it,steps)
!---------------------------------!
 implicit none

 real(8) :: energy,arate
 common/measurements/energy,arate

 integer :: it,steps

 energy=energy/dble(steps)
 arate=arate/dble(steps)
 open(1,file='bindata.dat',status='unknown',position='append')
 write(1,*)it,energy,arate
 close(1)

 end subroutine writebindata
!---------------------------!

!-----------------------!
 subroutine resetbindata
!-----------------------!
 implicit none

 real(8) :: energy,arate
 common/measurements/energy,arate

 energy=0.d0
 arate=0.d0

 end subroutine resetbindata
!---------------------------!

!-----------------------------------------------------!
 subroutine newcoordinates(l,delta,x0,y0,z0,x1,y1,z1)
!-----------------------------------------------------!
 implicit none

 real :: r,x0,y0,z0,x1,y1,z1,l,delta

 call random_number(r)
 x1=x0+delta*(r-0.5) 
 if (x1 < 0) x1=x1+l
 if (x1 > l) x1=x1-l
 call random_number(r)
 y1=y0+delta*(r-0.5)
 if (y1 < 0) y1=y1+l
 if (y1 > l) y1=y1-l
 call random_number(r)
 z1=z0+delta*(r-0.5)
 if (z1 < 0) z1=z1+l
 if (z1 > l) z1=z1-l
   
 end subroutine newcoordinates
!-----------------------------!

!--------------------------------------------!
 subroutine adjustedelta(steps,arate,delta,l)
!--------------------------------------------!
 implicit none

 integer :: steps
 real    :: arate,delta,l

 arate=arate/float(steps)
 if (arate < 0.4) delta=delta/1.5
 if (arate > 0.6 .and. delta < l/4.) delta=delta*1.5
 call openlog
 write(10,*)'Acceptance rate, delta : ',arate,delta
 call closelog
 arate=0.

 end subroutine adjustedelta
!---------------------------!

!-----------------------!
 subroutine initialstate
!-----------------------!
 use systemvariables
 implicit none

 integer :: i,j
 real    :: r,r2,dx,dy,dz

 n0=0
 do j=1,n
    1 call random_number(r)
    xyz(1,j)=r*l
    call random_number(r)
    xyz(2,j)=r*l
    call random_number(r)
    xyz(3,j)=r*l
    do i=1,j-1
       dx=abs(xyz(1,j)-xyz(1,i))
       dx=min(dx,l-dx)
       dy=abs(xyz(2,j)-xyz(2,i))
       dy=min(dy,l-dy)
       dz=abs(xyz(3,j)-xyz(3,i))
       dz=min(dz,l-dz)
       r2=dx**2+dy**2+dz**2
       if (r2 < r2mol1) then
          goto 1
       elseif (r2 < r2mol2) then
          n0=n0+1
       endif
    enddo
 enddo

 end subroutine initialstate
!---------------------------!

!----------------------!
 subroutine initrand(w)
!----------------------!
 implicit none

 integer :: i,w,s
 integer, allocatable :: seed(:)
 real    :: r

 call random_seed(s)
 allocate(seed(s))
 open(1,file='seed.in',status='old')
 read(1,*)seed(:)
 close(1)
 call random_seed(put=seed)
 if (w /= 0) then
   open(1,file='seed.in',status='replace')
   do i=1,s
      call random_number(r)
      write(1,*)abs(int((r-0.5)/0.23283064e-9))
   enddo
 endif

 end subroutine initrand
!-----------------------!

!------------------!
 subroutine openlog
!------------------!

 open(10,file='log.txt',status='unknown',position='append')

 end subroutine openlog
!----------------------!

!-------------------!
 subroutine closelog
!-------------------!

 close(10)

 end subroutine closelog
!-----------------------!
