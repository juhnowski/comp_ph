!-------------!
 module system
!-------------!
 implicit none

 integer :: l                        ! system length
 integer :: n                        ! number of spins (n=l*l)
 real(8) :: pflip(-4:4)              ! flip probabilities
 integer, allocatable :: spin(:)     ! spin array

 real(8) :: tt0                      ! measured time when magnetization reaches 0
 real(8), allocatable :: mag(:)      ! array for measred magnetization vs time

 end module system
!-----------------!

!-------------------------------------------------------------------!
!Metropolis-algorithm simulation of the two-dimensional Ising model.!
!Starts at a fully polarized state and computes the magnetization as!
!a function of simulation time or the time to reach m=0.
!-------------------------------------------------------------------!
 program ising2d
!---------------!
 use system; implicit none

 integer :: i,j,k,m,t,bins,reps,steps
 real(8) :: temp

 open(10,file='read.in',status='old')
 read(10,*)l,temp
 read(10,*)bins,reps,steps
 close(10)

 call initran(1)

 n=l*l
 do i=-4,4
    pflip(i)=exp(-2.d0*dble(i)/temp)
 enddo
 allocate (spin(0:n-1))
 if (steps/=0) allocate (mag(steps))
     
 do k=1,bins   
    if (steps==0) then
        tt0=0.d0
    else
        mag(:)=0.d0
    endif
    do j=1,reps
       m=n
       spin(:)=1
       if (steps==0) then
          i=0
          t=0
          do
             call mcstep(m,t)
             if (m==0) then
                tt0=tt0+dble(i)+dble(t)/dble(n)
                exit
             endif
             i=i+1
          enddo
       else
          t=1
          do i=1,steps
             call mcstep(m,t)
             mag(i)=mag(i)+dble(m)/dble(n)
          enddo
       endif
    enddo
    call writedata(reps,steps)
 enddo

 deallocate(spin)
 if (steps/=0) deallocate(mag)
     
 end program ising2d
!-------------------!

!--------------------------------------------!
!Carries out one Monte Carlo step, defined as! 
!n flip attempts of randomly selected spins. !
!--------------------------------------------!
 subroutine mcstep(m,t)
!----------------------!
 use system; implicit none

 integer :: i,m,t,s,x,y,s1,s2,s3,s4

 real(8), external :: ran

 do i=1,n
    s=int(ran()*n)
    x=mod(s,l); y=s/l
    s1=spin(mod(x+1,l)+y*l)
    s2=spin(x+mod(y+1,l)*l)
    s3=spin(mod(x-1+l,l)+y*l)
    s4=spin(x+mod(y-1+l,l)*l)
    if (ran()<pflip(spin(s)*(s1+s2+s3+s4))) then
        spin(s)=-spin(s)
        m=m+2*spin(s)
    endif
    if (t==0.and.m==0) then
       t=i
       return
    endif
 enddo

 end subroutine mcstep
!---------------------!

!-----------------------------------------------!
! Writes accumulated data to the file 'res.dat' !
!-----------------------------------------------!
 subroutine writedata(reps,steps)
!--------------------------------!
 use system; implicit none

 integer :: i,reps,steps
 
 open(10,file='res.dat',status='unknown',position='append')
 if (steps==0) then
    tt0=tt0/dble(reps)
    write(10,'(f16.8)')tt0
 else
    mag=mag/dble(reps)
    do i=1,steps
        write(10,'(i8,f14.8)')i,mag(i)       
    enddo
 endif
 close(10)

 end subroutine writedata
!------------------------!

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
