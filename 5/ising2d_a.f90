!----------------------!
 module systemvariables
 implicit none

 integer :: l                        ! system length
 integer :: n                        ! number of spins (n=l*l)
 real    :: pflip(-4:4)              ! flip probabilities
 integer, allocatable :: spin(:,:)   ! spin array

 end module systemvariables
!--------------------------!

!-------------------!
 module autocorrdata
 implicit none

 integer :: ntime
 integer :: count1
 integer :: count2      
 real(8), allocatable :: aobs(:)
 real(8), allocatable :: tobs(:,:)
 real(8), allocatable :: acor(:,:)

 end module autocorrdata
!-----------------------!

!-------------------------------------------------------------------!
!Metropolis-algorithm simulation of the two-dimensional Ising model.!
!This version calculates autocorrelation functions.                 !
!-------------------------------------------------------------------!
!Reads the following from a file 'ising.in':                        !
! Line 1: l, temp = system length, temperature (T/J).               !
! Line 2: bins, binsteps = number of bins, MC steps per bin).       !
!                          One bin is added for equilibration.      !
! Line 3: nt,tfreq = number of autocorrelations to be measured      !
!                    and the spacing (in MC steps) between them.    !
!Writes these bin averages to the file 'bindata.dat':               ! 
! <E>/N, <E**2>/N**2, <|M|>/N, <M**2>/N**2                          !
!These results should be processed by the program average.f90.      !
!Writes autocorrelation functions for E and |M| to 'auto.dat'.      !
!This data should be processed using the program autoaverage.f90.   !
!-------------------------------------------------------------------!
 program ising2d_a
!-----------------!
 use systemvariables
 implicit none

 integer :: i,j,bins,binsteps,ntime,tfreq
 real    :: temp

 open(1,file='ising.in',status='old')
 read(1,*)l,temp
 read(1,*)bins,binsteps
 read(1,*)ntime,tfreq
 close(1)

 call initialize(temp,ntime)
 do i=1,binsteps
    call mcstep
 end do
 do j=1,bins
    call resetdatasums
    do i=1,binsteps
       call mcstep
       call measure(i,tfreq)
    end do
    call writebindata(n,binsteps)
    call writeautodata
 end do
 call cleanup
     
!---------------------!
 end program ising2d_a
!---------------------!

!--------------------------------------------!
!Carries out one Monte Carlo step, defined as! 
!n flip attempts of randomly selected spins. !
!--------------------------------------------!
 subroutine mcstep
!-----------------!
 use systemvariables
 implicit none

 integer :: i,s,x,y,s1,s2,s3,s4
 real    :: r

 do i=1,n
    call random_number(r)
    s=int(r*n)
    x=mod(s,l); y=s/l
    s1=spin(mod(x+1,l),y)
    s2=spin(mod(x-1+l,l),y)
    s3=spin(x,mod(y+1,l))
    s4=spin(x,mod(y-1+l,l))
    call random_number(r)
    if (r<=pflip(spin(x,y)*(s1+s2+s3+s4))) spin(x,y)=-spin(x,y)
 end do

 end subroutine mcstep
!---------------------!

!--------------------------------------------------------!
!Measures the energy (enrg1), the squared energy (enrg2),! 
!the absolute value of the magnetization (magn1), and the!
!squared magnetization (magn2).                          !
!--------------------------------------------------------!
 subroutine measure(step,tfreq)
!------------------------------!
 use systemvariables
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2
 common/measurments/enrg1,enrg2,magn1,magn2

 integer :: x,y,x1,y1,e,m,step,tfreq
 real(8) :: obs(2)

 e=0
 m=0
 do y=0,l-1
    y1=mod(y+1,l)
    do x=0,l-1
       x1=mod(x+1,l)
       e=e-spin(x,y)*(spin(x1,y)+spin(x,y1))
       m=m+spin(x,y)
    end do
 end do
 m=abs(m)
 enrg1=enrg1+dble(e)
 enrg2=enrg2+dble(e)**2
 magn1=magn1+dble(m)
 magn2=magn2+dble(m)**2

 if (mod(step,tfreq)==0) then
    obs(1)=dble(e)/dble(n)
    obs(2)=dble(m)/dble(n)
    call autocorr(obs,2)
 end if

 end subroutine measure
!----------------------!

!-------------------------------------------------!
!Writes accumulated data to the file 'bindata.dat'!
!-------------------------------------------------!
 subroutine writebindata(n,steps)
!--------------------------------!
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2
 common/measurments/enrg1,enrg2,magn1,magn2

 integer :: n,steps
 
 open(1,file='bindata.dat',status='unknown',position='append')
 enrg1=enrg1/(dble(steps)*dble(n))
 enrg2=enrg2/(dble(steps)*dble(n)**2)
 magn1=magn1/(dble(steps)*dble(n))
 magn2=magn2/(dble(steps)*dble(n)**2)
 write(1,1)enrg1,enrg2,magn1,magn2
 1 format(4f18.12)
 close(1)

 end subroutine writebindata
!---------------------------!

!-------------------------------------------------------------------!
!Calculates the autocorrelation functions of nobs observables in the! 
!vector obs. A series of ntime measurements are stored in the array !
!tobs. Autocorrelations calculated on the basis of this data are    !
!accumulated in acor. The plain averages of the observables are     !
!accumulated in aobs.                                               !
!-------------------------------------------------------------------!
 subroutine autocorr(obs,nobs)
!-----------------------------!
 use autocorrdata
 implicit none

 integer :: j,nobs
 real(8) :: obs(nobs)

 if (count1<ntime) then
    count1=count1+1
    tobs(:,count1)=obs(:)
 else
    count2=count2+1
    aobs(:)=aobs(:)+obs(:)
    do j=1,ntime-1
       tobs(:,j)=tobs(:,j+1)
    end do
    tobs(:,ntime)=obs(:)
    do j=0,ntime-1
       acor(:,j)=acor(:,j)+tobs(:,1)*tobs(:,1+j)
    end do
 end if

 end subroutine autocorr
!-----------------------!

!----------------------------------------------------------!
!Writes the autocorrelation functions to the file auto.dat.!
!----------------------------------------------------------!
 subroutine writeautodata
!------------------------!
 use autocorrdata
 implicit none

 integer :: j

 open(2,file='auto.dat',status='unknown',position='append')
 write(2,1)aobs(:)/dble(count2)
 do j=0,ntime-1
    write(2,1)acor(:,j)/dble(count2)
 end do
 1 format(4f20.10)
 close(2)

 end subroutine writeautodata
!------------------------!

!--------------------------------------------!
!Sets the data accumulation variables to zero!
!--------------------------------------------!
 subroutine resetdatasums
!------------------------!
 use autocorrdata
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2
 common/measurments/enrg1,enrg2,magn1,magn2

 enrg1=0.d0
 enrg2=0.d0
 magn1=0.d0
 magn2=0.d0

 count1=0
 count2=0
 tobs=0.d0
 acor=0.d0
 aobs=0.d0

 end subroutine resetdatasums
!----------------------------!

!-------------------------------------------------------------------!
!Various initializations: system size n=l*l, flipping probabilities !
!pflip, random number generator, initial random spin state          !
!allocation of arrays used for evaluating autocorrelation functions.!
!-------------------------------------------------------------------!
 subroutine initialize(temp,nt)
!------------------------------!
 use systemvariables
 use autocorrdata
 implicit none

 integer :: i,j,nt,ns
 integer, allocatable :: seed(:)
 real    :: temp,r

 n=l*l
 pflip=0.
 do i=-4,4,2
    pflip(i)=exp(-i*2./temp)
 end do

 call random_seed(ns)
 allocate(seed(ns))
 open(1,file='seed.in',status='old')
 read(1,*)seed
 close(1)
 call random_seed(put=seed)

 allocate (spin(0:l-1,0:l-1))
 do j=0,l-1
 do i=0,l-1
    call random_number(r)
    spin(i,j)=2*int(2.*r)-1
 end do
 end do
 ntime=nt
 allocate(aobs(2))
 allocate(tobs(2,ntime))
 allocate(acor(2,0:ntime))

 end subroutine initialize
!-------------------------!

!------------------!
!Deallocates arrays!
!------------------!
 subroutine cleanup
!------------------!
 use systemvariables
 use autocorrdata
 
 deallocate(spin)
 deallocate(aobs)
 deallocate(tobs)
 deallocate(acor)

 end subroutine cleanup
!----------------------!
