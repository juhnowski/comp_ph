!---------------------------------------------------------------------------!
!Calculates the area of a circle with radius 1 using Monte Carlo integration!
!The final average is written to the file 'pi.dat'.                         !
!---------------------------------------------------------------------------!

!------------------------!
 program montecarlocircle
 implicit none

 integer, parameter :: nbin=20

 integer :: i,j,n
 real(8) :: x2,y2,av,er,sum!,rand

 !external :: rand
 
 write(*,'(a)',advance='no')'Number of samples per bin: ';read*,n
 call initran(1)

 av=0.d0
 er=0.d0
 do j=1,nbin
   sum=0.d0
   do i=1,n
     x2=(dble(ran()-0.5))**2
     y2=(dble(ran()-0.5))**2
     if (x2+y2 <= 0.25d0) sum=sum+1.d0
   end do
   sum=4.d0*sum/dble(n)
   av=av+sum
   er=er+sum**2
   write(*,'(a,I3,a,F11.6)')'Done bin ',j,'   Result of this bin: ',sum
 end do
 av=av/dble(nbin)
 er=er/dble(nbin)
 er=sqrt((er-av**2)/dble(nbin-1))
 open(1,file='pi.dat',status='replace')
 write(1,'(2F12.6)')av,er
 close(1)

 end program montecarlocircle
!----------------------------!

