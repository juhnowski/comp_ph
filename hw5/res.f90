!----------------!
 program averages
!------------------------------------------------------------------------------------!
! This program computes averages and error bars for binned data read from a file     ! 
! 'res.dat', which contains output of simulations done in PY502 assignment #5. Input !
! data are read from 'read.in', the same file as used in the simulation program:     !
!   l                                                                                !
!   i,i,st                                                                           !
! where 'l' is the system length, 'st' is the number of MC steps or 0, and 'i'       !
! is not used in this program. If st is not 0, res.dat should contain lines with     !
! t,m(t), where t (integer) is the time step t=1,...,st and m(t) (real) is the bin   !
! average for that time. An arbitrary number of bins can appear consecutively in     !
! the file. If st=0 the file contains rows of single real number. Output in          !
! either case is written to the file 'r.dat'.                                        !
!------------------------------------------------------------------------------------!
 implicit none

 integer :: i,l,st

 real(8), allocatable :: av(:,:)

 open(10,file='read.in',status='old')
 read(10,*)l
 read(10,*)i,i,st
 close(10)

 if (st/=0) then
    call average1(st)
 else
    call average2(l)
 endif

 end program averages
!--------------------!

!-----------------------!
 subroutine average1(st)
!-----------------------!
 implicit none

 integer :: i,j,st,nr
 real(8) :: mg

 real(8), allocatable :: av(:,:)

 allocate(av(2,st))

 nr=0
 av=0.d0
 open(10,file='res.dat',status='old')
 do
    do i=1,st
       read(10,*,end=10)j,mg
       av(1,i)=av(1,i)+mg
       av(2,i)=av(2,i)+mg**2
    enddo
    nr=nr+1
 enddo
 10 close(10)
 write(*,*)'Number of bins read: ',nr
 av=av/dble(nr)
 av(2,:)=sqrt((av(2,:)-av(1,:)**2)/dble(nr))

 open(10,file='r.dat')
 do i=1,st
    write(10,'(i8,2f14.8)')i,av(:,i)
 enddo
 close(10)

 deallocate(av)

 end subroutine average1
!-----------------------!

!----------------------!
 subroutine average2(l)
!----------------------!
 implicit none

 integer :: i,l,st,nr
 real(8) :: t,av(2)

 nr=0
 av=0.d0
 open(10,file='res.dat',status='old')
 do 
    read(10,*,end=10)t
    av(1)=av(1)+t
    av(2)=av(2)+t**2
    nr=nr+1
 enddo
 10 close(10)
 write(*,*)nr
 av=av/dble(nr)
 av(2)=sqrt((av(2)-av(1)**2)/dble(nr))

 open(10,file='r.dat')
 write(10,'(i8,2f16.4)')l,av(:)
 close(10)

 end subroutine average2
!-----------------------!
