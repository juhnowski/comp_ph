!---------------!
 program rantest
!-------------------------------------------------------------!
! Demonstration of the use of the random number function ran()
!-------------------------------------------------------------!
 implicit none

 integer :: i
 real(8) :: ran

 external :: ran

 call initran(1)  ! "1" means the seed file will be replaced; "0" means no change

do i=1,10
   write(*,*)i,ran(i)
 enddo

 end program rantest
!-------------------!
