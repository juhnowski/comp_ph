 implicit none

 integer :: i,size
 integer, allocatable :: seed(:)
 real :: r

 call random_seed(size)
 allocate (seed(size))
 write(*,*)'give ',size,' random seeds '
 read(*,*)seed
 write(*,*)

 call random_seed(put=seed)
 do i=1,10
    call random_number(r)
    write(*,*)r
 end do

 write(*,*)
 call random_seed(get=seed)
 write(*,*)seed

 end
