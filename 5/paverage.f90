!--------------------!
 program multiaverage
!--------------------!
 implicit none

 integer :: i,d,n,it,nt,bins
 real(8) :: data1(2),av,er,temp,tmax,dt
 real(8), allocatable :: data(:,:)

 open(1,file='param.in',status='old')
 read(1,*)n
 read(1,*)tmax,dt,nt
 close(1)

 open(1,file='e.dat',status='replace')
 open(2,file='a.dat',status='replace')

 do it=0,nt-1
    temp=tmax-it*dt

    open (10,file='bindata.dat',status='old')
    bins=0
    do 
       read(10,*,end=1)i
       if (i==it) bins=bins+1
    end do
    1 allocate(data(bins,2))
    print*,it,bins
    rewind(10)
    bins=0
    do        
       read(10,*,end=2)i,data1
       if (i==it) then
          bins=bins+1
          data(bins,:)=data1
       end if
    end do
    2 close(10)

    call averageanderror1(data(:,1),bins,av,er)
    write(1,3)temp,av,er    

    call averageanderror1(data(:,2),bins,av,er)
    write(2,3)temp,av,er    

    3 format(3f15.8)

    deallocate(data)

 end do

 close(1)
 close(2)

 end program multiaverage
!------------------------!

!--------------------------------------------!
 subroutine averageanderror1(data,bins,av,er)
!--------------------------------------------!
 implicit none

 integer :: i,bins
 real(8) :: data(bins),av,er

 av=sum(data)/dble(bins)
 er=sum(data**2)/dble(bins)
 er=sqrt((er-av**2)/dble(bins-1))

 end subroutine averageanderror1
!-------------------------------!
