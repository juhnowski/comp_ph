!--------------------!
 program multiaverage
!--------------------!
 implicit none

 integer :: i,d,n,it,nt,bins
 real(8) :: data1(5),av,er,avc,erc
 integer, allocatable :: l(:)
 real(8), allocatable :: temp(:)
 real(8), allocatable :: data(:,:)

 open(1,file='ising.in',status='old')
 read(1,*)d;  allocate(l(d))
 read(1,*)l;  n=product(l)
 read(1,*)i
 read(1,*)nt; allocate(temp(nt))
 do it=1,nt
    read(1,*)temp(it)
 enddo
 close(1)

 open(1,file='e.dat',status='replace')
 open(2,file='c.dat',status='replace')
 open(3,file='m1.dat',status='replace')
 open(4,file='m2.dat',status='replace')
 open(5,file='x.dat',status='replace')
 open(7,file='q.dat',status='replace')

 do it=1,nt

    open (10,file='bindata.dat',status='old')
    bins=0
    do 
       read(10,*,end=1)i
       if (i==it) bins=bins+1
    end do
    1 allocate(data(bins,5))
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
    write(1,3)temp(it),av,er    

    call averageanderror2(data(:,1),data(:,2),bins,n,temp(it),av,er)
    write(2,3)temp(it),av,er    

    call averageanderror1(data(:,3),bins,av,er)
    write(3,3)temp(it),av,er    

    call averageanderror1(data(:,4),bins,av,er)
    call averageanderror1(data(:,5),bins,avc,erc)
    write(4,3)temp(it),av,er,avc,erc    

    call averageanderror2(data(:,3),data(:,4),bins,n,temp(it),av,er)
    write(5,3)temp(it),av,er    

    call averageanderror3(data(:,3),data(:,4),bins,av,er)
    write(7,3)temp(it),av,er

    3 format(5f15.8)

    deallocate(data)

 end do

 close(1)
 close(2)
 close(3)
 close(4)
 close(5)
 close(6)

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

!----------------------------------------------------------!
 subroutine averageanderror2(data1,data2,bins,n,temp,av,er)
!----------------------------------------------------------!
 implicit none

 integer :: i,n,bins
 real(8) :: data1(bins),data2(bins),temp,diff,av,er

 av=sum(data2-data1**2)/dble(bins)
 er=sum((data2-data1**2)**2)/dble(bins)
 er=sqrt((er-av**2)/dble(bins-1))
 av=av*dble(n)/temp
 er=er*dble(n)/temp

 end subroutine averageanderror2
!-------------------------------!

!---------------------------------------------------!
 subroutine averageanderror3(data1,data2,bins,av,er)
!---------------------------------------------------!
 implicit none

 integer :: i,n,bins
 real(8) :: data1(bins),data2(bins),d1,d2,av,er

 d1=sum(data1)/dble(bins)
 d2=sum(data2)/dble(bins)
 av=d2/d1**2
 er=sqrt(sum((data2/data1**2-av)**2)/dble(bins-1))

 end subroutine averageanderror3
!-------------------------------!
