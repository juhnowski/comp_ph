!-------------------!
 program autoaverage
!-------------------!
 implicit none

 integer :: i,j,nt,dt,nt1,nt2,bins,seed(3)
 real(8) :: r,a1(2),tint(2),aint(2),eint(2)
 real(8), allocatable :: adata(:,:,:)
 real(8), allocatable :: abin(:,:)
 real(8), allocatable :: avt(:,:)
 real(8), allocatable :: ert(:,:)
 
 open(1,file='ising.in',status='old')
 read(1,*)nt 
 read(1,*)nt
 read(1,*)nt,dt
 close(1)
 seed(1)=355682928; seed(2)=6473671; seed(3)=8938927
 call random_seed(put=seed)

 open(1,file='auto.dat',status='old')
 bins=0
 do
    do i=0,nt
       read(1,*,end=1)a1(:)
    end do
    bins=bins+1
 end do
 1 write(*,*)'Number of bins: ',bins

 allocate(adata(2,0:nt,bins))
 allocate(abin(2,nt))
 allocate(avt(2,nt))
 allocate(ert(2,nt))

 rewind(1)
 do j=1,bins
    do i=0,nt
       read(1,*)adata(:,i,j)
    end do
 end do
 close(1)


!*** Average autocorrelation function
 avt(:,:)=0.d0
 do j=1,bins
    do i=1,nt
       avt(:,i)=avt(:,i)+adata(:,i,j)-adata(:,0,j)**2
    end do
 end do
 do i=nt,1,-1
    avt(:,i)=avt(:,i)/avt(:,1)
 end do

!*** Error bars of autocorrelation function
 ert(:,:)=0.d0
 do j=1,bins
    do i=1,nt
       abin(:,i)=adata(:,i,j)-adata(:,0,j)**2
    end do
    do i=nt,1,-1
       abin(:,i)=abin(:,i)/abin(:,1)
       ert(:,i)=ert(:,i)+(avt(:,i)-abin(:,i))**2
    end do
 end do
 ert(:,:)=sqrt(ert(:,:)/dble(bins*(bins-1)))

!*** Integrated autocorrelation times
 aint(:)=0.d0
 do i=1,nt
    if (ert(1,i)<avt(1,i)) then
       aint(1)=aint(1)+avt(1,i)
    else
       nt1=i
       exit
    end if
 end do
 do i=1,nt
    if (ert(2,i)<avt(2,i)) then
       aint(2)=aint(2)+avt(2,i)
    else
       nt2=i
       exit
    end if
 end do
 aint(:)=aint(:)-0.5d0

!*** Error bars of integrated autocorrelation times
 eint(:)=0.d0
 do j=1,bins
    a1(:)=0.d0
    do i=1,nt
       abin(:,i)=adata(:,i,j)-adata(:,0,j)**2
    end do
    do i=nt,1,-1
       abin(:,i)=abin(:,i)/abin(:,1)
       if (i<=nt1) a1(1)=a1(1)+abin(1,i)
       if (i<=nt2) a1(2)=a1(2)+abin(2,i)
    end do
    a1(:)=a1(:)-0.5d0
    eint(:)=eint(:)+(aint(:)-a1(:))**2
 end do
 eint(:)=sqrt(eint(:)/dble(bins*(bins-1)))

 aint(:)=dt*aint(:)
 eint(:)=dt*eint(:)
 avt(:,:)=dt*avt(:,:)
 ert(:,:)=dt*ert(:,:)
 open(1,file='aa.dat',status='replace')
 do i=1,nt
    write(1,2)(i-1)*dt,avt(1,i),ert(1,i),avt(2,i),ert(2,i)
    2 format(i6,' ',4f15.8)
 end do
 close(1)

 write(*,*)'Integrated autocorrelation times'
 write(*,3)' E   : ',aint(1),eint(1)
 write(*,3)' |M| : ',aint(2),eint(2)
 3 format(a,2F12.4)

 end program autoaverage
!-----------------------!


