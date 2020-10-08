 program hchain_temp

 implicit none

 integer :: i,j,k,p,z,n,s,nn,nt,nk,ne
 real(8) :: e0,ee,ss,w,t,dt

 real(8), allocatable :: et(:),ct(:),xt(:),zt(:)

 write(*,*)'System size N'
 read(*,*)nn
 write(*,*)'Number of temperatures, Delta-T'
 read(*,*)nt,dt

 allocate(et(nt)); et=0.d0
 allocate(ct(nt)); ct=0.d0
 allocate(xt(nt)); xt=0.d0
 allocate(zt(nt)); zt=0.d0

 open(10,file='eig.dat',status='old')
 e0=0.d0
 do
    read(10,*,end=10)k,p,z,ne
    do n=1,ne
       read(10,*)i,ee,ss
       if (ee<e0) e0=ee
    enddo
 enddo
 10 rewind(10)
 do
    read(10,*,end=20)k,p,z,ne
    if (k==0.or.k==nn/2) then
       nk=1
    else
       nk=2
    endif
    do n=1,ne
       read(10,*)i,ee,ss
       s=int(ss+1.d-6)
       if (abs(ss-dfloat(s))>1.d-8) then
          write(*,'(a,3i3,a,2f14.8)')' k,p,z : ',k,p,z,'    E, S : ',ee,ss
       endif
       do i=1,nt
          t=dt*i
          w=nk*(2*s+1)*dexp(-(ee-e0)/t)
          zt(i)=zt(i)+w
          et(i)=et(i)+w*ee
          ct(i)=ct(i)+w*ee**2
          xt(i)=xt(i)+w*ss*(ss+1.d0)/3.d0
       enddo
    enddo
 enddo
 20 close(10)

 open(10,file='t.dat',status='unknown')
 do i=1,nt
    t=dt*i
    et(i)=et(i)/zt(i)
    ct(i)=(ct(i)/zt(i)-et(i)**2)/(nn*t**2)
    xt(i)=(xt(i)/zt(i))/(nn*t)
    write(10,'(4f16.9)')t,et(i)/nn,ct(i),xt(i)
 enddo
 close(10)

 deallocate(et)
 deallocate(ct)
 deallocate(xt)
 deallocate(zt)

 end program hchain_temp
