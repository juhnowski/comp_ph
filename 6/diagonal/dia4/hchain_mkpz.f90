!--------------!
 module system
!--------------!
 implicit none

 integer :: nn
 integer :: nrep
 integer :: smax

 integer, allocatable :: repr(:)
 integer, allocatable :: type(:)
 integer, allocatable :: peri(:)
 integer, allocatable :: mtrf(:)
 integer, allocatable :: ntrf(:)

 real(8), allocatable :: mat(:,:)
 real(8), allocatable :: vec(:,:)
 real(8), allocatable :: enr(:)
 real(8), allocatable :: spn(:)

 real(8), allocatable :: ffun(:,:,:)
 real(8), allocatable :: gfun(:,:)

 end module system
!-----------------!

!===================!
 program hchain_mkpz
!===================!
 use system; implicit none

 integer :: k,p,z,p1,p2,rm

 open(10,file='read.in',status='old')
 read(10,*)nn,rm
 close(10)

 smax=2**nn-1 

 allocate(repr(rm))
 allocate(type(rm))
 allocate(peri(rm))
 allocate(mtrf(rm))
 allocate(ntrf(rm))

 do k=0,nn/2
    call fgfunctions(k)
    if (k==0.or.k==nn/2) then
       p1=-1; p2=+1
    else
       p1=+1; p2=+1
    endif
    do p=p2,p1,-2
    do z=+1,-1,-2
       call makebasis(k,p,z)
       if (nrep/=0) then
          allocate(mat(nrep,nrep))
          allocate(vec(nrep,nrep))
          allocate(enr(nrep))
          allocate(spn(nrep))
          call hamiltonian(p,z)
          call diagonalize(nrep,mat,vec,enr)
          call spinsquared(p,z)
          call transform(nrep,mat,vec,spn)
          spn(:)=0.5d0*abs(sqrt(1.d0+4.d0*spn(:))-1.d0)
          call writedata(k,p,z)          
          deallocate(mat)
          deallocate(vec)
          deallocate(enr)
          deallocate(spn)
       endif
    enddo
    enddo
 enddo

 deallocate(repr)
 deallocate(type)
 deallocate(peri)
 deallocate(mtrf)
 deallocate(ntrf)
 deallocate(ffun)
 deallocate(gfun)
 
 end program hchain_mkpz
!=======================!

!---------------------------!
 subroutine writedata(k,p,z)
!---------------------------!
 use system; implicit none

 integer :: i,k,p,z
 character (len=20) :: eig_filename
 
 write (eig_filename,"('eig_',i0,'_',i0,'_',i0,'_',i0,'.csv')") k,p,z,nrep

 open(10,file=eig_filename,status='unknown')
 WRITE(10,*) """X"",""Y"",""Z"""
 do i=1,nrep
    write(10,*)i,",",anint(enr(i)*10)/10,",",spn(i)
 end do
 close(10)

 !open(10,file='eig.dat',position='append')
 !write(10,10)k,p,z,nrep
 !10 format(3i5,i10)
 !do i=1,nrep
 !   write(10,20)i,enr(i),spn(i)
 !   20 format(i5,'  ',2f16.10)
 !enddo
 !close(10)

 !open(10,file='low.dat',position='append')
 !write(10,30)k,p,z,enr(1),spn(1),nrep
 !30 format(3i5,2f16.10,i10)
 !close(10)

 end subroutine writedata
!------------------------!

!---------------------------!
 subroutine hamiltonian(p,z)
!---------------------------!
 use system; implicit none

 integer :: i,j,a,b,p,z,l,q,g,ia,ib,sa,sb,na,nb
 real(8) :: matelement,ez

 mat(:,:)=0.d0
 do ia=1,nrep
    sa=repr(ia)     
    if (ia>1.and.sa==repr(ia-1)) then 
       cycle
    elseif (ia<nrep.and.sa==repr(ia+1)) then
       na=2
    else  
       na=1
    endif
    ez=0.d0
    do i=0,nn-1
       j=mod(i+1,nn)
       if (btest(sa,i).eqv.btest(sa,j)) then
          ez=ez+0.25d0
       else
          ez=ez-0.25d0
       endif
    enddo
    do a=ia,ia+na-1
       mat(a,a)=mat(a,a)+ez
    enddo  
    do i=0,nn-1
       j=mod(i+1,nn)
       if (btest(sa,i).neqv.btest(sa,j)) then
          if (btest(sa,i)) then
             sb=ibset(ibclr(sa,i),j)
          else
             sb=ibset(ibclr(sa,j),i)
          endif
          call representative(sb,b,l,q,g)
          call findstate(b,ib)
          if (ib>0) then
             if (ib>1.and.repr(ib)==repr(ib-1)) then
                ib=ib-1
                nb=2
             elseif (ib<nrep.and.repr(ib)==repr(ib+1)) then
                nb=2
            else
                nb=1
             endif
             do a=ia,ia+na-1
             do b=ib,ib+nb-1
                mat(a,b)=mat(a,b)+0.5d0*matelement(a,b,p,z,l,q,g)
             enddo
             enddo
          endif
       endif
    enddo
 enddo

 end subroutine hamiltonian
!--------------------------!

!---------------------------!
 subroutine spinsquared(p,z)
!---------------------------!
 use system; implicit none

 integer :: i,j,a,b,p,z,l,q,g,ia,ib,sa,sb,na,nb
 real(8) :: matelement

 mat(:,:)=0.d0
 do ia=1,nrep
    sa=repr(ia)    
    if (ia>1.and.sa==repr(ia-1)) then
       cycle
    elseif (ia<nrep.and.sa==repr(ia+1)) then
       na=2
    else  
       na=1
    endif
    do a=ia,ia+na-1
       mat(a,a)=mat(a,a)+0.5d0*dfloat(nn)
    enddo
    do i=0,nn-1
    do j=i+1,nn-1
       if (btest(sa,i).neqv.btest(sa,j)) then
          if (btest(sa,i)) then
             sb=ibset(ibclr(sa,i),j)
          else
             sb=ibset(ibclr(sa,j),i)
          endif
          call representative(sb,b,l,q,g)
          call findstate(b,ib)
          if (ib>0) then
             if (ib>1.and.repr(ib)==repr(ib-1)) then
                ib=ib-1
                nb=2
             elseif (ib<nrep.and.repr(ib)==repr(ib+1)) then
                nb=2
             else
                nb=1
             endif
             do a=ia,ia+na-1
             do b=ib,ib+nb-1
                mat(a,b)=mat(a,b)+matelement(a,b,p,z,l,q,g)
             enddo
             enddo
          endif
       endif
    enddo
    enddo
 enddo

 end subroutine spinsquared
!--------------------------!

!------------------------------------------!
 real(8) function matelement(a,b,p,z,l,q,g)
!------------------------------------------!
 use system; implicit none

 integer :: a,b,p,z,l,q,g,s,t,ca,cb

 ca=type(a)/2
 cb=type(b)/2
 s=2*mod(type(a),2)-1
 t=2*mod(type(b),2)-1
 matelement=dfloat(peri(a))/dfloat(peri(b))
 if (ca==2.or.ca==5) matelement=matelement/gfun(s*p,mtrf(a))
 if (ca==3.or.ca==5) matelement=matelement/gfun(z,ntrf(a))
 if (ca==4) matelement=matelement/gfun(s*p*z,mtrf(a))
 if (cb==2.or.cb==5) matelement=matelement*gfun(t*p,mtrf(b))
 if (cb==3.or.cb==5) matelement=matelement*gfun(z,ntrf(b))
 if (cb==4) matelement=matelement*gfun(t*p*z,mtrf(b))
 matelement=((s*p)**q)*(z**g)*dsqrt(matelement)
 if (cb==1.or.cb==3) then        
    matelement=matelement*ffun(t,s,l)
 elseif (cb==2.or.cb==5) then
    matelement=matelement*(ffun(t,s,l)+t*p*ffun(t,s,l-mtrf(b)))/gfun(t*p,mtrf(b))
 elseif (cb==4) then        
    matelement=matelement*(ffun(t,s,l)+t*p*z*ffun(t,s,l-mtrf(b)))/gfun(t*p*z,mtrf(b))
 endif
 
 end function matelement
!-----------------------!

!---------------------------!
 subroutine makebasis(k,p,z)
!---------------------------!
 use system; implicit none

 integer :: a,k,p,z,s,m,n,ca,ra,tp,tz,tpz
 logical :: pass

 nrep=0
 do a=0,2**nn-1
    call checkstate(a,k,ra,tp,tz,tpz,pass)
    if (pass) then 
       do s=-1,1,2
          if ((k==0.or.k==nn/2).and.s==-1) then
             pass=.false.
          elseif (tp==nn.and.tz==nn.and.tpz==nn) then
             ca=1
          elseif (tp/=nn.and.tz==nn) then
             ca=2; m=tp
             if (gfun(s*p,m)<1.d-8) pass=.false.
             if (s==-1.and.gfun(-s*p,m)>1.d-8) pass=.false.
          elseif (tp==nn.and.tz/=nn) then
             ca=3; n=tz
             if (gfun(z,n)<1.d-8) pass=.false.
          elseif (tp==nn.and.tz==nn) then      
             ca=4; m=tpz
             if (gfun(s*p*z,m)<1.d-8) pass=.false.
             if (s==-1.and.gfun(-s*p*z,m)>1.d-8) pass=.false.
          elseif (tp/=nn.and.tz/=nn) then
             ca=5; m=tp; n=tz
             if (gfun(z,n)*gfun(s*p,m)<1.d-8) pass=.false.
             if (s==-1.and.gfun(-s*p,m)>1.d-8) pass=.false.             
          endif
          if (pass) then
             nrep=nrep+1
             repr(nrep)=a
             type(nrep)=2*ca+(s+1)/2
             peri(nrep)=ra
             mtrf(nrep)=m
             ntrf(nrep)=n
          endif
          pass=.true.
       enddo
    endif
 enddo

 end subroutine makebasis
!------------------------!
   
!---------------------------------------------!
 subroutine checkstate(a0,k,ra,tp,tz,tpz,pass)
!---------------------------------------------!
 use system; implicit none

 integer :: i,k,t,m,a0,at,az,ra,tp,tz,tpz
 logical :: pass

 pass=.false.
 m=0
 do i=0,nn-1
    if (btest(a0,i)) m=m+1
 enddo
 if (m/=nn/2) return
 ra=nn
 tz=nn
 at=a0
 do t=1,nn-1
    if (btest(at,0)) then
       at=at/2
       at=ibset(at,nn-1)
    else
       at=at/2
    endif
    az=smax-at
    if (at<a0.or.az<a0) then
       return
    elseif (at==a0) then
       if (mod(k,nn/t)/=0) return
       ra=t
       exit
    elseif (az==a0) then
       tz=t
    endif
 enddo 
 tp=nn
 tpz=nn
 at=0
 do i=0,nn-1
    if (btest(a0,i)) at=ibset(at,nn-1-i)
 enddo
 az=smax-at
 do t=0,ra-1
    if (at<a0.or.az<a0) return
    if (at==a0) tp=t      
    if (az==a0) tpz=t      
    if (btest(at,0)) then
       at=at/2
       at=ibset(at,nn-1)
    else
       at=at/2
    endif
    az=smax-at
 enddo 
 pass=.true.

 end subroutine checkstate
!-------------------------!

!-------------------------------------!
 subroutine representative(a0,a,l,q,g)
!-------------------------------------!
 use system; implicit none

 integer :: i,a,t,l,q,g,a0,at,az

 at=a0; a=a0; l=0; q=0; g=0
 do t=1,nn-1
    if (btest(at,0)) then
       at=at/2; at=ibset(at,nn-1)
    else
       at=at/2
    endif
    if (at<a) then
       a=at; l=t; g=0
    endif
    az=smax-at
    if (az<a) then
       a=az; l=t; g=1
    endif
 enddo
 at=0
 do i=0,nn-1
    if (btest(a0,i)) at=ibset(at,nn-1-i)
 enddo
 do t=0,nn-1
    if (at<a) then
       a=at; l=t; q=1; g=0
    endif
    az=smax-at
    if (az<a) then
       a=az; l=t; q=1; g=1
    endif
    if (btest(at,0)) then
       at=at/2; at=ibset(at,nn-1)
    else
       at=at/2
    endif
 enddo 

 end subroutine representative
!-----------------------------!

!--------------------------!
 subroutine findstate(sa,a)
!--------------------------!
 use system; implicit none

 integer :: a,sa,amin,amax

 amin=1; amax=nrep
 do 
    a=amin+(amax-amin)/2
    if (sa<repr(a)) then
       amax=a-1
    elseif (sa>repr(a)) then
       amin=a+1
    else
       return
    endif
    if (amin>amax) then
       a=-1
       return
    endif
 enddo

 end subroutine findstate
!------------------------!

!-------------------------!
 subroutine fgfunctions(k)
!-------------------------!
 use system; implicit none

 real(8), parameter :: pi=3.14159265358979d0

 integer :: i,k
 real(8) :: kk

 kk=dfloat(2*k)*pi/dfloat(nn)
 if (.not.allocated(ffun)) then
    allocate(ffun(-1:1,-1:1,-nn:nn))
    allocate(gfun(-1:1,-nn:nn))
 endif
 do i=-nn,nn
    ffun(-1,-1,i)=+cos(kk*i)
    ffun(+1,+1,i)=+cos(kk*i)
    ffun(+1,-1,i)=+sin(kk*i)
    ffun(-1,+1,i)=-sin(kk*i)
    gfun(-1,i)=1.d0-cos(kk*i)
    gfun(+1,i)=1.d0+cos(kk*i)
 enddo

 end subroutine fgfunctions
!--------------------------!

!-------------------------------------!
 subroutine diagonalize(n,mat,vec,eig)
!-------------------------------------!
 implicit none

 integer :: n,info
 real(8) :: mat(n,n),vec(n,n),eig(n),work(n*(3+n/2))

 vec=mat
 call dsyev('V','U',n,vec,n,eig,work,n*(3+n/2),info)

 end subroutine diagonalize
!--------------------------!

!-----------------------------------!
 subroutine transform(n,mat,vec,dia)
!-----------------------------------!
 implicit none

 integer :: i,n
 real(8) :: mat(n,n),vec(n,n),dia(n)

 mat=matmul(mat,vec)
 mat=matmul(transpose(vec),mat)
 do i=1,n
    dia(i)=mat(i,i)
 enddo

 end subroutine transform
!------------------------!

