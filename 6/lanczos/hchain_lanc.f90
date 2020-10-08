!--------------!
 module system
!--------------!
 save

 integer :: nn
 integer :: nrep
 integer :: smax

 integer(4), allocatable :: repr(:)
 integer(1), allocatable :: peri(:)
 integer(2), allocatable :: mntr(:)

 real(8), allocatable :: sink(:)
 real(8), allocatable :: cosk(:)

 end module system
!-----------------!

!---------------------!
 module lanczosvectors
!---------------------!
 save

 real(8), allocatable :: p0(:)
 real(8), allocatable :: ff(:,:)
 real(8), allocatable :: aal(:)
 real(8), allocatable :: nnl(:)

 integer :: nham
 integer :: nval
 integer(1), allocatable :: nhab(:)
 integer(2), allocatable :: vhab(:)
 integer(4), allocatable :: bhab(:)
 real(8),    allocatable :: hval(:)

 end module lanczosvectors
!-------------------------!

!-------------------!
 program hchain_lanc
!-------------------!
 use system; use lanczosvectors; implicit none

 integer :: i,k,p,z,rm,hm,vm,ne,ns,nc,mlanc
 real(8) :: ee,ss

 real(8), allocatable :: eig(:)
 real(8), allocatable :: vec(:,:)

 open(10,file='read.in',status='old')
 read(10,*)nn,k,p,z
 read(10,*)mlanc,rm,hm,vm
 close(10)

 allocate(repr(rm))
 allocate(peri(rm))
 allocate(mntr(rm))
 allocate(nhab(rm))
 allocate(bhab(hm))
 allocate(vhab(hm))
 allocate(hval(vm))

 allocate(aal(0:mlanc-1))
 allocate(nnl(0:mlanc-1))
 allocate(eig(0:mlanc-1))
 allocate(vec(0:mlanc-1,0:mlanc-1))
 
 open(10,file='log.txt',position='append')
 write(10,*)'(k,p,z) = ',k,p,z
 close(10) 

 call trigfunctions(k)
 call makebasis(k,p,z)

 open(10,file='log.txt')
 write(10,*)'Basis size, max  ',nrep,rm
 close(10)
       
 call hamiltonian(p,z)

 open(10,file='log.txt',position='append')
 write(10,*)'nham, max = ',nham,size(vhab)
 write(10,*)'nval, max = ',nval,size(hval)
 close(10)

 allocate(p0(nrep))
 allocate(ff(nrep,0:mlanc))

 call initstate(nrep,1) 
 call lanczos1(mlanc,eig,vec)
 call diatridiag(mlanc,eig,vec)
 open(10,file='e.dat',position='append')
 do i=0,mlanc
    write(10,'(i6,f16.10)')i,eig(i)
 enddo
 close(10)
       
 deallocate(p0)
 deallocate(ff)
 deallocate(aal)
 deallocate(nnl)
 deallocate(eig)
 deallocate(vec)
 deallocate(nhab)
 deallocate(bhab)
 deallocate(vhab)
 deallocate(hval)
 deallocate(repr)
 deallocate(peri)
 deallocate(mntr)
 deallocate(sink)
 deallocate(cosk)

 end program hchain_lanc
!-----------------------!

!---------------------------!
 subroutine makebasis(k,p,z)
!---------------------------!
 use system; implicit none

 integer :: i,a,k,p,z,s,m,n,na,ca,ra,tp,tz,tpz
 logical :: pass

 na=0
 do i=0,nn-2
    na=ibset(na,i)
 enddo  
 smax=ibset(na,nn-1)
 nrep=0
 do i=0,na-1
    if (nn<32) then
       a=i
    else
       a=ibset(i,31)
    endif
    call checkstate(a,k,ra,tp,tz,tpz,pass)
    if (pass) then 
       do s=-1,1,2
          if ((k==0.or.k==nn/2).and.s==-1) then
             cycle
          elseif (tp==nn.and.tz==nn.and.tpz==nn) then
             ca=1
          elseif (tp/=nn.and.tz==nn) then
             ca=2; m=tp
             if ((1.d0+s*p*cosk(m))<1.d-8) cycle 
             if (s==-1.and.(1.d0-s*p*cosk(m))>1.d-8) cycle 
          elseif (tp==nn.and.tz/=nn) then
             ca=3; n=tz
             if (1.d0+z*cosk(n)<1.d-8) cycle 
          elseif (tp==nn.and.tz==nn) then      
             ca=4; m=tpz
             if ((1.d0+s*p*z*cosk(m))<1.d-8) cycle 
             if (s==-1.and.(1.d0-s*p*z*cosk(m))>1.d-8) cycle 
          elseif (tp/=nn.and.tz/=nn) then
             ca=5; m=tp; n=tz
             if ((1.d0+z*cosk(n))*(1.d0+s*p*cosk(m))<1.d-8) cycle 
             if (s==-1.and.(1.d0-s*p*cosk(m))>1.d-8) cycle 
          endif
          nrep=nrep+1
          repr(nrep)=a
          peri(nrep)=2*ra+(s+1)/2
          mntr(nrep)=m+n*(nn+1)+ca*(nn+1)**2
       enddo
    endif
    if (mod(i,1000000)==0) then
       open(10,file='log.txt')
       write(10,'(f10.6,i10)')dble(i)/dble(na),nrep
       close(10)
    endif
 enddo

 end subroutine makebasis
!------------------------!

!----------------------------------!
 subroutine lanczos1(mlanc,eig,vec)
!----------------------------------!
 use lanczosvectors; implicit none

 integer :: i,m,mlanc
 real(8) :: q1,q2,eig(0:mlanc-1),vec(0:mlanc-1,0:mlanc-1)

 ff(:,0)=p0
 call hoperation(0)
 aal(0)=dot_product(ff(:,0),ff(:,1))
 ff(:,1)=ff(:,1)-aal(0)*ff(:,0)
 nnl(1)=sqrt(dot_product(ff(:,1),ff(:,1)))
 ff(:,1)=ff(:,1)/nnl(1)

 do m=2,mlanc

    call hoperation(m-1)
    aal(m-1)=dot_product(ff(:,m-1),ff(:,m))
    ff(:,m)=ff(:,m)-aal(m-1)*ff(:,m-1)-nnl(m-1)*ff(:,m-2)
    nnl(m)=sqrt(dot_product(ff(:,m),ff(:,m)))
    ff(:,m)=ff(:,m)/nnl(m)

    do i=0,m-1
       q1=dot_product(ff(:,m),ff(:,i))
       q2=1.d0/sqrt(1.d0-q1**2)
       ff(:,m)=q2*(ff(:,m)-q1*ff(:,i))
    enddo

    call diatridiag(m,eig,vec)
    open(10,file='itr.dat',status='unknown',position='append')
    write(10,'(i5,5f16.10)')m,eig(0:4)
    close(10)

 end do

 end subroutine lanczos1
!-----------------------!

!------------------------!
 subroutine hoperation(m)
!------------------------!
 use lanczosvectors; implicit none

 integer :: m,i,a,b,c

 ff(:,m+1)=0.d0; i=0
 do a=1,size(p0)
    do c=1,nhab(a)
       i=i+1
       b=bhab(i)
       ff(b,m+1)=ff(b,m+1)+ff(a,m)*hval(vhab(i))
       if (b/=a) ff(a,m+1)=ff(a,m+1)+ff(b,m)*hval(vhab(i))
    enddo
 enddo

 end subroutine hoperation
!-------------------------!

!---------------------------!
 subroutine hamiltonian(p,z)
!---------------------------!
 use system; use lanczosvectors; implicit none

 integer :: i,j,a,b,c,p,z,l,q,g,v,ia,ib,sa,sb,bb,na,nb,nh(2)
 real(8) :: matelement,hh
 logical :: pass

 integer, allocatable :: bloc(:,:)
 integer, allocatable :: bval(:,:)
 allocate(bloc(2*(nn+1),2))
 allocate(bval(2*(nn+1),2))

 nham=0
 nval=0
 nhab(:)=0
 do ia=1,nrep
    sa=repr(ia)     
    if (ia>1.and.sa==repr(ia-1)) then 
       cycle
    elseif (ia<nrep.and.sa==repr(ia+1)) then
       na=2
    else  
       na=1
    endif
    hh=0.d0
    do i=0,nn-1
       j=mod(i+1,nn)
       if (btest(sa,i).eqv.btest(sa,j)) then
          hh=hh+0.25d0
       else
          hh=hh-0.25d0
       endif
    enddo
    do v=1,nval
       if (abs(hval(v)-hh)<1.d-8) exit
    enddo
    if (v>nval) then
       nval=nval+1
       hval(nval)=hh
    endif
    do a=ia,ia+na-1
       c=a-ia+1
       bloc(1,c)=a
       bval(1,c)=v
       nh(c)=1
    enddo
    do i=0,nn-1
       j=mod(i+1,nn)
       if (btest(sa,i).neqv.btest(sa,j)) then
          if (btest(sa,i)) then
             bb=ibset(ibclr(sa,i),j)
          else
             bb=ibset(ibclr(sa,j),i)
          endif
          call representative(bb,sb,l,q,g)
          call findstate(sb,ib)
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
                c=a-ia+1
                do b=ib,ib+nb-1
                   if (b<=a) then
                      hh=0.5d0*matelement(a,b,p,z,l,q,g)
                      do v=1,nval
                         if (abs(hval(v)-hh)<1.d-8) exit
                      enddo
                      if (v>nval) then
                         nval=nval+1
                         hval(nval)=hh
                      endif
                      nh(c)=nh(c)+1
                      bloc(nh(c),c)=b
                      bval(nh(c),c)=v
                   endif
                enddo
             enddo
          endif
       endif
    enddo   
    do a=ia,ia+na-1
       c=a-ia+1
       bhab(nham+1:nham+nh(c))=bloc(1:nh(c),c)         
       vhab(nham+1:nham+nh(c))=bval(1:nh(c),c)         
       nhab(a)=nh(c)
       nham=nham+nh(c)       
    enddo       
 enddo

 deallocate(bloc)
 deallocate(bval)

 end subroutine hamiltonian
!--------------------------!
   
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
    at=ishftc(at,-1,nn)
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
    at=ishftc(at,-1,nn)
    az=smax-at
 enddo 
 pass=.true.

 end subroutine checkstate
!-------------------------!

!-------------------------------------!
 subroutine representative(a0,a,l,q,g)
!-------------------------------------!
 use system; implicit none

 integer :: i,t,l,q,g,a,a0,at,az

 at=a0; a=a0; l=0; q=0; g=0
 do t=1,nn-1
    at=ishftc(at,-1,nn)
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
    at=ishftc(at,-1,nn)
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

!----------------------------!
 subroutine initstate(nrep,w)
!----------------------------!
  use lanczosvectors; use random; implicit none

 integer :: i,w,nrep

 call initran(w)
 do i=1,nrep
    p0(i)=ran()-0.5d0
 end do
 p0(:)=p0(:)/sqrt(dot_product(p0,p0))
 
 end subroutine initstate
!------------------------!

!---------------------------!
 subroutine trigfunctions(k)
!---------------------------!
 use system; implicit none

 real(8), parameter :: pi=3.14159265358979d0

 integer :: r,k
 real(8) :: kk

 kk=dfloat(2*k)*pi/dfloat(nn)
 if (.not.allocated(sink)) then
    allocate(sink(-nn:nn))
    allocate(cosk(-nn:nn))
 endif

 do r=-nn,nn
    sink(r)=sin(kk*r)
    cosk(r)=cos(kk*r)
 enddo

 end subroutine trigfunctions
!----------------------------!

!--------------------------------!
 subroutine diatridiag(n,eig,vec)
!--------------------------------!
 use lanczosvectors; implicit none

 integer :: n,inf
 real(8) :: d(n),e(n),eig(n),vec(n,n),work(max(1,2*n-2))

 d(1:n)=aal(0:n-1)
 e(1:n)=nnl(1:n)
 call dstev('V',n,d,e,vec,n,work,inf)
 eig=d

 end subroutine diatridiag
!-------------------------!

!------------------------------------------!
 real(8) function matelement(a,b,p,z,l,q,g)
!------------------------------------------!
 use system; implicit none

 integer :: a,b,p,z,l,q,g,s,t,ca,cb,ma,mb,na,nb
 real(8) :: norma,normb

 ca=mntr(a)/(nn+1)**2
 cb=mntr(b)/(nn+1)**2
 na=(mntr(a)-ca*(nn+1)**2)
 nb=(mntr(b)-cb*(nn+1)**2)
 ma=mod(na,nn+1)
 mb=mod(nb,nn+1)
 na=na/(nn+1)
 nb=nb/(nn+1)
 s=2*mod(peri(a),2_1)-1
 t=2*mod(peri(b),2_1)-1

 norma=1.d0/dble(peri(a)/2)
 if (ca==2.or.ca==5) norma=norma*(1.d0+s*p*cosk(ma))
 if (ca==3.or.ca==5) norma=norma*(1.d0+z*cosk(na))
 if (ca==4) norma=norma*(1.d0+s*p*z*cosk(ma))
 normb=1.d0/dble(peri(b)/2)
 if (cb==2.or.cb==5) normb=normb*(1.d0+t*p*cosk(mb))
 if (cb==3.or.cb==5) normb=normb*(1.d0+z*cosk(nb))
 if (cb==4) normb=normb*(1.d0+t*p*z*cosk(mb))
 matelement=((s*p)**q)*(z**g)*dsqrt(normb/norma)
 if (cb==1.or.cb==3) then        
    if (s==t) then
       matelement=matelement*cosk(l)
    else
       matelement=matelement*sink(l)*t
    endif
 elseif (cb==2.or.cb==5) then
    if (s==t) then
       matelement=matelement*(cosk(l)+t*p*cosk(l-mb))/(1.d0+t*p*cosk(mb))
    else
       matelement=matelement*(t*sink(l)+p*sink(l-mb))/(1.d0+t*p*cosk(mb))
    endif
 elseif (cb==4) then        
    if (s==t) then
       matelement=matelement*(cosk(l)+t*p*z*cosk(l-mb))/(1.d0+t*p*z*cosk(mb))
    else
       matelement=matelement*(t*sink(l)+p*z*sink(l-mb))/(1.d0+t*p*z*cosk(mb))
    endif
 endif

 end function matelement
!-----------------------!
