!--------------!
 module system
!--------------!
 implicit none

 integer :: nn
 integer :: nrep

 integer, allocatable :: repr(:)
 integer, allocatable :: peri(:)
 integer, allocatable :: mtrf(:)

 real(8), allocatable :: mat(:,:)
 real(8), allocatable :: vec(:,:)
 real(8), allocatable :: enr(:)
 real(8), allocatable :: spn(:)

 real(8), allocatable :: ffun(:,:,:)
 real(8), allocatable :: gfun(:,:)

 end module system
!-----------------!

!==================!
 program hchain_mkp
!==================!
 use system; implicit none

 integer :: i,k,p,nu,rm,p1,p2

 open(10,file='read.in',status='old')
 read(10,*)nn,rm
 close(10)

 print *, nn, rm

 allocate(repr(rm))
 allocate(peri(rm))
 allocate(mtrf(rm))

 do nu=0,nn/2
 do k=0,nn/2
    if (k==0.or.k==nn/2) then
       p1=-1; p2=+1
    else
       p1=+1; p2=+1
    endif
   call fgfunctions(k)
    do p=p1,p2,2
       call makebasis(nu,k,p)
       if (nrep/=0) then
          allocate(mat(nrep,nrep))
          allocate(vec(nrep,nrep))
          allocate(enr(nrep))
          allocate(spn(nrep))
          call hamiltonian(p)
          call diagonalize(nrep,mat,vec,enr)
          call spinsquared(nu,p) 
          call transform(nrep,mat,vec,spn)
          spn(:)=0.5d0*abs(sqrt(1.d0+4.d0*spn(:))-1.d0)
          deallocate(mat)
          deallocate(vec)
          deallocate(enr)
          deallocate(spn)
       endif
       call writedata(nu,k,p)          
    enddo
 enddo
 enddo

 deallocate(repr)
 deallocate(peri)
 deallocate(mtrf)
 deallocate(ffun)
 deallocate(gfun)
 
 end program hchain_mkp
!======================!

!----------------------------!
 subroutine writedata(nu,k,p)
!----------------------------!
 use system; implicit none

 integer :: i,k,p,nu

 character (len=20) :: eig_filename

 if (nrep>0) then 

 write (eig_filename,"('eig_',i0,'_',i0,'_',i0,'.csv')") nu,k,p

 open(10,file=eig_filename,status='unknown')
 WRITE(10,*) """X"",""Y"",""Z"""
 do i=1,nrep
    write(10,*)i-1,",",anint(enr(i)*10)/10,",",spn(i)
 end do
 close(10)

end if

! open(10,file='eig.dat',position='append')
! write(10,'(a,i2,a,i2,a,i2,a,i4)')'nu =',nu,', k = ',k,', p = ',p,',  nst =',nrep
! do i=1,nrep
!    write(10,'(i5,2f18.10)')i-1,enr(i),spn(i)
! enddo
! close(10)

 !open(10,file='low.dat',position='append')

 !if (nrep/=0) then 
 !  write(10,30)nu,k,p,enr(1),spn(1),nrep
 !end if

 !30 format(3i5,2f16.10,i10)
 !close(10)

 end subroutine writedata
!------------------------!

!-------------------------!
 subroutine hamiltonian(p)
!-------------------------!
 use system; implicit none

 integer :: i,j,a,b,p,r,l,q,ia,ib,sa,sb,na,nb
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
          sb=ieor(sa,2**i+2**j)
          call representative(sb,b,l,q)
          call findstate(b,ib)
          if (ib>=0) then
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
                mat(a,b)=mat(a,b)+0.5d0*matelement(a,b,p,l,q)
             enddo
             enddo
          endif
       endif
    enddo
 enddo

 end subroutine hamiltonian
!--------------------------!

!----------------------------!
 subroutine spinsquared(nu,p)
!----------------------------!
 use system; implicit none

 integer :: i,j,a,b,p,r,l,q,nu,ia,ib,sa,sb,na,nb
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
       mat(a,a)=mat(a,a)+(dfloat(nu)-dfloat(nn)/2.d0)**2+dfloat(nn)/2.d0
    enddo
    do i=0,nn-1
    do j=i+1,nn-1
       if (btest(sa,i).neqv.btest(sa,j)) then
          sb=ieor(sa,2**i+2**j)
          call representative(sb,b,l,q)
          call findstate(b,ib)
          if (ib>=0) then
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
                mat(a,b)=mat(a,b)+matelement(a,b,p,l,q)
             enddo            
             enddo
          endif
       endif
    enddo
    enddo
 enddo

 end subroutine spinsquared
!--------------------------!

!--------------------------------------!
 real(8) function matelement(a,b,p,l,q)
!--------------------------------------!
 use system; implicit none

 integer :: a,b,p,l,q,s,t

 s=abs(peri(a))/peri(a)
 t=abs(peri(b))/peri(b)
 matelement=abs(dfloat(peri(a))/dfloat(peri(b)))
 if (mtrf(a)/=nn) matelement=matelement/gfun(p*s,mtrf(a))
 if (mtrf(b)/=nn) matelement=matelement*gfun(p*t,mtrf(b))
 matelement=((p*s)**q)*dsqrt(matelement)
 if (mtrf(b)==nn) then
    matelement=matelement*ffun(t,s,l)
 else
    matelement=matelement*(ffun(t,s,l)+t*p*ffun(t,s,l-mtrf(b)))/gfun(t*p,mtrf(b))
 endif
 
 end function matelement
!-----------------------!

!----------------------------!
 subroutine makebasis(nu,k,p)
!----------------------------!
 use system; implicit none

 integer :: a,i,k,p,s,n1,nu,ra,tp
 logical :: pass

 nrep=0
 do a=0,2**nn-1
    call checkstate(a,nu,k,p,ra,tp,pass)
    if (pass) then
       do s=-1,1,2
          if ((k==0.or.k==nn/2).and.s==-1) then
             pass=.false.
          elseif (tp/=nn) then
             if (gfun(s*p,tp)<1.d-8.or.(s==-1.and.gfun(-s*p,tp)>1.d-8)) pass=.false.
          endif
          if (pass) then
             nrep=nrep+1
             repr(nrep)=a
             peri(nrep)=s*ra
             mtrf(nrep)=tp
          endif
          pass=.true.
       enddo
    endif
 enddo

 end subroutine makebasis
!------------------------!
   
!-------------------------------------------!
 subroutine checkstate(a0,nu,k,p,ra,tp,pass)
!-------------------------------------------!
 use system; implicit none

 integer :: i,t,k,p,nu,n1,ra,tp,a0,at
 logical :: pass

 pass=.false.
 n1=0
 do i=0,nn-1
    if (btest(a0,i)) n1=n1+1
 enddo
 if (n1/=nu) return
 ra=nn
 at=a0
 do t=1,nn-1
    if (btest(at,0)) then
       at=at/2
       at=ibset(at,nn-1)
    else
       at=at/2
    endif
    if (at<a0) then
       return
    elseif (at==a0) then
       if (mod(k,nn/t)/=0) return
       ra=t
       exit
    endif
 enddo 
 tp=nn
 at=0
 do i=0,nn-1
    if (btest(a0,i)) at=ibset(at,nn-1-i)
 enddo
 do t=0,ra-1
    if (at<a0) return
    if (at==a0) tp=t      
    if (btest(at,0)) then
       at=at/2
       at=ibset(at,nn-1)
    else
       at=at/2
    endif
 enddo 
 pass=.true.

 end subroutine checkstate
!-------------------------!

!-----------------------------------!
 subroutine representative(a0,a,l,q)
!-----------------------------------!
 use system; implicit none

 integer :: i,t,a,l,q,a0,at

 at=a0; a=a0; l=0; q=0
 do t=1,nn-1
    if (btest(at,0)) then
       at=at/2; at=ibset(at,nn-1)
    else
       at=at/2
    endif
    if (at<a) then
       a=at; l=t
    endif
 enddo
 at=0
 do i=0,nn-1
    if (btest(a0,i)) at=ibset(at,nn-1-i)
 enddo
 do t=0,nn-1
    if (at<a) then
       a=at; l=t; q=1
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
!-------------------------!

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
