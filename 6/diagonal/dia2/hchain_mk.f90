
!--------------!
 module system
!--------------!
 implicit none

 integer :: nn
 integer :: nrep

 integer, allocatable :: repr(:)
 integer, allocatable :: peri(:)

 real(8),    allocatable :: enr(:)
 real(8),    allocatable :: spn(:)

 complex(8), allocatable :: mat(:,:)
 complex(8), allocatable :: vec(:,:)
 complex(8), allocatable :: expk(:)

 end module system
!-----------------!

!=================!
 program hchain_mk
!=================!
 use system; implicit none

 integer :: i,k,nu,rm

 open(10,file='read.in',status='old')
 read(10,*)nn,rm
 close(10)
 
 allocate(repr(rm))
 allocate(peri(rm))

 do nu=0,nn/2
 do k=0,nn/2

    call expfunction(k)
    call makebasis(nu,k)

    if (nrep/=0) then

       allocate(mat(nrep,nrep))
       allocate(vec(nrep,nrep))
       allocate(enr(nrep))
       allocate(spn(nrep))

       call hamiltonian()
       call diagonalize(nrep,mat,vec,enr)

       call spinsquared(nu) 
       call transform(nrep,mat,vec,spn)
       spn(:)=0.5d0*abs(sqrt(1.d0+4.d0*spn(:))-1.d0)

       call writedata(nu,k)          

       deallocate(mat)
       deallocate(vec)
       deallocate(enr)
       deallocate(spn)

    endif

 enddo
 enddo

 deallocate(repr)
 deallocate(peri)
 deallocate(expk)
 
 end program hchain_mk
!=====================!

!--------------------------!
 subroutine writedata(nu,k)
!--------------------------!
 use system; implicit none

 integer :: i,k,nu

 character (len=20) :: eig_filename
 character (len=20) :: spn_filename
 character (len=20) :: nrep_filename

 write (eig_filename,"('eig_',i0,'_',i0,'_',i0,'.csv')") nu,k,nrep

 open(10,file=eig_filename,status='unknown')
 WRITE(10,*) """X"",""Y"",""Z"""
 do i=1,nrep
    write(10,*)i-1,",",anint(enr(i)*10)/10,",",spn(i)
 end do
 close(10)

 end subroutine writedata
!------------------------!

!------------------------!
 subroutine hamiltonian()
!------------------------!
 use system; implicit none

 integer :: i,j,a,b,l,sa,sb,bb

 mat(:,:)=0.d0
 do a=1,nrep
    sa=repr(a)    
    do i=0,nn-1
       j=mod(i+1,nn)
       if (btest(sa,i).eqv.btest(sa,j)) then
          mat(a,a)=mat(a,a)+0.25d0
       else
          mat(a,a)=mat(a,a)-0.25d0
          bb=ieor(sa,2**i+2**j)
          call representative(bb,sb,l)
          call findstate(sb,b)
          if (b>=0) then
             mat(a,b)=mat(a,b)+0.5d0*sqrt(dfloat(peri(a))/dfloat(peri(b)))*expk(l)
          endif
       endif
    enddo
 enddo

 end subroutine hamiltonian
!--------------------------!

!--------------------------!
 subroutine spinsquared(nu)
!--------------------------!
 use system; implicit none

 integer :: i,j,a,b,l,nu,sa,sb,bb

 mat(:,:)=0.d0
 do a=1,nrep
    sa=repr(a)    
    mat(a,a)=mat(a,a)+(dfloat(nu)-dfloat(nn)/2.d0)**2+dfloat(nn)/2.d0
    do i=0,nn-1
    do j=i+1,nn-1
       if (btest(sa,i).neqv.btest(sa,j)) then
          bb=ieor(sa,2**i+2**j)
          call representative(bb,sb,l)
          call findstate(sb,b)
          if (b>=0) mat(a,b)=mat(a,b)+sqrt(dfloat(peri(a))/dfloat(peri(b)))*expk(l)
       endif
    enddo
    enddo
 enddo

 end subroutine spinsquared
!--------------------------!

!--------------------------!
 subroutine makebasis(nu,k)
!--------------------------!
 use system; implicit none

 integer :: s,k,nu,ra
 logical :: pass

 nrep=0
 do s=0,2**nn-1
    call checkstate(s,nu,k,ra,pass)
    if (pass) then
       nrep=nrep+1
       repr(nrep)=s
       peri(nrep)=ra
    endif
 enddo

 end subroutine makebasis
!------------------------!
   
!--------------------------------------!
 subroutine checkstate(sa,nu,k,ra,pass)
!--------------------------------------!
 use system; implicit none

 integer :: i,t,k,nu,n1,ra,sa,at
 logical :: pass

 pass=.false.
 n1=0
 do i=0,nn-1
    if (btest(sa,i)) n1=n1+1
 enddo
 if (n1/=nu) return
 ra=nn
 at=sa
 do t=1,nn-1
    at=ishftc(at,-1,nn)
    if (at<sa) then
       return
    elseif (at==sa) then
       if (mod(k,nn/t)/=0) return
       ra=t
       exit
    endif
 enddo 
 pass=.true.

 end subroutine checkstate
!-------------------------!

!----------------------------------!
 subroutine representative(aa,sa,l)
!----------------------------------!
 use system; implicit none

 integer :: i,t,l,aa,sa,at

 sa=aa; at=aa; l=0
 do t=1,nn-1
    at=ishftc(at,-1,nn)
    if (at<sa) then
       sa=at; l=t
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
 subroutine expfunction(k)
!-------------------------!
 use system; implicit none

 real(8), parameter :: pi=3.14159265358979d0

 integer :: i,k
 real(8) :: kk

 kk=dfloat(2*k)*pi/dfloat(nn)
 if (.not.allocated(expk)) allocate(expk(-nn:nn))
 do i=-nn,nn
    expk(i)=exp(-(0,1)*kk*i)
 enddo

 end subroutine expfunction
!--------------------------!

!-------------------------------------!
 subroutine diagonalize(n,mat,vec,eig)
!-------------------------------------!
 implicit none

 integer :: n,ierr
 real(8) :: rz(n,n),iz(n,n),eig(n),fv1(n),fv2(n),fm1(2,n)
 complex(8) :: mat(n,n),vec(n,n)

 vec=mat
 call dch(n,n,real(vec),aimag(vec),eig,1,rz,iz,fv1,fv2,fm1,ierr)
 vec=rz+(0,1)*iz

 end subroutine diagonalize
!--------------------------!

!-----------------------------------!
 subroutine transform(n,mat,vec,dia)
!-----------------------------------!
 implicit none

 integer :: i,n
 real(8) :: dia(n)
 complex(8) :: mat(n,n),vec(n,n)

 mat=matmul(mat,vec)
 mat=matmul(conjg(transpose(vec)),mat)
 do i=1,n
    dia(i)=real(mat(i,i))
 enddo

 end subroutine transform
!------------------------!
