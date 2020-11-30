!-------------!
 module system
!-------------!

 integer :: nn
 integer :: nst

 integer, allocatable :: state(:)

 real(8), allocatable :: mat(:,:)
 real(8), allocatable :: vec(:,:)
 real(8), allocatable :: enr(:)
 real(8), allocatable :: spn(:)

 end module system
!-----------------!

!================!
 program hchain_m
!================!
 use system; implicit none

 integer :: nu,mst

 open(10,file='read.in',status='old')
 read(10,*)nn
 close(10)

 mst=1
 do nu=(nn+1)/2+1,nn
    mst=mst*nu
 enddo
 do nu=2,nn/2
    mst=mst/nu
 enddo
 
 allocate(state(mst))

 do nu=mod(nn,2),nn/2

    call makebasis(nu)

    allocate(mat(nst,nst))
    allocate(vec(nst,nst))
    allocate(enr(nst))
    allocate(spn(nst))

    call hamiltonian()
    call diagonalize(nst,mat,vec,enr)

    call spinsquared(nu)
    call transform(nst,mat,vec,spn)
    spn(:)=0.5d0*abs(sqrt(1.d0+4.d0*spn(:))-1.d0)

    call writedata(nu)

    deallocate(mat)
    deallocate(vec)
    deallocate(enr)
    deallocate(spn)

 enddo

 deallocate(state)

 end program hchain_m
!====================!

!------------------------!
 subroutine writedata(nu)
!------------------------!
 use system; implicit none

 integer :: i,nu
 
 character (len=20) :: filename
 write (filename,"('eig_',i0,'_',i0,'.csv')") nu,nst

 open(10,file=filename,status='unknown')

 WRITE(10,*) """X"",""Y"",""Z"""

 do i=1,nst
    write(10,*) i-1,",",anint(enr(i)*10)/10,",",spn(i)
 enddo
 close(10)

 end subroutine writedata
!------------------------!

!------------------------!
 subroutine makebasis(nu)
!------------------------!
 use system; implicit none

 integer :: i,s,n1,nu

 nst=0
 do s=0,2**nn-1
    n1=0
    do i=0,nn-1
       if (btest(s,i)) n1=n1+1
    enddo
    if (n1==nu) then
       nst=nst+1
       state(nst)=s
    endif
 enddo

 end subroutine makebasis
!------------------------!

!------------------------!
 subroutine hamiltonian()
!------------------------!
 use system; implicit none

 integer :: i,j,a,b,sa,sb

 mat(:,:)=0.d0
 do a=1,nst
    sa=state(a)
    do i=0,nn-1
       j=mod(i+1,nn)
       if (btest(sa,i).eqv.btest(sa,j)) then 
          mat(a,a)=mat(a,a)+0.25d0 
       else
          mat(a,a)=mat(a,a)-0.25d0                
          sb=ieor(sa,2**i+2**j)
          call findstate(sb,b)
          mat(a,b)=0.5d0    
       endif
    enddo
 enddo

 end subroutine hamiltonian
!--------------------------!

!--------------------------!
 subroutine spinsquared(nu)
!--------------------------!
 use system; implicit none

 integer :: i,j,a,b,sa,sb,nu
 

 mat(:,:)=0.d0
 do a=1,nst
    sa=state(a)
    mat(a,a)=(dfloat(nu)-dfloat(nn)/2.d0)**2+dfloat(nn)/2.d0
    do i=0,nn-1
    do j=i+1,nn-1
       if (btest(sa,i).neqv.btest(sa,j)) then 
          sb=ieor(sa,2**i+2**j)
          call findstate(sb,b)
          mat(a,b)=1.d0    
       endif
    enddo
    enddo
 enddo

 end subroutine spinsquared
!--------------------------!

!--------------------------!
 subroutine findstate(sa,a)
!--------------------------!
 use system; implicit none

 integer :: a,sa,amin,amax

 amin=1; amax=nst
 do 
    a=amin+(amax-amin)/2
    if (sa<state(a)) then
       amax=a-1
    elseif (sa>state(a)) then
       amin=a+1
    else
       return
    endif
 enddo

 end subroutine findstate
!-------------------------!

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
