!-------------!
 module system
!-------------!

 integer :: nn
 integer :: nst

 real(8), allocatable :: mat(:,:)
 real(8), allocatable :: vec(:,:)
 real(8), allocatable :: enr(:)
 real(8), allocatable :: spn(:)
 real(8), allocatable :: mag(:)

 end module system
!-----------------!

!================!
 program hchain_0
!================!
 use system; implicit none

 open(10,file='read.in',status='old')
 read(10,*)nn
 close(10)

 nst=2**nn
 allocate(mat(0:nst-1,0:nst-1))
 allocate(vec(0:nst-1,0:nst-1))
 allocate(enr(0:nst-1))
 allocate(spn(0:nst-1))
 allocate(mag(0:nst-1))

 call hamiltonian()
 call diagonalize(nst,mat,vec,enr)

 call spinsquared()
 call transform(nst,mat,vec,spn)
 spn(:)=0.5d0*abs(sqrt(1.d0+4.d0*spn(:))-1.d0)

 call magnetization()

 call writedata()

 deallocate(mat)
 deallocate(vec)
 deallocate(enr)
 deallocate(spn)
 deallocate(mag)

 end program hchain_0
!====================!

!----------------------!
 subroutine writedata()
!----------------------!
 use system; implicit none

 integer :: i

 open(10,file='eig.dat',status='unknown')
 do i=0,nst-1
    write(10,'(i5,3f18.10)')i,enr(i),spn(i),mag(i)
 enddo
 close(10)

 end subroutine writedata
!------------------------!

!------------------------!
 subroutine hamiltonian()
!------------------------!
 use system; implicit none

 integer :: i,j,a,b

 mat(:,:)=0.d0
 do a=0,nst-1
    do i=0,nn-1
       j=mod(i+1,nn)
       if (btest(a,i).eqv.btest(a,j)) then 
          mat(a,a)=mat(a,a)+0.25d0 
       else
          mat(a,a)=mat(a,a)-0.25d0                
          b=ieor(a,2**i+2**j)
          mat(a,b)=mat(a,b)+0.5d0    
       endif
    enddo
 enddo


 end subroutine hamiltonian
!--------------------------!

!------------------------!
 subroutine spinsquared()
!------------------------!
 use system; implicit none

 integer :: i,j,a,b,m

 mat(:,:)=0.d0
 do a=0,nst-1
    m=0
    do i=0,nn-1
       if (btest(a,i)) m=m+1
    enddo
    mat(a,a)=dfloat(m-nn/2)**2+0.5d0*dfloat(nn)
    do i=0,nn-1   
    do j=i+1,nn-1
       if (btest(a,i).neqv.btest(a,j)) then 
          b=ieor(a,2**i+2**j)
          mat(a,b)=mat(a,b)+1.d0    
       endif
    enddo
    enddo
 enddo

 end subroutine spinsquared
!--------------------------!

!--------------------------!
 subroutine magnetization()
!--------------------------!
 use system 
 implicit none

 integer :: i,j,a,b
 integer, allocatable :: mz(:)

 allocate(mz(0:nst-1))
 do a=0,nst-1
    mz(a)=0
    do i=0,nn-1
       if (btest(a,i)) mz(a)=mz(a)+1
    enddo
 enddo
 do a=0,nst-1
    mag(a)=0.d0
    do b=0,nst-1
       mag(a)=mag(a)+mz(b)*vec(b,a)**2
    enddo
 enddo
 mag(:)=(mag(:)-nn/2)*0.5d0
 deallocate(mz)

 end subroutine magnetization
!----------------------------!

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
