!---------------------------------------------------------!
!Diagonalization of a real symmetric tridiagonal matrix   !
!---------------------------------------------------------!
!input:  a(n)     = diagonal matrix elements              !
!        b(n)     = off-diagonal matrix elements          !
!        n        = size of the matrix                    !
!output: eig(n)   = eigenvalues in ascending order        !
!        vec(n,n) = eigenvectors                          !
!---------------------------------------------------------!
!Uses the LAPACK subroutine DSTEV                         !
!---------------------------------------------------------!
 subroutine diatri(a,b,eig,vec,n)
!-------------------------------!
 implicit none

 integer :: n,inf
 real(8) :: a(n),b(n),d(n),e(n),eig(n),vec(n,n),work(max(1,2*n-2))

 d=a
 e=b
 call dstev('V',n,d,e,vec,n,work,inf)
 eig=d

 end subroutine diatri
!---------------------!