program test
    implicit none
    integer :: nn
    integer :: nst
    real(8), allocatable :: mat(:,:)
    integer :: i,j,a,b

    nn=2
    nst=2**nn
    allocate(mat(0:nst-1,0:nst-1))

    mat(:,:)=0.d0
    do a=0,nst-1
       do i=0,nn-1
          j=mod(i+1,nn)
          print *, "a=",a, "i=", i, "j=" , j

          if (btest(a,i).eqv.btest(a,j)) then 
             mat(a,a)=mat(a,a)+0.25d0 
             print *, mat(a,a), btest(a,i) 
          else
             mat(a,a)=mat(a,a)-0.25d0           
             b=ieor(a,2**i+2**j)
             mat(a,b)=mat(a,b)+0.5d0  
             print *, mat(a,a), b, mat(a,b), btest(a,i), btest(a,j)      
          endif
       enddo
    enddo
    !print *, mat
end program test