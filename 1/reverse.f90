 subroutine reverse(word,n)

 implicit none

 integer :: i,n
 character(n) :: word,rword

 rword=''
 do i=1,n
   rword(i:i)=word(n-i+1:n-i+1)
 end do
 word=rword
  
 end subroutine reverse
