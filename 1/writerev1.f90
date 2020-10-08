 implicit none

 character(80) :: word

 print*,'Give a word'; read*,word
 call reverse
 print*,word

 contains

   subroutine reverse

   implicit none

   integer :: i,n
   character(80) :: rword

   rword=''
   n=len_trim(word)
   do i=1,n
     rword(i:i)=word(n-i+1:n-i+1)
   end do
   word=rword
  
   end subroutine reverse

 end

