 implicit none

 character(80) :: word1,word2

 print*,'Give two words'; read*,word1,word2
 call reverse(word1)
 call reverse(word2)
 print*,trim(word2),' ',trim(word1)

 contains

   subroutine reverse(word)

   implicit none

   integer :: i,n
   character(80) :: word,rword

   rword=''
   n=len_trim(word)
   do i=1,n
     rword(i:i)=word(n-i+1:n-i+1)
   end do
   word=rword
  
   end subroutine reverse

 end

