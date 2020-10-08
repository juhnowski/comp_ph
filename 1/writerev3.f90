 implicit none

 character(80) :: word1,word2

 print*,'Give two words'; read*,word1,word2
 call reverse(word1(1:len_trim(word1)),len_trim(word1))
 call reverse(word2(1:len_trim(word2)),len_trim(word2))
 print*,trim(word2),' ',trim(word1)

 end


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



