 implicit none

 character(80) :: word1,word2

 print*,'Give two words'; read*,word1,word2
 call reverse(word1,len_trim(word1))
 call reverse(word2,len_trim(word2))
 print*,trim(word2),' ',trim(word1)

 end
