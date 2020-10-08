 module powermodule

   integer, parameter :: nd=10, np=3
   real(8) :: xpowers(np,nd)

 contains

   subroutine initpowers(dx)

   implicit none

   integer :: i,j
   real(8) :: dx
 
   do i=1,nd
      xpowers(1,i)=dx*i
   enddo
   do i=2,np
      xpowers(i,:)=xpowers(1,:)**i
   enddo
    
   end subroutine initpowers

 end module powermodule
