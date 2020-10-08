!-------------!
 module random
!-------------!
 save

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64,oldseed
 
 contains

  !----------------------!
   real(8) function ran()
  !----------------------------------------------!
  ! 64-bit congruental generator                 !
  ! iran64=oran64*2862933555777941757+1013904243 !
  !----------------------------------------------!

   ran64=ran64*mul64+add64
   ran=0.5d0+dmu64*dble(ran64)

   end function ran
  !----------------!

  !---------------------!
   subroutine initran(w)
  !---------------------!
   integer(8) :: c1,c2,c3,irmax
   integer(4) :: w,nb,b
      
   irmax=2_8**31
   irmax=2*(irmax**2-1)+1
   mul64=2862933555777941757_8
   add64=1013904243_8
   dmu64=0.5d0/dble(irmax)

   if (w>=0) then
      open(10,file='seed.in',status='old')
      read(10,*)ran64
      close(10)
      oldseed=ran64
   else
      ran64=oldseed
   endif
   if (w>0) then
      open(10,file='seed.in',status='unknown')
      write(10,*)abs((ran64*mul64)/5+5265361)
      close(10)
   endif

   end subroutine initran
  !----------------------!

 end module random
!-----------------!
