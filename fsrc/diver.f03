subroutine diver()
   ! compute the right hand side for Poisson eq.
   !
   use myarrays_velo3
   !
   use scala3
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   !-----------------------------------------------------------------------
   !
   if(.not.potenziale)then
      do k=kparasta,kparaend
         do j=1,n2
            do i=1,n1
               !
               rhs(i,j,k)=(uc(i,j,k)-uc(i-1,j  ,k  )+ &
                  vc(i,j,k)-vc(i  ,j-1,k  )+ &
                  wc(i,j,k)-wc(i  ,j  ,k-1) )/dt
     
            !
            end do
         end do
      end do
   else
      do k=kparasta,kparaend
         do j=1,n2
            do i=1,n1
               !
               rhs(i,j,k)=0
            !
            end do
         end do
      end do
   end if
   !
   return
end
