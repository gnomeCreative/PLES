subroutine adams(ktime,rve,bapp)
   ! compute explicit term for momentum eq.
   !
   use myarrays_velo3, only: rhs
   use myarrays_metri3, only: f1,f2,f3,giac
   !
   use scala3, only: n1,n2,dt,kparasta,kparaend
      
   implicit none
   !
   !-----------------------------------------------------------------------
   ! array declaration
   integer :: i,j,k,ktime
   real :: rve(n1,n2,kparasta:kparaend)
   real :: bapp(n1,n2,kparasta:kparaend)

   !-----------------------------------------------------------------------
   ! compute convective term + explicit diffusive
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1

            rhs(i,j,k)=f1(i,j,k)-f1(i-1,j,k)+ &
               f2(i,j,k)-f2(i,j-1,k)+ &
               f3(i,j,k)-f3(i,j,k-1)+ &
               giac(i,j,k)*bapp(i,j,k)

         end do
      end do
   end do

   ! compute Adam-Bashforth part for momentum eq.
   if (ktime==1) then

      do k=kparasta,kparaend
         do j=1,n2
            do i=1,n1
               rve(i,j,k)=rhs(i,j,k)
            end do
         end do
      end do

   end if

   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1

            rhs(i,j,k)=dt*(1.5*rhs(i,j,k)-.5*rve(i,j,k)) / giac(i,j,k)

         end do
      end do
   end do

   ! compute explicit part at step n-1
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1

            rve(i,j,k)=f1(i,j,k)-f1(i-1,j,k)+ &
               f2(i,j,k)-f2(i,j-1,k)+ &
               f3(i,j,k)-f3(i,j,k-1)+ &
               giac(i,j,k)*bapp(i,j,k)

         end do
      end do
   end do

   return
end
