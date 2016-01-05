!***********************************************************************
subroutine adams(ktime,rve,bapp,kparasta,kparaend)
   !***********************************************************************
   ! compute explicit term for momentum eq.
   !
   use myarrays_velo3
   use myarrays_metri3
   !
   use scala3, only: jx,jy,n1,n2,dt
   !
   use mpi
      
   implicit none
   !
   !-----------------------------------------------------------------------
   ! array declaration
   integer i,j,k,ktime
   integer kparasta,kparaend
   real rve(n1,n2,kparasta:kparaend)
   real bapp(n1,n2,kparasta:kparaend)
   !
   !-----------------------------------------------------------------------
   ! compute convective term + explicit diffusive
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            !
            rhs(i,j,k)=f1(i,j,k)-f1(i-1,j,k)+ &
               f2(i,j,k)-f2(i,j-1,k)+ &
               f3(i,j,k)-f3(i,j,k-1)+ &
               giac(i,j,k)*bapp(i,j,k)
         !
         enddo
      enddo
   enddo
   !
   ! compute Adam-Bashforth part for momentum eq.
   !
   if (ktime.eq.1) then
      !
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               !
               rve(i,j,k)=rhs(i,j,k)
            !
            end do
         end do
      end do
   !
   end if
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            !
            rhs(i,j,k)=dt*(1.5*rhs(i,j,k)-.5*rve(i,j,k)) &
               /giac(i,j,k)
         !
         enddo
      enddo
   enddo
   !
   ! compute explicit part at step n-1
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            !
            rve(i,j,k)=f1(i,j,k)-f1(i-1,j,k)+ &
               f2(i,j,k)-f2(i,j-1,k)+ &
               f3(i,j,k)-f3(i,j,k-1)+ &
               giac(i,j,k)*bapp(i,j,k)
         !
         enddo
      enddo
   enddo
   !
   return
end
