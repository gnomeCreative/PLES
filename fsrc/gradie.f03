!***********************************************************************
subroutine gradie()
   !***********************************************************************
   ! compute pressure gradient at cell centroid *dt/J(-1) and computation
   ! of controvariant gradient on faces
   !
   use myarrays_metri3
   use myarrays_velo3, only: cgra1,cgra2,cgra3,fi,gra1,gra2,gra3
   use mysending, only: myid
   !
   use scala3, only: dt,n1,n2,kparasta,kparaend
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparastal,kparaendl
      
   real coef,coef_dt
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   !
   coef_dt = -0.5 * dt

   !     gradient for cartesian velocity
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
      
            coef = coef_dt/giac(i,j,k)
            !
            gra1(i,j,k)=csx(i,j,k)*(fi(i,j,k)+fi(i+1,j,k))- &
               csx(i-1,j,k)*(fi(i,j,k)+fi(i-1,j,k))+ &
               etx(i,j,k)*(fi(i,j,k)+fi(i,j+1,k))- &
               etx(i,j-1,k)*(fi(i,j,k)+fi(i,j-1,k))+ &
               ztx(i,j,k)*(fi(i,j,k)+fi(i,j,k+1))- &
               ztx(i,j,k-1)*(fi(i,j,k)+fi(i,j,k-1))
            !
            gra1(i,j,k)=coef*gra1(i,j,k)
            !
            gra2(i,j,k)=csy(i,j,k)*(fi(i,j,k)+fi(i+1,j,k))- &
               csy(i-1,j,k)*(fi(i,j,k)+fi(i-1,j,k))+ &
               ety(i,j,k)*(fi(i,j,k)+fi(i,j+1,k))- &
               ety(i,j-1,k)*(fi(i,j,k)+fi(i,j-1,k))+ &
               zty(i,j,k)*(fi(i,j,k)+fi(i,j,k+1))- &
               zty(i,j,k-1)*(fi(i,j,k)+fi(i,j,k-1))
            !
            gra2(i,j,k)=coef*gra2(i,j,k)
            !
            gra3(i,j,k)=csz(i,j,k)*(fi(i,j,k)+fi(i+1,j,k))- &
               csz(i-1,j,k)*(fi(i,j,k)+fi(i-1,j,k))+ &
               etz(i,j,k)*(fi(i,j,k)+fi(i,j+1,k))- &
               etz(i,j-1,k)*(fi(i,j,k)+fi(i,j-1,k))+ &
               ztz(i,j,k)*(fi(i,j,k)+fi(i,j,k+1))- &
               ztz(i,j,k-1)*(fi(i,j,k)+fi(i,j,k-1))
            !
            gra3(i,j,k)=coef*gra3(i,j,k)
         !
         end do
      end do
   end do
   !
   !     gradient for controvariant velocity
   do k=kparasta,kparaend
      do j=1,n2
         do i=0,n1
            !
            cgra1(i,j,k)=g11(i,j,k)*(fi(i+1,j,k)-fi(i,j,k))+ &
               .25*g12(i,j,k)*(fi(i+1,j+1,k)+fi(i,j+1,k) &
               -fi(i+1,j-1,k)-fi(i,j-1,k))+ &
               .25*g13(i,j,k)*(fi(i+1,j,k+1)+fi(i,j,k+1) &
               -fi(i+1,j,k-1)-fi(i,j,k-1))
         !
         end do
      end do
   end do
   !
   do k=kparasta,kparaend
      do j=0,n2
         do i=1,n1
            !
            cgra2(i,j,k)=g22(i,j,k)*(fi(i,j+1,k)-fi(i,j,k))+ &
               .25*g21(i,j,k)*(fi(i+1,j+1,k)+fi(i+1,j,k) &
               -fi(i-1,j+1,k)-fi(i-1,j,k))+ &
               .25*g23(i,j,k)*(fi(i,j+1,k+1)+fi(i,j,k+1) &
               -fi(i,j+1,k-1)-fi(i,j,k-1))
         !
         end do
      end do
   end do
   !
   if (myid==0) then
      kparastal=0
      kparaendl=kparaend
   else
      kparastal=kparasta
      kparaendl=kparaend
   endif

   do k=kparastal,kparaendl
      do j=1,n2
         do i=1,n1
            !
            cgra3(i,j,k)=g33(i,j,k)*(fi(i,j,k+1)-fi(i,j,k))+ &
               .25*g31(i,j,k)*(fi(i+1,j,k+1)+fi(i+1,j,k) &
               -fi(i-1,j,k+1)-fi(i-1,j,k))+ &
               .25*g32(i,j,k)*(fi(i,j+1,k+1)+fi(i,j+1,k) &
               -fi(i,j-1,k+1)-fi(i,j-1,k))
         !
         end do
      end do
   end do
   !
   return
end subroutine gradie
