!***********************************************************************
subroutine inverse_para2(myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! find the component starting from controvariant
   !
   use turbo_module
   !
   use scala3
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparasta,kparaend,myid,nproc

   real den
   !-----------------------------------------------------------------------
   !
   do k=kparasta,kparaend
      do j=1,jy    
         do i=1,jx
            !
            den=ap11(i,j,k)*ap22(i,j,k)*ap33(i,j,k)+ &
               ap12(i,j,k)*ap23(i,j,k)*ap31(i,j,k)+ &
               ap21(i,j,k)*ap32(i,j,k)*ap13(i,j,k)- &
               ap13(i,j,k)*ap22(i,j,k)*ap31(i,j,k)- &
               ap11(i,j,k)*ap23(i,j,k)*ap32(i,j,k)- &
               ap12(i,j,k)*ap21(i,j,k)*ap33(i,j,k)

            pp1(i,j,k)=(pc1(i,j,k)*ap22(i,j,k)*ap33(i,j,k)+ &
               pc2(i,j,k)*ap32(i,j,k)*ap13(i,j,k)+ &
               pc3(i,j,k)*ap12(i,j,k)*ap23(i,j,k)- &
               ap13(i,j,k)*ap22(i,j,k)*pc3(i,j,k)- &
               pc1(i,j,k)*ap23(i,j,k)*ap32(i,j,k)- &
               ap12(i,j,k)*pc2(i,j,k)*ap33(i,j,k))/den
            !
            pp2(i,j,k)=(ap11(i,j,k)*pc2(i,j,k)*ap33(i,j,k)+ &
               pc1(i,j,k)*ap23(i,j,k)*ap31(i,j,k)+ &
               ap21(i,j,k)*pc3(i,j,k)*ap13(i,j,k)- &
               ap13(i,j,k)*pc2(i,j,k)*ap31(i,j,k)- &
               ap11(i,j,k)*ap23(i,j,k)*pc3(i,j,k)- &
               pc1(i,j,k)*ap21(i,j,k)*ap33(i,j,k))/den
            !
            pp3(i,j,k)=(ap11(i,j,k)*ap22(i,j,k)*pc3(i,j,k)+ &
               ap12(i,j,k)*pc2(i,j,k)*ap31(i,j,k)+ &
               ap21(i,j,k)*ap32(i,j,k)*pc1(i,j,k)- &
               pc1(i,j,k)*ap22(i,j,k)*ap31(i,j,k)- &
               ap11(i,j,k)*pc2(i,j,k)*ap32(i,j,k)- &
               ap12(i,j,k)*ap21(i,j,k)*pc3(i,j,k))/den
         !
         end do
      end do
   end do
   !
   return
end
