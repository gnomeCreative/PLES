!***********************************************************************
subroutine buffer2d_par(rbuff,var,n,kest)
   !***********************************************************************
   !     put buff vector in the variable
   !
   use scala3

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,n,kest,m
   double precision var(0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   double precision rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   do j=0,jy+1
      do i=0,jx+1

         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j,kest)=rbuff(m)

      enddo
   enddo
   !
   return
end
