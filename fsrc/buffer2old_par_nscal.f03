!***********************************************************************
subroutine buffer2old_par_nscal(rbuff,var,n,kest,isc)
   !***********************************************************************
   ! put the reciving buff vector in the variables

   use mysending
   use scala3

   implicit none
   !----------------------------------------------------------------------
   !     array declaration
   real var(nscal,0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   real rbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m,isc
   !-----------------------------------------------------------------------
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(isc,i,j,kest) = rbuff(m)
      !
      enddo
   enddo
   !
   return
end
