!***********************************************************************
subroutine buffer1d_par(sbuff,var,n,kest,kparasta,kparaend)
   !***********************************************************************
   ! update vector sbuff to exchange variable
   !
   use scala3

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer kparasta,kparaend
   double precision var(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !0:n3+1)
   double precision sbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   do j=0,jy+1
      do i=0,jx+1

         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         sbuff(m)=var(i,j,kest)

      enddo
   enddo
   !
   return
end
