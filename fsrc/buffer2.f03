!***********************************************************************
subroutine buffer2(rbuff,var,n,kest,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   !     put data from buff in the variable
   !
   use scala3
   !
   use mpi

   implicit none
   !

   !-----------------------------------------------------------------------
   !     array declaration
   integer myid,nproc,kparasta,kparaend
   integer i,j,n,kest,m
      
   real var(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   real rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   !
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j,kest)=rbuff(m)
      !
      enddo
   enddo
   !
   return
end
