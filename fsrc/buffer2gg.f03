!***********************************************************************
subroutine buffer2gg(var,n,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! put the variables from the reciving buff in the variables

   use turbo3bis
   !
   use scala3
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer myid,nproc,kparasta,kparaend
   integer i,j,n,m
   real var(0:n1+1,0:n2+1)
   !-----------------------------------------------------------------------

   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j)=rbuff1(m)
      !
      end do
   end do
   !
   return
end
