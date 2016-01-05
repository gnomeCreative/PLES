!***********************************************************************
subroutine buffer1g(var,n,kest,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! update sbuff vector for sending variables
   !
   use turbo3bis
   !
   use scala3
   !
   use mpi

   implicit none
   !

   !-----------------------------------------------------------------------
   !     array declaration
   integer myid,nproc,kparasta,kparaend
   real var(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   !
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         sbuff1(m)=var(i,j,kest)
      !
      enddo
   enddo
   !
   return
end
