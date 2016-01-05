!***********************************************************************
subroutine diver(kparasta,kparaend,nproc,myid)
   !***********************************************************************
   ! compute the right hand side for Poisson eq.
   !
   use myarrays_velo3
   !
   use scala3
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr,myid,nproc
   integer ncolperproc,kparasta,kparaend,m
   integer i,j,k
   !-----------------------------------------------------------------------
   !
   if(potenziale .eq. 0)then
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               !
               rhs(i,j,k)=(uc(i,j,k)-uc(i-1,j  ,k  )+ &
                  vc(i,j,k)-vc(i  ,j-1,k  )+ &
                  wc(i,j,k)-wc(i  ,j  ,k-1) )/dt
     
            !
            enddo
         enddo
      enddo
   else
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               !
               rhs(i,j,k)=0
            !
            enddo
         enddo
      enddo
   end if
   !
   return
end
