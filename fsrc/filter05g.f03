!***********************************************************************
subroutine filter05g(p0ag,p0bg,pp2g,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on csi and eta for a product
   !
   use scala3
   use period
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparasta,kparaend
   integer kparastal,kparaendl
   integer myid,nproc
   real p0ag(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real p0bg(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered on csi
   real pp2g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp1g(i,j,k)=.125*(p0ag(i-1,j,k)*p0bg(i-1,j,k)+ &
               6.*p0ag(i,j,k)*p0bg(i,j,k)+ &
               p0ag(i+1,j,k)*p0bg(i+1,j,k))
         enddo
      enddo     
   enddo
   !
   !     next filter is not computed at the border
   !     so no periodicity applied
   !
   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=1+jp,jy-jp
         do i=1+ip,jx-ip
            pp2g(i,j,k)=.125*(pp1g(i,j-1,k)+6.*pp1g(i,j,k)+pp1g(i,j+1,k))
         enddo
      enddo     
   enddo
   !
   !     periodicity and filtering on zita are done in a next step
   !
   return
end
