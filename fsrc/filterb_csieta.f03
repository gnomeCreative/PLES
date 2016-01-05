!***********************************************************************
subroutine filterb_csieta(pp0g,pp2g,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter for scale similar model on csi and eta
   ! filtering coefficent 1/8 3/4 1/8
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
   !
   real pp0g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on csi
   real pp2g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp1g(i,j,k)=.125*(pp0g(i-1,j,k)+6.*pp0g(i,j,k)+pp0g(i+1,j,k))
         enddo
      enddo     
   enddo
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp2g(i,j,k)=.125*(pp1g(i,j-1,k)+6.*pp1g(i,j,k)+pp1g(i,j+1,k))
         enddo
      enddo
   enddo
   !
   !     filtering on zita is done in a next step
   !
   return
end
