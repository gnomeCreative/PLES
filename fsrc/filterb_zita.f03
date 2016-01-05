!***********************************************************************
subroutine filterb_zita(pp2,pp3,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on zita for scale similar model
   ! filtering coefficent 1/8, 6/8, 1/8
   !
   use scala3
   use period
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,myid,nproc
   integer kparasta,kparaend
   integer kparastal,kparaendl

   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on zita
   !-----------------------------------------------------------------------
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp3(i,j,k)=.125*(pp2(i,j,k-1)+6.*pp2(i,j,k)+pp2(i,j,k+1))
         enddo
      enddo     
   enddo
   !
   return
end
