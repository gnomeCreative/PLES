!***********************************************************************
subroutine filter06(pp2,pp3,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on zita
   !
   use scala3
   use period
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   integer i,j,k,myid,nproc
   integer kparasta,kparaend
   integer kparastal,kparaendl

   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on zita
   !-----------------------------------------------------------------------
   !
   ! imposed only non periodicity because periodicity has been already done
   ! passing the ghost plane
   !
   !     filter on zita
   !
   if(myid.eq.0)then
      kparastal=kparasta+kp
      kparaendl=kparaend
   else if (myid.eq.nproc-1) then
      kparastal=kparasta
      kparaendl=kparaend-kp
   else
      kparastal=kparasta
      kparaendl=kparaend
   endif
   !
   do k=kparastal,kparaendl
      do j=1+jp,jy-jp
         do i=1+ip,jx-ip
            pp3(i,j,k)=.125*(pp2(i,j,k-1)+6.*pp2(i,j,k)+pp2(i,j,k+1))
         enddo
      enddo     
   enddo
   !
   return
end
