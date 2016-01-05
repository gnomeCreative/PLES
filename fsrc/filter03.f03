!***********************************************************************
subroutine filter03(pp2,pp3,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   !     apply test filter on zita
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

   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !matrix to filter
   real pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !matrix filtered
   !-----------------------------------------------------------------------
   !
   !     impose only not periodicity case, because periodicity has been
   !     already done passing ghost plane
   !
   if((kp.eq.1).and.(myid.eq.0))then
      do j=1,jy
         do i=1,jx
            pp2(i,j,0)=2.*pp2(i,j,1)-pp2(i,j,2)
         enddo
      enddo
   else if((kp.eq.1).and.(myid.eq.nproc-1))then
      do j=1,jy
         do i=1,jx
            pp2(i,j,jz+1)=2.*pp2(i,j,jz)-pp2(i,j,jz-1)
         enddo
      enddo
   endif
   !
   !     filter on zita
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp3(i,j,k)=.25*(pp2(i,j,k-1)+2.*pp2(i,j,k)+pp2(i,j,k+1))
         enddo
      enddo
   enddo
   !
   return
end
