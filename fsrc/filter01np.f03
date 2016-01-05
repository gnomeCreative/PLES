!***********************************************************************
subroutine filter01np(myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! apply test filter on csi and eta
   use turbo2_data
   use turbo3bis
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
   !      pp0 matrix to filter
   !      pp1 matrix filtered on csi
   !      pp2 matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp1(i,j,k)=.25*(pp0(i-1,j,k)+2.*pp0(i,j,k)+pp0(i+1,j,k))
         enddo
      enddo     
   enddo
   !
   !     periodicity on eta (necessary for next filtering)
   !
   do k=kparasta,kparaend
      do i=1,jx
         pp1(i,0,k)=(1-jp)*pp1(i,jy,k) +  &
            jp*(2.*pp1(i,1,k)-pp1(i,2,k))
         pp1(i,jy+1,k)=(1-jp)*pp1(i,1,k) + &
            jp*(2.*pp1(i,jy,k)-pp1(i,jy-1,k))
      enddo
   enddo
   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp2(i,j,k)=.25*(pp1(i,j-1,k)+2.*pp1(i,j,k)+pp1(i,j+1,k))
         enddo
      enddo     
   enddo
   !
   !     periodicity and filtering on zita is done in a next step
   !
   return
end
