!***********************************************************************
subroutine filter02np(myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! apply test filter in csi and eta for a product
   !
   use turbo3bis
   use turbo2_data

   use scala3
   use period
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparasta,kparaend,myid,nproc
   !     p0a  matrix to filter
   !     p0b  matrix to filter
   !     pp1  product matrix filtered in csi
   !     pp2  product matrix filtered in eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp1(i,j,k)=.25*(p0a(i-1,j,k)*p0b(i-1,j,k)+ &
               2.*p0a(i,j,k)*p0b(i,j,k)+ &
               p0a(i+1,j,k)*p0b(i+1,j,k))
         enddo
      enddo     
   enddo
   !
   !     periodicity in eta (necessay for next filtering)
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
   !     periodicity and filtering in zita is done in a next step
   !
   return
end
