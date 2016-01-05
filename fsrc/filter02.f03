!***********************************************************************
subroutine filter02(p0a,p0b,pp2,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! apply test filter in csi and eta on a product
   !
   use scala3
   use period
   !
   use mpi

   implicit none
   !
   !-----------------------------------------------------------------------
   integer i,j,k
   integer kparasta,kparaend,myid,nproc

   real p0a(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real p0b(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered in csi
   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered in eta
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
