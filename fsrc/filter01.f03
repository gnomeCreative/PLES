!***********************************************************************
subroutine filter01(pp0,pp2,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! apply test filter on csi and eta
   !
   use scala3
   use period
   !
   use mpi

   implicit none
   !
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,myid,nproc
   integer kparasta,kparaend
   real pp0(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered in csi
   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered in eta
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
   !     periodicity in eta (necessary for next filtering)
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
   !     periodicity and filtering on zita will done later in a second step
   !
   return
end
