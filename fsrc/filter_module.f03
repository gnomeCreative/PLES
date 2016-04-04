module filter_module

   use scala3
   use period
   !
   use mpi

   implicit none

contains

!***********************************************************************
subroutine filter01(pp0,pp2,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! apply test filter on csi and eta
   !

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

!***********************************************************************
subroutine filter01np(myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! apply test filter on csi and eta
   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1
   !
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

!***********************************************************************
subroutine filter02np(myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! apply test filter in csi and eta for a product
   !
   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1

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

!***********************************************************************
subroutine filter03(pp2,pp3,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   !     apply test filter on zita
   !
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

!***********************************************************************
subroutine filter04b(pp0,pp2,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   !     apply test filter on csi and eta
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparasta,kparaend,myid,nproc

   real pp0(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on csi
   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filtered on csi
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            pp1(i,j,k)=.125*(pp0(i-1,j,k)+6.*pp0(i,j,k)+pp0(i+1,j,k))
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
            pp2(i,j,k)=.125*(pp1(i,j-1,k)+6.*pp1(i,j,k)+pp1(i,j+1,k))
         enddo
      enddo
   enddo
   !
   !     periodicity and filter in zita are done in a next step
   !
   return
end

!***********************************************************************
subroutine filter04g(pp0g,pp2g,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on csi and eta
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparasta,kparaend
   integer kparastal,kparaendl
   integer myid,nproc

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
      do j=1+jp,jy-jp
         do i=1+ip,jx-ip
            pp2g(i,j,k)=.125*(pp1g(i,j-1,k)+6.*pp1g(i,j,k)+pp1g(i,j+1,k))
         enddo
      enddo
   enddo
   !
   !     periodicty and filtering in zita are done in a next step
   !
   return
end

!***********************************************************************
subroutine filter05g(p0ag,p0bg,pp2g,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on csi and eta for a product
   !
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

!***********************************************************************
subroutine filter06(pp2,pp3,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on zita
   !
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

!***********************************************************************
subroutine filter06b(pp2,pp3,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on zita
   !

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,myid,nproc
   integer kparasta,kparaend

   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   ! matrix to filter
   real pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   ! matrix filtered on zita
   !-----------------------------------------------------------------------
   ! impose only not periodicity because periodicity has been already done
   ! passing ghost plane
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
            pp3(i,j,k)=.125*(pp2(i,j,k-1)+6.*pp2(i,j,k)+pp2(i,j,k+1))
         enddo
      enddo
   enddo
   !
   return
end

!***********************************************************************
subroutine filterb_csieta(pp0g,pp2g,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter for scale similar model on csi and eta
   ! filtering coefficent 1/8 3/4 1/8
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

!***********************************************************************
subroutine filterb_zita(pp2,pp3,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! apply test filter on zita for scale similar model
   ! filtering coefficent 1/8, 6/8, 1/8
   !
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

! BUFFERS ------------------------------------------------------------------

!***********************************************************************
subroutine buffer1d_par(sbuff,var,n,kest,kparasta,kparaend)
   !***********************************************************************
   ! update vector sbuff to exchange variable
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer kparasta,kparaend
   real var(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !0:n3+1)
   real sbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   do j=0,jy+1
      do i=0,jx+1

         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         sbuff(m)=var(i,j,kest)

      enddo
   enddo
   !
   return
end

!***********************************************************************
subroutine buffer1g(var,n,kest,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! update sbuff vector for sending variables
   !
   !use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3,
   use turbo_module, only: rbuff1, sbuff1

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer myid,nproc,kparasta,kparaend
   real var(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   !
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         sbuff1(m)=var(i,j,kest)
      !
      enddo
   enddo
   !
   return
end

!***********************************************************************
subroutine buffer1old_par_nscal(sbuff,var,n,kest,isc)
   !***********************************************************************
   ! update sbuff vector to send variables
   use mysending

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   real var(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real sbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m,isc
   !-----------------------------------------------------------------------
   !
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         sbuff(m)=var(isc,i,j,kest)
      !
      enddo
   enddo
   !
   return
end


!***********************************************************************
subroutine buffer1old_par(sbuff,var,n,kest)
   !***********************************************************************
   ! update sbuff vector to send variables
   use mysending

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   real var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real sbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   !
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         sbuff(m)=var(i,j,kest)
      !
      enddo
   enddo
   !
   return
end

!***********************************************************************
subroutine buffer2(rbuff,var,n,kest,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   !     put data from buff in the variable

   implicit none
   !

   !-----------------------------------------------------------------------
   !     array declaration
   integer myid,nproc,kparasta,kparaend
   integer i,j,n,kest,m

   real var(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   real rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   !
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j,kest)=rbuff(m)
      !
      enddo
   enddo
   !
   return
end

!***********************************************************************
subroutine buffer2d_par(rbuff,var,n,kest)
   !***********************************************************************
   !     put buff vector in the variable
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,n,kest,m
   real var(0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   real rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   do j=0,jy+1
      do i=0,jx+1

         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j,kest)=rbuff(m)

      enddo
   enddo
   !
   return
end

!**********************************************************************
subroutine buffer2g(rbuff,var,n,myid,nproc,kparasta,kparaend)
   !**********************************************************************
   !     put data from buff to the variable
   !
   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1


   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer myid,nproc,kparasta,kparaend
   integer i,j,n,m
   real var(0:n1+1,0:n2+1)
   real rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   !
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j)=rbuff(m)
      !
      end do
   end do
   !
   return
end

!***********************************************************************
subroutine buffer2gg(var,n,myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! put the variables from the reciving buff in the variables

   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer myid,nproc,kparasta,kparaend
   integer i,j,n,m
   real var(0:n1+1,0:n2+1)
   !-----------------------------------------------------------------------

   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j)=rbuff1(m)
      !
      end do
   end do
   !
   return
end

!***********************************************************************
subroutine buffer2old_par_nscal(rbuff,var,n,kest,isc)
   !***********************************************************************
   ! put the reciving buff vector in the variables

   use mysending

   implicit none
   !----------------------------------------------------------------------
   !     array declaration
   real var(nscal,0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   real rbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m,isc
   !-----------------------------------------------------------------------
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(isc,i,j,kest) = rbuff(m)
      !
      enddo
   enddo
   !
   return
end

!***********************************************************************
subroutine buffer2old_par(rbuff,var,n,kest)
   !***********************************************************************
   ! put the reciving buff vector in the variables

   use mysending

   implicit none
   !----------------------------------------------------------------------
   !     array declaration
   real var(0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   real rbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   do j=0,jy+1
      do i=0,jx+1
         !
         m=(n-1)*(jx+2)*(jy+2)+i+1+(jx+2)*j
         var(i,j,kest) = rbuff(m)
      !
      enddo
   enddo
   !
   return
end

!**********************************************************************
subroutine buffvect1d(sbuff,var,n)
   !**********************************************************************
   !
   ! aggiorna il vettore sbuffd per lo scambio di variabili
   implicit none
   !
   real var(n2)
   real sbuff((n1+2)*(n2+2)*40)
   integer j,n,m

   do j=1,jy

      m=(n-1)*jy+j
      sbuff(m)=var(j)

   enddo
   !
   return
end

!**********************************************************************
subroutine buffvect2d(rbuff,var,n)
   !**********************************************************************
   !
   ! aggiorna il vettore rbuffd per lo scambio di variabili

   use scala3

   implicit none
   !
   real var(n2)
   real rbuff((n1+2)*(n2+2)*40)
   integer j,n,m

   do j=1,jy

      m=(n-1)*jy+j
      var(j)=rbuff(m)

   enddo
   !
   return
end



end module filter_module
