module filter_module

   use scala3
   use mysending
   use period
   !
   use mpi

   implicit none

contains

subroutine filter01(pp0,pp2)
   ! apply test filter on csi and eta
   !

   implicit none
   !
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   real pp0(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered in csi
   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered in eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1(i,j,k)=.25*(pp0(i-1,j,k)+2.*pp0(i,j,k)+pp0(i+1,j,k))
         end do
      end do
   end do
   !
   !     periodicity in eta (necessary for next filtering)
   !
   do k=kparasta,kparaend
      do i=1,n1
         pp1(i,0,k)=2.*pp1(i,1,k)-pp1(i,2,k)
         pp1(i,n2+1,k)=2.*pp1(i,n2,k)-pp1(i,n2-1,k)
      end do
   end do

   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp2(i,j,k)=.25*(pp1(i,j-1,k)+2.*pp1(i,j,k)+pp1(i,j+1,k))
         end do
      end do
   end do
   !
   !     periodicity and filtering on zita will done later in a second step
   !
   return
end

subroutine filter01np()
   ! apply test filter on csi and eta
   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1
   !
   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   !      pp0 matrix to filter
   !      pp1 matrix filtered on csi
   !      pp2 matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1(i,j,k)=.25*(pp0(i-1,j,k)+2.*pp0(i,j,k)+pp0(i+1,j,k))
         end do
      end do
   end do
   !
   !     periodicity on eta (necessary for next filtering)
   !
   do k=kparasta,kparaend
      do i=1,n1
         pp1(i,0,k)=2.*pp1(i,1,k)-pp1(i,2,k)
         pp1(i,n2+1,k)=2.*pp1(i,n2,k)-pp1(i,n2-1,k)
      end do
   end do
   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp2(i,j,k)=.25*(pp1(i,j-1,k)+2.*pp1(i,j,k)+pp1(i,j+1,k))
         end do
      end do
   end do
   !
   !     periodicity and filtering on zita is done in a next step
   !
   return
end

subroutine filter02(p0a,p0b,pp2)
   ! apply test filter in csi and eta on a product
   !
   implicit none
   !
   !-----------------------------------------------------------------------
   integer i,j,k

   real p0a(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real p0b(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered in csi
   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered in eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1(i,j,k)=.25*(p0a(i-1,j,k)*p0b(i-1,j,k)+ 2.*p0a(i,j,k)*p0b(i,j,k)+ p0a(i+1,j,k)*p0b(i+1,j,k))
         end do
      end do
   end do
   !
   !     periodicity on eta (necessary for next filtering)
   !
   do k=kparasta,kparaend
      do i=1,n1
         pp1(i,0,k)=2.*pp1(i,1,k)-pp1(i,2,k)
         pp1(i,n2+1,k)=2.*pp1(i,n2,k)-pp1(i,n2-1,k)
      end do
   end do

   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp2(i,j,k)=.25*(pp1(i,j-1,k)+2.*pp1(i,j,k)+pp1(i,j+1,k))
         end do
      end do
   end do
   !
   !     periodicity and filtering on zita is done in a next step
   !
   return
end

subroutine filter02np()
   ! apply test filter in csi and eta for a product
   !
   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   !     p0a  matrix to filter
   !     p0b  matrix to filter
   !     pp1  product matrix filtered in csi
   !     pp2  product matrix filtered in eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1(i,j,k)=.25*(p0a(i-1,j,k)*p0b(i-1,j,k)+ 2.*p0a(i,j,k)*p0b(i,j,k)+ p0a(i+1,j,k)*p0b(i+1,j,k))
         end do
      end do
   end do
   !
   !     periodicity in eta (necessay for next filtering)
   !
   do k=kparasta,kparaend
      do i=1,n1
         pp1(i,0,k)=2.*pp1(i,1,k)-pp1(i,2,k)
         pp1(i,n2+1,k)=2.*pp1(i,n2,k)-pp1(i,n2-1,k)
      end do
   end do

   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp2(i,j,k)=.25*(pp1(i,j-1,k)+2.*pp1(i,j,k)+pp1(i,j+1,k))
         end do
      end do
   end do
   !
   !     periodicity and filtering in zita is done in a next step
   !
   return
end

subroutine filter03(pp2,pp3)
   !     apply test filter on zita
   !
   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k

   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !matrix to filter
   real pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !matrix filtered
   !-----------------------------------------------------------------------
   !
   !     impose only not periodicity case, because periodicity has been
   !     already done passing ghost plane
   !
   if((kp.eq.1).and.(myid.eq.0))then
      do j=1,n2
         do i=1,n1
            pp2(i,j,0)=2.*pp2(i,j,1)-pp2(i,j,2)
         end do
      end do
   else if((kp.eq.1).and.(myid.eq.nproc-1))then
      do j=1,n2
         do i=1,n1
            pp2(i,j,n3+1)=2.*pp2(i,j,n3)-pp2(i,j,n3-1)
         end do
      end do
   endif
   !
   !     filter on zita
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp3(i,j,k)=.25*(pp2(i,j,k-1)+2.*pp2(i,j,k)+pp2(i,j,k+1))
         end do
      end do
   end do
   !
   return
end

subroutine filter04b(pp0,pp2)
   !     apply test filter on csi and eta
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k

   real pp0(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on csi
   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filtered on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1(i,j,k)=.125*(pp0(i-1,j,k)+6.*pp0(i,j,k)+pp0(i+1,j,k))
         end do
      end do
   end do
   !
   !     periodicity in eta (necessary for next filtering)
   !
   do k=kparasta,kparaend
      do i=1,n1
         pp1(i,0,k)=2.*pp1(i,1,k)-pp1(i,2,k)
         pp1(i,n2+1,k)=2.*pp1(i,n2,k)-pp1(i,n2-1,k)
      end do
   end do
   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp2(i,j,k)=.125*(pp1(i,j-1,k)+6.*pp1(i,j,k)+pp1(i,j+1,k))
         end do
      end do
   end do
   !
   !     periodicity and filter in zita are done in a next step
   !
   return
end

subroutine filter04g(pp0g,pp2g)
   ! apply test filter on csi and eta
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparastal,kparaendl

   real pp0g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on csi
   real pp2g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1g(i,j,k)=.125*(pp0g(i-1,j,k)+6.*pp0g(i,j,k)+pp0g(i+1,j,k))
         end do
      end do
   end do
   !
   do k=kparasta,kparaend
      do j=2,n2-1
         do i=1+ip,n1-ip
            pp2g(i,j,k)=.125*(pp1g(i,j-1,k)+6.*pp1g(i,j,k)+pp1g(i,j+1,k))
         end do
      end do
   end do
   !
   !     periodicty and filtering in zita are done in a next step
   !
   return
end

subroutine filter05g(p0ag,p0bg,pp2g)
   ! apply test filter on csi and eta for a product
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparastal,kparaendl

   real p0ag(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real p0bg(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered on csi
   real pp2g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !product matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1g(i,j,k)=.125*(p0ag(i-1,j,k)*p0bg(i-1,j,k)+ 6.*p0ag(i,j,k)*p0bg(i,j,k)+ p0ag(i+1,j,k)*p0bg(i+1,j,k))
         end do
      end do
   end do
   !
   !     next filter is not computed at the border
   !     so no periodicity applied
   !
   !
   !     filter on eta
   !
   do k=kparasta,kparaend
      do j=2,n2-1
         do i=1+ip,n1-ip
            pp2g(i,j,k)=.125*(pp1g(i,j-1,k)+6.*pp1g(i,j,k)+pp1g(i,j+1,k))
         end do
      end do
   end do
   !
   !     periodicity and filtering on zita are done in a next step
   !
   return
end

subroutine filter06(pp2,pp3)
   ! apply test filter on zita
   !
   implicit none
   !-----------------------------------------------------------------------
   integer i,j,k
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
      do j=2,n2-1
         do i=1+ip,n1-ip
            pp3(i,j,k)=.125*(pp2(i,j,k-1)+6.*pp2(i,j,k)+pp2(i,j,k+1))
         end do
      end do
   end do
   !
   return
end

subroutine filter06b(pp2,pp3)
   ! apply test filter on zita
   !

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k

   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   ! matrix to filter
   real pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   ! matrix filtered on zita
   !-----------------------------------------------------------------------
   ! impose only not periodicity because periodicity has been already done
   ! passing ghost plane
   !
   if((kp.eq.1).and.(myid.eq.0))then
      do j=1,n2
         do i=1,n1
            pp2(i,j,0)=2.*pp2(i,j,1)-pp2(i,j,2)
         end do
      end do
   else if((kp.eq.1).and.(myid.eq.nproc-1))then
      do j=1,n2
         do i=1,n1
            pp2(i,j,n3+1)=2.*pp2(i,j,n3)-pp2(i,j,n3-1)
         end do
      end do
   endif
   !
   !     filter on zita
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp3(i,j,k)=.125*(pp2(i,j,k-1)+6.*pp2(i,j,k)+pp2(i,j,k+1))
         end do
      end do
   end do
   !
   return
end

subroutine filterb_csieta(pp0g,pp2g)
   ! apply test filter for scale similar model on csi and eta
   ! filtering coefficent 1/8 3/4 1/8
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kparastal,kparaendl
   !
   real pp0g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp1g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on csi
   real pp2g(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on eta
   !-----------------------------------------------------------------------
   !
   !     filter on csi
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp1g(i,j,k)=.125*(pp0g(i-1,j,k)+6.*pp0g(i,j,k)+pp0g(i+1,j,k))
         end do
      end do
   end do
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp2g(i,j,k)=.125*(pp1g(i,j-1,k)+6.*pp1g(i,j,k)+pp1g(i,j+1,k))
         end do
      end do
   end do
   !
   !     filtering on zita is done in a next step
   !
   return
end

subroutine filterb_zita(pp2,pp3)
   ! apply test filter on zita for scale similar model
   ! filtering coefficent 1/8, 6/8, 1/8
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,myid,nproc
   integer kparastal,kparaendl

   real pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix to filter
   real pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1)   !matrix filtered on zita
   !-----------------------------------------------------------------------
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=1,n1
            pp3(i,j,k)=.125*(pp2(i,j,k-1)+6.*pp2(i,j,k)+pp2(i,j,k+1))
         end do
      end do
   end do
   !
   return
end

! BUFFERS ------------------------------------------------------------------

subroutine buffer1d_par(sbuff,var,n,kest)
   ! update vector sbuff to exchange variable
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   real var(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !0:n3+1)
   real sbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   do j=0,n2+1
      do i=0,n1+1

         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         sbuff(m)=var(i,j,kest)

      end do
   end do
   !
   return
end

subroutine buffer1g(var,n,kest)
   ! update sbuff vector for sending variables
   !
   !use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3,
   use turbo_module, only: rbuff1, sbuff1

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   real,intent(in) :: var(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   integer,intent(in) :: n,kest
   integer i,j,m
   !-----------------------------------------------------------------------
   !
   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         sbuff1(m)=var(i,j,kest)
      !
      end do
   end do
   !
   return
end

subroutine buffer1old_par_nscal(sbuff,var,n,kest,isc)
   ! update sbuff vector to send variables

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   real var(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real sbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m,isc
   !-----------------------------------------------------------------------
   !
   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         sbuff(m)=var(isc,i,j,kest)
      !
      end do
   end do
   !
   return
end

subroutine buffer1old_par(sbuff,var,n,kest)
   ! update sbuff vector to send variables
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   real var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real sbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   !
   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         sbuff(m)=var(i,j,kest)
      !
      end do
   end do
   !
   return
end

subroutine buffer2(rbuff,var,n,kest)
   !     put data from buff in the variable

   implicit none
   !

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,n,kest,m

   real var(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   real rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   !
   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         var(i,j,kest)=rbuff(m)
      !
      end do
   end do
   !
   return
end

subroutine buffer2d_par(rbuff,var,n,kest)
   !     put buff vector in the variable
   !
   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,n,kest,m
   real var(0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   real rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   do j=0,n2+1
      do i=0,n1+1

         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         var(i,j,kest)=rbuff(m)

      end do
   end do
   !
   return
end

subroutine buffer2g(rbuff,var,n)
   !     put data from buff to the variable
   !
   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1


   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,n,m
   real var(0:n1+1,0:n2+1)
   real rbuff((n1+2)*(n2+2)*40)
   !-----------------------------------------------------------------------
   !
   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         var(i,j)=rbuff(m)
      !
      end do
   end do
   !
   return
end

subroutine buffer2gg(var,n)
   ! put the variables from the reciving buff in the variables

   use turbo_module, only: p0a, p0b, pp0, pp1, pp2, pp3, rbuff1, sbuff1

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,n,m
   real var(0:n1+1,0:n2+1)
   !-----------------------------------------------------------------------

   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         var(i,j)=rbuff1(m)
      !
      end do
   end do
   !
   return
end

subroutine buffer2old_par_nscal(rbuff,var,n,kest,isc)
   ! put the reciving buff vector in the variables

   implicit none
   !----------------------------------------------------------------------
   !     array declaration
   real var(nscal,0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   real rbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m,isc
   !-----------------------------------------------------------------------
   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         var(isc,i,j,kest) = rbuff(m)
      !
      end do
   end do
   !
   return
end

subroutine buffer2old_par(rbuff,var,n,kest)
   ! put the reciving buff vector in the variables
   implicit none
   !----------------------------------------------------------------------
   !     array declaration
   real var(0:n1+1,0:n2+1,kest:kest) !0:n3+1)
   real rbuff((n1+2)*(n2+2)*40)
   integer i,j,n,kest,m
   !-----------------------------------------------------------------------
   do j=0,n2+1
      do i=0,n1+1
         !
         m=(n-1)*(n1+2)*(n2+2)+i+1+(n1+2)*j
         var(i,j,kest) = rbuff(m)
      !
      end do
   end do
   !
   return
end

subroutine buffvect1d(sbuff,var,n)
   !
   ! aggiorna il vettore sbuffd per lo scambio di variabili
   implicit none
   !
   real var(n2)
   real sbuff((n1+2)*(n2+2)*40)
   integer j,n,m

   do j=1,n2

      m=(n-1)*n2+j
      sbuff(m)=var(j)

   end do
   !
   return
end

subroutine buffvect2d(rbuff,var,n)
   !
   ! aggiorna il vettore rbuffd per lo scambio di variabili

   implicit none
   !
   real var(n2)
   real rbuff((n1+2)*(n2+2)*40)
   integer j,n,m

   do j=1,n2

      m=(n-1)*n2+j
      var(j)=rbuff(m)

   end do
   !
   return
end



end module filter_module
