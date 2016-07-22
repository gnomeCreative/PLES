subroutine contra()
   ! compute intermediate controvariant velocity
   ! and boundary condition cs for pressure
   use myarrays_metri3
   use myarrays_velo3
   !
   use scala3
   use mysending
   use period
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,ii,jj,kk
   real uinter,vinter,winter
   !
   integer kparastal,kparaendl
   integer ierr,status(MPI_STATUS_SIZE)

   real,allocatable :: cs3col(:),cs3tot(:)
   real,allocatable :: cs4col(:),cs4tot(:)
   !
   !-----------------------------------------------------------------------
   ! compute UC and cs1, cs2 at sides left right
   !
   ! sides left and right uc* - uc^(n+1)
   do ii=1,ip
      !
      do k=kparasta,kparaend
         do j=1,n2
            !
            cs1(j,k)=u(0,j,k)*csx(0,j,k)+ &
               v(0,j,k)*csy(0,j,k)+ &
               w(0,j,k)*csz(0,j,k)-uc(0,j,k)


            if(potenziale)then
               cs1(j,k)=-uc(0,j,k)
            end if
            !
            uc(0,j,k)=u(0,j,k)*csx(0,j,k)+ &
               v(0,j,k)*csy(0,j,k)+ &
               w(0,j,k)*csz(0,j,k)

            !
            cs2(j,k)=u(n1+1,j,k)*csx(n1,j,k)+ &
               v(n1+1,j,k)*csy(n1,j,k)+ &
               w(n1+1,j,k)*csz(n1,j,k)-uc(n1,j,k)

            if(potenziale)then
               cs2(j,k)=-uc(n1,j,k)
            end if
            !
            uc(n1,j,k)=u(n1+1,j,k)*csx(n1,j,k)+ &
               v(n1+1,j,k)*csy(n1,j,k)+ &
               w(n1+1,j,k)*csz(n1,j,k)
         !
         end do
      end do
   !
   end do
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do j=1,n2
         do i=ip,n1-ip
            !
            uinter=.5*(u(i,j,k)+u(i+1,j,k))
            vinter=.5*(v(i,j,k)+v(i+1,j,k))
            winter=.5*(w(i,j,k)+w(i+1,j,k))

            uc(i,j,k)=uinter*csx(i,j,k)+ &
               vinter*csy(i,j,k)+ &
               winter*csz(i,j,k)
         !
         end do
      end do
   end do
   !
   !
   !-----------------------------------------------------------------------
   ! compute VC and cs3, cs4 at sides bottom and upper
   !
   ! sides bottom and upper vc* - vc^(n+1)
   !
   ! direction j is always periodic
      !
      do k=kparasta,kparaend
         do i=1,n1
            !
            cs3(i,k)=u(i,0,k)*etx(i,0,k)+ &
               v(i,0,k)*ety(i,0,k)+ &
               w(i,0,k)*etz(i,0,k)-vc(i,0,k)

            if(potenziale)then
               cs3(i,k)=-vc(i,0,k)
            end if

            !
            vc(i,0,k)=u(i,0,k)*etx(i,0,k)+ &
               v(i,0,k)*ety(i,0,k)+ &
               w(i,0,k)*etz(i,0,k)
            !
            cs4(i,k)=u(i,n2+1,k)*etx(i,n2,k)+ &
               v(i,n2+1,k)*ety(i,n2,k)+ &
               w(i,n2+1,k)*etz(i,n2,k)-vc(i,n2,k)

            if(potenziale)then
               cs4(i,k)=-vc(i,n2,k)
            end if

            vc(i,n2,k)=u(i,n2+1,k)*etx(i,n2,k)+ &
               v(i,n2+1,k)*ety(i,n2,k)+ &
               w(i,n2+1,k)*etz(i,n2,k)
         !
         end do
      end do
   !
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do j=1,n2-1
         do i=1,n1
            !
            uinter=.5*(u(i,j,k)+u(i,j+1,k))
            vinter=.5*(v(i,j,k)+v(i,j+1,k))
            winter=.5*(w(i,j,k)+w(i,j+1,k))
            !
            vc(i,j,k)=uinter*etx(i,j,k)+ &
               vinter*ety(i,j,k)+ &
               winter*etz(i,j,k)
         !
         end do
      end do
   end do

   !-----------------------------------------------------------------------
   ! compute WC and cs5, cs6 at sides front and back
   !
   ! sides front and back wc* - wc^(n+1)
   !
   do kk=1,kp
      !
      do j=1,n2
         do i=1,n1
            !
            if(myid.eq.0)then
               cs5(i,j)=u(i,j,0)*ztx(i,j,0)+ &
                  v(i,j,0)*zty(i,j,0)+ &
                  w(i,j,0)*ztz(i,j,0)-wc(i,j,0)
               if(potenziale)then
                  cs5(i,j)=-wc(i,j,0)
               end if

               !
               wc(i,j,0)=u(i,j,0)*ztx(i,j,0)+ &
                  v(i,j,0)*zty(i,j,0)+ &
                  w(i,j,0)*ztz(i,j,0)
            !
            else if(myid.eq.nproc-1)then
               cs6(i,j)=u(i,j,n3+1)*ztx(i,j,n3)+ &
                  v(i,j,n3+1)*zty(i,j,n3)+ &
                  w(i,j,n3+1)*ztz(i,j,n3)-wc(i,j,n3)

               if(potenziale)then
                  cs6(i,j)=-wc(i,j,n3)
               end if
               !
               wc(i,j,n3)=u(i,j,n3+1)*ztx(i,j,n3)+ &
                  v(i,j,n3+1)*zty(i,j,n3)+ &
                  w(i,j,n3+1)*ztz(i,j,n3)
   
            endif
         !
         end do
      end do
   !
   end do
   !
   ! into the field
   !
   if (myid.eq.0) then
      kparastal=kp
      kparaendl=kparaend
   else if (myid.eq.nproc-1) then
      kparastal=kparasta
      kparaendl=kparaend-kp
   else
      kparastal=kparasta
      kparaendl=kparaend
   endif

   do k=kparastal,kparaendl

      do j=1,n2
         do i=1,n1
            !
            uinter=.5*(u(i,j,k)+u(i,j,k+1))
            vinter=.5*(v(i,j,k)+v(i,j,k+1))
            winter=.5*(w(i,j,k)+w(i,j,k+1))
            !
            wc(i,j,k)=uinter*ztx(i,j,k)+ &
               vinter*zty(i,j,k)+ &
               winter*ztz(i,j,k)
         !
         end do
      end do

   end do


   ! subroutine diver needs wc(k-1) to compute rhs so I need to pass kparaend plane between proc
   call MPI_SENDRECV(wc(1,1,kparaend),n1*n2,MPI_REAL_SD,rightpem,51+myid,wc(1,1,kparasta-1),n1*n2, &
      MPI_REAL_SD,leftpem,50+myid,MPI_COMM_WORLD,status,ierr)

   return
end
