!***********************************************************************
subroutine contra(kparasta,kparaend,rightpe,leftpe,nproc,myid)
   !***********************************************************************
   ! compute intermediate controvariant velocity
   ! and boundary condition cs for pressure
   use myarrays_metri3
   use myarrays_velo3
   use mysending, only: MPI_REAL_SD
   !
   use scala3
   use period
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,ii,jj,kk
   real uinter,vinter,winter
   !
   integer ierr,myid,nproc
   integer ncolperproc,kparasta,kparaend,m
   integer kparastal,kparaendl
   integer status(MPI_STATUS_SIZE)
   integer rightpe,leftpe,rightpem,leftpem

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
         do j=1,jy
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
            cs2(j,k)=u(jx+1,j,k)*csx(jx,j,k)+ &
               v(jx+1,j,k)*csy(jx,j,k)+ &
               w(jx+1,j,k)*csz(jx,j,k)-uc(jx,j,k)

            if(potenziale)then
               cs2(j,k)=-uc(jx,j,k)
            end if
            !
            uc(jx,j,k)=u(jx+1,j,k)*csx(jx,j,k)+ &
               v(jx+1,j,k)*csy(jx,j,k)+ &
               w(jx+1,j,k)*csz(jx,j,k)
         !
         end do
      end do
   !
   enddo
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
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
   do jj=1,jp
      !
      do k=kparasta,kparaend
         do i=1,jx
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
            cs4(i,k)=u(i,jy+1,k)*etx(i,jy,k)+ &
               v(i,jy+1,k)*ety(i,jy,k)+ &
               w(i,jy+1,k)*etz(i,jy,k)-vc(i,jy,k)

            if(potenziale)then
               cs4(i,k)=-vc(i,jy,k)
            end if

            vc(i,jy,k)=u(i,jy+1,k)*etx(i,jy,k)+ &
               v(i,jy+1,k)*ety(i,jy,k)+ &
               w(i,jy+1,k)*etz(i,jy,k)
         !
         end do
      end do
   !
   enddo
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do j=jp,jy-jp
         do i=1,jx
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
      do j=1,jy
         do i=1,jx
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
               cs6(i,j)=u(i,j,jz+1)*ztx(i,j,jz)+ &
                  v(i,j,jz+1)*zty(i,j,jz)+ &
                  w(i,j,jz+1)*ztz(i,j,jz)-wc(i,j,jz)

               if(potenziale)then
                  cs6(i,j)=-wc(i,j,jz)
               end if
               !
               wc(i,j,jz)=u(i,j,jz+1)*ztx(i,j,jz)+ &
                  v(i,j,jz+1)*zty(i,j,jz)+ &
                  w(i,j,jz+1)*ztz(i,j,jz)
   
            endif
         !
         end do
      end do
   !
   enddo
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

      do j=1,jy
         do i=1,jx
            !
            uinter=.5*(u(i,j,k)+u(i,j,k+1))
            vinter=.5*(v(i,j,k)+v(i,j,k+1))
            winter=.5*(w(i,j,k)+w(i,j,k+1))
            !
            wc(i,j,k)=uinter*ztx(i,j,k)+ &
               vinter*zty(i,j,k)+ &
               winter*ztz(i,j,k)
         !
         enddo
      enddo

   enddo


   ! subroutine diver needs wc(k-1) to compute rhs
   ! so I need to pass kparaend plane  between proc
   !
   if(myid.eq.0)then
      leftpem=MPI_PROC_NULL
      rightpem=rightpe
   else if(myid.eq.nproc-1)then
      leftpem=leftpe
      rightpem=MPI_PROC_NULL
   else if((myid.ne.0).and.(myid.ne.nproc-1))then
      leftpem=leftpe
      rightpem=rightpe
   endif
   !
   call MPI_SENDRECV(wc(1,1,kparaend),jx*jy, &
      MPI_REAL_SD,rightpem,51+myid, &
      wc(1,1,kparasta-1),jx*jy, &
      MPI_REAL_SD,leftpem,50+myid, &
      MPI_COMM_WORLD,status,ierr)
   !
   return
end
