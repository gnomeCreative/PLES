!***********************************************************************
subroutine flucnesp(r,visualizzo,tau_wind,tipo,body)
   !***********************************************************************
   ! compute explicit diffusive term
   !
   ! NNI*G11*D(u,v,w)/D(csi)+
   ! NNI*G22*D(u,v,w)/D(eta)+
   ! NNI*G33*D(u,v,w)/D(zita)
   !
   use mysending
   use myarrays_metri3
   use myarrays2
   !
   use scala3
   use period
   use convex
   use tipologia
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr,status(MPI_STATUS_SIZE)
   integer m
   integer kparastal,kparaendl
   double precision bulk_loc,bulkn
   integer i,j,k,ii,jj,kk
   real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   integer visualizzo!,windyes
   integer req1,req2
   integer istatus(MPI_STATUS_SIZE)

   ! Giulia modificavento: !SANTIAGO FIXED
   !      real vel_tau(0:n1+1,kparasta-1:kparaend+1)
   real, optional, intent(in out) :: tau_wind(0:,kparasta-1:)!(0:n1+1,kparasta-1:kparaend+1)

   real rhow
   ! Giulia modificavento:

   ! giulia aggiungo questi per eliminare i flussi tra i solidi
   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   integer body
   !-----------------------------------------------------------------------
   ! term nni*g11*d/d(csi)
   !
   do ii=1,ip
      !
      !     side left
      do k=kparasta,kparaend
         do j=1,jy
            f1(0,j,k)=f1(0,j,k)+annit(0,j,k)*g11(0,j,k)* &
               (-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
         end do
      end do
      !
      !     side right
      do k=kparasta,kparaend
         do j=1,jy
            f1(jx,j,k)=f1(jx,j,k)+annit(jx+1,j,k)*g11(jx,j,k)* &
               (8.*r(jx+1,j,k)-9.*r(jx,j,k)+r(jx-1,j,k))/3.
         end do
      end do
   !
   enddo
   !
   !     into the field
   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
            !
            f1(i,j,k)=f1(i,j,k) &
               +.5*(annit(i,j,k)+annit(i+1,j,k))*g11(i,j,k) &
               *(r(i+1,j,k)-r(i,j,k))
         !
         enddo
      enddo
   enddo


   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
            !

            if(body.eq.1)then
               if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
               if(i.lt.jx)then
                  if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
               endif
            endif
         enddo
      enddo
   enddo

   !
   !-----------------------------------------------------------------------
   ! term nni*g22*d/d(eta)
   !
   do jj=1,jp
      !
      !       side bottom
      do k=kparasta,kparaend
         do i=1,jx

            f2(i,0,k)=f2(i,0,k)+annitV(i,0,k)*g22(i,0,k)* &
               (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
         end do
      end do
      !  c
      !  c     side upper
      if(present(tau_wind))then

         !  c Giulia modificavento:
         rhow = 1000.

         do k=kparasta,kparaend
            do i=1,jx

               !    c      f2(i,jy,k) = abs(vel_tau(i,k))*vel_tau(i,k)*areola4(i,k)
               f2(i,jy,k)=tau_wind(i,k)*areola4(i,k)/rhow

            end do
         end do
      !  c Giulia modificavento:
   
      else

         do k=kparasta,kparaend
            do i=1,jx

               f2(i,jy,k)=f2(i,jy,k) &
                  +annitV(i,jy+1,k)*g22(i,jy,k)* &
                  (8.*r(i,jy+1,k)-9.*r(i,jy,k)+r(i,jy-1,k))/3.
            end do
         end do

      end if
   !
   enddo
   !
   !     into the field
   do k=kparasta,kparaend
      do i=1,jx
         do j=jp,jy-jp
            !
            f2(i,j,k)=f2(i,j,k) &
               +.5*(annitV(i,j,k)+annitV(i,j+1,k)) &
               *g22(i,j,k)*(r(i,j+1,k)-r(i,j,k))
         !
         enddo
      enddo
   enddo

   do k=kparasta,kparaend
      do j=jp,jy-jp
         do i=1,jx
            if(body.eq.1)then
               if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
               if(j.lt.jy)then
                  if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
               endif
            endif
         enddo
      enddo
   enddo
   !
   !-----------------------------------------------------------------------
   ! term nni*g33d/d(zita)
   !
   do kk=1,kp
      !
      !     side back
      if (myid.eq.0) then

         do i=1,jx
            do j=1,jy
               f3(i,j,0)=f3(i,j,0)+annit(i,j,0)*g33(i,j,0)* &
                  (-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
            end do
         end do
  
      endif
      !
      !     side front
      if (myid.eq.nproc-1) then

         do i=1,jx
            do j=1,jy
               f3(i,j,jz)=f3(i,j,jz) &
                  +annit(i,j,jz+1)*g33(i,j,jz)* &
                  (8.*r(i,j,jz+1)-9.*r(i,j,jz)+r(i,j,jz-1))/3.
            end do
         end do
 
      endif
   !
   enddo
   !
   !     into the field
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

   do j=1,jy
      do i=1,jx
         do k=kparastal,kparaendl
            !
            f3(i,j,k)=f3(i,j,k) &
               +.5*(annit(i,j,k)+annit(i,j,k+1))*g33(i,j,k) &
               *(r(i,j,k+1)-r(i,j,k))
         !
         enddo
      enddo
   enddo

   do j=1,jy
      do i=1,jx
         do k=kparastal,kparaendl
            if(body.eq.1)then
               if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
               if(k.lt.jz)then
                  if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
               endif
            endif
         enddo
      enddo
   enddo

   !
   ! pass f3 at the border between procs

   if(myid.eq.0)then
      leftpem=MPI_PROC_NULL
      rightpem=rightpe
   else if(myid.eq.nproc-1)then
      leftpem=leftpe
      rightpem=MPI_PROC_NULL
   else
      leftpem=leftpe
      rightpem=rightpe
   endif


   if(rightpem /= MPI_PROC_NULL) then
      call MPI_SSEND(f3(1,1,kparaend),jx*jy,MPI_REAL_SD, &
         rightpem ,tagrs,MPI_COMM_WORLD,ierr)
   endif
   if(leftpem /= MPI_PROC_NULL) then
      call MPI_RECV(f3(1,1,kparasta-1),jx*jy,MPI_REAL_SD, &
         leftpem  ,taglr,MPI_COMM_WORLD,status,ierr)
   endif

   if(rightpem /= MPI_PROC_NULL) then
   !      call MPI_WAIT(req1,istatus,ierr)
   endif
   if(leftpem /= MPI_PROC_NULL) then
   !      call MPI_WAIT(req2,istatus,ierr)
   endif

   !     integral
   bulk_loc=0.
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            bulk_loc=bulk_loc+f3(i,j,k)-f3(i  ,j  ,k-1)+ &
               f2(i,j,k)-f2(i  ,j-1,k  )+ &
               f1(i,j,k)-f1(i-1,j  ,k  )
         end do
      end do
   end do

   ! make bulk known to all procs

   call MPI_ALLREDUCE(bulk_loc,bulkn,1,MPI_DOUBLE_PRECISION, &
      MPI_SUM, &
      MPI_COMM_WORLD,ierr)

   bulk=bulk+bulkn

   return
end
