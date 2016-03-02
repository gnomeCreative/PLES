!***********************************************************************
subroutine flucrhoesp(r,isc,tipo,body)
   !***********************************************************************
   ! compute explicit diffusive term
   !
   ! k*g11*D(rho)/D(csi) +
   ! k*g22*D(rho)/D(eta) +
   ! k*g33*D(rho)/D(zita)
   !
   ! the subroutine is for scalar eq.
   !
   use mysending
   use myarrays_metri3
   use myarrays_density
   !
   use scala3
   use period
   use tipologia
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr,status(MPI_STATUS_SIZE)
   integer m
   integer kparastal,kparaendl
   integer i,j,k,ii,jj,kk,isc
   real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)

   ! giulia aggiungo questi per eliminare i flussi tra i solidi
   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   integer body
   !-----------------------------------------------------------------------
   ! term nni*g11*d/d(csi)
   !
   do ii=1,ip
      !hicco perche' usa due cicli invece di uno???? ridurre
      ! analogamente per sotto
      !
      ! side left

      do k=kparasta,kparaend
         do j=1,jy
            f1(0,j,k)=f1(0,j,k)+akapt(isc,0,j,k)*g11(0,j,k)* &
               (-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
         end do
      end do
      !
      ! side right
      !
      do k=kparasta,kparaend
         do j=1,jy
            f1(jx,j,k)=f1(jx,j,k) &
               +akapt(isc,jx+1,j,k)*g11(jx,j,k)* &
               (8.*r(jx+1,j,k)-9.*r(jx,j,k)+r(jx-1,j,k))/3.
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
            f1(i,j,k)=f1(i,j,k) &
               +.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)* &
               (r(i+1,j,k)-r(i,j,k))
         !
         enddo
      enddo
   enddo

   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
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
      ! side bottom
      !
      do k=kparasta,kparaend
         do i=1,jx
            f2(i,0,k)=f2(i,0,k)+akaptV(isc,i,0,k)*g22(i,0,k)* &
               (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
         end do
      end do
      !
      ! side upper
      !
      do k=kparasta,kparaend
         do i=1,jx
            f2(i,jy,k)=f2(i,jy,k) &
               +akaptV(isc,i,jy+1,k)*g22(i,jy,k)* &
               (8.*r(i,jy+1,k)-9.*r(i,jy,k)+r(i,jy-1,k))/3.
         end do
      end do
   !
   enddo
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do i=1,jx
         do j=jp,jy-jp
            !
            f2(i,j,k)=f2(i,j,k) &
               +.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)* &
               (r(i,j+1,k)-r(i,j,k))
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
      ! side back
      !
      if (myid.eq.0) then

         do i=1,jx
            do j=1,jy
               f3(i,j,0)=f3(i,j,0)+akapt(isc,i,j,0)*g33(i,j,0)* &
                  (-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
            end do
         end do

      ! side front
      !
      else if (myid.eq.nproc-1) then

         do i=1,jx
            do j=1,jy
               f3(i,j,jz)=f3(i,j,jz) &
                  +akapt(isc,i,j,jz+1)*g33(i,j,jz)* &
                  (8.*r(i,j,jz+1)-9.*r(i,j,jz)+r(i,j,jz-1))/3.
            end do
         end do

      endif
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
   else if ((myid.ne.0).and.(myid.ne.nproc-1)) then
      kparastal=kparasta
      kparaendl=kparaend
   endif


   do j=1,jy
      do i=1,jx
         do k=kparastal,kparaendl
            !
            f3(i,j,k)=f3(i,j,k) &
               +.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))*g33(i,j,k)* &
               (r(i,j,k+1)-r(i,j,k))
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
   ! make border values f1, f2, f3 known to closer procs
   !
   call MPI_SENDRECV(f3(1,1,kparaend),jx*jy, &
      MPI_REAL_SD,rightpe,tagrs, &
      f3(1,1,kparasta-1),jx*jy, &
      MPI_REAL_SD,leftpe,taglr, &
      MPI_COMM_WORLD,status,ierr)
   !
   return
end
