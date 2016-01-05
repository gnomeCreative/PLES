!***********************************************************************
subroutine contourp(kparasta,kparaend,nproc,myid)
   !***********************************************************************
   ! compute cartesian velocity and controvariant in periodic cell at
   ! step n+1, at the corner computation at the end of the sub
   !
   use mysettings
   use myarrays_velo3
   !
   use scala3
   use period
   use tipologia
   use orl
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,ii,isc,kk
   !
   integer ierr,myid,nproc,status(MPI_STATUS_SIZE)
   integer kparasta,kparaend
   integer kparastal,kparaendl
   integer plantype!,lett
        
   real, allocatable:: rho(:,:,:)
   !-----------------------------------------------------------------------
   ! periodic cell in csi (also outside the boundary)
   !
   do i=1,1-ip

      do k=kparasta,kparaend
         do j=0,jy+1
            !
            u(0   ,j,k)=u(jx,j,k)
            v(0   ,j,k)=v(jx,j,k)
            w(0   ,j,k)=w(jx,j,k)
            u(jx+1,j,k)=u(1 ,j,k)
            v(jx+1,j,k)=v(1 ,j,k)
            w(jx+1,j,k)=w(1 ,j,k)
            do isc=1,nscal
               rhov(isc,0   ,j,k)=rhov(isc,jx,j,k)
               rhov(isc,jx+1,j,k)=rhov(isc,1 ,j,k)
            end do
         !
         enddo
      enddo

   enddo
   !
   ! periodic cell in eta (also outside the boundary)
   !
   do j=1,1-jp

      if (myid.eq.0) then
         kparastal=kp
         kparaendl=kparaend
      else if (myid.eq.nproc-1) then
         kparastal=kparasta
         kparaendl=kparaend+1-kp
      else
         kparastal=kparasta
         kparaendl=kparaend
      endif

      do i=ip,jx+1-ip
         do k=kparastal,kparaendl
            u(i,0   ,k)=u(i,jy,k)
            v(i,0   ,k)=v(i,jy,k)
            w(i,0   ,k)=w(i,jy,k)
            u(i,jy+1,k)=u(i,1 ,k)
            v(i,jy+1,k)=v(i,1 ,k)
            w(i,jy+1,k)=w(i,1 ,k)
            do isc=1,nscal
               rhov(isc,i,0   ,k)=rhov(isc,i,jy,k)
               rhov(isc,i,jy+1,k)=rhov(isc,i,1 ,k)
            end do
         enddo
      enddo

   enddo
   !
   ! periodic cell in zita (also outside the boundary)
   !
   do k=1,1-kp

      call MPI_TYPE_VECTOR(jy+2,jx,jx+2,MPI_REAL_SD,plantype,ierr)
      call MPI_TYPE_COMMIT(plantype,ierr)

      if (myid.eq.0) then
         call MPI_SENDRECV(u(1,0,1),1,plantype,nproc-1,12,u(1,0,0),1,plantype,nproc-1,11,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(v(1,0,1),1,plantype,nproc-1,14,v(1,0,0),1,plantype,nproc-1,13,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(w(1,0,1),1,plantype,nproc-1,16,w(1,0,0),1,plantype,nproc-1,15,MPI_COMM_WORLD,status,ierr)
      endif

      if (myid.eq.nproc-1) then

         call MPI_SENDRECV(u(1,0,jz),1,plantype,0,11,u(1,0,jz+1),1,plantype,0,12,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(v(1,0,jz),1,plantype,0,13,v(1,0,jz+1),1,plantype,0,14,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(w(1,0,jz),1,plantype,0,15,w(1,0,jz+1),1,plantype,0,16,MPI_COMM_WORLD,status,ierr)

      endif

      allocate(rho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))
     
      do isc=1,nscal
         do kk=kparasta-1,kparaend+1 !0,jz+1
            do j=0,jy+1
               do i=0,jx+1
                  rho(i,j,kk)=rhov(isc,i,j,kk)
               end do
            end do
         end do

         if (myid.eq.0) then
            call MPI_SENDRECV(rho(1,0,1),1,plantype,nproc-1,18+isc,rho(1,0,0),&
               1,plantype,nproc-1,17+isc,MPI_COMM_WORLD,status,ierr)
         end if
         if(myid.eq.nproc-1)then
            call MPI_SENDRECV(rho(1,0,jz),1,plantype,0,17+isc,rho(1,0,jz+1),1,plantype,0,18+isc,MPI_COMM_WORLD,status,ierr)
         end if

         do kk=kparasta-1,kparaend+1 !0,jz+1
            do j=0,jy+1
               do i=0,jx+1
                  rhov(isc,i,j,kk)=rho(i,j,kk)
               end do
            end do
         end do
      end do !nscal
      deallocate(rho)
      call MPI_TYPE_FREE(plantype,ierr)

   enddo
   !
   ! To calculate the next pressure at surface for the next iteration
   if(myid.eq.0)then
      kparastal=0
      kparaendl=kparaend
   else if (myid.eq.nproc-1) then
      kparastal=kparasta
      kparaendl=kparaend+1
   else
      kparastal=kparasta
      kparaendl=kparaend
   endif

   if(freesurface == 1)then !free surface ON.<<<<<<<<<
      if(myid.eq.0)then
         write(*,*)'Free Surface ON'
      end if
      do k=kparastal,kparaendl
         do i=0,jx+1
            next_prs(i,k)=((fi(i,jy+1,k)+fi(i,jy,k))*0.5) - (v(i,jy+1,k) * dt * bby)
         enddo
      enddo
   else
      if(myid.eq.0)then
         write(*,*)'Free Surface OFF'
      end if
   end if !free surface<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   return
end
