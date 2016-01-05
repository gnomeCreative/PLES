!***********************************************************************
subroutine filtro(r,xstart,xend,ystart,yend,zstart,zend)
   !***********************************************************************
   ! applica il test filter su csi ed eta
   !
   use mysending
   !
   use scala3
   use period
   use tipologia
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
          
   integer istatus,ierr
   integer req1,req2,req3,req4

   real  r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   real rf(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
     
   real r_destra,r_sinistra
   real r_sotto,r_sopra
   real r_indietro,r_avanti

   integer xstart,xend,ystart,yend,zstart,zend

   integer status(MPI_STATUS_SIZE)
   integer plantype
   !-----------------------------------------------------------------------
   !     initialize
   !
   do k=kparasta-1,kparaend+1
      do j=0,jy+1
         do i=0,jx+1
            rf(i,j,k)=0.
         enddo
      enddo
   enddo

   !     direction k

   ! if filter is applied first to k it is not necessary an intermediate sending
   ! for planes
   do k=zstart,zend !kparasta,kparaend
      do j=ystart,yend 
         do i=xstart,xend
      
            r_indietro=0.5*(r(i,j,k)+r(i,j,k-1))
            r_avanti=  0.5*(r(i,j,k)+r(i,j,k+1))

            if (k.eq.1) then
               r_indietro=r(i,j,0)
            end if

            if (k.eq.n3) then
               r_avanti=r(i,j,n3+1)
            end if
                
            rf(i,j,k)=(r_indietro+2.*r(i,j,k)+r_avanti)*0.25
    
         enddo
      enddo
   enddo

   do k=zstart,zend !kparasta,kparaend
      do j=ystart,yend 
         do i=xstart,xend
            r(i,j,k) = rf(i,j,k)
         enddo
      enddo
   enddo

   !     direction i
   do k=zstart,zend !kparasta,kparaend
      do j=ystart,yend 
         do i=xstart,xend
      
            r_sinistra=0.5*(r(i,j,k)+r(i-1,j,k))
            r_destra = 0.5*(r(i,j,k)+r(i+1,j,k))

            if (i.eq.1) then
               r_sinistra=r(0,j,k)
            end if

            if (i.eq.n1) then
               r_destra=r(n1+1,j,k)
            end if

            rf(i,j,k)=(r_sinistra+2.*r(i,j,k)+r_destra)*0.25
          
         enddo
      enddo
   enddo

   do k=zstart,zend !kparasta,kparaend
      do j=ystart,yend 
         do i=xstart,xend
            r(i,j,k) = rf(i,j,k)
         enddo
      enddo
   enddo
      
   !     direction j
   do k=zstart,zend !kparasta,kparaend
      do j=ystart,yend 
         do i=xstart,xend
      
      
            r_sotto= 0.5*(r(i,j,k)+r(i,j-1,k))
            r_sopra= 0.5*(r(i,j,k)+r(i,j+1,k))

            if (j.eq.1) then
               r_sotto=r(i,0,k)
            end if

            if (j.eq.n2) then
               r_sopra=r(i,n2+1,k)
            end if
      
            rf(i,j,k)=(r_sotto+2.*r(i,j,k)+r_sopra)*0.25

         enddo
      enddo
   enddo

   do k=zstart,zend !kparasta,kparaend
      do j=ystart,yend 
         do i=xstart,xend
            r(i,j,k) = rf(i,j,k)
         enddo
      enddo
   enddo
         
   !-----------------------------------------------------------------------
   !    send planes at the border
   !-----------------------------------------------------------------------

   if(myid.ne.0) then
      call MPI_SSEND(r(0,0,kparasta),(jx+2)*(jy+2), &
         MPI_REAL_SD,leftpem,tagls, &
         MPI_COMM_WORLD,ierr)
   !         call MPI_WAIT(req1,istatus,ierr)
   endif
      
   if(myid.ne.nproc-1) then
      call MPI_RECV(r(0,0,kparaend+1),(jx+2)*(jy+2), &
         MPI_REAL_SD,rightpem,tagrr, &
         MPI_COMM_WORLD,status,ierr)
   !         call MPI_WAIT(req2,istatus,ierr)
   endif

   if(myid.ne.nproc-1) then
      call MPI_SSEND(r(0,0,kparaend),(jx+2)*(jy+2), &
         MPI_REAL_SD,rightpem,tagrs, &
         MPI_COMM_WORLD,ierr)
   !         call MPI_WAIT(req3,istatus,ierr)
   endif
   if(myid.ne.0) then
      call MPI_RECV(r(0,0,kparasta-1),(jx+2)*(jy+2), &
         MPI_REAL_SD,leftpem,taglr, &
         MPI_COMM_WORLD,status,ierr)
   !         call MPI_WAIT(req4,istatus,ierr)
   endif
   !
   !     if periodic
   do k=1,1-kp

      call MPI_TYPE_VECTOR(jy+2,jx,jx+2,MPI_REAL_SD,plantype,ierr)
      call MPI_TYPE_COMMIT(plantype,ierr)


      if (myid.eq.nproc-1) then
         call MPI_SENDRECV(r(1,0,jz),1, &
            plantype,0,11, &
            r(1,0,jz+1),1, &
            plantype,0,12, &
            MPI_COMM_WORLD,status,ierr)


      endif

      if (myid.eq.0) then
         call MPI_SENDRECV(r(1,0,1),1, &
            plantype,nproc-1,12, &
            r(1,0,0),1, &
            plantype,nproc-1,11, &
            MPI_COMM_WORLD,status,ierr)

      endif

      call MPI_TYPE_FREE(plantype,ierr)

   enddo
   !hicco serve riapplicare periodicity per i etc...
   return
end
