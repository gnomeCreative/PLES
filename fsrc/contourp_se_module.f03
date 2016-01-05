module contourp_se_module

   use myarrays_velo3
   use mysending
   !
   use scala3
   use period
   use tipologia
   !
   use mpi

   implicit none

   private

   public :: contourp_se
contains

    !***********************************************************************
    subroutine contourp_se
        !***********************************************************************
        ! compute cartesian velocity and controvarian fluxes in periodic cells
        ! at time n+1
        ! at the corner the computation is made at the end
        !
        !
        !-----------------------------------------------------------------------
        ! arrays declaration

        integer i,j,k,kk,isc,err
        integer kpsta,kpend
        integer status(MPI_STATUS_SIZE),ierr, siz,MPI_UVW_TYPE
        !real, allocatable :: sendbbuf(:)
        !real, allocatable :: recvbbuf(:)
        !real, allocatable :: rr(:,:,:)
        !      real, allocatable, dimension(:,:) :: usendbuf,vsendbuf,wsendbuf,&
        !                                           urecvbuf,vrecvbuf,wrecvbuf
        real, dimension(1:jx,0:jy+1)  :: sendbuf,recvbuf
        !
        !-----------------------------------------------------------------------
           !allocate( sendbbuf( (nscal+3)*jx*(jy+2) ) )
           !allocate( recvbbuf( (nscal+3)*jx*(jy+2) ) )
        !call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        !call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
        siz = jx*(jy+2)
        !-----------------------------------------------------------------------
        ! set index for parallel loop in k
        if(myid.eq.0)then
            kpsta = kparasta+kp
            kpend = kparaend
        elseif(myid.eq.nproc-1)then
            kpsta = kparasta
            kpend = kparaend+1-kp
        else
            kpsta = kparasta
            kpend = kparaend
        end if
        !-----------------------------------------------------------------------
        !
        ! periodic cells in csi (also out of the domain)
        !
        do i=1,1-ip
            do k=kparasta,kparaend !1,jz
                do j=0,jy+1
                    !        face 1
                    u(0   ,j,k)=u(jx,j,k)
                    v(0   ,j,k)=v(jx,j,k)
                    w(0   ,j,k)=w(jx,j,k)
                    !        face 2
                    u(jx+1,j,k)=u(1 ,j,k)
                    v(jx+1,j,k)=v(1 ,j,k)
                    w(jx+1,j,k)=w(1 ,j,k)
                    !
                    do isc=1,nscal
                        rhov(isc,0   ,j,k)=rhov(isc,jx,j,k)
                        rhov(isc,jx+1,j,k)=rhov(isc,1 ,j,k)
                    end do
                !
                enddo
            enddo
        enddo
        !
        ! periodic cells in eta (also out of the domain)
        !
        do j=1,1-jp
            do i=ip,jx+1-ip
                do k=kpsta,kpend !kp,jz+1-kp
                    !        face 3
                    u(i,0   ,k)=u(i,jy,k)
                    v(i,0   ,k)=v(i,jy,k)
                    w(i,0   ,k)=w(i,jy,k)
                    !        face 4
                    u(i,jy+1,k)=u(i,1 ,k)
                    v(i,jy+1,k)=v(i,1 ,k)
                    w(i,jy+1,k)=w(i,1 ,k)

                    do isc=1,nscal
                        rhov(isc,i,0   ,k)=rhov(isc,i,jy,k)
                        rhov(isc,i,jy+1,k)=rhov(isc,i,1 ,k)
                    end do
                !
                enddo
            enddo
        enddo
        !
        ! periodic cells in zita (also out of the domain)
        !

        !allocate(recvbuf(1:jx, 0:jy+1),stat=err)
        !allocate(sendbuf(1:jx, 0:jy+1),stat=err)
        do k=1,1-kp
            !THIS WAS THE OLD WAY, KEPT FOR NOSTALGIA
            !#include "old_parallel.h"

            !call MPI_TYPE_VECTOR(jx, 1, jy+2, MPI_REAL_SD,MPI_UVW_TYPE, ierr)
            call MPI_TYPE_VECTOR(jy+2, jx, jx, MPI_REAL_SD,MPI_UVW_TYPE, ierr)
            call MPI_TYPE_COMMIT(MPI_UVW_TYPE, ierr)
            !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            if (myid==0) then

                !------------------------------------------------------------------------
                ! SEND U
                do i = 1, jx
                    do j = 0, jy+1
                        sendbuf(i,j  ) = u(i, j, 0)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,193,                 &
                    & MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND V
                do i = 1, jx
                    do j = 0, jy+1
                        sendbuf(i,j  ) = v(i, j, 0)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,194,                 &
                    & MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND W
                do i = 1, jx
                    do j = 0, jy+1
                        sendbuf(i,j  ) = w(i, j, 0)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,195,                 &
                    & MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND RHOV
                do isc = 1, nscal
                    do i = 1, jx
                        do j = 0, jy+1
                            sendbuf(i,j  ) = rhov(isc,i,j,0)
                        end do
                    end do
                    call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,196+isc,    &
                        & MPI_COMM_WORLD,ierr)
                end do

                !------------------------------------------------------------------------
                ! RECV U
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,nproc-1,182,           &
                    & MPI_COMM_WORLD,status,ierr)
                !call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,nproc-1,182,           &
                !   & MPI_COMM_WORLD,status,ierr)
                do i = 1, jx
                    do j = 0, jy + 1
                        u(i, j, 0) = recvbuf(i, j  )
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV V
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,nproc-1,183,           &
                    & MPI_COMM_WORLD,status,ierr)
                !call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,nproc-1,183,           &
                !   & MPI_COMM_WORLD,status,ierr)
                do i = 1, jx
                    do j = 0, jy + 1
                        v(i, j, 0) = recvbuf(i, j  )
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV W
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,nproc-1,184,           &
                    & MPI_COMM_WORLD,status,ierr)
                !call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,nproc-1,184,           &
                !   & MPI_COMM_WORLD,status,ierr)
                do i = 1, jx
                    do j = 0, jy + 1
                        w(i, j, 0) = recvbuf(i, j  )
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV RHOV
                do isc = 1, nscal
                    call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,nproc-1,185+isc,   &
                        & MPI_COMM_WORLD,status,ierr)
                    !call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,nproc-1,185+isc,   &
                    !   & MPI_COMM_WORLD,status,ierr)
                    do i = 1, jx
                        do j = 0, jy + 1
                            rhov(isc,i,j,0) = recvbuf(i,j  )
                        end do
                    end do
                end do

            !--------------------------------------------------------------------------------------
            else if (myid == nproc-1) then
                !--------------------------------------------------------------------------------------

                !------------------------------------------------------------------------
                ! RECV U
                !call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,0,193,                 &
                !   & MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,193,                 &
                    & MPI_COMM_WORLD,status,ierr)
                do i = 1, jx
                    do j = 0, jy + 1
                        u(i, j, jz+1) = recvbuf(i, j  )
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV V
                !call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,0,194,                 &
                !   & MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,194,                 &
                    & MPI_COMM_WORLD,status,ierr)
                do i = 1, jx
                    do j = 0, jy + 1
                        v(i, j, jz+1) = recvbuf(i, j  )
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV W
                !call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,0,195,                 &
                !   & MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,195,                 &
                    & MPI_COMM_WORLD,status,ierr)
                do i = 1, jx
                    do j = 0, jy + 1
                        w(i, j, jz+1) = recvbuf(i, j  )
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV  RHOV
                do isc = 1, nscal
                    call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,196+isc,    &
                        & MPI_COMM_WORLD,status,ierr)
                    !call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,0,196+isc,    &
                    !   & MPI_COMM_WORLD,status,ierr)
                    do i = 1, jx
                        do j = 0, jy + 1
                            rhov(isc,i,j,jz+1) = recvbuf(i,j  )
                        end do
                    end do
                end do

                !------------------------------------------------------------------------
                ! SEND U
                do i = 1, jx
                    do j = 0, jy+1
                        sendbuf(i,j  ) = u(i, j, jz+1)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,182,           &
                    & MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND V
                do i = 1, jx
                    do j = 0, jy+1
                        sendbuf(i,j  ) = v(i, j, jz+1)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,183,           &
                    & MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND W
                do i = 1, jx
                    do j = 0, jy+1
                        sendbuf(i,j  ) = w(i, j, jz+1)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,184,           &
                    & MPI_COMM_WORLD,ierr)


                !------------------------------------------------------------------------
                ! SEND RHOV
                do isc = 1, nscal
                    !if(.not.allocated(rhosendbuf)) then
                    !   call MPI_ABORT(MPI_COMM_WORLD, err, ierr)
                    !   stop "CARE CHIMBA ERROR"
                    !end if
                    do i = 1, jx
                        do j = 0, jy+1
                            sendbuf(i,j  ) = rhov(isc,i,j,jz+1)
                        end do
                    end do

                    call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,185+isc,   &
                        & MPI_COMM_WORLD,ierr)

                end do


            end if
            call MPI_TYPE_FREE(MPI_UVW_TYPE,ierr)
        enddo
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        ! at the corner (temporary version)
        !
        do i=1,1-ip
            do j=0,jy
                if(myid.eq.0)then
                    u(0   ,j,   0)=u(jx,j,   0)
                    u(jx+1,j,   0)=u(1 ,j,   0)

                    v(0   ,j,   0)=v(jx,j,   0)
                    v(jx+1,j,   0)=v(1 ,j,   0)

                    w(0   ,j,   0)=w(jx,j,   0)
                    w(jx+1,j,   0)=w(1 ,j,   0)

                    do isc=1,nscal
                        rhov(isc,0   ,j,   0)=rhov(isc,jx,j,   0)
                        rhov(isc,jx+1,j,   0)=rhov(isc,1 ,j,   0)
                    end do
                end if

                if(myid.eq.nproc-1)then
                    u(0   ,j,jz+1)=u(jx,j,jz+1)
                    u(jx+1,j,jz+1)=u(1 ,j,jz+1)

                    v(0   ,j,jz+1)=v(jx,j,jz+1)
                    v(jx+1,j,jz+1)=v(1 ,j,jz+1)

                    w(0   ,j,jz+1)=w(jx,j,jz+1)
                    w(jx+1,j,jz+1)=w(1 ,j,jz+1)

                    do isc=1,nscal
                        rhov(isc,0   ,j,jz+1)=rhov(isc,jx,j,jz+1)
                        rhov(isc,jx+1,j,jz+1)=rhov(isc,1 ,j,jz+1)
                    end do
                end if
            enddo
        enddo

    !      if(allocated(sendbbuf))deallocate(sendbbuf)
    !      if(allocated(recvbbuf))deallocate(recvbbuf)
    !
    end subroutine contourp_se

!   !***********************************************************************
!   subroutine contourp_se
!      !***********************************************************************
!      ! compute cartesian velocity and controvarian fluxes in periodic cells
!      ! at time n+1
!      ! at the corner the computation is made at the end
!      !
!      !
!      !-----------------------------------------------------------------------
!      ! arrays declaration
!
!      integer i,j,k,kk,isc
!      integer kpsta,kpend
!      integer status(MPI_STATUS_SIZE),ierr
!      real, allocatable :: sendbuf(:)
!      real, allocatable :: recvbuf(:)
!      real, allocatable :: rr(:,:,:)
!      !
!      !-----------------------------------------------------------------------
!      allocate( sendbuf( (nscal+3)*jx*(jy+2) ) )
!      allocate( recvbuf( (nscal+3)*jx*(jy+2) ) )
!      !-----------------------------------------------------------------------
!      ! set index for parallel loop in k
!      if(myid.eq.0)then
!         kpsta = kparasta+kp
!         kpend = kparaend
!      elseif(myid.eq.nproc-1)then
!         kpsta = kparasta
!         kpend = kparaend+1-kp
!      else
!         kpsta = kparasta
!         kpend = kparaend
!      end if
!      !-----------------------------------------------------------------------
!      !
!      ! periodic cells in csi (also out of the domain)
!      !
!      do i=1,1-ip
!         do k=kparasta,kparaend !1,jz
!            do j=0,jy+1
!               !        face 1
!               u(0   ,j,k)=u(jx,j,k)
!               v(0   ,j,k)=v(jx,j,k)
!               w(0   ,j,k)=w(jx,j,k)
!               !        face 2
!               u(jx+1,j,k)=u(1 ,j,k)
!               v(jx+1,j,k)=v(1 ,j,k)
!               w(jx+1,j,k)=w(1 ,j,k)
!               !
!               do isc=1,nscal
!                  rhov(isc,0   ,j,k)=rhov(isc,jx,j,k)
!                  rhov(isc,jx+1,j,k)=rhov(isc,1 ,j,k)
!               end do
!            !
!            enddo
!         enddo
!      enddo
!      !
!      ! periodic cells in eta (also out of the domain)
!      !
!      do j=1,1-jp
!         do i=ip,jx+1-ip
!            do k=kpsta,kpend !kp,jz+1-kp
!               !        face 3
!               u(i,0   ,k)=u(i,jy,k)
!               v(i,0   ,k)=v(i,jy,k)
!               w(i,0   ,k)=w(i,jy,k)
!               !        face 4
!               u(i,jy+1,k)=u(i,1 ,k)
!               v(i,jy+1,k)=v(i,1 ,k)
!               w(i,jy+1,k)=w(i,1 ,k)
!
!               do isc=1,nscal
!                  rhov(isc,i,0   ,k)=rhov(isc,i,jy,k)
!                  rhov(isc,i,jy+1,k)=rhov(isc,i,1 ,k)
!               end do
!            !
!            enddo
!         enddo
!      enddo
!      !
!      ! periodic cells in zita (also out of the domain)
!      !
!      do k=1,1-kp
!
!         !     needs to send plane jz and 1 for periodicity
!
!         !     initialize send/recv buffer
!         sendbuf=0.
!         recvbuf=0.
!
!         !     data upload in send buffer
!         !     for scalar  use a temporary array rr in order to
!         !     use the same subroutine
!
!         if (myid.eq.nproc-1) then
!            allocate(rr(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
!            rr = 0.
!            call create_sendbuf(sendbuf,u,1,jz)
!            call create_sendbuf(sendbuf,v,2,jz)
!            call create_sendbuf(sendbuf,w,3,jz)
!            do isc=1,nscal
!               do kk=kparasta,kparaend
!                  do j=0,jy+1
!                     do i=0,jx+1
!                        rr(i,j,kk)=rhov(isc,i,j,kk)
!                     end do
!                  end do
!               end do
!               call create_sendbuf(sendbuf,rr,3+isc,jz)
!               rr = 0.
!            end do
!          !deallocate(rr)
!         else if (myid.eq.0) then
!            allocate(rr(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
!            rr = 0.
!            call create_sendbuf(sendbuf,u,1,1 )
!            call create_sendbuf(sendbuf,v,2,1 )
!            call create_sendbuf(sendbuf,w,3,1 )
!
!            do isc=1,nscal
!               do kk=kparasta,kparaend
!                  do j=0,jy+1
!                     do i=0,jx+1
!                        rr(i,j,kk)=rhov(isc,i,j,kk)
!                     end do
!                  end do
!               end do
!               call create_sendbuf(sendbuf,rr,3+isc,1 )
!               rr = 0.
!            end do
!          !deallocate(rr)
!         endif
!         !     send the buffer so that 0 known plane jz and nproc-1 known plane 1
!         if (myid.eq.nproc-1) then
!            call MPI_SENDRECV(sendbuf(1),(nscal+3)*(jx)*(jy+2),MPI_REAL_SD,&
!               &           0,9001,recvbuf(1),(nscal+3)*(jx)*(jy+2),MPI_REAL_SD,&
!               &           0,8001,MPI_COMM_WORLD,status,ierr)
!
!         else if (myid.eq.0) then
!            call MPI_SENDRECV(sendbuf(1),(nscal+3)*(jx)*(jy+2),MPI_REAL_SD,&
!               &    nproc-1, 8001,recvbuf(1),(nscal+3)*(jx)*(jy+2),MPI_REAL_SD,&
!               &    nproc-1, 9001,MPI_COMM_WORLD,status,ierr)
!
!         endif
!
!
!         !     download data of recive buffer in temporary plane
!         !     for scalar  use a temporary array rr in order to
!         !     use the same subroutine
!         if(allocated(rr)) deallocate(rr)
!         if(myid.eq.0)then
!            allocate(rr(0:jx+1,0:jy+1,jz:jz))
!            call make_from_recvbuf(recvbuf,u_piano,1,jz)
!            call make_from_recvbuf(recvbuf,v_piano,2,jz)
!            call make_from_recvbuf(recvbuf,w_piano,3,jz)
!            do isc=1,nscal
!               call make_from_recvbuf(recvbuf,rr,3+isc,jz)
!               do j=0,jy+1
!                  do i=0,jx+1
!                     rhov_piano(isc,i,j,jz) = rr(i,j,jz)
!                  end do
!               end do
!               rr = 0.
!            end do
!           !deallocate(rr)
!         elseif(myid.eq.nproc-1)then
!            allocate(rr(0:jx+1,0:jy+1,1:1))
!            call make_from_recvbuf(recvbuf,u_piano,1,1)
!            call make_from_recvbuf(recvbuf,v_piano,2,1)
!            call make_from_recvbuf(recvbuf,w_piano,3,1)
!            do isc=1,nscal
!               call make_from_recvbuf(recvbuf,rr,3+isc,1)
!               do j=0,jy+1
!                  do i=0,jx+1
!                     rhov_piano(isc,i,j,1) = rr(i,j,1)
!                  end do
!               end do
!               rr = 0.
!            end do
!           !deallocate(rr)
!         endif
!         if(allocated(rr)) deallocate(rr)
!
!         !     apply periodicity
!         if(myid.eq.0)then
!            do i=1,jx
!               do j=0,jy+1
!                  !
!                  u(i,j,   0)=u_piano(i,j,jz)
!                  v(i,j,   0)=v_piano(i,j,jz)
!                  w(i,j,   0)=w_piano(i,j,jz)
!
!                  do isc=1,nscal
!                     rhov(isc,i,j,   0)=rhov_piano(isc,i,j,jz)
!                  end do
!               !
!               enddo
!            enddo
!         end if
!
!         if(myid.eq.nproc-1)then
!            do i=1,jx
!               do j=0,jy+1
!                  !
!                  u(i,j,jz+1)=u_piano(i,j, 1)
!                  v(i,j,jz+1)=v_piano(i,j, 1)
!                  w(i,j,jz+1)=w_piano(i,j, 1)
!
!                  do isc=1,nscal
!                     rhov(isc,i,j,jz+1)=rhov_piano(isc,i,j, 1)
!                  end do
!               !
!               enddo
!            enddo
!         end if
!
!      enddo
!      !
!      ! at the corner (temporary version)
!      !
!      do i=1,1-ip
!         do j=0,jy
!            if(myid.eq.0)then
!               u(0   ,j,   0)=u(jx,j,   0)
!               u(jx+1,j,   0)=u(1 ,j,   0)
!
!               v(0   ,j,   0)=v(jx,j,   0)
!               v(jx+1,j,   0)=v(1 ,j,   0)
!
!               w(0   ,j,   0)=w(jx,j,   0)
!               w(jx+1,j,   0)=w(1 ,j,   0)
!
!               do isc=1,nscal
!                  rhov(isc,0   ,j,   0)=rhov(isc,jx,j,   0)
!                  rhov(isc,jx+1,j,   0)=rhov(isc,1 ,j,   0)
!               end do
!            end if
!
!            if(myid.eq.nproc-1)then
!               u(0   ,j,jz+1)=u(jx,j,jz+1)
!               u(jx+1,j,jz+1)=u(1 ,j,jz+1)
!
!               v(0   ,j,jz+1)=v(jx,j,jz+1)
!               v(jx+1,j,jz+1)=v(1 ,j,jz+1)
!
!               w(0   ,j,jz+1)=w(jx,j,jz+1)
!               w(jx+1,j,jz+1)=w(1 ,j,jz+1)
!
!               do isc=1,nscal
!                  rhov(isc,0   ,j,jz+1)=rhov(isc,jx,j,jz+1)
!                  rhov(isc,jx+1,j,jz+1)=rhov(isc,1 ,j,jz+1)
!               end do
!            end if
!         enddo
!      enddo
!
!      deallocate(sendbuf)
!      deallocate(recvbuf)
!   !
!   end subroutine contourp_se
!
!
!
!   !***********************************************************************
!   subroutine create_sendbuf(sendbuf,var,n,kest)
!      !***********************************************************************
!      ! upload data in sending buffer
!      !
!      !use mysending
!      implicit none
!      !include 'scala3.h'
!      !-----------------------------------------------------------------------
!      integer i,j,m
!      integer, intent(in) :: n,kest
!      real, intent(in out) :: sendbuf(:)!( (nscal+3)*jx*(jy+2) )
!      real, intent(in) :: var(0:,0:,kparasta-deepl:)!(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
!      !-----------------------------------------------------------------------
!      do j=0,jy+1
!         do i=1,jx
!            m=(n-1)*(jx)*(jy+2)+i+(jx)*j
!            sendbuf(m)=var(i,j,kest)
!         enddo
!      enddo
!   !
!   end      subroutine create_sendbuf
!
!   !***********************************************************************
!   subroutine make_from_recvbuf(recvbuf,var,n,kest)
!      !***********************************************************************
!      ! download data from recived buffer
!      !
!      !use mysending
!      implicit none
!      !include 'scala3.h'
!      !-----------------------------------------------------------------------
!      integer, intent(in) :: n,kest
!      integer i,j,m
!      real recvbuf(:)!( (nscal+3)*jx*(jy+2) )
!      real var(0:,0:,kest:)!(0:n1+1,0:n2+1,kest:kest)
!      !-----------------------------------------------------------------------
!      do j=0,jy+1
!         do i=1,jx
!            m=(n-1)*(jx)*(jy+2)+i+(jx)*j
!            var(i,j,kest)=recvbuf(m)
!         end do
!      end do
!   !
!   end subroutine make_from_recvbuf

end module contourp_se_module
