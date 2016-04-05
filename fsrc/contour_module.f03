module contour_module

    use myarrays_velo3
    use mysending
    use orlansky_module
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    integer,bind(C) :: iboun1,iboun2
    integer,bind(C) :: iboun3,iboun4
    integer,bind(C) :: iboun5,iboun6

    integer,bind(C) :: ibodybuffer1,ibodybuffer2
    integer,bind(C) :: ibodybuffer3,ibodybuffer4
    integer,bind(C) :: ibodybuffer5,ibodybuffer6

    ! arrays for inflow
    real,allocatable :: up1(:,:),vp1(:,:),wp1(:,:),rhovp1(:,:,:)
    real,allocatable :: up2(:,:),vp2(:,:),wp2(:,:),rhovp2(:,:,:)
    real,allocatable :: up3(:,:),vp3(:,:),wp3(:,:),rhovp3(:,:,:)
    real,allocatable :: up4(:,:),vp4(:,:),wp4(:,:),rhovp4(:,:,:)
    real,allocatable :: up5(:,:),vp5(:,:),wp5(:,:),rhovp5(:,:,:)
    real,allocatable :: up6(:,:),vp6(:,:),wp6(:,:),rhovp6(:,:,:)

    real,allocatable,dimension(:,:,:)::rhovo1,rhovo2,rhovo5,rhovo6
    real,allocatable,dimension(:,:,:)::rhovn1,rhovn2,rhovn5,rhovn6

    real,allocatable,dimension(:,:) :: ucp1
    real,allocatable,dimension(:,:) :: ucp2
    real,allocatable,dimension(:,:) :: wcp5
    real,allocatable,dimension(:,:) :: wcp6

    real, allocatable :: tke1(:,:),tke2(:,:)
    real, allocatable :: tke5(:,:),tke6(:,:)
    real, allocatable :: tkepom1(:,:,:),tkepom2(:,:,:)
    real, allocatable :: tkepom5(:,:,:),tkepom6(:,:,:)

contains

    subroutine contourp_se()
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
        if (myid==0) then
            kpsta = kparasta+kp
            kpend = kparaend
        else if (myid==nproc-1) then
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
                end do
            end do
        end do
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
                end do
            end do
        end do
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
            else if (myid==nproc-1) then
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
                    !if (.not.allocated(rhosendbuf)) then
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
        end do
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        ! at the corner (temporary version)
        !
        do i=1,1-ip
            do j=0,jy
                if (myid==0) then
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

                if (myid==nproc-1) then
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
            end do
        end do

    !      if (allocated(sendbbuf))deallocate(sendbbuf)
    !      if (allocated(recvbbuf))deallocate(recvbbuf)
    !
    end subroutine contourp_se

    subroutine contour_se_nesting()
        !***********************************************************************
        ! boundary condition on cartesian velocities and scalars and
        ! for controvariant fluxes
        !
        use myarrays_metri3

        implicit none
        !
        !-----------------------------------------------------------------------
        ! arrays declaration
        !
        integer i,j,k
        integer ii,jj,kk
        integer isc
        !
        !-----------------------------------------------------------------------
        !     face 1 'sinistra' (for periodicity see vel_up)
        do ii=1,ip
            if (infout1==0) then  !inflow
                do k=kparasta,kparaend !1,jz
                    do j=1,jy
                        u(0,j,k)=up1(j,k)
                        v(0,j,k)=vp1(j,k)
                        w(0,j,k)=wp1(j,k)
                        do isc=1,nscal
                            rhov(isc,0,j,k)=rhovp1(isc,j,k)
                        end do
                        uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 1'
            else if (infout1==1) then  !outflow
                do k=kparasta,kparaend
                    do j=1,jy
                        u(0,j,k)=u(0,j,k)
                        v(0,j,k)=v(0,j,k)
                        w(0,j,k)=w(0,j,k)
                        do isc=1,nscal
                            rhov(isc,0,j,k)=rhov(isc,1,j,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 1'
            else if (infout1==2) then  !wall
                if (iboun1==0) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(0,j,k)=0.
                            v(0,j,k)=0.
                            w(0,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,0,j,k)=rhov(isc,1,j,k)
                            end do
                            uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                        end do
                    end do
                else if (iboun1==1) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(0,j,k)=0.
                            v(0,j,k)=v(1,j,k)
                            w(0,j,k)=w(1,j,k)
                            do isc=1,nscal
                                rhov(isc,0,j,k)=rhov(isc,1,j,k)
                            end do
                            uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                        end do
                    end do
                else if (iboun1==2) then
                    do k=kparasta,kparaend
                        do j=1,jy

                            u(0,j,k)= ucp1(j,k)*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                            !              v(0,j,k)= ucp1(j,k)*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                            w(0,j,k)= ucp1(j,k)*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                            v(0,j,k)=v(1,j,k)
                            !              w(0,j,k)=w(1,j,k)

                            uc(0,j,k)=ucp1(j,k)

                            do isc=1,nscal
                                if (potenziale==0) then
                                    rhov(isc,0,j,k)=rhovp1(isc,j,k)
                                else
                                    rhov(isc,0,j,k)=rhovo1(isc,j,k)
                                end if
                            end do

                        end do
                    end do
                end if
                if (myid==0)write(*,*)'solid wall face 1'
            end if

            !     face 2 'destra' (for periodicity see vel_up)
            if (infout2==0) then  !inflow
                do k=kparasta,kparaend
                    do j=1,jy
                        u(jx+1,j,k)=up2(j,k)
                        v(jx+1,j,k)=vp2(j,k)
                        w(jx+1,j,k)=wp2(j,k)
                        do isc=1,nscal
                            rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
                        end do
                        uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 2'
            else if (infout2==1) then  !outflow
                do k=kparasta,kparaend
                    do j=1,jy
                        u(jx+1,j,k)=u(jx+1,j,k)
                        v(jx+1,j,k)=v(jx+1,j,k)
                        w(jx+1,j,k)=w(jx+1,j,k)
                        do isc=1,nscal
                            rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 2'
            else if (infout2==2) then  !wall
                if (iboun2==0) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(jx+1,j,k)=0.
                            v(jx+1,j,k)=0.
                            w(jx+1,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                            end do
                            uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                        end do
                    end do
                else if (iboun2==1) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(jx+1,j,k)=0.
                            v(jx+1,j,k)=v(jx,j,k)
                            w(jx+1,j,k)=w(jx,j,k)
                            do isc=1,nscal
                                rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                            end do
                            uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                        end do
                    end do
                else if (iboun2==2) then
                    do k=kparasta,kparaend
                        do j=1,jy

                            u(jx+1,j,k)= ucp2(j,k)*csx(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                            !             v(jx+1,j,k)= ucp2(j,k)*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                            w(jx+1,j,k)= ucp2(j,k)*csz(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                            v(jx+1,j,k)=v(jx,j,k)
                            !             w(jx+1,j,k)=w(jx,j,k)

                            uc(jx,j,k)=ucp2(j,k)

                            do isc=1,nscal
                                if (potenziale==0) then
                                    rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
                                else
                                    rhov(isc,jx+1,j,k)=rhovo2(isc,j,k)
                                end if
                            end do

                        end do
                    end do
                end if
                if (myid==0)write(*,*)'solid wall face 2'
            end if
        end do   !end loop  ii=1,ip
        !
        !-----------------------------------------------------------------------
        !     face 3 'sotto' (for periodicity see vel_up)
        !
        do jj=1,jp
            if (infout3==0) then  !inflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,0,k)=up3(i,k)
                        v(i,0,k)=vp3(i,k)
                        w(i,0,k)=wp3(i,k)
                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhovp3(isc,i,k)
                        end do
                        vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 3'
            else if (infout3==1) then  !outflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,0,k)=u(i,0,k)
                        v(i,0,k)=v(i,0,k)
                        w(i,0,k)=w(i,0,k)
                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhov(isc,i,1,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 3'
            else if (infout3==2) then  !wall
                if (iboun3==0) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=0.
                            v(i,0,k)=0.
                            w(i,0,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                else if (iboun3==1) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=u(i,1,k)
                            v(i,0,k)=0.
                            w(i,0,k)=w(i,1,k)
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                else if (iboun3==2) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=0.
                            v(i,0,k)=0.
                            w(i,0,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                end if
                if (myid==0)write(*,*)'solid wall face 3'
            end if

            !     face 4 'sopra' (for periodicity see vel_up)
            if (infout4==0) then  !inflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,jy+1,k)=up4(i,k)
                        v(i,jy+1,k)=vp4(i,k)
                        w(i,jy+1,k)=wp4(i,k)
                        do isc=1,nscal
                            rhov(isc,i,jy+1,k)=rhovp4(isc,i,k)
                        end do
                        vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 4'
            else if (infout4==1) then  !outflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,jy+1,k)=u(i,jy+1,k)
                        v(i,jy+1,k)=v(i,jy+1,k)
                        w(i,jy+1,k)=w(i,jy+1,k)
                        do isc=1,nscal
                            rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 4'
            else if (infout4==2) then  !wall
                if (iboun4==0) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k)=0.
                            v(i,jy+1,k)=0.
                            w(i,jy+1,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                            vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                else if (iboun4==1) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k)=u(i,jy,k)
                            v(i,jy+1,k)=0.
                            w(i,jy+1,k)=w(i,jy,k)
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                            vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                else if (iboun4==2) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k)=0.
                            v(i,jy+1,k)=0.
                            w(i,jy+1,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                            vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                end if
                if (myid==0)write(*,*)'solid wall face 4'
            end if
        end do !end loop  jj=1,jp
        !
        !-----------------------------------------------------------------------
        !     face 5 'indietro' (for periodicity see vel_up)
        !
        do kk=1,kp
            if (myid==0) then
                if (infout5==0) then  !inflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,0)=up5(i,j)
                            v(i,j,0)=vp5(i,j)
                            w(i,j,0)=wp5(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,0)=rhovp5(isc,i,j)
                            end do
                            wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                        end do
                    end do
                    write(*,*)'inflow face 5'
                else if (infout5==1) then  !outflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,0)=u(i,j,0)
                            v(i,j,0)=v(i,j,0)
                            w(i,j,0)=w(i,j,0)
                            do isc=1,nscal
                                rhov(isc,i,j,0)=rhov(isc,i,j,1)
                            end do
                        end do
                    end do
                    write(*,*)'outflow face 5'
                else if (infout5==2) then  !wall
                    if (iboun5==0) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,0)=0.
                                v(i,j,0)=0.
                                w(i,j,0)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,0)=rhov(isc,i,j,1)
                                end do
                                wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                            end do
                        end do
                    else if (iboun5==1) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,0)=u(i,j,1)
                                v(i,j,0)=v(i,j,1)
                                w(i,j,0)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,0)=rhov(isc,i,j,1)
                                end do
                                wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                            end do
                        end do
                    else if (iboun5==2) then
                        do j=1,jy
                            do i=1,jx
                                !             u(i,j,0)=u(i,j,1)
                                v(i,j,0)=v(i,j,1)

                                u(i,j,0)= wcp5(i,j)*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                !             v(i,j,0)= wcp5(i,j)*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                w(i,j,0)= wcp5(i,j)*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                wc(i,j,0)=wcp5(i,j)

                                do isc=1,nscal
                                    if (potenziale==0) then
                                        rhov(isc,i,j,0)=rhovp5(isc,i,j)
                                    else
                                        rhov(isc,i,j,0)=rhovo5(isc,i,j)
                                    end if
                                end do

                            end do
                        end do
                    end if
                    write(*,*)'solid wall face 5'
                end if
            end if ! myid=0
            !................................................
            !     face 6 'avanti' (for periodicity see vel_up)
            !
            if (myid==nproc-1) then
                if (infout6==0) then  !inflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,jz+1)=up6(i,j)
                            v(i,j,jz+1)=vp6(i,j)
                            w(i,j,jz+1)=wp6(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)
                            end do
                            wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)

                        end do
                    end do
                    write(*,*)'inflow face 6'
                else if (infout6==1) then  !outflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,jz+1)=u(i,j,jz+1)
                            v(i,j,jz+1)=v(i,j,jz+1)
                            w(i,j,jz+1)=w(i,j,jz+1)
                            do isc=1,nscal
                                rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                            end do
                        end do
                    end do
                    write(*,*)'outflow face 6'
                else if (infout6==2) then  !wall
                    if (iboun6==0) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,jz+1)=0.
                                v(i,j,jz+1)=0.
                                w(i,j,jz+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                                end do
                                wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                            end do
                        end do
                    else if (iboun6==1) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,jz+1)=u(i,j,jz)
                                v(i,j,jz+1)=v(i,j,jz)
                                w(i,j,jz+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                                end do
                                wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                            end do
                        end do
                    else if (iboun6==2) then
                        do j=1,jy
                            do i=1,jx
                                !             u(i,j,jz+1)=u(i,j,jz)
                                v(i,j,jz+1)=v(i,j,jz)

                                u(i,j,jz+1)= wcp6(i,j)*ztx(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                !             v(i,j,jz+1)= wcp6(i,j)*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                w(i,j,jz+1)= wcp6(i,j)*ztz(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                wc(i,j,jz)=wcp6(i,j)

                                do isc=1,nscal
                                    if (potenziale==0) then
                                        rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)
                                    else
                                        rhov(isc,i,j,jz+1)=rhovo6(isc,i,j)
                                    end if
                                end do

                            end do
                        end do
                    end if
                    write(*,*)'solid wall face 6'
                end if
            end if !myid=nproc-1
        end do !end loop  kk=1,kp
        !-----------------------------------------------------------------------

        return

    end subroutine contour_se_nesting

    subroutine contour_se()
        !***********************************************************************
        ! boundary condition on cartesian velocities and scalars and
        ! for controvariant fluxes
        !
        use myarrays_metri3

        implicit none
        !
        !-----------------------------------------------------------------------
        ! arrays declaration
        !
        integer i,j,k
        integer ii,jj,kk
        integer isc
        !
        !-----------------------------------------------------------------------
        !     face 1 'sinistra' (for periodicity see vel_up)
        do ii=1,ip
            if (infout1==0) then  !inflow
                do k=kparasta,kparaend !1,jz
                    do j=1,jy
                        u(0,j,k)=up1(j,k)
                        v(0,j,k)=vp1(j,k)
                        w(0,j,k)=wp1(j,k)
                        do isc=1,nscal
                            rhov(isc,0,j,k)=rhovp1(isc,j,k)
                        end do
                        uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 1'
            else if (infout1==1) then  !outflow
                do k=kparasta,kparaend
                    do j=1,jy
                        u(0,j,k)=u(0,j,k)
                        v(0,j,k)=v(0,j,k)
                        w(0,j,k)=w(0,j,k)
                        do isc=1,nscal
                            rhov(isc,0,j,k)=rhov(isc,1,j,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 1'
            else if (infout1==2) then  !wall
                if (iboun1==0) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(0,j,k)=0.
                            v(0,j,k)=0.
                            w(0,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,0,j,k)=rhov(isc,1,j,k)
                            end do
                            uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                        end do
                    end do
                else if (iboun1==1) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(0,j,k)=0.
                            v(0,j,k)=v(1,j,k)
                            w(0,j,k)=w(1,j,k)
                            do isc=1,nscal
                                rhov(isc,0,j,k)=rhov(isc,1,j,k)
                            end do
                            uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                        end do
                    end do
                else if (iboun1==2) then

                end if
                if (myid==0)write(*,*)'solid wall face 1'
            end if

            !     face 2 'destra' (for periodicity see vel_up)
            if (infout2==0) then  !inflow
                do k=kparasta,kparaend
                    do j=1,jy
                        u(jx+1,j,k)=up2(j,k)
                        v(jx+1,j,k)=vp2(j,k)
                        w(jx+1,j,k)=wp2(j,k)
                        do isc=1,nscal
                            rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
                        end do
                        uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 2'
            else if (infout2==1) then  !outflow
                do k=kparasta,kparaend
                    do j=1,jy
                        u(jx+1,j,k)=u(jx+1,j,k)
                        v(jx+1,j,k)=v(jx+1,j,k)
                        w(jx+1,j,k)=w(jx+1,j,k)
                        do isc=1,nscal
                            rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 2'
            else if (infout2==2) then  !wall
                if (iboun2==0) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(jx+1,j,k)=0.
                            v(jx+1,j,k)=0.
                            w(jx+1,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                            end do
                            uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                        end do
                    end do
                else if (iboun2==1) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(jx+1,j,k)=0.
                            v(jx+1,j,k)=v(jx,j,k)
                            w(jx+1,j,k)=w(jx,j,k)
                            do isc=1,nscal
                                rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                            end do
                            uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                        end do
                    end do
                else if (iboun2==2) then

                end if
                if (myid==0)write(*,*)'solid wall face 2'
            end if
        end do   !end loop  ii=1,ip
        !
        !-----------------------------------------------------------------------
        !     face 3 'sotto' (for periodicity see vel_up)
        !
        do jj=1,jp
            if (infout3==0) then  !inflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,0,k)=up3(i,k)
                        v(i,0,k)=vp3(i,k)
                        w(i,0,k)=wp3(i,k)
                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhovp3(isc,i,k)
                        end do
                        vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 3'
            else if (infout3==1) then  !outflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,0,k)=u(i,0,k)
                        v(i,0,k)=v(i,0,k)
                        w(i,0,k)=w(i,0,k)
                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhov(isc,i,1,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 3'
            else if (infout3==2) then  !wall
                if (iboun3==0) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=0.
                            v(i,0,k)=0.
                            w(i,0,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                else if (iboun3==1) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=u(i,1,k)
                            v(i,0,k)=0.
                            w(i,0,k)=w(i,1,k)
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                else if (iboun3==2) then

                end if
                if (myid==0)write(*,*)'solid wall face 3'
            end if

            !     face 4 'sopra' (for periodicity see vel_up)
            if (infout4==0) then  !inflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,jy+1,k)=up4(i,k)
                        v(i,jy+1,k)=vp4(i,k)
                        w(i,jy+1,k)=wp4(i,k)
                        do isc=1,nscal
                            rhov(isc,i,jy+1,k)=rhovp4(isc,i,k)
                        end do
                        vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 4'
            else if (infout4==1) then  !outflow
                do k=kparasta,kparaend
                    do i=1,jx
                        u(i,jy+1,k)=u(i,jy+1,k)
                        v(i,jy+1,k)=v(i,jy+1,k)
                        w(i,jy+1,k)=w(i,jy+1,k)
                        do isc=1,nscal
                            rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 4'
            else if (infout4==2) then  !wall
                if (iboun4==0) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k)=1.
                            v(i,jy+1,k)=0.
                            w(i,jy+1,k)=1.
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                            vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k) +ety(i,jy,k)*v(i,jy+1,k) +etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                else if (iboun4==1) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k)=u(i,jy,k)
                            v(i,jy+1,k)=0.
                            w(i,jy+1,k)=w(i,jy,k)
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                            vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                else if (iboun4==2) then

                end if
                if (myid==0)write(*,*)'solid wall face 4'
            end if
        end do !end loop  jj=1,jp
        !
        !-----------------------------------------------------------------------
        !     face 5 'indietro' (for periodicity see vel_up)
        !
        do kk=1,kp
            if (myid==0) then
                if (infout5==0) then  !inflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,0)=up5(i,j)
                            v(i,j,0)=vp5(i,j)
                            w(i,j,0)=wp5(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,0)=rhovp5(isc,i,j)
                            end do
                            wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                        end do
                    end do
                    write(*,*)'inflow face 5'
                else if (infout5==1) then  !outflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,0)=u(i,j,0)
                            v(i,j,0)=v(i,j,0)
                            w(i,j,0)=w(i,j,0)
                            do isc=1,nscal
                                rhov(isc,i,j,0)=rhov(isc,i,j,1)
                            end do
                        end do
                    end do
                    write(*,*)'outflow face 5'
                else if (infout5==2) then  !wall
                    if (iboun5==0) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,0)=0.
                                v(i,j,0)=0.
                                w(i,j,0)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,0)=rhov(isc,i,j,1)
                                end do
                                wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                            end do
                        end do
                    else if (iboun5==1) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,0)=u(i,j,1)
                                v(i,j,0)=v(i,j,1)
                                w(i,j,0)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,0)=rhov(isc,i,j,1)
                                end do
                                wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                            end do
                        end do
                    else if (iboun5==2) then

                    end if
                    write(*,*)'solid wall face 5'
                end if
            end if ! myid=0
            !................................................
            !     face 6 'avanti' (for periodicity see vel_up)
            !
            if (myid==nproc-1) then
                if (infout6==0) then  !inflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,jz+1)=up6(i,j)
                            v(i,j,jz+1)=vp6(i,j)
                            w(i,j,jz+1)=wp6(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)
                            end do
                            wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)

                        end do
                    end do
                    write(*,*)'inflow face 6'
                else if (infout6==1) then  !outflow
                    do j=1,jy
                        do i=1,jx
                            u(i,j,jz+1)=u(i,j,jz+1)
                            v(i,j,jz+1)=v(i,j,jz+1)
                            w(i,j,jz+1)=w(i,j,jz+1)
                            do isc=1,nscal
                                rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                            end do
                        end do
                    end do
                    write(*,*)'outflow face 6'
                else if (infout6==2) then  !wall
                    if (iboun6==0) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,jz+1)=0.
                                v(i,j,jz+1)=0.
                                w(i,j,jz+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                                end do
                                wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                            end do
                        end do
                    else if (iboun6==1) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,jz+1)=u(i,j,jz)
                                v(i,j,jz+1)=v(i,j,jz)
                                w(i,j,jz+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                                end do
                                wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                            end do
                        end do
                    else if (iboun6==2) then

                    end if
                    write(*,*)'solid wall face 6'
                end if
            end if !myid=nproc-1
        end do !end loop  kk=1,kp
        !-----------------------------------------------------------------------

        return
    end subroutine contour_se

    subroutine contour()
        !***********************************************************************
        ! set boundary condition for cartesian velocity and controvariant flux
        !
        use mysettings, only: bby,freesurface,i_rest,windyes
        use myarrays_metri3
        use myarrays_WB

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,ii,jj,kk,isc
        integer ierr
        integer kparastal,kparaendl
        !      integer windyes,ktime,i_rest

        real theta(jx)
        real upp,delta,yc

        real delmassa1,delmassa2,delmassa3
        real delmassa4,delmassa5,delmassa6

        real segno_parete,segno_pre
        real segno_u,segno_v,segno_w
        real u_int,v_int,w_int
        real ucp_old,wcp_old
        !
        !-----------------------------------------------------------------------
        if (i_rest==3) then
            if (infout1 /= 0)iboun1 = 2
            if (infout2 /= 0)iboun2 = 2
            if (infout5 /= 0)iboun5 = 2
            if (infout6 /= 0)iboun6 = 2
        end if
        !-----------------------------------------------------------------------
        ! (for periodicity see vel_up)
        !chicco portare right and left dentro un blocco do unico
        !     left side
        do ii=1,ip
            !
            if (infout1==0) then        !inflow
                do k=kparasta,kparaend
                    do j=1,jy
                        u(0,j,k)=up1(j,k)
                        v(0,j,k)=vp1(j,k)
                        w(0,j,k)=wp1(j,k)
                        do isc=1,nscal
                            rhov(isc,0,j,k)=rhovp1(isc,j,k)
                        end do
                        uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                    end do
                end do

            else if (infout1==1) then   !outflow
                do k=kparasta,kparaend
                    do j=1,jy

                        uc(0,j,k) = uc1_orl(j,k)

                        delmassa1= uc(0,j,k)-(du_dx1(j,k)*csx(0,j,k)+dv_dx1(j,k)*csy(0,j,k)+dw_dx1(j,k)*csz(0,j,k))

                        u(0,j,k)= du_dx1(j,k) + delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                        v(0,j,k)= dv_dx1(j,k) + delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                        w(0,j,k)= dw_dx1(j,k) + delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                        do isc=1,nscal
                            rhov(isc,0,j,k)=rhov(isc,1,j,k)
                        end do

                    end do
                end do

            else if (infout1==2) then      ! wall
                if (iboun1==0) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(0,j,k)=0.
                            v(0,j,k)=0.
                            w(0,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,0,j,k)=rhov(isc,1,j,k)
                            end do
                            uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                        end do
                    end do
                else if (iboun1==1) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(0,j,k)=0.
                            v(0,j,k)=v(1,j,k)
                            w(0,j,k)=w(1,j,k)
                            do isc=1,nscal
                                rhov(isc,0,j,k)=rhov(isc,1,j,k)
                            end do
                            uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                        end do
                    end do
                else if (iboun1==2) then
                    do k=kparasta,kparaend
                        do j=1,jy

                            if ( ucp1(j,k) .gt. 0.) then

                                ucp_old = csx(0,j,k)*up1(j,k) + csy(0,j,k)*vp1(j,k) + csz(0,j,k)*wp1(j,k)

                                delmassa1 = ucp1(j,k) - ucp_old

                                u(0,j,k)= up1(j,k) + delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                !             v(0,j,k)=v(1,j,k)
                                v(0,j,k)= vp1(j,k) + delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                w(0,j,k)= wp1(j,k) + delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                            else

                                ucp_old = csx(0,j,k)*u(1,j,k)+ csy(0,j,k)*v(1,j,k)+ csz(0,j,k)*w(1,j,k)

                                delmassa1 = ucp1(j,k) - ucp_old


                                u(0,j,k)= u(1,j,k) + delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                !             v(0,j,k)=v(1,j,k)
                                v(0,j,k)= v(1,j,k) + delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                w(0,j,k)= w(1,j,k) + delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                            end if


                            uc(0,j,k)=ucp1(j,k)

                            do isc=1,nscal
                                if (potenziale==0) then
                                    rhov(isc,0,j,k)=rhovp1(isc,j,k)
                                    if (ucp1(j,k).lt.0.) rhov(isc,0,j,k)=rhov(isc,1,j,k)

                                else
                                    rhov(isc,0,j,k)=rhovo1(isc,j,k)
                                end if
                            end do

                        end do
                    end do
                end if
            end if

            !     right side
            if (infout2==0) then        !inflow
                do k=kparasta,kparaend
                    do j=1,jy

                        u(jx+1,j,k)=up2(j,k)
                        v(jx+1,j,k)=vp2(j,k)
                        w(jx+1,j,k)=wp2(j,k)
                        do isc=1,nscal
                            rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
                        end do
                        uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                    end do
                end do
            else if (infout2==1) then   !outflow
                do k=kparasta,kparaend
                    do j=1,jy

                        uc(jx,j,k) = uc2_orl(j,k)

                        delmassa2=uc(jx,j,k)-(du_dx2(j,k)*csx(jx,j,k)+dv_dx2(j,k)*csy(jx,j,k)+dw_dx2(j,k)*csz(jx,j,k))

                        u(jx+1,j,k)= du_dx2(j,k) + delmassa2*csx(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
                        v(jx+1,j,k)= dv_dx2(j,k) + delmassa2*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
                        w(jx+1,j,k)= dw_dx2(j,k) + delmassa2*csz(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
                        do isc=1,nscal
                            rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                        end do


                    end do
                end do

            else if (infout2==2) then      !wall
                if (iboun2==0) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(jx+1,j,k)=0.
                            v(jx+1,j,k)=0.
                            w(jx+1,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                            end do
                            uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                        end do
                    end do
                else if (iboun2==1) then
                    do k=kparasta,kparaend
                        do j=1,jy
                            u(jx+1,j,k)=0.
                            v(jx+1,j,k)=v(jx,j,k)
                            w(jx+1,j,k)=w(jx,j,k)
                            do isc=1,nscal
                                rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                            end do
                            uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
                        end do
                    end do
                else if (iboun2==2) then
                    do k=kparasta,kparaend
                        do j=1,jy

                            if (ucp2(j,k) .lt. 0.) then

                                ucp_old = csx(jx,j,k)*up2(j,k)+ csy(jx,j,k)*vp2(j,k)+ csz(jx,j,k)*wp2(j,k)

                                delmassa2 = ucp2(j,k) - ucp_old

                                u(jx+1,j,k)= up2(j,k) + delmassa2*csx(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                                !             v(jx+1,j,k)=v(1,j,k)
                                v(jx+1,j,k)= vp2(j,k) + delmassa2*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                                w(jx+1,j,k)= wp2(j,k) + delmassa2*csz(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                                uc(jx,j,k)=ucp2(j,k)

                            else

                                ucp_old = csx(jx,j,k)*u(jx,j,k)+ csy(jx,j,k)*v(jx,j,k) + csz(jx,j,k)*w(jx,j,k)

                                delmassa2 = ucp2(j,k) - ucp_old

                                u(jx+1,j,k)= u(jx,j,k) + delmassa2*csx(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                                !             v(jx+1,j,k)=v(1,j,k)
                                v(jx+1,j,k)= v(jx,j,k) + delmassa2*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                                w(jx+1,j,k)= w(jx,j,k) + delmassa2*csz(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                                uc(jx,j,k)=ucp2(j,k)

                            end if

                            do isc=1,nscal
                                if (potenziale==0) then
                                    rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
                                    if ( ucp2(j,k) .gt. 0.) then
                                        rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                                    end if

                                else
                                    rhov(isc,jx+1,j,k)=rhovo2(isc,j,k)
                                end if
                            end do

                        end do
                    end do
                end if
            end if

        end do   !end loop  ii=1,ip
        !-----------------------------------------------------------------------
        !
        !     bottom side
        do jj=1,jp

            if (infout3==0) then             !inflow
                do k=kparasta,kparaend
                    do i=1,jx

                        u(i,0,k)=up3(i,k)
                        v(i,0,k)=vp3(i,k)
                        w(i,0,k)=wp3(i,k)
                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhovp3(isc,i,k)
                        end do
                        vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                    end do
                end do
            else if (infout3==1) then        !outflow
                do k=kparasta,kparaend
                    do i=1,jx

                        vc(i,0,k) = vc3_orl(i,k)

                        delmassa3= vc(i,0,k)-(du_dy3(i,k)*etx(i,0,k)+dv_dy3(i,k)*ety(i,0,k)+dw_dy3(i,k)*etz(i,0,k))

                        u(i,0,k)=  du_dy3(i,k)+delmassa3*etx(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)
                        v(i,0,k)=  dv_dy3(i,k)+delmassa3*ety(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)
                        w(i,0,k)=  dw_dy3(i,k)+delmassa3*etz(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)

                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhov(isc,i,1,k)
                        end do


                    end do
                end do

            else if (infout3==2) then           !wall
                if (iboun3==0) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=0.
                            v(i,0,k)=0.
                            w(i,0,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                else if (iboun3==1) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=u(i,1,k)
                            v(i,0,k)=0.
                            w(i,0,k)=w(i,1,k)
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                else if (iboun3==2) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,0,k)=0.
                            v(i,0,k)=0.
                            w(i,0,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,0,k)=rhov(isc,i,1,k)
                            end do
                            vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
                        end do
                    end do
                end if
            end if
            !
            !     upper side
            if (freesurface==0) then !no freesurface<<<<<<<<<<<<<<<<<<<<<<<<<<

                if (infout4==0) then             !inflow
                    do k=kparasta,kparaend
                        do i=1,jx

                            u(i,jy+1,k)=up4(i,k)
                            v(i,jy+1,k)=vp4(i,k)
                            w(i,jy+1,k)=wp4(i,k)
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhovp4(isc,i,k)
                            end do
                            vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                else if (infout4==1) then        !outflow
                    do k=kparasta,kparaend
                        do i=1,jx

                            vc(i,jy,k) = vc4_orl(i,k)

                            delmassa4= vc(i,jy,k)-(du_dy4(i,k)*etx(i,jy,k)+dv_dy4(i,k)*ety(i,jy,k)+dw_dy4(i,k)*etz(i,jy,k))

                            u(i,jy+1,k)= du_dy4(i,k)+delmassa4*etx(i,jy,k)/(etx(i,jy,k)**2+ety(i,jy,k)**2+etz(i,jy,k)**2)
                            v(i,jy+1,k)= dv_dy4(i,k)+delmassa4*ety(i,jy,k)/(etx(i,jy,k)**2+ety(i,jy,k)**2+etz(i,jy,k)**2)
                            w(i,jy+1,k)= dw_dy4(i,k)+delmassa4*etz(i,jy,k)/(etx(i,jy,k)**2+ety(i,jy,k)**2+etz(i,jy,k)**2)

                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do



                        end do
                    end do

                ! windyes=1 option, if at the upper side there is a wind stress

                else if (infout4==2 .and. windyes==0) then ! wall without wind
                    if (iboun4==0) then
                        do k=kparasta,kparaend
                            do i=1,jx
                                u(i,jy+1,k)=0.
                                v(i,jy+1,k)=0.
                                w(i,jy+1,k)=0.
                                do isc=1,nscal
                                    rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                                end do
                                vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k) +ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                            end do
                        end do
                    else if (iboun4==1) then
                        do k=kparasta,kparaend
                            do i=1,jx
                                u(i,jy+1,k)=u(i,jy,k)
                                v(i,jy+1,k)=0.
                                w(i,jy+1,k)=w(i,jy,k)
                                do isc=1,nscal
                                    rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                                end do
                                vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                                if ((i==4).and.(k==17)) then
                                    write(*,*)'VC_contour:',vc(i,32,k)
                                end if
                            end do
                        end do
                    else if (iboun4==2) then
                        do k=kparasta,kparaend
                            do i=1,jx
                                u(i,jy+1,k)=0.
                                v(i,jy+1,k)=0.
                                w(i,jy+1,k)=0.
                                do isc=1,nscal
                                    rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                                end do
                                vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                            end do
                        end do
                    end if

                else if (infout4==2 .and. windyes==1) then ! wall with wind extrapolation from the interior

                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k) = 1.5*u(i,jy,k)-.5*u(i,jy-1,k)
                            v(i,jy+1,k)= 0.0
                            w(i,jy+1,k) = 1.5*w(i,jy,k)-.5*w(i,jy-1,k)
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                            vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                end if

            else if (freesurface==1) then !free surface ON.<<<<<<<<<
                if (windyes==0) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k)=u(i,jy,k)
                            v(i,jy+1,k)=v(i,jy,k)
                            w(i,jy+1,k)=w(i,jy,k)
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                        !           vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)
                        !     >               +ety(i,jy,k)*v(i,jy+1,k)
                        !     >               +etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                else if (windyes==1) then
                    do k=kparasta,kparaend
                        do i=1,jx
                            u(i,jy+1,k) = 1.5*u(i,jy,k)-.5*u(i,jy-1,k)
                            v(i,jy+1,k)= v(i,jy,k)
                            w(i,jy+1,k) = 1.5*w(i,jy,k)-.5*w(i,jy-1,k)
                            do isc=1,nscal
                                rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                            end do
                        !           vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)
                        !     >               +ety(i,jy,k)*v(i,jy+1,k)
                        !     >               +etz(i,jy,k)*w(i,jy+1,k)
                        end do
                    end do
                end if
            end if !free surface<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        end do !end loop  jj=1,jp
        !-----------------------------------------------------------------------
        !
        !
        do kk=1,kp
            !
            !     back side
            if (myid==0) then
                if (infout5==0) then            !inflow
                    do j=1,jy
                        do i=1,jx

                            u(i,j,0)=up5(i,j)
                            v(i,j,0)=vp5(i,j)
                            w(i,j,0)=wp5(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,0)=rhovp5(isc,i,j)
                            end do
                            wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)

                        end do
                    end do
                else if (infout5==1) then            !outflow
                    do j=1,jy
                        do i=1,jx

                            wc(i,j,0) = wc5_orl(i,j)

                            delmassa5= wc(i,j,0)-(du_dz5(i,j)*ztx(i,j,0)+dv_dz5(i,j)*zty(i,j,0)+dw_dz5(i,j)*ztz(i,j,0))

                            u(i,j,0)= du_dz5(i,j)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
                            v(i,j,0)= dv_dz5(i,j)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
                            w(i,j,0)= dw_dz5(i,j)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                            do isc=1,nscal
                                rhov(isc,i,j,0)=rhov(isc,i,j,1)
                            end do


                        end do
                    end do


                else if (infout5==2) then            !wall
                    if (iboun5==0) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,0)=0.
                                v(i,j,0)=0.
                                w(i,j,0)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,0)=rhov(isc,i,j,1)
                                end do
                                wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                            end do
                        end do
                    else if (iboun5==1) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,0)=u(i,j,1)
                                v(i,j,0)=v(i,j,1)
                                w(i,j,0)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,0)=rhov(isc,i,j,1)
                                end do
                                wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                            end do
                        end do
                    else if (iboun5==2) then
                        do j=1,jy
                            do i=1,jx

                                if ( wcp5(i,j) .gt. 0.) then

                                    wcp_old =ztx(i,j,0)*up5(i,j)+zty(i,j,0)*vp5(i,j)+ztz(i,j,0)*wp5(i,j)

                                    delmassa5 = wcp5(i,j) - wcp_old

                                    u(i,j,0)= up5(i,j)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    !             v(i,j,0)=v(i,j,1)
                                    v(i,j,0)= vp5(i,j)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    w(i,j,0)= wp5(i,j)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    wc(i,j,0)=wcp5(i,j)

                                else

                                    wcp_old =ztx(i,j,0)*u(i,j,1)+zty(i,j,0)*v(i,j,1)+ztz(i,j,0)*w(i,j,1)

                                    delmassa5 = wcp5(i,j) - wcp_old

                                    u(i,j,0)= u(i,j,1)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    !             v(i,j,0)=v(i,j,1)
                                    v(i,j,0)= v(i,j,1)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    w(i,j,0)= w(i,j,1)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    wc(i,j,0)=wcp5(i,j)

                                end if

                                do isc=1,nscal
                                    if (potenziale==0) then
                                        rhov(isc,i,j,0)=rhovp5(isc,i,j)
                                        if (wcp5(i,j).lt.0.) then
                                            rhov(isc,i,j,0)=rhov(isc,i,j,1)
                                        end if

                                    else
                                        rhov(isc,i,j,0)=rhovo5(isc,i,j)
                                    end if
                                end do

                            end do
                        end do
                    end if
                end if
            end if
            !
            !     front side
            if (myid==nproc-1) then
                if (infout6==0) then            !inflow
                    do j=1,jy
                        do i=1,jx

                            u(i,j,jz+1)=up6(i,j)
                            v(i,j,jz+1)=vp6(i,j)
                            w(i,j,jz+1)=wp6(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)
                            end do
                            wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)

                        end do
                    end do
                else if (infout6==1) then            !outflow
                    do j=1,jy
                        do i=1,jx

                            wc(i,j,jz) = wc6_orl(i,j)

                            delmassa6= wc(i,j,jz)-(du_dz6(i,j)*ztx(i,j,jz) +dv_dz6(i,j)*zty(i,j,jz)+dw_dz6(i,j)*ztz(i,j,jz))

                            u(i,j,jz+1)= du_dz6(i,j)+delmassa6*ztx(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)
                            v(i,j,jz+1)= dv_dz6(i,j)+delmassa6*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)
                            w(i,j,jz+1)= dw_dz6(i,j)+delmassa6*ztz(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                            do isc=1,nscal
                                rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                            end do

                        end do
                    end do


                else if (infout6==2) then            !wall
                    if (iboun6==0) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,jz+1)=0.
                                v(i,j,jz+1)=0.
                                w(i,j,jz+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                                end do
                                wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                            end do
                        end do
                    else if (iboun6==1) then
                        do j=1,jy
                            do i=1,jx
                                u(i,j,jz+1)=u(i,j,jz)
                                v(i,j,jz+1)=v(i,j,jz)
                                w(i,j,jz+1)=0.



                                do isc=1,nscal
                                    rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                                end do
                                wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                            end do
                        end do
                    else if (iboun6==2) then
                        do j=1,jy
                            do i=1,jx

                                if (wcp6(i,j) .lt. 0.) then
                                    wcp_old =ztx(i,j,jz)*up6(i,j)+zty(i,j,jz)*vp6(i,j)+ztz(i,j,jz)*wp6(i,j)

                                    delmassa6 = wcp6(i,j) - wcp_old

                                    u(i,j,jz+1)= up6(i,j)+delmassa6*ztx(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                    !             v(i,j,jz+1)=v(i,j,jz)
                                    v(i,j,jz+1)= vp6(i,j)+delmassa6*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                    w(i,j,jz+1)= wp6(i,j)+delmassa6*ztz(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                    wc(i,j,jz)=wcp6(i,j)

                                else

                                    wcp_old =ztx(i,j,jz)*u(i,j,jz)+zty(i,j,jz)*v(i,j,jz)+ztz(i,j,jz)*w(i,j,jz)

                                    delmassa6 = wcp6(i,j) - wcp_old

                                    u(i,j,jz+1)= u(i,j,jz)+delmassa6*ztx(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                    !             v(i,j,jz+1)=v(i,j,jz)
                                    v(i,j,jz+1)= v(i,j,jz)+delmassa6*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                    w(i,j,jz+1)= w(i,j,jz)+delmassa6*ztz(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                    wc(i,j,jz)=wcp6(i,j)

                                end if

                                do isc=1,nscal
                                    if (potenziale==0) then
                                        rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)

                                        if (wcp6(i,j).gt.0.) then
                                            rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                                        end if
                                    else
                                        rhov(isc,i,j,jz+1)= rhovo6(isc,i,j)
                                    end if
                                end do

                            end do
                        end do
                    end if
                end if
            end if

        end do !end loop  kk=1,kp
        !-----------------------------------------------------------------------
        !     set boundary condition at the corner

        if (myid==0) then
            kparastal=0
            kparaendl=kparaend
        else if (myid==nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend+1
        else
            kparastal=kparasta
            kparaendl=kparaend
        end if

        ! FIRST: corner with two sides like orlansky (infout=1)

        !     corner close to side 5
        if (myid==0) then

            !     corner 5/1
            if (infout5==1 .and. infout1==1 ) then
                do j=0,jy+1
                    u(0,j,0)=u(1,j,1)
                    v(0,j,0)=v(1,j,1)
                    w(0,j,0)=w(1,j,1)
                    do isc=1,nscal
                        rhov(isc,0,j,0)=rhov(isc,1,j,1)
                    end do
                end do
            end if  !corner 5/1

            !     corner 5/2
            if (infout5==1 .and. infout2==1 ) then
                do j=0,jy+1
                    u(jx+1,j,0)=u(jx,j,1)
                    v(jx+1,j,0)=v(jx,j,1)
                    w(jx+1,j,0)=w(jx,j,1)
                    do isc=1,nscal
                        rhov(isc,jx+1,j,0)=rhov(isc,jx,j,1)
                    end do
                end do
            end if  !corner 5/2

            !     corner 5/3
            if (infout5==1 .and. infout3==1 ) then
                do i=0,jx+1
                    u(i,0,0)=u(i,1,1)
                    v(i,0,0)=v(i,1,1)
                    w(i,0,0)=w(i,1,1)
                    do isc=1,nscal
                        rhov(isc,i,0,0)=rhov(isc,i,1,1)
                    end do
                end do
            end if  !corner 5/3

            !     corner 5/4
            if (infout5==1 .and. infout4==1 ) then
                do i=0,jx+1
                    u(i,jy+1,0)=u(i,jy,1)
                    v(i,jy+1,0)=v(i,jy,1)
                    w(i,jy+1,0)=w(i,jy,1)
                    do isc=1,nscal
                        rhov(isc,i,jy+1,0)=rhov(isc,i,jy,1)
                    end do
                end do
            end if  !corner 5/4

        end if !proc 0
        !-----------------------------------------------------------------------
        !     corner close to side 6
        if (myid==nproc-1) then
            !     corner 6/1
            if (infout6==1 .and. infout1==1 ) then
                do j=0,jy+1
                    u(0,j,jz+1)=u(1,j,jz)
                    v(0,j,jz+1)=v(1,j,jz)
                    w(0,j,jz+1)=w(1,j,jz)
                    do isc=1,nscal
                        rhov(isc,0,j,jz+1)=rhov(isc,1,j,jz)
                    end do
                end do
            end if  !corner 6/1

            !     corner 6/2
            if (infout6==1 .and. infout2==1 ) then
                do j=0,jy+1
                    u(jx+1,j,jz+1)=u(jx,j,jz)
                    v(jx+1,j,jz+1)=v(jx,j,jz)
                    w(jx+1,j,jz+1)=w(jx,j,jz)
                    do isc=1,nscal
                        rhov(isc,jx+1,j,jz+1)=rhov(isc,jx,j,jz)
                    end do
                end do
            end if  !corner 6/2

            !     corner 6/3
            if (infout6==1 .and. infout3==1 ) then
                do i=0,jx+1
                    u(i,0,jz+1)=u(i,1,jz)
                    v(i,0,jz+1)=v(i,1,jz)
                    w(i,0,jz+1)=w(i,1,jz)
                    do isc=1,nscal
                        rhov(isc,i,0,jz+1)=rhov(isc,i,1,jz)
                    end do
                end do
            end if  !corner 6/3

            !     corner 6/4
            if (infout6==1 .and. infout4==1 ) then
                do i=0,jx+1
                    u(i,jy+1,jz+1)=u(i,jy,jz)
                    v(i,jy+1,jz+1)=v(i,jy,jz)
                    w(i,jy+1,jz+1)=w(i,jy,jz)
                    do isc=1,nscal
                        rhov(isc,i,jy+1,jz+1)=rhov(isc,i,jy,jz)
                    end do
                end do
            end if  !corner 6/4

        end if ! nproc-1

        !-----------------------------------------------------------------------
        ! corner 1/3 e 1/4

        !     corner 1/3
        if (infout1==1 .and. infout3==1 ) then
            do k=kparastal,kparaendl
                u(0,0,k)=u(1,1,k)
                v(0,0,k)=v(1,1,k)
                w(0,0,k)=w(1,1,k)
                do isc=1,nscal
                    rhov(isc,0,0,k)=rhov(isc,1,1,k)
                end do
            end do
        end if  !corner 1/3

        !     corner 1/4
        if (infout1==1 .and. infout4==1 ) then
            do k=kparastal,kparaendl
                u(0,jy+1,k)=u(1,jy,k)
                v(0,jy+1,k)=v(1,jy,k)
                w(0,jy+1,k)=w(1,jy,k)
                do isc=1,nscal
                    rhov(isc,0,jy+1,k)=rhov(isc,1,jy,k)
                end do
            end do
        end if  !corner 1/4

        !-----------------------------------------------------------------------
        ! corner 2/3 e 2/4
        !     corner 2/3
        if (infout2==1 .and. infout3==1 ) then
            do k=kparastal,kparaendl
                u(jx+1,0,k)=u(jx,1,k)
                v(jx+1,0,k)=v(jx,1,k)
                w(jx+1,0,k)=w(jx,1,k)
                do isc=1,nscal
                    rhov(isc,jx+1,0,k)=rhov(isc,jx,1,k)
                end do
            end do
        end if  !corner 2/3

        !     corner 2/4
        if (infout2==1 .and. infout4==1 ) then
            do k=kparastal,kparaendl
                u(jx+1,jy+1,k)=u(jx,jy,k)
                v(jx+1,jy+1,k)=v(jx,jy,k)
                w(jx+1,jy+1,k)=w(jx,jy,k)
                do isc=1,nscal
                    rhov(isc,jx+1,jy+1,k)=rhov(isc,jx,jy,k)
                end do
            end do
        end if  !corner 2/4
        !-----------------------------------------------------------------------
        ! SECOND: case sides with inlow (infout=0)

        !     side 1
        if (infout1==0) then
            i=0
            do k=kparastal,kparaendl
                !     corner 1/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do

                !     corner 1/4
                u(i,jy+1,k)=u(i,jy,k)
                v(i,jy+1,k)=v(i,jy,k)
                w(i,jy+1,k)=w(i,jy,k)
                do isc=1,nscal
                    rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                end do
            end do

            do j=0,jy+1
                !     corner 1/5
                if (myid==0) then
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !     corner 1/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 1

        !-----------------------------------------------------------------------
        !     side 2
        if (infout2==0) then
            i=jx+1
            do k=kparastal,kparaendl
                !     corner 2/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do

                !     corner 2/4
                u(i,jy+1,k)=u(i,jy,k)
                v(i,jy+1,k)=v(i,jy,k)
                w(i,jy+1,k)=w(i,jy,k)
                do isc=1,nscal
                    rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                end do
            end do

            do j=0,jy+1
                if (myid==0) then
                    !     corner 2/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !     corner 2/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 2

        !-----------------------------------------------------------------------
        !     side 3
        if (infout3==0) then
            j=0
            do k=kparastal,kparaendl
                !     corner 3/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                !     corner 3/2
                u(jx+1,j,k)=u(jx,j,k)
                v(jx+1,j,k)=v(jx,j,k)
                w(jx+1,j,k)=w(jx,j,k)
                do isc=1,nscal
                    rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                end do
            end do

            do i=0,jx+1
                if (myid==0) then
                    !     corner 3/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !     corner 3/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 3

        !-----------------------------------------------------------------------
        !     side 4
        if (infout4==0) then
            j=jy+1
            do k=kparastal,kparaendl
                !     corner 4/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                !     corner 4/2
                u(jx+1,j,k)=u(jx,j,k)
                v(jx+1,j,k)=v(jx,j,k)
                w(jx+1,j,k)=w(jx,j,k)
                do isc=1,nscal
                    rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                end do
            end do

            do i=0,jx+1
                if (myid==0) then
                    !     corner 4/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !     corner 4/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 4

        !-----------------------------------------------------------------------
        !     side 5
        if (infout5==0) then
            if (myid==0) then
                k=0
                do j=0,jy+1
                    !     corner 5/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do

                    !     corner 5/2
                    u(jx+1,j,k)=u(jx,j,k)
                    v(jx+1,j,k)=v(jx,j,k)
                    w(jx+1,j,k)=w(jx,j,k)
                    do isc=1,nscal
                        rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                    end do
                end do

                do i=0,jx+1
                    !     corner 5/3
                    u(i,0,k)=u(i,1,k)
                    v(i,0,k)=v(i,1,k)
                    w(i,0,k)=w(i,1,k)
                    do isc=1,nscal
                        rhov(isc,i,0,k)=rhov(isc,i,1,k)
                    end do
                    !     corner 5/4
                    u(i,jy+1,k)=u(i,jy,k)
                    v(i,jy+1,k)=v(i,jy,k)
                    w(i,jy+1,k)=w(i,jy,k)
                    do isc=1,nscal
                        rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                    end do
                end do
            end if
        end if !side 5

        !-----------------------------------------------------------------------
        !     side 6
        if (infout6==0) then
            if (myid==nproc-1) then
                k=jz+1
                do j=0,jy+1
                    !     corner 6/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do
                    !     corner 6/2
                    u(jx+1,j,k)=u(jx,j,k)
                    v(jx+1,j,k)=v(jx,j,k)
                    w(jx+1,j,k)=w(jx,j,k)
                    do isc=1,nscal
                        rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                    end do
                end do

                do i=0,jx+1
                    if (myid==0) then
                        !     corner 6/3
                        u(i,0,k)=u(i,1,k)
                        v(i,0,k)=v(i,1,k)
                        w(i,0,k)=w(i,1,k)
                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhov(isc,i,1,k)
                        end do
                    else if (myid==nproc-1) then
                        !     corner 6/4
                        u(i,jy+1,k)=u(i,jy,k)
                        v(i,jy+1,k)=v(i,jy,k)
                        w(i,jy+1,k)=w(i,jy,k)
                        do isc=1,nscal
                            rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                        end do
                    end if
                end do
            end if
        end if !side 6
        !-----------------------------------------------------------------------
        ! THIRD: case wall side (infout=2)
        !
        !     side 1
        if (infout1==2) then
            i=0
            do k=kparastal,kparaendl
                !     corner 1/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do
                !     corner 1/4
                u(i,jy+1,k)=u(i,jy,k)
                v(i,jy+1,k)=v(i,jy,k)
                w(i,jy+1,k)=w(i,jy,k)
                do isc=1,nscal
                    rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                end do
            end do

            do j=0,jy+1
                if (myid==0) then
                    !       corner 1/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !       corner 1/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 1

        !-----------------------------------------------------------------------
        !     side 2
        if (infout2==2) then
            i=jx+1
            do k=kparastal,kparaendl
                !       corner 2/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do

                !       corner 2/4
                u(i,jy+1,k)=u(i,jy,k)
                v(i,jy+1,k)=v(i,jy,k)
                w(i,jy+1,k)=w(i,jy,k)
                do isc=1,nscal
                    rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                end do
            end do

            do j=0,jy+1
                if (myid==0) then
                    !       corner 2/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !       corner 2/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 2

        !-----------------------------------------------------------------------
        !     side 5
        if (infout5==2) then
            if (myid==0) then
                k=0
                do j=0,jy+1
                    !       corner 5/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do

                    !       corner 5/2
                    u(jx+1,j,k)=u(jx,j,k)
                    v(jx+1,j,k)=v(jx,j,k)
                    w(jx+1,j,k)=w(jx,j,k)
                    do isc=1,nscal
                        rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                    end do
                end do

                do i=0,jx+1
                    !       corner 5/3
                    u(i,0,k)=u(i,1,k)
                    v(i,0,k)=v(i,1,k)
                    w(i,0,k)=w(i,1,k)
                    do isc=1,nscal
                        rhov(isc,i,0,k)=rhov(isc,i,1,k)
                    end do

                    !       corner 5/4
                    u(i,jy+1,k)=u(i,jy,k)
                    v(i,jy+1,k)=v(i,jy,k)
                    w(i,jy+1,k)=w(i,jy,k)
                    do isc=1,nscal
                        rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                    end do
                end do
            end if
        end if !side 5

        !-----------------------------------------------------------------------
        !     side 6
        if (infout6==2) then
            if (myid==nproc-1) then
                k=jz+1
                do j=0,jy+1
                    !       corner 6/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do

                    !       corner 6/2
                    u(jx+1,j,k)=u(jx,j,k)
                    v(jx+1,j,k)=v(jx,j,k)
                    w(jx+1,j,k)=w(jx,j,k)
                    do isc=1,nscal
                        rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                    end do
                end do
                do i=0,jx+1
                    !       corner 6/3
                    u(i,0,k)=u(i,1,k)
                    v(i,0,k)=v(i,1,k)
                    w(i,0,k)=w(i,1,k)
                    do isc=1,nscal
                        rhov(isc,i,0,k)=rhov(isc,i,1,k)
                    end do
                    !       corner 6/4
                    u(i,jy+1,k)=u(i,jy,k)
                    v(i,jy+1,k)=v(i,jy,k)
                    w(i,jy+1,k)=w(i,jy,k)
                    do isc=1,nscal
                        rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                    end do
                end do
            end if
        end if !side 6


        !-----------------------------------------------------------------------
        !     side 3
        if (infout3==2) then
            j=0
            do k=kparastal,kparaendl
                !       corner 3/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                !       corner 3/2
                u(jx+1,j,k)=u(jx,j,k)
                v(jx+1,j,k)=v(jx,j,k)
                w(jx+1,j,k)=w(jx,j,k)
                do isc=1,nscal
                    rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                end do
            end do

            do i=0,jx+1
                if (myid==0) then
                    !       corner 3/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !       corner 3/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 3

        !-----------------------------------------------------------------------
        !     side 4
        if (infout4==2) then
            j=jy+1
            do k=kparastal,kparaendl
                !       corner 4/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                !       corner 4/2
                u(jx+1,j,k)=u(jx,j,k)
                v(jx+1,j,k)=v(jx,j,k)
                w(jx+1,j,k)=w(jx,j,k)
                do isc=1,nscal
                    rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                end do
            end do

            do i=0,jx+1
                if (myid==0) then
                    !       corner 4/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    !       corner 4/6
                    u(i,j,jz+1)=u(i,j,jz)
                    v(i,j,jz+1)=v(i,j,jz)
                    w(i,j,jz+1)=w(i,j,jz)
                    do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                    end do
                end if
            end do

        end if !side 4


        return
    end subroutine contour

    subroutine contourp()
        !***********************************************************************
        ! compute cartesian velocity and controvariant in periodic cell at
        ! step n+1, at the corner computation at the end of the sub
        !
        use mysettings, only: bby,freesurface,i_rest,windyes
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,ii,isc,kk
        !
        integer ierr,status(MPI_STATUS_SIZE)
        ! integer myid,nproc,
        !integer kparasta,kparaend
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
                end do
            end do

        end do
        !
        ! periodic cell in eta (also outside the boundary)
        !
        do j=1,1-jp

            if (myid==0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid==nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend+1-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            end if

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
                end do
            end do

        end do
        !
        ! periodic cell in zita (also outside the boundary)
        !
        do k=1,1-kp

            call MPI_TYPE_VECTOR(jy+2,jx,jx+2,MPI_REAL_SD,plantype,ierr)
            call MPI_TYPE_COMMIT(plantype,ierr)

            if (myid==0) then
                call MPI_SENDRECV(u(1,0,1),1,plantype,nproc-1,12,u(1,0,0),1,plantype,nproc-1,11,MPI_COMM_WORLD,status,ierr)

                call MPI_SENDRECV(v(1,0,1),1,plantype,nproc-1,14,v(1,0,0),1,plantype,nproc-1,13,MPI_COMM_WORLD,status,ierr)

                call MPI_SENDRECV(w(1,0,1),1,plantype,nproc-1,16,w(1,0,0),1,plantype,nproc-1,15,MPI_COMM_WORLD,status,ierr)
            end if

            if (myid==nproc-1) then

                call MPI_SENDRECV(u(1,0,jz),1,plantype,0,11,u(1,0,jz+1),1,plantype,0,12,MPI_COMM_WORLD,status,ierr)

                call MPI_SENDRECV(v(1,0,jz),1,plantype,0,13,v(1,0,jz+1),1,plantype,0,14,MPI_COMM_WORLD,status,ierr)

                call MPI_SENDRECV(w(1,0,jz),1,plantype,0,15,w(1,0,jz+1),1,plantype,0,16,MPI_COMM_WORLD,status,ierr)

            end if

            allocate(rho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))

            do isc=1,nscal
                do kk=kparasta-1,kparaend+1 !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            rho(i,j,kk)=rhov(isc,i,j,kk)
                        end do
                    end do
                end do

                if (myid==0) then
                    call MPI_SENDRECV(rho(1,0,1),1,plantype,nproc-1,18+isc,rho(1,0,0),&
                        1,plantype,nproc-1,17+isc,MPI_COMM_WORLD,status,ierr)
                end if
                if (myid==nproc-1) then
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

        end do
        !
        ! To calculate the next pressure at surface for the next iteration
        if (myid==0) then
            kparastal=0
            kparaendl=kparaend
        else if (myid==nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend+1
        else
            kparastal=kparasta
            kparaendl=kparaend
        end if

        if (freesurface==1) then !free surface ON.<<<<<<<<<
            if (myid==0) then
                write(*,*)'Free Surface ON'
            end if
            do k=kparastal,kparaendl
                do i=0,jx+1
                    next_prs(i,k)=((fi(i,jy+1,k)+fi(i,jy,k))*0.5) - (v(i,jy+1,k) * dt * bby)
                end do
            end do
        else
            if (myid==0) then
                write(*,*)'Free Surface OFF'
            end if
        end if !free surface<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        return
    end subroutine contourp

    subroutine initialize_contour()

            ! face 1
        allocate(up1(0:jy+1,0:jz+1))
        allocate(vp1(0:jy+1,0:jz+1))
        allocate(wp1(0:jy+1,0:jz+1))
        allocate(rhovp1(nscal,0:jy+1,0:jz+1))
        allocate(tke1(0:jy+1,kparasta-1:kparaend+1)) !face 1
        tke1 = 0.
        up1=0.
        vp1=0.
        wp1=0.
        rhovp1=0.

        ! face 2
        allocate(up2(0:jy+1,0:jz+1))
        allocate(vp2(0:jy+1,0:jz+1))
        allocate(wp2(0:jy+1,0:jz+1))
        allocate(rhovp2(nscal,0:jy+1,0:jz+1))
        allocate(tke2(0:jy+1,kparasta-1:kparaend+1)) !face 2
        tke2 = 0.
        up2=0.
        vp2=0.
        wp2=0.
        rhovp2=0.

        ! face 3
        allocate(up3(0:jx+1,0:jz+1))
        allocate(vp3(0:jx+1,0:jz+1))
        allocate(wp3(0:jx+1,0:jz+1))
        allocate(rhovp3(nscal,0:jx+1,0:jz+1))
        up3   = 0.
        vp3   = 0.
        wp3   = 0.
        rhovp3 = 0.

        ! face 4
        allocate(up4(0:jx+1,0:jz+1))
        allocate(vp4(0:jx+1,0:jz+1))
        allocate(wp4(0:jx+1,0:jz+1))
        allocate(rhovp4(nscal,0:jx+1,0:jz+1))
        up4   = 0.
        vp4   = 0.
        wp4   = 0.
        rhovp4 = 0.

        ! face 5
        allocate(up5(0:jx+1,0:jy+1))
        allocate(vp5(0:jx+1,0:jy+1))
        allocate(wp5(0:jx+1,0:jy+1))
        allocate(rhovp5(nscal,0:jx+1,0:jy+1))
        allocate(tke5(0:jx+1,0:jy+1))
        tke5 = 0.
        up5=0.
        vp5=0.
        wp5=0.
        rhovp5=0.

        ! face 6
        allocate(up6(0:jx+1,0:jy+1))
        allocate(vp6(0:jx+1,0:jy+1))
        allocate(wp6(0:jx+1,0:jy+1))
        allocate(rhovp6(nscal,0:jx+1,0:jy+1))
        allocate(tke6(0:jx+1,0:jy+1))
        tke6 = 0.
        up6=0.
        vp6=0.
        wp6=0.
        rhovp6=0.

    end subroutine initialize_contour

    subroutine condi()
        !************************************************************************
        ! set boundary conditions on du dv dw
        ! at faces 1 and 2 "sinistra" and "destra"
        ! as in Kim and Moin
        !
        implicit none
        !-----------------------------------------------------------------------
        ! array declaration
        integer i,j,k,ierr
        !-----------------------------------------------------------------------
        !
        do k=kparasta,kparaend
            do j=1,jy
                !     face 1 sinistra
                delu(0,j,k)= &
                    1.875*gra1(1,j,k)-1.25*gra1(2,j,k)+.375*gra1(3,j,k)
                delu(0,j,k)=-delu(0,j,k)
                !
                delv(0,j,k)= &
                    1.875*gra2(1,j,k)-1.25*gra2(2,j,k)+.375*gra2(3,j,k)
                delv(0,j,k)=-delv(0,j,k)
                !
                delw(0,j,k)= &
                    1.875*gra3(1,j,k)-1.25*gra3(2,j,k)+.375*gra3(3,j,k)
                delw(0,j,k)=-delw(0,j,k)
                !
                !     face 2 destra
                delu(jx+1,j,k)= &
                    .375*gra1(jx-2,j,k)-1.25*gra1(jx-1,j,k)+1.875*gra1(jx,j,k)
                delu(jx+1,j,k)=-delu(jx+1,j,k)
                !
                delv(jx+1,j,k)= &
                    .375*gra2(jx-2,j,k)-1.25*gra2(jx-1,j,k)+1.875*gra2(jx,j,k)
                delv(jx+1,j,k)=-delv(jx+1,j,k)
                !
                delw(jx+1,j,k)= &
                    .375*gra3(jx-2,j,k)-1.25*gra3(jx-1,j,k)+1.875*gra3(jx,j,k)
                delw(jx+1,j,k)=-delw(jx+1,j,k)
            !
            end do
        end do
        !
        !***********************************************************************
        ! set boundary conditions on du dv dw
        ! at faces 3 and 4 "sotto" and "sopra"
        ! as in Kim and Moin
        !
        !-----------------------------------------------------------------------
        do k=kparasta,kparaend
            do i=1,jx
                !     face 3 sotto
                delu(i,0,k)= &
                    1.875*gra1(i,1,k)-1.25*gra1(i,2,k)+.375*gra1(i,3,k)
                delu(i,0,k)=-delu(i,0,k)
                !
                delv(i,0,k)= &
                    1.875*gra2(i,1,k)-1.25*gra2(i,2,k)+.375*gra2(i,3,k)
                delv(i,0,k)=-delv(i,0,k)
                !
                delw(i,0,k)= &
                    1.875*gra3(i,1,k)-1.25*gra3(i,2,k)+.375*gra3(i,3,k)
                delw(i,0,k)=-delw(i,0,k)
                !
                !     face 4 sopra
                delu(i,jy+1,k)= &
                    .375*gra1(i,jy-2,k)-1.25*gra1(i,jy-1,k)+1.875*gra1(i,jy,k)
                delu(i,jy+1,k)=-delu(i,jy+1,k)
                !
                delv(i,jy+1,k)= &
                    .375*gra2(i,jy-2,k)-1.25*gra2(i,jy-1,k)+1.875*gra2(i,jy,k)
                delv(i,jy+1,k)=-delv(i,jy+1,k)
                !
                delw(i,jy+1,k)= &
                    .375*gra3(i,jy-2,k)-1.25*gra3(i,jy-1,k)+1.875*gra3(i,jy,k)
                delw(i,jy+1,k)=-delw(i,jy+1,k)
            !
            end do
        end do
        !
        !***********************************************************************
        ! set boundary conditions on du dv dw
        ! at faces 5 and 6 "avanti" and "indietro"
        ! as in Kim and Moin
        !
        !-----------------------------------------------------------------------
        !
        !     face 5 avanti
        if (myid.eq.0) then
            !
            do j=1,jy
                do i=1,jx
                    !
                    delu(i,j,0)= &
                        1.875*gra1(i,j,1)-1.25*gra1(i,j,2)+.375*gra1(i,j,3)
                    delu(i,j,0)=-delu(i,j,0)
                    !
                    delv(i,j,0)= &
                        1.875*gra2(i,j,1)-1.25*gra2(i,j,2)+.375*gra2(i,j,3)
                    delv(i,j,0)=-delv(i,j,0)
                    !
                    delw(i,j,0)= &
                        1.875*gra3(i,j,1)-1.25*gra3(i,j,2)+.375*gra3(i,j,3)
                    delw(i,j,0)=-delw(i,j,0)
                !
                end do
            end do

        end if
        !
        !     face 6 indietro
        if (myid.eq.nproc-1) then
            !
            do j=1,jy
                do i=1,jx
                    !
                    delu(i,j,jz+1)= &
                        .375*gra1(i,j,jz-2)-1.25*gra1(i,j,jz-1)+1.875*gra1(i,j,jz)
                    delu(i,j,jz+1)=-delu(i,j,jz+1)
                    !
                    delv(i,j,jz+1)= &
                        .375*gra2(i,j,jz-2)-1.25*gra2(i,j,jz-1)+1.875*gra2(i,j,jz)
                    delv(i,j,jz+1)=-delv(i,j,jz+1)
                    !
                    delw(i,j,jz+1)= &
                        .375*gra3(i,j,jz-2)-1.25*gra3(i,j,jz-1)+1.875*gra3(i,j,jz)
                    delw(i,j,jz+1)=-delw(i,j,jz+1)
                !
                end do
            end do
        !
        end if

        return

    end subroutine condi

end module contour_module
