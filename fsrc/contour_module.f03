module contour_module

    use myarrays_velo3
    use orlansky_module
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    integer(kind=c_int),bind(C) :: iboun1,iboun2
    integer(kind=c_int),bind(C) :: iboun3,iboun4
    integer(kind=c_int),bind(C) :: iboun5,iboun6


    logical(kind=c_bool),bind(C) :: ibodybuffer1,ibodybuffer2
    logical(kind=c_bool),bind(C) :: ibodybuffer3,ibodybuffer4
    logical(kind=c_bool),bind(C) :: ibodybuffer5,ibodybuffer6

    ! arrays for inflow
    real,allocatable :: up1(:,:),vp1(:,:),wp1(:,:),rhovp1(:,:,:)
    real,allocatable :: up2(:,:),vp2(:,:),wp2(:,:),rhovp2(:,:,:)
    real,allocatable :: up3(:,:),vp3(:,:),wp3(:,:),rhovp3(:,:,:)
    real,allocatable :: up4(:,:),vp4(:,:),wp4(:,:),rhovp4(:,:,:)
    real,allocatable :: up5(:,:),vp5(:,:),wp5(:,:),rhovp5(:,:,:)
    real,allocatable :: up6(:,:),vp6(:,:),wp6(:,:),rhovp6(:,:,:)

    real,allocatable,dimension(:,:,:) :: rhovo1,rhovo2,rhovo5,rhovo6
    real,allocatable,dimension(:,:,:) :: rhovn1,rhovn2,rhovn5,rhovn6

    real,allocatable,dimension(:,:) :: ucp1
    real,allocatable,dimension(:,:) :: ucp2
    real,allocatable,dimension(:,:) :: wcp5
    real,allocatable,dimension(:,:) :: wcp6

    real,allocatable :: tke1(:,:),tke2(:,:)
    real,allocatable :: tke5(:,:),tke6(:,:)
    real,allocatable :: tkepom1(:,:,:),tkepom2(:,:,:)
    real,allocatable :: tkepom5(:,:,:),tkepom6(:,:,:)

contains

    subroutine initialize_contour()

        ! face 1
        allocate(up1(0:n2+1,0:n3+1))
        allocate(vp1(0:n2+1,0:n3+1))
        allocate(wp1(0:n2+1,0:n3+1))
        allocate(rhovp1(nscal,0:n2+1,0:n3+1))
        allocate(tke1(0:n2+1,kparasta-1:kparaend+1)) !face 1
        tke1(:,:)=0.0
        up1(:,:)=0.0
        vp1(:,:)=0.0
        wp1(:,:)=0.0
        rhovp1(:,:,:)=0.0

        ! face 2
        allocate(up2(0:n2+1,0:n3+1))
        allocate(vp2(0:n2+1,0:n3+1))
        allocate(wp2(0:n2+1,0:n3+1))
        allocate(rhovp2(nscal,0:n2+1,0:n3+1))
        allocate(tke2(0:n2+1,kparasta-1:kparaend+1)) !face 2
        tke2(:,:)=0.0
        up2(:,:)=0.0
        vp2(:,:)=0.0
        wp2(:,:)=0.0
        rhovp2(:,:,:)=0.0

        ! face 3
        allocate(up3(0:n1+1,0:n3+1))
        allocate(vp3(0:n1+1,0:n3+1))
        allocate(wp3(0:n1+1,0:n3+1))
        allocate(rhovp3(nscal,0:n1+1,0:n3+1))
        up3(:,:)=0.0
        vp3(:,:)=0.0
        wp3(:,:)=0.0
        rhovp3(:,:,:)=0.0

        ! face 4
        allocate(up4(0:n1+1,0:n3+1))
        allocate(vp4(0:n1+1,0:n3+1))
        allocate(wp4(0:n1+1,0:n3+1))
        allocate(rhovp4(nscal,0:n1+1,0:n3+1))
        up4(:,:)=0.0
        vp4(:,:)=0.0
        wp4(:,:)=0.0
        rhovp4(:,:,:)=0.0

        ! face 5
        allocate(up5(0:n1+1,0:n2+1))
        allocate(vp5(0:n1+1,0:n2+1))
        allocate(wp5(0:n1+1,0:n2+1))
        allocate(rhovp5(nscal,0:n1+1,0:n2+1))
        allocate(tke5(0:n1+1,0:n2+1))
        tke5(:,:)=0.0
        up5(:,:)=0.0
        vp5(:,:)=0.0
        wp5(:,:)=0.0
        rhovp5(:,:,:)=0.0

        ! face 6
        allocate(up6(0:n1+1,0:n2+1))
        allocate(vp6(0:n1+1,0:n2+1))
        allocate(wp6(0:n1+1,0:n2+1))
        allocate(rhovp6(nscal,0:n1+1,0:n2+1))
        allocate(tke6(0:n1+1,0:n2+1))
        tke6(:,:)=0.0
        up6(:,:)=0.0
        vp6(:,:)=0.0
        wp6(:,:)=0.0
        rhovp6(:,:,:)=0.0

    end subroutine initialize_contour

    subroutine contourp_se()

        ! compute cartesian velocity and controvarian fluxes in periodic cells at time n+1
        ! at the corner the computation is made at the end
        !
        !-----------------------------------------------------------------------
        integer :: i,j,k,isc,err
        integer :: status(MPI_STATUS_SIZE),ierr,siz,MPI_UVW_TYPE
        real,dimension(1:n1,0:n2+1)  :: sendbuf,recvbuf
        !-----------------------------------------------------------------------

        siz=n1*(n2+2)

        ! periodic cells in csi (also out of the domain)
        if (ip==0) then
            do k=kparasta,kparaend
                do j=0,n2+1
                    ! face 1
                    u(0,j,k)=u(n1,j,k)
                    v(0,j,k)=v(n1,j,k)
                    w(0,j,k)=w(n1,j,k)
                    ! face 2
                    u(n1+1,j,k)=u(1,j,k)
                    v(n1+1,j,k)=v(1,j,k)
                    w(n1+1,j,k)=w(1,j,k)
                    !
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,n1,j,k)
                        rhov(isc,n1+1,j,k)=rhov(isc,1,j,k)
                    end do
                !
                end do
            end do
        end if


        ! periodic cells in zita (also out of the domain)
        if (kp==0) then

            call MPI_TYPE_VECTOR(n2+2,n1,n1,MPI_REAL_SD,MPI_UVW_TYPE,ierr)
            call MPI_TYPE_COMMIT(MPI_UVW_TYPE,ierr)
            if (myid==0) then

                !------------------------------------------------------------------------
                ! SEND U
                do i=1,n1
                    do j=0,n2+1
                        sendbuf(i,j)=u(i,j,0)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,193,MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND V
                do i=1,n1
                    do j=0,n2+1
                        sendbuf(i,j)=v(i,j,0)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,194,MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND W
                do i=1,n1
                    do j=0,n2+1
                        sendbuf(i,j)=w(i,j,0)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,195,MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND RHOV
                do isc=1,nscal
                    do i=1,n1
                        do j=0,n2+1
                            sendbuf(i,j)=rhov(isc,i,j,0)
                        end do
                    end do
                    call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,nproc-1,196+isc,MPI_COMM_WORLD,ierr)
                end do

                !------------------------------------------------------------------------
                ! RECV U
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,nproc-1,182,MPI_COMM_WORLD,status,ierr)
                do i=1,n1
                    do j=0,n2+1
                        u(i,j,0)=recvbuf(i,j)
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV V
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,nproc-1,183,MPI_COMM_WORLD,status,ierr)
                do i=1,n1
                    do j=0,n2+1
                        v(i,j,0)=recvbuf(i,j)
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV W
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,nproc-1,184,MPI_COMM_WORLD,status,ierr)
                do i=1,n1
                    do j=0,n2+1
                        w(i,j,0)=recvbuf(i,j)
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV RHOV
                do isc=1,nscal
                    call MPI_RECV(recvbuf(1,0),siz,MPI_REAL_SD,nproc-1,185+isc,MPI_COMM_WORLD,status,ierr)
                    do i=1,n1
                        do j=0,n2+1
                            rhov(isc,i,j,0)=recvbuf(i,j)
                        end do
                    end do
                end do

            !--------------------------------------------------------------------------------------
            else if (myid==nproc-1) then
                !--------------------------------------------------------------------------------------

                !------------------------------------------------------------------------
                ! RECV U
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,193,MPI_COMM_WORLD,status,ierr)
                do i=1,n1
                    do j=0,n2+1
                        u(i,j,n3+1)=recvbuf(i,j)
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV V
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,194,MPI_COMM_WORLD,status,ierr)
                do i=1,n1
                    do j=0,n2+1
                        v(i,j,n3+1)=recvbuf(i,j)
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV W
                call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,195,MPI_COMM_WORLD,status,ierr)
                do i=1,n1
                    do j=0,n2+1
                        w(i,j,n3+1)=recvbuf(i,j)
                    end do
                end do
                !------------------------------------------------------------------------
                ! RECV  RHOV
                do isc=1,nscal
                    call MPI_RECV(recvbuf(1,0),1,MPI_UVW_TYPE,0,196+isc,MPI_COMM_WORLD,status,ierr)
                    do i=1,n1
                        do j=0,n2+1
                            rhov(isc,i,j,n3+1)=recvbuf(i,j)
                        end do
                    end do
                end do

                !------------------------------------------------------------------------
                ! SEND U
                do i=1,n1
                    do j=0,n2+1
                        sendbuf(i,j)=u(i,j,n3+1)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,182,MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND V
                do i=1,n1
                    do j=0,n2+1
                        sendbuf(i,j)=v(i,j,n3+1)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,183,MPI_COMM_WORLD,ierr)
                !------------------------------------------------------------------------
                ! SEND W
                do i=1,n1
                    do j=0,n2+1
                        sendbuf(i,j)=w(i,j,n3+1)
                    end do
                end do
                call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,184,MPI_COMM_WORLD,ierr)


                !------------------------------------------------------------------------
                ! SEND RHOV
                do isc=1,nscal
                    !if (.not.allocated(rhosendbuf)) then
                    ! call MPI_ABORT(MPI_COMM_WORLD,err,ierr)
                    ! stop "CARE CHIMBA ERROR"
                    !end if
                    do i=1,n1
                        do j=0,n2+1
                            sendbuf(i,j)=rhov(isc,i,j,n3+1)
                        end do
                    end do

                    call MPI_SEND(sendbuf(1,0),1,MPI_UVW_TYPE,0,185+isc,MPI_COMM_WORLD,ierr)

                end do


            end if
            call MPI_TYPE_FREE(MPI_UVW_TYPE,ierr)
        end if

        ! at the corner (temporary version)
        !
        if (ip==0) then
            do j=0,n2
                if (myid==0) then
                    u(0,j,0)=u(n1,j,0)
                    u(n1+1,j,0)=u(1,j,0)

                    v(0,j,0)=v(n1,j,0)
                    v(n1+1,j,0)=v(1,j,0)

                    w(0,j,0)=w(n1,j,0)
                    w(n1+1,j,0)=w(1,j,0)

                    do isc=1,nscal
                        rhov(isc,0,j,0)=rhov(isc,n1,j,0)
                        rhov(isc,n1+1,j,0)=rhov(isc,1,j,0)
                    end do
                end if

                if (myid==nproc-1) then
                    u(0,j,n3+1)=u(n1,j,n3+1)
                    u(n1+1,j,n3+1)=u(1,j,n3+1)

                    v(0,j,n3+1)=v(n1,j,n3+1)
                    v(n1+1,j,n3+1)=v(1,j,n3+1)

                    w(0,j,n3+1)=w(n1,j,n3+1)
                    w(n1+1,j,n3+1)=w(1,j,n3+1)

                    do isc=1,nscal
                        rhov(isc,0,j,n3+1)=rhov(isc,n1,j,n3+1)
                        rhov(isc,n1+1,j,n3+1)=rhov(isc,1,j,n3+1)
                    end do
                end if
            end do
        end if

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
        ! face 1 'sinistra' (for periodicity see vel_up)
        do ii=1,ip
            if (infout1==0) then  !inflow
                do k=kparasta,kparaend !1,jz
                    do j=1,n2
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
                    do j=1,n2
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
                        do j=1,n2
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
                        do j=1,n2
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
                        do j=1,n2

                            u(0,j,k)=ucp1(j,k)*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                            ! v(0,j,k)=ucp1(j,k)*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                            w(0,j,k)=ucp1(j,k)*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                            v(0,j,k)=v(1,j,k)
                            ! w(0,j,k)=w(1,j,k)

                            uc(0,j,k)=ucp1(j,k)

                            do isc=1,nscal
                                if (.not.potenziale) then
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

            ! face 2 'destra' (for periodicity see vel_up)
            if (infout2==0) then  !inflow
                do k=kparasta,kparaend
                    do j=1,n2
                        u(n1+1,j,k)=up2(j,k)
                        v(n1+1,j,k)=vp2(j,k)
                        w(n1+1,j,k)=wp2(j,k)
                        do isc=1,nscal
                            rhov(isc,n1+1,j,k)=rhovp2(isc,j,k)
                        end do
                        uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 2'
            else if (infout2==1) then  !outflow
                do k=kparasta,kparaend
                    do j=1,n2
                        u(n1+1,j,k)=u(n1+1,j,k)
                        v(n1+1,j,k)=v(n1+1,j,k)
                        w(n1+1,j,k)=w(n1+1,j,k)
                        do isc=1,nscal
                            rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 2'
            else if (infout2==2) then  !wall
                if (iboun2==0) then
                    do k=kparasta,kparaend
                        do j=1,n2
                            u(n1+1,j,k)=0.
                            v(n1+1,j,k)=0.
                            w(n1+1,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                            end do
                            uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                        end do
                    end do
                else if (iboun2==1) then
                    do k=kparasta,kparaend
                        do j=1,n2
                            u(n1+1,j,k)=0.
                            v(n1+1,j,k)=v(n1,j,k)
                            w(n1+1,j,k)=w(n1,j,k)
                            do isc=1,nscal
                                rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                            end do
                            uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                        end do
                    end do
                else if (iboun2==2) then
                    do k=kparasta,kparaend
                        do j=1,n2

                            u(n1+1,j,k)=ucp2(j,k)*csx(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                            ! v(jx+1,j,k)=ucp2(j,k)*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                            w(n1+1,j,k)=ucp2(j,k)*csz(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                            v(n1+1,j,k)=v(n1,j,k)
                            ! w(jx+1,j,k)=w(jx,j,k)

                            uc(n1,j,k)=ucp2(j,k)

                            do isc=1,nscal
                                if (.not.potenziale) then
                                    rhov(isc,n1+1,j,k)=rhovp2(isc,j,k)
                                else
                                    rhov(isc,n1+1,j,k)=rhovo2(isc,j,k)
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
        ! face 3 'sotto' (for periodicity see vel_up)
        !
        ! direction 2 is always not periodic

        if (infout3==0) then  !inflow
            do k=kparasta,kparaend
                do i=1,n1
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
                do i=1,n1
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
                    do i=1,n1
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
                    do i=1,n1
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
                    do i=1,n1
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

        ! face 4 'sopra' (for periodicity see vel_up)
        if (infout4==0) then  !inflow
            do k=kparasta,kparaend
                do i=1,n1
                    u(i,n2+1,k)=up4(i,k)
                    v(i,n2+1,k)=vp4(i,k)
                    w(i,n2+1,k)=wp4(i,k)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,k)=rhovp4(isc,i,k)
                    end do
                    vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                end do
            end do
            if (myid==0)write(*,*)'inflow face 4'
        else if (infout4==1) then  !outflow
            do k=kparasta,kparaend
                do i=1,n1
                    u(i,n2+1,k)=u(i,n2+1,k)
                    v(i,n2+1,k)=v(i,n2+1,k)
                    w(i,n2+1,k)=w(i,n2+1,k)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                    end do
                end do
            end do
            if (myid==0)write(*,*)'outflow face 4'
        else if (infout4==2) then  !wall
            if (iboun4==0) then
                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=0.
                        v(i,n2+1,k)=0.
                        w(i,n2+1,k)=0.
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                        vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                    end do
                end do
            else if (iboun4==1) then
                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=u(i,n2,k)
                        v(i,n2+1,k)=0.
                        w(i,n2+1,k)=w(i,n2,k)
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                        vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                    end do
                end do
            else if (iboun4==2) then
                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=0.
                        v(i,n2+1,k)=0.
                        w(i,n2+1,k)=0.
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                        vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                    end do
                end do
            end if
            if (myid==0)write(*,*)'solid wall face 4'
        end if
        !
        !-----------------------------------------------------------------------
        ! face 5 'indietro' (for periodicity see vel_up)
        !
        do kk=1,kp
            if (myid==0) then
                if (infout5==0) then  !inflow
                    do j=1,n2
                        do i=1,n1
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
                    do j=1,n2
                        do i=1,n1
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
                        do j=1,n2
                            do i=1,n1
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
                        do j=1,n2
                            do i=1,n1
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
                        do j=1,n2
                            do i=1,n1
                                ! u(i,j,0)=u(i,j,1)
                                v(i,j,0)=v(i,j,1)

                                u(i,j,0)=wcp5(i,j)*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                ! v(i,j,0)=wcp5(i,j)*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                w(i,j,0)=wcp5(i,j)*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                wc(i,j,0)=wcp5(i,j)

                                do isc=1,nscal
                                    if (.not.potenziale) then
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
            ! face 6 'avanti' (for periodicity see vel_up)
            !
            if (myid==nproc-1) then
                if (infout6==0) then  !inflow
                    do j=1,n2
                        do i=1,n1
                            u(i,j,n3+1)=up6(i,j)
                            v(i,j,n3+1)=vp6(i,j)
                            w(i,j,n3+1)=wp6(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,n3+1)=rhovp6(isc,i,j)
                            end do
                            wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)

                        end do
                    end do
                    write(*,*)'inflow face 6'
                else if (infout6==1) then  !outflow
                    do j=1,n2
                        do i=1,n1
                            u(i,j,n3+1)=u(i,j,n3+1)
                            v(i,j,n3+1)=v(i,j,n3+1)
                            w(i,j,n3+1)=w(i,j,n3+1)
                            do isc=1,nscal
                                rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                            end do
                        end do
                    end do
                    write(*,*)'outflow face 6'
                else if (infout6==2) then  !wall
                    if (iboun6==0) then
                        do j=1,n2
                            do i=1,n1
                                u(i,j,n3+1)=0.
                                v(i,j,n3+1)=0.
                                w(i,j,n3+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                                end do
                                wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)
                            end do
                        end do
                    else if (iboun6==1) then
                        do j=1,n2
                            do i=1,n1
                                u(i,j,n3+1)=u(i,j,n3)
                                v(i,j,n3+1)=v(i,j,n3)
                                w(i,j,n3+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                                end do
                                wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)
                            end do
                        end do
                    else if (iboun6==2) then
                        do j=1,n2
                            do i=1,n1
                                ! u(i,j,jz+1)=u(i,j,jz)
                                v(i,j,n3+1)=v(i,j,n3)

                                u(i,j,n3+1)=wcp6(i,j)*ztx(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                ! v(i,j,jz+1)=wcp6(i,j)*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                                w(i,j,n3+1)=wcp6(i,j)*ztz(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                wc(i,j,n3)=wcp6(i,j)

                                do isc=1,nscal
                                    if (.not.potenziale) then
                                        rhov(isc,i,j,n3+1)=rhovp6(isc,i,j)
                                    else
                                        rhov(isc,i,j,n3+1)=rhovo6(isc,i,j)
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
        ! boundary condition on cartesian velocities and scalars and
        ! for controvariant fluxes
        !
        use myarrays_metri3

        implicit none

        !-----------------------------------------------------------------------
        ! arrays declaration
        !
        integer i,j,k
        integer ii,jj,kk
        integer isc
        !
        !-----------------------------------------------------------------------
        ! face 1 'sinistra' (for periodicity see vel_up)
        do ii=1,ip
            if (infout1==0) then  !inflow
                do k=kparasta,kparaend !1,jz
                    do j=1,n2
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
                    do j=1,n2
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
                        do j=1,n2
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
                        do j=1,n2
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

            ! face 2 'destra' (for periodicity see vel_up)
            if (infout2==0) then  !inflow
                do k=kparasta,kparaend
                    do j=1,n2
                        u(n1+1,j,k)=up2(j,k)
                        v(n1+1,j,k)=vp2(j,k)
                        w(n1+1,j,k)=wp2(j,k)
                        do isc=1,nscal
                            rhov(isc,n1+1,j,k)=rhovp2(isc,j,k)
                        end do
                        uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                    end do
                end do
                if (myid==0)write(*,*)'inflow face 2'
            else if (infout2==1) then  !outflow
                do k=kparasta,kparaend
                    do j=1,n2
                        u(n1+1,j,k)=u(n1+1,j,k)
                        v(n1+1,j,k)=v(n1+1,j,k)
                        w(n1+1,j,k)=w(n1+1,j,k)
                        do isc=1,nscal
                            rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                        end do
                    end do
                end do
                if (myid==0)write(*,*)'outflow face 2'
            else if (infout2==2) then  !wall
                if (iboun2==0) then
                    do k=kparasta,kparaend
                        do j=1,n2
                            u(n1+1,j,k)=0.
                            v(n1+1,j,k)=0.
                            w(n1+1,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                            end do
                            uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                        end do
                    end do
                else if (iboun2==1) then
                    do k=kparasta,kparaend
                        do j=1,n2
                            u(n1+1,j,k)=0.
                            v(n1+1,j,k)=v(n1,j,k)
                            w(n1+1,j,k)=w(n1,j,k)
                            do isc=1,nscal
                                rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                            end do
                            uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                        end do
                    end do
                else if (iboun2==2) then

                end if
                if (myid==0)write(*,*)'solid wall face 2'
            end if
        end do   !end loop  ii=1,ip
        !
        !-----------------------------------------------------------------------
        ! face 3 'sotto' (for periodicity see vel_up)
        ! direction 2 is always not periodic

        if (infout3==0) then  !inflow
            do k=kparasta,kparaend
                do i=1,n1
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
                do i=1,n1
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
                    do i=1,n1
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
                    do i=1,n1
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

        ! face 4 'sopra' (for periodicity see vel_up)
        if (infout4==0) then  !inflow
            do k=kparasta,kparaend
                do i=1,n1
                    u(i,n2+1,k)=up4(i,k)
                    v(i,n2+1,k)=vp4(i,k)
                    w(i,n2+1,k)=wp4(i,k)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,k)=rhovp4(isc,i,k)
                    end do
                    vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                end do
            end do
            if (myid==0)write(*,*)'inflow face 4'
        else if (infout4==1) then  !outflow
            do k=kparasta,kparaend
                do i=1,n1
                    u(i,n2+1,k)=u(i,n2+1,k)
                    v(i,n2+1,k)=v(i,n2+1,k)
                    w(i,n2+1,k)=w(i,n2+1,k)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                    end do
                end do
            end do
            if (myid==0)write(*,*)'outflow face 4'
        else if (infout4==2) then  !wall
            if (iboun4==0) then
                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=1.
                        v(i,n2+1,k)=0.
                        w(i,n2+1,k)=1.
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                        vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k) +ety(i,n2,k)*v(i,n2+1,k) +etz(i,n2,k)*w(i,n2+1,k)
                    end do
                end do
            else if (iboun4==1) then
                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=u(i,n2,k)
                        v(i,n2+1,k)=0.
                        w(i,n2+1,k)=w(i,n2,k)
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                        vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                    end do
                end do
            else if (iboun4==2) then

            end if
            if (myid==0)write(*,*)'solid wall face 4'
        end if
        !
        !-----------------------------------------------------------------------
        ! face 5 'indietro' (for periodicity see vel_up)
        !
        do kk=1,kp
            if (myid==0) then
                if (infout5==0) then  !inflow
                    do j=1,n2
                        do i=1,n1
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
                    do j=1,n2
                        do i=1,n1
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
                        do j=1,n2
                            do i=1,n1
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
                        do j=1,n2
                            do i=1,n1
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
            ! face 6 'avanti' (for periodicity see vel_up)
            !
            if (myid==nproc-1) then
                if (infout6==0) then  !inflow
                    do j=1,n2
                        do i=1,n1
                            u(i,j,n3+1)=up6(i,j)
                            v(i,j,n3+1)=vp6(i,j)
                            w(i,j,n3+1)=wp6(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,n3+1)=rhovp6(isc,i,j)
                            end do
                            wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)

                        end do
                    end do
                    write(*,*)'inflow face 6'
                else if (infout6==1) then  !outflow
                    do j=1,n2
                        do i=1,n1
                            u(i,j,n3+1)=u(i,j,n3+1)
                            v(i,j,n3+1)=v(i,j,n3+1)
                            w(i,j,n3+1)=w(i,j,n3+1)
                            do isc=1,nscal
                                rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                            end do
                        end do
                    end do
                    write(*,*)'outflow face 6'
                else if (infout6==2) then  !wall
                    if (iboun6==0) then
                        do j=1,n2
                            do i=1,n1
                                u(i,j,n3+1)=0.
                                v(i,j,n3+1)=0.
                                w(i,j,n3+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                                end do
                                wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)
                            end do
                        end do
                    else if (iboun6==1) then
                        do j=1,n2
                            do i=1,n1
                                u(i,j,n3+1)=u(i,j,n3)
                                v(i,j,n3+1)=v(i,j,n3)
                                w(i,j,n3+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                                end do
                                wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)
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
        ! set boundary condition for cartesian velocity and controvariant flux
        !
        use mysettings, only: freesurface,i_rest,windyes
        use myarrays_metri3, only: csx,csy,csz,etx,ety,etz,ztx,zty,ztz

        implicit none

        !-----------------------------------------------------------------------
        integer i,j,k,isc
        integer ierr
        integer kparastal,kparaendl
        real delmassa1,delmassa2,delmassa3
        real delmassa4,delmassa5,delmassa6
        real ucp_old,wcp_old
        !-----------------------------------------------------------------------

        if (i_rest==3) then
            if (infout1/=0) iboun1=2
            if (infout2/=0) iboun2=2
            if (infout5/=0) iboun5=2
            if (infout6/=0) iboun6=2
        end if

        ! (for periodicity see vel_up)
        !chicco portare right and left dentro un blocco do unico
        ! left side
        if (ip==1) then
            !
            if (infout1==0) then        !inflow
                do k=kparasta,kparaend
                    do j=1,n2
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
                    do j=1,n2

                        uc(0,j,k)=uc1_orl(j,k)

                        delmassa1=uc(0,j,k)-(du_dx1(j,k)*csx(0,j,k)+dv_dx1(j,k)*csy(0,j,k)+dw_dx1(j,k)*csz(0,j,k))

                        u(0,j,k)=du_dx1(j,k)+delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                        v(0,j,k)=dv_dx1(j,k)+delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                        w(0,j,k)=dw_dx1(j,k)+delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                        do isc=1,nscal
                            rhov(isc,0,j,k)=rhov(isc,1,j,k)
                        end do

                    end do
                end do

            else if (infout1==2) then      ! wall
                if (iboun1==0) then
                    do k=kparasta,kparaend
                        do j=1,n2
                            u(0,j,k)=0.0
                            v(0,j,k)=0.0
                            w(0,j,k)=0.0
                            do isc=1,nscal
                                rhov(isc,0,j,k)=rhov(isc,1,j,k)
                            end do
                            uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
                        end do
                    end do
                else if (iboun1==1) then
                    do k=kparasta,kparaend
                        do j=1,n2
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
                        do j=1,n2

                            if ( ucp1(j,k) .gt. 0.) then

                                ucp_old=csx(0,j,k)*up1(j,k)+csy(0,j,k)*vp1(j,k)+csz(0,j,k)*wp1(j,k)

                                delmassa1=ucp1(j,k)-ucp_old

                                u(0,j,k)=up1(j,k)+delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                ! v(0,j,k)=v(1,j,k)
                                v(0,j,k)=vp1(j,k)+delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                w(0,j,k)=wp1(j,k)+delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                            else

                                ucp_old=csx(0,j,k)*u(1,j,k)+ csy(0,j,k)*v(1,j,k)+ csz(0,j,k)*w(1,j,k)

                                delmassa1=ucp1(j,k)-ucp_old


                                u(0,j,k)=u(1,j,k)+delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                ! v(0,j,k)=v(1,j,k)
                                v(0,j,k)=v(1,j,k)+delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                                w(0,j,k)=w(1,j,k)+delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                            end if


                            uc(0,j,k)=ucp1(j,k)

                            do isc=1,nscal
                                if (.not.potenziale) then
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

            ! right side
            if (infout2==0) then        !inflow
                do k=kparasta,kparaend
                    do j=1,n2

                        u(n1+1,j,k)=up2(j,k)
                        v(n1+1,j,k)=vp2(j,k)
                        w(n1+1,j,k)=wp2(j,k)
                        do isc=1,nscal
                            rhov(isc,n1+1,j,k)=rhovp2(isc,j,k)
                        end do
                        uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                    end do
                end do
            else if (infout2==1) then   !outflow
                do k=kparasta,kparaend
                    do j=1,n2

                        uc(n1,j,k)=uc2_orl(j,k)

                        delmassa2=uc(n1,j,k)-(du_dx2(j,k)*csx(n1,j,k)+dv_dx2(j,k)*csy(n1,j,k)+dw_dx2(j,k)*csz(n1,j,k))

                        u(n1+1,j,k)=du_dx2(j,k)+delmassa2*csx(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)
                        v(n1+1,j,k)=dv_dx2(j,k)+delmassa2*csy(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)
                        w(n1+1,j,k)=dw_dx2(j,k)+delmassa2*csz(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)
                        do isc=1,nscal
                            rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                        end do


                    end do
                end do

            else if (infout2==2) then      !wall
                if (iboun2==0) then
                    do k=kparasta,kparaend
                        do j=1,n2
                            u(n1+1,j,k)=0.
                            v(n1+1,j,k)=0.
                            w(n1+1,j,k)=0.
                            do isc=1,nscal
                                rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                            end do
                            uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                        end do
                    end do
                else if (iboun2==1) then
                    do k=kparasta,kparaend
                        do j=1,n2
                            u(n1+1,j,k)=0.
                            v(n1+1,j,k)=v(n1,j,k)
                            w(n1+1,j,k)=w(n1,j,k)
                            do isc=1,nscal
                                rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                            end do
                            uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k)+csy(n1,j,k)*v(n1+1,j,k)+csz(n1,j,k)*w(n1+1,j,k)
                        end do
                    end do
                else if (iboun2==2) then
                    do k=kparasta,kparaend
                        do j=1,n2

                            if (ucp2(j,k) .lt. 0.) then

                                ucp_old=csx(n1,j,k)*up2(j,k)+ csy(n1,j,k)*vp2(j,k)+ csz(n1,j,k)*wp2(j,k)

                                delmassa2=ucp2(j,k)-ucp_old

                                u(n1+1,j,k)=up2(j,k)+delmassa2*csx(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                                ! v(jx+1,j,k)=v(1,j,k)
                                v(n1+1,j,k)=vp2(j,k)+delmassa2*csy(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                                w(n1+1,j,k)=wp2(j,k)+delmassa2*csz(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                                uc(n1,j,k)=ucp2(j,k)

                            else

                                ucp_old=csx(n1,j,k)*u(n1,j,k)+ csy(n1,j,k)*v(n1,j,k)+csz(n1,j,k)*w(n1,j,k)

                                delmassa2=ucp2(j,k)-ucp_old

                                u(n1+1,j,k)=u(n1,j,k)+delmassa2*csx(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                                ! v(jx+1,j,k)=v(1,j,k)
                                v(n1+1,j,k)=v(n1,j,k)+delmassa2*csy(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                                w(n1+1,j,k)=w(n1,j,k)+delmassa2*csz(n1,j,k)/(csx(n1,j,k)**2+csy(n1,j,k)**2+csz(n1,j,k)**2)

                                uc(n1,j,k)=ucp2(j,k)

                            end if

                            do isc=1,nscal
                                if (.not.potenziale) then
                                    rhov(isc,n1+1,j,k)=rhovp2(isc,j,k)
                                    if ( ucp2(j,k) .gt. 0.) then
                                        rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                                    end if

                                else
                                    rhov(isc,n1+1,j,k)=rhovo2(isc,j,k)
                                end if
                            end do

                        end do
                    end do
                end if
            end if

        end if   !end loop  ii=1,ip
        !-----------------------------------------------------------------------
        !
        ! bottom side
        ! direction 2 is always not periodic

        if (infout3==0) then             !inflow
            do k=kparasta,kparaend
                do i=1,n1

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
                do i=1,n1

                    vc(i,0,k)=vc3_orl(i,k)

                    delmassa3=vc(i,0,k)-(du_dy3(i,k)*etx(i,0,k)+dv_dy3(i,k)*ety(i,0,k)+dw_dy3(i,k)*etz(i,0,k))

                    u(i,0,k)=du_dy3(i,k)+delmassa3*etx(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)
                    v(i,0,k)=dv_dy3(i,k)+delmassa3*ety(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)
                    w(i,0,k)=dw_dy3(i,k)+delmassa3*etz(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)

                    do isc=1,nscal
                        rhov(isc,i,0,k)=rhov(isc,i,1,k)
                    end do


                end do
            end do

        else if (infout3==2) then           !wall
            if (iboun3==0) then
                do k=kparasta,kparaend
                    do i=1,n1
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
                    do i=1,n1
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
                    do i=1,n1
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
        ! upper side
        if (.not.freesurface) then !no freesurface<<<<<<<<<<<<<<<<<<<<<<<<<<

            if (infout4==0) then             !inflow
                do k=kparasta,kparaend
                    do i=1,n1

                        u(i,n2+1,k)=up4(i,k)
                        v(i,n2+1,k)=vp4(i,k)
                        w(i,n2+1,k)=wp4(i,k)
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhovp4(isc,i,k)
                        end do
                        vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                    end do
                end do
            else if (infout4==1) then        !outflow
                do k=kparasta,kparaend
                    do i=1,n1

                        vc(i,n2,k)=vc4_orl(i,k)

                        delmassa4=vc(i,n2,k)-(du_dy4(i,k)*etx(i,n2,k)+dv_dy4(i,k)*ety(i,n2,k)+dw_dy4(i,k)*etz(i,n2,k))

                        u(i,n2+1,k)=du_dy4(i,k)+delmassa4*etx(i,n2,k)/(etx(i,n2,k)**2+ety(i,n2,k)**2+etz(i,n2,k)**2)
                        v(i,n2+1,k)=dv_dy4(i,k)+delmassa4*ety(i,n2,k)/(etx(i,n2,k)**2+ety(i,n2,k)**2+etz(i,n2,k)**2)
                        w(i,n2+1,k)=dw_dy4(i,k)+delmassa4*etz(i,n2,k)/(etx(i,n2,k)**2+ety(i,n2,k)**2+etz(i,n2,k)**2)

                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do



                    end do
                end do

            ! windyes=1 option,if at the upper side there is a wind stress

            else if (infout4==2 .and. windyes==0) then ! wall without wind
                if (iboun4==0) then
                    do k=kparasta,kparaend
                        do i=1,n1
                            u(i,n2+1,k)=0.
                            v(i,n2+1,k)=0.
                            w(i,n2+1,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                            end do
                            vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k) +ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                        end do
                    end do
                else if (iboun4==1) then
                    do k=kparasta,kparaend
                        do i=1,n1
                            u(i,n2+1,k)=u(i,n2,k)
                            v(i,n2+1,k)=0.
                            w(i,n2+1,k)=w(i,n2,k)
                            do isc=1,nscal
                                rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                            end do
                            vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                            if ((i==4).and.(k==17)) then
                                write(*,*)'VC_contour:',vc(i,32,k)
                            end if
                        end do
                    end do
                else if (iboun4==2) then
                    do k=kparasta,kparaend
                        do i=1,n1
                            u(i,n2+1,k)=0.
                            v(i,n2+1,k)=0.
                            w(i,n2+1,k)=0.
                            do isc=1,nscal
                                rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                            end do
                            vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                        end do
                    end do
                end if

            else if (infout4==2 .and. windyes==1) then ! wall with wind extrapolation from the interior

                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=1.5*u(i,n2,k)-.5*u(i,n2-1,k)
                        v(i,n2+1,k)=0.0
                        w(i,n2+1,k)=1.5*w(i,n2,k)-.5*w(i,n2-1,k)
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                        vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k)+ety(i,n2,k)*v(i,n2+1,k)+etz(i,n2,k)*w(i,n2+1,k)
                    end do
                end do
            end if

        else if (freesurface) then !free surface ON.<<<<<<<<<
            if (windyes==0) then
                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=u(i,n2,k)
                        v(i,n2+1,k)=v(i,n2,k)
                        w(i,n2+1,k)=w(i,n2,k)
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                    ! vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)
                    ! >               +ety(i,jy,k)*v(i,jy+1,k)
                    ! >               +etz(i,jy,k)*w(i,jy+1,k)
                    end do
                end do
            else if (windyes==1) then
                do k=kparasta,kparaend
                    do i=1,n1
                        u(i,n2+1,k)=1.5*u(i,n2,k)-.5*u(i,n2-1,k)
                        v(i,n2+1,k)=v(i,n2,k)
                        w(i,n2+1,k)=1.5*w(i,n2,k)-.5*w(i,n2-1,k)
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                    ! vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)
                    ! >               +ety(i,jy,k)*v(i,jy+1,k)
                    ! >               +etz(i,jy,k)*w(i,jy+1,k)
                    end do
                end do
            end if
        end if !free surface<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !-----------------------------------------------------------------------
        !
        !
        if (kp==1) then
            !
            ! back side
            if (myid==0) then
                if (infout5==0) then            !inflow
                    do j=1,n2
                        do i=1,n1

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
                    do j=1,n2
                        do i=1,n1

                            wc(i,j,0)=wc5_orl(i,j)

                            delmassa5=wc(i,j,0)-(du_dz5(i,j)*ztx(i,j,0)+dv_dz5(i,j)*zty(i,j,0)+dw_dz5(i,j)*ztz(i,j,0))

                            u(i,j,0)=du_dz5(i,j)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
                            v(i,j,0)=dv_dz5(i,j)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
                            w(i,j,0)=dw_dz5(i,j)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                            do isc=1,nscal
                                rhov(isc,i,j,0)=rhov(isc,i,j,1)
                            end do


                        end do
                    end do


                else if (infout5==2) then            !wall
                    if (iboun5==0) then
                        do j=1,n2
                            do i=1,n1
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
                        do j=1,n2
                            do i=1,n1
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
                        do j=1,n2
                            do i=1,n1

                                if ( wcp5(i,j) .gt. 0.) then

                                    wcp_old=ztx(i,j,0)*up5(i,j)+zty(i,j,0)*vp5(i,j)+ztz(i,j,0)*wp5(i,j)

                                    delmassa5=wcp5(i,j)-wcp_old

                                    u(i,j,0)=up5(i,j)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    ! v(i,j,0)=v(i,j,1)
                                    v(i,j,0)=vp5(i,j)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    w(i,j,0)=wp5(i,j)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    wc(i,j,0)=wcp5(i,j)

                                else

                                    wcp_old=ztx(i,j,0)*u(i,j,1)+zty(i,j,0)*v(i,j,1)+ztz(i,j,0)*w(i,j,1)

                                    delmassa5=wcp5(i,j)-wcp_old

                                    u(i,j,0)=u(i,j,1)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    ! v(i,j,0)=v(i,j,1)
                                    v(i,j,0)=v(i,j,1)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    w(i,j,0)=w(i,j,1)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                                    wc(i,j,0)=wcp5(i,j)

                                end if

                                do isc=1,nscal
                                    if (.not.potenziale) then
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
            ! front side
            if (myid==nproc-1) then
                if (infout6==0) then            !inflow
                    do j=1,n2
                        do i=1,n1

                            u(i,j,n3+1)=up6(i,j)
                            v(i,j,n3+1)=vp6(i,j)
                            w(i,j,n3+1)=wp6(i,j)
                            do isc=1,nscal
                                rhov(isc,i,j,n3+1)=rhovp6(isc,i,j)
                            end do
                            wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)

                        end do
                    end do
                else if (infout6==1) then            !outflow
                    do j=1,n2
                        do i=1,n1

                            wc(i,j,n3)=wc6_orl(i,j)

                            delmassa6=wc(i,j,n3)-(du_dz6(i,j)*ztx(i,j,n3) +dv_dz6(i,j)*zty(i,j,n3)+dw_dz6(i,j)*ztz(i,j,n3))

                            u(i,j,n3+1)=du_dz6(i,j)+delmassa6*ztx(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)
                            v(i,j,n3+1)=dv_dz6(i,j)+delmassa6*zty(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)
                            w(i,j,n3+1)=dw_dz6(i,j)+delmassa6*ztz(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                            do isc=1,nscal
                                rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                            end do

                        end do
                    end do


                else if (infout6==2) then            !wall
                    if (iboun6==0) then
                        do j=1,n2
                            do i=1,n1
                                u(i,j,n3+1)=0.
                                v(i,j,n3+1)=0.
                                w(i,j,n3+1)=0.
                                do isc=1,nscal
                                    rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                                end do
                                wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)
                            end do
                        end do
                    else if (iboun6==1) then
                        do j=1,n2
                            do i=1,n1
                                u(i,j,n3+1)=u(i,j,n3)
                                v(i,j,n3+1)=v(i,j,n3)
                                w(i,j,n3+1)=0.



                                do isc=1,nscal
                                    rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                                end do
                                wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1)+zty(i,j,n3)*v(i,j,n3+1)+ztz(i,j,n3)*w(i,j,n3+1)
                            end do
                        end do
                    else if (iboun6==2) then
                        do j=1,n2
                            do i=1,n1

                                if (wcp6(i,j) .lt. 0.) then
                                    wcp_old=ztx(i,j,n3)*up6(i,j)+zty(i,j,n3)*vp6(i,j)+ztz(i,j,n3)*wp6(i,j)

                                    delmassa6=wcp6(i,j)-wcp_old

                                    u(i,j,n3+1)=up6(i,j)+delmassa6*ztx(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                    ! v(i,j,jz+1)=v(i,j,jz)
                                    v(i,j,n3+1)=vp6(i,j)+delmassa6*zty(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                    w(i,j,n3+1)=wp6(i,j)+delmassa6*ztz(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                    wc(i,j,n3)=wcp6(i,j)

                                else

                                    wcp_old=ztx(i,j,n3)*u(i,j,n3)+zty(i,j,n3)*v(i,j,n3)+ztz(i,j,n3)*w(i,j,n3)

                                    delmassa6=wcp6(i,j)-wcp_old

                                    u(i,j,n3+1)=u(i,j,n3)+delmassa6*ztx(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                    ! v(i,j,jz+1)=v(i,j,jz)
                                    v(i,j,n3+1)=v(i,j,n3)+delmassa6*zty(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                    w(i,j,n3+1)=w(i,j,n3)+delmassa6*ztz(i,j,n3)/(ztx(i,j,n3)**2+zty(i,j,n3)**2+ztz(i,j,n3)**2)

                                    wc(i,j,n3)=wcp6(i,j)

                                end if

                                do isc=1,nscal
                                    if (.not.potenziale) then
                                        rhov(isc,i,j,n3+1)=rhovp6(isc,i,j)

                                        if (wcp6(i,j).gt.0.) then
                                            rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                                        end if
                                    else
                                        rhov(isc,i,j,n3+1)=rhovo6(isc,i,j)
                                    end if
                                end do

                            end do
                        end do
                    end if
                end if
            end if

        end if
        !-----------------------------------------------------------------------
        ! set boundary condition at the corner

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

        ! corner close to side 5
        if (myid==0) then

            ! corner 5/1
            if (infout5==1 .and. infout1==1) then
                do j=0,n2+1
                    u(0,j,0)=u(1,j,1)
                    v(0,j,0)=v(1,j,1)
                    w(0,j,0)=w(1,j,1)
                    do isc=1,nscal
                        rhov(isc,0,j,0)=rhov(isc,1,j,1)
                    end do
                end do
            end if  !corner 5/1

            ! corner 5/2
            if (infout5==1 .and. infout2==1) then
                do j=0,n2+1
                    u(n1+1,j,0)=u(n1,j,1)
                    v(n1+1,j,0)=v(n1,j,1)
                    w(n1+1,j,0)=w(n1,j,1)
                    do isc=1,nscal
                        rhov(isc,n1+1,j,0)=rhov(isc,n1,j,1)
                    end do
                end do
            end if  !corner 5/2

            ! corner 5/3
            if (infout5==1 .and. infout3==1) then
                do i=0,n1+1
                    u(i,0,0)=u(i,1,1)
                    v(i,0,0)=v(i,1,1)
                    w(i,0,0)=w(i,1,1)
                    do isc=1,nscal
                        rhov(isc,i,0,0)=rhov(isc,i,1,1)
                    end do
                end do
            end if  !corner 5/3

            ! corner 5/4
            if (infout5==1 .and. infout4==1) then
                do i=0,n1+1
                    u(i,n2+1,0)=u(i,n2,1)
                    v(i,n2+1,0)=v(i,n2,1)
                    w(i,n2+1,0)=w(i,n2,1)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,0)=rhov(isc,i,n2,1)
                    end do
                end do
            end if  !corner 5/4

        end if !proc 0
        !-----------------------------------------------------------------------
        ! corner close to side 6
        if (myid==nproc-1) then
            ! corner 6/1
            if (infout6==1 .and. infout1==1) then
                do j=0,n2+1
                    u(0,j,n3+1)=u(1,j,n3)
                    v(0,j,n3+1)=v(1,j,n3)
                    w(0,j,n3+1)=w(1,j,n3)
                    do isc=1,nscal
                        rhov(isc,0,j,n3+1)=rhov(isc,1,j,n3)
                    end do
                end do
            end if  !corner 6/1

            ! corner 6/2
            if (infout6==1 .and. infout2==1) then
                do j=0,n2+1
                    u(n1+1,j,n3+1)=u(n1,j,n3)
                    v(n1+1,j,n3+1)=v(n1,j,n3)
                    w(n1+1,j,n3+1)=w(n1,j,n3)
                    do isc=1,nscal
                        rhov(isc,n1+1,j,n3+1)=rhov(isc,n1,j,n3)
                    end do
                end do
            end if  !corner 6/2

            ! corner 6/3
            if (infout6==1 .and. infout3==1) then
                do i=0,n1+1
                    u(i,0,n3+1)=u(i,1,n3)
                    v(i,0,n3+1)=v(i,1,n3)
                    w(i,0,n3+1)=w(i,1,n3)
                    do isc=1,nscal
                        rhov(isc,i,0,n3+1)=rhov(isc,i,1,n3)
                    end do
                end do
            end if  !corner 6/3

            ! corner 6/4
            if (infout6==1 .and. infout4==1) then
                do i=0,n1+1
                    u(i,n2+1,n3+1)=u(i,n2,n3)
                    v(i,n2+1,n3+1)=v(i,n2,n3)
                    w(i,n2+1,n3+1)=w(i,n2,n3)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,n3+1)=rhov(isc,i,n2,n3)
                    end do
                end do
            end if  !corner 6/4

        end if ! nproc-1

        !-----------------------------------------------------------------------
        ! corner 1/3 e 1/4

        ! corner 1/3
        if (infout1==1 .and. infout3==1) then
            do k=kparastal,kparaendl
                u(0,0,k)=u(1,1,k)
                v(0,0,k)=v(1,1,k)
                w(0,0,k)=w(1,1,k)
                do isc=1,nscal
                    rhov(isc,0,0,k)=rhov(isc,1,1,k)
                end do
            end do
        end if  !corner 1/3

        ! corner 1/4
        if (infout1==1 .and. infout4==1) then
            do k=kparastal,kparaendl
                u(0,n2+1,k)=u(1,n2,k)
                v(0,n2+1,k)=v(1,n2,k)
                w(0,n2+1,k)=w(1,n2,k)
                do isc=1,nscal
                    rhov(isc,0,n2+1,k)=rhov(isc,1,n2,k)
                end do
            end do
        end if  !corner 1/4

        !-----------------------------------------------------------------------
        ! corner 2/3 e 2/4
        ! corner 2/3
        if (infout2==1 .and. infout3==1) then
            do k=kparastal,kparaendl
                u(n1+1,0,k)=u(n1,1,k)
                v(n1+1,0,k)=v(n1,1,k)
                w(n1+1,0,k)=w(n1,1,k)
                do isc=1,nscal
                    rhov(isc,n1+1,0,k)=rhov(isc,n1,1,k)
                end do
            end do
        end if  !corner 2/3

        ! corner 2/4
        if (infout2==1 .and. infout4==1) then
            do k=kparastal,kparaendl
                u(n1+1,n2+1,k)=u(n1,n2,k)
                v(n1+1,n2+1,k)=v(n1,n2,k)
                w(n1+1,n2+1,k)=w(n1,n2,k)
                do isc=1,nscal
                    rhov(isc,n1+1,n2+1,k)=rhov(isc,n1,n2,k)
                end do
            end do
        end if  !corner 2/4
        !-----------------------------------------------------------------------
        ! SECOND: case sides with inlow (infout=0)

        ! side 1
        if (infout1==0) then
            i=0
            do k=kparastal,kparaendl
                ! corner 1/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do

                ! corner 1/4
                u(i,n2+1,k)=u(i,n2,k)
                v(i,n2+1,k)=v(i,n2,k)
                w(i,n2+1,k)=w(i,n2,k)
                do isc=1,nscal
                    rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                end do
            end do

            do j=0,n2+1
                ! corner 1/5
                if (myid==0) then
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 1/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 1

        !-----------------------------------------------------------------------
        ! side 2
        if (infout2==0) then
            i=n1+1
            do k=kparastal,kparaendl
                ! corner 2/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do

                ! corner 2/4
                u(i,n2+1,k)=u(i,n2,k)
                v(i,n2+1,k)=v(i,n2,k)
                w(i,n2+1,k)=w(i,n2,k)
                do isc=1,nscal
                    rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                end do
            end do

            do j=0,n2+1
                if (myid==0) then
                    ! corner 2/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 2/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 2

        !-----------------------------------------------------------------------
        ! side 3
        if (infout3==0) then
            j=0
            do k=kparastal,kparaendl
                ! corner 3/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                ! corner 3/2
                u(n1+1,j,k)=u(n1,j,k)
                v(n1+1,j,k)=v(n1,j,k)
                w(n1+1,j,k)=w(n1,j,k)
                do isc=1,nscal
                    rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                end do
            end do

            do i=0,n1+1
                if (myid==0) then
                    ! corner 3/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 3/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 3

        !-----------------------------------------------------------------------
        ! side 4
        if (infout4==0) then
            j=n2+1
            do k=kparastal,kparaendl
                ! corner 4/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                ! corner 4/2
                u(n1+1,j,k)=u(n1,j,k)
                v(n1+1,j,k)=v(n1,j,k)
                w(n1+1,j,k)=w(n1,j,k)
                do isc=1,nscal
                    rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                end do
            end do

            do i=0,n1+1
                if (myid==0) then
                    ! corner 4/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 4/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 4

        !-----------------------------------------------------------------------
        ! side 5
        if (infout5==0) then
            if (myid==0) then
                k=0
                do j=0,n2+1
                    ! corner 5/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do

                    ! corner 5/2
                    u(n1+1,j,k)=u(n1,j,k)
                    v(n1+1,j,k)=v(n1,j,k)
                    w(n1+1,j,k)=w(n1,j,k)
                    do isc=1,nscal
                        rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                    end do
                end do

                do i=0,n1+1
                    ! corner 5/3
                    u(i,0,k)=u(i,1,k)
                    v(i,0,k)=v(i,1,k)
                    w(i,0,k)=w(i,1,k)
                    do isc=1,nscal
                        rhov(isc,i,0,k)=rhov(isc,i,1,k)
                    end do
                    ! corner 5/4
                    u(i,n2+1,k)=u(i,n2,k)
                    v(i,n2+1,k)=v(i,n2,k)
                    w(i,n2+1,k)=w(i,n2,k)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                    end do
                end do
            end if
        end if !side 5

        !-----------------------------------------------------------------------
        ! side 6
        if (infout6==0) then
            if (myid==nproc-1) then
                k=n3+1
                do j=0,n2+1
                    ! corner 6/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do
                    ! corner 6/2
                    u(n1+1,j,k)=u(n1,j,k)
                    v(n1+1,j,k)=v(n1,j,k)
                    w(n1+1,j,k)=w(n1,j,k)
                    do isc=1,nscal
                        rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                    end do
                end do

                do i=0,n1+1
                    if (myid==0) then
                        ! corner 6/3
                        u(i,0,k)=u(i,1,k)
                        v(i,0,k)=v(i,1,k)
                        w(i,0,k)=w(i,1,k)
                        do isc=1,nscal
                            rhov(isc,i,0,k)=rhov(isc,i,1,k)
                        end do
                    else if (myid==nproc-1) then
                        ! corner 6/4
                        u(i,n2+1,k)=u(i,n2,k)
                        v(i,n2+1,k)=v(i,n2,k)
                        w(i,n2+1,k)=w(i,n2,k)
                        do isc=1,nscal
                            rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                        end do
                    end if
                end do
            end if
        end if !side 6
        !-----------------------------------------------------------------------
        ! THIRD: case wall side (infout=2)
        !
        ! side 1
        if (infout1==2) then
            i=0
            do k=kparastal,kparaendl
                ! corner 1/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do
                ! corner 1/4
                u(i,n2+1,k)=u(i,n2,k)
                v(i,n2+1,k)=v(i,n2,k)
                w(i,n2+1,k)=w(i,n2,k)
                do isc=1,nscal
                    rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                end do
            end do

            do j=0,n2+1
                if (myid==0) then
                    ! corner 1/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 1/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 1

        !-----------------------------------------------------------------------
        ! side 2
        if (infout2==2) then
            i=n1+1
            do k=kparastal,kparaendl
                ! corner 2/3
                u(i,0,k)=u(i,1,k)
                v(i,0,k)=v(i,1,k)
                w(i,0,k)=w(i,1,k)
                do isc=1,nscal
                    rhov(isc,i,0,k)=rhov(isc,i,1,k)
                end do

                ! corner 2/4
                u(i,n2+1,k)=u(i,n2,k)
                v(i,n2+1,k)=v(i,n2,k)
                w(i,n2+1,k)=w(i,n2,k)
                do isc=1,nscal
                    rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                end do
            end do

            do j=0,n2+1
                if (myid==0) then
                    ! corner 2/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 2/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 2

        !-----------------------------------------------------------------------
        ! side 5
        if (infout5==2) then
            if (myid==0) then
                k=0
                do j=0,n2+1
                    ! corner 5/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do

                    ! corner 5/2
                    u(n1+1,j,k)=u(n1,j,k)
                    v(n1+1,j,k)=v(n1,j,k)
                    w(n1+1,j,k)=w(n1,j,k)
                    do isc=1,nscal
                        rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                    end do
                end do

                do i=0,n1+1
                    ! corner 5/3
                    u(i,0,k)=u(i,1,k)
                    v(i,0,k)=v(i,1,k)
                    w(i,0,k)=w(i,1,k)
                    do isc=1,nscal
                        rhov(isc,i,0,k)=rhov(isc,i,1,k)
                    end do

                    ! corner 5/4
                    u(i,n2+1,k)=u(i,n2,k)
                    v(i,n2+1,k)=v(i,n2,k)
                    w(i,n2+1,k)=w(i,n2,k)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                    end do
                end do
            end if
        end if !side 5

        !-----------------------------------------------------------------------
        ! side 6
        if (infout6==2) then
            if (myid==nproc-1) then
                k=n3+1
                do j=0,n2+1
                    ! corner 6/1
                    u(0,j,k)=u(1,j,k)
                    v(0,j,k)=v(1,j,k)
                    w(0,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,1,j,k)
                    end do

                    ! corner 6/2
                    u(n1+1,j,k)=u(n1,j,k)
                    v(n1+1,j,k)=v(n1,j,k)
                    w(n1+1,j,k)=w(n1,j,k)
                    do isc=1,nscal
                        rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                    end do
                end do
                do i=0,n1+1
                    ! corner 6/3
                    u(i,0,k)=u(i,1,k)
                    v(i,0,k)=v(i,1,k)
                    w(i,0,k)=w(i,1,k)
                    do isc=1,nscal
                        rhov(isc,i,0,k)=rhov(isc,i,1,k)
                    end do
                    ! corner 6/4
                    u(i,n2+1,k)=u(i,n2,k)
                    v(i,n2+1,k)=v(i,n2,k)
                    w(i,n2+1,k)=w(i,n2,k)
                    do isc=1,nscal
                        rhov(isc,i,n2+1,k)=rhov(isc,i,n2,k)
                    end do
                end do
            end if
        end if !side 6


        !-----------------------------------------------------------------------
        ! side 3
        if (infout3==2) then
            j=0
            do k=kparastal,kparaendl
                ! corner 3/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                ! corner 3/2
                u(n1+1,j,k)=u(n1,j,k)
                v(n1+1,j,k)=v(n1,j,k)
                w(n1+1,j,k)=w(n1,j,k)
                do isc=1,nscal
                    rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                end do
            end do

            do i=0,n1+1
                if (myid==0) then
                    ! corner 3/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 3/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 3

        !-----------------------------------------------------------------------
        ! side 4
        if (infout4==2) then
            j=n2+1
            do k=kparastal,kparaendl
                ! corner 4/1
                u(0,j,k)=u(1,j,k)
                v(0,j,k)=v(1,j,k)
                w(0,j,k)=w(1,j,k)
                do isc=1,nscal
                    rhov(isc,0,j,k)=rhov(isc,1,j,k)
                end do

                ! corner 4/2
                u(n1+1,j,k)=u(n1,j,k)
                v(n1+1,j,k)=v(n1,j,k)
                w(n1+1,j,k)=w(n1,j,k)
                do isc=1,nscal
                    rhov(isc,n1+1,j,k)=rhov(isc,n1,j,k)
                end do
            end do

            do i=0,n1+1
                if (myid==0) then
                    ! corner 4/5
                    u(i,j,0)=u(i,j,1)
                    v(i,j,0)=v(i,j,1)
                    w(i,j,0)=w(i,j,1)
                    do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                    end do
                else if (myid==nproc-1) then
                    ! corner 4/6
                    u(i,j,n3+1)=u(i,j,n3)
                    v(i,j,n3+1)=v(i,j,n3)
                    w(i,j,n3+1)=w(i,j,n3)
                    do isc=1,nscal
                        rhov(isc,i,j,n3+1)=rhov(isc,i,j,n3)
                    end do
                end if
            end do

        end if !side 4

        return

    end subroutine contour

    subroutine contourp()
        ! compute cartesian velocity and controvariant in periodic cell at
        ! step n+1,at the corner computation at the end of the sub

        use mysettings, only: bby,freesurface,i_rest,windyes

        implicit none

        !-----------------------------------------------------------------------
        ! array declaration
        integer i,j,k,ii,isc,kk
        integer ierr,status(MPI_STATUS_SIZE)
        integer kparastal,kparaendl
        integer plantype
        real,allocatable:: rho(:,:,:)
        !-----------------------------------------------------------------------

        ! periodic cell in csi (also outside the boundary)
        if (ip==0) then
            do k=kparasta,kparaend
                do j=0,n2+1

                    u(0,j,k)=u(n1,j,k)
                    v(0,j,k)=v(n1,j,k)
                    w(0,j,k)=w(n1,j,k)
                    u(n1+1,j,k)=u(1,j,k)
                    v(n1+1,j,k)=v(1,j,k)
                    w(n1+1,j,k)=w(1,j,k)
                    do isc=1,nscal
                        rhov(isc,0,j,k)=rhov(isc,n1,j,k)
                        rhov(isc,n1+1,j,k)=rhov(isc,1,j,k)
                    end do

                end do
            end do
        end if

        ! direction 2 is always not periodic

        ! periodic cell in zita (also outside the boundary)
        if (kp==0) then

            call MPI_TYPE_VECTOR(n2+2,n1,n1+2,MPI_REAL_SD,plantype,ierr)
            call MPI_TYPE_COMMIT(plantype,ierr)

            if (myid==0) then
                call MPI_SENDRECV(u(1,0,1),1,plantype,nproc-1,12,u(1,0,0),1,plantype,nproc-1,11,MPI_COMM_WORLD,status,ierr)
                call MPI_SENDRECV(v(1,0,1),1,plantype,nproc-1,14,v(1,0,0),1,plantype,nproc-1,13,MPI_COMM_WORLD,status,ierr)
                call MPI_SENDRECV(w(1,0,1),1,plantype,nproc-1,16,w(1,0,0),1,plantype,nproc-1,15,MPI_COMM_WORLD,status,ierr)
            end if

            if (myid==nproc-1) then
                call MPI_SENDRECV(u(1,0,n3),1,plantype,0,11,u(1,0,n3+1),1,plantype,0,12,MPI_COMM_WORLD,status,ierr)
                call MPI_SENDRECV(v(1,0,n3),1,plantype,0,13,v(1,0,n3+1),1,plantype,0,14,MPI_COMM_WORLD,status,ierr)
                call MPI_SENDRECV(w(1,0,n3),1,plantype,0,15,w(1,0,n3+1),1,plantype,0,16,MPI_COMM_WORLD,status,ierr)
            end if

            allocate(rho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))

            do isc=1,nscal
                do kk=kparasta-1,kparaend+1
                    do j=0,n2+1
                        do i=0,n1+1
                            rho(i,j,kk)=rhov(isc,i,j,kk)
                        end do
                    end do
                end do

                if (myid==0) then
                    call MPI_SENDRECV(rho(1,0,1),1,plantype,nproc-1,18+isc,rho(1,0,0), &
                        1,plantype,nproc-1,17+isc,MPI_COMM_WORLD,status,ierr)
                end if
                if (myid==nproc-1) then
                    call MPI_SENDRECV(rho(1,0,n3),1,plantype,0,17+isc,rho(1,0,n3+1),1,plantype,0,18+isc,MPI_COMM_WORLD,status,ierr)
                end if

                do kk=kparasta-1,kparaend+1
                    do j=0,n2+1
                        do i=0,n1+1
                            rhov(isc,i,j,kk)=rho(i,j,kk)
                        end do
                    end do
                end do
            end do !nscal
            deallocate(rho)
            call MPI_TYPE_FREE(plantype,ierr)

        end if
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

        if (freesurface) then !free surface ON.<<<<<<<<<
            if (myid==0) then
                write(*,*)'Free Surface ON'
            end if
            do k=kparastal,kparaendl
                do i=0,n1+1
                    next_prs(i,k)=((fi(i,n2+1,k)+fi(i,n2,k))*0.5)-(v(i,n2+1,k)*dt*bby)
                end do
            end do
        else
            if (myid==0) then
                write(*,*)'Free Surface OFF'
            end if
        end if !free surface<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        return
    end subroutine contourp

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
            do j=1,n2
                ! face 1 sinistra
                delu(0,j,k)=1.875*gra1(1,j,k)-1.25*gra1(2,j,k)+.375*gra1(3,j,k)
                delu(0,j,k)=-delu(0,j,k)
                !
                delv(0,j,k)=1.875*gra2(1,j,k)-1.25*gra2(2,j,k)+.375*gra2(3,j,k)
                delv(0,j,k)=-delv(0,j,k)
                !
                delw(0,j,k)=1.875*gra3(1,j,k)-1.25*gra3(2,j,k)+.375*gra3(3,j,k)
                delw(0,j,k)=-delw(0,j,k)
                !
                ! face 2 destra
                delu(n1+1,j,k)=.375*gra1(n1-2,j,k)-1.25*gra1(n1-1,j,k)+1.875*gra1(n1,j,k)
                delu(n1+1,j,k)=-delu(n1+1,j,k)
                !
                delv(n1+1,j,k)=.375*gra2(n1-2,j,k)-1.25*gra2(n1-1,j,k)+1.875*gra2(n1,j,k)
                delv(n1+1,j,k)=-delv(n1+1,j,k)
                !
                delw(n1+1,j,k)=.375*gra3(n1-2,j,k)-1.25*gra3(n1-1,j,k)+1.875*gra3(n1,j,k)
                delw(n1+1,j,k)=-delw(n1+1,j,k)
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
            do i=1,n1
                ! face 3 sotto
                delu(i,0,k)=1.875*gra1(i,1,k)-1.25*gra1(i,2,k)+.375*gra1(i,3,k)
                delu(i,0,k)=-delu(i,0,k)
                !
                delv(i,0,k)=1.875*gra2(i,1,k)-1.25*gra2(i,2,k)+.375*gra2(i,3,k)
                delv(i,0,k)=-delv(i,0,k)
                !
                delw(i,0,k)=1.875*gra3(i,1,k)-1.25*gra3(i,2,k)+.375*gra3(i,3,k)
                delw(i,0,k)=-delw(i,0,k)
                !
                ! face 4 sopra
                delu(i,n2+1,k)=.375*gra1(i,n2-2,k)-1.25*gra1(i,n2-1,k)+1.875*gra1(i,n2,k)
                delu(i,n2+1,k)=-delu(i,n2+1,k)
                !
                delv(i,n2+1,k)=.375*gra2(i,n2-2,k)-1.25*gra2(i,n2-1,k)+1.875*gra2(i,n2,k)
                delv(i,n2+1,k)=-delv(i,n2+1,k)
                !
                delw(i,n2+1,k)=.375*gra3(i,n2-2,k)-1.25*gra3(i,n2-1,k)+1.875*gra3(i,n2,k)
                delw(i,n2+1,k)=-delw(i,n2+1,k)
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
        ! face 5 avanti
        if (myid==0) then
            !
            do j=1,n2
                do i=1,n1
                    !
                    delu(i,j,0)=1.875*gra1(i,j,1)-1.25*gra1(i,j,2)+.375*gra1(i,j,3)
                    delu(i,j,0)=-delu(i,j,0)
                    !
                    delv(i,j,0)=1.875*gra2(i,j,1)-1.25*gra2(i,j,2)+.375*gra2(i,j,3)
                    delv(i,j,0)=-delv(i,j,0)
                    !
                    delw(i,j,0)=1.875*gra3(i,j,1)-1.25*gra3(i,j,2)+.375*gra3(i,j,3)
                    delw(i,j,0)=-delw(i,j,0)
                !
                end do
            end do

        end if
        !
        ! face 6 indietro
        if (myid==nproc-1) then
            !
            do j=1,n2
                do i=1,n1
                    !
                    delu(i,j,n3+1)=.375*gra1(i,j,n3-2)-1.25*gra1(i,j,n3-1)+1.875*gra1(i,j,n3)
                    delu(i,j,n3+1)=-delu(i,j,n3+1)
                    !
                    delv(i,j,n3+1)=.375*gra2(i,j,n3-2)-1.25*gra2(i,j,n3-1)+1.875*gra2(i,j,n3)
                    delv(i,j,n3+1)=-delv(i,j,n3+1)
                    !
                    delw(i,j,n3+1)=.375*gra3(i,j,n3-2)-1.25*gra3(i,j,n3-1)+1.875*gra3(i,j,n3)
                    delw(i,j,n3+1)=-delw(i,j,n3+1)
                !
                end do
            end do
        !
        end if

        return

    end subroutine condi

    subroutine update()

        ! update for intermediate velocity into the flow field
        ! at the walls parabolic extrapolation

        use scala3
        use period
        use mysending


        implicit none

        !-----------------------------------------------------------------------
        integer :: i,j,k
        !-----------------------------------------------------------------------

        ! into the field
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    u(i,j,k)=u(i,j,k)+delu(i,j,k)
                    v(i,j,k)=v(i,j,k)+delv(i,j,k)
                    w(i,j,k)=w(i,j,k)+delw(i,j,k)
                end do
            end do
        end do

        ! extrapolation left and right sides (1 and 2)
        !
        if (ip==1) then
            call extrap_from_nodes(u,u,1)
            call extrap_from_nodes(v,v,1)
            call extrap_from_nodes(w,w,1)
            call extrap_from_nodes(u,u,2)
            call extrap_from_nodes(v,v,2)
            call extrap_from_nodes(w,w,2)
        else
            do k=kparasta,kparaend
                do j=1,n2
                    u(0,j,k)=u(n1,j,k)
                    v(0,j,k)=v(n1,j,k)
                    w(0,j,k)=w(n1,j,k)
                    u(n1+1,j,k)=u(1,j,k)
                    v(n1+1,j,k)=v(1,j,k)
                    w(n1+1,j,k)=w(1,j,k)
                end do
            end do
        end if

        ! extrapolation on bottom and upper sides (3 and 4)
        ! direction 2 is always not periodic

        call extrap_from_nodes(u,u,3)
        call extrap_from_nodes(v,v,3)
        call extrap_from_nodes(w,w,3)
        call extrap_from_nodes(u,u,4)
        call extrap_from_nodes(v,v,4)
        call extrap_from_nodes(w,w,4)


        ! extrapolation on front and back sides (5 and 6)
        if (kp==1) then
            call extrap_from_nodes(u,u,5)
            call extrap_from_nodes(v,v,5)
            call extrap_from_nodes(w,w,5)
            call extrap_from_nodes(u,u,6)
            call extrap_from_nodes(v,v,6)
            call extrap_from_nodes(w,w,6)
        else

            call periodic_exchange_z(u)
            call periodic_exchange_z(v)
            call periodic_exchange_z(w)

        end if
        !
        ! every proc needs to know values in kparaend+1 to compute
        ! controvariant velocity in contra

        call border_exchange_centroids(u,1)
        call border_exchange_centroids(v,1)
        call border_exchange_centroids(w,1)

        return

    end subroutine update

    subroutine extrap_from_nodes(var_out,var_in,direc)

        ! performs parabolic interpoliation o a variable on the wall
        ! the direction of interpolation is given by the variable direc,which means:
        ! 1: positive x direction
        ! 2: negative x direction
        ! 3: positive y direction
        ! 4: negative y direction
        ! 5: positive z direction
        ! 6: negative z direction
        ! walls are assumed to be located at i=0 and i=jx+1 (same for other directions)

        implicit none

        !-----------------------------------------------------------------------
        real,intent(inout) :: var_out(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer,intent(in) :: direc
        real,intent(in) :: var_in(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        real,parameter :: c1=1.875
        real,parameter :: c2=-1.250
        real,parameter :: c3=0.375
        integer :: i,j,k
        !-----------------------------------------------------------------------

        select case (direc)
            case (1)
                do k=kparasta,kparaend
                    do j=1,n2
                        var_out(0,j,k)=c1*var_in(1,j,k)+c2*var_in(2,j,k)+c3*var_in(3,j,k)
                    end do
                end do
            case (2)
                do k=kparasta,kparaend
                    do j=1,n2
                        var_out(n1+1,j,k)=c3*var_in(n1-2,j,k)+c2*var_in(n1-1,j,k)+c1*var_in(n1,j,k)
                    end do
                end do
            case (3)
                do k=kparasta,kparaend
                    do i=1,n1
                        var_out(i,0,k)=c1*var_in(i,1,k)+c2*var_in(i,2,k)+c3*var_in(i,3,k)
                    end do
                end do
            case (4)
                do k=kparasta,kparaend
                    do i=1,n1
                        var_out(i,n2+1,k)=c3*var_in(i,n2-2,k)+c2*var_in(i,n2-1,k)+c1*var_in(i,n2,k)
                    end do
                end do
            case (5)
                if (myid==0) then
                    do j=1,n2
                        do i=1,n1
                            var_out(i,j,0)=c1*var_in(i,j,1)+c2*var_in(i,j,2)+c3*var_in(i,j,3)
                        end do
                    end do
                end if
            case (6)
                if (myid==nproc-1) then
                    do j=1,n2
                        do i=1,n1
                            var_out(i,j,n3+1)=c3*var_in(i,j,n3-2)+c2*var_in(i,j,n3-1)+c1*var_in(i,j,n3)
                        end do
                    end do
                end if
        end select
        return

    end subroutine extrap_from_nodes

    subroutine extrap_from_centroids(var_out,var_in,direc)

        ! performs parabolic interpoliation o a variable on the wall
        ! the direction of interpolation is given by the variable direc,which means:
        ! 1: positive x direction
        ! 2: negative x direction
        ! 3: positive y direction
        ! 4: negative y direction
        ! 5: positive z direction
        ! 6: negative z direction
        ! walls are assumed to be located at i=0 and i=jx+1 (same for other directions)

        implicit none

        !-----------------------------------------------------------------------
        real,intent(inout) :: var_out(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer,intent(in) :: direc
        real,intent(in) :: var_in(n1,n2,kparasta-1:kparaend+1)
        !-----------------------------------------------------------------------
        real,parameter :: c1=1.875
        real,parameter :: c2=-1.250
        real,parameter :: c3=0.375
        integer :: i,j,k
        !-----------------------------------------------------------------------

        select case (direc)
            case (1)
                do k=kparasta,kparaend
                    do j=1,n2
                        var_out(0,j,k)=c1*var_in(1,j,k)+c2*var_in(2,j,k)+c3*var_in(3,j,k)
                    end do
                end do
            case (2)
                do k=kparasta,kparaend
                    do j=1,n2
                        var_out(n1+1,j,k)=c3*var_in(n1-2,j,k)+c2*var_in(n1-1,j,k)+c1*var_in(n1,j,k)
                    end do
                end do
            case (3)
                do k=kparasta,kparaend
                    do i=1,n1
                        var_out(i,0,k)=c1*var_in(i,1,k)+c2*var_in(i,2,k)+c3*var_in(i,3,k)
                    end do
                end do
            case (4)
                do k=kparasta,kparaend
                    do i=1,n1
                        var_out(i,n2+1,k)=c3*var_in(i,n2-2,k)+c2*var_in(i,n2-1,k)+c1*var_in(i,n2,k)
                    end do
                end do
            case (5)
                if (myid==0) then
                    do j=1,n2
                        do i=1,n1
                            var_out(i,j,0)=c1*var_in(i,j,1)+c2*var_in(i,j,2)+c3*var_in(i,j,3)
                        end do
                    end do
                end if
            case (6)
                if (myid==nproc-1) then
                    do j=1,n2
                        do i=1,n1
                            var_out(i,j,n3+1)=c3*var_in(i,j,n3-2)+c2*var_in(i,j,n3-1)+c1*var_in(i,j,n3)
                        end do
                    end do
                end if
        end select
        return

    end subroutine extrap_from_centroids

end module contour_module
