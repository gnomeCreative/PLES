module myarrays_velo3

    ! contains arrays like velocity etc

    implicit none
    !-----------------------------------------------------------------------
    real,allocatable :: u(:,:,:),v(:,:,:),w(:,:,:),rhov(:,:,:,:)
    real,allocatable :: uc(:,:,:),vc(:,:,:),wc(:,:,:)
    real,allocatable :: rhs(:,:,:),fi(:,:,:),next_prs(:,:)

    real,allocatable ::  u_piano(:,:,:),v_piano(:,:,:),w_piano(:,:,:)
    real,allocatable ::  rhov_piano(:,:,:,:)
    real,allocatable ::  uc_piano(:,:,:)
    real,allocatable ::  vc_piano(:,:,:)
    real,allocatable ::  wc_piano(:,:,:)

    real,allocatable ::  cgra1(:,:,:),cgra2(:,:,:),cgra3(:,:,:)
    real,allocatable ::  gra1(:,:,:), gra2(:,:,:), gra3(:,:,:)

    real,allocatable ::  gra1_appoggio(:,:,:)
    real,allocatable ::  gra2_appoggio(:,:,:)
    real,allocatable ::  gra3_appoggio(:,:,:)

    real,allocatable :: delu(:,:,:),delv(:,:,:),delw(:,:,:)
	
    ! array for pressure
    real,allocatable :: cs1(:,:),cs2(:,:)
    real,allocatable :: cs3(:,:),cs4(:,:)
    real,allocatable :: cs5(:,:),cs6(:,:)

    real,allocatable :: uc1_orl(:,:),uc2_orl(:,:)
    real,allocatable :: vc3_orl(:,:),vc4_orl(:,:)
    real,allocatable :: wc5_orl(:,:),wc6_orl(:,:)

contains

    subroutine communication_velocity()

        use mysettings, only: insc
        use mysending
        use scala3, only: jx,jy,jz

        use mpi

        implicit none

        integer :: ierr
        integer :: status(MPI_STATUS_SIZE)

        ! send u
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(u(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !    quick
            if (insc==1) then
                call MPI_SSEND(u(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            end if
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(u(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            !    quick
            if (insc==1) then
                call MPI_RECV(u(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            end if
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(u(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(u(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if
        ! send v
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(v(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !    quick
            if (insc==1) then
                call MPI_SSEND(v(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            end if
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(v(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            !    quick
            if (insc==1) then
                call MPI_RECV(v(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            end if
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(v(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(v(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if
        ! send w
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(w(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !    quick
            if (insc==1) then
                call MPI_SSEND(w(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            end if
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(w(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            !    quick
            if (insc==1) then
                call MPI_RECV(w(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            end if
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(w(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(w(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if

    end subroutine communication_velocity

    subroutine communication_pressure()

        use mysending
        use scala3, only: jx,jy,jz

        use mpi

        implicit none

        integer :: ierr
        integer :: status(MPI_STATUS_SIZE)

        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(fi(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(fi(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(fi(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(fi(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if

    end subroutine communication_pressure

    subroutine communication_velpiano()

        use mysending
        use scala3, only: jx,jy,jz

        use mpi

        implicit none

        integer :: ierr,status(MPI_STATUS_SIZE)

        if (myid==nproc-1) then
            call MPI_SSEND(u(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
        else if (myid==0) then
            call MPI_RECV(u_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
        end if
        if (myid==0) then
            call MPI_SSEND(u(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
        end if
        if (myid==nproc-1) then
            call MPI_RECV(u_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
        end if

        if (myid==nproc-1) then
            call MPI_SSEND(v(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
        else if (myid==0) then
            call MPI_RECV(v_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
        end if
        if (myid==0) then
            call MPI_SSEND(v(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
        end if
        if (myid==nproc-1) then
            call MPI_RECV(v_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
        end if

        if (myid==nproc-1) then
            call MPI_SSEND(w(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
        else if (myid==0) then
            call MPI_RECV(w_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
        end if
        if (myid==0) then
            call MPI_SSEND(w(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
        end if
        if (myid==nproc-1) then
            call MPI_RECV(w_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
        end if

    end subroutine communication_velpiano

    subroutine average_periodicfluxes()

        use mysending
        use scala3, only: jx,jy,jz
        use period

        use mpi

        implicit none

        integer i,j,k,ii,jj,kk
        integer :: ierr,status(MPI_STATUS_SIZE)


        ! average on fluxes in periodicity direction
        do ii=1,1-ip
            do k=kparasta,kparaend
                do j=1,jy
                    uc(0,j,k)=.5*(uc(0,j,k)+uc(jx,j,k))
                    uc(jx,j,k)=uc(0,j,k)
                end do
            end do
        end do

        do jj=1,1-jp
            do k=kparasta,kparaend
                do i=1,jx
                    vc(i,0,k)=.5*(vc(i,0,k)+vc(i,jy,k))
                    vc(i,jy,k)=vc(i,0,k)
                end do
            end do
        end do

        ! send wc(i,j,jz) to P0 in wc_piano of myid=0
        do kk=1,1-kp
            if (myid==nproc-1) then
                call MPI_SSEND(wc(1,1,jz),jx*jy,MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
            else if (myid==0) then
                call MPI_RECV(wc_piano(1,1,jz),jx*jy,MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
            end if
            if (myid==0) then
                do i=1,jx
                    do j=1,jy
                        wc(i,j,0)=.5*(wc(i,j,0)+wc_piano(i,j,jz))
                    end do
                end do
                call MPI_SSEND(wc(1,1,0),jx*jy,MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(wc(1,1,jz),jx*jy,MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
            end if
        end do

    end subroutine average_periodicfluxes
end module myarrays_velo3





