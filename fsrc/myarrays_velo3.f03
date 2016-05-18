module myarrays_velo3

    ! contains arrays like velocity etc

    implicit none

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

    ! gradients
    real,allocatable :: delu(:,:,:),delv(:,:,:),delw(:,:,:)
    real,allocatable :: delrho(:,:,:),delrhov(:,:,:,:)

    ! eddy diffusivity
    real,allocatable :: akapt(:,:,:,:),akaptV(:,:,:,:)
    real,allocatable :: akapt_piano(:,:,:,:),akaptV_piano(:,:,:,:)
	
    ! array for pressure
    real,allocatable :: cs1(:,:),cs2(:,:)
    real,allocatable :: cs3(:,:),cs4(:,:)
    real,allocatable :: cs5(:,:),cs6(:,:)

    real,allocatable :: uc1_orl(:,:),uc2_orl(:,:)
    real,allocatable :: vc3_orl(:,:),vc4_orl(:,:)
    real,allocatable :: wc5_orl(:,:),wc6_orl(:,:)


contains

    subroutine iniz(f1ve,f2ve,f3ve,bcsi,beta,bzet)

        use mysending
        use scala3, only: n1,n2,n3,nscal

        use mpi

        ! matrixes allocation and initialization

        implicit none

        !-----------------------------------------------------------------------
        real,intent(out) :: f1ve(n1,n2,kparasta:kparaend)
        real,intent(out) :: f2ve(n1,n2,kparasta:kparaend)
        real,intent(out) :: f3ve(n1,n2,kparasta:kparaend)

        real,intent(out) :: bcsi(n1,n2,kparasta:kparaend)
        real,intent(out) :: beta(n1,n2,kparasta:kparaend)
        real,intent(out) :: bzet(n1,n2,kparasta:kparaend)
        !-----------------------------------------------------------------------
        integer :: i,j,k
        integer :: kpsta_alloc,kpend_alloc
        !-----------------------------------------------------------------------

        allocate(delu(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(delv(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(delw(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        allocate(delrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(delrhov(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))


        allocate(akaptV(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(akapt (nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        akaptV=0.0
        akapt=0.0
        if (myid.eq.0) then
            allocate( akapt_piano(nscal,0:n1+1,0:n2+1,n3:n3))
            allocate(akaptV_piano(nscal,0:n1+1,0:n2+1,n3:n3))
            akapt_piano  = 0.0
            akaptV_piano = 0.0
        else if (myid.eq. nproc-1) then
            allocate( akapt_piano(nscal,0:n1+1,0:n2+1,1:1))
            allocate(akaptV_piano(nscal,0:n1+1,0:n2+1,1:1))
            akapt_piano  = 0.0
            akaptV_piano = 0.0
        else
            allocate( akapt_piano(nscal,0:n1+1,0:n2+1,1:1))
            allocate(akaptV_piano(nscal,0:n1+1,0:n2+1,1:1))
            akapt_piano  = 0.0
            akaptV_piano = 0.0
        end if
        !-----------------------------------------------------------------------
        !    if (imoist==1) then
        !        allocate(tpotm(1:n2))
        !        allocate(qm(1:n2))
        !    end if
        !-----------------------------------------------------------------------
        !
        kpsta_alloc=kparasta
        kpend_alloc=kparaend
        if (myid==0) kpsta_alloc = kparasta-1
        !
        !-----------------------------------------------------------------------
        allocate(cs1(n2,n3))
        allocate(cs2(n2,n3))
        allocate(cs3(n1,n3))
        allocate(cs4(n1,n3))
        allocate(cs5(n1,n2))
        allocate(cs6(n1,n2))

        cs1=0.0
        cs2=0.0
        cs3=0.0
        cs4=0.0
        cs5=0.0
        cs6=0.0

        !-----------------------------------------------------------------------
        allocate(uc1_orl(n2,n3))
        allocate(uc2_orl(n2,n3))
        allocate(vc3_orl(n1,n3))
        allocate(vc4_orl(n1,n3))
        allocate(wc5_orl(n1,n2))
        allocate(wc6_orl(n1,n2))

        uc1_orl=0.0
        uc2_orl=0.0
        vc3_orl=0.0
        vc4_orl=0.0
        wc5_orl=0.0
        wc6_orl=0.0

        !-----------------------------------------------------------------------c
        allocate(rhs(1:n1,1:n2,kparasta:kparaend))

        allocate(u(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(v(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(w(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(fi(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(next_prs(0:n1+1,kparasta-deepl:kparaend+deepr))
        allocate(rhov(1:nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))

        allocate( uc(0:n1,1:n2,kparasta  :kparaend))
        allocate( vc(1:n1,0:n2,kparasta  :kparaend))
        allocate( wc(1:n1,1:n2,kparasta-1:kparaend))


        u = 0.0
        v = 0.0
        w = 0.0
        fi = 0.0
        next_prs=0.0
        rhov = 0.0
        uc = 0.0
        vc = 0.0
        wc = 0.0
        rhs = 0.0

        if (myid.eq.0) then
            allocate(u_piano(0:n1+1,0:n2+1,n3:n3))
            allocate(v_piano(0:n1+1,0:n2+1,n3:n3))
            allocate(w_piano(0:n1+1,0:n2+1,n3:n3))
            allocate(rhov_piano(1:nscal,0:n1+1,0:n2+1,n3:n3))

            u_piano=0.0
            v_piano=0.0
            w_piano=0.0

            allocate(wc_piano(1:n1,1:n2,n3:n3))

            wc_piano=0.0

        else if (myid.eq.nproc-1) then

            allocate(u_piano(0:n1+1,0:n2+1,1:1))
            allocate(v_piano(0:n1+1,0:n2+1,1:1))
            allocate(w_piano(0:n1+1,0:n2+1,1:1))
            allocate(rhov_piano(1:nscal,0:n1+1,0:n2+1,1:1))

            u_piano=0.0
            v_piano=0.0
            w_piano=0.0

            allocate(wc_piano(1:n1,1:n2,0:0))

            wc_piano=0.0

        end if
        !
        !-----------------------------------------------------------------------
        !
        if (myid.eq.0) then
            allocate(gra1_appoggio(n1,n2,n3:n3))
            allocate(gra2_appoggio(n1,n2,n3:n3))
            allocate(gra3_appoggio(n1,n2,n3:n3))

            do i=1,n1
                do j=1,n2
                    do k=n3,n3
                        gra1_appoggio(i,j,k)=0.
                        gra2_appoggio(i,j,k)=0.
                        gra3_appoggio(i,j,k)=0.
                    end do
                end do
            end do
        else if (myid.eq.nproc-1) then
            allocate(gra1_appoggio(n1,n2,1:1))
            allocate(gra2_appoggio(n1,n2,1:1))
            allocate(gra3_appoggio(n1,n2,1:1))

            do i=1,n1
                do j=1,n2
                    do k=1,1
                        gra1_appoggio(i,j,k)=0.
                        gra2_appoggio(i,j,k)=0.
                        gra3_appoggio(i,j,k)=0.
                    end do
                end do
            end do
        end if

        allocate(cgra1(0:n1,  n2,kparasta-1:kparaend+1))
        allocate(cgra2(  n1,0:n2,kparasta-1:kparaend+1))
        allocate(cgra3(  n1,  n2,kparasta-1:kparaend+1))

        allocate(gra1(n1,n2,kparasta-1:kparaend+1))
        allocate(gra2(n1,n2,kparasta-1:kparaend+1))
        allocate(gra3(n1,n2,kparasta-1:kparaend+1))

        do k=kparasta-1,kparaend+1
            do j=1,n2
                do i=0,n1
                    cgra1(i,j,k) = 0.0
                end do
            end do
        end do

        do k=kparasta-1,kparaend+1
            do j=0,n2
                do i=1,n1
                    cgra2(i,j,k) = 0.0
                end do
            end do
        end do

        do k=kparasta-1,kparaend+1
            do j=1,n2
                do i=1,n1
                    cgra3(i,j,k) = 0.0
                end do
            end do
        end do

        do k=kparasta-1,kparaend+1
            do j=1,n2
                do i=1,n1
                    !
                    gra1(i,j,k)=0.0
                    gra2(i,j,k)=0.0
                    gra3(i,j,k)=0.0
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    !
                    f1ve(i,j,k)=0.0
                    f2ve(i,j,k)=0.0
                    f3ve(i,j,k)=0.0
                    bcsi(i,j,k)=0.0
                    beta(i,j,k)=0.0
                    bzet(i,j,k)=0.0
                !
                end do
            end do
        end do

        ! matrix defined in the field and out
        !
        do  i=0,n1+1
            do  j=0,n2+1
                do  k=kparasta-1,kparaend+1
                    delu(i,j,k)=0.0
                    delv(i,j,k)=0.0
                    delw(i,j,k)=0.0
                end do
            end do
        end do

        !
        return
    end subroutine iniz


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
            call MPI_SEND(u(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !    quick
            if (insc==1) then
                call MPI_SEND(u(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
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
            call MPI_SEND(u(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(u(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if
        ! send v
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SEND(v(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !    quick
            if (insc==1) then
                call MPI_SEND(v(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
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
            call MPI_SEND(v(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(v(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if
        ! send w
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SEND(w(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !    quick
            if (insc==1) then
                call MPI_SEND(w(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
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
            call MPI_SEND(w(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
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
            call MPI_SEND(fi(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(fi(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
        end if
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SEND(fi(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
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
            call MPI_SEND(u(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
        else if (myid==0) then
            call MPI_RECV(u_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
        end if
        if (myid==0) then
            call MPI_SEND(u(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
        end if
        if (myid==nproc-1) then
            call MPI_RECV(u_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
        end if

        if (myid==nproc-1) then
            call MPI_SEND(v(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
        else if (myid==0) then
            call MPI_RECV(v_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
        end if
        if (myid==0) then
            call MPI_SEND(v(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
        end if
        if (myid==nproc-1) then
            call MPI_RECV(v_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
        end if

        if (myid==nproc-1) then
            call MPI_SEND(w(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
        else if (myid==0) then
            call MPI_RECV(w_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
        end if
        if (myid==0) then
            call MPI_SEND(w(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
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
        integer :: ierr,status(MPI_STATUS_SIZE),req


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
        if (kp==0) then
            if (myid==nproc-1) then
                call MPI_SEND(wc(1,1,jz),jx*jy,MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
            else if (myid==0) then
                call MPI_RECV(wc_piano(1,1,jz),jx*jy,MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
            end if
            if (myid==0) then
                do i=1,jx
                    do j=1,jy
                        wc(i,j,0)=.5*(wc(i,j,0)+wc_piano(i,j,jz))
                    end do
                end do
                call MPI_SEND(wc(1,1,0),jx*jy,MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(wc(1,1,jz),jx*jy,MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
            end if
        end if


    end subroutine average_periodicfluxes

end module myarrays_velo3





