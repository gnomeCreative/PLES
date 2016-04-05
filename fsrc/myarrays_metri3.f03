!***********************************************************************
module myarrays_metri3
    !***********************************************************************
    ! coordinates x,y,z
    ! 19 metric terms, 6 for each face plus the giacobian
    ! fluxes
    ! eddy viscosity
    ! NOTE better to move annit and fluxes in an other module

    implicit none
    !-----------------------------------------------------------------------
    real,allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
    real,allocatable :: csx(:,:,:),csy(:,:,:),csz(:,:,:)
    real,allocatable :: etx(:,:,:),ety(:,:,:),etz(:,:,:)
    real,allocatable :: ztx(:,:,:),zty(:,:,:),ztz(:,:,:)

    real,allocatable :: xcd(:,:,:),ycd(:,:,:),zcd(:,:,:)

    real,allocatable :: g11(:,:,:),g12(:,:,:),g13(:,:,:)
    real,allocatable :: g21(:,:,:),g22(:,:,:),g23(:,:,:)
    real,allocatable :: g31(:,:,:),g32(:,:,:),g33(:,:,:)
      
    real,allocatable :: giac(:,:,:)

    real,allocatable :: f1(:,:,:),f2(:,:,:),f3(:,:,:)
    real,allocatable :: annit(:,:,:),annitV(:,:,:)
    real,allocatable :: annit_piano(:,:,:),annitV_piano(:,:,:)

    !CORIOLIS OIL
    real,allocatable :: g_co11(:,:),g_co12(:,:),g_co13(:,:)
    real,allocatable :: g_co31(:,:),g_co32(:,:),g_co33(:,:)

    real :: box_vol,solid_vol,fluid_vol

contains

subroutine init_metrica(nlevmultimax,nlevel)
    ! initialize metrica

    use mysending
    use scala3

    implicit none

    !-----------------------------------------------------------------------
    !     arrays declaration

    integer,intent(inout) :: nlevel
    integer,intent(out) :: nlevmultimax

    integer i,j,k
    integer var_piani,livello

    !-----------------------------------------------------------------------
    ! griglia

    !    qui decido la profondita' di allocazione della griglia
    !    per ogni processore in funzione della subroutine mul_met
    !    questo perche' in relazione al nlevel di multigrid la
    !    allocazione locale di griglia deve essere maggiore
    var_piani=ncolperproc
    livello=1
    do while(var_piani.ne.1.and.livello.ne.5)
        var_piani=int(var_piani/2)
        livello=livello+1
    end do
    livello=livello-1
    !      write(*,*)myid,'multigrid level for each proc',livello
    if(livello.ne.1)then
        deep_mul=2*2**(livello-2) ! allocation deep
    !        write(*,*)'grid allocation deep',deep_mul
    else
        if(myid==0)then
            write(*,*)'PROBLEM on multigrid level for proc'
        end if
        stop
    end if

    if(livello .lt. nlevel)then
        if(myid==0)then
            write(*,*)'PROBLEM on multigrid level for proc'
            write(*,*)'level reduced from ',nlevel,' to ',livello
        end if
        nlevel = livello
        nlevmultimax = nlevel
    end if

    !-----------------------------------------------------------------------
    !     coordinates allocation
    if(myid .eq. 0)then
        allocate (x(-8:n1+8,-8:n2+8,-8:kparaend+deep_mul+1))
        allocate (y(-8:n1+8,-8:n2+8,-8:kparaend+deep_mul+1))
        allocate (z(-8:n1+8,-8:n2+8,-8:kparaend+deep_mul+1))
    elseif(myid.eq.nproc-1)then
        allocate (x(-8:n1+8,-8:n2+8,kparasta-deep_mul-1:kparaend+8))
        allocate (y(-8:n1+8,-8:n2+8,kparasta-deep_mul-1:kparaend+8))
        allocate (z(-8:n1+8,-8:n2+8,kparasta-deep_mul-1:kparaend+8))
    else
        allocate(x(-8:n1+8,-8:n2+8, &
            kparasta-deep_mul-1:kparaend+deep_mul+1))
        allocate(y(-8:n1+8,-8:n2+8, &
            kparasta-deep_mul-1:kparaend+deep_mul+1))
        allocate(z(-8:n1+8,-8:n2+8, &
            kparasta-deep_mul-1:kparaend+deep_mul+1))
    end if

    x = 0.
    y = 0.
    z = 0.

    if(myid.eq.0)write(*,*)myid,'coordinates allocation done'
    !-----------------------------------------------------------------------
    !     metric terms allocation


    allocate (csx(0:n1,1:n2,kparasta:kparaend))
    allocate (csy(0:n1,1:n2,kparasta:kparaend))
    allocate (csz(0:n1,1:n2,kparasta:kparaend))
    allocate (g11(0:n1,1:n2,kparasta:kparaend))
    allocate (g12(0:n1,1:n2,kparasta:kparaend))
    allocate (g13(0:n1,1:n2,kparasta:kparaend))

    allocate (etx(1:n1,0:n2,kparasta:kparaend))
    allocate (ety(1:n1,0:n2,kparasta:kparaend))
    allocate (etz(1:n1,0:n2,kparasta:kparaend))
    allocate (g21(1:n1,0:n2,kparasta:kparaend))
    allocate (g22(1:n1,0:n2,kparasta:kparaend))
    allocate (g23(1:n1,0:n2,kparasta:kparaend))

    allocate (ztx(1:n1,1:n2,kparasta-1:kparaend)) !NOTA: ksta per myid=0 vale 0
    allocate (zty(1:n1,1:n2,kparasta-1:kparaend))
    allocate (ztz(1:n1,1:n2,kparasta-1:kparaend))
    allocate (g31(1:n1,1:n2,kparasta-1:kparaend))
    allocate (g32(1:n1,1:n2,kparasta-1:kparaend))
    allocate (g33(1:n1,1:n2,kparasta-1:kparaend))

    allocate (giac(1:n1,1:n2,kparasta:kparaend))

    csx = 0.
    csy = 0.
    csz = 0.

    etx = 0.
    ety = 0.
    etz = 0.

    ztx = 0.
    zty = 0.
    ztz = 0.

    g11 = 0.
    g12 = 0.
    g13 = 0.

    g21 = 0.
    g22 = 0.
    g23 = 0.

    g31 = 0.
    g32 = 0.
    g33 = 0.

    giac = 0.

    !CORIOLISOIL
    allocate (g_co11(0:n1,kparasta:kparaend))
    allocate (g_co12(0:n1,kparasta:kparaend))
    allocate (g_co13(0:n1,kparasta:kparaend))

    allocate (g_co31(1:n1,kparasta-1:kparaend))
    allocate (g_co32(1:n1,kparasta-1:kparaend))
    allocate (g_co33(1:n1,kparasta-1:kparaend))

    g_co11 = 0.
    g_co12 = 0.
    g_co13 = 0.

    g_co31 = 0.
    g_co32 = 0.
    g_co33 = 0.

    if(myid.eq.0)write(*,*)myid,'metric allocation done'
    !-----------------------------------------------------------------------
    !     fluxes allocation

    allocate (f1(0:n1,1:n2,kparasta  :kparaend))
    allocate (f2(1:n1,0:n2,kparasta  :kparaend))
    allocate (f3(1:n1,1:n2,kparasta-1:kparaend))

    f1 = 0.
    f2 = 0.
    f3 = 0.

    if(myid.eq.0)write(*,*)myid,'fluxes allocation done'
    !-----------------------------------------------------------------------

    ! viscosita' turbolenta

    allocate ( annit(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
    allocate (annitV(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))

    annit  = 0.
    annitV = 0.

    if(myid.eq.0)then
        allocate( annit_piano(0:n1+1,0:n2+1,n3:n3))
        allocate(annitV_piano(0:n1+1,0:n2+1,n3:n3))
        annit_piano  = 0.
        annitV_piano = 0.
    elseif(myid.eq. nproc-1)then
        allocate( annit_piano(0:n1+1,0:n2+1,1:1))
        allocate(annitV_piano(0:n1+1,0:n2+1,1:1))
        annit_piano  = 0.
        annitV_piano = 0.
    end if

    if(myid.eq.0)write(*,*)myid,'eddy viscosity allocation done'

    return

end subroutine init_metrica

    subroutine communication_viscosity()

        use mysending
        use scala3
        use mpi

        implicit none

        integer :: ierror
        integer :: status(MPI_STATUS_SIZE)

        ! send to left
        if (myid/=0) then
            call MPI_SSEND(annit(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(annit(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if
        ! send to right
        if (myid/=nproc-1) then
            call MPI_SSEND(annit(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(annit(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if
        !-----------------------------------------------------------------------
        ! send to left
        if (myid/=0) then
            call MPI_SSEND(annitV(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(annitV(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if
        ! send to right
        if (myid/=nproc-1) then
            call MPI_SSEND(annitV(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(annitV(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

    end subroutine communication_viscosity

    subroutine compute_centroids() !jx,jy,jz,myid,nproc,kparasta,kparaend

            use mysending
    use scala3
    use mpi
        implicit none

        !integer,intent(in) :: jx,jy,jz
        !integer,intent(in) :: myid,nproc,kparasta,kparaend
        !-----------------------------------------------------------------------
        integer i,j,k
        integer ksta,kend
        integer iprintgrid,ighost,ibound

        integer req1,req2,req3,req4,req5,req6
        integer ierr,istatus,status(MPI_STATUS_SIZE)
        !-----------------------------------------------------------------------
        iprintgrid =1
        ighost = 1

        if(myid==0)then
            ksta = kparasta
            kend = kparaend + 1
        elseif(myid==nproc-1)then
            ksta = kparasta - 1
            kend = kparaend
        else
            ksta = kparasta - 1
            kend = kparaend + 1
        end if

        ! compute centroids

        allocate(xcd(0:jx+1,0:jy+1,ksta-1:kend+1)) ! fix the periodicity!!!!!!
        allocate(ycd(0:jx+1,0:jy+1,ksta-1:kend+1))
        allocate(zcd(0:jx+1,0:jy+1,ksta-1:kend+1))


        !-----------------------------------------------------------------------
        ! centroids
        do k=ksta,kend  !1,jz
            do j=1,jy
                do i=1,jx

                    xcd(i,j,k)=0.125*(x(i,j,k)+x(i,j,k-1)+x(i-1,j,k-1)+x(i-1,j,k) &
                        +x(i,j-1,k)+x(i,j-1,k-1)+x(i-1,j-1,k-1)+x(i-1,j-1,k))

                    ycd(i,j,k)=0.125*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k) &
                        +y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k))

                    zcd(i,j,k)=0.125*(z(i,j,k)+z(i,j,k-1)+z(i-1,j,k-1)+z(i-1,j,k) &
                        +z(i,j-1,k)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j-1,k))

                end do
            end do
        end do

        ! side 1 and 2
        do k=ksta,kend  !1,jz
            do j=1,jy

                i=0

                xcd(i,j,k)=0.25*(x(i,j,k)+x(i,j,k-1)+x(i,j-1,k-1)+x(i,j-1,k))
                ycd(i,j,k)=0.25*(y(i,j,k)+y(i,j,k-1)+y(i,j-1,k-1)+y(i,j-1,k))
                zcd(i,j,k)=0.25*(z(i,j,k)+z(i,j,k-1)+z(i,j-1,k-1)+z(i,j-1,k))

                i=jx

                xcd(i+1,j,k)=0.25*(x(i,j,k)+x(i,j,k-1)+x(i,j-1,k-1)+x(i,j-1,k))
                ycd(i+1,j,k)=0.25*(y(i,j,k)+y(i,j,k-1)+y(i,j-1,k-1)+y(i,j-1,k))
                zcd(i+1,j,k)=0.25*(z(i,j,k)+z(i,j,k-1)+z(i,j-1,k-1)+z(i,j-1,k))

            end do
        end do

        do i=1,jx
            j=0
            k=ksta-1 !0
            xcd(i,j,k)=0.5*(x(i,j,k)+x(i-1,j,k))
            ycd(i,j,k)=0.5*(y(i,j,k)+y(i-1,j,k))
            zcd(i,j,k)=0.5*(z(i,j,k)+z(i-1,j,k))

            j=0
            k=kend !jz
            xcd(i,j,k+1)=0.5*(x(i,j,k)+x(i-1,j,k))
            ycd(i,j,k+1)=0.5*(y(i,j,k)+y(i-1,j,k))
            zcd(i,j,k+1)=0.5*(z(i,j,k)+z(i-1,j,k))

            j=jy
            k=kend !jz
            xcd(i,j+1,k+1)=0.5*(x(i,j,k)+x(i-1,j,k))
            ycd(i,j+1,k+1)=0.5*(y(i,j,k)+y(i-1,j,k))
            zcd(i,j+1,k+1)=0.5*(z(i,j,k)+z(i-1,j,k))

            j=jy
            k=ksta-1 !0
            xcd(i,j+1,k)=0.5*(x(i,j,k)+x(i-1,j,k))
            ycd(i,j+1,k)=0.5*(y(i,j,k)+y(i-1,j,k))
            zcd(i,j+1,k)=0.5*(z(i,j,k)+z(i-1,j,k))
        end do
        !...........................................................
        ! side 3 and 4 plus corners
        do k=ksta,kend !1,jz
            do i=1,jx

                j=0

                xcd(i,j,k)=0.25*(x(i,j,k)+x(i,j,k-1)+x(i-1,j,k-1)+x(i-1,j,k))
                ycd(i,j,k)=0.25*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k))
                zcd(i,j,k)=0.25*(z(i,j,k)+z(i,j,k-1)+z(i-1,j,k-1)+z(i-1,j,k))

                j=jy

                xcd(i,j+1,k)=0.25*(x(i,j,k)+x(i,j,k-1)+x(i-1,j,k-1)+x(i-1,j,k))
                ycd(i,j+1,k)=0.25*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k))
                zcd(i,j+1,k)=0.25*(z(i,j,k)+z(i,j,k-1)+z(i-1,j,k-1)+z(i-1,j,k))

            end do
        end do

        do j=1,jy

            k=ksta-1 !0
            i=0
            xcd(i,j,k)=0.5*(x(i,j,k)+x(i,j-1,k))
            ycd(i,j,k)=0.5*(y(i,j,k)+y(i,j-1,k))
            zcd(i,j,k)=0.5*(z(i,j,k)+z(i,j-1,k))

            k=ksta-1 !0
            i=jx
            xcd(i+1,j,k)=0.5*(x(i,j,k)+x(i,j-1,k))
            ycd(i+1,j,k)=0.5*(y(i,j,k)+y(i,j-1,k))
            zcd(i+1,j,k)=0.5*(z(i,j,k)+z(i,j-1,k))

            k=kend !jz
            i=jx
            xcd(i+1,j,k+1)=0.5*(x(i,j,k)+x(i,j-1,k))
            ycd(i+1,j,k+1)=0.5*(y(i,j,k)+y(i,j-1,k))
            zcd(i+1,j,k+1)=0.5*(z(i,j,k)+z(i,j-1,k))

            k=kend !jz
            i=0
            xcd(i,j,k+1)=0.5*(x(i,j,k)+x(i,j-1,k))
            ycd(i,j,k+1)=0.5*(y(i,j,k)+y(i,j-1,k))
            zcd(i,j,k+1)=0.5*(z(i,j,k)+z(i,j-1,k))

        end do
        !..........................................................
        ! parete 5 e 6 e spigoli
        do j=1,jy
            do i=1,jx

                k=ksta-1 !0

                xcd(i,j,k)=0.25*(x(i,j,k)+x(i,j-1,k)+x(i-1,j-1,k)+x(i-1,j,k))
                ycd(i,j,k)=0.25*(y(i,j,k)+y(i,j-1,k)+y(i-1,j-1,k)+y(i-1,j,k))
                zcd(i,j,k)=0.25*(z(i,j,k)+z(i,j-1,k)+z(i-1,j-1,k)+z(i-1,j,k))

                k=kend !jz

                xcd(i,j,k+1)=0.25*(x(i,j,k)+x(i,j-1,k)+x(i-1,j-1,k)+x(i-1,j,k))
                ycd(i,j,k+1)=0.25*(y(i,j,k)+y(i,j-1,k)+y(i-1,j-1,k)+y(i-1,j,k))
                zcd(i,j,k+1)=0.25*(z(i,j,k)+z(i,j-1,k)+z(i-1,j-1,k)+z(i-1,j,k))

            end do
        end do

        do k=ksta,kend !1,jz
            j=0
            i=0
            xcd(i,j,k)=0.5*(x(i,j,k)+x(i,j,k-1))
            ycd(i,j,k)=0.5*(y(i,j,k)+y(i,j,k-1))
            zcd(i,j,k)=0.5*(z(i,j,k)+z(i,j,k-1))

            j=0
            i=jx
            xcd(i+1,j,k)=0.5*(x(i,j,k)+x(i,j,k-1))
            ycd(i+1,j,k)=0.5*(y(i,j,k)+y(i,j,k-1))
            zcd(i+1,j,k)=0.5*(z(i,j,k)+z(i,j,k-1))

            j=jy
            i=jx
            xcd(i+1,j+1,k)=0.5*(x(i,j,k)+x(i,j,k-1))
            ycd(i+1,j+1,k)=0.5*(y(i,j,k)+y(i,j,k-1))
            zcd(i+1,j+1,k)=0.5*(z(i,j,k)+z(i,j,k-1))

            j=jy
            i=0
            xcd(i,j+1,k)=0.5*(x(i,j,k)+x(i,j,k-1))
            ycd(i,j+1,k)=0.5*(y(i,j,k)+y(i,j,k-1))
            zcd(i,j+1,k)=0.5*(z(i,j,k)+z(i,j,k-1))
        end do


        i=0
        j=0
        k=ksta-1 !0
        xcd(i,j,k)=x(i,j,k)
        ycd(i,j,k)=y(i,j,k)
        zcd(i,j,k)=z(i,j,k)

        i=0
        j=0
        k=kend !jz
        xcd(i,j,k+1)=x(i,j,k)
        ycd(i,j,k+1)=y(i,j,k)
        zcd(i,j,k+1)=z(i,j,k)

        i=0
        j=jy
        k=kend !jz
        xcd(i,j+1,k+1)=x(i,j,k)
        ycd(i,j+1,k+1)=y(i,j,k)
        zcd(i,j+1,k+1)=z(i,j,k)

        i=0
        j=jy
        k=ksta-1 !0
        xcd(i,j+1,k)=x(i,j,k)
        ycd(i,j+1,k)=y(i,j,k)
        zcd(i,j+1,k)=z(i,j,k)
        !.....................
        i=jx
        j=0
        k=ksta-1 !0
        xcd(i+1,j,k)=x(i,j,k)
        ycd(i+1,j,k)=y(i,j,k)
        zcd(i+1,j,k)=z(i,j,k)

        i=jx
        j=0
        k=kend !jz

        xcd(i+1,j,k+1)=x(i,j,k)
        ycd(i+1,j,k+1)=y(i,j,k)
        zcd(i+1,j,k+1)=z(i,j,k)

        i=jx
        j=jy
        k=kend !jz
        xcd(i+1,j+1,k+1)=x(i,j,k)
        ycd(i+1,j+1,k+1)=y(i,j,k)
        zcd(i+1,j+1,k+1)=z(i,j,k)

        i=jx
        j=jy
        k=ksta-1 !0
        xcd(i+1,j+1,k)=x(i,j,k)
        ycd(i+1,j+1,k)=y(i,j,k)
        zcd(i+1,j+1,k)=z(i,j,k)


        !     if(iprintgrid == 1)then
        !         write(450+myid,*)'zone f=point i=', &
        !             jx+2,' j=',jy+2,'k=',kend-ksta+3
        !         do k=ksta-1,kend+1
        !             do j=0,jy+1
        !                 do i=0,jx+1
        !                     write(450+myid,*)xcd(i,j,k),ycd(i,j,k),zcd(i,j,k)
        !                 end do
        !             end do
        !         end do
        !     end if
        !
        !-----------------------------------------------------------------------
        !     construct ghost cells

        if(ighost==1)then
            do k=ksta-1,kend+1
                do j=0,jy+1
                    !        side 1
                    xcd(0,j,k) = 2.*xcd(0,j,k)-xcd(1,j,k)
                    ycd(0,j,k) = 2.*ycd(0,j,k)-ycd(1,j,k)
                    zcd(0,j,k) = 2.*zcd(0,j,k)-zcd(1,j,k)

                    !        side 2
                    xcd(jx+1,j,k) = 2.*xcd(jx+1,j,k)-xcd(jx,j,k)
                    ycd(jx+1,j,k) = 2.*ycd(jx+1,j,k)-ycd(jx,j,k)
                    zcd(jx+1,j,k) = 2.*zcd(jx+1,j,k)-zcd(jx,j,k)
                end do
            end do


            do k=ksta-1,kend+1 !parasta,kparaend
                do i=0,jx+1
                    !        side 3
                    xcd(i,0,k) = 2.*xcd(i,0,k)-xcd(i,1,k)
                    ycd(i,0,k) = 2.*ycd(i,0,k)-ycd(i,1,k)
                    zcd(i,0,k) = 2.*zcd(i,0,k)-zcd(i,1,k)
                    !        side 4
                    xcd(i,jy+1,k) = 2.*xcd(i,jy+1,k)-xcd(i,jy,k)
                    ycd(i,jy+1,k) = 2.*ycd(i,jy+1,k)-ycd(i,jy,k)
                    zcd(i,jy+1,k) = 2.*zcd(i,jy+1,k)-zcd(i,jy,k)
                end do
            end do

            if(myid==0)then
                do j=0,jy+1
                    do i=0,jx+1
                        !        side 5
                        xcd(i,j,0) = 2.*xcd(i,j,0)-xcd(i,j,1)
                        ycd(i,j,0) = 2.*ycd(i,j,0)-ycd(i,j,1)
                        zcd(i,j,0) = 2.*zcd(i,j,0)-zcd(i,j,1)
                    end do
                end do
            end if

            if(myid==nproc-1)then
                do j=0,jy+1
                    do i=0,jx+1
                        !        side 6
                        xcd(i,j,jz+1) = 2.*xcd(i,j,jz+1)-xcd(i,j,jz)
                        ycd(i,j,jz+1) = 2.*ycd(i,j,jz+1)-ycd(i,j,jz)
                        zcd(i,j,jz+1) = 2.*zcd(i,j,jz+1)-zcd(i,j,jz)
                    end do
                end do
            end if


        !         if(iprintgrid == 1)then
        !             write(550+myid,*)'zone f=point i=', &
        !                 jx+2,' j=',jy+2,'k=',kend-ksta+3
        !             do k=ksta-1,kend+1
        !                 do j=0,jy+1
        !                     do i=0,jx+1
        !                         write(550+myid,*)xcd(i,j,k),ycd(i,j,k),zcd(i,j,k),0,0,0,0,0,0,0
        !                     end do
        !                 end do
        !             end do
        !         end if

        end if
        !-----------------------------------------------------------------------
        !     for periodicity I need to know the grid length in z
        !     if(myid==0 .or. myid==nproc-1)then
        !         allocate(x_bound(0:jx,0:jy))
        !         allocate(y_bound(0:jx,0:jy))
        !         allocate(z_bound(0:jx,0:jy))
        !         x_bound = 0.
        !         y_bound = 0.
        !         z_bound = 0.
        !     end if
        !     ibound = 1
        !     if(ibound == 1)then
        !         if(myid.eq.nproc-1)then
        !             call MPI_SSEND(x(0,0,jz),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 0,1001,MPI_COMM_WORLD,ierr)
        !             !         call MPI_WAIT(req1,istatus,ierr)
        !
        !             call MPI_SSEND(y(0,0,jz),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 0,1002,MPI_COMM_WORLD,ierr)
        !             !         call MPI_WAIT(req2,istatus,ierr)
        !
        !             call MPI_SSEND(z(0,0,jz),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 0,1003,MPI_COMM_WORLD,ierr)
        !         !         call MPI_WAIT(req3,istatus,ierr)
        !         elseif(myid.eq.0)then
        !             call MPI_RECV(x_bound(0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 nproc-1,1001,MPI_COMM_WORLD,status,ierr)
        !             !         call MPI_WAIT(req4,istatus,ierr)
        !
        !             call MPI_RECV(y_bound(0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 nproc-1,1002,MPI_COMM_WORLD,status,ierr)
        !             !         call MPI_WAIT(req5,istatus,ierr)
        !
        !             call MPI_RECV(z_bound(0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 nproc-1,1003,MPI_COMM_WORLD,status,ierr)
        !         !         call MPI_WAIT(req6,istatus,ierr)
        !         endif
        !
        !
        !
        !         if(myid.eq.0)then
        !             call MPI_SSEND(x(0,0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 nproc-1,2001,MPI_COMM_WORLD,ierr)
        !             !         call MPI_WAIT(req1,istatus,ierr)
        !
        !             call MPI_SSEND(y(0,0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 nproc-1,2002,MPI_COMM_WORLD,ierr)
        !             !         call MPI_WAIT(req2,istatus,ierr)
        !
        !             call MPI_SSEND(z(0,0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 nproc-1,2003,MPI_COMM_WORLD,ierr)
        !         !         call MPI_WAIT(req3,istatus,ierr)
        !         endif
        !         if(myid.eq.nproc-1)then
        !             call MPI_RECV(x_bound(0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 0,2001,MPI_COMM_WORLD,status,ierr)
        !             !         call MPI_WAIT(req4,istatus,ierr)
        !
        !             call MPI_RECV(y_bound(0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 0,2002,MPI_COMM_WORLD,status,ierr)
        !             !         call MPI_WAIT(req5,istatus,ierr)
        !
        !             call MPI_RECV(z_bound(0,0),(jx+1)*(jy+1),MPI_REAL_SD, &
        !                 0,2003,MPI_COMM_WORLD,status,ierr)
        !         !         call MPI_WAIT(req6,istatus,ierr)
        !         endif
        !     end if
        !-----------------------------------------------------------------------
        return
    end
                       
    subroutine compute_volume(tipo) !n1,n2,deepl,deepr,myid,kparasta,kparaend

    use mysending
    use scala3
    use mpi

        implicit none

        !integer,intent(in) :: n1,n2,deepl,deepr,myid,kparasta,kparaend
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        ! for volume coputation
        integer :: i,j,k
        integer :: ierr

        real :: vol_loc,fluid_loc,solid_loc

            !-----------------------------------------------------------------------
        !     compute the volume once
        vol_loc=0.
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    vol_loc=vol_loc+giac(i,j,k)
                    if (tipo(i,j,k)/=0) then
                        fluid_loc=fluid_loc+giac(i,j,k)
                    end if
                    if (tipo(i,j,k)==0) then
                        solid_loc=solid_loc+giac(i,j,k)
                    end if
                end do
            end do
        end do
        !
        !     reduce operation on vol
        !
        box_vol=0.0
        call MPI_REDUCE(vol_loc,box_vol,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        solid_vol=0.0
        call MPI_REDUCE(solid_loc,solid_vol,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        fluid_vol=0.0
        call MPI_REDUCE(fluid_loc,fluid_vol,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (myid == 0) then
            write(*,*)'box volume: ',box_vol
            write(*,*)'(fluid: ',fluid_vol,'; solid: ',solid_vol,')'
        end if

    end subroutine compute_volume


!***********************************************************************
end module myarrays_metri3
!***********************************************************************
