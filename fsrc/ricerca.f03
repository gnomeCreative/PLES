module ricerca_module
    !-----------------------------------------------------------------------------------------------
    !***********************************************************************************************
    !-----------------------------------------------------------------------------------------------
    ! release 21-11-2009
    ! Federico Roman
    ! modified by Alessandro Leonardi, starting December 2015

    ! Pre-processing program to use the Immersed Boundary Method as in
    ! Roman,Napoli,Armenio,Milici 2008

    ! GRID: non staggered
    ! IMMERSED BODY: must be discretize with a triangular mesh,
    ! besides the mesh must belong to a close volume.

    ! INPUT: "gri3dp_in.dat"                grid
    !        "Coordinate_IMMB.inp"              mesh point coordinates
    !    "Indici_piani_IMMB.inp"            mesh triangle indices
    !        "Coordinate//n body//_IMMB.inp"        mesh point coordinatess multi body
    !    "Indici//n body//_piani_IMMB.inp"      mesh triangle indices multi body
    !    "tipologia.dat"                phase type if input=1
    !    "tipologia//n body//.dat"          phase type if input=2
    !    "restart_solide//100+iteration//.dat"      restart the phase search, restrat = 1

    ! OUTPUT: "Celle_IB_indici.inp"             IB and V indices i,j,k              NS solver input
    !     "Celle_IB_distanze.inp"           PP-V distances in x y z                 "
    !     "Celle_Bloccate_indici.inp"           solid nodes indices i j k               "
    !     "distanze_interpolazioni.inp"         PP-IB and IB-IP distances, IP coordinates       "
    !
    !     "nodi_solidi.dat"             solid node coordinates              to visualize the results
    !     "nodi_parete.dat"             IP         coordinates                         "
    !     "nodi_ib.dat"                 IB    node coordinates                         "
    !     "nodi_proiezione.dat"             PP         coordinates                         "
    !     "nodi_vicini.dat"             V     node coordinates                         "


    ! CHECK IP
    !    fort.1000 IP on plane    ---> ok
    !    fort.1100 IP on plane with iterative box procedure  ---> ok
    !    fort.1200 IP on line     ---> ok, but maybe calibration is needed, increase calibro_iter
    !                           it will increase fort.1100 and fort.1200 will be reduce
    !    fort.1300 IP on vertex     ---> this is not good but can happen
    !    fort.1400 IP no found   ---> there is something wrong
    !
    ! CHECK PP
    !
    !-----------------------------------------------------------------------
    !***********************************************************************
    !-----------------------------------------------------------------------
    use :: mysending
    use :: trilinear
    use :: scala3
    use :: myarrays_metri3
    use :: myarrays_ibm
    use :: mpi


    use,intrinsic :: iso_c_binding

    implicit none


contains

    !---------------------------------------------------------------------------------------------------------------------------------------------
    !*********************************************************************************************************************************************
    !---------------------------------------------------------------------------------------------------------------------------------------------

    subroutine do_ricerca(tipo)

        implicit none

        integer,intent(out) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        ! number of computational nodes in x y z
        !integer :: n1,n2,n3
        integer :: iinizio, ifine
        integer :: jinizio, jfine
        integer :: kinizio, kfine

        integer :: ierr

        real(8) :: time1,time2,time

        !---------------------------------------------------------------------------------------------------------------------------------------------
        !*********************************************************************************************************************************************
        !---------------------------------------------------------------------------------------------------------------------------------------------
        ! CHECK INPUT
        !        if(miro_solide == 0)then
        !            iinizio = 1
        !            ifine   = n1
        !
        !            jinizio = 1
        !            jfine   = n2
        !
        !            kinizio = 1
        !            kfine   = n3
        !        else
        !            istop = 0
        !            if(iinizio .lt. 1) istop = 1
        !            if(iinizio .lt. 1) istop = 1
        !            if(iinizio .lt. 1) istop = 1
        !            if(ifine   .gt. n1)istop = 1
        !            if(jfine   .gt. n2)istop = 1
        !            if(kfine   .gt. n3)istop = 1
        !
        !            if( (ifine-iinizio).lt.0 ) istop = 1
        !            if( (jfine-jinizio).lt.0 ) istop = 1
        !            if( (kfine-kinizio).lt.0 ) istop = 1
        !
        !            do i=1,istop
        !                if(myid==0)then
        !                    write(*,*)'error in solid phase limit, see Aricerca.inp'
        !                    write(*,*)'search limited to:',iinizio, ifine,'while in x the domain is',1,n1
        !                    write(*,*)'search limited to:',jinizio, jfine,'while in y the domain is',1,n2
        !                    write(*,*)'search limited to:',kinizio, kfine,'while in z the domain is',1,n3
        !                end if
        !                stop
        !            end do
        !        end if
        !
        !        do i=1,miro_solide
        !            if(myid==0)then
        !                write(*,*)'PAY ATTENTION'
        !                write(*,*)'search limited to:',iinizio, ifine,'while in x the domain is',1,n1
        !                write(*,*)'search limited to:',jinizio, jfine,'while in y the domain is',1,n2
        !                write(*,*)'search limited to:',kinizio, kfine,'while in z the domain is',1,n3
        !            end if
        !        end do


        !        MN = 0
        !        MP = 0
        !
        !        open (11,file='Coordinate_IMMB.inp')
        !        read (11,*) MN
        !        close(11)
        !
        !        open (12,file='Indici_piani_IMMB.inp')
        !        read (12,*) MP
        !        close(12)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid==0)then
            write(*,*)' ---------------------------------------'
            write(*,*)'| SECTION 1: phases search --> started'
            write(*,*)' ---------------------------------------'
            write(*,*)' '
        end if

        call cpu_time(time1)
        call phase_search(tipo)
        call cpu_time(time2)
        time = time2 - time1

        if(myid==0)then
            write(*,*)' ---------------------------------------'
            write(*,*)'| SECTION 1: phases search --> finished'
            write(*,*)'| Time expired: ',time
            write(*,*)' ---------------------------------------'
            write(*,*)' '
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid==0)then
            write(*,*)' ---------------------------------------'
            write(*,*)'| SECTION 2: IP search --> started'
            write(*,*)' ---------------------------------------'
            write(*,*)' '
        end if

        call cpu_time(time1)
        call interface_search(tipo)
        call cpu_time(time2)
        time = time2 - time1

        if(myid==0)then
            write(*,*)' ---------------------------------------'
            write(*,*)'| SECTION 2: IP search --> finished'
            write(*,*)'| Time expired: ',time
            write(*,*)' ---------------------------------------'
            write(*,*)' '
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid==0)then
            write(*,*)' ---------------------------------------'
            write(*,*)'| SECTION 3: PP search --> started'
            write(*,*)' ---------------------------------------'
            write(*,*)' '
        end if

        call cpu_time(time1)
        call projection_search(tipo)
        call cpu_time(time2)
        time = time2 - time1

        if(myid==0)then
            write(*,*)' ---------------------------------------'
            write(*,*)'| SECTION 3: PP search --> finished'
            write(*,*)'| Time expired: ',time
            write(*,*)' ---------------------------------------'
            write(*,*)' '
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


        call prepare_communication()

        return

    end subroutine do_ricerca


    subroutine phase_search(tipo)

        implicit none

        integer,intent(out) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        ! ----------------------------------------------------
        ! For solid/fluid/IBM research

        real(8) :: isoparameter
        integer :: solid_count
        integer :: counterTot,reducedCounterTot
        integer :: counterFluid,contatore,num_solide_real
        integer :: reducedCounterFluid

        integer,parameter :: maxSolidIndex=6

        integer,allocatable :: neighbor(:,:,:,:,:)
        logical,allocatable :: isValid(:,:,:,:)
        integer :: i,j,k,n
        integer :: ni,nj,nk
        integer :: indexHere
        integer :: sx,sy,sz
        integer :: ierr
        integer :: ksta,kend

        if(myid==0)then
            !write(*,*)'Number of points: ',MN
            !write(*,*)'Number of triangle: ',MP
            !write(*,*)'the domain dimensions are:',n1,n2,n3
        end if

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

        allocate(neighbor(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,6,3))
        allocate(isValid(0:n1+1,0:n2+1,kparasta-1:kparaend+1,6))

        !allocate(Indici_Boundary(MP,3))
        allocate(solidIndexSize(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(solidIndex(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,maxSolidindex))

        !---------------------------------------------------------------------------------------------------------------------------------------------
        !*********************************************************************************************************************************************
        !---------------------------------------------------------------------------------------------------------------------------------------------
        ! SECTION 2: PHASES

        ! determine if a cell is solid, fluid or ib starting from zero or from a restart file
        ! (0,1,2 MEANS fluid, solid, ib)


        !        ista = iinizio
        !        iend = ifine
        !        jsta = jinizio
        !        jend = jfine
        !        ksta = kinizio
        !        kend = kfine

        write(*,*) 'Proc ',myid,'scanning ',kparasta-deepl,' to ',kparaend+deepr


        ! total number of cells
        counterTot=n1*n2*(kparaend-kparasta+1)
        call MPI_ALLREDUCE(counterTot,reducedCounterTot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! set everything fluid
        tipo=2 !0

        ! reset solid indices
        solidIndex=0
        solidIndexSize=0

        num_solide=0
        contatore=0
        if (myid==0) then
            write(*,*) 'Looking for solid nodes'
        end if
        do k=kparasta-deepl,kparaend+deepr
            do j=0,n2+1
                do i=0,n1+1
                    do n=1,totParticles
                        if (tipo(i,j,k)/=0) then !1
                            isoparameter=(xcd(i,j,k)-x_sphere(n))**2 + &
                                (ycd(i,j,k)-y_sphere(n))**2 + &
                                (zcd(i,j,k)-z_sphere(n))**2

                                if (isoparameter<radius2(n)) then

                                    tipo(i,j,k)=0 !1
                                    solidIndex(i,j,k,1)=n

                                    if (i>=1 .and. i<=n1 .and. j>=1 .and. j<=n2) then
                                        num_solide=num_solide+1
                                        if (k > kparaend .or. k < kparasta) then
                                            contatore=contatore+1
                                        end if
                                    end if
                                end if
                            end if
                        end do
                end do
            end do
        end do

        num_solide_real=num_solide-contatore !without border cells

        if (myid==0) then
            write(*,*) 'Building up neighbor matrix'
        end if
        ! compute neighbors
        isValid=.false.
        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    do n=1,6
                        neighbor(i,j,k,n,1)=i
                        neighbor(i,j,k,n,2)=j
                        neighbor(i,j,k,n,3)=k
                    end do
                    neighbor(i,j,k,1,1)=i+1
                    neighbor(i,j,k,2,1)=i-1
                    neighbor(i,j,k,3,2)=j+1
                    neighbor(i,j,k,4,2)=j-1
                    neighbor(i,j,k,5,3)=k+1
                    neighbor(i,j,k,6,3)=k-1
                    if (i+1<=n1+1) isValid(i,j,k,1)=.true.
                    if (i-1>=0) isValid(i,j,k,2)=.true.
                    if (j+1<=n2+1) isValid(i,j,k,3)=.true.
                    if (j-1>=0) isValid(i,j,k,4)=.true.
                    if (k+1<=n3+1) isValid(i,j,k,5)=.true.
                    if (k-1>=0) isValid(i,j,k,6)=.true.
                end do
            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! set boundaries
        if (myid==0) then
            write(*,*) 'Looking for IB nodes'
        end if
        num_ib=0
        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    ! check if it's fluid
                    if (tipo(i,j,k)==2) then !0
                        ! check if neighbors are solid. in this case, turn the node into IB
                        do n = 1,6
                            if (isValid(i,j,k,n)) then
                                ni=neighbor(i,j,k,n,1)
                                nj=neighbor(i,j,k,n,2)
                                nk=neighbor(i,j,k,n,3)
                                if (tipo(ni,nj,nk)==0) then !
                                    ! update solid neighbors
                                    indexHere=solidIndex(ni,nj,nk,1)
                                    solidIndexSize(i,j,k)=solidIndexSize(i,j,k)+1
                                    solidIndex(i,j,k,solidIndexSize(i,j,k))=indexHere
                                    ! turn it into IB only if that has not happened before
                                    ! (e.g. with another neighbor)
                                    if (tipo(i,j,k)/=1) then !2
                                        tipo(i,j,k)=1 !2
                                        if (k>=kparasta .and. k<=kparaend .and. &
                                            i>=1 .and. i<=n1 .and. j>=1 .and. j<=n2) then
                                            num_ib=num_ib+1
                                        end if
                                    end if

                                end if
                            end if
                        end do
                    end if
                end do
            end do
        end do

        MN=num_ib

        ! check number of fluid cells
        counterFluid=0
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    if (tipo(i,j,k)==2) then !0
                        counterFluid=counterFluid+1
                    end if
                end do
            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        call MPI_ALLREDUCE(num_solide,numero_celle_bloccate,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(num_solide_real,numero_celle_bloccate_real,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(counterFluid,reducedCounterFluid,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(num_ib,numero_celle_IB,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! checking for mistakes
        if (myid==0) then
            write(*,*) 'Total cells = ',reducedCounterTot
            write(*,*) '   Total solid = ',numero_celle_bloccate_real
            write(*,*) '   Total IB = ',numero_celle_IB
            write(*,*) '   Total fluid = ',reducedCounterFluid


            if (reducedCounterTot/=numero_celle_bloccate_real+numero_celle_IB+reducedCounterFluid) then
                write(*,*) 'Error in assigning type'
                call exit(0)
            end if
        end if

        ! now we can allocate

        allocate(indici_CELLE_IB(num_ib,6))
        allocate(indici_celle_bloccate(num_solide,3))
        allocate(distanze_CELLE_IB(num_ib,3))
        allocate(dist_pp_ib(num_ib))
        allocate(dist_ib_parete(num_ib))
        allocate(proiezioni(num_ib,3))
        allocate(ustar(num_ib),pressure_ib(num_ib),shear_ib(num_ib,3))

        allocate(tricoef(num_ib,4))
        allocate(trind(num_ib,4,3))

        allocate(rot(num_ib,3,3))
        allocate(rot_inverse(num_ib,3,3))
        allocate(r_solid(nscal,num_solide))

        solid_count=0
        do k=kparasta-deepl,kparaend+deepr
            do j=1,n2
                do i=1,n1
                    if (tipo(i,j,k)==0) then !1

                        solid_count=solid_count+1

                        indici_celle_bloccate(solid_count,1)=i
                        indici_celle_bloccate(solid_count,2)=j
                        indici_celle_bloccate(solid_count,3)=k

                    end if
                end do
            end do
        end do

        if (solid_count/=num_solide) then
            write(*,*) 'Proc: ',myid,' solid_count=',solid_count,', num_solide=',num_solide
            call exit(0)
        end if

        MP=num_solide

        return

    end subroutine phase_search

    subroutine interface_search(tipo)

        implicit none

        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        ! normal to surface vector components
        real(8) :: normalVectorXHere, normalVectorYHere, normalVectorZHere
        real(8),allocatable :: normalVectorX(:,:,:), normalVectorY(:,:,:), normalVectorZ(:,:,:)
        ! IP_BP distance
        real(8),allocatable :: surfaceDistance(:,:,:)
        real(8) :: distanceNormHere
        integer :: fileUnit

        integer :: i,j,k,n
        integer :: sx,sy,sz
        integer :: indexHere, solidIndexSizeHere
        integer :: ierr

        !---------------------------------------------------------------------------------------------------------------------------------------------
        ! SECTION 2: search the projection point IP on the body, through the normal from the IB to the body


        allocate(nodo_vicino_x_array(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(nodo_vicino_y_array(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(nodo_vicino_z_array(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(normalVectorX(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(normalVectorY(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(normalVectorZ(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(surfaceDistance(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))

        nodo_vicino_x_array=0.0
        nodo_vicino_y_array=0.0
        nodo_vicino_z_array=0.0
        normalVectorX=0.0
        normalVectorY=0.0
        normalVectorZ=0.0
        surfaceDistance=0.0

        do k=kparasta-deepl,kparaend+deepr
            do j=0,n2+1
                do i=0,n1+1
                    if(tipo(i,j,k)==1)then !2
                        solidIndexSizeHere=solidIndexSize(i,j,k)
                        ! we start from the position of the node itself, then we shift it later
                        nodo_vicino_x_array(i,j,k)=xcd(i,j,k)
                        nodo_vicino_y_array(i,j,k)=ycd(i,j,k)
                        nodo_vicino_z_array(i,j,k)=zcd(i,j,k)
                        ! update normal and surface distance (one shift for every solid enighbor)
                        if (solidIndexSizeHere<=0) then
                            write(*,*) 'Problem: negative index size'
                            call exit(0)
                        end if
                        do n=1,solidIndexSizeHere
                            ! compute vector normal to surface
                            indexHere=solidIndex(i,j,k,n)
                            if (indexHere<=0 .or. indexHere>totParticles) then
                                write(*,*) 'Problem: solid index: ',indexHere
                                call exit(0)
                            end if
                            normalVectorXHere=xcd(i,j,k)-x_sphere(indexHere)
                            normalVectorYHere=ycd(i,j,k)-y_sphere(indexHere)
                            normalVectorZHere=zcd(i,j,k)-z_sphere(indexHere)
                            ! normalize it
                            distanceNormHere=sqrt(normalVectorXHere**2.0+normalVectorYHere**2.0+normalVectorZHere**2.0)
                            if (distanceNormHere==0.0) then
                                write(*,*) 'XP=',xcd(i,j,k),'; XS=',x_sphere(indexHere),'; DX=',normalVectorXHere
                            end if
                            normalVectorXHere=normalVectorXHere/distanceNormHere
                            normalVectorYHere=normalVectorYHere/distanceNormHere
                            normalVectorZHere=normalVectorZHere/distanceNormHere
                            ! shift normal
                            normalVectorX(i,j,k)=normalVectorX(i,j,k)+normalVectorXHere/solidIndexSizeHere
                            normalVectorY(i,j,k)=normalVectorY(i,j,k)+normalVectorYHere/solidIndexSizeHere
                            normalVectorZ(i,j,k)=normalVectorZ(i,j,k)+normalVectorZHere/solidIndexSizeHere
                            ! update distance from surface
                            surfaceDistance(i,j,k)=surfaceDistance(i,j,k)+(distanceNormHere-radius(indexHere))/solidIndexSizeHere
                            if (surfaceDistance(i,j,k)<=0) then
                                write(*,*) 'Problem: negative surface distance: ',surfaceDistance(i,j,k),&
                                    '; solidIndexSizeHere = ',solidIndexSizeHere,'; n = ',n
                                write(*,*) 'sphere = ',indexHere,'; radius = ',radius(indexHere),'; norm = ',distanceNormHere
                                call exit(0)
                            end if
                        end do
                        ! shift interface point IP
                        nodo_vicino_x_array(i,j,k)=nodo_vicino_x_array(i,j,k)-normalVectorX(i,j,k)*surfaceDistance(i,j,k)
                        nodo_vicino_y_array(i,j,k)=nodo_vicino_y_array(i,j,k)-normalVectorY(i,j,k)*surfaceDistance(i,j,k)
                        nodo_vicino_z_array(i,j,k)=nodo_vicino_z_array(i,j,k)-normalVectorZ(i,j,k)*surfaceDistance(i,j,k)

                    end if
                end do
            end do
        end do

    end subroutine interface_search

    subroutine projection_search(tipo)

        implicit none

        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        !

        real(8), allocatable :: vicini_x(:),vicini_y(:),vicini_z(:)

        real(8), allocatable :: proiezioni_x(:),proiezioni_y(:),proiezioni_z(:)

        real(8), allocatable :: parete_x(:),parete_y(:),parete_z(:)

        real(8), allocatable :: distanza_1(:),distanza_2(:)

        integer, allocatable :: ib_ind_i(:),ib_ind_j(:),ib_ind_k(:)
        integer, allocatable :: v_ind_i(:),v_ind_j(:),v_ind_k(:)

        real(8), allocatable :: dist_dir_x(:),dist_dir_y(:),dist_dir_z(:)

        real(8), allocatable :: rotore1(:),rotore2(:),rotore3(:), &
            rotore4(:),rotore5(:),rotore6(:), &
            rotore7(:),rotore8(:),rotore9(:)

        real(8), allocatable :: rotore1_inv(:),rotore2_inv(:),rotore3_inv(:), &
            rotore4_inv(:),rotore5_inv(:),rotore6_inv(:), &
            rotore7_inv(:),rotore8_inv(:),rotore9_inv(:)

        real(8), allocatable :: var_real1(:),var_real2(:),var_real3(:), &
            var_real4(:),var_real5(:),var_real6(:), &
            var_real7(:),var_real8(:),var_real9(:)

        real(8), allocatable :: var_real1_inv(:),var_real2_inv(:),var_real3_inv(:), &
            var_real4_inv(:),var_real5_inv(:),var_real6_inv(:), &
            var_real7_inv(:),var_real8_inv(:),var_real9_inv(:)

        integer, allocatable :: var_int1(:),var_int2(:),var_int3(:), &
            var_int4(:),var_int5(:),var_int6(:)


        real(8) trilinear_coef(4)
        real(8), allocatable :: trilinear_coef_col1(:),trilinear_coef_col2(:), &
            trilinear_coef_col3(:),trilinear_coef_col4(:)

        integer, allocatable :: trilinear_index_col1(:),trilinear_index_col2(:), &
            trilinear_index_col3(:),trilinear_index_col4(:)

        integer, allocatable :: trilinear_index(:,:,:)

        ! SECTION 3: calculate the mirror point PP, on the other side respect to IB of the IP
        !                                          | solid phase
        !                PP -------- IB -------- IP| solid phase
        !                                          | solid phase
        !          then determine the closest fluid node V to PP

        !     STEP 1: if PP is outside the grid I interpolate it on the boundary

        !     STEP 2: compute the closest fluid node V to PP, ic,jc,kc will be its indices

        !     STEP 3: compute the "real" colsest node to PP considering also the IB

        !     STEP 4: compute how man2 IB nodes (contatore_stencil) stay in the V node stencil

        !     STEP 5: if there are too man2 IB nodes in the V stencil I search for a new V
        !             NOTE: the number of allowed IB point in the stencil is defined as input parameter

        !     STEP 6: put the PP inside the V cell, then PP-IB =/ IB_IP
        !             NOTE: this step is activated by input parameter "shift_proiezione"

        !------------------------------------------------------------------------

        integer, parameter ::  iread_linea = 0
        real(8), parameter ::  calibro_iter = 3.5
        real(8), parameter ::  calibro_iter_inc = 0.2
        real(8), parameter ::  calibro_punto = 100.
        real(8), parameter ::  calibro_punto_inc = 5.
        ! parameter


        !MN,MP matrix dimension to store mesh: points and triangle
        !integer :: MN,MP

        !input=0 start from zero to search the solid phase
        !     =1 start from restart file: "tipologia.dat" , solid phase is known
        !     =2 start from restart file for n bodies "tipologia//n bodies//.dat", solid phase is known
        integer,parameter :: input=0
        integer,parameter :: corpi=1
        integer,parameter :: restart_solid=0
        integer,parameter :: restart_solid_num=0
        integer,parameter :: irestart_ip=0
        integer,parameter :: kprint=90
        integer,parameter :: miro_solide=0 ! consider switching this back on, it's actually cool

        !do not change this parameter, interpolate PP node on boundary if it go outside the grid,
        integer, parameter :: correzione_parete=1
        !do not change this parameter, search V node in a faster way
        integer, parameter :: direzione=1
        !number of IB points allowed to be present in the stencil around the V node
        !integer, parameter :: num_ib_stencil=3
        !1 la proiezione viene avvicinata al nodo vicino, 0 distanze uguali
        integer, parameter :: shift_proiezione=1
        !1 wall, 0 periodic or inflow put to solid an ib cell close to the boundary and solid
        integer, parameter :: ip=1,jp=1,kp=1


        real(8), parameter :: tol=1.0D-14
        real(8), parameter :: tolC=1.0D-14
        !---------------------------------------------------------------------------------
        ! matrix
        !real(8),allocatable :: X(:,:,:,:)       ! grid
        !real(8),allocatable :: xcd(:,:,:),ycd(:,:,:),zcd(:,:,:)  ! grid cell center
        !real(8),allocatable :: Coord_Boundary(:,:)    ! mesh points coordinates x y z

        ! mesh triangle indicies, every index refers to a mesh point
        !integer,allocatable :: Indici_Boundary(:,:)

        !integer,allocatable :: vertice_comune(:)
        !integer,allocatable :: vertice_comune_prov(:)

        !index
        integer :: i,j,k
        ! V indices
        integer :: ic,jc,kc

        ! search area for V
        integer :: ista,iend
        integer :: jsta,jend
        integer :: ksta,kend

        ! active procedure
        logical :: trovato_vicino

        ! incremental
        integer :: contatore_mesh_rozza
        ! cell center coordinates
        real(8) :: sx,sy,sz
        ! V coordinates
        ! real(8) :: mx,my,mz
        !PP coordinates
        real(8) :: nodo_proiezione_x, nodo_proiezione_y, nodo_proiezione_z
        !IP coordinates
        real(8) :: nodo_vicino_x, nodo_vicino_y, nodo_vicino_z

        !grid nodes
        !integer :: n1,n2,n3
        !integer :: Numero_Celle_IB,Numero_Celle_Bloccate

        ! distance between PP and V
        real(8) :: distanza_x, distanza_y,distanza_z
        real(8) :: distanza_np_ib,distanza_ib_parete
        !grid length
        !real(8) :: h1,h2,h3

        !rotation
        real(8) :: rot_ric(3,3),rot_inverse_ric(3,3)

        !integer :: lll
        !integer,allocatable :: linea(:,:)
        integer :: istop

        !character*15 :: FMT1
        !character*8 :: grid_format

        integer,parameter :: formatted=0
        !character*150 commento

        integer ib_count

        !integer :: ib_count

        real(8) r1,r2,r3,r4,r5,r6,r7,r8,r9
        real(8) r1_in,r2_in,r3_in,r4_in,r5_in,r6_in,r7_in,r8_in,r9_in
        integer count_rot
        integer ntri,check_box_solid
        real(8) count_trilinear_coef
        integer icount_tri
        character*60 filename1,filename2

        integer pp_found
        integer iref,jref,kref
        real(8) x0,y0,z0,FF,GG,HH
        logical trovato_pp
        real(8) val_trilinear

        integer,parameter :: icheck_print=0

        integer :: ierr

        !-----------------------------------------------------------------------
        ! PROJECTION POINT PP SEARCHIG PROCEDURE
        !-----------------------------------------------------------------------

        contatore_mesh_rozza = 0

        allocate(trilinear_coef_col1(n1*n2*n3/nproc))
        allocate(trilinear_coef_col2(n1*n2*n3/nproc))
        allocate(trilinear_coef_col3(n1*n2*n3/nproc))
        allocate(trilinear_coef_col4(n1*n2*n3/nproc))

        allocate(trilinear_index(n1*n2*n3/nproc,4,3))

        trilinear_coef_col1 = -999999.
        trilinear_coef_col2 = -999999.
        trilinear_coef_col3 = -999999.
        trilinear_coef_col4 = -999999.


        trilinear_index = -999999

        trilinear_coef = -999999.

        !-------------------------------------------------------------------------------
        icount_tri = 0
        ib_count = 1
        do k = kparasta,kparaend !1,n3
            if(myid==0)write(*,*)'internal loop PP:',k,'of',n3/nproc
            do j = 1,n2
                do i = 1,n1
                    if(tipo(i,j,k) == 1)then !2

                        trovato_vicino = .false.

                        sx = xcd(i,j,k)
                        sy = ycd(i,j,k)
                        sz = zcd(i,j,k)

                        nodo_vicino_x = nodo_vicino_x_array(i,j,k)
                        nodo_vicino_y = nodo_vicino_y_array(i,j,k)
                        nodo_vicino_z = nodo_vicino_z_array(i,j,k)


                        !compute the distance on the normal
                        distanza_ib_parete = sqrt( (nodo_vicino_x-sx)**2. &
                            + (nodo_vicino_y-sy)**2. &
                            + (nodo_vicino_z-sz)**2. )



                        !----------------------------------------------------------------------------------
                        !search closest fluid node to PP
                        do ntri=1,4 ! initialize current box to IB index
                            trilinear_index(ib_count,ntri,1)=i
                            trilinear_index(ib_count,ntri,2)=j
                            trilinear_index(ib_count,ntri,3)=k

                            trilinear_coef(ntri) = 0.
                        end do


                        !----------------------------------------------------------------------------------
                                  !explicit line IB-IP
                        ic = i
                        jc = j
                        kc = k


                        call line_exp2par_3d(sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            FF,GG,HH,x0,y0,z0)
                        !                    if(icheck_print == 1)then
                        !                        write(3500+myid,*)'ib',sx,sy,sz
                        !                        write(3500+myid,*)'ip',nodo_vicino_x,nodo_vicino_y,nodo_vicino_z
                        !                        write(3500+myid,*)'retta',FF,GG,HH,x0,y0,z0
                        !                        write(3500+myid,*)'indice ref',i,j,k
                        !                    end if

                        ! AAA da chiamare una volta sola fuori dal loop
                        call triangles_mesh


                        !AAA se trovato lo metto qui false e non dentro poi quando e' trovato non passa nelle altre
                        ! sub ma stare attenti alla parte finale della sub
                        pp_found = 0
                        iref = i-1
                        call pp_on_meshbox_diri(iref,j,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if(trovato_pp)pp_found = pp_found + 1
                        if(icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if(myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        iref = i+1
                        call pp_on_meshbox_diri(iref,j,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if(trovato_pp)pp_found = pp_found + 1
                        if(icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if(myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        jref = j-1
                        call pp_on_meshbox_dirj(jref,i,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if(trovato_pp)pp_found = pp_found + 1
                        if(icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if(myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        jref = j+1
                        call pp_on_meshbox_dirj(jref,i,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if(trovato_pp)pp_found = pp_found + 1
                        if(icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if(myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        kref = k-1
                        call pp_on_meshbox_dirk(kref,i,j,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if(trovato_pp)pp_found = pp_found + 1
                        if(icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if(myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        kref = k+1
                        call pp_on_meshbox_dirk(kref,i,j,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if(trovato_pp)pp_found = pp_found + 1
                        if(icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                             !if(myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        if(pp_found == 1)then
                            call find_trilinear_coef2(ib_count,trilinear_index,trilinear_coef, &
                                nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                                n1,n2,n3,icheck_print,i,j,k)

                            distanza_np_ib = sqrt( (nodo_proiezione_x-sx)**2. &
                                + (nodo_proiezione_y-sy)**2. &
                                + (nodo_proiezione_z-sz)**2. )

                            trilinear_coef_col1(ib_count) = trilinear_coef(1)
                            trilinear_coef_col2(ib_count) = trilinear_coef(2)
                            trilinear_coef_col3(ib_count) = trilinear_coef(3)
                            trilinear_coef_col4(ib_count) = trilinear_coef(4)

                            ! v node
                            val_trilinear = 0.
                            do ntri=1,4
                                if(trilinear_coef(ntri) > val_trilinear)then
                                    val_trilinear = trilinear_coef(ntri)

                                    ic = trilinear_index(ib_count,ntri,1)
                                    jc = trilinear_index(ib_count,ntri,2)
                                    kc = trilinear_index(ib_count,ntri,3)

                                    if(ic == n1+1)ic = n1
                                    if(ic == 0   )ic = 1

                                    if(jc == n2+1)jc = n2
                                    if(jc == 0   )jc = 1

                                    if(kc == n3+1)kc = n3
                                    if(kc == 0   )kc = 1
                                    !write(*,*)ntri,ib_count,'ntri',i,j,k,ic,jc,kc

                                    distanza_x=nodo_proiezione_x-xcd(ic,jc,kc)
                                    distanza_y=nodo_proiezione_y-ycd(ic,jc,kc)
                                    distanza_z=nodo_proiezione_z-zcd(ic,jc,kc)
                                end if
                            end do

                            !trilinear coef for solid = 0

                            do ntri = 1,4
                                if(trilinear_index(ib_count,ntri,2).le. n2)then
                                    if(trilinear_index(ib_count,ntri,2).ge. 1 )then
                                        if(trilinear_index(ib_count,ntri,1).le. n1)then
                                            if(trilinear_index(ib_count,ntri,1).ge. 1 )then
                                                if(trilinear_index(ib_count,ntri,3).le. n3)then
                                                    if(trilinear_index(ib_count,ntri,3).ge. 1 )then
                                                        if( tipo(trilinear_index(ib_count,ntri,1), &
                                                            trilinear_index(ib_count,ntri,2), &
                                                            trilinear_index(ib_count,ntri,3)) == 0)then !1

                                                            trilinear_coef(ntri) = 0.

                                                        end if
                                                    end if
                                                end if
                                            end if
                                        end if
                                    end if
                                end if
                            end do


                        else
                            write(*,*)'ID=',myid,', problem on ib',i,j,k,'in finding PP',pp_found, &
                            'vicino=(',sx,sy,sz,')'
                        end if




                        count_trilinear_coef = 0.
                        do ntri=1,4
                            count_trilinear_coef = count_trilinear_coef + trilinear_coef(ntri)
                        end do
                        if(count_trilinear_coef == 0.)then
                            icount_tri = icount_tri + 1
                           !write(*,*)'coef nulli',i,j,k
                        end if

                        if(icheck_print == 1)then
                            write(6700+myid,*)'ib',i,j,k,'found',pp_found,ib_count
                            write(6700+myid,*)trilinear_index(ib_count,1,1), &
                                trilinear_index(ib_count,1,2), &
                                trilinear_index(ib_count,1,3)
                            write(6700+myid,*)trilinear_index(ib_count,2,1), &
                                trilinear_index(ib_count,2,2), &
                                trilinear_index(ib_count,2,3)
                            write(6700+myid,*)trilinear_index(ib_count,3,1), &
                                trilinear_index(ib_count,3,2), &
                                trilinear_index(ib_count,3,3)
                            write(6700+myid,*)trilinear_index(ib_count,4,1), &
                                trilinear_index(ib_count,4,2), &
                                trilinear_index(ib_count,4,3)
                            write(6700+myid,*)(trilinear_coef(ntri),ntri=1,4)



                            write(3000,*)nodo_vicino_x,nodo_vicino_y,nodo_vicino_z
                            write(3000,*)sx,sy,sz
                            write(3000,*)nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z
                        end if

                        if(icheck_print == 1)then
                            write(3100,*)xcd(trilinear_index(ib_count,1,1), &
                                trilinear_index(ib_count,1,2), &
                                trilinear_index(ib_count,1,3)), &
                                ycd(trilinear_index(ib_count,1,1), &
                                trilinear_index(ib_count,1,2), &
                                trilinear_index(ib_count,1,3)), &
                                zcd(trilinear_index(ib_count,1,1), &
                                trilinear_index(ib_count,1,2), &
                                trilinear_index(ib_count,1,3))
                            write(3100,*)xcd(trilinear_index(ib_count,2,1), &
                                trilinear_index(ib_count,2,2), &
                                trilinear_index(ib_count,2,3)), &
                                ycd(trilinear_index(ib_count,2,1), &
                                trilinear_index(ib_count,2,2), &
                                trilinear_index(ib_count,2,3)), &
                                zcd(trilinear_index(ib_count,2,1), &
                                trilinear_index(ib_count,2,2), &
                                trilinear_index(ib_count,2,3))
                            write(3100,*)xcd(trilinear_index(ib_count,3,1), &
                                trilinear_index(ib_count,3,2), &
                                trilinear_index(ib_count,3,3)), &
                                ycd(trilinear_index(ib_count,3,1), &
                                trilinear_index(ib_count,3,2), &
                                trilinear_index(ib_count,3,3)), &
                                zcd(trilinear_index(ib_count,3,1), &
                                trilinear_index(ib_count,3,2), &
                                trilinear_index(ib_count,3,3))
                            write(3100,*)xcd(trilinear_index(ib_count,4,1), &
                                trilinear_index(ib_count,4,2), &
                                trilinear_index(ib_count,4,3)), &
                                ycd(trilinear_index(ib_count,4,1), &
                                trilinear_index(ib_count,4,2), &
                                trilinear_index(ib_count,4,3)), &
                                zcd(trilinear_index(ib_count,4,1), &
                                trilinear_index(ib_count,4,2), &
                                trilinear_index(ib_count,4,3))
                        end if

                        !----------------------------------------------------------------------------------
                        !----------------------------------------------------------------------------------
                        !----------------------------------------------------------------------------------
                                  ! I check how man2 IB nodes (contatore_stencil) stay in the V node stencil
                        !     contatore_stencil=0
                        !     do ii=ic-1,ic+1,2
                        !        if(ii.ne.0 .and. ii.ne.(n1+1))then
                        !          if(tipo(ii,jc,kc)==2) contatore_stencil = contatore_stencil+1
                        !        end if
                        !!    end do
                        !     do jj=jc-1,jc+1,2
                        !        if(jj.ne.0 .and. jj.ne.(n2+1))then
                        !          if(tipo(ic,jj,kc)==2) contatore_stencil = contatore_stencil+1
                        !        end if
                        !     end do
                        !     do kk=kc-1,kc+1,2
                        !!       if(kk.ne.0 .and. kk.ne.(n3+1))then
                        !          if(tipo(ic,jc,kk)==2) contatore_stencil = contatore_stencil+1
                        !!       end if
                        !!    end do


                        !----------------------------------------------------------------------------------
                        ! store data in vector

                        !v node

                        indici_CELLE_IB(ib_count,1)=i
                        indici_CELLE_IB(ib_count,2)=j
                        indici_CELLE_IB(ib_count,3)=k
                        indici_CELLE_IB(ib_count,4)=ic
                        indici_CELLE_IB(ib_count,5)=jc
                        indici_CELLE_IB(ib_count,6)=kc

                        distanze_CELLE_IB(ib_count,1)=distanza_x
                        distanze_CELLE_IB(ib_count,2)=distanza_y
                        distanze_CELLE_IB(ib_count,3)=distanza_z

                        dist_pp_ib(ib_count)=distanza_np_ib
                        dist_ib_parete(ib_count)=distanza_ib_parete

                        proiezioni(ib_count,1)=nodo_vicino_x
                        proiezioni(ib_count,2)=nodo_vicino_y
                        proiezioni(ib_count,3)=nodo_vicino_z

                        !if (myid==0) then
!                        write(*,*) 'IB:',ib_count, &
!                        ' (',indici_CELLE_IB(ib_count,1), &
!                        ',',indici_CELLE_IB(ib_count,2), &
!                        ',',indici_CELLE_IB(ib_count,3), &
!                        '), V:(',indici_CELLE_IB(ib_count,4), &
!                        ',',indici_CELLE_IB(ib_count,4), &
!                        ',',indici_CELLE_IB(ib_count,5), &
!                        ')'
                        !end if
                        !-----------------------------------

                        !vicini_x(ib_count) = xcd(ic,jc,kc)
                        !vicini_y(ib_count) = ycd(ic,jc,kc)
                        !vicini_z(ib_count) = zcd(ic,jc,kc)

                        !proiezioni_x(ib_count)=nodo_proiezione_x
                        !proiezioni_y(ib_count)=nodo_proiezione_y
                        !proiezioni_z(ib_count)=nodo_proiezione_z

                        !parete_x(ib_count) = nodo_vicino_x
                        !parete_y(ib_count) = nodo_vicino_y
                        !parete_z(ib_count) = nodo_vicino_z

                        !distanza_1(ib_count) = distanza_np_ib
                        !distanza_2(ib_count) = distanza_ib_parete

                        !ib_ind_i(ib_count) = i
                        !ib_ind_j(ib_count) = j
                        !ib_ind_k(ib_count) = k

                        !v_ind_i(ib_count) = ic
                        !v_ind_j(ib_count) = jc
                        !v_ind_k(ib_count) = kc

                        !dist_dir_x(ib_count)=distanza_x
                        !dist_dir_y(ib_count)=distanza_y
                        !dist_dir_z(ib_count)=distanza_z

                        !----------------------------------------------------------------------------------
                        ! ROTAZION MATRIX
                        call rotation2(sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z,rot_ric,rot_inverse_ric) !,alfa,beta,gamma,vertice)

                        ! rotation matrix (to construct tangential and normal velocity)
                        rot(ib_count,1,1)=rot_ric(1,1)
                        rot(ib_count,1,2)=rot_ric(1,2)
                        rot(ib_count,1,3)=rot_ric(1,3)

                        rot(ib_count,2,1)=rot_ric(2,1)
                        rot(ib_count,2,2)=rot_ric(2,2)
                        rot(ib_count,2,3)=rot_ric(2,3)

                        rot(ib_count,3,1)=rot_ric(3,1)
                        rot(ib_count,3,2)=rot_ric(3,2)
                        rot(ib_count,3,3)=rot_ric(3,3)


                        !         inverse rotation matrix
                        rot_inverse(ib_count,1,1)=rot_inverse_ric(1,1)
                        rot_inverse(ib_count,1,2)=rot_inverse_ric(1,2)
                        rot_inverse(ib_count,1,3)=rot_inverse_ric(1,3)

                        rot_inverse(ib_count,2,1)=rot_inverse_ric(2,1)
                        rot_inverse(ib_count,2,2)=rot_inverse_ric(2,2)
                        rot_inverse(ib_count,2,3)=rot_inverse_ric(2,3)

                        rot_inverse(ib_count,3,1)=rot_inverse_ric(3,1)
                        rot_inverse(ib_count,3,2)=rot_inverse_ric(3,2)
                        rot_inverse(ib_count,3,3)=rot_inverse_ric(3,3)


                        !                        rotore1(ib_count)   = rot_ric(1,1)
                        !                        rotore2(ib_count)   = rot_ric(1,2)
                        !                        rotore3(ib_count)   = rot_ric(1,3)
                        !                        rotore4(ib_count)   = rot_ric(2,1)
                        !                        rotore5(ib_count)   = rot_ric(2,2)
                        !                        rotore6(ib_count)   = rot_ric(2,3)
                        !                        rotore7(ib_count)   = rot_ric(3,1)
                        !                        rotore8(ib_count)   = rot_ric(3,2)
                        !                        rotore9(ib_count)   = rot_ric(3,3)
                        !
                        !                        rotore1_inv(ib_count)     = rot_inverse_ric(1,1)
                        !                        rotore2_inv(ib_count)     = rot_inverse_ric(1,2)
                        !                        rotore3_inv(ib_count)     = rot_inverse_ric(1,3)
                        !                        rotore4_inv(ib_count)     = rot_inverse_ric(2,1)
                        !                        rotore5_inv(ib_count)     = rot_inverse_ric(2,2)
                        !                        rotore6_inv(ib_count)     = rot_inverse_ric(2,3)
                        !                        rotore7_inv(ib_count)     = rot_inverse_ric(3,1)
                        !                        rotore8_inv(ib_count)     = rot_inverse_ric(3,2)
                        !                        rotore9_inv(ib_count)     = rot_inverse_ric(3,3)

                        !----------------------------------------------------------------------------------

                        !first control, solid inside the box?
                        check_box_solid = 0
                        do ntri = 1,4
                                                    !write(*,*)    trilinear_index(ib_count,ntri,1), &
                                                    !    trilinear_index(ib_count,ntri,2), &
                                                    !    trilinear_index(ib_count,ntri,3)

                            if(trilinear_index(ib_count,ntri,2).le. n2)then
                                if(trilinear_index(ib_count,ntri,1).le. n1)then
                                    if(trilinear_index(ib_count,ntri,1).ge. 1)then
                                        if(trilinear_index(ib_count,ntri,3).le. n3)then
                                            if(trilinear_index(ib_count,ntri,3).ge. 1)then


                                                if( tipo(trilinear_index(ib_count,ntri,1), &
                                                    trilinear_index(ib_count,ntri,2), &
                                                    trilinear_index(ib_count,ntri,3)) &
                                                    == 0)then !1

                                                    check_box_solid = check_box_solid + 1
                                                    write(4500+myid,*)trilinear_index(ib_count,ntri,1), &
                                                        trilinear_index(ib_count,ntri,2), &
                                                        trilinear_index(ib_count,ntri,3)
                                                    write(4500+myid,*)trilinear_coef(ntri)

                                                    if(trilinear_coef(ntri) .gt. 1.d-7)then

                                                        write(4600+myid,*)trilinear_index(ib_count,ntri,1), &
                                                            trilinear_index(ib_count,ntri,2), &
                                                            trilinear_index(ib_count,ntri,3)
                                                        write(4600+myid,*)trilinear_coef(ntri)
                                                        write(4600+myid,*)i,j,k,ic,jc,kc
                                                        write(4600+myid,*)' '
                                                    end if


                                                end if

                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end do
                        if(check_box_solid .gt. 0)then
                            write(4500+myid,*)check_box_solid
                            write(4500+myid,*)i,j,k,ic,jc,kc
                            write(4500+myid,*)'--------'
                        end if

                        ! trilinear coefficent
                        tricoef(ib_count,1) = trilinear_coef_col1(ib_count)
                        tricoef(ib_count,2) = trilinear_coef_col2(ib_count)
                        tricoef(ib_count,3) = trilinear_coef_col3(ib_count)
                        tricoef(ib_count,4) = trilinear_coef_col4(ib_count)


                        ! trilinear index
                        trind(ib_count,1,1) = trilinear_index(ib_count,1,1)
                        trind(ib_count,2,1) = trilinear_index(ib_count,2,1)
                        trind(ib_count,3,1) = trilinear_index(ib_count,3,1)
                        trind(ib_count,4,1) = trilinear_index(ib_count,4,1)


                        trind(ib_count,1,2) = trilinear_index(ib_count,1,2)
                        trind(ib_count,2,2) = trilinear_index(ib_count,2,2)
                        trind(ib_count,3,2) = trilinear_index(ib_count,3,2)
                        trind(ib_count,4,2) = trilinear_index(ib_count,4,2)


                        trind(ib_count,1,3) = trilinear_index(ib_count,1,3)
                        trind(ib_count,2,3) = trilinear_index(ib_count,2,3)
                        trind(ib_count,3,3) = trilinear_index(ib_count,3,3)
                        trind(ib_count,4,3) = trilinear_index(ib_count,4,3)

                        ib_count = ib_count + 1

                    endif  !if tipo
                end do
            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(myid==0)then

            write(*,*)' -----------------------------------------'
            write(*,*)'| SECTION 3: geometry search --> finished |'
            write(*,*)' -----------------------------------------'

            if(contatore_mesh_rozza/=0 )then
                if(int(Numero_Celle_IB/contatore_mesh_rozza)<100)then
                    write(*,*)'too man2 projection-point results, maybe your mesh is not good'
                    write(*,*)'SUGGESTION 1: use fort.1200 to see which zones give problems'
                    write(*,*)'SUGGESTION 2: try to make structured mesh in that area'
                    write(*,*)'SUGGESTION 3: put more mesh points on edge'
                end if
            end if
        end if

        ! deallocate data not used an2more
        !deallocate(X)
        !deallocate(Coord_Boundary,Indici_Boundary,vertice_comune)
        deallocate(nodo_vicino_x_array,nodo_vicino_y_array,nodo_vicino_z_array)
        !deallocate(tipo)
        !deallocate(vertice_comune_prov)
        !deallocate(linea)
        !deallocate(num_punti,num_triangoli)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        return

    end subroutine projection_search

    subroutine prepare_communication()

        implicit none

        integer :: i,j,k,kp,l,np
        integer :: ierror
        integer :: status(MPI_STATUS_SIZE)


        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !     PREPARE COMUNICATION BETWEEN PROCS FOR IB
        allocate( stencil_left_snd(2*jx*jy,3))
        allocate(stencil_right_snd(2*jx*jy,3))
        allocate( stencil_left_rcv(2*jx*jy,3))
        allocate(stencil_right_rcv(2*jx*jy,3))

        allocate(tipo_spedito(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        stencil_left_rcv  = 0
        stencil_right_rcv = 0
        stencil_left_snd  = 0
        stencil_right_snd = 0

        ! my stencil requires node of the close proc
        ! therefore kparasta-1 and kparaend+1

        num_left_snd = 0
        num_right_snd = 0
        num_right_rcv = 0
        num_left_rcv = 0

        tipo_spedito = 0

        do l=1,num_ib
            do np = 1,4
                i=trind(l,np,1)
                j=trind(l,np,2)
                k=trind(l,np,3)


                if(tipo_spedito(i,j,k)==0)then

                    if(k==kparasta-1)then
                        num_left_rcv = num_left_rcv + 1
                        stencil_left_rcv(num_left_rcv,1) = i
                        stencil_left_rcv(num_left_rcv,2) = j
                        stencil_left_rcv(num_left_rcv,3) = k

                        tipo_spedito(i,j,k) = 1
                    end if


                    if(k==kparaend+1)then
                        num_right_rcv = num_right_rcv + 1
                        stencil_right_rcv(num_right_rcv,1) = i
                        stencil_right_rcv(num_right_rcv,2) = j
                        stencil_right_rcv(num_right_rcv,3) = k

                        tipo_spedito(i,j,k) = 1
                    end if

                end if ! tipo_spedito

            end do
        end do

        !      write(*,*)myid,'num left recv',num_left_rcv
        !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !      write(*,*)myid,'num right recv',num_right_rcv
        !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !     COMUNICATION OF INDEX AND COUNTER
        !
        !     what I send in stencil_left_send will be recived in stencil_right_recive
        !     in closer left proc
        !
        !     and the same for the right
        !
        !     what I send in stencil_right_send will be recived in stencil_left_recive
        !     in closer right proc
        !
        !     step 1, comunication of counter, how much I send/recive


        !     comunicate left the point left has to send right
        if(myid.ne.0)then
            call MPI_SSEND(num_left_rcv,1,MPI_INTEGER, &
                leftpe,tagls,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (ireq1,status,ierror)
        end if
        if(myid.ne.nproc-1)then
            call MPI_RECV(num_right_snd,1,MPI_INTEGER, &
                rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (ireq2,status,ierror)
        end if

        !     comunicate right the point right has to send left
        if(myid.ne.nproc-1)then
            call MPI_SSEND(num_right_rcv,1,MPI_INTEGER, &
                rightpe,tagrs,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (ireq3,status,ierror)
        end if
        if(myid.ne.0)then
            call MPI_RECV(num_left_snd,1,MPI_INTEGER, &
                leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (ireq4,status,ierror)
        end if

        !     if periodic
        if(kp==0)then
            if(myid.eq.0)then
                call MPI_SSEND(num_left_rcv,1,MPI_INTEGER, &
                    nproc-1,tagls,MPI_COMM_WORLD,ierror)
            !          call MPI_WAIT (ireq1,status,ierror)
            end if
            if(myid.eq.nproc-1)then
                call MPI_RECV(num_right_snd,1,MPI_INTEGER, &
                    0,tagrr,MPI_COMM_WORLD,status,ierror)
            !          call MPI_WAIT (ireq2,status,ierror)
            end if

            if(myid.eq.nproc-1)then
                call MPI_SSEND(num_right_rcv,1,MPI_INTEGER, &
                    0,tagrs,MPI_COMM_WORLD,ierror)
            !          call MPI_WAIT (ireq3,status,ierror)
            end if
            if(myid.eq.0)then
                call MPI_RECV(num_left_snd,1,MPI_INTEGER, &
                    nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            !          call MPI_WAIT (ireq4,status,ierror)
            end if
        else
            if(myid.eq.0)then
                num_left_snd = 0
                num_left_rcv = 0
            end if
            if(myid.eq.nproc-1)then
                num_right_rcv = 0
                num_right_snd = 0
            end if
        end if

        !     I know the counter! now I send the values

        !     send left the indeces left has to send right
        if(myid.ne.0)then
            call MPI_SSEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER, &
                leftpe,tagls,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (ireq1,status,ierror)
        end if
        if(myid.ne.nproc-1)then
            call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
                rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (ireq2,status,ierror)
        end if

        !     send right the indeces right has to send left
        if(myid.ne.nproc-1)then
            call MPI_SSEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER, &
                rightpe,tagrs,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (ireq3,status,ierror)
        end if
        if(myid.ne.0)then
            call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
                leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (ireq4,status,ierror)
        end if

        !     if periodic
        if(kp==0)then

            !       left to 0 is nproc-1
            if(myid.eq.0)then
                call MPI_SSEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,   &
                    nproc-1,tagls,MPI_COMM_WORLD,ierror)
            !          call MPI_WAIT (ireq1,status,ierror)
            end if
            if(myid.eq.nproc-1)then
                call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
                    0,tagrr,MPI_COMM_WORLD,status,ierror)
            !          call MPI_WAIT (ireq2,status,ierror)
            end if

            ! right to nproc-1 is 0
            if(myid.eq.nproc-1)then
                call MPI_SSEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER, &
                    0,tagrs,MPI_COMM_WORLD,ierror)
            !          call MPI_WAIT (ireq3,status,ierror)
            end if
            if(myid.eq.0)then
                call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
                    nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            !          call MPI_WAIT (ireq4,status,ierror)
            end if

            !       now myid=0 and myid=nproc-1 know the index they will recive,
            !       but they need to save
            !       on the border, so I change the k

            if(myid == nproc-1)stencil_right_snd(:,3)=jz
            if(myid == 0)      stencil_left_snd(:,3)=1

        end if

        return

    end subroutine prepare_communication

    subroutine readGeometry()

        use mysending

        implicit none

        character*16,parameter :: geometryFileName='bedParticles.dat'
        integer,parameter :: geometryFileID=1234

        integer :: n
        integer :: dummyInteger
        real :: dummyReal

        !-------------------------------------------------
        if (myid==0) then
            write(*,*) 'Reading particle file ',geometryFileName
        end if

        open(geometryFileID,file=geometryFileName)
        read(geometryFileID,*) totParticles

        allocate(sphereIndex(totParticles))
        allocate(radius(totParticles))
        allocate(radius2(totParticles))
        allocate(x_sphere(totParticles))
        allocate(y_sphere(totParticles))
        allocate(z_sphere(totParticles))
        allocate(x_force_sphere(totParticles))
        allocate(y_force_sphere(totParticles))
        allocate(z_force_sphere(totParticles))
        allocate(surface_sphere(totParticles))

        do n=1,totParticles
            read(geometryFileID,*) sphereIndex(n),dummyInteger, &
                radius(n),x_sphere(n),y_sphere(n),z_sphere(n), &
                dummyReal,dummyReal,dummyReal
            radius2(n)=radius(n)**2.0
            surface_sphere(n)=4.0*3.1428*radius(n)*radius(n)
        end do

        close(geometryFileID)

        !        if (myid==0) then
        !            write(*,*) 'Particles:'
        !            do n=1,totParticles
        !                write(*,*) radius(n),x_sphere(n),y_sphere(n),z_sphere(n)
        !            end do
        !        end if

        if (myid==0) then
            write(*,*) "Read particle file. Tot particles: ",totParticles
        end if

        return

    end subroutine readGeometry

end module ricerca_module


