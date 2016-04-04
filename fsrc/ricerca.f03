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

    use,intrinsic :: iso_c_binding
    use :: mysending
    use :: trilinear
    use :: scala3
    use :: myarrays_metri3
    use :: myarrays_ibm
    use :: mpi

    implicit none

    private

    integer,allocatable,public :: tipo(:,:,:),tipo2(:,:,:)
    integer,allocatable :: neighbor(:,:,:,:,:)
    logical,allocatable :: isValid(:,:,:,:)
    integer,parameter :: maxSolidIndex=6

    public :: set_ibm
    private :: particle_ricerca,phase_search,interface_search,projection_search,prepare_communication,rotationAle

contains

    subroutine set_ibm(ktime)

        implicit none

        integer,intent(in) :: ktime

        integer itiposta,itipoend
        integer jtiposta,jtipoend
        integer ktiposta,ktipoend
        integer i,j,k,ii,jj,kk

        ! -------------------------------------------------------------------------
        ! allocation: what needs to be allocated depends on whether bodyforce or particles are active

        ! any case tipo is needed
        if (ktime==0) then
            allocate(tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            allocate(tipo2(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            tipo(:,:,:)=0
            tipo2(:,:,:)=2
        end if



        ! now check if the ibm has to be applied
        if (bodyforce) then

            if (particles) then

                ! these are allocated at every time step to manage particles
                ! most of them are just debugging facilities, consider removing
                if (.not.allocated(solidIndex)) then
                    allocate(solidIndexSize(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(solidIndex(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,maxSolidIndex))
                    allocate(nodo_vicino_x_array(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(nodo_vicino_y_array(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(nodo_vicino_z_array(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(normalVectorX(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(normalVectorY(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(normalVectorZ(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(surfaceDistance(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
                    allocate(neighbor(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,6,3))
                    allocate(isValid(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,6))
                end if

                ! do the ricerca cycle (to be changed A LOT!)
                call particle_ricerca(ktime)

            else if (ktime==0) then
                ! this is the old ibm style, with loading of data from external files

                open(155,file='Celle_IB_indici.inp',status='old')
                open(156,file='Celle_Bloccate_Indici.inp',status='old')

                read(155,*)MN
                read(156,*)MP

                close(155)
                close(156)

                allocate(indici_CELLE_IB(MN,6))
                allocate(indici_celle_bloccate(MP,3))
                allocate(distanze_CELLE_IB(MN,3))
                allocate(dist_pp_ib(MN))
                allocate(dist_ib_parete(MN))
                allocate(proiezioni(MN,3))
                allocate(ustar(MN))

                allocate(tricoef(MN,4))
                allocate(trind(MN,4,3))


                allocate(rot(MN,3,3))
                allocate(rot_inverse(MN,3,3))

                call carico_immb()

                allocate(r_solid(nscal,num_solide))

                if (myid==0) then
                    write(*,*)'CHECK: call carico_immb----> OK'
                    write(*,*)'values for MN and MP: ',MN,MP
                end if

            end if


        else
            ! if no ibm set all the cells to fluid
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        tipo(i,j,k)=2 ! ATTENTION HERE, MODIFIED!!!!!!!!!!!
                    end do
                end do
            end do
        end if


        ! recompute tipo2 if first iteration or if particles
        if (ktime==0 .or. particles) then
            do k=kparasta,kparaend !1,jz
                do j=1,jy
                    do i=1,jx

                        itiposta=1
                        itipoend=1
                        jtiposta=1
                        jtipoend=1
                        ktiposta=1
                        ktipoend=1

                        if (i==1) itiposta=0
                        if (j==1) jtiposta=0
                        if (k==1) ktiposta=0

                        if (i==jx) itipoend=0
                        if (j==jy) jtipoend=0
                        if (k==jz) ktipoend=0

                        if (tipo(i,j,k)==2) then
                            do kk=k-ktiposta,k+ktipoend
                                do jj=j-jtiposta,j+jtipoend
                                    do ii=i-itiposta,i+itipoend

                                        if (tipo(ii,jj,kk)==1) then
                                            tipo2(i,j,k)=0    ! fluid close to ib
                                        end if

                                    end do
                                end do
                            end do
                        end if

                        if (tipo(i,j,k)==0) tipo2(i,j,k)=0
                        if (tipo(i,j,k)==1) tipo2(i,j,k)=0

                    end do
                end do
            end do
        end if

    end subroutine set_ibm

    subroutine particle_ricerca(ktime)

        implicit none

        integer,intent(in) :: ktime

        integer :: ierr

        real :: time1,time2,time

        ! ----------------------------------------------------

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call cpu_time(time1)
        call phase_search(ktime)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call cpu_time(time2)
        time = time2 - time1
        if (myid==0) then
            write(*,*)'| IBM 1: Phase search --> finished, time expired: ',time
        end if

        call cpu_time(time1)
        call interface_search()
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call cpu_time(time2)
        time = time2 - time1

        if (myid==0) then
            write(*,*)'| IBM 2: Projection search --> finished, time expired: ',time
        end if

        call cpu_time(time1)
        call projection_search()
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call cpu_time(time2)
        time = time2 - time1

        if (myid==0) then
            write(*,*)'| IBM 3: Trilinear search --> finished, time expired: ',time
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        call prepare_communication()

        return

    end subroutine particle_ricerca

    subroutine phase_search(ktime)

        use :: particle_module, only: totParticles,spherePosition,sphereRadius2

        implicit none

        integer,intent(in) :: ktime

        ! ----------------------------------------------------
        ! For solid/fluid/IBM research

        real :: isoparameter
        integer :: solid_count
        integer :: counterTot,reducedCounterTot
        integer :: counterFluid,contatore,num_solide_real
        integer :: reducedCounterFluid
        integer :: i,j,k,l,p,n,ni,nj,nk
        integer :: indexHere
        integer :: ierr

        ! SECTION 2: PHASES

        ! determine if a cell is solid, fluid or ib starting from zero or from a restart file
        ! (0,1,2 MEANS: fluid, solid, ib)


        !write(*,*) 'Proc ',myid,'scanning ',kparasta-deepl,' to ',kparaend+deepr


        ! total number of cells
        counterTot=n1*n2*(kparaend-kparasta+1)
        call MPI_ALLREDUCE(counterTot,reducedCounterTot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! set everything fluid
        tipo(:,:,:)=2 !0

        ! reset solid indices
        solidIndex(:,:,:,:)=0
        solidIndexSize(:,:,:)=0

        num_solide=0
        contatore=0
        if (myid==0) then
            write(*,*) 'Looking for solid nodes'
        end if
        do k=kparasta-deepl,kparaend+deepr
            do j=0,n2+1
                do i=0,n1+1
                    do p=1,totParticles
                        if (tipo(i,j,k)/=0) then !1
                            isoparameter=(xcd(i,j,k)-spherePosition(1,p))**2 + &
                                (ycd(i,j,k)-spherePosition(2,p))**2 + &
                                (zcd(i,j,k)-spherePosition(3,p))**2

                            if (isoparameter<sphereRadius2(p)) then


                                tipo(i,j,k)=0 !1
                                solidIndex(i,j,k,1)=p

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
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
        if (allocated(indici_CELLE_IB)) then
            deallocate(indici_CELLE_IB,indici_celle_bloccate)
            deallocate(distanze_CELLE_IB,dist_pp_ib,dist_ib_parete,proiezioni)
            deallocate(ustar,pressure_ib,shear_ib,caso_ib)
            deallocate(tricoef,trind,rot,rot_inverse,r_solid)
        end if

        allocate(indici_CELLE_IB(num_ib,6))
        allocate(indici_celle_bloccate(num_solide,3))
        allocate(distanze_CELLE_IB(num_ib,3))
        allocate(dist_pp_ib(num_ib))
        allocate(dist_ib_parete(num_ib))
        allocate(proiezioni(num_ib,3))
        allocate(ustar(num_ib),pressure_ib(num_ib,3),shear_ib(num_ib,3),caso_ib(num_ib))

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

    subroutine interface_search()

        use :: particle_module, only: totParticles,spherePosition,sphereRadius

        implicit none

        ! normal to surface vector components
        real :: normalVectorXHere, normalVectorYHere, normalVectorZHere
        ! IP_BP distance
        real :: distanceNormHere
        integer :: i,j,k,n
        integer :: indexHere, solidIndexSizeHere

        !---------------------------------------------------------------------------------------------------------------------------------------------
        ! SECTION 2: search the projection point IP on the body, through the normal from the IB to the body
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
                    if (tipo(i,j,k)==1) then !2
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
                                write(*,*) 'Problem: solid index=',indexHere
                                call exit(0)
                            end if
                            normalVectorXHere=xcd(i,j,k)-spherePosition(1,indexHere)
                            normalVectorYHere=ycd(i,j,k)-spherePosition(2,indexHere)
                            normalVectorZHere=zcd(i,j,k)-spherePosition(3,indexHere)
                            ! normalize it
                            distanceNormHere=sqrt(normalVectorXHere**2.0+normalVectorYHere**2.0+normalVectorZHere**2.0)
                            if (distanceNormHere==0.0) then
                                write(*,*) 'XP=',xcd(i,j,k),'; XS=',spherePosition(1,indexHere),'; DX=',normalVectorXHere
                            end if
                            normalVectorXHere=normalVectorXHere/distanceNormHere
                            normalVectorYHere=normalVectorYHere/distanceNormHere
                            normalVectorZHere=normalVectorZHere/distanceNormHere
                            ! shift normal
                            normalVectorX(i,j,k)=normalVectorX(i,j,k)+normalVectorXHere/solidIndexSizeHere
                            normalVectorY(i,j,k)=normalVectorY(i,j,k)+normalVectorYHere/solidIndexSizeHere
                            normalVectorZ(i,j,k)=normalVectorZ(i,j,k)+normalVectorZHere/solidIndexSizeHere
                            ! update distance from surface
                            surfaceDistance(i,j,k)=surfaceDistance(i,j,k)+ &
                                (distanceNormHere-sphereRadius(indexHere))/solidIndexSizeHere
                            if (surfaceDistance(i,j,k)<=0) then
                                write(*,*) 'Problem: negative surface distance: ',surfaceDistance(i,j,k),&
                                    '; solidIndexSizeHere = ',solidIndexSizeHere,'; n = ',n
                                write(*,*) 'sphere = ',indexHere,'; radius = ',sphereRadius(indexHere),'; norm = ',distanceNormHere
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

    subroutine projection_search()

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
        use mpi

        implicit none

        real,allocatable :: trilinear_coef_col1(:),trilinear_coef_col2(:),trilinear_coef_col3(:),trilinear_coef_col4(:)
        real trilinear_coef(4)

        integer :: ierr
        integer, allocatable :: trilinear_index(:,:,:)

        integer, parameter ::  iread_linea = 0
        real, parameter ::  calibro_iter = 3.5
        real, parameter ::  calibro_iter_inc = 0.2
        real, parameter ::  calibro_punto = 100.
        real, parameter ::  calibro_punto_inc = 5.
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


        real,parameter :: tol=1.0e-14
        real,parameter :: tolC=1.0e-14
        !---------------------------------------------------------------------------------

        !index
        integer :: i,j,k
        ! V indices
        integer :: ic,jc,kc

        ! active procedure
        logical :: trovato_vicino

        ! incremental
        integer :: contatore_mesh_rozza
        ! cell center coordinates
        real :: sx,sy,sz
        ! V coordinates
        ! real :: mx,my,mz
        !PP coordinates
        real :: nodo_proiezione_x, nodo_proiezione_y, nodo_proiezione_z
        !IP coordinates
        real :: nodo_vicino_x, nodo_vicino_y, nodo_vicino_z

        ! distance between PP and V
        real :: distanza_x, distanza_y,distanza_z
        real :: distanza_np_ib,distanza_ib_parete

        !rotation
        real :: rot_ric(3,3),rot_inverse_ric(3,3)

        integer,parameter :: formatted=0

        integer ib_count

        integer ntri,check_box_solid
        real count_trilinear_coef
        integer icount_tri

        integer pp_found
        integer iref,jref,kref
        real x0,y0,z0,FF,GG,HH
        logical trovato_pp
        real val_trilinear

        integer,parameter :: icheck_print=0

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
            !if (myid==0)write(*,*)'internal loop PP:',k,'of',n3/nproc
            do j = 1,n2
                do i = 1,n1
                    if (tipo(i,j,k) == 1) then !2

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


                        ! AAA da chiamare una volta sola fuori dal loop
                        call triangles_mesh


                        !AAA se trovato lo metto qui false e non dentro poi quando e' trovato non passa nelle altre
                        ! sub ma stare attenti alla parte finale della sub
                        pp_found = 0
                        iref = i-1
                        call pp_on_meshbox_diri(iref,j,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if (trovato_pp)pp_found = pp_found + 1
                        if (icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if (myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        iref = i+1
                        call pp_on_meshbox_diri(iref,j,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if (trovato_pp)pp_found = pp_found + 1
                        if (icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if (myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        jref = j-1
                        call pp_on_meshbox_dirj(jref,i,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if (trovato_pp)pp_found = pp_found + 1
                        if (icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if (myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        jref = j+1
                        call pp_on_meshbox_dirj(jref,i,k,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if (trovato_pp)pp_found = pp_found + 1
                        if (icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if (myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        kref = k-1
                        call pp_on_meshbox_dirk(kref,i,j,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if (trovato_pp)pp_found = pp_found + 1
                        if (icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                        !if (myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        kref = k+1
                        call pp_on_meshbox_dirk(kref,i,j,sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z, &
                            nodo_proiezione_x,nodo_proiezione_y,nodo_proiezione_z, &
                            ib_count,trovato_pp,trilinear_index,FF,GG,HH,x0,y0,z0,n1,n2,n3,icheck_print)
                        if (trovato_pp)pp_found = pp_found + 1
                        if (icheck_print == 1)write(6700+myid,*)'trovato',trovato_pp
                             !if (myid==1)write(*,*)'trovato',trovato_pp,i,j,k

                        if (pp_found == 1) then
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
                                if (trilinear_coef(ntri) > val_trilinear) then
                                    val_trilinear = trilinear_coef(ntri)

                                    ic = trilinear_index(ib_count,ntri,1)
                                    jc = trilinear_index(ib_count,ntri,2)
                                    kc = trilinear_index(ib_count,ntri,3)

                                    if (ic == n1+1)ic = n1
                                    if (ic == 0   )ic = 1

                                    if (jc == n2+1)jc = n2
                                    if (jc == 0   )jc = 1

                                    if (kc == n3+1)kc = n3
                                    if (kc == 0   )kc = 1
                                    !write(*,*)ntri,ib_count,'ntri',i,j,k,ic,jc,kc

                                    distanza_x=nodo_proiezione_x-xcd(ic,jc,kc)
                                    distanza_y=nodo_proiezione_y-ycd(ic,jc,kc)
                                    distanza_z=nodo_proiezione_z-zcd(ic,jc,kc)
                                end if
                            end do

                            !trilinear coef for solid = 0

                            do ntri = 1,4
                                if (trilinear_index(ib_count,ntri,2).le. n2) then
                                    if (trilinear_index(ib_count,ntri,2).ge. 1 ) then
                                        if (trilinear_index(ib_count,ntri,1).le. n1) then
                                            if (trilinear_index(ib_count,ntri,1).ge. 1 ) then
                                                if (trilinear_index(ib_count,ntri,3).le. n3) then
                                                    if (trilinear_index(ib_count,ntri,3).ge. 1 ) then
                                                        if ( tipo(trilinear_index(ib_count,ntri,1), &
                                                            trilinear_index(ib_count,ntri,2), &
                                                            trilinear_index(ib_count,ntri,3)) == 0) then !1

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


                            !write(*,'(A,I0.1,A,I0.1,I0.1,I0.1,A,I0.1,A,ES11.2E3,ES11.2E3,ES11.2E3,A)') &
                             !   'ID = ',myid,', problem on ib ',i,j,k,' in finding PP (',pp_found, &
                             !   ') vicino=(',sx,sy,sz,')'
                        end if

                        count_trilinear_coef = 0.
                        do ntri=1,4
                            count_trilinear_coef = count_trilinear_coef + trilinear_coef(ntri)
                        end do
                        if (count_trilinear_coef == 0.) then
                            icount_tri = icount_tri + 1
                           !write(*,*)'coef nulli',i,j,k
                        end if

                        if (icheck_print == 1) then
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

                        if (icheck_print == 1) then
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

                        !----------------------------------------------------------------------------------
                        ! ALESSANDRO : THIS IS THE OLD STYLE. NEW STYLE IS IN INTERFACE_SEARCH
                        ! ROTAZION MATRIX
                        !                        call rotation2(sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z,rot_ric,rot_inverse_ric) !,alfa,beta,gamma,vertice)
                        !
                        !                        ! rotation matrix (to construct tangential and normal velocity)
                        !                        rot(ib_count,1,1)=rot_ric(1,1)
                        !                        rot(ib_count,1,2)=rot_ric(1,2)
                        !                        rot(ib_count,1,3)=rot_ric(1,3)
                        !
                        !                        rot(ib_count,2,1)=rot_ric(2,1)
                        !                        rot(ib_count,2,2)=rot_ric(2,2)
                        !                        rot(ib_count,2,3)=rot_ric(2,3)
                        !
                        !                        rot(ib_count,3,1)=rot_ric(3,1)
                        !                        rot(ib_count,3,2)=rot_ric(3,2)
                        !                        rot(ib_count,3,3)=rot_ric(3,3)
                        !
                        !                        ! inverse rotation matrix
                        !                        rot_inverse(ib_count,1,1)=rot_inverse_ric(1,1)
                        !                        rot_inverse(ib_count,1,2)=rot_inverse_ric(1,2)
                        !                        rot_inverse(ib_count,1,3)=rot_inverse_ric(1,3)
                        !
                        !                        rot_inverse(ib_count,2,1)=rot_inverse_ric(2,1)
                        !                        rot_inverse(ib_count,2,2)=rot_inverse_ric(2,2)
                        !                        rot_inverse(ib_count,2,3)=rot_inverse_ric(2,3)
                        !
                        !                        rot_inverse(ib_count,3,1)=rot_inverse_ric(3,1)
                        !                        rot_inverse(ib_count,3,2)=rot_inverse_ric(3,2)
                        !                        rot_inverse(ib_count,3,3)=rot_inverse_ric(3,3)

                        call rotation2(sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z,rot_ric,rot_inverse_ric) !,alfa,beta,gamma,vertice)
                        !call rotationAle(sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z,rot_ric,rot_inverse_ric) !,alfa,beta,gamma,vertice)

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


                        ! inverse rotation matrix
                        rot_inverse(ib_count,1,1)=rot_inverse_ric(1,1)
                        rot_inverse(ib_count,1,2)=rot_inverse_ric(1,2)
                        rot_inverse(ib_count,1,3)=rot_inverse_ric(1,3)

                        rot_inverse(ib_count,2,1)=rot_inverse_ric(2,1)
                        rot_inverse(ib_count,2,2)=rot_inverse_ric(2,2)
                        rot_inverse(ib_count,2,3)=rot_inverse_ric(2,3)

                        rot_inverse(ib_count,3,1)=rot_inverse_ric(3,1)
                        rot_inverse(ib_count,3,2)=rot_inverse_ric(3,2)
                        rot_inverse(ib_count,3,3)=rot_inverse_ric(3,3)


                        !----------------------------------------------------------------------------------

                        !first control, solid inside the box?
                        check_box_solid = 0
                        do ntri = 1,4
                                                    !write(*,*)    trilinear_index(ib_count,ntri,1), &
                                                    !    trilinear_index(ib_count,ntri,2), &
                                                    !    trilinear_index(ib_count,ntri,3)

                            if (trilinear_index(ib_count,ntri,2).le. n2) then
                                if (trilinear_index(ib_count,ntri,1).le. n1) then
                                    if (trilinear_index(ib_count,ntri,1).ge. 1) then
                                        if (trilinear_index(ib_count,ntri,3).le. n3) then
                                            if (trilinear_index(ib_count,ntri,3).ge. 1) then


                                                if ( tipo(trilinear_index(ib_count,ntri,1), &
                                                    trilinear_index(ib_count,ntri,2), &
                                                    trilinear_index(ib_count,ntri,3)) &
                                                    == 0) then !1

                                                    check_box_solid = check_box_solid + 1
                                                    write(4500+myid,*)trilinear_index(ib_count,ntri,1), &
                                                        trilinear_index(ib_count,ntri,2), &
                                                        trilinear_index(ib_count,ntri,3)
                                                    write(4500+myid,*)trilinear_coef(ntri)

                                                    if (trilinear_coef(ntri) .gt. 1.d-7) then

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
                        if (check_box_solid .gt. 0) then
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

        if (myid==0) then

            !write(*,*)' -----------------------------------------'
            !write(*,*)'| SECTION 3: geometry search --> finished |'
            !write(*,*)' -----------------------------------------'

            if (contatore_mesh_rozza/=0 ) then
                if (int(Numero_Celle_IB/contatore_mesh_rozza)<100) then
                    write(*,*)'too many projection-point results, maybe your mesh is not good'
                    write(*,*)'SUGGESTION 1: use fort.1200 to see which zones give problems'
                    write(*,*)'SUGGESTION 2: try to make structured mesh in that area'
                    write(*,*)'SUGGESTION 3: put more mesh points on edge'
                end if
            end if
        end if

        ! deallocate data not used an2more
        !deallocate(X)
        !deallocate(Coord_Boundary,Indici_Boundary,vertice_comune)
        !deallocate(nodo_vicino_x_array,nodo_vicino_y_array,nodo_vicino_z_array)
        !deallocate(tipo)
        !deallocate(vertice_comune_prov)
        !deallocate(linea)
        !deallocate(num_punti,num_triangoli)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        return

    end subroutine projection_search

    subroutine prepare_communication()

        use period, only: kp

        implicit none

        integer :: i,j,k,l,np
        integer :: ierror
        integer :: status(MPI_STATUS_SIZE)

        !     PREPARE COMUNICATION BETWEEN PROCS FOR IB

        if (allocated(stencil_left_snd)) then
            deallocate(stencil_left_snd,stencil_right_snd)
            deallocate(stencil_left_rcv,stencil_right_rcv)
            deallocate(tipo_spedito)
        end if

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

                if (tipo_spedito(i,j,k)==0) then

                    if (k==kparasta-1) then
                        num_left_rcv = num_left_rcv + 1
                        stencil_left_rcv(num_left_rcv,1) = i
                        stencil_left_rcv(num_left_rcv,2) = j
                        stencil_left_rcv(num_left_rcv,3) = k

                        tipo_spedito(i,j,k) = 1
                    end if


                    if (k==kparaend+1) then
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
        if (myid/=0) then
            call MPI_SSEND(num_left_rcv,1,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(num_right_snd,1,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     comunicate right the point right has to send left
        if (myid/=nproc-1) then
            call MPI_SSEND(num_right_rcv,1,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(num_left_snd,1,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if (kp==0) then
            if (myid==0) then
                call MPI_SSEND(num_left_rcv,1,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(num_right_snd,1,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            if (myid==nproc-1) then
                call MPI_SSEND(num_right_rcv,1,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
            end if
            if (myid==0) then
                call MPI_RECV(num_left_snd,1,MPI_INTEGER,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            end if
        else
            if (myid==0) then
                num_left_snd = 0
                num_left_rcv = 0
            end if
            if (myid==nproc-1) then
                num_right_rcv = 0
                num_right_snd = 0
            end if
        end if

        !     I know the counter! now I send the values

        !     send left the indeces left has to send right
        if (myid/=0) then
            call MPI_SSEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     send right the indeces right has to send left
        if (myid/=nproc-1) then
            call MPI_SSEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if (kp==0) then

            !       left to 0 is nproc-1
            if (myid==0) then
                call MPI_SSEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            ! right to nproc-1 is 0
            if (myid==nproc-1) then
                call MPI_SSEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
            end if
            if (myid==0) then
                call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            end if

            !       now myid=0 and myid=nproc-1 know the index they will recive,
            !       but they need to save
            !       on the border, so I change the k

            if (myid == nproc-1)stencil_right_snd(:,3)=jz
            if (myid == 0)      stencil_left_snd(:,3)=1

        end if

        return

    end subroutine prepare_communication

    subroutine rotationAle(sx,sy,sz,ipx,ipy,ipz,rot,rot_inverse)
        !-----------------------------------------------------------------------
        ! This is specific for the IBM, the complexity of rotation2 is reduced

        implicit none

        real,intent(in) :: sx,sy,sz
        real,intent(in) :: ipx,ipy,ipz

        real,intent(out) :: rot(:,:),rot_inverse(:,:)

        ! declaration
        ! normal vector
        real :: nx,ny,nz
        ! tangent vector 1
        real :: tx,ty,tz
        ! tangent vector 2 (binormal)
        real :: bx,by,bz

        real :: norm_t,norm_n,norm_b

        real,parameter :: pi = acos(-1.)
        real,parameter :: toll=1.0e-5

        !-----------------------------------------------------------------------
        ! normal vector n
        nx=sx-ipx
        ny=sy-ipy
        nz=sz-ipz

        ! norm of n
        norm_n=sqrt(nx*nx+ny*ny+nz*nz)

        ! normalization of n
        nx=nx/norm_n
        ny=ny/norm_n
        nz=nz/norm_n

        ! recompute norm
        norm_n=sqrt(nx*nx+ny*ny+nz*nz)

        if (norm_n<1.0-toll .or. norm_n>1.0+toll) then
            write(*,*) 'Very bad error, norm fucked up'
        end if

        ! we try to align t1 or t2 to the x axes

        if (abs(nx)>toll) then
            tz=1/sqrt(1+(nz/nx)**2)
            tx=-1.0*tz*nz/nx
            ty=0.0
        else if (abs(nz)>toll) then
            ty=1/sqrt(1+(ny/nz)**2)
            tz=-1.0*ty*ny/nz
            tx=0.0
        else if (abs(ny)>toll) then
            tx=1/sqrt(1+(nx/ny)**2)
            ty=-1.0*tx*nx/ny
            tz=0.0
        else
            write(*,*) 'Problem in rotation routine, t computation. n_x=',nx,' ny=',ny,' nz=',nz
        end if

        ! norm of t
        norm_t=sqrt(tx*tx+ty*ty+tz*tz)

        ! compute secont tangent vector with cross product
        bx=ny*tz-nz*ty
        by=nz*tx-nx*tz
        bz=nx*ty-ny*tx

        ! norm of t
        norm_b=sqrt(bx*bx+by*by+bz*bz)

        if (norm_b<1.0-toll .or. norm_b>1.0+toll .or. &
            norm_t<1.0-toll .or. norm_t>1.0+toll) then
            write(*,*) 'Problem in rotation routine, norm_b=',norm_b,' norm_t=',norm_t
        end if

        ! rotation matrix construction
        rot(1,1)=tx
        rot(1,2)=ty
        rot(1,3)=tz

        rot(2,1)=bx
        rot(2,2)=by
        rot(2,3)=bz

        rot(3,1)=nx
        rot(3,2)=ny
        rot(3,3)=nz

        ! inverse is just the transpose
        rot_inverse(1,1)=rot(1,1)
        rot_inverse(1,2)=rot(2,1)
        rot_inverse(1,3)=rot(3,1)

        rot_inverse(2,1)=rot(1,2)
        rot_inverse(2,2)=rot(2,2)
        rot_inverse(2,3)=rot(3,2)

        rot_inverse(3,1)=rot(1,3)
        rot_inverse(3,2)=rot(2,3)
        rot_inverse(3,3)=rot(3,3)

        return

    end subroutine rotationAle

end module ricerca_module

