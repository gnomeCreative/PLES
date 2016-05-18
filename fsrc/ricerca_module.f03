module ricerca_module

    ! release 21-11-2009
    ! Federico Roman
    ! modified by Alessandro Leonardi, starting December 2015

    ! Pre-processing program to use the Immersed Boundary Method as in
    ! Roman,Napoli,Armenio,Milici 2008

    ! GRID: non staggered


    use,intrinsic :: iso_c_binding
    use mysending
    use scala3
    use myarrays_metri3
    use ibm_module
    use mpi

    implicit none

    private

    integer,allocatable,public :: tipo(:,:,:),tipo2(:,:,:)
    ! Identifier for the solid object the IB belongs to
    integer,allocatable :: solidIndex(:,:,:,:)
    ! Identifier for the solid object the IB belongs to
    integer,allocatable :: solidIndexSize(:,:,:)
    integer,allocatable :: neighbor(:,:,:,:,:)
    logical,allocatable :: isValid(:,:,:,:)
    integer,parameter :: maxSolidIndex=6

    public :: set_ibm

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
            if (bodyforce .and. particles) then
                allocate(neighbor(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,6,3))
                allocate(isValid(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,6))
                allocate(solidIndex(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,maxSolidIndex))
                allocate(solidIndexSize(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            end if
        end if

        ! now check if the ibm has to be applied
        if (bodyforce) then

            if (particles) then

                ! do the ricerca cycle
                call particle_ricerca(ktime)

            else if (ktime==0) then
                ! this is the old ibm style, with loading of data from external files

                open(155,file='Celle_IB_indici.inp',status='old')
                open(156,file='Celle_Bloccate_Indici.inp',status='old')

                read(155,*)num_ib
                read(156,*)num_solide

                close(155)
                close(156)

                allocate(indici_CELLE_IB(num_ib,6))
                allocate(indici_celle_bloccate(num_solide,3))
                !allocate(distanze_CELLE_IB(MN,3))
                allocate(dist_pp_ib(num_ib))
                allocate(dist_ib_parete(num_ib))
                allocate(proiezioni(num_ib,3))
                allocate(ustar(num_ib))

                allocate(tricoef(num_ib,4))
                allocate(trind(num_ib,4,3))


                allocate(rot(num_ib,3,3))
                allocate(rot_inverse(num_ib,3,3))

                call carico_immb(tipo)

                !allocate(r_solid(nscal,num_solide))

                if (myid==0) then
                    write(*,*)'CHECK: call carico_immb----> OK'
                    write(*,*)'values for num_ib and MP: ',num_ib,num_solide
                end if

            end if


        else
            ! if no ibm set all the cells to fluid
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        tipo(i,j,k)=2
                    end do
                end do
            end do
        end if


        ! compute tipo2 if first iteration OR if particles
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
        integer :: solid_count,ib_count
        integer :: counterTot,reducedCounterTot
        integer :: counterFluid,contatore,num_solide_real
        integer :: reducedCounterFluid
        integer :: i,j,k,p,m,ni,nj,nk
        integer :: indexHere
        integer :: ierr

        ! SECTION 2: PHASES

        ! determine if a cell is solid, fluid or ib starting from zero or from a restart file
        ! (0,1,2 MEANS: fluid, solid, ib)


        !write(*,*) 'Proc ',myid,'scanning ',kparasta-deepl,' to ',kparaend+deepr

        ! compute neighbor matrix, only at first time step
        if (ktime==0) then
            if (myid==0) then
                write(*,*) 'Building up neighbor matrix'
            end if

            ! compute neighbors
            isValid=.false.
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        do m=1,6
                            neighbor(i,j,k,m,1)=i
                            neighbor(i,j,k,m,2)=j
                            neighbor(i,j,k,m,3)=k
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
        end if


        ! total number of cells
        counterTot=n1*n2*(kparaend-kparasta+1)
        call MPI_ALLREDUCE(counterTot,reducedCounterTot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! set everything fluid
        tipo(:,:,:)=2 !0

        ! reset solid indices
        solidIndex(:,:,:,:)=0 ! consider -1
        solidIndexSize(:,:,:)=0 ! consider -1

        num_solide=0
        contatore=0
        if (myid==0) then
            write(*,*) 'Looking for solid nodes'
        end if
        ! nodes inside a particle are treated as solid
        do k=kparasta-deepl,kparaend+deepr
            do j=0,n2+1
                do i=0,n1+1
                    do p=1,totParticles
                        if (tipo(i,j,k)/=0) then
                            isoparameter=(xcd(i,j,k)-spherePosition(1,p))**2 + &
                                (ycd(i,j,k)-spherePosition(2,p))**2 + &
                                (zcd(i,j,k)-spherePosition(3,p))**2

                            if (isoparameter<sphereRadius2(p)) then

                                tipo(i,j,k)=0
                                solidIndex(i,j,k,1)=p

                                if (i>=1 .and. i<=n1 .and. j>=1 .and. j<=n2) then
                                    num_solide=num_solide+1
                                    if (k>kparaend .or. k<kparasta) then
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
                        do m = 1,6
                            if (isValid(i,j,k,m)) then
                                ni=neighbor(i,j,k,m,1)
                                nj=neighbor(i,j,k,m,2)
                                nk=neighbor(i,j,k,m,3)
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

        call MPI_ALLREDUCE(num_solide,numero_celle_bloccate,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(num_solide_real,numero_celle_bloccate_real,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(counterFluid,reducedCounterFluid,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(num_ib,numero_celle_IB,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

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

        ! now we can allocate variables that are defined in this subroutine
        if (allocated(indici_CELLE_IB)) then
            deallocate(indici_CELLE_IB,indici_celle_bloccate)
            deallocate(position_ib)
            deallocate(indexsize_ib,index_ib)
        end if

        allocate(indici_CELLE_IB(num_ib,6))
        allocate(indici_celle_bloccate(num_solide,3))
        allocate(position_ib(num_ib,3))
        allocate(index_ib(num_ib,maxSolidIndex),indexsize_ib(num_ib))

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

        ib_count=0
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1

                    if (tipo(i,j,k)==1) then

                        ib_count=ib_count+1

                        indici_CELLE_IB(ib_count,1)=i
                        indici_CELLE_IB(ib_count,2)=j
                        indici_CELLE_IB(ib_count,3)=k

                        position_ib(ib_count,1)=xcd(i,j,k)
                        position_ib(ib_count,2)=ycd(i,j,k)
                        position_ib(ib_count,3)=zcd(i,j,k)

                        indexsize_ib(ib_count)=solidIndexSize(i,j,k)

                        do m=1,indexsize_ib(ib_count)
                            index_ib(ib_count,m)=solidIndex(i,j,k,m)
                        end do

                    end if

                end do
            end do
        end do

        ! check if everything went ok so far
        if (solid_count/=num_solide) then
            write(*,*) 'PROBLEM Proc: ',myid,' solid_count=',solid_count,', num_solide=',num_solide
            call exit(0)
        end if
        if (ib_count/=num_ib) then
            write(*,*) 'PROBLEM Proc: ',myid,' ib_count=',ib_count,', num_ib=',num_ib
            call exit(0)
        end if

        !MP=num_solide

        return

    end subroutine phase_search

    subroutine interface_search()

        use particle_module, only: totParticles,spherePosition,sphereRadius,sphereVelocity,sphereSpin

        implicit none

        ! normal to surface vector components
        real :: normalHere(3)
        real :: surfVelHere(3)
        real :: surfDistHere
        real :: leverHere(3)
        ! IP_BP distance
        real :: distanceNormHere
        ! iterators
        integer :: m,l
        integer :: indexHere,solidIndexSizeHere

        ! ----------------------------------------------------------------------------

        ! allocate variables that will be defined here
        if (allocated(proiezioni)) then
            !deallocate(distanze_CELLE_IB)
            deallocate(dist_ib_parete)
            deallocate(normalVector)
            deallocate(proiezioni)
            deallocate(surfVel)
        end if

        !allocate(distanze_CELLE_IB(num_ib,3))
        allocate(dist_ib_parete(num_ib))
        allocate(normalVector(num_ib,3))
        allocate(proiezioni(num_ib,3))
        allocate(surfVel(num_ib,3))

        dist_ib_parete(:)=0.0
        normalVector(:,:)=0.0
        proiezioni(:,:)=0.0
        surfVel(:,:)=0.0

        !---------------------------------------------------------------------------------------------------------------------------------------------
        ! SECTION 2: search the projection point IP on the body, through the normal from the IB to the body


        do l=1,num_ib

            solidIndexSizeHere=indexsize_ib(l)

            ! we start from the position of the node itself, then we shift it later
            proiezioni(l,:)=position_ib(l,:)
            if (solidIndexSizeHere<=0) then
                write(*,*) 'Problem: negative index size'
                call exit(0)
            end if

            ! iteratively compute normal and surface distance, averaging on the number of neighbors [solidIndex(l,m)]
            do m=1,solidIndexSizeHere
                ! compute vector normal to surface
                indexHere=index_ib(l,m)
                ! check if there is a problem with the indices
!                if (indexHere<=0 .or. indexHere>totParticles) then
!                    write(*,*) 'Problem: solid index=',indexHere
!                    call exit(0)
!                end if
                normalHere(:)=position_ib(l,:)-spherePosition(:,indexHere)
                ! normalize it
                distanceNormHere=norm2(normalHere)
!                if (distanceNormHere==0.0) then
!                    write(*,*) 'XP=',position_ib(l,1),'; XS=',spherePosition(1,indexHere),'; DX=',normalVectorHere(1)
!                end if
                normalHere(:)=normalHere(:)/distanceNormHere
                ! compute surface distance
                surfDistHere=distanceNormHere-sphereRadius(indexHere)
                ! compute surface velocity
                leverHere=surfDistHere*normalHere(:)
                surfVelHere(:)=sphereVelocity(:,indexHere)+cross(leverHere,sphereSpin(:,indexHere))
                ! shift normal
                normalVector(l,:)=normalVector(l,:)+normalHere(:)/solidIndexSizeHere
                ! update distance from surface
                dist_ib_parete(l)=dist_ib_parete(l)+surfDistHere/solidIndexSizeHere
                ! update surface velocity
                surfVel(l,:)=surfVel(l,:)+surfVelHere(:)/solidIndexSizeHere
                ! check if everything is all right
!                if (dist_ib_parete(l)<=0.0) then
!                    write(*,*) 'Problem: negative surface distance: ',dist_ib_parete(l),&
!                        '; solidIndexSizeHere = ',solidIndexSizeHere,'; m = ',m
!                    write(*,*) 'sphere = ',indexHere,'; radius = ',sphereRadius(indexHere),'; norm = ',distanceNormHere
!                    call exit(0)
!                end if
            end do

            ! shift interface point IP
            proiezioni(l,:)=proiezioni(l,:)-normalVector(l,:)*dist_ib_parete(l)

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

        use trilinear
        use geometricRoutines

        implicit none

        real :: trilinear_coef(4)

        !index
        integer :: l,i,j,k
        ! V indices
        integer :: v_indices(3)
        integer :: ic,jc,kc
        ! Ip,IB,PP coordinates
        real :: ib_position(3),ip_position(3),pp_position(3)

        !rotation
        real :: rot_ric(3,3),rot_inverse_ric(3,3)

        integer :: ntri
        real :: count_trilinear_coef

        integer :: iref,jref,kref
        real :: x0,y0,z0,FF,GG,HH
        logical :: trovato_pp
        real :: val_trilinear

        integer,parameter :: icheck_print=0

        ! ----------------------------------------------------------------------
        ! Allocate variables that will be defined here
        if (allocated(tricoef)) then
            deallocate(dist_pp_ib)
            deallocate(tricoef)
            deallocate(trind)
            deallocate(rot)
            deallocate(rot_inverse)
        end if

        ! move to next step
        allocate(dist_pp_ib(num_ib))
        allocate(tricoef(num_ib,4))
        allocate(trind(num_ib,4,3))
        allocate(rot(num_ib,3,3))
        allocate(rot_inverse(num_ib,3,3))

        tricoef(:,:)=-999999.0
        trind(:,:,:)=-999999

        !-----------------------------------------------------------------------
        ! PROJECTION POINT PP SEARCHIG PROCEDURE
        !-----------------------------------------------------------------------

        ! initialization of node coordination for trilinear
        call triangles_mesh()

        do l=1,num_ib

            ! position of current IB point
            ib_position(:)=position_ib(l,:)

            ! position of current IP point
            ip_position(:)=proiezioni(l,:)

            i=indici_CELLE_IB(l,1)
            j=indici_CELLE_IB(l,2)
            k=indici_CELLE_IB(l,3)

            !----------------------------------------------------------------------------------
            !search closest fluid node to PP

            ! initialize current box to IB index, and coefficients to zero
            v_indices(1)=i
            v_indices(2)=j
            v_indices(3)=k
            trilinear_coef(:)=0.0

            ! coefficients of line passing through IP and  IB
            call line_exp2par_3d(ib_position(1),ib_position(2),ib_position(3), &
                ip_position(1),ip_position(2),ip_position(3),FF,GG,HH,x0,y0,z0)

            ! start search for trilinear box
            trovato_pp=.false.

            if (.not.trovato_pp) then
                iref = i-1
                call pp_on_meshbox(1,iref,j,k,ib_position,ip_position, &
                    pp_position,l,trovato_pp,trind,FF,GG,HH,x0,y0,z0,num_ib)
            end if

            if (.not.trovato_pp) then
                iref = i+1
                call pp_on_meshbox(1,iref,j,k,ib_position,ip_position, &
                    pp_position,l,trovato_pp,trind,FF,GG,HH,x0,y0,z0,num_ib)
            end if

            if (.not.trovato_pp) then
                jref = j-1
                call pp_on_meshbox(2,i,jref,k,ib_position,ip_position, &
                    pp_position,l,trovato_pp,trind,FF,GG,HH,x0,y0,z0,num_ib)
            end if

            if (.not.trovato_pp) then
                jref = j+1
                call pp_on_meshbox(2,i,jref,k,ib_position,ip_position, &
                    pp_position,l,trovato_pp,trind,FF,GG,HH,x0,y0,z0,num_ib)
            end if

            if (.not.trovato_pp) then
                kref = k-1
                call pp_on_meshbox(3,i,j,kref,ib_position,ip_position, &
                    pp_position,l,trovato_pp,trind,FF,GG,HH,x0,y0,z0,num_ib)
            end if

            if (.not.trovato_pp) then
                kref = k+1
                call pp_on_meshbox(3,i,j,kref,ib_position,ip_position, &
                    pp_position,l,trovato_pp,trind,FF,GG,HH,x0,y0,z0,num_ib)
            end if

            if (trovato_pp) then ! Alessandro: it used to be pp_found == 1 (why??)

                call find_trilinear_coef2(l,trind,trilinear_coef, &
                    pp_position,icheck_print,i,j,k,num_ib)

                dist_pp_ib(l)=norm2(pp_position(:)-ib_position(:))

                tricoef(l,:)=trilinear_coef(:)

                ! v node
                val_trilinear=0.0
                do ntri=1,4
                    if (trilinear_coef(ntri)>val_trilinear) then

                        val_trilinear=trilinear_coef(ntri)

                        v_indices=trind(l,ntri,:)


                        if (v_indices(1)==n1+1) v_indices(1)=n1
                        if (v_indices(1)==0) v_indices(1)=1

                        if (v_indices(2)==n2+1) v_indices(2)=n2
                        if (v_indices(2)==0) v_indices(2)=1

                        if (v_indices(3)==n3+1) v_indices(3)=n3
                        if (v_indices(3)==0) kc=1
                        !write(*,*)ntri,ib_count,'ntri',i,j,k,ic,jc,kc

                        !distanze_CELLE_IB(l,1)=nodo_proiezione(1)-xcd(ic,jc,kc)
                        !distanze_CELLE_IB(l,2)=nodo_proiezione(2)-ycd(ic,jc,kc)
                        !distanze_CELLE_IB(l,3)=nodo_proiezione(3)-zcd(ic,jc,kc)

                    end if
                end do

                ! turn off coefficient if trilinear node is solid
                do ntri=1,4
                    if (trind(l,ntri,2)<=n2 .and. trind(l,ntri,2)>=1) then
                        if (trind(l,ntri,1)<=n1 .and. trind(l,ntri,1)>=1) then
                            if (trind(l,ntri,3)<=n3 .and. trind(l,ntri,3)>=1) then
                                if (tipo(trind(l,ntri,1),trind(l,ntri,2),trind(l,ntri,3))==0) then

                                    trilinear_coef(ntri)=0.0

                                end if
                            end if
                        end if
                    end if
                end do

            else

                write(*,'(A,I0.1,A,3I0.1,A)') 'ID = ',myid,', problem on ib ',i,j,k,' in finding PP'

            end if

            ! check whether all trilinear nodes are invalid
            count_trilinear_coef = 0.0
            do ntri=1,4
                count_trilinear_coef = count_trilinear_coef + trilinear_coef(ntri)
            end do
            if (count_trilinear_coef == 0.0) then
                write(*,'(A,3I0.1)') 'coef nulli',i,j,k
            end if

            !----------------------------------------------------------------------------------
            ! store data in vector

            !v node
!            indici_CELLE_IB(l,:)=ic
!            indici_CELLE_IB(l,5)=jc
!            indici_CELLE_IB(l,6)=kc

            indici_CELLE_IB(l,4)=v_indices(1)
            indici_CELLE_IB(l,5)=v_indices(2)
            indici_CELLE_IB(l,6)=v_indices(3)

            ! PP-IB distance
            !distanze_CELLE_IB(l,1)=distanza_x
            !distanze_CELLE_IB(l,2)=distanza_y
            !distanze_CELLE_IB(l,3)=distanza_z

            !dist_pp_ib(l)=distanza_np_ib

            !----------------------------------------------------------------------------------
            ! ROTAZION MATRIX (to construct tangential and normal velocity)
            call rotation2(ib_position,ip_position,rot_ric,rot_inverse_ric)
            !call rotationAle(sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z,rot_ric,rot_inverse_ric)

            ! forward rotation matrix
            rot(l,:,:)=rot_ric(:,:)

            ! inverse rotation matrix
            rot_inverse(l,:,:)=rot_inverse_ric(:,:)

        end do

        if (myid==0) then

            !write(*,*)' -----------------------------------------'
            !write(*,*)'| SECTION 3: geometry search --> finished |'
            !write(*,*)' -----------------------------------------'

        end if

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
            call MPI_SEND(num_left_rcv,1,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(num_right_snd,1,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     comunicate right the point right has to send left
        if (myid/=nproc-1) then
            call MPI_SEND(num_right_rcv,1,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(num_left_snd,1,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if (kp==0) then
            if (myid==0) then
                call MPI_SEND(num_left_rcv,1,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(num_right_snd,1,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            if (myid==nproc-1) then
                call MPI_SEND(num_right_rcv,1,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
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
            call MPI_SEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     send right the indeces right has to send left
        if (myid/=nproc-1) then
            call MPI_SEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if (kp==0) then

            !       left to 0 is nproc-1
            if (myid==0) then
                call MPI_SEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            ! right to nproc-1 is 0
            if (myid==nproc-1) then
                call MPI_SEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
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


end module ricerca_module

