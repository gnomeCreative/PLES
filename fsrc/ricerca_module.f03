module ricerca_module

    ! release 21-11-2009
    ! Federico Roman
    ! modified by Alessandro Leonardi, starting December 2015

    ! Pre-processing program to use the Immersed Boundary Method as in
    ! Roman,Napoli,Armenio,Milici 2008

    ! GRID: non staggered

    use,intrinsic :: iso_c_binding

    use mysending
    use myarrays_metri3
    use ibm_module
    use scala3

    use mpi

    implicit none

    private

    integer,parameter :: max_solid_index=6
    integer,parameter :: number_of_directions=6

    integer,allocatable,public :: tipo(:,:,:),tipo2(:,:,:)
    ! Identifier for the solid object the IB belongs to
    integer,allocatable :: solid_index(:,:,:,:)
    ! Identifier for the solid object the IB belongs to
    integer,allocatable :: solid_index_size(:,:,:)
    integer,allocatable :: neighbor(:,:,:,:,:)
    logical,allocatable :: is_valid(:,:,:,:)

    public :: set_ibm

contains

    subroutine set_ibm(ktime)

        use particle_module

        integer,intent(in) :: ktime

        integer :: itiposta,itipoend
        integer :: jtiposta,jtipoend
        integer :: ktiposta,ktipoend
        integer :: i,j,k,ii,jj,kk

        ! -------------------------------------------------------------------------
        ! allocation: what needs to be allocated depends on whether bodyforce or particles are active

        ! in any case tipo and tipo2 are needed, allocate it
        if (ktime==0) then
            allocate(tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            allocate(tipo2(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            tipo(:,:,:)=0
            tipo2(:,:,:)=2
        end if

        ! now check if the ibm has to be applied
        if (bodyforce) then
            ! check if ibm comes from particles or from external files
            if (particles) then
                ! do the ricerca cycle
                if (update_ibm .or. ktime==0) then

                    call particle_ricerca(ktime)

                end if

                ! if particle_ricerca will not be executed anymore, we can get rid of arrays

                if (.not.update_ibm .and. ktime==0) then

                    deallocate(neighbor)
                    deallocate(is_valid)
                    deallocate(solid_index)
                    deallocate(solid_index_size)

                end if
            else
                ! this is the old ibm style, with loading of data from external files
                if (ktime==0) then
                    call carico_immb()
                end if
            end if
        else
            ! if no ibm set all the cells to fluid
            tipo(:,:,:)=2
        end if


        ! compute tipo2 if first iteration OR if particles
        if (ktime==0) then
            do k=kparasta,kparaend !1,jz
                do j=1,n2
                    do i=1,n1

                        itiposta=1
                        itipoend=1
                        jtiposta=1
                        jtipoend=1
                        ktiposta=1
                        ktipoend=1

                        if (i==1) itiposta=0
                        if (j==1) jtiposta=0
                        if (k==1) ktiposta=0

                        if (i==n1) itipoend=0
                        if (j==n2) jtipoend=0
                        if (k==n3) ktipoend=0

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

    subroutine carico_immb()

        ! read ibm input file
        use period

        !-----------------------------------------------------------------------
        integer,parameter :: celleib_file=20
        integer,parameter :: cellebloccate_file=22
        integer,parameter :: distanze_file=23
        integer,parameter :: rotazione_file=24
        integer,parameter :: trilinearibm_file=25
        integer,parameter :: trilineari_file=26
        integer,parameter :: trilinearj_file=27
        integer,parameter :: trilineark_file=28

        integer i,j,k,in,jn,kn

        integer status(MPI_STATUS_SIZE),ierror

        integer contatore,num_solide_real
        integer ib_totali,solide_totali

        real :: dist_pp_ib_here,dist_ib_parete_here
        real,dimension(3) :: ip_here
        real,dimension(3,3) :: rot_here,irot_here
        real,dimension(4) :: tricoef_here
        integer,dimension(4,3) :: triindex_here

        integer :: l,np

        allocate(tipo_spedito(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        !-----------------------------------------------------------------------


        if (myid==0) then
            WRITE(*,*)' '
            write(*,*)'*****************************************'
            write(*,*)'     LOAD IBM'
        end if

        ! initialization
        tipo(:,:,:)=2

        ! read number of IB and solid nodes
        open(celleib_file,file='Celle_IB_indici.inp',status='old')
        open(cellebloccate_file,file='Celle_Bloccate_Indici.inp',status='old')
        read(celleib_file,*) numero_celle_IB
        read(cellebloccate_file,*)num_solide
        close(celleib_file)
        close(cellebloccate_file)


        if (myid==0) then
            write(*,*)'matrix ---> OK'
        end if

        ! read input files for ibm

        open(celleib_file,file='Celle_IB_indici.inp',status='old')
        open(distanze_file,file='distanze_interpolazioni.inp',status='old')
        open(rotazione_file,file='rotazione.inp',status='old')

        open(trilinearibm_file,file='trilinear_ibm.inp',status='old')
        open(trilineari_file,file='trilinear_i.inp',status='old')
        open(trilinearj_file,file='trilinear_j.inp',status='old')
        open(trilineark_file,file='trilinear_k.inp',status='old')

        read(celleib_file,*) numero_celle_IB

        num_ib=0

        allocate(indici_CELLE_IB(6,numero_celle_IB))
        allocate(dist_pp_ib(numero_celle_IB))
        allocate(dist_ib_parete(numero_celle_IB))
        allocate(proiezioni(3,numero_celle_IB))
        allocate(ustar(numero_celle_IB))
        allocate(tricoef(4,numero_celle_IB))
        allocate(trind(4,3,numero_celle_IB))
        allocate(rot(3,3,numero_celle_IB))
        allocate(rot_inverse(3,3,numero_celle_IB))

        indici_CELLE_IB(:,:)=0

        do l=1,numero_celle_IB

            read(celleib_file,*) i,j,k,in,jn,kn

            read(distanze_file,'(5E15.8)') dist_pp_ib_here,dist_ib_parete_here, &
                ip_here(1),ip_here(2),ip_here(3)

            read(rotazione_file,'(9E15.8)') &
                rot_here(1,1),rot_here(1,2),rot_here(1,3), &
                rot_here(2,1),rot_here(2,2),rot_here(3,3), &
                rot_here(3,1),rot_here(3,2),rot_here(3,3)


            read(rotazione_file,'(9E15.8)') &
                irot_here(1,1),irot_here(1,2),irot_here(1,3), &
                irot_here(2,1),irot_here(2,2),irot_here(2,3), &
                irot_here(3,1),irot_here(3,2),irot_here(3,3)

            read(trilinearibm_file,'(4E15.8)') &
                tricoef_here(1),tricoef_here(2),tricoef_here(3),tricoef_here(4)
                !tri1,tri2,tri3,tri4

            read(trilineari_file,'(4(I8,1X))') &
                triindex_here(1,1),triindex_here(2,1),triindex_here(3,1),triindex_here(4,1)

            read(trilinearj_file,'(4(I8,1X))') &
                triindex_here(1,2),triindex_here(2,2),triindex_here(3,2),triindex_here(4,2)

            read(trilineark_file,'(4(I8,1X))') &
                triindex_here(1,3),triindex_here(2,3),triindex_here(3,3),triindex_here(4,3)

            if ((k>=kparasta-deepl).and.(k<=kparaend+deepr)) then
                tipo(i,j,k)=1
            end if

            if ((k>=kparasta).and.(k<=kparaend)) then

                tipo(i,j,k)=1

                num_ib=num_ib+1

                indici_CELLE_IB(:,num_ib)=(/i,j,k,in,jn,kn/)

                dist_pp_ib(num_ib)=dist_pp_ib_here
                dist_ib_parete(num_ib)=dist_ib_parete_here

                proiezioni(:,num_ib)=ip_here(:)

                ! rotation matrix (to construct tangential and normal velocity)
                rot(:,:,num_ib)=rot_here(:,:)

                ! inverse rotation matrix
                rot_inverse(:,:,num_ib)=irot_here(:,:)

                ! trilinear coefficent
                tricoef(:,num_ib)=tricoef_here(:)

                ! trilinear index
                trind(:,:,num_ib)=triindex_here(:,:)

            end if
        end do

        close(celleib_file)
        close(distanze_file)
        close(rotazione_file)
        close(trilinearibm_file)
        close(trilineari_file)
        close(trilinearj_file)
        close(trilineark_file)

        write(*,*)myid,', number of ib points: ',num_ib
        !.......................................................................


        ! solid points
        open(cellebloccate_file,file='Celle_Bloccate_Indici.inp',status='old')

        read(cellebloccate_file,*) numero_celle_bloccate

        allocate(indici_celle_bloccate(3,numero_celle_bloccate))

        indici_celle_bloccate(:,:)=0

        num_solide=0
        contatore=0
        do l=1,numero_celle_bloccate

            read(cellebloccate_file,*)i,j,k

            if ((k>=kparasta-deepl).and.(k<=kparaend+deepr)) then

                tipo(i,j,k)=0

                num_solide=num_solide+1

                indici_celle_bloccate(1,num_solide)=i
                indici_celle_bloccate(2,num_solide)=j
                indici_celle_bloccate(3,num_solide)=k

                if (k>kparaend .or. k<kparasta) then
                    contatore=contatore+1
                end if
            end if
        end do

        close(cellebloccate_file)

        ! velocity is always zero if immersed body comes from file
        allocate(surfvel_ib(3,num_ib))
        allocate(solidvel_ib(3,num_solide))

        surfvel_ib(:,:)=0.0
        solidvel_ib(:,:)=0.0


        if (myid==0) then
            write(*,*)'CHECK: call carico_immb----> OK'
            write(*,*)'values for num_ib and num_solide: ',num_ib,num_solide
        end if

        write(*,*)'border solid',contatore

        !without border cells
        num_solide_real=num_solide-contatore

        write(*,*)myid,'number of solid cells: ',num_solide

        !-----------------------------------------------------------------------
        !     check number of ibm for each proc

        call MPI_REDUCE(num_ib,ib_totali,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
        call MPI_REDUCE(num_solide_real,solide_totali,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)

        if (myid==0) then
            write(*,*)'--------------------------------------------'
            if (ib_totali==numero_celle_IB) then
                write(*,*)'check IB for each proc --> OK'
            else
                write(*,*)'check IB for each proc --> NO'
            end if

            if (solide_totali==numero_celle_bloccate) then
                write(*,*)'check solid cells for each proc --> OK'
            else
                write(*,*)'check solid cells for each proc --> NO'
            end if

            write(*,*)'--------------------------------------------'
            write(*,*)'read input file for Immersed Boundaries'
            write(*,*)'total number of IB',numero_celle_IB
            write(*,*)'total number of solid',numero_celle_bloccate
            write(*,*)'--------------------------------------------'
        end if

        write(*,*)'number IB proc --->',num_ib
        write(*,*)'number solid cells proc --->',num_solide
        write(*,*)'--------------------------------------------'

        if (myid==0) then
            write(*,*)'       LOAD IBM finished'
            write(*,*)'*****************************************'
            write (*,*)' '
        end if

        !-----------------------------------------------------------------------
        !     PREPARE COMUNICATION BETWEEN PROCS FOR IB
        allocate( stencil_left_snd(2*n1*n2,3))
        allocate(stencil_right_snd(2*n1*n2,3))
        allocate( stencil_left_rcv(2*n1*n2,3))
        allocate(stencil_right_rcv(2*n1*n2,3))

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
                i=trind(np,1,l)
                j=trind(np,2,l)
                k=trind(np,3,l)


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
            call MPI_SEND(stencil_left_rcv(1,1),3*2*n1*n2,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(stencil_right_snd(1,1),3*2*n1*n2,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     send right the indeces right has to send left
        if (myid/=nproc-1) then
            call MPI_SEND(stencil_right_rcv(1,1),3*2*n1*n2,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(stencil_left_snd(1,1),3*2*n1*n2,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if (kp==0) then

            !       left to 0 is nproc-1
            if (myid==0) then
                call MPI_SEND(stencil_left_rcv(1,1),3*2*n1*n2,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(stencil_right_snd(1,1),3*2*n1*n2,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            ! right to nproc-1 is 0
            if (myid==nproc-1) then
                call MPI_SEND(stencil_right_rcv(1,1),3*2*n1*n2,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
            end if
            if (myid==0) then
                call MPI_RECV(stencil_left_snd(1,1),3*2*n1*n2,MPI_INTEGER,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            end if

            !       now myid=0 and myid=nproc-1 know the index they will recive,
            !       but they need to save
            !       on the border, so I change the k

            if (myid == nproc-1) stencil_right_snd(:,3)=n3
            if (myid == 0) stencil_left_snd(:,3)=1

        end if

        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------

        return

    end subroutine carico_immb

    subroutine particle_ricerca(ktime)

        ! ----------------------------------------------------
        integer,intent(in) :: ktime
        ! ----------------------------------------------------
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

        use :: particle_module, only: num_part_tot,sphereIndex,is_inside,is_inprocessor

        integer,intent(in) :: ktime

        ! ----------------------------------------------------
        ! For solid/fluid/IBM research

        real :: isoparameter
        integer :: solid_count,ib_count,doublers_count
        integer :: counterTot,reducedCounterTot
        integer :: num_solide_real,numero_celle_bloccate_real
        integer :: num_doublers,num_doublers_tot
        integer :: counterFluid,contatore
        integer :: reducedCounterFluid
        integer :: i,j,k,p,m,ni,nj,nk
        integer :: indexHere
        integer :: ierr

        ! SECTION 2: PHASES

        ! determine if a cell is solid, fluid or ib starting from zero or from a restart file

        !write(*,*) 'Proc ',myid,'scanning ',kparasta-deepl,' to ',kparaend+deepr

        ! compute neighbor matrix, only at first time step
        if (ktime==0) then
            if (myid==0) then
                write(*,*) 'Building up neighbor matrix'
            end if

            allocate(neighbor(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,number_of_directions,3))
            allocate(is_valid(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,number_of_directions))
            allocate(solid_index(max_solid_index,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            allocate(solid_index_size(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))

            ! compute neighbors
            is_valid=.false.
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        do m=1,number_of_directions
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
                        if (i+1<=n1+1) is_valid(i,j,k,1)=.true.
                        if (i-1>=0) is_valid(i,j,k,2)=.true.
                        if (j+1<=n2+1) is_valid(i,j,k,3)=.true.
                        if (j-1>=0) is_valid(i,j,k,4)=.true.
                        if (k+1<=n3+1) is_valid(i,j,k,5)=.true.
                        if (k-1>=0) is_valid(i,j,k,6)=.true.
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
        solid_index(:,:,:,:)=0 ! consider -1
        solid_index_size(:,:,:)=0 ! consider -1

        ! PART 1 - LOOKING FOR SOLID NODES -----------------------------------------------------------------

        num_solide=0
        contatore=0
        if (myid==0) then
            write(*,*) 'Looking for solid nodes'
        end if
        ! nodes inside a particle are treated as solid
        do p=1,num_part_tot

            indexHere=sphereIndex(p)
            if (is_inprocessor(indexHere)) then
            do k=kparasta-deepl,kparaend+deepr
                do j=0,n2+1
                    do i=0,n1+1
                        if (tipo(i,j,k)/=0) then

                            if (is_inside(centroid(:,i,j,k),indexHere)) then

                                tipo(i,j,k)=0
                                solid_index(1,i,j,k)=indexHere

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
            end if
        end do

        num_solide_real=num_solide-contatore !without border cells

        ! now we can allocate variables that are defined in this subroutine
        if (allocated(indici_celle_bloccate)) then
            deallocate(indici_celle_bloccate)
        end if

        allocate(indici_celle_bloccate(3,num_solide))

        solid_count=0
        do k=kparasta-deepl,kparaend+deepr
            do j=1,n2
                do i=1,n1
                    if (tipo(i,j,k)==0) then !1

                        solid_count=solid_count+1

                        indici_celle_bloccate(1,solid_count)=i
                        indici_celle_bloccate(2,solid_count)=j
                        indici_celle_bloccate(3,solid_count)=k

                    end if
                end do
            end do
        end do

        ! check if everything went ok so far
        if (solid_count/=num_solide) then
            write(*,*) 'PROBLEM Proc: ',myid,' solid_count=',solid_count,', num_solide=',num_solide
            call exit(0)
        end if

        ! PART 2 - LOOKING FOR IB NODES -----------------------------------------------------------------

        ! set boundaries
        if (myid==0) then
            write(*,*) 'Looking for IB nodes'
        end if

        num_ib=0
        num_doublers=0
        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    ! check if it's fluid
                    if (tipo(i,j,k)==2) then !0
                        ! check if neighbors are solid. in this case, turn the node into IB
                        do m=1,number_of_directions
                            if (is_valid(i,j,k,m)) then
                                ni=neighbor(i,j,k,m,1)
                                nj=neighbor(i,j,k,m,2)
                                nk=neighbor(i,j,k,m,3)
                                if (tipo(ni,nj,nk)==0) then
                                    if (k>=kparasta .and. k<=kparaend .and. &
                                            i>=1 .and. i<=n1 .and. j>=1 .and. j<=n2) then
                                        num_ib=num_ib+1
                                    end if
                                    ! update solid neighbors
                                    solid_index_size(i,j,k)=solid_index_size(i,j,k)+1
                                    solid_index(solid_index_size(i,j,k),i,j,k)=solid_index(1,ni,nj,nk)
                                    ! turn it into IB only if that has not happened before
                                    ! (e.g. with another neighbor)
                                    if (tipo(i,j,k)/=1) then
                                        tipo(i,j,k)=1
                                    else
                                        if (k>=kparasta .and. k<=kparaend .and. &
                                            i>=1 .and. i<=n1 .and. j>=1 .and. j<=n2) then
                                            num_doublers=num_doublers+1
                                        end if
                                    end if
                                end if
                            end if
                        end do
                    end if
                end do
            end do
        end do

        ! now we can allocate variables that are defined in this subroutine
        if (allocated(indici_CELLE_IB)) then
            deallocate(indici_CELLE_IB)
            deallocate(position_ib)
            deallocate(indexsize_ib,index_ib)
        end if

        allocate(indici_CELLE_IB(6,num_ib))
        allocate(position_ib(3,num_ib))
        allocate(index_ib(num_ib),indexsize_ib(num_ib))

        ib_count=0
        doublers_count=0
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    if (tipo(i,j,k)==1) then

                        ! add one ib_point for each solid neighbor
                        do m=1,solid_index_size(i,j,k)

                            ib_count=ib_count+1

                            indici_CELLE_IB(1,ib_count)=i
                            indici_CELLE_IB(2,ib_count)=j
                            indici_CELLE_IB(3,ib_count)=k

                            position_ib(:,ib_count)=centroid(:,i,j,k)

                            index_ib(ib_count)=solid_index(m,i,j,k)

                            indexsize_ib(ib_count)=solid_index_size(i,j,k)

                        end do

                        if (solid_index_size(i,j,k)>1) then

                            doublers_count=doublers_count+solid_index_size(i,j,k)-1

                        end if
                    end if
                end do
            end do
        end do

        ! check if everything went ok so far
        if (ib_count/=num_ib) then
            write(*,*) 'PROBLEM Proc: ',myid,' ib_count=',ib_count,', num_ib=',num_ib
            call exit(0)
        end if
        if (doublers_count/=num_doublers) then
            write(*,*) 'PROBLEM Proc: ',myid,' doublers_count=',doublers_count,', num_doublers=',num_doublers
            call exit(0)
        end if
        ! PART 3 - CHECK IF EVERYTHING IS ALRIGHT -----------------------------------------------------------------

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
        call MPI_ALLREDUCE(num_doublers,num_doublers_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(num_solide_real,numero_celle_bloccate_real,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(counterFluid,reducedCounterFluid,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(num_ib,numero_celle_IB,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! checking for mistakes
        if (myid==0) then
            write(*,*) 'Total cells = ',reducedCounterTot
            write(*,*) '   Total solid = ',numero_celle_bloccate_real
            write(*,*) '   Total IB = ',numero_celle_IB,' of which doublers = ',num_doublers_tot
            write(*,*) '   Total fluid = ',reducedCounterFluid


            if (reducedCounterTot/=numero_celle_bloccate_real+numero_celle_IB+reducedCounterFluid-num_doublers_tot) then
                write(*,*) 'Error in assigning type'
                call exit(0)
            end if
        end if

        return

    end subroutine phase_search

    subroutine interface_search()

        use particle_module

        ! vectors for solid node velocity
        real,dimension(3) :: solidpoint_here
        real,dimension(3) :: solidvel_here
        ! iterators
        integer :: i,j,k,l
        integer :: index_here,solidindex_size_here

        ! ----------------------------------------------------------------------------

        ! allocate variables that will be defined here
        if (allocated(proiezioni)) then
            deallocate(dist_ib_parete)
            deallocate(surfnormal_ib)
            deallocate(proiezioni)
            deallocate(surfvel_ib)
            deallocate(solidvel_ib)
        end if

        allocate(dist_ib_parete(num_ib))
        allocate(surfnormal_ib(3,num_ib))
        allocate(proiezioni(3,num_ib))
        allocate(surfvel_ib(3,num_ib))
        allocate(solidvel_ib(3,num_solide))

        dist_ib_parete(:)=0.0
        surfnormal_ib(:,:)=0.0
        proiezioni(:,:)=0.0
        surfvel_ib(:,:)=0.0
        solidvel_ib(:,:)=0.0

        !---------------------------------------------------------------------------------------------------------------------------------------------
        ! SECTION 2: search the projection point IP on the body, through the normal from the IB to the body
        do l=1,num_ib

            !solidindex_size_here=indexsize_ib(l)

            !if (solidindex_size_here<=0) then
            !   write(*,*) 'Problem: negative index size'
            !    call exit(0)
            !end if

            ! particle index
            index_here=index_ib(l)

            ! normal
            surfnormal_ib(:,l)=surface_normal(position_ib(:,l),index_here)
            ! distance from surface
            dist_ib_parete(l)=norm2(surface_distance(position_ib(:,l),index_here))
            ! compute surface point
            proiezioni(:,l)=position_ib(:,l)-surfnormal_ib(:,l)*dist_ib_parete(l)
            ! surface velocity
            surfvel_ib(:,l)=point_velocity(proiezioni(:,l),index_here)
        end do

        ! loop also on solid point to determine the particle velocity at that position
        do l=1,num_solide

            i=indici_celle_bloccate(1,l)
            j=indici_celle_bloccate(2,l)
            k=indici_celle_bloccate(3,l)

            solidpoint_here(:)=centroid(:,i,j,k)

            ! for solid points the particle index is always at position 1
            index_here=solid_index(1,i,j,k)

            ! velocty of the particle at that position
            solidvel_here(:)=point_velocity(solidpoint_here,index_here)

            solidvel_ib(:,l)=solidvel_here(:)

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

        real :: trilinear_coef(4)

        !index
        integer :: i,j,k,l
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
        logical :: coef_nulli

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
        allocate(tricoef(4,num_ib))
        allocate(trind(4,3,num_ib))
        allocate(rot(3,3,num_ib))
        allocate(rot_inverse(3,3,num_ib))

        tricoef(:,:)=-999999.0
        trind(:,:,:)=-999999

        !-----------------------------------------------------------------------
        ! PROJECTION POINT PP SEARCHIG PROCEDURE
        !-----------------------------------------------------------------------

        ! initialization of node coordination for trilinear
        call triangles_mesh()

        coef_nulli=.false.

        do l=1,num_ib

            ! position of current IB point
            ib_position(:)=position_ib(:,l)

            i=indici_CELLE_IB(1,l)
            j=indici_CELLE_IB(2,l)
            k=indici_CELLE_IB(3,l)

            !solidindex_size_here=indexsize_ib(l)

            ! position of current IP point
            ip_position(:)=proiezioni(:,l)


            !----------------------------------------------------------------------------------
            !search closest fluid node to PP

            ! initialize current box to IB index, and coefficients to zero
            v_indices(:)=(/i,j,k/)
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

                tricoef(:,l)=trilinear_coef(:)

                ! v node
                val_trilinear=0.0
                do ntri=1,4
                    if (trilinear_coef(ntri)>val_trilinear) then

                        val_trilinear=trilinear_coef(ntri)

                        v_indices=trind(ntri,:,l)


                        if (v_indices(1)==n1+1) v_indices(1)=n1
                        if (v_indices(1)==0) v_indices(1)=1

                        if (v_indices(2)==n2+1) v_indices(2)=n2
                        if (v_indices(2)==0) v_indices(2)=1

                        if (v_indices(3)==n3+1) v_indices(3)=n3
                        if (v_indices(3)==0) kc=1


                    end if
                end do

                ! turn off coefficient if trilinear node is solid
                do ntri=1,4
                    if (trind(ntri,2,l)<=n2 .and. trind(ntri,2,l)>=1) then
                        if (trind(ntri,1,l)<=n1 .and. trind(ntri,1,l)>=1) then
                            if (trind(ntri,3,l)<=n3 .and. trind(ntri,3,l)>=1) then
                                if (tipo(trind(ntri,1,l),trind(ntri,2,l),trind(ntri,3,l))==0) then

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
            if (count_trilinear_coef==0.0 .and. .not.coef_nulli) then
                coef_nulli=.true.
                write(*,'(A,3I0.1)') 'coef nulli',i,j,k
            end if

            !----------------------------------------------------------------------------------
            ! store data in vector

            indici_CELLE_IB(4,l)=v_indices(1)
            indici_CELLE_IB(5,l)=v_indices(2)
            indici_CELLE_IB(6,l)=v_indices(3)

            !----------------------------------------------------------------------------------
            ! ROTAZION MATRIX (to construct tangential and normal velocity)
            call rotation2(ib_position,ip_position,rot_ric,rot_inverse_ric)
            !call rotationAle(sx,sy,sz,nodo_vicino_x,nodo_vicino_y,nodo_vicino_z,rot_ric,rot_inverse_ric)

            ! forward rotation matrix
            rot(:,:,l)=rot_ric(:,:)

            ! inverse rotation matrix
            rot_inverse(:,:,l)=rot_inverse_ric(:,:)

        end do

        if (coef_nulli) then
                write(*,'(A)') 'WARNING: coef nulli'
            end if

        if (myid==0) then

            !write(*,*)' -----------------------------------------'
            !write(*,*)'| SECTION 3: geometry search --> finished |'
            !write(*,*)' -----------------------------------------'

        end if

        return

    end subroutine projection_search

    subroutine prepare_communication()

        use period, only: kp

        integer :: i,j,k,l,np
        integer :: ierror
        integer :: status(MPI_STATUS_SIZE)

        !     PREPARE COMUNICATION BETWEEN PROCS FOR IB

        if (allocated(stencil_left_snd)) then
            deallocate(stencil_left_snd,stencil_right_snd)
            deallocate(stencil_left_rcv,stencil_right_rcv)
            deallocate(tipo_spedito)
        end if

        allocate( stencil_left_snd(2*n1*n2,3))
        allocate(stencil_right_snd(2*n1*n2,3))
        allocate( stencil_left_rcv(2*n1*n2,3))
        allocate(stencil_right_rcv(2*n1*n2,3))

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
                i=trind(np,1,l)
                j=trind(np,2,l)
                k=trind(np,3,l)

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
            call MPI_SEND(stencil_left_rcv(1,1),3*2*n1*n2,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(stencil_right_snd(1,1),3*2*n1*n2,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     send right the indeces right has to send left
        if (myid/=nproc-1) then
            call MPI_SEND(stencil_right_rcv(1,1),3*2*n1*n2,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(stencil_left_snd(1,1),3*2*n1*n2,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if (kp==0) then

            !       left to 0 is nproc-1
            if (myid==0) then
                call MPI_SEND(stencil_left_rcv(1,1),3*2*n1*n2,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if (myid==nproc-1) then
                call MPI_RECV(stencil_right_snd(1,1),3*2*n1*n2,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            ! right to nproc-1 is 0
            if (myid==nproc-1) then
                call MPI_SEND(stencil_right_rcv(1,1),3*2*n1*n2,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
            end if
            if (myid==0) then
                call MPI_RECV(stencil_left_snd(1,1),3*2*n1*n2,MPI_INTEGER,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            end if

            !       now myid=0 and myid=nproc-1 know the index they will recive,
            !       but they need to save
            !       on the border, so I change the k

            if (myid == nproc-1)stencil_right_snd(:,3)=n3
            if (myid == 0)      stencil_left_snd(:,3)=1

        end if

        return

    end subroutine prepare_communication


end module ricerca_module

