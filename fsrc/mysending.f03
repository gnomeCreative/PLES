module mysending

    use,intrinsic :: iso_c_binding

    ! contains info for the message passaging


    use mpi

    integer :: MPI_REAL_SD

    integer :: tagls,taglr,tagrs,tagrr
    integer :: leftpe,rightpe
    integer :: leftpem,rightpem
    integer,bind(C) :: myid,nproc
    integer :: iparasta,iparaend
    integer :: ncolperproc

    integer :: deepgl,deepgr
    integer :: deep_mul

contains

    subroutine init_parallel(kgridparasta,kgridparaend,nlevmultimax,bodyforce,insc)


        use scala3, only: n3,kparasta,kparaend,deepl,deepr

        implicit none

        !-----------------------------------------------------------------------
        integer,intent(inout) :: kgridparasta,kgridparaend,nlevmultimax
        integer,intent(in) :: insc
        logical(kind=c_bool),intent(in) :: bodyforce
        !-----------------------------------------------------------------------
        integer :: ierr,ierror
        ! ------------------------------------------------
        ! depending on the setting in imput file MPI_REAL_SD was set toassume the type MPI_REAL4 or MPI_REAL8
        ! ~~~~~~~~ FAKE! only real8 actually works ~~~~~~~~~~
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,8,MPI_REAL_SD,ierr)

        ! definition of :
        ! - the domain decomposition between the processor
        ! - processor recognize left and right processor
        ! - tags for sending
        !
        ! jz must be a multiple of nproc and like nproc*2^n
        ! to allow multigrid
        ! jx must be multiple of nproc to allow transpose procedure
        ncolperproc=int(n3/nproc)

        kparasta=(myid*ncolperproc+1)
        kparaend=((myid+1)*ncolperproc)

        ! myid=0 has kgridparasta=0 because the number of points are odd
        if (myid==0) then
            kgridparasta=kparasta-1
        else
            kgridparasta=kparasta
        end if

        if (myid==nproc-1) then
            kgridparaend=kparaend+1
        else
            kgridparaend=kparaend
        end if

        ! recognize processors dx and sn for each proc
        ! sx of myid=0 is nproc-1
        ! dx of myid=nproc-1 is 0

        leftpe=myid-1
        rightpe=myid+1
        if (leftpe==-1) then
            leftpe=nproc-1
        end if
        if (rightpe==nproc) then
            rightpe=0
        end if

        !tag send to leftpe
        tagls=100+myid
        !tag recv from leftpe
        taglr=110+myid-1
        !tag send to rightpe
        tagrs=110+myid
        !tag recv from rightpe
        tagrr=100+myid+1

        if (myid==0) then
            tagls=tagls+nproc
            taglr=taglr+nproc
        end if

        if (myid==0) then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if (myid==nproc-1) then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else if ((myid/=0).and.(myid/=nproc-1)) then
            leftpem=leftpe
            rightpem=rightpe
        end if

        ! DEPTH OF GHOST LAYERS -----------------------------------
        ! how many plane to allocate less than kparasta --> deepl
        ! ow many plane to allocate less than kparaend --> deepr
        deepl = 1
        deepr = 1

        ! allocation needs one more plane left for quick
        if (insc==1 .and. myid/= nproc-1) then
            deepr = 2
        end if

        if (insc==2 .and. myid/= nproc-1) then
            deepr = 2
        end if


        ! ibm stencil for taylor needs two plane on closer proc
        if (bodyforce) then
            deepl = 2
            deepr = 2
            if (myid==0)       deepl=1
            if (myid==nproc-1) deepr=1
        end if

        ! deep grid right -> deepgr
        ! deep grid left -> deepgl
        ! values in mysending
        deepgr = 2
        deepgl = 2

        ! CHECKING CONDITIONS ON PROCESSORS ARE MET
        if (mod(n3,nproc)/= 0 .or. mod(n3,(2**nlevmultimax)) /= 0) then
            call MPI_ABORT(MPI_COMM_WORLD,ierr,ierror)
            error stop 'ERROR: NUM. PROC. INCORRECT'
        end if


    end subroutine init_parallel

    subroutine border_exchange_centroids(var,depth)

        use scala3, only: n1,n2,kparasta,kparaend,deepl,deepr

        implicit none

        ! ------------------------------------------------
        real,intent(inout) :: var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer,intent(in) :: depth
        ! ------------------------------------------------
        integer :: ierr,status(MPI_STATUS_SIZE)
        integer :: block_size1,block_size2
        ! ------------------------------------------------

        block_size1=(n1+2)*(n2+2)
        block_size2=(n1+2)*(n2+2)*depth

        if (leftpem/=MPI_PROC_NULL) then
            call MPI_SEND(var(0,0,kparasta),block_size2,MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
        end if
        if (rightpem/=MPI_PROC_NULL) then
            call MPI_RECV(var(0,0,kparaend+1),block_size2,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (rightpem/=MPI_PROC_NULL) then
            call MPI_SEND(var(0,0,kparaend),block_size1,MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem/=MPI_PROC_NULL) then
            call MPI_RECV(var(0,0,kparasta-1),block_size1,MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


    end subroutine border_exchange_centroids

    subroutine border_exchange_nodes(var)

        use scala3, only: n1,n2,kparasta,kparaend

        implicit none

        ! ------------------------------------------------
        real,intent(inout) :: var(n1,n2,kparasta-1:kparaend+1)
        ! ------------------------------------------------
        integer :: ierr,status(MPI_STATUS_SIZE)
        integer :: block_size
        ! ------------------------------------------------

        block_size=n1*n2

        if (leftpem/=MPI_PROC_NULL) then
            call MPI_SEND(var(1,1,kparasta),block_size,MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
        end if
        if (rightpem/=MPI_PROC_NULL) then
            call MPI_RECV(var(1,1,kparaend+1),block_size,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (leftpem/=MPI_PROC_NULL) then
            call MPI_RECV(var(1,1,kparasta-1),block_size,MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if
        if (rightpem/=MPI_PROC_NULL) then
            call MPI_SEND(var(1,1,kparaend),block_size,MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


    end subroutine border_exchange_nodes

    subroutine periodic_exchange_z(var)

        use scala3, only: n1,n2,n3,kparasta,kparaend,deepl,deepr

        implicit none

        ! ------------------------------------------------
        real,intent(inout) :: var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        ! ------------------------------------------------
        integer :: ierr,status(MPI_STATUS_SIZE)
        ! ------------------------------------------------

        if (myid==0) then
            call MPI_SEND(var(0,0,1),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,12,MPI_COMM_WORLD,ierr)
            call MPI_RECV(var(0,0,0),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,22,MPI_COMM_WORLD,status,ierr)
        end if
        if (myid==nproc-1) then
            call MPI_RECV(var(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,0,12,MPI_COMM_WORLD,status,ierr)
            call MPI_SEND(var(0,0,n3),(n1+2)*(n2+2),MPI_REAL_SD,0,22,MPI_COMM_WORLD,ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


    end subroutine periodic_exchange_z

    subroutine periodic_exchange_x(var)

        use scala3, only: n1,n2,kparasta,kparaend,deepl,deepr

        implicit none

        ! ------------------------------------------------
        real,intent(inout) :: var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        ! ------------------------------------------------
        integer :: ierr,status(MPI_STATUS_SIZE)
        ! ------------------------------------------------

        var(0,:,:)=var(n1,:,:)
        var(n1+1,:,:)=var(1,:,:)

    end subroutine periodic_exchange_x

    subroutine var_extension(var)

        use scala3, only: n1,n2,n3,kparasta,kparaend,deepl,deepr

        implicit none

        ! ------------------------------------------------
        real,intent(inout) :: var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        ! ------------------------------------------------

        ! sides left and right
        var(0,:,:)=var(1,:,:)
        var(n1+1,:,:)=var(n1,:,:)

        ! sides back and front
        if (myid==0) then
            var(:,:,0)=var(:,:,1)
        else if (myid==nproc-1) then
            var(:,:,n3+1)=var(:,:,n3)
        end if

    end subroutine var_extension

    subroutine var_complete_exchange(var)

        use scala3, only: n1,n2,kparasta,kparaend,deepl,deepr
        use period, only: ip,kp

        implicit none

        ! ------------------------------------------------
        real,intent(inout) :: var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        ! ------------------------------------------------

        ! extend beyond border, useful in case of no periodicity
        call var_extension(var)

        ! border between processors
        call border_exchange_centroids(var,1)

        ! paralle boundaries
        if (ip==0) then
            call periodic_exchange_x(var)
        end if

        if (kp==0) then
            call periodic_exchange_z(var)
        end if

    end subroutine var_complete_exchange



end module mysending
