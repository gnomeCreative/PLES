!***********************************************************************
module output_module
    !***********************************************************************
    !     contains arrays and settings for printing procedure

    use iso_c_binding

    use mysending, only: kparaend,kparasta,myid,nproc,deepl,deepr
    !
    use scala3, only: jx,jy,jz,nscal,re,n1,n2,n3,dt
            !
    use print

    implicit none

    ! flags
    integer,bind(C) :: print_iter_or_time,i_printfile,iformat_newres
    logical,bind(C) :: i_paraview,i_newres,i_medietempo
    logical,bind(C) :: i_cumulative,iformat_grid
    character(len=500) :: result_folder
    !character,bind(C),dimension(*) :: result_folder_name,simulation_name


    logical :: iscrivo

    ! information output files
    integer,parameter :: info_run_file=11
    integer,parameter :: new_res_file=14
    integer,parameter :: medietempo_file=14
    integer,parameter :: paraview_file=17
    integer,parameter :: paraview_piece_file=18
    integer,parameter :: turbo_file=29
    integer,parameter :: subdis_file=30
    integer,parameter :: subscl_file=31
    integer,parameter :: subrho_file=33
    integer,parameter :: bulkvel_file=34
    integer,parameter :: encheck_file=13
    integer,parameter :: drag_file=54
    integer,parameter :: maxvel_file=53
    integer,parameter :: sonde_file_base=2000

    ! formats
    character(len=12) :: string_newres_format
    character(len=12) :: string_grid_format
    character(len=20) :: fmt_newres

    ! file names
    character(len=500) :: filename_newres,filename_medietempo,filename_paraview
    character(len=500) :: filename_movie,filename_paraview_piece
    !
    character(len=500) :: filename_inforun,filename_turbo,filename_subdis
    character(len=500) :: filename_subscl,filename_subrho,filename_drag
    character(len=500) :: filename_bulkvel,filename_maxvel,filename_encheck
    character(len=500) :: flename_sondefilebase

    character*120,allocatable   :: paraview_piecename(:)

    ! folders to store results (see output_init)
    !character(len=23) :: result_folder ! length = length(result_folder_name)+length('results/')
    character(len=34) :: new_res_folder
    character(len=36) :: paraview_folder
    character(len=38) :: medietempo_folder

    character(len=4) :: char_myid

    real,bind(C) :: i_time
    real ti_start
      
    ! variables for printing
    real, allocatable :: umed(:,:,:,:) !(1:n1,1:n2,kparasta:kparaend,4)
    real, allocatable :: uutau(:,:,:,:) !(1:n1,1:n2,kparasta:kparaend,7)
    real :: deltatempo

    real,allocatable :: print_u(:,:,:),print_v(:,:,:),print_w(:,:,:)
      
    real,allocatable :: print_uc(:,:,:),print_vc(:,:,:),print_wc(:,:,:)
      
    real,allocatable :: print_fi(:,:,:)
    real,allocatable :: print_rhov(:,:,:,:)
    real,allocatable :: print_annit(:,:,:)
    real,allocatable :: print_annitv(:,:,:)
    real,allocatable :: print_akapt(:,:,:)
    real,allocatable :: print_akaptv(:,:,:)
       
    real,allocatable :: rhocol(:),rhotot(:)
    real,allocatable :: ficol(:),fitot(:)
    real,allocatable :: ucol(:),utot(:)
    real,allocatable :: vcol(:),vtot(:)
    real,allocatable :: wcol(:),wtot(:)
    real,allocatable :: uccol(:),uctot(:)
    real,allocatable :: vccol(:),vctot(:)
    real,allocatable :: wccol(:),wctot(:)
    real,allocatable :: annitcol(:),annittot(:)
    real,allocatable :: akaptcol(:),akapttot(:)
    real,allocatable :: upncol(:),upntot(:)
    real,allocatable :: vpncol(:),vpntot(:)
    real,allocatable :: wpncol(:),wpntot(:)
    real,allocatable :: updcol(:),updtot(:)
    real,allocatable :: vpdcol(:),vpdtot(:)
    real,allocatable :: wpdcol(:),wpdtot(:)
    real,allocatable :: uptcol(:),upttot(:)
    real,allocatable :: vptcol(:),vpttot(:)
    real,allocatable :: wptcol(:),wpttot(:)
    real,allocatable :: upqcol(:),upqtot(:)
    real,allocatable :: vpqcol(:),vpqtot(:)
    real,allocatable :: wpqcol(:),wpqtot(:)

contains

    !**********************************************************************
    subroutine create_output_folder()

        use mysettings, only: lagr,nsonde

        use mpi

        implicit none

        character(len=12),parameter :: paraview_folder_name = 'paraviewData'
        character(len=10),parameter :: new_res_folder_name = 'newResData'
        character(len=14),parameter :: medietempo_folder_name = 'medietempoData'
        !character(len=7),parameter  :: global_result_folder_name = 'results'
        !character(len=8)            :: date
        !character(len=6)            :: time
        !character(len=15)           :: result_folder_name ! length = date+time+1

        ! for sonde
        character*3  identificosonda
        character*20 filesonda

        integer isonde
        integer ierr
        integer ichar

        !-----------------------------------------------------------------------
        ! check date and time for the result folder creation
        !call date_and_time(date,time)

        if (myid == 0) then

            write(*,*) 'Initialize output step'

            if (lagr == 0) then
                write(*,*) 'Warning: in output files no subrid section has been implemented'
            end if

            ! global result folder
            !write(*,*) 'Check if global result folder exist'
            !call system('mkdir '//trim(global_result_folder_name))

        end if

        ! result folder
        !result_folder_name=date//'_'//time
        !result_folder=global_result_folder_name//'/'//result_folder_name
        !if (myid == 0) then
        !    write(*,*) 'Create result folder (format=date_time)'
        !    call system('mkdir '//trim(result_folder))
        !end if

        ! paraview folder
        if (i_paraview) then
            paraview_folder=trim(result_folder)//'/'//paraview_folder_name
            if (myid == 0) then
                call system('mkdir '//paraview_folder)
                write(*,*) 'Created folder ',paraview_folder
            end if

        end if

        ! newRes folder
        if (i_newres) then
            new_res_folder=trim(result_folder)//'/'//new_res_folder_name
            if (myid == 0) then
                call system('mkdir '//new_res_folder)
                write(*,*) 'Created folder ',new_res_folder
            end if
        end if

        ! medietempo folder
        if (i_medietempo) then
            medietempo_folder=trim(result_folder)//'/'//medietempo_folder_name
            if (myid == 0) then
                call system('mkdir '//medietempo_folder)
                write(*,*) 'Created folder ',medietempo_folder
            end if

        end if

        ! Files for information output -------------------------------------
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        filename_inforun=trim(result_folder)//'/'//'info_run.txt'
        open(info_run_file,file=filename_inforun,action='write',status='replace')

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     open output files

        if (myid==0) then
            if (lagr==0) then
                filename_turbo=trim(result_folder)//'/'//'turbo.out'
                filename_subdis=trim(result_folder)//'/'//'subdis.out'
                filename_subscl=trim(result_folder)//'/'//'subscl.out'
                filename_subrho=trim(result_folder)//'/'//'subrho.out'
                open(turbo_file,file=filename_turbo,action='write',status='replace')
                open(subdis_file,file=filename_subdis,action='write',status='replace')
                open(subscl_file,file=filename_subscl,action='write',status='replace')
                open(subrho_file,file=filename_subrho,action='write',status='replace')
            end if

            filename_bulkvel=trim(result_folder)//'/'//'bulk_velocity.dat'
            open(bulkvel_file,file=filename_bulkvel,action='write',status='replace')
            filename_maxvel=trim(result_folder)//'/'//'max_velocity.dat'
            open(maxvel_file,file=filename_maxvel,action='write',status='replace')
            filename_encheck=trim(result_folder)//'/'//'encheck.dat'
            open(encheck_file,file=filename_encheck,action='write',status='replace')
            filename_drag=trim(result_folder)//'/'//'drag.dat'
            open(drag_file,file=filename_drag,action='write',status='replace')

            !  sonde files
            do isonde=1,nsonde
                write(identificosonda,'(i3.3)')isonde
                filesonda = trim(result_folder)//'nsonda'//identificosonda//'.dat'
                open(sonde_file_base+isonde,file=filesonda,action='write',status='unknown')
            end do

        end if


    end subroutine create_output_folder

    !***********************************************************************
    subroutine output_init()

        use mysettings, only: lagr,nsonde

        use mpi

        implicit none

        integer ierr
        integer ichar

        !-----------------------------------------------------------------------
        ! check date and time for the result folder creation

        !-----------------------------------------------------------------------
        ! Setting for input/output formats
        string_newres_format="10e18.10"
        string_grid_format="10e18.10"


        allocate(paraview_piecename(0:nproc-1))

        ! medietempo allocation
        if (i_medietempo) then

            ! allocation
            allocate( umed(1:n1,1:n2,kparasta:kparaend,4+nscal))
            allocate(uutau(1:n1,1:n2,kparasta:kparaend,6+nscal+nscal))
            ! initialization
            deltatempo=0.
            umed  = 0.
            uutau = 0.

        end if

        ! string form of processor identifier
        write (char_myid,'(i4)') myid
        do ichar = 1,len(char_myid)
            if (char_myid(ichar:ichar)==' ') char_myid(ichar:ichar)='0'
        end do
        !rid_myid   =len_trim(char_myid)

        ! if I write with format recognize it
        ! set the format
        if (iformat_newres == 1) then
            string_newres_format = adjustl(string_newres_format)  ! with format

            fmt_newres = '('//trim(string_newres_format)//')'
            if (myid == 0) write(*,*) 'New_res file format: ',fmt_newres
        end if

    end subroutine output_init

    !***********************************************************************
    subroutine output_step(count_print_time,dbbx,tipo)

        !-----------------------------------------------------------------------
        !++++++++++  WRITE OUTPUT  +++++++++++++++++++++++++++++++++++++++++++++

        use parti, only: alx, alz, ti
        use mysettings, only: bbx, niter
        use convex, only: bulk
        use myarrays_ibm, only: bodyupdate,bodyforce
        !
        use mpi

        implicit none

        ! arguments
        integer,intent(inout) :: count_print_time
        real,intent(in) :: dbbx
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        ! local variables
        integer :: next_print
        integer :: ierr
        real :: time_evaluation
        logical :: writeprocedure
        !-----------------------------------------------------------------------
        ! print restart file in the form new_res

        ! check if it's time to print an output file
        ! initial guess: no
        writeprocedure = .false.

        ! check condition on reached number of iterations
        if (print_iter_or_time ==0 .and. ktime==i_print*(ktime/i_print)) then
            writeprocedure = .true.

        ! check condition on reached time
        else if (print_iter_or_time == 1 .and. &
            ( ((ti-ti_start)-count_print_time*i_time)>= i_time .or. &
            abs(((ti-ti_start)-count_print_time*i_time)-i_time )<10e-7 ) ) then

            writeprocedure = .true.

            count_print_time = count_print_time + 1
        end if


        ! if it's time to write output, start prodecure
        if (writeprocedure) then

            !-----------------------------------------------------------------------
            if (myid == 0) write(*,*) 'OUTPUT STEP BEGINS ---------------------------------'

                ! set the folder for the output and the filename
            call write_preparation(ti)

            ! write output new_res_form without mpi procedure
            if ((i_newres .and. i_printfile == 0) .or. (i_newres .and. i_printfile == 2)) then
                ! collect data to print
                call prepare_printdata
                ! and produce new_res
                call write_newres(ti,bbx,dbbx,tipo)
            end if

            ! print paraview
            if (i_paraview) then
                call write_paraview(tipo)
            end if

            ! print medietempo
            if (i_medietempo) then
                call write_medietempo()
            end if


            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

                !
            !-----------------------------------------------------------------------
            !
            if (myid == 0) write(*,*) 'OUTPUT STEP ENDS -----------------------------------'

            !     check if I can write another new_res before ending the number of iterations
            if (print_iter_or_time == 1) then

                ! time_evaluation = (niter-ktime)*dt
                ! time_evaluation =  niter*dt - ti_start

                time_evaluation = (ti-ti_start) + (niter-ktime)*dt
                next_print = time_evaluation/real(i_time)
                if (myid == 0) write(*,*) 'time evaluation ',time_evaluation,next_print,count_print_time

                if (next_print <= count_print_time .and. ktime/=niter) then

                    !         if (time_evaluation < i_time) then

                    if (myid ==0) then
                        write(*,*)'remaining iteration not sufficent'
                        write(*,*)'to write another new_res'
                        write(*,*)'run is stopped'
                    end if

                    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

                    stop
                end if
            end if


        end if !write procedure

        ! print data from sonde
        call write_tracers(ti)

        ! print kinetic energy
        call print_energy(tipo)

        ! print bulk velocity
        call print_velocity(tipo)

        ! print bulk velocity
        if (bodyforce==1 .and. bodyupdate) then
            call print_drag()
        end if

    end subroutine

    !***********************************************************************
    subroutine write_preparation(ti)
        ! this subroutine prepares the file names and the folder

        use mpi

        implicit none

        ! file names
        character(len=7),parameter :: filename_header= 'new_res'
        character(len=9),parameter :: filemedie_header='meantime'
        character(len=4),parameter  :: paraview_header ='para'

        ! arguments
        real,intent(in) :: ti

        !-----------------------------------------------------------------------
        !  local variables declaration
        integer ichar,iproc
        integer ierr

        character(len=6)  :: char_ikiter
        character(len=13) :: char_iktime
        character(len=4)  :: char_piece

        !-------------------------------------------------------------------------

        ! set the folder for the output
        ! print for iterations
        if (print_iter_or_time == 0) then

            write (char_ikiter,'(i6)') ktime
            do ichar = 1,len(char_ikiter)
                if (char_ikiter(ichar:ichar)==' ') char_ikiter(ichar:ichar)='0'
            end do

        ! print for time
        else if (print_iter_or_time == 1) then

            write(char_iktime,'(1e13.7)') ti ! INT(ti-ti_start)
            write(*,*) char_iktime

            do ichar = 1,len(char_iktime)
                if (char_iktime(ichar:ichar)==' ') char_iktime(ichar:ichar)='0'
            end do


        end if

        ! flag for telling every processor if it must write
        iscrivo = .false.

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !-----------------------------------------------------------------------
        ! setting up names

        ! 1: set the base filename for new_res, paraview and medetempo
        filename_newres=trim(new_res_folder)//'/'//filename_header
        filename_medietempo=trim(medietempo_folder)//'/'//filemedie_header
        filename_paraview=trim(paraview_folder)//'/'//paraview_header
        filename_paraview_piece=trim(paraview_folder)//'/'//paraview_header
        do iproc=0,nproc-1
            paraview_piecename(iproc)=paraview_header
        end do

        ! 2: add iteration or time ID
        if (print_iter_or_time == 0) then
            ! add interation number
            filename_newres=trim(filename_newres)//'_iter'//char_ikiter
            filename_medietempo=trim(filename_medietempo)//'_iter'//char_ikiter
            filename_paraview=trim(filename_paraview)//'_iter'//char_ikiter
            filename_paraview_piece=trim(filename_paraview_piece)//'_iter'//char_ikiter
            do iproc=0,nproc-1
                paraview_piecename(iproc)=trim(paraview_piecename(iproc))//'_iter'//char_ikiter
            end do
        else if (print_iter_or_time == 1) then
            ! add time
            filename_newres=trim(filename_newres)//'_time'//char_iktime//'sec'
            filename_medietempo=trim(filename_medietempo)//'_time'//char_iktime//'sec'
            filename_paraview=trim(filename_paraview)//'_time'//char_iktime//'sec'
            filename_paraview_piece=trim(filename_paraview_piece)//'_time'//char_iktime//'sec'
            do iproc=0,nproc-1
                paraview_piecename(iproc)=trim(paraview_piecename(iproc))//'_time'//char_iktime//'sec'
            end do
        end if

        ! 3: add processor ID
        ! new_res can be parallel or not
        if (i_printfile == 0 .or. i_printfile == 1) then
            filename_newres=trim(filename_newres)//'.dat'
        else if (i_printfile == 2) then
            ! add processor id
            filename_newres=trim(filename_newres)//'_proc'//char_myid//'.dat'
        end if
        ! medietempo is always parallel
        filename_medietempo=trim(filename_medietempo)//'_proc'//char_myid//'.dat'
        ! paraview is also always parallel (but with header)
        filename_paraview_piece=trim(filename_paraview_piece)//'_proc'//char_myid//'.vti'
        filename_paraview=trim(filename_paraview)//'.pvti'

        ! this is the LOCAL paraview piece name, necessary for
        do iproc=0,nproc-1
            write (char_piece,'(i4.1)') iproc
            do ichar = 1,len(char_piece)
                if (char_piece(ichar:ichar)==' ') char_piece(ichar:ichar)='0'
            end do
            paraview_piecename(iproc)=trim(paraview_piecename(iproc))//'_proc'//char_piece//'.vti'
        end do

        !-----------------------------------------------------------------------
        ! who writes the new_res??
        if (i_printfile == 0 .or. i_printfile == 1) then
            ! only proc 0 writes
            if (myid == 0) iscrivo=.true.
        else if (i_printfile == 2) then
            ! all procs write
            iscrivo=.true.
        end if

        !-----------------------------------------------------------------------
        ! set output file names for paraview
        ! also for paraview all procs write, and in addition there is one header file
        ! that is handled only by processor 0

    !        if (i_paraview) then
    !            ! filename for iter
    !            if (print_iter_or_time == 0) then
    !                if (myid==0) then
    !                    filename_paraview=paraview_folder//'/'// &
    !                        paraview_header//char_ikiter//'.pvti'
    !                end if
    !                do iproc=0,nproc-1
    !                    write (char_piece,'(i4.1)') iproc
    !                    do ichar = 1,len(char_piece)
    !                        if (char_piece(ichar:ichar)==' ') char_piece(ichar:ichar)='0'
    !                    end do
    !                    paraview_piecename(iproc)=paraview_header//char_ikiter//'_'//char_piece//'.vti'
    !                end do
    !                filename_paraview_piece=paraview_folder//'/'// &
    !                    paraview_piecename(myid)
    !            ! filename for time
    !            else if (print_iter_or_time == 1) then
    !                if (myid==0) then
    !                    filename_paraview=paraview_folder//'/'// &
    !                        paraview_header//char_iktime//'.pvti'
    !                end if
    !                do iproc=0,nproc-1
    !                    write (char_piece,'(i4.1)') char_piece
    !                    do ichar = 1,len(char_piece)
    !                        if (char_piece(ichar:ichar)==' ') char_piece(ichar:ichar)='0'
    !                    end do
    !                    paraview_piecename(iproc)=paraview_header//char_iktime//'_'//char_piece//'.vti'
    !                end do
    !                filename_paraview_piece=paraview_folder//'/'// &
    !                    paraview_piecename(myid)
    !            end if
    !
    !        end if
        !-----------------------------------------------------------------------

    end subroutine write_preparation

    !***********************************************************************
    subroutine prepare_printdata()
        !***********************************************************************
        !     in these sub the variable that must be print are given to temp
        !     arrays like print_u.
        !     if iprint_file = 0  print_u is allocated on all the domain
        !     all procs known the flow field and just one will print
        !     if iprint_file = 1 print_u is allocated locally and it is used
        !     instead of directly u just for simplicity in the writing subroutine
        !
        use myarrays_velo3
        use myarrays_metri3
        use myarrays_density
        !
        use tipologia
        !
        use mpi

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k
        integer m
        integer ierr
        integer isc

        !integer requ1,requ2,reqv1,reqv2,reqw1,reqw2
        !integer reqannit1,reqannit2,reqakapt1,reqakapt2
        !integer reqr1,reqr2,reqf1,reqf2
        integer istatus(MPI_STATUS_SIZE)

        real, allocatable :: rho(:,:,:)
        !-----------------------------------------------------------------------
        if (i_printfile == 0) then

            allocate(print_u(0:jx+1,0:jy+1,0:jz+1))
            allocate(print_v(0:jx+1,0:jy+1,0:jz+1))
            allocate(print_w(0:jx+1,0:jy+1,0:jz+1))

            allocate(print_uc(0:jx,jy,jz))
            allocate(print_vc(jx,0:jy,jz))
            allocate(print_wc(jx,jy,0:jz))

            allocate(print_fi(0:jx+1,0:jy+1,0:jz+1))

            allocate(print_annit(0:jx+1,0:jy+1,0:jz+1))
            allocate(print_annitV(0:jx+1,0:jy+1,0:jz+1))
            allocate(print_akapt(0:jx+1,0:jy+1,0:jz+1))
            allocate(print_akaptV(0:jx+1,0:jy+1,0:jz+1))

            allocate(print_rhov(nscal,0:jx+1,0:jy+1,0:jz+1))

            !      startoutputtime=MPI_WTIME()
            !
            ! adesso sarebbe da mettere gli ALLGATHER - basterebbe
            ! comunicare solo a P0 che e' l'unico che scrive l'ouptut
            ! devo passare a P0 anche k=jz+1 di u v w fi rho
            !

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_annit(i,j,k)=annit(i,j,k)
                    end do
                end do
            end do

            allocate (annitcol((jx+2)*(jy+2)*jz/nproc))
            allocate (annittot((jx+2)*(jy+2)*jz))

            do k=kparasta,kparaend
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                        annitcol(m)=annit(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER( &
                annitcol(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                annittot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                        print_annit(i,j,k)=annittot(m)
                    end do
                end do
            end do

            deallocate(annitcol,annittot)

            if (myid==nproc-1) then
                call MPI_SSEND(print_annit(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,0,1061, &
                    MPI_COMM_WORLD,ierr)
            !      call MPI_WAIT(reqannit1,istatus,ierr)
            else if (myid==0) then
                call MPI_RECV(print_annit(0,0,jz+1),(jx+2)*(jy+2),MPI_REAL_SD, &
                    nproc-1,1061, &
                    MPI_COMM_WORLD,istatus,ierr)
            !      call MPI_WAIT(reqannit2,istatus,ierr)
            end if


            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_annitV(i,j,k)=annitV(i,j,k)
                    end do
                end do
            end do

            allocate (annitcol((jx+2)*(jy+2)*jz/nproc))
            allocate (annittot((jx+2)*(jy+2)*jz))

            do k=kparasta,kparaend
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                        annitcol(m)=annitV(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER( &
                annitcol(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                annittot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                        print_annitV(i,j,k)=annittot(m)
                    end do
                end do
            end do

            deallocate(annitcol,annittot)

            if (myid==nproc-1) then
                call MPI_SSEND(print_annitV(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,0,1061, &
                    MPI_COMM_WORLD,ierr)
            !      call MPI_WAIT(reqannit1,istatus,ierr)
            else if (myid==0) then
                call MPI_RECV(print_annitV(0,0,jz+1),(jx+2)*(jy+2),MPI_REAL_SD, &
                    nproc-1,1061, &
                    MPI_COMM_WORLD,istatus,ierr)
            !      call MPI_WAIT(reqannit2,istatus,ierr)
            end if

            !.................................................................

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_akapt(i,j,k)=akapt(1,i,j,k)
                    end do
                end do
            end do


            allocate (akaptcol((jx+2)*(jy+2)*jz/nproc))
            allocate (akapttot((jx+2)*(jy+2)*jz))

            do k=kparasta,kparaend
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                        akaptcol(m)=akapt(1,i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER( &
                akaptcol(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                akapttot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                        print_akapt(i,j,k)=akapttot(m)
                    end do
                end do
            end do

            deallocate(akaptcol,akapttot)

            if (myid==nproc-1) then
                call MPI_SSEND(print_akapt(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,0,1071, &
                    MPI_COMM_WORLD,ierr)
            !      call MPI_WAIT(reqakapt1,istatus,ierr)
            else if (myid==0) then
                call MPI_RECV(print_akapt(0,0,jz+1),(jx+2)*(jy+2),MPI_REAL_SD, &
                    nproc-1,1071, &
                    MPI_COMM_WORLD,istatus,ierr)
            !      call MPI_WAIT(reqakapt2,istatus,ierr)
            end if

            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_u(i,j,k)=u(i,j,k)
                    end do
                end do
            end do


            allocate (ucol((jx+2)*(jy+2)*jz/nproc))
            allocate (utot((jx+2)*(jy+2)*jz))

            do k=kparasta,kparaend
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                        ucol(m)=u(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(ucol(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                utot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                        print_u(i,j,k)=utot(m)
                    end do
                end do
            end do

            deallocate(ucol,utot)

            if (myid==nproc-1) then
                call MPI_SSEND(print_u(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,0,1011, &
                    MPI_COMM_WORLD,ierr)
            !      call MPI_WAIT(requ1,istatus,ierr)
            else if (myid==0) then
                call MPI_RECV(print_u(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,nproc-1,1011, &
                    MPI_COMM_WORLD,istatus,ierr)
            !      call MPI_WAIT(requ2,istatus,ierr)
            end if

            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_v(i,j,k)=v(i,j,k)
                    end do
                end do
            end do

            allocate (vcol((jx+2)*(jy+2)*jz/nproc))
            allocate (vtot((jx+2)*(jy+2)*jz))

            do k=kparasta,kparaend
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                        vcol(m)=v(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(vcol(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                vtot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                        print_v(i,j,k)=vtot(m)
                    end do
                end do
            end do

            deallocate(vcol,vtot)

            if (myid==nproc-1) then
                call MPI_SSEND(print_v(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,0,1021, &
                    MPI_COMM_WORLD,ierr)
            !      call MPI_WAIT(reqv1,istatus,ierr)
            else if (myid==0) then
                call MPI_RECV(print_v(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,nproc-1,1021, &
                    MPI_COMM_WORLD,istatus,ierr)
            !      call MPI_WAIT(reqv2,istatus,ierr)
            end if
            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_w(i,j,k)=w(i,j,k)
                    end do
                end do
            end do


            allocate (wcol((jx+2)*(jy+2)*jz/nproc))
            allocate (wtot((jx+2)*(jy+2)*jz))

            do k=kparasta,kparaend
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                        wcol(m)=w(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(wcol(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                wtot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                        print_w(i,j,k)=wtot(m)
                    end do
                end do
            end do

            deallocate(wcol,wtot)

            if (myid==nproc-1) then
                call MPI_SSEND(print_w(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,0,1031, &
                    MPI_COMM_WORLD,ierr)
            !      call MPI_WAIT(reqw1,istatus,ierr)
            else if (myid==0) then
                call MPI_RECV(print_w(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,nproc-1,1031, &
                    MPI_COMM_WORLD,istatus,ierr)
            !      call MPI_WAIT(reqw2,istatus,ierr)
            end if

            !.......................................................................
            allocate(rho(0:n1+1,0:n2+1,0:n3+1))

            do isc=1,nscal

                do k=kparasta-1,kparaend+1 !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            print_rhov(isc,i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do

                do k=kparasta-1,kparaend+1 !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            rho(i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do

                allocate (rhocol((jx+2)*(jy+2)*jz/nproc))
                allocate (rhotot((jx+2)*(jy+2)*jz))

                do k=kparasta,kparaend
                    do j=0,jy+1
                        do i=0,jx+1
                            m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                            rhocol(m)=rho(i,j,k)
                        end do
                    end do
                end do

                call MPI_ALLGATHER(rhocol(1),(jx+2)*(jy+2)*(jz/nproc), &
                    MPI_REAL_SD, &
                    rhotot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                    MPI_COMM_WORLD,ierr)

                do k=1,jz
                    do j=0,jy+1
                        do i=0,jx+1
                            m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                            rho(i,j,k)=rhotot(m)
                        end do
                    end do
                end do

                deallocate(rhocol,rhotot)

                if (myid==nproc-1) then
                    call MPI_SSEND(rho(0,0,jz+1),(jx+2)*(jy+2),MPI_REAL_SD,0,1041, &
                        MPI_COMM_WORLD,ierr)
                !      call MPI_WAIT(reqr1,istatus,ierr)
                else if (myid==0) then
                    call MPI_RECV(rho(0,0,jz+1),(jx+2)*(jy+2), &
                        MPI_REAL_SD,nproc-1,1041, &
                        MPI_COMM_WORLD,istatus,ierr)
                !      call MPI_WAIT(reqr2,istatus,ierr)
                end if

                do k=0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            print_rhov(isc,i,j,k)=rho(i,j,k)
                        end do
                    end do
                end do

            end do !nscal
            deallocate(rho)
            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_fi(i,j,k)=fi(i,j,k)
                    end do
                end do
            end do


            allocate (ficol((jx+2)*(jy+2)*jz/nproc))
            allocate (fitot((jx+2)*(jy+2)*jz))

            do k=kparasta,kparaend
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-kparasta))
                        ficol(m)=fi(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(ficol(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                fitot(1),(jx+2)*(jy+2)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy+1
                    do i=0,jx+1
                        m=i+1+(jx+2)*(j+(jy+2)*(k-1))
                        print_fi(i,j,k)=fitot(m)
                    end do
                end do
            end do

            deallocate(ficol,fitot)

            if (myid==nproc-1) then
                call MPI_SSEND(print_fi(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,0,1051, &
                    MPI_COMM_WORLD,ierr)
            !      call MPI_WAIT(reqf1,istatus,ierr)
            else if (myid==0) then
                call MPI_RECV(print_fi(0,0,jz+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,nproc-1,1051, &
                    MPI_COMM_WORLD,istatus,ierr)
            !      call MPI_WAIT(reqf2,istatus,ierr)
            end if

            !.......................................................................
            do k=kparasta,kparaend
                do j=1,jy
                    do i=0,jx
                        print_uc(i,j,k)=uc(i,j,k)
                    end do
                end do
            end do



            allocate (uccol((jx+1)*jy*jz/nproc))
            allocate (uctot((jx+1)*jy*jz))

            do k=kparasta,kparaend
                do j=1,jy
                    do i=0,jx
                        m=i+1+(jx+1)*(j-1+jy*(k-kparasta))
                        uccol(m)=uc(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(uccol(1),(jx+1)*jy*(jz/nproc),MPI_REAL_SD, &
                uctot(1),(jx+1)*jy*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=1,jy
                    do i=0,jx
                        m=i+1+(jx+1)*(j-1+jy*(k-1))
                        print_uc(i,j,k)=uctot(m)
                    end do
                end do
            end do

            deallocate(uccol,uctot)

            !.......................................................................
            do k=kparasta,kparaend
                do j=0,jy
                    do i=1,jx
                        print_vc(i,j,k)=vc(i,j,k)
                    end do
                end do
            end do

            allocate (vccol(jx*(jy+1)*jz/nproc))
            allocate (vctot(jx*(jy+1)*jz))

            do k=kparasta,kparaend
                do j=0,jy
                    do i=1,jx
                        m=i+jx*(j+(jy+1)*(k-kparasta))
                        vccol(m)=vc(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(vccol(1),jx*(jy+1)*(jz/nproc),MPI_REAL_SD, &
                vctot(1),jx*(jy+1)*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=0,jy
                    do i=1,jx
                        m=i+jx*(j+(jy+1)*(k-1))
                        print_vc(i,j,k)=vctot(m)
                    end do
                end do
            end do

            deallocate(vccol,vctot)
            !.......................................................................

            do k=kparasta-1,kparaend
                do j=1,jy
                    do i=1,jx
                        print_wc(i,j,k)=wc(i,j,k)
                    end do
                end do
            end do

            allocate (wccol(jx*jy*jz/nproc))
            allocate (wctot(jx*jz*jy))

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        m=i+jx*(j-1+jy*(k-kparasta))
                        wccol(m)=wc(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(wccol(1),jx*jy*(jz/nproc),MPI_REAL_SD, &
                wctot(1),jx*jy*(jz/nproc),MPI_REAL_SD, &
                MPI_COMM_WORLD,ierr)

            do k=1,jz
                do j=1,jy
                    do i=1,jx
                        m=i+jx*((j-1)+jy*(k-1))
                        print_wc(i,j,k)=wctot(m)
                    end do
                end do
            end do

            deallocate(wccol,wctot)


        !-----------------------------------------------------------------------
        !     move variable to print_u etc if each procs are going to print

        else if (i_printfile == 2) then

            allocate(print_u(0:jx+1,0:jy+1,kparasta-1:kparaend+1))
            allocate(print_v(0:jx+1,0:jy+1,kparasta-1:kparaend+1))
            allocate(print_w(0:jx+1,0:jy+1,kparasta-1:kparaend+1))

            allocate(    print_fi(0:jx+1,0:jy+1,kparasta-1:kparaend+1))
            allocate( print_annit(0:jx+1,0:jy+1,kparasta-1:kparaend+1))
            allocate(print_annitV(0:jx+1,0:jy+1,kparasta-1:kparaend+1))
            allocate( print_akapt(0:jx+1,0:jy+1,kparasta-1:kparaend+1))
            allocate(print_akaptV(0:jx+1,0:jy+1,kparasta-1:kparaend+1))

            allocate(print_rhov(nscal,0:jx+1,0:jy+1,kparasta-1:kparaend+1))

            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        print_u(i,j,k)  =  u(i,j,k)
                        print_v(i,j,k)  =  v(i,j,k)
                        print_w(i,j,k)  =  w(i,j,k)
                        print_fi(i,j,k) = fi(i,j,k)

                        print_annit(i,j,k) =   annit(i,j,k)
                        print_annitV(i,j,k) =  annitV(i,j,k)
                        !          print_akapt(i,j,k) =   akapt(1,i,j,k)
                        !         print_akaptV(i,j,k) =  akaptV(1,i,j,k)

                        do isc = 1,nscal
                            print_rhov(isc,i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do
            end do


            allocate(print_uc(0:jx,jy,kparasta:kparaend))
            allocate(print_vc(jx,0:jy,kparasta:kparaend))

            do k=kparasta,kparaend
                do j=1,jy
                    do i=0,jx
                        print_uc(i,j,k) = uc(i,j,k)
                    end do
                end do
            end do

            do k=kparasta,kparaend
                do j=0,jy
                    do i=1,jx
                        print_vc(i,j,k) = vc(i,j,k)
                    end do
                end do
            end do

            allocate(print_wc(jx,jy,kparasta-1:kparaend))
            do k=kparasta-1,kparaend
                do j=1,jy
                    do i=1,jx
                        print_wc(i,j,k) = wc(i,j,k)
                    end do
                end do
            end do

        end if ! if on i_printfile


        return

    end subroutine prepare_printdata

    !***********************************************************************
    subroutine write_newres(ti,bbx,dbbx,tipo)

        !***********************************************************************
        !     this sub write the new_res_form which contains the flow field for
        !     restart the computation or for analysis.
        !     depending on i_printfile only proc0 or all procs write their file
        !     depending on iformat_newres the output is written as
        !     iformat_newres = 0 ---> *
        !     iformat_newres = 1 ---> use a format from Agenerale.in
        !     iformat_newres = 2 ---> binary
        ! Paraview always uses standard format

        implicit none

        ! arguments
        real,intent(in) :: ti,bbx,dbbx
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! local variables declaration
        integer kstawrite,kendwrite
        integer i,j,k,n,isc
        !integer ierr

        !-----------------------------------------------------------------------
        !     start to write new_res file (only proc that have iscrivo true write)
        if (iscrivo) then

            write(*,*) 'Write newres on file: ',filename_newres

            ! with no format
            if (iformat_newres == 0) then

                open(new_res_file,file=filename_newres,status='unknown')
                write (new_res_file,*) nscal
                write (new_res_file,*) ti
                write (new_res_file,*) bbx,dbbx

                !  write the cartesian component
                if (i_printfile == 0) then
                    kstawrite = 0
                    kendwrite = jz+1
                end if
                if (i_printfile == 2) then
                    kstawrite = kparasta
                    kendwrite = kparaend
                    if (myid == 0) kstawrite = kparasta -1
                    if (myid == nproc-1) kendwrite = kparaend +1
                end if

                do k=kstawrite,kendwrite         !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            write (new_res_file,*) print_u(i,j,k)
                            write (new_res_file,*) print_v(i,j,k)
                            write (new_res_file,*) print_w(i,j,k)
                            write (new_res_file,*) print_fi(i,j,k)
                            do isc = 1,nscal
                                write (new_res_file,*) print_rhov(isc,i,j,k)
                            end do
                        end do
                    end do
                end do

                ! write the controvariant component
                if (i_printfile == 0) then
                    kstawrite = 1
                    kendwrite = jz
                end if
                if (i_printfile == 2) then
                    kstawrite = kparasta
                    kendwrite = kparaend
                end if


                do k=kstawrite,kendwrite !1,jz
                    do j=1,jy
                        do i=0,jx
                            write(new_res_file,*)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,jz
                    do j=0,jy
                        do i=1,jx
                            write(new_res_file,*)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid == 0) kstawrite = 0
                do k=kstawrite,kendwrite  !0,jz
                    do j=1,jy
                        do i=1,jx
                            write(new_res_file,*)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(new_res_file)

            !     file with annit
            !         write(rsti,'(i6)')itin
            !         open(new_res_file,file='annit'//prti,status='unknown')

            !         write(stringprint)ti
            !         do k=kstawrite,kendwrite !0,jz+1
            !         do j=0,jy+1
            !         do i=0,jx+1
            !            write(stringprint)print_annit(i,j,k),print_annitV(i,j,k)
            !         end do
            !         end do
            !         end do
            !         close(new_res_file)

            !-----------------------------------------------------------------------
            ! with format fmt_newres
            else if (iformat_newres == 1) then
                open(new_res_file,file=trim(filename_newres),status='unknown')
                write(new_res_file,*) nscal
                write(new_res_file,*)ti
                write(new_res_file,*)bbx,dbbx

                ! write the cartesiano component
                if (i_printfile == 0) then
                    kstawrite = 0
                    kendwrite = jz+1
                end if
                if (i_printfile == 2) then
                    kstawrite = kparasta
                    kendwrite = kparaend
                    if (myid == 0)kstawrite = kparasta -1
                    if (myid == nproc-1)kendwrite = kparaend +1
                end if

                do k=kstawrite,kendwrite         !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            write(new_res_file,fmt_newres)print_u(i,j,k), &
                                print_v(i,j,k), &
                                print_w(i,j,k), &
                                print_fi(i,j,k), &
                                (print_rhov(isc,i,j,k),isc=1,nscal)
                        end do
                    end do
                end do

                !        write the controvariant component
                if (i_printfile == 0) then
                    kstawrite = 1
                    kendwrite = jz
                end if
                if (i_printfile == 2) then
                    kstawrite = kparasta
                    kendwrite = kparaend
                end if


                do k=kstawrite,kendwrite !1,jz
                    do j=1,jy
                        do i=0,jx
                            write(new_res_file,fmt_newres)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,jz
                    do j=0,jy
                        do i=1,jx
                            write(new_res_file,fmt_newres)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid == 0) kstawrite = 0
                do k=kstawrite,kendwrite  !0,jz
                    do j=1,jy
                        do i=1,jx
                            write(new_res_file,fmt_newres)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(new_res_file)


            !     file with annit
            !         write(rsti,'(i6)')itin
            !         open(new_res_file,file='annit'//prti,status='unknown')

            !         write(stringprint)ti
            !         do k=kstawrite,kendwrite !0,jz+1
            !         do j=0,jy+1
            !         do i=0,jx+1
            !            write(stringprint)print_annit(i,j,k),print_annitV(i,j,k)
            !         end do
            !         end do
            !         end do
            !         close(new_res_file)

            !-----------------------------------------------------------------------
            !     binary
            else if (iformat_newres == 2) then

                open(new_res_file,file=filename_newres,status='unknown',form='unformatted')

                write(new_res_file)nscal
                write(new_res_file)ti
                write(new_res_file)bbx,dbbx

                !        write the cartesiano component
                if (i_printfile == 0) then
                    kstawrite = 0
                    kendwrite = jz+1
                end if
                if (i_printfile == 2) then
                    kstawrite = kparasta
                    kendwrite = kparaend
                    if (myid == 0)kstawrite = kparasta -1
                    if (myid == nproc-1)kendwrite = kparaend +1
                end if

                do k=kstawrite,kendwrite         !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            write(new_res_file)print_u(i,j,k), &
                                print_v(i,j,k), &
                                print_w(i,j,k), &
                                print_fi(i,j,k), &
                                (print_rhov(isc,i,j,k),isc=1,nscal)
                        end do
                    end do
                end do

                !        write the controvariant component
                if (i_printfile == 0) then
                    kstawrite = 1
                    kendwrite = jz
                end if
                if (i_printfile == 2) then
                    kstawrite = kparasta
                    kendwrite = kparaend
                end if


                do k=kstawrite,kendwrite !1,jz
                    do j=1,jy
                        do i=0,jx
                            write(new_res_file)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,jz
                    do j=0,jy
                        do i=1,jx
                            write(new_res_file)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid == 0) kstawrite = 0
                do k=kstawrite,kendwrite  !0,jz
                    do j=1,jy
                        do i=1,jx
                            write(new_res_file)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(new_res_file)

            !     file with annit
            !         write(rsti,'(i6)')itin
            !         open(new_res_file,file='annit'//prti,status='unknown')

            !         write(stringprint)ti
            !         do k=kstawrite,kendwrite !0,jz+1
            !         do j=0,jy+1
            !         do i=0,jx+1
            !            write(stringprint)print_annit(i,j,k),print_annitV(i,j,k)
            !         end do
            !         end do
            !         end do
            !         close(new_res_file)

            end if ! on the format
        end if ! if iscrivo


        !-----------------------------------------------------------------------
        !     these arrays have been allocated in prepare_printdata subroutine
        !
        deallocate(print_u)
        deallocate(print_v)
        deallocate(print_w)
        deallocate(print_uc)
        deallocate(print_vc)
        deallocate(print_wc)
        deallocate(print_fi)
        deallocate(print_annit)
        deallocate(print_annitv)
        deallocate(print_akapt)
        deallocate(print_akaptv)
        deallocate(print_rhov)

    end subroutine write_newres

    !***********************************************************************
    subroutine write_medietempo()

        use mpi

        implicit none

        integer ifmt_mean
        character*2 char_fmt
        character*12 fmt_mean

        integer :: ierr
        integer :: i,j,k,n

        !-----------------------------------------------------------------------
        ! write the medietempo for time statistics

        ifmt_mean = 4 + nscal + 6 + nscal + nscal
        write(char_fmt,'(i2)')ifmt_mean
        fmt_mean = '('//char_fmt//'e18.10)'

        open(medietempo_file,file=trim(filename_medietempo),status='new',action='write')
        write(*,*) 'Write medietempo file: ',trim(filename_medietempo)
        write(medietempo_file,*)nscal,nproc,re
        write(medietempo_file,*)jx,jy,jz
        write(medietempo_file,*)deltatempo
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    !        write(medietempo_file,1001)umed(i,j,k,1),umed(i,j,k,2),
                    !     >         umed(i,j,k,3),umed(i,j,k,4),
                    !     >         umed(i,j,k,5), uutau(i,j,k,1),
                    !     >         uutau(i,j,k,2),uutau(i,j,k,3),
                    !     >         uutau(i,j,k,4),uutau(i,j,k,5),
                    !     >         uutau(i,j,k,6),uutau(i,j,k,7)


                    write(medietempo_file,fmt_mean)(umed(i,j,k,n),n=1,4+nscal), &       !u,v,w,fi,r_n
                        (uutau(i,j,k,n),n=1,6+nscal+nscal)  !uu_ii, uu_ij, rr_n, rv_n



                end do
            end do
        end do
        close(medietempo_file)

        !      end if



        !1001    format(12e18.10)

        !     the data are used only for average in the time window and not globally
        if (.not.i_cumulative) then
            deltatempo = 0.
            umed = 0.
            uutau = 0.
        end if

    end subroutine write_medietempo

    !***********************************************************************
    subroutine write_tracers(ti)

        use mysettings, only: sonde,nsonde

        implicit none

        real, intent(in) :: ti
        !
        integer :: isc,isonde

        !use myarrays
        !-----------------------------------------------------------------------
        !     print tracers
        !
        do isonde = 1,nsonde
            if (sonde(3,isonde)>=kparasta.and.sonde(3,isonde)<=kparaend) then
                write(sonde_file_base+isonde,'(10e18.10)') ti, &
                    print_u(sonde(1,isonde),sonde(2,isonde),sonde(3,isonde)), &
                    print_v(sonde(1,isonde),sonde(2,isonde),sonde(3,isonde)), &
                    print_w(sonde(1,isonde),sonde(2,isonde),sonde(3,isonde)), &
                    print_fi(sonde(1,isonde),sonde(2,isonde),sonde(3,isonde)), &
                    (print_rhov(isc,sonde(1,isonde),sonde(2,isonde),sonde(3,isonde)), &
                    isc=1,nscal)
            end if
        end do

    !1000    format(10e18.10)

    end subroutine write_tracers

    !***********************************************************************
    subroutine write_paraview(tipo)

        use myarrays_metri3, only: x,y,z
        use myarrays_velo3, only: u,v,w,fi,rhov
        use myarrays_ibm, only: bodyforce
        !use myarrays_levelset2, only: phi0, phi_temp, fext

        use mpi

        implicit none

        !
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! local variables declaration
        integer piecefile
        integer i,j,k,isc,iproc
        integer z_piece_begin(0:nproc-1),z_piece_end(0:nproc-1)
        integer x_piece_begin,x_piece_end
        integer y_piece_begin,y_piece_end
        integer ierr
        ! we need kparasta and kparaend saved at the root processor
        ! for this, a gather is needed

        !piece_begin(myid)=kparasta
        !piece_end(myid)=kparaend
        call MPI_GATHER(kparasta-1,1,MPI_INTEGER,z_piece_begin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(kparaend+1,1,MPI_INTEGER,z_piece_end,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        x_piece_begin=0
        x_piece_end=jx+1
        y_piece_begin=0
        y_piece_end=jy+1
        z_piece_begin(myid)=kparasta-1
        z_piece_end(myid)=kparaend+1

        ! header file is written only by root processor
        if (myid == 0) then
            write(*,*) 'Write Paraview header file ',trim(filename_paraview)
            open(paraview_file,file=trim(filename_paraview),status='unknown')
            write(paraview_file,'(A)')'<?xml version="1.0"?>'
            write(paraview_file,'(A)')'<VTKFile type="PImageData" version="0.1">'
            write(paraview_file,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <PImageData WholeExtent="',&
                x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(0),' ',z_piece_end(nproc-1),&
                '" GhostLevel="1" Origin="0.0 0.0 0.0" Spacing="1.0 1.0 1.0">' ! spacing is missing
            write(paraview_file,'(A)')'    <PPointData>'
            write(paraview_file,'(A)')'      <PDataArray type="Float64" Name="warp" NumberOfComponents="3"/>'
            !            write(paraview_file,'(A)')'      <PDataArray type="Int32" Name="y"/>'
            !            write(paraview_file,'(A)')'      <PDataArray type="Int32" Name="z"/>'
            write(paraview_file,'(A)')'      <PDataArray type="Float64" Name="vel" NumberOfComponents="3"/>'
            write(paraview_file,'(A)')'      <PDataArray type="Float64" Name="pres"/>'
            !            write(paraview_file,'(A)')'      <PDataArray type="Float64" Name="phi_temp"/>'
            !            write(paraview_file,'(A)')'      <PDataArray type="Float64" Name="fext"/>'
            !            write(paraview_file,'(A)')'      <PDataArray type="Float64" Name="phi0"/>'
            !if (bodyforce == 1) then
            write(paraview_file,'(A)')'      <PDataArray type="Int32" Name="tipo"/>'
            !end if
            !            if (nscal>0) then
            !                do isc = 1,nscal
            !                    write(paraview_file,'(A,I0.1,A)')'      <PDataArray type="Float64" Name="scal',isc,'"/>'
            !                end do
            !            end if
            write(paraview_file,'(A)')'    </PPointData>'
            write(paraview_file,'(A)')'    <PCellData>'
            write(paraview_file,'(A)')'    </PCellData>'
            !write(paraview_file,'(A,I0.1,A,I0.1,A)')'    <Piece Extent="0 ',jx+1,' 0 ',jy+1,' ',piece_begin,' ',piece_end,'" Source="para_iter000010_0.vti"/>'
            do iproc=0,nproc-1
                write(paraview_file,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,A,A)') '    <Piece Extent="',&
                    x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',&
                    z_piece_begin(iproc),' ',z_piece_end(iproc),&
                    '" Source="',trim(paraview_piecename(iproc)),'">'
                write(paraview_file,'(A)')'    </Piece>'
            end do

            write(paraview_file,'(A)')'  </PImageData>'
            write(paraview_file,'(A)')'</VTKFile>'

            close(paraview_file)
        end if

        ! every processor writes a piece file
        piecefile=paraview_piece_file+myid
        write(*,*) 'Write Paraview piece file ',trim(paraview_piecename(myid))
        open(piecefile,file=trim(filename_paraview_piece),status='unknown')
        write(piecefile,'(A)')'<?xml version="1.0"?>'
        write(piecefile,'(A)')'<VTKFile type="ImageData" version="0.1">'
        write(piecefile,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <ImageData WholeExtent="', &
            x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
            '">'
        write(piecefile,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <Piece Extent="',&
            x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
            '">'
        write(piecefile,'(A)')'    <PointData>'
        write(piecefile,'(A)')'      <DataArray type="Float64" Name="warp" NumberOfComponents="3">'
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    write(piecefile,'(ES14.7,A,ES14.7,A,ES14.7,A)',advance='no') &
                        x(i,j,k)-1.0*i,' ',y(i,j,k)-1.0*j,' ',z(i,j,k)-1.0*k,' '
                end do
            end do
        end do
        write(piecefile,*)
        write(piecefile,'(A)')'      </DataArray>'
        write(piecefile,'(A)')'      <DataArray type="Float64" Name="vel" NumberOfComponents="3">'
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    write(piecefile,'(ES14.7,A,ES14.7,A,ES14.7,A)',advance='no') &
                        u(i,j,k),' ',v(i,j,k),' ',w(i,j,k),' '
                end do
            end do
        end do
        write(piecefile,*)
        write(piecefile,'(A)')'      </DataArray>'
        write(piecefile,'(A)')'      <DataArray type="Float64" Name="pres">'
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    write(piecefile,'(ES14.7,A)',advance='no') fi(i,j,k),' '
                end do
            end do
        end do
        write(piecefile,*)
        write(piecefile,'(A)')'      </DataArray>'
        !            write(piecefile,'(A)')'      <DataArray type="Float64" Name="phi_temp">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(piecefile,'(ES14.7,A)',advance='no') phi_temp(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(piecefile,*)
        !            write(piecefile,'(A)')'      </DataArray>'
        !            write(piecefile,'(A)')'      <DataArray type="Float64" Name="Fext">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(piecefile,'(ES14.7,A)',advance='no') fext(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(piecefile,*)
        !            write(piecefile,'(A)')'      </DataArray>'
        !            write(piecefile,'(A)')'      <DataArray type="Float64" Name="phi0">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(piecefile,'(ES14.7,A)',advance='no') phi0(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(piecefile,*)
        !            write(piecefile,'(A)')'      </DataArray>'
        !if (bodyforce == 1) then
        write(piecefile,'(A)')'      <DataArray type="Int32" Name="tipo">'
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    write(piecefile,'(I0.1,A)',advance='no') tipo(i,j,k),' '
                end do
            end do
        end do
        write(piecefile,'(A)')'      </DataArray>'
        !            write(fileunit,'(A)')'      <DataArray type="Float64" Name="surfNorm" NumberOfComponents="3">'
        !            do k = 1,nz
        !                do j = 1,ny
        !                    do i = 1,nx
        !                        write(fileunit,'(ES14.7,A,ES14.7,A,ES14.7,A)',advance='no') &
        !                            normalVectorX(i,j,k),' ',normalVectorY(i,j,k),' ',normalVectorZ(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(fileunit,*)
        !            write(fileunit,'(A)')'      </DataArray>'
            !end if
        !            if (nscal>0) then
        !                do isc = 1,nscal
        !                write(piecefile,'(A,I0.1,A)')'      <DataArray type="Float64" Name="scal',isc,'">'
        !                do k=kparasta-1,kparaend+1
        !                    do j=0,jy+1
        !                        do i=0,jx+1
        !                            write(piecefile,'(ES14.7,A)',advance='no') rhov(isc,i,j,k),' '
        !                        end do
        !                    end do
        !                end do
        !                write(piecefile,*)
        !            end do
        !        end if
        !        write(piecefile,'(A)')'      </DataArray>'
        write(piecefile,'(A)')'    </PointData>'
        write(piecefile,'(A)')'    <CellData>'
        write(piecefile,'(A)')'    </CellData>'
        write(piecefile,'(A)')'    </Piece>'
        write(piecefile,'(A)')'  </ImageData>'
        write(piecefile,'(A)')'</VTKFile>'

        close(piecefile)

    end subroutine write_paraview

    !***********************************************************************
    subroutine output_finalize()

        use mysettings, only: lagr,nsonde

        implicit none

        integer :: isonde

        if (myid==0) then

            if (lagr==0) then
                close(turbo_file)   ! turbo.out
                close(subdis_file)
                close(subscl_file)
                close(subrho_file)
            end if
            close(bulkvel_file)   ! bulk velocity nel tempo
            close(maxvel_file)
            close(encheck_file)
            !      close(35)   ! eddy viscosity
            !close(50)

            do isonde = 1,nsonde
                close(sonde_file_base+isonde)
            end do

        end if

    end subroutine output_finalize

    !***********************************************************************
    subroutine print_energy(tipo)

        use myarrays_metri3, only: fluid_vol,giac
        use myarrays_velo3, only: u,v,w
        use parti, only: ti

        use tipologia
        use mpi

        implicit none

        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        ! for energy coputation
        real :: enck,enck_loc
        integer :: i,j,k
        integer :: ierr

        !-----------------------------------------------------------------------
        !     print kinetic energy in time
        !
        enck_loc=0.

        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    if (tipo(i,j,k)/=0) then
                        enck_loc=enck_loc+.5*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))*giac(i,j,k)
                    end if
                end do
            end do
        end do

        ! reduce operation on enck and vol
        enck = 0.
        call MPI_REDUCE(enck_loc,enck,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (myid==0) then
            enck=enck/fluid_vol
            open(encheck_file,file=filename_encheck,status='OLD',action='WRITE',position='APPEND')
            write(encheck_file,*) ti,enck
            close(encheck_file)
        end if
    end subroutine print_energy

    !***********************************************************************
    subroutine print_velocity(tipo)

        use myarrays_metri3, only: fluid_vol,giac
        use myarrays_velo3, only: u,v,w
        use parti, only: ti

        use tipologia
        use mpi

        implicit none

        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        real :: u_bulk,u_bulk_tot,v_bulk,v_bulk_tot,w_bulk,w_bulk_tot
        real :: u_max,u_max_tot,v_max,v_max_tot,w_max,w_max_tot
        integer :: i,j,k
        integer :: ierr

        ! bulk velocity
        u_bulk=0.
        u_bulk_tot=0.
        v_bulk=0.
        v_bulk_tot=0.
        w_bulk=0.
        w_bulk_tot=0.

        ! max velocity
        u_max=0.
        u_max_tot=0.
        v_max=0.
        v_max_tot=0.
        w_max=0.
        w_max_tot=0.

        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    if (tipo(i,j,k)/=0) then
                        u_bulk=u_bulk+u(i,j,k)*giac(i,j,k)
                        v_bulk=v_bulk+v(i,j,k)*giac(i,j,k)
                        w_bulk=w_bulk+w(i,j,k)*giac(i,j,k)

                        u_max=max(u_max,abs(u(i,j,k)))
                        v_max=max(v_max,abs(v(i,j,k)))
                        w_max=max(w_max,abs(w(i,j,k)))
                    end if
                end do
            end do
        end do

        ! u_bulk sum between procs
        call MPI_REDUCE(u_bulk,u_bulk_tot,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(v_bulk,v_bulk_tot,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(w_bulk,w_bulk_tot,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        ! u_max reduce between procs
        call MPI_REDUCE(u_max,u_max_tot,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(v_max,v_max_tot,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(w_max,w_max_tot,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)

        if (myid==0) then

            u_bulk_tot=u_bulk_tot/fluid_vol

            write(*,*)'u_bulk:',u_bulk_tot,', u_max:',u_max_tot

            open(bulkvel_file,file=filename_bulkvel,status='old',action='write',position='append')
            write(bulkvel_file,*) ti,u_bulk_tot,v_bulk_tot,w_bulk_tot
            close(bulkvel_file)

            open(maxvel_file,file=filename_maxvel,status='old',action='write',position='append')
            write(maxvel_file,*) ti,u_max_tot,v_max_tot,w_max_tot
            close(maxvel_file)
        end if

    end subroutine print_velocity

    !***********************************************************************
    subroutine print_drag()

        use myarrays_ibm, only: x_force_sphere,y_force_sphere,z_force_sphere
        use parti, only: ti

        implicit none

        real :: u_bulk,u_bulk_tot
        integer :: i,j,k
        integer :: ierr

        if (myid==0) then

            write(*,*)'u_bulk:',u_bulk_tot

            open(drag_file,file=filename_drag,status='old',action='write',position='append')
            write(drag_file,*) ti,x_force_sphere(1),y_force_sphere(1),z_force_sphere(1)
            close(drag_file)
        end if

    end subroutine print_drag

    !***********************************************************************
    subroutine update_medie()

        use myarrays_velo3, only: u,v,w,fi,rhov

        implicit none

        integer :: i,j,k,n

        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    umed(i,j,k,1) = umed(i,j,k,1) +      u(i,j,k)*dt          ! u media * Tempo totale
                    umed(i,j,k,2) = umed(i,j,k,2) +      v(i,j,k)*dt          ! v media * Tempo totale
                    umed(i,j,k,3) = umed(i,j,k,3) +      w(i,j,k)*dt          ! w media * Tempo totale
                    umed(i,j,k,4) = umed(i,j,k,4) +      fi(i,j,k)*dt

                    uutau(i,j,k,1)=uutau(i,j,k,1)+u(i,j,k)*u(i,j,k)*dt    ! u*u    * Tempo totale
                    uutau(i,j,k,2)=uutau(i,j,k,2)+u(i,j,k)*v(i,j,k)*dt    ! u*v    * Tempo totale
                    uutau(i,j,k,3)=uutau(i,j,k,3)+u(i,j,k)*w(i,j,k)*dt    ! u*w    * Tempo totale
                    uutau(i,j,k,4)=uutau(i,j,k,4)+v(i,j,k)*v(i,j,k)*dt    ! v*v    * Tempo totale
                    uutau(i,j,k,5)=uutau(i,j,k,5)+v(i,j,k)*w(i,j,k)*dt    ! v*w    * Tempo totale
                    uutau(i,j,k,6)=uutau(i,j,k,6)+w(i,j,k)*w(i,j,k)*dt    ! w*w    * Tempo totale


                    do n=1,nscal
                        umed(i,j,k,4+n) = umed(i,j,k,4+n) + rhov(n,i,j,k)*dt          ! r media * Tempo totale

                        uutau(i,j,k,6+n)= uutau(i,j,k,6+n) + rhov(n,i,j,k)*rhov(n,i,j,k)*dt  ! r*r

                        uutau(i,j,k,6+nscal+n)= uutau(i,j,k,6+nscal+n) + rhov(n,i,j,k)*v(i,j,k)*dt  ! r*v   * Tempo totale

                    end do

                end do
            end do
        end do

        !     update the total time
        deltatempo = deltatempo + dt

        !     print the data
    !        if (ktime==niter) then
    !            if (myid<10) then
    !                write(idproc,'(i1)')myid
    !                filepntempo='medietempo0'//idproc//'.dat'
    !            else
    !                write(idproc2,'(i2)')myid
    !                filepntempo='medietempo'//idproc2//'.dat'
    !            end if
    !            open(1000,file=filepntempo,status='unknown')
    !            write(1000,1001)deltatempo
    !            do k=kparasta,kparaend
    !                do j=1,jy
    !                    do i=1,jx
    !                        write(1000,1001)umed(i,j,k,1),umed(i,j,k,2),
    !                        >                    umed(i,j,k,3),umed(i,j,k,4),
    !                        >                    umed(i,j,k,5), uutau(i,j,k,1),
    !                        >                    uutau(i,j,k,2),uutau(i,j,k,3),
    !                        >                    uutau(i,j,k,4),uutau(i,j,k,5),
    !                        >                    uutau(i,j,k,6),uutau(i,j,k,7)
    !                    end do
    !                end do
    !            end do
    !            close(1000)
    !        end if
    !1001    format(12e18.10)


    end subroutine update_medie

    !***********************************************************************
    subroutine set_simulation_folder(simulation_folder_string) bind(C,name='set_simulation_folder')

        use, intrinsic :: iso_c_binding

        implicit none

        character(kind=c_char), dimension(*), intent(IN) :: simulation_folder_string
        integer :: i,l


        do l=1,500
            if (simulation_folder_string(l) == C_NULL_CHAR) exit
            result_folder(l:l)=simulation_folder_string(l)
        end do
        result_folder(l:500)=" "


        print *, 'In Fortran:'
        print *, 'simulation folder: ', trim(result_folder)
    end subroutine set_simulation_folder

end module output_module

