!***********************************************************************
module output_module
    !***********************************************************************
    !     contains arrays and settings for printing procedure

    use mysending, only: kparaend,kparasta,myid,nproc,deepl,deepr
    !
    use scala3, only: jx,jy,jz,nscal,re,n1,n2,n3,dt
            !
    use print

    implicit none

    ! flags
    integer print_iter_or_time
    integer i_printfile
    integer i_paraview
    integer i_cumulative
    integer iscrivo
    integer iformat_newres
    integer iformat_grid
    integer ifolder

    ! information output files
    integer,parameter  :: info_run_file=11

    ! formats
    character*12 :: string_newres_format
    character*50 :: new_res_step_folder
    character*12 :: string_grid_format
    character*20 :: fmt_newres

    ! file names
    character*120               :: filename_newres,filename_medietempo,filename_paraview,filename_movie,filename_paraview_piece
    character*120,allocatable   :: paraview_piecename(:)


    ! folders to store results (see output_init)
    character(len=23) :: result_folder ! length = length(result_folder_name)+length('results/')
    character(len=34) :: new_res_folder
    character(len=36) :: paraview_folder

    character(len=4) :: char_myid

    real i_time
    real ti_start
      
    ! variables for printing
    real, allocatable :: umed(:,:,:,:) !(1:n1,1:n2,kparasta:kparaend,4)
    real, allocatable :: uutau(:,:,:,:) !(1:n1,1:n2,kparasta:kparaend,7)
    real deltatempo

    real,allocatable :: print_u(:,:,:)
    real,allocatable :: print_v(:,:,:)
    real,allocatable :: print_w(:,:,:)
      
    real,allocatable :: print_uc(:,:,:)
    real,allocatable :: print_vc(:,:,:)
    real,allocatable :: print_wc(:,:,:)
      
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

    !***********************************************************************
    subroutine output_init()

        use mysettings, only: lagr

        use mpi

        implicit none

        character(len=12),parameter :: paraview_folder_name = 'paraviewData'
        character(len=10),parameter :: new_res_folder_name = 'newResData'
        character(len=7),parameter  :: global_result_folder_name = 'results'
        character(len=8)            :: date
        character(len=6)            :: time
        character(len=15)           :: result_folder_name ! length = date+time+1

        integer ierr
        integer ichar

        !-----------------------------------------------------------------------
        ! check date and time for the result folder creation
        call date_and_time(date,time)


        if (myid == 0) then

            write(*,*) 'Initialize output step'

            if (lagr == 0) then
                write(*,*) 'Warning: in output files no subrid section has been implemented'
            end if

            ! global result folder
            write(*,*) 'Check if global result folder exist'
            call system('mkdir '//trim(global_result_folder_name))

        end if

        ! result folder
        result_folder_name=date//'_'//time
        result_folder=global_result_folder_name//'/'//result_folder_name
        if (myid == 0) then
            write(*,*) 'Create result folder (format=date_time)'
            call system('mkdir '//trim(result_folder))
        end if

        ! paraview folder
        if (i_paraview == 1) then
            paraview_folder=result_folder//'/'//paraview_folder_name
            if (myid == 0) then
                call system('mkdir '//paraview_folder)
                write(*,*) 'Created folder ',paraview_folder
            end if

        end if

        ! newRes folder
        if (1 == 1) then
            new_res_folder=result_folder//'/'//new_res_folder_name
            if (myid == 0) then
                call system('mkdir '//new_res_folder)
                write(*,*) 'Created folder ',new_res_folder
            end if
        end if


        if (i_paraview == 1) then
            if (myid == 0) write(*,*) 'Preparing paraview output split into ',nproc,' pieces'
            allocate(paraview_piecename(0:nproc-1))
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


        ! Files for information output -------------------------------------
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        open(info_run_file,file=trim(result_folder)//'/'//'info_run.txt',status='unknown')

    end subroutine output_init

    !***********************************************************************
    subroutine output_step(count_print_time,dbbx,tipo)
        !-----------------------------------------------------------------------
        !++++++++++  WRITE OUTPUT  +++++++++++++++++++++++++++++++++++++++++++++

        use parti, only: alx, alz, ti
        use mysettings, only: bbx, niter
        use convex, only: bulk
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

            ! set the folder for the output and the filename
            call write_preparation(ti)

            ! write output new_res_form without mpi procedure
            if (i_printfile == 0 .or. i_printfile == 2 .or. i_paraview==1) then

                ! collect data to print
                call prepare_printdata

                ! print data on file
                call write_output(ti,bbx,dbbx,tipo)

            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            !     check if I can write another new_res before ending the
            !     number of iteration
            if (print_iter_or_time == 1) then
                !         time_evaluation = (niter-ktime)*dt

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

    end subroutine

    !***********************************************************************
    subroutine write_preparation(ti)
        ! this subroutine prepares the file names and the folder

        use mpi

        implicit none

        ! file names
        character(len=12),parameter :: filename_header= 'new_res_iter'
        character(len=18),parameter :: filemedie_header='meantime_iter_proc'
        character(len=9),parameter  :: paraview_header ='para_iter'

        ! arguments
        real,intent(in) :: ti

        !-----------------------------------------------------------------------
        !  local variables declaration
        integer ichar,iproc
        integer ierr

        character(len=6)  :: char_ikiter
        character(len=13) :: char_iktime
        character(len=4)  :: char_piece


        !-----------------------------------------------------------------------
        if (myid == 0) write(*,*) 'OUTPUT STEP BEGIN ---------------------------'

        ! set the folder for the output

        ! print for iterations
        if (print_iter_or_time == 0) then

            write (char_ikiter,'(i6)') ktime
            do ichar = 1,len(char_ikiter)
                if (char_ikiter(ichar:ichar)==' ') char_ikiter(ichar:ichar)='0'
            end do

            if (ifolder == 1) then
                new_res_step_folder = new_res_folder//'/RUN_ITER_'//char_ikiter
            else
                new_res_step_folder= new_res_folder

            end if

        ! print for time
        else if (print_iter_or_time == 1) then

            write(char_iktime,'(1e13.7)') ti ! INT(ti-ti_start)
            write(*,*) char_iktime

            do ichar = 1,len(char_iktime)
                if (char_iktime(ichar:ichar)==' ') char_iktime(ichar:ichar)='0'
            end do

            if (ifolder == 1) then
                new_res_step_folder = new_res_folder//'/RUN_TIME_'//char_iktime//'sec'
            else
                new_res_step_folder= new_res_folder
            end if

        end if

        !rid_folder = len_trim(new_res_step_folder)
        if (ifolder == 1 .and. myid == 0) then
            call system('mkdir '//trim(new_res_step_folder))
            write(*,*) 'Created folder "',trim(new_res_step_folder),'"'
        end if

        ! flag for telling every processor if it must write
        iscrivo = 0

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !-----------------------------------------------------------------------
        ! set the output filename for new_res

        ! one proc write
        if (i_printfile == 0 .or. i_printfile == 1) then
            ! filename for iter
            if (print_iter_or_time == 0) then
                filename_newres=trim(new_res_step_folder)//'/'// &
                    filename_header//char_ikiter//'.dat'
            end if
            ! filename for time
            if (print_iter_or_time == 1) then
                filename_newres=trim(new_res_step_folder)//'/'// &
                    filename_header//char_iktime//'sec.dat'
            end if
            ! only procs 0 write
            if (myid == 0) iscrivo = 1

        end if

        ! all procs write
        if (i_printfile == 2) then
            ! filename for iter
            if (print_iter_or_time == 0) then
                filename_newres=trim(new_res_step_folder)//'/'// &
                    filename_header//char_ikiter//'_proc'//char_myid//'.dat'
            end if
            ! filename for time
            if (print_iter_or_time == 1) then
                filename_newres=trim(new_res_step_folder)//'/'// &
                    filename_header//char_iktime//'sec_proc'//char_myid//'.dat'
            end if
            ! all procs write
            iscrivo = 1
        end if


        !-----------------------------------------------------------------------
        ! set the output filename for medietempo
        filename_medietempo = trim(new_res_step_folder)//'/'// &
            filemedie_header//char_myid//'.dat'

        !-----------------------------------------------------------------------
        ! set output file names for paraview
        ! also for paraview all procs write, and in addition there is one header file
        ! that is handled only by processor 0
        if (i_paraview == 1) then
            ! filename for iter
            if (print_iter_or_time == 0) then
                if (myid==0) then
                    filename_paraview=paraview_folder//'/'// &
                        paraview_header//char_ikiter//'.pvti'
                end if
                do iproc=0,nproc-1
                    write (char_piece,'(i4.1)') iproc
                    do ichar = 1,len(char_piece)
                        if (char_piece(ichar:ichar)==' ') char_piece(ichar:ichar)='0'
                    end do
                    paraview_piecename(iproc)=paraview_header//char_ikiter//'_'//char_piece//'.vti'
                end do
                filename_paraview_piece=paraview_folder//'/'// &
                    paraview_piecename(myid)
            ! filename for time
            else if (print_iter_or_time == 1) then
                if (myid==0) then
                    filename_paraview=paraview_folder//'/'// &
                        paraview_header//char_iktime//'.pvti'
                end if
                do iproc=0,nproc-1
                    write (char_piece,'(i4.1)') char_piece
                    do ichar = 1,len(char_piece)
                        if (char_piece(ichar:ichar)==' ') char_piece(ichar:ichar)='0'
                    end do
                    paraview_piecename(iproc)=paraview_header//char_iktime//'_'//char_piece//'.vti'
                end do
                filename_paraview_piece=paraview_folder//'/'// &
                    paraview_piecename(myid)
            end if

        end if



        !-----------------------------------------------------------------------


    end subroutine write_preparation

    !***********************************************************************
    subroutine prepare_printdata
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
    subroutine write_output(ti,bbx,dbbx,tipo)

        !***********************************************************************
        !     this sub write the new_res_form which contains the flow field for
        !     restart the computation or for analysis.
        !     depending on i_printfile only proc0 or all procs write their file
        !     depending on iformat_newres the output is written as
        !     iformat_newres = 0 ---> *
        !     iformat_newres = 1 ---> use a format from Agenerale.in
        !     iformat_newres = 2 ---> binary
        ! Paraview always usus standard format

        implicit none

        ! arguments
        real,intent(in) :: ti,bbx,dbbx
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! local variables declaration
        integer fileunit
        integer kstawrite,kendwrite
        integer i,j,k,n,isc
        !integer ierr

        integer ifmt_mean
        character*2 char_fmt
        character*12 fmt_mean

        fileunit= 14
        !-----------------------------------------------------------------------
        !     start to write new_res file (only proc that have iscrivo=1 write)
        if (iscrivo == 1) then

            write(*,*) 'Write newres on file: ',filename_newres

            ! with no format
            if (iformat_newres == 0) then

                open(fileunit,file=filename_newres,status='unknown')
                write (fileunit,*) nscal
                write (fileunit,*) ti
                write (fileunit,*) bbx,dbbx

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
                            write (fileunit,*) print_u(i,j,k)
                            write (fileunit,*) print_v(i,j,k)
                            write (fileunit,*) print_w(i,j,k)
                            write (fileunit,*) print_fi(i,j,k)
                            do isc = 1,nscal
                                write (fileunit,*) print_rhov(isc,i,j,k)
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
                            write(fileunit,*)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,jz
                    do j=0,jy
                        do i=1,jx
                            write(fileunit,*)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid == 0) kstawrite = 0
                do k=kstawrite,kendwrite  !0,jz
                    do j=1,jy
                        do i=1,jx
                            write(fileunit,*)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(fileunit)

            !     file with annit
            !         write(rsti,'(i6)')itin
            !         open(fileunit,file='annit'//prti,status='unknown')

            !         write(stringprint)ti
            !         do k=kstawrite,kendwrite !0,jz+1
            !         do j=0,jy+1
            !         do i=0,jx+1
            !            write(stringprint)print_annit(i,j,k),print_annitV(i,j,k)
            !         end do
            !         end do
            !         end do
            !         close(fileunit)

            !-----------------------------------------------------------------------
            ! with format fmt_newres
            else if (iformat_newres == 1) then
                open(fileunit,file=trim(filename_newres),status='unknown')
                write(fileunit,*) nscal
                write(fileunit,*)ti
                write(fileunit,*)bbx,dbbx

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
                            write(fileunit,fmt_newres)print_u(i,j,k), &
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
                            write(fileunit,fmt_newres)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,jz
                    do j=0,jy
                        do i=1,jx
                            write(fileunit,fmt_newres)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid == 0) kstawrite = 0
                do k=kstawrite,kendwrite  !0,jz
                    do j=1,jy
                        do i=1,jx
                            write(fileunit,fmt_newres)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(fileunit)


            !     file with annit
            !         write(rsti,'(i6)')itin
            !         open(fileunit,file='annit'//prti,status='unknown')

            !         write(stringprint)ti
            !         do k=kstawrite,kendwrite !0,jz+1
            !         do j=0,jy+1
            !         do i=0,jx+1
            !            write(stringprint)print_annit(i,j,k),print_annitV(i,j,k)
            !         end do
            !         end do
            !         end do
            !         close(fileunit)

            !-----------------------------------------------------------------------
            !     binary
            else if (iformat_newres == 2) then

                open(fileunit,file=filename_newres,status='unknown',form='unformatted')

                write(fileunit)nscal
                write(fileunit)ti
                write(fileunit)bbx,dbbx

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
                            write(fileunit)print_u(i,j,k), &
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
                            write(fileunit)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,jz
                    do j=0,jy
                        do i=1,jx
                            write(fileunit)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid == 0) kstawrite = 0
                do k=kstawrite,kendwrite  !0,jz
                    do j=1,jy
                        do i=1,jx
                            write(fileunit)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(fileunit)

            !     file with annit
            !         write(rsti,'(i6)')itin
            !         open(fileunit,file='annit'//prti,status='unknown')

            !         write(stringprint)ti
            !         do k=kstawrite,kendwrite !0,jz+1
            !         do j=0,jy+1
            !         do i=0,jx+1
            !            write(stringprint)print_annit(i,j,k),print_annitV(i,j,k)
            !         end do
            !         end do
            !         end do
            !         close(fileunit)

            end if ! on the format
        end if ! if iscrivo


        if (i_paraview==1) then
            call write_paraview_files(tipo)
        end if


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
        !
        !-----------------------------------------------------------------------
        ! write the medietempo for time statistics

        ifmt_mean = 4 + nscal + 6 + nscal + nscal
        write(char_fmt,'(i2)')ifmt_mean
        fmt_mean = '('//char_fmt//'e18.10)'

        fileunit = 15

        open(fileunit,file=trim(filename_medietempo),status='unknown')
        write(fileunit,*)nscal,nproc,re
        write(fileunit,*)jx,jy,jz
        write(fileunit,*)deltatempo
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    !        write(fileunit,1001)umed(i,j,k,1),umed(i,j,k,2),
                    !     >         umed(i,j,k,3),umed(i,j,k,4),
                    !     >         umed(i,j,k,5), uutau(i,j,k,1),
                    !     >         uutau(i,j,k,2),uutau(i,j,k,3),
                    !     >         uutau(i,j,k,4),uutau(i,j,k,5),
                    !     >         uutau(i,j,k,6),uutau(i,j,k,7)


                    write(fileunit,fmt_mean)(umed(i,j,k,n),n=1,4+nscal), &       !u,v,w,fi,r_n
                        (uutau(i,j,k,n),n=1,6+nscal+nscal)  !uu_ii, uu_ij, rr_n, rv_n



                end do
            end do
        end do
        close(fileunit)

        !      end if



        !1001    format(12e18.10)

        !     the data are used only for average in the time window and not globally
        if (i_cumulative == 0) then
            deltatempo = 0.
            umed = 0.
            uutau = 0.
        end if

        !
        !-----------------------------------------------------------------------
        !
        if (myid == 0) write(*,*) 'OUTPUT STEP END ---------------------------'

    end subroutine write_output

    !***********************************************************************
    subroutine write_paraview_files(tipo)

        use myarrays_metri3, only: x,y,z
        use myarrays_ibm, only: bodyforce
        !use myarrays_levelset2, only: phi0, phi_temp, fext

        use mpi

        implicit none

        !
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! local variables declaration
        integer fileunit
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
            fileunit=17
            write(*,*) 'Write Paraview header file ',trim(filename_paraview)
            open(fileunit,file=trim(filename_paraview),status='unknown')
            write(fileunit,'(A)')'<?xml version="1.0"?>'
            write(fileunit,'(A)')'<VTKFile type="PImageData" version="0.1">'
            write(fileunit,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <PImageData WholeExtent="',&
                x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(0),' ',z_piece_end(nproc-1),&
                '" GhostLevel="1" Origin="0.0 0.0 0.0" Spacing="1.0 1.0 1.0">' ! spacing is missing
            write(fileunit,'(A)')'    <PPointData>'
!            write(fileunit,'(A)')'      <PDataArray type="Float64" Name="coord" NumberOfComponents="3">'
!            write(fileunit,'(A)')'      </PDataArray>'
            !            write(fileunit,'(A)')'      <PDataArray type="Int32" Name="y">'
            !            write(fileunit,'(A)')'      </PDataArray>'
            !            write(fileunit,'(A)')'      <PDataArray type="Int32" Name="z">'
            !            write(fileunit,'(A)')'      </PDataArray>'
            write(fileunit,'(A)')'      <PDataArray type="Float64" Name="vel" NumberOfComponents="3">'
            write(fileunit,'(A)')'      </PDataArray>'
            write(fileunit,'(A)')'      <PDataArray type="Float64" Name="pres">'
            write(fileunit,'(A)')'      </PDataArray>'
            !            write(fileunit,'(A)')'      <PDataArray type="Float64" Name="phi_temp">'
            !            write(fileunit,'(A)')'      </PDataArray>'
            !            write(fileunit,'(A)')'      <PDataArray type="Float64" Name="fext">'
            !            write(fileunit,'(A)')'      </PDataArray>'
            !            write(fileunit,'(A)')'      <PDataArray type="Float64" Name="phi0">'
            !            write(fileunit,'(A)')'      </PDataArray>'
            if (bodyforce == 1) then
                write(fileunit,'(A)')'      <PDataArray type="Int32" Name="tipo">'
                write(fileunit,'(A)')'      </PDataArray>'
            end if
            if (nscal>0) then
                do isc = 1,nscal
                    write(fileunit,'(A,I0.1,A)')'      <PDataArray type="Float64" Name="scal',isc,'">'
                    write(fileunit,'(A)')'      </PDataArray>'
                end do
            end if
            write(fileunit,'(A)')'    </PPointData>'
            write(fileunit,'(A)')'    <PCellData>'
            write(fileunit,'(A)')'    </PCellData>'
            !write(fileunit,'(A,I0.1,A,I0.1,A)')'    <Piece Extent="0 ',jx+1,' 0 ',jy+1,' ',piece_begin,' ',piece_end,'" Source="para_iter000010_0.vti"/>'
            do iproc=0,nproc-1
                write(fileunit,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,A,A)') '    <Piece Extent="',&
                    x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',&
                    z_piece_begin(iproc),' ',z_piece_end(iproc),&
                    '" Source="',trim(paraview_piecename(iproc)),'">'
                write(fileunit,'(A)')'    </Piece>'
            end do

            write(fileunit,'(A)')'  </PImageData>'
            write(fileunit,'(A)')'</VTKFile>'

            close(fileunit)
        end if

        ! every processor writes a piece file
        fileunit=18+myid
        write(*,*) 'Write Paraview piece file ',trim(filename_paraview_piece)
        open(fileunit,file=trim(filename_paraview_piece),status='unknown')
        write(fileunit,'(A)')'<?xml version="1.0"?>'
        write(fileunit,'(A)')'<VTKFile type="ImageData" version="0.1">'
        write(fileunit,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <ImageData WholeExtent="', &
            x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
            '">'
        write(fileunit,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <Piece Extent="',&
            x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
            '">'
        write(fileunit,'(A)')'    <PointData>'
!        write(fileunit,'(A)')'      <DataArray type="Float64" Name="coord" NumberOfComponents="3">'
!        do k=kparasta-1,kparaend+1
!            do j=0,jy+1
!                do i=0,jx+1
!                    write(fileunit,'(ES14.7,A,ES14.7,A,ES14.7,A)',advance='no') &
!                        x(i,j,k),' ',y(i,j,k),' ',z(i,j,k),' '
!                end do
!            end do
!        end do
!        write(fileunit,*)
!        write(fileunit,'(A)')'      </DataArray>'
        !            write(fileunit,'(A)')'      <DataArray type="Float64" Name="y">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(fileunit,'(I0.1,A)',advance='no') y(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(fileunit,*)
        !            write(fileunit,'(A)')'      </DataArray>'
        !            write(fileunit,'(A)')'      <DataArray type="Float64" Name="z">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(fileunit,'(I0.1,A)',advance='no') z(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(fileunit,*)
        !            write(fileunit,'(A)')'      </DataArray>'
        write(fileunit,'(A)')'      <DataArray type="Float64" Name="vel" NumberOfComponents="3">'
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    write(fileunit,'(ES14.7,A,ES14.7,A,ES14.7,A)',advance='no') &
                        print_u(i,j,k),' ',print_v(i,j,k),' ',print_w(i,j,k),' '
                end do
            end do
        end do
        write(fileunit,*)
        write(fileunit,'(A)')'      </DataArray>'
        write(fileunit,'(A)')'      <DataArray type="Float64" Name="pres">'
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    write(fileunit,'(ES14.7,A)',advance='no') print_fi(i,j,k),' '
                end do
            end do
        end do
        write(fileunit,*)
        write(fileunit,'(A)')'      </DataArray>'
        !            write(fileunit,'(A)')'      <DataArray type="Float64" Name="phi_temp">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(fileunit,'(ES14.7,A)',advance='no') phi_temp(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(fileunit,*)
        !            write(fileunit,'(A)')'      </DataArray>'
        !            write(fileunit,'(A)')'      <DataArray type="Float64" Name="Fext">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(fileunit,'(ES14.7,A)',advance='no') fext(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(fileunit,*)
        !            write(fileunit,'(A)')'      </DataArray>'
        !            write(fileunit,'(A)')'      <DataArray type="Float64" Name="phi0">'
        !            do k=kparasta-1,kparaend+1
        !                do j=0,jy+1
        !                    do i=0,jx+1
        !                        write(fileunit,'(ES14.7,A)',advance='no') phi0(i,j,k),' '
        !                    end do
        !                end do
        !            end do
        !            write(fileunit,*)
        !            write(fileunit,'(A)')'      </DataArray>'
        if (bodyforce == 1) then
            write(fileunit,'(A)')'      <DataArray type="Int32" Name="tipo">'
            do k=kparasta-1,kparaend+1
                do j=0,jy+1
                    do i=0,jx+1
                        write(fileunit,'(I0.1,A)',advance='no') tipo(i,j,k),' '
                    end do
                end do
            end do
            write(fileunit,'(A)')'      </DataArray>'
        end if
        if (nscal>0) then
            do isc = 1,nscal
                write(fileunit,'(A,I0.1,A)')'      <DataArray type="Float64" Name="scal',isc,'">'
                do k=kparasta-1,kparaend+1
                    do j=0,jy+1
                        do i=0,jx+1
                            write(fileunit,'(ES14.7,A)',advance='no') print_rhov(isc,i,j,k),' '
                        end do
                    end do
                end do
                write(fileunit,*)
            end do
        end if
        write(fileunit,'(A)')'      </DataArray>'
        write(fileunit,'(A)')'    </PointData>'
        write(fileunit,'(A)')'    <CellData>'
        write(fileunit,'(A)')'    </CellData>'
        write(fileunit,'(A)')'    </Piece>'
        write(fileunit,'(A)')'  </ImageData>'
        write(fileunit,'(A)')'</VTKFile>'

        close(fileunit)

    end subroutine write_paraview_files

end module output_module
