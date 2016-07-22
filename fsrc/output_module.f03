module output_module

    !     contains arrays and settings for printing procedure

    use iso_c_binding

    use mysending, only: myid,nproc,MPI_REAL_SD
    use scala3, only: n1,n2,n3,nscal,re,dt,kparaend,kparasta,deepl,deepr

    use mpi

    implicit none

    private


    ! public elements: -------------------------------------------------------
    integer(kind=c_int),bind(C),public :: i_print,ktime
    real,public :: ti

    ! flags
    integer(kind=c_int),bind(C),public :: print_iter_or_time,i_printfile,iformat_newres
    logical(kind=c_bool),bind(C),public :: i_newres,i_medietempo
    logical(kind=c_bool),bind(C),public :: paraview_aver,paraview_inst,paraview_compact
    logical(kind=c_bool),bind(C),public :: i_cumulative
    real(kind=c_double),bind(C),public :: i_time

    character(len=500),public :: result_folder
    character(len=500),public :: grid_file
    character(len=500),public :: restart_file

    ! formats
    character(len=500),public :: string_newres_format

    ! piani & sonde
    integer(kind=c_int),bind(C),public :: npiani,nsonde
    type(c_ptr),bind(C),public :: c_sondeindexi,c_sondeindexj,c_sondeindexk,c_piani

    ! scalar dimensionless numbers
    type(c_ptr),bind(C),public :: c_pran,c_prsc

    public :: create_output_folder,initialize_output,initialize_input
    public :: update_medie,output_step,output_finalize,get_cpp_strings,print_info

    ! private elements: -------------------------------------------------------
    logical :: iscrivo

    ! information output files
    integer,parameter,public :: info_run_file=11
    integer,parameter :: encheck_file=13
    integer,parameter :: new_res_file=14
    integer,parameter :: medietempo_file=15
    integer,parameter :: paraview_aver_file=17
    integer,parameter :: paraview_inst_file=18
    integer,parameter :: paraview_aver_piece_file_base=1400
    integer,parameter :: paraview_inst_piece_file_base=1700
    integer,parameter :: turbo_file=129
    integer,parameter :: subdis_file=130
    integer,parameter :: subscl_file=131
    integer,parameter :: subrho_file=133
    integer,parameter :: bulkvel_file=134
    integer,parameter :: wallstress_file=155
    integer,parameter :: particle_drag_file=154
    integer,parameter :: particle_vel_file=156
    integer,parameter :: maxvel_file=153
    integer,parameter :: sonde_file_base=2000

    ! formats
    character(len=20) :: fmt_newres

    ! file names
    character(len=500) :: filename_newres,filename_medietempo
    character(len=500) :: filename_paraview_aver,filename_paraview_inst
    character(len=500) :: filename_paraview_aver_piece,filename_paraview_inst_piece
    character(len=500) :: filename_inforun,filename_particle_drag,filename_particle_vel
    character(len=500) :: filename_bulkvel,filename_maxvel,filename_encheck,filename_wallstress

    character*120,allocatable :: paraview_inst_piecename(:)
    character*120,allocatable :: paraview_aver_piecename(:)

    ! folders to store results (see output_init)
    character(len=500) :: new_res_folder
    character(len=500) :: paraview_aver_folder
    character(len=500) :: paraview_inst_folder
    character(len=500) :: medietempo_folder

    character(len=4) :: char_myid

    real ti_start
    integer count_print_time

    ! piani & sonde
    integer,pointer :: sonde(:,:),piani(:)
    integer,pointer :: sondeindexi(:),sondeindexj(:),sondeindexk(:)

      
    ! variables for printing
    real, allocatable :: umed(:,:,:,:),uutau(:,:,:,:)
    real :: deltatempo

    real,allocatable :: print_u(:,:,:),print_v(:,:,:),print_w(:,:,:)
    real,allocatable :: print_uc(:,:,:),print_vc(:,:,:),print_wc(:,:,:)
    real,allocatable :: print_fi(:,:,:),print_rhov(:,:,:,:)

    real,allocatable :: rhocol(:),rhotot(:),ficol(:),fitot(:)
    real,allocatable :: ucol(:),utot(:),vcol(:),vtot(:),wcol(:),wtot(:)
    real,allocatable :: uccol(:),uctot(:),vccol(:),vctot(:),wccol(:),wctot(:)

contains

    subroutine create_output_folder()

        use wallmodel_module, only: wfp3,wfp4
        use mysettings, only: bodyforce

        character(len=12),parameter :: paraview_aver_folder_name='paraviewAverData'
        character(len=12),parameter :: paraview_inst_folder_name='paraviewInstData'
        character(len=10),parameter :: new_res_folder_name='newResData'
        character(len=14),parameter :: medietempo_folder_name='medietempoData'

        ! for sonde
        character(len=3) :: identificosonda
        character(len=20) :: filesonda

        integer :: isonde

        !-----------------------------------------------------------------------
        ! check date and time for the result folder creation
        !call date_and_time(date,time)

        if (myid==0) then
            write(*,*) 'Initialize output step'
        end if


        ! paraview folder, average
        if (paraview_aver) then
            paraview_aver_folder=trim(result_folder)//'/'//paraview_aver_folder_name
            if (myid==0) then
                call system('mkdir '//trim(paraview_aver_folder))
                write(*,*) 'Created folder ',trim(paraview_aver_folder)
            end if

        end if

        ! paraview folder, instantaneous
        if (paraview_inst) then
            paraview_inst_folder=trim(result_folder)//'/'//paraview_inst_folder_name
            if (myid==0) then
                call system('mkdir '//trim(paraview_inst_folder))
                write(*,*) 'Created folder ',trim(paraview_inst_folder)
            end if

        end if

        ! newRes folder
        if (i_newres) then
            new_res_folder=trim(result_folder)//'/'//new_res_folder_name
            if (myid==0) then
                call system('mkdir '//trim(new_res_folder))
                write(*,*) 'Created folder ',trim(new_res_folder)
            end if
        end if

        ! medietempo folder
        if (i_medietempo) then
            medietempo_folder=trim(result_folder)//'/'//medietempo_folder_name
            if (myid==0) then
                call system('mkdir '//trim(medietempo_folder))
                write(*,*) 'Created folder ',trim(medietempo_folder)
            end if

        end if

        ! Files for information output -------------------------------------
        filename_inforun=trim(result_folder)//'/'//'info_run.txt'
        open(info_run_file,file=filename_inforun,action='write',status='replace')

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     open output files

        if (myid==0) then

            filename_bulkvel=trim(result_folder)//'/'//'bulk_velocity.dat'
            open(bulkvel_file,file=filename_bulkvel,action='write',status='replace')
            filename_maxvel=trim(result_folder)//'/'//'max_velocity.dat'
            open(maxvel_file,file=filename_maxvel,action='write',status='replace')
            filename_encheck=trim(result_folder)//'/'//'encheck.dat'
            open(encheck_file,file=filename_encheck,action='write',status='replace')
            filename_particle_drag=trim(result_folder)//'/'//'particle_drag.dat'
            open(particle_drag_file,file=filename_particle_drag,action='write',status='replace')
            filename_particle_vel=trim(result_folder)//'/'//'particle_vel.dat'
            open(particle_vel_file,file=filename_particle_vel,action='write',status='replace')

            !  sonde files
            do isonde=1,nsonde
                write(identificosonda,'(i3.3)')isonde
                filesonda=trim(result_folder)//'nsonda'//identificosonda//'.dat'
                open(sonde_file_base+isonde,file=filesonda,action='write',status='unknown')
            end do

            if (wfp4.or.wfp3.or.bodyforce) then
                filename_wallstress=trim(result_folder)//'/'//'wall_stress.dat'
                open(wallstress_file,file=filename_wallstress,action='write',status='replace')
            end if

        end if


    end subroutine create_output_folder

    subroutine initialize_output()

        integer :: ichar

        ! check date and time for the result folder creation

        allocate(paraview_aver_piecename(0:nproc-1))
        allocate(paraview_inst_piecename(0:nproc-1))

        ! medietempo allocation
        if (i_medietempo .or. paraview_aver) then

            ! allocation
            allocate(umed(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,4+nscal))
            allocate(uutau(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,6+nscal+nscal))
            ! initialization
            deltatempo=0.0
            umed(:,:,:,:)=0.0
            uutau(:,:,:,:)=0.0

        end if

        ! string form of processor identifier
        write (char_myid,'(i4)') myid
        do ichar=1,len(char_myid)
            if (char_myid(ichar:ichar)==' ') char_myid(ichar:ichar)='0'
        end do
        !rid_myid=len_trim(char_myid)


        count_print_time=0
        ti_start=ti

    end subroutine initialize_output

    subroutine print_info()

        use mysending

        integer :: iproc

        !-----------------------------------------------------------------------
        ! print information about domain decomposition
        if (myid==0) then
            write(*,*)'----------------------------------------'
            write(info_run_file,*)'----------------------------------------'

            write(*,*)'number of procs: ',nproc
            write(info_run_file,*)'number of procs: ',nproc
        end if

        write(*,*)'I am proc',myid, 'of', nproc,'procs'
        write(info_run_file,*)'I am proc',myid, 'of', nproc,'procs'
        !
        !-----------------------------------------------------------------------
        ! MPI_REAL_SD assumes MPI_REAL8
        if (myid==0) then
            write(*,*)'DOUBLE PRECISION'
            write(info_run_file,*)'DOUBLE PRECISION'

            write(*,*)'n col. per PE: ',ncolperproc
            write(info_run_file,*)'n col. per PE: ',ncolperproc
        end if

        do iproc=0,nproc-1
            if (myid==iproc) then
                write(*,*)'PE: ',myid,'kparasta=', kparasta,' kparaend=',kparaend
                write(*,*)'right ',rightpe,'left ',leftpe
                write(info_run_file,*)'PE: ',myid,'kparasta=', kparasta,' kparaend=',kparaend
                write(info_run_file,*)'right ',rightpe,'left ',leftpe
            end if
        end do

    end subroutine print_info

    subroutine initialize_input()

        use mysettings, only: rich,pran,prsc

        integer :: i,ierr

        ! variable allocation for piani and sonde and convert (see above)
        allocate(piani(npiani))
        call C_F_POINTER(c_piani,piani,[npiani])
        allocate(sonde(3,nsonde))
        allocate(sondeindexi(nsonde),sondeindexj(nsonde),sondeindexk(nsonde))
        call C_F_POINTER(c_sondeindexi,sondeindexi,[nsonde])
        call C_F_POINTER(c_sondeindexj,sondeindexj,[nsonde])
        call C_F_POINTER(c_sondeindexk,sondeindexk,[nsonde])

        !-----------------------------------------------------------------------
        ! variable allocation for the scalar equations
        allocate(pran(nscal))
        allocate(prsc(nscal))

        ! convert values from c++format to fortran format
        call C_F_POINTER(c_pran,pran,[nscal])
        call C_F_POINTER(c_prsc,prsc,[nscal])

        if (myid==0) then
            write (*,*) "Reynolds=",re
            write (*,*) "Richardson=",rich
            do i=1,nscal
                write(*,*) "Scalar ",i," Prandtl=",pran(i),' Prscr=',prsc(i)
            end do

            write (*,*) "Total piani=",npiani
            do i=1,npiani
                write(*,*) "Piano number ",i," x=",piani(i)
            end do

            write (*,*) "Total sonde=",nsonde
            do i=1,nsonde
                write(*,*) "Sonda number ",i,"point:(",sondeindexi(i),",",sondeindexj(i),",",sondeindexk(i),")"
            end do
        end if

    end subroutine initialize_input

    subroutine output_step(tipo)

        ! write output
        use wallmodel_module, only: wfp3,wfp4
        use mysettings, only: bbx,particles,bodyforce
        !
        use mpi

        ! arguments
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        ! local variables
        integer :: ierr
        logical :: writeprocedure

        !-----------------------------------------------------------------------
        ! print restart file in the form new_res

        ! check if it's time to print an output file
        ! initial guess: no
        writeprocedure=.false.

        ! check condition on reached number of iterations
        if (print_iter_or_time==0 .and. ktime==i_print*(ktime/i_print)) then

            writeprocedure=.true.

        ! check condition on reached time
        else if (print_iter_or_time==1 .and. &
            ( ((ti-ti_start)-count_print_time*i_time)>=i_time .or. &
            abs(((ti-ti_start)-count_print_time*i_time)-i_time)<10e-7) ) then

            writeprocedure=.true.
            count_print_time=count_print_time+1

        end if


        ! if it's time to write output, start prodecure
        if (writeprocedure) then

            if (myid==0) write(*,*) 'OUTPUT STEP BEGINS ---------------------------------'

            ! set the folder for the output and the filename
            call write_preparation()

            ! write output new_res_form without mpi procedure
            if ((i_newres .and. i_printfile==0) .or. (i_newres .and. i_printfile==2)) then
                ! collect data to print
                call prepare_printdata
                ! and produce new_res
                call write_newres(bbx)
            end if

            ! print paraview
            if (paraview_aver .or. paraview_inst) then
                call write_paraview(tipo)
            end if

            ! print medietempo
            if (i_medietempo) then
                call write_medietempo()
            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (myid==0) write(*,*) 'OUTPUT STEP ENDS -----------------------------------'

            ! check if I can write another new_res before ending the number of iterations
        !            if (print_iter_or_time==1) then
        !
        !                ! time_evaluation=(niter-ktime)*dt
        !                ! time_evaluation=niter*dt - ti_start
        !
        !                time_evaluation=(ti-ti_start)+(niter-ktime)*dt
        !                next_print=time_evaluation/real(i_time)
        !                if (myid==0) write(*,*) 'time evaluation ',time_evaluation,next_print,count_print_time
        !
        !                if (next_print <=count_print_time .and. ktime/=niter) then
        !
        !                    if (myid==0) then
        !                        write(*,*)'remaining iteration not sufficent'
        !                        write(*,*)'to write another new_res'
        !                        write(*,*)'run is stopped'
        !                    end if
        !
        !                    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        !                    stop
        !                end if
        !            end if


        end if !write procedure

        ! print data from sonde
        call write_tracers()

        ! print kinetic energy
        call print_energy(tipo)

        ! print bulk velocity
        call print_velocity(tipo)

        ! print average drag force on particles
        if (bodyforce.and.particles) then
            call print_particle_drag()
            call print_particle_velocity()
        end if

        ! print bulk velocity
        if (wfp3.or.wfp4.or.bodyforce) then
            call print_wallstress()
        end if

        ! print priano for inflow
        call print_inflow_planes(ktime)


    end subroutine

    subroutine write_preparation()
        ! this subroutine prepares the file names and the folder

        !-----------------------------------------------------------------------
        ! file names
        character(len=7),parameter :: filename_header='new_res'
        character(len=9),parameter :: filemedie_header='meantime'
        character(len=8),parameter :: paraview_aver_header='paraAver'
        character(len=8),parameter :: paraview_inst_header='paraInst'

        integer :: ichar,iproc
        character(len=6)  :: char_ikiter
        character(len=13) :: char_iktime
        character(len=4)  :: char_piece
        !-------------------------------------------------------------------------

        ! set the folder for the output

        ! print for iterations
        if (print_iter_or_time==0) then

            write (char_ikiter,'(i6)') ktime
            do ichar=1,len(char_ikiter)
                if (char_ikiter(ichar:ichar)==' ') char_ikiter(ichar:ichar)='0'
            end do

        ! print for time
        else if (print_iter_or_time==1) then

            write(char_iktime,'(1e13.7)') ti ! INT(ti-ti_start)
            write(*,*) char_iktime

            do ichar=1,len(char_iktime)
                if (char_iktime(ichar:ichar)==' ') char_iktime(ichar:ichar)='0'
            end do


        end if

        ! flag for telling every processor if it must write
        iscrivo=.false.

        !-----------------------------------------------------------------------
        ! setting up names

        ! 1: set the base filename for new_res, paraview and medetempo
        filename_newres=trim(new_res_folder)//'/'//filename_header
        filename_medietempo=trim(medietempo_folder)//'/'//filemedie_header
        filename_paraview_aver=trim(paraview_aver_folder)//'/'//paraview_aver_header
        filename_paraview_inst=trim(paraview_inst_folder)//'/'//paraview_inst_header
        filename_paraview_aver_piece=trim(paraview_aver_folder)//'/'//paraview_aver_header
        filename_paraview_inst_piece=trim(paraview_inst_folder)//'/'//paraview_inst_header
        do iproc=0,nproc-1
            paraview_aver_piecename(iproc)=paraview_aver_header
            paraview_inst_piecename(iproc)=paraview_inst_header
        end do

        ! 2: add iteration or time ID
        if (print_iter_or_time==0) then
            ! add interation number
            filename_newres=trim(filename_newres)//'_iter'//char_ikiter
            filename_medietempo=trim(filename_medietempo)//'_iter'//char_ikiter
            filename_paraview_aver=trim(filename_paraview_aver)//'_iter'//char_ikiter
            filename_paraview_inst=trim(filename_paraview_inst)//'_iter'//char_ikiter
            filename_paraview_aver_piece=trim(filename_paraview_aver_piece)//'_iter'//char_ikiter
            filename_paraview_inst_piece=trim(filename_paraview_inst_piece)//'_iter'//char_ikiter
            do iproc=0,nproc-1
                paraview_aver_piecename(iproc)=trim(paraview_aver_piecename(iproc))//'_iter'//char_ikiter
                paraview_inst_piecename(iproc)=trim(paraview_inst_piecename(iproc))//'_iter'//char_ikiter
            end do
        else if (print_iter_or_time==1) then
            ! add time
            filename_newres=trim(filename_newres)//'_time'//char_iktime//'sec'
            filename_medietempo=trim(filename_medietempo)//'_time'//char_iktime//'sec'
            filename_paraview_aver=trim(filename_paraview_aver)//'_time'//char_iktime//'sec'
            filename_paraview_inst=trim(filename_paraview_inst)//'_time'//char_iktime//'sec'
            filename_paraview_aver_piece=trim(filename_paraview_aver_piece)//'_time'//char_iktime//'sec'
            filename_paraview_inst_piece=trim(filename_paraview_inst_piece)//'_time'//char_iktime//'sec'
            do iproc=0,nproc-1
                paraview_aver_piecename(iproc)=trim(paraview_aver_piecename(iproc))//'_time'//char_iktime//'sec'
                paraview_inst_piecename(iproc)=trim(paraview_inst_piecename(iproc))//'_time'//char_iktime//'sec'
            end do
        end if

        ! 3: add processor ID
        ! new_res can be parallel or not
        if (i_printfile==0 .or. i_printfile==1) then
            filename_newres=trim(filename_newres)//'.dat'
        else if (i_printfile==2) then
            ! add processor id
            filename_newres=trim(filename_newres)//'_proc'//char_myid//'.dat'
        end if
        ! medietempo is always parallel
        filename_medietempo=trim(filename_medietempo)//'_proc'//char_myid//'.dat'
        ! paraview is also always parallel (but with header)
        filename_paraview_aver_piece=trim(filename_paraview_aver_piece)//'_proc'//char_myid//'.vti'
        filename_paraview_inst_piece=trim(filename_paraview_inst_piece)//'_proc'//char_myid//'.vti'
        filename_paraview_inst=trim(filename_paraview_inst)//'.pvti'
        filename_paraview_aver=trim(filename_paraview_aver)//'.pvti'

        ! this is the LOCAL paraview piece name, necessary for
        do iproc=0,nproc-1
            write (char_piece,'(i4.1)') iproc
            do ichar=1,len(char_piece)
                if (char_piece(ichar:ichar)==' ') char_piece(ichar:ichar)='0'
            end do
            paraview_aver_piecename(iproc)=trim(paraview_aver_piecename(iproc))//'_proc'//char_piece//'.vti'
            paraview_inst_piecename(iproc)=trim(paraview_inst_piecename(iproc))//'_proc'//char_piece//'.vti'
        end do

        !-----------------------------------------------------------------------
        ! who writes the new_res??
        if (i_printfile==0 .or. i_printfile==1) then
            ! only proc 0 writes
            if (myid==0) iscrivo=.true.
        else if (i_printfile==2) then
            ! all procs write
            iscrivo=.true.
        end if

    end subroutine write_preparation

    subroutine write_paraview(tipo)

        use myarrays_metri3, only: centroid
        use myarrays_velo3, only: u,v,w,fi,rhov
        use particle_module, only: fluidPressureForce,fluidShearForce,forceCaso, &
            fluidMomentumForce,fluidParticleVel,print_fields
        use mysettings, only: attiva_scal,bodyforce,particles
        use turbo_module, only: q_crit
        use mysending

        !-----------------------------------------------------------------------
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        integer :: piecefile
        integer :: i,j,k,isc,iproc
        integer :: z_piece_begin(0:nproc-1),z_piece_end(0:nproc-1)
        integer :: x_piece_begin,x_piece_end
        integer :: y_piece_begin,y_piece_end
        integer :: ierr
        !-----------------------------------------------------------------------


        ! we need kparasta and kparaend saved at the root processor
        ! for this, a gather is needed
        call MPI_GATHER(kparasta-1,1,MPI_INTEGER,z_piece_begin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(kparaend+1,1,MPI_INTEGER,z_piece_end,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        x_piece_begin=0
        x_piece_end=n1+1
        y_piece_begin=0
        y_piece_end=n2+1
        z_piece_begin(myid)=kparasta-1
        z_piece_end(myid)=kparaend+1

        ! header file is written only by root processor

        if (paraview_inst) then

            ! PART 1: WRITE TOP PART OF FILES -------------

            if (myid==0) then

                write(*,*) 'Write instantaneous Paraview header file ',trim(filename_paraview_inst)
                open(paraview_inst_file,file=trim(filename_paraview_inst),status='unknown')
                write(paraview_inst_file,'(A)')'<?xml version="1.0"?>'
                write(paraview_inst_file,'(A)')'<VTKFile type="PImageData" version="0.1">'
                write(paraview_inst_file,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <PImageData WholeExtent="',&
                    x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(0),' ',z_piece_end(nproc-1),&
                    '" GhostLevel="1" Origin="0.0 0.0 0.0" Spacing="1.0 1.0 1.0">' ! spacing is missing
                write(paraview_inst_file,'(A)')'    <PPointData>'
                write(paraview_inst_file,'(A)')'      <PDataArray type="Float64" Name="vel" NumberOfComponents="3"/>'
                write(paraview_inst_file,'(A)')'      <PDataArray type="Float64" Name="pres"/>'
                if (.not.paraview_compact) then
                    write(paraview_inst_file,'(A)')'      <PDataArray type="Float64" Name="warp" NumberOfComponents="3"/>'

                    write(paraview_inst_file,'(A)')'      <PDataArray type="Float64" Name="q_crit"/>'
                    if (bodyforce) then
                        write(paraview_inst_file,'(A)')'      <PDataArray type="Int32" Name="tipo"/>'
                        if (particles .and. print_fields) then
                            !write(paraview_inst_file,'(A)')'      <PDataArray type="Float64" Name="normal" NumberOfComponents="3"/>'
                            write(paraview_inst_file,'(A)')'      <PDataArray type="Int32" Name="caso"/>'
                            write(paraview_inst_file,'(A)') &
                            '      <PDataArray type="Float64" Name="fluidPressureForce" NumberOfComponents="3"/>'
                            write(paraview_inst_file,'(A)') &
                            '      <PDataArray type="Float64" Name="fluidShearForce" NumberOfComponents="3"/>'
                            write(paraview_inst_file,'(A)') &
                            '      <PDataArray type="Float64" Name="particleIndex"/>'
                        !                    write(paraview_inst_file,'(A)')'      <PDataArray type="Float64" Name="fluidMomentumForce" NumberOfComponents="3"/>'
                        !                    write(paraview_inst_file,'(A)')'      <PDataArray type="Float64" Name="fluidParticleVel" NumberOfComponents="3"/>'
                        end if
                    end if
                    if (attiva_scal) then
                        do isc=1,nscal
                            write(paraview_inst_file,'(A,I0.1,A)')'      <PDataArray type="Float64" Name="scal',isc,'"/>'
                        end do
                    end if
                end if
                write(paraview_inst_file,'(A)')'    </PPointData>'
                write(paraview_inst_file,'(A)')'    <PCellData>'
                write(paraview_inst_file,'(A)')'    </PCellData>'
                do iproc=0,nproc-1
                    write(paraview_inst_file,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,A,A)') '    <Piece Extent="',&
                        x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',&
                        z_piece_begin(iproc),' ',z_piece_end(iproc),&
                        '" Source="',trim(paraview_inst_piecename(iproc)),'">'
                    write(paraview_inst_file,'(A)')'    </Piece>'
                end do

                write(paraview_inst_file,'(A)')'  </PImageData>'
                write(paraview_inst_file,'(A)')'</VTKFile>'

                close(paraview_inst_file)

            end if

            ! PART 2: ADD VARIABLES -------------

            ! q criterion nees completion of ghost cells
            call var_complete_exchange(q_crit)

            ! every processor writes a piece file
            piecefile=paraview_inst_piece_file_base+myid
            write(*,*) 'Write instantaneous Paraview piece file ',trim(paraview_inst_piecename(myid))
            open(piecefile,file=trim(filename_paraview_inst_piece),status='unknown')
            write(piecefile,'(A)')'<?xml version="1.0"?>'
            write(piecefile,'(A)')'<VTKFile type="ImageData" version="0.1">'
            write(piecefile,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <ImageData WholeExtent="', &
                x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
                '">'
            write(piecefile,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <Piece Extent="',&
                x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
                '">'
            write(piecefile,'(A)')'    <PointData>'
            write(piecefile,'(A)')'      <DataArray type="Float64" Name="vel" NumberOfComponents="3">'
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                            u(i,j,k),' ',v(i,j,k),' ',w(i,j,k),' '
                    end do
                end do
            end do
            write(piecefile,*)
            write(piecefile,'(A)')'      </DataArray>'
            write(piecefile,'(A)')'      <DataArray type="Float64" Name="pres">'
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        write(piecefile,'(ES14.6E3,A)',advance='no') fi(i,j,k),' '
                    end do
                end do
            end do
            write(piecefile,*)
            write(piecefile,'(A)')'      </DataArray>'
            if (.not.paraview_compact) then
                write(piecefile,'(A)')'      <DataArray type="Float64" Name="warp" NumberOfComponents="3">'
                do k=kparasta-1,kparaend+1
                    do j=0,n2+1
                        do i=0,n1+1
                            write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                                centroid(1,i,j,k)-1.0*real(i),' ', &
                                centroid(2,i,j,k)-1.0*real(j),' ', &
                                centroid(3,i,j,k)-1.0*real(k),' '
                        end do
                    end do
                end do
                write(piecefile,*)
                write(piecefile,'(A)')'      </DataArray>'
                write(piecefile,'(A)')'      <DataArray type="Float64" Name="q_crit">'
                do k=kparasta-1,kparaend+1
                    do j=0,n2+1
                        do i=0,n1+1
                            write(piecefile,'(ES14.6E3,A)',advance='no') q_crit(i,j,k),' '
                        end do
                    end do
                end do
                write(piecefile,*)
                write(piecefile,'(A)')'      </DataArray>'
                if (bodyforce) then
                    write(piecefile,'(A)')'      <DataArray type="Int32" Name="tipo">'
                    do k=kparasta-1,kparaend+1
                        do j=0,n2+1
                            do i=0,n1+1
                                write(piecefile,'(I0.1,A)',advance='no') tipo(i,j,k),' '
                            end do
                        end do
                    end do
                    write(piecefile,*)
                    write(piecefile,'(A)')'      </DataArray>'
                    if (particles .and. print_fields) then
                                        !            write(piecefile,'(A)')'      <DataArray type="Float64" Name="normal" NumberOfComponents="3">'
                    !            do k=kparasta-1,kparaend+1
                    !                do j=0,n2+1
                    !                    do i=0,n1+1
                    !                        write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                    !                            normalVectorX(i,j,k),' ',normalVectorY(i,j,k),' ',normalVectorZ(i,j,k),' '
                    !                    end do
                    !                end do
                    !            end do
                    !            write(piecefile,*)
                    !            write(piecefile,'(A)')'      </DataArray>'
                    write(piecefile,'(A)')'      <DataArray type="Int32" Name="caso">'
                    do k=kparasta-1,kparaend+1
                        do j=0,n2+1
                            do i=0,n1+1
                                write(piecefile,'(I0.1,A)',advance='no') forceCaso(i,j,k),' '
                            end do
                        end do
                    end do
                    write(piecefile,*)
                    write(piecefile,'(A)')'      </DataArray>'
                        write(piecefile,'(A)')'      <DataArray type="Float64" Name="fluidPressureForce" NumberOfComponents="3">'
                        do k=kparasta-1,kparaend+1
                            do j=0,n2+1
                                do i=0,n1+1
                                    write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                                        fluidPressureForce(1,i,j,k),' ', &
                                        fluidPressureForce(2,i,j,k),' ', &
                                        fluidPressureForce(3,i,j,k),' '
                                end do
                            end do
                        end do
                        write(piecefile,*)
                        write(piecefile,'(A)')'      </DataArray>'
                        write(piecefile,'(A)')'      <DataArray type="Float64" Name="fluidShearForce" NumberOfComponents="3">'
                        do k=kparasta-1,kparaend+1
                            do j=0,n2+1
                                do i=0,n1+1
                                    write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                                        fluidShearForce(1,i,j,k),' ', &
                                        fluidShearForce(2,i,j,k),' ', &
                                        fluidShearForce(3,i,j,k),' '
                                end do
                            end do
                        end do
                        write(piecefile,*)
                        write(piecefile,'(A)')'      </DataArray>'
                    !                write(piecefile,'(A)')'      <DataArray type="Float64" Name="fluidMomentumForce" NumberOfComponents="3">'
                    !                do k=kparasta-1,kparaend+1
                    !                    do j=0,n2+1
                    !                        do i=0,n1+1
                    !                            write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                    !                                fluidMomentumForce(1,i,j,k),' ',fluidMomentumForce(2,i,j,k),' ',fluidMomentumForce(3,i,j,k),' '
                    !                        end do
                    !                    end do
                    !                end do
                    !                write(piecefile,*)
                    !                write(piecefile,'(A)')'      </DataArray>'
                    !                write(piecefile,'(A)')'      <DataArray type="Float64" Name="fluidParticleVel" NumberOfComponents="3">'
                    !                do k=kparasta-1,kparaend+1
                    !                    do j=0,n2+1
                    !                        do i=0,n1+1
                    !                            write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                    !                                fluidParticleVel(1,i,j,k),' ',fluidParticleVel(2,i,j,k),' ',fluidParticleVel(3,i,j,k),' '
                    !                        end do
                    !                    end do
                    !                end do
                    !                write(piecefile,*)
                    !                write(piecefile,'(A)')'      </DataArray>'
                    end if
                end if
                if (attiva_scal) then
                    do isc=1,nscal
                        write(piecefile,'(A,I0.1,A)')'      <DataArray type="Float64" Name="scal',isc,'">'
                        do k=kparasta-1,kparaend+1
                            do j=0,n2+1
                                do i=0,n1+1
                                    write(piecefile,'(ES14.6E3,A)',advance='no') rhov(isc,i,j,k),' '
                                end do
                            end do
                        end do
                        write(piecefile,*)
                        write(piecefile,'(A)')'      </DataArray>'
                    end do
                end if
            end if
            write(piecefile,'(A)')'    </PointData>'
            write(piecefile,'(A)')'    <CellData>'
            write(piecefile,'(A)')'    </CellData>'
            write(piecefile,'(A)')'    </Piece>'
            write(piecefile,'(A)')'  </ImageData>'
            write(piecefile,'(A)')'</VTKFile>'

            close(piecefile)

        end if

        if (paraview_aver) then

            ! PART 1: WRITE TOP PART OF FILES -------------

            if (myid==0) then

                write(*,*) 'Write averaged Paraview header file ',trim(filename_paraview_aver)
                open(paraview_aver_file,file=trim(filename_paraview_aver),status='unknown')
                write(paraview_aver_file,'(A)')'<?xml version="1.0"?>'
                write(paraview_aver_file,'(A)')'<VTKFile type="PImageData" version="0.1">'
                write(paraview_aver_file,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <PImageData WholeExtent="',&
                    x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(0),' ',z_piece_end(nproc-1),&
                    '" GhostLevel="1" Origin="0.0 0.0 0.0" Spacing="1.0 1.0 1.0">' ! spacing is missing
                write(paraview_aver_file,'(A)')'    <PPointData>'
                write(paraview_aver_file,'(A)')'      <PDataArray type="Float64" Name="av_vel" NumberOfComponents="3"/>'
                write(paraview_aver_file,'(A)')'      <PDataArray type="Float64" Name="av_pres"/>'
                write(paraview_aver_file,'(A)')'      <PDataArray type="Float64" Name="eddy_vel" NumberOfComponents="6"/>'
                if (.not.paraview_compact) then
                    write(paraview_aver_file,'(A)')'      <PDataArray type="Float64" Name="warp" NumberOfComponents="3"/>'
                    if (attiva_scal) then
                        do isc=1,nscal
                            write(paraview_aver_file,'(A,I0.1,A)')'      <PDataArray type="Float64" Name="scal',isc,'"/>'
                        end do
                    end if
                end if
                write(paraview_aver_file,'(A)')'    </PPointData>'
                write(paraview_aver_file,'(A)')'    <PCellData>'
                write(paraview_aver_file,'(A)')'    </PCellData>'
                do iproc=0,nproc-1
                    write(paraview_aver_file,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,A,A)') '    <Piece Extent="',&
                        x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',&
                        z_piece_begin(iproc),' ',z_piece_end(iproc),&
                        '" Source="',trim(paraview_aver_piecename(iproc)),'">'
                    write(paraview_aver_file,'(A)')'    </Piece>'
                end do

                write(paraview_aver_file,'(A)')'  </PImageData>'
                write(paraview_aver_file,'(A)')'</VTKFile>'

                close(paraview_aver_file)

            end if

            ! PART 2: ADD VARIABLES -------------

            ! every processor writes a piece file
            piecefile=paraview_aver_piece_file_base+myid
            write(*,*) 'Write averaged Paraview piece file ',trim(paraview_aver_piecename(myid))
            open(piecefile,file=trim(filename_paraview_aver_piece),status='unknown')
            write(piecefile,'(A)')'<?xml version="1.0"?>'
            write(piecefile,'(A)')'<VTKFile type="ImageData" version="0.1">'
            write(piecefile,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <ImageData WholeExtent="', &
                x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
                '">'
            write(piecefile,'(A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A,I0.1,A)')'  <Piece Extent="',&
                x_piece_begin,' ',x_piece_end,' ',y_piece_begin,' ',y_piece_end,' ',z_piece_begin(myid),' ',z_piece_end(myid),&
                '">'
            write(piecefile,'(A)')'    <PointData>'
            write(piecefile,'(A)')'      <DataArray type="Float64" Name="av_vel" NumberOfComponents="3">'
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                            umed(i,j,k,1)/deltatempo,' ',umed(i,j,k,2)/deltatempo,' ',umed(i,j,k,3)/deltatempo,' '
                    end do
                end do
            end do
            write(piecefile,*)
            write(piecefile,'(A)')'      </DataArray>'

            write(piecefile,'(A)')'      <DataArray type="Float64" Name="av_pres">'
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        write(piecefile,'(ES14.6E3,A)',advance='no') umed(i,j,k,4)/deltatempo,' '
                    end do
                end do
            end do
            write(piecefile,*)
            write(piecefile,'(A)')'      </DataArray>'
            write(piecefile,'(A)')'      <DataArray type="Float64" Name="eddy_vel" NumberOfComponents="6">'
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A,ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                            sqrt(abs( uutau(i,j,k,1)/deltatempo-umed(i,j,k,1)**2.0/deltatempo**2.0 )),' ', & !urms
                            sqrt(abs( uutau(i,j,k,4)/deltatempo-umed(i,j,k,2)**2.0/deltatempo**2.0 )),' ', & !vrms
                            sqrt(abs( uutau(i,j,k,6)/deltatempo-umed(i,j,k,3)**2.0/deltatempo**2.0 )),' ', & !wrms
                            uutau(i,j,k,2)/deltatempo-umed(i,j,k,1)*umed(i,j,k,2)/deltatempo**2.0,' ', & !uv
                            uutau(i,j,k,3)/deltatempo-umed(i,j,k,1)*umed(i,j,k,3)/deltatempo**2.0,' ', & !uw
                            uutau(i,j,k,5)/deltatempo-umed(i,j,k,2)*umed(i,j,k,3)/deltatempo**2.0,' ' !vw
                    end do
                end do
            end do
            write(piecefile,*)
            write(piecefile,'(A)')'      </DataArray>'
            if (.not.paraview_compact) then

                write(piecefile,'(A)')'      <DataArray type="Float64" Name="warp" NumberOfComponents="3">'
                do k=kparasta-1,kparaend+1
                    do j=0,n2+1
                        do i=0,n1+1
                            write(piecefile,'(ES14.6E3,A,ES14.6E3,A,ES14.6E3,A)',advance='no') &
                                centroid(1,i,j,k)-1.0*real(i),' ',&
                                centroid(2,i,j,k)-1.0*real(j),' ',&
                                centroid(3,i,j,k)-1.0*real(k),' '
                        end do
                    end do
                end do
                write(piecefile,*)
                write(piecefile,'(A)')'      </DataArray>'
                if (attiva_scal) then
                    do isc=1,nscal
                        write(piecefile,'(A,I0.1,A)')'      <DataArray type="Float64" Name="scal',isc,'">'
                        do k=kparasta-1,kparaend+1
                            do j=0,n2+1
                                do i=0,n1+1
                                    write(piecefile,'(ES14.6E3,A)',advance='no') rhov(isc,i,j,k),' '
                                end do
                            end do
                        end do
                        write(piecefile,*)
                        write(piecefile,'(A)')'      </DataArray>'
                    end do
                end if
            end if
            write(piecefile,'(A)')'    </PointData>'
            write(piecefile,'(A)')'    <CellData>'
            write(piecefile,'(A)')'    </CellData>'
            write(piecefile,'(A)')'    </Piece>'
            write(piecefile,'(A)')'  </ImageData>'
            write(piecefile,'(A)')'</VTKFile>'

            close(piecefile)

        end if


    end subroutine write_paraview

    subroutine write_newres(bbx)

        ! this sub write the new_res_form which contains the flow field for restart the computation or for analysis.
        ! depending on i_printfile only proc0 or all procs write their file depending on iformat_newres the output is written as
        !  iformat_newres=0 --->*
        !  iformat_newres=1 ---> use a format from Agenerale.in
        !  iformat_newres=2 ---> binary

        ! arguments
        real,intent(in) :: bbx

        !-----------------------------------------------------------------------
        ! local variables declaration
        integer :: kstawrite,kendwrite
        integer :: i,j,k,isc
        real,parameter :: dumpvalue=0.0

        !-----------------------------------------------------------------------
        !     start to write new_res file (only proc that have iscrivo true write)
        if (iscrivo) then

            write(*,*) 'Write newres on file: ',filename_newres

            ! with no format
            if (iformat_newres==0) then

                open(new_res_file,file=filename_newres,status='unknown')
                write (new_res_file,*) nscal
                write (new_res_file,*) ti
                write (new_res_file,*) bbx,dumpvalue

                kstawrite=0
                kendwrite=0

                !  write the cartesian component
                if (i_printfile==0) then
                    kstawrite=0
                    kendwrite=n3+1
                else if (i_printfile==2) then
                    kstawrite=kparasta
                    kendwrite=kparaend
                    if (myid==0) kstawrite=kparasta -1
                    if (myid==nproc-1) kendwrite=kparaend+1
                end if

                do k=kstawrite,kendwrite         !0,n3+1
                    do j=0,n2+1
                        do i=0,n1+1
                            write (new_res_file,*) print_u(i,j,k)
                            write (new_res_file,*) print_v(i,j,k)
                            write (new_res_file,*) print_w(i,j,k)
                            write (new_res_file,*) print_fi(i,j,k)
                            do isc=1,nscal
                                write (new_res_file,*) print_rhov(isc,i,j,k)
                            end do
                        end do
                    end do
                end do

                ! write the controvariant component
                if (i_printfile==0) then
                    kstawrite=1
                    kendwrite=n3
                end if
                if (i_printfile==2) then
                    kstawrite=kparasta
                    kendwrite=kparaend
                end if


                do k=kstawrite,kendwrite !1,n3
                    do j=1,n2
                        do i=0,n1
                            write(new_res_file,*)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,n3
                    do j=0,n2
                        do i=1,n1
                            write(new_res_file,*)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid==0) kstawrite=0
                do k=kstawrite,kendwrite  !0,n3
                    do j=1,n2
                        do i=1,n1
                            write(new_res_file,*)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(new_res_file)


            !-----------------------------------------------------------------------
            ! with format fmt_newres
            else if (iformat_newres==1) then
                open(new_res_file,file=trim(filename_newres),status='unknown')
                write(new_res_file,*) nscal
                write(new_res_file,*) ti
                write(new_res_file,*) bbx,dumpvalue

                kstawrite=0
                kendwrite=0

                ! write the cartesiano component
                if (i_printfile==0) then
                    kstawrite=0
                    kendwrite=n3+1
                end if
                if (i_printfile==2) then
                    kstawrite=kparasta
                    kendwrite=kparaend
                    if (myid==0)kstawrite=kparasta -1
                    if (myid==nproc-1)kendwrite=kparaend+1
                end if

                do k=kstawrite,kendwrite         !0,n3+1
                    do j=0,n2+1
                        do i=0,n1+1
                            write (new_res_file,fmt_newres) print_u(i,j,k)
                            write (new_res_file,fmt_newres) print_v(i,j,k)
                            write (new_res_file,fmt_newres) print_w(i,j,k)
                            write (new_res_file,fmt_newres) print_fi(i,j,k)
                            do isc=1,nscal
                                write (new_res_file,fmt_newres) print_rhov(isc,i,j,k)
                            end do
                        end do
                    end do
                end do

                !        write the controvariant component
                if (i_printfile==0) then
                    kstawrite=1
                    kendwrite=n3
                end if
                if (i_printfile==2) then
                    kstawrite=kparasta
                    kendwrite=kparaend
                end if


                do k=kstawrite,kendwrite !1,n3
                    do j=1,n2
                        do i=0,n1
                            write(new_res_file,fmt_newres)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,n3
                    do j=0,n2
                        do i=1,n1
                            write(new_res_file,fmt_newres)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid==0) kstawrite=0
                do k=kstawrite,kendwrite  !0,n3
                    do j=1,n2
                        do i=1,n1
                            write(new_res_file,fmt_newres)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(new_res_file)


            !-----------------------------------------------------------------------
            !     binary
            else if (iformat_newres==2) then

                open(new_res_file,file=filename_newres,status='unknown',form='unformatted')

                write(new_res_file)nscal
                write(new_res_file)ti
                write(new_res_file)bbx,dumpvalue

                kstawrite=0
                kendwrite=0

                ! write the cartesiano component
                if (i_printfile==0) then
                    kstawrite=0
                    kendwrite=n3+1
                end if
                if (i_printfile==2) then
                    kstawrite=kparasta
                    kendwrite=kparaend
                    if (myid==0)kstawrite=kparasta -1
                    if (myid==nproc-1)kendwrite=kparaend+1
                end if

                do k=kstawrite,kendwrite  !0,n3+1
                    do j=0,n2+1
                        do i=0,n1+1
                            write(new_res_file)print_u(i,j,k), &
                                print_v(i,j,k), &
                                print_w(i,j,k), &
                                print_fi(i,j,k), &
                                (print_rhov(isc,i,j,k),isc=1,nscal)
                        end do
                    end do
                end do

                !        write the controvariant component
                if (i_printfile==0) then
                    kstawrite=1
                    kendwrite=n3
                end if
                if (i_printfile==2) then
                    kstawrite=kparasta
                    kendwrite=kparaend
                end if


                do k=kstawrite,kendwrite !1,n3
                    do j=1,n2
                        do i=0,n1
                            write(new_res_file)print_uc(i,j,k)
                        end do
                    end do
                end do

                do k=kstawrite,kendwrite  !1,n3
                    do j=0,n2
                        do i=1,n1
                            write(new_res_file)print_vc(i,j,k)
                        end do
                    end do
                end do

                if (myid==0) kstawrite=0
                do k=kstawrite,kendwrite  !0,n3
                    do j=1,n2
                        do i=1,n1
                            write(new_res_file)print_wc(i,j,k)
                        end do
                    end do
                end do

                close(new_res_file)


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
        deallocate(print_rhov)

    end subroutine write_newres

    subroutine write_medietempo()

        ! write the medietempo for time statistics

        !-----------------------------------------------------------------------
        integer :: ifmt_mean
        character(len=2) :: char_fmt
        character(len=12) :: fmt_mean
        integer :: i,j,k,n
        !-----------------------------------------------------------------------

        ifmt_mean=4+nscal+6+nscal+nscal
        write(char_fmt,'(i2)')ifmt_mean
        fmt_mean='('//char_fmt//'e18.10)'

        open(medietempo_file,file=trim(filename_medietempo),status='new',action='write')
        write(*,*) 'Write medietempo file: ',trim(filename_medietempo)
        write(medietempo_file,*)nscal,nproc,re
        write(medietempo_file,*)n1,n2,n3
        write(medietempo_file,*)deltatempo

        do k=kparasta-deepl,kparaend+deepr
            do j=0,n2+1
                do i=0,n1+1

                    write(medietempo_file,fmt_mean)(umed(i,j,k,n),n=1,4+nscal), &       !u,v,w,fi,r_n
                        (uutau(i,j,k,n),n=1,6+nscal+nscal)  !uu_ii, uu_ij, rr_n, rv_n

                end do
            end do
        end do

        close(medietempo_file)

        !     the data are used only for average in the time window and not globally
        if (.not.i_cumulative) then
            deltatempo=0.0
            umed(:,:,:,:)=0.0
            uutau(:,:,:,:)=0.0
        end if

    end subroutine write_medietempo

    subroutine write_tracers()

        !-----------------------------------------------------------------------
        integer :: isc,isonde
        !-----------------------------------------------------------------------

        !     print tracers
        do isonde=1,nsonde
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

    end subroutine write_tracers

    subroutine output_finalize()


        integer :: isonde

        if (myid==0) then

            close(bulkvel_file)   ! bulk velocity nel tempo
            close(maxvel_file)
            close(encheck_file)
            !      close(35)   ! eddy viscosity
            !close(50)

            do isonde=1,nsonde
                close(sonde_file_base+isonde)
            end do

        end if

    end subroutine output_finalize

    subroutine print_energy(tipo)

        use myarrays_metri3, only: fluid_vol,giac
        use myarrays_velo3, only: u,v,w

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
            do j=1,n2
                do i=1,n1
                    if (tipo(i,j,k)/=0) then
                        enck_loc=enck_loc+.5*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))*giac(i,j,k)
                    end if
                end do
            end do
        end do

        ! reduce operation on enck and vol
        enck=0.
        call MPI_REDUCE(enck_loc,enck,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (myid==0) then
            enck=enck/fluid_vol
            open(encheck_file,file=filename_encheck,status='OLD',action='WRITE',position='APPEND')
            write(encheck_file,*) ti,enck
            close(encheck_file)
        end if
    end subroutine print_energy

    subroutine print_velocity(tipo)

        use myarrays_metri3, only: fluid_vol,giac
        use myarrays_velo3, only: u,v,w,fi

        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        real :: u_bulk,u_bulk_tot,v_bulk,v_bulk_tot,w_bulk,w_bulk_tot
        real :: u_max,u_max_tot,v_max,v_max_tot,w_max,w_max_tot
        real :: fi_bulk,fi_bulk_tot
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

        ! bulk computational pressure
        fi_bulk=0.
        fi_bulk_tot=0.


        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    if (tipo(i,j,k)/=0) then
                        u_bulk=u_bulk+u(i,j,k)*giac(i,j,k)
                        v_bulk=v_bulk+v(i,j,k)*giac(i,j,k)
                        w_bulk=w_bulk+w(i,j,k)*giac(i,j,k)
                        fi_bulk=fi_bulk+fi(i,j,k)*giac(i,j,k)

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
        call MPI_REDUCE(fi_bulk,fi_bulk_tot,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! u_max max between procs
        call MPI_REDUCE(u_max,u_max_tot,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(v_max,v_max_tot,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(w_max,w_max_tot,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (myid==0) then

            u_bulk_tot=u_bulk_tot/fluid_vol
            v_bulk_tot=v_bulk_tot/fluid_vol
            w_bulk_tot=w_bulk_tot/fluid_vol
            fi_bulk_tot=fi_bulk_tot/fluid_vol

            write(*,*)'u_bulk:',u_bulk_tot,'fi_bulk:',fi_bulk_tot,', u_max:',u_max_tot

            open(bulkvel_file,file=filename_bulkvel,status='old',action='write',position='append')
            write(bulkvel_file,*) ti,u_bulk_tot,v_bulk_tot,w_bulk_tot
            close(bulkvel_file)

            open(maxvel_file,file=filename_maxvel,status='old',action='write',position='append')
            write(maxvel_file,*) ti,u_max_tot,v_max_tot,w_max_tot
            close(maxvel_file)

        end if

    end subroutine print_velocity

    subroutine print_particle_drag()

        use particle_module, only: sphereShearForce,spherePressureForce,sphereMomentumForce,num_part_tot

        real,dimension(3) :: averageShearForce,averagePressureForce,averageMomentumForce
        integer :: p

        ! compute average drag
        averageShearForce(:)=0.0
        averagePressureForce(:)=0.0
        averageMomentumForce(:)=0.0


        do p=1,num_part_tot
            averageShearForce(:)=averageShearForce(:)+sphereShearForce(:,p)/real(num_part_tot)
            averagePressureForce(:)=averagePressureForce(:)+spherePressureForce(:,p)/real(num_part_tot)
            averageMomentumForce(:)=averageMomentumForce(:)+sphereMomentumForce(:,p)/real(num_part_tot)
        end do

        if (myid==0) then

            open(particle_drag_file,file=filename_particle_drag,status='old',action='write',position='append')
            write(*,*) 'Forces written on drag file :'
            write(*,*) ti,averageShearForce(1),averageShearForce(2),averageShearForce(3), &
                averagePressureForce(1),averagePressureForce(2),averagePressureForce(3), &
                averageMomentumForce(1),averageMomentumForce(2),averageMomentumForce(3)
            write(particle_drag_file,*) ti,averageShearForce(1),averageShearForce(2),averageShearForce(3), &
                averagePressureForce(1),averagePressureForce(2),averagePressureForce(3), &
                averageMomentumForce(1),averageMomentumForce(2),averageMomentumForce(3)
            close(particle_drag_file)

        end if

    end subroutine print_particle_drag

    subroutine print_particle_velocity()

        use particle_module, only: sphereVelocity,sphereMoves,spherePosition,num_part_tot,border_left,border_right
        use ibm_module, only: update_ibm

        ! ---------------------------------------------------------------
        real,dimension(3) :: average_velocity_loc,average_velocity_tot
        integer :: p,counter
        integer :: ierr
        ! ---------------------------------------------------------------

        if (update_ibm) then
            ! compute average drag
            average_velocity_loc(:)=0.0
            average_velocity_tot(:)=0.0
            counter=0

            do p=1,num_part_tot
                if (sphereMoves(p)) then
                    if (spherePosition(3,p)>border_left .and. spherePosition(3,p)<border_right) then
                        counter=counter+1
                        average_velocity_loc(:)=average_velocity_loc(:)+sphereVelocity(:,p)
                    end if
                end if
            end do

            average_velocity_loc(:)=average_velocity_loc(:)/real(counter)

            call MPI_REDUCE(average_velocity_loc,average_velocity_tot,3,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (myid==0) then

                average_velocity_tot=average_velocity_tot/real(nproc)

                open(particle_vel_file,file=filename_particle_vel,status='old',action='write',position='append')
                write(*,*) 'Average particle velocity:'
                write(*,*) ti,average_velocity_tot(1),average_velocity_tot(2),average_velocity_tot(3)
                write(particle_vel_file,*) ti,average_velocity_tot(1),average_velocity_tot(2),average_velocity_tot(3)
                close(particle_vel_file)

            end if

        end if

    end subroutine print_particle_velocity

    subroutine update_medie()

        use myarrays_velo3, only: u,v,w,fi,rhov

        integer :: i,j,k,n

        do k=kparasta-deepl,kparaend+deepr
            do j=0,n2+1
                do i=0,n1+1
                    umed(i,j,k,1)=umed(i,j,k,1)+u(i,j,k)*dt          ! u media*Tempo totale
                    umed(i,j,k,2)=umed(i,j,k,2)+v(i,j,k)*dt          ! v media*Tempo totale
                    umed(i,j,k,3)=umed(i,j,k,3)+w(i,j,k)*dt          ! w media*Tempo totale
                    umed(i,j,k,4)=umed(i,j,k,4)+fi(i,j,k)*dt

                    uutau(i,j,k,1)=uutau(i,j,k,1)+u(i,j,k)*u(i,j,k)*dt    ! u*u*Tempo totale
                    uutau(i,j,k,2)=uutau(i,j,k,2)+u(i,j,k)*v(i,j,k)*dt    ! u*v*Tempo totale
                    uutau(i,j,k,3)=uutau(i,j,k,3)+u(i,j,k)*w(i,j,k)*dt    ! u*w*Tempo totale
                    uutau(i,j,k,4)=uutau(i,j,k,4)+v(i,j,k)*v(i,j,k)*dt    ! v*v*Tempo totale
                    uutau(i,j,k,5)=uutau(i,j,k,5)+v(i,j,k)*w(i,j,k)*dt    ! v*w*Tempo totale
                    uutau(i,j,k,6)=uutau(i,j,k,6)+w(i,j,k)*w(i,j,k)*dt    ! w*w*Tempo totale


                    do n=1,nscal
                        umed(i,j,k,4+n)=umed(i,j,k,4+n)+rhov(n,i,j,k)*dt          ! r media*Tempo totale

                        uutau(i,j,k,6+n)=uutau(i,j,k,6+n)+rhov(n,i,j,k)*rhov(n,i,j,k)*dt  ! r*r

                        uutau(i,j,k,6+nscal+n)=uutau(i,j,k,6+nscal+n)+rhov(n,i,j,k)*v(i,j,k)*dt  ! r*v*Tempo totale

                    end do

                end do
            end do
        end do

        ! update the total time
        deltatempo=deltatempo+dt



    end subroutine update_medie

    subroutine get_cpp_strings(simulation_folder_string,grid_file_string,restart_file_string, &
        newresformat_file_string) bind(C,name='get_cpp_strings')

        use, intrinsic :: iso_c_binding

        character(kind=c_char),dimension(500),intent(IN) :: simulation_folder_string
        character(kind=c_char),dimension(500),intent(IN) :: grid_file_string
        character(kind=c_char),dimension(500),intent(IN) :: restart_file_string
        character(kind=c_char),dimension(500),intent(IN) :: newresformat_file_string

        integer :: l

        ! ----------------------------------------------------------------------

        do l=1,500
            if (simulation_folder_string(l)==C_NULL_CHAR) exit
            result_folder(l:l)=simulation_folder_string(l)
        end do
        result_folder(l:500)=" "

        do l=1,500
            if (grid_file_string(l)==C_NULL_CHAR) exit
            grid_file(l:l)=grid_file_string(l)
        end do
        grid_file(l:500)=" "

        do l=1,500
            if (restart_file_string(l)==C_NULL_CHAR) exit
            restart_file(l:l)=restart_file_string(l)
        end do
        restart_file(l:500)=" "

        ! if I write with format recognize it: set the format
        if (iformat_newres==1) then
            do l=1,500
                if (newresformat_file_string(l)==C_NULL_CHAR) exit
                string_newres_format(l:l)=newresformat_file_string(l)
            end do
            string_newres_format(l:500)=" "

            fmt_newres='('//trim(adjustl(string_newres_format))//')'

            !if (myid==0) write(*,*) 'newresformat_file_string: ',newresformat_file_string
            !if (myid==0) write(*,*) 'string_newres_format: ',string_newres_format
            !if (myid==0) write(*,*) 'New_res file format: ',fmt_newres

        end if

        if (myid==0) then
            print*, 'simulation folder: ', trim(result_folder)
            print*, 'grid file: ', trim(grid_file)
            print*, 'restart file: ', trim(restart_file)
        end if

    end subroutine get_cpp_strings

    subroutine print_inflow_planes(ktime)

        use mysettings, only: inf,niter
        use myarrays_velo3, only: u,v,w,rhov

        integer,intent(in) :: ktime
        !
        integer ipiani
        character*30 filepiano
        character*2  idproc2
        character*1  idpiano,idproc
        integer i,j,k,isc
        ! ----------------------------------------------

        if (inf) then
            !       open file
            if (ktime==1) then
                do ipiani=1,npiani
                    if (myid<10) then
                        write(idpiano,'(i1)')ipiani
                        write(idproc ,'(i1)')myid
                        filepiano='npiano'//idpiano//'processore0'//idproc//'.dat'
                        open(5000+ipiani*100+myid,file=filepiano, status='unknown')
                    else
                        write(idpiano,'(i1)')ipiani
                        write(idproc2,'(i2)')myid
                        filepiano='npiano'//idpiano//'processore'//idproc2//'.dat'
                        open(5000+ipiani*100+myid,file=filepiano, status='unknown')
                    end if
                end do
            end if

            do ipiani=1,npiani
                i=piani(ipiani)
                write(5000+ipiani*100+myid,*)ktime
                do k=kparasta,kparaend
                    do j=1,n2
                        write(5000+ipiani*100+myid,100)u(i,j,k),v(i,j,k),w(i,j,k),(rhov(isc,i,j,k),isc=1,nscal)
                    end do
                end do
            end do

            !       close file
            if (ktime==niter) then
                do ipiani=1,npiani
                    close(5000+ipiani*100+myid)
                end do
            end if
        end if

100     format(10e13.5)

    end subroutine print_inflow_planes

    subroutine print_wallstress()

        use wallmodel_module, only: att_mod_par,u_t
        use ibm_module, only: num_ib,shear_ib,pressure_ib,momentum_ib,ustar,indici_CELLE_IB
        use mysettings, only: bodyforce
        use myarrays_metri3, only: ref_length,ref_area,giac

        !-----------------------------------------------------------------------
        integer :: l
        ! tau and ustar for the walls
        real :: u_t_sotto,u_t_sopra
        real :: tau_sotto,tau_sopra
        ! tau in case of ibm
        real :: tautot,tau_loc
        real :: sheartot,shear_loc
        real :: pressuretot,pressure_loc
        real :: momentumtot,momentum_loc
        real :: ustartot,ustar_loc
        real :: refLengthHere,refAreaHere,refVolumeHere

        integer :: i,k,i0,j0,k0
        integer :: ierr
        !-----------------------------------------------------------------------

        if (bodyforce) then

            ! total stress
            shear_loc=0.0
            pressure_loc=0.0
            momentum_loc=0.0
            ustar_loc=0.0
            do l=1,num_ib

                ! index ib
                i0=indici_CELLE_IB(1,l)
                j0=indici_CELLE_IB(2,l)
                k0=indici_CELLE_IB(3,l)

                ! reference geometric feature
                refLengthHere=ref_length(i0,j0,k0)
                refAreaHere=ref_area(i0,j0,k0)
                refVolumeHere=giac(i0,j0,k0)

                ! update integral quantities
                shear_loc=shear_loc+shear_ib(1,l)*refAreaHere
                pressure_loc=pressure_loc+pressure_ib(3,l)*refAreaHere
                momentum_loc=momentum_loc+momentum_ib(1,l)*refVolumeHere
                ustar_loc=ustar_loc+ustar(l)*refAreaHere
                tau_loc=tau_loc+ustar(l)*ustar(l)

            end do


            call MPI_ALLREDUCE(tau_loc,tautot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(shear_loc,sheartot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(pressure_loc,pressuretot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(momentum_loc,momentumtot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ustar_loc,ustartot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            !tautot=tautot/(2.0*real(n1*n3))

            if (myid==0) then


                write(*,*)'Av. IBM tau',tautot/real(n1)/real(n3)
                write(*,*)'IBM ustar',ustartot
                write(*,*)'IBM shear_ib',sheartot
                write(*,*)'IBM pressure_ib',pressuretot
                write(*,*)'IBM momentum_ib',momentumtot

                write(*,*)'-----------------------------------------------'

                open(wallstress_file,file=filename_wallstress,status='old',action='write',position='append')
                write(wallstress_file,*) ti,ustartot,sheartot,pressuretot,momentumtot
                close(wallstress_file)
            end if

        else

            u_t_sotto=0.0
            u_t_sopra=0.0
            do k=kparasta,kparaend
                do i=1,n1
                    if (att_mod_par(i,1,k)) then
                        u_t_sotto=u_t_sotto+u_t(i,1,k)*u_t(i,1,k)
                    end if
                    if (att_mod_par(i,2,k)) then
                        u_t_sopra=u_t_sopra+u_t(i,2,k)*u_t(i,2,k)
                    end if

                end do
            end do

            call MPI_ALLREDUCE(u_t_sotto,tau_sotto,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(u_t_sopra,tau_sopra,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (myid==0) then

                tau_sotto=tau_sotto/real(n1)/real(n3)
                tau_sopra=tau_sopra/real(n1)/real(n3)

                write(*,*)'tau',tau_sotto,tau_sopra

                open(wallstress_file,file=filename_wallstress,status='old',action='write',position='append')
                write(wallstress_file,*) ti,tau_sotto,tau_sopra
                close(wallstress_file)

            end if
        end if

    end subroutine print_wallstress

    subroutine prepare_printdata()

        ! in these sub the variable that must be print are given to temparrays like print_u.
        ! if iprint_file=0  print_u is allocated on all the domain all procs known the flow field and just one will print
        ! if iprint_file=1 print_u is allocated locally and it is used instead of directly u just for simplicity in the writing subroutine

        use myarrays_velo3, only: fi,rhov,u,uc,v,vc,w,wc,akapt
        use myarrays_metri3, only: annit, annitV

        use mpi

        !-----------------------------------------------------------------------
        !     array declaration
        integer :: i,j,k,m,ierr,isc

        !integer requ1,requ2,reqv1,reqv2,reqw1,reqw2
        !integer reqannit1,reqannit2,reqakapt1,reqakapt2
        !integer reqr1,reqr2,reqf1,reqf2
        integer :: istatus(MPI_STATUS_SIZE)

        real, allocatable :: rho(:,:,:)
        logical,parameter :: short=.true.
        !-----------------------------------------------------------------------


        if (i_printfile==0) then

            allocate(print_u(0:n1+1,0:n2+1,0:n3+1))
            allocate(print_v(0:n1+1,0:n2+1,0:n3+1))
            allocate(print_w(0:n1+1,0:n2+1,0:n3+1))

            allocate(print_uc(0:n1,n2,n3))
            allocate(print_vc(n1,0:n2,n3))
            allocate(print_wc(n1,n2,0:n3))

            allocate(print_fi(0:n1+1,0:n2+1,0:n3+1))

            allocate(print_rhov(nscal,0:n1+1,0:n2+1,0:n3+1))

            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        print_u(i,j,k)=u(i,j,k)
                    end do
                end do
            end do


            allocate (ucol((n1+2)*(n2+2)*n3/nproc))
            allocate (utot((n1+2)*(n2+2)*n3))

            do k=kparasta,kparaend
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-kparasta))
                        ucol(m)=u(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(ucol(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD, &
                utot(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

            do k=1,n3
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-1))
                        print_u(i,j,k)=utot(m)
                    end do
                end do
            end do

            deallocate(ucol,utot)

            if (myid==nproc-1) then
                call MPI_SEND(print_u(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,0,1011,MPI_COMM_WORLD,ierr)
            else if (myid==0) then
                call MPI_RECV(print_u(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,1011,MPI_COMM_WORLD,istatus,ierr)
            end if

            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        print_v(i,j,k)=v(i,j,k)
                    end do
                end do
            end do

            allocate (vcol((n1+2)*(n2+2)*n3/nproc))
            allocate (vtot((n1+2)*(n2+2)*n3))

            do k=kparasta,kparaend
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-kparasta))
                        vcol(m)=v(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(vcol(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD, &
                vtot(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

            do k=1,n3
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-1))
                        print_v(i,j,k)=vtot(m)
                    end do
                end do
            end do

            deallocate(vcol,vtot)

            if (myid==nproc-1) then
                call MPI_SEND(print_v(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,0,1021,MPI_COMM_WORLD,ierr)
            else if (myid==0) then
                call MPI_RECV(print_v(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,1021,MPI_COMM_WORLD,istatus,ierr)
            end if
            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        print_w(i,j,k)=w(i,j,k)
                    end do
                end do
            end do


            allocate (wcol((n1+2)*(n2+2)*n3/nproc))
            allocate (wtot((n1+2)*(n2+2)*n3))

            do k=kparasta,kparaend
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-kparasta))
                        wcol(m)=w(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(wcol(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD,&
                wtot(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

            do k=1,n3
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-1))
                        print_w(i,j,k)=wtot(m)
                    end do
                end do
            end do

            deallocate(wcol,wtot)

            if (myid==nproc-1) then
                call MPI_SEND(print_w(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,0,1031,MPI_COMM_WORLD,ierr)
            else if (myid==0) then
                call MPI_RECV(print_w(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,1031,MPI_COMM_WORLD,istatus,ierr)
            end if

            !.......................................................................
            allocate(rho(0:n1+1,0:n2+1,0:n3+1))

            do isc=1,nscal

                do k=kparasta-1,kparaend+1 !0,n3+1
                    do j=0,n2+1
                        do i=0,n1+1
                            print_rhov(isc,i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do

                do k=kparasta-1,kparaend+1 !0,n3+1
                    do j=0,n2+1
                        do i=0,n1+1
                            rho(i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do

                allocate (rhocol((n1+2)*(n2+2)*n3/nproc))
                allocate (rhotot((n1+2)*(n2+2)*n3))

                do k=kparasta,kparaend
                    do j=0,n2+1
                        do i=0,n1+1
                            m=i+1+(n1+2)*(j+(n2+2)*(k-kparasta))
                            rhocol(m)=rho(i,j,k)
                        end do
                    end do
                end do

                call MPI_ALLGATHER(rhocol(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD, &
                    rhotot(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

                do k=1,n3
                    do j=0,n2+1
                        do i=0,n1+1
                            m=i+1+(n1+2)*(j+(n2+2)*(k-1))
                            rho(i,j,k)=rhotot(m)
                        end do
                    end do
                end do

                deallocate(rhocol,rhotot)

                if (myid==nproc-1) then
                    call MPI_SEND(rho(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,0,1041,MPI_COMM_WORLD,ierr)
                else if (myid==0) then
                    call MPI_RECV(rho(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,1041,MPI_COMM_WORLD,istatus,ierr)
                end if

                do k=0,n3+1
                    do j=0,n2+1
                        do i=0,n1+1
                            print_rhov(isc,i,j,k)=rho(i,j,k)
                        end do
                    end do
                end do

            end do !nscal
            deallocate(rho)
            !.......................................................................

            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        print_fi(i,j,k)=fi(i,j,k)
                    end do
                end do
            end do


            allocate (ficol((n1+2)*(n2+2)*n3/nproc))
            allocate (fitot((n1+2)*(n2+2)*n3))

            do k=kparasta,kparaend
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-kparasta))
                        ficol(m)=fi(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(ficol(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD, &
                fitot(1),(n1+2)*(n2+2)*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

            do k=1,n3
                do j=0,n2+1
                    do i=0,n1+1
                        m=i+1+(n1+2)*(j+(n2+2)*(k-1))
                        print_fi(i,j,k)=fitot(m)
                    end do
                end do
            end do

            deallocate(ficol,fitot)

            if (myid==nproc-1) then
                call MPI_SEND(print_fi(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,0,1051,MPI_COMM_WORLD,ierr)
            else if (myid==0) then
                call MPI_RECV(print_fi(0,0,n3+1),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,1051,MPI_COMM_WORLD,istatus,ierr)
            end if

            !.......................................................................
            do k=kparasta,kparaend
                do j=1,n2
                    do i=0,n1
                        print_uc(i,j,k)=uc(i,j,k)
                    end do
                end do
            end do



            allocate (uccol((n1+1)*n2*n3/nproc))
            allocate (uctot((n1+1)*n2*n3))

            do k=kparasta,kparaend
                do j=1,n2
                    do i=0,n1
                        m=i+1+(n1+1)*(j-1+n2*(k-kparasta))
                        uccol(m)=uc(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(uccol(1),(n1+1)*n2*(n3/nproc),MPI_REAL_SD, &
                uctot(1),(n1+1)*n2*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

            do k=1,n3
                do j=1,n2
                    do i=0,n1
                        m=i+1+(n1+1)*(j-1+n2*(k-1))
                        print_uc(i,j,k)=uctot(m)
                    end do
                end do
            end do

            deallocate(uccol,uctot)

            !.......................................................................
            do k=kparasta,kparaend
                do j=0,n2
                    do i=1,n1
                        print_vc(i,j,k)=vc(i,j,k)
                    end do
                end do
            end do

            allocate (vccol(n1*(n2+1)*n3/nproc))
            allocate (vctot(n1*(n2+1)*n3))

            do k=kparasta,kparaend
                do j=0,n2
                    do i=1,n1
                        m=i+n1*(j+(n2+1)*(k-kparasta))
                        vccol(m)=vc(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(vccol(1),n1*(n2+1)*(n3/nproc),MPI_REAL_SD, &
                vctot(1),n1*(n2+1)*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

            do k=1,n3
                do j=0,n2
                    do i=1,n1
                        m=i+n1*(j+(n2+1)*(k-1))
                        print_vc(i,j,k)=vctot(m)
                    end do
                end do
            end do

            deallocate(vccol,vctot)
            !.......................................................................

            do k=kparasta-1,kparaend
                do j=1,n2
                    do i=1,n1
                        print_wc(i,j,k)=wc(i,j,k)
                    end do
                end do
            end do

            allocate (wccol(n1*n2*n3/nproc))
            allocate (wctot(n1*n3*n2))

            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        m=i+n1*(j-1+n2*(k-kparasta))
                        wccol(m)=wc(i,j,k)
                    end do
                end do
            end do

            call MPI_ALLGATHER(wccol(1),n1*n2*(n3/nproc),MPI_REAL_SD, &
                wctot(1),n1*n2*(n3/nproc),MPI_REAL_SD,MPI_COMM_WORLD,ierr)

            do k=1,n3
                do j=1,n2
                    do i=1,n1
                        m=i+n1*((j-1)+n2*(k-1))
                        print_wc(i,j,k)=wctot(m)
                    end do
                end do
            end do

            deallocate(wccol,wctot)


        !-----------------------------------------------------------------------
        !     move variable to print_u etc if each procs are going to print

        else if (i_printfile==2) then

            allocate(print_u(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(print_v(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(print_w(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate(print_fi(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate(print_rhov(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        print_u(i,j,k)=u(i,j,k)
                        print_v(i,j,k)=v(i,j,k)
                        print_w(i,j,k)=w(i,j,k)
                        print_fi(i,j,k)=fi(i,j,k)

                        do isc=1,nscal
                            print_rhov(isc,i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do
            end do


            allocate(print_uc(0:n1,n2,kparasta:kparaend))
            allocate(print_vc(n1,0:n2,kparasta:kparaend))

            do k=kparasta,kparaend
                do j=1,n2
                    do i=0,n1
                        print_uc(i,j,k)=uc(i,j,k)
                    end do
                end do
            end do

            do k=kparasta,kparaend
                do j=0,n2
                    do i=1,n1
                        print_vc(i,j,k)=vc(i,j,k)
                    end do
                end do
            end do

            allocate(print_wc(n1,n2,kparasta-1:kparaend))
            do k=kparasta-1,kparaend
                do j=1,n2
                    do i=1,n1
                        print_wc(i,j,k)=wc(i,j,k)
                    end do
                end do
            end do

        end if ! if on i_printfile


        return

    end subroutine prepare_printdata

end module output_module











