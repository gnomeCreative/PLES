module strati

        !
        ! GENERALIZED COORDINATE DYNAMIC LAGRANGIAN MIXED SGS MODEL
        !
        ! ######################################################################
        ! parallel version MPI
        ! (APRIL 2011, Roman F.)
        ! ######################################################################
        ! modified by Alessandro Leonardi, starting October 2015
        ! ######################################################################
        !
        ! NavierStokes solver
        ! Kim and Moin scheme with generalized coordinates
        ! central scheme for convective term with quick option
        ! sor and line sor + multigrid for pressure
        ! esplicit / semimplicit time scheme AB or AB+CN
        !
        ! les model dynamic or static, isotropic or anisotropic
        ! check lagrangian model
        ! check scale similar part
        !
        ! INPUT:
        !  grid: gri3dp_in.dat (no format)
        !  parameters: Agenerale.in
        !              Aboundary.in
        !              Apianisonde.in
        !              Afiltraggio.in
        ! depending on parameter settings: inflow plane "piano1.dat"
        !
        !
        !              for IBM: Celle_IB_indici.inp
        !                       Celle_IB_distanze.inp
        !                       Celle_Bloccate_Indici.inp
        !                       distanze_interpolazioni.inp
        !                       rotazione.inp
        ! OUTPUT:
        !  new_res : the flow field
        !  medietempo : to make time statistics
        !
        !
        !-----------------------------------------------------------------------

    use,intrinsic :: iso_c_binding

    ! MODULE AND COMMON AREA
    use :: mysettings           ! simulation settings when not in the include
    use :: turbo_module         ! module for turbulence model
    use :: myarrays_WB          ! for wave breaking
    use :: myarrays_LC          ! for langmuir circulation
    use :: myarrays_ibm         ! for immersed boundary
    use :: myarrays_velo3
    use :: myarrays_metri3
    use :: myarrays_density
    use :: mysending
    use :: myarrays_wallmodel

    !-------------------------------------------------------------------------
    ! EXTERN SUBROUTINES TO MODULES
    !-------------------------------------------------------------------------
    use :: multigrid_module
    use :: contour_module
    use :: flu_module
    use :: jord_module
    use :: output_module
    use :: ricerca_module
    use :: buffer_bodyforce_module
    use :: inflow_module
    use :: nesting_module
    use :: filtro_module
    use :: orlansky_module
    use :: tridag_module
    !------------------------------------------------------------------------
    use :: scala3 !domain dimension + re + dt
    use :: subgrid
    use :: period ! periodicity
    use :: velpar
    !------------------------------------------------------------------------
    use :: mpi

    implicit none

    private

    !-----------------------------------------------------------------------
    ! ARRAY AND VARIABLES DECLARATION
    ! index
    integer :: kgridparasta,kgridparaend
    integer :: kpsta,kpend

    ! pressure gradient rhs^n-1 etc
    real,allocatable :: delrho(:,:,:)
    real,allocatable :: delrhov(:,:,:,:)
    real,allocatable :: bcsi(:,:,:),beta(:,:,:),bzet(:,:,:)
    real,allocatable :: f1ve(:,:,:),f2ve(:,:,:),f3ve(:,:,:)
    real,allocatable :: rho(:,:,:)
    real,allocatable :: fdve(:,:,:,:)

    ! for message passing
    integer :: ierr,status(MPI_STATUS_SIZE)
           
    ! coef species decay
    real,allocatable :: kdeg(:)
           
    ! for transposed tridiag and approximate factorization
    real,allocatable :: g33_tr(:,:,:), giac_tr(:,:,:)
    real,allocatable :: aaa(:,:),rh(:)
    real,allocatable :: aa(:),bb(:),cc(:)

    ! multigrid info
    integer  :: nlevel
    integer :: jxc(0:4),jyc(0:4),jzc(0:4)

    ! for timing to see code performance
    integer startc, endc, ratc
    real elapsed_time,start_cput,end_cput,elapsed_cput
    real :: starttime,endtime
    real :: ticks

    real, allocatable :: rho_piano(:,:,:)

    public :: les_initialize,les_core,les_finalize,les_domain_size
    private :: init_parallel,read_grid

contains

    subroutine les_initialize() bind (C, name="les_initialize")

        implicit none

        ! for reading grid and restart
        real :: val_u,val_v,val_w
        real,allocatable :: val_rhov(:)

        ! filename for wind, consider removing
        character*60 filename

        integer :: i,j,k,kk,l,isc

        ! grid has beeen read outside

        ! initialize output files and folders
        call create_output_folder()

        ! initialize parallelization variables
        call init_parallel(kgridparasta,kgridparaend,nlevmultimax)

        call initialize_input()

        !-----------------------------------------------------------------------
        ticks = MPI_WTICK() ! FOR SCALING WTIME RESULTS
        !-----------------------------------------------------------------------

        ! INITIALIZATION

        allocate(f1ve(n1,n2,kparasta:kparaend))
        allocate(f2ve(n1,n2,kparasta:kparaend))
        allocate(f3ve(n1,n2,kparasta:kparaend))

        allocate(bcsi(n1,n2,kparasta:kparaend))
        allocate(beta(n1,n2,kparasta:kparaend))
        allocate(bzet(n1,n2,kparasta:kparaend))

        allocate(delu(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(delv(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(delw(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        allocate(delrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(delrhov(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        !
        call iniz(f1ve,f2ve,f3ve,bcsi,beta,bzet)
        ti=0.

        ! species decay initialized to zero
        allocate(kdeg(nscal))
        do i=1,nscal
            kdeg(i) = 0.
        end do

        ! temp matrix to allocate periodicity in k
        call indy(nlevel,jxc,jyc,jzc)

        call init_metrica(nlevmultimax,nlevel)

        call read_grid(grid_file)
        !-----------------------------------------------------------------------
        !  read IBM data
        call set_ibm(ktime)

        !-----------------------------------------------------------------------
        ! compute total area for each sides and cell
        call intialize_myarrays2()
        call facce()

        !-----------------------------------------------------------------------
        ! allocate plane at sides
        call initialize_contour()

        !-----------------------------------------------------------------------
        ! initialization for delrhov: needed if attiva_scal = 0
        delrhov=0.

        ! variables allocation for ORLANSKY and INFLOW
        call initialize_Orlansky()

        ! variables allocation for VELPAR
        call initialize_velpar()

        ! variables allocation for SUBGRID
        call initialize_subgrid()

        !-----------------------------------------------------------------------
        ! no idea about what these are used for...
        allocate(aaa(3,n1+n2+n3),rh(n1+n2+n3))
        allocate(aa(n1+n2+n3),bb(n1+n2+n3),cc(n1+n2+n3))

        !-----------------------------------------------------------------------
        ! Check the filtering (FILTRAGGIO) parameters
        call initialize_filtro()

        !-----------------------------------------------------------------------
        ! ***** RESTART *****
        !-----------------------------------------------------------------------
        !
        if (i_rest==1) then  !start from previous solution old_res_form
            call restart(restart_file,ti)
        else if (i_rest==2) then ! start from dns interpolation
            kpsta = kparasta - deepl
            kpend = kparaend + deepr
            if (myid== 0) kpsta = 0
            if (myid== nproc-1) kpend = jz+1

            allocate(val_rhov(nscal))
            do k=0,jz+1
                do j=0,jy+1
                    do i=0,jx+1
                        read(12,string_newres_format)val_u !u(i,j,k)
                        read(12,string_newres_format)val_v !v(i,j,k)
                        read(12,string_newres_format)val_w !w(i,j,k)
                        do isc=1,nscal
                            read(12,string_newres_format)val_rhov(isc) !rhov(isc,i,j,k)
                        end do

                        if (k>= kpsta .and. k<=kpend ) then
                            u(i,j,k)   = val_u
                            v(i,j,k)   = val_v
                            w(i,j,k)   = val_w
                            do isc=1,nscal
                                rhov(isc,i,j,k) = val_rhov(isc)
                            end do
                        end if
                    end do
                end do
            end do
            deallocate(val_rhov)
        else if (i_rest==3) then ! start with nesting

            !    read side 1
            if (infout1 /=0) then
                open(81,file='parete1.dat',status='old')
            end if
            !    read side 2
            if (infout2 /=0) then
                open(82,file='parete2.dat',status='old')
            end if

            !    read side 5
            if (infout5 /=0) then
                open(85,file='parete5.dat',status='old')
            end if
            !    read side 6
            if (infout6 /=0) then
                open(86,file='parete6.dat',status='old')
            end if

            call aree_parziali()

            call prepare_nesting(ti,81,82,83,84,85,86)

        end if !i_rest


        ! eddy viscosity initialization to molecular value
        do k=kparasta-1,kparaend+1 !0,jz+1
            do j=0,jy+1
                do i=0,jx+1
                    annit(i,j,k)=1./re
                    annitV(i,j,k)=1./re
                    do isc=1,nscal
                        akapt(isc,i,j,k)=1./re/pran(isc)
                        akaptV(isc,i,j,k)=1./re/pran(isc)
                    end do
                end do
            end do
        end do

        if (myid==0 .or. myid==nproc-1) then
            do isc=1,nscal
                akapt_piano  =1./re/pran(isc)
                akaptV_piano =1./re/pran(isc)
            end do
        end if

        !-----------------------------------------------------------------------
        ! compute metric terms as in Zang Street Koseff and index at the walls
        !
        call metrica()
        call mul_met(nlevel,jxc,jyc,jzc)
        call wall(nlevel,jxc,jyc,jzc)

        !-----------------------------------------------------------------------
        ! if semimplict allocation for tridiag of g33 and giac

        iparasta=(myid* int(n1/nproc)  +1)
        iparaend=((myid+1)* int(n1/nproc))

        allocate(giac_tr(n3,n2,iparasta:iparaend))
        giac_tr = 0.
        allocate(g33_tr(0:n3,n2,iparasta:iparaend))
        g33_tr = 0.

        call set_transpose_implicit(g33_tr,giac_tr)

        !-----------------------------------------------------------------------
        ! compute initial contravariant flux if i_rest==2 dns-interpolation
        if (i_rest==2) then
            call contrin()
        end if
        ! for nesting
        if (i_rest==3.and.potenziale==1) then
            call contrin_pot()
        else if (i_rest==3.and.potenziale==0) then
            call contrin_lat()
        end if


        !-----------------------------------------------------------------------
        ! read inflow files and sett index for orlansky
        if (lett/=0) then
            call aree_parziali()

            call inflow()
        end if

        !-----------------------------------------------------------------------
        ! for nesting: redistribution of mass on controvariant fluxes
        ! to obtain a divergence free flow
        if (i_rest==3) then
            call aree_parziali()

            if (freesurface==0) then !if freesurface off
                if (myid==0) then
                    write(*,*)'freesurface is off and entering redistribuzione.'
                end if
                call redistribuzione()
            else if (freesurface==1) then !if freesurface is on
                if (myid==0) then
                    write(*,*)'free surface is on. redistribuzione skipped.'
                end if
            end if !if freesurface on/off
        end if

        !-----------------------------------------------------------------------
        ! boundary conditions
        if (i_rest==3) then
            call contour_se_nesting()
        else
            call contour_se()
        end if

        ! boundary conditions for periodicity
        call contourp_se()

        ! compute cartesian velocity and controvariant
        call update()
        !

        !-----------------------------------------------------------------------
        !
        !if (myid==0 .and.lagr==0) then
        !    write(29,*) niter/i_print
        !end if
        !
        !-----------------------------------------------------------------------
        !
        ! SPT (Single Processor Timing)
        ! to check code performance
        !
        !  call setrteopts('cpu_time_type=total_alltime')
        call cpu_time(start_cput)
        !  write (*,*)myid, 'start_cput= ',start_cput
        call system_clock(startc,ratc)
        !  write (*,*)myid, 'startc= ',startc,'ratc= ',ratc
        starttime=MPI_WTIME()

        ! variables allocation for turbulence model (turbo_statico)
        call initialize_turbo()

        !
        !-----------------------------------------------------------------------
        ! variables allocation for WB and LC
        !
        if (windyes==1) then

            allocate(u_wind(0:n1+1,kparasta-1:kparaend+1))
            allocate(w_wind(0:n1+1,kparasta-1:kparaend+1))
            ! Giulia modificavento:  alloca tauu_att
            allocate(tauu_att(0:n1+1,kparasta-1:kparaend+1))
            allocate(tauw_att(0:n1+1,kparasta-1:kparaend+1))
            ! Giulia modificavento:
            allocate(u_att(0:n1+1,kparasta-1:kparaend+1))
            allocate(v_att(0:n1+1,kparasta-1:kparaend+1))
            allocate(w_att(0:n1+1,kparasta-1:kparaend+1))
            !-----------------------------------------------------------------------
            ! Andrea: alloco e leggo i cf
            allocate(cf(1:n1,kparasta:kparaend))
            INQUIRE(FILE="cf.dat", EXIST=cf_exists)
            if (cf_exists.eqv. .true.) then
                filename='cf.dat'
                if (myid==0)write(*,*)'open the file: ',filename
                open(2,file=filename,status='old')
                if (myid==0)write(*,*)'Reading wind damping coefficients '
                do k=1,jz
                    do i=1,jx
                        read(2,*)varcf
                        if (k>=kparasta.and.k<=kparaend)cf(i,k)=varcf
                    end do
                end do
                close(2)
            !   write(*,*)'Done '
            else
                do k=kparasta,kparaend
                    do i=1,jx
                        cf(i,k)=1.
                    end do
                end do
            end if
            !-----------------------------------------------------------------------
            do k=kparasta-1,kparaend+1
                do i=0,jx+1
                    v_att(i,k)=0.
                end do
            end do

        end if

        if ((windyes==1.or.wavebk==1)) then ! .or. imoist==1

            allocate(Fx(0:n1+1,kparasta:kparaend))
            allocate(Fz(0:n1+1,kparasta:kparaend))

        end if

        allocate(vortx(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(vorty(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(vortz(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(u_drift(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(w_drift(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        allocate(ucs(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(wcs(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        if (attiva_scal==1) then
            allocate(fdve(nscal,n1,n2,kparasta:kparaend))
        end if


        ! initialization

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    vortx(i,j,k)   = 0.
                    vorty(i,j,k)   = 0.
                    vortz(i,j,k)   = 0.
                    u_drift(i,j,k) = 0.
                    w_drift(i,j,k) = 0.
                    ucs(i,j,k)     = 0.
                    wcs(i,j,k)     = 0.
                end do
            end do
        end do


        !-----------------------------------------------------------------------
        !***********************************************************************
        !-----------------------------------------------------------------------

        if (myid==0)write(*,*)myid,'end allocation'

        !-----------------------------------------------------------------------
        ! allocation for wall function

        ! comment: if no wall function I put to zero the control for each wall
        !     this is necessary if one forget to change accordingly Aboundary.in

        call initialize_wallmodel(coef_wall)

        !  nesting
        termina=.false.
        !

        !-----------------------------------------------------------------------
        ! plane comunication for periodicity in k

        do kk=1,1-kp

            call communication_velpiano()

            !..............................................................................
            allocate(rho(0:jx+1,0:jy+1,kparasta-deepl:kparaend+deepr))
            if (myid==0) then
                allocate(rho_piano(0:jx+1,0:jy+1,n3:n3))
            else if (myid==nproc-1) then
                allocate(rho_piano(0:jx+1,0:jy+1,1:1))
            end if

            do isc=1,nscal

                rho(:,:,:)=0.0

                if (myid==0 .or. myid==nproc-1) then
                    rho_piano=0.0
                end if

                do k=kparasta-deepl,kparaend+deepr
                    do j=0,jy+1
                        do i=0,jx+1
                            rho(i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do


                if (myid==nproc-1) then
                    call MPI_SSEND(rho(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
                else if (myid==0) then
                    call MPI_RECV(rho_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
                end if
                if (myid==0) then
                    call MPI_SSEND(rho(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
                end if
                if (myid==nproc-1) then
                    call MPI_RECV(rho_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
                end if


                if (myid==0) then
                    do j=0,jy+1
                        do i=0,jx+1
                            rhov_piano(isc,i,j,jz)=rho_piano(i,j,jz)
                        end do
                    end do
                else if (myid==0) then
                    do j=0,jy+1
                        do i=0,jx+1
                            rhov_piano(isc,i,j,1)=rho_piano(i,j,1)
                        end do
                    end do
                end if

            end do  ! isc

            deallocate(rho)

            if (myid==0) then
                deallocate(rho_piano)
            else if (myid==nproc-1) then
                deallocate(rho_piano)
            end if


        end do
        !

        ! THIS IS THE LEVEL SET
        !call iniz_levelset()

        ! intialize output data
        call initialize_output()

        ! compute volume
        call compute_volume(tipo) !n1,n2,deepl,deepr,myid,kparasta,kparaend

    end subroutine les_initialize

    subroutine les_core() bind ( C, name="les_core" )

        implicit none

        real :: startturbotime,endturbotime
        real :: startmultitime,endmultitime
        real :: starteqstime,endeqstime
        real :: startitertime,enditertime
        real :: mpitimestart,mpitimeend
        real :: startrhoghost,endrhoghost
        real :: startrhoa,endrhoa,startrhob,endrhob
        real :: starteq1a,endeq1a,starteq1b,endeq1b
        real :: starteq2a,endeq2a,starteq2b,endeq2b
        real :: starteq3a,endeq3a,starteq3b,endeq3b
        real :: startdiv,enddiv
        real :: startgrad,endgrad

        real :: divint,divint_loc

        real :: l_x,l_y,l_z,l_f,coef_annit

        ! iterators ad indices
        integer :: i,j,k,ii,jj,kk,i0,j0,k0,l,isc

        ! for clipping, consider removing
        real :: rhomin,rhomax

        integer :: iq1

        !***********************************************************************
        ! CYCLE
        !***********************************************************************

        startitertime=MPI_WTIME()

        !-----------------------------------------------------------------------
        ! update position of IBM nodes
        call set_ibm(ktime)

        !-----------------------------------------------------------------------
        ! boundary conditions on du, dv, dw
        call condi()

        ! boundary conditions for periodicity
        !call contourp_se()
        !call contour()

        ! average on fluxes in periodicity direction
        call average_periodicfluxes()

        !-----------------------------------------------------------------------
        ! compute the utau for wall model
        if (coef_wall==1 .and. potenziale==0) then
            call wall_function_bodyfitted(ktime,niter,tipo,i_rest)
        end if

        !-----------------------------------------------------------------------
        ! computation of turbulent viscosity and diffusion
        startturbotime=MPI_WTIME()
        call execute_turbo(ktime,i_print,i_rest,in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kpstamg,kpendmg)
        endturbotime=MPI_WTIME()
        if (myid==0 .and. nsgs/=0) then
            print*,myid,'turbo procedure: ',endturbotime-startturbotime
        end if

        !-----------------------------------------------------------------------
        !  IBM CORRECTION
        if (bodyforce .and. coef_wall>=1 .and. i_rest/=0) then
            if (ktime==1) then

                ! call correggi_ib (...)
                !correggo_rho = 0
                !correggo_delu = 0
                ipressione_ibm = 0
                call correggi_ib(ktime,tipo)

            end if

            do l=1,num_solide
                i=indici_celle_bloccate(l,1)
                j=indici_celle_bloccate(l,2)
                k=indici_celle_bloccate(l,3)
                annit(i,j,k)=1./re
                annitV(i,j,k)=1./re
            end do

            do l=1,num_ib
                i0=indici_CELLE_IB(l,1) !ib
                j0=indici_CELLE_IB(l,2)
                k0=indici_CELLE_IB(l,3)

                i=indici_CELLE_IB(l,4) !v
                j=indici_CELLE_IB(l,5)
                k=indici_CELLE_IB(l,6)


                !      if (ktime ==1 .and. i_rest/=0) then
                !       fornisco un primo valore per la ustar
                !       calcolo la velocita' tangente all' IB
                !        call vel_tangente(i,j,k,MN,MP,proiezioni,
                ! >      u(i,j,k),v(i,j,k),w(i,j,k),vtan,alfa,l)
                !        call wernerwengle(l,MN,vtan,dist_ib_parete,
                ! >                                dist_pp_ib,ustar)
                !      end if

                l_x =abs(0.25*(x(i,j,k)+x(i,j,k-1)+x(i,j-1,k-1)+x(i,j-1,k)) &
                    -0.25*(x(i-1,j,k)+x(i-1,j,k-1)+x(i-1,j-1,k-1)+x(i-1,j-1,k)))
        
                l_y =abs(0.25*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k)) &
                    -.25*(y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k)))

                l_z =abs(0.25*(z(i,j,k)+z(i,j-1,k)+z(i-1,j-1,k)+z(i-1,j,k)) &
                    -.25*(z(i,j,k-1)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j,k-1)))
     
                l_f = min(l_x,l_y,l_z)

                l_f = dist_ib_parete(l)+0.5*l_f

                if (ustar(l)>0) then
                    coef_annit = l_f / dist_ib_parete(l)
                else
                    coef_annit = 0.
                end if

                if (dist_ib_parete(l)*re*ustar(l)<11) then
                    coef_annit = 0.
                end if

                annit(i0,j0,k0)=coef_annit*0.41*ustar(l)*dist_ib_parete(l)

                if (annit(i0,j0,k0)<1./re) then
                    annit(i0,j0,k0)=1./re
                end if

                annitV(i0,j0,k0)=annit(i0,j0,k0)

            end do

        end if !end bodyforce = 1

        !-----------------------------------------------------------------------
        !  nesting
        !  set for ti the value of the first POM time
        !  if (ktime==1.and.i_rest==3) then
        !  ti=ti_pom_old
        !  end if

        ! compute dt at the first iteration
        if (ktime == 1) then
            dt_start = dt
            call courant(ktime)
        end if
      
        ! update time
        ti=ti+dt
        if (myid==0) then
            write(*,*)'----> ktime, tempo  ',ktime,ti,'----------------------'
        end if

        !-----------------------------------------------------------------------
        if (i_rest==3 .and. potenziale==0) then
            call nesting(ti,freesurface)
        else if (i_rest==3 .and. potenziale==1) then
            call redistribuzione()
        end if
 
        !-----------------------------------------------------------------------
        ! viscosity communication
        call communication_viscosity()

        !-----------------------------------------------------------------------
        ! read wind file
        if (windyes==1) then
            call leggivento(kparasta,kparaend,myid,nproc)
            if (langyes==1) then
                !      langmuir circulation
                call vorticitag(myid,nproc,kparasta,kparaend)
                call drift(myid,nproc,kparasta,kparaend)
            end if
        end if

        !-----------------------------------------------------------------------
        ! compute bodyforce acting on the fluid
        call fmassa(bcsi,beta,bzet,ktime)

        starteqstime=MPI_WTIME()

        !-----------------------------------------------------------------------
        ! nesting: generate disturbance on the inflow
        if (ibb==1) then
            call buffer_bodyforce(ti,ktime,bcsi,beta,bzet)
        end if

        !-----------------------------------------------------------------------
        ! SCALAR EQUATION
        !-----------------------------------------------------------------------

        ! if attiva_scal=0, scalar eq. solution is bypassed
        if (attiva_scal==0. .or. potenziale==1) then
      
            if (myid==0) then
                print*,myid,'no density equation'
            end if

        else

            startrhoa=MPI_WTIME()

            allocate(rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            ! cycle on n scalars
            do isc=1,nscal

                do k=kparasta-deepl,kparaend+deepr !0,n3+1
                    do j=0,n2+1
                        do i=0,n1+1
                            rho(i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do

                ! call vislam(pran,akapt)         !coefficienti di diffusione laminare

                call mixrho_para(inmodrho,rho) ! scale similar part for the model

                call flud1(uc,cgra1,rho,akapt,insc,isc,tipo2)  ! expl. term in R11 !tipo2
                call flud2(vc,cgra2,rho,akapt,insc,isc,tipo2)  ! expl. term in R22 !tipo2
                call flud3(wc,cgra3,rho,akapt,insc,isc,tipo2)  ! expl. term in R33 !tipo2

                if (espl==1) then
                    call flucrhoesp(rho,isc,tipo)  ! expl. term Crank-Nicolson

                    call ada_rho(ktime,fdve,isc,rho,kdeg(isc))   ! Adams-Bashforth
       
                    do k=kparasta,kparaend !1,jz
                        do j=1,jy
                            do i=1,jx
                                delrho(i,j,k)=rhs(i,j,k)
                            end do
                        end do
                    end do
                else

                    call ada_rho(ktime,fdve,isc,rho,kdeg(isc))   ! Adams-Bashforth

                    call flucrho(ti,rho,akaptV,akapt,isc,tipo)  ! expl. term Crank-Nicolson

                    call rhs1_rho(kparasta,kparaend)    ! right hand side scalar eq.

                    if (bodyforce) then

                        do k=kparasta,kparaend
                            do j=1,jy
                                do i=1,jx
                                    if (tipo(i,j,k)==1) then
                                        rhs(i,j,k) = 0.
                                        fdve(isc,i,j,k) = 0.
                                    end if
                                end do
                            end do
                        end do

                    end if

                    !-----------------------------------------------------------------------
                    ! approximate factorization
                    !
                    ! upload csi
                    do k=kparasta,kparaend
                        do j=1,jy
                            call coed1(j,k,delrho,aaa,rh,kparasta,kparaend,isc) !coefficent construction upload csi
                            do ii=1,jx
                                aa(ii)=aaa(1,ii)
                                bb(ii)=aaa(2,ii)
                                cc(ii)=aaa(3,ii)
                            end do
                            do ii=1,1-ip
                                call triper(aa,bb,cc,rh,jx-1)
                            end do
                            do ii=1,ip
                                call tridag(aa,bb,cc,rh,jx)
                            end do
                            do i=1,jx
                                delrho(i,j,k)=rh(i)          ! put out in delrho
                            end do
                        end do
                    end do

                    ! upload eta
                    do k=kparasta,kparaend
                        do i=1,jx
                            call coed2(i,k,delrho,aaa,rh,kparasta,kparaend,isc)! coefficent construction upload eta
                            do ii=1,jy
                                aa(ii)=aaa(1,ii)
                                bb(ii)=aaa(2,ii)
                                cc(ii)=aaa(3,ii)
                            end do
                            do jj=1,1-jp
                                call triper(aa,bb,cc,rh,jy-1)
                            end do
                            do jj=1,jp
                                call tridag(aa,bb,cc,rh,jy)
                            end do
                            do j=1,jy
                                delrho(i,j,k)=rh(j)          ! put out in delrho
                            end do
                        end do
                    end do

                    endrhoa=MPI_WTIME()
                    !
                    ! new subroutine to solve the third part of approximate factorization
                    ! with transposed method
                    !
                    startrhob=MPI_WTIME()

                    call tridiag_trasp_para_rho(akapt,g33_tr,giac_tr,delrho,akapt_piano,pran,isc)
     
                end if

                ! clipping
                if (myid == 0)write(*,*)'CLIPPING'
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=1,jx
                            delrhov(isc,i,j,k)=delrho(i,j,k)
                            rho(i,j,k)=rho(i,j,k)+delrho(i,j,k)

!                            if (rhov(1,i,j,k)<0.0) then
!                                rhov(1,i,j,k)=0.0
!                                end if
!                            if (isc==1) then
!                                rhomax=9.5
!                                rhomin=6.5
!                            else if (isc==2) then
!                                rhomax=40.00
!                                rhomin=0.00
!                            end if
!
!                            if (rho(i,j,k) < rhomin) then
!                                rho(i,j,k) = rhomin
!                            end if
!
!                            if (rho(i,j,k) > rhomax) then
!                                rho(i,j,k) = rhomax
!                            end if
                        !   if (rho(i,j,k) < 0.) then
                        !  rho(i,j,k) = 0.
                        !   end if

                        !   if (rho(i,j,k) > 1.) then
                        !  rho(i,j,k) = 1.
                        !   end if


                        end do
                    end do
                end do

                endrhob=MPI_WTIME()

                !
                ! distribution of closer plane between procs
                !
                startrhoghost=MPI_WTIME()

                if (leftpem /= MPI_PROC_NULL) then
                    call MPI_SSEND(rho(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
                    !   quick
                    if (insc==1) then
                        call MPI_SSEND(rho(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
                    end if
                end if
      
                if (rightpem /= MPI_PROC_NULL) then
                    call MPI_RECV(rho(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
                    !   quick
                    if (insc==1) then
                        call MPI_RECV(rho(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
                    end if
                end if
      
                if (rightpem /= MPI_PROC_NULL) then
                    call MPI_SSEND(rho(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
                end if
                if (leftpem /= MPI_PROC_NULL) then
                    call MPI_RECV(rho(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
                end if

                ! put on rhov the computed values
                do k=kparasta-deepl,kparaend+deepr !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            rhov(isc,i,j,k)=rho(i,j,k)
                        end do
                    end do
                end do

            end do !cyclo on scalar
      
            deallocate(rho)

        end if  !attiva_scal

        endrhoghost=MPI_WTIME()
        !
        !--------------------------------------------------------------------
        !                       MOMENTUM EQUATION
        !--------------------------------------------------------------------
        !
        ! predictor step,
        ! solution of the equation for u, v ,w to determine the
        ! intermediate flow field
        !
        !
        !-----------------------------------------------------------------------
        !                         FIRST EQUATION
        !-----------------------------------------------------------------------
        starteq1a=MPI_WTIME()
        !
        ! compute convective and diffusive explicit terms in first eq.
        iq1=1
        call mix_para(inmod,iq1,ktime,i_print,lagr)

        call jord1(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)
        call flu_turbo()

        if (langyes==1) then
            call langmuir2(myid,nproc,kparasta,kparaend)
        end if

        call flux1(uc,cgra1,u,insc,tipo)     ! explicit term in F21
        call flux2(vc,cgra2,u,insc,tipo)     ! explicit term in F12
        call flux3(wc,cgra3,u,insc,tipo)     ! explicit term in F13

        if (espl==1) then
            if (windyes==1) then
                call flucnesp(u,visualizzo,tipo,tauu_att) ! expl. term Crank-Nicolson
            else
                call flucnesp(u,visualizzo,tipo) ! expl. term Crank-Nicolson
            end if
            call adams(ktime,f1ve,bcsi,kparasta,kparaend)     ! Adams-Bashforth

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        delu(i,j,k)=rhs(i,j,k)
                    end do
                end do
            end do

        else

            call adams(ktime,f1ve,bcsi,kparasta,kparaend)     ! Adams-Bashforth
            !
            !  wall model is on in flucn only if wfp3=1
            !    tauwz(i,j,k) = tauw(i,j,k)*u(i,j,k)/(u(i,j,k)+w(i,j,k))
            eseguo34=0
            if (wfp3==1.or.wfp4==1) then
                eseguo34 = 1
            end if
            if (windyes==1) then
                call flucn(u,visualizzo,coef_wall,ti,tipo,tauu_att)
            else
                call flucn(u,visualizzo,coef_wall,ti,tipo)
            end if
            eseguo34=0
            call rhs1_rho(kparasta,kparaend)   ! right hand side of momentum eq.

            !
            !-----------------------------------------------------------------------
            !  approximate factorization
            !
            !  first eq., upload csi
            do k=kparasta,kparaend
                do j=1,jy
                    call coef1_par(j,k,delu,aaa,rh,kparasta,kparaend)  !coefficent construction upload csi
                    do ii=1,jx
                        aa(ii)=aaa(1,ii)
                        bb(ii)=aaa(2,ii)
                        cc(ii)=aaa(3,ii)
                    end do
                    do ii=1,1-ip
                        call triper(aa,bb,cc,rh,jx-1)
                    end do
                    do ii=1,ip
                        call tridag(aa,bb,cc,rh,jx)
                    end do
                    do i=1,jx
                        delu(i,j,k)=rh(i)        ! put out in delu
                    end do
                end do
            end do

            !
            !  first eq., upload eta
            do k=kparasta,kparaend
                do i=1,jx
                    call coef2_par(i,k,delu,aaa,rh,kparasta,kparaend)  !coefficent construction upload eta
                    do ii=1,jy
                        aa(ii)=aaa(1,ii)
                        bb(ii)=aaa(2,ii)
                        cc(ii)=aaa(3,ii)
                    end do
                    do jj=1,1-jp
                        call triper(aa,bb,cc,rh,jy-1)
                    end do
                    do jj=1,jp
                        call tridag(aa,bb,cc,rh,jy)
                    end do
                    do j=1,jy
                        delu(i,j,k)=rh(j)          ! put out in delu
                    end do
                end do
            end do

            endeq1a=MPI_WTIME()

            starteq1b=MPI_WTIME()

            call tridiag_trasp_para(annit,g33_tr,giac_tr,delu,ktime)!annit_piano,


            endeq1b=MPI_WTIME()

        end if



        !-----------------------------------------------------------------------
        !                         SECOND EQUATION
        !-----------------------------------------------------------------------

        starteq2a=MPI_WTIME()
        ! compute convective and diffusive explicit terms in second eq.
        !
        iq1=2
        call mix_para(inmod,iq1,ktime,i_print,lagr)

        call jord2(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)

        call flu_turbo()

        if (langyes==1) then
            call langmuir2(myid,nproc,kparasta,kparaend)
        end if

        call flux1(uc,cgra1,v,insc,tipo)     ! explicit term in F21
        call flux2(vc,cgra2,v,insc,tipo)     ! explicit term in F22
        call flux3(wc,cgra3,v,insc,tipo)     ! explicit term in F23

        if (espl==1) then
            if (windyes==1) then
                call flucnesp(v,visualizzo,tipo,v_att)
            else
                call flucnesp(v,visualizzo,tipo)
            end if
            call adams(ktime,f2ve,beta,kparasta,kparaend)    ! Adams-Bashforth

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        delv(i,j,k)=rhs(i,j,k)
                    end do
                end do
            end do
        else

            call adams(ktime,f2ve,beta,kparasta,kparaend)    ! Adams-Bashforth

            eseguo34 = 0
            if (windyes==1) then
                call flucn(v,visualizzo,coef_wall,ti,tipo,v_att)  ! expl. term Crank-Nicolson
            else
                call flucn(v,visualizzo,coef_wall,ti,tipo)  ! expl. term Crank-Nicolson
            end if
       
            call rhs1_rho(kparasta,kparaend)  ! right hand side of momentum eq.
            !
            !-----------------------------------------------------------------------
            ! approximate factorization
            !
            ! second eq., upload csi
            do k=kparasta,kparaend
                do j=1,jy
                    call coef1_par(j,k,delv,aaa,rh,kparasta,kparaend)  !coefficent construction upload csi
                    do ii=1,jx
                        aa(ii)=aaa(1,ii)
                        bb(ii)=aaa(2,ii)
                        cc(ii)=aaa(3,ii)
                    end do
                    do ii=1,1-ip
                        call triper(aa,bb,cc,rh,jx-1)
                    end do
                    do ii=1,ip
                        call tridag(aa,bb,cc,rh,jx)
                    end do
                    do i=1,jx
                        delv(i,j,k)=rh(i)          !put out in delv
                    end do
                end do
            end do
            !
            ! second eq., upload eta
            do k=kparasta,kparaend
                do i=1,jx
                    call coef2_par(i,k,delv,aaa,rh,kparasta,kparaend)  !coefficent construction upload eta
                    do ii=1,jy
                        aa(ii)=aaa(1,ii)
                        bb(ii)=aaa(2,ii)
                        cc(ii)=aaa(3,ii)
                    end do
                    do jj=1,1-jp
                        call triper(aa,bb,cc,rh,jy-1)
                    end do
                    do jj=1,jp
                        call tridag(aa,bb,cc,rh,jy)
                    end do
                    do j=1,jy
                        delv(i,j,k)=rh(j)          ! put out in delv
                    end do
                end do
            end do
            !
            endeq2a=MPI_WTIME()

            starteq2b=MPI_WTIME()

            call tridiag_trasp_para(annit,g33_tr,giac_tr,delv,ktime)!annit_piano,

            endeq2b=MPI_WTIME()

        end if

        !-----------------------------------------------------------------------
        !                         THIRD EQUATION
        !-----------------------------------------------------------------------
        starteq3a=MPI_WTIME()
      
        ! compute convective and diffusive explicit terms in third eq.
        !
        iq1=3
        call mix_para(inmod,iq1,ktime,i_print,lagr)
     
        call jord3(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)

        call flu_turbo()

        if (langyes==1) then
            call langmuir2(myid,nproc,kparasta,kparaend)
        end if

        call flux1(uc,cgra1,w,insc,tipo)     ! explicit term in F31
        call flux2(vc,cgra2,w,insc,tipo)     ! explicit term in F32
        call flux3(wc,cgra3,w,insc,tipo)     ! explicit term in F33

        if (espl==1) then
            if (windyes==1) then
                call flucnesp(w,visualizzo,tipo,tauw_att) ! expl. term in Crank-Nicolson
            else
                call flucnesp(w,visualizzo,tipo) ! expl. term in Crank-Nicolson
            end if

            call adams(ktime,f3ve,bzet,kparasta,kparaend)    ! Adams-Bashforth

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        delw(i,j,k)=rhs(i,j,k)
                    end do
                end do
            end do
        else

            call adams(ktime,f3ve,bzet,kparasta,kparaend)    ! Adams-Bashforth

            !  wall model is on in flucn only if wfp3=1
            !      tauwz(i,j,k) = tauw(i,j,k)*w(i,j,k)/(u(i,j,k)+w(i,j,k))
            eseguo34=0
            if (wfp3==1.or.wfp4==1) then
                eseguo34 = 2
            end if
            if (windyes==1) then
                call flucn(w,visualizzo,coef_wall,ti,tipo,tauw_att)
            else
                call flucn(w,visualizzo,coef_wall,ti,tipo)
            end if
            eseguo34=0

            call rhs1_rho(kparasta,kparaend)  ! right hand side of momentum eq.

            !-----------------------------------------------------------------------
            !  approximate factorization
            !
            !  third eq., upload csi
            do k=kparasta,kparaend
                do j=1,jy
                    call coef1_par(j,k,delw,aaa,rh,kparasta,kparaend)  !coefficent construction upload csi
                    do ii=1,jx
                        aa(ii)=aaa(1,ii)
                        bb(ii)=aaa(2,ii)
                        cc(ii)=aaa(3,ii)
                    end do
                    do ii=1,1-ip
                        call triper(aa,bb,cc,rh,jx-1)
                    end do
                    do ii=1,ip
                        call tridag(aa,bb,cc,rh,jx)
                    end do
                    do i=1,jx
                        delw(i,j,k)=rh(i)          !put out in delw
                    end do
                end do
            end do
            !
            !  third eq., upload eta
            do k=kparasta,kparaend
                do i=1,jx
                    call coef2_par(i,k,delw,aaa,rh,kparasta,kparaend)  !coefficent construction upload eta
                    do ii=1,jy
                        aa(ii)=aaa(1,ii)
                        bb(ii)=aaa(2,ii)
                        cc(ii)=aaa(3,ii)
                    end do
                    do jj=1,1-jp
                        call triper(aa,bb,cc,rh,jy-1)
                    end do
                    do jj=1,jp
                        call tridag(aa,bb,cc,rh,jy)
                    end do
                    do j=1,jy
                        delw(i,j,k)=rh(j)          ! put out in delw
                    end do
                end do
            end do
            !
            endeq3a=MPI_WTIME()

            starteq3b=MPI_WTIME()

            call tridiag_trasp_para(annit,g33_tr,giac_tr,delw,ktime)!annit_piano,

            endeq3b=MPI_WTIME()

        end if

        !-----------------------------------------------------------------------
        ! apply orlansky boundary condition

        startdiv=MPI_WTIME()
        !
        if (lett/=0) then
            call orlansky_generale()
        end if
        !
        !-----------------------------------------------------------------------
        ! filtering procedure for the flow field
        call filtering_procedure(ktime)

        !-----------------------------------------------------------------------
        ! update velocity u^ = u+du from time step n to intermediate one
        ! on wall found with parabolic extrapolation
        call update()
        !
        !-----------------------------------------------------------------------
        ! correction on IBM

!        if (bodyforce==1 .and. potenziale==0) then
!            if (attiva_scal==1) then
!                correggo_rho=1
!            else
!                correggo_rho=0
!            end if
!
!            correggo_delu=1
!        !         call correggi_ib (...)
!
!        end if
        !-----------------------------------------------------------------------
        ! compute intermediate controvariant component
        !
        if (lett/=0 .and. i_rest/=3) then
            call contra_infout(ktime)
        else
            call contra(kparasta,kparaend,rightpe,leftpe,nproc,myid)
        end if

        call mass_balance()

        !-----------------------------------------------------------------------
        ! COMPUTE INTERMEDIATE DIVERGENCE
        !
        call diver(kparasta,kparaend,nproc,myid)

        divint_loc=0.
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    divint_loc=divint_loc+rhs(i,j,k)
                end do
            end do
        end do
        !
        ! sum local contribution to divint
        ! with MPI_REDUCE only myid=0 knows the value
        !
        divint=0.
        call MPI_REDUCE(divint_loc,divint,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        !
        ! I know divint
        !
        if (myid==0) then
            write(*,*)'divint ',divint
        end if
        !
        enddiv=MPI_WTIME()
        endeqstime=MPI_WTIME()
        !
        !-----------------------------------------------------------------------
        ! COMPUTE COMPUTATIONAL PRESSURE
        !
        startmultitime=MPI_WTIME()

        ! there used to be also more possibilities: all removed except multi
        call multi(eps,ficycle,nlevel,jxc,jyc,jzc,islor,tipo,ktime,freesurface)

        call communication_pressure()
                                                                                     
        endmultitime=MPI_WTIME()

        if (myid==0) then
            print*,myid,'mlt ',endmultitime-startmultitime
        end if

        !-----------------------------------------------------------------------
        ! compute pressure gradient
        startgrad=MPI_WTIME()
        call gradie()

        !-----------------------------------------------------------------------
        ! compute cartesian velocity and controvariant
        call vel_up(bodyforce,freesurface)

        !-----------------------------------------------------------------------
        ! correction on IBM
        if (bodyforce) then
            if (potenziale==1) then
                ! solid cell
                do l=1,num_solide

                    i=indici_celle_bloccate(l,1)
                    j=indici_celle_bloccate(l,2)
                    k=indici_celle_bloccate(l,3)

                    u(i,j,k)=0.
                    v(i,j,k)=0.
                    w(i,j,k)=0.

                end do
                ! Giulia boundary condition before correggi_ib
                call contour() !contour(kparasta,kparaend,nproc,myid,windyes,ktime,i_rest)
                !
                call contourp() !contourp(kparasta,kparaend,nproc,myid,lett)
            else !potenziale
                !correggo_rho=0
                !correggo_delu=0

                ! call correggi_ib (...)

                !----------------------------------------------------------------------
                ! boundary condition

                call contour() !contour(kparasta,kparaend,nproc,myid,windyes,ktime,i_rest)
                !
                call contourp() !contourp(kparasta,kparaend,nproc,myid,lett)

                ipressione_ibm = 0
                call correggi_ib(ktime,tipo)

            end if !potenziale

        else !bodyforce

            call contour() !call contour(kparasta,kparaend,nproc,myid,windyes,ktime,i_rest)
            !
            call contourp() !contourp(kparasta,kparaend,nproc,myid,lett)

        end if

        !----------------------------------------------------ant 21Jan TEST
        ! giulia contour va prima di correggi_ib
        !  call contour(kparasta,kparaend,nproc,myid)
        !
        !  call contourp(kparasta,kparaend,nproc,myid)
        !----------------------------------------------------ant 21Jan TEST

        !-----------------------------------------------------------------------
        ! compute max courant or dt
        call courant(ktime)

        !-----------------------------------------------------------------------
        ! print planes for inflow
        ! --> moved to output_step

        !-----------------------------------------------------------------------
        ! inflow: read the next planes of data
        call read_inflow(ktime)

        ! update averages and Reynolds stresses
        if (i_medietempo) then
            call update_medie()
        end if

        !-----------------------------------------------------------------------
        ! pass the ghost cells for u,v,w between procs
        !
        call communication_velocity()
                                                                                     
        !-----------------------------------------------------------------------
        ! COMPUTE THE DIVERGENCE
        !-----------------------------------------------------------------------
        call check_divergence(tipo)

        endgrad=MPI_WTIME()

        !-----------------------------------------------------------------------
        mpitimestart=MPI_WTIME()
        mpitimeend=MPI_WTIME()

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        enditertime=MPI_WTIME()

        if (myid == 0) then
            print '(I2,A,E15.8)',myid,'itr ',(enditertime-startitertime)/ticks
            write(*,*)'                   '
        end if

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! END CYCLE
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !-----------------------------------------------------------------------
        ! write output
        call output_step(tipo)

    end subroutine les_core

    subroutine les_finalize() bind ( C, name="les_finalize" )

        implicit none

        real :: resolution

        call output_finalize()

        endtime=MPI_WTIME()
        resolution=MPI_WTICK()
        if (myid==0) then
            print*,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
            print*,myid,'elapsed=',endtime-starttime,'resolution=',resolution
            print*,'------------------------------'
        end if

        deallocate(u,v,w,uc,vc,wc,fi,rhov)
        deallocate(annit,annitV,akapt,akaptV)
        deallocate(x,y,z)
        deallocate(csx,csy,csz,etx,ety,etz,ztx,zty,ztz)
        deallocate(g11,g12,g13,g21,g22,g23,g31,g32,g33)
        deallocate(giac)

        deallocate(areola1,areola2,areola3,areola4,areola5,areola6)
        deallocate(gg221,gg231,gg331,gg112,gg222,gg122,gg113,gg333,gg133)

        ! deallocate
        if (inmod==1.or.inmodrho==1) then

            deallocate (m21,m22,m23,m33)
            deallocate (ucof,vcof,wcof)
            deallocate (uucof,vvcof,wwcof)
            deallocate (l11,l12,l13,l21,l22,l23,l31,l32,l33)
            deallocate (ass11,ass12,ass13,ass21,ass22,ass23,ass31,ass32,ass33)
            deallocate (piano1,piano2,piano3,piano4,piano5,piano6,piano7,piano8,piano9,piano10)
            deallocate (piano11,piano12,piano13,piano14,piano15,piano16,piano17,piano18,piano19,piano20)
            deallocate (piano21,piano22,piano23,piano24,piano25,piano26,piano27,piano28,piano29,piano30)
            deallocate (uco,vco,wco)
            deallocate (uuco,uvco,uwco,vuco,vvco,vwco,wuco,wvco,wwco)

        end if
        !
        if (inmod==1) then
            deallocate (lmf11,lmf12,lmf13,lmf21,lmf22,lmf23,lmf31,lmf32,lmf33)
            deallocate (uf,vf,wf)
            deallocate (uvcof,uwcof,vucof,vwcof,wucof,wvcof)
            deallocate (m11m,m12m,m13m,m21m,m22m,m23m,m31m,m32m,m33m)

        end if

        if (inmodrho==1) then

            deallocate (rhof,rhofl)

        end if

        deallocate (smod,smodV,smodH)
        deallocate (apcsx,apcsy,apcsz,apetx,apety,apetz,apztx,apzty)
        deallocate (rbuff1,sbuff1,pp0)

        ! matrici che servono in inverse_para2
        deallocate (pp1,pp2,pp3,pc1,pc2,pc3,p0a,p0b)
        deallocate (ap11,ap12,ap13,ap21,ap22,ap23,ap31,ap32,ap33)


        !----------------------------------------------
        !**********************************************
        !----------------------------------------------
        ! matrix deallocation

        if (wavebk==1) then
            deallocate(u_wind,w_wind)
            deallocate(u_att,w_att)
            ! Giulia modificavento:
            deallocate(tauu_att,tauw_att)
        ! Giulia modificavento:
        end if

        if ((wavebk==1.and.windyes==1)) then ! .or. imoist==1
            deallocate(Fx,Fz)
        end if

        deallocate(vortx,vorty,vortz)
        deallocate(u_drift,w_drift)

        deallocate(ucs,wcs)
      
        deallocate(pran,prsc)
      
        !-----------------------------------------------------
        !*****************************************************
        !-----------------------------------------------------

        !call MPI_FINALIZE(ierr)

        call cpu_time(end_cput)
        !ccc      write (*,*)myid, 'end_cput= ',end_cput
        !
        call system_clock(endc,ratc)
        !ccc      write (*,*)myid, 'endc= ',endc,'ratc= ',ratc
        elapsed_time=real(endc-startc)/real(ratc)
        elapsed_cput=end_cput-start_cput
      
        !  do iproc=0,nproc-1
        write (*,*)myid, 'elapsed_time= ',elapsed_time,' seconds'
        write (*,*)myid, 'elapsed_cput= ',elapsed_cput,' seconds'
    !  end do
  
    !ccc      close(1000)

    end subroutine les_finalize

    subroutine init_parallel(kgridparasta,kgridparaend,nlevmultimax)

        use output_module, only: info_run_file
        use mysending

        implicit none

        integer ierr, ierror
        integer iproc
        integer,intent(inout) :: kgridparasta,kgridparaend,nlevmultimax

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (myid==0) then
            write(*,*)'----------------------------------------'
            write(info_run_file,*)'----------------------------------------'
        end if

        !-----------------------------------------------------------------------
        if (myid==0) then
            write(*,*)'number of procs: ',nproc
            write(info_run_file,*)'number of procs: ',nproc
        end if
        write(*,*)'I am proc',myid, 'of', nproc,'procs'
        write(info_run_file,*)'I am proc',myid, 'of', nproc,'procs'
        !
        !-----------------------------------------------------------------------
        ! depending on the setting in scala3.h MPI_REAL_SD assumes
        ! the type MPI_REAL4 or MPI_REAL8
        if (single_or_double == 1) then
            call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,8,MPI_REAL_SD,ierr)
            if (myid==0) then
                write(*,*)'DOUBLE PRECISION'
                write(info_run_file,*)'DOUBLE PRECISION'
            end if
        else
            call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,4,MPI_REAL_SD,ierr)
            if (myid==0) then
                write(*,*)'SINGLE PRECISION'
                write(info_run_file,*)'SINGLE PRECISION'
            end if
        end if
        !
        !-----------------------------------------------------------------------
        ! definition of :
        !   - the domain decomposition between the processor
        !   - processor recognize left and right processor
        !   - tags for sending
        !
        !   jz must be a multiple of nproc and like nproc*2^n
        !   to allow multigrid
        !   jx must be multiple of nproc to allow transpose procedure
        ncolperproc=int(n3/nproc)

        if (myid==0) then
            write(* ,*)'n col. per PE: ',ncolperproc
            write(info_run_file,*)'n col. per PE: ',ncolperproc
        end if

        kparasta=(myid*ncolperproc+1)
        kparaend=((myid+1)*ncolperproc)

        !
        ! myid=0 has kgridparasta=0 because the number of points are odd
        !
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

        do iproc=0,nproc-1
            if (myid==iproc) then
                write(*,*)'PE: ',myid,'kparasta= ', kparasta,' kparaend= ',kparaend
                write(*,*)'right ',rightpe,'left ',leftpe
                write(info_run_file,*)'PE: ',myid,'kparasta= ', kparasta,' kparaend= ',kparaend
                write(info_run_file,*)'right ',rightpe,'left ',leftpe
            end if
        end do


        tagls=100+myid        !tag send to leftpe
        taglr=110+myid-1      !tag recv from leftpe
        tagrs=110+myid        !tag send to rightpe
        tagrr=100+myid+1      !tag recv from rightpe

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

        !
        !-----------------------------------------------------------------------
        ! CHECKING CONDITIONS ON PROCESSORS ARE MET
        !
        if (mod(n3,nproc)/= 0 .or. mod(n3,(2**nlevmultimax)) /= 0) then

            ! call MPI_ABORT(ierr)
            call MPI_ABORT(MPI_COMM_WORLD,ierr,ierror)
            error stop 'ERROR: NUM. PROC. INCORRECT'
        end if


        !-----------------------------------------------------------------------
        ! check on conflicts
!    if (imoist==1 .and. nscal < 2) then
!        if (myid==0) then
!            write(*,*)'there is a conflict between the moisture\n'// &
!                'procedure and the number of scalars,\n'// &
!                'be sure nscal>=2 in scala3.h'
!            write(info_run_file,*)'there is a conflict between the moisture\n'// &
!                'procedure and the number of scalars,\n'// &
!                'be sure nscal>=2 in scala3.h'
!        end if
!
!        stop
!    end if

        ! turn off some features
        if (potenziale==1) then
            if (myid==0) then
                write(*,*)'TURN OFF att_wm_sgs and coef_wall'
                write(info_run_file,*)'TURN OFF att_wm_sgs and coef_wall'
            end if
            att_wm_sgs = 0
            coef_wall = 0
        end if

        if (coef_wall==0) then
            if (myid==0) then
                write(*,*)'TURN OFF att_wm_sgs'
                write(info_run_file,*)'TURN OFF att_wm_sgs'
            end if
            att_wm_sgs = 0
        end if

    end subroutine init_parallel

    subroutine read_grid(grid_file)

        implicit none

        character(len=500),intent(in) :: grid_file

        ! domain size
        integer,parameter :: grid_file_id=12

        integer :: kini,kfin

        real coor_x,coor_y,coor_z
        real,allocatable :: xx1(:,:,:),yy1(:,:,:),zz1(:,:,:)
        real,allocatable :: xx2(:,:,:),yy2(:,:,:),zz2(:,:,:)

        integer :: i,j,k

        !-----------------------------------------------------------------------
        ! read the grid

        allocate(xx1(-8:n1+8,-8:n2+8,n3-8:n3))
        allocate(yy1(-8:n1+8,-8:n2+8,n3-8:n3))
        allocate(zz1(-8:n1+8,-8:n2+8,n3-8:n3))

        allocate(xx2(-8:n1+8,-8:n2+8,0:8))
        allocate(yy2(-8:n1+8,-8:n2+8,0:8))
        allocate(zz2(-8:n1+8,-8:n2+8,0:8))

        open(grid_file_id,file=trim(grid_file),status='old')
        !
        read(grid_file_id,*)alx,aly,alz
        read(grid_file_id,*)jx,jy,jz

        if (jx/=n1 .and. myid==0) then
            write(info_run_file,*)'MISMATCH GRID/CODE DIMENSION ON x',jx,n1
            write(info_run_file,*)'check gri3dp_in.dat and scala3.h in the code'
            write(*,*)'read info_run.txt'
            stop
        end if
        if (jy/=n2 .and. myid==0) then
            write(info_run_file,*)'MISMATCH GRID/CODE DIMENSION ON y',jy,n2
            write(info_run_file,*)'check gri3dp_in.dat and scala3.h in the code'
            write(*,*)'read info_run.txt'
            stop
        end if
        if (jz/=n3 .and. myid==0) then
            write(info_run_file,*)'MISMATCH GRID/CODE DIMENSION ON z',jz,n3
            write(info_run_file,*)'check gri3dp_in.dat and scala3.h in the code'
            write(*,*)'read info_run.txt'
            stop
        end if


        if (myid==0) then
            kini = kparasta - 1
            kfin = kparaend + deep_mul +1 !deepgr
        else if (myid==nproc-1) then
            kini = kparasta - deep_mul -1 !deepgl
            kfin = kparaend
        else
            kini = kparasta - deep_mul -1 !deepgl
            kfin = kparaend + deep_mul +1 !deepgr
        end if

        write(info_run_file,*)myid,'READ GRID BETWEEN',kini,kfin,kparasta,kparaend
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        do k=0,jz
            do j=0,jy
                do i=0,jx
                    read(grid_file_id,*)coor_x
                    read(grid_file_id,*)coor_y
                    read(grid_file_id,*)coor_z

                    if (k>=kini .and. k<=kfin) then
                        x(i,j,k) = coor_x
                        y(i,j,k) = coor_y
                        z(i,j,k) = coor_z
                    end if

                    if (myid==0) then
                        if (k>=(n3-8).and.k<=n3) then
                            xx1(i,j,k) = coor_x
                            yy1(i,j,k) = coor_y
                            zz1(i,j,k) = coor_z
                        end if
                    end if

                    if (myid==nproc-1) then
                        if (k>=0.and.k<=8) then
                            xx2(i,j,k) = coor_x
                            yy2(i,j,k) = coor_y
                            zz2(i,j,k) = coor_z
                        end if
                    end if

                end do
            end do
        end do

        ! close the grid file
        close(12)

        if (myid==0) then
            write(*,*)'CHECK: grid read ---> OK'
        end if

        ! compute centroid position
        if(myid==0)write(*,*)'Compute centroids'
        call compute_centroids() !jx,jy,jz,myid,nproc,kparasta,kparaend

        ! periodic cells
        call cellep(xx1,yy1,zz1,xx2,yy2,zz2,kgridparasta)

        if (myid==0) then
            write(*,*) 'Grid dimension: ',jx,jy,jz
            write(*,*) 'Array dimension: ',n1,n2,n3

            write(*,*)'CHECK: bodyforce --->',bodyforce
        end if

        deallocate(xx1,yy1,zz1,xx2,yy2,zz2)

    end subroutine read_grid

    subroutine les_domain_size()  bind (C, name="les_domain_size")

        implicit none

        !character(len=500),intent(in) :: grid_file

        ! domain size
        integer,parameter :: grid_file_id=12

        ! open grid file and read first two lines
        open(grid_file_id,file=trim(grid_file),status='old')
        !
        read(grid_file_id,*)alx,aly,alz
        read(grid_file_id,*)jx,jy,jz

        ! close the grid file
        close(12)

        ! set also this shit equal
        n1=jx
        n2=jy
        n3=jz

        if (myid==0) then
            write(*,*) 'Grid dimension: ',jx,jy,jz
        end if

    end subroutine les_domain_size

    subroutine mass_balance()

        implicit none

        real :: massa1,massa2,massa3,massa4,massa5,massa6
        real :: massa1tot,massa2tot,massa3tot,massa4tot,massa5tot,massa6tot
        real :: bilancio

        integer :: i,j,k

            !-----------------------------------------------------------------------
        ! check mass on sides 1 and 2
        massa1=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa1=massa1+uc(0,j,k)
            end do
        end do

        massa2=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa2=massa2-uc(jx,j,k)
            end do
        end do

        ! check mass on sides 3 and 4
        massa3=0.
        do k=kparasta,kparaend
            do i=1,jx
                massa3=massa3+vc(i,0,k)
            end do
        end do

        massa4=0.
        do k=kparasta,kparaend
            do i=1,jx
                massa4=massa4-vc(i,jy,k)
            end do
        end do

        ! check mass on sides 5 and 6
        massa5=0.
        if (myid==0) then
            do j=1,jy
                do i=1,jx
                    massa5=massa5+wc(i,j,0)
                end do
            end do
        end if

        massa6=0.
        if (myid==nproc-1) then
            do j=1,jy
                do i=1,jx
                    massa6=massa6-wc(i,j,jz)
                end do
            end do
        end if

        ! from local to total quantities
        call MPI_ALLREDUCE(massa1,massa1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa2,massa2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa3,massa3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa4,massa4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa5,massa5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa6,massa6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

        if (myid==0) then
            bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot
            write(*,*)'mass balance: ',bilancio
        end if

        !-----------------------------------------------------------------------
        ! check cs

        ! check mass on sides 1 and 2
        massa1=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa1=massa1+cs1(j,k)
            end do
        end do

        massa2=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa2=massa2-cs2(j,k)
            end do
        end do

        ! check mass on sides 3 and 4
        massa3=0.
        do k=kparasta,kparaend
            do i=1,jx
                massa3=massa3+cs3(i,k)
            end do
        end do

        massa4=0.
        do k=kparasta,kparaend
            do i=1,jx
                massa4=massa4-cs4(i,k)
            end do
        end do

        ! check mass on sides 5 and 6
        massa5=0.
        if (myid==0) then
            do j=1,jy
                do i=1,jx
                    massa5=massa5+cs5(i,j)
                end do
            end do
        end if

        massa6=0.
        if (myid==nproc-1) then
            do j=1,jy
                do i=1,jx
                    massa6=massa6-cs6(i,j)
                end do
            end do
        end if

        ! from local to total quantities
        call MPI_ALLREDUCE(massa1,massa1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa2,massa2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa3,massa3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa4,massa4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa5,massa5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa6,massa6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

        if (myid==0) then
            bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot
            write(*,*)'cs balance: ',bilancio
        end if

    end subroutine mass_balance

end module strati
