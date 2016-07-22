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
    !      Aboundary.in
    !      Apianisonde.in
    !      Afiltraggio.in
    !
    !  for IBM: Celle_IB_indici.inp
    !       Celle_IB_distanze.inp
    !       Celle_Bloccate_Indici.inp
    !       distanze_interpolazioni.inp
    !       rotazione.inp
    ! OUTPUT:
    ! new_res : the flow field
    ! medietempo : to make time statistics
    !
    !-----------------------------------------------------------------------

    use,intrinsic :: iso_c_binding

    use mysettings           ! simulation settings when not in the include
    use turbo_module         ! module for turbulence model
    use wind_module
    use myarrays_velo3
    use myarrays_metri3
    use mysending
    use wallmodel_module
    use ibm_module
    use multigrid_module
    use contour_module
    use flu_module
    use jord_module
    use output_module
    use ricerca_module
    use buffer_bodyforce_module
    use inflow_module
    use nesting_module
    use filtro_module
    use orlansky_module
    use tridag_module
    use particle_module
    !------------------------------------------------------------------------
    use scala3 !domain dimension + re + dt
    use subgrid
    use period ! periodicity
    use velpar
    !------------------------------------------------------------------------
    use mpi

    implicit none

    private

    !-----------------------------------------------------------------------
    ! ARRAY AND VARIABLES DECLARATION
    ! index
    integer :: kgridparasta,kgridparaend
    integer :: kpsta,kpend

    real,allocatable :: bcsi(:,:,:),beta(:,:,:),bzet(:,:,:)
    real,allocatable :: f1ve(:,:,:),f2ve(:,:,:),f3ve(:,:,:)
    real,allocatable :: rho(:,:,:)
    real,allocatable :: rho_piano(:,:,:)
    real,allocatable :: fdve(:,:,:,:)

    ! for message passing
    integer :: ierr,status(MPI_STATUS_SIZE)
           
    ! coef species decay
    real,allocatable :: kdeg(:)
           
    ! multigrid info
    integer  :: nlevel
    integer :: jxc(0:4),jyc(0:4),jzc(0:4)

    ! for timing to see code performance
    real :: starttime,endtime
    real :: ticks

    public :: les_initialize_part1,les_initialize_part2,les_core,les_finalize,les_domain_size

contains

    subroutine les_initialize_part1() bind (C, name="les_initialize_part1")

        !-----------------------------------------------------------------------

        ! initialize output files and folders
        call create_output_folder()

        ! initialize parallelization variables
        call init_parallel(kgridparasta,kgridparaend,nlevmultimax,bodyforce,insc)

        ! process imput variables that need processing
        call initialize_input()

        ! print info about domain decomposition, processor, etc...
        call print_info()

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

        call iniz(f1ve,f2ve,f3ve,bcsi,beta,bzet)

        ti=0.0

        ! species decay initialized to zero
        allocate(kdeg(nscal))
        kdeg(:) = 0.

        ! temp matrix to allocate periodicity in k
        call indy(nlevel,jxc,jyc,jzc)

        call init_metrica(nlevmultimax,nlevel)

        call read_grid(grid_file)

    end subroutine les_initialize_part1

    subroutine les_initialize_part2() bind (C, name="les_initialize_part2")

        use mysending
        use output_module

                !-----------------------------------------------------------------------
        real :: ustar_prov,speed ! DELETE!!!

        ! for reading grid and restart
        real :: val_u,val_v,val_w
        real,allocatable :: val_rhov(:)

        ! filename for wind, consider removing
        character*60 filename

        integer :: i,j,k,kk,isc

        !-----------------------------------------------------------------------
        !  read/produce IBM data
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
            if (myid== nproc-1) kpend = n3+1

            allocate(val_rhov(nscal))
            do k=0,n3+1
                do j=0,n2+1
                    do i=0,n1+1
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
            do j=0,n2+1
                do i=0,n1+1
                    annit(i,j,k)=1./re
                    annitV(i,j,k)=1./re
                    do isc=1,nscal
                        akapt(isc,i,j,k)=1./re/pran(isc)
                        akaptV(isc,i,j,k)=1./re/pran(isc)
                    end do
                end do
            end do
        end do

        !-----------------------------------------------------------------------
        ! compute metric terms as in Zang Street Koseff and index at the walls
        !
        call metrica()
        call mul_met(nlevel,jxc,jyc,jzc)
        call wall(nlevel,jxc,jyc,jzc)

        ! variables allocation and initialization for approximate factorization
        call initialize_tridiag()

        !-----------------------------------------------------------------------
        ! compute initial contravariant flux if i_rest==2 dns-interpolation
        if (i_rest==2) then
            call contrin()
        end if
        ! for nesting
        if (i_rest==3.and.potenziale) then
            call contrin_pot()
        else if (i_rest==3.and.potenziale) then
            call contrin_lat()
        end if


        !-----------------------------------------------------------------------
        ! read inflow files and sett index for orlansky
        if (lett) then
            call aree_parziali()

            call inflow()
        end if

        !-----------------------------------------------------------------------
        ! for nesting: redistribution of mass on controvariant fluxes
        ! to obtain a divergence free flow
        if (i_rest==3) then
            call aree_parziali()
            if (freesurface) then
                if (myid==0) then
                    write(*,*) 'free surface is on. redistribuzione skipped.'
                end if
            else !if freesurface off
                if (myid==0) then
                    write(*,*) 'freesurface is off and entering redistribuzione.'
                end if
                call redistribuzione()
            end if
        end if

        !-----------------------------------------------------------------------
        ! boundary conditions
        if (myid==0) then
            write(*,*)'call contour_se'
        end if
        if (i_rest==3) then
            call contour_se_nesting()
        else
            call contour_se()
        end if

        ! boundary conditions for periodicity
        if (myid==0) then
            write(*,*)'call contourp_se'
        end if
        call contourp_se()


        ! compute cartesian velocity and controvariant
        if (myid==0) then
            write(*,*)'call update'
        end if
        call update()

        !-----------------------------------------------------------------------
        !
        ! SPT (Single Processor Timing)
        ! to check code performance
        starttime=MPI_WTIME()

        ! variables allocation for turbulence model (turbo_statico)
        call initialize_turbo()

        !
        !-----------------------------------------------------------------------
        ! variables allocation for WB and LC
        !
        if (windyes==1) then

            allocate(tauu_att(0:n1+1,kparasta-1:kparaend+1))
            allocate(tauw_att(0:n1+1,kparasta-1:kparaend+1))
            !-----------------------------------------------------------------------
            ! Andrea: alloco e leggo i cf
            allocate(cf(1:n1,kparasta:kparaend))
            allocate(v_att(0:n1+1,kparasta-1:kparaend+1))

            INQUIRE(FILE="cf.dat", EXIST=cf_exists)
            if (cf_exists.eqv. .true.) then
                filename='cf.dat'
                if (myid==0)write(*,*)'open the file: ',filename
                open(2,file=filename,status='old')
                if (myid==0)write(*,*)'Reading wind damping coefficients '
                do k=1,n3
                    do i=1,n1
                        read(2,*)varcf
                        if (k>=kparasta.and.k<=kparaend)cf(i,k)=varcf
                    end do
                end do
                close(2)
            !   write(*,*)'Done '
            else
                do k=kparasta,kparaend
                    do i=1,n1
                        cf(i,k)=1.
                    end do
                end do
            end if
            !-----------------------------------------------------------------------
            do k=kparasta-1,kparaend+1
                do i=0,n1+1
                    v_att(i,k)=0.
                end do
            end do

        end if

        if (attiva_scal) then
            allocate(fdve(nscal,n1,n2,kparasta:kparaend))
        end if


        !-----------------------------------------------------------------------
        !***********************************************************************
        !-----------------------------------------------------------------------

        if (myid==0) write(*,*) myid,'end allocation'

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
            allocate(rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            if (myid==0) then
                allocate(rho_piano(0:n1+1,0:n2+1,n3:n3))
            else if (myid==nproc-1) then
                allocate(rho_piano(0:n1+1,0:n2+1,1:1))
            end if

            do isc=1,nscal

                rho(:,:,:)=0.0

                if (myid==0 .or. myid==nproc-1) then
                    rho_piano=0.0
                end if

                do k=kparasta-deepl,kparaend+deepr
                    do j=0,n2+1
                        do i=0,n1+1
                            rho(i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do


                if (myid==nproc-1) then
                    call MPI_SEND(rho(0,0,n3),(n1+2)*(n2+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
                else if (myid==0) then
                    call MPI_RECV(rho_piano(0,0,n3),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
                end if
                if (myid==0) then
                    call MPI_SEND(rho(0,0,1),(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
                end if
                if (myid==nproc-1) then
                    call MPI_RECV(rho_piano(0,0,1),(n1+2)*(n2+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
                end if


                if (myid==0) then
                    do j=0,n2+1
                        do i=0,n1+1
                            rhov_piano(isc,i,j,n3)=rho_piano(i,j,n3)
                        end do
                    end do
                else if (myid==0) then
                    do j=0,n2+1
                        do i=0,n1+1
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

        if (myid==0) then
            write(*,*) '**********************************************************'
            write(*,*) '**********************************************************'
            write(*,*) '**********************************************************'
            write(*,*) 'Wall function test:'
            ustar_prov=ustar_lawofwall(0.01,1.e6)
            write(*,*) 'log ustar=',ustar_prov
            speed=speed_lawofwall(0.01,ustar_prov)
            write(*,*) 'log speed=',speed
            ustar_prov=ustar_wernerwengle(0.01,1.e6)
            write(*,*) 'ww ustar=',ustar_prov
            speed=speed_wernerwengle(0.01,ustar_prov)
            write(*,*) 'ww speed=',speed

            write(*,*) '**********************************************************'
            write(*,*) '**********************************************************'
            write(*,*) '**********************************************************'
        end if

    end subroutine les_initialize_part2

    subroutine les_core() bind ( C, name="les_core" )

        implicit none

        real :: iter_time_start,iter_time_finish
        real :: ricerca_time_start,ricerca_time_finish
        real :: bc_time_start,bc_time_finish
        real :: ibm_time_start,ibm_time_finish
        real :: turbo_time_start,turbo_time_finish
        real :: scalar_time_start,scalar_time_finish
        real :: momentum_time_start,momentum_time_finish
        real :: pressure_time_start,pressure_time_finish

        real :: divint,divint_loc

        real :: l_x,l_y,l_z,l_f,coef_annit

        ! iterators ad indices
        integer :: i,j,k,i0,j0,k0,l,isc

        ! for clipping, consider removing
        !real :: rhomin,rhomax

        integer :: iq1
        logical :: scalar


        !-----------------------------------------------------------------------

        iter_time_start=MPI_WTIME()

        !-----------------------------------------------------------------------
        ! update position of IBM nodes if necessary

        if (particles) then
            ricerca_time_start=MPI_WTIME()
            call set_ibm(ktime)
            ricerca_time_finish=MPI_WTIME()
            if (myid==0) then
                write(*,*) 'Ricerca time: ',ricerca_time_finish-ricerca_time_start
            end if
        end if

        !-----------------------------------------------------------------------
        ! boundary conditions on du, dv, dw
        bc_time_start=MPI_WTIME()
        call condi()

        ! boundary conditions for periodicity
        !call contourp_se()
        !call contour()

        ! average on fluxes in periodicity direction
        call average_periodicfluxes()

        !-----------------------------------------------------------------------
        ! compute the utau for wall model
        if (coef_wall==1 .and. .not.potenziale) then
            call correggi_walls(ktime,tipo,i_rest)
            !call wall_function_bodyfitted(ktime,niter,tipo,i_rest)
        end if
        bc_time_finish=MPI_WTIME()
        if (myid==0) then
            write(*,*) 'BC time: ',bc_time_finish-bc_time_start
        end if


        !-----------------------------------------------------------------------
        ! computation of turbulent viscosity and diffusion
        turbo_time_start=MPI_WTIME()
        call execute_turbo(ktime,i_rest,in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kpstamg,kpendmg)
        turbo_time_finish=MPI_WTIME()
        if (myid==0 .and. nsgs/=0) then
            write(*,*) 'Turbo time: ',turbo_time_finish-turbo_time_start
        end if

        !-----------------------------------------------------------------------
        !  IBM CORRECTION
        if (bodyforce .and. coef_wall>=1 .and. i_rest/=0) then

            if (ktime==1) then

                call correggi_ib(ktime,tipo)

            end if

            do l=1,num_solide
                i=indici_celle_bloccate(1,l)
                j=indici_celle_bloccate(2,l)
                k=indici_celle_bloccate(3,l)
                annit(i,j,k)=1./re
                annitV(i,j,k)=1./re
            end do

            do l=1,num_ib
!                i0=indici_CELLE_IB(1,l) !ib
!                j0=indici_CELLE_IB(2,l)
!                k0=indici_CELLE_IB(3,l)
!
!                i=indici_CELLE_IB(4,l) !v
!                j=indici_CELLE_IB(5,l)
!                k=indici_CELLE_IB(6,l)
!
!
!                !      if (ktime ==1 .and. i_rest/=0) then
!                !       fornisco un primo valore per la ustar
!                !       calcolo la velocita' tangente all' IB
!                !        call vel_tangente(i,j,k,MN,MP,proiezioni,
!                ! >      u(i,j,k),v(i,j,k),w(i,j,k),vtan,alfa,l)
!                !        call wernerwengle(l,MN,vtan,dist_ib_parete,
!                ! >                                dist_pp_ib,ustar)
!                !      end if
!
!                l_x =abs(0.25*(x(i,j,k)+x(i,j,k-1)+x(i,j-1,k-1)+x(i,j-1,k)) &
!                    -0.25*(x(i-1,j,k)+x(i-1,j,k-1)+x(i-1,j-1,k-1)+x(i-1,j-1,k)))
!
!                l_y =abs(0.25*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k)) &
!                    -.25*(y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k)))
!
!                l_z =abs(0.25*(z(i,j,k)+z(i,j-1,k)+z(i-1,j-1,k)+z(i-1,j,k)) &
!                    -.25*(z(i,j,k-1)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j,k-1)))
!
!                l_f = min(l_x,l_y,l_z)
!
!                l_f = dist_ib_parete(l)+0.5*l_f
!
!                if (ustar(l)>0) then
!                    coef_annit = l_f / dist_ib_parete(l)
!                else
!                    coef_annit = 0.
!                end if
!
!                if (dist_ib_parete(l)*re*ustar(l)<11) then
!                    coef_annit = 0.
!                end if
!
!                annit(i0,j0,k0)=coef_annit*0.41*ustar(l)*dist_ib_parete(l)
!
!                if (annit(i0,j0,k0)<1./re) then
!                    annit(i0,j0,k0)=1./re
!                end if
!
!                annitV(i0,j0,k0)=annit(i0,j0,k0)

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
        if (i_rest==3 .and. .not.potenziale) then
            call nesting(ti)
        else if (i_rest==3 .and. potenziale) then
            call redistribuzione()
        end if
 
        !-----------------------------------------------------------------------
        ! viscosity communication
        call communication_viscosity()

        !-----------------------------------------------------------------------
        ! read wind file
        if (windyes==1) then
            call leggivento()
        end if

        !-----------------------------------------------------------------------
        ! compute bodyforce acting on the fluid
        call fmassa(bcsi,beta,bzet,ktime)

        !-----------------------------------------------------------------------
        ! nesting: generate disturbance on the inflow
        if (ibb) then
            call buffer_bodyforce(ti,ktime,bcsi,beta,bzet)
        end if

        !-----------------------------------------------------------------------
        ! SCALAR EQUATION
        !-----------------------------------------------------------------------

        ! if attiva_scal=0, scalar eq. solution is bypassed
        if (.not.attiva_scal .or. potenziale) then
      
            if (myid==0) then
                write(*,*) 'no density equation'
            end if

        else

            scalar_time_start=MPI_WTIME()

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

                call mixrho_para(rho) ! scale similar part for the model

                call flud1(uc,cgra1,rho,isc,tipo2)  ! expl. term in R11 !tipo2
                call flud2(vc,cgra2,rho,isc,tipo2)  ! expl. term in R22 !tipo2
                call flud3(wc,cgra3,rho,isc,tipo2)  ! expl. term in R33 !tipo2

                if (espl==1) then
                    call flucrhoesp(rho,isc,tipo)  ! expl. term Crank-Nicolson

                    call ada_rho(ktime,fdve,isc,rho,kdeg(isc))   ! Adams-Bashforth
       
                    do k=kparasta,kparaend !1,jz
                        do j=1,n2
                            do i=1,n1
                                delrho(i,j,k)=rhs(i,j,k)
                            end do
                        end do
                    end do
                else

                    call ada_rho(ktime,fdve,isc,rho,kdeg(isc))   ! Adams-Bashforth

                    call flucrho(rho,isc,tipo)  ! expl. term Crank-Nicolson

                    call rhs1_rho(kparasta,kparaend)    ! right hand side scalar eq.

                    if (bodyforce) then

                        do k=kparasta,kparaend
                            do j=1,n2
                                do i=1,n1
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
                    scalar=.true.
                    call factorization(scalar,ktime,delrho,isc)
     
                end if

                ! clipping
                if (myid == 0)write(*,*)'CLIPPING'
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=1,n1
                            delrhov(isc,i,j,k)=delrho(i,j,k)
                            rho(i,j,k)=rho(i,j,k)+delrho(i,j,k)
                            rho(i,j,k)=1.0

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

                ! distribution of closer plane between procs
                if (leftpem /= MPI_PROC_NULL) then
                    call MPI_SEND(rho(0,0,kparasta),(n1+2)*(n2+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
                    !   quick
                    if (insc==1) then
                        call MPI_SEND(rho(0,0,kparasta+1),(n1+2)*(n2+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
                    end if
                end if
      
                if (rightpem /= MPI_PROC_NULL) then
                    call MPI_RECV(rho(0,0,kparaend+1),(n1+2)*(n2+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
                    !   quick
                    if (insc==1) then
                        call MPI_RECV(rho(0,0,kparaend+2),(n1+2)*(n2+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
                    end if
                end if
      
                if (rightpem /= MPI_PROC_NULL) then
                    call MPI_SEND(rho(0,0,kparaend),(n1+2)*(n2+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
                end if
                if (leftpem /= MPI_PROC_NULL) then
                    call MPI_RECV(rho(0,0,kparasta-1),(n1+2)*(n2+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
                end if

                ! put on rhov the computed values
                do k=kparasta-deepl,kparaend+deepr !0,jz+1
                    do j=0,n2+1
                        do i=0,n1+1
                            rhov(isc,i,j,k)=rho(i,j,k)
                        end do
                    end do
                end do

            end do !cyclo on scalar
      
            deallocate(rho)

            scalar_time_finish=MPI_WTIME()

            if (myid==0) then
                write(*,*) 'Scalar equation time: ',scalar_time_finish-scalar_time_start
            end if

        end if  !attiva_scal

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
        momentum_time_start=MPI_WTIME()
        !
        ! compute convective and diffusive explicit terms in first eq.
        iq1=1
        call mix_para(iq1)

        call jord1(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)
        call flu_turbo()

        call flux1(uc,cgra1,u,tipo)     ! explicit term in F21
        call flux2(vc,cgra2,u,tipo)     ! explicit term in F12
        call flux3(wc,cgra3,u,tipo)     ! explicit term in F13


        if (espl==1) then
            if (windyes==1) then
                !call flucnesp(u,tipo,tauu_att) ! expl. term Crank-Nicolson
                call flucn(u,espl,coef_wall,tipo,tauu_att) ! expl. term Crank-Nicolson
            else
                !call flucnesp(u,tipo) ! expl. term Crank-Nicolson
                call flucn(u,espl,coef_wall,tipo) ! expl. term Crank-Nicolson
            end if
            call adams(ktime,f1ve,bcsi)     ! Adams-Bashforth

            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        delu(i,j,k)=rhs(i,j,k)
                    end do
                end do
            end do

        else

            call adams(ktime,f1ve,bcsi)     ! Adams-Bashforth

            !  wall model is on in flucn only if wfp3=1
            !    tauwz(i,j,k) = tauw(i,j,k)*u(i,j,k)/(u(i,j,k)+w(i,j,k))
            eseguo34=0
            if (wfp3.or.wfp4) then
                eseguo34 = 1
            end if
            if (windyes==1) then
                !call flucn(u,coef_wall,tipo,tauu_att)
                call flucn(u,espl,coef_wall,tipo,tauu_att)
            else
                call flucn(u,espl,coef_wall,tipo)
                !call flucn(u,coef_wall,tipo)
            end if
            eseguo34=0
            call rhs1_rho(kparasta,kparaend)   ! right hand side of momentum eq.

            scalar=.false.
            call factorization(scalar,ktime,delu,isc)

        end if


        !-----------------------------------------------------------------------
        !                         SECOND EQUATION
        !-----------------------------------------------------------------------

        ! compute convective and diffusive explicit terms in second eq.
        !
        iq1=2
        call mix_para(iq1)

        call jord2(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)

        call flu_turbo()

        call flux1(uc,cgra1,v,tipo)     ! explicit term in F21
        call flux2(vc,cgra2,v,tipo)     ! explicit term in F22
        call flux3(wc,cgra3,v,tipo)     ! explicit term in F23

        if (espl==1) then
            ! expl. term in Crank-Nicolson
            if (windyes==1) then
                !call flucnesp(v,tipo,v_att)
                call flucn(v,espl,coef_wall,tipo,v_att)
            else
                !call flucnesp(v,tipo)
                call flucn(v,espl,coef_wall,tipo)
            end if
            call adams(ktime,f2ve,beta)    ! Adams-Bashforth

            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        delv(i,j,k)=rhs(i,j,k)
                    end do
                end do
            end do
        else

            call adams(ktime,f2ve,beta)    ! Adams-Bashforth

            eseguo34=0
            ! expl. term Crank-Nicolson
            if (windyes==1) then
                !call flucn(v,coef_wall,tipo,v_att)
                call flucn(v,espl,coef_wall,tipo,v_att)
            else
                !call flucn(v,coef_wall,tipo)
                call flucn(v,espl,coef_wall,tipo)
            end if
       
            call rhs1_rho(kparasta,kparaend)  ! right hand side of momentum eq.

            scalar=.false.
            call factorization(scalar,ktime,delv,isc)

        end if

        !-----------------------------------------------------------------------
        !                         THIRD EQUATION
        !-----------------------------------------------------------------------
      
        ! compute convective and diffusive explicit terms in third eq.
        !
        iq1=3
        call mix_para(iq1)
     
        call jord3(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)

        call flu_turbo()

        call flux1(uc,cgra1,w,tipo)     ! explicit term in F31
        call flux2(vc,cgra2,w,tipo)     ! explicit term in F32
        call flux3(wc,cgra3,w,tipo)     ! explicit term in F33

        if (espl==1) then
            ! expl. term in Crank-Nicolson
            if (windyes==1) then
                !call flucnesp(w,tipo,tauw_att)
                call flucn(w,espl,coef_wall,tipo,tauw_att)
            else
                !call flucnesp(w,tipo)
                call flucn(w,espl,coef_wall,tipo)
            end if

            call adams(ktime,f3ve,bzet)    ! Adams-Bashforth

            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        delw(i,j,k)=rhs(i,j,k)
                    end do
                end do
            end do
        else

            call adams(ktime,f3ve,bzet)    ! Adams-Bashforth

            !  wall model is on in flucn only if wfp3=1
            !      tauwz(i,j,k) = tauw(i,j,k)*w(i,j,k)/(u(i,j,k)+w(i,j,k))
            eseguo34=0
            if (wfp3.or.wfp4) then
                eseguo34=2
            end if
            ! expl. term in Crank-Nicolson
            if (windyes==1) then
                !call flucn(w,coef_wall,tipo,tauw_att)
                call flucn(w,espl,coef_wall,tipo,tauw_att)
            else
                !call flucn(w,coef_wall,tipo)
                call flucn(w,espl,coef_wall,tipo)
            end if
            eseguo34=0

            call rhs1_rho(kparasta,kparaend)  ! right hand side of momentum eq.

            scalar=.false.
            call factorization(scalar,ktime,delw,isc)

        end if
        momentum_time_finish=MPI_WTIME()

        if (myid==0) then
            write(*,*) 'Momentum equation time: ',momentum_time_finish-momentum_time_start
        end if



        !-----------------------------------------------------------------------
        ! apply orlansky boundary condition
        if (lett) then
            call orlansky_generale(giac)
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
        if (lett .and. i_rest/=3) then
            call contra_infout(ktime)
        else
            call contra()
        end if

        call mass_balance()

        !-----------------------------------------------------------------------
        ! COMPUTE INTERMEDIATE DIVERGENCE
        !
        call diver()

        divint_loc=0.
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    divint_loc=divint_loc+rhs(i,j,k)
                end do
            end do
        end do
        !
        ! sum local contribution to divint
        ! with MPI_REDUCE only myid=0 knows the value
        !
        divint=0.0
        call MPI_REDUCE(divint_loc,divint,1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        ! I know divint
        !
        if (myid==0) then
            write(*,*)'divint ',divint
        end if


        !-----------------------------------------------------------------------
        ! COMPUTE COMPUTATIONAL PRESSURE
        !
        pressure_time_start=MPI_WTIME()

        ! there used to be also more possibilities: all removed except multi
        call multi(eps,ficycle,nlevel,jxc,jyc,jzc,tipo,ktime)

        call communication_pressure()
                                                                                     
        pressure_time_finish=MPI_WTIME()

        if (myid==0) then
            write(*,*) 'Multigrid time ',pressure_time_finish-pressure_time_start
        end if

        !-----------------------------------------------------------------------
        ! compute pressure gradient
        call gradie()

        !-----------------------------------------------------------------------
        ! compute cartesian velocity and controvariant
        call vel_up()

        !-----------------------------------------------------------------------
        ! correction on IBM
        if (bodyforce) then
            if (potenziale) then
                ! solid cell
                do l=1,num_solide

                    i=indici_celle_bloccate(1,l)
                    j=indici_celle_bloccate(2,l)
                    k=indici_celle_bloccate(3,l)

                    u(i,j,k)=0.
                    v(i,j,k)=0.
                    w(i,j,k)=0.

                end do
                ! Giulia boundary condition before correggi_ib
                call contour()
                call contourp()
            else

                !----------------------------------------------------------------------
                ! boundary condition

                call contour()
                call contourp()

                ricerca_time_start=MPI_WTIME()

                call correggi_ib(ktime,tipo)

                ! compute forces on the spheres
                if (particles) then
                    call compute_sphere_forces()
                end if

                ricerca_time_finish=MPI_WTIME()
                if (myid==0) then
                    write(*,*) 'IBM time: ',ibm_time_finish-ibm_time_start
                end if

            end if

        else

            call contour()
            call contourp()

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
        ! inflow: read the next planes of data
        call read_inflow(ktime)

        ! update averages and Reynolds stresses
        if (i_medietempo .or. paraview_aver) then
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

        iter_time_finish=MPI_WTIME()

        if (myid == 0) then
            write(*,*) 'Iteration time ',iter_time_finish-iter_time_start
        end if

        !-----------------------------------------------------------------------
        ! write output
        call output_step(tipo)
    !1234 continue
    end subroutine les_core

    subroutine les_finalize() bind ( C, name="les_finalize" )

        implicit none

        real :: resolution

        ! ---------------------------------------------------------

        call output_finalize()

        endtime=MPI_WTIME()
        resolution=MPI_WTICK()
        if (myid==0) then
            write(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
            write(*,*) 'elapsed=',endtime-starttime,'resolution=',resolution
            write(*,*) '------------------------------'
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
        if (inmod.or.inmodrho) then

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
        if (inmod) then
            deallocate (lmf11,lmf12,lmf13,lmf21,lmf22,lmf23,lmf31,lmf32,lmf33)
            deallocate (uf,vf,wf)
            deallocate (uvcof,uwcof,vucof,vwcof,wucof,wvcof)
            deallocate (m11m,m12m,m13m,m21m,m22m,m23m,m31m,m32m,m33m)

        end if

        if (inmodrho) then

            deallocate (rhof,rhofl)

        end if

        deallocate (smod,smodV,smodH,q_crit)
        deallocate (apcsx,apcsy,apcsz,apetx,apety,apetz,apztx,apzty)
        deallocate (rbuff1,sbuff1,pp0)

        ! matrici che servono in inverse_para2
        deallocate (pp1,pp2,pp3,pc1,pc2,pc3,p0a,p0b)
        deallocate (ap11,ap12,ap13,ap21,ap22,ap23,ap31,ap32,ap33)


        !----------------------------------------------
        !**********************************************
        !----------------------------------------------
        ! matrix deallocation

        if (windyes==1) then
            ! Giulia modificavento:
            deallocate(tauu_att,tauw_att)
        end if

        deallocate(pran,prsc)


    end subroutine les_finalize

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
        read(grid_file_id,*)n1,n2,n3

        if (n1/=n1 .and. myid==0) then
            write(info_run_file,*)'MISMATCH GRID/CODE DIMENSION ON x',n1,n1
            write(info_run_file,*)'check gri3dp_in.dat and scala3.h in the code'
            write(*,*)'read info_run.txt'
            stop
        end if
        if (n2/=n2 .and. myid==0) then
            write(info_run_file,*)'MISMATCH GRID/CODE DIMENSION ON y',n2,n2
            write(info_run_file,*)'check gri3dp_in.dat and scala3.h in the code'
            write(*,*)'read info_run.txt'
            stop
        end if
        if (n3/=n3 .and. myid==0) then
            write(info_run_file,*)'MISMATCH GRID/CODE DIMENSION ON z',n3,n3
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

        do k=0,n3
            do j=0,n2
                do i=0,n1
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
            write(*,*) 'Grid dimension: ',n1,n2,n3
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
        read(grid_file_id,*)n1,n2,n3

        ! close the grid file
        close(12)

        if (myid==0) then
            write(*,*) 'Grid dimension: ',n1,n2,n3
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
            do j=1,n2
                massa1=massa1+uc(0,j,k)
            end do
        end do

        massa2=0.
        do k=kparasta,kparaend
            do j=1,n2
                massa2=massa2-uc(n1,j,k)
            end do
        end do

        ! check mass on sides 3 and 4
        massa3=0.
        do k=kparasta,kparaend
            do i=1,n1
                massa3=massa3+vc(i,0,k)
            end do
        end do

        massa4=0.
        do k=kparasta,kparaend
            do i=1,n1
                massa4=massa4-vc(i,n2,k)
            end do
        end do

        ! check mass on sides 5 and 6
        massa5=0.
        if (myid==0) then
            do j=1,n2
                do i=1,n1
                    massa5=massa5+wc(i,j,0)
                end do
            end do
        end if

        massa6=0.
        if (myid==nproc-1) then
            do j=1,n2
                do i=1,n1
                    massa6=massa6-wc(i,j,n3)
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
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (myid==0) then
            bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot
            write(*,*)'mass balance: ',bilancio
        end if

        !-----------------------------------------------------------------------
        ! check cs

        ! check mass on sides 1 and 2
        massa1=0.
        do k=kparasta,kparaend
            do j=1,n2
                massa1=massa1+cs1(j,k)
            end do
        end do

        massa2=0.
        do k=kparasta,kparaend
            do j=1,n2
                massa2=massa2-cs2(j,k)
            end do
        end do

        ! check mass on sides 3 and 4
        massa3=0.
        do k=kparasta,kparaend
            do i=1,n1
                massa3=massa3+cs3(i,k)
            end do
        end do

        massa4=0.
        do k=kparasta,kparaend
            do i=1,n1
                massa4=massa4-cs4(i,k)
            end do
        end do

        ! check mass on sides 5 and 6
        massa5=0.
        if (myid==0) then
            do j=1,n2
                do i=1,n1
                    massa5=massa5+cs5(i,j)
                end do
            end do
        end if

        massa6=0.
        if (myid==nproc-1) then
            do j=1,n2
                do i=1,n1
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
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (myid==0) then
            bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot
            write(*,*)'cs balance: ',bilancio
        end if

    end subroutine mass_balance

end module strati
