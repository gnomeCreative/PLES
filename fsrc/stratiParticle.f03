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
        !      grid: gri3dp_in.dat (no format)
        !      parameters: Agenerale.in
        !                  Aboundary.in
        !                  Apianisonde.in
        !                  Afiltraggio.in
        !     depending on parameter settings: inflow plane "piano1.dat"
        !
        !
        !                  for IBM: Celle_IB_indici.inp
        !                           Celle_IB_distanze.inp
        !                           Celle_Bloccate_Indici.inp
        !                           distanze_interpolazioni.inp
        !                           rotazione.inp
        ! OUTPUT:
        !      new_res : the flow field
        !      medietempo : to make time statistics
        !
        !
        !-----------------------------------------------------------------------

    use iso_c_binding
    ! MODULE AND COMMON AREA
    use mysettings           ! simulation settings when not in the include

    use turbo2_data          ! module for turbulence model
    use turbo3bis
    use myarrays2
    use myarrays_WB          ! for wave breaking
    use myarrays_LC          ! for langmuir circulation
    use myarrays_cor         ! for coriolis
    use myarrays_ibm         ! for immersed boundary
    use myarrays_velo3
    use myarrays_metri3
    use myarrays_density
    use mysending
    use myarrays_wallmodel
    use myarrays_nesting
    use myarrays_buffer_bodyforce
    use myarrays_moisture

    !-------------------------------------------------------------------------
    ! EXTERN SUBROUTINES TO MODULES
    !-------------------------------------------------------------------------
    use multigrid_module
    use contour_module
    use flucn_module
    use jord_module
    use output_module
    use ricerca_module

    !------------------------------------------------------------------------
    use scala3 !domain dimension + re + dt
    use subgrid
    use period !periodicity
    use convex
    use print
    use mpi ! to initialize MPI
    use tipologia! for the data type MPI_REAL_SD
    use orl ! for orlansky condition
    use velpar
    ! only used here, consider removing
    use parti
    use parete

    implicit none


    !-----------------------------------------------------------------------
    ! ARRAY AND VARIABLES DECLARATION (not in common area)
    !
    ! index
    integer i,j,k,i0,j0,k0,ii,jj,kk
    integer l,n,iproc,mmm
    integer iq1,isc
    integer nlevel
    integer jxc(0:4),jyc(0:4),jzc(0:4)
    integer kpstamg(0:4),kpendmg(0:4),kfacepstamg(0:4)
    integer kgridparasta,kgridparaend
    integer kpsta,kpend,kini,kfin
            
    ! quantities
    double precision divint,divint_loc
    real dbbx,dbby,dbbz
    real um_loc,um,umm1
    !double precision um_loc,um,umm1
    real massa1,massa2,massa3,massa4,massa5,massa6
    real massa1tot,massa2tot
    real massa3tot,massa4tot
    real massa5tot,massa6tot
    real bilancio
    real area1,area2,area3,area4,area5,area6
                  
    ! pressure gradient rhs^n-1 etc
    real,allocatable :: delrho(:,:,:)
    real,allocatable :: delrhov(:,:,:,:)
    real,allocatable :: bcsi(:,:,:),beta(:,:,:),bzet(:,:,:)
    real,allocatable :: f1ve(:,:,:),f2ve(:,:,:),f3ve(:,:,:)
    real,allocatable :: rho(:,:,:)
    real,allocatable :: fdve(:,:,:,:)

    ! for message passing
    integer ierr,ierror
    integer ncolprocmg(0:4)
    integer kparastam,kparaendm
    integer status(MPI_STATUS_SIZE)
    !integer istatus(MPI_STATUS_SIZE)
           
    !integer :: correggo_rho,correggo_delu

    integer, allocatable :: tipo(:,:,:),tipo2(:,:,:)

    integer itiposta,itipoend
    integer jtiposta,jtipoend
    integer ktiposta,ktipoend

    real l_x,l_y,l_z,l_f,coef_annit
    
    ! coef species decay
    real,allocatable :: kdeg(:)
           
    ! for reading grid and restart
    real val_u,val_v,val_w
    real, allocatable :: val_rhov(:)

    ! for transposed tridiag and approximate factorization
    real,allocatable :: g33_tr(:,:,:), giac_tr(:,:,:)
    real,allocatable :: aaa(:,:),rh(:)
    real,allocatable :: aa(:),bb(:),cc(:)

    ! for planes and tracers
    integer ipiani
    integer pianoturbo
    character*30 filepiano
    character*2  identificosonda2
    character*2  idproc2
    character*1  idpiano,idproc

    ! for filtering
    integer loopfiltro

    ! for timing to see code performance
    integer startc, endc, ratc
    real elapsed_time,start_cput,end_cput,elapsed_cput
    double precision starttime,endtime,resolution
    double precision startturbotime,endturbotime
    double precision startmultitime,endmultitime
    double precision starteqstime,endeqstime
    double precision startitertime,enditertime
    double precision mpitimestart,mpitimeend
    double precision startrhoghost,endrhoghost
    double precision startrhoa,endrhoa
    double precision startrhob,endrhob
    double precision starteq1a,endeq1a
    double precision starteq1b,endeq1b
    double precision starteq2a,endeq2a
    double precision starteq2b,endeq2b
    double precision starteq3a,endeq3a
    double precision starteq3b,endeq3b
    double precision startdiv,enddiv
    double precision startgrad,endgrad
    double precision start_part,end_part
    double precision :: ticks

    ! definizione matrice T per temperatura reale e per la densita'
    ! chicco AAA allocare tra kparasta e kparaend e mettere in modulo
    ! real tpotm(1:n2), qm(1:n2)
    real tm_tot, qm_tot
    real tm_loc, qm_loc
    real dist_mean
      
    character*60 filename
    real rhomin,rhomax
    real, allocatable :: rho_piano(:,:,:)
      
    integer count_print_time

contains

    subroutine les_initialize() bind ( C, name="les_initialize" )

        implicit none

        !-----------------------------------------------------------------------
        ! initialize output files and folders
        call create_output_folder()

        ! initialize parallelization variables
        call init_parallel(kgridparasta,kgridparaend,nlevmultimax)

        ! intialize output data
        call output_init()

        !-----------------------------------------------------------------------
        ticks = MPI_WTICK() ! FOR SCALING WTIME RESULTS
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        ! read imput data
        !call read_simulation_setting()
        !
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

        allocate(tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(tipo2(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        tipo  = 0
        tipo2 = 0
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

        if (myid==0)write(*,*)'NLEVEL: ',nlevel

        !do n=0,nlevel
        !    jxc(n) = jxc(n)
        !    jyc(n) = jyc(n)
        !    jzc(n) = jzc(n)
        !end do

        call iniz_metrica(nlevel)

        call read_grid()

        !-----------------------------------------------------------------------
        !  read IBM data

        if (bodyforce>=1) then

            if (bodyupdate) then

                ! read geometry input
                call readGeometry()

                ! do the ricerca cycle (to be changed A LOT!)
                call do_ricerca(tipo)


            else

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

                call carico_immb(tipo)
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

        tipo2 = 2
        do k=kparasta,kparaend !1,jz
            do j=1,jy
                do i=1,jx

                    itiposta = 1
                    itipoend = 1
                    jtiposta = 1
                    jtipoend = 1
                    ktiposta = 1
                    ktipoend = 1

                    if (i==1)itiposta=0
                    if (j==1)jtiposta=0
                    if (k==1)ktiposta=0

                    if (i==jx)itipoend=0
                    if (j==jy)jtipoend=0
                    if (k==jz)ktipoend=0


                    if (tipo(i,j,k)==2) then
                        do kk=k-ktiposta,k+ktipoend
                            do jj=j-jtiposta,j+jtipoend
                                do ii=i-itiposta,i+itipoend

                                    if (  tipo(ii,jj,kk)==1 ) then
                                        tipo2(i,j,k) = 0    ! fluid close to ib
                                    end if

                                end do
                            end do
                        end do
                    end if

                    if (tipo(i,j,k)==0)tipo2(i,j,k)=0
                    if (tipo(i,j,k)==1)tipo2(i,j,k)=0

                end do
            end do
        end do

        !-----------------------------------------------------------------------
        ! compute total area for each sides and cell

        allocate(areola1(0:n2+1,kparasta-1:kparaend+1))
        allocate(areola2(0:n2+1,kparasta-1:kparaend+1))
        allocate(areola3(0:n1+1,kparasta-1:kparaend+1))
        allocate(areola4(0:n1+1,kparasta-1:kparaend+1))
        allocate(areola5(0:n1+1,0:n2+1))
        allocate(areola6(0:n1+1,0:n2+1))
        !
        allocate(gg221(0:n2+1,kparasta-1:kparaend+1))
        allocate(gg331(0:n2+1,kparasta-1:kparaend+1))
        allocate(gg231(0:n2+1,kparasta-1:kparaend+1))
        allocate(gg112(0:n1+1,0:n2+1))
        allocate(gg222(0:n1+1,0:n2+1))
        allocate(gg122(0:n1+1,0:n2+1))
        allocate(gg113(0:n1+1,kparasta-1:kparaend+1))
        allocate(gg133(0:n1+1,kparasta-1:kparaend+1))
        allocate(gg333(0:n1+1,kparasta-1:kparaend+1))
        !

        call facce(myid,nproc,kparasta,kparaend,area1,area2,area3,area4,area5,area6)
        !-----------------------------------------------------------------------
        ! allocate plane at sides

        ! face 1
        allocate(up1(0:jy+1,0:jz+1))
        allocate(vp1(0:jy+1,0:jz+1))
        allocate(wp1(0:jy+1,0:jz+1))
        allocate(rhovp1(nscal,0:jy+1,0:jz+1))
        allocate(tke1(0:jy+1,kparasta-1:kparaend+1)) !face 1
        tke1 = 0.
        up1=0.
        vp1=0.
        wp1=0.
        rhovp1=0.

        ! face 2
        allocate(up2(0:jy+1,0:jz+1))
        allocate(vp2(0:jy+1,0:jz+1))
        allocate(wp2(0:jy+1,0:jz+1))
        allocate(rhovp2(nscal,0:jy+1,0:jz+1))
        allocate(tke2(0:jy+1,kparasta-1:kparaend+1)) !face 2
        tke2 = 0.
        up2=0.
        vp2=0.
        wp2=0.
        rhovp2=0.

        ! face 3
        allocate(up3(0:jx+1,0:jz+1))
        allocate(vp3(0:jx+1,0:jz+1))
        allocate(wp3(0:jx+1,0:jz+1))
        allocate(rhovp3(nscal,0:jx+1,0:jz+1))
        up3   = 0.
        vp3   = 0.
        wp3   = 0.
        rhovp3 = 0.

        ! face 4
        allocate(up4(0:jx+1,0:jz+1))
        allocate(vp4(0:jx+1,0:jz+1))
        allocate(wp4(0:jx+1,0:jz+1))
        allocate(rhovp4(nscal,0:jx+1,0:jz+1))
        up4   = 0.
        vp4   = 0.
        wp4   = 0.
        rhovp4 = 0.

        ! face 5
        allocate(up5(0:jx+1,0:jy+1))
        allocate(vp5(0:jx+1,0:jy+1))
        allocate(wp5(0:jx+1,0:jy+1))
        allocate(rhovp5(nscal,0:jx+1,0:jy+1))
        allocate(tke5(0:jx+1,0:jy+1))
        tke5 = 0.
        up5=0.
        vp5=0.
        wp5=0.
        rhovp5=0.

        ! face 6
        allocate(up6(0:jx+1,0:jy+1))
        allocate(vp6(0:jx+1,0:jy+1))
        allocate(wp6(0:jx+1,0:jy+1))
        allocate(rhovp6(nscal,0:jx+1,0:jy+1))
        allocate(tke6(0:jx+1,0:jy+1))
        tke6 = 0.
        up6=0.
        vp6=0.
        wp6=0.
        rhovp6=0.

        !-----------------------------------------------------------------------
        ! initialization for delrhov: needed if attiva_scal = 0
        delrhov=0.

        !-----------------------------------------------------------------------
        ! variables allocation for ORLANSKY and INFLOW

        allocate(du_dx1(n2,n3),dv_dx1(n2,n3),dw_dx1(n2,n3))
        allocate(du_dx2(n2,n3),dv_dx2(n2,n3),dw_dx2(n2,n3))

        allocate(du_dy3(n1,n3),dv_dy3(n1,n3),dw_dy3(n1,n3))
        allocate(du_dy4(n1,n3),dv_dy4(n1,n3),dw_dy4(n1,n3))

        allocate(du_dz5(n1,n2),dv_dz5(n1,n2),dw_dz5(n1,n2))
        allocate(du_dz6(n1,n2),dv_dz6(n1,n2),dw_dz6(n1,n2))

        allocate(drho_dx1(nscal,n2,n3),drho_dx2(nscal,n2,n3))
        allocate(drho_dy3(nscal,n1,n3),drho_dy4(nscal,n1,n3))
        allocate(drho_dz5(nscal,n1,n2),drho_dz6(nscal,n1,n2))

        allocate(index_out1(0:n2+1,0:n3+1),index_out2(0:n2+1,0:n3+1))
        allocate(index_out3(0:n1+1,0:n3+1),index_out4(0:n1+1,0:n3+1))
        allocate(index_out5(0:n1+1,0:n2+1),index_out6(0:n1+1,0:n2+1))

        allocate(index_rho1(0:n2+1,0:n3+1),index_rho2(0:n2+1,0:n3+1))
        allocate(index_rho3(0:n1+1,0:n3+1),index_rho4(0:n1+1,0:n3+1))
        allocate(index_rho5(0:n1+1,0:n2+1),index_rho6(0:n1+1,0:n2+1))

        index_out1=1
        index_out2=1

        index_rho1=1
        index_rho2=1

        index_out3=1
        index_out4=1

        index_rho3=1
        index_rho4=1

        index_out5=1
        index_out6=1

        index_rho5=1
        index_rho6=1

        !-----------------------------------------------------------------------
        ! variables allocation for VELPAR

        allocate(usn(n2,n3),vsn(n2,n3),wsn(n2,n3))
        allocate(udx(n2,n3),vdx(n2,n3),wdx(n2,n3))
        allocate(rhosn(nscal,n2,n3),rhodx(nscal,n2,n3))
        !
        allocate(usp(n1,n3),vsp(n1,n3),wsp(n1,n3))
        allocate(ust(n1,n3),vst(n1,n3),wst(n1,n3))
        allocate(rhosp(nscal,n1,n3),rhost(nscal,n1,n3))
        !
        allocate(uav(n1,n2),vav(n1,n2),wav(n1,n2))
        allocate(uin(n1,n2),vin(n1,n2),win(n1,n2))
        allocate(rhoav(nscal,n1,n2),rhoin(nscal,n1,n2))


        !-----------------------------------------------------------------------
        ! variables allocation for SUBGRID

        allocate(sub(n2),sub11(n2),sub22(n2),sub33(n2))
        allocate(sub12(n2),sub13(n2),sub23(n2))
        allocate(sus(n2),sus11(n2),sus22(n2),sus33(n2))
        allocate(sus12(n2),sus13(n2),sus23(n2))
        allocate(subrho11(n2),subrho22(n2),subrho33(n2))
        allocate(susrho11(n2),susrho22(n2),susrho33(n2))
        allocate(c11(n2),c22(n2),c33(n2))

        !-----------------------------------------------------------------------
        ! no idea about what these are used for...

        allocate(aaa(3,n1+n2+n3),rh(n1+n2+n3))
        allocate(aa(n1+n2+n3),bb(n1+n2+n3),cc(n1+n2+n3))

        !-----------------------------------------------------------------------
        ! variable allocation for the scalar equations
        allocate(pran(nscal))
        allocate(prsc(nscal))

        ! convert values from c++ format to fortran format
        call C_F_POINTER(c_pran,pran,[nscal])
        call C_F_POINTER(c_prsc,prsc,[nscal])

        ! variable allocation for piani and sonde and convert (see above)
        allocate(piani(npiani))
        call C_F_POINTER(c_piani,piani,[npiani])
        allocate(sonde(3,nsonde))
        allocate(sondeindexi(nsonde),sondeindexj(nsonde),sondeindexk(nsonde))
        call C_F_POINTER(c_sondeindexi,sondeindexi,[nsonde])
        call C_F_POINTER(c_sondeindexj,sondeindexj,[nsonde])
        call C_F_POINTER(c_sondeindexk,sondeindexk,[nsonde])

        if (myid==0) then
            write (*,*) "Total piani = ",npiani
            do i=1,npiani
                write(*,*) "Piano number ",i," x=",piani(i)
            end do

            write (*,*) "Total sonde = ",nsonde
            do i=1,nsonde
                write(*,*) "Sonda number ",i,"point:(",sondeindexi(i),",",sondeindexj(i),",",sondeindexk(i),")"
            end do
        end if

        !-----------------------------------------------------------------------
        ! Check the filtering (FILTRAGGIO) parameters

        if(ifiltro .eq. 1)then

            if(myid .eq. 0)then
                write(*,*)'PAY ATTENTION YOU ARE FILTERING THE FLOW FIELD'
                write(info_run_file,*)'FILTERING ON',ifiltro
                if(xend.gt.n1 .or. yend.gt.n2 .or.zend.gt.n3 &
                    .or. xstart.lt.1 .or. ystart.lt.1 .or. zstart.lt.1)then
                    write(*,*)'FILTERING AREA TOO LARGE'
                    write(*,*)'OUT OF BOUNDS'
                    write(*,*)'xend:',xend,',n1:',n1
                    write(*,*)'yend:',yend,',n2:',n2
                    write(*,*)'zend:',zend,',n3:',n3
                    write(*,*)'xend:',xstart,1
                    write(*,*)'yend:',ystart,1
                    write(*,*)'zend:',zstart,1

                    write(info_run_file,*)'FILTERING AREA TOO LARGE'
                    write(info_run_file,*)'OUT OF BOUNDS'
                    write(info_run_file,*)'xend:',xend,',n1:',n1
                    write(info_run_file,*)'yend:',yend,',n2:',n2
                    write(info_run_file,*)'zend:',zend,',n3:',n3
                    write(info_run_file,*)'xend:',xstart,1
                    write(info_run_file,*)'yend:',ystart,1
                    write(info_run_file,*)'zend:',zstart,1

                    stop

                end if
            end if

            if(zstart .le. kparasta)then
                zstart=kparasta
            end if
            if(zend .ge. kparaend)then
                zend = kparaend
            endif

            if(zstart .gt. kparaend .or. zend .lt. kparasta)then
                zstart = kparaend
                zend =   kparaend -1
                write(*,*)'proc: ',myid,'no filtering',zstart,zend
                write(info_run_file,*)'proc: ',myid,'no filtering',zstart,zend
            else
                write(*,*) 'proc:',myid,'filter in k between',zstart,'and',zend
                write(info_run_file,*)'proc:',myid,'filter in k between',zstart,'and',zend
            end if

        else
            if(myid .eq. 0)then
                write(*,*)'FILTERING OFF',ifiltro
                write(info_run_file,*)'FILTERING OFF',ifiltro
            end if
        end if

        !-----------------------------------------------------------------------
        ! ***** RESTART *****
        !-----------------------------------------------------------------------
        !
        if (i_rest==1) then  !start from previous solution old_res_form
            call restart(ti,dbbx)
        elseif (i_rest==2) then ! start from dns interpolation
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
        elseif (i_rest==3) then ! start with nesting

            !        read side 1
            if (infout1 /=0) then
                open(81,file='parete1.dat',status='old')
            end if
            !        read side 2
            if (infout2 /=0) then
                open(82,file='parete2.dat',status='old')
            end if

            if (imoist==1) then
                !        read side 3
                if (infout3 /=0) then
                    open(83,file='BCbottom.dat',status='old')
                end if
                !        read side 4
                if (infout4 /=0) then
                    open(84,file='BCupper.dat',status='old')
                end if
            end if

            !        read side 5
            if (infout5 /=0) then
                open(85,file='parete5.dat',status='old')
            end if
            !        read side 6
            if (infout6 /=0) then
                open(86,file='parete6.dat',status='old')
            end if

            call aree_parziali(myid,nproc,lett,i_rest,area1,area2,area3,area4,area5,area6,kparasta,kparaend)

            call prepare_nesting(ti,dbbx,81,82,83,84,85,86)

        endif !i_rest
        !     close the grid file
        close(12)


        !
        ! eddy viscosity initialization to molecular value
        !
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
        ! index computation for multigrid
        !
        !      call indy(jx,jy,jz,nlevel,jxc,jyc,jzc)
        !      do n=0,nlevel
        !      jxc(n) = jxc(n)
        !      jyc(n) = jyc(n)
        !      jzc(n) = jzc(n)
        !      end do
        !
        !-----------------------------------------------------------------------
        ! compute metric terms as in Zang Street Koseff and index at the walls
        !
        call metrica()
        call mul_met(nlevel,jxc,jyc,jzc)
        call wall(nlevel,jxc,jyc,jzc)
        !compute volume
        call compute_volume(n1,n2,deepl,deepr,myid,kparasta,kparaend,tipo)
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
        !
        if (i_rest==2) then
            call contrin()
        end if
        !     for nesting
        if (i_rest==3.and.potenziale==1) then
            call contrin_pot()
        elseif (i_rest==3.and.potenziale==0) then
            call contrin_lat()
        end if


        !-----------------------------------------------------------------------
        ! read inflow files and sett index for orlansky
        if (lett/=0) then
            call aree_parziali(myid,nproc,lett,i_rest,area1,area2,area3,area4,area5,area6,kparasta,kparaend)

            call inflow(myid,nproc,lett,bodyforce,area1,area2,area3,area4,area5,area6,kparasta,kparaend)
        end if

        !-----------------------------------------------------------------------
        ! for nesting: redistribution of mass on controvariant fluxes
        ! to obtain a divergence free flow
        if (i_rest==3) then
            call aree_parziali(myid,nproc,lett,i_rest,area1,area2,area3,area4,area5,area6,kparasta,kparaend)

            if (freesurface==0) then !if freesurface off
                if (myid==0) then
                    write(*,*)'freesurface is off and entering redistribuzione.'
                end if
                call redistribuzione(bodyforce,area1,area2,area5,area6,kparasta,kparaend,myid,nproc)
            elseif (freesurface==1) then !if freesurface is on
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
        !endif
        !
        !-----------------------------------------------------------------------
        !
        ! SPT (Single Processor Timing)
        ! to check code performance
        !
        !      call setrteopts('cpu_time_type=total_alltime')
        call cpu_time(start_cput)
        !      write (*,*)myid, 'start_cput= ',start_cput
        call system_clock(startc,ratc)
        !      write (*,*)myid, 'startc= ',startc,'ratc= ',ratc
        starttime=MPI_WTIME()

        !-----------------------------------------------------------------------
        !***********************************************************************
        !-----------------------------------------------------------------------
        ! variables allocation for turbulence model (turbo_statico)
        !
        if (myid==0)write(*,*)myid,'start allocation for turbo_statico'

        if (inmod==1 .or. inmodrho==1 .or. nsgs>=2) then

            allocate (m21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (ucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (l11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (ass11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (piano4(0:n1+1,0:n2+1))
            allocate (piano5(0:n1+1,0:n2+1))
            allocate (piano6(0:n1+1,0:n2+1))
            allocate (piano7(0:n1+1,0:n2+1))
            allocate (piano8(0:n1+1,0:n2+1))
            allocate (piano9(0:n1+1,0:n2+1))
            allocate (piano10(0:n1+1,0:n2+1))
            allocate (piano11(0:n1+1,0:n2+1))
            allocate (piano12(0:n1+1,0:n2+1))
            allocate (piano13(0:n1+1,0:n2+1))
            allocate (piano14(0:n1+1,0:n2+1))
            allocate (piano15(0:n1+1,0:n2+1))
            allocate (piano16(0:n1+1,0:n2+1))
            allocate (piano17(0:n1+1,0:n2+1))
            allocate (piano18(0:n1+1,0:n2+1))
            allocate (piano19(0:n1+1,0:n2+1))
            allocate (piano20(0:n1+1,0:n2+1))
            allocate (piano21(0:n1+1,0:n2+1))
            allocate (piano22(0:n1+1,0:n2+1))

            allocate (uco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uuco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (uvco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (uwco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vuco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vvco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vwco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wuco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wvco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wwco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        end if

        if (inmod==1 .or. nsgs >=2) then

            allocate (lmf11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (uwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (m11m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m12m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m13m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m21m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m22m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m23m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m31m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m32m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m33m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        end if

        if (inmodrho==1 .or. nsgs>=2) then
            allocate (rhof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhofl(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        end if
        !
        allocate (piano1(0:n1+1,0:n2+1))
        allocate (piano2(0:n1+1,0:n2+1))
        allocate (piano3(0:n1+1,0:n2+1))
        !
        !     now turbo_statico matrix
        !
        allocate(smod(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(smodV(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(smodH(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        !

        allocate (apcsx(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apcsy(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apcsz(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apetx(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apety(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apetz(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apztx(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apzty(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apztz(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        !
        allocate (pp0(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pc1(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pc2(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pc3(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (p0a(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (p0b(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        allocate (rbuff1((n1+2)*(n2+2)*40))
        allocate (sbuff1((n1+2)*(n2+2)*40))

        !-----------------------------------------------------------------------
        !***********************************************************************
        !-----------------------------------------------------------------------
        !     allocation for dynamic procedure

        if (nsgs>=2) then


            allocate (s11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s12f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s13f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s21f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s23f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s31f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s32f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods12f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods13f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods21f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods23f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods31f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods32f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (rho11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodrho11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodrho22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodrho33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (m11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (mrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (mrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (mrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (amrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (amrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (amrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))


            allocate (rhoucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhovcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhowcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !

            allocate (lrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (sgs11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgsrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgsrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgsrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !
            allocate (   appo1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))
            allocate (   appo2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))
            allocate (appo1rho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))
            allocate (appo2rho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))

            if (myid==0) then
                allocate(   appo1_piano(0:n1+1,0:n2+1,jz:jz)) !0:n3+1))
                allocate(   appo2_piano(0:n1+1,0:n2+1,jz:jz)) !0:n3+1))
                allocate(appo1rho_piano(0:n1+1,0:n2+1,jz:jz)) !0:n3+1))
                allocate(appo2rho_piano(0:n1+1,0:n2+1,jz:jz)) !0:n3+1))
            elseif (myid==nproc-1) then
                allocate(   appo1_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
                allocate(   appo2_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
                allocate(appo1rho_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
                allocate(appo2rho_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
            end if

            allocate(   alalm(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(   alamm(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(alalmrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(alammrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1))


            if (inmod==1) then

                allocate (assrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (assrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (assrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (al11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (alrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (alrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (alrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (lmfb11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

                allocate (rhofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (ucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (uucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (uvcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (uwcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vvcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vwcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wvcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wwcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhoucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhovcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhowcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

                allocate (ufucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (ufvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (ufwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vfucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vfvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vfwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wfucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wfvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wfwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhofucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhofvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhofwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (bl11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (blrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (blrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (blrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            end if

            !      allocate (piano22(0:n1+1,0:n2+1))
            allocate (piano23(0:n1+1,0:n2+1))
            allocate (piano24(0:n1+1,0:n2+1))
            allocate (piano25(0:n1+1,0:n2+1))
            allocate (piano26(0:n1+1,0:n2+1))
            allocate (piano27(0:n1+1,0:n2+1))
            allocate (piano28(0:n1+1,0:n2+1))
            allocate (piano29(0:n1+1,0:n2+1))
            allocate (piano30(0:n1+1,0:n2+1))
            !
            !
            !     allocation for turbo3bis:
            !
            allocate (s11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !
            !      allocate (rhofl(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho11(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho22(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho33(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhouco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhovco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhowco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        end if ! dynamic

        !
        !-----------------------------------------------------------------------
        !***********************************************************************
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
            !       write(*,*)'Done '
            else
                do k=kparasta,kparaend
                    do i=1,jx
                        cf(i,k)=1.
                    end do
                end do
            endif
            !-----------------------------------------------------------------------
            do k=kparasta-1,kparaend+1
                do i=0,jx+1
                    v_att(i,k)=0.
                end do
            end do

        end if

        if ((windyes==1.or.wavebk==1) .or. imoist==1) then

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

        if (indm==0) then
            dbbx=0.
            dbby=0.
            dbbz=0.
        end if

        !-----------------------------------------------------------------------
        ! allocation for wall function

        ! comment: if no wall function I put to zero the control for each wall
        !         this is necessary if one forget to change accordingly Aboundary.in

        if (coef_wall == 0) then
            wfp1 = 0
            wfp2 = 0
            wfp3 = 0
            wfp4 = 0
            wfp5 = 0
            wfp6 = 0
            att_wm_sgs=0
        elseif (coef_wall==1) then
            allocate(att_mod_par(1:jx,2,kparasta:kparaend))
            allocate(        u_t(1:jx,2,kparasta:kparaend))
            allocate(  utangente(1:jx,2,kparasta:kparaend))
            allocate(punto_wfp3(3,3,jx,kparasta:kparaend))
            allocate(punto_wfp4(3,3,jx,kparasta:kparaend))
            att_mod_par=0
            utangente=1.
            u_t=1.

            if (wfp3==1) then
                att_mod_par(:,1,:)=1 ! first index=1 : face 3
            end if

            if (wfp4==1) then
                att_mod_par(:,2,:)=1 ! second index=2 : face 4
            end if

        end if

        !      nesting
        termina=.false.
        !

        !-----------------------------------------------------------------------
        ! sett kparasta and kparaend for each multigrid level
        ! k parallel start multigrid -> kpstamg
        ! k parallel end multigrid   -> kpstamg
        do n=1,nlevel
            ncolprocmg(n) = jzc(n)/nproc
            kpstamg(n) = myid*ncolprocmg(n)+1
            kpendmg(n) = (myid+1)*ncolprocmg(n)
            !        only myid=0 has kfacepstamg = 0
            kfacepstamg(n) = kpstamg(n)
            if (myid==0)kfacepstamg(n) = kpstamg(n) -1
        end do
        !
        !-----------------------------------------------------------------------
        ! plane comunication for periodicity in k

        do kk=1,1-kp
            if (myid==nproc-1) then
                call MPI_SSEND(u(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req1,istatus,ierr)
            elseif (myid==0) then
                call MPI_RECV(u_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req2,istatus,ierr)
            endif
            if (myid==0) then
                call MPI_SSEND(u(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req3,istatus,ierr)
            endif
            if (myid==nproc-1) then
                call MPI_RECV(u_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req4,istatus,ierr)
            endif

            if (myid==nproc-1) then
                call MPI_SSEND(v(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req1,istatus,ierr)
            elseif (myid==0) then
                call MPI_RECV(v_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req2,istatus,ierr)
            endif
            if (myid==0) then
                call MPI_SSEND(v(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req3,istatus,ierr)
            endif
            if (myid==nproc-1) then
                call MPI_RECV(v_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req4,istatus,ierr)
            endif


            if (myid==nproc-1) then
                call MPI_SSEND(w(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req1,istatus,ierr)
            elseif (myid==0) then
                call MPI_RECV(w_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req2,istatus,ierr)
            endif
            if (myid==0) then
                call MPI_SSEND(w(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req3,istatus,ierr)
            endif
            if (myid==nproc-1) then
                call MPI_RECV(w_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req4,istatus,ierr)
            endif


            !..............................................................................
            allocate(rho(0:jx+1,0:jy+1,kparasta-deepl:kparaend+deepr))
            if (myid==0) then
                allocate(rho_piano(0:jx+1,0:jy+1,n3:n3))
            elseif (myid==nproc-1) then
                allocate(rho_piano(0:jx+1,0:jy+1,1:1))
            end if

            do isc=1,nscal

                rho  = 0.
                if (myid==0 .or. myid==nproc-1)rho_piano=0.

                do k=kparasta-deepl,kparaend+deepr
                    do j=0,jy+1
                        do i=0,jx+1
                            rho(i,j,k)=rhov(isc,i,j,k)
                        end do
                    end do
                end do


                if (myid==nproc-1) then
                    call MPI_SSEND(rho(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
                !            call MPI_WAIT(req1,istatus,ierr)
                elseif (myid==0) then
                    call MPI_RECV(rho_piano(0,0,jz),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
                !            call MPI_WAIT(req2,istatus,ierr)
                endif
                if (myid==0) then
                    call MPI_SSEND(rho(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
                !            call MPI_WAIT(req3,istatus,ierr)
                endif
                if (myid==nproc-1) then
                    call MPI_RECV(rho_piano(0,0,1),(jx+2)*(jy+2),MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
                !            call MPI_WAIT(req4,istatus,ierr)
                endif


                if (myid==0) then
                    do j=0,jy+1
                        do i=0,jx+1
                            rhov_piano(isc,i,j,jz)=rho_piano(i,j,jz)
                        end do
                    end do
                elseif (myid==0) then
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
            elseif (myid==nproc-1) then
                deallocate(rho_piano)
            end if


        enddo
        !

        ! THIS IS THE LEVEL SET
        !call iniz_levelset()

        count_print_time = 0
        ti_start = ti



    end subroutine les_initialize

    subroutine les_core() bind ( C, name="les_core" )

        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !***********************************************************************
        !     CYCLE
        !***********************************************************************
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------

        !
        startitertime=MPI_WTIME()

        !-----------------------------------------------------------------------
        !++++++++++  WRITE OUTPUT  +++++++++++++++++++++++++++++++++++++++++++++
        !-----------------------------------------------------------------------
        call output_step(count_print_time,dbbx,tipo)

        !
        !-----------------------------------------------------------------------
        ! boundary conditions on du, dv, dw
        !
        call condi1()
        call condi2()
        call condi3()

        ! boundary conditions for periodicity
        !call contourp_se()
        !call contour()

        !
        ! average on fluxes in periodicity direction
        !
        do ii=1,1-ip
            do k=kparasta,kparaend
                do j=1,jy
                    uc(0,j,k)=.5*(uc(0,j,k)+uc(jx,j,k))
                    uc(jx,j,k)=uc(0,j,k)
                end do
            end do
        end do

        do jj=1,1-jp
            do k=kparasta,kparaend
                do i=1,jx
                    vc(i,0,k)=.5*(vc(i,0,k)+vc(i,jy,k))
                    vc(i,jy,k)=vc(i,0,k)
                end do
            end do
        end do
      
        ! send wc(i,j,jz) to P0 in wc_piano of myid=0

        do kk=1,1-kp
            if (myid==nproc-1) then

                call MPI_SSEND(wc(1,1,jz),jx*jy,MPI_REAL_SD,0,1001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req1,istatus,ierr)

            elseif (myid==0) then

                call MPI_RECV(wc_piano(1,1,jz),jx*jy,MPI_REAL_SD,nproc-1,1001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req2,istatus,ierr)

            endif
            !
            if (myid==0) then
                do i=1,jx
                    do j=1,jy
                        wc(i,j,0)=.5*(wc(i,j,0)+wc_piano(i,j,jz))
                    end do
                end do

                call MPI_SSEND(wc(1,1,0),jx*jy,MPI_REAL_SD,nproc-1,2001,MPI_COMM_WORLD,ierr)
            !            call MPI_WAIT(req3,istatus,ierr)
            endif

            if (myid==nproc-1) then
                call MPI_RECV(wc(1,1,jz),jx*jy,MPI_REAL_SD,0,2001,MPI_COMM_WORLD,status,ierr)
            !            call MPI_WAIT(req4,istatus,ierr)
            endif
        enddo

        !-----------------------------------------------------------------------
        !     compute the utau for wall model
        if (coef_wall==1 .and. potenziale==0) then
            call wall_function_bodyfitted(ktime,niter,tipo,i_rest)
        end if

        !-----------------------------------------------------------------------


        ! AAA this external loop is a temporary solution to have the eddy diffusivity computed
        ! for each scalar without reynolds analogy


        if (nsgs==1) then
            !       compute the eddy viscosity and diffusivity with smagorinksy
            !       model with fixed constant
            startturbotime=MPI_WTIME()
            call turbo_statico(ktime,i_print,in_dx1,in_sn1,&
                in_sp1,in_st1,in_av1,in_in1,isotropo,kpstamg,kpendmg)
            endturbotime=MPI_WTIME()
            if (myid==0) then
                print*,myid,'trb static',endturbotime-startturbotime
            endif
        elseif (nsgs==2) then

            !       compute the eddy viscosity and diffusivity with smagorinksy
            !       model with dynamic procedure for the constant
            startturbotime=MPI_WTIME()
            call turbo_dinamico(inmod,inmodrho,ktime,i_print,i_rest,in_dx1,in_sn1,in_sp1, &
                in_st1,in_av1,in_in1,isotropo,kpstamg,kpendmg)
            endturbotime=MPI_WTIME()
            if (myid==0) then
                print*,myid,'trb dynamic',endturbotime-startturbotime
            endif

        elseif (nsgs==3) then
            !       compute the eddy viscosity and diffusivity with smagorinksy
            !       model with lagrangian procedure for the constant
            startturbotime=MPI_WTIME()
            call turbo_lagrangian(inmod,inmodrho,ktime,i_print,i_rest,in_dx1, &
                in_sn1,in_sp1,in_st1,in_av1,in_in1,isotropo,kpstamg,kpendmg)
            endturbotime=MPI_WTIME()
            if (myid==0) then
                print*,myid,'trb dynamic',endturbotime-startturbotime
            endif
        elseif (nsgs==0) then
            !
            !       if DNS:
            do k=kparasta-1,kparaend+1 !0,jz+1
                do j=0,jy+1
                    do i=0,jx+1
                        !          eddy viscosity
                        annit(i,j,k) =1./re
                        annitV(i,j,k)=1./re
       
                        !          eddy diffusivity
                        do isc=1,nscal
                            akapt(isc,i,j,k) =1./re/pran(isc)
                            akaptV(isc,i,j,k)=1./re/pran(isc)
                        end do
                    enddo
                enddo
            enddo
            if (myid==0 .or. myid==nproc-1) then
                annit_piano  = 1./re
                annitV_piano = 1./re
            end if

            if (myid==0 .or. myid==nproc-1) then
                do isc=1,nscal
                    akapt_piano  = 1./re/pran(isc)
                    akaptV_piano = 1./re/pran(isc)
                end do
            end if
            inmod    = 0
            inmodrho = 0

        endif

        !-----------------------------------------------------------------------
        !      IBM CORRECTION

        if ( bodyforce==1 .and. coef_wall >=1 .and. ktime ==1 .and. i_rest /=0 ) then
          
            ! call correggi_ib (...)
            !correggo_rho = 0
            !correggo_delu = 0
            ipressione_ibm = 0
            call correggi_ib(ktime,tipo)
     
     
        end if

        if (bodyforce==1 .and. coef_wall >=1 .and. i_rest/=0) then

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


                !          if (ktime ==1 .and. i_rest/=0) then
                !           fornisco un primo valore per la ustar
                !           calcolo la velocita' tangente all' IB
                !            call vel_tangente(i,j,k,MN,MP,proiezioni,
                !     >      u(i,j,k),v(i,j,k),w(i,j,k),vtan,alfa,l)
                !            call wernerwengle(l,MN,vtan,dist_ib_parete,
                !     >                                dist_pp_ib,ustar)
                !          end if

                l_x =abs( &
                    .25*(x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
                    +x(i  ,j-1,k-1)+x(i  ,j-1,k  )) &
                    -.25*(x(i-1,j  ,k  )+x(i-1,j  ,k-1) &
                    +x(i-1,j-1,k-1)+x(i-1,j-1,k  )))
     
        
                l_y =abs( &
                    .25*(y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
                    +y(i-1,j  ,k-1)+y(i-1,j  ,k  )) &
                    -.25*(y(i  ,j-1,k  )+y(i  ,j-1,k-1) &
                    +y(i-1,j-1,k-1)+y(i-1,j-1,k  )))

                l_z =abs( &
                    .25*(z(i  ,j  ,k  )+z(i  ,j-1,k  ) &
                    +z(i-1,j-1,k  )+z(i-1,j  ,k  )) &
                    -.25*(z(i  ,j  ,k-1)+z(i  ,j-1,k-1) &
                    +z(i-1,j-1,k-1)+z(i-1,j  ,k-1)))
     
                l_f = min(l_x,l_y,l_z)

                l_f = dist_ib_parete(l) + 0.5 * l_f


                if (ustar(l)>0) then
                    coef_annit = l_f / dist_ib_parete(l)
                else
                    coef_annit = 0.
                end if

                if (dist_ib_parete(l)*re*ustar(l) <11) then
                    coef_annit = 0.
                end if


                annit(i0,j0,k0)=coef_annit*.41*ustar(l)*dist_ib_parete(l)

                if (annit(i0,j0,k0)<1./re) then
                    annit(i0,j0,k0)=1./re
                end if


                annitV(i0,j0,k0)=annit(i0,j0,k0)

            end do

        end if !end bodyforce = 1



        !-----------------------------------------------------------------------
        !  nesting
        !  set for ti the value of the first POM time
        !      if (ktime==1.and.i_rest==3) then
        !      ti=ti_pom_old
        !      end if

        ! compute dt at the first iteration
        if (ktime == 1) then
            dt_start = dt
            call courant(ind_cou,cou,kparasta,kparaend,nproc,myid,espl,i_rest,ktime)

        end if
      
        ! update time
        ti=ti+dt
        if (myid==0) then
            write(*,*)'----> ktime, tempo  ',ktime,ti,'----------------------'
            !write(*,*)'bbx, dbbx  ',bbx,dbbx
            !write(*,*)'bby, dbby  ',bby,dbby
        endif
        !-----------------------------------------------------------------------

        if (i_rest==3 .and. potenziale==0) then
            call nesting(bodyforce,ti,area1,area2,area5,area6,kparasta,kparaend,myid,nproc,freesurface)

        elseif (i_rest==3 .and. potenziale==1) then
            call redistribuzione(bodyforce,area1,area2,area5,area6,kparasta,kparaend,myid,nproc)
        end if
 
        !-----------------------------------------------------------------------
        !***********************************************************************
        !-----------------------------------------------------------------------
        ! VISCOSITY COMUNICATION

        !     send to left
        if (myid/=0) then
            call MPI_SSEND(annit(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,tagls,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (req1,status,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(annit(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (req2,status,ierror)
        end if

        !     send to right
        if (myid/=nproc-1) then
            call MPI_SSEND(annit(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (req5,status,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(annit(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (req6,status,ierror)
        end if

        !-----------------------------------------------------------------------
        !     send to left
        if (myid/=0) then
            call MPI_SSEND(annitV(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,tagls,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (req1,status,ierror)
        end if
        if (myid/=nproc-1) then
            call MPI_RECV(annitV(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (req2,status,ierror)
        end if

        !     send to right
        if (myid/=nproc-1) then
            call MPI_SSEND(annitV(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (req5,status,ierror)
        end if
        if (myid/=0) then
            call MPI_RECV(annitV(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        !      call MPI_WAIT (req6,status,ierror)
        end if

        !
        !-----------------------------------------------------------------------
        !     read wind file

        if (windyes==1) then
            call leggivento(kparasta,kparaend,myid,nproc)
            if (langyes==1) then
                !          langmuir circulation
                call vorticitag(myid,nproc,kparasta,kparaend)
                call drift(myid,nproc,kparasta,kparaend)
            end if
        end if
        !
        !-----------------------------------------------------------------------
        ! compute bodyforce acting on the fluid
        !
        call fmassa(dbbx,dbby,dbbz,rich,bbx,bby,bbz,bcsi,beta,bzet, &
            myid,nproc,kparasta,kparaend,ti,langyes,wavebk,ktime)
        !
        starteqstime=MPI_WTIME()

        !-----------------------------------------------------------------------
        ! nesting: generate disturbance on the inflow
        if (ibb==1) then
            call buffer_bodyforce(ti,ktime,bcsi,beta,bzet)
        end if
 

        !-----------------------------------------------------------------------
        ! SCALAR EQUATION
        !-----------------------------------------------------------------------
        !
        ! if attiva_scal=0, scalar eq. solution is bypassed
        if (attiva_scal==0. .or. potenziale==1) then
      
            if (myid==0) then
                print*,myid,'no density equation'
            endif

        else

            !.....................................................................
            ! Check gauss_random
            if (imoist==1) then
                tm_tot=0
                tm_loc = 0
                dist_mean = 0
                do k=kparasta,kparaend
                    do i=1,jx
                        tm_loc=tm_loc+Fx(i,k)
                    end do
                end do
                ! ora faccio la somma di tutte le Tm_loc
                call MPI_ALLREDUCE(tm_loc,tm_tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
                dist_mean = tm_tot/real(jx*jz)
                if (myid==0) then
                    write(*,*)'Disturbance mean = ', dist_mean
                end if

                do k=kparasta,kparaend
                    do i=1,jx
                        Fx(i,k) = Fx(i,k) - dist_mean
                    end do
                end do

                tm_tot = 0
                tm_loc = 0
                do k=kparasta,kparaend
                    do i=1,jx
                        tm_loc = tm_loc + Fx(i,k)
                    end do
                end do
                ! ora faccio la somma di tutte le Tm_loc
                call MPI_ALLREDUCE(tm_loc,tm_tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
                dist_mean=tm_tot/real(jx*jz)
                if (myid==0) then
                    write(*,*)'Disturbance mean = ', dist_mean
                end if
            end if
            !.......................................................

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

                call flud1(uc,cgra1,rho,akapt,insc,isc,tipo2)  ! expl. term in R11
                call flud2(vc,cgra2,rho,akapt,insc,isc,tipo2)  ! expl. term in R22
                call flud3(wc,cgra3,rho,akapt,insc,isc,tipo2)  ! expl. term in R33

                if (espl==1) then
                    call flucrhoesp(rho,isc,tipo,bodyforce)  ! expl. term Crank-Nicolson

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

                    call flucrho(ti,rho,akaptV,akapt,isc,tipo,bodyforce)  ! expl. term Crank-Nicolson

                    call rhs1_rho(kparasta,kparaend)    ! right hand side scalar eq.


                    if (bodyforce==1) then

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
                    !
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
                    !
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
                            if (isc==1) then
                                rhomax=9.5
                                rhomin=6.5
                            elseif (isc==2) then
                                rhomax=40.00
                                rhomin=0.00
                            endif

                            if (rho(i,j,k) < rhomin) then
                                rho(i,j,k) = rhomin
                            end if
      
                            if (rho(i,j,k) > rhomax) then
                                rho(i,j,k) = rhomax
                            end if
                        !   if (rho(i,j,k) < 0.) then
                        !      rho(i,j,k) = 0.
                        !   end if

                        !   if (rho(i,j,k) > 1.) then
                        !      rho(i,j,k) = 1.
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
                    !       quick
                    if (insc==1) then
                        call MPI_SSEND(rho(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
                    end if
                endif
      
                if (rightpem /= MPI_PROC_NULL) then
                    call MPI_RECV(rho(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
                    !       quick
                    if (insc==1) then
                        call MPI_RECV(rho(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
                    end if
                endif
      
                if (rightpem /= MPI_PROC_NULL) then
                    call MPI_SSEND(rho(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
                endif
                if (leftpem /= MPI_PROC_NULL) then
                    call MPI_RECV(rho(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
                endif

                !     put on rhov the computed values
                do k=kparasta-deepl,kparaend+deepr !0,jz+1
                    do j=0,jy+1
                        do i=0,jx+1
                            rhov(isc,i,j,k)=rho(i,j,k)
                        end do
                    end do
                end do

            end do !cyclo on scalar
      
            deallocate(rho)

        endif  !attiva_scal

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
        call mix_para(inmod,iq1,kparasta,kparaend,rightpe,leftpe,tagls,taglr,tagrs,tagrr, &
            nproc,myid,rightpem,leftpem,ktime,i_print,lagr)


        if (myid==0) then
            kparastam=kparasta+kp
            kparaendm=kparaend
        else if (myid==nproc-1) then
            kparastam=kparasta
            kparaendm=kparaend-kp
        else
            kparastam=kparasta
            kparaendm=kparaend
        endif

        call jord1(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1, &
            kparasta,kparaend,rightpe,leftpe,tagls,taglr,tagrs,tagrr,rightpem,leftpem)

        call flu_turbo(kparasta,kparaend)

        if (langyes==1) then
            call langmuir2(myid,nproc,kparasta,kparaend)
        end if

        call flux1(uc,cgra1,u,insc,tipo,bodyforce)     ! explicit term in F21
        call flux2(vc,cgra2,u,insc,tipo,bodyforce)     ! explicit term in F12
        call flux3(wc,cgra3,u,insc,tipo,bodyforce)     ! explicit term in F13

        if (espl==1) then
            if (windyes==1) then
                call flucnesp(u,visualizzo,tauu_att,tipo,bodyforce) ! expl. term Crank-Nicolson
            else
                call flucnesp(u,visualizzo,tipo,bodyforce) ! expl. term Crank-Nicolson
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
            !      wall model is on in flucn only if wfp3=1
            !    tauwz(i,j,k) = tauw(i,j,k)*u(i,j,k)/(u(i,j,k)+w(i,j,k))
            eseguo34=0
            if (wfp3==1.or.wfp4==1) then
                eseguo34 = 1
            end if
            if (windyes==1) then
                call flucn(u,visualizzo,coef_wall,ti,tipo,bodyforce,tauu_att)
            else
                call flucn(u,visualizzo,coef_wall,ti,tipo,bodyforce)
            end if
            eseguo34=0
            call rhs1_rho(kparasta,kparaend)   ! right hand side of momentum eq.

            !
            !-----------------------------------------------------------------------
            !      approximate factorization
            !
            !      first eq., upload csi
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
            !      first eq., upload eta
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
        !
        !-----------------------------------------------------------------------
        !                         SECOND EQUATION
        !-----------------------------------------------------------------------

        starteq2a=MPI_WTIME()
        ! compute convective and diffusive explicit terms in second eq.
        !
        iq1=2
        call mix_para(inmod,iq1,kparasta,kparaend,rightpe,leftpe,tagls,taglr,tagrs,tagrr, &
            nproc,myid,rightpem,leftpem,ktime,i_print,lagr)

        call jord2(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1, &
            kparasta,kparaend,rightpe,leftpe,tagls,taglr,tagrs,tagrr,rightpem,leftpem)

        call flu_turbo(kparasta,kparaend)

        if (langyes==1) then
            call langmuir2(myid,nproc,kparasta,kparaend)
        end if

        call flux1(uc,cgra1,v,insc,tipo,bodyforce)     ! explicit term in F21
        call flux2(vc,cgra2,v,insc,tipo,bodyforce)     ! explicit term in F22
        call flux3(wc,cgra3,v,insc,tipo,bodyforce)     ! explicit term in F23

        if (espl==1) then
            if (windyes==1) then
                call flucnesp(v,visualizzo,v_att,tipo,bodyforce)
            else
                call flucnesp(v,visualizzo,tipo,bodyforce)
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
                call flucn(v,visualizzo,coef_wall,ti,tipo,bodyforce,v_att)  ! expl. term Crank-Nicolson
            else
                call flucn(v,visualizzo,coef_wall,ti,tipo,bodyforce)  ! expl. term Crank-Nicolson
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
      
        !     compute convective and diffusive explicit terms in third eq.
        !
        iq1=3
        call mix_para(inmod,iq1,kparasta,kparaend,rightpe,leftpe,tagls,taglr,tagrs,tagrr, &
            nproc,myid,rightpem,leftpem,ktime,i_print,lagr)
     
        call jord3(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kparasta,kparaend,rightpe,leftpe,&
            tagls,taglr,tagrs,tagrr,rightpem,leftpem)

        call flu_turbo(kparasta,kparaend)

        if (langyes==1) then
            call langmuir2(myid,nproc,kparasta,kparaend)
        end if

        call flux1(uc,cgra1,w,insc,tipo,bodyforce)     ! explicit term in F31
        call flux2(vc,cgra2,w,insc,tipo,bodyforce)     ! explicit term in F32
        call flux3(wc,cgra3,w,insc,tipo,bodyforce)     ! explicit term in F33

        if (espl==1) then
            if (windyes==1) then
                call flucnesp(w,visualizzo,tauw_att,tipo,bodyforce) ! expl. term in Crank-Nicolson
            else
                call flucnesp(w,visualizzo,tipo,bodyforce) ! expl. term in Crank-Nicolson
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

            !      wall model is on in flucn only if wfp3=1
            !          tauwz(i,j,k) = tauw(i,j,k)*w(i,j,k)/(u(i,j,k)+w(i,j,k))
            eseguo34=0
            if (wfp3==1.or.wfp4==1) then
                eseguo34 = 2
            end if
            if (windyes==1) then
                call flucn(w,visualizzo,coef_wall,ti,tipo,bodyforce,tauw_att)
            else
                call flucn(w,visualizzo,coef_wall,ti,tipo,bodyforce)
            end if
            eseguo34=0

            call rhs1_rho(kparasta,kparaend)  ! right hand side of momentum eq.

            !-----------------------------------------------------------------------
            !      approximate factorization
            !
            !      third eq., upload csi
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
            !      third eq., upload eta
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
        !
        !-----------------------------------------------------------------------
        !     apply orlansky boundary condition

        startdiv=MPI_WTIME()
        !
        if (lett/=0) then
            call orlansky_generale(ktime)
        end if
        !
        !-----------------------------------------------------------------------
        !     filtering procedure for the flow field
        if (ifiltro == 1 .and. ktime==nfiltro*(ktime/nfiltro)) then
      
            do loopfiltro=1,filtrou
                call filtro(u,xstart,xend,ystart,yend,zstart,zend)
                if (myid==0) then
                    write(*,*)'call filtering for u'
                end if
            end do

            do loopfiltro=1,filtrov
                call filtro(v,xstart,xend,ystart,yend,zstart,zend)
                if (myid==0) then
                    write(*,*)'call filtering for v'
                end if
            end do

      
            do loopfiltro=1,filtrow
                call filtro(w,xstart,xend,ystart,yend,zstart,zend)
                if (myid==0) then
                    write(*,*)'call filtering for w'
                end if
            end do

      
            do loopfiltro=1,filtrorho
                allocate(rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)) !0:n3+1))
                do isc=1,nscal
                    do k=kparasta-1,kparaend+1
                        do j=0,n2+1
                            do i=0,n1+1
                                rho(i,j,k)=rhov(isc,i,j,k)
                            end do
                        end do
                    end do

                    call filtro(rho,xstart,xend,ystart,yend,zstart,zend)
                    if (myid==0) then
                        write(*,*)'call filtering for rho'
                    end if
        
                    do k=kparasta-1,kparaend+1
                        do j=0,n2+1
                            do i=0,n1+1
                                rhov(isc,i,j,k)=rho(i,j,k)
                            end do
                        end do
                    end do
                end do
                deallocate(rho)
            end do
      
            do loopfiltro=1,filtrofi
                call filtro(fi,xstart,xend,ystart,yend,zstart,zend)
                if (myid==0) then
                    write(*,*)'call filtering for fi'
                end if
            end do
      
        end if
        !-----------------------------------------------------------------------
        !     update velocity u^ = u+du from time step n to intermediate one
        !     on wall found with parabolic extrapolation
        call update()
        !
        !-----------------------------------------------------------------------
        !     correction on IBM

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
        !     compute intermediate controvariant component
        !
        if (lett/=0 .and. i_rest/=3) then
            call contra_infout(kparasta,kparaend,rightpe,leftpe, &
                nproc,myid,area1,area2,area3,area4,area5,area6,ktime)
        else
            call contra(kparasta,kparaend,rightpe,leftpe,nproc,myid)
        endif
        !-----------------------------------------------------------------------
        !     check mass on sides 1 and 2
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

        !     check mass on sides 3 and 4
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

        !     check mass on sides 5 and 6
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

        !     from local to total quantities
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

        !     check mass on sides 1 and 2
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

        !     check mass on sides 3 and 4
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

        !     check mass on sides 5 and 6
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

        !     from local to total quantities
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

        !
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
        endif
        !
        enddiv=MPI_WTIME()
        endeqstime=MPI_WTIME()
        !
        !-----------------------------------------------------------------------
        ! COMPUTE COMPUTATIONAL PRESSURE
        !
        startmultitime=MPI_WTIME()

        if (jpos==0) then
        !                call poiss(omega,eps,myid,nproc,kparasta,kparaend, &
        !                    tagls,taglr,tagrs,tagrr,leftpe,rightpe)
        else
            if (ipress_cart==0) then
                call multi(eps,ficycle,nlevel,jxc,jyc,jzc, &
                    kparasta,kparaend,myid,nproc,rightpe,leftpe, &
                    tagls,taglr,tagrs,tagrr,islor, &
                    bodyforce,bodypressure,tipo,deepl,deepr, &
                    ktime,freesurface,ti)
            else
                ! call multi_cart(nlevmultimax,eps,omega,ficycle,nlevel, &
            !                        jxc,jyc,jzc, &
            !                        kparasta,kparaend,myid,nproc, &
            !                        rightpe,leftpe, &
            !                        tagls,taglr,tagrs,tagrr,islor, &
            !                        bodyforce,bodypressure,tipo,deepl,deepr)
            end if
        end if


        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(fi(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
        endif
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(fi(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
        endif
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(fi(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        endif
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(fi(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        endif
                                                                                     
        endmultitime=MPI_WTIME()

        if (myid==0) then
            print*,myid,'mlt ',endmultitime-startmultitime
        endif
        !
        !-----------------------------------------------------------------------
        !     compute pressure gradient
        startgrad=MPI_WTIME()
        call gradie(kparasta,kparaend)
        !
        !-----------------------------------------------------------------------
        !     compute cartesian velocity and controvariant
        call vel_up(kparasta,kparaend,rightpe,leftpe,bodyforce,freesurface)
        !
        !-----------------------------------------------------------------------
        !     correction on IBM
        if (bodyforce==1) then
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
                !     boundary condition

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

        !---------------------------------------------------------------------

        !
        !-----------------------------------------------------------------------
        !     case constant mass
        !     update the pressure gradient in csi
        if (indm==1) then
            !        compute u mean from Uc
            um_loc=0.
            i=jx/2
            do k=kparasta,kparaend
                do j=1,jy
                    if (tipo(i,j,k)/=0) then
                        um_loc=um_loc+uc(i,j,k)/alz/aly
                    end if
                end do
            end do

            !        sum of local um
            call MPI_ALLREDUCE(um_loc,um,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

            !        all procs know um

            if (ktime==1 .and. i_rest/=0) then
                !          leave dbbx unchanged, and um as reference for the next step
                umm1=um
            else
                !          compute delta bbx
                if (myid==0)write(*,*)'myid',myid,'um',um,umm1-um
                dbbx=1.0/dt*0.5*(umm1 - um) !(0.5*umm1+0.5*umm-um) questo era di anna
                dbby=0.
                dbbz=0.
                if (myid==0) then
                    write(444,*)dbbx
                    write(445,*)um
                end if
            end if
        !
        elseif (indm==0) then
            dbbx=0.
            dbby=0.
            dbbz=0.
        endif
        !

        !----------------------------------------------------ant 21Jan TEST
        ! giulia contour va prima di correggi_ib
        !      call contour(kparasta,kparaend,nproc,myid)
        !
        !      call contourp(kparasta,kparaend,nproc,myid)
        !----------------------------------------------------ant 21Jan TEST

        !-----------------------------------------------------------------------
        !     compute max courant or dt
        !
        call courant(ind_cou,cou,kparasta,kparaend,nproc,myid,espl,i_rest,ktime)


        !
        !-----------------------------------------------------------------------
        !     print planes for inflow

        if (inf==1) then
            !       open file
            if (ktime == 1) then
                do ipiani = 1,npiani
                    if (myid<10) then
                        write(idpiano,'(i1)')ipiani
                        write(idproc ,'(i1)')myid
                        filepiano = 'npiano'//idpiano//'processore0'//idproc//'.dat'
                        open(5000+ipiani*100+myid,file=filepiano, status='unknown')
                    else
                        write(idpiano,'(i1)')ipiani
                        write(idproc2,'(i2)')myid
                        filepiano = 'npiano'//idpiano//'processore'//idproc2//'.dat'
                        open(5000+ipiani*100+myid,file=filepiano, status='unknown')
                    end if
                end do
            end if

            do ipiani = 1,npiani
                i = piani(ipiani)
                write(5000+ipiani*100+myid,*)ktime
                do k=kparasta,kparaend
                    do j=1,jy
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

        !-----------------------------------------------------------------------
        !     inflow: read the next planes of data

        if (lett /= 0) then
            !     side 1
            do ipiani=1,1-infout1
                if (npn1 > 1 .and. ktime /= niter) then
                    read(10011,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 1'
                    end if
                    do k=1,jz
                        do j=1,jy
                            read(10011,100)up1(j,k),vp1(j,k),wp1(j,k),(rhovp1(isc,j,k),isc=1,nscal)

                            up1(j,k)=up1(j,k)*real(index_out1(j,k))
                            vp1(j,k)=vp1(j,k)*real(index_out1(j,k))
                            wp1(j,k)=wp1(j,k)*real(index_out1(j,k))
                            do isc=1,nscal
                                rhovp1(isc,j,k)=rhovp1(isc,j,k)*real(index_out1(j,k))
                            end do
                        end do
                    end do
                elseif (ktime == niter) then
                    close(10011)
                end if
            end do
            !     side 2
            do ipiani=1,1-infout2
                if (npn2 > 1 .and. ktime /= niter) then
                    read(10012,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 2'
                    end if
                    do k=1,jz
                        do j=1,jy
                            read(10012,100)up2(j,k),vp2(j,k),wp2(j,k),(rhovp2(isc,j,k),isc=1,nscal)

                            up2(j,k)=up2(j,k)*real(index_out2(j,k))
                            vp2(j,k)=vp2(j,k)*real(index_out2(j,k))
                            wp2(j,k)=wp2(j,k)*real(index_out2(j,k))
                            do isc=1,nscal
                                rhovp2(isc,j,k)=rhovp2(isc,j,k)*real(index_out2(j,k))
                            end do
                        end do
                    end do
                elseif (ktime == niter) then
                    close(10012)
                end if
            end do
            !     side 3
            do ipiani=1,1-infout3
                if (npn3 > 1 .and. ktime /= niter) then
                    read(10013,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 3'
                    end if
                    do k=1,jz
                        do i=1,jx
                            read(10013,100)up3(i,k),vp3(i,k),wp3(i,k),(rhovp3(isc,i,k),isc=1,nscal)

                            up3(i,k)=up3(i,k)*real(index_out3(i,k))
                            vp3(i,k)=vp3(i,k)*real(index_out3(i,k))
                            wp3(i,k)=wp3(i,k)*real(index_out3(i,k))
                            do isc=1,nscal
                                rhovp3(isc,i,k)=rhovp3(isc,i,k)*real(index_out3(i,k))
                            end do

                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10013)
                end if
            end do
            !     side 4
            do ipiani=1,1-infout4
                if (npn4 > 1 .and. ktime /= niter) then
                    read(10014,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 4'
                    end if
                    do k=1,jz
                        do i=1,jx
                            read(10014,100)up4(i,k),vp4(i,k),wp4(i,k),(rhovp4(isc,i,k),isc=1,nscal)

                            up4(i,k)=up4(i,k)*real(index_out4(i,k))
                            vp4(i,k)=vp4(i,k)*real(index_out4(i,k))
                            wp4(i,k)=wp4(i,k)*real(index_out4(i,k))
                            do isc=1,nscal
                                rhovp4(isc,i,k)=rhovp4(isc,i,k)*real(index_out4(i,k))
                            end do

                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10014)
                end if
            end do
            !     side 5
            do ipiani=1,1-infout5
                if (npn5 > 1 .and. ktime /= niter) then
                    read(10015,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 5'
                    end if
                    do j=1,jy
                        do i=1,jx
                            read(10015,100)up5(i,j),vp5(i,j),wp5(i,j),(rhovp5(isc,i,j),isc=1,nscal)

                            up5(i,j)=up5(i,j)*real(index_out5(i,j))
                            vp5(i,j)=vp5(i,j)*real(index_out5(i,j))
                            wp5(i,j)=wp5(i,j)*real(index_out5(i,j))
                            do isc=1,nscal
                                rhovp5(isc,i,j)=rhovp5(isc,i,j)*real(index_out5(i,j))
                            end do

                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10015)
                end if
            end do
            !     side 6
            do ipiani=1,1-infout6
                if (npn6 > 1 .and. ktime /= niter) then
                    read(10016,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 6'
                    end if
                    do j=1,jy
                        do i=1,jx
                            read(10016,100)up6(i,j),vp6(i,j),wp6(i,j),(rhovp6(isc,i,j),isc=1,nscal)

                            up6(i,j)=up6(i,j)*real(index_out6(i,j))
                            vp6(i,j)=vp6(i,j)*real(index_out6(i,j))
                            wp6(i,j)=wp6(i,j)*real(index_out6(i,j))
                            do isc=1,nscal
                                rhovp6(isc,i,j)=rhovp6(isc,i,j)*real(index_out6(i,j))
                            end do
                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10016)
                end if
            end do

100         format(10e13.5)

        end if


            ! THIS IS THE LEVEL SET
            !call Extension_velocities()


        !-----------------------------------------------------------------------
        ! Compute the mean temperature Theta and Humidity Q at the surface
        if (imoist==1) then
            tm_loc = 0.0
            qm_loc = 0.0
            tm_tot = 0.0
            qm_tot = 0.0

            do k=kparasta,kparaend
                do i=1,jx
                    tm_loc=tm_loc+rhov(1,i,1,k)
                    qm_loc=qm_loc+rhov(2,i,1,k)
                end do
            end do
            !     sum on Tm_loc
            call MPI_ALLREDUCE(tm_loc,tm_tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            tm_tot = tm_tot/real(jx*jz)
            !     sum on qm_loc
            call MPI_ALLREDUCE(qm_loc,qm_tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            qm_tot = qm_tot/real(jx*jz)
            if (myid==0) then
                write(*,*)'Tm ',tm_tot,' Qm ',qm_tot
            end if

        end if

        ! update averages and Reynolds stresses
        if (i_medietempo) then
            call update_medie()
        end if

        !-----------------------------------------------------------------------
        !     pass the ghost cells for u,v,w between procs
        !
        !     send u
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(u(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !        quick
            if (insc==1) then
                call MPI_SSEND(u(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            end if
        endif
      
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(u(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            !        quick
            if (insc==1) then
                call MPI_RECV(u(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            end if
        endif
      
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(u(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        endif
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(u(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        endif
                                                                                     
        if (leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            !      call MPI_WAIT(req4,istatus,ierr)
            if (insc==1) then
            !             call MPI_WAIT(req5,istatus,ierr)
            end if
        endif
        if (rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            !      call MPI_WAIT(req3,istatus,ierr)
            if (insc==1) then
            !            call MPI_WAIT(req6,istatus,ierr)
            end if
        endif

        !     send v
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(v(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !        quick
            if (insc==1) then
                call MPI_SSEND(v(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            end if
        endif
      
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(v(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            !        quick
            if (insc==1) then
                call MPI_RECV(v(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            end if
        endif
      
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(v(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        endif
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(v(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        endif
                                                                                     
        if (leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            !      call MPI_WAIT(req4,istatus,ierr)
            if (insc==1) then
            !            call MPI_WAIT(req5,istatus,ierr)
            end if
        endif
        if (rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            !      call MPI_WAIT(req3,istatus,ierr)
            if (insc==1) then
            !            call MPI_WAIT(req6,istatus,ierr)
            end if
        endif

        !     send w
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(w(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            !        quick
            if (insc==1) then
                call MPI_SSEND(w(0,0,kparasta+1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
            end if
        endif
      
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(w(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            !        quick
            if (insc==1) then
                call MPI_RECV(w(0,0,kparaend+2),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
            end if
        endif
      
        if (rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(w(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        endif
        if (leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(w(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        endif
                                                                                     
        if (leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            !      call MPI_WAIT(req4,istatus,ierr)
            if (insc==1) then
            !           call MPI_WAIT(req5,istatus,ierr)
            end if
        endif
        if (rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            !      call MPI_WAIT(req3,istatus,ierr)
            if (insc==1) then
            !            call MPI_WAIT(req6,istatus,ierr)
            end if
        endif

        !-----------------------------------------------------------------------
        !     COMPUTE THE DIVERGENCE
        !-----------------------------------------------------------------------
        call check_divergence(tipo)

        endgrad=MPI_WTIME()
        !
        !-----------------------------------------------------------------------
        !
        mpitimestart=MPI_WTIME()

        mpitimeend=MPI_WTIME()



        !-----------------------------------------------------------------------


        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        enditertime=MPI_WTIME()

        if (myid == 0) then
            print '(i2,a,e15.8)',myid,'itr ',(enditertime-startitertime)/ticks
            write(*,*)'                   '
        endif

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     END CYCLE
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    end subroutine les_core

    subroutine les_finalize() bind ( C, name="les_finalize" )

        call output_finalize()

        endtime=MPI_WTIME()
        resolution=MPI_WTICK()
        do iproc=0,nproc-1
            if (myid==iproc) then
                print*,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
                print*,myid,'elapsed=',endtime-starttime,'resolution=',resolution
                print*,'------------------------------'
            endif
        enddo

        !-----------------------------------------------------------------------
        !***********************************************************************
        !-----------------------------------------------------------------------

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

        if ((wavebk==1.and.windyes==1) .or. imoist==1) then
            deallocate(Fx,Fz)
        end if

        deallocate(vortx,vorty,vortz)
        deallocate(u_drift,w_drift)

        deallocate(ucs,wcs)
      
        deallocate(pran,prsc)
      
        if (imoist==1) then
            deallocate(tpotm)
            deallocate(qm)
            deallocate(tauw3nest)
            deallocate(scalarflux3nest)
        !   deallocate(heatflux3nest)
        !   deallocate(vapflux3nest)
        end if

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
      
        !      do iproc=0,nproc-1
        write (*,*)myid, 'elapsed_time= ',elapsed_time,' seconds'
        write (*,*)myid, 'elapsed_cput= ',elapsed_cput,' seconds'
    !      enddo
  
    !ccc      close(1000)

    !...carlo       call deallocate_grids()
    end subroutine les_finalize

    subroutine init_parallel(kgridparasta,kgridparaend,nlevmultimax)

        use scala3
        use tipologia
        use mpi
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
        !     depending on the setting in scala3.h MPI_REAL_SD assumes
        !     the type MPI_REAL4 or MPI_REAL8
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
        endif

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
        endif
        if (rightpe==nproc) then
            rightpe=0
        endif

        do iproc=0,nproc-1
            if (myid==iproc) then
                write(*,*)'PE: ',myid,'kparasta= ', kparasta,' kparaend= ',kparaend
                write(*,*)'right ',rightpe,'left ',leftpe
                write(info_run_file,*)'PE: ',myid,'kparasta= ', kparasta,' kparaend= ',kparaend
                write(info_run_file,*)'right ',rightpe,'left ',leftpe
            endif
        enddo


        tagls=100+myid        !tag send to leftpe
        taglr=110+myid-1      !tag recv from leftpe
        tagrs=110+myid        !tag send to rightpe
        tagrr=100+myid+1      !tag recv from rightpe

        if (myid==0) then
            tagls=tagls+nproc
            taglr=taglr+nproc
        endif

        if (myid==0) then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if (myid==nproc-1) then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else if ((myid/=0).and.(myid/=nproc-1)) then
            leftpem=leftpe
            rightpem=rightpe
        endif

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
        if (bodyforce == 1) then
            deepl = 2
            deepr = 2
            if (myid==0)       deepl=1
            if (myid==nproc-1) deepr=1
        end if

        !     deep grid right -> deepgr
        !     deep grid left -> deepgl
        !     values in mysending
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
        !     check on conflicts
        if (imoist==1 .and. nscal < 2) then
            if (myid==0) then
                write(*,*)'there is a conflict between the moisture\n'// &
                    'procedure and the number of scalars,\n'// &
                    'be sure nscal>=2 in scala3.h'
                write(info_run_file,*)'there is a conflict between the moisture\n'// &
                    'procedure and the number of scalars,\n'// &
                    'be sure nscal>=2 in scala3.h'
            end if

            stop
        end if

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

    subroutine read_grid()

        real coor_x,coor_y,coor_z
        real,allocatable :: xx1(:,:,:),yy1(:,:,:),zz1(:,:,:)
        real,allocatable :: xx2(:,:,:),yy2(:,:,:),zz2(:,:,:)

        !-----------------------------------------------------------------------
        ! read the grid

        allocate(xx1(-8:n1+8,-8:n2+8,n3-8:n3))
        allocate(yy1(-8:n1+8,-8:n2+8,n3-8:n3))
        allocate(zz1(-8:n1+8,-8:n2+8,n3-8:n3))

        allocate(xx2(-8:n1+8,-8:n2+8,0:8))
        allocate(yy2(-8:n1+8,-8:n2+8,0:8))
        allocate(zz2(-8:n1+8,-8:n2+8,0:8))

        open(12,file='gri3dp_in.dat',status='old')
        !
        read(12,*)alx,aly,alz
        read(12,*)jx,jy,jz

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
        elseif (myid==nproc-1) then
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
                    read(12,*)coor_x
                    read(12,*)coor_y
                    read(12,*)coor_z

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
        if (myid==0) then
            write(*,*)'CHECK: grid read ---> OK'
        end if

        !-----------------------------------------------------------------------
        ! compute centroid position
        if(myid==0)write(*,*)'Compute centroids'
        call compute_centroids(jx,jy,jz,myid,nproc,kparasta,kparaend)

        !-----------------------------------------------------------------------
        ! periodic cells
        !
        call cellep(xx1,yy1,zz1,xx2,yy2,zz2,kgridparasta)
        deallocate(xx1,yy1,zz1,xx2,yy2,zz2)
        !
        !-----------------------------------------------------------------------
        !
        if (myid==0) then
            write(*,*) 'Grid dimension: ',jx,jy,jz
            write(*,*) 'Array dimension: ',n1,n2,n3

            write(*,*)'CHECK: bodyforce --->',bodyforce
        end if

    end subroutine read_grid

    !subroutine convertstring(instr,)

end module strati





