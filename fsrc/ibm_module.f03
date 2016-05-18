module ibm_module

    ! for immersed boundary method

    use,intrinsic :: iso_c_binding

    use mysettings, only: particles,bodyforce
    use myarrays_velo3, only: fi,rhov,u,v,w
    use wallmodel_module

    ! --------------------------------------------------

    integer(kind=c_int),bind(C) :: num_iter
    !MN,MP matrix dimension to store mesh: points and triangle
    !integer :: MN,MP
    integer :: numero_celle_IB,numero_celle_bloccate,numero_celle_bloccate_real
    integer :: num_ib,num_solide

    integer,allocatable :: indici_CELLE_IB(:,:)
    integer,allocatable :: indici_celle_bloccate(:,:)

    !real,allocatable :: distanze_CELLE_IB(:,:) ! apparently useless

    real,allocatable :: dist_pp_ib(:)
    real,allocatable :: dist_ib_parete(:)
    real,allocatable :: ustar(:)
    integer,allocatable :: caso_ib(:)
    real,allocatable :: proiezioni(:,:)
    real,allocatable :: surfVel(:,:)
      
    ! array for rotation with eulerian angles c$
    ! PROBLEMATIC, left hand convention
    real, allocatable :: rot(:,:,:)
    real, allocatable :: rot_inverse(:,:,:)
      
    ! array for trilinear
    real, allocatable :: tricoef(:,:)
    integer, allocatable :: trind(:,:,:)
      
    ! buffer
    real, allocatable :: sbuff_ibm(:)
    real, allocatable :: rbuff_ibm(:)

    ! to send the ib stencil
    integer :: num_left_snd,num_right_snd
    integer :: num_left_rcv,num_right_rcv
    integer, allocatable :: stencil_left_snd(:,:)
    integer, allocatable :: stencil_left_rcv(:,:)
    integer, allocatable :: stencil_right_snd(:,:)
    integer, allocatable :: stencil_right_rcv(:,:)
    integer, allocatable :: tipo_spedito(:,:,:)

    ! to send the solid cells
    integer :: numsolid_left_snd,numsolid_right_snd
    integer :: numsolid_left_rcv,numsolid_right_rcv
    integer, allocatable :: solid_left_snd(:,:)
    integer, allocatable :: solid_left_rcv(:,:)
    integer, allocatable :: solid_right_snd(:,:)
    integer, allocatable :: solid_right_rcv(:,:)

    integer :: ipressione_ibm

    ! ----------------------------------------------------
    ! For particles
    real,allocatable :: pressure_ib(:,:),shear_ib(:,:),momentum_ib(:,:)
    real,allocatable :: position_ib(:,:)
    integer,allocatable :: index_ib(:,:),indexsize_ib(:)
    ! move this to the subroutine in ricerca
    real,allocatable :: normalVector(:,:)

contains

    subroutine correggi_ib(ktime,tipo)

        ! Roman et al. 2008 Computer and Fluids
        ! Roman et al. 2009 Physics of Fluids
        !
        ! subroutine to correct velocity with Immersed Boundary Method (IBM)
        ! correction to zero for velocity and pressure inside the body
        !
        ! Interpolation on PP point with Taylor series
        !
        ! IB velocity reconstruction with linear or log profile
        !-----------------------------------------------------------------------
        ! used to be:
        ! correggi_ib(ktime,coef_wall,visualizzo,integrale,correggo_rho,correggo_delu,delrhov,tipo,z0,ti,nfinale)
        !
        use mysending
        use scala3
        use period
        !
        use mpi

        implicit none

        !-----------------------------------------------------------------------
        ! input stuff
        integer,intent(in) :: ktime
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        !-----------------------------------------------------------------------
        ! local variables
        integer :: kiter
        integer :: i0,j0,k0
        integer :: i,j,k,l
        integer :: coincide_check !,caso,contatore_star
        logical :: node_is_ok(num_ib)

        ! check for errors
        integer :: ierr
        real :: errore,errore_max,errore_max_loc
        real :: variazione1,variazione2,variazione3

        ! velocity at ib points
        real,dimension(3) :: vel_ib
        real,dimension(3) :: vel_ib_old
        real,allocatable :: vel_ib_array(:,:)
        ! velocity values before ib_procedure
        real,allocatable :: vel_ib_oldarray(:,:)

        !-----------------------------------------------------------------------
        ! check call subroutine

        if (myid==0) then
            write(*,*)'-----------------------------------------------'
            write(*,*)'       IMMERSED BOUNDARY'

            if (ktime==1) then
                write(*,*)'IB proc0',num_ib,'su',numero_celle_IB
                write(*,*)'solide proc0',num_solide,'on',numero_celle_bloccate
            end if
        end if

        allocate(vel_ib_array(num_ib,3))
        allocate(vel_ib_oldarray(num_ib,3))


        ! intialize stuff
        if (allocated(ustar)) then
            deallocate(ustar)
        end if
        allocate(ustar(num_ib))
        ustar(:)=0.0

        ! intialize ib properties specific for particles
        !if (particles) then
            if (allocated(shear_ib)) then
                deallocate(shear_ib)
                deallocate(pressure_ib)
                deallocate(momentum_ib)
                deallocate(caso_ib)
            end if
            allocate(shear_ib(num_ib,3))
            allocate(pressure_ib(num_ib,3))
            allocate(momentum_ib(num_ib,3))
            allocate(caso_ib(num_ib))
            shear_ib(:,:)=0.0
            pressure_ib(:,:)=0.0
            momentum_ib(:,:)=0.0
            caso_ib(:)=5
        !end if

        !-----------------------------------------------------------------------
        !     solid border
        do j=1,jy
            do i=1,jx
                k=kparasta-1
                if (tipo(i,j,k)==0) then
                    u(i,j,k)=0.0
                    v(i,j,k)=0.0
                    w(i,j,k)=0.0
                end if
                k=kparaend+1
                if (tipo(i,j,k)==0) then
                    u(i,j,k)=0.0
                    v(i,j,k)=0.0
                    w(i,j,k)=0.0
                end if
            end do
        end do

        !-----------------------------------------------------------------------
        ! correction to solid for IB without V

        do l=1,num_ib

            ! index ib
            i0=indici_CELLE_IB(l,1)
            j0=indici_CELLE_IB(l,2)
            k0=indici_CELLE_IB(l,3)

            ! index V
            i=indici_CELLE_IB(l,4)
            j=indici_CELLE_IB(l,5)
            k=indici_CELLE_IB(l,6)

            coincide_check=(i0-i)*(i0-i)+(j0-j)*(j0-j)+(k0-k)*(k0-k)

            if (coincide_check==0) then
                u(i0,j0,k0)=0.0
                v(i0,j0,k0)=0.0
                w(i0,j0,k0)=0.0
                node_is_ok(l)=.false.
            else
                node_is_ok(l)=.true.
            end if

        end do

        !-----------------------------------------------------------------------
        ! CORRECTION ON IB
        !
        ! initialize
        do l=1,num_ib

            i0=indici_CELLE_IB(l,1)
            j0=indici_CELLE_IB(l,2)
            k0=indici_CELLE_IB(l,3)

            vel_ib_array(l,1)=u(i0,j0,k0)
            vel_ib_array(l,2)=v(i0,j0,k0)
            vel_ib_array(l,3)=w(i0,j0,k0)

        end do

        ! this is the value before the IBM cycle, used for computing momentum exchange
        vel_ib_oldarray(:,:)=vel_ib_array(:,:)

        !-----------------------------------------------------------------------
        ! iterative procedure, necessary if the stencil has IB nodes
        !
        ! chicco mettere errore adimensionale !!!!
        errore_max=1.1e-7 !-5 !1.1d-7
        kiter=1

        ! do while(errore_max.gt.1.0d-7.and.kiter.le.num_iter)
        do while (kiter<num_iter)

            errore_max=1.0e-7
            errore_max_loc=errore_max
            errore=errore_max

            do l=1,num_ib

                i0=indici_CELLE_IB(l,1)
                j0=indici_CELLE_IB(l,2)
                k0=indici_CELLE_IB(l,3)

                ! get velocity from last iteration step
                u(i0,j0,k0)=vel_ib_array(l,1)
                v(i0,j0,k0)=vel_ib_array(l,2)
                w(i0,j0,k0)=vel_ib_array(l,3)

            end do

            ! exchange planes between procs:
            ! - to compute derivative for Taylor
            ! - to coorect IB values in the iterative procedure
            call passo_ibm()

            do l=1,num_ib

                if (node_is_ok(l)) then !only for IB /= V

                    ! index IB
                    i0=indici_CELLE_IB(l,1)
                    j0=indici_CELLE_IB(l,2)
                    k0=indici_CELLE_IB(l,3)

                    vel_ib_old(:)=vel_ib_oldarray(l,:)

                    errore=0.0

                    ! compute u_ib and ustar
                    call compute_u_ib(kiter,ktime,l,i0,j0,k0,vel_ib_old,vel_ib)

                    ! chicco mettere errore adimensionale !!!!
                    ! error at ib
                    variazione1=vel_ib(1)-u(i0,j0,k0)
                    variazione2=vel_ib(2)-v(i0,j0,k0)
                    variazione3=vel_ib(3)-w(i0,j0,k0)

                    errore=max(abs(variazione1),abs(variazione2),abs(variazione3))

                    ! update velocity at IB
                    vel_ib_array(l,:)=vel_ib(:)

                    ! update error
                    if (errore>=errore_max) then
                        errore_max=errore
                        errore_max_loc=errore
                    end if

                end if

            end do !fine loop su celle ib

            ! ----------------------------------------------------------------

            call MPI_ALLREDUCE(errore_max_loc,errore_max,1,MPI_REAL_SD,MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (myid==0) then
                write(*,*)'***',kiter,'errore ciclo IBM',errore_max,'***'
            end if

            kiter=kiter+1

        end do ! end loop correzione

        do l=1,num_solide

            i=indici_celle_bloccate(l,1)
            j=indici_celle_bloccate(l,2)
            k=indici_celle_bloccate(l,3)

            u(i,j,k)=0.0
            v(i,j,k)=0.0
            w(i,j,k)=0.0

            ! giulia
            if (i==1) then
                u(i-1,j,k)=0.0
                v(i-1,j,k)=0.0
                w(i-1,j,k)=0.0
            else if (i==jx) then
                u(i+1,j,k)=0.0
                v(i+1,j,k)=0.0
                w(i+1,j,k)=0.0
            end if
            if (j==1) then
                u(i,j-1,k)=0.0
                v(i,j-1,k)=0.0
                w(i,j-1,k)=0.0
            else if (j==jy) then
                u(i,j+1,k)=0.0
                v(i,j+1,k)=0.0
                w(i,j+1,k)=0.0
            end if
            if (myid==0.and.k==1) then
                u(i,j,k-1)=0.0
                v(i,j,k-1)=0.0
                w(i,j,k-1)=0.0
            else if (myid==(nproc-1).and.k==jz) then
                u(i,j,k+1)=0.0
                v(i,j,k+1)=0.0
                w(i,j,k+1)=0.0
            end if

        !        fi(i,j,k)=0.
        end do

        !-----------------------------------------------------------------------
        !     solid border
        do j=1,jy
            do i=1,jx
                k = kparasta -1
                if (tipo(i,j,k)==0) then
                    u(i,j,k)=0.
                    v(i,j,k)=0.
                    w(i,j,k)=0.
                end if
                k = kparaend +1
                if (tipo(i,j,k)==0) then
                    u(i,j,k)=0.
                    v(i,j,k)=0.
                    w(i,j,k)=0.
                end if
            end do
        end do
        !
        !-----------------------------------------------------------------------
        !     pass data
        call passo_ibm()



        !-----------------------------------------------------------------------
        deallocate(vel_ib_array)
        deallocate(vel_ib_oldarray)
        !-----------------------------------------------------------------------

        return
    end subroutine correggi_ib

    subroutine compute_u_ib(kiter,ktime,l,i0,j0,k0,vel_ib_old,vel_ib)

        use mysettings, only: coef_wall
        use scala3, only: re,dt

        implicit none

        integer,intent(in) :: l,i0,j0,k0
        integer,intent(in) :: kiter,ktime
        real,dimension(3),intent(in) :: vel_ib_old
        ! output of the function: new IB point velocity
        real,dimension(3),intent(out) :: vel_ib

        !-----------------------------------------------------------------------
        ! IP_PP distance
        real :: distanza,fattore_distanza
        ! pressure values
        real :: f_pp,f_ib,fi_pro
        integer :: caso
        !real :: pro1,pro2,pro3
        real :: local_shear1,local_shear2
        ! velocities at PP - Projection Point
        real,dimension(3) :: vel_pp,vloc_pp
        real :: vtan_pp
        ! velocities at IP - Interface Point
        real,dimension(3) :: vel_ip,vloc_ip
        real :: vtan_ip
        ! IB relative velocities (IB/IP - Immersed-Boundary point relative to Interface point)
        real,dimension(3) :: vel_ib_rel,vloc_ib_rel
        real :: vtan_ib_rel
        ! PP relative velocities (PP/IP - Projection point relative to Interface point)
        real,dimension(3) :: vel_pp_rel,vloc_pp_rel
        real :: vtan_pp_rel
        ! shear velocity
        real :: ustarHere

        ! ---------------------------------------------------------------------

        ! distance IP-PP
        distanza=dist_ib_parete(l)+dist_pp_ib(l)
        fattore_distanza=dist_ib_parete(l)/distanza

        ! velocity at interpolation points (PP, V?)
        vel_pp(1)=trilinear_interpolation(u,l)
        vel_pp(2)=trilinear_interpolation(v,l)
        vel_pp(3)=trilinear_interpolation(w,l)

        ! velocity at IB points (initial, global)
        vel_ib(1)=u(i0,j0,k0)
        vel_ib(2)=v(i0,j0,k0)
        vel_ib(3)=w(i0,j0,k0)

        ! velocity at PP (local)
        vloc_pp(:)=global2local(vel_pp(:),l)
        vtan_pp=sqrt(vloc_pp(1)**2.0+vloc_pp(2)**2.0)

        ! if the surface is moving (i.e. there are particles) the law of the wall does not apply in
        ! absolute terms, but relative to the velocity of the surface (here surfVel)

        ! compute velocities (global)= at IP  (/=0 only if particles are active)
        if (particles) then
            vel_ip(:)=surfVel(l,:)
        else
            vel_ip(:)=(/0.0,0.0,0.0/)
        end if

        ! velocity at IP (local)
        vloc_ip(:)=global2local(vel_ip(:),l)

        vtan_ip=sqrt(vloc_ip(1)**2.0+vloc_ip(2)**2.0)

        ! velocities at PP relative to IP (global)
        vel_pp_rel(:)=vel_pp(:)-vel_ip(:)

        ! velocities at PP relative to IP (local)
        vloc_pp_rel(:)=vloc_pp(:)-vloc_ip(:)
        vtan_pp_rel=sqrt(vloc_pp_rel(1)**2.0+vloc_pp_rel(2)**2.0)

        ! determine case
        ! (0 = no vel; 1 = linear profile; 2 = full law)
        if (coef_wall==0) then
            ! linear profile
            caso=1
        else if ((dist_ib_parete(l)<1.e-9).or.(vtan_pp_rel<1.e-9)) then
            ! values too small: zero velocity
            caso=0
        else
            ! full law of the wall
            caso=2
        end if

        !-----------------------------------------------------------------------
        ! update velocity at IB

        if (caso==0) then
            ! zero velocity
            vel_ib(:)=0.0
            ustarHere=0.0

        else if (caso==1) then
            ! linear profile (in relative form)

            !ustar(l)=sqrt(abs(uip)/(distanza*re)) ! bug???
            !ustar(l)=sqrt(abs(vtan_pp/(distanza*re)))
            ustarHere=ustar_linprofile(distanza,vtan_pp_rel)

            ! linear profile
            vel_ib_rel(:)=fattore_distanza*vel_pp_rel(:)+vel_ip(:)

            vel_ib(:)=vel_ib_rel(:)+vel_ip(:)

        else if (caso==2) then
            ! law of the wall (with respect to wall velocity)

            ! find the shear velocity with the law of the wall
            ! procedure is different depending on whether we have an initial value for ustar
            if (kiter==1 .and. ktime==1) then
                ustarHere=ustar_lawofwall(distanza,vtan_pp_rel)
            else
                ustarHere=ustar_lawofwall(distanza,vtan_pp_rel,ustar(l))
            end if

            ! I know ustar which is fix for the normal
            vtan_ib_rel=speed_lawofwall(dist_ib_parete(l),ustarHere)

            ! log profile

            ! normal and tangential component at IB in the local frame of reference
            ! here they suppose same behavior for tangential and normal velocity (questionable)
            vloc_ib_rel(:)=vloc_pp_rel(:)*(vtan_ib_rel/vtan_pp_rel)

            ! construct the cartesian component in the general frame of reference
            vel_ib_rel(:)=local2global(vloc_ib_rel(:),l)

            vel_ib(:)=vel_ib_rel(:)+vel_ip(:)

        end if ! caso

        ! finally we determine ustar
        ustar(l)=ustarHere

        ! other things necessary for particles
        !if (particles) then

            caso_ib(l)=caso

            ! compute components of shear stress (conserve sign)
            local_shear1=sign((vloc_pp_rel(1)*(ustarHere/vtan_pp_rel))**2,vloc_pp_rel(1))
            local_shear2=sign((vloc_pp_rel(2)*(ustarHere/vtan_pp_rel))**2,vloc_pp_rel(2))

            f_pp=trilinear_interpolation(fi,l)

            f_ib=fi(i0,j0,k0)

            ! this assumes a linear behavior of pressure along the normal line
            fi_pro=((dist_ib_parete(l)+dist_pp_ib(l))*f_ib-f_pp*dist_ib_parete(l))/dist_pp_ib(l)

            ! force on solid objects are made up of two components: pressure and shear (2,3 INVERTED)
            shear_ib(l,1)=local_shear1*rot_inverse(l,1,1)+local_shear2*rot_inverse(l,1,2)
            shear_ib(l,3)=local_shear1*rot_inverse(l,2,1)+local_shear2*rot_inverse(l,2,2)
            shear_ib(l,2)=local_shear1*rot_inverse(l,3,1)+local_shear2*rot_inverse(l,3,2)

            ! pressure has a minus sign because the normal is directed outwards (2,3 INVERTED)
            pressure_ib(l,1)=-1.0*fi_pro*rot_inverse(l,1,3)
            pressure_ib(l,3)=-1.0*fi_pro*rot_inverse(l,2,3)
            pressure_ib(l,2)=-1.0*fi_pro*rot_inverse(l,3,3)

            momentum_ib(l,:)=(vel_ib_old(:)-vel_ib(:))/dt

        !   the following is valid if rotationAle is used!!!!
        ! force on solid objects are made up of two components: pressure and shear
        !   fluidShearForce(i0,j0,k0,1)=local_shear1*rot_inverse(l,1,1) &
        !                +local_shear2*rot_inverse(l,1,2)
        !   fluidShearForce(i0,j0,k0,2)=local_shear1*rot_inverse(l,2,1) &
        !                +local_shear2*rot_inverse(l,2,2)
        !   fluidShearForce(i0,j0,k0,3)=local_shear1*rot_inverse(l,3,1) &
        !                +local_shear2*rot_inverse(l,3,2)
        !   ! pressure has a minus sign because the normal is directed outwards
        !   fluidPressureForce(i0,j0,k0,1)=-1.0*fi_pro*rot_inverse(l,1,3)
        !   fluidPressureForce(i0,j0,k0,2)=-1.0*fi_pro*rot_inverse(l,2,3)
        !   fluidPressureForce(i0,j0,k0,3)=-1.0*fi_pro*rot_inverse(l,3,3)
        !
        !    shear_ib(l,1)=fluidShearForce(i0,j0,k0,1)
        !    shear_ib(l,2)=fluidShearForce(i0,j0,k0,2)
        !    shear_ib(l,3)=fluidShearForce(i0,j0,k0,3)
        !
        !    pressure_ib(l,1)=fluidPressureForce(i0,j0,k0,1)
        !    pressure_ib(l,2)=fluidPressureForce(i0,j0,k0,2)
        !    pressure_ib(l,3)=fluidPressureForce(i0,j0,k0,3)
        !end if

        return

    end subroutine compute_u_ib

    subroutine passo_ibm()

        ! used to be:
        ! passo_ibm(correggo_delu)

        !***********************************************************************
        ! send IB at the border between the procs
        use mysending
        use scala3
        use period
        !
        use mpi


        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,l,isc
        integer iper,kper,jper
        integer ireq1,ireq2,ireq3,ireq4
        integer status(MPI_STATUS_SIZE),ierror
        !integer correggo_delu

        logical icheck

        !-----------------------------------------------------------------------
        icheck = .false.
        !-----------------------------------------------------------------------
        if(icheck)then
            write(900+myid,*)myid,num_left_snd,num_right_snd,num_left_rcv,num_right_rcv
        end if

        !     left send, right recive
        if(num_left_snd .ne.0)then
            allocate(sbuff_ibm( (3+nscal)*num_left_snd ))
            sbuff_ibm = 0.
        end if
        if(num_right_rcv .ne.0)then
            allocate(rbuff_ibm( (3+nscal)*num_right_rcv))
            rbuff_ibm = 0.
        end if

        if(num_left_snd .ne.0)then
            if(myid.ne.0)then
                do l=1,num_left_snd
                    i = stencil_left_snd(l,1)
                    j = stencil_left_snd(l,2)
                    k = stencil_left_snd(l,3)

                    sbuff_ibm(l               ) = u(i,j,k)
                    sbuff_ibm(l+  num_left_snd) = v(i,j,k)
                    sbuff_ibm(l+2*num_left_snd) = w(i,j,k)
                    do isc=1,nscal
                        sbuff_ibm(l+3*num_left_snd &
                            +(isc-1)*num_left_snd) &
                            = rhov(isc,i,j,k)
                    end do
                end do

                call MPI_SEND(sbuff_ibm(1),(3+nscal)*num_left_snd,MPI_REAL_SD, &
                    leftpe,tagls,MPI_COMM_WORLD,ierror)
            !      call MPI_WAIT (ireq1,status,ierror)
            end if
        end if

        if(num_right_rcv.ne.0)then
            if(myid.ne.nproc-1)then
                call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_right_rcv,MPI_REAL_SD, &
                    rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
                !      call MPI_WAIT (ireq2,status,ierror)

                do l=1,num_right_rcv
                    i = stencil_right_rcv(l,1)
                    j = stencil_right_rcv(l,2)
                    k = stencil_right_rcv(l,3)

                    u(i,j,k) = rbuff_ibm(l                )
                    v(i,j,k) = rbuff_ibm(l+  num_right_rcv)
                    w(i,j,k) = rbuff_ibm(l+2*num_right_rcv)
                    do isc=1,nscal
                        rhov(isc,i,j,k)= rbuff_ibm(l+3*num_right_rcv  &
                            +(isc-1)*num_right_rcv)
                    end do
                end do
            end if
        end if

        ! now periodicity
        do kper = 1,1-kp
            if(num_left_snd .ne.0)then
                if(myid.eq.0)then
                    do l=1,num_left_snd
                        i = stencil_left_snd(l,1)
                        j = stencil_left_snd(l,2)
                        k = stencil_left_snd(l,3)

                        sbuff_ibm(l               ) = u(i,j,k)
                        sbuff_ibm(l+  num_left_snd) = v(i,j,k)
                        sbuff_ibm(l+2*num_left_snd) = w(i,j,k)
                        do isc=1,nscal
                            sbuff_ibm(l+3*num_left_snd+(isc-1)*num_left_snd) &
                                = rhov(isc,i,j,k)
                        end do
                    end do

                    call MPI_SEND(sbuff_ibm(1),(3+nscal)*num_left_snd, &
                        MPI_REAL_SD,nproc-1,tagls,MPI_COMM_WORLD,ierror)
                !          call MPI_WAIT (ireq1,status,ierror)
                end if
            end if

            if(num_right_rcv .ne.0)then
                if(myid.eq.nproc-1)then
                    call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_right_rcv, &
                        MPI_REAL_SD,0,tagrr,MPI_COMM_WORLD,status,ierror)
                    !          call MPI_WAIT (ireq2,status,ierror)


                    ! chicco AAA nno puo' andare bene qui con la modifica applicata!
                    ! deve inserire questi dati non in u(i,j,k) ma in u_piano(i,j)
                    ! oppure direttamente su u(i,j,0)

                    do l=1,num_right_rcv
                        i = stencil_right_rcv(l,1)
                        j = stencil_right_rcv(l,2)
                        k = stencil_right_rcv(l,3)

                        u(i,j,k) = rbuff_ibm(l                )
                        v(i,j,k) = rbuff_ibm(l+  num_right_rcv)
                        w(i,j,k) = rbuff_ibm(l+2*num_right_rcv)
                        do isc=1,nscal
                            rhov(isc,i,j,k)  &
                                = rbuff_ibm(l+3*num_right_rcv+(isc-1)*num_right_rcv)
                        end do
                    end do
                end if
            end if
        end do


        if(num_left_snd .ne.0)then
            deallocate(sbuff_ibm)
        end if
        if(num_right_rcv .ne.0)then
            deallocate(rbuff_ibm)
        end if


        !-----------------------------------------------------------------------
        !     right send, left recive
        !
        if(num_right_snd .ne.0)then
            allocate(sbuff_ibm( (3+nscal)*num_right_snd ))
        end if
        if(num_left_rcv.ne.0)then
            allocate(rbuff_ibm( (3+nscal)*num_left_rcv))
        end if

        if(num_right_snd .ne.0)then
            if(myid.ne.nproc-1)then
                do l=1,num_right_snd
                    i = stencil_right_snd(l,1)
                    j = stencil_right_snd(l,2)
                    k = stencil_right_snd(l,3)

                    sbuff_ibm(l                ) = u(i,j,k)
                    sbuff_ibm(l+  num_right_snd) = v(i,j,k)
                    sbuff_ibm(l+2*num_right_snd) = w(i,j,k)
                    do isc=1,nscal
                        sbuff_ibm(l+3*num_right_snd+(isc-1)*num_right_snd) &
                            =rhov(isc,i,j,k)
                    end do
                end do

                call MPI_SEND(sbuff_ibm(1),(3+nscal)*num_right_snd,MPI_REAL_SD, &
                    rightpe,tagrs,MPI_COMM_WORLD,ierror)
            !      call MPI_WAIT (ireq3,status,ierror)
            end if
        end if

        if(num_left_rcv.ne.0)then
            if(myid.ne.0)then
                call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_left_rcv,MPI_REAL_SD, &
                    leftpe,taglr,MPI_COMM_WORLD,status,ierror)
                !      call MPI_WAIT (ireq4,status,ierror)

                do l=1,num_left_rcv
                    i = stencil_left_rcv(l,1)
                    j = stencil_left_rcv(l,2)
                    k = stencil_left_rcv(l,3)

                    u(i,j,k) = rbuff_ibm(l               )
                    v(i,j,k) = rbuff_ibm(l+  num_left_rcv)
                    w(i,j,k) = rbuff_ibm(l+2*num_left_rcv)
                    do isc=1,nscal
                        rhov(isc,i,j,k)  &
                            = rbuff_ibm(l+3*num_left_rcv+(isc-1)*num_left_rcv)
                    end do
                end do
            end if
        end if

        !     now periodicity
        do kper = 1,1-kp
            if(num_right_snd .ne.0)then
                if(myid.eq.nproc-1)then
                    do l=1,num_right_snd
                        i = stencil_right_snd(l,1)
                        j = stencil_right_snd(l,2)
                        k = stencil_right_snd(l,3)

                        sbuff_ibm(l                ) = u(i,j,k)
                        sbuff_ibm(l+  num_right_snd) = v(i,j,k)
                        sbuff_ibm(l+2*num_right_snd) = w(i,j,k)
                        do isc=1,nscal
                            sbuff_ibm(l+3*num_right_snd+(isc-1)*num_right_snd) &
                                =rhov(isc,i,j,k)
                        end do
                    end do

                    call MPI_SEND(sbuff_ibm(1),(3+nscal)*num_right_snd, &
                        MPI_REAL_SD,0,tagrs,MPI_COMM_WORLD,ierror)
                !           call MPI_WAIT (ireq3,status,ierror)
                end if
            end if

            if(num_left_rcv.ne.0)then
                if(myid.eq.0)then
                    call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_left_rcv, &
                        MPI_REAL_SD,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
                    !           call MPI_WAIT (ireq4,status,ierror)

                    do l=1,num_left_rcv
                        i = stencil_left_rcv(l,1)
                        j = stencil_left_rcv(l,2)
                        k = stencil_left_rcv(l,3)

                        u(i,j,k) = rbuff_ibm(l               )
                        v(i,j,k) = rbuff_ibm(l+  num_left_rcv)
                        w(i,j,k) = rbuff_ibm(l+2*num_left_rcv)
                        do isc=1,nscal
                            rhov(isc,i,j,k)  &
                                = rbuff_ibm(l+3*num_left_rcv+(isc-1)*num_left_rcv)
                        end do
                    end do
                end if
            end if
        end do

        if(num_right_snd .ne.0)then
            deallocate(sbuff_ibm)
        end if
        if(num_left_rcv.ne.0)then
            deallocate(rbuff_ibm)
        end if

        !-----------------------------------------------------------------------
        !     i direction, impose periodicity or wall
        !
        !     parabolic extrapolation y=ax2+bx+c
        if(ip==1)then

        !        if(correggo_delu==1)then
        !            do k=kparasta,kparaend
        !                do j=1,jy
        !                    !
        !                    u(0,j,k)=1.875*u(1,j,k) &
        !                        -1.25*u(2,j,k) &
        !                        +.375*u(3,j,k)
        !                    v(0,j,k)=1.875*v(1,j,k) &
        !                        -1.25*v(2,j,k) &
        !                        +.375*v(3,j,k)
        !                    w(0,j,k)=1.875*w(1,j,k) &
        !                        -1.25*w(2,j,k) &
        !                        +.375*w(3,j,k)
        !                    !
        !                    u(jx+1,j,k)=.375*u(jx-2,j,k) &
        !                        -1.25*u(jx-1,j,k) &
        !                        +1.875*u(jx  ,j,k)
        !                    !
        !                    v(jx+1,j,k)=.375*v(jx-2,j,k) &
        !                        -1.25*v(jx-1,j,k) &
        !                        +1.875*v(jx  ,j,k)
        !                    !
        !                    w(jx+1,j,k)=.375*w(jx-2,j,k) &
        !                        -1.25*w(jx-1,j,k) &
        !                        +1.875*w(jx  ,j,k)
        !                !
        !                end do
        !            end do
        !        end if
        !
        else
            !
            !     periodicity
            !
            do k=kparasta,kparaend
                do j=1,jy

                    u(0   ,j,k)=u(jx,j,k)
                    v(0   ,j,k)=v(jx,j,k)
                    w(0   ,j,k)=w(jx,j,k)
                    u(jx+1,j,k)=u(1 ,j,k)
                    v(jx+1,j,k)=v(1 ,j,k)
                    w(jx+1,j,k)=w(1 ,j,k)
                !
                end do
            end do

        end if
        !
        !-----------------------------------------------------------------------------
        !     j direction, impose periodicity or wall
        !
        !     parabolic extrapolation y=ax2+bx+c
        if(jp==1)then
            !
        !        if(correggo_delu==1)then
        !            do k=kparasta,kparaend
        !                do i=1,jx
        !
        !                    u(i,0,k)=1.875*u(i,1,k) &
        !                        -1.25*u(i,2,k) &
        !                        +.375*u(i,3,k)
        !                    v(i,0,k)=1.875*v(i,1,k) &
        !                        -1.25*v(i,2,k) &
        !                        +.375*v(i,3,k)
        !                    w(i,0,k)=1.875*w(i,1,k) &
        !                        -1.25*w(i,2,k) &
        !                        +.375*w(i,3,k)
        !                    !
        !                    u(i,jy+1,k)=.375*u(i,jy-2,k) &
        !                        -1.25*u(i,jy-1,k) &
        !                        +1.875*u(i,jy  ,k)
        !                    !
        !                    v(i,jy+1,k)=.375*v(i,jy-2,k) &
        !                        -1.25*v(i,jy-1,k) &
        !                        +1.875*v(i,jy  ,k)
        !                    !
        !                    w(i,jy+1,k)=.375*w(i,jy-2,k) &
        !                        -1.25*w(i,jy-1,k) &
        !                        +1.875*w(i,jy  ,k)
        !                !
        !                end do
        !            end do
        !        end if

        else
            !
            !     periodicity on eta
            do k=kparasta,kparaend
                do i=1,jx

                    u(i,   0,k)=u(i,jy,k)
                    v(i,   0,k)=v(i,jy,k)
                    w(i,   0,k)=w(i,jy,k)
                    u(i,jy+1,k)=u(i, 1,k)
                    v(i,jy+1,k)=v(i, 1,k)
                    w(i,jy+1,k)=w(i, 1,k)
                !
                end do
            end do

        end if
        !
        !-----------------------------------------------------------------------------
        !    k direction, impose periodicity or wall
        !
        ! extrapolation on front and back sides (5 and 6)
        if(kp==1)then

        !        if(correggo_delu==1)then
        !            if (myid.eq.0) then
        !                do j=1,jy
        !                    do i=1,jx
        !                        !
        !                        u(i,j,0)=1.875*u(i,j,1) &
        !                            -1.25*u(i,j,2) &
        !                            +.375*u(i,j,3)
        !                        v(i,j,0)=1.875*v(i,j,1) &
        !                            -1.25*v(i,j,2) &
        !                            +.375*v(i,j,3)
        !                        w(i,j,0)=1.875*w(i,j,1) &
        !                            -1.25*w(i,j,2) &
        !                            +.375*w(i,j,3)
        !                    end do
        !                end do
        !            endif
        !            !
        !            if (myid.eq.nproc-1) then
        !                do j=1,jy
        !                    do i=1,jx
        !
        !                        u(i,j,jz+1)=.375*u(i,j,jz-2) &
        !                            -1.25*u(i,j,jz-1) &
        !                            +1.875*u(i,j,jz  )
        !                        !
        !                        v(i,j,jz+1)=.375*v(i,j,jz-2) &
        !                            -1.25*v(i,j,jz-1) &
        !                            +1.875*v(i,j,jz  )
        !                        !
        !                        w(i,j,jz+1)=.375*w(i,j,jz-2) &
        !                            -1.25*w(i,j,jz-1) &
        !                            +1.875*w(i,j,jz  )
        !                    !
        !                    enddo
        !                enddo
        !            end if
        !        end if

        end if

        ! the periodicity in k has been already applied with the message passing
        !-----------------------------------------------------------------------------
        if(icheck)then

            k=kparasta
            do j=1,jy
                do i=1,jx
                    write(4500+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do

            k=kparasta-1
            do j=1,jy
                do i=1,jx
                    write(4600+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do

            k=kparaend
            do j=1,jy
                do i=1,jx
                    write(4700+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do

            k=kparaend+1
            do j=1,jy
                do i=1,jx
                    write(4800+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do


        end if
        !-----------------------------------------------------------------------------
        return
    end

    subroutine carico_immb(tipo)
        !***********************************************************************
        !     read ibm input file
        use mysending
        !
        use scala3
        use period
        !
        use mpi

        implicit none

        !-----------------------------------------------------------------------
        integer,intent(out) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        !     array declaration
        integer l,i,j,k,in,jn,kn
        integer iloop,jloop,kloop

        integer status(MPI_STATUS_SIZE),ierror

        integer contatore,num_solide_real
        integer ib_totali,solide_totali

        real a1,a2,a3
        real dist1,dist2

        real rot11,rot12,rot13
        real rot21,rot22,rot23
        real rot31,rot32,rot33

        real irot11,irot12,irot13
        real irot21,irot22,irot23
        real irot31,irot32,irot33

        real tri1,tri2,tri3,tri4
        integer i_tri1,i_tri2,i_tri3,i_tri4
        integer j_tri1,j_tri2,j_tri3,j_tri4
        integer k_tri1,k_tri2,k_tri3,k_tri4

        integer np

        integer :: ierr
        !-----------------------------------------------------------------------
        allocate(tipo_spedito(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        !-----------------------------------------------------------------------
        if(myid.eq.0)then
            WRITE(*,*)' '
            write(*,*)'*****************************************'
            write(*,*)'     LOAD IBM'
        end if
        !-----------------------------------------------------------------------
        !     initialization
        do k=kparasta-deepl,kparaend+deepr !0,jz+1
            do j=0,jy+1
                do i=0,jx+1
                    tipo(i,j,k)=2
                end do
            end do
        end do

        do i=1,num_ib
            do j=1,6
                indici_CELLE_IB(i,j)=0
            end do
        end do

        !   do i=1,MN
        !      do j=1,3
        !         distanze_CELLE_IB(i,j)=0
        !      end do
        !   end do

        do i=1,num_solide
            do j=1,3
                indici_celle_bloccate(i,j)=0
            end do
        end do

        if(myid.eq.0)then
            write(*,*)'matrix ---> OK'
        end if

        !-----------------------------------------------------------------------
        !     read input files for ibm


        !     ib point
        open (20,file='Celle_IB_indici.inp',status='old')
        !open (21,file='Celle_IB_distanze.inp',status='old')
        open (23,file='distanze_interpolazioni.inp',status='old')
        open (24,file='rotazione.inp',status='old')

        open (25,file='trilinear_ibm.inp',status='old')
        open (26,file='trilinear_i.inp',status='old')
        open (27,file='trilinear_j.inp',status='old')
        open (28,file='trilinear_k.inp',status='old')

        read (20,*) numero_celle_IB

        num_ib=0

        do l=1,numero_celle_IB

            read(20,*)i,j,k,in,jn,kn

            !read(21,310)dist_x,dist_y,dist_z


            read(23,320)dist1,dist2,a1,a2,a3

            read(24,350)rot11,rot12,rot13, &
                rot21,rot22,rot23, &
                rot31,rot32,rot33

            read(24,350)irot11,irot12,irot13, &
                irot21,irot22,irot23, &
                irot31,irot32,irot33

            read(25,360)tri1,tri2,tri3,tri4

            read(26,370)i_tri1,i_tri2,i_tri3,i_tri4

            read(27,370)j_tri1,j_tri2,j_tri3,j_tri4

            read(28,370)k_tri1,k_tri2,k_tri3,k_tri4

            if((k.ge.kparasta-deepl).and.(k.le.kparaend+deepr))tipo(i,j,k)=1

            if((k.ge.kparasta).and.(k.le.kparaend))then

                tipo(i,j,k)=1

                num_ib=num_ib+1

                indici_CELLE_IB(num_ib,1)=i
                indici_CELLE_IB(num_ib,2)=j
                indici_CELLE_IB(num_ib,3)=k
                indici_CELLE_IB(num_ib,4)=in
                indici_CELLE_IB(num_ib,5)=jn
                indici_CELLE_IB(num_ib,6)=kn

                !         distanze_CELLE_IB(num_ib,1)=dist_x
                !         distanze_CELLE_IB(num_ib,2)=dist_y
                !         distanze_CELLE_IB(num_ib,3)=dist_z

                dist_pp_ib(num_ib)=dist1
                dist_ib_parete(num_ib)=dist2

                proiezioni(num_ib,1)=a1
                proiezioni(num_ib,2)=a2
                proiezioni(num_ib,3)=a3


                !         rotation matrix (to construct tangential and normal velocity)
                rot(num_ib,1,1)=rot11
                rot(num_ib,1,2)=rot12
                rot(num_ib,1,3)=rot13

                rot(num_ib,2,1)=rot21
                rot(num_ib,2,2)=rot22
                rot(num_ib,2,3)=rot23

                rot(num_ib,3,1)=rot31
                rot(num_ib,3,2)=rot32
                rot(num_ib,3,3)=rot33


                !         inverse rotation matrix
                rot_inverse(num_ib,1,1)=irot11
                rot_inverse(num_ib,1,2)=irot12
                rot_inverse(num_ib,1,3)=irot13

                rot_inverse(num_ib,2,1)=irot21
                rot_inverse(num_ib,2,2)=irot22
                rot_inverse(num_ib,2,3)=irot23

                rot_inverse(num_ib,3,1)=irot31
                rot_inverse(num_ib,3,2)=irot32
                rot_inverse(num_ib,3,3)=irot33

                !         trilinear coefficent
                tricoef(num_ib,1)=tri1
                tricoef(num_ib,2)=tri2
                tricoef(num_ib,3)=tri3
                tricoef(num_ib,4)=tri4


                !         trilinear index
                trind(num_ib,1,1)=i_tri1
                trind(num_ib,2,1)=i_tri2
                trind(num_ib,3,1)=i_tri3
                trind(num_ib,4,1)=i_tri4


                trind(num_ib,1,2)=j_tri1
                trind(num_ib,2,2)=j_tri2
                trind(num_ib,3,2)=j_tri3
                trind(num_ib,4,2)=j_tri4


                trind(num_ib,1,3)=k_tri1
                trind(num_ib,2,3)=k_tri2
                trind(num_ib,3,3)=k_tri3
                trind(num_ib,4,3)=k_tri4


            end if
        end do

        close(20)
        close(21)
        close(23)
        close(24)

        close(25)
        close(26)
        close(27)
        close(28)



300     format(6(i4,1x))
310     format(3e15.8)
320     format (5e15.8)
350     format(9e15.8)
360     format(4e15.8)     ! tri_ibm
370     format(4(i8,1x))

        write(*,*)myid,', number of ib points: ',num_ib !200+myid
        !.......................................................................


        !     solid points
        open (22,file='Celle_Bloccate_Indici.inp',status='old')

        read (22,*) numero_celle_bloccate

        num_solide=0
        contatore=0
        do l=1,numero_celle_bloccate

            read(22,*)i,j,k

            if((k.ge.kparasta-deepl).and.(k.le.kparaend+deepr))then

                tipo(i,j,k)=0

                num_solide=num_solide+1

                indici_celle_bloccate(num_solide,1)=i
                indici_celle_bloccate(num_solide,2)=j
                indici_celle_bloccate(num_solide,3)=k

                if(k.gt.kparaend .or. k.lt.kparasta)then
                    contatore=contatore+1
                end if
            end if
        end do

        close (22)
        !      if(myid.eq.0)then
        !      write(*,*)'read solid cells ---> OK'
        !      end if
        write(*,*)'border solid',contatore

        num_solide_real=num_solide-contatore !without border cells

        write(*,*)myid,'number of solid cells: ',num_solide

        !-----------------------------------------------------------------------
        !     check number of ibm for each proc

        call MPI_REDUCE(num_ib,ib_totali,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
        call MPI_REDUCE(num_solide_real,solide_totali,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(myid.eq.0)then
            write(*,*)'--------------------------------------------'
            if(ib_totali.eq.numero_celle_IB)then
                write(*,*)'check IB for each proc --> OK'
            else
                write(*,*)'check IB for each proc --> NO'
            end if

            if(solide_totali.eq.numero_celle_bloccate)then
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

        if(myid.eq.0)then
            write(*,*)'       LOAD IBM finished'
            write(*,*)'*****************************************'
            write (*,*)' '
        end if

        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !     PREPARE COMUNICATION BETWEEN PROCS FOR IB
        allocate( stencil_left_snd(2*jx*jy,3))
        allocate(stencil_right_snd(2*jx*jy,3))
        allocate( stencil_left_rcv(2*jx*jy,3))
        allocate(stencil_right_rcv(2*jx*jy,3))

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


                if(tipo_spedito(i,j,k)==0)then

                    if(k==kparasta-1)then
                        num_left_rcv = num_left_rcv + 1
                        stencil_left_rcv(num_left_rcv,1) = i
                        stencil_left_rcv(num_left_rcv,2) = j
                        stencil_left_rcv(num_left_rcv,3) = k

                        tipo_spedito(i,j,k) = 1
                    end if


                    if(k==kparaend+1)then
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
        if(myid.ne.0)then
            call MPI_SEND(num_left_rcv,1,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if(myid.ne.nproc-1)then
            call MPI_RECV(num_right_snd,1,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     comunicate right the point right has to send left
        if(myid.ne.nproc-1)then
            call MPI_SEND(num_right_rcv,1,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if(myid.ne.0)then
            call MPI_RECV(num_left_snd,1,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if(kp==0)then
            if(myid.eq.0)then
                call MPI_SEND(num_left_rcv,1,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if(myid.eq.nproc-1)then
                call MPI_RECV(num_right_snd,1,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            if(myid.eq.nproc-1)then
                call MPI_SEND(num_right_rcv,1,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
            end if
            if(myid.eq.0)then
                call MPI_RECV(num_left_snd,1,MPI_INTEGER,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            end if
        else
            if(myid.eq.0)then
                num_left_snd = 0
                num_left_rcv = 0
            end if
            if(myid.eq.nproc-1)then
                num_right_rcv = 0
                num_right_snd = 0
            end if
        end if

        !     I know the counter! now I send the values

        !     send left the indeces left has to send right
        if(myid.ne.0)then
            call MPI_SEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,leftpe,tagls,MPI_COMM_WORLD,ierror)
        end if
        if(myid.ne.nproc-1)then
            call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER,rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
        end if

        !     send right the indeces right has to send left
        if(myid.ne.nproc-1)then
            call MPI_SEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER,rightpe,tagrs,MPI_COMM_WORLD,ierror)
        end if
        if(myid.ne.0)then
            call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER,leftpe,taglr,MPI_COMM_WORLD,status,ierror)
        end if

        !     if periodic
        if(kp==0)then

            !       left to 0 is nproc-1
            if(myid.eq.0)then
                call MPI_SEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            end if
            if(myid.eq.nproc-1)then
                call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER,0,tagrr,MPI_COMM_WORLD,status,ierror)
            end if

            ! right to nproc-1 is 0
            if(myid.eq.nproc-1)then
                call MPI_SEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER,0,tagrs,MPI_COMM_WORLD,ierror)
            end if
            if(myid.eq.0)then
                call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
            end if

            !       now myid=0 and myid=nproc-1 know the index they will recive,
            !       but they need to save
            !       on the border, so I change the k

            if (myid == nproc-1) stencil_right_snd(:,3)=jz
            if (myid == 0) stencil_left_snd(:,3)=1

        end if

        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------

        return

    end subroutine carico_immb

    function global2local(var,l)

        implicit none

        real,intent(in) :: var(3)
        integer,intent(in) :: l
        real,dimension(3) :: global2local

        global2local(1)=var(1)*rot(l,1,1)+var(3)*rot(l,1,2)+var(2)*rot(l,1,3)
        global2local(2)=var(1)*rot(l,2,1)+var(3)*rot(l,2,2)+var(2)*rot(l,2,3)
        global2local(3)=var(1)*rot(l,3,1)+var(3)*rot(l,3,2)+var(2)*rot(l,3,3)

        return

    end function global2local

    function local2global(var,l)

        implicit none

        real,intent(in) :: var(3)
        integer,intent(in) :: l
        real,dimension(3) :: local2global

        local2global(1)=var(1)*rot_inverse(l,1,1)+var(2)*rot_inverse(l,1,2)+var(3)*rot_inverse(l,1,3)
        local2global(2)=var(1)*rot_inverse(l,3,1)+var(2)*rot_inverse(l,3,2)+var(3)*rot_inverse(l,3,3)
        local2global(3)=var(1)*rot_inverse(l,2,1)+var(2)*rot_inverse(l,2,2)+var(3)*rot_inverse(l,2,3)

        return

    end function local2global

    function trilinear_interpolation(var,l)

        use scala3, only: n1,n2
        use mysending, only: kparasta,kparaend,deepl,deepr

        implicit none

        real :: trilinear_interpolation
        real,intent(in) :: var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer,intent(in) :: l
        !real,intent(out) :: trilinear_interpolation

        trilinear_interpolation= &
            tricoef(l,1)*var(trind(l,1,1),trind(l,1,2),trind(l,1,3)) &
            +tricoef(l,2)*var(trind(l,2,1),trind(l,2,2),trind(l,2,3)) &
            +tricoef(l,3)*var(trind(l,3,1),trind(l,3,2),trind(l,3,3)) &
            +tricoef(l,4)*var(trind(l,4,1),trind(l,4,2),trind(l,4,3))

        return

    end function trilinear_interpolation

end module ibm_module

!***********************************************************************





