module ibm_module

    ! for immersed boundary method

    use,intrinsic :: iso_c_binding

    use mysettings, only: particles,bodyforce
    use myarrays_velo3, only: fi,rhov,u,v,w
    use wallmodel_module
    use scala3

    implicit none

    ! --------------------------------------------------

    ! trigger for updating ibm (i.e. ricerca for particles) at every time step
    logical,public :: update_ibm

    ! IBM iterations, from imput file
    integer(kind=c_int),bind(C) :: num_iter
    ! number of IB points (= number of PP, V and IP points), all processors
    integer :: numero_celle_IB
    ! number of IB points (= number of PP, V and IP points), local processor
    integer :: num_ib
    ! number of solid points, all processors
    integer :: numero_celle_bloccate
    ! number of solid points, local processor
    integer :: num_solide

    ! indices of IB points and V points (6,num_ib)
    integer,allocatable :: indici_CELLE_IB(:,:)
    ! indices of solid points (3,num_solide)
    integer,allocatable :: indici_celle_bloccate(:,:)

    ! distance between projection points (PP) and IB points
    real,allocatable :: dist_pp_ib(:)
    ! distance between projection points (IP, on the surface) and IB points (num_ib)
    real,allocatable :: dist_ib_parete(:)
    ! shear velocity at IB points
    real,allocatable :: ustar(:)
    ! used locally wall model for every IB point, determined in correggi
    integer,allocatable :: caso_ib(:)
    ! location in space of IP points (3,num_ib)
    real,allocatable :: proiezioni(:,:)
    ! velocity of the immersed solid wall at IP points (3,num_ib)
    real,allocatable :: surfvel_ib(:,:)
    ! velocity of the immersed body at solid points (3,num_solide)
    real,allocatable :: solidvel_ib(:,:)
      
    ! array for rotation with eulerian angles (direct and inverse)
    ! PROBLEMATIC, left hand convention
    real,allocatable :: rot(:,:,:),rot_inverse(:,:,:)
      
    ! trilinear interpolation: interpolation points
    integer,allocatable :: trind(:,:,:)
    ! trilinear interpolation: interpolation coefficients
    real,allocatable :: tricoef(:,:)
      
    ! buffer
    real,allocatable :: sbuff_ibm(:)
    real,allocatable :: rbuff_ibm(:)

    ! to send the ib stencil
    integer :: num_left_snd,num_right_snd
    integer :: num_left_rcv,num_right_rcv
    integer,allocatable :: stencil_left_snd(:,:)
    integer,allocatable :: stencil_left_rcv(:,:)
    integer,allocatable :: stencil_right_snd(:,:)
    integer,allocatable :: stencil_right_rcv(:,:)
    integer,allocatable :: tipo_spedito(:,:,:)

    ! to send the solid cells
    integer :: numsolid_left_snd,numsolid_right_snd
    integer :: numsolid_left_rcv,numsolid_right_rcv
    integer,allocatable :: solid_left_snd(:,:)
    integer,allocatable :: solid_left_rcv(:,:)
    integer,allocatable :: solid_right_snd(:,:)
    integer,allocatable :: solid_right_rcv(:,:)

    ! ----------------------------------------------------
    ! For particles
    real,allocatable :: pressure_ib(:,:),shear_ib(:,:),momentum_ib(:,:)
    real,allocatable :: position_ib(:,:)
    integer,allocatable :: index_ib(:),indexsize_ib(:)
    ! move this to the subroutine in ricerca (3,num_ib)
    real,allocatable :: surfnormal_ib(:,:)

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
        use period
        !
        use mpi

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

        ! velocity at ib points
        real,dimension(3) :: vel_ib
        real,dimension(3) :: vel_ib_old
        real,dimension(3) :: vel_solido
        real,dimension(3) :: delta_vel
        ! velocity values at the end of ib_procedure
        real,allocatable :: vel_ib_final(:,:)
        ! intermediate velocity values at the end of every iteration
        real,allocatable :: vel_ib_iter(:,:)
        ! velocity values before ib_procedure
        real,allocatable :: vel_ib_oldarray(:,:)

        ! to check conditions
        integer :: caso0_loc,caso1_loc,caso2_loc,notok_loc
        integer :: caso0_tot,caso1_tot,caso2_tot,notok_tot


        !-----------------------------------------------------------------------
        ! check call subroutine

        if (myid==0) then
            write(*,*)'-----------------------------------------------'
            write(*,*)'       IMMERSED BOUNDARY ~ CORREGGI '

            if (ktime==1) then
                write(*,*)'IB proc0',num_ib,'su',numero_celle_IB
                write(*,*)'solide proc0',num_solide,'on',numero_celle_bloccate
            end if
        end if

        allocate(vel_ib_final(3,num_ib))
        allocate(vel_ib_iter(3,num_ib))
        allocate(vel_ib_oldarray(3,num_ib))


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
        allocate(shear_ib(3,num_ib))
        allocate(pressure_ib(3,num_ib))
        allocate(momentum_ib(3,num_ib))
        allocate(caso_ib(num_ib))
        shear_ib(:,:)=0.0
        pressure_ib(:,:)=0.0
        momentum_ib(:,:)=0.0
        caso_ib(:)=-10
        !end if

        !-----------------------------------------------------------------------
        !     solid border
        !        do j=1,jy
        !            do i=1,jx
        !                k=kparasta-1
        !                if (tipo(i,j,k)==0) then
        !                    vel_solido(:)=point_velocity(centroid(:,i,j,k),solid_index(1,i,j,k))
        !                    u(i,j,k)=vel_solido(1)
        !            v(i,j,k)=vel_solido(2)
        !            w(i,j,k)=vel_solido(3)
        !                end if
        !                k=kparaend+1
        !                if (tipo(i,j,k)==0) then
        !                vel_solido(:)=point_velocity(centroid(:,i,j,k),solid_index(1,i,j,k))
        !                    u(i,j,k)=vel_solido(1)
        !            v(i,j,k)=vel_solido(2)
        !            w(i,j,k)=vel_solido(3)
        !                end if
        !            end do
        !        end do


        !-----------------------------------------------------------------------
        ! correction to solid for IB without V

        do l=1,num_ib

            ! index ib
            i0=indici_CELLE_IB(1,l)
            j0=indici_CELLE_IB(2,l)
            k0=indici_CELLE_IB(3,l)

            ! index V
            i=indici_CELLE_IB(4,l)
            j=indici_CELLE_IB(5,l)
            k=indici_CELLE_IB(6,l)

            coincide_check=(i0-i)*(i0-i)+(j0-j)*(j0-j)+(k0-k)*(k0-k)

            if (coincide_check==0) then
                u(i0,j0,k0)=0.0
                v(i0,j0,k0)=0.0
                w(i0,j0,k0)=0.0
                node_is_ok(l)=.false.
                !write(*,*) 'node ',l,' is not ok'
            else
                node_is_ok(l)=.true.
            end if

        end do

        !-----------------------------------------------------------------------
        ! CORRECTION ON IB
        !
        ! initialize
        do l=1,num_ib

            i0=indici_CELLE_IB(1,l)
            j0=indici_CELLE_IB(2,l)
            k0=indici_CELLE_IB(3,l)

            !this is the value before the IBM cycle, used for computing momentum exchange
            vel_ib_oldarray(1,l)=u(i0,j0,k0)
            vel_ib_oldarray(2,l)=v(i0,j0,k0)
            vel_ib_oldarray(3,l)=w(i0,j0,k0)

        end do

        ! first iteration: take old value
        vel_ib_iter(:,:)=vel_ib_oldarray(:,:)

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

            ! exchange planes between procs:
            ! - to compute derivative for Taylor
            ! - to coorect IB values in the iterative procedure
            call passo_ibm()

            do l=1,num_ib

                if (node_is_ok(l)) then !only for IB /= V


                    vel_ib_old(:)=vel_ib_oldarray(:,l)
                    vel_ib(:)=0.0

                    ! compute u_ib and ustar
                    call compute_u_ib(kiter,ktime,l,vel_ib_old,vel_ib,tipo)

                    ! update velocity at IB
                    vel_ib_final(:,l)=vel_ib(:)

                    ! error at ib
                    delta_vel(:)=vel_ib_final(:,l)-vel_ib_iter(:,l)
                    errore=max(abs(delta_vel(1)),abs(delta_vel(2)),abs(delta_vel(3)))

                    ! update error
                    if (errore>=errore_max) then
                        errore_max=errore
                        errore_max_loc=errore
                    end if

                end if

            end do !fine loop su celle ib

            ! reset velocity at IB points
            ! get velocity at IB points from last iteration step
            do l=1,num_ib

                i0=indici_CELLE_IB(1,l)
                j0=indici_CELLE_IB(2,l)
                k0=indici_CELLE_IB(3,l)

                u(i0,j0,k0)=0.0
                v(i0,j0,k0)=0.0
                w(i0,j0,k0)=0.0

            end do

            ! get velocity at IB points from last iteration step
            do l=1,num_ib

                i0=indici_CELLE_IB(1,l)
                j0=indici_CELLE_IB(2,l)
                k0=indici_CELLE_IB(3,l)

                u(i0,j0,k0)=vel_ib_final(1,l)/real(indexsize_ib(l))
                v(i0,j0,k0)=vel_ib_final(2,l)/real(indexsize_ib(l))
                w(i0,j0,k0)=vel_ib_final(3,l)/real(indexsize_ib(l))

            end do

            ! for next time step
            vel_ib_iter(:,:)=vel_ib_final(:,:)

            ! ----------------------------------------------------------------

            call MPI_ALLREDUCE(errore_max_loc,errore_max,1,MPI_REAL_SD,MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (myid==0) then
                write(*,*)'***',kiter,'errore ciclo IBM',errore_max,'***'
            end if

            kiter=kiter+1

        end do ! end loop correzione


        do l=1,num_solide

            i=indici_celle_bloccate(1,l)
            j=indici_celle_bloccate(2,l)
            k=indici_celle_bloccate(3,l)


            vel_solido=solidvel_ib(:,l)

            u(i,j,k)=vel_solido(1)
            v(i,j,k)=vel_solido(2)
            w(i,j,k)=vel_solido(3)

            ! set also pressure to zero
            !fi(i,j,k)=0.0

            ! giulia
            if (i==1) then
                u(i-1,j,k)=vel_solido(1)
                v(i-1,j,k)=vel_solido(2)
                w(i-1,j,k)=vel_solido(3)
            else if (i==n1) then
                u(i+1,j,k)=vel_solido(1)
                v(i+1,j,k)=vel_solido(2)
                w(i+1,j,k)=vel_solido(3)
            end if
            if (j==1) then
                u(i,j-1,k)=vel_solido(1)
                v(i,j-1,k)=vel_solido(2)
                w(i,j-1,k)=vel_solido(3)
            else if (j==n2) then
                u(i,j+1,k)=vel_solido(1)
                v(i,j+1,k)=vel_solido(2)
                w(i,j+1,k)=vel_solido(3)
            end if
            if (myid==0.and.k==1) then
                u(i,j,k-1)=vel_solido(1)
                v(i,j,k-1)=vel_solido(2)
                w(i,j,k-1)=vel_solido(3)
            else if (myid==(nproc-1).and.k==n3) then
                u(i,j,k+1)=vel_solido(1)
                v(i,j,k+1)=vel_solido(2)
                w(i,j,k+1)=vel_solido(3)
            end if

        !        fi(i,j,k)=0.
        end do

        !-----------------------------------------------------------------------
        !     solid border
        !        do j=1,jy
        !            do i=1,jx
        !                k=kparasta-1
        !                if (tipo(i,j,k)==0) then
        !                    vel_solido(:)=point_velocity(centroid(:,i,j,k),solid_index(1,i,j,k))
        !                    u(i,j,k)=vel_solido(1)
        !            v(i,j,k)=vel_solido(2)
        !            w(i,j,k)=vel_solido(3)
        !                end if
        !                k=kparaend+1
        !                if (tipo(i,j,k)==0) then
        !                vel_solido(:)=point_velocity(centroid(:,i,j,k),solid_index(1,i,j,k))
        !                    u(i,j,k)=vel_solido(1)
        !            v(i,j,k)=vel_solido(2)
        !            w(i,j,k)=vel_solido(3)
        !                end if
        !            end do
        !        end do

        !     pass data
        call passo_ibm()

        deallocate(vel_ib_final)
        deallocate(vel_ib_iter)
        deallocate(vel_ib_oldarray)

                ! check conditions
        caso0_loc=0
        caso1_loc=0
        caso2_loc=0
        notok_loc=0
        caso0_tot=0
        caso1_tot=0
        caso2_tot=0
        notok_tot=0
        do l=1,num_ib
            if (caso_ib(l)==0) then
                caso0_loc=caso0_loc+1
            else if (caso_ib(l)==1) then
                caso1_loc=caso1_loc+1
            else if (caso_ib(l)==2) then
                caso2_loc=caso2_loc+1
            end if
            if (.not.node_is_ok(l)) then
                notok_loc=notok_loc+1
            end if
        end do
        call MPI_ALLREDUCE(caso0_loc,caso0_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(caso1_loc,caso1_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(caso2_loc,caso2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(notok_loc,notok_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        if (myid==0) then
            write(*,*) 'statistics'
            write(*,*) '   not ok: ',notok_tot
            write(*,*) '   caso=0: ',caso0_tot
            write(*,*) '   caso=1: ',caso1_tot
            write(*,*) '   caso=2: ',caso2_tot
            write(*,*)'-----------------------------------------------'
        end if

        return

    end subroutine correggi_ib

    subroutine compute_u_ib(kiter,ktime,l,vel_ib_old,vel_ib,tipo)

        use mysettings, only: coef_wall

        integer,intent(in) :: l
        integer,intent(in) :: kiter,ktime
        real,dimension(3),intent(in) :: vel_ib_old
        ! output of the function: new IB point velocity
        real,dimension(3),intent(out) :: vel_ib
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        !-----------------------------------------------------------------------
        integer :: i0,j0,k0
        ! IP_PP distance
        real :: distanza,fattore_distanza
        ! pressure values
        real :: f_pp,f_ib,fi_pro
        integer :: caso
        !real :: pro1,pro2,pro3
        real :: local_shear1,local_shear2
        real :: local_ustar1,local_ustar2
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

        ! index IB
        i0=indici_CELLE_IB(1,l)
        j0=indici_CELLE_IB(2,l)
        k0=indici_CELLE_IB(3,l)

        ! distance IP-PP
        distanza=dist_ib_parete(l)+dist_pp_ib(l)
        fattore_distanza=dist_ib_parete(l)/distanza

        ! velocity at interpolation points (PP, V?)
        vel_pp(1)=trilinear_interpolation(u,l,tipo)
        vel_pp(2)=trilinear_interpolation(v,l,tipo)
        vel_pp(3)=trilinear_interpolation(w,l,tipo)

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
            vel_ip(:)=surfvel_ib(:,l)
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
        else if ((dist_ib_parete(l)<1.0e-9).or.(vtan_pp_rel<1.0e-9)) then
            ! values too small: zero velocity
            caso=0
            !write(*,*) ' vel_ib=(',vel_ib(1),' ',vel_ib(2),' ',vel_ib(3),')'
            !write(*,*) 'caso=0, l=',l
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

        if (caso==2) then

            ! find the two components of ustar (ustar is aligned with vtan_rel)
            local_ustar1=vloc_pp_rel(1)*ustarHere/vtan_pp_rel
            local_ustar2=vloc_pp_rel(2)*ustarHere/vtan_pp_rel

            ! compute components of shear stress (conserve sign) tau=rho*ustar^2
            local_shear1=sign((local_ustar1)**2,local_ustar1)
            local_shear2=sign((local_ustar2)**2,local_ustar2)

            ! force on solid objects are made up of two components: pressure and shear (2,3 INVERTED)
            shear_ib(1,l)=local_shear1*rot_inverse(1,1,l)+local_shear2*rot_inverse(1,2,l)
            shear_ib(3,l)=local_shear1*rot_inverse(2,1,l)+local_shear2*rot_inverse(2,2,l)
            shear_ib(2,l)=local_shear1*rot_inverse(3,1,l)+local_shear2*rot_inverse(3,2,l)

        end if



        ! pressure is independent of caso
        f_ib=fi(i0,j0,k0)

        ! OLD: this assumes a linear behavior of pressure along the normal line
        !f_pp=trilinear_interpolation(fi,l)
        !fi_pro=((dist_ib_parete(l)+dist_pp_ib(l))*f_ib-f_pp*dist_ib_parete(l))/dist_pp_ib(l)

        ! NEW: we take the pressure value at the IB node
        fi_pro=f_ib

        ! pressure has a minus sign because the normal is directed outwards (2,3 INVERTED)
        pressure_ib(1,l)=-1.0*fi_pro*rot_inverse(1,3,l)
        pressure_ib(3,l)=-1.0*fi_pro*rot_inverse(2,3,l)
        pressure_ib(2,l)=-1.0*fi_pro*rot_inverse(3,3,l)

        momentum_ib(:,l)=(vel_ib_old(:)-vel_ib(:))/dt

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
        use period
        !
        use mpi

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
            end if
        end if

        if(num_right_rcv.ne.0)then
            if(myid.ne.nproc-1)then
                call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_right_rcv,MPI_REAL_SD, &
                    rightpe,tagrr,MPI_COMM_WORLD,status,ierror)

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
            end if
        end if

        if(num_left_rcv.ne.0)then
            if(myid.ne.0)then
                call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_left_rcv,MPI_REAL_SD, &
                    leftpe,taglr,MPI_COMM_WORLD,status,ierror)

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
                end if
            end if

            if(num_left_rcv.ne.0)then
                if(myid.eq.0)then
                    call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_left_rcv, &
                        MPI_REAL_SD,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)

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
                do j=1,n2

                    u(0   ,j,k)=u(n1,j,k)
                    v(0   ,j,k)=v(n1,j,k)
                    w(0   ,j,k)=w(n1,j,k)
                    u(n1+1,j,k)=u(1 ,j,k)
                    v(n1+1,j,k)=v(1 ,j,k)
                    w(n1+1,j,k)=w(1 ,j,k)
                !
                end do
            end do

        end if
        !
        !-----------------------------------------------------------------------------
        !     j direction, impose periodicity or wall
        !
        !     parabolic extrapolation y=ax2+bx+c
        ! direction 2 is always not periodic
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
            do j=1,n2
                do i=1,n1
                    write(4500+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do

            k=kparasta-1
            do j=1,n2
                do i=1,n1
                    write(4600+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do

            k=kparaend
            do j=1,n2
                do i=1,n1
                    write(4700+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do

            k=kparaend+1
            do j=1,n2
                do i=1,n1
                    write(4800+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do


        end if
        !-----------------------------------------------------------------------------
        return
    end

    function global2local(var,l)

        real,intent(in) :: var(3)
        integer,intent(in) :: l
        real,dimension(3) :: global2local

        global2local(1)=var(1)*rot(1,1,l)+var(3)*rot(1,2,l)+var(2)*rot(1,3,l)
        global2local(2)=var(1)*rot(2,1,l)+var(3)*rot(2,2,l)+var(2)*rot(2,3,l)
        global2local(3)=var(1)*rot(3,1,l)+var(3)*rot(3,2,l)+var(2)*rot(3,3,l)

        return

    end function global2local

    function local2global(var,l)

        real,intent(in) :: var(3)
        integer,intent(in) :: l
        real,dimension(3) :: local2global

        local2global(1)=var(1)*rot_inverse(1,1,l)+var(2)*rot_inverse(1,2,l)+var(3)*rot_inverse(1,3,l)
        local2global(2)=var(1)*rot_inverse(3,1,l)+var(2)*rot_inverse(3,2,l)+var(3)*rot_inverse(3,3,l)
        local2global(3)=var(1)*rot_inverse(2,1,l)+var(2)*rot_inverse(2,2,l)+var(3)*rot_inverse(2,3,l)

        return

    end function local2global

    function trilinear_interpolation(var,l,tipo)

        real :: trilinear_interpolation
        real,intent(in) :: var(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer,intent(in) :: l
        !real,intent(out) :: trilinear_interpolation
        integer :: a

        trilinear_interpolation= &
            tricoef(1,l)*var(trind(1,1,l),trind(1,2,l),trind(1,3,l)) &
            +tricoef(2,l)*var(trind(2,1,l),trind(2,2,l),trind(2,3,l)) &
            +tricoef(3,l)*var(trind(3,1,l),trind(3,2,l),trind(3,3,l)) &
            +tricoef(4,l)*var(trind(4,1,l),trind(4,2,l),trind(4,3,l))

!            if (tipo(trind(1,1,l),trind(1,2,l),trind(1,3,l))==0 .and. &
!                tipo(trind(2,1,l),trind(2,2,l),trind(2,3,l))==0 .and. &
!                tipo(trind(3,1,l),trind(3,2,l),trind(3,3,l))==0 .and. &
!                tipo(trind(4,1,l),trind(4,2,l),trind(4,3,l))==0 .and. &
!                var(trind(1,1,l),trind(1,2,l),trind(1,3,l))==0.0 .and. &
!                var(trind(2,1,l),trind(2,2,l),trind(2,3,l))==0.0 .and. &
!                var(trind(3,1,l),trind(3,2,l),trind(3,3,l))==0.0 .and. &
!                var(trind(4,1,l),trind(4,2,l),trind(4,3,l))==0.0) then
!
!                write(*,*) 'both problems(0) ',l
!
!            else if (tipo(trind(1,1,l),trind(1,2,l),trind(1,3,l))==1 .and. &
!                tipo(trind(2,1,l),trind(2,2,l),trind(2,3,l))==1 .and. &
!                tipo(trind(3,1,l),trind(3,2,l),trind(3,3,l))==1 .and. &
!                tipo(trind(4,1,l),trind(4,2,l),trind(4,3,l))==1 .and. &
!                var(trind(1,1,l),trind(1,2,l),trind(1,3,l))==0.0 .and. &
!                var(trind(2,1,l),trind(2,2,l),trind(2,3,l))==0.0 .and. &
!                var(trind(3,1,l),trind(3,2,l),trind(3,3,l))==0.0 .and. &
!                var(trind(4,1,l),trind(4,2,l),trind(4,3,l))==0.0) then
!
!            write(*,*) 'both problems(1) ',l
!
!                else
!
!            if (tipo(trind(1,1,l),trind(1,2,l),trind(1,3,l))==0 .and. &
!                tipo(trind(2,1,l),trind(2,2,l),trind(2,3,l))==0 .and. &
!                tipo(trind(3,1,l),trind(3,2,l),trind(3,3,l))==0 .and. &
!                tipo(trind(4,1,l),trind(4,2,l),trind(4,3,l))==0) then
!                write(*,*) 'tipo problem(0) ',l
!            end if
!
!            if (tipo(trind(1,1,l),trind(1,2,l),trind(1,3,l))==1 .and. &
!                tipo(trind(2,1,l),trind(2,2,l),trind(2,3,l))==1 .and. &
!                tipo(trind(3,1,l),trind(3,2,l),trind(3,3,l))==1 .and. &
!                tipo(trind(4,1,l),trind(4,2,l),trind(4,3,l))==1) then
!                write(*,*) 'tipo problem(0) ',l
!            end if
!
!            if (var(trind(1,1,l),trind(1,2,l),trind(1,3,l))==0.0 .and. &
!                var(trind(2,1,l),trind(2,2,l),trind(2,3,l))==0.0 .and. &
!                var(trind(3,1,l),trind(3,2,l),trind(3,3,l))==0.0 .and. &
!                var(trind(4,1,l),trind(4,2,l),trind(4,3,l))==0.0) then
!                write(*,*) 'var problem ',l
!            end if
!
!            end if




        return

    end function trilinear_interpolation

end module ibm_module

!***********************************************************************





