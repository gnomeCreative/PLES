module particle_module

    use,intrinsic :: iso_c_binding

    use ibm_module
    use mysending

    use mpi

    implicit none

    integer :: num_part_loc,num_part_tot

    real,allocatable :: sphereRadius(:),sphereRadius2(:),sphereSurface(:)
    real,allocatable :: spherePosition(:,:)
    real,allocatable :: sphereVelocity(:,:)
    real,allocatable :: sphereSpin(:,:)
    integer,allocatable :: sphereIndex(:)
    logical,allocatable :: sphereMoves(:)
    real,allocatable :: sphereShearForce(:,:),spherePressureForce(:,:),sphereMomentumForce(:,:)
    ! field vector for visualization
    real,allocatable :: fluidShearForce(:,:,:,:),fluidPressureForce(:,:,:,:),fluidMomentumForce(:,:,:,:)
    real,allocatable :: fluidParticleVel(:,:,:,:)
    integer,allocatable :: forceCaso(:,:,:)

    real :: border_left,border_right
    real :: deep_border_left,deep_border_right

    logical,parameter :: print_fields=.false.

contains

    subroutine compute_sphere_forces()

        use myarrays_metri3, only: giac,ref_length,ref_area
        use scala3, only: dt

        ! for force geometric computation
        ! IP_IB distance
        real :: refDistance
        ! characteristics of fluid element
        real :: refLengthHere,refAreaHere,volumeHere,refVolume
        real,parameter :: sqrt2=sqrt(2.0)
        real,parameter :: sqrt2d2=sqrt(2.0)/2.0

        ! for output
        integer :: l,m
        integer :: ierr
        integer :: solidIndexHere,p,indexSizeHere,particleIndex
        integer :: ibCounter
        real,dimension(3) :: forceShearHere,forceShearTot
        real,dimension(3) :: forcePressureHere,forcePressureTot
        real,dimension(3) :: forceMomentumHere,forceMomentumTot
        !
        integer :: i0,j0,k0

        ! remove this
        if (print_fields) then
            call compute_sphere_force_fields()
        end if

        ! deallocate previous arrays if necessary
        if (allocated(sphereShearForce)) then
            deallocate(sphereShearForce,spherePressureForce,sphereMomentumForce)
        end if

        allocate(sphereShearForce(3,num_part_tot))
        allocate(spherePressureForce(3,num_part_tot))
        allocate(sphereMomentumForce(3,num_part_tot))


        ! Force on particles
        sphereShearForce(:,:)=0.0
        spherePressureForce(:,:)=0.0
        sphereMomentumForce(:,:)=0.0

        do p=1,num_part_tot

            particleIndex=sphereIndex(p)


            forceShearHere(:)=0.0
            forcePressureHere(:)=0.0
            forceMomentumHere(:)=0.0
            forceShearTot(:)=0.0
            forcePressureTot(:)=0.0
            forceMomentumTot(:)=0.0

            ibCounter=0

            do l=1,num_ib

                ! index ib
                i0=indici_CELLE_IB(1,l)
                j0=indici_CELLE_IB(2,l)
                k0=indici_CELLE_IB(3,l)

                solidIndexHere=index_ib(l)

                if (solidIndexHere==particleIndex) then

                    ! momentum needs an estimation of the reference volume of the IB node
                    ! for ref_distance->0, ref_volume=1/2*cell_volume

                    ! distance IP-IB
                    refDistance=dist_ib_parete(l)

                    ! refence length of fluid element
                    refLengthHere=ref_length(i0,j0,k0)
                    refAreaHere=ref_area(i0,j0,k0)
                    volumeHere=giac(i0,j0,k0)

                    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    ! CHECK THIS !!!!!!!!!!!!!!!!!!!!!!
                    if (refDistance<refLengthHere) then
                        refVolume=0.5*volumeHere*(1.0+refDistance/refLengthHere)
                    else
                        refVolume=volumeHere*(1.0+refDistance/refLengthHere)
                    end if

                    !write(*,*) 'VOL=',refVolume

                    forceShearHere(:)=forceShearHere(:)+shear_ib(:,l)*refAreaHere
                    forcePressureHere(:)=forcePressureHere(:)+pressure_ib(:,l)*refAreaHere
                    forceMomentumHere(:)=forceMomentumHere(:)+momentum_ib(:,l)*refVolume

                    ibCounter=ibCounter+1

                end if

            end do

            call MPI_ALLREDUCE(forceShearHere(1),forceShearTot(1),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(forceShearHere(2),forceShearTot(2),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(forceShearHere(3),forceShearTot(3),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

            call MPI_ALLREDUCE(forcePressureHere(1),forcePressureTot(1),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(forcePressureHere(2),forcePressureTot(2),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(forcePressureHere(3),forcePressureTot(3),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

            call MPI_ALLREDUCE(forceMomentumHere(1),forceMomentumTot(1),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(forceMomentumHere(2),forceMomentumTot(2),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(forceMomentumHere(3),forceMomentumTot(3),1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            sphereShearForce(:,particleIndex)=forceShearTot(:)
            spherePressureForce(:,particleIndex)=forcePressureTot(:)
            sphereMomentumForce(:,particleIndex)=forceMomentumTot(:)


        end do


    end subroutine compute_sphere_forces

    subroutine compute_sphere_force_fields()

        use scala3

        integer :: i,j,k,l

        if (.not.allocated(fluidShearForce)) then
            allocate(fluidPressureForce(3,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            allocate(fluidShearForce(3,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            allocate(fluidMomentumForce(3,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            allocate(fluidParticleVel(3,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            allocate(forceCaso(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        end if

            fluidShearForce(:,:,:,:)=0.0
            fluidPressureForce(:,:,:,:)=0.0
            fluidMomentumForce(:,:,:,:)=0.0
            fluidParticleVel(:,:,:,:)=0.0
            forceCaso(:,:,:)=-10

        do l=1,num_ib

            i=indici_CELLE_IB(1,l)
            j=indici_CELLE_IB(2,l)
            k=indici_CELLE_IB(3,l)

            fluidShearForce(:,i,j,k)=shear_ib(:,l)

            fluidPressureForce(:,i,j,k)=pressure_ib(:,l)

            fluidMomentumForce(:,i,j,k)=momentum_ib(:,l)

            fluidParticleVel(:,i,j,k)=surfvel_ib(:,l)

            forceCaso(i,j,k)=caso_ib(l)

        end do


    end subroutine compute_sphere_force_fields

    subroutine pass_geometry(totpart,pposx,pposy,pposz,pvelx,pvely,pvelz,pspinx,pspiny,pspinz,prad,pmoves) &
        bind (C, name="pass_geometry")

        use myarrays_metri3, only: z
        use scala3

        integer(c_int),intent(in) :: totpart
        real(c_double),intent(in) :: pposx(totpart),pposy(totpart),pposz(totpart)
        real(c_double),intent(in) :: pvelx(totpart),pvely(totpart),pvelz(totpart)
        real(c_double),intent(in) :: pspinx(totpart),pspiny(totpart),pspinz(totpart)
        real(c_double),intent(in) :: prad(totpart)
        logical(c_bool),intent(in) :: pmoves(totpart)
        integer :: p
        integer :: moving_counter,ierr
        real :: max_radius
        !-------------------------------------------------


        ! compute max radius for determination of borders between processors
        max_radius=maxval(prad)
        write(*,*) 'Max radius ',max_radius
        border_left=z(1,1,kparasta-1)
        border_right=z(1,1,kparaend+1)
        deep_border_left=border_left-1.5*max_radius
        deep_border_right=border_right+1.5*max_radius
        write(*,*) 'Proc ',myid,' left=',deep_border_left,' right=',deep_border_right

        ! total number of particles
        num_part_tot=totpart

        ! deallocate previous arrays if necessary
        if (allocated(sphereIndex)) then
            deallocate(sphereIndex,sphereMoves)
            deallocate(sphereRadius,sphereRadius2)
            deallocate(spherePosition)
            deallocate(sphereVelocity)
            deallocate(sphereSpin)
            deallocate(sphereSurface)
        end if

        allocate(sphereIndex(num_part_tot))
        allocate(sphereMoves(num_part_tot))
        allocate(sphereRadius(num_part_tot),sphereRadius2(num_part_tot))
        allocate(spherePosition(3,num_part_tot))
        allocate(sphereVelocity(3,num_part_tot))
        allocate(sphereSpin(3,num_part_tot))
        allocate(sphereSurface(num_part_tot))

        moving_counter=0

        do p=1,num_part_tot
            !write(*,*) 'Read',n
            sphereIndex(p)=p
            sphereMoves(p)=pmoves(p)
            if (sphereMoves(p)) then
                moving_counter=moving_counter+1
            end if
            ! directly taken variables
            spherePosition(:,p)=(/pposx(p),pposy(p),pposz(p)/)
            sphereVelocity(:,p)=(/pvelx(p),pvely(p),pvelz(p)/)
            sphereSpin(:,p)=(/pspinx(p),pspiny(p),pspinz(p)/)
            sphereRadius(p)=prad(p)
            ! computed variables
            sphereRadius2(p)=sphereRadius(p)**2.0
            sphereSurface(p)=4.0*3.1428*sphereRadius(p)*sphereRadius(p)
        end do

        if (moving_counter>0) then
            update_ibm=.true.
        end if

!        if (myid==0) then
!            write(*,'(A,I0.1,A)') ' LES: Got geometry data of ',totParticles,' particle(s) from DEM:'
!            do p=1,totParticles
!                write(*,'(I0.1,ES11.2E3,ES11.2E3,ES11.2E3,ES11.2E3,A,L1)') &
!                    sphereIndex(p),sphereRadius(p),spherePosition(1,p),spherePosition(2,p),spherePosition(3,p),' ',sphereMoves(p)
!            end do
!        end if

        return

    end subroutine pass_geometry

    subroutine pass_forces(sforx,sfory,sforz,pforx,pfory,pforz,mforx,mfory,mforz) bind (C, name="pass_forces")

        real(c_double),intent(out) :: sforx(num_part_tot),sfory(num_part_tot),sforz(num_part_tot)
        real(c_double),intent(out) :: pforx(num_part_tot),pfory(num_part_tot),pforz(num_part_tot)
        real(c_double),intent(out) :: mforx(num_part_tot),mfory(num_part_tot),mforz(num_part_tot)

        integer :: p

        if (allocated(sphereShearForce)) then
            !-------------------------------------------------
            if (myid==1) then
                write(*,'(A,I0.1,A)') 'LES: Send force data of ',num_part_tot,' particle(s) to DEM '
            end if

            do p=1,num_part_tot
                sforx(p)=sphereShearForce(1,p)
                sfory(p)=sphereShearForce(2,p)
                sforz(p)=sphereShearForce(3,p)
                pforx(p)=spherePressureForce(1,p)
                pfory(p)=spherePressureForce(2,p)
                pforz(p)=spherePressureForce(3,p)
                mforx(p)=sphereMomentumForce(1,p)
                mfory(p)=sphereMomentumForce(2,p)
                mforz(p)=sphereMomentumForce(3,p)
            end do

!            if (myid==1) then
!                write(*,*) 'Forces:'
!                do p=1,totParticles
!                    write(*,*) sforx(p),sfory(p),sforz(p),&
!                        pforx(p),pfory(p),pforz(p),&
!                        mforx(p),mfory(p),mforz(p)
!                end do
!            end if

        else

            do p=1,num_part_tot
                sforx(p)=0.0
                sfory(p)=0.0
                sforz(p)=0.0
                pforx(p)=0.0
                pfory(p)=0.0
                pforz(p)=0.0
                mforx(p)=0.0
                mfory(p)=0.0
                mforz(p)=0.0
            end do

            if (myid==0) then
                write(*,*) 'No forces from LES, values resetted'
            end if

        end if

        return

    end subroutine pass_forces

    function is_inside(point,p_index) result (check)

        use scala3, only: dot

        ! ----------------------------------------------------
        logical :: check
        real,dimension(3),intent(in) :: point
        integer,intent(in) :: p_index
        ! ----------------------------------------------------
        real :: iso
        ! ----------------------------------------------------

        iso=dot(point(:)-spherePosition(:,p_index),point(:)-spherePosition(:,p_index))

        if (iso>sphereRadius2(p_index)) then
            check=.false.
            return
        else
            check=.true.
            return
        end if

    end function is_inside


    function is_inprocessor(p_index) result (check)

        ! ----------------------------------------------------
        logical :: check
        integer,intent(in) :: p_index
        ! ----------------------------------------------------


        if (spherePosition(3,p_index)>deep_border_left .and. spherePosition(3,p_index)<deep_border_right) then
            check=.true.
            return
        else
            check=.false.
            return
        end if

    end function is_inprocessor


    function point_velocity(point,p_index) result (velocity)

        use scala3, only: cross

        ! ----------------------------------------------------
        real,dimension(3) :: velocity
        real,dimension(3),intent(in) :: point
        integer,intent(in) :: p_index
        ! ----------------------------------------------------
        real,dimension(3) :: lever,center_distance,normal
        real :: norm_center_distance
        ! ----------------------------------------------------

        lever(:)=point(:)-spherePosition(:,p_index)

        velocity(:)=sphereVelocity(:,p_index)+cross(lever,sphereSpin(:,p_index))

    end function point_velocity

    function surface_point(normal,p_index) result (point)

        ! find direction vector passing by the point and pointing outward from
        ! particle center

        ! ----------------------------------------------------
        real,dimension(3) :: point
        real,dimension(3),intent(in) :: normal
        integer,intent(in) :: p_index
        ! ----------------------------------------------------

        point(:)=spherePosition(:,p_index)+normal(:)*sphereRadius(p_index)

    end function surface_point

    function surface_normal(point,p_index) result (normal)

        ! find direction vector passing by the point and pointing outward from
        ! particle center

        ! ----------------------------------------------------
        real,dimension(3) :: normal
        real,dimension(3),intent(in) :: point
        integer,intent(in) :: p_index
        ! ----------------------------------------------------

        normal(:)=point(:)-spherePosition(:,p_index)

        normal(:)=normal(:)/norm2(normal)

    end function surface_normal

    function surface_distance(point,p_index) result (distance)

        ! ----------------------------------------------------
        real,dimension(3) :: distance
        real,dimension(3),intent(in) :: point
        integer,intent(in) :: p_index
        ! ----------------------------------------------------

        ! ----------------------------------------------------

        distance(:)=point(:)-spherePosition(:,p_index)-surface_normal(point,p_index)*sphereRadius(p_index)

    end function surface_distance

end module particle_module
