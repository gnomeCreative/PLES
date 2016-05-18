module particle_module

    use,intrinsic :: iso_c_binding
    use ibm_module
    use mysending
    !
    use mpi

    implicit none

    private

    integer,public :: totParticles

    real,allocatable,public :: sphereRadius(:),sphereRadius2(:),sphereSurface(:)
    real,allocatable,public :: spherePosition(:,:)
    real,allocatable,public :: sphereVelocity(:,:)
    real,allocatable,public :: sphereSpin(:,:)
    integer,allocatable :: sphereIndex(:)
    logical,allocatable :: sphereMoves(:)
    real,allocatable,public :: sphereShearForce(:,:),spherePressureForce(:,:),sphereMomentumForce(:,:)
    ! field vector for visualization
    real,allocatable,public :: fluidShearForce(:,:,:,:),fluidPressureForce(:,:,:,:),fluidMomentumForce(:,:,:,:)
    real,allocatable,public :: fluidParticleVel(:,:,:,:)
    integer,allocatable,public :: forceCaso(:,:,:)


    public :: pass_geometry,pass_forces,compute_sphere_forces,compute_sphere_force_fields

contains


    subroutine compute_sphere_forces()

        use myarrays_metri3, only: giac,ref_length
        use scala3, only: dt

        implicit none

        ! for force geometric computation
        ! IP_IB distance
        real :: refDistance
        ! characteristics of fluid element
        real :: refLengthHere,volumeHere,refVolume
        real,parameter :: sqrt2=sqrt(2.0)
        real,parameter :: sqrt2d2=sqrt(2.0)/2.0

        ! for output
        integer :: l,m
        integer :: ierr
        integer :: solidIndexHere,p,indexHere,indexSizeHere
        integer :: ibCounter
        integer :: MPI_3_REAL
        real,dimension(3) :: forceShearHere,forceShearTot
        real,dimension(3) :: forcePressureHere,forcePressureTot
        real,dimension(3) :: forceMomentumHere,forceMomentumTot

        !
        integer :: i0,j0,k0

        ! remove this
        call compute_sphere_force_fields()

        ! deallocate previous arrays if necessary
        if (allocated(sphereShearForce)) then
            deallocate(sphereShearForce,spherePressureForce,sphereMomentumForce)
        end if

        allocate(sphereShearForce(3,totParticles))
        allocate(spherePressureForce(3,totParticles))
        allocate(sphereMomentumForce(3,totParticles))


        ! Force on particles
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !call MPI_TYPE_CONTIGUOUS(3,MPI_REAL_SD,MPI_3_REAL,ierr)

        sphereShearForce(:,:)=0.0
        spherePressureForce(:,:)=0.0
        sphereMomentumForce(:,:)=0.0

        do p=1,totParticles


            forceShearHere(:)=0.0
            forcePressureHere(:)=0.0
            forceMomentumHere(:)=0.0
            forceShearTot(:)=0.0
            forcePressureTot(:)=0.0
            forceMomentumTot(:)=0.0

            indexHere=sphereIndex(p)

            ibCounter=0

            do l=1,num_ib

                ! index ib
                i0=indici_CELLE_IB(l,1)
                j0=indici_CELLE_IB(l,2)
                k0=indici_CELLE_IB(l,3)

                do m=1,indexsize_ib(l)

                    solidIndexHere=index_ib(l,m)

                    if (solidIndexHere==indexHere) then

                        ! momentum needs an estimation of the reference volume of the IB node
                        ! for ref_distance->0, ref_volume=1/2*cell_volume

                        ! distance IP-IB
                        refDistance=dist_ib_parete(l)

                        ! refence length of fluid element
                        refLengthHere=ref_length(i0,j0,k0)
                        volumeHere=giac(i0,j0,k0)

                        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        ! CHECK THIS !!!!!!!!!!!!!!!!!!!!!!
                        if (refDistance<refLengthHere) then
                            refVolume=0.5*volumeHere*(1.0+refDistance/refLengthHere)
                        else
                            refVolume=volumeHere*(1.0+refDistance/refLengthHere)
                        end if

                        !write(*,*) 'VOL=',refVolume

                        forceShearHere(:)=forceShearHere(:)+shear_ib(l,:)*refLengthHere**2.0
                        forcePressureHere(:)=forcePressureHere(:)+pressure_ib(l,:)*refLengthHere**2.0
                        forceMomentumHere(:)=forceMomentumHere(:)+momentum_ib(l,:)*refVolume/dt

                        ibCounter=ibCounter+1

                    end if

                end do

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

            sphereShearForce(:,p)=forceShearTot(:)
            spherePressureForce(:,p)=forcePressureTot(:)
            sphereMomentumForce(:,p)=forceMomentumTot(:)


        end do


    end subroutine compute_sphere_forces

    subroutine compute_sphere_force_fields()

        use scala3

        implicit none

        integer :: i0,j0,k0,l

        if (allocated(fluidShearForce)) then
                    ! Alessandro: interaction force --------------------------------
            fluidShearForce(:,:,:,:)=0.0
            fluidPressureForce(:,:,:,:)=0.0
            fluidMomentumForce(:,:,:,:)=0.0
            fluidParticleVel(:,:,:,:)=0.0
            forceCaso(:,:,:)=5
        else
            allocate(fluidPressureForce(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,3))
            allocate(fluidShearForce(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,3))
            allocate(fluidMomentumForce(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,3))
            allocate(fluidParticleVel(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,3))
            allocate(forceCaso(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            fluidShearForce(:,:,:,:)=0.0
            fluidPressureForce(:,:,:,:)=0.0
            fluidMomentumForce(:,:,:,:)=0.0
            fluidParticleVel(:,:,:,:)=0.0
            forceCaso(:,:,:)=5
        end if

        do l=1,num_ib

            i0=indici_CELLE_IB(l,1)
            j0=indici_CELLE_IB(l,2)
            k0=indici_CELLE_IB(l,3)

            fluidShearForce(i0,j0,k0,:)=shear_ib(l,:)

            fluidPressureForce(i0,j0,k0,:)=pressure_ib(l,:)

            fluidMomentumForce(i0,j0,k0,:)=momentum_ib(l,:)

            fluidParticleVel(i0,j0,k0,:)=surfVel(l,:)

            forceCaso(i0,j0,k0)=caso_ib(l)

        end do

    end subroutine compute_sphere_force_fields

    subroutine pass_geometry(totpart,pposx,pposy,pposz,pvelx,pvely,pvelz,pspinx,pspiny,pspinz,prad,pmoves) &
        bind (C, name="pass_geometry")

        implicit none

        integer(c_int),intent(in) :: totpart
        real(c_double),intent(in) :: pposx(totpart),pposy(totpart),pposz(totpart)
        real(c_double),intent(in) :: pvelx(totpart),pvely(totpart),pvelz(totpart)
        real(c_double),intent(in) :: pspinx(totpart),pspiny(totpart),pspinz(totpart)
        real(c_double),intent(in) :: prad(totpart)
        logical(c_bool),intent(in) :: pmoves(totpart)
        integer :: p

        !-------------------------------------------------
        totParticles=totpart

        ! deallocate previous arrays if necessary
        if (allocated(sphereIndex)) then
            deallocate(sphereIndex,sphereMoves)
            deallocate(sphereRadius,sphereRadius2)
            deallocate(spherePosition)
            deallocate(sphereVelocity)
            deallocate(sphereSpin)
            deallocate(sphereSurface)
        end if

        allocate(sphereIndex(totParticles),sphereMoves(totParticles))
        allocate(sphereRadius(totParticles),sphereRadius2(totParticles))
        allocate(spherePosition(3,totParticles))
        allocate(sphereVelocity(3,totParticles))
        allocate(sphereSpin(3,totParticles))
        allocate(sphereSurface(totParticles))

        do p=1,totParticles
            !write(*,*) 'Read',n
            sphereIndex(p)=p
            sphereMoves(p)=pmoves(p)
            spherePosition(1,p)=pposx(p)
            spherePosition(2,p)=pposy(p)
            spherePosition(3,p)=pposz(p)
            sphereVelocity(1,p)=pvelx(p)
            sphereVelocity(2,p)=pvely(p)
            sphereVelocity(3,p)=pvelz(p)
            sphereSpin(1,p)=pspinx(p)
            sphereSpin(2,p)=pspiny(p)
            sphereSpin(3,p)=pspinz(p)
            sphereRadius(p)=prad(p)
            sphereRadius2(p)=sphereRadius(p)**2.0
            sphereSurface(p)=4.0*3.1428*sphereRadius(p)*sphereRadius(p)
        end do

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

        implicit none

        real(c_double),intent(out) :: sforx(totParticles),sfory(totParticles),sforz(totParticles)
        real(c_double),intent(out) :: pforx(totParticles),pfory(totParticles),pforz(totParticles)
        real(c_double),intent(out) :: mforx(totParticles),mfory(totParticles),mforz(totParticles)

        integer :: p

        if (allocated(sphereShearForce)) then
            !-------------------------------------------------
            if (myid==1) then
                write(*,'(A,I0.1,A)') 'LES: Send force data of ',totParticles,' particle(s) to DEM '
            end if

            do p=1,totParticles
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

            do p=1,totParticles
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


end module particle_module
