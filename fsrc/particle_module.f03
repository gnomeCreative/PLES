module particle_module

    use,intrinsic :: iso_c_binding

    implicit none

    private

    integer,public :: totParticles
    integer :: totMoving,totFixed

    real,allocatable,public :: sphereRadius(:),sphereRadius2(:),sphereSurface(:)
    real,allocatable,public :: spherePosition(:,:)
    real,allocatable,public :: sphereShearForce(:,:)
    real,allocatable,public :: spherePressureForce(:,:)
    real,allocatable,public :: fluidShearForce(:,:,:,:),fluidPressureForce(:,:,:,:)
    integer,allocatable,public :: forceCaso(:,:,:)
    integer,allocatable :: sphereIndex(:)
    logical,allocatable :: sphereMoves(:)

    public :: pass_geometry,pass_forces,compute_sphere_forces,compute_sphere_force_fields

contains


    subroutine compute_sphere_forces()

        use myarrays_ibm
        use mysending

        use tipologia
        use mpi

        implicit none

        ! for output
        integer :: i0,j0,k0,l
        integer :: ierr
        integer :: solidIndexHere,p,indexHere
        integer :: ibCounter
        integer :: MPI_3_REAL
        real,dimension(3) :: forceShearHere,forcePressureHere,forceShearTot,forcePressureTot

        ! remove this
        call compute_sphere_force_fields()

        ! Force on particles
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        call MPI_TYPE_CONTIGUOUS(3,MPI_REAL_SD,MPI_3_REAL,ierr)

        sphereShearForce(:,:)=0.0
        spherePressureForce(:,:)=0.0

        do p=1,totParticles

            forceShearHere(:)=0.0
            forcePressureHere(:)=0.0
            forceShearTot(:)=0.0
            forcePressureTot(:)=0.0

            indexHere=sphereIndex(p)

            ibCounter=0

            do l=1,num_ib

                i0=indici_CELLE_IB(l,1)
                j0=indici_CELLE_IB(l,2)
                k0=indici_CELLE_IB(l,3)

                solidIndexHere=solidIndex(i0,j0,k0,1)

                if (solidIndexHere==indexHere) then

                    forceShearHere(:)=forceShearHere(:)+shear_ib(l,:)
                    forcePressureHere(:)=forcePressureHere(:)+pressure_ib(l,:)

                    ibCounter=ibCounter+1

                end if

            end do

            if (ibCounter>0) then

                forceShearHere(:)=forceShearHere(:)/real(ibCounter)
                forcePressureHere(:)=forcePressureHere(:)/real(ibCounter)

            end if

            call MPI_REDUCE(forceShearHere(1),forceShearTot(1),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(forceShearHere(2),forceShearTot(2),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(forceShearHere(3),forceShearTot(3),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)

            call MPI_REDUCE(forcePressureHere(1),forcePressureTot(1),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(forcePressureHere(2),forcePressureTot(2),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(forcePressureHere(3),forcePressureTot(3),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)

            sphereShearForce(:,p)=forceShearTot(:)
            spherePressureForce(:,p)=forcePressureTot(:)


        end do

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        do p=1,totParticles

            sphereShearForce(:,p)=sphereShearForce(:,p)*sphereSurface(p)

            spherePressureForce(:,p)=spherePressureForce(:,p)*sphereSurface(p)

            if (myid==0) then
                !write(*,*) 'Tot particles: ',totParticles
                write(*,*) '--> Shear force = (' &
                    ,sphereShearForce(1,p),',',sphereShearForce(2,p),',',sphereShearForce(3,p),')'
                write(*,*) '--> Pressure force = (' &
                    ,spherePressureForce(1,p),',',spherePressureForce(2,p),',',spherePressureForce(3,p),')'
            end if
        end do


    end subroutine compute_sphere_forces

    subroutine compute_sphere_force_fields()

        use myarrays_ibm
        use mysending
        use scala3

        implicit none

        integer :: i0,j0,k0,l

        if (allocated(fluidShearForce)) then
                    ! Alessandro: interaction force --------------------------------
            fluidShearForce(:,:,:,:)=0.0
            fluidPressureForce(:,:,:,:)=0.0
            forceCaso(:,:,:)=5
        else
            allocate(fluidPressureForce(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,3))
            allocate(fluidShearForce(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr,3))
            allocate(forceCaso(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
            fluidShearForce(:,:,:,:)=0.0
            fluidPressureForce(:,:,:,:)=0.0
            forceCaso(:,:,:)=5
        end if

        do l=1,num_ib

            i0=indici_CELLE_IB(l,1)
            j0=indici_CELLE_IB(l,2)
            k0=indici_CELLE_IB(l,3)

            fluidShearForce(i0,j0,k0,:)=shear_ib(l,:)

            fluidPressureForce(i0,j0,k0,:)=pressure_ib(l,:)

            forceCaso(i0,j0,k0)=caso_ib(l)

        end do

end subroutine compute_sphere_force_fields

    subroutine pass_geometry(totpart,pposx,pposy,pposz,prad,pmoves) bind (C, name="pass_geometry")

        use mysending

        implicit none

        integer(c_int),intent(in) :: totpart
        real(c_double),intent(in) :: pposx(totpart),pposy(totpart),pposz(totpart),prad(totpart)
        logical(c_bool),intent(in) :: pmoves(totpart)
        integer :: p

        !-------------------------------------------------
        totParticles=totpart

        if (myid==0) then
            write(*,'(A,I0.1,A)') 'Get geometry data of ',totParticles,' particles from C++ '
        end if


        ! deallocate previous arrays if necessary
        if (allocated(sphereIndex)) then
            deallocate(sphereIndex,sphereMoves)
            deallocate(sphereRadius,sphereRadius2)
            deallocate(spherePosition)
            deallocate(sphereShearForce,spherePressureForce)
            deallocate(sphereSurface)
        end if

        allocate(sphereIndex(totParticles),sphereMoves(totParticles))
        allocate(sphereRadius(totParticles),sphereRadius2(totParticles))
        allocate(spherePosition(3,totParticles))
        allocate(sphereShearForce(3,totParticles))
        allocate(spherePressureForce(3,totParticles))
        allocate(sphereSurface(totParticles))

        do p=1,totParticles
            !write(*,*) 'Read',n
            sphereIndex=p
            sphereMoves(p)=pmoves(p)
            spherePosition(1,p)=pposx(p)
            spherePosition(2,p)=pposy(p)
            spherePosition(3,p)=pposz(p)
            sphereRadius(p)=prad(p)
            sphereRadius2(p)=sphereRadius(p)**2.0
            sphereSurface(p)=4.0*3.1428*sphereRadius(p)*sphereRadius(p)
        end do

        if (myid==0) then
            write(*,*) 'Particles:'
            do p=1,totParticles
                write(*,'(ES11.2E3,ES11.2E3,ES11.2E3,ES11.2E3,A,L1)') &
                    sphereRadius(p),spherePosition(1,p),spherePosition(2,p),spherePosition(3,p),' ',sphereMoves(p)
            end do
        end if

        return

    end subroutine pass_geometry

    subroutine pass_forces(sforx,sfory,sforz,pforx,pfory,pforz) bind (C, name="pass_forces")

        use mysending

        implicit none

        real(c_double),intent(out) :: sforx(totParticles),sfory(totParticles),sforz(totParticles)
        real(c_double),intent(out) :: pforx(totParticles),pfory(totParticles),pforz(totParticles)

        integer :: p

        !-------------------------------------------------
        if (myid==0) then
            write(*,'(A,I0.1,A)') 'Send force data of ',totParticles,' particles to C++ '
        end if

        do p=1,totParticles
            sforx(p)=sphereShearForce(1,p)
            sfory(p)=sphereShearForce(2,p)
            sforz(p)=sphereShearForce(3,p)
            pforx(p)=spherePressureForce(1,p)
            pfory(p)=spherePressureForce(2,p)
            pforz(p)=spherePressureForce(3,p)
        end do

        if (myid==0) then
            write(*,*) 'Forces:'
            do p=1,totParticles
                write(*,*) sforx(p),sfory(p),sforz(p),pforx(p),pfory(p),pforz(p)
            end do
        end if

        return

    end subroutine pass_forces


end module particle_module
