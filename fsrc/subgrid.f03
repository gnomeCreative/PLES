module subgrid

    !
    ! common per variabili subgrid cartesiani eddy-viscosity
    ! viene dichiarata anche la media, max e rms della costante
    ! del modello - tutti valori mediati sui piani di omogeneità xz
    !
    real,allocatable :: sub(:),sub11(:),sub22(:),sub33(:)
    real,allocatable ::         sub12(:),sub13(:),sub23(:)
    real,allocatable :: sus(:),sus11(:),sus22(:),sus33(:)
    real,allocatable ::         sus12(:),sus13(:),sus23(:)
    real,allocatable :: subrho11(:),subrho22(:),subrho33(:)
    real,allocatable :: susrho11(:),susrho22(:),susrho33(:)
    real,allocatable :: c11(:),c22(:),c33(:)

contains

    subroutine initialize_subgrid()

        use scala3, only: n2

        implicit none

        allocate(sub(n2),sub11(n2),sub22(n2),sub33(n2))
        allocate(sub12(n2),sub13(n2),sub23(n2))
        allocate(sus(n2),sus11(n2),sus22(n2),sus33(n2))
        allocate(sus12(n2),sus13(n2),sus23(n2))
        allocate(subrho11(n2),subrho22(n2),subrho33(n2))
        allocate(susrho11(n2),susrho22(n2),susrho33(n2))
        allocate(c11(n2),c22(n2),c33(n2))

    end subroutine initialize_subgrid

end module subgrid
