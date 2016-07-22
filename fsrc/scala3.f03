module scala3

    use,intrinsic :: iso_c_binding
    ! SONO DEFINITI:
    ! i parametri n1,n2,n3 massimi di dimensione
    ! il numero di scalari che si considerano nscal
    ! il numero di celle in x,y,z
    !

    ! -----------------------------------------------------
    integer(kind=c_int),bind(C) :: n1,n2,n3
    integer(kind=c_int),bind(C) :: nscal
    logical(kind=c_bool),bind(C) :: potenziale
    real(kind=c_double),bind(C) :: dt,dt_start,re
    ! -----------------------------------------------------
    integer :: kparasta,kparaend
    integer :: deepl,deepr
    ! -----------------------------------------------------

    ! IMPORTANTE se doppia compilare con -qrealsize=8 o flag equivalente
    ! per ridefinire i real in real*8

contains

    function cross(a,b)

        real,dimension(3) :: cross
        real,dimension(3),intent(in) :: a,b

        cross(1)=a(2)*b(3)-a(3)*b(2)
        cross(2)=a(3)*b(1)-a(1)*b(3)
        cross(3)=a(1)*b(2)-a(2)*b(1)

    end function cross

    function dot(a,b)

        real :: dot
        real,dimension(3),intent(in) :: a,b

        dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

    end function dot

end module scala3

