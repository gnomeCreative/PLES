module wind_module

    !     for wavebreaking
    !-----------------------------------------------------------------------
    use iso_c_binding
    !
    real :: betaw
    real :: a_wind
    ! Giulia modificavento: variabili nuove
    real, allocatable :: tauu_att(:,:)
    real, allocatable :: tauw_att(:,:)
    real, allocatable :: v_att(:,:)
    logical :: cf_exists
    real, allocatable :: cf(:,:)
    real :: varcf
    real :: c10
contains

    subroutine leggivento

        ! wind from file is read, file written as wind amplitude and angle
        ! in radians
        use scala3, only: n1,kparasta,kparaend

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k
        !integer kparasta,kparaend,myid,nproc
        integer ktime
        real rhoa,rhow,pi,ang

        real cf_i(1:n1,kparasta:kparaend)
        !-----------------------------------------------------------------------
        !
        pi = acos(-1.)

        rhoa = 1.300  ! air density
        rhow = 1000.  ! water density

        a_wind = 3.45
        ang= 45.
        betaw =(-ang+270)*pi/180. ! wind angle, meteorological convention: clockwise, 0° N axis.
        !if(myid.eq.0)write(*,*)'wind',a_wind,ang

        !------------------------------Wind computation-------

        do k=kparasta,kparaend
            do i=1,n1
                cf_i(i,k) = cf(i,k)
            enddo
        enddo

        !     read file and compute x and z component
        !      open(500,file='dativento.dat',status='unknown')
        !      do k=1,n3
        !      do i=1,n1
                !read(500,*)a,betaw
        !    betaw = betaw !(-betaw+270)*pi/180. ! wind angle, meteorological convention: clockwise, 0° N axis.
        !        if (k.ge.kparasta.and.k.le.kparaend) then
        !         u_wind(i,k) =  a*cos(betaw)
        !         w_wind(i,k) =  a*sin(betaw)
        !         end if
        !       end do
        !      end do
        !      close(500)

        !     compute friction velocity

        ! Giulia modificavento:
        ! Giulia: correction wind stress
        !Wu: Wind-stress Coefficients Over Sea Surface...
        ! C_10=(.8+.065*U_10)*.001
        ! C_10=|tau|/rho_a U_10^2

        c10=(.8+.065*a_wind)*.001
        do k=kparasta-1,kparaend+1
            do i=1,n1
                tauu_att(i,k) = a_wind**2*cos(betaw)*c10*rhoa
                tauw_att(i,k) = a_wind**2*sin(betaw)*c10*rhoa
            !        u_att(i,k) = a*cos(betaw)*(sqrt(c10*rhoa/rhow))
            !        w_att(i,k) = a*sin(betaw)*(sqrt(c10*rhoa/rhow))
            !        c10  = (.8+.065*abs(u_wind(i,k)) )*.001
            !        write(*,*) myid,i,k, c10, u_wind
            !        u_att(i,k) = u_wind(i,k)*(sqrt(c10*rhoa/rhow))
            !        c10  = (.8+.065*abs(w_wind(i,k)) )*.001
            !        w_att(i,k) = w_wind(i,k)*(sqrt(c10*rhoa/rhow))
            end do
        end do

    end subroutine leggivento

end module wind_module
