!***********************************************************************
subroutine fmassa(bcsi,beta,bzet,ktime)
    !***********************************************************************
    ! compute covariant component of bodyforce
    !
    use iso_c_binding
    use myarrays_velo3, only: u,v,w,rhov
    !
    use scala3
    use mysending, only: myid,nproc,kparasta,kparaend,MPI_REAL_SD
    use mysettings, only: bbx,bby,bbz,rich,latitude
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer,intent(in) :: ktime
    real,intent(out) :: bcsi(n1,n2,kparasta:kparaend)
    real,intent(out) :: beta(n1,n2,kparasta:kparaend)
    real,intent(out) :: bzet(n1,n2,kparasta:kparaend)

    integer ierr,i,j,k
    real ang,smooth
    real :: omega1,omega2,omega3
      
    real,parameter :: earth_rotation = 0.000073
    real,parameter :: pi = acos(-1.)

    !-----------------------------------------------------------------------
    ang = LATITUDE*pi/180.0
    if (abs(LATITUDE)<1.0e-10) then
        omega1 = 0.0
        omega2 = 0.0
        omega3 = 0.0
    else
        omega1 = 0.0
        omega2 = 2.0*earth_rotation*cos(ang)
        omega3 = 2.0*earth_rotation*sin(ang)
    end if
    !-----------------------------------------------------------------------
    ! omega1, omega2, omega3 are the three component of the rotation vector
    ! must be defined in Agenerale.in

    smooth=1.0

    if (ktime<50) then
        smooth=(real(ktime)/50.0)
    end if

    if (myid==0) then
        write(*,*)'smooth',smooth
    end if
      
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx

                bcsi(i,j,k) = bbx+omega2*w(i,j,k)-omega3*v(i,j,k)

                beta(i,j,k) = bby+omega3*u(i,j,k)-omega1*w(i,j,k)-rich*rhov(1,i,j,k) !smooth*

                bzet(i,j,k) = bbz+omega1*v(i,j,k)-omega2* u(i,j,k)

            end do
        end do
    end do

    return
end
