subroutine rhs1_rho()

    ! compute the right hand side for the eq.

    use myarrays_velo3
    use myarrays_metri3

    use scala3

    implicit none

    !-----------------------------------------------------------------------
    integer :: i,j,k
    real :: espl
    !-----------------------------------------------------------------------


    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1

                espl=f1(i,j,k)-f1(i-1,j,k) &
                    +f2(i,j,k)-f2(i,j-1,k) &
                    +f3(i,j,k)-f3(i,j,k-1)

                rhs(i,j,k)=rhs(i,j,k)+dt*espl/giac(i,j,k)

            end do
        end do
    end do

    return

end
