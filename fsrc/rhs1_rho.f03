!***********************************************************************
subroutine rhs1_rho(kparasta,kparaend)
    !***********************************************************************
    ! compute the right hand side for the eq.
    !
    use myarrays_velo3
    use myarrays_metri3
    !
    use scala3
    !
    use mpi

    implicit none
    !-----------------------------------------------------------------------
    !     array declaration
    integer ierr,myid,nproc
    integer ncolperproc,kparasta,kparaend,m
    !
    integer i,j,k
    real espl
    !-----------------------------------------------------------------------

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                !
                espl=f1(i,j,k)-f1(i-1,j,k) &
                    +f2(i,j,k)-f2(i,j-1,k) &
                    +f3(i,j,k)-f3(i,j,k-1)
                !
                rhs(i,j,k)=rhs(i,j,k)+dt*espl/giac(i,j,k)
            !
            enddo
        enddo
    enddo
    !
    return
end
