!***********************************************************************
subroutine ada_rho(ktime,rven,isc,rho,kdeg)
    !***********************************************************************
    ! compute explicit term for scalar equation
    !
    use mysending
    use myarrays_metri3
    use myarrays_velo3
    !
    use scala3
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer i,j,k,ktime,isc
    real rven(nscal,n1,n2,kparasta:kparaend)
    real rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
    real kdeg
    !
    integer ierr
    integer ncolperproc,m
    !-----------------------------------------------------------------------
    ! compute convective term + explicit diffusive term
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                !
                rhs(i,j,k)=f1(i,j,k)-f1(i-1,j,k)+ &
                    f2(i,j,k)-f2(i,j-1,k)+&
                    f3(i,j,k)-f3(i,j,k-1) &
                    -kdeg*rho(i,j,k)*giac(i,j,k)  !species decay, kdeg=0 default
            !
            enddo
        enddo
    enddo
    !
    ! compute Adam-Bashforth part for scalar equation
    !
    if (ktime.eq.1) then
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    !
                    rven(isc,i,j,k)=rhs(i,j,k)
                !
                end do
            end do
        end do
    !
    end if
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                !
                rhs(i,j,k)=dt*(1.5*rhs(i,j,k)-.5*rven(isc,i,j,k)) &
                    /giac(i,j,k)
            !
            end do
        end do
    end do
    !
    !
    ! compute explict part at step n-1
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                !
                rven(isc,i,j,k)=f1(i,j,k)-f1(i-1,j,k)+ &
                    f2(i,j,k)-f2(i,j-1,k)+ &
                    f3(i,j,k)-f3(i,j,k-1) &
                    -kdeg*rho(i,j,k)*giac(i,j,k)  !species decay, kdeg=0 default
            !
            enddo
        enddo
    enddo
    !
    return
end
