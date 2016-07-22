subroutine ada_rho(ktime,rven,isc,rho,kdeg)
    ! compute explicit term for scalar equation

    use myarrays_metri3, only: f1,f2,f3,giac
    use myarrays_velo3, only: rhs
    use scala3

    implicit none

    !-----------------------------------------------------------------------
    integer i,j,k,ktime,isc
    real rven(nscal,n1,n2,kparasta:kparaend)
    real rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
    real kdeg
    !-----------------------------------------------------------------------

    ! compute convective term + explicit diffusive term
    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                rhs(i,j,k)=f1(i,j,k)-f1(i-1,j,k)+f2(i,j,k)-f2(i,j-1,k)+&
                    f3(i,j,k)-f3(i,j,k-1)-kdeg*rho(i,j,k)*giac(i,j,k)  !species decay, kdeg=0 default
            end do
        end do
    end do

    ! compute Adam-Bashforth part for scalar equation
    if (ktime==1) then
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    rven(isc,i,j,k)=rhs(i,j,k)
                end do
            end do
        end do
    end if

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                rhs(i,j,k)=dt*(1.5*rhs(i,j,k)-.5*rven(isc,i,j,k))/giac(i,j,k)
            end do
        end do
    end do

    ! compute explict part at step n-1
    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                rven(isc,i,j,k)=f1(i,j,k)-f1(i-1,j,k)+f2(i,j,k)-f2(i,j-1,k)+ &
                    f3(i,j,k)-f3(i,j,k-1)-kdeg*rho(i,j,k)*giac(i,j,k)  !species decay, kdeg=0 default
            end do
        end do
    end do

    return
end
