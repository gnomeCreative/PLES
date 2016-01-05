!***********************************************************************
subroutine tridag(aa,bb,cc,ff,n)
    !***********************************************************************
    ! solution for a tridiagonal system
    implicit none
    !-----------------------------------------------------------------------
    !     array declaration
    integer n1,n2,k,k1,n,n1n
    real aa(*),bb(*),cc(*),ff(*)
    !-----------------------------------------------------------------------
    n1=1
    bb(n1)=1./bb(n1)
    aa(n1)=ff(n1)*bb(n1)
    n2=n1+1
    n1n=n1+n
    do k=n2,n
        k1=k-1
        cc(k1)=cc(k1)*bb(k1)
        bb(k)=bb(k)-aa(k)*cc(k1)
        bb(k)=1./bb(k)
        aa(k)=(ff(k)-aa(k)*aa(k1))*bb(k)
    end do
    !
    ! back substitution
    !
    ff(n)=aa(n)
    do k1=n2,n
        k=n1n-k1
        ff(k)=aa(k)-cc(k)*ff(k+1)
    end do
    return
end
