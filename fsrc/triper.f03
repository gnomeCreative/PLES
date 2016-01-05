!***********************************************************************
subroutine triper(aa,bb,cc,ff,n)
    !***********************************************************************
    ! subeoutine from book ...
    !
    ! aa,bb,cc coefficent vectors: aa(i)=a(1,i)
    !                              bb(i)=a(2,i)
    ! the periodic term up/right is in aa(N1)
    ! the periodic term bottom/left is in cc(N+1)
    ! ff is the vector with the known terms, then transformed
    ! in the solution vector.
    ! aa,bb,cc,ff,gam2 are vectors with dimension N+2 to be specified in the
    ! calling program
    !
    use scala3

    implicit none
    !-----------------------------------------------------------------------
    !     array declaration
    integer k,k1,k2,n,np1,np2,n1n
    real zaa
    real aa(*),bb(*),cc(*),ff(*)
    real gam2(n1+n2+n3)
    !-----------------------------------------------------------------------
    np1=1
    bb(np1)=1./bb(np1)
    gam2(np1)=-aa(np1)*bb(np1)
    aa(np1)=ff(np1)*bb(np1)
    np2=np1+1
    n1n=np1+n
    do k=np2,n
        k1=k-1
        cc(k1)=cc(k1)*bb(k1)
        bb(k)=bb(k)-aa(k)*cc(k1)
        bb(k)=1./bb(k)
        gam2(k)=-aa(k)*gam2(k1)*bb(k)
        aa(k)=(ff(k)-aa(k)*aa(k1))*bb(k)
    end do
    gam2(n)=gam2(n)-cc(n)*bb(n)
    !
    ! back substitution
    !
    ff(n)=aa(n)
    bb(n)=gam2(n)
    do k1=np2,n
        k=n1n-k1
        k2=k+1
        ff(k)=aa(k)-cc(k)*ff(k2)
        bb(k)=gam2(k)-cc(k)*bb(k2)
    end do
    !
    k1=n+1
    zaa=ff(k1)-cc(k1)*ff(np1)-aa(k1)*ff(n)
    zaa=zaa/(bb(k1)+aa(k1)*bb(n)+cc(k1)*bb(np1))
    ff(k1)=zaa
    do k=np1,n
        ff(k)=ff(k)+bb(k)*zaa
    end do    !
    ff(n+2)=ff(np1)
    return
end
