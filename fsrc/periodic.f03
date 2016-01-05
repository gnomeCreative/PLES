!***********************************************************************
subroutine periodic(r1,r2,r3,myid,nproc,kparasta,kparaend)
    !***********************************************************************
    ! set periodic boundary condition on model matrix
    !
    use scala3
    use period
    !
    use mpi

    implicit none
    !
    !-----------------------------------------------------------------------
    !     array declaration
    integer i,j,k
    integer ii,jj,kk
    integer myid,nproc,kparasta,kparaend

    real r1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
    real r2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
    real r3(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
    !-----------------------------------------------------------------------
    !
    !      periodicity on csi
    !
    do ii=1,1-ip
        do k=kparasta,kparaend
            do j=1,jy
                r1(0,j,k)=r1(jx,j,k)
                r2(0,j,k)=r2(jx,j,k)
                r3(0,j,k)=r3(jx,j,k)
                !
                r1(jx+1,j,k)=r1(1,j,k)
                r2(jx+1,j,k)=r2(1,j,k)
                r3(jx+1,j,k)=r3(1,j,k)
            enddo
        enddo
    enddo
    !
    !      periodicity on eta
    !
    do jj=1,1-jp
        do k=kparasta,kparaend
            do i=1,jx
                r1(i,0,k)=r1(i,jy,k)
                r2(i,0,k)=r2(i,jy,k)
                r3(i,0,k)=r3(i,jy,k)
                !
                r1(i,jy+1,k)=r1(i,1,k)
                r2(i,jy+1,k)=r2(i,1,k)
                r3(i,jy+1,k)=r3(i,1,k)
            enddo
        enddo
    enddo
    !
    !      periodicity on zeta made in a next step
    !
    return
end
