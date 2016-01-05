!***********************************************************************
subroutine restrict_sndrcv(n,i1,j1,k1,i2,j2,k2, &
    jxc,jyc,jzc,rh,rh1, &
    kstamg,kendmg)
    !***********************************************************************
    ! compute residual on lower level grid
    ! restriction operation with average value
    !
    use mpi

    implicit none
    !
    !-----------------------------------------------------------------------
    !     array declaration
    integer ierr,myid,nproc,status
    integer kstamg(4),kendmg(4),m
    !
    integer i,j,k,id,jd,kd,n,i1,j1,k1,i2,j2,k2
    integer jxc(0:4),jyc(0:4),jzc(0:4)
    real  rh(i1,j1,kstamg(n):kendmg(n)) !k1)
    real rh1(i2,j2,kstamg(n+1):kendmg(n+1)) !k2)

    real, allocatable :: rh1col(:),rh1tot(:)
      
    real inv_8
    !-----------------------------------------------------------------------
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    !-----------------------------------------------------------------------
    inv_8 = 1./8.

    do k=kstamg(n+1),kendmg(n+1)
        do j=1,jyc(n+1)
            do i=1,jxc(n+1)
                !
                id=2*i
                jd=2*j
                kd=2*k
                !
                rh1(i,j,k)=rh(id  ,jd  ,kd  ) &
                    +rh(id-1,jd  ,kd  ) &
                    +rh(id  ,jd-1,kd  ) &
                    +rh(id-1,jd-1,kd  ) &
                    +rh(id  ,jd  ,kd-1) &
                    +rh(id-1,jd  ,kd-1) &
                    +rh(id  ,jd-1,kd-1) &
                    +rh(id-1,jd-1,kd-1)
                !
                rh1(i,j,k)=-rh1(i,j,k)*inv_8
            !
            enddo
        enddo
    enddo
    !
    return
end
