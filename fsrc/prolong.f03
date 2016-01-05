!***********************************************************************
subroutine prolong(n,i1,j1,k1,i2,j2,k2,jxc,jyc,jzc,prf,pr, &
    kstamg,kendmg,rightpe,leftpe, &
    tagls,taglr,tagrs,tagrr)
    !***********************************************************************
    ! prolongation operation from coarse grid to fine one with bilinear
    ! interpolation
    !
    use tipologia
    !
    use mpi

    implicit none
    !

    !-----------------------------------------------------------------------
    !     array declaration
    integer kstamg(4),kendmg(4)
    integer i,j,k,i1,j1,k1,i2,j2,k2,n,if,jf,kf
    integer jxc(0:4),jyc(0:4),jzc(0:4)
    real a1,a2,a3,a4,a5,a6,a7,a8
    real inv_64
    real prf(0:i1+1,0:j1+1,kstamg(n  )-1:kendmg(n  )+1) !0:k1+1)
    real  pr(0:i2+1,0:j2+1,kstamg(n+1)-1:kendmg(n+1)+1) !0:k2+1)

    integer kstamg0(4)
    integer nproc,myid,ierr
    integer leftpe,rightpe,ktime
    integer leftpem,rightpem,status(MPI_STATUS_SIZE)
    integer tagls,taglr,tagrs,tagrr
    integer lll
    integer plantypef
    integer req1,req2,req3,req4
    integer istatus(MPI_STATUS_SIZE)
    !-----------------------------------------------------------------------
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    !-----------------------------------------------------------------------
    !
    !     prf pressure on fine grid
    !     pr  pressure on coarse grid
    !
    call MPI_TYPE_VECTOR(jyc(n)+1,jxc(n)+1,jxc(n)+2, &
        MPI_REAL_SD,plantypef,ierr)
    call MPI_TYPE_COMMIT(plantypef,ierr)

    if(myid.eq.0)then
        leftpem=MPI_PROC_NULL
        rightpem=rightpe
        kstamg0(n+1)=0
    else if(myid.eq.nproc-1)then
        leftpem=leftpe
        rightpem=MPI_PROC_NULL
        kstamg0(n+1)=kstamg(n+1)
    else if((myid.ne.0).and.(myid.ne.nproc-1))then
        leftpem=leftpe
        rightpem=rightpe
        kstamg0(n+1)=kstamg(n+1)
    endif


    inv_64 = 1./64.

    do k=kstamg0(n+1),kendmg(n+1)
        do j=0,jyc(n+1)
            do i=0,jxc(n+1)
                !
                if=2*i
                jf=2*j
                kf=2*k
                !
                a1=27.
                a2=9.
                a3=9.
                a4=9.
                a5=3.
                a6=3.
                a7=3.
                a8=1.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64


                if=2*i+1
                jf=2*j
                kf=2*k
                !
                a1=9.
                a2=27.
                a3=3.
                a4=3.
                a5=9.
                a6=9.
                a7=1.
                a8=3.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64
                !
                !
                if=2*i
                jf=2*j+1
                kf=2*k
                !
                a1=9.
                a2=3.
                a3=27.
                a4=3.
                a5=9.
                a6=1.
                a7=9.
                a8=3.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64
                !
                !
                if=2*i
                jf=2*j
                kf=2*k+1
                !
                a1=9.
                a2=3.
                a3=3.
                a4=27.
                a5=1.
                a6=9.
                a7=9.
                a8=3.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64
                !
                !
                if=2*i+1
                jf=2*j+1
                kf=2*k
                !
                a1=3.
                a2=9.
                a3=9.
                a4=1.
                a5=27.
                a6=3.
                a7=3.
                a8=9.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64
                !
                !
                if=2*i+1
                jf=2*j
                kf=2*k+1
                !
                a1=3.
                a2=9.
                a3=1.
                a4=9.
                a5=3.
                a6=27.
                a7=3.
                a8=9.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64
                !
                !
                if=2*i+1
                jf=2*j+1
                kf=2*k+1
                !
                a1=1.
                a2=3.
                a3=3.
                a4=3.
                a5=9.
                a6=9.
                a7=9.
                a8=27.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64
                !
                !
                if=2*i
                jf=2*j+1
                kf=2*k+1
                !
                a1=3.
                a2=1.
                a3=9.
                a4=9.
                a5=3.
                a6=3.
                a7=27.
                a8=9.
                !
                prf(if,jf,kf)=prf(if,jf,kf)+ &
                    (  a1*pr(i  ,j  ,k  )+ &
                    a2*pr(i+1,j  ,k  )+ &
                    a3*pr(i  ,j+1,k  )+ &
                    a4*pr(i  ,j  ,k+1)+ &
                    a5*pr(i+1,j+1,k  )+ &
                    a6*pr(i+1,j  ,k+1)+ &
                    a7*pr(i  ,j+1,k+1)+ &
                    a8*pr(i+1,j+1,k+1) )*inv_64
            !
            enddo
        enddo
    enddo
    !
    ! before exchange I need to collocate in the proper position
    ! the plane computed by each procs which belong to its right
    !
    if(rightpem /= MPI_PROC_NULL) then
        call MPI_SSEND(prf(0,0,kendmg(n)+1),1, &
            plantypef,rightpem,tagrs, &
            MPI_COMM_WORLD,ierr)
    endif
    if(leftpem /= MPI_PROC_NULL) then
        call MPI_RECV(prf(0,0,kstamg(n)),1, &
            plantypef,leftpem,taglr, &
            MPI_COMM_WORLD,istatus,ierr)
    endif

    if(rightpem /= MPI_PROC_NULL) then
    !      call MPI_WAIT(req1,istatus,ierr)
    endif
    if(leftpem /= MPI_PROC_NULL) then
    !      call MPI_WAIT(req2,istatus,ierr)
    endif

    !

    if(rightpem /= MPI_PROC_NULL) then
        call MPI_SSEND(prf(0,0,kendmg(n)),1, &
            plantypef,rightpem,tagrs, &
            MPI_COMM_WORLD,ierr)
    endif
    if(leftpem /= MPI_PROC_NULL) then
        call MPI_RECV(prf(0,0,kstamg(n)-1),1, &
            plantypef,leftpem,taglr, &
            MPI_COMM_WORLD,istatus,ierr)
    endif
                                                                                     
    if(rightpem /= MPI_PROC_NULL) then
    !      call MPI_WAIT(req3,istatus,ierr)
    endif
    if(leftpem /= MPI_PROC_NULL) then
    !      call MPI_WAIT(req4,istatus,ierr)
    endif


    call MPI_TYPE_FREE(plantypef,ierr)


    return
end
