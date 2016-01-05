!***********************************************************************
subroutine resid_sndrcv_cart(n,i1,j1,k1,jxc,jyc,jzc,pr,rh1,rh, &
    cs1,cs2,cs3,cs4,cs5,cs6,       &
    r11,r12,r13,r21,r22,r23,r31,r32,r33, &
    i_dx,i_sn,i_sp,i_st,i_av,i_in, &
    kstamg,kendmg, &
    bodyforce,bodypressure,kparasta,kparaend)
    !***********************************************************************
    ! compute residuals on current grid
    !
    use scala3
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer kstamg(4),kendmg(4)
    integer kparasta,kparaend
    integer i,j,k,n,i1,j1,k1,ipot
    integer i_dx(i1,j1,kstamg(n):kendmg(n)) !k1)
    integer i_sn(i1,j1,kstamg(n):kendmg(n)) !k1)
    integer i_sp(i1,j1,kstamg(n):kendmg(n)) !k1)
    integer i_st(i1,j1,kstamg(n):kendmg(n)) !k1)
    integer i_av(i1,j1,kstamg(n):kendmg(n)) !k1)
    integer i_in(i1,j1,kstamg(n):kendmg(n)) !k1)
    integer jxc(0:4),jyc(0:4),jzc(0:4)
    !
    real res_av,res_ind,res_sot,res_sop,res_sn,res_dx
    real pot,inv_pot,inv_dt,coef
    real r11(0:i1,  j1,  kstamg(n):kendmg(n))
    real r12(0:i1,  j1,  kstamg(n):kendmg(n))
    real r13(0:i1,  j1,  kstamg(n):kendmg(n))
    real r21(i1  ,0:j1,  kstamg(n):kendmg(n))
    real r22(i1  ,0:j1,  kstamg(n):kendmg(n))
    real r23(i1  ,0:j1,  kstamg(n):kendmg(n))
    real r31(i1  ,  j1,kstamg(n)-1:kendmg(n))
    real r32(i1  ,  j1,kstamg(n)-1:kendmg(n))
    real r33(i1  ,  j1,kstamg(n)-1:kendmg(n))
    real pr(0:i1+1,0:j1+1,kstamg(n)-1:kendmg(n)+1) !0:k1+1)
    real rh(i1,j1,kstamg(n):kendmg(n))  !k1)
    real rh1(i1,j1,kstamg(n):kendmg(n)) !k1)
      
    real cs1(n2,n3),cs2(n2,n3)
    real cs3(n1,n3),cs4(n1,n3)
    real cs5(n1,n2),cs6(n1,n2)
    real an

    !
    integer ierr,myid,nproc,status
    integer ncolperproc,m
    !
    !      integer tipo(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
    integer bodyforce,ibodyforce
    integer bodypressure,ibodypressure
    integer ilivello,itipo_ib,itipo_solida
    integer ibloop,solidaloop

    !-----------------------------------------------------------------------
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    !-----------------------------------------------------------------------

    ipot=2**(n-1)
    pot=float(ipot)

    inv_pot = 1./pot

      
    coef = 0.25

    if(n.eq.1)then
        an=1.
    else
        an=0.
    end if

    inv_dt  = an/dt
    !-----------------------------------------------------------------------
    !
    do k=kstamg(n),kendmg(n)
        do j=1,jyc(n)
            do i=1,jxc(n)
                !
                !      residual right
                res_dx=i_dx(i,j,k)*(  &
                    r11(i,j,k)*(pr(i+1,j,k)-pr(i,j,k)))*inv_pot &
                    +(1.-i_dx(i,j,k))*cs2(j,k)*inv_dt
    
                !
                !      left residual
                res_sn=i_sn(i,j,k)*( &
                    r11(i-1,j,k)*(pr(i,j,k)-pr(i-1,j,k)))*inv_pot &
                    +(1.-i_sn(i,j,k))*cs1(j,k)*inv_dt
                !
                !      upper residual
                res_sop=i_sp(i,j,k)*( &
                    r22(i,j,k)*(pr(i,j+1,k)-pr(i,j,k)))*inv_pot &
                    +(1.-i_sp(i,j,k))*cs4(i,k)*inv_dt
                !
                !     bottom residual
                res_sot=i_st(i,j,k)*( &
                    r22(i,j-1,k)*(pr(i,j,k)-pr(i,j-1,k)))*inv_pot &
                    +(1.-i_st(i,j,k))*cs3(i,k)*inv_dt
                !
                !     front residual
                res_av=i_av(i,j,k)*( &
                    r33(i,j,k)*(pr(i,j,k+1)-pr(i,j,k)))*inv_pot &
                    +(1.-i_av(i,j,k))*cs6(i,j)*inv_dt
                !
                !     back residual
                res_ind=i_in(i,j,k)*( &
                    r33(i,j,k-1)*(pr(i,j,k)-pr(i,j,k-1)))*inv_pot &
                    +(1.-i_in(i,j,k))*cs5(i,j)*inv_dt
                !
                rh1(i,j,k)= &
                    (res_dx  - res_sn  +  &
                    res_sop - res_sot + &
                    res_av  - res_ind )*inv_pot - rh(i,j,k)
            !
            enddo
        enddo
    enddo
    !-----------------------------------------------------------------------


    return
end
