!***********************************************************************
subroutine iniz(f1ve,f2ve,f3ve,bcsi,beta,bzet)
    !***********************************************************************
    !
    !     matrixes allocation and initialization
    !
    !-----------------------------------------------------------------------
    !
    use myarrays_velo3
    use myarrays_density
    use myarrays_moisture
    use mysending
    use scala3

    implicit none
    !      include 'metri3.h'
    !      include 'velo3.h'
    !
    !-----------------------------------------------------------------------
    !     array declaration
    !
    integer i,j,k
    integer kpsta_alloc,kpend_alloc
      
    real f1ve(n1,n2,kparasta:kparaend)
    real f2ve(n1,n2,kparasta:kparaend)
    real f3ve(n1,n2,kparasta:kparaend)
  
    real bcsi(n1,n2,kparasta:kparaend)
    real beta(n1,n2,kparasta:kparaend)
    real bzet(n1,n2,kparasta:kparaend)
    !-----------------------------------------------------------------------

    allocate(akaptV(nscal,0:n1+1,0:n2+1, &
        kparasta-deepl:kparaend+deepr))
    allocate(akapt (nscal,0:n1+1,0:n2+1, &
        kparasta-deepl:kparaend+deepr))
    akaptV = 0.
    akapt = 0.
    if(myid.eq.0)then
        allocate( akapt_piano(nscal,0:n1+1,0:n2+1,n3:n3))
        allocate(akaptV_piano(nscal,0:n1+1,0:n2+1,n3:n3))
        akapt_piano  = 0.
        akaptV_piano = 0.
    elseif(myid.eq. nproc-1)then
        allocate( akapt_piano(nscal,0:n1+1,0:n2+1,1:1))
        allocate(akaptV_piano(nscal,0:n1+1,0:n2+1,1:1))
        akapt_piano  = 0.
        akaptV_piano = 0.
    end if
    !-----------------------------------------------------------------------
    if(imoist==1)then
        allocate(tpotm(1:n2))
        allocate(qm(1:n2))
    end if
    !-----------------------------------------------------------------------
    !
    kpsta_alloc = kparasta
    kpend_alloc = kparaend
    if(myid .eq. 0) kpsta_alloc = kparasta-1
    !
    !-----------------------------------------------------------------------
    allocate(cs1(n2,n3))
    allocate(cs2(n2,n3))
    allocate(cs3(n1,n3))
    allocate(cs4(n1,n3))
    allocate(cs5(n1,n2))
    allocate(cs6(n1,n2))
      
    cs1 = 0.0
    cs2 = 0.0
    cs3 = 0.0
    cs4 = 0.0
    cs5 = 0.0
    cs6 = 0.0
      
    !-----------------------------------------------------------------------
    allocate(uc1_orl(n2,n3))
    allocate(uc2_orl(n2,n3))
    allocate(vc3_orl(n1,n3))
    allocate(vc4_orl(n1,n3))
    allocate(wc5_orl(n1,n2))
    allocate(wc6_orl(n1,n2))
      
    uc1_orl = 0.0
    uc2_orl = 0.0
    vc3_orl = 0.0
    vc4_orl = 0.0
    wc5_orl = 0.0
    wc6_orl = 0.0
      
    !-----------------------------------------------------------------------c
    allocate(rhs(1:n1,1:n2,kparasta:kparaend))

    allocate( u(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
    allocate( v(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
    allocate( w(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
    allocate(fi(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
    allocate(next_prs(0:n1+1,kparasta-deepl:kparaend+deepr))
    allocate(rhov(1:nscal,0:n1+1,0:n2+1, &
        kparasta-deepl:kparaend+deepr))
       
    allocate( uc(0:n1,1:n2,kparasta  :kparaend))
    allocate( vc(1:n1,0:n2,kparasta  :kparaend))
    allocate( wc(1:n1,1:n2,kparasta-1:kparaend))
 

    u = 0.
    v = 0.
    w = 0.
    fi = 0.
    next_prs=0.
    rhov = 0.
    uc = 0.
    vc = 0.
    wc = 0.
    rhs = 0.
       
    if(myid.eq.0)then
        allocate(u_piano(0:n1+1,0:n2+1,n3:n3))
        allocate(v_piano(0:n1+1,0:n2+1,n3:n3))
        allocate(w_piano(0:n1+1,0:n2+1,n3:n3))
        allocate(rhov_piano(1:nscal,0:n1+1,0:n2+1,n3:n3))
	
        u_piano=0.
        v_piano=0.
        w_piano=0.

        allocate(wc_piano(1:n1,1:n2,n3:n3))
		 
        wc_piano=0.

    elseif(myid.eq.nproc-1)then
       
        allocate(u_piano(0:n1+1,0:n2+1,1:1))
        allocate(v_piano(0:n1+1,0:n2+1,1:1))
        allocate(w_piano(0:n1+1,0:n2+1,1:1))
        allocate(rhov_piano(1:nscal,0:n1+1,0:n2+1,1:1))
	
        u_piano=0.
        v_piano=0.
        w_piano=0.

        allocate(wc_piano(1:n1,1:n2,0:0))
		 
        wc_piano=0.
     	
    end if
    !
    !-----------------------------------------------------------------------
    !
    if(myid.eq.0)then
        allocate(gra1_appoggio(n1,n2,n3:n3))
        allocate(gra2_appoggio(n1,n2,n3:n3))
        allocate(gra3_appoggio(n1,n2,n3:n3))
	
        do i=1,n1
            do j=1,n2
                do k=n3,n3
                    gra1_appoggio(i,j,k)=0.
                    gra2_appoggio(i,j,k)=0.
                    gra3_appoggio(i,j,k)=0.
                end do
            end do
        end do
    elseif(myid.eq.nproc-1)then
        allocate(gra1_appoggio(n1,n2,1:1))
        allocate(gra2_appoggio(n1,n2,1:1))
        allocate(gra3_appoggio(n1,n2,1:1))
	
        do i=1,n1
            do j=1,n2
                do k=1,1
                    gra1_appoggio(i,j,k)=0.
                    gra2_appoggio(i,j,k)=0.
                    gra3_appoggio(i,j,k)=0.
                end do
            end do
        end do
    end if

    allocate(cgra1(0:n1,  n2,kparasta-1:kparaend+1))
    allocate(cgra2(  n1,0:n2,kparasta-1:kparaend+1))
    allocate(cgra3(  n1,  n2,kparasta-1:kparaend+1))

    allocate(gra1(n1,n2,kparasta-1:kparaend+1))
    allocate(gra2(n1,n2,kparasta-1:kparaend+1))
    allocate(gra3(n1,n2,kparasta-1:kparaend+1))

    do k=kparasta-1,kparaend+1
        do j=1,n2
            do i=0,n1
                cgra1(i,j,k) = 0.
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2
            do i=1,n1
                cgra2(i,j,k) = 0.
            end do
        end do
    end do
              
    do k=kparasta-1,kparaend+1
        do j=1,n2
            do i=1,n1
                cgra3(i,j,k) = 0.
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=1,n2
            do i=1,n1
                !
                gra1(i,j,k)=0.
                gra2(i,j,k)=0.
                gra3(i,j,k)=0.
            !
            end do
        end do
    end do

    ! matrix defined at the centroid
    !
    !       do 1 i=1,n1
    !       do 1 j=1,n2
    !       do 1 k=1,n3
    !
    !       giac(i,j,k)=0.
    !
    !  1    continue


    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                !
                f1ve(i,j,k)=0.
                f2ve(i,j,k)=0.
                f3ve(i,j,k)=0.
                bcsi(i,j,k)=0.
                beta(i,j,k)=0.
                bzet(i,j,k)=0.
            !
            end do
        end do
    end do

    !
    ! matrix defined on points
    !
    !       do 2 i=0,n1
    !       do 2 j=0,n2
    !       do 2 k=0,n3
    !
    !       x(i,j,k)=0.
    !       y(i,j,k)=0.
    !       z(i,j,k)=0.
    !
    !  2    continue
    !
    ! matrix on face with constant csi
    !
    do j=1,n2
        do k=1,n3
            !
            !       cs1(j,k)=0.
            !       cs2(j,k)=0.
            !
            do i=0,n1
            !
            !       csx(i,j,k)=0.
            !       csy(i,j,k)=0.
            !       csz(i,j,k)=0.
            !       g11(i,j,k)=0.
            !       g12(i,j,k)=0.
            !       g13(i,j,k)=0.
            !
            !       uc(i,j,k)=0.
            !       f1(i,j,k)=0.
            !
            END DO
        END DO
    END DO
    !
    ! matrix on face with constant eta
    !
    do i=1,n1
        do k=1,n3
            !
            !       cs3(i,k)=0.
            !       cs4(i,k)=0.
            !
            do j=0,n2
            !
            !       etx(i,j,k)=0.
            !       ety(i,j,k)=0.
            !       etz(i,j,k)=0.
            !       g21(i,j,k)=0.
            !       g22(i,j,k)=0.
            !       g23(i,j,k)=0.
            !
            !       vc(i,j,k)=0.
            !       f2(i,j,k)=0.
            !
            END DO
        END DO
    END DO
    !
    ! matrix on face with constant zeta
    !
    do i=1,n1
        do j=1,n2
            !
            !       cs5(i,j)=0.
            !       cs6(i,j)=0.
            !
            do k=0,n3
            !
            !       ztx(i,j,k)=0.
            !       zty(i,j,k)=0.
            !       ztz(i,j,k)=0.
            !       g31(i,j,k)=0.
            !       g32(i,j,k)=0.
            !       g33(i,j,k)=0.
            !
            !       wc(i,j,k)=0.
            !       f3(i,j,k)=0.
            !
            END DO
        END DO
    END DO
    !
    ! matrix defined in the field and out
    !
    do  i=0,n1+1
        do  j=0,n2+1
            do  k=kparasta-1,kparaend+1
                delu(i,j,k)=0.
                delv(i,j,k)=0.
                delw(i,j,k)=0.
            end do
        end do
    end do
    !

    do i=0,n1+1
        do j=0,n2+1
            do k=0,n3+1
            !
            !       rhov(:,i,j,k)=0.0
            !       u(i,j,k)=0.
            !       v(i,j,k)=0.
            !       w(i,j,k)=0.
            !       fi(i,j,k)=0.
            !
            !       annit(i,j,k)=0.
            !
            END DO
        END DO
    END DO
    !
    return
end
