module ts_module

    use turbo2_data
    use turbo3bis
    use myarrays_wallmodel

    use mysending
    use myarrays_metri3
    use myarrays_velo3
    use myarrays_density
    !
    use scala3
    use subgrid
    use period
    use tipologia
    use velpar
    !
    use mpi

    implicit none
    private
    !
    public :: turbo_statico
contains
    !***********************************************************************
    subroutine turbo_statico(ktime,i_print, &
        in_dx1,in_sn1,in_sp1, &
        in_st1,in_av1,in_in1, &
        isotropo,kpstamg,kpendmg)
        !***********************************************************************
        ! compute the eddy viscosity and diffusivity with Smagorinsky
        ! with fixed constant
        !
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,i_print
        integer lagr
        integer ktime,kper,lll,m

        integer kpstamg(0:4),kpendmg(0:4)
        integer in_dx1(n1  ,n2  ,kpstamg(1):kpendmg(1)) !    (n1,n2,n3)
        integer in_sn1(n1  ,n2  ,kpstamg(1):kpendmg(1)) !    (n1,n2,n3)
        integer in_sp1(n1  ,n2  ,kpstamg(1):kpendmg(1)) !    (n1,n2,n3)
        integer in_st1(n1  ,n2  ,kpstamg(1):kpendmg(1)) !    (n1,n2,n3)
        integer in_av1(n1  ,n2  ,kpstamg(1):kpendmg(1)) !    (n1,n2,n3)
        integer in_in1(n1  ,n2  ,kpstamg(1):kpendmg(1)) !    (n1,n2,n3)
        !
        integer ierr
        integer ncolperproc
        integer status(MPI_STATUS_SIZE)
        integer req1,req2,req3,req4
        integer istatus(MPI_STATUS_SIZE)
        integer kparastal,kparaendl

        integer debugg,ttt
        !
        real sbuff((n1+2)*(n2+2)*40)
        real rbuff((n1+2)*(n2+2)*40)
        !
        real somma,delt
        real grg
        !
        real xcsi,ycsi,zcsi
        real xeta,yeta,zeta
        real xzet,yzet,zzet
        real dudx,dudy,dudz
        real dvdx,dvdy,dvdz
        real dwdx,dwdy,dwdz
        real drhodx,drhody,drhodz
        real u1,u2,u3,u4,u5,u6
        real v1,v2,v3,v4,v5,v6
        real w1,w2,w3,w4,w5,w6
        real s1,s2,s3,s4,s5,s6
        real rho1,rho2,rho3,rho4,rho5,rho6
        real,allocatable :: s1_loc(:),s2_loc(:),s3_loc(:)
        real,allocatable :: s4_loc(:),s5_loc(:),s6_loc(:)
        real,allocatable :: s1_tot(:),s2_tot(:),s3_tot(:)
        real,allocatable :: s4_tot(:),s5_tot(:),s6_tot(:)

        real,allocatable :: s1rho_loc(:),s2rho_loc(:),s3rho_loc(:)
        real,allocatable :: s1rho_tot(:),s2rho_tot(:),s3rho_tot(:)

        real,allocatable :: grg_loc(:)
        real,allocatable :: grg_tot(:)

        integer pp
        !
        real tp,ti
        double precision costante
        real lhor,lhorx,lhorz,lver
        real costanteV,costanteH
        !
        integer iso,isotropo

        integer iii
        integer isc,cont

        real, allocatable :: P_akapt(:,:,:),P_akaptV(:,:,:)
        !-----------------------------------------------------------------------
        !     for debugging
        lagr=1.
        debugg = 200
        if (ktime.eq.1) then
            ttt = 0
        end if
        if (ktime.eq.2) then
            ttt = 10000
        end if
        !-----------------------------------------------------------------------
        ! velocity extrapolation on sides i=0 and i=jx+1
        !
        !     not periodic
        do kper=1,ip
            !
            do k=kparasta,kparaend
                do j=1,jy
                    !
                    usn(j,k)=  2.*  u(0,j,k)-  u(1,j,k)
                    vsn(j,k)=  2.*  v(0,j,k)-  v(1,j,k)
                    wsn(j,k)=  2.*  w(0,j,k)-  w(1,j,k)
                    rhosn(1,j,k)=  2.*rhov(1,0,j,k)-rhov(1,1,j,k)
                    !
                    udx(j,k)=  2.*  u(jx+1,j,k)-  u(jx,j,k)
                    vdx(j,k)=  2.*  v(jx+1,j,k)-  v(jx,j,k)
                    wdx(j,k)=  2.*  w(jx+1,j,k)-  w(jx,j,k)
                    rhodx(1,j,k)=  2.*rhov(1,jx+1,j,k)-rhov(1,jx,j,k)
                !
                end do
            end do
        end do
        !
        !     periodic
        do kper=1,1-ip
            do k=kparasta,kparaend
                do j=1,jy
                    !
                    usn(j,k)=u(jx,j,k)
                    vsn(j,k)=v(jx,j,k)
                    wsn(j,k)=w(jx,j,k)
                    rhosn(1,j,k)=rhov(1,jx,j,k)
                    !
                    udx(j,k)=u(1,j,k)
                    vdx(j,k)=v(1,j,k)
                    wdx(j,k)=w(1,j,k)
                    rhodx(1,j,k)=rhov(1,1,j,k)
                !
                end do
            end do
        end do
        !
        ! extrapolation for velocity on sides j=0 j=jy+1
        !
        !     not periodic
        do kper=1,jp
            !
            do k=kparasta,kparaend
                do i=1,jx
                    !
                    ust(i,k)=  2.*  u(i,0,k)-  u(i,1,k)
                    vst(i,k)=  2.*  v(i,0,k)-  v(i,1,k)
                    wst(i,k)=  2.*  w(i,0,k)-  w(i,1,k)
                    rhost(1,i,k)=  2.*rhov(1,i,0,k)-rhov(1,i,1,k)
                    !
                    usp(i,k)=  2.*  u(i,jy+1,k)-  u(i,jy,k)
                    vsp(i,k)=  2.*  v(i,jy+1,k)-  v(i,jy,k)
                    wsp(i,k)=  2.*  w(i,jy+1,k)-  w(i,jy,k)
                    rhosp(1,i,k)=  2.*rhov(1,i,jy+1,k)-rhov(1,i,jy,k)
                !
                end do
            end do
        end do
        !
        !     periodic
        do kper=1,1-jp
            !
            do k=kparasta,kparaend
                do i=1,jx
                    !
                    ust(i,k)=u(i,jy,k)
                    vst(i,k)=v(i,jy,k)
                    wst(i,k)=w(i,jy,k)
                    rhost(1,i,k)=rhov(1,i,jy,k)
                    !
                    usp(i,k)=u(i,1,k)
                    vsp(i,k)=v(i,1,k)
                    wsp(i,k)=w(i,1,k)
                    rhosp(1,i,k)=rhov(1,i,1,k)
                !
                end do
            end do
        end do
        !
        ! extrapolation for velocity on sides k=0 k=jz+1
        !
        !     not periodic
        do kper=1,kp
            !
            do j=1,jy
                do i=1,jx

                    if(myid.eq.0)then
                        uin(i,j)=  2.*  u(i,j,0)-  u(i,j,1)
                        vin(i,j)=  2.*  v(i,j,0)-  v(i,j,1)
                        win(i,j)=  2.*  w(i,j,0)-  w(i,j,1)
                        rhoin(1,i,j)= 2.*rhov(1,i,j,0)-rhov(1,i,j,1)

                    else if(myid.eq.nproc-1)then
                        uav(i,j)=  2.*  u(i,j,jz+1)-  u(i,j,jz)
                        vav(i,j)=  2.*  v(i,j,jz+1)-  v(i,j,jz)
                        wav(i,j)=  2.*  w(i,j,jz+1)-  w(i,j,jz)
                        rhoav(1,i,j)=2.*rhov(1,i,j,jz+1)-rhov(1,i,j,jz)
                    endif
                !
                end do
            end do
        end do
        !
        !     periodic
        !     planes for P0 and Pn-1 already known from the previous time step
        do kper=1,1-kp
            !
            do j=1,jy
                do i=1,jx

                    uin(i,j)=0.
                    vin(i,j)=0.
                    win(i,j)=0.
                    rhoin(1,i,j)=0.

                    uav(i,j)=0.
                    vav(i,j)=0.
                    wav(i,j)=0.
                    rhoav(1,i,j)=0.
                    !
                    if(myid.eq.0)then
                        uin(i,j)=u_piano(i,j,jz)
                        vin(i,j)=v_piano(i,j,jz)
                        win(i,j)=w_piano(i,j,jz)
                        rhoin(1,i,j)=rhov_piano(1,i,j,jz)
                    !
                    else if(myid.eq.nproc-1)then
                        uav(i,j)=u_piano(i,j,1)
                        vav(i,j)=v_piano(i,j,1)
                        wav(i,j)=w_piano(i,j,1)
                        rhoav(1,i,j)=rhov_piano(1,i,j,1)
                    endif
                !
                end do
            end do
        end do
        !-----------------------------------------------------------------------
        ! START COMPUTATION EDDY VISCOSITY AND DIFFUSIVITY
        !
        somma=float(jx)*float(jz)

        !     initialization
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    annit(i,j,k)=1./re
                    annitV(i,j,k)=1./re
                    do isc=1,nscal
                        akapt(isc,i,j,k) =1./re/pran(isc)
                        akaptV(isc,i,j,k)=1./re/pran(isc)
                    end do
                end do
            end do
        end do
        if(myid==0 .or. myid==nproc-1) then
            annit_piano  = 1./re
            annitV_piano = 1./re
        end if

        if(myid==0)then
            do j=1,jy
                do i=1,jx
                    do isc=1,nscal
                        akapt_piano(isc,i,j,jz)  = 1./re/pran(isc)
                        akaptV_piano(isc,i,j,jz) = 1./re/pran(isc)
                    end do
                end do
            end do
        end if

        if(myid==nproc-1)then
            do j=1,jy
                do i=1,jx
                    do isc=1,nscal
                        akapt_piano(isc,i,j,1)  = 1./re/pran(isc)
                        akaptV_piano(isc,i,j,1) = 1./re/pran(isc)
                    end do
                end do
            end do
        end if
        !
        !--------------------------------------------------------------------
        ! start computation of the controvariant form of the strian rate
        ! tensor Sij, its module is needed to compute the eddy viscosity
        ! with Smagorinsky model

        if(lagr .eq. 0) then
            ! compute dissipative strain from Smagorinsky model, averaged on
            ! horizontal plane
            allocate(s1_loc(jy))
            allocate(s2_loc(jy))
            allocate(s3_loc(jy))
            allocate(s4_loc(jy))
            allocate(s5_loc(jy))
            allocate(s6_loc(jy))
            allocate(s1rho_loc(jy))
            allocate(s2rho_loc(jy))
            allocate(s3rho_loc(jy))

            allocate(s1_tot(jy))
            allocate(s2_tot(jy))
            allocate(s3_tot(jy))
            allocate(s4_tot(jy))
            allocate(s5_tot(jy))
            allocate(s6_tot(jy))
            allocate(s1rho_tot(jy))
            allocate(s2rho_tot(jy))
            allocate(s3rho_tot(jy))

            allocate(grg_loc(jy))
            allocate(grg_tot(jy))
            s1_loc=0.
            s2_loc=0.
            s3_loc=0.
            s4_loc=0.
            s5_loc=0.
            s6_loc=0.

            s1rho_loc=0.
            s2rho_loc=0.
            s3rho_loc=0.

            grg_loc=0.
        end if


        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx

                    !
                    apcsx(i,j,k)=0.5*(csx(i,j,k)+csx(i-1,j,k))
                    apcsy(i,j,k)=0.5*(csy(i,j,k)+csy(i-1,j,k))
                    apcsz(i,j,k)=0.5*(csz(i,j,k)+csz(i-1,j,k))
                    !
                    apetx(i,j,k)=0.5*(etx(i,j,k)+etx(i,j-1,k))
                    apety(i,j,k)=0.5*(ety(i,j,k)+ety(i,j-1,k))
                    apetz(i,j,k)=0.5*(etz(i,j,k)+etz(i,j-1,k))
                    !
                    apztx(i,j,k)=0.5*(ztx(i,j,k)+ztx(i,j,k-1))
                    apzty(i,j,k)=0.5*(zty(i,j,k)+zty(i,j,k-1))
                    apztz(i,j,k)=0.5*(ztz(i,j,k)+ztz(i,j,k-1))
                    !
                    u2=in_dx1(i,j,k)*u(i+1,j,k) &
                        +(1-in_dx1(i,j,k))*udx(j,k)
                    u1=in_sn1(i,j,k)*u(i-1,j,k) &
                        +(1-in_sn1(i,j,k))*usn(j,k)
                    !
                    u4=in_sp1(i,j,k)*u(i,j+1,k) &
                        +(1-in_sp1(i,j,k))*usp(i,k)
                    u3=in_st1(i,j,k)*u(i,j-1,k) &
                        +(1-in_st1(i,j,k))*ust(i,k)
                    !
                    u6=in_av1(i,j,k)*u(i,j,k+1) &
                        +(1-in_av1(i,j,k))*uav(i,j)
                    u5=in_in1(i,j,k)*u(i,j,k-1) &
                        +(1-in_in1(i,j,k))*uin(i,j)
                    !
                    v2=in_dx1(i,j,k)*v(i+1,j,k) &
                        +(1-in_dx1(i,j,k))*vdx(j,k)
                    v1=in_sn1(i,j,k)*v(i-1,j,k) &
                        +(1-in_sn1(i,j,k))*vsn(j,k)
                    !
                    v4=in_sp1(i,j,k)*v(i,j+1,k) &
                        +(1-in_sp1(i,j,k))*vsp(i,k)
                    v3=in_st1(i,j,k)*v(i,j-1,k) &
                        +(1-in_st1(i,j,k))*vst(i,k)
                    !
                    v6=in_av1(i,j,k)*v(i,j,k+1) &
                        +(1-in_av1(i,j,k))*vav(i,j)
                    v5=in_in1(i,j,k)*v(i,j,k-1) &
                        +(1-in_in1(i,j,k))*vin(i,j)
                    !
                    w2=in_dx1(i,j,k)*w(i+1,j,k) &
                        +(1-in_dx1(i,j,k))*wdx(j,k)
                    w1=in_sn1(i,j,k)*w(i-1,j,k) &
                        +(1-in_sn1(i,j,k))*wsn(j,k)
                    !
                    w4=in_sp1(i,j,k)*w(i,j+1,k) &
                        +(1-in_sp1(i,j,k))*wsp(i,k)
                    w3=in_st1(i,j,k)*w(i,j-1,k) &
                        +(1-in_st1(i,j,k))*wst(i,k)
                    !
                    w6=in_av1(i,j,k)*w(i,j,k+1) &
                        +(1-in_av1(i,j,k))*wav(i,j)
                    w5=in_in1(i,j,k)*w(i,j,k-1) &
                        +(1-in_in1(i,j,k))*win(i,j)
                    !
                    rho2=in_dx1(i,j,k)*rhov(1,i+1,j,k) &
                        +(1-in_dx1(i,j,k))*rhodx(1,j,k)
                    rho1=in_sn1(i,j,k)*rhov(1,i-1,j,k) &
                        +(1-in_sn1(i,j,k))*rhosn(1,j,k)
                    !
                    rho4=in_sp1(i,j,k)*rhov(1,i,j+1,k) &
                        +(1-in_sp1(i,j,k))*rhosp(1,i,k)
                    rho3=in_st1(i,j,k)*rhov(1,i,j-1,k) &
                        +(1-in_st1(i,j,k))*rhost(1,i,k)
                    !
                    rho6=in_av1(i,j,k)*rhov(1,i,j,k+1) &
                        +(1-in_av1(i,j,k))*rhoav(1,i,j)
                    rho5=in_in1(i,j,k)*rhov(1,i,j,k-1) &
                        +(1-in_in1(i,j,k))*rhoin(1,i,j)


                    dudx=(.5*(u2-u1)*apcsx(i,j,k)+ &
                        .5*(u4-u3)*apetx(i,j,k)+ &
                        .5*(u6-u5)*apztx(i,j,k))/giac(i,j,k)
                    !
                    dudy=(.5*(u2-u1)*apcsy(i,j,k)+ &
                        .5*(u4-u3)*apety(i,j,k)+ &
                        .5*(u6-u5)*apzty(i,j,k))/giac(i,j,k)
                    !
                    dudz=(.5*(u2-u1)*apcsz(i,j,k)+ &
                        .5*(u4-u3)*apetz(i,j,k)+ &
                        .5*(u6-u5)*apztz(i,j,k))/giac(i,j,k)
                    !
                    dvdx=(.5*(v2-v1)*apcsx(i,j,k)+ &
                        .5*(v4-v3)*apetx(i,j,k)+ &
                        .5*(v6-v5)*apztx(i,j,k))/giac(i,j,k)
                    !
                    dvdy=(.5*(v2-v1)*apcsy(i,j,k)+ &
                        .5*(v4-v3)*apety(i,j,k)+ &
                        .5*(v6-v5)*apzty(i,j,k))/giac(i,j,k)
                    !
                    dvdz=(.5*(v2-v1)*apcsz(i,j,k)+ &
                        .5*(v4-v3)*apetz(i,j,k)+ &
                        .5*(v6-v5)*apztz(i,j,k))/giac(i,j,k)
                    !
                    dwdx=(.5*(w2-w1)*apcsx(i,j,k)+ &
                        .5*(w4-w3)*apetx(i,j,k)+ &
                        .5*(w6-w5)*apztx(i,j,k))/giac(i,j,k)
                    !
                    dwdy=(.5*(w2-w1)*apcsy(i,j,k)+ &
                        .5*(w4-w3)*apety(i,j,k)+ &
                        .5*(w6-w5)*apzty(i,j,k))/giac(i,j,k)
                    !
                    dwdz=(.5*(w2-w1)*apcsz(i,j,k)+ &
                        .5*(w4-w3)*apetz(i,j,k)+ &
                        .5*(w6-w5)*apztz(i,j,k))/giac(i,j,k)


                    drhodx=(.5*(rho2-rho1)*apcsx(i,j,k)+ &
                        .5*(rho4-rho3)*apetx(i,j,k)+ &
                        .5*(rho6-rho5)*apztx(i,j,k))/giac(i,j,k)
                    !
                    drhody=(.5*(rho2-rho1)*apcsy(i,j,k)+ &
                        .5*(rho4-rho3)*apety(i,j,k)+ &
                        .5*(rho6-rho5)*apzty(i,j,k))/giac(i,j,k)
                    !
                    drhodz=(.5*(rho2-rho1)*apcsz(i,j,k)+ &
                        .5*(rho4-rho3)*apetz(i,j,k)+ &
                        .5*(rho6-rho5)*apztz(i,j,k))/giac(i,j,k)

                    !-----------------------------------------------------------------------
                    ! for wall modeling, dudy and dwdy computation considering the coarse grid

                    !hicco mettere anche in funzione di coeff wall
                    do lll=1,att_wm_sgs
                        if(j==1)then
                            do iii=1,att_mod_par(i,1,k)
                                dudy =(u_t(i,1,k)/(0.41*punto_wfp3(3,1,i,k))) &
                                    *u(i,j,k)/utangente(i,1,k)
                                dwdy =(u_t(i,1,k)/(0.41*punto_wfp3(3,1,i,k))) &
                                    *w(i,j,k)/utangente(i,1,k)
                            enddo
                        elseif(j==jy) then
                            do iii=1,att_mod_par(i,2,k)
                                dudy =(u_t(i,2,k)/(0.41*punto_wfp4(3,1,i,k))) &
                                    *u(i,j,k)/utangente(i,2,k)
                                dwdy =(u_t(i,2,k)/(0.41*punto_wfp4(3,1,i,k))) &
                                    *w(i,j,k)/utangente(i,2,k)
                            enddo
                        endif
                    enddo
                    !-----------------------------------------------------------------------
                    !
                    s1=dudx
                    s2=.5*(dudy+dvdx)
                    s3=.5*(dudz+dwdx)
                    s4=dvdy
                    s5=.5*(dvdz+dwdy)
                    s6=dwdz
                    !
                    smod(i,j,k)=sqrt(2.*s1*s1+ &
                        2.*s4*s4+ &
                        2.*s6*s6+ &
                        4.*s2*s2+ &
                        4.*s3*s3+ &
                        4.*s5*s5)

                    smodV(i,j,k)=sqrt(2.*s4*s4+ &
                        4.*s2*s2+ &
                        4.*s5*s5)

                    smodH(i,j,k)=sqrt(2.*s1*s1+ &
                        2.*s6*s6+ &
                        4.*s3*s3)
                    !
                    if(ktime.eq.i_print*(ktime/i_print)) then

                        do pp=1,1-lagr
                            !
                            grg=cost*cost*giac(i,j,k)**(2./3.)*smod(i,j,k)

                            s1=grg*s1
                            s2=grg*s2
                            s3=grg*s3
                            s4=grg*s4
                            s5=grg*s5
                            s6=grg*s6


                            s1_loc(j)=s1_loc(j)+s1
                            s2_loc(j)=s2_loc(j)+s2
                            s3_loc(j)=s3_loc(j)+s3
                            s4_loc(j)=s4_loc(j)+s4
                            s5_loc(j)=s5_loc(j)+s5
                            s6_loc(j)=s6_loc(j)+s6

                            s1rho_loc(j)=s1rho_loc(j)+0.5*grg*drhodx
                            s2rho_loc(j)=s2rho_loc(j)+0.5*grg*drhody
                            s3rho_loc(j)=s3rho_loc(j)+0.5*grg*drhodz

                            grg_loc(j)=grg_loc(j)+grg
                        end do
                    end if

                enddo
            enddo
        enddo

        if(ktime.eq.i_print*(ktime/i_print)) then
            do pp=1,1-lagr
                call MPI_REDUCE(s1_loc(1),s1_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(s2_loc(1),s2_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(s3_loc(1),s3_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(s4_loc(1),s4_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(s5_loc(1),s5_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(s6_loc(1),s6_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)

                call MPI_REDUCE(s1rho_loc(1),s1rho_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(s2rho_loc(1),s2rho_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(s3rho_loc(1),s3rho_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(grg_loc(1),grg_tot(1),jy,MPI_REAL_SD, &
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                if( myid.eq.0)then
                    s1_tot=s1_tot/(float(jx)*float(jz))
                    s2_tot=s2_tot/(float(jx)*float(jz))
                    s3_tot=s3_tot/(float(jx)*float(jz))
                    s4_tot=s4_tot/(float(jx)*float(jz))
                    s5_tot=s5_tot/(float(jx)*float(jz))
                    s6_tot=s6_tot/(float(jx)*float(jz))

                    s1rho_tot=s1rho_tot/(float(jx)*float(jz))
                    s2rho_tot=s2rho_tot/(float(jx)*float(jz))
                    s3rho_tot=s3rho_tot/(float(jx)*float(jz))

                    grg_tot=grg_tot/(float(jx)*float(jz))
                    do j=1,jy
                        if(j==1)then
                            grg_tot(1)=grg_tot(2)+((y(i,1,k)-y(i,2,k))/ &
                                (y(i,3,k)-y(i,2,k))) &
                                *(grg_tot(3)-grg_tot(2))
                        end if
                    end do
                end if
            end do
        end if !iprint
6666    format(8e14.4)
        !
        ! extrapolation on sides i=0 i=jx+1
        ! periodic/ not periodic
        do k=kparasta,kparaend
            do j=1,jy
                !
                smod(0,j,k) =(1-ip)* smod(jx,j,k) + &
                    ip*(2.* smod(1 ,j,k)- smod(2,j,k))
                smodV(0,j,k)=(1-ip)*smodV(jx,j,k) + &
                    ip*(2.*smodV(1 ,j,k)-smodV(2,j,k))
                smodH(0,j,k)=(1-ip)*smodH(jx,j,k) + &
                    ip*(2.*smodH(1 ,j,k)-smodH(2,j,k))


                smod(jx+1,j,k) =(1-ip)* smod(1 ,j,k) + &
                    ip*(2.* smod(jx,j,k)- smod(jx-1,j,k))
                smodV(jx+1,j,k)=(1-ip)*smodV(1 ,j,k) + &
                    ip*(2.*smodV(jx,j,k)-smodV(jx-1,j,k))
                smodH(jx+1,j,k)=(1-ip)*smodH(1 ,j,k) + &
                    ip*(2.*smodH(jx,j,k)-smodH(jx-1,j,k))

            end do
        end do
        !
        ! extrapolation on sides j=0 j=jy+1
        ! periodic/ not periodic
        do k=kparasta,kparaend
            do i=1,jx
                !
                smod(i,0,k) =(1-jp)* smod(i,jy,k) + &
                    jp*(2.*smod(i,1,k)-smod(i,2,k))
                smodV(i,0,k)=(1-jp)*smodV(i,jy,k) + &
                    jp*(2.*smodV(i,1,k)-smodV(i,2,k))
                smodH(i,0,k)=(1-jp)*smodH(i,jy,k) + &
                    jp*(2.*smodH(i,1,k)-smodH(i,2,k))

                !
                smod(i,jy+1,k) =(1-jp)*smod(i,1,k) + &
                    jp*(2.*smod(i,jy,k)-smod(i,jy-1,k))
                smodV(i,jy+1,k)=(1-jp)*smodV(i,1,k) + &
                    jp*(2.*smodV(i,jy,k)-smodV(i,jy-1,k))
                smodH(i,jy+1,k)=(1-jp)*smodH(i,1,k) + &
                    jp*(2.*smodH(i,jy,k)-smodH(i,jy-1,k))

            end do
        end do

        ! now periodicity/not periodicity on front and back
        ! including ghost cells
        ! extrapolation for velocity on sides k=0 and k=jz+1
        !
        ! I make a vector buff which contains all the planes jx*jy
        ! for exchange between P0 and Pn-1

        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        if (myid.eq.nproc-1) then

            call buffer1g(smod ,1,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(smodV,2,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(smodH,3,jz,myid,nproc,kparasta,kparaend)

        else if (myid.eq.0) then

            call buffer1g(smod ,1,1,myid,nproc,kparasta,kparaend)
            call buffer1g(smodV,2,1,myid,nproc,kparasta,kparaend)
            call buffer1g(smodH,3,1,myid,nproc,kparasta,kparaend)

        endif
        !
        ! exchange so P0 knows k=jz and Pn-1 knows j=1
        !

        if (myid.eq.nproc-1) then

            call MPI_SENDRECV(sbuff1(1),(3)*(jx+2)*(jy+2),MPI_REAL_SD, &
                0,9001, &
                rbuff1(1),(3)*(jx+2)*(jy+2),MPI_REAL_SD, &
                0,8001, &
                MPI_COMM_WORLD,status,ierr)

        !if(ierr == MPI_SUCCESS) print *, "va bene", myid
        else if (myid.eq.0) then

            call MPI_SENDRECV(sbuff1(1),(3)*(jx+2)*(jy+2),MPI_REAL_SD, &
                nproc-1, 8001, &
                rbuff1(1),(3)*(jx+2)*(jy+2),MPI_REAL_SD, &
                nproc-1, 9001, &
                MPI_COMM_WORLD,status,ierr)

        !if(ierr == MPI_SUCCESS) print *, "va bene", myid
        endif
        ! put the values on the correct place
        if (myid.eq.0) then
            call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)

        else if (myid.eq.nproc-1) then
            call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)


        endif
        !
        ! now P0 knows k=jz
        ! and Pn-1 knows k=1
        !
        if(myid.eq.0)then

            do j=1,jy
                do i=1,jx

                    smod(i,j,0)=(1-kp)*piano1(i,j) + &
                        kp*(2.*smod(i,j,1)-smod(i,j,2))
                    smodV(i,j,0)=(1-kp)*piano2(i,j) + &
                        kp*(2.*smodV(i,j,1)-smodV(i,j,2))
                    smodH(i,j,0)=(1-kp)*piano3(i,j) + &
                        kp*(2.*smodH(i,j,1)-smodH(i,j,2))

                end do
            end do

        endif
        !
        if(myid.eq.nproc-1)then

            do j=1,jy
                do i=1,jx

                    smod(i,j,jz+1)=(1-kp)*piano1(i,j) + &
                        kp*(2.*smod(i,j,jz)-smod(i,j,jz-1))
                    smodV(i,j,jz+1)=(1-kp)*piano2(i,j) + &
                        kp*(2.*smodV(i,j,jz)-smodV(i,j,jz-1))
                    smodH(i,j,jz+1)=(1-kp)*piano3(i,j) + &
                        kp*(2.*smodH(i,j,jz)-smodH(i,j,jz-1))

                end do
            end do

        endif


        if(isotropo == 1)then
            !     isotropic model, one eddy viscosity

            costante = cost * cost
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        !
                        annit(i,j,k)=annit(i,j,k) &
                            +giac(i,j,k)**(2./3.)*costante*smod(i,j,k)
                        annitV(i,j,k)=annit(i,j,k)
                        !
                        !          prandtl turbolent = 0.5 ==>  Kt=2*Vt
                        do isc=1,nscal
                            akapt(isc,i,j,k)  =akapt(isc,i,j,k) &
                                +prsc(isc)* giac(i,j,k)**(2./3.)*costante*smod(i,j,k)
                            akaptV(isc,i,j,k) =akapt(isc,i,j,k)

                        end do
                    !
                    enddo
                enddo
            enddo
        end if


        if(isotropo == 0)then
            !     anisotropic model, two eddy viscosity
            costanteH = costH * costH
            costanteV = costV * costV

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx

                        lhorx=.25*(x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
                            +x(i  ,j-1,k-1)+x(i  ,j-1,k  )) &
                            -.25*(x(i-1,j  ,k  )+x(i-1,j  ,k-1) &
                            +x(i-1,j-1,k-1)+x(i-1,j-1,k  ))

                        lhorz=.25*(z(i  ,j  ,k  )+z(i  ,j-1,k  ) &
                            +z(i-1,j-1,k  )+z(i-1,j  ,k  )) &
                            -.25*(z(i  ,j  ,k-1)+z(i  ,j-1,k-1) &
                            +z(i-1,j-1,k-1)+z(i-1,j  ,k-1))

                        !           lhor = lhorx**2. + lhorz**2. !length square
                        lhor = lhorx*lhorx + lhorz*lhorz !length square

                        lver =.25*(y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
                            +y(i-1,j  ,k-1)+y(i-1,j  ,k  )) &
                            -.25*(y(i  ,j-1,k  )+y(i  ,j-1,k-1) &
                            +y(i-1,j-1,k-1)+y(i-1,j-1,k  ))

                        lver = lver*lver !**2.


                        annit(i,j,k) =annit(i,j,k) &
                            + lhor*costanteH*smodH(i,j,k)

                        annitV(i,j,k)=annitV(i,j,k) &
                            + lver*costanteV*smodV(i,j,k)

                        !          prandtl turbolent = 0.5 ==>  Kt=2*Vt
                        do isc=1,nscal


                            akapt(isc,i,j,k) =akapt(isc,i,j,k) &
                                +prsc(isc)* lhor*costanteH*smodH(i,j,k)


                            !          prandtl turbolent = 0.5 ==>  Kt=2*Vt
                            akaptV(isc,i,j,k)=akaptV(isc,i,j,k) &
                                +prsc(isc)* lver*costanteV*smodV(i,j,k)

                        end do
                    !
                    enddo
                enddo
            enddo
        end if
        !
        ! add eddy viscosity and diffusivity on periodic sides
        !
        ! sides left and right
        do k=kparasta,kparaend
            do j=1,jy
                annit(0,j,k)   =(1-ip)*annit(jx,j,k) &
                    +ip*annit(1,j,k)
                annitV(0,j,k)   =(1-ip)*annitV(jx,j,k) &
                    +ip*annitV(1,j,k)


                annit(jx+1,j,k)=(1-ip)*annit(1,j,k) &
                    +ip*annit(jx,j,k)
                annitV(jx+1,j,k)=(1-ip)*annitV(1,j,k) &
                    +ip*annitV(jx,j,k)

                do isc=1,nscal
                    akapt(isc,0,j,k)   =(1-ip)*akapt(isc,jx,j,k) &
                        +ip*akapt(isc,1,j,k)
                    akaptV(isc,0,j,k)   =(1-ip)*akaptV(isc,jx,j,k) &
                        +ip*akaptV(isc,1,j,k)
                    !
                    akapt(isc,jx+1,j,k)=(1-ip)*akapt(isc,1,j,k) &
                        +ip*akapt(isc,jx,j,k)
                    akaptV(isc,jx+1,j,k)=(1-ip)*akaptV(isc,1,j,k) &
                        +ip*akaptV(isc,jx,j,k)
                end do
            end do
        end do
        !
        ! sides back and front
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        if (myid.eq.nproc-1) then

            call buffer1old_par(sbuff,annit ,1,jz)
            call buffer1old_par(sbuff,annitV,2,jz)
            cont = 2
            do isc=1,nscal
                cont = cont + 1
                call buffer1old_par_nscal(sbuff,akapt ,cont,jz,isc)
                cont = cont + 1
                call buffer1old_par_nscal(sbuff,akaptV,cont,jz,isc)
            end do

        else if (myid.eq.0) then

            call buffer1old_par(sbuff,annit ,1,1)
            call buffer1old_par(sbuff,annitV,2,1)
            cont = 2
            do isc= 1,nscal
                cont = cont+1
                call buffer1old_par_nscal(sbuff,akapt ,cont,1,isc)
                cont = cont+1
                call buffer1old_par_nscal(sbuff,akaptV,cont,1,isc)
            end do

        endif

        cont = 2+2*nscal
        !
        ! exchange the planes so P0 knows k=jz
        ! and Pn-1 knows k=1
        !
        if (myid.eq.nproc-1) then

            call MPI_SENDRECV(sbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,0,3001, &
                rbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,0,2001, &
                MPI_COMM_WORLD,status,ierr)

        !if(ierr == MPI_SUCCESS) print *, "va bene2", myid
        else if (myid.eq.0) then

            call MPI_SENDRECV( &
                sbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001, &
                rbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,3001, &
                MPI_COMM_WORLD,status,ierr)

        !if(ierr == MPI_SUCCESS) print *, "va bene2", myid
        endif

        if (myid.eq.0) then

            call buffer2old_par(rbuff,annit_piano ,1,jz)
            call buffer2old_par(rbuff,annitV_piano,2,jz)
            cont = 2
            do isc=1,nscal
                cont = cont + 1
                call buffer2old_par_nscal(rbuff,akapt_piano ,cont,jz,isc)
                cont = cont + 1
                call buffer2old_par_nscal(rbuff,akaptV_piano,cont,jz,isc)
            end do

        else if (myid.eq.nproc-1) then

            call buffer2old_par(rbuff,annit_piano ,1,1)
            call buffer2old_par(rbuff,annitV_piano,2,1)
            cont = 2
            do isc=1,nscal
                cont = cont + 1
                call buffer2old_par_nscal(rbuff,akapt_piano ,cont,1,isc)
                cont = cont + 1
                call buffer2old_par_nscal(rbuff,akaptV_piano,cont,1,isc)
            end do

        endif
        !
        ! now P0 knows k=jz
        ! and Pn-1 knows k=1
        !
        if(myid.eq.0)then

            do j=1,jy
                do i=1,jx
                    annit(i,j,0) =(1-kp)*annit_piano(i,j,jz)+kp*annit(i,j,1)
                    annitV(i,j,0)=(1-kp)*annitV_piano(i,j,jz)+kp*annitV(i,j,1)
                    do isc=1,nscal
                        akapt(isc,i,j,0) =(1-kp)*akapt_piano(isc,i,j,jz) &
                            +kp*akapt(isc,i,j,1)
                        akaptV(isc,i,j,0)=(1-kp)*akaptV_piano(isc,i,j,jz) &
                            +kp*akaptV(isc,i,j,1)
                    end do
                end do
            end do

        else if(myid.eq.nproc-1)then

            do j=1,jy
                do i=1,jx
                    annit(i,j,jz+1)=(1-kp)*annit_piano(i,j,1)+kp*annit(i,j,jz)
                    annitV(i,j,jz+1)=(1-kp)*annitV_piano(i,j,1) &
                        +kp*annitV(i,j,jz)
                    do isc=1,nscal
                        akapt(isc,i,j,jz+1)=(1-kp)*akapt_piano(isc,i,j,1) &
                            +kp*akapt(isc,i,j,jz)


                        akaptV(isc,i,j,jz+1)=(1-kp)*akaptV_piano(isc,i,j,1) &
                            +kp*akaptV(isc,i,j,jz)
                    end do
                end do
            end do

        endif
        !
        ! exchange closer planes between procs
        !
        !     annit
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(annit(0,0,kparasta),(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,tagls, &
                MPI_COMM_WORLD,ierr)
        !if(ierr == MPI_SUCCESS) print *, "va bene3", myid
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(annit(0,0,kparaend+1),(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrr, &
                MPI_COMM_WORLD,status,ierr)
        !if(ierr == MPI_SUCCESS) print *, "va bene3", myid
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(annit(0,0,kparaend),(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrs, &
                MPI_COMM_WORLD,ierr)
        !if(ierr == MPI_SUCCESS) print *, "va bene3", myid
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(annit(0,0,kparasta-1),(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,taglr, &
                MPI_COMM_WORLD,status,ierr)
        !if(ierr == MPI_SUCCESS) print *, "va bene3", myid
        endif

        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        !      call MPI_WAIT(req4,istatus,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        !      call MPI_WAIT(req3,istatus,ierr)
        endif

        !     annit and annitV
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(annitV(0,0,kparasta),(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,tagls, &
                MPI_COMM_WORLD,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(annitV(0,0,kparaend+1),(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrr, &
                MPI_COMM_WORLD,status,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(annitV(0,0,kparaend),(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrs, &
                MPI_COMM_WORLD,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(annitV(0,0,kparasta-1),(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,taglr, &
                MPI_COMM_WORLD,status,ierr)
        endif

        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        !      call MPI_WAIT(req4,istatus,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        !      call MPI_WAIT(req3,istatus,ierr)
        endif
        !
        !-----------------------------------------------------------------------
        !
        !     annitV and akaptV vertical component
        !
        allocate (P_akapt(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate (P_akaptV(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))

        do isc=1,nscal

            P_akapt = 0.
            P_akaptV= 0.

            do k=kparasta-deepl,kparasta+deepr
                do j=0,n2+1
                    do i=0,n1+1
                        P_akapt(i,j,k) = akapt(isc,i,j,k)
                        P_akaptV(i,j,k)=akaptV(isc,i,j,k)
                    end do
                end do
            end do

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SSEND(P_akapt(0,0,kparasta),(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(P_akapt(0,0,kparaend+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SSEND(P_akapt(0,0,kparaend),(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(P_akapt(0,0,kparasta-1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            !      call MPI_WAIT(req4,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            !      call MPI_WAIT(req3,istatus,ierr)
            endif


            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SSEND(P_akaptV(0,0,kparasta),(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(P_akaptV(0,0,kparaend+1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SSEND(P_akaptV(0,0,kparaend),(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(P_akaptV(0,0,kparasta-1),(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            !      call MPI_WAIT(req4,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            !      call MPI_WAIT(req3,istatus,ierr)
            endif

            do k=kparasta-deepl,kparasta+deepr
                do j=0,n2+1
                    do i=0,n1+1
                        akapt(isc,i,j,k) =  P_akapt(i,j,k)
                        akaptV(isc,i,j,k) = P_akaptV(i,j,k)
                    end do
                end do
            end do

        end do !cycle on nscal

        deallocate(P_akapt)
        deallocate(P_akaptV)

    end subroutine turbo_statico
end module ts_module

