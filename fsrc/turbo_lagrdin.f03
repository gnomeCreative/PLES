subroutine turbo_lagrdin(ktime,in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kpstamg,kpendmg)

    use filter_module
    use turbo_module
    use mysending
    use mysettings, only: cost,costH,costV,inmod,inmodrho,nsgs,isotropo,pran,prsc,re_analogy

    use myarrays_velo3
    use myarrays_metri3
    !
    use subgrid
    use period
    use velpar
    use scala3
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    ! parameters
    integer,intent(in) :: ktime
    integer,intent(in) :: kpstamg(0:4),kpendmg(0:4)
    integer,intent(in) :: in_dx1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_sn1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_sp1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_st1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_av1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_in1(n1,n2,kpstamg(1):kpendmg(1))

    !-----------------------------------------------------------------------
    !     variable declaration
    integer i,j,k
    integer kper,lll,m

    integer ierr
    integer status(MPI_STATUS_SIZE)
    integer kparastal,kparaendl

    integer debugg,ttt
    !
    real sbuff((n1+2)*(n2+2)*40)
    real rbuff((n1+2)*(n2+2)*40)
    real sbuffd((n1+2)*(n2+2)*40)
    real rbuffd((n1+2)*(n2+2)*40)
    !
    real cb(n2)
    real cbrho(n2)
    real c1,c2
    real somma,delt,psp,alfa,alfa3
    real grg,grgrho,costrho,rzer
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
    !
    real c(n2),crho(n2)
    real den(n1,n2,kparasta:kparaend)
    real num(n1,n2,kparasta:kparaend)
    real denrho(n1,n2,kparasta:kparaend)
    real numrho(n1,n2,kparasta:kparaend)
    real nummed_loc(n2),denmed_loc(n2)
    real numrhomed_loc(n2),denrhomed_loc(n2)
    real nummed(n2),denmed(n2)
    real numrhomed(n2),denrhomed(n2)
    !
    integer medio_in_k ! dinamico only
    real nummed_z(n1,n2),denmed_z(n1,n2) ! dinamico only
    real nummed_loc_z(n1,n2),denmed_loc_z(n1,n2) ! dinamico only
    real numrhomed_z(n1,n2),denrhomed_z(n1,n2) ! dinamico only
    real numrhomed_loc_z(n1,n2),denrhomed_loc_z(n1,n2) ! dinamico only
    real linea_z ! dinamico only
    !
    real sub11_loc(n2),sub12_loc(n2),sub13_loc(n2)
    real sub22_loc(n2),sub23_loc(n2),sub33_loc(n2)
    real sub_loc(n2)
    real subrho_loc11(n2),subrho_loc22(n2), &
        subrho_loc33(n2)
    real sus_loc11(n2),sus_loc12(n2),sus_loc13(n2)
    real sus_loc22(n2),sus_loc23(n2),sus_loc33(n2)
    real sus_loc(n2)
    real susrho_loc11(n2),susrho_loc22(n2), &
        susrho_loc33(n2)
    !
    real a1nu(n2),a1dn(n2),a1nurho(n2),a1dnrho(n2)
    !
    real tp,ti

    real lhor,lhorx,lhorz,lver
    !
    integer iso
    real exp2over3

    integer isc,cont
    real, allocatable :: P_akapt(:,:,:)
    real,allocatable :: rho(:,:,:)
    integer i_isc ! only lagrangian
    !-----------------------------------------------------------------------
    ! parte lagrangiana
    real,allocatable :: x1(:,:,:),y1(:,:,:),z1(:,:,:)
    real,allocatable :: eps(:,:,:),eps_rho(:,:,:),vvf(:,:,:,:)

    real base,tem,supres
    real base_rho,tem_rho
    !-----------------------------------------------------------------------

    if (nsgs==2) then
        medio_in_k = 0   ! if=0 average in i and k; if =1 average in k
    else if (nsgs==3) then

        allocate(x1(n1,n2,kparasta:kparaend))
        allocate(y1(n1,n2,kparasta:kparaend))
        allocate(z1(n1,n2,kparasta:kparaend))
        allocate(eps(n1,n2,kparasta:kparaend))
        allocate(eps_rho(n1,n2,kparasta:kparaend))
        allocate(vvf(n1,n2,kparasta:kparaend,4))

    end if


    allocate(rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
    rho = 0

    do k=kparasta-1,kparaend+1
        do j=0,jy+1
            do i=0,jx+1
                rho(i,j,k) = rhov(1,i,j,k)
            end do
        end do
    end do

    debugg = 200
    if (ktime.eq.1) then
        ttt = 0
    end if
    if (ktime.eq.2) then
        ttt = 10000
    end if
    !-----------------------------------------------------------------------
    ! velocity extrapolation on face 1 and 2, i=0 and i=jx+1
    !
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
                rhosn(1,j,k)=  2.*     rho(0,j,k)-     rho(1,j,k)
                !
                udx(j,k)=  2.*  u(jx+1,j,k)-  u(jx,j,k)
                vdx(j,k)=  2.*  v(jx+1,j,k)-  v(jx,j,k)
                wdx(j,k)=  2.*  w(jx+1,j,k)-  w(jx,j,k)
                rhodx(1,j,k)=  2.*     rho(jx+1,j,k)-     rho(jx,j,k)
            !
            end do
        end do
    end do
    !
    !     periodic
    do kper=1,1-ip
        !
        do k=kparasta,kparaend
            do j=1,jy
                !
                usn(j,k)=u(jx,j,k)
                vsn(j,k)=v(jx,j,k)
                wsn(j,k)=w(jx,j,k)
                rhosn(1,j,k)=rho(jx,j,k)
                !
                udx(j,k)=u(1,j,k)
                vdx(j,k)=v(1,j,k)
                wdx(j,k)=w(1,j,k)
                rhodx(1,j,k)=rho(1,j,k)
            !
            end do
        end do
    end do
    !
    ! velocity extrapolation on face 3 and 4, j=0 and j=jy+1
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
                rhost(1,i,k)=2.*rho(i,0,k)-rho(i,1,k)
                !
                usp(i,k)=  2.*  u(i,jy+1,k)-  u(i,jy,k)
                vsp(i,k)=  2.*  v(i,jy+1,k)-  v(i,jy,k)
                wsp(i,k)=  2.*  w(i,jy+1,k)-  w(i,jy,k)
                rhosp(1,i,k)=2.*rho(i,jy+1,k)-rho(i,jy,k)
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
                rhost(1,i,k)=rho(i,jy,k)
                !
                usp(i,k)=u(i,1,k)
                vsp(i,k)=v(i,1,k)
                wsp(i,k)=w(i,1,k)
                rhosp(1,i,k)=rho(i,1,k)
            !
            end do
        end do
    end do
    !
    ! velocity extrapolation on face 5 and 6, k=0 and k=jz+1
    !
    !     not periodic
    do kper=1,kp
        !
        do j=1,jy
            do i=1,jx
                !
                if(myid.eq.0)then
                    uin(i,j)=  2.*  u(i,j,0)-  u(i,j,1)
                    vin(i,j)=  2.*  v(i,j,0)-  v(i,j,1)
                    win(i,j)=  2.*  w(i,j,0)-  w(i,j,1)
                    rhoin(1,i,j)=2.*rho(i,j,0)-rho(i,j,1)
                !
                else if(myid.eq.nproc-1)then
                    uav(i,j)=  2.*  u(i,j,jz+1)-  u(i,j,jz)
                    vav(i,j)=  2.*  v(i,j,jz+1)-  v(i,j,jz)
                    wav(i,j)=  2.*  w(i,j,jz+1)-  w(i,j,jz)
                    rhoav(1,i,j)=2.*rho(i,j,jz+1)-rho(i,j,jz)
                endif
            !
            end do
        end do
    end do
    !
    !     periodic
    !     transfer data between P0 and Pn-1
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
    !     START COMPUTATION EDDY VISCOSITY AND DIFFUSIVITY
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

    if(myid==0 .or. myid==nproc-1)then
        annit_piano  = 1./re
        annitV_piano = 1./re
        do isc=1,nscal
            akapt_piano  = 1./re/pran(isc)
            akaptV_piano = 1./re/pran(isc)
        end do
    end if
    !
    if (inmod) then
        do j=1,jy
            cb(j)=1.0       ! 1 mixed , 0 smagorinski
        end do
    else
        do j=1,jy
            cb(j)=0.0       ! 1 mixed , 0 smagorinski
        end do
    end if

    if (inmodrho) then
        do j=1,jy
            cbrho(j)=1.0 ! 1 mixed , 0 smagorinski
        end do
    else
        do j=1,jy
            cbrho(j)=0.0 ! 1 mixed , 0 smagorinski
        end do
    end if

    cb(1)=0.
    cb(jy)=0.
    cbrho(1)=0.
    cbrho(jy)=0.
    !
    psp=1./3.
    alfa=sqrt(6.)*sqrt(6.)*sqrt(6.)
    alfa3=alfa**(2./3.)
    !
    !-----------------------------------------------------------------------
    ! start computation of controvariant form of strain rate tensor Sij
    ! its module is needed to find the eddy viscosity with Smagorinsky model
    !
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
                rho2=in_dx1(i,j,k)*rho(i+1,j,k) &
                    +(1-in_dx1(i,j,k))*rhodx(1,j,k)
                rho1=in_sn1(i,j,k)*rho(i-1,j,k) &
                    +(1-in_sn1(i,j,k))*rhosn(1,j,k)
                !
                rho4=in_sp1(i,j,k)*rho(i,j+1,k) &
                    +(1-in_sp1(i,j,k))*rhosp(1,i,k)
                rho3=in_st1(i,j,k)*rho(i,j-1,k) &
                    +(1-in_st1(i,j,k))*rhost(1,i,k)
                !
                rho6=in_av1(i,j,k)*rho(i,j,k+1) &
                    +(1-in_av1(i,j,k))*rhoav(1,i,j)
                rho5=in_in1(i,j,k)*rho(i,j,k-1) &
                    +(1-in_in1(i,j,k))*rhoin(1,i,j)
                !

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
                !
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
                !
                s1=dudx
                s2=.5*(dudy+dvdx)
                s3=.5*(dudz+dwdx)
                s4=dvdy
                s5=.5*(dvdz+dwdy)
                s6=dwdz
                !
                ! controvariant strain rate tensor elements
                !
                s11(i,j,k)=apcsx(i,j,k)*s1+ &
                    apcsy(i,j,k)*s2+ &
                    apcsz(i,j,k)*s3
                !
                s12(i,j,k)=apetx(i,j,k)*s1+ &
                    apety(i,j,k)*s2+ &
                    apetz(i,j,k)*s3
                !
                s13(i,j,k)=apztx(i,j,k)*s1+ &
                    apzty(i,j,k)*s2+ &
                    apztz(i,j,k)*s3
                !
                s21(i,j,k)=apcsx(i,j,k)*s2+ &
                    apcsy(i,j,k)*s4+ &
                    apcsz(i,j,k)*s5
                !
                s22(i,j,k)=apetx(i,j,k)*s2+ &
                    apety(i,j,k)*s4+ &
                    apetz(i,j,k)*s5
                !
                s23(i,j,k)=apztx(i,j,k)*s2+ &
                    apzty(i,j,k)*s4+ &
                    apztz(i,j,k)*s5
                !
                s31(i,j,k)=apcsx(i,j,k)*s3+ &
                    apcsy(i,j,k)*s5+ &
                    apcsz(i,j,k)*s6
                !
                s32(i,j,k)=apetx(i,j,k)*s3+ &
                    apety(i,j,k)*s5+ &
                    apetz(i,j,k)*s6
                !
                s33(i,j,k)=apztx(i,j,k)*s3+ &
                    apzty(i,j,k)*s5+ &
                    apztz(i,j,k)*s6
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
                rho11(1,i,j,k)=apcsx(i,j,k)*drhodx+ &
                    apcsy(i,j,k)*drhody+ &
                    apcsz(i,j,k)*drhodz
                !
                rho22(1,i,j,k)=apetx(i,j,k)*drhodx+ &
                    apety(i,j,k)*drhody+ &
                    apetz(i,j,k)*drhodz
                !
                rho33(1,i,j,k)=apztx(i,j,k)*drhodx+ &
                    apzty(i,j,k)*drhody+ &
                    apztz(i,j,k)*drhodz

                !
                ! compute controvariant component at cell centroid
                ! uco(i,j,k), it is a flux
                !
                uco(i,j,k)= apcsx(i,j,k)*u(i,j,k)+ &
                    apcsy(i,j,k)*v(i,j,k)+ &
                    apcsz(i,j,k)*w(i,j,k)

                !
                vco(i,j,k)= apetx(i,j,k)*u(i,j,k)+ &
                    apety(i,j,k)*v(i,j,k)+ &
                    apetz(i,j,k)*w(i,j,k)
                !
                wco(i,j,k)= apztx(i,j,k)*u(i,j,k)+ &
                    apzty(i,j,k)*v(i,j,k)+ &
                    apztz(i,j,k)*w(i,j,k)
            !
            enddo
        enddo
    enddo

    !
    !
    ! velocity extrapolation on face 1 and 2, i=0 and i=jx+1
    !
    do k=kparasta,kparaend
        do j=1,jy
            !
            s11(0,j,k)=real(1-ip)*s11(jx,j,k) +  &
                real(ip)*(2.*s11(1,j,k)-s11(2,j,k))
            s12(0,j,k)=real(1-ip)*s12(jx,j,k) +  &
                real(ip)*(2.*s12(1,j,k)-s12(2,j,k))
            s13(0,j,k)=real(1-ip)*s13(jx,j,k) +  &
                real(ip)*(2.*s13(1,j,k)-s13(2,j,k))
            s21(0,j,k)=real(1-ip)*s21(jx,j,k) +  &
                real(ip)*(2.*s21(1,j,k)-s21(2,j,k))
            s22(0,j,k)=real(1-ip)*s22(jx,j,k) +  &
                real(ip)*(2.*s22(1,j,k)-s22(2,j,k))
            s23(0,j,k)=real(1-ip)*s23(jx,j,k) +  &
                real(ip)*(2.*s23(1,j,k)-s23(2,j,k))
            s31(0,j,k)=real(1-ip)*s31(jx,j,k) +  &
                real(ip)*(2.*s31(1,j,k)-s31(2,j,k))
            s32(0,j,k)=real(1-ip)*s32(jx,j,k) +  &
                real(ip)*(2.*s32(1,j,k)-s32(2,j,k))
            s33(0,j,k)=real(1-ip)*s33(jx,j,k) +  &
                real(ip)*(2.*s33(1,j,k)-s33(2,j,k))
            rho11(1,0,j,k)=real(1-ip)*rho11(1,jx,j,k) + &
                real(ip)*(2.*rho11(1,1,j,k)-rho11(1,2,j,k))
            rho22(1,0,j,k)=real(1-ip)*rho22(1,jx,j,k) + &
                real(ip)*(2.*rho22(1,1,j,k)-rho22(1,2,j,k))
            rho33(1,0,j,k)=real(1-ip)*rho33(1,jx,j,k) + &
                real(ip)*(2.*rho33(1,1,j,k)-rho33(1,2,j,k))
            !
            smod(0,j,k)=real(1-ip)*smod(jx,j,k) +  &
                real(ip)*(2.*smod(1,j,k)-smod(2,j,k))
            smodV(0,j,k)=real(1-ip)*smodV(jx,j,k) +  &
                real(ip)*(2.*smodV(1,j,k)-smodV(2,j,k))
            smodH(0,j,k)=real(1-ip)*smodH(jx,j,k) +  &
                real(ip)*(2.*smodH(1,j,k)-smodH(2,j,k))

            uco(0,j,k)=real(1-ip)*uco(jx,j,k) +  &
                real(ip)*(2.*uco(1,j,k)-uco(2,j,k))
            vco(0,j,k)=real(1-ip)*vco(jx,j,k) +  &
                real(ip)*(2.*vco(1,j,k)-vco(2,j,k))
            wco(0,j,k)=real(1-ip)*wco(jx,j,k) +  &
                real(ip)*(2.*wco(1,j,k)-wco(2,j,k))
            !
            ! define border metric terms
            !
            apcsx(0,j,k)=real(1-ip)*apcsx(jx,j,k) + &
                real(ip)*(2.*apcsx(1,j,k)-apcsx(2,j,k))
            apcsy(0,j,k)=real(1-ip)*apcsy(jx,j,k) + &
                real(ip)*(2.*apcsy(1,j,k)-apcsy(2,j,k))
            apcsz(0,j,k)=real(1-ip)*apcsz(jx,j,k) + &
                real(ip)*(2.*apcsz(1,j,k)-apcsz(2,j,k))
            apetx(0,j,k)=real(1-ip)*apetx(jx,j,k) + &
                real(ip)*(2.*apetx(1,j,k)-apetx(2,j,k))
            apety(0,j,k)=real(1-ip)*apety(jx,j,k) + &
                real(ip)*(2.*apety(1,j,k)-apety(2,j,k))
            apetz(0,j,k)=real(1-ip)*apetz(jx,j,k) + &
                real(ip)*(2.*apetz(1,j,k)-apetz(2,j,k))
            apztx(0,j,k)=real(1-ip)*apztx(jx,j,k) + &
                real(ip)*(2.*apztx(1,j,k)-apztx(2,j,k))
            apzty(0,j,k)=real(1-ip)*apzty(jx,j,k) + &
                real(ip)*(2.*apzty(1,j,k)-apzty(2,j,k))
            apztz(0,j,k)=real(1-ip)*apztz(jx,j,k) + &
                real(ip)*(2.*apztz(1,j,k)-apztz(2,j,k))
            !
            apcsx(jx+1,j,k)=real(1-ip)*apcsx(1,j,k) + &
                real(ip)*(2.*apcsx(jx,j,k)-apcsx(jx-1,j,k))
            apcsy(jx+1,j,k)=real(1-ip)*apcsy(1,j,k) + &
                real(ip)*(2.*apcsy(jx,j,k)-apcsy(jx-1,j,k))
            apcsz(jx+1,j,k)=real(1-ip)*apcsz(1,j,k) + &
                real(ip)*(2.*apcsz(jx,j,k)-apcsz(jx-1,j,k))
            apetx(jx+1,j,k)=real(1-ip)*apetx(1,j,k) + &
                real(ip)*(2.*apetx(jx,j,k)-apetx(jx-1,j,k))
            apety(jx+1,j,k)=real(1-ip)*apety(1,j,k) + &
                real(ip)*(2.*apety(jx,j,k)-apety(jx-1,j,k))
            apetz(jx+1,j,k)=real(1-ip)*apetz(1,j,k) + &
                real(ip)*(2.*apetz(jx,j,k)-apetz(jx-1,j,k))
            apztx(jx+1,j,k)=real(1-ip)*apztx(1,j,k) + &
                real(ip)*(2.*apztx(jx,j,k)-apztx(jx-1,j,k))
            apzty(jx+1,j,k)=real(1-ip)*apzty(1,j,k) + &
                real(ip)*(2.*apzty(jx,j,k)-apzty(jx-1,j,k))
            apztz(jx+1,j,k)=real(1-ip)*apztz(1,j,k) + &
                real(ip)*(2.*apztz(jx,j,k)-apztz(jx-1,j,k))
            !
            !
            s11(jx+1,j,k)=real(1-ip)*s11(1,j,k) +  &
                real(ip)*(2.*s11(jx,j,k)-s11(jx-1,j,k))
            s12(jx+1,j,k)=real(1-ip)*s12(1,j,k) +  &
                real(ip)*(2.*s12(jx,j,k)-s12(jx-1,j,k))
            s13(jx+1,j,k)=real(1-ip)*s13(1,j,k) +  &
                real(ip)*(2.*s13(jx,j,k)-s13(jx-1,j,k))
            s21(jx+1,j,k)=real(1-ip)*s21(1,j,k) +  &
                real(ip)*(2.*s21(jx,j,k)-s21(jx-1,j,k))
            s22(jx+1,j,k)=real(1-ip)*s22(1,j,k) +  &
                real(ip)*(2.*s22(jx,j,k)-s22(jx-1,j,k))
            s23(jx+1,j,k)=real(1-ip)*s23(1,j,k) +  &
                real(ip)*(2.*s23(jx,j,k)-s23(jx-1,j,k))
            s31(jx+1,j,k)=real(1-ip)*s31(1,j,k) +  &
                real(ip)*(2.*s31(jx,j,k)-s31(jx-1,j,k))
            s32(jx+1,j,k)=real(1-ip)*s32(1,j,k) +  &
                real(ip)*(2.*s32(jx,j,k)-s32(jx-1,j,k))
            s33(jx+1,j,k)=real(1-ip)*s33(1,j,k) +  &
                real(ip)*(2.*s33(jx,j,k)-s33(jx-1,j,k))
            rho11(1,jx+1,j,k)=real(1-ip)*rho11(1,1,j,k) +  &
                real(ip)*(2.*rho11(1,jx,j,k)-rho11(1,jx-1,j,k))
            rho22(1,jx+1,j,k)=real(1-ip)*rho22(1,1,j,k) +  &
                real(ip)*(2.*rho22(1,jx,j,k)-rho22(1,jx-1,j,k))
            rho33(1,jx+1,j,k)=real(1-ip)*rho33(1,1,j,k) +  &
                real(ip)*(2.*rho33(1,jx,j,k)-rho33(1,jx-1,j,k))
            !
            smod(jx+1,j,k)=real(1-ip)*smod(1,j,k) +  &
                real(ip)*(2.*smod(jx,j,k)-smod(jx-1,j,k))
            smodV(jx+1,j,k)=real(1-ip)*smodV(1,j,k) +  &
                real(ip)*(2.*smodV(jx,j,k)-smodV(jx-1,j,k))
            smodH(jx+1,j,k)=real(1-ip)*smodH(1,j,k) +  &
                real(ip)*(2.*smodH(jx,j,k)-smodH(jx-1,j,k))

            uco(jx+1,j,k)=real(1-ip)*uco(1,j,k) +  &
                real(ip)*(2.*uco(jx,j,k)-uco(jx-1,j,k))
            vco(jx+1,j,k)=real(1-ip)*vco(1,j,k) +  &
                real(ip)*(2.*vco(jx,j,k)-vco(jx-1,j,k))
            wco(jx+1,j,k)=real(1-ip)*wco(1,j,k) +  &
                real(ip)*(2.*wco(jx,j,k)-wco(jx-1,j,k))
        !
        end do
    end do
    !
    !
    ! velocity extrapolation on face 3 and 4, j=0 and j=jy+1
    !
    do k=kparasta,kparaend
        do i=1,jx
            !
            s11(i,0,k)=real(1-jp)*s11(i,jy,k) +  &
                real(jp)*(2.*s11(i,1,k)-s11(i,2,k))
            s12(i,0,k)=real(1-jp)*s12(i,jy,k) +  &
                real(jp)*(2.*s12(i,1,k)-s12(i,2,k))
            s13(i,0,k)=real(1-jp)*s13(i,jy,k) +  &
                real(jp)*(2.*s13(i,1,k)-s13(i,2,k))
            s21(i,0,k)=real(1-jp)*s21(i,jy,k) +  &
                real(jp)*(2.*s21(i,1,k)-s21(i,2,k))
            s22(i,0,k)=real(1-jp)*s22(i,jy,k) +  &
                real(jp)*(2.*s22(i,1,k)-s22(i,2,k))
            s23(i,0,k)=real(1-jp)*s23(i,jy,k) +  &
                real(jp)*(2.*s23(i,1,k)-s23(i,2,k))
            s31(i,0,k)=real(1-jp)*s31(i,jy,k) +  &
                real(jp)*(2.*s31(i,1,k)-s31(i,2,k))
            s32(i,0,k)=real(1-jp)*s32(i,jy,k) +  &
                real(jp)*(2.*s32(i,1,k)-s32(i,2,k))
            s33(i,0,k)=real(1-jp)*s33(i,jy,k) +  &
                real(jp)*(2.*s33(i,1,k)-s33(i,2,k))
            rho11(1,i,0,k)=real(1-jp)*rho11(1,i,jy,k) + &
                real(jp)*(2.*rho11(1,i,1,k)-rho11(1,i,2,k))
            rho22(1,i,0,k)=real(1-jp)*rho22(1,i,jy,k) + &
                real(jp)*(2.*rho22(1,i,1,k)-rho22(1,i,2,k))
            rho33(1,i,0,k)=real(1-jp)*rho33(1,i,jy,k) + &
                real(jp)*(2.*rho33(1,i,1,k)-rho33(1,i,2,k))
            !
            smod(i,0,k)=real(1-jp)*smod(i,jy,k) +  &
                real(jp)*(2.*smod(i,1,k)-smod(i,2,k))
            smodV(i,0,k)=real(1-jp)*smodV(i,jy,k) +  &
                real(jp)*(2.*smodV(i,1,k)-smodV(i,2,k))
            smodH(i,0,k)=real(1-jp)*smodH(i,jy,k) +  &
                real(jp)*(2.*smodH(i,1,k)-smodH(i,2,k))

            uco(i,0,k)=real(1-jp)*uco(i,jy,k) +  &
                real(jp)*(2.*uco(i,1,k)-uco(i,2,k))
            vco(i,0,k)=real(1-jp)*vco(i,jy,k) +  &
                real(jp)*(2.*vco(i,1,k)-vco(i,2,k))
            wco(i,0,k)=real(1-jp)*wco(i,jy,k) +  &
                real(jp)*(2.*wco(i,1,k)-wco(i,2,k))
            !
            ! define border metric term
            !
            apcsx(i,0,k)=real(1-jp)*apcsx(i,jy,k) + &
                real(jp)*(2.*apcsx(i,1,k)-apcsx(i,2,k))
            apcsy(i,0,k)=real(1-jp)*apcsy(i,jy,k) + &
                real(jp)*(2.*apcsy(i,1,k)-apcsy(i,2,k))
            apcsz(i,0,k)=real(1-jp)*apcsz(i,jy,k) + &
                real(jp)*(2.*apcsz(i,1,k)-apcsz(i,2,k))
            apetx(i,0,k)=real(1-jp)*apetx(i,jy,k) + &
                real(jp)*(2.*apetx(i,1,k)-apetx(i,2,k))
            apety(i,0,k)=real(1-jp)*apety(i,jy,k) + &
                real(jp)*(2.*apety(i,1,k)-apety(i,2,k))
            apetz(i,0,k)=real(1-jp)*apetz(i,jy,k) + &
                real(jp)*(2.*apetz(i,1,k)-apetz(i,2,k))
            apztx(i,0,k)=real(1-jp)*apztx(i,jy,k) + &
                real(jp)*(2.*apztx(i,1,k)-apztx(i,2,k))
            apzty(i,0,k)=real(1-jp)*apzty(i,jy,k) + &
                real(jp)*(2.*apzty(i,1,k)-apzty(i,2,k))
            apztz(i,0,k)=real(1-jp)*apztz(i,jy,k) + &
                real(jp)*(2.*apztz(i,1,k)-apztz(i,2,k))
            !
            apcsx(i,jy+1,k)=real(1-jp)*apcsx(i,1,k) + &
                real(jp)*(2.*apcsx(i,jy,k)-apcsx(i,jy-1,k))
            apcsy(i,jy+1,k)=real(1-jp)*apcsy(i,1,k) + &
                real(jp)*(2.*apcsy(i,jy,k)-apcsy(i,jy-1,k))
            apcsz(i,jy+1,k)=real(1-jp)*apcsz(i,1,k) + &
                real(jp)*(2.*apcsz(i,jy,k)-apcsz(i,jy-1,k))
            apetx(i,jy+1,k)=real(1-jp)*apetx(i,1,k) + &
                real(jp)*(2.*apetz(i,jy,k)-apetx(i,jy-1,k))
            apety(i,jy+1,k)=real(1-jp)*apety(i,1,k) + &
                real(jp)*(2.*apety(i,jy,k)-apety(i,jy-1,k))
            apetz(i,jy+1,k)=real(1-jp)*apetz(i,1,k) + &
                real(jp)*(2.*apetz(i,jy,k)-apetz(i,jy-1,k))
            apztx(i,jy+1,k)=real(1-jp)*apztx(i,1,k) + &
                real(jp)*(2.*apztx(i,jy,k)-apztx(i,jy-1,k))
            apzty(i,jy+1,k)=real(1-jp)*apzty(i,1,k) + &
                real(jp)*(2.*apzty(i,jy,k)-apzty(i,jy-1,k))
            apztz(i,jy+1,k)=real(1-jp)*apztz(i,1,k) + &
                real(jp)*(2.*apztz(i,jy,k)-apztz(i,jy-1,k))
            !
            !
            s11(i,jy+1,k)=real(1-jp)*s11(i,1,k) +  &
                real(jp)*(2.*s11(i,jy,k)-s11(i,jy-1,k))
            s12(i,jy+1,k)=real(1-jp)*s12(i,1,k) +  &
                real(jp)*(2.*s12(i,jy,k)-s12(i,jy-1,k))
            s13(i,jy+1,k)=real(1-jp)*s13(i,1,k) +  &
                real(jp)*(2.*s13(i,jy,k)-s13(i,jy-1,k))
            s21(i,jy+1,k)=real(1-jp)*s21(i,1,k) +  &
                real(jp)*(2.*s21(i,jy,k)-s21(i,jy-1,k))
            s22(i,jy+1,k)=real(1-jp)*s22(i,1,k) +  &
                real(jp)*(2.*s22(i,jy,k)-s22(i,jy-1,k))
            s23(i,jy+1,k)=real(1-jp)*s23(i,1,k) +  &
                real(jp)*(2.*s23(i,jy,k)-s23(i,jy-1,k))
            s31(i,jy+1,k)=real(1-jp)*s31(i,1,k) +  &
                real(jp)*(2.*s31(i,jy,k)-s31(i,jy-1,k))
            s32(i,jy+1,k)=real(1-jp)*s32(i,1,k) +  &
                real(jp)*(2.*s32(i,jy,k)-s32(i,jy-1,k))
            s33(i,jy+1,k)=real(1-jp)*s33(i,1,k) +  &
                real(jp)*(2.*s33(i,jy,k)-s33(i,jy-1,k))
            rho11(1,i,jy+1,k)=real(1-jp)*rho11(1,i,1,k) +  &
                real(jp)*(2.*rho11(1,i,jy,k)-rho11(1,i,jy-1,k))
            rho22(1,i,jy+1,k)=real(1-jp)*rho22(1,i,1,k) +  &
                real(jp)*(2.*rho22(1,i,jy,k)-rho22(1,i,jy-1,k))
            rho33(1,i,jy+1,k)=real(1-jp)*rho33(1,i,1,k) +  &
                real(jp)*(2.*rho33(1,i,jy,k)-rho33(1,i,jy-1,k))
            !
            smod(i,jy+1,k)=real(1-jp)*smod(i,1,k) +  &
                real(jp)*(2.*smod(i,jy,k)-smod(i,jy-1,k))
            smodV(i,jy+1,k)=real(1-jp)*smodV(i,1,k) +  &
                real(jp)*(2.*smodV(i,jy,k)-smodV(i,jy-1,k))
            smodH(i,jy+1,k)=real(1-jp)*smodH(i,1,k) +  &
                real(jp)*(2.*smodH(i,jy,k)-smodH(i,jy-1,k))

            uco(i,jy+1,k)=real(1-jp)*uco(i,1,k) +  &
                real(jp)*(2.*uco(i,jy,k)-uco(i,jy-1,k))
            vco(i,jy+1,k)=real(1-jp)*vco(i,1,k) +  &
                real(jp)*(2.*vco(i,jy,k)-vco(i,jy-1,k))
            wco(i,jy+1,k)=real(1-jp)*wco(i,1,k) +  &
                real(jp)*(2.*wco(i,jy,k)-wco(i,jy-1,k))
        !
        end do
    end do

    ! now I need to define periodicity for k direction
    !
    ! velocity extrapolation on face 5 and 6, k=0 and k=jz+1
    !
    ! buildd a buff vector which contains all jx*jy planes
    ! for information exchange between P0 and Pn-1
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    if (myid.eq.nproc-1) then

        call buffer1g(s11,1,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s12,2,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s13,3,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s21,4,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s22,5,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s23,6,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s31,7,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s32,8,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(s33,9,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(rho11,10,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(rho22,11,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(rho33,12,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(smod,13,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(uco,14,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(vco,15,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(wco,16,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apcsx,17,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apcsy,18,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apcsz,19,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apetx,20,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apety,21,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apetz,22,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apztx,23,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apzty,24,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(apztz,25,jz,myid,nproc,kparasta,kparaend)
        !     for anisotropic cell
        call buffer1g(smodV,26,jz,myid,nproc,kparasta,kparaend)
        call buffer1g(smodH,27,jz,myid,nproc,kparasta,kparaend)

    else if (myid.eq.0) then

        call buffer1g(s11,1,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s12,2,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s13,3,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s21,4,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s22,5,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s23,6,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s31,7,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s32,8,1,myid,nproc,kparasta,kparaend)
        call buffer1g(s33,9,1,myid,nproc,kparasta,kparaend)
        call buffer1g(rho11,10,1,myid,nproc,kparasta,kparaend)
        call buffer1g(rho22,11,1,myid,nproc,kparasta,kparaend)
        call buffer1g(rho33,12,1,myid,nproc,kparasta,kparaend)
        call buffer1g(smod,13,1,myid,nproc,kparasta,kparaend)
        call buffer1g(uco,14,1,myid,nproc,kparasta,kparaend)
        call buffer1g(vco,15,1,myid,nproc,kparasta,kparaend)
        call buffer1g(wco,16,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apcsx,17,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apcsy,18,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apcsz,19,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apetx,20,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apety,21,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apetz,22,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apztx,23,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apzty,24,1,myid,nproc,kparasta,kparaend)
        call buffer1g(apztz,25,1,myid,nproc,kparasta,kparaend)
        !     for anisotropic cell
        call buffer1g(smodV,26,1,myid,nproc,kparasta,kparaend)
        call buffer1g(smodH,27,1,myid,nproc,kparasta,kparaend)

    endif
    !
    ! send procedure so that Po knows k=jz and Pn-1 knows k=1

    if (myid.eq.nproc-1) then

        call MPI_SENDRECV(sbuff1(1),27*(jx+2)*(jy+2),MPI_REAL_SD,0,9001, &
            rbuff1(1),27*(jx+2)*(jy+2),MPI_REAL_SD,0,8001, &
            MPI_COMM_WORLD,status,ierr)

    else if (myid.eq.0) then

        call MPI_SENDRECV(sbuff1(1),27*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1, &
            8001, &
            rbuff1(1),27*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1, &
            9001, &
            MPI_COMM_WORLD,status,ierr)

    endif

    if (myid.eq.0) then

        call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano4,4,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano5,5,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano6,6,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano7,7,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano8,8,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano9,9,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano10,10,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano11,11,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano12,12,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano13,13,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano14,14,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano15,15,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano16,16,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano17,17,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano18,18,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano19,19,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano20,20,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano21,21,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano22,22,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano23,23,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano24,24,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano25,25,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano26,26,myid,nproc, &
            kparasta,kparaend)
        call buffer2gg(piano27,27,myid,nproc, &
            kparasta,kparaend)


    else if (myid.eq.nproc-1) then

        call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano4,4,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano5,5,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano6,6,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano7,7,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano8,8,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano9,9,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano10,10,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano11,11,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano12,12,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano13,13,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano14,14,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano15,15,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano16,16,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano17,17,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano18,18,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano19,19,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano20,20,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano21,21,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano22,22,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano23,23,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano24,24,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano25,25,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano26,26,myid,nproc,kparasta,kparaend)
        call buffer2gg(piano27,27,myid,nproc,kparasta,kparaend)

    endif

    ! now P0   knows the plane k=jz
    !     Pn-1 knows the plane k=1

    if(myid.eq.0)then

        do j=1,jy
            do i=1,jx

                s11(i,j,0)=real(1-kp)*piano1(i,j) +  &
                    real(kp)*(2.*s11(i,j,1)-s11(i,j,2))
                s12(i,j,0)=real(1-kp)*piano2(i,j) +  &
                    real(kp)*(2.*s12(i,j,1)-s12(i,j,2))
                s13(i,j,0)=real(1-kp)*piano3(i,j) +  &
                    real(kp)*(2.*s13(i,j,1)-s13(i,j,2))
                s21(i,j,0)=real(1-kp)*piano4(i,j) +  &
                    real(kp)*(2.*s21(i,j,1)-s21(i,j,2))
                s22(i,j,0)=real(1-kp)*piano5(i,j) +  &
                    real(kp)*(2.*s22(i,j,1)-s22(i,j,2))
                s23(i,j,0)=real(1-kp)*piano6(i,j) +  &
                    real(kp)*(2.*s23(i,j,1)-s23(i,j,2))
                s31(i,j,0)=real(1-kp)*piano7(i,j) +  &
                    real(kp)*(2.*s31(i,j,1)-s31(i,j,2))
                s32(i,j,0)=real(1-kp)*piano8(i,j) +  &
                    real(kp)*(2.*s32(i,j,1)-s32(i,j,2))
                s33(i,j,0)=real(1-kp)*piano9(i,j) +  &
                    real(kp)*(2.*s33(i,j,1)-s33(i,j,2))
                rho11(1,i,j,0)=real(1-kp)*piano10(i,j) + &
                    real(kp)*(2.*rho11(1,i,j,1)-rho11(1,i,j,2))
                rho22(1,i,j,0)=real(1-kp)*piano11(i,j) + &
                    real(kp)*(2.*rho22(1,i,j,1)-rho22(1,i,j,2))
                rho33(1,i,j,0)=real(1-kp)*piano12(i,j) + &
                    real(kp)*(2.*rho33(1,i,j,1)-rho33(1,i,j,2))

                smod(i,j,0)=real(1-kp)*piano13(i,j) +  &
                    real(kp)*(2.*smod(i,j,1)-smod(i,j,2))
                smodV(i,j,0)=real(1-kp)*piano26(i,j) +  &
                    real(kp)*(2.*smodV(i,j,1)-smodV(i,j,2))
                smodH(i,j,0)=real(1-kp)*piano27(i,j) +  &
                    real(kp)*(2.*smodH(i,j,1)-smodH(i,j,2))



                uco(i,j,0)=real(1-kp)*piano14(i,j) +  &
                    real(kp)*(2.*uco(i,j,1)-uco(i,j,2))
                vco(i,j,0)=real(1-kp)*piano15(i,j) +  &
                    real(kp)*(2.*vco(i,j,1)-vco(i,j,2))
                wco(i,j,0)=real(1-kp)*piano16(i,j) +  &
                    real(kp)*(2.*wco(i,j,1)-wco(i,j,2))
                !
                ! define border metric term
                !
                apcsx(i,j,0)=real(1-kp)*piano17(i,j) + &
                    real(kp)*(2.*apcsx(i,j,1)-apcsx(i,j,2))
                apcsy(i,j,0)=real(1-kp)*piano18(i,j) + &
                    real(kp)*(2.*apcsy(i,j,1)-apcsy(i,j,2))
                apcsz(i,j,0)=real(1-kp)*piano19(i,j) + &
                    real(kp)*(2.*apcsz(i,j,1)-apcsz(i,j,2))
                apetx(i,j,0)=real(1-kp)*piano20(i,j) + &
                    real(kp)*(2.*apetx(i,j,1)-apetx(i,j,2))
                apety(i,j,0)=real(1-kp)*piano21(i,j) + &
                    real(kp)*(2.*apety(i,j,1)-apety(i,j,2))
                apetz(i,j,0)=real(1-kp)*piano22(i,j) + &
                    real(kp)*(2.*apetz(i,j,1)-apetz(i,j,2))
                apztx(i,j,0)=real(1-kp)*piano23(i,j) + &
                    real(kp)*(2.*apztx(i,j,1)-apztx(i,j,2))
                apzty(i,j,0)=real(1-kp)*piano24(i,j) + &
                    real(kp)*(2.*apzty(i,j,1)-apzty(i,j,2))
                apztz(i,j,0)=real(1-kp)*piano25(i,j) + &
                    real(kp)*(2.*apztz(i,j,1)-apztz(i,j,2))


            end do
        end do

    endif
    !
    !
    if(myid.eq.nproc-1)then

        do j=1,jy
            do i=1,jx


                s11(i,j,jz+1)=real(1-kp)*piano1(i,j) +  &
                    real(kp)*(2.*s11(i,j,jz)-s11(i,j,jz-1))
                s12(i,j,jz+1)=real(1-kp)*piano2(i,j) +  &
                    real(kp)*(2.*s12(i,j,jz)-s12(i,j,jz-1))
                s13(i,j,jz+1)=real(1-kp)*piano3(i,j) +  &
                    real(kp)*(2.*s13(i,j,jz)-s13(i,j,jz-1))
                s21(i,j,jz+1)=real(1-kp)*piano4(i,j) +  &
                    real(kp)*(2.*s21(i,j,jz)-s21(i,j,jz-1))
                s22(i,j,jz+1)=real(1-kp)*piano5(i,j) +  &
                    real(kp)*(2.*s22(i,j,jz)-s22(i,j,jz-1))
                s23(i,j,jz+1)=real(1-kp)*piano6(i,j) +  &
                    real(kp)*(2.*s23(i,j,jz)-s23(i,j,jz-1))
                s31(i,j,jz+1)=real(1-kp)*piano7(i,j) +  &
                    real(kp)*(2.*s31(i,j,jz)-s31(i,j,jz-1))
                s32(i,j,jz+1)=real(1-kp)*piano8(i,j) +  &
                    real(kp)*(2.*s32(i,j,jz)-s32(i,j,jz-1))
                s33(i,j,jz+1)=real(1-kp)*piano9(i,j) +  &
                    real(kp)*(2.*s33(i,j,jz)-s33(i,j,jz-1))
                rho11(1,i,j,jz+1)=real(1-kp)*piano10(i,j) + &
                    real(kp)*(2.*rho11(1,i,j,jz)-rho11(1,i,j,jz-1))
                rho22(1,i,j,jz+1)=real(1-kp)*piano11(i,j) + &
                    real(kp)*(2.*rho22(1,i,j,jz)-rho22(1,i,j,jz-1))
                rho33(1,i,j,jz+1)=real(1-kp)*piano12(i,j) + &
                    real(kp)*(2.*rho33(1,i,j,jz)-rho33(1,i,j,jz-1))
                !

                smod(i,j,jz+1)=real(1-kp)*piano13(i,j) +  &
                    real(kp)*(2.*smod(i,j,jz)-smod(i,j,jz-1))
                smodV(i,j,jz+1)=real(1-kp)*piano26(i,j) +  &
                    real(kp)*(2.*smodV(i,j,jz)-smodV(i,j,jz-1))
                smodH(i,j,jz+1)=real(1-kp)*piano27(i,j) +  &
                    real(kp)*(2.*smodH(i,j,jz)-smodH(i,j,jz-1))

                uco(i,j,jz+1)=real(1-kp)*piano14(i,j) +  &
                    real(kp)*(2.*uco(i,j,jz)-uco(i,j,jz-1))
                vco(i,j,jz+1)=real(1-kp)*piano15(i,j) +  &
                    real(kp)*(2.*vco(i,j,jz)-vco(i,j,jz-1))
                wco(i,j,jz+1)=real(1-kp)*piano16(i,j) +  &
                    real(kp)*(2.*wco(i,j,jz)-wco(i,j,jz-1))
                !
                ! define metric border term
                !
                apcsx(i,j,jz+1)=real(1-kp)*piano17(i,j) + &
                    real(kp)*(2.*apcsx(i,j,jz)-apcsx(i,j,jz-1))
                apcsy(i,j,jz+1)=real(1-kp)*piano18(i,j) + &
                    real(kp)*(2.*apcsy(i,j,jz)-apcsy(i,j,jz-1))
                apcsz(i,j,jz+1)=real(1-kp)*piano19(i,j) + &
                    real(kp)*(2.*apcsz(i,j,jz)-apcsz(i,j,jz-1))
                apetx(i,j,jz+1)=real(1-kp)*piano20(i,j) + &
                    real(kp)*(2.*apetx(i,j,jz)-apetx(i,j,jz-1))
                apety(i,j,jz+1)=real(1-kp)*piano21(i,j) + &
                    real(kp)*(2.*apety(i,j,jz)-apety(i,j,jz-1))
                apetz(i,j,jz+1)=real(1-kp)*piano22(i,j) + &
                    real(kp)*(2.*piano23(i,j)-apetz(i,j,jz-1))
                apztx(i,j,jz+1)=real(1-kp)*piano23(i,j) + &
                    real(kp)*(2.*apztx(i,j,jz)-apztx(i,j,jz-1))
                apzty(i,j,jz+1)=real(1-kp)*piano24(i,j) + &
                    real(kp)*(2.*apzty(i,j,jz)-apzty(i,j,jz-1))
                apztz(i,j,jz+1)=real(1-kp)*piano25(i,j) + &
                    real(kp)*(2.*apztz(i,j,jz)-apztz(i,j,jz-1))

            end do
        end do

    endif
    !
    ! check
    !
    if (debugg.eq.1) then

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(3000+myid+ttt,*)apcsx(i,j,k),apcsy(i,j,k),apcsz(i,j,k)
                    write(3005+myid+ttt,*)apetx(i,j,k),apety(i,j,k),apetz(i,j,k)
                    write(3010+myid+ttt,*)apztx(i,j,k),apzty(i,j,k),apztz(i,j,k)
                    write(3015+myid+ttt,*)uco(i,j,k),vco(i,j,k),wco(i,j,k)
                    write(3020+myid+ttt,*)s11(i,j,k),s12(i,j,k),s13(i,j,k)
                    write(3025+myid+ttt,*)s21(i,j,k),s22(i,j,k),s23(i,j,k)
                    write(3030+myid+ttt,*)s31(i,j,k),s32(i,j,k),s33(i,j,k)
                    write(3035+myid+ttt,*)rho11(1,i,j,k),rho22(1,i,j,k),rho33(1,i,j,k)
                end do
            end do
        end do

    end if

    !-----------------------------------------------------------------------
    !***********************************************************************
    !-----------------------------------------------------------------------
    ! start filtering procedure of all component previously computed
    ! (also for cross product), to compute the matrix Mij and Lij (cartesian????),
    ! they are necessary to find the model constant C for the subgrid stress
    ! (SGS) model
    !
    ! Mij = Df^2 * |S|f * Sijf - D^2 * (|S|Sij)f
    !                             (D=2*giac^-1/3 e Df=sqrt(6)*D)
    !
    ! Lij = (ui * uj)f - uif * ujf
    !
    ! I need to compute the controvariant form of these quantities
    !
    ! Mik = Df^2 * |S|f * Sikf - D^2 * (|S|Sik)f
    !
    ! Lik = (ui * Uk)f - uif * Ukf
    !
    ! same computation is done for density
    !

    call filter01(s11,m11,myid,nproc,kparasta,kparaend)
    call filter01(s12,m12,myid,nproc,kparasta,kparaend)
    call filter01(s13,m13,myid,nproc,kparasta,kparaend)
    call filter01(s21,m21,myid,nproc,kparasta,kparaend)
    call filter01(s22,m22,myid,nproc,kparasta,kparaend)
    call filter01(s23,m23,myid,nproc,kparasta,kparaend)
    call filter01(s31,m31,myid,nproc,kparasta,kparaend)
    call filter01(s32,m32,myid,nproc,kparasta,kparaend)
    call filter01(s33,m33,myid,nproc,kparasta,kparaend)
    call filter01(rho11,mrho11,myid,nproc,kparasta,kparaend)
    call filter01(rho22,mrho22,myid,nproc,kparasta,kparaend)
    call filter01(rho33,mrho33,myid,nproc,kparasta,kparaend)

    call filter02(s11,smod,l11,myid,nproc,kparasta,kparaend)
    call filter02(s12,smod,l12,myid,nproc,kparasta,kparaend)
    call filter02(s13,smod,l13,myid,nproc,kparasta,kparaend)
    call filter02(s21,smod,l21,myid,nproc,kparasta,kparaend)
    call filter02(s22,smod,l22,myid,nproc,kparasta,kparaend)
    call filter02(s23,smod,l23,myid,nproc,kparasta,kparaend)
    call filter02(s31,smod,l31,myid,nproc,kparasta,kparaend)
    call filter02(s32,smod,l32,myid,nproc,kparasta,kparaend)
    call filter02(s33,smod,l33,myid,nproc,kparasta,kparaend)
    call filter02(rho11,smod,lrho11,myid,nproc,kparasta,kparaend)
    call filter02(rho22,smod,lrho22,myid,nproc,kparasta,kparaend)
    call filter02(rho33,smod,lrho33,myid,nproc,kparasta,kparaend)


    !     check
    if (debugg.eq.1) then

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(3100+myid+ttt,*)m11(i,j,k),m12(i,j,k),m13(i,j,k)
                    write(3105+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                    write(3110+myid+ttt,*)m31(i,j,k),m32(i,j,k),m33(i,j,k)
                    write(3115+myid+ttt,*)mrho11(i,j,k),mrho22(i,j,k),mrho33(i,j,k)
                    write(3120+myid+ttt,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                    write(3125+myid+ttt,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                    write(3130+myid+ttt,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                    write(3135+myid+ttt,*)lrho11(i,j,k),lrho22(i,j,k),lrho33(i,j,k)
                end do
            end do
        end do

    end if

    !
    ! now I need to transfer the "ghost" plane to include periodicity
    ! in order to complete the filtering procedure in zita (k+1 and k-1)
    !
    ! first kparasta
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    call buffer1g(m11,1,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m12,2,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m13,3,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m21,4,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m22,5,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m23,6,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m31,7,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m32,8,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m33,9,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(mrho11,10,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(mrho22,11,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(mrho33,12,kparasta,myid,nproc, &
        kparasta,kparaend)

    call buffer1g(l11,13,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l12,14,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l13,15,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l21,16,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l22,17,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l23,18,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l31,19,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l32,20,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l33,21,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho11,22,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho22,23,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho33,24,kparasta,myid,nproc, &
        kparasta,kparaend)
    !
    ! if I put pe and not pem I have implicit periodicity on k
    !
    if (kp.eq.0) then

        call MPI_SENDRECV(sbuff1(1),24*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,tagls, &
            rbuff1(1),24*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrr, &
            MPI_COMM_WORLD,status,ierr)


    else if (kp.eq.1) then


        if(leftpem /= MPI_PROC_NULL) then

            call MPI_SEND(sbuff1(1),24*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,tagls, &
                MPI_COMM_WORLD,ierr)

        endif

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),24*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        endif

        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif

    endif

    call buffer2(rbuff1,m11,1,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m12,2,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m13,3,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m21,4,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m22,5,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m23,6,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m31,7,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m32,8,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m33,9,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,mrho11,10,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,mrho22,11,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,mrho33,12,kparaend+1,myid,nproc, &
        kparasta,kparaend)

    call buffer2(rbuff1,l11,13,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l12,14,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l13,15,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l21,16,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l22,17,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l23,18,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l31,19,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l32,20,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l33,21,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho11,22,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho22,23,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho33,24,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    !
    ! now kparaend
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    call buffer1g(m11,1,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m12,2,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m13,3,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m21,4,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m22,5,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m23,6,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m31,7,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m32,8,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m33,9,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(mrho11,10,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(mrho22,11,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(mrho33,12,kparaend,myid,nproc, &
        kparasta,kparaend)

    call buffer1g(l11,13,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l12,14,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l13,15,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l21,16,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l22,17,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l23,18,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l31,19,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l32,20,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l33,21,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho11,22,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho22,23,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho33,24,kparaend,myid,nproc, &
        kparasta,kparaend)

    if (kp.eq.0) then

        call MPI_SENDRECV(sbuff1(1),24*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrs, &
            rbuff1(1),24*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)

    else if (kp.eq.1) then

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SEND(sbuff1(1),24*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrs, &
                MPI_COMM_WORLD,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),24*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,taglr, &
                MPI_COMM_WORLD,status,ierr)
        endif

        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif


    endif

    call buffer2(rbuff1,m11,1,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m12,2,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m13,3,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m21,4,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m22,5,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m23,6,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m31,7,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m32,8,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m33,9,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,mrho11,10,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,mrho22,11,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,mrho33,12,kparasta-1,myid,nproc, &
        kparasta,kparaend)

    call buffer2(rbuff1,l11,13,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l12,14,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l13,15,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l21,16,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l22,17,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l23,18,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l31,19,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l32,20,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l33,21,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho11,22,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho22,23,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho33,24,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    !
    ! periodicity in k (implemented in filter03)
    ! now I complete the filtering procedure in zita
    !
    call filter03(m11,s11f,kparasta,kparaend,myid,nproc)
    call filter03(m12,s12f,kparasta,kparaend,myid,nproc)
    call filter03(m13,s13f,kparasta,kparaend,myid,nproc)
    call filter03(m21,s21f,kparasta,kparaend,myid,nproc)
    call filter03(m22,s22f,kparasta,kparaend,myid,nproc)
    call filter03(m23,s23f,kparasta,kparaend,myid,nproc)
    call filter03(m31,s31f,kparasta,kparaend,myid,nproc)
    call filter03(m32,s32f,kparasta,kparaend,myid,nproc)
    call filter03(m33,s33f,kparasta,kparaend,myid,nproc)
    call filter03(mrho11,rho11f,kparasta,kparaend,myid,nproc)
    call filter03(mrho22,rho22f,kparasta,kparaend,myid,nproc)
    call filter03(mrho33,rho33f,kparasta,kparaend,myid,nproc)

    call filter03(l11,smods11f,kparasta,kparaend,myid,nproc)
    call filter03(l12,smods12f,kparasta,kparaend,myid,nproc)
    call filter03(l13,smods13f,kparasta,kparaend,myid,nproc)
    call filter03(l21,smods21f,kparasta,kparaend,myid,nproc)
    call filter03(l22,smods22f,kparasta,kparaend,myid,nproc)
    call filter03(l23,smods23f,kparasta,kparaend,myid,nproc)
    call filter03(l31,smods31f,kparasta,kparaend,myid,nproc)
    call filter03(l32,smods32f,kparasta,kparaend,myid,nproc)
    call filter03(l33,smods33f,kparasta,kparaend,myid,nproc)
    call filter03(lrho11,smodrho11f,kparasta,kparaend,myid,nproc)
    call filter03(lrho22,smodrho22f,kparasta,kparaend,myid,nproc)
    call filter03(lrho33,smodrho33f,kparasta,kparaend,myid,nproc)



    ! check
    if (debugg.eq.1) then

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(3200+myid+ttt,*)s11f(i,j,k),s12f(i,j,k),s13f(i,j,k)
                    write(3205+myid+ttt,*)s21f(i,j,k),s22f(i,j,k),s23f(i,j,k)
                    write(3210+myid+ttt,*)s31f(i,j,k),s32f(i,j,k),s33f(i,j,k)
                    write(3215+myid+ttt,*)rho11f(i,j,k),rho22f(i,j,k),rho33f(i,j,k)
                    write(3220+myid+ttt,*)smods11f(i,j,k),smods12f(i,j,k), &
                        smods13f(i,j,k)
                    write(3225+myid+ttt,*)smods21f(i,j,k),smods22f(i,j,k), &
                        smods23f(i,j,k)
                    write(3230+myid+ttt,*)smods31f(i,j,k),smods32f(i,j,k), &
                        smods33f(i,j,k)
                    write(3235+myid+ttt,*)smodrho11f(i,j,k),smodrho22f(i,j,k), &
                        smodrho33f(i,j,k)
                end do
            end do
        end do

    end if

    ! still need to compute |S| filtered (smodf)
    ! not yet done because it needs all filtered components of Sik

    !
    ! METRIC FILTERING
    !
    call filter01(apcsx,m11,myid,nproc,kparasta,kparaend)
    call filter01(apcsy,m12,myid,nproc,kparasta,kparaend)
    call filter01(apcsz,m13,myid,nproc,kparasta,kparaend)
    call filter01(apetx,m21,myid,nproc,kparasta,kparaend)
    call filter01(apety,m22,myid,nproc,kparasta,kparaend)
    call filter01(apetz,m23,myid,nproc,kparasta,kparaend)
    call filter01(apztx,m31,myid,nproc,kparasta,kparaend)
    call filter01(apzty,m32,myid,nproc,kparasta,kparaend)
    call filter01(apztz,m33,myid,nproc,kparasta,kparaend)
    !
    ! now I need to exchange the "ghost" plane to include periodicity
    ! to complete the filtering procedure on metric in zita(k+1 and k-1)

    !
    ! first kparasta
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    call buffer1g(m11,1,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m12,2,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m13,3,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m21,4,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m22,5,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m23,6,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m31,7,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m32,8,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m33,9,kparasta,myid,nproc,kparasta,kparaend)
    !
    ! if I put pe and not pem I implicitly have periodicity in k
    !
    if (kp.eq.0) then

        call MPI_SENDRECV(sbuff1(1),9*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,tagls, &
            rbuff1(1),9*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrr, &
            MPI_COMM_WORLD,status,ierr)

    else if (kp.eq.1) then

        if(leftpem /= MPI_PROC_NULL) then
            call MPI_SEND(sbuff1(1),9*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,tagls, &
                MPI_COMM_WORLD,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),9*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrr, &
                MPI_COMM_WORLD,status,ierr)
        endif

        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif


    endif


    call buffer2(rbuff1,m11,1,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m12,2,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m13,3,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m21,4,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m22,5,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m23,6,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m31,7,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m32,8,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m33,9,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    !
    ! now kparaend
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    call buffer1g(m11,1,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m12,2,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m13,3,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m21,4,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m22,5,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m23,6,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m31,7,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m32,8,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m33,9,kparaend,myid,nproc,kparasta,kparaend)

    if (kp.eq.0) then

        call MPI_SENDRECV(sbuff1(1),9*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrs, &
            rbuff1(1),9*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)

    else if (kp.eq.1) then

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SEND(sbuff1(1),9*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrs, &
                MPI_COMM_WORLD,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),9*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,taglr, &
                MPI_COMM_WORLD,status,ierr)
        endif

        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif


    endif

    call buffer2(rbuff1,m11,1,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m12,2,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m13,3,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m21,4,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m22,5,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m23,6,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m31,7,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m32,8,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m33,9,kparasta-1,myid,nproc, &
        kparasta,kparaend)

    if (debugg.eq.1) then

        do k=kparasta-1,kparaend+1
            do j=1,n2
                do i=1,n1
                    write(3300+myid+ttt,*)m11(i,j,k),m12(i,j,k),m13(i,j,k)
                    write(3305+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                    write(3310+myid+ttt,*)m31(i,j,k),m32(i,j,k),m33(i,j,k)
                end do
            end do
        end do

    end if

    !
    ! periodicity in k (implemented in filter03)
    ! now I complete the filtering on zita
    !
    call filter03(m11,lmf11,kparasta,kparaend,myid,nproc)
    call filter03(m12,lmf12,kparasta,kparaend,myid,nproc)
    call filter03(m13,lmf13,kparasta,kparaend,myid,nproc)
    call filter03(m21,lmf21,kparasta,kparaend,myid,nproc)
    call filter03(m22,lmf22,kparasta,kparaend,myid,nproc)
    call filter03(m23,lmf23,kparasta,kparaend,myid,nproc)
    call filter03(m31,lmf31,kparasta,kparaend,myid,nproc)
    call filter03(m32,lmf32,kparasta,kparaend,myid,nproc)
    call filter03(m33,lmf33,kparasta,kparaend,myid,nproc)

    !
    ! I find the filtered cartesian component from the controvariant
    ! then I compute |S| filtered as a function of theese last
    !

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = s11f(i,j,k)
                pc2(i,j,k)  = s12f(i,j,k)
                pc3(i,j,k)  = s13f(i,j,k)
                ap11(i,j,k) = lmf11(i,j,k)
                ap12(i,j,k) = lmf12(i,j,k)
                ap13(i,j,k) = lmf13(i,j,k)
                ap21(i,j,k) = lmf21(i,j,k)
                ap22(i,j,k) = lmf22(i,j,k)
                ap23(i,j,k) = lmf23(i,j,k)
                ap31(i,j,k) = lmf31(i,j,k)
                ap32(i,j,k) = lmf32(i,j,k)
                ap33(i,j,k) = lmf33(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil11(i,j,k) = pp1(i,j,k)
                fil12(i,j,k) = pp2(i,j,k)
                fil13(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = s21f(i,j,k)
                pc2(i,j,k)  = s22f(i,j,k)
                pc3(i,j,k)  = s23f(i,j,k)
                ap11(i,j,k) = lmf11(i,j,k)
                ap12(i,j,k) = lmf12(i,j,k)
                ap13(i,j,k) = lmf13(i,j,k)
                ap21(i,j,k) = lmf21(i,j,k)
                ap22(i,j,k) = lmf22(i,j,k)
                ap23(i,j,k) = lmf23(i,j,k)
                ap31(i,j,k) = lmf31(i,j,k)
                ap32(i,j,k) = lmf32(i,j,k)
                ap33(i,j,k) = lmf33(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil21(i,j,k) = pp1(i,j,k)
                fil22(i,j,k) = pp2(i,j,k)
                fil23(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do


    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = s31f(i,j,k)
                pc2(i,j,k)  = s32f(i,j,k)
                pc3(i,j,k)  = s33f(i,j,k)
                ap11(i,j,k) = lmf11(i,j,k)
                ap12(i,j,k) = lmf12(i,j,k)
                ap13(i,j,k) = lmf13(i,j,k)
                ap21(i,j,k) = lmf21(i,j,k)
                ap22(i,j,k) = lmf22(i,j,k)
                ap23(i,j,k) = lmf23(i,j,k)
                ap31(i,j,k) = lmf31(i,j,k)
                ap32(i,j,k) = lmf32(i,j,k)
                ap33(i,j,k) = lmf33(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil31(i,j,k) = pp1(i,j,k)
                fil32(i,j,k) = pp2(i,j,k)
                fil33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                smodf(i,j,k)=sqrt(2.*fil11(i,j,k)*fil11(i,j,k)+ &
                    2.*fil12(i,j,k)*fil12(i,j,k)+ &
                    2.*fil13(i,j,k)*fil13(i,j,k)+ &
                    2.*fil21(i,j,k)*fil21(i,j,k)+ &
                    2.*fil22(i,j,k)*fil22(i,j,k)+ &
                    2.*fil23(i,j,k)*fil23(i,j,k)+ &
                    2.*fil31(i,j,k)*fil31(i,j,k)+ &
                    2.*fil32(i,j,k)*fil32(i,j,k)+ &
                    2.*fil33(i,j,k)*fil33(i,j,k))
            enddo
        enddo
    enddo

    if (debugg.eq.1) then
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    write(80+myid+ttt,*)i,j,k,smodf(i,j,k)
                enddo
            enddo
        enddo
    end if

    ! compute the first part of Mik  ( Df^2 * |S|f * Sikf )

    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx

                m11(i,j,k)=alfa3*smodf(i,j,k)*s11f(i,j,k)
                m12(i,j,k)=alfa3*smodf(i,j,k)*s12f(i,j,k)
                m13(i,j,k)=alfa3*smodf(i,j,k)*s13f(i,j,k)
                m21(i,j,k)=alfa3*smodf(i,j,k)*s21f(i,j,k)
                m22(i,j,k)=alfa3*smodf(i,j,k)*s22f(i,j,k)
                m23(i,j,k)=alfa3*smodf(i,j,k)*s23f(i,j,k)
                m31(i,j,k)=alfa3*smodf(i,j,k)*s31f(i,j,k)
                m32(i,j,k)=alfa3*smodf(i,j,k)*s32f(i,j,k)
                m33(i,j,k)=alfa3*smodf(i,j,k)*s33f(i,j,k)
                mrho11(i,j,k)=alfa3*smodf(i,j,k)*rho11f(i,j,k)
                mrho22(i,j,k)=alfa3*smodf(i,j,k)*rho22f(i,j,k)
                mrho33(i,j,k)=alfa3*smodf(i,j,k)*rho33f(i,j,k)

            enddo
        enddo
    enddo
    !
    ! now compute the total Mik (multiplied by -1)
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx

                delt=4.*giac(i,j,k)**(2./3.)
                am11(i,j,k)=-(m11(i,j,k)-smods11f(i,j,k))*delt
                am12(i,j,k)=-(m12(i,j,k)-smods12f(i,j,k))*delt
                am13(i,j,k)=-(m13(i,j,k)-smods13f(i,j,k))*delt
                am21(i,j,k)=-(m21(i,j,k)-smods21f(i,j,k))*delt
                am22(i,j,k)=-(m22(i,j,k)-smods22f(i,j,k))*delt
                am23(i,j,k)=-(m23(i,j,k)-smods23f(i,j,k))*delt
                am31(i,j,k)=-(m31(i,j,k)-smods31f(i,j,k))*delt
                am32(i,j,k)=-(m32(i,j,k)-smods32f(i,j,k))*delt
                am33(i,j,k)=-(m33(i,j,k)-smods33f(i,j,k))*delt
                amrho11(i,j,k)=-(mrho11(i,j,k)-smodrho11f(i,j,k))*delt
                amrho22(i,j,k)=-(mrho22(i,j,k)-smodrho22f(i,j,k))*delt
                amrho33(i,j,k)=-(mrho33(i,j,k)-smodrho33f(i,j,k))*delt

            enddo
        enddo
    enddo
    !
    ! now I need to compute the cartesian component of Mij
    ! starting from the controvariant term Mik
    !
    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = am11(i,j,k)
                pc2(i,j,k)  = am12(i,j,k)
                pc3(i,j,k)  = am13(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil11(i,j,k) = pp1(i,j,k)
                fil12(i,j,k) = pp2(i,j,k)
                fil13(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do


    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = am21(i,j,k)
                pc2(i,j,k)  = am22(i,j,k)
                pc3(i,j,k)  = am23(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil21(i,j,k) = pp1(i,j,k)
                fil22(i,j,k) = pp2(i,j,k)
                fil23(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = am31(i,j,k)
                pc2(i,j,k)  = am32(i,j,k)
                pc3(i,j,k)  = am33(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil31(i,j,k) = pp1(i,j,k)
                fil32(i,j,k) = pp2(i,j,k)
                fil33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do


    !     check
    if (debugg.eq.1) then
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(3320+myid+ttt,*)fil11(i,j,k),fil12(i,j,k),fil13(i,j,k)
                    write(3325+myid+ttt,*)fil21(i,j,k),fil22(i,j,k),fil23(i,j,k)
                    write(3330+myid+ttt,*)fil31(i,j,k),fil32(i,j,k),fil33(i,j,k)
                end do
            end do
        end do
    end if


    !-----------------------------------------------------------------------
    !       Lik
    !-----------------------------------------------------------------------
    !
    ! now I move to the terms of Lik
    ! as a function of ui and UK
    !
    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                pp0(i,j,k) = rho(i,j,k)
            end do
        end do
    end do

    call filter01np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                m33(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                pp0(i,j,k) =uco(i,j,k)
            end do
        end do
    end do

    call filter01np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                m21(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                pp0(i,j,k) = vco(i,j,k)
            end do
        end do
    end do

    call filter01np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                m22(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                pp0(i,j,k) = wco(i,j,k)
            end do
        end do
    end do

    call filter01np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                m23(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do
    !
    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = u(i,j,k)
                p0b(i,j,k) = uco(i,j,k)
            end do
        end do
    end do

    !     (uU)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l11(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = u(i,j,k)
                p0b(i,j,k) = vco(i,j,k)
            end do
        end do
    end do

    !     (uV)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l12(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do


    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = u(i,j,k)
                p0b(i,j,k) = wco(i,j,k)
            end do
        end do
    end do

    !     (uW)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l13(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = v(i,j,k)
                p0b(i,j,k) = uco(i,j,k)
            end do
        end do
    end do

    !     (vU)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l21(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = v(i,j,k)
                p0b(i,j,k) = vco(i,j,k)
            end do
        end do
    end do

    !     (vV)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l22(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do


    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = v(i,j,k)
                p0b(i,j,k) = wco(i,j,k)
            end do
        end do
    end do

    !     (vW)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l23(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = w(i,j,k)
                p0b(i,j,k) = uco(i,j,k)
            end do
        end do
    end do

    !     (wU)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l31(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do


    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = w(i,j,k)
                p0b(i,j,k) = vco(i,j,k)
            end do
        end do
    end do

    !     (wV)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l32(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = w(i,j,k)
                p0b(i,j,k) = wco(i,j,k)
            end do
        end do
    end do

    !     (wW)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                l33(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = rho(i,j,k)
                p0b(i,j,k) = uco(i,j,k)
            end do
        end do
    end do

    !     (rhoU)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                lrho11(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = rho(i,j,k)
                p0b(i,j,k) = vco(i,j,k)
            end do
        end do
    end do

    !     (rhoV)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                lrho22(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do


    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                p0a(i,j,k) = rho(i,j,k)
                p0b(i,j,k) = wco(i,j,k)
            end do
        end do
    end do


    !     (rhoW)f
    call filter02np(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                lrho33(i,j,k) = pp2(i,j,k)
            end do
        end do
    end do

    if (debugg.eq.1) then
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(5000+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                    write(5005+myid+ttt,*)m33(i,j,k)
                    write(5010+myid+ttt,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                    write(5015+myid+ttt,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                    write(5020+myid+ttt,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                    write(5025+myid+ttt,*)lrho11(i,j,k),lrho22(i,j,k),lrho33(i,j,k)
                end do
            end do
        end do
    end if

    !-----------------------------------------------------------------------
    ! now I need to transfer the "ghost" planes to include periodicity
    ! to complete the filtering in zita (k+1 and k-1)
    !
    ! first kparasta
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    call buffer1g(m33,1,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m21,2,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m22,3,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(m23,4,kparasta,myid,nproc,kparasta,kparaend)

    call buffer1g(l11,5,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(l12,6,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(l13,7,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(l21,8,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(l22,9,kparasta,myid,nproc,kparasta,kparaend)
    call buffer1g(l23,10,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l31,11,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l32,12,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l33,13,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho11,14,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho22,15,kparasta,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho33,16,kparasta,myid,nproc, &
        kparasta,kparaend)

    if (kp.eq.0) then

        call MPI_SENDRECV(sbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,tagls, &
            rbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrr, &
            MPI_COMM_WORLD,status,ierr)

    else if (kp.eq.1) then

        if(leftpem /= MPI_PROC_NULL) then
            call MPI_SEND(sbuff1(1),16*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,tagls, &
                MPI_COMM_WORLD,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),16*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrr, &
                MPI_COMM_WORLD,status,ierr)
        endif

        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif


    endif

    call buffer2(rbuff1,m33,1,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m21,2,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m22,3,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m23,4,kparaend+1,myid,nproc, &
        kparasta,kparaend)

    call buffer2(rbuff1,l11,5,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l12,6,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l13,7,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l21,8,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l22,9,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l23,10,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l31,11,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l32,12,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l33,13,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho11,14,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho22,15,kparaend+1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho33,16,kparaend+1,myid,nproc, &
        kparasta,kparaend)

    !
    ! now kparaend
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    call buffer1g(m33,1,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m21,2,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m22,3,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(m23,4,kparaend,myid,nproc,kparasta,kparaend)

    call buffer1g(l11,5,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(l12,6,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(l13,7,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(l21,8,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(l22,9,kparaend,myid,nproc,kparasta,kparaend)
    call buffer1g(l23,10,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l31,11,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l32,12,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(l33,13,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho11,14,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho22,15,kparaend,myid,nproc, &
        kparasta,kparaend)
    call buffer1g(lrho33,16,kparaend,myid,nproc, &
        kparasta,kparaend)

    if (kp.eq.0) then

        call MPI_SENDRECV(sbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrs, &
            rbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)

    else if (kp.eq.1) then

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SEND(sbuff1(1),16*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrs, &
                MPI_COMM_WORLD,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),16*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,taglr, &
                MPI_COMM_WORLD,status,ierr)
        endif

        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif


    endif

    call buffer2(rbuff1,m33,1,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m21,2,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m22,3,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,m23,4,kparasta-1,myid,nproc, &
        kparasta,kparaend)

    call buffer2(rbuff1,l11,5,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l12,6,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l13,7,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l21,8,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l22,9,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l23,10,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l31,11,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l32,12,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,l33,13,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho11,14,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho22,15,kparasta-1,myid,nproc, &
        kparasta,kparaend)
    call buffer2(rbuff1,lrho33,16,kparasta-1,myid,nproc, &
        kparasta,kparaend)

    ! now I can complete the filtering in zita
    !
    call filter03(m33,rhof,kparasta,kparaend,myid,nproc)
    !c      call filter03(m11,uf)
    !c      call filter03(m12,vf)
    !c      call filter03(m13,wf)
    call filter03(m21,ucof,kparasta,kparaend,myid,nproc)
    call filter03(m22,vcof,kparasta,kparaend,myid,nproc)
    call filter03(m23,wcof,kparasta,kparaend,myid,nproc)

    call filter03(l11,uucof,kparasta,kparaend,myid,nproc)
    call filter03(l12,uvcof,kparasta,kparaend,myid,nproc)
    call filter03(l13,uwcof,kparasta,kparaend,myid,nproc)
    call filter03(l21,vucof,kparasta,kparaend,myid,nproc)
    call filter03(l22,vvcof,kparasta,kparaend,myid,nproc)
    call filter03(l23,vwcof,kparasta,kparaend,myid,nproc)
    call filter03(l31,wucof,kparasta,kparaend,myid,nproc)
    call filter03(l32,wvcof,kparasta,kparaend,myid,nproc)
    call filter03(l33,wwcof,kparasta,kparaend,myid,nproc)
    call filter03(lrho11,rhoucof,kparasta,kparaend,myid,nproc)
    call filter03(lrho22,rhovcof,kparasta,kparaend,myid,nproc)
    call filter03(lrho33,rhowcof,kparasta,kparaend,myid,nproc)

    if (debugg.eq.1) then
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(3350+myid+ttt,*)uucof(i,j,k),uvcof(i,j,k),uwcof(i,j,k)
                    write(3355+myid+ttt,*)vucof(i,j,k),vvcof(i,j,k),vwcof(i,j,k)
                    write(3360+myid+ttt,*)wucof(i,j,k),wvcof(i,j,k),wwcof(i,j,k)
                    write(3365+myid+ttt,*)rhoucof(i,j,k),rhoucof(i,j,k), &
                        rhowcof(i,j,k)
                    write(3370+myid+ttt,*)lrho11(i,j,k),lrho22(i,j,k),lrho33(i,j,k)
                    write(3375+myid+ttt,*)rhoucof(i,j,k),rhovcof(i,j,k), &
                        rhowcof(i,j,k)
                    write(3380+myid+ttt,*)rhof(i,j,k)
                    write(3385+myid+ttt,*)ucof(i,j,k),vcof(i,j,k),wcof(i,j,k)
                end do
            end do
        end do
    end if
    !
    ! as in the original code: I derive uf vf and wf by inversion
    ! of filtered controvariant components
    !
    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = ucof(i,j,k)
                pc2(i,j,k)  = vcof(i,j,k)
                pc3(i,j,k)  = wcof(i,j,k)
                ap11(i,j,k) = lmf11(i,j,k)
                ap12(i,j,k) = lmf12(i,j,k)
                ap13(i,j,k) = lmf13(i,j,k)
                ap21(i,j,k) = lmf21(i,j,k)
                ap22(i,j,k) = lmf22(i,j,k)
                ap23(i,j,k) = lmf23(i,j,k)
                ap31(i,j,k) = lmf31(i,j,k)
                ap32(i,j,k) = lmf32(i,j,k)
                ap33(i,j,k) = lmf33(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                uf(i,j,k) = pp1(i,j,k)
                vf(i,j,k) = pp2(i,j,k)
                wf(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    !
    ! now I compute Lik
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx

                l11(i,j,k)=uucof(i,j,k)-uf(i,j,k)*ucof(i,j,k)     !(uU)f-ufUf
                l12(i,j,k)=vucof(i,j,k)-vf(i,j,k)*ucof(i,j,k)     !(vU)f-vfUf
                l13(i,j,k)=wucof(i,j,k)-wf(i,j,k)*ucof(i,j,k)     !(wU)f-wfUf
                l21(i,j,k)=uvcof(i,j,k)-uf(i,j,k)*vcof(i,j,k)     !(uV)f-ufVf
                l22(i,j,k)=vvcof(i,j,k)-vf(i,j,k)*vcof(i,j,k)     !(vV)f-vfVf
                l23(i,j,k)=wvcof(i,j,k)-wf(i,j,k)*vcof(i,j,k)     !(wV)f-wfVf
                l31(i,j,k)=uwcof(i,j,k)-uf(i,j,k)*wcof(i,j,k)     !(uW)f-ufWf
                l32(i,j,k)=vwcof(i,j,k)-vf(i,j,k)*wcof(i,j,k)     !(vW)f-vfWf
                l33(i,j,k)=wwcof(i,j,k)-wf(i,j,k)*wcof(i,j,k)     !(wW)f-wfWf
                lrho11(i,j,k)=rhoucof(i,j,k)-rhof(i,j,k)*ucof(i,j,k)  !(rhoU)f-rhofUf
                lrho22(i,j,k)=rhovcof(i,j,k)-rhof(i,j,k)*vcof(i,j,k)  !(rhoV)f-rhofVf
                lrho33(i,j,k)=rhowcof(i,j,k)-rhof(i,j,k)*wcof(i,j,k)  !(rhoW)f-rhofWf

            enddo
        enddo
    enddo
    !
    ! now I compute the cartesian product MijMij and LijMij
    ! from the controvariant term Mik and Lik
    ! (pay attention: in the original code Lik are computed with a
    ! transpose way, therefore L12=(vU)f-vfUf  instead of  (uV)f-ufVf
    ! in the computation of cartesian component I consider this aspect,
    ! in fact the 3 terms are transpose)
    !

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = l11(i,j,k)
                pc2(i,j,k)  = l21(i,j,k)
                pc3(i,j,k)  = l31(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fipp11(i,j,k) = pp1(i,j,k)
                fipp12(i,j,k) = pp2(i,j,k)
                fipp13(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = l12(i,j,k)
                pc2(i,j,k)  = l22(i,j,k)
                pc3(i,j,k)  = l32(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fipp21(i,j,k) = pp1(i,j,k)
                fipp22(i,j,k) = pp2(i,j,k)
                fipp23(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = l13(i,j,k)
                pc2(i,j,k)  = l23(i,j,k)
                pc3(i,j,k)  = l33(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fipp31(i,j,k) = pp1(i,j,k)
                fipp32(i,j,k) = pp2(i,j,k)
                fipp33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                den(i,j,k)=fil11(i,j,k)*fil11(i,j,k)+fil12(i,j,k)*fil12(i,j,k)+ &
                    fil13(i,j,k)*fil13(i,j,k)+fil21(i,j,k)*fil21(i,j,k)+ &
                    fil22(i,j,k)*fil22(i,j,k)+fil23(i,j,k)*fil23(i,j,k)+ &
                    fil31(i,j,k)*fil31(i,j,k)+fil32(i,j,k)*fil32(i,j,k)+ &
                    fil33(i,j,k)*fil33(i,j,k)
                num(i,j,k)=fil11(i,j,k)*fipp11(i,j,k)+ &
                    fil12(i,j,k)*fipp12(i,j,k)+ &
                    fil13(i,j,k)*fipp13(i,j,k)+ &
                    fil21(i,j,k)*fipp21(i,j,k)+ &
                    fil22(i,j,k)*fipp22(i,j,k)+ &
                    fil23(i,j,k)*fipp23(i,j,k)+ &
                    fil31(i,j,k)*fipp31(i,j,k)+ &
                    fil32(i,j,k)*fipp32(i,j,k)+ &
                    fil33(i,j,k)*fipp33(i,j,k)
            enddo
        enddo
    enddo
    !
    ! same thing for rho terms
    ! I find the physical component rhoMij and rhoLij
    ! from the contravariant ones to compute rhoMijrhoMij and rhoMijrhoLij
    !

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = amrho11(i,j,k)
                pc2(i,j,k)  = amrho22(i,j,k)
                pc3(i,j,k)  = amrho33(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil11(i,j,k) = pp1(i,j,k)
                fil22(i,j,k) = pp2(i,j,k)
                fil33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = lrho11(i,j,k)
                pc2(i,j,k)  = lrho22(i,j,k)
                pc3(i,j,k)  = lrho33(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fipp11(i,j,k) = pp1(i,j,k)
                fipp22(i,j,k) = pp2(i,j,k)
                fipp33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do


    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                denrho(i,j,k)=fil11(i,j,k)*fil11(i,j,k)+ &
                    fil22(i,j,k)*fil22(i,j,k)+ &
                    fil33(i,j,k)*fil33(i,j,k)
                numrho(i,j,k)=fil11(i,j,k)*fipp11(i,j,k)+ &
                    fil22(i,j,k)*fipp22(i,j,k)+ &
                    fil33(i,j,k)*fipp33(i,j,k)
            enddo
        enddo
    enddo
    !
    !-----------------------------------------------------------------------
    !               SCALE SIMILAR PART
    !-----------------------------------------------------------------------
    !
    ! second term of numerator: Nij defined as the difference between two
    ! terms Bij and Aij:
    !
    ! Bij = (uiF*ujF)fF - uiFfF*ujFfF
    !
    ! Aij = (ui*uj)fF - (uif*ujf)F
    !
    ! where f --> filter bar defined in filter04, filter05 and filter06
    !       F --> filter hat defined previously in (f01, f02, f03)
    !
    ! as in the previous case, the terms are computed as product of
    ! cartesian and controvariant so all ui are therefore UK and through
    ! inversion I obtain the filtered cartesian term --> Aik and Bik
    !
    ! a similar approach can be used for density
    !
    !-----------------------------------------------------------------------
    !
    ! Aik term (not computed at j=1 and j=jy)
    ! bar filtering for metric

    if(inmod)then

        call filter04g(uco,m21,kparasta,kparaend,myid,nproc)
        call filter04g(vco,m22,kparasta,kparaend,myid,nproc)
        call filter04g(wco,m23,kparasta,kparaend,myid,nproc)

        if (debugg.eq.1) then
            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        write(415+myid+ttt,*)m11(i,j,k),m12(i,j,k),m13(i,j,k)
                        write(420+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                        write(425+myid+ttt,*)m31(i,j,k),m32(i,j,k),m33(i,j,k)
                    end do
                end do
            end do
        end if

        !
        ! I need to treat a 3D matrix, the density
        !
        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = rho(i,j,k)
                end do
            end do
        end do

        call filter04g(p0a,m33,kparasta,kparaend,myid,nproc)

        if (debugg.eq.1) then

            do k=kparasta-1,kparaend+1
                do j=0,n2+1
                    do i=0,n1+1
                        write(540+myid+ttt,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                        write(545+myid+ttt,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                        write(550+myid+ttt,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                        write(555+myid+ttt,*)v(i,j,k), &
                            uco(i,j,k),vco(i,j,k),wco(i,j,k)
                    end do
                end do
            end do
        end if


        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = u(i,j,k)
                    p0b(i,j,k) = uco(i,j,k)
                end do
            end do
        end do

        !     (uU)f
        call filter05g(p0a,p0b,l11,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = u(i,j,k)
                    p0b(i,j,k) = vco(i,j,k)
                end do
            end do
        end do

        !     (uV)f
        call filter05g(p0a,p0b,l12,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = u(i,j,k)
                    p0b(i,j,k) = wco(i,j,k)
                end do
            end do
        end do

        !     (uW)f
        call filter05g(p0a,p0b,l13,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = v(i,j,k)
                    p0b(i,j,k) = uco(i,j,k)
                end do
            end do
        end do

        !     (vU)f
        call filter05g(p0a,p0b,l21,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = v(i,j,k)
                    p0b(i,j,k) = vco(i,j,k)
                end do
            end do
        end do

        !     (vV)f
        call filter05g(p0a,p0b,l22,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = v(i,j,k)
                    p0b(i,j,k) = wco(i,j,k)
                end do
            end do
        end do

        !     (vW)f
        call filter05g(p0a,p0b,l23,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = w(i,j,k)
                    p0b(i,j,k) = uco(i,j,k)
                end do
            end do
        end do

        !     (wU)f
        call filter05g(p0a,p0b,l31,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = w(i,j,k)
                    p0b(i,j,k) = vco(i,j,k)
                end do
            end do
        end do

        !     (wV)f
        call filter05g(p0a,p0b,l32,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = w(i,j,k)
                    p0b(i,j,k) = wco(i,j,k)
                end do
            end do
        end do

        !     (wW)f
        call filter05g(p0a,p0b,l33,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = rho(i,j,k)
                    p0b(i,j,k) = uco(i,j,k)
                end do
            end do
        end do

        !     (rhoU)f
        call filter05g(p0a,p0b,lrho11,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = rho(i,j,k)
                    p0b(i,j,k) = vco(i,j,k)
                end do
            end do
        end do

        !     (rhoV)f
        call filter05g(p0a,p0b,lrho22,kparasta,kparaend,myid,nproc)

        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = rho(i,j,k)
                    p0b(i,j,k) = wco(i,j,k)
                end do
            end do
        end do

        !     (rhoW)f
        call filter05g(p0a,p0b,lrho33,kparasta,kparaend,myid,nproc)

        !
        call filter04g(apcsx,m11m,kparasta,kparaend,myid,nproc)
        call filter04g(apcsy,m12m,kparasta,kparaend,myid,nproc)
        call filter04g(apcsz,m13m,kparasta,kparaend,myid,nproc)
        call filter04g(apetx,m21m,kparasta,kparaend,myid,nproc)
        call filter04g(apety,m22m,kparasta,kparaend,myid,nproc)
        call filter04g(apetz,m23m,kparasta,kparaend,myid,nproc)
        call filter04g(apztx,m31m,kparasta,kparaend,myid,nproc)
        call filter04g(apzty,m32m,kparasta,kparaend,myid,nproc)
        call filter04g(apztz,m33m,kparasta,kparaend,myid,nproc)

        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(1500+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k),m33(i,j,k)
                        write(1505+myid+ttt,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)

                        write(1510+myid+ttt,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)

                        write(1515+myid+ttt,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                        write(1520+myid+ttt,*)lrho11(i,j,k),lrho22(i,j,k),lrho33(i,j,k)
                        write(1525+myid+ttt,*)m11m(i,j,k),m12m(i,j,k),m13m(i,j,k)
                        write(1530+myid+ttt,*)m21m(i,j,k),m22m(i,j,k),m23m(i,j,k)
                        write(1535+myid+ttt,*)m31m(i,j,k),m32m(i,j,k),m33m(i,j,k)
                    end do
                end do
            end do
        end if

        ! now i need to transfer the "ghost" plane, to impose periodicity
        ! and complete the filtering in zita (k+1 and k-1)
        !
        ! first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m21,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m33,4,kparasta,myid,nproc,kparasta,kparaend)
        !
        call buffer1g(l11,5,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l12,6,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l13,7,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l21,8,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l22,9,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l23,10,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l31,11,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l32,12,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l33,13,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(lrho11,14,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho22,15,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho33,16,kparasta,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer1g(m11m,17,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,18,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,19,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,20,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,21,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,22,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,23,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,24,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,25,kparasta,myid,nproc, &
            kparasta,kparaend)
        !
        ! if I put pe and not pem I have implicitly periodicity in k
        ! and values in k=0 and k=jz+1 are defined in filter06
        !
        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),25*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),25*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),25*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),25*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif
        !
        call buffer2(rbuff1,m21,1,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,2,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,3,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33,4,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer2(rbuff1,l11,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l12,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l13,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l21,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l22,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l23,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l31,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l32,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l33,13,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho11,14,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho22,15,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho33,16,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer2(rbuff1,m11m,17,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,18,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,19,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,20,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,21,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,22,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,23,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,24,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,25,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo
        !
        call buffer1g(m21,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m33,4,kparaend,myid,nproc,kparasta,kparaend)
        !
        call buffer1g(l11,5,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l12,6,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l13,7,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l21,8,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l22,9,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l23,10,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(l31,11,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(l32,12,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(l33,13,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho11,14,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho22,15,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho33,16,kparaend,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer1g(m11m,17,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,18,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,19,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,20,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,21,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,22,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,23,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,24,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,25,kparaend,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),25*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),25*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),25*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),25*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m21,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer2(rbuff1,l11,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l12,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l13,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l21,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l22,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l23,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l31,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l32,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l33,13,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho11,14,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho22,15,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho33,16,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer2(rbuff1,m11m,17,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,18,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,19,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,20,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,21,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,22,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,23,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,24,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,25,kparasta-1,myid,nproc, &
            kparasta,kparaend)

        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(4000+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k), &
                            m33(i,j,k)
                        write(4005+myid+ttt,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                        write(4010+myid+ttt,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                        write(4015+myid+ttt,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                        write(4020+myid+ttt,*)lrho11(i,j,k),lrho22(i,j,k), &
                            lrho33(i,j,k)
                        write(4025+myid+ttt,*)m11m(i,j,k),m12m(i,j,k),m13m(i,j,k)
                        write(4030+myid+ttt,*)m21m(i,j,k),m22m(i,j,k),m23m(i,j,k)
                        write(4035+myid+ttt,*)m31m(i,j,k),m32m(i,j,k),m33m(i,j,k)
                    end do
                end do
            end do
        end if
        !
        ! complete the filtering in zita
        !
        call filter06(m21,ucof,kparasta,kparaend,myid,nproc)
        call filter06(m22,vcof,kparasta,kparaend,myid,nproc)
        call filter06(m23,wcof,kparasta,kparaend,myid,nproc)
        call filter06(m33,rhof,kparasta,kparaend,myid,nproc)
        !
        call filter06(l11,uucof,kparasta,kparaend,myid,nproc)
        call filter06(l12,uvcof,kparasta,kparaend,myid,nproc)
        call filter06(l13,uwcof,kparasta,kparaend,myid,nproc)
        call filter06(l21,vucof,kparasta,kparaend,myid,nproc)
        call filter06(l22,vvcof,kparasta,kparaend,myid,nproc)
        call filter06(l23,vwcof,kparasta,kparaend,myid,nproc)
        call filter06(l31,wucof,kparasta,kparaend,myid,nproc)
        call filter06(l32,wvcof,kparasta,kparaend,myid,nproc)
        call filter06(l33,wwcof,kparasta,kparaend,myid,nproc)
        call filter06(lrho11,rhoucof,kparasta,kparaend,myid,nproc)
        call filter06(lrho22,rhovcof,kparasta,kparaend,myid,nproc)
        call filter06(lrho33,rhowcof,kparasta,kparaend,myid,nproc)
        !
        call filter06(m11m,lmf11,kparasta,kparaend,myid,nproc)
        call filter06(m12m,lmf12,kparasta,kparaend,myid,nproc)
        call filter06(m13m,lmf13,kparasta,kparaend,myid,nproc)
        call filter06(m21m,lmf21,kparasta,kparaend,myid,nproc)
        call filter06(m22m,lmf22,kparasta,kparaend,myid,nproc)
        call filter06(m23m,lmf23,kparasta,kparaend,myid,nproc)
        call filter06(m31m,lmf31,kparasta,kparaend,myid,nproc)
        call filter06(m32m,lmf32,kparasta,kparaend,myid,nproc)
        call filter06(m33m,lmf33,kparasta,kparaend,myid,nproc)

        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(4060+myid+ttt,*)ucof(i,j,k),vcof(i,j,k),wcof(i,j,k), &
                            rhof(i,j,k)
                        write(4065+myid+ttt,*)uucof(i,j,k),uvcof(i,j,k),uwcof(i,j,k)
                        write(4070+myid+ttt,*)vucof(i,j,k),vvcof(i,j,k),vwcof(i,j,k)
                        write(4075+myid+ttt,*)wucof(i,j,k),wvcof(i,j,k),wwcof(i,j,k)
                        write(4080+myid+ttt,*)rhoucof(i,j,k),rhovcof(i,j,k), &
                            rhowcof(i,j,k)
                        write(4085+myid+ttt,*)lmf11(i,j,k),lmf12(i,j,k),lmf13(i,j,k)
                        write(4090+myid+ttt,*)lmf21(i,j,k),lmf22(i,j,k),lmf23(i,j,k)
                        write(4095+myid+ttt,*)lmf31(i,j,k),lmf32(i,j,k),lmf33(i,j,k)

                    end do
                end do
            end do
        end if
        !
        ! cartesian component derived from the controvariant ones
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = ucof(i,j,k)
                    pc2(i,j,k)  = vcof(i,j,k)
                    pc3(i,j,k)  = wcof(i,j,k)
                    ap11(i,j,k) = lmf11(i,j,k)
                    ap12(i,j,k) = lmf12(i,j,k)
                    ap13(i,j,k) = lmf13(i,j,k)
                    ap21(i,j,k) = lmf21(i,j,k)
                    ap22(i,j,k) = lmf22(i,j,k)
                    ap23(i,j,k) = lmf23(i,j,k)
                    ap31(i,j,k) = lmf31(i,j,k)
                    ap32(i,j,k) = lmf32(i,j,k)
                    ap33(i,j,k) = lmf33(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    uf(i,j,k) = pp1(i,j,k)
                    vf(i,j,k) = pp2(i,j,k)
                    wf(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        !
        ! Leonard term Aij
        !
        do k=kparasta,kparaend
            !cc      do j=2,jy-1
            do j=1,jy
                do i=1,jx
                    ass11(i,j,k)=uucof(i,j,k)-uf(i,j,k)*ucof(i,j,k)       !uU
                    ass12(i,j,k)=vucof(i,j,k)-vf(i,j,k)*ucof(i,j,k)       !vU
                    ass13(i,j,k)=wucof(i,j,k)-wf(i,j,k)*ucof(i,j,k)       !wU
                    ass21(i,j,k)=uvcof(i,j,k)-uf(i,j,k)*vcof(i,j,k)       !uV
                    ass22(i,j,k)=vvcof(i,j,k)-vf(i,j,k)*vcof(i,j,k)       !vV
                    ass23(i,j,k)=wvcof(i,j,k)-wf(i,j,k)*vcof(i,j,k)       !wV
                    ass31(i,j,k)=uwcof(i,j,k)-uf(i,j,k)*wcof(i,j,k)       !uW
                    ass32(i,j,k)=vwcof(i,j,k)-vf(i,j,k)*wcof(i,j,k)       !vW
                    ass33(i,j,k)=wwcof(i,j,k)-wf(i,j,k)*wcof(i,j,k)       !wW
                    assrho11(i,j,k)=rhoucof(i,j,k)-rhof(i,j,k)*ucof(i,j,k)   !rhoU
                    assrho22(i,j,k)=rhovcof(i,j,k)-rhof(i,j,k)*vcof(i,j,k)   !rhoV
                    assrho33(i,j,k)=rhowcof(i,j,k)-rhof(i,j,k)*wcof(i,j,k)   !rhoW
                enddo
            enddo
        enddo


        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        write(8875+myid+ttt,*)ass11(i,j,k),ass12(i,j,k),ass13(i,j,k)
                        write(8880+myid+ttt,*)ass21(i,j,k),ass22(i,j,k),ass23(i,j,k)
                        write(8885+myid+ttt,*)ass31(i,j,k),ass32(i,j,k),ass33(i,j,k)
                        write(8890+myid+ttt,*)assrho11(i,j,k),assrho22(i,j,k), &
                            assrho33(i,j,k)
                    enddo
                enddo
            enddo
        end if

        ! periodicity for Aik and Arhoik is done in the subroutine
        ! only in csi and eta, while on zita is done one time at the
        ! end
        !
        call periodic(ass11,ass12,ass13)
        call periodic(ass21,ass22,ass23)
        call periodic(ass31,ass32,ass33)
        call periodic(assrho11,assrho22,assrho33)


        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(4040+myid+ttt,*)ass11(i,j,k),ass12(i,j,k),ass13(i,j,k)
                        write(4045+myid+ttt,*)ass21(i,j,k),ass22(i,j,k),ass23(i,j,k)
                        write(4050+myid+ttt,*)ass31(i,j,k),ass32(i,j,k),ass33(i,j,k)
                        write(4055+myid+ttt,*)assrho11(i,j,k),assrho22(i,j,k), &
                            assrho33(i,j,k)
                    end do
                end do
            end do
        end if

        !
        ! periodicity in zita
        !
        if(kp.eq.0)then
            !
            ! here I transfer planes k=0 and k=jz+1
            ! so that P0 knows k=jz and Pn-1 knows k=1
            !
            do m=1,40*(jx+2)*(jy+2)
                sbuff(m)=0.
                rbuff(m)=0.
            enddo

            if (myid.eq.nproc-1) then

                call buffer1g(ass11,1,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass12,2,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass13,3,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass21,4,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass22,5,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass23,6,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass31,7,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass32,8,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ass33,9,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(assrho11,10,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(assrho22,11,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(assrho33,12,jz,myid,nproc,kparasta,kparaend)

            else if (myid.eq.0) then

                call buffer1g(ass11,1,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass12,2,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass13,3,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass21,4,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass22,5,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass23,6,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass31,7,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass32,8,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ass33,9,1,myid,nproc,kparasta,kparaend)
                call buffer1g(assrho11,10,1,myid,nproc,kparasta,kparaend)
                call buffer1g(assrho22,11,1,myid,nproc,kparasta,kparaend)
                call buffer1g(assrho33,12,1,myid,nproc,kparasta,kparaend)

            endif
            !
            ! now the exchange so
            ! P0 knows k=jz
            ! Pn-1 knows k=1

            if (myid.eq.nproc-1) then

                call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,0,901, &
                    rbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,0,801, &
                    MPI_COMM_WORLD,status,ierr)

            else if (myid.eq.0) then

                call MPI_SENDRECV( &
                    sbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,801, &
                    rbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,901, &
                    MPI_COMM_WORLD,status,ierr)

            endif

            if (myid.eq.0) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)

            else if (myid.eq.nproc-1) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)

            endif

            ! now P0   knows the plane k=jz
            ! and Pn-1 knows the plane k=1

            if(myid.eq.0)then
                do j=1,jy
                    do i=1,jx

                        ass11(i,j,0)    = piano1(i,j)
                        ass12(i,j,0)    = piano2(i,j)
                        ass13(i,j,0)    = piano3(i,j)
                        ass21(i,j,0)    = piano4(i,j)
                        ass22(i,j,0)    = piano5(i,j)
                        ass23(i,j,0)    = piano6(i,j)
                        ass31(i,j,0)    = piano7(i,j)
                        ass32(i,j,0)    = piano8(i,j)
                        ass33(i,j,0)    = piano9(i,j)
                        assrho11(i,j,0) = piano10(i,j)
                        assrho22(i,j,0) = piano11(i,j)
                        assrho33(i,j,0) = piano12(i,j)

                    enddo
                enddo
            endif
            !
            if(myid.eq.nproc-1)then
                do j=1,jy
                    do i=1,jx

                        ass11(i,j,jz+1)    = piano1(i,j)
                        ass12(i,j,jz+1)    = piano2(i,j)
                        ass13(i,j,jz+1)    = piano3(i,j)
                        ass21(i,j,jz+1)    = piano4(i,j)
                        ass22(i,j,jz+1)    = piano5(i,j)
                        ass23(i,j,jz+1)    = piano6(i,j)
                        ass31(i,j,jz+1)    = piano7(i,j)
                        ass32(i,j,jz+1)    = piano8(i,j)
                        ass33(i,j,jz+1)    = piano9(i,j)
                        assrho11(i,j,jz+1) = piano10(i,j)
                        assrho22(i,j,jz+1) = piano11(i,j)
                        assrho33(i,j,jz+1) = piano12(i,j)

                    enddo
                enddo
            endif
        !
        endif

        !
        ! now final filtering on Aik
        !
        call filter01(ass11,m11,myid,nproc,kparasta,kparaend)
        call filter01(ass12,m12,myid,nproc,kparasta,kparaend)
        call filter01(ass13,m13,myid,nproc,kparasta,kparaend)
        call filter01(ass21,m21,myid,nproc,kparasta,kparaend)
        call filter01(ass22,m22,myid,nproc,kparasta,kparaend)
        call filter01(ass23,m23,myid,nproc,kparasta,kparaend)
        call filter01(ass31,m31,myid,nproc,kparasta,kparaend)
        call filter01(ass32,m32,myid,nproc,kparasta,kparaend)
        call filter01(ass33,m33,myid,nproc,kparasta,kparaend)
        call filter01(assrho11,mrho11,myid,nproc,kparasta,kparaend)
        call filter01(assrho22,mrho22,myid,nproc,kparasta,kparaend)
        call filter01(assrho33,mrho33,myid,nproc,kparasta,kparaend)
        !
        ! now I need to transfer the "ghost" plane to
        ! include periodicity to complete the filtering in zita (k+1 and k-1)
        !
        ! first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m11,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m12,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m13,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,4,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,5,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,6,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m31,7,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m32,8,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m33,9,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(mrho11,10,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(mrho22,11,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(mrho33,12,kparasta,myid,nproc, &
            kparasta,kparaend)
        !
        ! if I put pe and not pem, implicitly I have periodicity in k
        !
        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m11,1,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12,2,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13,3,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21,4,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,mrho11,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,mrho22,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,mrho33,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m11,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m12,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m13,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,4,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,5,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,6,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m31,7,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m32,8,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m33,9,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(mrho11,10,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(mrho22,11,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(mrho33,12,kparaend,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then


            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m11,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,mrho11,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,mrho22,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,mrho33,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)


        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(4100+myid+ttt,*)m11(i,j,k),m12(i,j,k),m13(i,j,k)
                        write(4105+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                        write(4110+myid+ttt,*)m31(i,j,k),m32(i,j,k),m33(i,j,k)
                        write(4115+myid+ttt,*)mrho11(i,j,k),mrho22(i,j,k),mrho33(i,j,k)
                    end do
                end do
            end do
        end if
        !
        ! periodicity in k (implemented in filter03)
        !
        ! now I complete the filtering in zita
        !

        call filter03(m11,al11,kparasta,kparaend,myid,nproc)
        call filter03(m12,al12,kparasta,kparaend,myid,nproc)
        call filter03(m13,al13,kparasta,kparaend,myid,nproc)
        call filter03(m21,al21,kparasta,kparaend,myid,nproc)
        call filter03(m22,al22,kparasta,kparaend,myid,nproc)
        call filter03(m23,al23,kparasta,kparaend,myid,nproc)
        call filter03(m31,al31,kparasta,kparaend,myid,nproc)
        call filter03(m32,al32,kparasta,kparaend,myid,nproc)
        call filter03(m33,al33,kparasta,kparaend,myid,nproc)
        call filter03(mrho11,alrho11,kparasta,kparaend,myid,nproc)
        call filter03(mrho22,alrho22,kparasta,kparaend,myid,nproc)
        call filter03(mrho33,alrho33,kparasta,kparaend,myid,nproc)

        !
        ! store to compute SGS stress
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    sgs11(i,j,k)=ass11(i,j,k)
                    sgs12(i,j,k)=ass21(i,j,k)
                    sgs13(i,j,k)=ass31(i,j,k)
                    sgs21(i,j,k)=ass12(i,j,k)
                    sgs22(i,j,k)=ass22(i,j,k)
                    sgs23(i,j,k)=ass32(i,j,k)
                    sgs31(i,j,k)=ass13(i,j,k)
                    sgs32(i,j,k)=ass23(i,j,k)
                    sgs33(i,j,k)=ass33(i,j,k)
                    sgsrho11(i,j,k)=assrho11(i,j,k)
                    sgsrho22(i,j,k)=assrho22(i,j,k)
                    sgsrho33(i,j,k)=assrho33(i,j,k)
                enddo
            enddo
        enddo

        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(9000+myid+ttt,*)sgs11(i,j,k),sgs12(i,j,k),sgs13(i,j,k)
                        write(9005+myid+ttt,*)sgs21(i,j,k),sgs22(i,j,k),sgs23(i,j,k)
                        write(9010+myid+ttt,*)sgs31(i,j,k),sgs32(i,j,k),sgs33(i,j,k)
                        write(9015+myid+ttt,*)sgsrho11(i,j,k),sgsrho22(i,j,k), &
                            sgsrho33(i,j,k)
                    end do
                end do
            end do
        end if

        !
        ! compute the term Bik - first part
        !
        ! (observation: defining in a correct way the filtering
        ! to compute Aik and using different matrix like ucofb etc
        ! - already dclared - the original ucof could be used again
        ! to compute Lik !!)
        !
        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = rho(i,j,k)
                end do
            end do
        end do

        call filter01(p0a,m33,myid,nproc,kparasta,kparaend)


        call filter01(uco,m21,myid,nproc,kparasta,kparaend)
        call filter01(vco,m22,myid,nproc,kparasta,kparaend)
        call filter01(wco,m23,myid,nproc,kparasta,kparaend)
        !
        call filter01(apcsx,m11m,myid,nproc,kparasta,kparaend)
        call filter01(apcsy,m12m,myid,nproc,kparasta,kparaend)
        call filter01(apcsz,m13m,myid,nproc,kparasta,kparaend)
        call filter01(apetx,m21m,myid,nproc,kparasta,kparaend)
        call filter01(apety,m22m,myid,nproc,kparasta,kparaend)
        call filter01(apetz,m23m,myid,nproc,kparasta,kparaend)
        call filter01(apztx,m31m,myid,nproc,kparasta,kparaend)
        call filter01(apzty,m32m,myid,nproc,kparasta,kparaend)
        call filter01(apztz,m33m,myid,nproc,kparasta,kparaend)

        !
        ! now I transfer the "ghost" plane to include periodicity
        ! to complete the filtering on metric in zita (k+1 and k-1)
        !
        ! first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m33,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,4,kparasta,myid,nproc,kparasta,kparaend)

        call buffer1g(m11m,5,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,6,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,7,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,8,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,9,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,10,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,11,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,12,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,13,kparasta,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m33,1,kparaend+1,myid,nproc,kparasta,kparaend)
        call buffer2(rbuff1,m21,2,kparaend+1,myid,nproc,kparasta,kparaend)
        call buffer2(rbuff1,m22,3,kparaend+1,myid,nproc,kparasta,kparaend)
        call buffer2(rbuff1,m23,4,kparaend+1,myid,nproc,kparasta,kparaend)

        call buffer2(rbuff1,m11m,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m33,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,4,kparaend,myid,nproc,kparasta,kparaend)

        call buffer1g(m11m,5,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m12m,6,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m13m,7,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m21m,8,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m22m,9,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m23m,10,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m31m,11,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m32m,12,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m33m,13,kparaend,myid,nproc,kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then


            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m33,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)

        call buffer2(rbuff1,m11m,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now complete filtering in zita
        !
        call filter03(m33,rhof,kparasta,kparaend,myid,nproc)
        call filter03(m21,ucof,kparasta,kparaend,myid,nproc)
        call filter03(m22,vcof,kparasta,kparaend,myid,nproc)
        call filter03(m23,wcof,kparasta,kparaend,myid,nproc)

        call filter03(m11m,lmf11,kparasta,kparaend,myid,nproc)
        call filter03(m12m,lmf12,kparasta,kparaend,myid,nproc)
        call filter03(m13m,lmf13,kparasta,kparaend,myid,nproc)
        call filter03(m21m,lmf21,kparasta,kparaend,myid,nproc)
        call filter03(m22m,lmf22,kparasta,kparaend,myid,nproc)
        call filter03(m23m,lmf23,kparasta,kparaend,myid,nproc)
        call filter03(m31m,lmf31,kparasta,kparaend,myid,nproc)
        call filter03(m32m,lmf32,kparasta,kparaend,myid,nproc)
        call filter03(m33m,lmf33,kparasta,kparaend,myid,nproc)
        !
        ! in the original code define values in j=0 and j=jy+1
        ! for the necessary quantities
        !
        do k=kparasta,kparaend
            do i=1,jx
                ucof(i,0,k)=uco(i,0,k)
                vcof(i,0,k)=vco(i,0,k)
                wcof(i,0,k)=wco(i,0,k)
                rhof(i,0,k)=rho(i,0,k)
                lmf11(i,0,k)=apcsx(i,0,k)
                lmf12(i,0,k)=apcsy(i,0,k)
                lmf13(i,0,k)=apcsz(i,0,k)
                lmf21(i,0,k)=apetx(i,0,k)
                lmf22(i,0,k)=apety(i,0,k)
                lmf23(i,0,k)=apetz(i,0,k)
                lmf31(i,0,k)=apztx(i,0,k)
                lmf32(i,0,k)=apzty(i,0,k)
                lmf33(i,0,k)=apztz(i,0,k)
                !
                ucof(i,jy+1,k)=uco(i,jy+1,k)
                vcof(i,jy+1,k)=vco(i,jy+1,k)
                wcof(i,jy+1,k)=wco(i,jy+1,k)
                rhof(i,jy+1,k)=rho(i,jy+1,k)
                lmf11(i,jy+1,k)=apcsx(i,jy+1,k)
                lmf12(i,jy+1,k)=apcsy(i,jy+1,k)
                lmf13(i,jy+1,k)=apcsz(i,jy+1,k)
                lmf21(i,jy+1,k)=apetx(i,jy+1,k)
                lmf22(i,jy+1,k)=apety(i,jy+1,k)
                lmf23(i,jy+1,k)=apetz(i,jy+1,k)
                lmf31(i,jy+1,k)=apztx(i,jy+1,k)
                lmf32(i,jy+1,k)=apzty(i,jy+1,k)
                lmf33(i,jy+1,k)=apztz(i,jy+1,k)
            enddo
        enddo




        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5030+myid+ttt,*)ucof(i,j,k),vcof(i,j,k),wcof(i,j,k)
                        write(5035+myid+ttt,*)rhof(i,j,k)
                        write(5040+myid+ttt,*)lmf11(i,j,k),lmf12(i,j,k),lmf13(i,j,k)
                        write(5045+myid+ttt,*)lmf21(i,j,k),lmf22(i,j,k),lmf23(i,j,k)
                        write(5050+myid+ttt,*)lmf31(i,j,k),lmf32(i,j,k),lmf33(i,j,k)
                    end do
                end do
            end do
        end if

        !
        ! periodicity (necessary for the next filtering) is done only
        ! in csi and eta in the subroutine, while in zita is done once at the end
        !
        call periodic(ucof,vcof,wcof)
        call periodic(rhof,rhof,rhof)
        call periodic(lmf11,lmf12,lmf13)
        call periodic(lmf21,lmf22,lmf23)
        call periodic(lmf31,lmf32,lmf33)
        !
        ! periodicity in zita
        !
        if(kp.eq.0)then
            !
            ! here I need to transfer the planes k=0 and k=jz+1
            ! so taht P0 knows k=jz and Pn-1 knows k=1
            !
            do m=1,40*(jx+2)*(jy+2)
                sbuff(m)=0.
                rbuff(m)=0.
            enddo

            if (myid.eq.nproc-1) then

                call buffer1g(ucof,1,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vcof,2,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wcof,3,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhof,4,jz,myid,nproc,kparasta,kparaend)

                call buffer1g(lmf11,5,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf12,6,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf13,7,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf21,8,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf22,9,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf23,10,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf31,11,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf32,12,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf33,13,jz,myid,nproc,kparasta,kparaend)

            else if (myid.eq.0) then

                call buffer1g(ucof,1,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vcof,2,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wcof,3,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhof,4,1,myid,nproc,kparasta,kparaend)

                call buffer1g(lmf11,5,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf12,6,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf13,7,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf21,8,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf22,9,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf23,10,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf31,11,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf32,12,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmf33,13,1,myid,nproc,kparasta,kparaend)

            endif
            !
            ! now the exchange so
            ! P0 knows k=jz
            ! Pn-1 knows k=1
            !
            if (myid.eq.nproc-1) then

                call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,0,91, &
                    rbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,0,81, &
                    MPI_COMM_WORLD,status,ierr)

            else if (myid.eq.0) then

                call MPI_SENDRECV( &
                    sbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,81, &
                    rbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,91, &
                    MPI_COMM_WORLD,status,ierr)

            endif

            if (myid.eq.0) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano13,13,myid,nproc,kparasta,kparaend)


            else if (myid.eq.nproc-1) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano13,13,myid,nproc,kparasta,kparaend)

            endif

            ! now P0 knows the plane k=jz
            ! and Pn-1 knows the plane k=1
            !

            if (myid.eq.0) then

                do j=1,jy
                    do i=1,jx
                        ucof(i,j,0)  = piano1(i,j)
                        vcof(i,j,0)  = piano2(i,j)
                        wcof(i,j,0)  = piano3(i,j)
                        rhof(i,j,0)  = piano4(i,j)
                        lmf11(i,j,0) = piano5(i,j)
                        lmf12(i,j,0) = piano6(i,j)
                        lmf13(i,j,0) = piano7(i,j)
                        lmf21(i,j,0) = piano8(i,j)
                        lmf22(i,j,0) = piano9(i,j)
                        lmf23(i,j,0) = piano10(i,j)
                        lmf31(i,j,0) = piano11(i,j)
                        lmf32(i,j,0) = piano12(i,j)
                        lmf33(i,j,0) = piano13(i,j)
                    end do
                end do

            end if

            if (myid.eq.nproc-1) then

                do j=1,jy
                    do i=1,jx
                        ucof(i,j,jz+1)  = piano1(i,j)
                        vcof(i,j,jz+1)  = piano2(i,j)
                        wcof(i,j,jz+1)  = piano3(i,j)
                        rhof(i,j,jz+1)  = piano4(i,j)
                        lmf11(i,j,jz+1) = piano5(i,j)
                        lmf12(i,j,jz+1) = piano6(i,j)
                        lmf13(i,j,jz+1) = piano7(i,j)
                        lmf21(i,j,jz+1) = piano8(i,j)
                        lmf22(i,j,jz+1) = piano9(i,j)
                        lmf23(i,j,jz+1) = piano10(i,j)
                        lmf31(i,j,jz+1) = piano11(i,j)
                        lmf32(i,j,jz+1) = piano12(i,j)
                        lmf33(i,j,jz+1) = piano13(i,j)
                    end do
                end do

            end if
        !
        endif


        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5060+myid+ttt,*)ucof(i,j,k),vcof(i,j,k),wcof(i,j,k)
                        write(5065+myid+ttt,*)rhof(i,j,k)
                        write(5070+myid+ttt,*)lmf11(i,j,k),lmf12(i,j,k),lmf13(i,j,k)
                        write(5075+myid+ttt,*)lmf21(i,j,k),lmf22(i,j,k),lmf23(i,j,k)
                        write(5080+myid+ttt,*)lmf31(i,j,k),lmf32(i,j,k),lmf33(i,j,k)
                    end do
                end do
            end do
        end if
        !
        ! now bar filteringon velocity and density already filtered with hat filter
        ! (observation: here also redefining the computation, the already defined soubroutine
        ! can be used again, without the b flag)
        !
        call filter04b(ucof,m21,myid,nproc,kparasta,kparaend)
        call filter04b(vcof,m22,myid,nproc,kparasta,kparaend)
        call filter04b(wcof,m23,myid,nproc,kparasta,kparaend)
        call filter04b(rhof,m33,myid,nproc,kparasta,kparaend)
        !
        call filter04b(lmf11,m11m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf12,m12m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf13,m13m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf21,m21m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf22,m22m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf23,m23m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf31,m31m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf32,m32m,myid,nproc,kparasta,kparaend)
        call filter04b(lmf33,m33m,myid,nproc,kparasta,kparaend)


        !
        ! now I transfer the "ghost" plane to include periodicity
        ! to complete the filtering on metric in zita (k+1 and k-1)
        !
        !  first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m21,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m33,4,kparasta,myid,nproc,kparasta,kparaend)
        !
        call buffer1g(m11m,5,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m12m,6,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m13m,7,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m21m,8,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m22m,9,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m23m,10,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,11,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,12,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,13,kparasta,myid,nproc, &
            kparasta,kparaend)
        !
        ! if I pute pe and not pem the periodicity is implicitly included in k
        ! and values in k=0 and k=jz+1 are defined in filter06
        !
        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif
        !
        call buffer2(rbuff1,m21,1,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,2,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,3,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33,4,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer2(rbuff1,m11m,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo
        !
        call buffer1g(m21,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m33,4,kparaend,myid,nproc,kparasta,kparaend)
        !
        call buffer1g(m11m,5,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,6,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,7,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,8,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,9,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,10,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,11,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,12,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,13,kparaend,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m21,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        !
        call buffer2(rbuff1,m11m,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparasta-1,myid,nproc, &
            kparasta,kparaend)


        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5140+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                        write(5145+myid+ttt,*)m33(i,j,k)
                        write(5150+myid+ttt,*)m11m(i,j,k),m12m(i,j,k),m13m(i,j,k)
                        write(5155+myid+ttt,*)m21m(i,j,k),m22m(i,j,k),m23m(i,j,k)
                        write(5160+myid+ttt,*)m31m(i,j,k),m32m(i,j,k),m33m(i,j,k)
                    end do
                end do
            end do
        end if

        !
        ! complete filtering in zita
        !
        call filter06b(m21,ucofb,kparasta,kparaend,myid,nproc)
        call filter06b(m22,vcofb,kparasta,kparaend,myid,nproc)
        call filter06b(m23,wcofb,kparasta,kparaend,myid,nproc)
        call filter06b(m33,rhofb,kparasta,kparaend,myid,nproc)
        !
        call filter06b(m11m,lmfb11,kparasta,kparaend,myid,nproc)
        call filter06b(m12m,lmfb12,kparasta,kparaend,myid,nproc)
        call filter06b(m13m,lmfb13,kparasta,kparaend,myid,nproc)
        call filter06b(m21m,lmfb21,kparasta,kparaend,myid,nproc)
        call filter06b(m22m,lmfb22,kparasta,kparaend,myid,nproc)
        call filter06b(m23m,lmfb23,kparasta,kparaend,myid,nproc)
        call filter06b(m31m,lmfb31,kparasta,kparaend,myid,nproc)
        call filter06b(m32m,lmfb32,kparasta,kparaend,myid,nproc)
        call filter06b(m33m,lmfb33,kparasta,kparaend,myid,nproc)

        !
        ! periodicity necessary for next filtering in the subroutine is done
        ! only for csi and eta, while zita is done once at the end
        !
        call periodic(ucofb,vcofb,wcofb)
        call periodic(rhofb,rhofb,rhofb)
        call periodic(lmfb11,lmfb12,lmfb13)
        call periodic(lmfb21,lmfb22,lmfb23)
        call periodic(lmfb31,lmfb32,lmfb33)
        !
        ! periodicity in zita
        !
        if(kp.eq.0)then
            !
            ! here I tranfer planes k=0 and k=jz+1
            ! so taht P0 knows k=jz and Pn-1 knows k=1
            !
            do m=1,40*(jx+2)*(jy+2)
                sbuff(m)=0.
                rbuff(m)=0.
            enddo

            if (myid.eq.nproc-1) then

                call buffer1g(ucofb,1,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vcofb,2,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wcofb,3,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofb,4,jz,myid,nproc,kparasta,kparaend)

                call buffer1g(lmfb11,5,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb12,6,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb13,7,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb21,8,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb22,9,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb23,10,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb31,11,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb32,12,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb33,13,jz,myid,nproc,kparasta,kparaend)

            else if (myid.eq.0) then

                call buffer1g(ucofb,1,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vcofb,2,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wcofb,3,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofb,4,1,myid,nproc,kparasta,kparaend)

                call buffer1g(lmfb11,5,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb12,6,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb13,7,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb21,8,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb22,9,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb23,10,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb31,11,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb32,12,1,myid,nproc,kparasta,kparaend)
                call buffer1g(lmfb33,13,1,myid,nproc,kparasta,kparaend)

            endif
            !
            ! exchange so
            ! P0 knows k=jz
            ! Pn-1 knows k=1
            !
            if (myid.eq.nproc-1) then

                call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,0,91, &
                    rbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,0,81, &
                    MPI_COMM_WORLD,status,ierr)

            else if (myid.eq.0) then

                call MPI_SENDRECV( &
                    sbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,81, &
                    rbuff1(1),13*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,91, &
                    MPI_COMM_WORLD,status,ierr)

            endif

            if (myid.eq.0) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)

                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano13,13,myid,nproc,kparasta,kparaend)

            else if (myid.eq.nproc-1) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)

                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano13,13,myid,nproc,kparasta,kparaend)

            endif

            ! now P0 knows k=jz plane
            ! and Pn-1 knows k=1 plane

            if(myid.eq.0)then

                do j=1,jy
                    do i=1,jx
                        ucofb(i,j,0)  = piano1(i,j)
                        vcofb(i,j,0)  = piano2(i,j)
                        wcofb(i,j,0)  = piano3(i,j)
                        rhofb(i,j,0)  = piano4(i,j)
                        lmfb11(i,j,0) = piano5(i,j)
                        lmfb12(i,j,0) = piano6(i,j)
                        lmfb13(i,j,0) = piano7(i,j)
                        lmfb21(i,j,0) = piano8(i,j)
                        lmfb22(i,j,0) = piano9(i,j)
                        lmfb23(i,j,0) = piano10(i,j)
                        lmfb31(i,j,0) = piano11(i,j)
                        lmfb32(i,j,0) = piano12(i,j)
                        lmfb33(i,j,0) = piano13(i,j)
                    end do
                end do

            endif


            if (myid.eq.nproc-1) then

                do j=1,jy
                    do i=1,jx
                        ucofb(i,j,jz+1)  = piano1(i,j)
                        vcofb(i,j,jz+1)  = piano2(i,j)
                        wcofb(i,j,jz+1)  = piano3(i,j)
                        rhofb(i,j,jz+1)  = piano4(i,j)
                        lmfb11(i,j,jz+1) = piano5(i,j)
                        lmfb12(i,j,jz+1) = piano6(i,j)
                        lmfb13(i,j,jz+1) = piano7(i,j)
                        lmfb21(i,j,jz+1) = piano8(i,j)
                        lmfb22(i,j,jz+1) = piano9(i,j)
                        lmfb23(i,j,jz+1) = piano10(i,j)
                        lmfb31(i,j,jz+1) = piano11(i,j)
                        lmfb32(i,j,jz+1) = piano12(i,j)
                        lmfb33(i,j,jz+1) = piano13(i,j)
                    end do
                end do

            endif

        endif
        !
        ! hat filtering on Uik and rho (last)
        !
        call filter01(rhofb,m33,myid,nproc,kparasta,kparaend)
        call filter01(ucofb,m21,myid,nproc,kparasta,kparaend)
        call filter01(vcofb,m22,myid,nproc,kparasta,kparaend)
        call filter01(wcofb,m23,myid,nproc,kparasta,kparaend)
        !
        call filter01(lmfb11,m11m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb12,m12m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb13,m13m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb21,m21m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb22,m22m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb23,m23m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb31,m31m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb32,m32m,myid,nproc,kparasta,kparaend)
        call filter01(lmfb33,m33m,myid,nproc,kparasta,kparaend)
        !
        ! now I need to transfer the "ghost" planes
        ! to impose periodicity to complete the filtering
        ! on metrica in zita (k+1 and k-1)
        !
        ! first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m33,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,4,kparasta,myid,nproc,kparasta,kparaend)

        call buffer1g(m11m,5,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,6,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,7,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,8,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,9,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,10,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,11,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,12,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,13,kparasta,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m33,1,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21,2,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,3,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,4,kparaend+1,myid,nproc, &
            kparasta,kparaend)

        call buffer2(rbuff1,m11m,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m33,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,4,kparaend,myid,nproc,kparasta,kparaend)

        call buffer1g(m11m,5,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,6,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,7,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,8,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,9,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,10,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,11,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,12,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,13,kparaend,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then


            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m33,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)

        call buffer2(rbuff1,m11m,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparasta-1,myid,nproc, &
            kparasta,kparaend)



        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5200+myid+ttt,*)m33(i,j,k)
                        write(5205+myid+ttt,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                        write(5210+myid+ttt,*)m11m(i,j,k),m12m(i,j,k),m13m(i,j,k)
                        write(5215+myid+ttt,*)m21m(i,j,k),m22m(i,j,k),m23m(i,j,k)
                        write(5220+myid+ttt,*)m31m(i,j,k),m32m(i,j,k),m33m(i,j,k)
                    end do
                end do
            end do
        end if
        !
        ! compete filtering in zita
        !
        call filter03(m33,rhof,kparasta,kparaend,myid,nproc)
        call filter03(m21,ucof,kparasta,kparaend,myid,nproc)
        call filter03(m22,vcof,kparasta,kparaend,myid,nproc)
        call filter03(m23,wcof,kparasta,kparaend,myid,nproc)

        call filter03(m11m,lmf11,kparasta,kparaend,myid,nproc)
        call filter03(m12m,lmf12,kparasta,kparaend,myid,nproc)
        call filter03(m13m,lmf13,kparasta,kparaend,myid,nproc)
        call filter03(m21m,lmf21,kparasta,kparaend,myid,nproc)
        call filter03(m22m,lmf22,kparasta,kparaend,myid,nproc)
        call filter03(m23m,lmf23,kparasta,kparaend,myid,nproc)
        call filter03(m31m,lmf31,kparasta,kparaend,myid,nproc)
        call filter03(m32m,lmf32,kparasta,kparaend,myid,nproc)
        call filter03(m33m,lmf33,kparasta,kparaend,myid,nproc)
        !
        ! cartesian component derived from controvariant
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = ucof(i,j,k)
                    pc2(i,j,k)  = vcof(i,j,k)
                    pc3(i,j,k)  = wcof(i,j,k)
                    ap11(i,j,k) = lmf11(i,j,k)
                    ap12(i,j,k) = lmf12(i,j,k)
                    ap13(i,j,k) = lmf13(i,j,k)
                    ap21(i,j,k) = lmf21(i,j,k)
                    ap22(i,j,k) = lmf22(i,j,k)
                    ap23(i,j,k) = lmf23(i,j,k)
                    ap31(i,j,k) = lmf31(i,j,k)
                    ap32(i,j,k) = lmf32(i,j,k)
                    ap33(i,j,k) = lmf33(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    uf(i,j,k) = pp1(i,j,k)
                    vf(i,j,k) = pp2(i,j,k)
                    wf(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        !
        ! now I have (Uik)FfF and (uij)FfF, I make the product
        ! to obtain the first erm of Bik
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    bl11(i,j,k)=uf(i,j,k)*ucof(i,j,k)       !uU
                    bl12(i,j,k)=vf(i,j,k)*ucof(i,j,k)       !vU
                    bl13(i,j,k)=wf(i,j,k)*ucof(i,j,k)       !wU
                    bl21(i,j,k)=uf(i,j,k)*vcof(i,j,k)       !uV
                    bl22(i,j,k)=vf(i,j,k)*vcof(i,j,k)       !vV
                    bl23(i,j,k)=wf(i,j,k)*vcof(i,j,k)       !wV
                    bl31(i,j,k)=uf(i,j,k)*wcof(i,j,k)       !uW
                    bl32(i,j,k)=vf(i,j,k)*wcof(i,j,k)       !vW
                    bl33(i,j,k)=wf(i,j,k)*wcof(i,j,k)       !wW
                    blrho11(i,j,k)=rhof(i,j,k)*ucof(i,j,k)   !rhoU
                    blrho22(i,j,k)=rhof(i,j,k)*vcof(i,j,k)   !rhoV
                    blrho33(i,j,k)=rhof(i,j,k)*wcof(i,j,k)   !rhoW
                enddo
            enddo
        enddo
        !
        ! the second part of Bik
        !
        do k=kparasta-1,kparaend+1
            do j=0,n2+1
                do i=0,n1+1
                    p0a(i,j,k) = rho(i,j,k)
                end do
            end do
        end do


        call filter01(p0a,m33,myid,nproc,kparasta,kparaend)

        call filter01(uco,m21,myid,nproc,kparasta,kparaend)
        call filter01(vco,m22,myid,nproc,kparasta,kparaend)
        call filter01(wco,m23,myid,nproc,kparasta,kparaend)
        !
        call filter01(apcsx,m11m,myid,nproc,kparasta,kparaend)
        call filter01(apcsy,m12m,myid,nproc,kparasta,kparaend)
        call filter01(apcsz,m13m,myid,nproc,kparasta,kparaend)
        call filter01(apetx,m21m,myid,nproc,kparasta,kparaend)
        call filter01(apety,m22m,myid,nproc,kparasta,kparaend)
        call filter01(apetz,m23m,myid,nproc,kparasta,kparaend)
        call filter01(apztx,m31m,myid,nproc,kparasta,kparaend)
        call filter01(apzty,m32m,myid,nproc,kparasta,kparaend)
        call filter01(apztz,m33m,myid,nproc,kparasta,kparaend)
        !
        ! I need to exchange the "ghost" planes to impose periodicity
        ! to complete filtering in zita on metrica (k+1 and k-1)
        !
        ! first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m33,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,4,kparasta,myid,nproc,kparasta,kparaend)

        call buffer1g(m11m,5,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,6,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,7,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,8,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,9,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,10,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,11,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,12,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,13,kparasta,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m33,1,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21,2,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,3,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,4,kparaend+1,myid,nproc, &
            kparasta,kparaend)

        call buffer2(rbuff1,m11m,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(m33,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m21,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m22,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(m23,4,kparaend,myid,nproc,kparasta,kparaend)

        call buffer1g(m11m,5,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m12m,6,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m13m,7,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m21m,8,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m22m,9,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m23m,10,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m31m,11,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m32m,12,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(m33m,13,kparaend,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),13*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),13*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,m33,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)

        call buffer2(rbuff1,m11m,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m12m,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m13m,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m21m,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m22m,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m23m,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m31m,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m32m,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,m33m,13,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        !
        ! complete the filtering in zita
        !
        call filter03(m33,rhof,kparasta,kparaend,myid,nproc)
        call filter03(m21,ucof,kparasta,kparaend,myid,nproc)
        call filter03(m22,vcof,kparasta,kparaend,myid,nproc)
        call filter03(m23,wcof,kparasta,kparaend,myid,nproc)

        call filter03(m11m,lmf11,kparasta,kparaend,myid,nproc)
        call filter03(m12m,lmf12,kparasta,kparaend,myid,nproc)
        call filter03(m13m,lmf13,kparasta,kparaend,myid,nproc)
        call filter03(m21m,lmf21,kparasta,kparaend,myid,nproc)
        call filter03(m22m,lmf22,kparasta,kparaend,myid,nproc)
        call filter03(m23m,lmf23,kparasta,kparaend,myid,nproc)
        call filter03(m31m,lmf31,kparasta,kparaend,myid,nproc)
        call filter03(m32m,lmf32,kparasta,kparaend,myid,nproc)
        call filter03(m33m,lmf33,kparasta,kparaend,myid,nproc)


        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5225+myid+ttt,*)rhof(i,j,k),ucof(i,j,k),vcof(i,j,k), &
                            wcof(i,j,k)
                        write(5230+myid+ttt,*)lmf11(i,j,k),lmf12(i,j,k),lmf13(i,j,k)
                        write(5235+myid+ttt,*)lmf21(i,j,k),lmf22(i,j,k),lmf23(i,j,k)
                        write(5240+myid+ttt,*)lmf31(i,j,k),lmf32(i,j,k),lmf33(i,j,k)
                    end do
                end do
            end do
        end if
        !
        ! derive the filtered cartesian component
        ! from the controvariant ones
        ! (I need this to make the product that
        ! will be filtered to obtain the second term of Bik)
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = ucof(i,j,k)
                    pc2(i,j,k)  = vcof(i,j,k)
                    pc3(i,j,k)  = wcof(i,j,k)
                    ap11(i,j,k) = lmf11(i,j,k)
                    ap12(i,j,k) = lmf12(i,j,k)
                    ap13(i,j,k) = lmf13(i,j,k)
                    ap21(i,j,k) = lmf21(i,j,k)
                    ap22(i,j,k) = lmf22(i,j,k)
                    ap23(i,j,k) = lmf23(i,j,k)
                    ap31(i,j,k) = lmf31(i,j,k)
                    ap32(i,j,k) = lmf32(i,j,k)
                    ap33(i,j,k) = lmf33(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    uf(i,j,k) = pp1(i,j,k)
                    vf(i,j,k) = pp2(i,j,k)
                    wf(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        !
        ! the product that will be filtered
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    ufucof(i,j,k)=uf(i,j,k)*ucof(i,j,k)
                    vfucof(i,j,k)=vf(i,j,k)*ucof(i,j,k)
                    wfucof(i,j,k)=wf(i,j,k)*ucof(i,j,k)
                    ufvcof(i,j,k)=uf(i,j,k)*vcof(i,j,k)
                    vfvcof(i,j,k)=vf(i,j,k)*vcof(i,j,k)
                    wfvcof(i,j,k)=wf(i,j,k)*vcof(i,j,k)
                    ufwcof(i,j,k)=uf(i,j,k)*wcof(i,j,k)
                    vfwcof(i,j,k)=vf(i,j,k)*wcof(i,j,k)
                    wfwcof(i,j,k)=wf(i,j,k)*wcof(i,j,k)
                    rhofucof(i,j,k)=rhof(i,j,k)*ucof(i,j,k)
                    rhofvcof(i,j,k)=rhof(i,j,k)*vcof(i,j,k)
                    rhofwcof(i,j,k)=rhof(i,j,k)*wcof(i,j,k)
                enddo
            enddo
        enddo

        !
        ! periodicity (necessary for the next filtering) is done
        ! for csi and eta, for zita it is done at the end just one time
        !
        call periodic(ufucof,vfucof,wfucof)
        call periodic(ufvcof,vfvcof,wfvcof)
        call periodic(ufwcof,vfwcof,wfwcof)
        call periodic(rhofucof,rhofvcof,rhofwcof)
        !
        ! periodicity in zita
        !
        if(kp.eq.0)then
            !
            ! I need to transfer planes in k=0 and k=jz+1
            ! so P0 knows k=jz and Pn-1 knows k=1
            !
            do m=1,40*(jx+2)*(jy+2)
                sbuff(m)=0.
                rbuff(m)=0.
            enddo

            if (myid.eq.nproc-1) then

                call buffer1g(ufucof,1,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vfucof,2,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wfucof,3,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ufvcof,4,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vfvcof,5,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wfvcof,6,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(ufwcof,7,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vfwcof,8,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wfwcof,9,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofucof,10,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofvcof,11,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofwcof,12,jz,myid,nproc,kparasta,kparaend)

            else if (myid.eq.0) then

                call buffer1g(ufucof,1,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vfucof,2,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wfucof,3,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ufvcof,4,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vfvcof,5,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wfvcof,6,1,myid,nproc,kparasta,kparaend)
                call buffer1g(ufwcof,7,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vfwcof,8,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wfwcof,9,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofucof,10,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofvcof,11,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhofwcof,12,1,myid,nproc,kparasta,kparaend)

            endif
            !
            ! now the exchange so that
            ! P0 knows k=jz
            ! Pn-1 knows k=1
            !
            if (myid.eq.nproc-1) then

                call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,0,91, &
                    rbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,0,81, &
                    MPI_COMM_WORLD,status,ierr)

            else if (myid.eq.0) then

                call MPI_SENDRECV( &
                    sbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,81, &
                    rbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,91, &
                    MPI_COMM_WORLD,status,ierr)

            endif

            if (myid.eq.0) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)


            else if (myid.eq.nproc-1) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)

            endif
            !
            ! now P0 knows planes k=jz
            ! and Pn-1 knows k=1
            !

            if (myid.eq.0) then

                do j=1,jy
                    do i=1,jx
                        ufucof(i,j,0)   = piano1(i,j)
                        vfucof(i,j,0)   = piano2(i,j)
                        wfucof(i,j,0)   = piano3(i,j)
                        ufvcof(i,j,0)   = piano4(i,j)
                        vfvcof(i,j,0)   = piano5(i,j)
                        wfvcof(i,j,0)   = piano6(i,j)
                        ufwcof(i,j,0)   = piano7(i,j)
                        vfwcof(i,j,0)   = piano8(i,j)
                        wfwcof(i,j,0)   = piano9(i,j)
                        rhofucof(i,j,0) = piano10(i,j)
                        rhofvcof(i,j,0) = piano11(i,j)
                        rhofwcof(i,j,0) = piano12(i,j)
                    end do
                end do

            endif

            if (myid.eq.nproc-1) then

                do j=1,jy
                    do i=1,jx
                        ufucof(i,j,jz+1)   = piano1(i,j)
                        vfucof(i,j,jz+1)   = piano2(i,j)
                        wfucof(i,j,jz+1)   = piano3(i,j)
                        ufvcof(i,j,jz+1)   = piano4(i,j)
                        vfvcof(i,j,jz+1)   = piano5(i,j)
                        wfvcof(i,j,jz+1)   = piano6(i,j)
                        ufwcof(i,j,jz+1)   = piano7(i,j)
                        vfwcof(i,j,jz+1)   = piano8(i,j)
                        wfwcof(i,j,jz+1)   = piano9(i,j)
                        rhofucof(i,j,jz+1) = piano10(i,j)
                        rhofvcof(i,j,jz+1) = piano11(i,j)
                        rhofwcof(i,j,jz+1) = piano12(i,j)
                    end do
                end do

            end if

        !
        endif
        !
        ! compute the filtered product, I need them to
        ! compute the second term of Bik
        !
        call filter04b(ufucof,l11,myid,nproc,kparasta,kparaend)   !(uU)f
        call filter04b(ufvcof,l12,myid,nproc,kparasta,kparaend)   !(uV)f
        call filter04b(ufwcof,l13,myid,nproc,kparasta,kparaend)   !(uW)f
        call filter04b(vfucof,l21,myid,nproc,kparasta,kparaend)   !(vU)f
        call filter04b(vfvcof,l22,myid,nproc,kparasta,kparaend)   !(vV)f
        call filter04b(vfwcof,l23,myid,nproc,kparasta,kparaend)   !(vW)f
        call filter04b(wfucof,l31,myid,nproc,kparasta,kparaend)   !(wU)f
        call filter04b(wfvcof,l32,myid,nproc,kparasta,kparaend)   !(wV)f
        call filter04b(wfwcof,l33,myid,nproc,kparasta,kparaend)   !(wW)f
        call filter04b(rhofucof,lrho11,myid,nproc,kparasta,kparaend)   !(rhoU)f
        call filter04b(rhofvcof,lrho22,myid,nproc,kparasta,kparaend)   !(rhoV)f
        call filter04b(rhofwcof,lrho33,myid,nproc,kparasta,kparaend)   !(rhoW)f

        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5250+myid+ttt,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                        write(5255+myid+ttt,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                        write(5260+myid+ttt,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                        write(5265+myid+ttt,*)lrho11(i,j,k),lrho22(i,j,k),lrho33(i,j,k)
                    end do
                end do
            end do
        end if

        !
        ! I need to exchange the "ghost" planes to impose periodicity
        ! to complete filtering in zita on metrica (k+1 and k-1)
        !
        ! first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(l11,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l12,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l13,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l21,4,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l22,5,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l23,6,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l31,7,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l32,8,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l33,9,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(lrho11,10,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho22,11,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho33,12,kparasta,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,l11,1,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l12,2,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l13,3,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l21,4,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l22,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l23,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l31,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l32,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l33,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho11,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho22,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho33,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(l11,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l12,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l13,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l21,4,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l22,5,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l23,6,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l31,7,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l32,8,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l33,9,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(lrho11,10,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho22,11,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho33,12,kparaend,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,l11,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l12,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l13,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l21,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l22,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l23,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l31,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l32,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l33,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho11,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho22,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho33,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)



        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5270+myid+ttt,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                        write(5275+myid+ttt,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                        write(5280+myid+ttt,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                        write(5285+myid+ttt,*)lrho11(i,j,k),lrho22(i,j,k),lrho33(i,j,k)
                    end do
                end do
            end do
        end if
        !
        ! complete filtering in zita
        !
        call filter06b(l11,uucof,kparasta,kparaend,myid,nproc)
        call filter06b(l12,uvcof,kparasta,kparaend,myid,nproc)
        call filter06b(l13,uwcof,kparasta,kparaend,myid,nproc)
        call filter06b(l21,vucof,kparasta,kparaend,myid,nproc)
        call filter06b(l22,vvcof,kparasta,kparaend,myid,nproc)
        call filter06b(l23,vwcof,kparasta,kparaend,myid,nproc)
        call filter06b(l31,wucof,kparasta,kparaend,myid,nproc)
        call filter06b(l32,wvcof,kparasta,kparaend,myid,nproc)
        call filter06b(l33,wwcof,kparasta,kparaend,myid,nproc)
        call filter06b(lrho11,rhoucof,kparasta,kparaend,myid,nproc)
        call filter06b(lrho22,rhovcof,kparasta,kparaend,myid,nproc)
        call filter06b(lrho33,rhowcof,kparasta,kparaend,myid,nproc)
        !
        ! periodicity (necessary for the next filtering) is done
        ! for csi and eta, for zita it is done at the end just one time
        !
        call periodic(uucof,uvcof,uwcof)
        call periodic(vucof,vvcof,vwcof)
        call periodic(wucof,wvcof,wwcof)
        call periodic(rhoucof,rhovcof,rhowcof)

        !
        ! periodicity in zita
        !
        if(kp.eq.0)then
            !
            ! I need to transfer planes in k=0 and k=jz+1
            ! so P0 knows k=jz and Pn-1 knows k=1
            !
            do m=1,40*(jx+2)*(jy+2)
                sbuff(m)=0.
                rbuff(m)=0.
            enddo

            if (myid.eq.nproc-1) then

                call buffer1g(uucof,1,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vucof,2,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wucof,3,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(uvcof,4,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vvcof,5,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wvcof,6,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(uwcof,7,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(vwcof,8,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(wwcof,9,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhoucof,10,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhovcof,11,jz,myid,nproc,kparasta,kparaend)
                call buffer1g(rhowcof,12,jz,myid,nproc,kparasta,kparaend)

            else if (myid.eq.0) then

                call buffer1g(uucof,1,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vucof,2,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wucof,3,1,myid,nproc,kparasta,kparaend)
                call buffer1g(uvcof,4,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vvcof,5,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wvcof,6,1,myid,nproc,kparasta,kparaend)
                call buffer1g(uwcof,7,1,myid,nproc,kparasta,kparaend)
                call buffer1g(vwcof,8,1,myid,nproc,kparasta,kparaend)
                call buffer1g(wwcof,9,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhoucof,10,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhovcof,11,1,myid,nproc,kparasta,kparaend)
                call buffer1g(rhowcof,12,1,myid,nproc,kparasta,kparaend)

            endif
            !
            ! now the exchange so that
            ! P0 knows k=jz
            ! Pn-1 knows k=1
            !
            if (myid.eq.nproc-1) then

                call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,0,91, &
                    rbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,0,81, &
                    MPI_COMM_WORLD,status,ierr)

            else if (myid.eq.0) then

                call MPI_SENDRECV( &
                    sbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,81, &
                    rbuff1(1),12*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,91, &
                    MPI_COMM_WORLD,status,ierr)

            endif

            if (myid.eq.0) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)

            else if (myid.eq.nproc-1) then

                call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano4,4,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano5,5,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano6,6,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano7,7,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano8,8,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano9,9,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano10,10,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano11,11,myid,nproc,kparasta,kparaend)
                call buffer2g(rbuff1,piano12,12,myid,nproc,kparasta,kparaend)

            endif
            !
            ! now P0 knows planes k=jz
            ! and Pn-1 knows k=1
            !
            if (myid.eq.0) then

                do j=1,jy
                    do i=1,jx
                        uucof(i,j,0)   = piano1(i,j)
                        vucof(i,j,0)   = piano2(i,j)
                        wucof(i,j,0)   = piano3(i,j)
                        uvcof(i,j,0)   = piano4(i,j)
                        vvcof(i,j,0)   = piano5(i,j)
                        wvcof(i,j,0)   = piano6(i,j)
                        uwcof(i,j,0)   = piano7(i,j)
                        vwcof(i,j,0)   = piano8(i,j)
                        wwcof(i,j,0)   = piano9(i,j)
                        rhoucof(i,j,0) = piano10(i,j)
                        rhovcof(i,j,0) = piano11(i,j)
                        rhowcof(i,j,0) = piano12(i,j)
                    end do
                end do

            end if

            if (myid.eq.nproc-1) then

                do j=1,jy
                    do i=1,jx

                        uucof(i,j,jz+1)   = piano1(i,j)
                        vucof(i,j,jz+1)   = piano2(i,j)
                        wucof(i,j,jz+1)   = piano3(i,j)
                        uvcof(i,j,jz+1)   = piano4(i,j)
                        vvcof(i,j,jz+1)   = piano5(i,j)
                        wvcof(i,j,jz+1)   = piano6(i,j)
                        uwcof(i,j,jz+1)   = piano7(i,j)
                        vwcof(i,j,jz+1)   = piano8(i,j)
                        wwcof(i,j,jz+1)   = piano9(i,j)
                        rhoucof(i,j,jz+1) = piano10(i,j)
                        rhovcof(i,j,jz+1) = piano11(i,j)
                        rhowcof(i,j,jz+1) = piano12(i,j)

                    end do
                end do

            end if

        !
        endif


        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5300+myid+ttt,*)uucof(i,j,k),vucof(i,j,k),wucof(i,j,k)
                        write(5305+myid+ttt,*)uvcof(i,j,k),vvcof(i,j,k),wvcof(i,j,k)
                        write(5310+myid+ttt,*)uwcof(i,j,k),vwcof(i,j,k),wwcof(i,j,k)
                        write(5315+myid+ttt,*)rhoucof(i,j,k),rhovcof(i,j,k), &
                            rhowcof(i,j,k)
                    end do
                end do
            end do
        end if

        !
        ! last filtering hat
        !
        call filter01(uucof,l11,myid,nproc,kparasta,kparaend)
        call filter01(uvcof,l12,myid,nproc,kparasta,kparaend)
        call filter01(uwcof,l13,myid,nproc,kparasta,kparaend)
        call filter01(vucof,l21,myid,nproc,kparasta,kparaend)
        call filter01(vvcof,l22,myid,nproc,kparasta,kparaend)
        call filter01(vwcof,l23,myid,nproc,kparasta,kparaend)
        call filter01(wucof,l31,myid,nproc,kparasta,kparaend)
        call filter01(wvcof,l32,myid,nproc,kparasta,kparaend)
        call filter01(wwcof,l33,myid,nproc,kparasta,kparaend)
        call filter01(rhoucof,lrho11,myid,nproc,kparasta,kparaend)
        call filter01(rhovcof,lrho22,myid,nproc,kparasta,kparaend)
        call filter01(rhowcof,lrho33,myid,nproc,kparasta,kparaend)
        !
        ! I need to exchange the "ghost" planes to impose periodicity
        ! to complete filtering in zita on metrica (k+1 and k-1)
        !
        ! first kparasta
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(l11,1,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l12,2,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l13,3,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l21,4,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l22,5,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l23,6,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l31,7,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l32,8,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(l33,9,kparasta,myid,nproc,kparasta,kparaend)
        call buffer1g(lrho11,10,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho22,11,kparasta,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho33,12,kparasta,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,tagls, &
                rbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(leftpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,l11,1,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l12,2,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l13,3,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l21,4,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l22,5,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l23,6,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l31,7,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l32,8,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l33,9,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho11,10,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho22,11,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho33,12,kparaend+1,myid,nproc, &
            kparasta,kparaend)
        !
        ! now kparaend
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
        enddo

        call buffer1g(l11,1,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l12,2,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l13,3,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l21,4,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l22,5,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l23,6,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l31,7,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l32,8,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(l33,9,kparaend,myid,nproc,kparasta,kparaend)
        call buffer1g(lrho11,10,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho22,11,kparaend,myid,nproc, &
            kparasta,kparaend)
        call buffer1g(lrho33,12,kparaend,myid,nproc, &
            kparasta,kparaend)

        if (kp.eq.0) then

            call MPI_SENDRECV(sbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpe,tagrs, &
                rbuff1(1),12*(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpe,taglr, &
                MPI_COMM_WORLD,status,ierr)

        else if (kp.eq.1) then

            if(rightpem /= MPI_PROC_NULL) then
                call MPI_SEND(sbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
                call MPI_RECV(rbuff1(1),12*(jx+2)*(jy+2), &
                    MPI_REAL_SD,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req1,istatus,ierr)
            endif
            if(leftpem /= MPI_PROC_NULL) then
            !      call MPI_WAIT(req2,istatus,ierr)
            endif


        endif

        call buffer2(rbuff1,l11,1,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l12,2,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l13,3,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l21,4,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l22,5,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l23,6,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l31,7,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l32,8,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,l33,9,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho11,10,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho22,11,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        call buffer2(rbuff1,lrho33,12,kparasta-1,myid,nproc, &
            kparasta,kparaend)
        !
        ! complete filtering in zita
        !
        call filter03(l11,uucofb,kparasta,kparaend,myid,nproc)
        call filter03(l12,uvcofb,kparasta,kparaend,myid,nproc)
        call filter03(l13,uwcofb,kparasta,kparaend,myid,nproc)
        call filter03(l21,vucofb,kparasta,kparaend,myid,nproc)
        call filter03(l22,vvcofb,kparasta,kparaend,myid,nproc)
        call filter03(l23,vwcofb,kparasta,kparaend,myid,nproc)
        call filter03(l31,wucofb,kparasta,kparaend,myid,nproc)
        call filter03(l32,wvcofb,kparasta,kparaend,myid,nproc)
        call filter03(l33,wwcofb,kparasta,kparaend,myid,nproc)
        call filter03(lrho11,rhoucofb,kparasta,kparaend,myid,nproc)
        call filter03(lrho22,rhovcofb,kparasta,kparaend,myid,nproc)
        call filter03(lrho33,rhowcofb,kparasta,kparaend,myid,nproc)


        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    bl11(i,j,k)=uucofb(i,j,k)-bl11(i,j,k)
                    bl12(i,j,k)=vucofb(i,j,k)-bl12(i,j,k)
                    bl13(i,j,k)=wucofb(i,j,k)-bl13(i,j,k)
                    bl21(i,j,k)=uvcofb(i,j,k)-bl21(i,j,k)
                    bl22(i,j,k)=vvcofb(i,j,k)-bl22(i,j,k)
                    bl23(i,j,k)=wvcofb(i,j,k)-bl23(i,j,k)
                    bl31(i,j,k)=uwcofb(i,j,k)-bl31(i,j,k)
                    bl32(i,j,k)=vwcofb(i,j,k)-bl32(i,j,k)
                    bl33(i,j,k)=wwcofb(i,j,k)-bl33(i,j,k)
                    blrho11(i,j,k)=rhoucofb(i,j,k)-blrho11(i,j,k)
                    blrho22(i,j,k)=rhovcofb(i,j,k)-blrho22(i,j,k)
                    blrho33(i,j,k)=rhowcofb(i,j,k)-blrho33(i,j,k)
                enddo
            enddo
        enddo

        !
        ! difference between Bik-Aik
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    al11(i,j,k)=bl11(i,j,k)-al11(i,j,k)
                    al12(i,j,k)=bl12(i,j,k)-al12(i,j,k)
                    al13(i,j,k)=bl13(i,j,k)-al13(i,j,k)
                    al21(i,j,k)=bl21(i,j,k)-al21(i,j,k)
                    al22(i,j,k)=bl22(i,j,k)-al22(i,j,k)
                    al23(i,j,k)=bl23(i,j,k)-al23(i,j,k)
                    al31(i,j,k)=bl31(i,j,k)-al31(i,j,k)
                    al32(i,j,k)=bl32(i,j,k)-al32(i,j,k)
                    al33(i,j,k)=bl33(i,j,k)-al33(i,j,k)
                    alrho11(i,j,k)=blrho11(i,j,k)-alrho11(i,j,k)
                    alrho22(i,j,k)=blrho22(i,j,k)-alrho22(i,j,k)
                    alrho33(i,j,k)=blrho33(i,j,k)-alrho33(i,j,k)
                enddo
            enddo
        enddo

        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5320+myid+ttt,*)al11(i,j,k),al12(i,j,k),al13(i,j,k)
                        write(5325+myid+ttt,*)al21(i,j,k),al22(i,j,k),al23(i,j,k)
                        write(5330+myid+ttt,*)al31(i,j,k),al32(i,j,k),al33(i,j,k)
                        write(5335+myid+ttt,*)alrho11(i,j,k),alrho22(i,j,k), &
                            alrho33(i,j,k)
                    end do
                end do
            end do
        end if
        !
        !-----------------------------------------------------------------------
        !     END SCALE SIMILAR
        !-----------------------------------------------------------------------
        !
        ! recompute the cartesian component Mij starting from
        ! the controvariant already computed
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = am11(i,j,k)
                    pc2(i,j,k)  = am12(i,j,k)
                    pc3(i,j,k)  = am13(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil11(i,j,k) = pp1(i,j,k)
                    fil12(i,j,k) = pp2(i,j,k)
                    fil13(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = am21(i,j,k)
                    pc2(i,j,k)  = am22(i,j,k)
                    pc3(i,j,k)  = am23(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil21(i,j,k) = pp1(i,j,k)
                    fil22(i,j,k) = pp2(i,j,k)
                    fil23(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do


        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = am31(i,j,k)
                    pc2(i,j,k)  = am32(i,j,k)
                    pc3(i,j,k)  = am33(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil31(i,j,k) = pp1(i,j,k)
                    fil32(i,j,k) = pp2(i,j,k)
                    fil33(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        !
        ! cartesian component Nij from the difference Bik-Aik
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = al11(i,j,k)
                    pc2(i,j,k)  = al21(i,j,k)
                    pc3(i,j,k)  = al31(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fipp11(i,j,k) = pp1(i,j,k)
                    fipp12(i,j,k) = pp2(i,j,k)
                    fipp13(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = al12(i,j,k)
                    pc2(i,j,k)  = al22(i,j,k)
                    pc3(i,j,k)  = al32(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fipp21(i,j,k) = pp1(i,j,k)
                    fipp22(i,j,k) = pp2(i,j,k)
                    fipp23(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = al13(i,j,k)
                    pc2(i,j,k)  = al23(i,j,k)
                    pc3(i,j,k)  = al33(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fipp31(i,j,k) = pp1(i,j,k)
                    fipp32(i,j,k) = pp2(i,j,k)
                    fipp33(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do
        !
        ! term LijMij-NijMij (first term is already in num)
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    num(i,j,k)=num(i,j,k)-cb(j)* &
                        (fil11(i,j,k)*fipp11(i,j,k)+fil12(i,j,k)*fipp12(i,j,k)+ &
                        fil13(i,j,k)*fipp13(i,j,k)+fil21(i,j,k)*fipp21(i,j,k)+ &
                        fil22(i,j,k)*fipp22(i,j,k)+fil32(i,j,k)*fipp32(i,j,k)+ &
                        fil31(i,j,k)*fipp31(i,j,k)+fil23(i,j,k)*fipp23(i,j,k)+ &
                        fil33(i,j,k)*fipp33(i,j,k))
                enddo
            enddo
        enddo

        if (debugg.eq.1) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5340+myid+ttt,*)num(i,j,k)
                    end do
                end do
            end do
        end if
        !
        ! recompute the cartesian component rhoMij starting
        ! from the controvariant already computed
        !
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = amrho11(i,j,k)
                    pc2(i,j,k)  = amrho22(i,j,k)
                    pc3(i,j,k)  = amrho33(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil11(i,j,k) = pp1(i,j,k)
                    fil22(i,j,k) = pp2(i,j,k)
                    fil33(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do
        !
        ! cartesian component rhoNij from the difference Bik-Aik
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = alrho11(i,j,k)
                    pc2(i,j,k)  = alrho22(i,j,k)
                    pc3(i,j,k)  = alrho33(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fipp11(i,j,k) = pp1(i,j,k)
                    fipp22(i,j,k) = pp2(i,j,k)
                    fipp33(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        !
        ! term LijMij-NijMij for density (first term already in numrho)
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    numrho(i,j,k)=numrho(i,j,k)-cbrho(j)* &
                        (fil11(i,j,k)*fipp11(i,j,k)+ &
                        fil22(i,j,k)*fipp22(i,j,k)+ &
                        fil33(i,j,k)*fipp33(i,j,k))
                enddo
            enddo
        enddo


    end if ! end first block scale similar inmod.eq.1
    !
    !-----------------------------------------------------------------------
    !     END SCALE SIMILAR
    !-----------------------------------------------------------------------
    !
    !
    !
    !-----------------------------------------------------------------------
    !     COMPUTE THE CONSTANT C Crho LAGRANGIAN
    !-----------------------------------------------------------------------
    !
    ! compute the constants C and Crho that in the eddy viscosity (smagorinsky)
    ! model are necessary to:
    ! - tau ij ---> sgs stress eddy viscosity and diffusivity
    ! - eddy viscosity and diffusivity
    !
    ! first numerator and denominator average
    ! on homogeneous plane (x and z)
    ! (plane averaged model)
    ! or just z is defined at the beginning of the subroutine

    ! this part is different for dinamico and lagrangian
    if (nsgs==3) then
        !
        if(ktime==1)then
            do j=1,jy

                nummed_loc(j)=0.
                denmed_loc(j)=0.
                numrhomed_loc(j)=0.
                denrhomed_loc(j)=0.

                do k=kparasta,kparaend
                    do i=1,jx
                        nummed_loc(j)=nummed_loc(j)+num(i,j,k)
                        denmed_loc(j)=denmed_loc(j)+den(i,j,k)
                        numrhomed_loc(j)=numrhomed_loc(j)+numrho(i,j,k)
                        denrhomed_loc(j)=denrhomed_loc(j)+denrho(i,j,k)
                    enddo
                enddo

            enddo
            !
            ! now local sum and send the result to each proc
            ! for every plane
            !
            do m=1,40*(jx+2)*(jy+2)
                sbuffd(m)=0.
                rbuffd(m)=0.
            enddo

            call buffvect1d(sbuffd,nummed_loc,1)
            call buffvect1d(sbuffd,denmed_loc,2)
            call buffvect1d(sbuffd,numrhomed_loc,3)
            call buffvect1d(sbuffd,denrhomed_loc,4)

            call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),4*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            call buffvect2d(rbuffd,nummed,1)
            call buffvect2d(rbuffd,denmed,2)
            call buffvect2d(rbuffd,numrhomed,3)
            call buffvect2d(rbuffd,denrhomed,4)
            !
            ! finally the mean value, sum divided by the number of points
            ! on j plane (known by all procs)
            !
            rzer=0.
            do j=1,jy

                nummed(j)=nummed(j)/somma
                denmed(j)=denmed(j)/somma
                numrhomed(j)=numrhomed(j)/somma
                denrhomed(j)=denrhomed(j)/somma

                ! check on density denominator ---> see grg and grgrho
                !
                ! then the constant (already multiplied by -1)

                do k=kparasta,kparaend
                    do i=1,jx


                        alalm(i,j,k)=nummed(j)
                        alamm(i,j,k)=denmed(j)

                        alalmrho(i,j,k)=numrhomed(j)
                        alammrho(i,j,k)=denrhomed(j)


                       !appo1(i,j,k)=nummed(j)
                       !appo2(i,j,k)=denmed(j)
                       !appo1rho(i,j,k)=numrhomed(j)
                       !appo2rho(i,j,k)=denrhomed(j)
                    enddo
                enddo

            enddo


        else

            if(leftpem/=MPI_PROC_NULL) then
                call MPI_SEND(alalm(1,1,kparasta),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_RECV(alalm(1,1,kparaend+1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem/=MPI_PROC_NULL) then
                call MPI_SEND(alamm(1,1,kparasta),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_RECV(alamm(1,1,kparaend+1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_SEND(alalm(1,1,kparaend),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(leftpem/=MPI_PROC_NULL) then
                call MPI_RECV(alalm(1,1,kparasta-1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_SEND(alamm(1,1,kparaend),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(leftpem/=MPI_PROC_NULL) then
                call MPI_RECV(alamm(1,1,kparasta-1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif
            !-----------------------------
            if(leftpem/=MPI_PROC_NULL) then
                call MPI_SEND(alalmrho(1,1,kparasta),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_RECV(alalmrho(1,1,kparaend+1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(leftpem/=MPI_PROC_NULL) then
                call MPI_SEND(alammrho(1,1,kparasta),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,tagls, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_RECV(alammrho(1,1,kparaend+1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_SEND(alalmrho(1,1,kparaend),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(leftpem/=MPI_PROC_NULL) then
                call MPI_RECV(alalmrho(1,1,kparasta-1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            if(rightpem/=MPI_PROC_NULL) then
                call MPI_SEND(alammrho(1,1,kparaend),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,rightpem,tagrs, &
                    MPI_COMM_WORLD,ierr)
            endif

            if(leftpem/=MPI_PROC_NULL) then
                call MPI_RECV(alammrho(1,1,kparasta-1),((jx+2)*(jy+2)), &
                    MPI_DOUBLE_PRECISION,leftpem,taglr, &
                    MPI_COMM_WORLD,status,ierr)
            endif

            supres=1.e-30
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        base=8.*alalm(i,j,k)*alamm(i,j,k)
                        base=max(base,supres)
                        tem=1.5*2.*giac(i,j,k)**(2./3.)/base**.125
                        eps(i,j,k)=(dt/tem)/(1.+dt/tem)

                        base_rho=8.*alalmrho(i,j,k)*alammrho(i,j,k)
                        base_rho=max(base_rho,supres)
                        tem_rho=1.5*2.*giac(i,j,k)**(2./3.)/base_rho**.125
                        eps_rho(i,j,k)=(dt/tem_rho)/(1.+dt/tem_rho)


                        x1(i,j,k)=-uco(i,j,k)*dt/giac(i,j,k)
                        y1(i,j,k)=-vco(i,j,k)*dt/giac(i,j,k)
                        z1(i,j,k)=-wco(i,j,k)*dt/giac(i,j,k)
                    enddo
                enddo
            enddo

            call interp3(x1,y1,z1,vvf,kparasta,kparaend,alalm,alamm, &
                alalmrho,alammrho)
            !
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        appo1(i,j,k)=eps(i,j,k)*num(i,j,k) &
                            +(1.-eps(i,j,k))*vvf(i,j,k,1)
                        appo1(i,j,k)=.5*(appo1(i,j,k)+abs(appo1(i,j,k)))
                        appo2(i,j,k)=eps(i,j,k)*den(i,j,k) &
                            +(1.-eps(i,j,k))*vvf(i,j,k,2)

                        appo1rho(i,j,k)=eps_rho(i,j,k)*numrho(i,j,k) &
                            +(1.-eps_rho(i,j,k))*vvf(i,j,k,3)
                        appo1rho(i,j,k)=.5*(appo1rho(i,j,k)+abs(appo1rho(i,j,k)))
                        appo2rho(i,j,k)=eps_rho(i,j,k)*denrho(i,j,k) &
                            +(1.-eps_rho(i,j,k))*vvf(i,j,k,4)

                    !       appo1rho(i,j,k)=alalmrho(i,j,k)
                    !       appo2rho(i,j,k)=alamm(i,j,k)
                    enddo
                enddo
            enddo

        end if

    else if (nsgs==2) then ! DINAMICO

            !------------------------------------------------------------------
        !     mean in i and k
        !------------------------------------------------------------------
        if(medio_in_k .eq. 0)then

            do j=1,jy

                nummed_loc(j)=0.
                denmed_loc(j)=0.
                numrhomed_loc(j)=0.
                denrhomed_loc(j)=0.

                do k=kparasta,kparaend
                    do i=1,jx
                        nummed_loc(j)=nummed_loc(j)+num(i,j,k)
                        denmed_loc(j)=denmed_loc(j)+den(i,j,k)
                        numrhomed_loc(j)=numrhomed_loc(j)+numrho(i,j,k)
                        denrhomed_loc(j)=denrhomed_loc(j)+denrho(i,j,k)
                    enddo
                enddo

            enddo
            !
            ! now local sum and send the result to each proc
            ! for every plane
            !
            do m=1,40*(jx+2)*(jy+2)
                sbuffd(m)=0.
                rbuffd(m)=0.
            enddo

            call buffvect1d(sbuffd,nummed_loc,1)
            call buffvect1d(sbuffd,denmed_loc,2)
            call buffvect1d(sbuffd,numrhomed_loc,3)
            call buffvect1d(sbuffd,denrhomed_loc,4)

            call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),4*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            call buffvect2d(rbuffd,nummed,1)
            call buffvect2d(rbuffd,denmed,2)
            call buffvect2d(rbuffd,numrhomed,3)
            call buffvect2d(rbuffd,denrhomed,4)
            !
            ! finally the mean value, sum divided by the number of points
            ! on j plane (known by all procs)
            !
            rzer=0.
            do j=1,jy

                nummed(j)=nummed(j)/somma
                denmed(j)=denmed(j)/somma
                numrhomed(j)=numrhomed(j)/somma
                denrhomed(j)=denrhomed(j)/somma

                ! check on density denominator ---> see grg and grgrho
                !
                ! then the constant (already multiplied by -1)

                do k=kparasta,kparaend
                    do i=1,jx
                        appo1(i,j,k)=nummed(j)
                        appo2(i,j,k)=denmed(j)
                        appo1rho(i,j,k)=numrhomed(j)
                        appo2rho(i,j,k)=denrhomed(j)
                    enddo
                enddo

            enddo

        end if

        !------------------------------------------------------------------
        !     mean in k only
        !------------------------------------------------------------------
        if(medio_in_k .eq. 1)then

            do j=1,jy
                do i=1,jx

                    nummed_loc_z(i,j)=0.
                    denmed_loc_z(i,j)=0.
                    numrhomed_loc_z(i,j)=0.
                    denrhomed_loc_z(i,j)=0.

                    do k=kparasta,kparaend

                        nummed_loc_z(i,j)=    nummed_loc_z(i,j)+    num(i,j,k)
                        denmed_loc_z(i,j)=    denmed_loc_z(i,j)+    den(i,j,k)
                        numrhomed_loc_z(i,j)= numrhomed_loc_z(i,j)+ numrho(i,j,k)
                        denrhomed_loc_z(i,j)= denrhomed_loc_z(i,j)+ denrho(i,j,k)
                    enddo
                enddo
            end do

            !     per la media in i e k
            nummed_loc=0.
            denmed_loc=0.
            numrhomed_loc=0.
            denrhomed_loc=0.


            do m=1,40*(jx+2)*(jy+2)
                sbuffd(m)=0.
                rbuffd(m)=0.
            enddo

            do i=1,jx
                do j=1,jy
                    nummed_loc(j)   = nummed_loc_z(i,j)
                end do
                call buffvect1d(sbuffd,nummed_loc,i)
            end do

            do i=1,jx
                do j=1,jy
                    denmed_loc(j)   = denmed_loc_z(i,j)
                end do
                call buffvect1d(sbuffd,denmed_loc,jx+i)
            end do

            do i=1,jx
                do j=1,jy
                    numrhomed_loc(j)= numrhomed_loc_z(i,j)
                end do
                call buffvect1d(sbuffd,numrhomed_loc,jx*2+i)
            end do

            do i=1,jx
                do j=1,jy
                    denrhomed_loc(j)= denrhomed_loc_z(i,j)
                end do
                call buffvect1d(sbuffd,denrhomed_loc,jx*3+i)
            end do

            call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),4*jx*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            do i=1,jx
                call buffvect2d(rbuffd,nummed,i)
                do j=1,jy
                    nummed_z(i,j)   = nummed(j)
                end do
            end do

            do i=1,jx
                call buffvect2d(rbuffd,denmed,jx+i)
                do j=1,jy
                    denmed_z(i,j)   = denmed(j)
                end do
            end do

            do i=1,jx
                call buffvect2d(rbuffd,numrhomed,jx*2+i)
                do j=1,jy
                    numrhomed_z(i,j)   = numrhomed(j)
                end do
            end do

            do i=1,jx
                call buffvect2d(rbuffd,denrhomed,jx*3+i)
                do j=1,jy
                    denrhomed_z(i,j)   = denrhomed(j)
                end do
            end do

            linea_z = real(jz)

            rzer=0.
            do j=1,jy
                do i=1,jx

                    nummed_z(i,j) =    nummed_z(i,j)/linea_z !/float(jz)
                    denmed_z(i,j) =    denmed_z(i,j)/linea_z !/float(jz)
                    numrhomed_z(i,j) = numrhomed_z(i,j)/linea_z !/float(jz)
                    denrhomed_z(i,j) = denrhomed_z(i,j)/linea_z !/float(jz)

                    do k=kparasta,kparaend
                        appo1(i,j,k)=nummed_z(i,j)
                        appo2(i,j,k)=denmed_z(i,j)
                        appo1rho(i,j,k)=numrhomed_z(i,j)
                        appo2rho(i,j,k)=denrhomed_z(i,j)
                    enddo
                enddo
            enddo

        end if
    end if




    !
    !------------------------------------------------------------------
    !
    ! compute sgs (scale similar) cartesian stress momentum
    !
    if(inmod)then

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = sgs11(i,j,k)
                    pc2(i,j,k)  = sgs12(i,j,k)
                    pc3(i,j,k)  = sgs13(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil11(i,j,k) = pp1(i,j,k)
                    fil12(i,j,k) = pp2(i,j,k)
                    fil13(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = sgs21(i,j,k)
                    pc2(i,j,k)  = sgs22(i,j,k)
                    pc3(i,j,k)  = sgs23(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil21(i,j,k) = pp1(i,j,k)
                    fil22(i,j,k) = pp2(i,j,k)
                    fil23(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = sgs31(i,j,k)
                    pc2(i,j,k)  = sgs32(i,j,k)
                    pc3(i,j,k)  = sgs33(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil31(i,j,k) = pp1(i,j,k)
                    fil32(i,j,k) = pp2(i,j,k)
                    fil33(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do
        !
        ! average on homogeneous plane
        !
        do j=1,jy

            sus_loc11(j)=0.
            sus_loc12(j)=0.
            sus_loc13(j)=0.
            sus_loc22(j)=0.
            sus_loc23(j)=0.
            sus_loc33(j)=0.

            do k=kparasta,kparaend
                do i=1,jx
                    sus_loc11(j)=sus_loc11(j)+fil11(i,j,k)
                    sus_loc22(j)=sus_loc22(j)+fil22(i,j,k)
                    sus_loc33(j)=sus_loc33(j)+fil33(i,j,k)
                    sus_loc12(j)=sus_loc12(j)+.5*(fil12(i,j,k)+fil21(i,j,k))
                    sus_loc13(j)=sus_loc13(j)+.5*(fil13(i,j,k)+fil31(i,j,k))
                    sus_loc23(j)=sus_loc23(j)+.5*(fil23(i,j,k)+fil32(i,j,k))
                enddo
            enddo

        enddo


        if (debugg.eq.1) then
            do j=1,n2
                write(9020+myid+ttt,*)sus_loc11(j),sus_loc22(j),sus_loc33(j)
                write(9025+myid+ttt,*)sus_loc12(j),sus_loc13(j),sus_loc23(j)
            end do

            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
                        write(5345+myid+ttt,*)fil11(i,j,k),fil22(i,j,k),fil33(i,j,k), &
                            fil13(i,j,k),fil31(i,j,k)
                        write(5350+myid+ttt,*)fil12(i,j,k),fil21(i,j,k),fil32(i,j,k), &
                            fil23(i,j,k)
                    end do
                end do
            end do

        end if
        !
        ! sum over local sum and send the result to all procs
        ! for each plane
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuffd(m)=0.
            rbuffd(m)=0.
        enddo

        call buffvect1d(sbuffd,sus_loc11,1)
        call buffvect1d(sbuffd,sus_loc12,2)
        call buffvect1d(sbuffd,sus_loc13,3)
        call buffvect1d(sbuffd,sus_loc22,4)
        call buffvect1d(sbuffd,sus_loc23,5)
        call buffvect1d(sbuffd,sus_loc33,6)

        call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),6*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        call buffvect2d(rbuffd,sus11,1)
        call buffvect2d(rbuffd,sus12,2)
        call buffvect2d(rbuffd,sus13,3)
        call buffvect2d(rbuffd,sus22,4)
        call buffvect2d(rbuffd,sus23,5)
        call buffvect2d(rbuffd,sus33,6)
        !
        ! finally the mean value, the sum divided by the number of points
        ! on each plane j (known by all procs
        !
        do j=1,jy
            sus11(j)=cb(j)*sus11(j)/somma
            sus22(j)=cb(j)*sus22(j)/somma
            sus33(j)=cb(j)*sus33(j)/somma
            sus12(j)=cb(j)*sus12(j)/somma
            sus13(j)=cb(j)*sus13(j)/somma
            sus23(j)=cb(j)*sus23(j)/somma
        enddo
        !
        ! compute dissipation scale similar tau_ij*s_ij
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    sgs11(i,j,k)=s11(i,j,k)
                    sgs12(i,j,k)=s12(i,j,k)
                    sgs13(i,j,k)=s13(i,j,k)
                    sgs21(i,j,k)=s21(i,j,k)
                    sgs22(i,j,k)=s22(i,j,k)
                    sgs23(i,j,k)=s23(i,j,k)
                    sgs31(i,j,k)=s31(i,j,k)
                    sgs32(i,j,k)=s32(i,j,k)
                    sgs33(i,j,k)=s33(i,j,k)
                enddo
            enddo
        enddo
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = sgs11(i,j,k)
                    pc2(i,j,k)  = sgs12(i,j,k)
                    pc3(i,j,k)  = sgs13(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fipp11(i,j,k) = pp1(i,j,k)
                    fipp12(i,j,k) = pp2(i,j,k)
                    fipp13(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = sgs21(i,j,k)
                    pc2(i,j,k)  = sgs22(i,j,k)
                    pc3(i,j,k)  = sgs23(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fipp21(i,j,k) = pp1(i,j,k)
                    fipp22(i,j,k) = pp2(i,j,k)
                    fipp23(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = sgs31(i,j,k)
                    pc2(i,j,k)  = sgs32(i,j,k)
                    pc3(i,j,k)  = sgs33(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fipp31(i,j,k) = pp1(i,j,k)
                    fipp32(i,j,k) = pp2(i,j,k)
                    fipp33(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        do j=1,jy
            sus_loc(j)=0.

            do k=kparasta,kparaend
                do i=1,jx
                    sus_loc(j)=sus_loc(j)+ &
                        fil11(i,j,k)*fipp11(i,j,k)+ &
                        fil12(i,j,k)*fipp12(i,j,k)+ &
                        fil13(i,j,k)*fipp13(i,j,k)+ &
                        fil21(i,j,k)*fipp21(i,j,k)+ &
                        fil22(i,j,k)*fipp22(i,j,k)+ &
                        fil23(i,j,k)*fipp23(i,j,k)+ &
                        fil31(i,j,k)*fipp31(i,j,k)+ &
                        fil32(i,j,k)*fipp32(i,j,k)+ &
                        fil33(i,j,k)*fipp33(i,j,k)
                enddo
            enddo

        enddo

        if (debugg.eq.1) then
            do j=1,jy
                write(5355+myid+ttt,*)sus_loc(j)
            end do
        end if

        !
        ! now sum of local sum and send the result to all procs
        ! for each plane
        !
        call MPI_ALLREDUCE(sus_loc(1),sus(1),jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        ! finally the mean value, the sum divided by the number of points
        ! for each plane (known by all the procs)
        !
        do j=1,jy
            sus(j)=cb(j)*sus(j)/somma
        enddo

        !
        ! compute sgs (scale similar) cartesian stress density
        !

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    pc1(i,j,k)  = sgsrho11(i,j,k)
                    pc2(i,j,k)  = sgsrho22(i,j,k)
                    pc3(i,j,k)  = sgsrho33(i,j,k)
                    ap11(i,j,k) = apcsx(i,j,k)
                    ap12(i,j,k) = apcsy(i,j,k)
                    ap13(i,j,k) = apcsz(i,j,k)
                    ap21(i,j,k) = apetx(i,j,k)
                    ap22(i,j,k) = apety(i,j,k)
                    ap23(i,j,k) = apetz(i,j,k)
                    ap31(i,j,k) = apztx(i,j,k)
                    ap32(i,j,k) = apzty(i,j,k)
                    ap33(i,j,k) = apztz(i,j,k)
                end do
            end do
        end do

        call inverse_para2(myid,nproc,kparasta,kparaend)

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    fil11(i,j,k) = pp1(i,j,k)
                    fil22(i,j,k) = pp2(i,j,k)
                    fil33(i,j,k) = pp3(i,j,k)
                end do
            end do
        end do

        !
        ! average oo homogeneous planes
        !
        do j=1,jy
            susrho_loc11(j)=0.
            susrho_loc22(j)=0.
            susrho_loc33(j)=0.

            do k=kparasta,kparaend
                do i=1,jx
                    susrho_loc11(j)=susrho_loc11(j)+fil11(i,j,k)
                    susrho_loc22(j)=susrho_loc22(j)+fil22(i,j,k)
                    susrho_loc33(j)=susrho_loc33(j)+fil33(i,j,k)
                enddo
            enddo

        enddo
        !
        ! now sum over local sum and send the results to all procs
        ! for each plane
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuffd(m)=0.
            rbuffd(m)=0.
        enddo

        call buffvect1d(sbuffd,susrho_loc11,1)
        call buffvect1d(sbuffd,susrho_loc22,2)
        call buffvect1d(sbuffd,susrho_loc33,3)

        call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),3*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        call buffvect2d(rbuffd,susrho11,1)
        call buffvect2d(rbuffd,susrho22,2)
        call buffvect2d(rbuffd,susrho33,3)

        !
        ! finally the mean value, the sum divide by the number of points
        ! for each j-plane (known by all procs)
        !
        do j=1,jy
            susrho11(j)=cbrho(j)*susrho11(j)/somma
            susrho22(j)=cbrho(j)*susrho22(j)/somma
            susrho33(j)=cbrho(j)*susrho33(j)/somma
        enddo

    end if  !if on inmod .eq. 1
    !
    !-----------------------------------------------------------------------
    !
    ! controvariant stress tau_ik
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                !     check on denominator
                if(appo2(i,j,k).eq.0.0)then
                    grg=0.
                else
                    grg=-appo1(i,j,k)/appo2(i,j,k)*4.*giac(i,j,k)**(2./3.) &
                        *smod(i,j,k)
                endif
                if(appo2rho(i,j,k).eq.0.0)then
                    grgrho=0.
                else
                    grgrho=-appo1rho(i,j,k)/appo2rho(i,j,k) &
                        *4.*giac(i,j,k)**(2./3.)*smod(i,j,k)
                endif
                !     end check on denominator

                sgs11(i,j,k)=grg*s11(i,j,k)
                sgs12(i,j,k)=grg*s12(i,j,k)
                sgs13(i,j,k)=grg*s13(i,j,k)
                sgs21(i,j,k)=grg*s21(i,j,k)
                sgs22(i,j,k)=grg*s22(i,j,k)
                sgs23(i,j,k)=grg*s23(i,j,k)
                sgs31(i,j,k)=grg*s31(i,j,k)
                sgs32(i,j,k)=grg*s32(i,j,k)
                sgs33(i,j,k)=grg*s33(i,j,k)
                sgsrho11(i,j,k)=grgrho*rho11(1,i,j,k)
                sgsrho22(i,j,k)=grgrho*rho22(1,i,j,k)
                sgsrho33(i,j,k)=grgrho*rho33(1,i,j,k)

            enddo
        enddo
    enddo
    !
    ! cartesian stress tau_ij
    !

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = sgs11(i,j,k)
                pc2(i,j,k)  = sgs12(i,j,k)
                pc3(i,j,k)  = sgs13(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil11(i,j,k) = pp1(i,j,k)
                fil12(i,j,k) = pp2(i,j,k)
                fil13(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = sgs21(i,j,k)
                pc2(i,j,k)  = sgs22(i,j,k)
                pc3(i,j,k)  = sgs23(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil21(i,j,k) = pp1(i,j,k)
                fil22(i,j,k) = pp2(i,j,k)
                fil23(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = sgs31(i,j,k)
                pc2(i,j,k)  = sgs32(i,j,k)
                pc3(i,j,k)  = sgs33(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil31(i,j,k) = pp1(i,j,k)
                fil32(i,j,k) = pp2(i,j,k)
                fil33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    if (debugg.eq.1) then

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(5360+myid+ttt,*)fil11(i,j,k),fil12(i,j,k),fil13(i,j,k)
                    write(5365+myid+ttt,*)fil21(i,j,k),fil22(i,j,k),fil23(i,j,k)
                    write(5370+myid+ttt,*)fil31(i,j,k),fil32(i,j,k),fil33(i,j,k)
                end do
            end do
        end do

    end if
    !
    ! average on homogeneous plane
    !
    do j=1,jy

        sub11_loc(j)=0.
        sub12_loc(j)=0.
        sub13_loc(j)=0.
        sub22_loc(j)=0.
        sub23_loc(j)=0.
        sub33_loc(j)=0.

        do k=kparasta,kparaend
            do i=1,jx

                sub11_loc(j)=sub11_loc(j)+fil11(i,j,k)
                sub22_loc(j)=sub22_loc(j)+fil22(i,j,k)
                sub33_loc(j)=sub33_loc(j)+fil33(i,j,k)
                sub12_loc(j)=sub12_loc(j)+0.5*(fil12(i,j,k)+fil21(i,j,k))
                sub13_loc(j)=sub13_loc(j)+0.5*(fil13(i,j,k)+fil31(i,j,k))
                sub23_loc(j)=sub23_loc(j)+0.5*(fil23(i,j,k)+fil32(i,j,k))

            enddo
        enddo

    enddo
    !
    ! now sum of local sum and the result is send to all procs
    ! for each plane
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuffd(m)=0.
        rbuffd(m)=0.
    enddo

    call buffvect1d(sbuffd,sub11_loc,1)
    call buffvect1d(sbuffd,sub22_loc,2)
    call buffvect1d(sbuffd,sub33_loc,3)
    call buffvect1d(sbuffd,sub12_loc,4)
    call buffvect1d(sbuffd,sub13_loc,5)
    call buffvect1d(sbuffd,sub23_loc,6)

    call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),6*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call buffvect2d(rbuffd,sub11,1)
    call buffvect2d(rbuffd,sub22,2)
    call buffvect2d(rbuffd,sub33,3)
    call buffvect2d(rbuffd,sub12,4)
    call buffvect2d(rbuffd,sub13,5)
    call buffvect2d(rbuffd,sub23,6)
    !
    ! finally the mean value, the sum divided by the number
    ! of points on each j-plane (known by all procs)
    !
    do j=1,jy

        sub11(j)=sub11(j)/somma
        sub22(j)=sub22(j)/somma
        sub33(j)=sub33(j)/somma
        sub12(j)=sub12(j)/somma
        sub13(j)=sub13(j)/somma
        sub23(j)=sub23(j)/somma

    enddo
    !
    ! now computee the eddy viscosity dissipation defined
    ! as the contraction of tau_ij*s_ij, I need to derive
    ! the cartesian component from the controvariant one
    !
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                sgs11(i,j,k)=s11(i,j,k)
                sgs12(i,j,k)=s12(i,j,k)
                sgs13(i,j,k)=s13(i,j,k)
                sgs21(i,j,k)=s21(i,j,k)
                sgs22(i,j,k)=s22(i,j,k)
                sgs23(i,j,k)=s23(i,j,k)
                sgs31(i,j,k)=s31(i,j,k)
                sgs32(i,j,k)=s32(i,j,k)
                sgs33(i,j,k)=s33(i,j,k)
            end do
        end do
    end do
    !
    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = sgs11(i,j,k)
                pc2(i,j,k)  = sgs12(i,j,k)
                pc3(i,j,k)  = sgs13(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fipp11(i,j,k) = pp1(i,j,k)
                fipp12(i,j,k) = pp2(i,j,k)
                fipp13(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = sgs21(i,j,k)
                pc2(i,j,k)  = sgs22(i,j,k)
                pc3(i,j,k)  = sgs23(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fipp21(i,j,k) = pp1(i,j,k)
                fipp22(i,j,k) = pp2(i,j,k)
                fipp23(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = sgs31(i,j,k)
                pc2(i,j,k)  = sgs32(i,j,k)
                pc3(i,j,k)  = sgs33(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fipp31(i,j,k) = pp1(i,j,k)
                fipp32(i,j,k) = pp2(i,j,k)
                fipp33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    if (debugg.eq.1) then

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(5375+myid+ttt,*)fipp11(i,j,k),fipp12(i,j,k),fipp13(i,j,k)
                    write(5380+myid+ttt,*)fipp21(i,j,k),fipp22(i,j,k),fipp23(i,j,k)
                    write(5385+myid+ttt,*)fipp31(i,j,k),fipp32(i,j,k),fipp33(i,j,k)
                end do
            end do
        end do

    end if

    do j=1,jy
        sub_loc(j)=0.

        do k=kparasta,kparaend
            do i=1,jx
                sub_loc(j)=sub_loc(j)+fil11(i,j,k)*fipp11(i,j,k)+ &
                    fil12(i,j,k)*fipp12(i,j,k)+ &
                    fil13(i,j,k)*fipp13(i,j,k)+ &
                    fil21(i,j,k)*fipp21(i,j,k)+ &
                    fil22(i,j,k)*fipp22(i,j,k)+ &
                    fil23(i,j,k)*fipp23(i,j,k)+ &
                    fil31(i,j,k)*fipp31(i,j,k)+ &
                    fil32(i,j,k)*fipp32(i,j,k)+ &
                    fil33(i,j,k)*fipp33(i,j,k)
            enddo
        enddo

    enddo
    !
    ! now sum over local sum and send the result to all procs
    ! for  each plane
    !
    call MPI_ALLREDUCE(sub_loc(1),sub(1),jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !
    ! finally the mean value, the sum divided by the number
    ! of point for each plane (known by all prc)
    !
    do j=1,jy
        sub(j)=sub(j)/somma
    enddo
    !
    ! to conclude, computation of sgs cartesian
    ! eddy for density
    !
    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                pc1(i,j,k)  = sgsrho11(i,j,k)
                pc2(i,j,k)  = sgsrho22(i,j,k)
                pc3(i,j,k)  = sgsrho33(i,j,k)
                ap11(i,j,k) = apcsx(i,j,k)
                ap12(i,j,k) = apcsy(i,j,k)
                ap13(i,j,k) = apcsz(i,j,k)
                ap21(i,j,k) = apetx(i,j,k)
                ap22(i,j,k) = apety(i,j,k)
                ap23(i,j,k) = apetz(i,j,k)
                ap31(i,j,k) = apztx(i,j,k)
                ap32(i,j,k) = apzty(i,j,k)
                ap33(i,j,k) = apztz(i,j,k)
            end do
        end do
    end do

    call inverse_para2(myid,nproc,kparasta,kparaend)

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1
                fil11(i,j,k) = pp1(i,j,k)
                fil22(i,j,k) = pp2(i,j,k)
                fil33(i,j,k) = pp3(i,j,k)
            end do
        end do
    end do

    if (debugg.eq.1) then
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(5390+myid+ttt,*)fil11(i,j,k),fil22(i,j,k),fil33(i,j,k)
                end do
            end do
        end do
    end if
    !
    do j=1,jy
        subrho_loc11(j)=0.
        subrho_loc22(j)=0.
        subrho_loc33(j)=0.

        do k=kparasta,kparaend
            do i=1,jx
                subrho_loc11(j)=subrho_loc11(j)+fil11(i,j,k)
                subrho_loc22(j)=subrho_loc22(j)+fil22(i,j,k)
                subrho_loc33(j)=subrho_loc33(j)+fil33(i,j,k)
            enddo
        enddo

    enddo
    !
    ! now sum over local sum and send the result to all procs
    ! for each plane
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuffd(m)=0.
        rbuffd(m)=0.
    enddo

    call buffvect1d(sbuffd,subrho_loc11,1)
    call buffvect1d(sbuffd,subrho_loc22,2)
    call buffvect1d(sbuffd,subrho_loc33,3)

    call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),3*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call buffvect2d(rbuffd,subrho11,1)
    call buffvect2d(rbuffd,subrho22,2)
    call buffvect2d(rbuffd,subrho33,3)
    !
    ! finally the mean value, the sum divided by the number
    ! of points for each j-plane (known from all procs)
    !
    do j=1,jy
        subrho11(j)=subrho11(j)/somma
        subrho22(j)=subrho22(j)/somma
        subrho33(j)=subrho33(j)/somma
    enddo
    !
    !-----------------------------------------------------------------------
    !***********************************************************************
    !-----------------------------------------------------------------------
    ! eddy viscosity and diffusivity computation

    if (myid.eq.0) then
        kparastal=0
        kparaendl=kparaend
    else if (myid.eq.nproc-1) then
        kparastal=kparasta
        kparaendl=jz+1
    else
        kparastal=kparasta
        kparaendl=kparaend
    endif

    !
    if (ktime.gt.0) then
        !
        !     extrapolation for appo1 and appo2 on j=1 and j=jy
        !
        do k=kparasta,kparaend
            do i=1,jx
                appo1(i,1,k)=0.
                appo1rho(i,1,k)=0.
                appo1(i,jy,k)=0.
                appo1rho(i,jy,k)=0.
                !
                appo2(i,1,k)=appo2(i,2,k)                      ! move to wall
                appo2(i,jy,k)=appo2(i,jy-1,k)                  ! move to wall
                appo2rho(i,1,k)=appo2rho(i,2,k)                ! move to wall
                appo2rho(i,jy,k)=appo2rho(i,jy-1,k)            ! move to wall
            enddo
        enddo
        !
        rzer=0.

        !-----------------------------------------------------------------------
        if(isotropo)then
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx

                        !        check on denominator
                        if(abs(appo2(i,j,k)).lt. 1.d-10)then
                            cost=0.
                        else
                            cost=.5*appo1(i,j,k)/appo2(i,j,k)
                        endif
                        if(abs(appo2rho(i,j,k)).lt. 1.d-10)then
                            costrho=0.
                        else
                            costrho=.5*appo1rho(i,j,k)/appo2rho(i,j,k)
                        endif
                        !        end check on denominator
                        !
                        cost=max(cost,rzer)
                        costrho=max(costrho,rzer)

                        !        pay attention, consider that:
                        !          smagorinsky cost**2
                        !          dynamic s.  cost

                        annit(i,j,k)=annit(i,j,k) &
                            +giac(i,j,k)**(2./3.)*4.*cost*smod(i,j,k)
                        annitV(i,j,k)=annit(i,j,k)

                        !        eddy diffusivity
                        if(re_analogy)then
                            do isc=1,nscal
                                akapt(isc,i,j,k)  = annit(i,j,k)/prsc(isc)
                                akaptV(isc,i,j,k) = annitV(i,j,k)/prsc(isc)
                            end do
                        else
                            do isc=1,nscal
                                akapt(isc,i,j,k)=akapt(isc,i,j,k) &
                                    +giac(i,j,k)**(2./3.)*4.*costrho*smod(i,j,k)
                                akaptV(isc,i,j,k)=akapt(isc,i,j,k)
                            end do
                        end if

                    enddo
                enddo
            enddo
        end if

        !-----------------------------------------------------------------------

        if(.not.isotropo)then
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx

                        !        check on denominator
                        if(abs(appo2(i,j,k)).lt. 1.d-10)then
                            cost=0.
                        else
                            cost=.5*appo1(i,j,k)/appo2(i,j,k)
                        endif
                        if(abs(appo2rho(i,j,k)).lt. 1.d-10)then
                            costrho=0.
                        else
                            costrho=.5*appo1rho(i,j,k)/appo2rho(i,j,k)
                        endif
                        !        end check on denominator
                        !
                        cost=max(cost,rzer)
                        costrho=max(costrho,rzer)

                        !       ------------- anisotropic model --------------------------------

                        lhorx=.25*(x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
                            +x(i  ,j-1,k-1)+x(i  ,j-1,k  )) &
                            -.25*(x(i-1,j  ,k  )+x(i-1,j  ,k-1) &
                            +x(i-1,j-1,k-1)+x(i-1,j-1,k  ))

                        lhorz=.25*(z(i  ,j  ,k  )+z(i  ,j-1,k  ) &
                            +z(i-1,j-1,k  )+z(i-1,j  ,k  )) &
                            -.25*(z(i  ,j  ,k-1)+z(i  ,j-1,k-1) &
                            +z(i-1,j-1,k-1)+z(i-1,j  ,k-1))

                        lhor = lhorx*lhorx + lhorz*lhorz !length square


                        lver =.25*(y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
                            +y(i-1,j  ,k-1)+y(i-1,j  ,k  )) &
                            -.25*(y(i  ,j-1,k  )+y(i  ,j-1,k-1) &
                            +y(i-1,j-1,k-1)+y(i-1,j-1,k  ))

                        lver = lver*lver

                        annit(i,j,k) =annit(i,j,k) + lhor*costH*smodH(i,j,k)

                        annitV(i,j,k)=annitV(i,j,k) + lver*costV*smodV(i,j,k)

                        !        ------------------------------------------------
                        !        eddy diffusivity
                        if(re_analogy)then
                            do isc=1,nscal
                                akapt(isc,i,j,k)  = prsc(isc)*annit(i,j,k)
                                akaptV(isc,i,j,k) = prsc(isc)*annitV(i,j,k)
                            end do
                        else
                            do isc=1,nscal
                                akapt(isc,i,j,k)=akapt(isc,i,j,k)+ &
                                    giac(i,j,k)**(2./3.)*4.*costrho*smod(i,j,k)
                                akaptV(isc,i,j,k)=akaptV(isc,i,j,k)+ &
                                    giac(i,j,k)**(2./3.)*4.*costrho*smod(i,j,k)

                            end do
                        end if
                    !
                    enddo
                enddo
            enddo
        end if

        !-----------------------------------------------------------------------
        !
        ! periodicity for appo1 and appo2 (imposed over x and z)
        !
        ! generalized periodicity on x
        !
        do k=kparasta,kparaend
            do j=1,jy
                appo1(0,j,k)=real(1-ip)*appo1(jx,j,k)+real(ip)*appo1(1,j,k)
                appo2(0,j,k)=real(1-ip)*appo2(jx,j,k)+real(ip)*appo2(1,j,k)
                appo1rho(0,j,k)=real(1-ip)*appo1rho(jx,j,k)+real(ip)*appo1rho(1,j,k)
                appo2rho(0,j,k)=real(1-ip)*appo2rho(jx,j,k)+real(ip)*appo2rho(1,j,k)
                !
                appo1(jx+1,j,k)=real(1-ip)*appo1(1,j,k)+real(ip)*appo1(jx,j,k)
                appo2(jx+1,j,k)=real(1-ip)*appo2(1,j,k)+real(ip)*appo2(jx,j,k)
                appo1rho(jx+1,j,k)=real(1-ip)*appo1rho(1,j,k)+real(ip)*appo1rho(jx,j,k)
                appo2rho(jx+1,j,k)=real(1-ip)*appo2rho(1,j,k)+real(ip)*appo2rho(jx,j,k)
            enddo
        enddo
        !
        !
        ! generalized periodicity on z
        !
        do m=1,40*(jx+2)*(jy+2)
            sbuffd(m)=0.
            rbuffd(m)=0.
        enddo

        if (myid.eq.nproc-1) then

            call buffer1d_par(sbuffd,appo1   ,1,jz,kparasta,kparaend)
            call buffer1d_par(sbuffd,appo2   ,2,jz,kparasta,kparaend)
            call buffer1d_par(sbuffd,appo1rho,3,jz,kparasta,kparaend)
            call buffer1d_par(sbuffd,appo2rho,4,jz,kparasta,kparaend)

        else if (myid.eq.0) then

            call buffer1d_par(sbuffd,appo1   ,1,1,kparasta,kparaend)
            call buffer1d_par(sbuffd,appo2   ,2,1,kparasta,kparaend)
            call buffer1d_par(sbuffd,appo1rho,3,1,kparasta,kparaend)
            call buffer1d_par(sbuffd,appo2rho,4,1,kparasta,kparaend)

        endif
        !
        ! now I exchange, so that P0 knows k=jz and Pn-1 knows k=1
        !
        if (myid.eq.nproc-1) then

            call MPI_SENDRECV(sbuffd(1),4*(jx+2)*(jy+2), &
                MPI_DOUBLE_PRECISION,0,5001, &
                rbuffd(1),4*(jx+2)*(jy+2), &
                MPI_DOUBLE_PRECISION,0,4001, &
                MPI_COMM_WORLD,status,ierr)

        else if (myid.eq.0) then

            call MPI_SENDRECV(sbuffd(1),4*(jx+2)*(jy+2), &
                MPI_DOUBLE_PRECISION,nproc-1,4001, &
                rbuffd(1),4*(jx+2)*(jy+2), &
                MPI_DOUBLE_PRECISION,nproc-1,5001, &
                MPI_COMM_WORLD,status,ierr)

        endif

        if (myid.eq.0) then

            call buffer2d_par(rbuffd,appo1_piano   ,1,jz)
            call buffer2d_par(rbuffd,appo2_piano   ,2,jz)
            call buffer2d_par(rbuffd,appo1rho_piano,3,jz)
            call buffer2d_par(rbuffd,appo2rho_piano,4,jz)

        else if (myid.eq.nproc-1) then

            call buffer2d_par(rbuffd,appo1_piano   ,1,1)
            call buffer2d_par(rbuffd,appo2_piano   ,2,1)
            call buffer2d_par(rbuffd,appo1rho_piano,3,1)
            call buffer2d_par(rbuffd,appo2rho_piano,4,1)

        endif
        !
        ! now P0 knows the plane k=jz
        ! and Pn-1 knows the plane k=1
        !
        if(myid.eq.0)then

            do j=1,jy
                do i=1,jx
                    appo1(i,j,0)=real(1-kp)*   appo1_piano(i,j,jz)+   real(kp)*appo1(i,j,1)
                    appo2(i,j,0)=real(1-kp)*   appo2_piano(i,j,jz)+   real(kp)*appo2(i,j,1)
                    appo1rho(i,j,0)=real(1-kp)*appo1rho_piano(i,j,jz)+real(kp)*appo1rho(i,j,1)
                    appo2rho(i,j,0)=real(1-kp)*appo2rho_piano(i,j,jz)+real(kp)*appo2rho(i,j,1)
                enddo
            enddo

        else if(myid.eq.nproc-1)then

            do j=1,jy
                do i=1,jx
                    appo1(i,j,jz+1)=real(1-kp)*appo1_piano(i,j,1)+real(kp)*appo1(i,j,jz)
                    appo2(i,j,jz+1)=real(1-kp)*appo2_piano(i,j,1)+real(kp)*appo2(i,j,jz)
                    appo1rho(i,j,jz+1)=real(1-kp)*appo1rho_piano(i,j,1) &
                        +real(kp)*appo1rho(i,j,jz)
                    appo2rho(i,j,jz+1)=real(1-kp)*appo2rho_piano(i,j,1) &
                        +real(kp)*appo2rho(i,j,jz)
                enddo
            enddo

        endif
        !
        ! put appo1 and 2 in alalm and alamm
        !
        if(myid.eq.0)then
            kparastal=0
            kparaendl=kparaend
        else if (myid.eq.nproc-1) then
            kparastal=kparasta
            kparaendl=jz+1
        else
            kparastal=kparasta
            kparaendl=kparaend
        endif

        do k=kparastal,kparaendl
            do j=1,jy
                do i=0,jx+1
                    alalm(i,j,k)=appo1(i,j,k)
                    alamm(i,j,k)=appo2(i,j,k)
                    alalmrho(i,j,k)=appo1rho(i,j,k)
                    alammrho(i,j,k)=appo2rho(i,j,k)
                enddo
            enddo
        enddo
    !
    end if ! close if ktime
    !
    !
    ! compute alalm and alamm on ghost cell bottom and upper
    ! if periodic ????
    !
    do k=kparastal,kparaendl
        do i=0,jx+1
            alalm(i,0,k)=0.
            alamm(i,0,k)=alamm(i,1,k)
            alalmrho(i,0,k)=0.
            alammrho(i,0,k)=alammrho(i,1,k)
            !
            alalm(i,jy+1,k)=0.
            alamm(i,jy+1,k)=alamm(i,jy,k)
            alalmrho(i,jy+1,k)=0.
            alammrho(i,jy+1,k)=alammrho(i,jy,k)
        enddo
    enddo
    !
    ! add eddy viscosity and diffusivity on periodic boundary
    ! left and right side
    !
    do k=kparasta,kparaend
        do j=1,jy
            annit(0,j,k)   =real(1-ip)*annit(jx,j,k)+real(ip)*annit(1,j,k)
            annit(jx+1,j,k)=real(1-ip)*annit(1,j,k)+real(ip)*annit(jx,j,k)
            !
            do isc=1,nscal
                akapt(isc,0,j,k)   =real(1-ip)*akapt(isc,jx,j,k) &
                    +real(ip)*akapt(isc,1,j,k)
                akapt(isc,jx+1,j,k)=real(1-ip)*akapt(isc,1,j,k) &
                    +real(ip)*akapt(isc,jx,j,k)
            end do
        end do
    end do
    !
    ! side front and back
    !
    do m=1,40*(jx+2)*(jy+2)
        sbuff(m)=0.
        rbuff(m)=0.
    enddo

    if (myid.eq.nproc-1) then

        call buffer1old_par(sbuff,annit,1,jz)
        cont = 1
        do isc=1,nscal
            cont = cont+1
            call buffer1old_par_nscal(sbuff,akapt,cont,jz,isc)
        end do

    else if (myid.eq.0) then

        call buffer1old_par(sbuff,annit,1,1)
        cont = 1
        do isc=1,nscal
            cont = cont + 1
            call buffer1old_par_nscal(sbuff,akapt,cont,1,isc)
        end do

    endif

    cont = 1+nscal
    !
    ! now the exchange so that P0 knows k=jz
    ! and Pn-1 knows k=1
    !
    if (myid.eq.nproc-1) then

        call MPI_SENDRECV(sbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,0,3001, &
            rbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,0,2001, &
            MPI_COMM_WORLD,status,ierr)

    else if (myid.eq.0) then

        call MPI_SENDRECV( &
            sbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,2001, &
            rbuff(1),cont*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,3001, &
            MPI_COMM_WORLD,status,ierr)

    endif

    if (myid.eq.0) then

        call buffer2old_par(rbuff,annit_piano,1,jz)
        cont=1
        do isc=1,nscal
            cont = cont + 1
            call buffer2old_par_nscal(rbuff,akapt_piano,cont,jz,isc)
        end do

    else if (myid.eq.nproc-1) then

        call buffer2old_par(rbuff,annit_piano,1,1)
        cont = 1
        do isc=1,nscal
            cont = cont + 1
            call buffer2old_par_nscal(rbuff,akapt_piano,cont,1,isc)
        end do

    endif
    !
    ! now P0 knows the plane k=jz
    ! and Pn-1 knows the plane k=1
    !
    if(myid.eq.0)then

        do j=1,jy
            do i=1,jx
                annit(i,j,0)=real(1-kp)*annit_piano(i,j,jz)+real(kp)*annit(i,j,1)
                do isc=1,nscal
                    akapt(isc,i,j,0)=real(1-kp)*akapt_piano(isc,i,j,jz) &
                        +real(kp)*akapt(isc,i,j,1)
                end do
            end do
        end do

    else if(myid.eq.nproc-1)then

        do j=1,jy
            do i=1,jx
                annit(i,j,jz+1)=real(1-kp)*annit_piano(i,j,1)+real(kp)*annit(i,j,jz)
                do isc=1,nscal
                    akapt(isc,i,j,jz+1)=real(1-kp)*akapt_piano(isc,i,j,1) &
                        +real(kp)*akapt(isc,i,j,jz)
                end do
            end do
        end do

    endif
    !
    ! pass the closer plane
    !
    if(leftpem /= MPI_PROC_NULL) then
        call MPI_SEND(annit(0,0,kparasta),(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpem,tagls, &
            MPI_COMM_WORLD,ierr)
    endif
    if(rightpem /= MPI_PROC_NULL) then
        call MPI_RECV(annit(0,0,kparaend+1),(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpem,tagrr, &
            MPI_COMM_WORLD,status,ierr)
    endif
    if(rightpem /= MPI_PROC_NULL) then
        call MPI_SEND(annit(0,0,kparaend),(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpem,tagrs, &
            MPI_COMM_WORLD,ierr)
    endif
    if(leftpem /= MPI_PROC_NULL) then
        call MPI_RECV(annit(0,0,kparasta-1),(jx+2)*(jy+2), &
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

    !     pass akapt for nscal
    allocate (P_akapt(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))

    do isc=1,nscal

        do k=kparasta-deepl,kparasta+deepr
            do j=0,n2+1
                do i=0,n1+1
                    P_akapt(i,j,k) = akapt(isc,i,j,k)
                end do
            end do
        end do

        if(leftpem /= MPI_PROC_NULL) then
            call MPI_SEND(P_akapt(0,0,kparasta),(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,tagls, &
                MPI_COMM_WORLD,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(P_akapt(0,0,kparaend+1),(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrr, &
                MPI_COMM_WORLD,status,ierr)
        endif
        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SEND(P_akapt(0,0,kparaend),(jx+2)*(jy+2), &
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

        do k=kparasta-deepl,kparasta+deepr
            do j=0,n2+1
                do i=0,n1+1
                    akapt(isc,i,j,k) =  P_akapt(i,j,k)
                end do
            end do
        end do

    end do !cycle on nscal

    deallocate(P_akapt)




    if (debugg.eq.2) then


        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1
                    write(95+myid+ttt,*)uco(i,j,k),vco(i,j,k),wco(i,j,k)
                    write(100+myid+ttt,*)u(i,j,k),v(i,j,k),w(i,j,k)
                    write(105+myid+ttt,*)rho(i,j,k)
                    write(110+myid+ttt,*)annit(i,j,k),akapt(1,i,j,k)
                    write(115+myid+ttt,*)subrho11(j),subrho22(j),subrho33(j)
                    write(120+myid+ttt,*)susrho11(j),susrho22(j),susrho33(j)
                end do
            end do
        end do

        !
        do j=1,n2
            write(125+myid+ttt,*)sus_loc(j)
            write(130+myid+ttt,*)sus_loc11(j),sus_loc12(j),sus_loc13(j)
            write(135+myid+ttt,*)sus_loc22(j),sus_loc23(j),sus_loc33(j)
            write(140+myid+ttt,*)sub_loc(j)
            write(145+myid+ttt,*)subrho_loc11(j),subrho_loc22(j), &
                subrho_loc33(j)
        end do


    end if


    deallocate(rho)

    if (nsgs==3) then
        deallocate(x1)
        deallocate(y1)
        deallocate(z1)
        deallocate(eps)
        deallocate(eps_rho)
        deallocate(vvf)
    end if

    return
end subroutine turbo_lagrdin


!***********************************************************************
subroutine interp3(xx1,yy1,zz1,vvf,kparasta,kparaend,alalm,alamm,alalmrho,alammrho)
   !***********************************************************************
   ! interpolazione delle funzioni num e den per modello
   ! lagrangiano in spazio computazionale
   !
   use myarrays_metri3
   use scala3

   implicit none
   !-----------------------------------------------------------------------
   !     variable declaration
   integer i,j,k,inp,jnp,knp,ist,jst,kst,ii
   integer kparasta,kparaend

   real alalm(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   real alamm(0:n1+1,0:n2+1,kparasta-1:kparaend+1)

   real alalmrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   real alammrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   !
   real dx1,dx2,dy1,dy2,dz1,dz2
   real epsi,ax,ay,az
   real vv1(4),vv2(4),vv3(4),vv4(4),vv5(4),vv6(4),vv7(4),vv8(4)
   real vv13(4),vv24(4),vv57(4),vv68(4)
   real vv1324(4),vv5768(4)
   real vvf(n1,n2,kparasta:kparaend,4)
   real xx1(n1,n2,kparasta:kparaend)
   real yy1(n1,n2,kparasta:kparaend)
   real zz1(n1,n2,kparasta:kparaend)
   !-----------------------------------------------------------------------
   !
   !     start cycling on plane
   !
   epsi=1.e-12
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            !
            ax=abs(xx1(i,j,k))
            ist=int((xx1(i,j,k)+ax)/(2.*ax+epsi)+.5)
            inp=i-1+ist
            !
            ay=abs(yy1(i,j,k))
            jst=int((yy1(i,j,k)+ay)/(2.*ay+epsi)+.5)
            jnp=j-1+jst
            !
            az=abs(zz1(i,j,k))
            kst=int((zz1(i,j,k)+az)/(2.*az+epsi)+.5)
            knp=k-1+kst
            !
            vv1(1)=alalm(inp,jnp,knp)
            vv1(2)=alamm(inp,jnp,knp)
            vv1(3)=alalmrho(inp,jnp,knp)
            vv1(4)=alammrho(inp,jnp,knp)
            !
            vv2(1)=alalm(inp+1,jnp,knp)
            vv2(2)=alamm(inp+1,jnp,knp)
            vv2(3)=alalmrho(inp+1,jnp,knp)
            vv2(4)=alammrho(inp+1,jnp,knp)

            !
            vv3(1)=alalm(inp,jnp,knp+1)
            vv3(2)=alamm(inp,jnp,knp+1)
            vv3(3)=alalmrho(inp,jnp,knp+1)
            vv3(4)=alammrho(inp,jnp,knp+1)

            !
            vv4(1)=alalm(inp+1,jnp,knp+1)
            vv4(2)=alamm(inp+1,jnp,knp+1)
            vv4(3)=alalmrho(inp+1,jnp,knp+1)
            vv4(4)=alammrho(inp+1,jnp,knp+1)
            !
            vv5(1)=alalm(inp,jnp+1,knp)
            vv5(2)=alamm(inp,jnp+1,knp)
            vv5(3)=alalmrho(inp,jnp+1,knp)
            vv5(4)=alammrho(inp,jnp+1,knp)
            !
            vv6(1)=alalm(inp+1,jnp+1,knp)
            vv6(2)=alamm(inp+1,jnp+1,knp)
            vv6(3)=alalmrho(inp+1,jnp+1,knp)
            vv6(4)=alammrho(inp+1,jnp+1,knp)
            !
            vv7(1)=alalm(inp,jnp+1,knp+1)
            vv7(2)=alamm(inp,jnp+1,knp+1)
            vv7(3)=alalmrho(inp,jnp+1,knp+1)
            vv7(4)=alammrho(inp,jnp+1,knp+1)
            !
            vv8(1)=alalm(inp+1,jnp+1,knp+1)
            vv8(2)=alamm(inp+1,jnp+1,knp+1)
            vv8(3)=alalmrho(inp+1,jnp+1,knp+1)
            vv8(4)=alammrho(inp+1,jnp+1,knp+1)
            !
            ! trova distanze in spazio computazionale
            !
            !      if(xx1(i,j,k).gt.0.) then
            !      dx1=xx1(i,j,k)
            !      dx2=1-xx1(i,j,k)
            !      else
            !      dx1=1.+xx1(i,j,k)
            !      dx2=-xx1(i,j,k)
            !      end if
            !
            dx1=    xx1(i,j,k) *ist + (1.+xx1(i,j,k))*(1-ist)
            dx2=(1.-xx1(i,j,k))*ist + (  -xx1(i,j,k))*(1-ist)
            !
            !      if(yy1(i,j,k).gt.0.) then
            !      dy1=yy1(i,j,k)
            !      dy2=1-yy1(i,j,k)
            !      else
            !      dy1=1.+yy1(i,j,k)
            !      dy2=-yy1(i,j,k)
            !      end if
            !
            dy1=    yy1(i,j,k) *jst + (1.+yy1(i,j,k))*(1-jst)
            dy2=(1.-yy1(i,j,k))*jst + (  -yy1(i,j,k))*(1-jst)
            !
            !      if(zz1(i,j,k).gt.0.) then
            !      dz1=zz1(i,j,k)
            !      dz2=1.-zz1(i,j,k)
            !      else
            !      dz1=1.+zz1(i,j,k)
            !      dz2=-zz1(i,j,k)
            !      end if
            !
            dz1=    zz1(i,j,k) *kst + (1.+zz1(i,j,k))*(1-kst)
            dz2=(1.-zz1(i,j,k))*kst + (  -zz1(i,j,k))*(1-kst)
            !
            ! parte interpolazione
            ! interseca il volume sul piano z=zp
            !

            ! interpola valori delle funzioni sui punti
            !
            do ii=1,4
               vv13(ii)=vv1(ii)*dz2+vv3(ii)*dz1
               vv24(ii)=vv2(ii)*dz2+vv4(ii)*dz1
               vv57(ii)=vv5(ii)*dz2+vv7(ii)*dz1
               vv68(ii)=vv6(ii)*dz2+vv8(ii)*dz1

            end do

            !
            ! adesso interseco il piano con linea x=xp
            !
            do ii=1,4
               vv5768(ii)=vv57(ii)*dx2+vv68(ii)*dx1
               vv1324(ii)=vv13(ii)*dx2+vv24(ii)*dx1
            end do
            !
            ! infine interpola sul punto
            !
            do ii=1,4
               vvf(i,j,k,ii)=vv1324(ii)*dy2+vv5768(ii)*dy1
            end do
         !
         enddo
      enddo
   enddo
   !
   return
end



