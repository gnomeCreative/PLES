subroutine turbo_statico(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kpstamg,kpendmg)

    ! compute the eddy viscosity and diffusivity with Smagorinsky
    ! with fixed constant

    use turbo_module
    use filter_module, only: buffer1g,buffer1old_par,buffer1old_par_nscal, &
        buffer2gg,buffer2old_par,buffer2old_par_nscal
    use wallmodel_module, only: att_mod_par,att_wm_sgs,u_t,utangente,wf_distance

    use mysending
    use mysettings, only: cost,costH,costV,pran,prsc,isotropo
    use myarrays_metri3, only: annit,annitV,giac,ref_area,x,y,z
    use myarrays_velo3

    use scala3
    use subgrid
    use period
    use velpar

    use mpi

    implicit none
    !-----------------------------------------------------------------------
    integer,intent(in) :: kpstamg(0:4),kpendmg(0:4)
    integer,intent(in) :: in_dx1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_sn1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_sp1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_st1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_av1(n1,n2,kpstamg(1):kpendmg(1))
    integer,intent(in) :: in_in1(n1,n2,kpstamg(1):kpendmg(1))

    !-----------------------------------------------------------------------
    integer :: i,j,k,isc

    real :: dudx,dudy,dudz
    real :: dvdx,dvdy,dvdz
    real :: dwdx,dwdy,dwdz
    real :: u1,u2,u3,u4,u5,u6
    real :: v1,v2,v3,v4,v5,v6
    real :: w1,w2,w3,w4,w5,w6
    real :: s1,s2,s3,s4,s5,s6

    real :: costante
    real :: lhor,lhorx,lhorz,lver
    real :: costanteV,costanteH

    real, allocatable :: P_akapt(:,:,:),P_akaptV(:,:,:)
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! velocity extrapolation on sides i=0 and i=jx+1

    ! not periodic
    if (ip==1) then
        do k=kparasta,kparaend
            do j=1,n2

                usn(j,k)=2.0*u(0,j,k)-u(1,j,k)
                vsn(j,k)=2.0*v(0,j,k)-v(1,j,k)
                wsn(j,k)=2.0*w(0,j,k)-w(1,j,k)

                udx(j,k)=2.0*u(n1+1,j,k)-u(n1,j,k)
                vdx(j,k)=2.0*v(n1+1,j,k)-v(n1,j,k)
                wdx(j,k)=2.0*w(n1+1,j,k)-w(n1,j,k)

            end do
        end do
    ! periodic
    else if (ip==0) then
        do k=kparasta,kparaend
            do j=1,n2

                usn(j,k)=u(n1,j,k)
                vsn(j,k)=v(n1,j,k)
                wsn(j,k)=w(n1,j,k)

                udx(j,k)=u(1,j,k)
                vdx(j,k)=v(1,j,k)
                wdx(j,k)=w(1,j,k)

            end do
        end do
    end if

    ! extrapolation for velocity on sides j=0 j=jy+1

    ! direction j not periodic
    do k=kparasta,kparaend
        do i=1,n1

            ust(i,k)=2.0*u(i,0,k)-u(i,1,k)
            vst(i,k)=2.0*v(i,0,k)-v(i,1,k)
            wst(i,k)=2.0*w(i,0,k)-w(i,1,k)

            usp(i,k)=2.0*u(i,n2+1,k)-u(i,n2,k)
            vsp(i,k)=2.0*v(i,n2+1,k)-v(i,n2,k)
            wsp(i,k)=2.0*w(i,n2+1,k)-w(i,n2,k)

        end do
    end do


    ! extrapolation for velocity on sides k=0 k=jz+1

    ! not periodic
    if (kp==1) then
        do j=1,n2
            do i=1,n1
                if (myid==0) then

                    uin(i,j)=2.0*u(i,j,0)-u(i,j,1)
                    vin(i,j)=2.0*v(i,j,0)-v(i,j,1)
                    win(i,j)=2.0*w(i,j,0)-w(i,j,1)

                else if (myid==nproc-1) then

                    uav(i,j)=2.0*u(i,j,n3+1)-u(i,j,n3)
                    vav(i,j)=2.0*v(i,j,n3+1)-v(i,j,n3)
                    wav(i,j)=2.0*w(i,j,n3+1)-w(i,j,n3)

                end if
            end do
        end do
    ! periodic
    ! planes for P0 and Pn-1 already known from the previous time step
    else if (kp==0) then

        do j=1,n2
            do i=1,n1

                uin(i,j)=0.0
                vin(i,j)=0.0
                win(i,j)=0.0

                uav(i,j)=0.0
                vav(i,j)=0.0
                wav(i,j)=0.0

                if (myid==0) then
                    uin(i,j)=u_piano(i,j,n3)
                    vin(i,j)=v_piano(i,j,n3)
                    win(i,j)=w_piano(i,j,n3)

                else if (myid==nproc-1) then
                    uav(i,j)=u_piano(i,j,1)
                    vav(i,j)=v_piano(i,j,1)
                    wav(i,j)=w_piano(i,j,1)
                end if

            end do
        end do

    end if

    !-----------------------------------------------------------------------
    ! START COMPUTATION EDDY VISCOSITY AND DIFFUSIVITY

    ! initialization
    do k=kparasta-1,kparaend+1
        do j=0,n2+1
            do i=0,n1+1
                annit(i,j,k)=1.0/re
                annitV(i,j,k)=1.0/re
                do isc=1,nscal
                    akapt(isc,i,j,k)=1.0/re/pran(isc)
                    akaptV(isc,i,j,k)=1.0/re/pran(isc)
                end do
            end do
        end do
    end do

    !--------------------------------------------------------------------
    ! start computation of the controvariant form of the strain rate
    ! tensor Sij, its module is needed to compute the eddy viscosity
    ! with Smagorinsky model

    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1

                u2=in_dx1(i,j,k)*u(i+1,j,k)+(1-in_dx1(i,j,k))*udx(j,k)
                u1=in_sn1(i,j,k)*u(i-1,j,k)+(1-in_sn1(i,j,k))*usn(j,k)
                u4=in_sp1(i,j,k)*u(i,j+1,k)+(1-in_sp1(i,j,k))*usp(i,k)
                u3=in_st1(i,j,k)*u(i,j-1,k)+(1-in_st1(i,j,k))*ust(i,k)
                u6=in_av1(i,j,k)*u(i,j,k+1)+(1-in_av1(i,j,k))*uav(i,j)
                u5=in_in1(i,j,k)*u(i,j,k-1)+(1-in_in1(i,j,k))*uin(i,j)

                v2=in_dx1(i,j,k)*v(i+1,j,k)+(1-in_dx1(i,j,k))*vdx(j,k)
                v1=in_sn1(i,j,k)*v(i-1,j,k)+(1-in_sn1(i,j,k))*vsn(j,k)
                v4=in_sp1(i,j,k)*v(i,j+1,k)+(1-in_sp1(i,j,k))*vsp(i,k)
                v3=in_st1(i,j,k)*v(i,j-1,k)+(1-in_st1(i,j,k))*vst(i,k)
                v6=in_av1(i,j,k)*v(i,j,k+1)+(1-in_av1(i,j,k))*vav(i,j)
                v5=in_in1(i,j,k)*v(i,j,k-1)+(1-in_in1(i,j,k))*vin(i,j)

                w2=in_dx1(i,j,k)*w(i+1,j,k)+(1-in_dx1(i,j,k))*wdx(j,k)
                w1=in_sn1(i,j,k)*w(i-1,j,k)+(1-in_sn1(i,j,k))*wsn(j,k)
                w4=in_sp1(i,j,k)*w(i,j+1,k)+(1-in_sp1(i,j,k))*wsp(i,k)
                w3=in_st1(i,j,k)*w(i,j-1,k)+(1-in_st1(i,j,k))*wst(i,k)
                w6=in_av1(i,j,k)*w(i,j,k+1)+(1-in_av1(i,j,k))*wav(i,j)
                w5=in_in1(i,j,k)*w(i,j,k-1)+(1-in_in1(i,j,k))*win(i,j)

                dudx=(0.5*(u2-u1)*apcsx(i,j,k)+0.5*(u4-u3)*apetx(i,j,k)+0.5*(u6-u5)*apztx(i,j,k))/giac(i,j,k)
                dudy=(0.5*(u2-u1)*apcsy(i,j,k)+0.5*(u4-u3)*apety(i,j,k)+0.5*(u6-u5)*apzty(i,j,k))/giac(i,j,k)
                dudz=(0.5*(u2-u1)*apcsz(i,j,k)+0.5*(u4-u3)*apetz(i,j,k)+0.5*(u6-u5)*apztz(i,j,k))/giac(i,j,k)
                dvdx=(0.5*(v2-v1)*apcsx(i,j,k)+0.5*(v4-v3)*apetx(i,j,k)+0.5*(v6-v5)*apztx(i,j,k))/giac(i,j,k)
                dvdy=(0.5*(v2-v1)*apcsy(i,j,k)+0.5*(v4-v3)*apety(i,j,k)+0.5*(v6-v5)*apzty(i,j,k))/giac(i,j,k)
                dvdz=(0.5*(v2-v1)*apcsz(i,j,k)+0.5*(v4-v3)*apetz(i,j,k)+0.5*(v6-v5)*apztz(i,j,k))/giac(i,j,k)
                dwdx=(0.5*(w2-w1)*apcsx(i,j,k)+0.5*(w4-w3)*apetx(i,j,k)+0.5*(w6-w5)*apztx(i,j,k))/giac(i,j,k)
                dwdy=(0.5*(w2-w1)*apcsy(i,j,k)+0.5*(w4-w3)*apety(i,j,k)+0.5*(w6-w5)*apzty(i,j,k))/giac(i,j,k)
                dwdz=(0.5*(w2-w1)*apcsz(i,j,k)+0.5*(w4-w3)*apetz(i,j,k)+0.5*(w6-w5)*apztz(i,j,k))/giac(i,j,k)

                !-----------------------------------------------------------------------
                ! for wall modeling, dudy and dwdy computation considering the coarse grid

                ! chicco mettere anche in funzione di coeff wall
                if (att_wm_sgs) then
                    if (j==1) then
                        if (att_mod_par(i,1,k)) then
                            dudy=(u_t(i,1,k)/(0.41*wf_distance(i,k,1)))*u(i,j,k)/utangente(i,1,k)
                            dwdy=(u_t(i,1,k)/(0.41*wf_distance(i,k,1)))*w(i,j,k)/utangente(i,1,k)
                        end if
                    else if (j==n2) then
                        if (att_mod_par(i,2,k)) then
                            dudy=(u_t(i,2,k)/(0.41*wf_distance(i,k,2)))*u(i,j,k)/utangente(i,2,k)
                            dwdy=(u_t(i,2,k)/(0.41*wf_distance(i,k,2)))*w(i,j,k)/utangente(i,2,k)
                        end if
                    end if
                end if
                !-----------------------------------------------------------------------

                s1=dudx
                s2=0.5*(dudy+dvdx)
                s3=0.5*(dudz+dwdx)
                s4=dvdy
                s5=0.5*(dvdz+dwdy)
                s6=dwdz

                smod(i,j,k)=sqrt(2.0*s1*s1+2.0*s4*s4+2.0*s6*s6+4.0*s2*s2+4.0*s3*s3+4.0*s5*s5)

                smodV(i,j,k)=sqrt(2.0*s4*s4+4.0*s2*s2+4.0*s5*s5)

                smodH(i,j,k)=sqrt(2.0*s1*s1+2.0*s6*s6+4.0*s3*s3)

                q_crit(i,j,k)=-0.5*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)-dudy*dvdx-dwdx*dudz-dvdz*dwdy

            end do
        end do
    end do
    if (isotropo) then
        ! isotropic model, one eddy viscosity

        costante=cost*cost
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1

                    annit(i,j,k)=annit(i,j,k)+ref_area(i,j,k)*costante*smod(i,j,k)
                    annitV(i,j,k)=annit(i,j,k)

                    ! prandtl turbolent=0.5 ==> Kt=2*Vt
                    do isc=1,nscal
                        akapt(isc,i,j,k)=akapt(isc,i,j,k)+prsc(isc)*ref_area(i,j,k)*costante*smod(i,j,k)
                        akaptV(isc,i,j,k)=akapt(isc,i,j,k)
                    end do

                end do
            end do
        end do
    else
        ! anisotropic model, two eddy viscosity
        costanteH=costH*costH
        costanteV=costV*costV

        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1

                    lhorx=.25*(x(i,j,k)+x(i,j,k-1)+x(i,j-1,k-1)+x(i,j-1,k)) &
                        -.25*(x(i-1,j,k)+x(i-1,j,k-1)+x(i-1,j-1,k-1)+x(i-1,j-1,k))

                    lhorz=.25*(z(i,j,k)+z(i,j-1,k)+z(i-1,j-1,k)+z(i-1,j,k)) &
                        -.25*(z(i,j,k-1)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j,k-1))

                    ! lhor=lhorx**2.+lhorz**2. !length square
                    lhor=lhorx*lhorx+lhorz*lhorz !length square

                    lver=.25*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k)) &
                        -.25*(y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k))

                    lver=lver*lver !**2.


                    annit(i,j,k)=annit(i,j,k)+lhor*costanteH*smodH(i,j,k)
                    annitV(i,j,k)=annitV(i,j,k)+lver*costanteV*smodV(i,j,k)

                    ! prandtl turbolent=0.5==>  Kt=2*Vt
                    do isc=1,nscal

                        akapt(isc,i,j,k)=akapt(isc,i,j,k)+prsc(isc)*lhor*costanteH*smodH(i,j,k)

                        ! prandtl turbolent=0.5==>  Kt=2*Vt
                        akaptV(isc,i,j,k)=akaptV(isc,i,j,k)+prsc(isc)*lver*costanteV*smodV(i,j,k)

                    end do

                end do
            end do
        end do
    end if

    ! add eddy viscosity and diffusivity on periodic sides
    call var_complete_exchange(annit)
    call var_complete_exchange(annitV)

    do isc=1,nscal

        allocate(P_akapt(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
        allocate(P_akaptV(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))

        P_akapt(:,:,:)=akapt(isc,:,:,:)
        P_akaptV(:,:,:)=akaptV(isc,:,:,:)

        call var_complete_exchange(P_akapt)
        call var_complete_exchange(P_akaptV)

        akapt(isc,:,:,:)=P_akapt(:,:,:)
        akaptV(isc,:,:,:)=P_akaptV(:,:,:)

        deallocate(P_akapt,P_akaptV)

    end do

end subroutine turbo_statico

