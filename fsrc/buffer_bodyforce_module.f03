module buffer_bodyforce_module

    use mysending
    use myarrays_velo3
    use contour_module
    use nesting_module
    use orlansky_module
    !
    use scala3
    use tipologia
    !
    use mpi

    implicit none

    !     face 1
    real, allocatable :: noise_f1_xtime(:,:),old_noise_f1_xtime(:)
    real, allocatable :: noise_f1_ytime(:,:),old_noise_f1_ytime(:)
    real, allocatable :: noise_f1_ztime(:,:),old_noise_f1_ztime(:)
    real, allocatable :: noise_f1_xspace(:,:)
    real, allocatable :: noise_f1_yspace(:,:)
    real, allocatable :: noise_f1_zspace(:,:)
    real, allocatable :: bcsi_f1(:,:,:)
    real, allocatable :: beta_f1(:,:,:)
    real, allocatable :: bzet_f1(:,:,:)

    !     face 2
    real, allocatable :: noise_f2_xtime(:,:),old_noise_f2_xtime(:)
    real, allocatable :: noise_f2_ytime(:,:),old_noise_f2_ytime(:)
    real, allocatable :: noise_f2_ztime(:,:),old_noise_f2_ztime(:)
    real, allocatable :: noise_f2_xspace(:,:)
    real, allocatable :: noise_f2_yspace(:,:)
    real, allocatable :: noise_f2_zspace(:,:)
    real, allocatable :: bcsi_f2(:,:,:)
    real, allocatable :: beta_f2(:,:,:)
    real, allocatable :: bzet_f2(:,:,:)

    !     face 3
    real, allocatable :: noise_f5_xtime(:,:),old_noise_f5_xtime(:)
    real, allocatable :: noise_f5_ytime(:,:),old_noise_f5_ytime(:)
    real, allocatable :: noise_f5_ztime(:,:),old_noise_f5_ztime(:)
    real, allocatable :: noise_f5_xspace(:,:)
    real, allocatable :: noise_f5_yspace(:,:)
    real, allocatable :: noise_f5_zspace(:,:)
    real, allocatable :: bcsi_f5(:,:,:)
    real, allocatable :: beta_f5(:,:,:)
    real, allocatable :: bzet_f5(:,:,:)
      
    !     face 4
    real, allocatable :: noise_f6_xtime(:,:),old_noise_f6_xtime(:)
    real, allocatable :: noise_f6_ytime(:,:),old_noise_f6_ytime(:)
    real, allocatable :: noise_f6_ztime(:,:),old_noise_f6_ztime(:)
    real, allocatable :: noise_f6_xspace(:,:)
    real, allocatable :: noise_f6_yspace(:,:)
    real, allocatable :: noise_f6_zspace(:,:)
    real, allocatable :: bcsi_f6(:,:,:)
    real, allocatable :: beta_f6(:,:,:)
    real, allocatable :: bzet_f6(:,:,:)
      
    real,bind(C) :: corr_factor
      
    integer,bind(C) :: ispon,kspon
    integer ntpom
      
    public :: buffer_bodyforce
    private :: sponge_dumping,sponge_dumping_visco
    private :: disturbance_face1,disturbance_face2,disturbance_face5,disturbance_face6
    private :: creo_xyz_f1,creo_xyz_f2,creo_xyz_f5,creo_xyz_f6

contains


    subroutine buffer_bodyforce(ti,ktime,bcsi,beta,bzet)

        ! generates a coloured noise to start turbulence from tke information
        ! tke comes from a RANS model, this sub is used for nesting procedure

        use mysending
        use myarrays_metri3

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer ktime
        integer nti,ntime

        integer i,j,k,n
        integer cont
        !      integer ispon,kspon
        integer istampo

        real cx,cy,cz
        real L_d
        real bcsi(n1,n2,kparasta:kparaend)
        real beta(n1,n2,kparasta:kparaend)
        real bzet(n1,n2,kparasta:kparaend)


        real apply_dist1(n1,kparasta:kparaend)
        real apply_dist2(n1,kparasta:kparaend)
        real apply_dist5(n1,n2)
        real apply_dist6(n1,n2)
        real mass_sign

        real, allocatable :: var_prov(:,:,:),var_prov2(:,:,:)
        real, allocatable :: local_dump(:,:,:)

        real v1,v2,ti
        !
        !-----------------------------------------------------------------------
        ! some parameter
        !      n_ti_pom = 1
        !
        !-----------------------------------------------------------------------
        !for tkepom
        if(ktime .eq.1)then
            ntpom = 1
        end if
        if(myid==0)write(*,*)'BUFFER TURBULENCE on'
        if(ibodybuffer1==1)write(*,*)'side 1 on'
        if(ibodybuffer2==1)write(*,*)'side 2 on'
        if(ibodybuffer5==1)write(*,*)'side 5 on'
        if(ibodybuffer6==1)write(*,*)'side 6 on'



        if(ti .gt. ti_pom_new)ntpom = ntpom + 1
        !-----------------------------------------------------------------------

        !     scaling factor for forcing term
        cx = 0.1
        cy = 0.1
        cz = 0.1

        !     time of colored disturbance
        nti = 100

        !     check sponge region dimension in i and k directions
        if(kspon .gt. jz/nproc .and. myid.eq.0)then
            write(*,*)'NESTING PROBLEM'
            write(*,*)'program not suited to have a sponge region'
            write(*,*)'larger than the processor z allocation'
            stop
        end if
        !
        !-----------------------------------------------------------------------
        !     allocation and initialization
        if(ktime .eq. 1)then
            !        allocate(tkepom1(jy,kparasta:kparaend,n_ti_pom)) !face 1
            !        allocate(tkepom2(jy,kparasta:kparaend,n_ti_pom)) !face 2
            !        allocate(tkepom5(jx,jy               ,n_ti_pom)) !face 3
            !        allocate(tkepom6(jx,jy               ,n_ti_pom)) !face 4
            ! tkepom1 = 0.
            ! tkepom2 = 0.
            ! tkepom5 = 0.
            ! tkepom6 = 0.

            allocate(bcsi_f1(ispon,jy,kparasta:kparaend))
            allocate(beta_f1(ispon,jy,kparasta:kparaend))
            allocate(bzet_f1(ispon,jy,kparasta:kparaend))
            bcsi_f1 = 0.
            beta_f1 = 0.
            bzet_f1 = 0.

            allocate(bcsi_f2(jx-ispon+1:jx,jy,kparasta:kparaend))
            allocate(beta_f2(jx-ispon+1:jx,jy,kparasta:kparaend))
            allocate(bzet_f2(jx-ispon+1:jx,jy,kparasta:kparaend))
            bcsi_f2 = 0.
            beta_f2 = 0.
            bzet_f2 = 0.

            allocate(bcsi_f5(jx,jy,kspon))
            allocate(beta_f5(jx,jy,kspon))
            allocate(bzet_f5(jx,jy,kspon))
            bcsi_f5 = 0.
            beta_f5 = 0.
            bzet_f5 = 0.

            allocate(bcsi_f6(jx,jy,jz-kspon+1:jz))
            allocate(beta_f6(jx,jy,jz-kspon+1:jz))
            allocate(bzet_f6(jx,jy,jz-kspon+1:jz))
            bcsi_f6 = 0.
            beta_f6 = 0.
            bzet_f6 = 0.
        end if
        !
        !-----------------------------------------------------------------------
        !     DSITURBANCE
        !     on face 1
        !      if(infout1 /= 0 .or. ibodybuffer1==1)then
        if( ibodybuffer1==1)then
            call disturbance_face1(ti,ktime,cx,cy,cz,nti)
        end if

        !     on face 2
        if(ibodybuffer2==1)then
            call disturbance_face2(ti,ktime,cx,cy,cz,nti)
        end if

        !     on face 5
        if(ibodybuffer5==1)then
            if(myid==0)then
                call disturbance_face5(ti,ktime,cx,cy,cz,nti)
            end if
        end if

        !     on face 6
        if(ibodybuffer6==1)then
            if(myid==nproc-1)then
                call disturbance_face6(ti,ktime,cx,cy,cz,nti)
            end if
        end if
        !
        !-----------------------------------------------------------------------
        !     SPOONGE DUMPING

        !     face 1
        istampo = 1
        allocate(local_dump(ispon,n2,kparasta:kparaend))
        L_d = (x(ispon,1,kparasta)-x(0,1,kparasta))
        !      if(myid.eq.0)write(*,*)'L_d face 1 ',L_d
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,ispon
                    local_dump(i,j,k) = x(i,j,k) - x(0,j,k)
                    local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
                end do
            end do
        end do
        !quando faccia la f2 dovrebbe bastare invertire la local_dump?
        call sponge_dumping(ispon,n2,n3/nproc,local_dump,bcsi_f1,beta_f1,bzet_f1,istampo)
        deallocate(local_dump)



        !     face 2
        istampo = 1
        allocate(local_dump(n1+1-ispon:n1,n2,kparasta:kparaend))
        L_d = (x(jx,1,kparasta)-x(jx-ispon,1,kparasta))
        !      if(myid.eq.0)write(*,*)'L_d face 2 ',L_d
        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx  !1,ispon
                    local_dump(i,j,k) = x(jx,j,k) - x(i-1,j,k)
                    local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
                end do
            end do
        end do
        !quando faccia la f2 dovrebbe bastare invertire la local_dump?
        call sponge_dumping(ispon,n2,n3/nproc,local_dump,bcsi_f2,beta_f2,bzet_f2,istampo)
        deallocate(local_dump)

        !     face 5
        if(myid.eq.0)then
            istampo = 1
            allocate(local_dump(n1,n2,kspon))
            L_d = (z(1,1,kspon)-z(1,1,0))
            !      if(myid.eq.0)write(*,*)'L_d face 5 ',L_d
            do k=1,kspon !kparasta,kparaend
                do j=1,jy
                    do i=1,jx  !ispon
                        local_dump(i,j,k) = z(i,j,k) - z(i,j,0)
                        local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
                    !    if(i.eq.1 .and. j.eq.10)then
                    !       write(452,*)k,local_dump(i,j,k)
                    !    end if
                    end do
                end do
            end do

            call sponge_dumping(n1,n2,kspon,local_dump,bcsi_f5,beta_f5,bzet_f5,istampo)
            deallocate(local_dump)
        end if



        !     face 6
        if(myid.eq.nproc-1)then
            istampo = 1
            allocate(local_dump(n1,n2,n3+1-kspon:n3))
            L_d = (z(1,1,jz)-z(1,1,jz-kspon))
            !      if(myid.eq.nproc-1)write(*,*)'L_d face 6 ',L_d
            do k=jz+1-kspon,jz !1,kspon !kparasta,kparaend
                do j=1,jy
                    do i=1,jx  !ispon
                        local_dump(i,j,k) = z(i,j,jz) - z(i,j,k-1)
                        local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
                    !    if(i.eq.1 .and. j.eq.10)then
                    !       write(453,*)k,local_dump(i,j,k)
                    !    end if

                    end do
                end do
            end do

            call sponge_dumping(n1,n2,kspon,local_dump,bcsi_f6,beta_f6,bzet_f6,istampo)
            deallocate(local_dump)
        end if


        !----------------------------------------------------------------
        allocate( var_prov(jx+1-ispon:jx,jy,kparasta:kparaend))
        allocate(var_prov2(jx+1-ispon:jx,jy,kparasta:kparaend))
        allocate(local_dump(n1+1-ispon:n1,n2,kparasta:kparaend))

        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx  !1,ispon
                    var_prov(i,j,k) = annit(i,j,k)
                    var_prov2(i,j,k) = annitV(i,j,k)
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx  !1,ispon
                    !         local_dump(i,j,k) = x(jx,j,k) - x(i-1,j,k)
                    local_dump(i,j,k) = (x(jx,j,k) - x(i-1,j,k) ) -  L_d
                    local_dump(i,j,k) = 1. - local_dump(i,j,k) * L_d
                end do
            end do
        end do

        call sponge_dumping_visco(ispon,n2,n3/nproc,local_dump,var_prov,var_prov,var_prov2,istampo)

        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx  !1,ispon
                !          annit(i,j,k) =  var_prov(i,j,k)
                !         annitV(i,j,k) = var_prov2(i,j,k)
                end do
            end do
        end do
        deallocate(var_prov)
        deallocate(var_prov2)

        deallocate(local_dump)


        !-----------------------------------------------------------------------
        !     CORNER CORRECTION BETWEEN SPONGE REGION IN i AND k
        !-----------------------------------------------------------------------
        !     CHECK IF DISTURBANCE MUST BE APPLIED
        apply_dist1 = 0.
        apply_dist2 = 0.
        do k=kparasta,kparaend
            do j=1,jy
                         !side 1
                ! uc>0 inflow; uc<0 outflow
                !inflow  mass_sign = 1
                !outflow mass_sign =-1
                mass_sign = INT(sign(1.,uc(0,j,k)))

                !inflow  apply_dist = 1.
                !outflow apply_dist = 0.
                apply_dist1(j,k) = 0.5*(1+mass_sign)
                !--------------------------------------

                         !side 2
                ! uc<0 inflow; uc>0 outflow
                !inflow  mass_sign = 1
                !outflow mass_sign =-1
                mass_sign = - INT(sign(1.,uc(jx,j,k)))

                !inflow  apply_dist = 1.
                !outflow apply_dist = 0.
                apply_dist2(j,k) = 0.5*(1+mass_sign)
            end do
        end do

        if(myid==0)then
            apply_dist5 = 0.
            do j=1,jy
                do i=1,jx
                             !side 5
                    ! wc>0 inflow; wc<0 outflow
                    !inflow  mass_sign = 1
                    !outflow mass_sign =-1
                    mass_sign = INT(sign(1.,wc(i,j,0)))

                    !inflow  apply_dist = 1.
                    !outflow apply_dist = 0.
                    apply_dist5(i,j) = 0.5*(1+mass_sign)
                end do
            end do
        end if

        if(myid==nproc-1)then
            apply_dist6 = 0.
            do j=1,jy
                do i=1,jx
                             !side 6
                    ! wc<0 inflow; wc>0 outflow
                    !inflow  mass_sign = 1
                    !outflow mass_sign =-1
                    mass_sign = - INT(sign(1.,wc(i,j,jz)))

                    !inflow  apply_dist = 1.
                    !outflow apply_dist = 0.
                    apply_dist6(i,j) = 0.5*(1+mass_sign)
                end do
            end do
        end if
        !
        !-----------------------------------------------------------------------
        !     GIVE THE VALUE TO THE BODYFORCE

        !     sponge in front of face 1
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,ispon

                    bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f1(i,j,k)*apply_dist1(j,k)
                    beta(i,j,k) = beta(i,j,k) + beta_f1(i,j,k)*apply_dist1(j,k)
                    bzet(i,j,k) = bzet(i,j,k) + bzet_f1(i,j,k)*apply_dist1(j,k)

                end do
            end do
        end do

        !     sponge in front of face 2
        do k=kparasta,kparaend
            do j=1,jy
                do i=jx-ispon+1,jx

                    bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f2(i,j,k)*apply_dist2(j,k)
                    beta(i,j,k) = beta(i,j,k) + beta_f2(i,j,k)*apply_dist2(j,k)
                    bzet(i,j,k) = bzet(i,j,k) + bzet_f2(i,j,k)*apply_dist2(j,k)

                end do
            end do
        end do

        !     sponge in front of face 5
        if(myid.eq.0)then
            do k=1,kspon
                do j=1,jy
                    do i=1,jx

                        bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f5(i,j,k)*apply_dist5(i,j)
                        beta(i,j,k) = beta(i,j,k) + beta_f5(i,j,k)*apply_dist5(i,j)
                        bzet(i,j,k) = bzet(i,j,k) + bzet_f5(i,j,k)*apply_dist5(i,j)

                    end do
                end do
            end do
        end if

        !     sponge in front of face 6
        if(myid.eq.nproc-1)then
            do k=jz-kspon+1,jz
                do j=1,jy
                    do i=1,jx

                        bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f6(i,j,k)*apply_dist6(i,j)
                        beta(i,j,k) = beta(i,j,k) + beta_f6(i,j,k)*apply_dist6(i,j)
                        bzet(i,j,k) = bzet(i,j,k) + bzet_f6(i,j,k)*apply_dist6(i,j)

                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------
        return
    end

    subroutine sponge_dumping(nx,ny,nz,local_dump,var_csi,var_eta,var_zet,istampo)

        !      I use the form:
        !            (u^(n+1)-u^(n))/dt = - sigma(u^(n)-u_desired)
        !      this means that sigma at i=end_domain must be of order 1/dt
        !      I use an exponential form for sigma:
        !            sigma = e^beta*g(x)
        !      with g(x)= L - (Xend - X) with L the length of the sponge region
        !      so at the domain end g(x) = L
        !      beta is chosen so that sigma has order 1/dt

        !       L_sponge = x(ispon)-x(0)

        !       angbeta = (1/L_sponge)*log(1/dt)
        integer i,j,k
        integer nx,ny,nz
        integer istampo

        real local_dump(nx,ny,nz)

        real var_csi(nx,ny,nz)
        real var_eta(nx,ny,nz)
        real var_zet(nx,ny,nz)

        real var_csi_damp(nx,ny,nz)
        real var_eta_damp(nx,ny,nz)
        real var_zet_damp(nx,ny,nz)

        real esponente,sigma,Acoef

        Acoef = 1.

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    esponente = -.5*( 3.5 * local_dump(i,j,k) )**2.
                    sigma = Acoef * 2.7183**esponente
                    !   if(j.eq.1 .and. k.eq.1 .and. istampo.eq.1)then
                    !     write(600,*)i,esponente,sigma
                    !   end if
                    var_csi_damp(i,j,k) = var_csi(i,j,k)-sigma*(var_csi(i,j,k)-0.)
                    var_eta_damp(i,j,k) = var_eta(i,j,k)-sigma*(var_eta(i,j,k)-0.)
                    var_zet_damp(i,j,k) = var_zet(i,j,k)-sigma*(var_zet(i,j,k)-0.)

                end do
            end do
        end do

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    var_csi(i,j,k) = var_csi_damp(i,j,k)
                    var_eta(i,j,k) = var_eta_damp(i,j,k)
                    var_zet(i,j,k) = var_zet_damp(i,j,k)
                end do
            end do
        end do

        return
    end

    subroutine sponge_dumping_visco(nx,ny,nz,local_dump,var_csi,var_eta,var_zet,istampo)

        !      I use the form:
        !            (u^(n+1)-u^(n))/dt = - sigma(u^(n)-u_desired)
        !      this means that sigma at i=end_domain must be of order 1/dt
        !      I use an exponential form for sigma:
        !            sigma = e^beta*g(x)
        !      with g(x)= L - (Xend - X) with L the length of the sponge region
        !      so at the domain end g(x) = L
        !      beta is chosen so that sigma has order 1/dt

        !       L_sponge = x(ispon)-x(0)

        !       angbeta = (1/L_sponge)*log(1/dt)
        integer i,j,k
        integer nx,ny,nz
        integer istampo

        real local_dump(nx,ny,nz)

        real var_csi(nx,ny,nz)
        real var_eta(nx,ny,nz)
        real var_zet(nx,ny,nz)

        real var_csi_damp(nx,ny,nz)
        real var_eta_damp(nx,ny,nz)
        real var_zet_damp(nx,ny,nz)

        real esponente,sigma,Acoef

        Acoef = 1.

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    esponente = -.5*( 3.5 * local_dump(i,j,k) )**2.
                    sigma = Acoef * 2.7183**esponente
                    !   if(j.eq.1 .and. k.eq.1 .and. istampo.eq.1)then
                    !     write(600,*)i,esponente,sigma
                    !   end if



                    var_csi_damp(i,j,k) = var_csi(i,j,k)-sigma*(var_csi(i,j,k)-4.*var_csi(i,j,k))
                    var_eta_damp(i,j,k) = var_eta(i,j,k)-sigma*(var_eta(i,j,k)-4.*var_eta(i,j,k))
                    var_zet_damp(i,j,k) = var_zet(i,j,k)-sigma*(var_zet(i,j,k)-4.*var_zet(i,j,k))



                end do
            end do
        end do

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    var_csi(i,j,k) = var_csi_damp(i,j,k)
                    var_eta(i,j,k) = var_eta_damp(i,j,k)
                    var_zet(i,j,k) = var_zet_damp(i,j,k)
                end do
            end do
        end do

        return
    end

    subroutine disturbance_face1(ti,ktime,cx,cy,cz,nti)

        ! generate disturbance for face 1 for nesting procedure
        use contour_module, only: up1,vp1,wp1

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer ktime,ntime
        integer i,j,k,n
        integer nti
        !      integer ispon
        integer istampo
        integer cont

        integer nreal,nstep

        real corr_factor_sp
        real cx,cy,cz

        real uprime,vprime,wprime

        real coef_x(ispon,n2,kparasta:kparaend)
        real coef_y(ispon,n2,kparasta:kparaend)
        real coef_z(ispon,n2,kparasta:kparaend)

        real ti,dt1,dt2,dt_tot
        real tke,tke_new
        integer cont_mean
        real meanx,meany,meanz
        real dev_x,dev_y,dev_z
        real var_x,var_y,var_z
        real velrms
        real segno_u,segno_v,segno_w
        real v_scale

        real a1,a2,a3,a_den
        real inv_dt

        real ke,ke_ref
        !-----------------------------------------------------------------------
        inv_dt = 1./dt
        !-----------------------------------------------------------------------
        ! allocation and initialization
        if(ktime .eq. 1)then

            corr_factor_sp = 0.

            !       face 1
            !       in x in space and time
            nreal = jy*jz/nproc
            nstep = nti
            allocate(noise_f1_xtime(nreal,-1:nstep*2))
            allocate(old_noise_f1_xtime(nreal))
            noise_f1_xtime      = 0.
            old_noise_f1_xtime  = 0.

            nreal = jy*jz/nproc
            nstep = ispon
            allocate(noise_f1_xspace(nreal,-1:nstep*2))
            noise_f1_xspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f1_xspace,istampo,myid)

            !       in y in space and time
            nreal = ispon*jz/nproc
            nstep = nti
            allocate(noise_f1_ytime(nreal,-1:nstep*2))
            allocate(old_noise_f1_ytime(nreal))
            noise_f1_ytime     = 0.
            old_noise_f1_ytime = 0.

            nreal = ispon*jz/nproc
            nstep = jy
            allocate(noise_f1_yspace(nreal,-1:nstep*2))
            noise_f1_yspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f1_yspace,istampo,myid)

            !       in z in space and time
            nreal = jy*ispon
            nstep = nti
            allocate(noise_f1_ztime(nreal,-1:nstep*2))
            allocate(old_noise_f1_ztime(nreal))
            noise_f1_ztime     = 0.
            old_noise_f1_ztime = 0.

            nreal = jy*ispon
            nstep = jz/nproc
            allocate(noise_f1_zspace(nreal,-1:nstep*2))
            noise_f1_zspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f1_zspace,istampo,myid)
        end if
        !-----------------------------------------------------------------------
        !     to have a colored noise connected to the previous number generation
        !
        if(ktime .gt.1 .and. ktime .eq. nti*(ktime/nti))then
            ntime = nti
            if(myid.eq.0)write(*,*)'TIME',ktime,ntime
            cont = 0
            do k=kparasta,kparaend
                do j=1,jy
                    cont = cont + 1
                    old_noise_f1_xtime(cont) = noise_f1_xtime(cont,ntime)
                end do
            end do

            cont = 0
            do k=kparasta,kparaend
                do i=1,ispon
                    cont = cont + 1
                    old_noise_f1_ytime(cont) = noise_f1_ytime(cont,ntime)
                end do
            end do

            cont = 0
            do j=1,jy
                do i=1,ispon
                    cont = cont + 1
                    old_noise_f1_ztime(cont) = noise_f1_ztime(cont,ntime)
                end do
            end do
        end if


        !-----------------------------------------------------------------------
        !     DISTURBANCE GENERATION IN TIME
        if(ktime .eq. 1 .or. ktime.eq.nti*(ktime/nti))then

            !      corr_factor = 100.

            !     in x in time
            nreal = jy*jz/nproc
            nstep = nti
            noise_f1_xtime=0.
            call genero_random(corr_factor,nreal,nstep,noise_f1_xtime,istampo,myid)

            !     in y in time
            nreal = ispon*jz/nproc
            nstep = nti
            noise_f1_ytime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f1_ytime,istampo,myid)

            !     in z in time
            nreal = jy*ispon
            nstep = nti
            noise_f1_ztime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f1_ztime,istampo,myid)

        end if  !if ktime

        !-----------------------------------------------------------------------
        !     shift for disturbance
        if(ktime .gt. 1 .and. ktime .eq. nti*(ktime/nti))then
            do ntime = 1,nti*2
                cont = 0
                do k=kparasta,kparaend
                    do j=1,jy
                        cont = cont + 1
                        noise_f1_xtime(cont,ntime)=noise_f1_xtime(cont,ntime)+(old_noise_f1_xtime(cont)-noise_f1_xtime(cont,1))
                    end do
                end do

                cont = 0
                do k=kparasta,kparaend
                    do i=1,ispon
                        cont = cont + 1
                        noise_f1_ytime(cont,ntime)=noise_f1_ytime(cont,ntime)+(old_noise_f1_ytime(cont)-noise_f1_ytime(cont,1))
                    end do
                end do

                cont = 0
                do j=1,jy
                    do i=1,ispon
                        cont = cont + 1
                        noise_f1_ztime(cont,ntime)=noise_f1_ztime(cont,ntime)+(old_noise_f1_ztime(cont)-noise_f1_ztime(cont,1))
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------
        !     noise in space and time in the sponge region

        ntime = ktime -(ktime/nti)*nti
        call creo_xyz_f1(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)


        meanx = 0.
        meany = 0.
        meanz = 0.
        cont_mean = 0
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,ispon
                    cont_mean = cont_mean + 1
                    meanx = meanx + coef_x(i,j,k)
                    meany = meany + coef_y(i,j,k)
                    meanz = meanz + coef_z(i,j,k)
                end do
            end do
        end do
        meanx = meanx/real(cont_mean)
        meany = meany/real(cont_mean)
        meanz = meanz/real(cont_mean)
        !      if(myid ==0)write(*,*)'MEAN:',meanx,meany,meanz
        dev_x = 0.
        dev_y = 0.
        dev_z = 0.
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,ispon
                    coef_x(i,j,k) = coef_x(i,j,k)-meanx
                    coef_y(i,j,k) = coef_y(i,j,k)-meany
                    coef_z(i,j,k) = coef_z(i,j,k)-meanz

                    dev_x = dev_x + coef_x(i,j,k)*coef_x(i,j,k)
                    dev_y = dev_y + coef_y(i,j,k)*coef_y(i,j,k)
                    dev_z = dev_z + coef_z(i,j,k)*coef_z(i,j,k)
                end do
            end do
        end do

        dev_x = sqrt(dev_x/real(cont_mean))
        dev_y = sqrt(dev_y/real(cont_mean))
        dev_z = sqrt(dev_z/real(cont_mean))
        var_x = 0.05
        var_y = 0.05
        var_z = 0.05
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,ispon
                    coef_x(i,j,k) = coef_x(i,j,k)*var_x/dev_x
                    coef_y(i,j,k) = coef_y(i,j,k)*var_y/dev_y
                    coef_z(i,j,k) = coef_z(i,j,k)*var_z/dev_z
                end do
            end do
        end do
        !-----------------------------------------------------------------------
        !     for interpolation in time of tke
        !      dt1 = ti - ti_pom_old
        !      dt2 = ti_pom_new - ti
        !      dt_tot = 1./(dt1+dt2)

        bcsi_f1 = 0.
        beta_f1 = 0.
        bzet_f1 = 0.


        do k=kparasta,kparaend
            do j=1,jy
                do i=1,ispon

                    if(index_out1(j,k) .gt. 0.01)then



                        !fluctuation
                        uprime = u(i,j,k)-up1(j,k)
                        vprime = v(i,j,k)-vp1(j,k)
                        wprime = w(i,j,k)-wp1(j,k)

                        ! coef to scale tke for each direction a1+a2+a3 = 1
                        a_den=1./( up1(j,k)*up1(j,k)+vp1(j,k)*vp1(j,k)+wp1(j,k)*wp1(j,k) )

                        a1 = 1. !up1(j,k)*up1(j,k) * a_den
                        a2 = 1. !vp1(j,k)*vp1(j,k) * a_den
                        a3 = 1. !wp1(j,k)*wp1(j,k) * a_den

                                 !kinetic energy
                        ke_ref =   up1(j,k)*up1(j,k)+vp1(j,k)*vp1(j,k)+wp1(j,k)*wp1(j,k)

                        ke =   u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)

                        if(ke .gt. 1.3*ke_ref)cycle

                        ! left hand side ( TKE(n+1)-TKE(n) ) / Delta T
                        tke = tke1(j,k)

                        tke_new = .5*(  uprime*uprime + vprime*vprime + wprime*wprime )




                        tke = abs((tke_new - tke)*inv_dt)

                        ! fluctuation velocity scale
                        velrms = sqrt(2.*tke)

                        !         v_scale = up1(j,k)*up1(j,k) + vp1(j,k)*vp1(j,k) + wp1(j,k)*wp1(j,k)

                        !         v_scale = 0.01*sqrt(v_scale)

                        !         velrms = v_scale

                        segno_u = sign(1.,up1(j,k))
                        segno_v = sign(1.,vp1(j,k))
                        segno_w = sign(1.,wp1(j,k))

                        segno_u = sign(1.,uprime)
                        segno_v = sign(1.,vprime)
                        segno_w = sign(1.,wprime)

                        !     write(400+myid,*)ktime,i,j,k,up1(j,k),vp1(j,k),wp1(j,k),velrms,
                        !     >          uprime,vprime,wprime,tke/uprime,coef_x(i,j,k)
                        !         if(coef_x(i,j,k).lt.0.)write(420+myid,*)coef_x(i,j,k)
                        !         if(coef_x(i,j,k).gt.0.)write(430+myid,*)coef_x(i,j,k)

                        if(abs(uprime) .lt. 0.1*velrms )then
                            uprime = 0.1*segno_u*velrms
                        end if
                        if(abs(vprime) .lt. 0.1*velrms )then
                            vprime = 0.1*segno_v*velrms
                        end if
                        if(abs(wprime) .lt. 0.1*velrms )then
                            wprime = 0.1*segno_w*velrms
                        end if


                        ! bodyforce
                        bcsi_f1(i,j,k)=cx*(a1*tke/abs(uprime))*coef_x(i,j,k)

                        beta_f1(i,j,k)=cy*(a2*tke/abs(vprime))*coef_y(i,j,k)

                        bzet_f1(i,j,k)=cz*(a3*tke/abs(wprime))*coef_z(i,j,k)

                        bcsi_f1(i,j,k)=cx*(a1*tke/abs(velrms))*coef_x(i,j,k)

                        beta_f1(i,j,k)=cy*(a2*tke/abs(velrms))*coef_y(i,j,k)

                        bzet_f1(i,j,k)=cz*(a3*tke/abs(velrms))*coef_z(i,j,k)

                        write(500+myid,*)bcsi_f1(i,j,k)*uprime,(tke_new-tke)*inv_dt

                    end if

                end do
            end do
        end do
9200    format(21e18.10) ! there was a label '9200' here

        return
    end

    subroutine creo_xyz_f1(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)

        implicit none
        !-----------------------------------------------------------------------
        integer i,j,k,cont,ntime
        integer n1,n2,n3,kparasta,kparaend
        integer myid

        real coef_x(ispon,n2,kparasta:kparaend)
        real coef_y(ispon,n2,kparasta:kparaend)
        real coef_z(ispon,n2,kparasta:kparaend)
        !-----------------------------------------------------------------------
        !     in x
        do i=1,ispon
            cont = 0
            do k=kparasta,kparaend
                do j=1,n2
                    cont = cont + 1
                    coef_x(i,j,k) = noise_f1_xtime(cont,ntime)+noise_f1_xspace(cont,i)
                end do
            end do
        end do

        !     in y
        do j=1,n2
            cont = 0
            do i=1,ispon
                do k=kparasta,kparaend
                    cont = cont + 1
                    coef_y(i,j,k) = noise_f1_ytime(cont,ntime)+noise_f1_yspace(cont,j)
                end do
            end do
        end do

        !     in z
        do k=kparasta,kparaend
            cont = 0
            do i=1,ispon
                do j=1,n2
                    cont = cont + 1
                    coef_z(i,j,k) = noise_f1_ztime(cont,ntime)+noise_f1_zspace(cont,k-kparasta+1)
                end do
            end do
        end do

        return
    end

    subroutine disturbance_face2(ti,ktime,cx,cy,cz,nti)

        !     generate disturbance for face 2 for nesting procedure
        use contour_module, only: up2,vp2,wp2

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer ktime,ntime
        integer i,j,k,n
        integer nti
        !      integer ispon
        integer istampo
        integer cont


        integer nreal,nstep

        real corr_factor_sp
        real cx,cy,cz

        real uprime,vprime,wprime

        real coef_x(n1+1-ispon:n1,n2,kparasta:kparaend)
        real coef_y(n1+1-ispon:n1,n2,kparasta:kparaend)
        real coef_z(n1+1-ispon:n1,n2,kparasta:kparaend)

        real ti,dt1,dt2,dt_tot
        real tke,tke_new

        integer cont_mean
        real meanx,meany,meanz
        real dev_x,dev_y,dev_z
        real var_x,var_y,var_z
        real velrms
        real segno_u,segno_v,segno_w
        real v_scale

        real a1,a2,a3,a_den
        real inv_dt

        real ke,ke_ref
        !-----------------------------------------------------------------------
        inv_dt = 1./dt
        !-----------------------------------------------------------------------
        ! allocation and initialization
        if(ktime .eq. 1)then

            corr_factor_sp = 0.

            !       face 1
            !       in x in space and time
            nreal = jy*jz/nproc
            nstep = nti
            allocate(noise_f2_xtime(nreal,-1:nstep*2))
            allocate(old_noise_f2_xtime(nreal))
            noise_f2_xtime      = 0.
            old_noise_f2_xtime  = 0.

            nreal = jy*jz/nproc
            nstep = ispon
            allocate(noise_f2_xspace(nreal,-1:nstep*2))
            noise_f2_xspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f2_xspace,istampo,myid)

            !       in y in space and time
            nreal = ispon*jz/nproc
            nstep = nti
            allocate(noise_f2_ytime(nreal,-1:nstep*2))
            allocate(old_noise_f2_ytime(nreal))
            noise_f2_ytime     = 0.
            old_noise_f2_ytime = 0.

            nreal = ispon*jz/nproc
            nstep = jy
            allocate(noise_f2_yspace(nreal,-1:nstep*2))
            noise_f2_yspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f2_yspace,istampo,myid)

            !       in z in space and time
            nreal = jy*ispon
            nstep = nti
            allocate(noise_f2_ztime(nreal,-1:nstep*2))
            allocate(old_noise_f2_ztime(nreal))
            noise_f2_ztime     = 0.
            old_noise_f2_ztime = 0.

            nreal = jy*ispon
            nstep = jz/nproc
            allocate(noise_f2_zspace(nreal,-1:nstep*2))
            noise_f2_zspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f2_zspace,istampo,myid)
        end if
        !-----------------------------------------------------------------------
        !     to have a colored noise connected to the previous number generation
        !
        if(ktime .gt.1 .and. ktime .eq. nti*(ktime/nti))then
            ntime = nti
            if(myid.eq.0)write(*,*)'TIME',ktime,ntime
            cont = 0
            do k=kparasta,kparaend
                do j=1,jy
                    cont = cont + 1
                    old_noise_f2_xtime(cont) = noise_f2_xtime(cont,ntime)
                end do
            end do

            cont = 0
            do k=kparasta,kparaend
                do i=1,ispon
                    cont = cont + 1
                    old_noise_f2_ytime(cont) = noise_f2_ytime(cont,ntime)
                end do
            end do

            cont = 0
            do j=1,jy
                do i=1,ispon
                    cont = cont + 1
                    old_noise_f2_ztime(cont) = noise_f2_ztime(cont,ntime)
                end do
            end do
        end if


        !-----------------------------------------------------------------------
        !     DISTURBANCE GENERATION IN TIME
        if(ktime .eq. 1 .or. ktime.eq.nti*(ktime/nti))then

            !      corr_factor = 100.

            !     in x in time
            nreal = jy*jz/nproc
            nstep = nti
            noise_f2_xtime=0.
            call genero_random(corr_factor,nreal,nstep,noise_f2_xtime,istampo,myid)

            !     in y in time
            nreal = ispon*jz/nproc
            nstep = nti
            noise_f2_ytime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f2_ytime,istampo,myid)

            !     in z in time
            nreal = jy*ispon
            nstep = nti
            noise_f2_ztime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f2_ztime,istampo,myid)

        end if  !if ktime

        !-----------------------------------------------------------------------
        !     shift for disturbance
        if(ktime .gt. 1 .and. ktime .eq. nti*(ktime/nti))then
            do ntime = 1,nti*2
                cont = 0
                do k=kparasta,kparaend
                    do j=1,jy
                        cont = cont + 1
                        noise_f2_xtime(cont,ntime)=noise_f2_xtime(cont,ntime)+(old_noise_f2_xtime(cont)-noise_f2_xtime(cont,1))
                    end do
                end do

                cont = 0
                do k=kparasta,kparaend
                    do i=1,ispon
                        cont = cont + 1
                        noise_f2_ytime(cont,ntime)=noise_f2_ytime(cont,ntime)+(old_noise_f2_ytime(cont)-noise_f2_ytime(cont,1))
                    end do
                end do

                cont = 0
                do j=1,jy
                    do i=1,ispon
                        cont = cont + 1
                        noise_f2_ztime(cont,ntime)=noise_f2_ztime(cont,ntime)+(old_noise_f2_ztime(cont)-noise_f2_ztime(cont,1))
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------
        !     noise in space and time in the sponge region

        ntime = ktime -(ktime/nti)*nti
        call creo_xyz_f2(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)


        meanx = 0.
        meany = 0.
        meanz = 0.
        cont_mean = 0
        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx
                    cont_mean = cont_mean + 1
                    meanx = meanx + coef_x(i,j,k)
                    meany = meany + coef_y(i,j,k)
                    meanz = meanz + coef_z(i,j,k)
                end do
            end do
        end do
        meanx = meanx/real(cont_mean)
        meany = meany/real(cont_mean)
        meanz = meanz/real(cont_mean)
        !      if(myid ==0)write(*,*)'MEAN:',meanx,meany,meanz
        dev_x = 0.
        dev_y = 0.
        dev_z = 0.
        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx
                    coef_x(i,j,k) = coef_x(i,j,k)-meanx
                    coef_y(i,j,k) = coef_y(i,j,k)-meany
                    coef_z(i,j,k) = coef_z(i,j,k)-meanz

                    dev_x = dev_x + coef_x(i,j,k)*coef_x(i,j,k)
                    dev_y = dev_y + coef_y(i,j,k)*coef_y(i,j,k)
                    dev_z = dev_z + coef_z(i,j,k)*coef_z(i,j,k)
                end do
            end do
        end do

        dev_x = sqrt(dev_x/real(cont_mean))
        dev_y = sqrt(dev_y/real(cont_mean))
        dev_z = sqrt(dev_z/real(cont_mean))
        var_x = 0.05
        var_y = 0.05
        var_z = 0.05
        !      if(myid ==0)write(*,*)'DEV:',dev_x,dev_y,dev_z
        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx
                    coef_x(i,j,k) = coef_x(i,j,k)*var_x/dev_x
                    coef_y(i,j,k) = coef_y(i,j,k)*var_y/dev_y
                    coef_z(i,j,k) = coef_z(i,j,k)*var_z/dev_z
                end do
            end do
        end do
        !-----------------------------------------------------------------------
        !     for interpolation in time of tke
        !      dt1 = ti - ti_pom_old
        !      dt2 = ti_pom_new - ti
        !      dt_tot = 1./(dt1+dt2)

        do k=kparasta,kparaend
            do j=1,jy
                do i=jx+1-ispon,jx  !1,ispon
                    if(index_out2(j,k) .gt. 0.01)then

                        !fluctuation
                        uprime = u(i,j,k)-up2(j,k)
                        vprime = v(i,j,k)-vp2(j,k)
                        wprime = w(i,j,k)-wp2(j,k)

                        ! coef to scale tke for each direction a1+a2+a3 = 1
                        a_den=1./( up2(j,k)*up2(j,k)+vp2(j,k)*vp2(j,k)+wp2(j,k)*wp2(j,k) )

                        a1 = up2(j,k)*up2(j,k) * a_den
                        a2 = vp2(j,k)*vp2(j,k) * a_den
                        a3 = wp2(j,k)*wp2(j,k) * a_den

                                 !kinetic energy
                        ke_ref =   up2(j,k)*up2(j,k)+vp2(j,k)*vp2(j,k)+wp2(j,k)*wp2(j,k)

                        ke =   u(i,j,k)*u(i,j,k) +v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)

                        if(ke .gt. 1.3*ke_ref)cycle


                        ! left hand side ( TKE(n+1)-TKE(n) ) / Delta T
                        tke = tke2(j,k)

                        tke_new = .5*(  uprime*uprime + vprime*vprime + wprime*wprime )

                        tke = abs((tke_new - tke)*inv_dt)

                        ! fluctuation velocity scale
                        velrms = sqrt(2.*tke_new)

                        !         v_scale = up2(j,k)*up2(j,k) + vp2(j,k)*vp2(j,k) + wp2(j,k)*wp2(j,k)

                        !         v_scale = 0.01*sqrt(v_scale)

                        !         velrms = v_scale

                        segno_u = sign(1.,up2(j,k))
                        segno_v = sign(1.,vp2(j,k))
                        segno_w = sign(1.,wp2(j,k))

                        segno_u = sign(1.,uprime)
                        segno_v = sign(1.,vprime)
                        segno_w = sign(1.,wprime)


                        if(abs(uprime) .lt. 0.1*velrms )then
                            uprime = 0.1*segno_u*velrms
                        end if
                        if(abs(vprime) .lt. 0.1*velrms )then
                            vprime = 0.1*segno_v*velrms
                        end if
                        if(abs(wprime) .lt. 0.1*velrms )then
                            wprime = 0.1*segno_w*velrms
                        end if


                        ! bodyforce
                        bcsi_f2(i,j,k)=cx*(a1*tke/abs(uprime))*coef_x(i,j,k)

                        beta_f2(i,j,k)=cy*(a2*tke/abs(vprime))*coef_y(i,j,k)

                        bzet_f2(i,j,k)=cz*(a3*tke/abs(wprime))*coef_z(i,j,k)


                    end if
                end do
            end do
        end do
9200    format(21e18.10)

        return
    end

    subroutine creo_xyz_f2(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)

        implicit none
        !-----------------------------------------------------------------------
        integer i,j,k,cont,icont,ntime
        integer n1,n2,n3,kparasta,kparaend
        integer myid

        real coef_x(n1+1-ispon:n1,n2,kparasta:kparaend)
        real coef_y(n1+1-ispon:n1,n2,kparasta:kparaend)
        real coef_z(n1+1-ispon:n1,n2,kparasta:kparaend)
        !-----------------------------------------------------------------------
        !     in x
        icont = 0
        do i=n1+1-ispon,n1 !1,ispon
            cont = 0
            icont = icont + 1
            do k=kparasta,kparaend
                do j=1,n2
                    cont = cont + 1
                    coef_x(i,j,k) = noise_f2_xtime(cont,ntime)+noise_f2_xspace(cont,icont)
                end do
            end do
        end do

        !     in y
        do j=1,n2
            cont = 0
            do i=n1+1-ispon,n1 !1,ispon
                do k=kparasta,kparaend
                    cont = cont + 1
                    coef_y(i,j,k) = noise_f2_ytime(cont,ntime)+noise_f2_yspace(cont,j)
                end do
            end do
        end do

        !     in z
        do k=kparasta,kparaend
            cont = 0
            do i=n1+1-ispon,n1 !1,ispon
                do j=1,n2
                    cont = cont + 1
                    coef_z(i,j,k) = noise_f2_ztime(cont,ntime)+noise_f2_zspace(cont,k-kparasta+1)
                end do
            end do
        end do

        return
    end

    subroutine disturbance_face5(ti,ktime,cx,cy,cz,nti)
        !     generate disturbance for face 5 for nesting procedure
        use contour_module, only: up5,vp5,wp5
        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer ktime,ntime
        integer i,j,k,n
        integer nti
        !      integer kspon
        integer istampo
        integer cont

        integer nreal,nstep

        real corr_factor_sp
        real cx,cy,cz

        real uprime,vprime,wprime

        real coef_x(n1,n2,kspon)
        real coef_y(n1,n2,kspon)
        real coef_z(n1,n2,kspon)

        real ti,dt1,dt2,dt_tot
        real tke,tke_new
        integer cont_mean
        real meanx,meany,meanz
        real dev_x,dev_y,dev_z
        real var_x,var_y,var_z
        real velrms
        real segno_u,segno_v,segno_w
        real v_scale

        real a1,a2,a3,a_den
        real inv_dt

        real ke,ke_ref

        integer icount
        !-----------------------------------------------------------------------
        inv_dt = 1./dt
        !-----------------------------------------------------------------------
        ! allocation and initialization
        if(ktime .eq. 1)then

            corr_factor_sp = 0.

            !       face 5
            !       in x in space and time
            nreal = jy*kspon   !jy*jz/nproc
            nstep = nti
            allocate(noise_f5_xtime(nreal,-1:nstep*2))
            allocate(old_noise_f5_xtime(nreal))
            noise_f5_xtime      = 0.
            old_noise_f5_xtime  = 0.

            nreal = jy*kspon   !jy*jz/nproc
            nstep = jx         !ispon
            allocate(noise_f5_xspace(nreal,-1:nstep*2))
            noise_f5_xspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f5_xspace,istampo,myid)

            !       in y in space and time
            nreal = jx*kspon   !ispon*jz/nproc
            nstep = nti
            allocate(noise_f5_ytime(nreal,-1:nstep*2))
            allocate(old_noise_f5_ytime(nreal))
            noise_f5_ytime     = 0.
            old_noise_f5_ytime = 0.

            nreal = jx*kspon !ispon*jz/nproc
            nstep = jy
            allocate(noise_f5_yspace(nreal,-1:nstep*2))
            noise_f5_yspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f5_yspace,istampo,myid)

            !       in z in space and time
            nreal = jy*jx  !jy*ispon
            nstep = nti
            allocate(noise_f5_ztime(nreal,-1:nstep*2))
            allocate(old_noise_f5_ztime(nreal))
            noise_f5_ztime     = 0.
            old_noise_f5_ztime = 0.

            nreal = jy*jx !jy*ispon
            nstep = kspon !jz/nproc
            allocate(noise_f5_zspace(nreal,-1:nstep*2))
            noise_f5_zspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f5_zspace,istampo,myid)
        end if
        !-----------------------------------------------------------------------
        !     to have a colored noise connected to the previous number generation
        !
        if(ktime .gt.1 .and. ktime .eq. nti*(ktime/nti))then
            ntime = nti
            if(myid.eq.0)write(*,*)'TIME',ktime,ntime
            cont = 0
            do k=1,kspon !kparasta,kparaend
                do j=1,jy
                    cont = cont + 1
                    old_noise_f5_xtime(cont) = noise_f5_xtime(cont,ntime)
                end do
            end do

            cont = 0
            do k=1,kspon !kparasta,kparaend
                do i=1,jx  !ispon
                    cont = cont + 1
                    old_noise_f5_ytime(cont) = noise_f5_ytime(cont,ntime)
                end do
            end do

            cont = 0
            do j=1,jy
                do i=1,jx  !ispon
                    cont = cont + 1
                    old_noise_f5_ztime(cont) = noise_f5_ztime(cont,ntime)
                end do
            end do
        end if


        !-----------------------------------------------------------------------
        !     DISTURBANCE GENERATION IN TIME
        if(ktime .eq. 1 .or. ktime.eq.nti*(ktime/nti))then

            !      corr_factor = 100.

            !     in x in time
            nreal = jy*kspon !jy*jz/nproc
            nstep = nti
            noise_f5_xtime=0.
            call genero_random(corr_factor,nreal,nstep,noise_f5_xtime,istampo,myid)

            !     in y in time
            nreal = jx*kspon !ispon*jz/nproc
            nstep = nti
            noise_f5_ytime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f5_ytime,istampo,myid)

            !     in z in time
            nreal = jy*jx  !jy*ispon
            nstep = nti
            noise_f5_ztime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f5_ztime,istampo,myid)

        end if  !if ktime

        !-----------------------------------------------------------------------
        !     shift for disturbance
        if(ktime .gt. 1 .and. ktime .eq. nti*(ktime/nti))then
            do ntime = 1,nti*2
                cont = 0
                do k=1,kspon !kparasta,kparaend
                    do j=1,jy
                        cont = cont + 1
                        noise_f5_xtime(cont,ntime)=noise_f5_xtime(cont,ntime)+(old_noise_f5_xtime(cont)-noise_f5_xtime(cont,1))
                    end do
                end do

                cont = 0
                do k=1,kspon !kparasta,kparaend
                    do i=1,jx    !ispon
                        cont = cont + 1
                        noise_f5_ytime(cont,ntime)=noise_f5_ytime(cont,ntime)+(old_noise_f5_ytime(cont)-noise_f5_ytime(cont,1))
                    end do
                end do

                cont = 0
                do j=1,jy
                    do i=1,jx  !ispon
                        cont = cont + 1
                        noise_f5_ztime(cont,ntime)=noise_f5_ztime(cont,ntime)+(old_noise_f5_ztime(cont)-noise_f5_ztime(cont,1))
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------
        !     noise in space and time in the sponge region

        ntime = ktime -(ktime/nti)*nti
        call creo_xyz_f5(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)

        meanx = 0.
        meany = 0.
        meanz = 0.
        cont_mean = 0
        do k=1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    cont_mean = cont_mean + 1
                    meanx = meanx + coef_x(i,j,k)
                    meany = meany + coef_y(i,j,k)
                    meanz = meanz + coef_z(i,j,k)
                end do
            end do
        end do
        meanx = meanx/real(cont_mean)
        meany = meany/real(cont_mean)
        meanz = meanz/real(cont_mean)
        !      if(myid ==0)write(*,*)'MEAN:',meanx,meany,meanz
        dev_x = 0.
        dev_y = 0.
        dev_z = 0.
        do k=1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    coef_x(i,j,k) = coef_x(i,j,k)-meanx
                    coef_y(i,j,k) = coef_y(i,j,k)-meany
                    coef_z(i,j,k) = coef_z(i,j,k)-meanz

                    dev_x = dev_x + coef_x(i,j,k)*coef_x(i,j,k)
                    dev_y = dev_y + coef_y(i,j,k)*coef_y(i,j,k)
                    dev_z = dev_z + coef_z(i,j,k)*coef_z(i,j,k)
                end do
            end do
        end do

        dev_x = sqrt(dev_x/real(cont_mean))
        dev_y = sqrt(dev_y/real(cont_mean))
        dev_z = sqrt(dev_z/real(cont_mean))
        var_x = 0.05
        var_y = 0.05
        var_z = 0.05
        do k=1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    coef_x(i,j,k) = coef_x(i,j,k)*var_x/dev_x
                    coef_y(i,j,k) = coef_y(i,j,k)*var_y/dev_y
                    coef_z(i,j,k) = coef_z(i,j,k)*var_z/dev_z
                end do
            end do
        end do
        !-----------------------------------------------------------------------
        !     for interpolation in time of tke
        !      dt1 = ti - ti_pom_old
        !      dt2 = ti_pom_new - ti
        !      dt_tot = 1./(dt1+dt2)
        icount = 0
        do k=1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    if(index_out5(i,j) .gt. 0.01)then

                        !fluctuation
                        uprime = u(i,j,k)-up5(i,j)
                        vprime = v(i,j,k)-vp5(i,j)
                        wprime = w(i,j,k)-wp5(i,j)

                        ! coef to scale tke for each direction a1+a2+a3 = 1
                        a_den=1./( up5(i,j)*up5(i,j)+vp5(i,j)*vp5(i,j)+wp5(i,j)*wp5(i,j) )

                        a1 = up5(i,j)*up5(i,j) * a_den
                        a2 = vp5(i,j)*vp5(i,j) * a_den
                        a3 = wp5(i,j)*wp5(i,j) * a_den

                                 !kinetic energy
                        ke_ref =   up5(i,j)*up5(i,j)+vp5(i,j)*vp5(i,j)+wp5(i,j)*wp5(i,j)

                        ke =   u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)

                        if(ke .gt. 1.3*ke_ref)cycle

                        ! left hand side ( TKE(n+1)-TKE(n) ) / Delta T
                        tke = tke5(i,j)

                        tke_new = .5*(  uprime*uprime + vprime*vprime + wprime*wprime )

                        tke = abs((tke_new - tke)*inv_dt)

                        ! fluctuation velocity scale
                        velrms = sqrt(2.*tke_new)

                        !         v_scale = up5(i,j)*up5(i,j) + vp5(i,j)*vp5(i,j) + wp5(i,j)*wp5(i,j)

                        !         v_scale = 0.01*sqrt(v_scale)

                        !         velrms = v_scale

                        segno_u = sign(1.,up5(i,j))
                        segno_v = sign(1.,vp5(i,j))
                        segno_w = sign(1.,wp5(i,j))

                        segno_u = sign(1.,uprime)
                        segno_v = sign(1.,vprime)
                        segno_w = sign(1.,wprime)


                        if(abs(uprime) .lt. 0.1*velrms )then
                            uprime = 0.1*segno_u*velrms
                            icount = icount + 1
                        end if
                        if(abs(vprime) .lt. 0.1*velrms )then
                            vprime = 0.1*segno_v*velrms
                        end if
                        if(abs(wprime) .lt. 0.1*velrms )then
                            wprime = 0.1*segno_w*velrms
                        end if


                        ! bodyforce
                        bcsi_f5(i,j,k)=cx*(a1*tke/abs(uprime))*coef_x(i,j,k)

                        beta_f5(i,j,k)=cy*(a2*tke/abs(vprime))*coef_y(i,j,k)

                        bzet_f5(i,j,k)=cz*(a3*tke/abs(wprime))*coef_z(i,j,k)

                    end if
                end do
            end do
        end do
9200    format(21e18.10)
        write(*,*)'COUNT',icount

        return
    end

    subroutine creo_xyz_f5(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)

        implicit none
        !-----------------------------------------------------------------------
        integer i,j,k,cont,ntime
        integer n1,n2,n3,kparasta,kparaend
        integer myid

        real coef_x(n1,n2,kspon)
        real coef_y(n1,n2,kspon)
        real coef_z(n1,n2,kspon)
        !-----------------------------------------------------------------------
        !     in x
        do i=1,n1   !ispon
            cont = 0
            do k=1,kspon !kparasta,kparaend
                do j=1,n2
                    cont = cont + 1
                    coef_x(i,j,k) = noise_f5_xtime(cont,ntime)+noise_f5_xspace(cont,i)
                end do
            end do
        end do

        !     in y
        do j=1,n2
            cont = 0
            do i=1,n1 !ispon
                do k=1,kspon !kparasta,kparaend
                    cont = cont + 1
                    coef_y(i,j,k) = noise_f5_ytime(cont,ntime)+noise_f5_yspace(cont,j)
                end do
            end do
        end do

        !     in z
        do k=1,kspon !kparasta,kparaend
            cont = 0
            do i=1,n1  !ispon
                do j=1,n2
                    cont = cont + 1
                    coef_z(i,j,k) = noise_f5_ztime(cont,ntime)+noise_f5_zspace(cont,k) !-kparasta+1)
                end do
            end do
        end do

        return
    end

    subroutine disturbance_face6(ti,ktime,cx,cy,cz,nti)

        use contour_module, only: up6,vp6,wp6

        !     generate disturbance for face 6 for nesting procedure

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer ktime,ntime
        integer i,j,k,n
        integer nti
        !      integer kspon
        integer istampo
        integer cont

        integer nreal,nstep

        real corr_factor_sp
        real cx,cy,cz

        real uprime,vprime,wprime

        real coef_x(n1,n2,n3+1-kspon:n3)
        real coef_y(n1,n2,n3+1-kspon:n3)
        real coef_z(n1,n2,n3+1-kspon:n3)

        real ti,dt1,dt2,dt_tot
        real tke,tke_new
        integer cont_mean
        real meanx,meany,meanz
        real dev_x,dev_y,dev_z
        real var_x,var_y,var_z
        real velrms
        real segno_u,segno_v,segno_w
        real v_scale

        real a1,a2,a3,a_den
        real inv_dt

        real ke,ke_ref
        !-----------------------------------------------------------------------
        inv_dt = 1./dt
        !-----------------------------------------------------------------------
        ! allocation and initialization
        if(ktime .eq. 1)then

            corr_factor_sp = 0.

            !       face 5
            !       in x in space and time
            nreal = jy*kspon   !jy*jz/nproc
            nstep = nti
            allocate(noise_f6_xtime(nreal,-1:nstep*2))
            allocate(old_noise_f6_xtime(nreal))
            noise_f6_xtime      = 0.
            old_noise_f6_xtime  = 0.

            nreal = jy*kspon   !jy*jz/nproc
            nstep = jx         !ispon
            allocate(noise_f6_xspace(nreal,-1:nstep*2))
            noise_f6_xspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f6_xspace,istampo,myid)

            !       in y in space and time
            nreal = jx*kspon   !ispon*jz/nproc
            nstep = nti
            allocate(noise_f6_ytime(nreal,-1:nstep*2))
            allocate(old_noise_f6_ytime(nreal))
            noise_f6_ytime     = 0.
            old_noise_f6_ytime = 0.

            nreal = jx*kspon !ispon*jz/nproc
            nstep = jy
            allocate(noise_f6_yspace(nreal,-1:nstep*2))
            noise_f6_yspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f6_yspace,istampo,myid)

            !       in z in space and time
            nreal = jy*jx  !jy*ispon
            nstep = nti
            allocate(noise_f6_ztime(nreal,-1:nstep*2))
            allocate(old_noise_f6_ztime(nreal))
            noise_f6_ztime     = 0.
            old_noise_f6_ztime = 0.

            nreal = jy*jx !jy*ispon
            nstep = kspon !jz/nproc
            allocate(noise_f6_zspace(nreal,-1:nstep*2))
            noise_f6_zspace=0.
            istampo = 0
            call genero_random(corr_factor_sp,nreal,nstep,noise_f6_zspace,istampo,myid)
        end if
        !-----------------------------------------------------------------------
        !     to have a colored noise connected to the previous number generation
        !
        if(ktime .gt.1 .and. ktime .eq. nti*(ktime/nti))then
            ntime = nti
            if(myid.eq.0)write(*,*)'TIME',ktime,ntime
            cont = 0
            do k=1,kspon !kparasta,kparaend
                do j=1,jy
                    cont = cont + 1
                    old_noise_f6_xtime(cont) = noise_f6_xtime(cont,ntime)
                end do
            end do

            cont = 0
            do k=1,kspon !kparasta,kparaend
                do i=1,jx  !ispon
                    cont = cont + 1
                    old_noise_f6_ytime(cont) = noise_f6_ytime(cont,ntime)
                end do
            end do

            cont = 0
            do j=1,jy
                do i=1,jx  !ispon
                    cont = cont + 1
                    old_noise_f6_ztime(cont) = noise_f6_ztime(cont,ntime)
                end do
            end do
        end if


        !-----------------------------------------------------------------------
        !     DISTURBANCE GENERATION IN TIME
        if(ktime .eq. 1 .or. ktime.eq.nti*(ktime/nti))then

            !      corr_factor = 100.

            !     in x in time
            nreal = jy*kspon !jy*jz/nproc
            nstep = nti
            noise_f6_xtime=0.
            call genero_random(corr_factor,nreal,nstep,noise_f6_xtime,istampo,myid)

            !     in y in time
            nreal = jx*kspon !ispon*jz/nproc
            nstep = nti
            noise_f6_ytime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f6_ytime,istampo,myid)

            !     in z in time
            nreal = jy*jx  !jy*ispon
            nstep = nti
            noise_f6_ztime=0.
            istampo = 0
            call genero_random(corr_factor,nreal,nstep,noise_f6_ztime,istampo,myid)

        end if  !if ktime

        !-----------------------------------------------------------------------
        !     shift for disturbance
        if(ktime .gt. 1 .and. ktime .eq. nti*(ktime/nti))then
            do ntime = 1,nti*2
                cont = 0
                do k=1,kspon !kparasta,kparaend
                    do j=1,jy
                        cont = cont + 1
                        noise_f6_xtime(cont,ntime)=noise_f6_xtime(cont,ntime)+(old_noise_f6_xtime(cont)-noise_f6_xtime(cont,1))
                    end do
                end do

                cont = 0
                do k=1,kspon !kparasta,kparaend
                    do i=1,jx    !ispon
                        cont = cont + 1
                        noise_f6_ytime(cont,ntime)=noise_f6_ytime(cont,ntime)+(old_noise_f6_ytime(cont)-noise_f6_ytime(cont,1))
                    end do
                end do

                cont = 0
                do j=1,jy
                    do i=1,jx  !ispon
                        cont = cont + 1
                        noise_f6_ztime(cont,ntime)=noise_f6_ztime(cont,ntime)+(old_noise_f6_ztime(cont)-noise_f6_ztime(cont,1))
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------
        !     noise in space and time in the sponge region

        ntime = ktime -(ktime/nti)*nti
        call creo_xyz_f6(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)

        meanx = 0.
        meany = 0.
        meanz = 0.
        cont_mean = 0
        do k=jz+1-kspon,jz !1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    cont_mean = cont_mean + 1
                    meanx = meanx + coef_x(i,j,k)
                    meany = meany + coef_y(i,j,k)
                    meanz = meanz + coef_z(i,j,k)
                end do
            end do
        end do
        meanx = meanx/real(cont_mean)
        meany = meany/real(cont_mean)
        meanz = meanz/real(cont_mean)
        !      if(myid ==0)write(*,*)'MEAN:',meanx,meany,meanz
        dev_x = 0.
        dev_y = 0.
        dev_z = 0.
        do k=jz+1-kspon,jz !1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    coef_x(i,j,k) = coef_x(i,j,k)-meanx
                    coef_y(i,j,k) = coef_y(i,j,k)-meany
                    coef_z(i,j,k) = coef_z(i,j,k)-meanz

                    dev_x = dev_x + coef_x(i,j,k)*coef_x(i,j,k)
                    dev_y = dev_y + coef_y(i,j,k)*coef_y(i,j,k)
                    dev_z = dev_z + coef_z(i,j,k)*coef_z(i,j,k)
                end do
            end do
        end do

        dev_x = sqrt(dev_x/real(cont_mean))
        dev_y = sqrt(dev_y/real(cont_mean))
        dev_z = sqrt(dev_z/real(cont_mean))
        var_x = 0.05
        var_y = 0.05
        var_z = 0.05
        do k=jz+1-kspon,jz !1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    coef_x(i,j,k) = coef_x(i,j,k)*var_x/dev_x
                    coef_y(i,j,k) = coef_y(i,j,k)*var_y/dev_y
                    coef_z(i,j,k) = coef_z(i,j,k)*var_z/dev_z
                end do
            end do
        end do
        !-----------------------------------------------------------------------
        !     for interpolation in time of tke
        !      dt1 = ti - ti_pom_old
        !      dt2 = ti_pom_new - ti
        !      dt_tot = 1./(dt1+dt2)

        do k=jz+1-kspon,jz !1,kspon !kparasta,kparaend
            do j=1,jy
                do i=1,jx    !ispon
                    if(index_out6(i,j) .gt. 0.01)then

                        !fluctuation
                        uprime = u(i,j,k)-up6(i,j)
                        vprime = v(i,j,k)-vp6(i,j)
                        wprime = w(i,j,k)-wp6(i,j)

                        ! coef to scale tke for each direction a1+a2+a3 = 1
                        a_den=1./( up6(i,j)*up6(i,j)+vp6(i,j)*vp6(i,j)+wp6(i,j)*wp6(i,j) )

                        a1 = up6(i,j)*up6(i,j) * a_den
                        a2 = vp6(i,j)*vp6(i,j) * a_den
                        a3 = wp6(i,j)*wp6(i,j) * a_den

                                 !kinetic energy
                        ke_ref =   up6(i,j)*up6(i,j)+vp6(i,j)*vp6(i,j)+wp6(i,j)*wp6(i,j)

                        ke =   u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)

                        if(ke .gt. 1.3*ke_ref)cycle

                        ! left hand side ( TKE(n+1)-TKE(n) ) / Delta T
                        tke = tke6(i,j)

                        tke_new = .5*(  uprime*uprime + vprime*vprime + wprime*wprime )

                        tke = abs((tke_new - tke)*inv_dt)

                        ! fluctuation velocity scale
                        velrms = sqrt(2.*tke_new)

                        !         v_scale = up6(i,j)*up6(i,j) + vp6(i,j)*vp6(i,j) + wp6(i,j)*wp6(i,j)

                        !         v_scale = 0.01*sqrt(v_scale)

                        !         velrms = v_scale

                        segno_u = sign(1.,up6(i,j))
                        segno_v = sign(1.,vp6(i,j))
                        segno_w = sign(1.,wp6(i,j))

                        segno_u = sign(1.,uprime)
                        segno_v = sign(1.,vprime)
                        segno_w = sign(1.,wprime)


                        if(abs(uprime) .lt. 0.1*velrms )then
                            uprime = 0.1*segno_u*velrms
                        end if
                        if(abs(vprime) .lt. 0.1*velrms )then
                            vprime = 0.1*segno_v*velrms
                        end if
                        if(abs(wprime) .lt. 0.1*velrms )then
                            wprime = 0.1*segno_w*velrms
                        end if


                        ! bodyforce
                        bcsi_f6(i,j,k)=cx*(a1*tke/abs(uprime))*coef_x(i,j,k)

                        beta_f6(i,j,k)=cy*(a2*tke/abs(vprime))*coef_y(i,j,k)

                        bzet_f6(i,j,k)=cz*(a3*tke/abs(wprime))*coef_z(i,j,k)


                    end if
                end do
            end do
        end do
9200    format(21e18.10)

        return
    end

    subroutine creo_xyz_f6(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)

        implicit none
        !-----------------------------------------------------------------------
        integer i,j,k,cont,ntime
        integer n1,n2,n3,kparasta,kparaend
        integer myid

        real coef_x(n1,n2,n3+1-kspon:n3)
        real coef_y(n1,n2,n3+1-kspon:n3)
        real coef_z(n1,n2,n3+1-kspon:n3)
        !-----------------------------------------------------------------------
        !     in x
        do i=1,n1   !ispon
            cont = 0
            do k=n3+1-kspon,n3 !1,kspon !kparasta,kparaend
                do j=1,n2
                    cont = cont + 1
                    coef_x(i,j,k) = noise_f6_xtime(cont,ntime)+noise_f6_xspace(cont,i)
                end do
            end do
        end do

        !     in y
        do j=1,n2
            cont = 0
            do i=1,n1 !ispon
                do k=n3+1-kspon,n3 !1,kspon !kparasta,kparaend
                    cont = cont + 1
                    coef_y(i,j,k) = noise_f6_ytime(cont,ntime)+noise_f6_yspace(cont,j)
                end do
            end do
        end do

        !     in z
        do k=n3+1-kspon,n3 !1,kspon !kparasta,kparaend
            cont = 0
            do i=1,n1  !ispon
                do j=1,n2
                    cont = cont + 1
                    coef_z(i,j,k) = noise_f6_ztime(cont,ntime)+noise_f6_zspace(cont,k-(n3-kspon)) !kparasta+1)
                end do
            end do
        end do

        return
    end

    subroutine genero_random(corr_factor,nreal,nstep,disturbance,istampo,myid)
        !***********************************************************************
        !
        implicit none
        ! program to illustrate the colored Gaussian Noise generator cgaussA
        ! The routine must be initialized with cgaus0A and calls a flat distribution
        ! random number generator available with most compilers or you can write your
        ! own. Here we used the routine RAN1 from Numerical Recipes 2nd Edition, by
        ! Press, Teukolsky, Vetterling, and Flannery.

        ! It now uses the F90 intrinsic subroutine RANDOM_NUMBER.

        ! The White Guassian noise generator GASDEV from Numerical Recipes was
        ! adapted to produce Colored Gaussian noise. The basic equations for this
        ! computation are presented in the article by
        ! Fox et al., Physical Review A vol.38(1988) page 5938.
        ! This code was [originally] compiled and tested with Microsoft Powerstation.

        ! It was modified by Walt Brainerd to be standard Fortran and
        ! compiled on NAGWare F90.

        integer i,j
        integer istampo,myid
        integer nreal,nstep,npts,idly
        real, allocatable :: eps(:,:),sum(:)
        real smean
        real deviation,variance
        real std,mean,dt !cgaussA  -> Alessandro
        real cortim
        real corr_factor
        real disturbance(nreal,-1:nstep*2)
        real start_time,end_time
        ! get input parameters (typical values shown)
        !c        open(1,file='cnoise.dat')
        !c        read(1,*)nreal             !number of realizations=1000
        !c        read(1,*)nstep             !max delay in corr. func=10
        !c        read(1,*)dt                !time step size=.5
        !c        read(1,*)cortim            !corr. time in the same units as DT=5

        !         nreal = 10
        !     nstep = 1000
        dt = 1.
        !     corr_factor = 0.
        cortim = corr_factor*dt

        allocate(eps(nreal,-1:nstep*2))
        allocate(sum(nreal))
        eps = 0.
        sum = 0.

        ! initialize
        !         call cpu_time(start_time)
        call cgaus0A(dt,cortim)
        !     call cpu_time(end_time)
        !     write(*,*)'TIME'
        !     write(*,'(1e18.10)')end_time-start_time

        ! store several series of Gaussian noise values in array EPS.
        !        call cpu_time(start_time)

        do i=1,nreal
            mean = 0.
            do j=1,nstep*2 !0,nstep*2
                eps(i,j) = cgaussA()
                mean = mean + eps(i,j)
            end do
            mean = mean/float(nstep*2)  !float(nstep*2+1)

            do j=1,nstep*2 !0,nstep*2
                eps(i,j) = eps(i,j) - mean
            end do
        end do
        mean = 0.


        do i=1,nreal             !realizations
            sum(i) = 0.
            do j=1,nstep*2  !0,nstep*2          !time delays
                !cccc          eps(i,j)=cgaussA()
                !          write(100+i,*)j,eps(i,j)
                sum(i) = sum(i) + eps(i,j)
            enddo
            sum(i)=sum(i)/float(nstep*2)   !float(nstep*2+1)

            deviation = 0.
            do j=1,nstep*2  !0,nstep*2          !time delays
                deviation = deviation + (eps(i,j)-sum(i))*(eps(i,j)-sum(i))
            enddo
            deviation = deviation/float(nstep*2) !float(nstep*2+1)
            deviation = sqrt(deviation)

            variance = 1
            variance = variance/deviation
            do j=1,nstep*2 !0,nstep*2
                eps(i,j)=eps(i,j)*variance
            end do
        enddo


        ! calculate the autocorrelation function in variable MEAN.
        !        call cpu_time(start_time)
        npts=nstep*nreal
        do idly=1,nstep !0,nstep
            mean=0.
            std=0.
            do i=1,nreal
                do j=1,nstep !0,nstep
                    mean=mean+real(eps(i,j)*eps(i,j+idly))
                enddo
            enddo
            mean=mean/real(npts)
            smean=sngl(mean)          !single precision speeds up calculations

            ! calculate the error in autocorrelation function in variable STD.
            do i=1,nreal
                do j=1,nstep !0,nstep
                    std=std+real((eps(i,j)*eps(i,j+idly)-smean)**2.)
                enddo
            enddo
            std=sqrt(std)/real(npts) !dble(npts-1.)
            write(200+myid,*)idly,mean,std            !output results
        enddo
        !    call cpu_time(end_time)
        !    write(*,*)'TIME'
        !    write(*,'(1e18.10)')end_time-start_time


        !       storage
        do i=1,nreal
            do j=0,nstep*2
                disturbance(i,j) = eps(i,j)
            end do
        end do


        deallocate(eps)
        deallocate(sum)
        !        end
        return
    end subroutine genero_random

    subroutine cgaus0A(dt,cortim)
        !==========================================================================
        ! initialize the RNG's
        ! and set the color of gaussian noise
        ! DT is the time step used in whatever process the colored Gaussian noise
        !   is used.
        ! CORTIM is correlation time in the same units as time step DT.
        ! WHITE=.true. means generate white gaussian noise which happens when
        !   CORTIM=0. This flag is used in cgaussA.
        ! Here we use the flat distribution RAN1 also taken from Numerical Recipe
        ! but any other good flat distribution random number generator will do.

        !        double precision ran1,cape,dt,l1me2,cgaussA
        real cape,dt,l1me2 !,cgaussA -> Alessandro
        real cortim,x
        logical white
        common /color/l1me2,cape,white
        if(cortim.eq.0.)then
            white=.true.
            l1me2=-2.000                        !white noise
            cape=0.0
        else
            white=.false.
            cape=exp(-dt/real(cortim))
            !parameter needed in cgaussA
            l1me2=-(real(1.)-cape*cape)*real(2./cortim)
        endif
        !        idum=-1
        !        x=ran1(idum)            !initialize flat rng
        x=cgaussA()            !initialize cgaussA value
        return
    end subroutine cgaus0A

    real function cgaussA()

        !==========================================================================
        ! Program to produce exponentially correlated colored (Gaussian) noise.
        ! based on Fox et al Physical Review A vol.38(1988)5938 and
        ! modification of GASDEV from Numerical Recipes for Fortran(2nd ed.pg279)

        ! CAPE is capital E in the article by Fox et. al.
        ! PREV is the previous value of cgaussA used in the next iteration
        ! L1ME2 is the main parameters causing colored noise in Fox et al
        !       and represents (lamda*(1-E)^2). Ditto for H in that article.

        ! routine is illustrated in Double Precision in case it is needed in this
        ! mode, otherwise all Double Precision variables maybe changed to REAL
        ! but the corresponding changes must be made to cgaus0A and the calling
        ! programs.


        Implicit none

        !      INTEGER idum,iset
        INTEGER iset
        logical white
        !      double precision  fac,gset,rsq,v1,v2,ran1,l1me2,h,cape
        real  fac,gset,rsq,v1,v2,l1me2,h,cape,prev
        common /color/l1me2,cape,white

        SAVE iset,gset,prev
        DATA iset/0/
        DATA prev/0.0d0/

        if (iset.eq.0) then
            !1       v1=2.*ran1(idum)-1.
1           call random_number(v1)
            v1=2.*v1-1
            !        v2=2.*ran1(idum)-1.
            call random_number(v2)
            v2=2.*v2-1
            rsq=v1**2.+v2**2.
            if(rsq>=1. .or. rsq==0.)goto 1
            !took out sqrt(2) vs eq(28) Fox etal
            fac=sqrt(l1me2*log(rsq)/rsq)
            gset=v1*fac
            h=v2*fac
            iset=1
        else
            h=gset
            iset=0
        endif

        if(white)then  !please note that the time step vs its sqrt
            cgaussA=h      !in integration is previously set in PARAM
        else
            cgaussA=prev*cape+h
            prev=cgaussA
        endif

        return
    end function cgaussA

end module buffer_bodyforce_module
