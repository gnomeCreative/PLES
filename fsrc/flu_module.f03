module flu_module

    use mysending
    use myarrays_metri3
    use myarrays_velo3
    !use wallmodel_module, only: att_mod_par, eseguo34, u_t, utangente
    use mysettings, only: bodyforce,insc
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    !real :: :: bulk

contains

    subroutine flucn(r,espl,coef_wall,tipo,tau_wind)
        ! compute explicit diffusive term
        !
        ! NNI*G11*D(u,v,w)/D(csi)+
        ! NNI*G22*D(u,v,w)/D(eta)+
        ! NNI*G33*D(u,v,w)/D(zita)
        !
        !-----------------------------------------------------------------------

        use inflow_module, only: areola3,areola4
        use wallmodel_module, only: att_mod_par,eseguo34,u_t,utangente

        implicit none

        !-----------------------------------------------------------------------
        real,intent(in) :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        real,optional,intent(in) :: tau_wind(0:,kparasta-1:)
        integer,intent(in) :: coef_wall
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer,intent(in) :: espl
        !-----------------------------------------------------------------------
        integer :: ierr,status(MPI_STATUS_SIZE)
        integer :: kparastal,kparaendl
        integer :: i,j,k
        real :: rhow
        real :: es
        logical,allocatable :: prov(:,:,:)
        !-----------------------------------------------------------------------

        ! explicit integration requires an additional term, f1=f1+...
        ! thi is enforced by using the following switch
        es=real(espl)


        ! term nni*g11*d/d(csi)
        !
        if (ip==1) then
            !
            !     side left
            do k=kparasta,kparaend
                do j=1,n2
                    f1(0,j,k)=es*f1(0,j,k)+annit(0,j,k)*g11(0,j,k)*(-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
                end do
            end do
            !
            !     side right
            do k=kparasta,kparaend
                do j=1,n2
                    f1(n1,j,k)=es*f1(n1,j,k)+annit(n1+1,j,k)*g11(n1,j,k)*(8.*r(n1+1,j,k)-9.*r(n1,j,k)+r(n1-1,j,k))/3.
                end do
            end do
        !
        end if
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !
                    f1(i,j,k)=es*f1(i,j,k)+.5*(annit(i,j,k)+annit(i+1,j,k))*g11(i,j,k)*(r(i+1,j,k)-r(i,j,k))
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !
                    if (bodyforce) then
                        if (tipo(i,j,k)==0) f1(i,j,k)=0.
                        if (i<n1) then
                            if (tipo(i+1,j,k)==0) f1(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do

        !
        !-----------------------------------------------------------------------
        ! term nni*g22*d/d(eta)
        !
        ! direction 2 is always not periodic
            !
            !     side bottom
            !
            !     for direction 2, allthough coef_wall=1, the wall model is switched off
            if (eseguo34==0.and.coef_wall==1) then
                allocate(prov(n1,2,kparasta:kparaend))
                prov(:,:,:)=att_mod_par(:,:,:)
                att_mod_par(:,:,:)=.false.
            end if

            !     wall model off
            if (coef_wall==0) then
                do k=kparasta,kparaend
                    do i=1,n1
                        f2(i,0,k)=es*f2(i,0,k)+annitV(i,0,k)*g22(i,0,k)*(-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                    end do
                end do
            !     wall model on
            else if (coef_wall==1) then
                do k=kparasta,kparaend
                    do i=1,n1
                        if (att_mod_par(i,1,k)) then
                            f2(i,0,k)=es*f2(i,0,k)+u_t(i,1,k)*u_t(i,1,k)*areola3(i,k)* &
                                ((2-eseguo34)*u(i,1,k)-(1-eseguo34)*w(i,1,k))/utangente(i,1,k)
                        else
                            f2(i,0,k)=es*f2(i,0,k)+annitV(i,0,k)*g22(i,0,k)*(-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                        end if

                    end do
                end do
            end if

            !         if (coef_wall==1 .and. imoist==1 .and. eseguo34/=0) then
            !            call readtau(ti,r)
            !         end if

            !
            ! side upper
            !
            if (present(tau_wind)) then

                ! Giulia modificavento:
                rhow=1000.  ! SMELLS LIKE MAGIC
                do k=kparasta,kparaend
                    do i=1,n1
                        !          f2(i,jy,k)=abs(vel_tau(i,k))*vel_tau(i,k)*areola4(i,k)
                        f2(i,n2,k)=tau_wind(i,k)*areola4(i,k)/rhow
                    end do
                end do
            ! Giulia modificavento:

            else

                !     wall model off
                if (coef_wall==0) then
                    do k=kparasta,kparaend
                        do i=1,n1
                            f2(i,n2,k)=es*f2(i,n2,k)+annitV(i,n2+1,k)*g22(i,n2,k)*&
                                (8.*r(i,n2+1,k)-9.*r(i,n2,k)+r(i,n2-1,k))/3.
                        end do
                    end do

                !      wall model on
                elseif (coef_wall==1) then
                    do k=kparasta,kparaend
                        do i=1,n1
                            if (att_mod_par(i,2,k)) then
                                f2(i,n2,k)=es*f2(i,n2,k)-u_t(i,2,k)*u_t(i,2,k)*areola4(i,k)*&
                                    ((2-eseguo34)*u(i,n2,k) &
                                    -(1-eseguo34)*w(i,n2,k))/utangente(i,2,k)
                            else
                                f2(i,n2,k)=es*f2(i,n2,k)+annitV(i,n2+1,k)*g22(i,n2,k)*&
                                    (8.*r(i,n2+1,k)-9.*r(i,n2,k)+r(i,n2-1,k))/3.

                            end if
                        end do
                    end do
                end if



            end if  !end if windyes
            !
                    ! wall model on but switched for direction2
            if (eseguo34==0.and.coef_wall==1) then
                att_mod_par(:,:,:)=prov(:,:,:)
                deallocate(prov)
            end if
        !
        !     into the field
        do k=kparasta,kparaend
            do i=1,n1
                do j=1,n2-1
                    !
                    f2(i,j,k)=es*f2(i,j,k)+.5*(annitV(i,j,k)+annitV(i,j+1,k))&
                        &              *g22(i,j,k)*(r(i,j+1,k)-r(i,j,k))
                !
                end do
            end do
        end do
        !
        do k=kparasta,kparaend
            do j=1,n2-1
                do i=1,n1
                    !

                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f2(i,j,k)=0.
                        if (j<n2) then
                            if (tipo(i,j+1,k)==0)f2(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! term nni*g33d/d(zita)
        !
        if (kp==1) then
            !
            !     side back
            if (myid==0) then

                do i=1,n1
                    do j=1,n2
                        f3(i,j,0)=es*f3(i,j,0)+annit(i,j,0)*g33(i,j,0)*&
                            &     (-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
                    end do
                end do

            end if
            !
            !     side front
            if (myid==nproc-1) then

                do i=1,n1
                    do j=1,n2
                        f3(i,j,n3)=es*f3(i,j,n3)+annit(i,j,n3+1)*g33(i,j,n3)*&
                            &   (8.*r(i,j,n3+1)-9.*r(i,j,n3)+r(i,j,n3-1))/3.
                    end do
                end do

            end if
        !
        end if
        !
        !     into the field
        if (myid==0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid==nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else
            kparastal=kparasta
            kparaendl=kparaend
        end if

        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    !
                    f3(i,j,k)=es*f3(i,j,k)+.5*(annit(i,j,k)+annit(i,j,k+1))*g33(i,j,k)*(r(i,j,k+1)-r(i,j,k))
                !
                end do
            end do
        end do

        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    if (bodyforce) then
                        if (tipo(i,j,k)==0) f3(i,j,k)=0.
                        if (k<n3) then
                            if (tipo(i,j,k+1)==0) f3(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do
        !
        ! pass f3 at the border to all procs
        !
        if (myid==0) then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if (myid==nproc-1) then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else
            leftpem=leftpe
            rightpem=rightpe
        end if

        if (rightpem/=MPI_PROC_NULL) then
            call MPI_SEND(f3(1,1,kparaend),n1*n2,MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem/=MPI_PROC_NULL) then
            call MPI_RECV(f3(1,1,kparasta-1),n1*n2,MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if

    !        !     integral
    !        bulk_loc=0.
    !        !
    !        do k=kparasta,kparaend
    !            do j=1,jy
    !                do i=1,jx
    !                    bulk_loc=bulk_loc+f3(i,j,k)-f3(i  ,j  ,k-1)+&
    !                        &                  f2(i,j,k)-f2(i  ,j-1,k  )+&
    !                        &                  f1(i,j,k)-f1(i-1,j  ,k  )
    !                end do
    !            end do
    !        end do
    !
    !        ! make the value known to all procs
    !
    !        call MPI_ALLREDUCE(bulk_loc,bulkn,1,MPI_REAL_SD,&
    !            &                   MPI_SUM,&
    !            &                   MPI_COMM_WORLD,ierr)
    !
    !        bulk=bulk+bulkn

    end subroutine flucn

    subroutine flu_turbo()
        !***********************************************************************
        ! compute explicit terms for turbulence model like
        ! d/dcsi(Uf*annit), periodic version
        !
        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer :: i,j,k
        !
        integer :: kparaendp

        !-----------------------------------------------------------------------
        ! compute controvariant flux Uf (jordan notation)
        !
        !     side left and right, periodic
        do k=kparasta,kparaend
            do j=1,n2
                !
                cgra1(0,j,k)=csx(0,j,k)*.5*(gra1(n1,j,k)+gra1(1,j,k))+ &
                    csy(0,j,k)*.5*(gra2(n1,j,k)+gra2(1,j,k))+ &
                    csz(0,j,k)*.5*(gra3(n1,j,k)+gra3(1,j,k))- &
                    cgra1(0,j,k)
                cgra1(0,j,k)=(1-ip)*cgra1(0,j,k)

                cgra1(n1,j,k)=csx(n1,j,k)*.5*(gra1(n1,j,k)+gra1(1,j,k))+ &
                    csy(n1,j,k)*.5*(gra2(n1,j,k)+gra2(1,j,k))+ &
                    csz(n1,j,k)*.5*(gra3(n1,j,k)+gra3(1,j,k))- &
                    cgra1(n1,j,k)
                cgra1(n1,j,k)=(1-ip)*cgra1(n1,j,k)
            !
            end do
        end do
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,n2
                do i=1,n1-1
                    !
                    cgra1(i,j,k)=csx(i,j,k)*.5*(gra1(i+1,j,k)+gra1(i,j,k))+ &
                        csy(i,j,k)*.5*(gra2(i+1,j,k)+gra2(i,j,k))+ &
                        csz(i,j,k)*.5*(gra3(i+1,j,k)+gra3(i,j,k))- &
                        cgra1(i,j,k)
                !
                end do
            end do
        end do
        !
        ! compute controvariant flux Vf (jordan notation)
        !
        !     sides bottom and upper, not periodic
        do k=kparasta,kparaend
            do i=1,n1
                cgra2(i,0,k)=0.0
                cgra2(i,n2,k)=0.0
            end do
        end do

        !     into the field
        do k=kparasta,kparaend
            do j=1,n2-1
                do i=1,n1
                    !
                    cgra2(i,j,k)=etx(i,j,k)*.5*(gra1(i,j+1,k)+gra1(i,j,k))+ &
                        ety(i,j,k)*.5*(gra2(i,j+1,k)+gra2(i,j,k))+ &
                        etz(i,j,k)*.5*(gra3(i,j+1,k)+gra3(i,j,k))- &
                        cgra2(i,j,k)
                !
                end do
            end do
        end do
        !
        ! compute controvariant flux Wf (jordan notation)
        !
        do j=1,n2
            do i=1,n1
                !
                if (myid==0) then

                    cgra3(i,j,0 )=ztx(i,j,0 ) *.5*(gra1(i,j,1)+gra1_appoggio(i,j,n3))+ &
                        zty(i,j,0 ) *.5*(gra2(i,j,1)+gra2_appoggio(i,j,n3))+ &
                        ztz(i,j,0 ) *.5*(gra3(i,j,1)+gra3_appoggio(i,j,n3))- &
                        cgra3(i,j,0)
                    cgra3(i,j,0 )=(1-kp)*cgra3(i,j,0 )

                else if (myid==nproc-1) then

                    cgra3(i,j,n3)=ztx(i,j,n3) *.5*(gra1_appoggio(i,j,1)+gra1(i,j,n3))+ &
                        zty(i,j,n3) *.5*(gra2_appoggio(i,j,1)+gra2(i,j,n3))+ &
                        ztz(i,j,n3) *.5*(gra3_appoggio(i,j,1)+gra3(i,j,n3))- &
                        cgra3(i,j,n3)
                    cgra3(i,j,n3)=(1-kp)*cgra3(i,j,n3)

                end if
            !
            end do
        end do
        !
        !     into the field
        if (myid==nproc-1) then
            kparaendp=kparaend-1
        else
            kparaendp=kparaend
        end if
        do k=kparasta,kparaendp
            do j=1,n2
                do i=1,n1
                    !
                    cgra3(i,j,k)=ztx(i,j,k)*.5*(gra1(i,j,k+1)+gra1(i,j,k))+ &
                        zty(i,j,k)*.5*(gra2(i,j,k+1)+gra2(i,j,k))+ &
                        ztz(i,j,k)*.5*(gra3(i,j,k+1)+gra3(i,j,k))- &
                        cgra3(i,j,k)
                !
                end do
            end do
        end do

        !hicco ??? come mai questo vecchio commento senza poi il passaggio tra procs
        ! rendo visibili cgra1 cgra2 cgra3 a tutti i PEs

        return
    end subroutine flu_turbo

    subroutine flucrho(r,isc,tipo)
        !***********************************************************************
        ! compute diffusive explicit terms
        !
        ! k*g11*D(rho)/D(csi) +
        ! k*g22*D(rho)/D(eta) +
        ! k*g33*D(rho)/D(zita)
        !
        ! the subroutine is for scalar eq.
        !
        !use mysettings, only: coef_wall
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer :: ierr,status(MPI_STATUS_SIZE)
        integer :: kparastal,kparaendl
        integer :: i,j,k,ii,jj,kk,isc

        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real :: akaptV(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real :: akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)

        !real :: f2_loc,f2_tot,f2_mean

        ! giulia aggiungo questi per eliminare i flussi tra i solidi
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! term nni*g11*d/d(csi)
        !
        do ii=1,ip
            !
            !hicco usare un solo ciclo!!!!

            !     side left
            do k=kparasta,kparaend
                do j=1,n2
                    f1(0,j,k)=akapt(isc,0,j,k)*g11(0,j,k)* (-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
                end do
            end do
            !
            !     side right
            do k=kparasta,kparaend
                do j=1,n2
                    f1(n1,j,k)=akapt(isc,n1+1,j,k)*g11(n1,j,k)* (8.*r(n1+1,j,k)-9.*r(n1,j,k)+r(n1-1,j,k))/3.
                end do
            end do
        !
        end do
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !
                    f1(i,j,k)=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)* (r(i+1,j,k)-r(i,j,k))
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !

                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f1(i,j,k)=0.
                        if (i<n1) then
                            if (tipo(i+1,j,k)==0)f1(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do


        !
        !-----------------------------------------------------------------------
        ! term nni*g22*d/d(eta)
        !
        ! direction 2 is always not periodic
            !
            !      if (coef_wall==1 .and. imoist==1) then
            !         call readflux(ti,isc)
            !      else

            !     side bottom
            do k=kparasta,kparaend
                do i=1,n1
                    f2(i,0,k)=akaptV(isc,i,0,k)*g22(i,0,k)* (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                end do
            end do

            !    end if
            !
            !     side upper
            do k=kparasta,kparaend
                do i=1,n1
                    f2(i,n2,k)=akaptV(isc,i,n2+1,k)*g22(i,n2,k)* (8.*r(i,n2+1,k)-9.*r(i,n2,k)+r(i,n2-1,k))/3.
                end do
            end do

        !hicco      end if
        !
        !
        !     into the field
        do k=kparasta,kparaend
            do i=1,n1
                do j=1,n2-1
                    !
                    f2(i,j,k)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)* (r(i,j+1,k)-r(i,j,k))
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2-1
                do i=1,n1
                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f2(i,j,k)=0.
                        if (j<n2) then
                            if (tipo(i,j+1,k)==0)f2(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do



        !
        !-----------------------------------------------------------------------
        ! term nni*g33d/d(zita)
        !
        do kk=1,kp
            !
            !     side back
            if (myid==0) then

                do i=1,n1
                    do j=1,n2
                        f3(i,j,0)=akapt(isc,i,j,0)*g33(i,j,0)*(-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
                    end do
                end do
            !
            !     side front
            else if (myid==nproc-1) then

                do i=1,n1
                    do j=1,n2
                        f3(i,j,n3)=akapt(isc,i,j,n3+1)*g33(i,j,n3)*(8.*r(i,j,n3+1)-9.*r(i,j,n3)+r(i,j,n3-1))/3.
                    end do
                end do

            end if
        !
        end do
        !
        !     into the field
        if (myid==0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid==nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else if ((myid/=0).and.(myid/=nproc-1)) then
            kparastal=kparasta
            kparaendl=kparaend
        end if


        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    !
                    f3(i,j,k)=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))*g33(i,j,k)* (r(i,j,k+1)-r(i,j,k))
                !
                end do
            end do
        end do

        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f3(i,j,k)=0.
                        if (k<n3) then
                            if (tipo(i,j,k+1)==0)f3(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do
        !
        ! make border values f1, f2, f3 known to all procs

        call MPI_SENDRECV(f3(1,1,kparaend),n1*n2, &
            MPI_REAL_SD,rightpe,tagrs, &
            f3(1,1,kparasta-1),n1*n2, &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)

        !-----------------------------------------------------------------------
        !   check for the fluxes
        !    if (imoist==1) then
        !        !     Compute the mean fluxes at bottom
        !        f2_loc=0.
        !        f2_tot=0.
        !        f2_mean=0.
        !        do k=kparasta,kparaend
        !            do i=1,jx
        !                f2_loc=f2_loc+f2(i,0,k)
        !            end do
        !        end do
        !        !     sum on f2_loc
        !        call MPI_ALLREDUCE(f2_loc,f2_tot,1,MPI_REAL_SD, &
        !            MPI_SUM,MPI_COMM_WORLD,ierr)
        !        f2_mean=f2_tot/real(jx*jz)
        !
        !        if (myid==0) then
        !            write(*,*)'scalar ',isc,' f2 mean bot: ',f2_mean
        !        end if
        !
        !        !     Compute the mean fluxes at top
        !        f2_loc=0.
        !        f2_tot=0.
        !        f2_mean=0.
        !        do k=kparasta,kparaend
        !            do i=1,jx
        !                f2_loc=f2_loc+f2(i,jy,k)
        !            end do
        !        end do
        !        !     sum on f2_loc
        !        call MPI_ALLREDUCE(f2_loc,f2_tot,1,MPI_REAL_SD, &
        !            MPI_SUM,MPI_COMM_WORLD,ierr)
        !        f2_mean=f2_tot/real(jx*jz)
        !
        !        if (myid==0) then
        !            write(*,*)'scalar ',isc,' f2 mean top: ',f2_mean
        !        end if
        !
        !    end if

        return
    end subroutine flucrho

    subroutine flucrhoesp(r,isc,tipo)
        !***********************************************************************
        ! compute explicit diffusive term
        !
        ! k*g11*D(rho)/D(csi) +
        ! k*g22*D(rho)/D(eta) +
        ! k*g33*D(rho)/D(zita)
        !
        ! the subroutine is for scalar eq.
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer :: ierr,status(MPI_STATUS_SIZE)
        integer :: kparastal,kparaendl
        integer :: i,j,k,ii,jj,kk,isc
        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)

        ! giulia aggiungo questi per eliminare i flussi tra i solidi
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! term nni*g11*d/d(csi)
        !
        do ii=1,ip
            !hicco perche' usa due cicli invece di uno???? ridurre
            ! analogamente per sotto
            !
            ! side left

            do k=kparasta,kparaend
                do j=1,n2
                    f1(0,j,k)=f1(0,j,k)+akapt(isc,0,j,k)*g11(0,j,k)*(-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
                end do
            end do
            !
            ! side right
            !
            do k=kparasta,kparaend
                do j=1,n2
                    f1(n1,j,k)=f1(n1,j,k)+akapt(isc,n1+1,j,k)*g11(n1,j,k)*(8.*r(n1+1,j,k)-9.*r(n1,j,k)+r(n1-1,j,k))/3.
                end do
            end do
        !
        end do
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !
                    f1(i,j,k)=f1(i,j,k)+.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)*(r(i+1,j,k)-r(i,j,k))
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f1(i,j,k)=0.
                        if (i<n1) then
                            if (tipo(i+1,j,k)==0)f1(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do

        !
        !-----------------------------------------------------------------------
        ! term nni*g22*d/d(eta)
        !
        ! direction 2 is always not periodic
            !
            ! side bottom
            !
            do k=kparasta,kparaend
                do i=1,n1
                    f2(i,0,k)=f2(i,0,k)+akaptV(isc,i,0,k)*g22(i,0,k)*(-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                end do
            end do
            !
            ! side upper
            !
            do k=kparasta,kparaend
                do i=1,n1
                    f2(i,n2,k)=f2(i,n2,k)+akaptV(isc,i,n2+1,k)*g22(i,n2,k)*(8.*r(i,n2+1,k)-9.*r(i,n2,k)+r(i,n2-1,k))/3.
                end do
            end do
        !
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do i=1,n1
                do j=1,n2-1
                    !
                    f2(i,j,k)=f2(i,j,k)+.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)*(r(i,j+1,k)-r(i,j,k))
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2-1
                do i=1,n1
                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f2(i,j,k)=0.
                        if (j<n2) then
                            if (tipo(i,j+1,k)==0)f2(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! term nni*g33d/d(zita)
        !
        do kk=1,kp
            !
            ! side back
            !
            if (myid==0) then

                do i=1,n1
                    do j=1,n2
                        f3(i,j,0)=f3(i,j,0)+akapt(isc,i,j,0)*g33(i,j,0)*(-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
                    end do
                end do

            ! side front
            !
            else if (myid==nproc-1) then

                do i=1,n1
                    do j=1,n2
                        f3(i,j,n3)=f3(i,j,n3)+akapt(isc,i,j,n3+1)*g33(i,j,n3)*(8.*r(i,j,n3+1)-9.*r(i,j,n3)+r(i,j,n3-1))/3.
                    end do
                end do

            end if
        !
        end do
        !
        ! into the field
        !
        if (myid==0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid==nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else if ((myid/=0).and.(myid/=nproc-1)) then
            kparastal=kparasta
            kparaendl=kparaend
        end if


        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    !
                    f3(i,j,k)=f3(i,j,k)+.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))*g33(i,j,k)*(r(i,j,k+1)-r(i,j,k))
                !
                end do
            end do
        end do

        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f3(i,j,k)=0.
                        if (k<n3) then
                            if (tipo(i,j,k+1)==0)f3(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do


        !
        ! make border values f1, f2, f3 known to closer procs
        !
        call MPI_SENDRECV(f3(1,1,kparaend),n1*n2, &
            MPI_REAL_SD,rightpe,tagrs, &
            f3(1,1,kparasta-1),n1*n2, &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)
        !
        return
    end subroutine flucrhoesp

    subroutine flud1(rc,cgra1,r,isc,tipo)
        !***********************************************************************
        ! compute explicit flux on csi component for scalar equation
        !
        ! convective uc*(rho) with centered scheme or quick
        ! diffusive akapt*g12*d(rho)/d(eta)
        ! diffusive akapt*g13*d(rho)/d(zita)
        !
        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer :: isc
        integer :: kparastak,kparaendk
        integer :: i,j,k,is,iss,jj,kk,ii
        !
        real :: ak,r0,r1,r2,fg
        real :: rc(0:n1,n2,kparasta:kparaend) !n3)
        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real :: cgra1(0:n1,n2,kparasta-1:kparaend+1)
        real :: ravanti,rindietro1,rindietro2
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !      integer :: tipo(0:n1+1,0:n2+1,0:n3+1) c7
        logical :: found_solid
        !-----------------------------------------------------------------------
        ! periodicity not in eta
        !
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Uc*rho:
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell

        !     sides left and right
        do i=1,ip

            do k=kparasta,kparaend
                do j=1,n2
                    !
                    f1(0,j,k)=rc(0,j,k)*r(0,j,k)-cgra1(0,j,k)
                    f1(n1,j,k)=rc(n1,j,k)*r(n1+1,j,k)-cgra1(n1,j,k)
                !
                end do
            end do

        end do
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !
                    f1(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
                !
                end do
            end do
        end do
        !
        ! quick
        !
        if (insc==1) then
            !     sides 1 and 2
            do k=kparasta,kparaend
                do j=1,n2

                    i=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(2,j,k)
                        rindietro1=r(1,j,k)
                        rindietro2=ip*(2.*r(0,j,k)-r(1,j,k))+(1.-ip)*r(0,j,k)

                        if (tipo(i,j,k)==1) then !ib cell
                            if (tipo(i+1,j,k)==0) then !solido davanti
                                f1(i,j,k)=0.
                            elseif (tipo(i-1,j,k)==0) then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-ravanti +2.*rindietro1-rindietro2)
                            end if !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo

                    else
                        ravanti=r(1,j,k)
                        rindietro1=r(2,j,k)
                        rindietro2=r(3,j,k)

                        if (tipo(i+1,j,k)==1) then
                            if (tipo(i,j,k)==0) then !solido indietro
                                f1(i,j,k)=0.
                            elseif (tipo(i+2,j,k)==0) then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else !solido di lato
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*( -ravanti +2.*rindietro1-rindietro2)
                            end if !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo
                    end if

                    !         if (tipo(i,j,k)==2) then
                    !    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    !         end if

                    i=n1-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(n1  ,j,k)
                        rindietro1=r(n1-1,j,k)
                        rindietro2=r(n1-2,j,k)

                        if (tipo(i,j,k)==1) then !ib cell
                            if (tipo(i+1,j,k)==0) then !solido davanti
                                f1(i,j,k)=0.
                            elseif (tipo(i-1,j,k)==0) then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-ravanti +2.*rindietro1-rindietro2)
                            end if !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo
                    else
                        ravanti=r(n1-1,j,k)
                        rindietro1=r(n1  ,j,k)
                        rindietro2=ip*(2.*r(n1+1,j,k)-r(n1,j,k))+(1.-ip)*r(n1+1,j,k)

                        if (tipo(i+1,j,k)==1) then
                            if (tipo(i,j,k)==0) then !solido indietro
                                f1(i,j,k)=0.
                            elseif (tipo(i+2,j,k)==0) then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else !solido di lato
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*( -ravanti +2.*rindietro1-rindietro2)
                            end if !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo
                    end if

                !         if (tipo(i,j,k)==2) then
                !   f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                !     >              *.125*( -ravanti +2.*rindietro1-rindietro2 )
                !         end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=1,n2
                    do i=2,n1-2
                        !Giulia non comprende tutti i casi, da un lato è quick dall'altro dc
                        !      if (tipo(i,j,k)==2) then
                        !     if (rc(i,j,k)>0.) then
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                        !      else
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                        !      end if
                        !         end if
                        !
                        if (rc(i,j,k)>0.) then

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i+1,j,k)==0) then !solido davanti
                                    f1(i,j,k)=0.
                                elseif (tipo(i-1,j,k)==0) then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                            end if !tipo

                        else !rc(i,j,k)<=0.

                            if (tipo(i+1,j,k)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f1(i,j,k)=0.
                                elseif (tipo(i+2,j,k)==0) then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else !solido di lato
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                            end if !tipo

                        end if !rc
                    end do
                end do
            end do


        ! QUICK MODIFICATO /SMART
        elseif (insc==2) then

            ! GIULIA UCS-CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,n2
                        !
                        f1(0,j,k)=rc(0,j,k)*r(0,j,k)-cgra1(0,j,k)
                        f1(n1,j,k)=rc(n1,j,k)*r(n1+1,j,k)-cgra1(n1,j,k)
                    !
                    end do
                end do

            end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2
                    do i=ip,n1-ip
                        !
                        f1(i,j,k)=-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,n2
                    i=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(2,j,k)
                        rindietro1=r(1,j,k)
                        rindietro2=ip*(2.*r(0,j,k)-r(1,j,k)) &
                            +(1.-ip)*r(0,j,k)
                    else
                        ravanti=r(1,j,k)
                        rindietro1=r(2,j,k)
                        rindietro2=r(3,j,k)
                    end if
                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)
                    !if (den<0.0000000001.and.den>=0) den=0.0000000001
                    !if (den>-0.0000000001.and.den<0) den=-0.000000001
                    !        if (den==0) then
                    !        phic(i,j,k)=1.
                    !   else
                    !   phic(i,j,k)=(num /den)
                    !        end if
                    !   !if (phic(i,j,k)>1.or.phic(i,j,k)<0)phic(i,j,k)=-0.1
                    !
                    !    if (phic(i,j,k)>=0.17.and.phic(i,j,k)<=0.83) then !1/6<phic<4/5
                    !       f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                    !     >           *.125*( 3*ravanti +6*rindietro1-rindietro2 )
                    !   else if (phic(i,j,k)>0.and.phic(i,j,k)<0.17) then !0<phic<1/6
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                    !   else if (phic(i,j,k)>0.83.and.phic(i,j,k)<1) then !4/5<phic<1
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                    !   else
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1
                    !   end if!phic
                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(3*ravanti+6*rindietro1-rindietro2)


                    i=n1-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(n1   ,j,k)
                        rindietro1=r(n1-1,j, k)
                        rindietro2=r(n1-2,j, k)
                    else
                        ravanti=r(n1-1,j,k)
                        rindietro1=r(n1  ,j,k)
                        rindietro2=ip*(2.*r(n1+1,j,k)-r(n1,j,k)) &
                            +(1.-ip)*r(n1+1,j,k)
                    end if

                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)
                    !   !if (den<0.0000000001.and.den>=0) den=0.0000000001
                    !   !if (den>-0.0000000001.and.den<0) den=-0.000000001
                    !        if (den==0) then
                    !        phic(i,j,k)=1.
                    !   else
                    !   phic(i,j,k)=(num /den)
                    !        end if
                    !   !if (phic(i,j,k)>1.or.phic(i,j,k)<0)phic(i,j,k)=-0.1

                    !   if (phic(i,j,k)>=0.17.and.phic(i,j,k)<=0.83) then
                    !       f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                    !     >           *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )
                    !   else if (phic(i,j,k)>0.and.phic(i,j,k)<0.17) then
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                    !   else if (phic(i,j,k)>0.83.and.phic(i,j,k)<1) then
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                    !   else
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1
                    !   end if!phic

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(3*ravanti+6*rindietro1-rindietro2)
                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                            else  !(rc(i,j,k)<=0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                                    .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                            end if
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i+1,j,k)==0) then !solido davanti

                                        f1(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    elseif (tipo(i-1,j,k)==0) then

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                            *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i-1,j,k)==1) then !la cella prima è un ib
                                    !       if (tipo(i-2,j,k)==0) then !perché considero anche questa?
                                    !              rindietro2=r(i,j,k)

                                    !       end if
                                    !         end if

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                end if !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i+1,j,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f1(i,j,k)=0.
                                    elseif (tipo(i+2,j,k)==0) then !solido avanti
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                    else !solido di lato
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                                            .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))

                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib



                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if (tipo(i+1,j,k)==1) then
                                    !       if (tipo(i+2,j,k)==0) then
                                    !             rindietro2=r(i+1,j,k)
                                    !       end if
                                    !       end if

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                                        .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        !    end if !SMART /quick modificato upwind

        elseif (insc==3) then !upwind

            ! GIULIA UCS-CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,n2
                        !
                        f1(0,j,k)=rc(0,j,k)*r(0,j,k)-cgra1(0,j,k)
                        f1(n1,j,k)=rc(n1,j,k)*r(n1+1,j,k)-cgra1(n1,j,k)
                    !
                    end do
                end do

            end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2
                    do i=ip,n1-ip
                        !
                        f1(i,j,k)=-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,n2
                    i=1
                    if (rc(i,j,k)>0.) then
                        rindietro1=r(1,j,k)
                    else
                        rindietro1=r(2,j,k)

                    end if
                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1



                    i=n1-1
                    if (rc(i,j,k)>0.) then
                        rindietro1=r(n1-1,j, k)
                    else
                        rindietro1=r(n1  ,j,k)
                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k)<=0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                            end if
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i+1,j,k)==0) then !solido davanti

                                        f1(i,j,k)=0.

                                    !       elseif (tipo(i-1,j,k)==0) then
                                    !
                                    !               f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione o dietro

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i-1,j,k)==1) then !la cella prima è un ib
                                    !       if (tipo(i-2,j,k)==0) then !perché considero anche questa?
                                    !              rindietro2=r(i,j,k)

                                    !       end if
                                    !         end if

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                                end if !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i+1,j,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f1(i,j,k)=0.
                                        !        elseif (tipo(i+2,j,k)==0) then !solido avanti
                                        !              f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                        !               else !solido di lato
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib



                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if (tipo(i+1,j,k)==1) then
                                    !       if (tipo(i+2,j,k)==0) then
                                    !             rindietro2=r(i+1,j,k)
                                    !       end if
                                    !       end if

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        end if !quick





        !      end if
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G12*Drho/D(ETA)
        !
        !     sides left and right
        is=0
        iss=0
        do i=1,2*ip
            do k=kparasta,kparaend
                !
                do j=2,n2-1
                    !
                    fg=.5*(r(iss,j+1,k)-r(iss,j-1,k))
                    !             f1(is,j,k)=-f1(is,j,k)
                    !     >                 +akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    !
                    if (tipo(iss,j-1,k)==0.or.tipo(iss,j+1,k)==0) then
                        f1(is,j,k)=-f1(is,j,k)
                    !         elseif (tipo(iss,j+1,k)==0) then
                    !           f1(is,j,k)=-f1(is,j,k)
                    else
                        f1(is,j,k)=-f1(is,j,k) &
                            +akapt(isc,iss,j,k)*g12(is,j,k)*fg

                    end if

                end do

                !
                ! direction 2 is always not periodic
                    !
                    !        check derivative at wall bottom and upper
                    !
                    j=1
                    fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j+1,k)-r(iss,j+2,k))
                    !         f1(is,j,k)=-f1(is,j,k)
                    !     >                +akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    !
                    if (tipo(iss,j+1,k)==0.or.tipo(iss,j+2,k)==0) then
                        f1(is,j,k)=-f1(is,j,k)
                    !         elseif (tipo(iss,j+2,k)==0) then
                    !           f1(is,j,k)=-f1(is,j,k)
                    else
                        f1(is,j,k)=-f1(is,j,k) &
                            +akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    end if


                    j=n2
                    fg=.5*(3.*r(iss,n2,k)-4.*r(iss,n2-1,k)+r(iss,n2-2,k))
                    !         f1(is,j,k)=-f1(is,j,k)
                    !     >                +akapt(isc,iss,j,k)*g12(is,j,k)*fg

                    if (tipo(iss,j-1,k)==0.or.tipo(iss,j-2,k)==0) then
                        f1(is,j,k)=-f1(is,j,k)
                    !         elseif (tipo(iss,j-2,k)==0) then
                    !           f1(is,j,k)=-f1(is,j,k)
                    else
                        f1(is,j,k)=-f1(is,j,k) &
                            +akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    end if
                !
            !
            end do
            !
            is=n1
            iss=n1+1
        !
        end do
        !
        ! inside the field
        !
        do k=kparasta,kparaend
            do i=ip,n1-ip
                !
                do j=2,n2-1
                    !
                    r1=.5*(r(i,j-1,k)+r(i+1,j-1,k))
                    r2=.5*(r(i,j+1,k)+r(i+1,j+1,k))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))


                    found_solid=.false.
                    do kk=k,k
                        do jj=j-1,j+1
                            do ii=i,i+1
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f1(i,j,k)=-f1(i,j,k)
                    else
                        f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    end if

                !
                end do
                !
                !     check derivative on sides bottom and upper
                !
                ! direction 2 is always not periodic
                    !
                    j=1
                    r0=.5*(r(i,j  ,k)+r(i+1,j  ,k))
                    r1=.5*(r(i,j+1,k)+r(i+1,j+1,k))
                    r2=.5*(r(i,j+2,k)+r(i+1,j+2,k))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
                    !      f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg

                    found_solid=.false.
                    do kk=k,k
                        do jj=j,j+2
                            do ii=i,i+1
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f1(i,j,k)=-f1(i,j,k)
                    else
                        f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    end if

                    !
                    j=n2
                    r0=.5*(r(i,n2  ,k)+r(i+1,n2  ,k))
                    r1=.5*(r(i,n2-1,k)+r(i+1,n2-1,k))
                    r2=.5*(r(i,n2-2,k)+r(i+1,n2-2,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
                    !      f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    !
                    found_solid=.false.
                    do kk=k,k
                        do jj=j-2,j
                            do ii=i,i+1
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f1(i,j,k)=-f1(i,j,k)
                    else
                        f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    end if
            !
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G13*Drho/D(ZITA)
        !
        ! sides left and right
        !
        ! define the limit dependeing on periodicity in z
        !
        if (myid==0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid==nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else if ((myid/=0).and.(myid/=nproc-1)) then
            kparastak=kparasta
            kparaendk=kparaend
        end if
        !
        is=0
        iss=0
        do i=1,2*ip
            !
            do j=1,n2
                !
                do k=kparastak,kparaendk
                    !
                    fg=.5*(r(iss,j,k+1)-r(iss,j,k-1))
                    !       f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
                    !

                    if (tipo(iss,j,k-1)==0.or.tipo(iss,j,k+1)==0) then
                        f1(is,j,k)=f1(is,j,k)
                    !         elseif (tipo(iss,j,k+1)==0) then
                    !           f1(is,j,k)=f1(is,j,k)
                    else
                        f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
                    end if
                end do
                !
                do kk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j,k+1)-r(iss,j,k+2))
                        !      f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg

                        if (tipo(iss,j,k+1)==0.or.tipo(iss,j,k+2)==0) then
                            f1(is,j,k)=f1(is,j,k)
                        !         elseif (tipo(iss,j,k+2)==0) then
                        !           f1(is,j,k)=f1(is,j,k)
                        else
                            f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
                        end if

                    end if
                    !
                    if (myid==nproc-1) then
                        k=n3
                        fg=.5*(3.*r(iss,j,n3)-4.*r(iss,j,n3-1)+r(iss,j,n3-2))
                        !      f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg

                        if (tipo(iss,j,k-1)==0.or.tipo(iss,j,k-2)==0) then
                            f1(is,j,k)=f1(is,j,k)
                        !         elseif (tipo(iss,j,k-2)==0) then
                        !           f1(is,j,k)=f1(is,j,k)
                        else
                            f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
                        end if
                    end if
                !
                end do
            !
            end do
            !
            is=n1
            iss=n1+1
        !
        end do
        !
        ! into the field
        !

        do j=1,n2
            do i=ip,n1-ip
                !
                do k=kparastak,kparaendk
                    !
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
                    r1=.5*(r(i,j,k-1)+r(i+1,j,k-1))
                    r2=.5*(r(i,j,k+1)+r(i+1,j,k+1))
                    fg=.5*(r2-r1)

                    !
                    found_solid=.false.
                    do kk=k-1,k+1
                        do jj=j,j
                            do ii=i,i+1
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f1(i,j,k)=f1(i,j,k)
                    else
                        f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
                    end if




                end do
                !
                if (kp==1) then
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        r0=.5*(r(i,j,k  )+r(i+1,j,k  ))
                        r1=.5*(r(i,j,k+1)+r(i+1,j,k+1))
                        r2=.5*(r(i,j,k+2)+r(i+1,j,k+2))
                        fg=.5*(-3.*r0+4.*r1-r2)
                        ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
                        !      f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg


                        found_solid=.false.
                        do kk=k,k+2
                            do jj=j,j
                                do ii=i,i+1
                                    if (tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if (found_solid) then
                            f1(i,j,k)=f1(i,j,k)
                        else
                            f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
                        end if



                    end if
                    !
                    if (myid==nproc-1) then

                        k=n3
                        r0=.5*(r(i,j,n3  )+r(i+1,j,n3  ))
                        r1=.5*(r(i,j,n3-1)+r(i+1,j,n3-1))
                        r2=.5*(r(i,j,n3-2)+r(i+1,j,n3-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
                        !      f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg

                        found_solid=.false.
                        do kk=k-2,k
                            do jj=j,j
                                do ii=i,i+1
                                    if (tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if (found_solid) then
                            f1(i,j,k)=f1(i,j,k)
                        else
                            f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
                        end if

                    end if
                !
                end if
            !
            end do
        end do


        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !

                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f1(i,j,k)=0.
                        if (i<n1) then
                            if (tipo(i+1,j,k)==0)f1(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do
        !
        return
    end subroutine flud1

    subroutine flud2(rc,cgra2,r,isc,tipo)
        !***********************************************************************
        ! compute explicit flux on eta component for scalar equation
        !
        ! convective vc*(rho) with centered scheme or quick
        ! diffusive akapt*g21*d(rho)/d(csi)
        ! diffusive akapt*g23*d(rho)/d(zita)
        !
        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer :: isc
        integer :: kparastak,kparaendk
        integer :: i,j,k,js,jss,ii,kk,jj

        real :: ak,r0,r1,r2,fg
        real :: rc(n1,0:n2,kparasta:kparaend) !n3)
        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real :: cgra2(n1,0:n2,kparasta-1:kparaend+1)
        real :: ravanti,rindietro1,rindietro2
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        logical :: found_solid
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Vc*rho:
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell

        !     sides bottom and upper
        ! direction 2 is always not periodic
            !
            do k=kparasta,kparaend
                do i=1,n1
                    !
                    f2(i,0,k)=rc(i,0,k)*r(i,0,k)-cgra2(i,0,k)
                    f2(i,n2,k)=rc(i,n2,k)*r(i,n2+1,k)-cgra2(i,n2,k)
                !
                end do
            end do

        !
        ! into the field
        !
        do k=kparasta,kparaend
            do i=1,n1
                do j=1,n2-1
                    !
                    f2(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j+1,k))-cgra2(i,j,k)
                !
                end do
            end do
        end do
        !
        ! quick
        !
        if (insc==1) then

            !     sides 3 and 4
            do k=kparasta,kparaend
                do i=1,n1

                    j=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,2,k)
                        rindietro1=r(i,1,k)
                        rindietro2=2.*r(i,0,k)-r(i,1,k)

                        if (tipo(i,j,k)==1) then !ib cell
                            if (tipo(i,j+1,k)==0) then !solido davanti
                                f2(i,j,k)=0.
                            elseif (tipo(i,j-1,k)==0) then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2)

                            end if !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo
                    else
                        ravanti=r(i,1,k)
                        rindietro1=r(i,2,k)
                        rindietro2=r(i,3,k)

                        if (tipo(i,j+1,k)==1) then
                            if (tipo(i,j,k)==0) then !solido indietro
                                f2(i,j,k)=0.
                            elseif (tipo(i,j+2,k)==0) then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)

                            else !solido di lato
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2)


                            end if !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo
                    end if

                    !         if (tipo(i,j,k)==2) then
                    !    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    !         end if

                    j=n2-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,n2  ,k)
                        rindietro1=r(i,n2-1,k)
                        rindietro2=r(i,n2-2,k)

                        if (tipo(i,j,k)==1) then !ib cell
                            if (tipo(i,j+1,k)==0) then !solido davanti
                                f2(i,j,k)=0.
                            elseif (tipo(i,j-1,k)==0) then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2)

                            end if !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo
                    else
                        ravanti=r(i,n2-1,k)
                        rindietro1=r(i,n2  ,k)
                        rindietro2=2.*r(i,n2+1,k)-r(i,n2,k)

                        if (tipo(i,j+1,k)==1) then
                            if (tipo(i,j,k)==0) then !solido indietro
                                f2(i,j,k)=0.
                            elseif (tipo(i,j+2,k)==0) then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)

                            else !solido di lato
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2)


                            end if !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1-rindietro2)
                        end if !tipo
                    end if
                !         if (tipo(i,j,k)==2) then
                !   f2(i,j,k)=f2(i,j,k)+rc(i,j,k)
                !     >              *.125*( -ravanti +2.*rindietro1-rindietro2 )
                !         end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=2,n2-2
                    do i=1,n1
                        !
                        !         if (tipo(i,j,k)==2) then
                        !      if (rc(i,j,k)>0.) then
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                        !      else
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                        !      end if
                        !      end if
                        !

                        if (rc(i,j,k)>0.) then

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i,j+1,k)==0) then !solido davanti
                                    f2(i,j,k)=0.
                                elseif (tipo(i,j-1,k)==0) then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                            end if !tipo

                        else !rc(i,j,k)<=0.
                            if (tipo(i,j+1,k)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f2(i,j,k)=0.
                                elseif (tipo(i,j+2,k)==0) then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)

                                else !solido di lato
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))


                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                            end if !tipo

                        end if !rc
                    !
                    end do
                end do
            end do

        ! smart/quick modificato
        elseif (insc==2) then

            ! GIULIA UCS-CGRA
            !     side left and right
            ! direction 2 is always not periodic
                !
                do k=kparasta,kparaend
                    do i=1,n1
                        !
                        f2(i,0,k)=rc(i,0,k)*r(i,0,k)-cgra2(i,0,k)
                        f2(i,n2,k)=rc(i,n2,k)*r(i,n2+1,k)-cgra2(i,n2,k)
                    !
                    end do
                end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2-1
                    do i=1,n1
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,n1
                    j=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,2,k)
                        rindietro1=r(i,1,k)
                        rindietro2=2.*r(i,0,k)-r(i,1,k)
                    else
                        ravanti=r(i,1,k)
                        rindietro1=r(i,2,k)
                        rindietro2=r(i,3,k)
                    end if
                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)
                    !   !if (den<0.0000000001.and.den>=0) den=0.0000000001
                    !   !if (den>-0.0000000001.and.den<0) den=-0.000000001
                    !        if (den==0) then
                    !        phic(i,j,k)=1.
                    !   else
                    !   phic(i,j,k)=(num /den)
                    !        end if
                    !   !if (phic(i,j,k)>1.or.phic(i,j,k)<0)phic(i,j,k)=-0.1

                    !    if (phic(i,j,k)>=0.17.and.phic(i,j,k)<=0.83) then !1/6<phic<4/5
                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )
                    !   else if (phic(i,j,k)>0.and.phic(i,j,k)<0.17) then !0<phic<1/6
                    !        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                    !   else if (phic(i,j,k)>0.83.and.phic(i,j,k)<1) then !4/5<phic<1
                    !        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                    !   else
                    !        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1
                    !   end if!phic

                    j=n2-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,n2  , k)
                        rindietro1=r(i,n2-1, k)
                        rindietro2=r(i,n2-2, k)
                    else
                        ravanti=r(i,n2-1,k)
                        rindietro1=r(i,n2  ,k)
                        rindietro2=2.*r(i,n2+1,k)-r(i,n2,k)
                    end if

                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)

                    !        !if (den<0.0000000001.and.den>=0) den=0.0000000001
                    !   !if (den>-0.0000000001.and.den<0) den=-0.000000001
                    !        if (den==0) then
                    !        phic(i,j,k)=1.
                    !   else
                    !   phic(i,j,k)=(num /den)
                    !        end if
                    !   !if (phic(i,j,k)>1.or.phic(i,j,k)<0)phic(i,j,k)=-0.1

                    !   if (phic(i,j,k)>=0.17.and.phic(i,j,k)<=0.83) then
                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )
                !  else if (phic(i,j,k)>0.and.phic(i,j,k)<0.17) then
                !       f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                !  else if (phic(i,j,k)>0.83.and.phic(i,j,k)<1) then
                !       f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                !  else
                !       f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1
                !  end if!phic
                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                            else  !(rc(i,j,k)<=0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k)) !my index +1

                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j+1,k)==0) then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo
                                    elseif (tipo(i,j-1,k)==0) then

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione


                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                                end if   !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j+1,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f2(i,j,k)=0.
                                    elseif (tipo(i,j+2,k)==0) then !solido avanti
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                    else !solido di lato

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))
                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        !    end if !SMART/QUICK modificato upwind

        ! smart/quick modificato
        elseif (insc==3) then

            ! GIULIA UCS-CGRA
            !     side left and right
            ! direction 2 is always not periodic
                !
                do k=kparasta,kparaend
                    do i=1,n1
                        !
                        f2(i,0,k)=rc(i,0,k)*r(i,0,k)-cgra2(i,0,k)
                        f2(i,n2,k)=rc(i,n2,k)*r(i,n2+1,k)-cgra2(i,n2,k)
                    !
                    end do
                end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2-1
                    do i=1,n1
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,n1
                    j=1
                    if (rc(i,j,k)>0.) then

                        rindietro1=r(i,1,k)

                    else

                        rindietro1=r(i,2,k)

                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1


                    j=n2-1
                    if (rc(i,j,k)>0.) then

                        rindietro1=r(i,n2-1, k)

                    else

                        rindietro1=r(i,n2  ,k)

                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)


                            else  !(rc(i,j,k)<=0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j+1,k)==0) then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo
                                    !       elseif (tipo(i,j-1,k)==0) then

                                    !               f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione


                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                end if   !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j+1,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f2(i,j,k)=0.
                                    else
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        end if !quick
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G21*Drho/D(CSI)
        !
        !     sides bottom and upper
        js=0
        jss=0
        do j=1,2
            !
            do k=kparasta,kparaend
                !
                do  i=1+ip,n1-ip
                    !
                    fg=.5*(r(i+1,jss,k)-r(i-1,jss,k))
                    !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
                    !

                    if (tipo(i-1,jss,k)==0.or.tipo(i+1,jss,k)==0) then
                        f2(i,js,k)=-f2(i,js,k)
                    !         elseif (tipo(i+1,jss,k)==0) then
                    !           f2(i,js,k)=-f2(i,js,k)
                    else
                        f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
                    end if

                end do
                !
                do ii=1,ip
                    !
                    !     check derivative on sides left and right
                    !
                    i=1
                    fg=.5*(-3.*r(i,jss,k)+4.*r(i+1,jss,k)-r(i+2,jss,k))
                    !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg

                    if (tipo(i+1,jss,k)==0.or.tipo(i+2,jss,k)==0) then
                        f2(i,js,k)=-f2(i,js,k)
                    !         elseif (tipo(i+2,jss,k)==0) then
                    !           f2(i,js,k)=-f2(i,js,k)
                    else
                        f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
                    end if
                    !
                    i=n1
                    fg=.5*(3.*r(n1,jss,k)-4.*r(n1-1,jss,k)+r(n1-2,jss,k))
                    !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg

                    if (tipo(i-1,jss,k)==0.or.tipo(i-2,jss,k)==0) then
                        f2(i,js,k)=-f2(i,js,k)
                    !         elseif (tipo(i-2,jss,k)==0) then
                    !           f2(i,js,k)=-f2(i,js,k)
                    else
                        f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
                    end if
                !
                end do
            !
            end do
            !
            js=n2
            jss=n2+1
        !
        end do
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=1,n2-1
                !
                do i=1+ip,n1-ip
                    !
                    r1=.5*(r(i-1,j,k)+r(i-1,j+1,k))
                    r2=.5*(r(i+1,j,k)+r(i+1,j+1,k))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))


                    found_solid=.false.
                    do kk=k,k
                        do jj=j,j+1
                            do ii=i-1,i+1
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f2(i,j,k)=-f2(i,j,k)
                    else
                        f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    end if


                end do
                !
                if (ip==1) then
                    !
                    ! check derivative on sides left and right
                    !
                    i=1
                    r0=.5*(r(i  ,j,k)+r(i  ,j+1,k))
                    r1=.5*(r(i+1,j,k)+r(i+1,j+1,k))
                    r2=.5*(r(i+2,j,k)+r(i+2,j+1,k))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
                    !      f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg

                    found_solid=.false.
                    do kk=k,k
                        do jj=j,j+1
                            do ii=i,i+2
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f2(i,j,k)=-f2(i,j,k)
                    else
                        f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    end if
                    !
                    i=n1
                    r0=.5*(r(n1  ,j,k)+r(n1  ,j+1,k))
                    r1=.5*(r(n1-1,j,k)+r(n1-1,j+1,k))
                    r2=.5*(r(n1-2,j,k)+r(n1-2,j+1,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
                    !      f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    !
                    found_solid=.false.
                    do kk=k,k
                        do jj=j,j+1
                            do ii=i-2,i
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f2(i,j,k)=-f2(i,j,k)
                    else
                        f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    end if

                end if
            !
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G23*Drho/D(ZITA)
        !
        ! sides bottom and upper
        !
        if (myid==0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid==nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else
            kparastak=kparasta
            kparaendk=kparaend
        end if

        js=0
        jss=0
        do j=1,2
            !
            do i=1,n1
                !
                do k=kparastak,kparaendk
                    !
                    fg=.5*(r(i,jss,k+1)-r(i,jss,k-1))
                    !         f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

                    if (tipo(i,jss,k-1)==0.or.tipo(i+1,jss,k+1)==0) then
                        f2(i,js,k)=f2(i,js,k)
                    !         elseif (tipo(i+1,jss,k+1)==0) then
                    !           f2(i,js,k)=f2(i,js,k)
                    else
                        f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
                    end if
                !
                end do
                !
                do kk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        fg=.5*(-3.*r(i,jss,k)+4.*r(i,jss,k+1)-r(i,jss,k+2))
                        !      f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

                        if (tipo(i,jss,k+1)==0.or.tipo(i+1,jss,k+2)==0) then
                            f2(i,js,k)=f2(i,js,k)
                        !         elseif (tipo(i+1,jss,k+2)==0) then
                        !           f2(i,js,k)=f2(i,js,k)
                        else
                            f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
                        end if

                    end if
                    !
                    if (myid==nproc-1) then

                        k=n3
                        fg=.5*(3.*r(i,jss,n3)-4.*r(i,jss,n3-1)+r(i,jss,n3-2))
                        !      f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

                        if (tipo(i,jss,k-1)==0.or.tipo(i+1,jss,k-2)==0) then
                            f2(i,js,k)=f2(i,js,k)
                        !         elseif (tipo(i+1,jss,k-2)==0) then
                        !           f2(i,js,k)=f2(i,js,k)
                        else
                            f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
                        end if

                    end if
                !
                end do
            !
            end do
            !
            js=n2
            jss=n2+1
        !
        end do
        !
        ! into the field
        !
        do j=1,n2-1
            do i=1,n1
                !
                do k=kparastak,kparaendk
                    !
                    r1=.5*(r(i,j,k-1)+r(i,j+1,k-1))
                    r2=.5*(r(i,j,k+1)+r(i,j+1,k+1))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
                    !
                    found_solid=.false.
                    do kk=k-1,k+1
                        do jj=j,j+1
                            do ii=i,i
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f2(i,j,k)=f2(i,j,k)
                    else
                        f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
                    end if



                end do
                !
                if (kp==1) then
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        r0=.5*(r(i,j,k  )+r(i,j+1,k  ))
                        r1=.5*(r(i,j,k+1)+r(i,j+1,k+1))
                        r2=.5*(r(i,j,k+2)+r(i,j+1,k+2))
                        fg=.5*(-3.*r0+4.*r1-r2)
                        ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
                        !      f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg

                        found_solid=.false.
                        do kk=k,k+2
                            do jj=j,j+1
                                do ii=i,i
                                    if (tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if (found_solid) then
                            f2(i,j,k)=f2(i,j,k)
                        else
                            f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
                        end if


                    end if
                    !
                    if (myid==nproc-1) then

                        k=n3
                        r0=.5*(r(i,j,n3  )+r(i,j+1,n3  ))
                        r1=.5*(r(i,j,n3-1)+r(i,j+1,n3-1))
                        r2=.5*(r(i,j,n3-2)+r(i,j+1,n3-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
                        !      f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg

                        found_solid=.false.
                        do kk=k-2,k
                            do jj=j,j+1
                                do ii=i,i
                                    if (tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if (found_solid) then
                            f2(i,j,k)=f2(i,j,k)
                        else
                            f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
                        end if
                    end if
                !
                end if
            !
            end do
        end do
        !


        do k=kparasta,kparaend
            do j=1,n2-1
                do i=1,n1
                    !

                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f2(i,j,k)=0.
                        if (j<n2) then
                            if (tipo(i,j+1,k)==0)f2(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do

        return
    end subroutine flud2

    subroutine flud3(rc,cgra3,r,isc,tipo)
        !***********************************************************************
        ! compute explicit flux on zita component for scalar equation
        !
        ! convective wc*(rho) with centered scheme or quick
        ! diffusive akapt*g31*d(rho)/d(csi)
        ! diffusive akapt*g32*d(rho)/d(eta)       zita)
        !
        !hicco riguardare la dicitura dei commenti per i termini diffusivi, non mi quadra!!!

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer :: ierr
        integer :: isc
        integer :: kparastal,kparaendl
        integer :: kparastall,kparaendll
        integer :: i,j,k,ks,kss,ii,jj,kk
        !
        real :: ak,r0,r1,r2,fg
        real :: rc(n1,n2,kparasta-1:kparaend) !0:n3)
        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real :: cgra3(n1,n2,kparasta-1:kparaend+1)
        real :: ravanti,rindietro1,rindietro2
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        integer :: status(MPI_STATUS_SIZE)
        logical :: found_solid
        !
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Wc*rho:
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell
        !
        ! sides back and front
        do k=1,kp

            do j=1,n2
                do i=1,n1
                    !
                    if (myid==0) then
                        f3(i,j,0)=rc(i,j,0 )*r(i,j,0   )-cgra3(i,j,0 )
                    end if

                    if (myid==nproc-1) then
                        f3(i,j,n3)=rc(i,j,n3)*r(i,j,n3+1)-cgra3(i,j,n3)
                    end if
                !
                end do
            end do

        end do
        !
        ! into the field
        !
        if (myid==0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid==nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else
            kparastal=kparasta
            kparaendl=kparaend
        end if

        do k=kparastal,kparaendl
            do j=1,n2
                do i=1,n1
                    !
                    f3(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
                !
                end do
            end do
        end do

        ! quick
        !
        if (insc==1) then

            if (myid==0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid==nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            end if

            !     sides 5 and 6
            if (myid==0) then
                do j=1,n2
                    do i=1,n1

                        k=1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,2)
                            rindietro1=r(i,j,1)
                            rindietro2=kp*(2.*r(i,j,0)-r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i,j,k+1)==0) then !solido davanti
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k-1)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1-rindietro2 )

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2 )
                            end if !tipo
                        else
                            ravanti=r(i,j,1)
                            rindietro1=r(i,j,2)
                            rindietro2=r(i,j,3)

                            if (tipo(i,j,k+1)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k+2)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1-rindietro2)
                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2)
                            end if !tipo
                        end if

                    !         if (tipo(i,j,k)==2) then
                    !    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    !         end if
                    end do
                end do
            end if

            if (myid==nproc-1) then
                do j=1,n2
                    do i=1,n1

                        k=n3-1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,n3  )
                            rindietro1=r(i,j,n3-1)
                            rindietro2=r(i,j,n3-2)

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i,j,k+1)==0) then !solido davanti
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k-1)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1-rindietro2 )

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2 )
                            end if !tipo
                        else
                            ravanti=r(i,j,n3-1)
                            rindietro1=r(i,j,n3  )
                            rindietro2=kp*(2.*r(i,j,n3+1)-r(i,j,n3)) &
                                +(1.-kp)*r(i,j,n3+1)


                            if (tipo(i,j,k+1)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k+2)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1-rindietro2)
                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1-rindietro2)
                            end if !tipo

                        end if

                    !         if (tipo(i,j,k)==2) then
                    !    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    !         end if
                    end do
                end do
            end if

            !     into the field

            do k=kparastall,kparaendll
                do j=1,n2
                    do i=1,n1
                        !
                        !      if (tipo(i,j,k)==2) then
                        !      if (rc(i,j,k)>0.) then
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                        !      else
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                        !      end if
                        !      end if

                        if (rc(i,j,k)>0.) then

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i,j,k+1)==0) then !solido davanti
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k-1)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                            end if !tipo

                        else !rc(i,j,k)<=0.

                            if (tipo(i,j,k+1)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k+2)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                            end if !tipo

                        end if !rc
                    !
                    end do
                end do
            end do


        ! quick modificato
        elseif (insc==2) then

            ! GIULIA UCS-CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,n2
                    do i=1,n1
                        !
                        if (myid==0) then
                            f3(i,j,0)=rc(i,j,0 )*r(i,j,0   )-cgra3(i,j,0 )
                        end if

                        if (myid==nproc-1) then
                            f3(i,j,n3)=rc(i,j,n3)*r(i,j,n3+1)-cgra3(i,j,n3)
                        end if
                    !
                    end do
                end do
            !
            end do

            if (myid==0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid==nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            end if


            do k=kparastal,kparaendl
                do j=1,n2
                    do i=1,n1
                        !
                        f3(i,j,k)=-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid==0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid==nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            end if

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if (myid==0) then
                do j=1,n2
                    do i=1,n1
                        k=1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,2)
                            rindietro1=r(i,j,1)
                            rindietro2=kp*(2.*r(i,j,0)-r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti=r(i,j,1)
                            rindietro1=r(i,j,2)
                            rindietro2=r(i,j,3)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*(3*ravanti+6*rindietro1-rindietro2)

                    end do
                end do

            end if !myid==0

            if (myid==nproc-1) then
                do j=1,n2
                    do i=1,n1

                        k=n3-1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,n3  )
                            rindietro1=r(i,j,n3-1)
                            rindietro2=r(i,j,n3-2)
                        else
                            ravanti=r(i,j,n3-1)
                            rindietro1=r(i,j,n3  )
                            rindietro2=kp*(2.*r(i,j,n3+1)-r(i,j,n3)) &
                                +(1.-kp)*r(i,j,n3+1)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*(3*ravanti+6*rindietro1-rindietro2)
                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                            else  !(rc(i,j,k)<=0.)
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(3*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j,k+1)==0) then !solido davanti

                                        f3(i,j,k)=0.

                                    elseif (tipo(i,j,k-1)==0) then !solido dietro

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else !solido lati
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                            *.125*(3*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))
                                    end if !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k-1)==1) then !la cella prima è un ib
                                    !       if (tipo(i,j,k-2)==0) then !perché considero anche questa?
                                    !         rindietro2=r(i,j,k)
                                    !       end if
                                    !         end if !fine k-1 IB

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                        *.125*(3*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))


                                end if !tip02

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j,k+1)==1) then !ib
                                    if (tipo(i,j,k)==0) then !solido avanti
                                        f3(i,j,k)=0.
                                    elseif (tipo(i,j,k+2)==0) then !solido indietro
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)
                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                            .125*(3*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

                                    end if !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k+1)==1) then
                                    !       if (tipo(i,j,k+2)==0) then
                                    !       rindietro2=r(i,j,k+1)
                                    !       end if
                                    !         end if


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(3*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))
                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        !    end if !SMART
        !-----SMART

        ! quick modificato
        elseif (insc==3) then

            ! GIULIA UCS-CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,n2
                    do i=1,n1
                        !
                        if (myid==0) then
                            f3(i,j,0)=rc(i,j,0 )*r(i,j,0   )-cgra3(i,j,0 )
                        end if

                        if (myid==nproc-1) then
                            f3(i,j,n3)=rc(i,j,n3)*r(i,j,n3+1)-cgra3(i,j,n3)
                        end if
                    !
                    end do
                end do
            !
            end do

            if (myid==0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid==nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            end if


            do k=kparastal,kparaendl
                do j=1,n2
                    do i=1,n1
                        !
                        f3(i,j,k)=-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid==0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid==nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            end if

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if (myid==0) then
                do j=1,n2
                    do i=1,n1
                        k=1
                        if (rc(i,j,k)>0.) then

                            rindietro1=r(i,j,1)

                        else

                            rindietro1=r(i,j,2)

                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1

                    end do
                end do

            end if !myid==0

            if (myid==nproc-1) then
                do j=1,n2
                    do i=1,n1

                        k=n3-1
                        if (rc(i,j,k)>0.) then

                            rindietro1=r(i,j,n3-1)

                        else

                            rindietro1=r(i,j,n3  )

                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1
                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k)<=0.)
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j,k+1)==0) then !solido davanti

                                        f3(i,j,k)=0.


                                    else !solido lati
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    end if !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k-1)==1) then !la cella prima è un ib
                                    !       if (tipo(i,j,k-2)==0) then !perché considero anche questa?
                                    !         rindietro2=r(i,j,k)
                                    !       end if
                                    !         end if !fine k-1 IB

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)


                                end if !tip02

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j,k+1)==1) then !ib
                                    if (tipo(i,j,k)==0) then !solido avanti
                                        f3(i,j,k)=0.

                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

                                    end if !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k+1)==1) then
                                    !       if (tipo(i,j,k+2)==0) then
                                    !       rindietro2=r(i,j,k+1)
                                    !       end if
                                    !         end if


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)
                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        end if !SMART
        !-----SMART


        !----------------------------------------------------------------------
        ! DIFFUSIVE TERM  NNI*G31*Drho/D(CSI)
        !
        !     sides back and front
        ks=0
        kss=0
        do k=1,2*kp

            if (myid==(k-1)*(nproc-1)) then

                do j=1,n2
                    !
                    do i=1+ip,n1-ip
                        !
                        fg=.5*(r(i+1,j,kss)-r(i-1,j,kss))
                        !         f3(i,j,ks)=-f3(i,j,ks)
                        !     >                 +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        !

                        if (tipo(i-1,j,kss)==0) then
                            f3(i,j,ks)=-f3(i,j,ks)
                        elseif (tipo(i+1,j,kss)==0) then
                            f3(i,j,ks)=-f3(i,j,ks)
                        else
                            f3(i,j,ks)=-f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        end if
                    end do
                    !
                    do ii=1,ip
                        !
                        !     check derivative on sides left and right
                        i=1
                        fg=.5*(-3.*r(i,j,kss)+4.*r(i+1,j,kss)-r(i+2,j,kss))
                        !      f3(i,j,ks)=-f3(i,j,ks)
                        !      >              +akapt(isc,i,j,kss)*g31(i,j,ks)*fg

                        if (tipo(i+1,j,kss)==0.or.tipo(i+2,j,kss)==0) then
                            f3(i,j,ks)=-f3(i,j,ks)
                        !         elseif (tipo(i+2,j,kss)==0) then
                        !           f3(i,j,ks)=-f3(i,j,ks)
                        else
                            f3(i,j,ks)=-f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        end if
                        !
                        i=n1
                        fg=.5*(3.*r(n1,j,kss)-4.*r(n1-1,j,kss)+r(n1-2,j,kss))
                        !      f3(i,j,ks)=-f3(i,j,ks)
                        !     >               +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        !
                        if (tipo(i-1,j,kss)==0.or.tipo(i-2,j,kss)==0) then
                            f3(i,j,ks)=-f3(i,j,ks)
                        !         elseif (tipo(i-2,j,kss)==0) then
                        !           f3(i,j,ks)=-f3(i,j,ks)
                        else
                            f3(i,j,ks)=-f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        end if
                    end do
                !
                end do
            !
            end if

            ks=n3
            kss=n3+1
        !
        end do
        !
        ! into the field
        !
        do k=kparastal,kparaendl
            do j=1,n2
                !
                do i=1+ip,n1-ip
                    !
                    r1=.5*(r(i-1,j,k)+r(i-1,j,k+1))
                    r2=.5*(r(i+1,j,k)+r(i+1,j,k+1))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))


                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j,j
                            do ii=i-1,i+1
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f3(i,j,k)=-f3(i,j,k)
                    else
                        f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg
                    end if



                !
                end do
                !
                if (ip==1) then
                    !
                    ! check derivative on sides left and right
                    !
                    i=1
                    r0=.5*(r(i  ,j,k)+r(i  ,j,k+1))
                    r1=.5*(r(i+1,j,k)+r(i+1,j,k+1))
                    r2=.5*(r(i+2,j,k)+r(i+2,j,k+1))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))
                    !      f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg

                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j,j
                            do ii=i,i+2
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f3(i,j,k)=-f3(i,j,k)
                    else
                        f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg
                    end if
                    !
                    i=n1
                    r0=.5*(r(n1  ,j,k)+r(n1  ,j,k+1))
                    r1=.5*(r(n1-1,j,k)+r(n1-1,j,k+1))
                    r2=.5*(r(n1-2,j,k)+r(n1-2,j,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))
                    !      f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg

                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j,j
                            do ii=i-2,i
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f3(i,j,k)=-f3(i,j,k)
                    else
                        f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg
                    end if
                !
                end if
            !
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G32*Drho/D(eta)   !ZITA)
        !
        !     sides back and front
        ks=0
        kss=0
        do k=1,2*kp

            if (myid==(k-1)*(nproc-1)) then

                do i=1,n1
                    !
                    do j=2,n2-1
                        !
                        fg=.5*(r(i,j+1,kss)-r(i,j-1,kss))
                        !         f3(i,j,ks)=f3(i,j,ks)
                        !     >                 +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        !
                        if (tipo(i,j-1,kss)==0.or.tipo(i,j+1,kss)==0) then
                            f3(i,j,ks)=f3(i,j,ks)
                        !         elseif (tipo(i,j+1,kss)==0) then
                        !           f3(i,j,ks)=f3(i,j,ks)
                        else
                            f3(i,j,ks)=f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        end if
                    end do
                    !
                    ! direction 2 is always not periodic
                        !
                        ! check derivative on sides back and front
                        !
                        j=1
                        fg=.5*(-3.*r(i,j,kss)+4.*r(i,j+1,kss)-r(i,j+2,kss))
                        !      f3(i,j,ks)=f3(i,j,ks)
                        !     >              +akapt(isc,i,j,kss)*g32(i,j,ks)*fg

                        if (tipo(i,j+1,kss)==0.or.tipo(i,j+2,kss)==0) then
                            f3(i,j,ks)=f3(i,j,ks)
                        !         elseif (tipo(i,j+2,kss)==0) then
                        !           f3(i,j,ks)=f3(i,j,ks)
                        else
                            f3(i,j,ks)=f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        end if
                        !
                        j=n2
                        fg=.5*(3.*r(i,n2,kss)-4.*r(i,n2-1,kss)+r(i,n2-2,kss))
                        !      f3(i,j,ks)=f3(i,j,ks)
                        !     >              +akapt(isc,i,j,kss)*g32(i,j,ks)*fg

                        if (tipo(i,j-1,kss)==0.or.tipo(i,j-2,kss)==0) then
                            f3(i,j,ks)=f3(i,j,ks)
                        !         elseif (tipo(i,j-2,kss)==0) then
                        !           f3(i,j,ks)=f3(i,j,ks)
                        else
                            f3(i,j,ks)=f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        end if
                    !
                end do
            !
            end if

            ks=n3
            kss=n3+1
        !
        end do
        !
        ! into the field
        !

        do k=kparastal,kparaendl
            do i=1,n1
                !
                do j=2,n2-1
                    !
                    r1=.5*(r(i,j-1,k)+r(i,j-1,k+1))
                    r2=.5*(r(i,j+1,k)+r(i,j+1,k+1))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))


                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j-1,j+1
                            do ii=i,i
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f3(i,j,k)=f3(i,j,k)
                    else
                        f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    end if



                !
                end do
                !
                ! direction 2 is always not periodic
                    !
                    ! check derivative on sides back and front
                    !
                    j=1
                    r0=.5*(r(i,j  ,k)+r(i,j  ,k+1))
                    r1=.5*(r(i,j+1,k)+r(i,j+1,k+1))
                    r2=.5*(r(i,j+2,k)+r(i,j+2,k+1))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))
                    !      f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg

                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j,j+2
                            do ii=i,i
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f3(i,j,k)=f3(i,j,k)
                    else
                        f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    end if

                    !
                    j=n2
                    r0=.5*(r(i,n2  ,k)+r(i,n2  ,k+1))
                    r1=.5*(r(i,n2-1,k)+r(i,n2-1,k+1))
                    r2=.5*(r(i,n2-2,k)+r(i,n2-2,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))
                    !      f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    !

                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j-2,j
                            do ii=i,i
                                if (tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if (found_solid) then
                        f3(i,j,k)=f3(i,j,k)
                    else
                        f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    end if
            !
            end do
        end do


        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    !
                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f3(i,j,k)=0.
                        if (k<n3) then
                            if (tipo(i,j,k+1)==0)f3(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do


        !     pass f3 at the border between procs
        if (myid==0) then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if (myid==nproc-1) then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else
            leftpem=leftpe
            rightpem=rightpe
        end if

        if (rightpem/=MPI_PROC_NULL) then
            call MPI_SEND(f3(1,1,kparaend),n1*n2,MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem/=MPI_PROC_NULL) then
            call MPI_RECV(f3(1,1,kparasta-1),n1*n2,MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if

        return
    end subroutine flud3

    subroutine flux1(rc,cgra1,r,tipo)
        !***********************************************************************
        ! compute explict flux on csi component:
        !
        ! convective  uc*(u,v,w)  with cenetered scheme or quick
        ! diffusive   nni*g12*d(u,v,w)/d(eta)
        ! diffusive   nni*g13*d(u,v,w)/d(zita)
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer :: kparastak,kparaendk
        !real bulk_loc
        !
        integer :: i,j,k,jj,kk,is,iss
        real :: fg,r0,r1,r2,an
        real :: rc(0:n1,1:n2,kparasta  :kparaend)
        real :: cgra1(0:n1,n2,kparasta-1:kparaend+1)
        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        real :: ravanti,rindietro1,rindietro2
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Uc*u
        ! implemented with central scheme or with quick depending on setting
        ! in point closer to wall it uses uf as ghost cell
        !
        !
        !     side left and right
        do i=1,ip
            !
            do k=kparasta,kparaend
                do j=1,n2
                    !
                    f1(0,j,k)=rc(0,j,k)*r(0,j,k)-cgra1(0,j,k)
                    f1(n1,j,k)=rc(n1,j,k)*r(n1+1,j,k)-cgra1(n1,j,k)
                !
                end do
            end do

        end do
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !
                    f1(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
                !
                end do
            end do
        end do
        !
        !
        !
        ! quick
        !
        if (insc==1) then
            !     sides 1 and 2
            do k=kparasta,kparaend
                do j=1,n2

                    i=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(2,j,k)
                        rindietro1=r(1,j,k)
                        rindietro2=ip*(2.*r(0,j,k)-r(1,j,k)) &
                            +(1.-ip)*r(0,j,k)
                    else
                        ravanti=r(1,j,k)
                        rindietro1=r(2,j,k)
                        rindietro2=r(3,j,k)
                    end if

                    if (tipo(i,j,k)==2) then
                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    end if

                    i=n1-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(n1  ,j,k)
                        rindietro1=r(n1-1,j,k)
                        rindietro2=r(n1-2,j,k)
                    else
                        ravanti=r(n1-1,j,k)
                        rindietro1=r(n1  ,j,k)
                        rindietro2=ip*(2.*r(n1+1,j,k)-r(n1,j,k)) &
                            +(1.-ip)*r(n1+1,j,k)
                    end if

                    if (tipo(i,j,k)==2) then
                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=1,n2
                    do i=2,n1-2
                        !Giulia non comprende tutti i casi, da un lato è quick dall'altro dc
                        !      if (tipo(i,j,k)==2) then
                        !     if (rc(i,j,k)>0.) then
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                        !      else
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                        !      end if
                        !         end if
                        !
                        if (rc(i,j,k)>0.) then

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i+1,j,k)==0) then !solido davanti
                                    f1(i,j,k)=0.
                                elseif (tipo(i-1,j,k)==0) then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                            end if !tipo

                        else !rc(i,j,k)<=0.

                            if (tipo(i+1,j,k)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f1(i,j,k)=0.
                                elseif (tipo(i+2,j,k)==0) then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else !solido di lato
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                            end if !tipo

                        end if !rc

                    end do
                end do
            end do
        ! SMART/QUICK modificato upwind
        elseif (insc==2) then

            ! GIULIA UCS-CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,n2
                        !
                        f1(0,j,k)=rc(0,j,k)*r(0,j,k)-cgra1(0,j,k)
                        f1(n1,j,k)=rc(n1,j,k)*r(n1+1,j,k)-cgra1(n1,j,k)
                    !
                    end do
                end do

            end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2
                    do i=ip,n1-ip
                        !
                        f1(i,j,k)=-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,n2
                    i=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(2,j,k)
                        rindietro1=r(1,j,k)
                        rindietro2=ip*(2.*r(0,j,k)-r(1,j,k)) &
                            +(1.-ip)*r(0,j,k)
                    else
                        ravanti=r(1,j,k)
                        rindietro1=r(2,j,k)
                        rindietro2=r(3,j,k)
                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )


                    i=n1-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(n1   ,j,k)
                        rindietro1=r(n1-1,j, k)
                        rindietro2=r(n1-2,j, k)
                    else
                        ravanti=r(n1-1,j,k)
                        rindietro1=r(n1  ,j,k)
                        rindietro2=ip*(2.*r(n1+1,j,k)-r(n1,j,k)) &
                            +(1.-ip)*r(n1+1,j,k)
                    end if


                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))



                            else  !(rc(i,j,k)<=0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k)) !my index +1

                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo



                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i+1,j,k)==0) then !solido davanti

                                        f1(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    elseif (tipo(i-1,j,k)==0) then

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)


                                    else ! solido sta nell'altra direzione

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib
                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i-1,j,k)==1) then !la cella prima è un ib
                                    !       if (tipo(i-2,j,k)==0) then !perché considero anche questa?
                                    !              rindietro2=r(i,j,k)

                                    !       end if
                                    !         end if

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                end if !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i+1,j,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f1(i,j,k)=0.
                                    elseif (tipo(i+2,j,k)==0) then !solido avanti
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                    else !solido di lato

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))

                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if (tipo(i+1,j,k)==1) then
                                    !       if (tipo(i+2,j,k)==0) then
                                    !             rindietro2=r(i+1,j,k)
                                    !       end if
                                    !       end if


                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        !    end if !SMART/quick modificato upwind
        ! SMART/QUICK modificato upwind
        elseif (insc==3) then

            ! GIULIA UCS-CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,n2
                        !
                        f1(0,j,k)=rc(0,j,k)*r(0,j,k)-cgra1(0,j,k)
                        f1(n1,j,k)=rc(n1,j,k)*r(n1+1,j,k)-cgra1(n1,j,k)
                    !
                    end do
                end do

            end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2
                    do i=ip,n1-ip
                        !
                        f1(i,j,k)=-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,n2
                    i=1
                    if (rc(i,j,k)>0.) then
                        rindietro1=r(1,j,k)

                    else
                        rindietro1=r(2,j,k)
                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1


                    i=n1-1
                    if (rc(i,j,k)>0.) then
                        rindietro1=r(n1-1,j, k)

                    else

                        rindietro1=r(n1  ,j,k)

                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k)<=0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=1,n2
                        do i=2,n1-2
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i+1,j,k)==0) then !solido davanti

                                        f1(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    else ! solido sta nell'altra direzione

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)


                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib
                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i-1,j,k)==1) then !la cella prima è un ib
                                    !       if (tipo(i-2,j,k)==0) then !perché considero anche questa?
                                    !              rindietro2=r(i,j,k)

                                    !       end if
                                    !         end if

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                                end if !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i+1,j,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f1(i,j,k)=0.
                                    !       elseif (tipo(i+2,j,k)==0) then !solido avanti
                                    !             f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                    else !solido di lato

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if (tipo(i+1,j,k)==1) then
                                    !       if (tipo(i+2,j,k)==0) then
                                    !             rindietro2=r(i+1,j,k)
                                    !       end if
                                    !       end if


                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        end if !SMART/quick modificato upwind

        !
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G12*DU/D(ETA)
        !
        !     sides left and right
        is=0
        iss=0

        do i=1,2*ip
            !
            do k=kparasta,kparaend
                !
                do j=2,n2-1
                    !
                    fg=.5*(r(iss,j+1,k)-r(iss,j-1,k))
                    f1(is,j,k)=-f1(is,j,k) &
                        +annit(iss,j,k)*g12(is,j,k)*fg
                !
                end do
                !
                ! direction 2 is always not periodic
                    !
                    !        check derivative at wall bottom and upper
                    !
                    j=1
                    fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j+1,k)-r(iss,j+2,k))
                    f1(is,j,k)=-f1(is,j,k) &
                        +annit(iss,j,k)*g12(is,j,k)*fg
                    !
                    j=n2
                    fg=.5*(3.*r(iss,n2,k)-4.*r(iss,n2-1,k)+r(iss,n2-2,k))
                    f1(is,j,k)=-f1(is,j,k) &
                        +annit(iss,j,k)*g12(is,j,k)*fg
                !
            !
            end do
            !
            is=n1
            iss=n1+1
        !
        end do
        !
        !
        !     into the field
        do k=kparasta,kparaend
            do i=ip,n1-ip
                !
                do j=2,n2-1
                    !
                    r1=.5*(r(i,j-1,k)+r(i+1,j-1,k))
                    r2=.5*(r(i,j+1,k)+r(i+1,j+1,k))
                    fg=.5*(r2-r1)
                    an=.5*(annit(i,j,k)+annit(i+1,j,k))
                    f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
                !
                end do
                !
                !     check derivative on sides bottom and upper
                !
                ! direction 2 is always not periodic
                    !
                    j=1
                    r0=.5*(r(i,j  ,k)+r(i+1,j  ,k))
                    r1=.5*(r(i,j+1,k)+r(i+1,j+1,k))
                    r2=.5*(r(i,j+2,k)+r(i+1,j+2,k))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    an=.5*(annit(i,j,k)+annit(i+1,j,k))
                    f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
                    !
                    j=n2
                    r0=.5*(r(i,n2  ,k)+r(i+1,n2  ,k))
                    r1=.5*(r(i,n2-1,k)+r(i+1,n2-1,k))
                    r2=.5*(r(i,n2-2,k)+r(i+1,n2-2,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i+1,j,k))
                    f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
                !
            !
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G13*DU/D(ZITA)
        !
        ! sides left and right
        !
        ! define the limit depending on periodicity in z
        !
        if (myid==0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid==nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else
            kparastak=kparasta
            kparaendk=kparaend
        end if

        is=0
        iss=0

        do i=1,2*ip
            !
            do j=1,n2
                !
                do k=kparastak,kparaendk
                    !
                    fg=.5*(r(iss,j,k+1)-r(iss,j,k-1))
                    f1(is,j,k)=f1(is,j,k) &
                        +annit(iss,j,k)*g13(is,j,k)*fg
                !
                end do
                !
                do kk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j,k+1)-r(iss,j,k+2))
                        f1(is,j,k)=f1(is,j,k) &
                            +annit(iss,j,k)*g13(is,j,k)*fg

                    end if
                    !
                    if (myid==nproc-1) then

                        k=n3
                        fg=.5*(3.*r(iss,j,n3)-4.*r(iss,j,n3-1)+r(iss,j,n3-2))
                        f1(is,j,k)=f1(is,j,k) &
                            +annit(iss,j,k)*g13(is,j,k)*fg

                    end if
                !
                end do
            !
            end do
            !
            is=n1
            iss=n1+1
        !
        end do


        !
        ! into the field
        !
        do j=1,n2
            do i=ip,n1-ip
                do k=kparastak,kparaendk
                    !
                    !
                    an=.5*(annit(i,j,k)+annit(i+1,j,k))
                    r1=.5*(r(i,j,k-1)+r(i+1,j,k-1))
                    r2=.5*(r(i,j,k+1)+r(i+1,j,k+1))
                    fg=.5*(r2-r1)
                    f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg
                !

                end do
                !

                do kk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        r0=.5*(r(i,j,k  )+r(i+1,j,k  ))
                        r1=.5*(r(i,j,k+1)+r(i+1,j,k+1))
                        r2=.5*(r(i,j,k+2)+r(i+1,j,k+2))
                        fg=.5*(-3.*r0+4.*r1-r2)
                        an=.5*(annit(i,j,k)+annit(i+1,j,k))
                        f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg

                    end if
                    !
                    if (myid==nproc-1) then

                        k=n3
                        r0=.5*(r(i,j,n3  )+r(i+1,j,n3  ))
                        r1=.5*(r(i,j,n3-1)+r(i+1,j,n3-1))
                        r2=.5*(r(i,j,n3-2)+r(i+1,j,n3-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        an=.5*(annit(i,j,k)+annit(i+1,j,k))
                        f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg

                    end if
                !
                end do

            end do

        end do


        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    !

                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f1(i,j,k)=0.
                        if (i<n1) then
                            if (tipo(i+1,j,k)==0)f1(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do

        !
        !
        !        ! integral of convective plus diffusive off diagonal
        !        !
        !        bulk_loc=0.
        !        do k=kparasta,kparaend
        !            do j=1,jy
        !                do i=1,jx
        !                    bulk_loc=bulk_loc+f1(i,j,k)-f1(i-1,j,k)
        !                end do
        !            end do
        !        end do
        !        !
        !        ! make the value known to all procs to compute the total
        !        !
        !        call MPI_ALLREDUCE(bulk_loc,bulk,1,MPI_REAL_SD, &
        !            MPI_SUM, &
        !            MPI_COMM_WORLD,ierr)
        !
        !        ! now bulk is known to all procs
        !        !
        return
    end subroutine flux1

    subroutine flux2(rc,cgra2,r,tipo)
        !***********************************************************************
        ! compute explicit flux on eta component
        !
        ! convective  vc*(u,v,w)   with central scheme or quick
        ! diffusive   nni*g21*d(u,v,w)/d(csi)
        ! diffusive   nni*g23*d(u,v,w)/d(zita)
        !
        use wallmodel_module, only: wfp3, wfp4

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer :: kparastak,kparaendk
        !real bulk_loc,bulk2

        integer :: i,j,k,ii,kk,js,jss
        real :: fg,r0,r1,r2,an
        real :: rc(n1,0:n2,kparasta:kparaend)
        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)    !0:n3+1)
        real :: cgra2(n1,0:n2,kparasta-1:kparaend+1)
        real :: ravanti,rindietro1,rindietro2
        logical :: iwfp
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Vc*u
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell

        !     sides bottom and upper
        ! direction 2 is always not periodic
            !
            do k=kparasta,kparaend
                do i=1,n1
                    !
                    f2(i,0,k)=rc(i,0,k)*r(i,0,k)-cgra2(i,0,k)
                    f2(i,n2,k)=rc(i,n2,k)*r(i,n2+1,k)-cgra2(i,n2,k)
                !
                end do
            end do

        !     into the field
        do k=kparasta,kparaend
            do j=1,n2-1
                do i=1,n1
                    !
                    f2(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j+1,k))-cgra2(i,j,k)
                !
                end do
            end do
        end do
        !
        ! quick
        !
        if (insc==1) then

            !     sides 3 and 4
            do k=kparasta,kparaend
                do i=1,n1

                    j=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,2,k)
                        rindietro1=r(i,1,k)
                        rindietro2=2.*r(i,0,k)-r(i,1,k)
                    else
                        ravanti=r(i,1,k)
                        rindietro1=r(i,2,k)
                        rindietro2=r(i,3,k)
                    end if

                    if (tipo(i,j,k)==2) then
                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    end if

                    j=n2-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,n2  ,k)
                        rindietro1=r(i,n2-1,k)
                        rindietro2=r(i,n2-2,k)
                    else
                        ravanti=r(i,n2-1,k)
                        rindietro1=r(i,n2  ,k)
                        rindietro2=2.*r(i,n2+1,k)-r(i,n2,k)
                    end if

                    if (tipo(i,j,k)==2) then
                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1-rindietro2 )
                    end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=2,n2-2
                    do i=1,n1
                        !
                        !         if (tipo(i,j,k)==2) then
                        !      if (rc(i,j,k)>0.) then
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                        !      else
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                        !      end if
                        !      end if
                        !

                        if (rc(i,j,k)>0.) then

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i,j+1,k)==0) then !solido davanti
                                    f2(i,j,k)=0.
                                elseif (tipo(i,j-1,k)==0) then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                            end if !tipo

                        else !rc(i,j,k)<=0.
                            if (tipo(i,j+1,k)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f2(i,j,k)=0.
                                elseif (tipo(i,j+2,k)==0) then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)

                                else !solido di lato
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))


                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                            end if !tipo

                        end if !rc
                    end do
                end do
            end do

        ! quick upwind
        elseif (insc==2) then

            ! GIULIA UCS-CGRA
            !     side left and right
            ! direction 2 is always not periodic
                !
                do k=kparasta,kparaend
                    do i=1,n1
                        !
                        f2(i,0,k)=rc(i,0,k)*r(i,0,k)-cgra2(i,0,k)
                        f2(i,n2,k)=rc(i,n2,k)*r(i,n2+1,k)-cgra2(i,n2,k)
                    !
                    end do
                end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2-1
                    do i=1,n1
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,n1
                    j=1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,2,k)
                        rindietro1=r(i,1,k)
                        rindietro2=2.*r(i,0,k)-r(i,1,k)
                    else
                        ravanti=r(i,1,k)
                        rindietro1=r(i,2,k)
                        rindietro2=r(i,3,k)
                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )

                    j=n2-1
                    if (rc(i,j,k)>0.) then
                        ravanti=r(i,n2  , k)
                        rindietro1=r(i,n2-1, k)
                        rindietro2=r(i,n2-2, k)
                    else
                        ravanti=r(i,n2-1,k)
                        rindietro1=r(i,n2  ,k)
                        rindietro2=2.*r(i,n2+1,k)-r(i,n2,k)
                    end if


                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                            else  !(rc(i,j,k)<=0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k)) !my index +1


                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo


                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j+1,k)==0) then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo
                                    elseif (tipo(i,j-1,k)==0) then

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i-1,j,k)==1) then !la cella prima è un ib
                                    !       if (tipo(i-2,j,k)==0) then !perché considero anche questa?
                                    !              rindietro2=r(i,j,k)

                                    !       end if
                                    !         end if


                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                                end if !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j+1,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f2(i,j,k)=0.
                                    elseif (tipo(i,j+2,k)==0) then !solido avanti
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                    else !solido di lato

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if (tipo(i+1,j,k)==1) then
                                    !       if (tipo(i+2,j,k)==0) then
                                    !             rindietro2=r(i+1,j,k)
                                    !       end if
                                    !       end if


                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        !    end if !SMART/   quick modificato upwind

        ! quick upwind
        elseif (insc==3) then

            ! GIULIA UCS-CGRA
            !     side left and right
            ! direction 2 is always not periodic
                !
                do k=kparasta,kparaend
                    do i=1,n1
                        !
                        f2(i,0,k)=rc(i,0,k)*r(i,0,k)-cgra2(i,0,k)
                        f2(i,n2,k)=rc(i,n2,k)*r(i,n2+1,k)-cgra2(i,n2,k)
                    !
                    end do
                end do
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,n2-1
                    do i=1,n1
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,n1
                    j=1
                    if (rc(i,j,k)>0.) then

                        rindietro1=r(i,1,k)

                    else

                        rindietro1=r(i,2,k)

                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

                    j=n2-1
                    if (rc(i,j,k)>0.) then

                        rindietro1=r(i,n2-1, k)

                    else

                        rindietro1=r(i,n2  ,k)

                    end if


                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)


                            else  !(rc(i,j,k)<=0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)


                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=2,n2-2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j+1,k)==0) then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    else ! solido sta nell'altra direzione

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    end if !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i-1,j,k)==1) then !la cella prima è un ib
                                    !       if (tipo(i-2,j,k)==0) then !perché considero anche questa?
                                    !              rindietro2=r(i,j,k)

                                    !       end if
                                    !         end if

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                                end if !tipo

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j+1,k)==1) then
                                    if (tipo(i,j,k)==0) then !solido indietro
                                        f2(i,j,k)=0.
                                    elseif (tipo(i,j+2,k)==0) then !solido avanti
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                    else !solido di lato

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                                    end if !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if (tipo(i+1,j,k)==1) then
                                    !       if (tipo(i+2,j,k)==0) then
                                    !             rindietro2=r(i+1,j,k)
                                    !       end if
                                    !       end if

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        end if !SMART/    quick modificato upwind
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G21*DU/D(CSI)
        !
        !     sides bottom and upper
        !
        js=0
        jss=0
        iwfp=wfp3

        do j=1,2
            !
            do k=kparasta,kparaend
                !

                do  i=1+ip,n1-ip
                    fg=.5*(r(i+1,jss,k)-r(i-1,jss,k))
                    !         WALL MODEL : all in flucn
                    if (iwfp) then
                        f2(i,js,k)=-f2(i,js,k)
                    else
                        f2(i,js,k)=-f2(i,js,k) &
                            +annit(i,jss,k)*g21(i,js,k)*fg
                    end if
                end do

                !
                do ii=1,ip
                    !     check derivative on sides left and right
                    i=1
                    fg=.5*(-3.*r(i,jss,k)+4.*r(i+1,jss,k)-r(i+2,jss,k))
                    f2(i,js,k)=-f2(i,js,k) &
                        +annit(i,jss,k)*g21(i,js,k)*fg
                    i=n1
                    fg=.5*(3.*r(n1,jss,k)-4.*r(n1-1,jss,k)+r(n1-2,jss,k))
                    f2(i,js,k)=-f2(i,js,k) &
                        +annit(i,jss,k)*g21(i,js,k)*fg
                end do


            !
            end do
            !
            js=n2
            jss=n2+1
            iwfp=wfp4
        !
        end do
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=1,n2-1
                !
                do i=1+ip,n1-ip
                    !
                    r1=.5*(r(i-1,j,k)+r(i-1,j+1,k))
                    r2=.5*(r(i+1,j,k)+r(i+1,j+1,k))
                    fg=.5*(r2-r1)
                    an=.5*(annit(i,j,k)+annit(i,j+1,k))
                    f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*an*fg
                !
                end do
                !
                do ii=1,ip
                    !
                    ! check derivative on sides left and right
                    !
                    i=1
                    r0=.5*(r(i  ,j,k)+r(i  ,j+1,k))
                    r1=.5*(r(i+1,j,k)+r(i+1,j+1,k))
                    r2=.5*(r(i+2,j,k)+r(i+2,j+1,k))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    an=.5*(annit(i,j,k)+annit(i,j+1,k))
                    f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*an*fg
                    !
                    i=n1
                    r0=.5*(r(n1  ,j,k)+r(n1  ,j+1,k))
                    r1=.5*(r(n1-1,j,k)+r(n1-1,j+1,k))
                    r2=.5*(r(n1-2,j,k)+r(n1-2,j+1,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i,j+1,k))
                    f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*an*fg
                !
                end do
            !
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G23*DU/D(ZITA)
        !
        ! sides bottom and upper
        !
        ! define computation limit depending on periodicity in z
        !
        if (myid==0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid==nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else
            kparastak=kparasta
            kparaendk=kparaend
        end if

        js=0
        jss=0
        iwfp=wfp3
        do j=1,2
            !
            do i=1,n1
                !
                do k=kparastak,kparaendk
                    !
                    fg=.5*(r(i,jss,k+1)-r(i,jss,k-1))

                    !           WALL MODEL: all in flucn
                    if (iwfp) then
                        f2(i,js,k)=f2(i,js,k)
                    else
                        f2(i,js,k)=f2(i,js,k) &
                            +annit(i,jss,k)*g23(i,js,k)*fg
                    end if
                !
                end do
                !

                do kk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        fg=.5*(-3.*r(i,jss,k)+4.*r(i,jss,k+1)-r(i,jss,k+2))
                        f2(i,js,k)=f2(i,js,k) &
                            +annit(i,jss,k)*g23(i,js,k)*fg

                    end if
                    !
                    if (myid==nproc-1) then

                        k=n3
                        fg=.5*(3.*r(i,jss,n3)-4.*r(i,jss,n3-1)+r(i,jss,n3-2))
                        f2(i,js,k)=f2(i,js,k) &
                            +annit(i,jss,k)*g23(i,js,k)*fg

                    end if
                !
                end do
            !
            end do
            !
            js=n2
            jss=n2+1
            iwfp=wfp4
        !
        end do
        !
        ! into the field
        !
        do j=1,n2-1
            do i=1,n1
                !
                do k=kparastak,kparaendk
                    !
                    r1=.5*(r(i,j,k-1)+r(i,j+1,k-1))
                    r2=.5*(r(i,j,k+1)+r(i,j+1,k+1))
                    fg=.5*(r2-r1)
                    an=.5*(annit(i,j,k)+annit(i,j+1,k))
                    f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg
                !
                end do
                !
                do kk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid==0) then

                        k=1
                        r0=.5*(r(i,j,k  )+r(i,j+1,k  ))
                        r1=.5*(r(i,j,k+1)+r(i,j+1,k+1))
                        r2=.5*(r(i,j,k+2)+r(i,j+1,k+2))
                        fg=.5*(-3.*r0+4.*r1-r2)
                        an=.5*(annit(i,j,k)+annit(i,j+1,k))
                        f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg

                    end if
                    !
                    if (myid==nproc-1) then

                        k=n3
                        r0=.5*(r(i,j,n3  )+r(i,j+1,n3  ))
                        r1=.5*(r(i,j,n3-1)+r(i,j+1,n3-1))
                        r2=.5*(r(i,j,n3-2)+r(i,j+1,n3-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        an=.5*(annit(i,j,k)+annit(i,j+1,k))
                        f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg

                    end if
                !
                end do
            !
            end do
        end do


        do k=kparasta,kparaend
            do j=1,n2-1
                do i=1,n1
                    !

                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f2(i,j,k)=0.
                        if (j<n2) then
                            if (tipo(i,j+1,k)==0)f2(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do


        !
        !        ! integral on f2
        !        bulk_loc=0.
        !        !
        !        do k=kparasta,kparaend
        !            do j=1,jy
        !                do i=1,jx
        !                    bulk_loc=bulk_loc+f2(i,j,k)-f2(i,j-1,k)
        !                end do
        !            end do
        !        end do
        !
        !        ! make the value known to all procs
        !        call MPI_ALLREDUCE(bulk_loc,bulk2,1,MPI_REAL_SD, &
        !            MPI_SUM, &
        !            MPI_COMM_WORLD,ierr)
        !
        !        ! now bulk is known by all procs
        !
        !        bulk=bulk+bulk2

        return
    end subroutine flux2

    subroutine flux3(rc,cgra3,r,tipo)
        !***********************************************************************
        ! compute explicit flux on component zita
        !
        ! convective  wc*(u,v,w)   with central schema or quick
        ! diffusive   nni*g31*d(u,v,w)/d(csi)
        ! diffusive   nni*g32*d(u,v,w)/d(zita)
        !
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer :: ierr
        integer :: kparastal,kparaendl
        integer :: kparastall,kparaendll
        !real bulk_loc,bulk3
        integer :: i,j,k,ii,jj,ks,kss
        real :: fg,r0,r1,r2,an
        real :: rc(n1,n2,kparasta-1:kparaend) !0:n3)
        real :: r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real :: cgra3(n1,n2,kparasta-1:kparaend+1)
        real :: ravanti,rindietro1,rindietro2
        integer :: status(MPI_STATUS_SIZE)
        integer :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM  Wc*u:
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell
        !
        !     sides back and front
        do k=1,kp

            do j=1,n2
                do i=1,n1
                    !
                    if (myid==0) then
                        f3(i,j,0)=rc(i,j,0)*r(i,j,0   )-cgra3(i,j,0 )
                    end if

                    if (myid==nproc-1) then
                        f3(i,j,n3)=rc(i,j,n3)*r(i,j,n3+1)-cgra3(i,j,n3)
                    end if
                !
                end do
            end do
        !
        end do
        !
        !     into the field
        !
        if (myid==0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid==nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else
            kparastal=kparasta
            kparaendl=kparaend
        end if


        do k=kparastal,kparaendl
            do j=1,n2
                do i=1,n1
                    !
                    f3(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
                !
                end do
            end do
        end do

        ! quick
        if (insc==1) then

            if (myid==0) then
                kparastall=kparasta +2
                kparaendll=kparaend
            else if (myid==nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            end if

            !     sides 5 and 6

            if (myid==0) then
                do j=1,n2
                    do i=1,n1

                        k=1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,2)
                            rindietro1=r(i,j,1)
                            rindietro2=kp*(2.*r(i,j,0)-r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti=r(i,j,1)
                            rindietro1=r(i,j,2)
                            rindietro2=r(i,j,3)
                        end if

                        if (tipo(i,j,k)==2) then
                            f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                *.125*( -ravanti +2.*rindietro1-rindietro2 )
                        end if

                    end do
                end do
            end if

            if (myid==nproc-1) then
                do j=1,n2
                    do i=1,n1

                        k=n3-1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,n3  )
                            rindietro1=r(i,j,n3-1)
                            rindietro2=r(i,j,n3-2)
                        else
                            ravanti=r(i,j,n3-1)
                            rindietro1=r(i,j,n3  )
                            rindietro2=kp*(2.*r(i,j,n3+1)-r(i,j,n3)) &
                                +(1.-kp)*r(i,j,n3+1)
                        end if

                        if (tipo(i,j,k)==2) then
                            f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                *.125*( -ravanti +2.*rindietro1-rindietro2 )
                        end if
                    end do
                end do
            end if

            !     into the field
            do k=kparastall,kparaendll
                do j=1,n2
                    do i=1,n1
                        !
                        !      if (tipo(i,j,k)==2) then
                        !      if (rc(i,j,k)>0.) then
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                        !      else
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                        !      end if
                        !      end if

                        if (rc(i,j,k)>0.) then

                            if (tipo(i,j,k)==1) then !ib cell
                                if (tipo(i,j,k+1)==0) then !solido davanti
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k-1)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))

                                end if !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                            end if !tipo

                        else !rc(i,j,k)<=0.

                            if (tipo(i,j,k+1)==1) then
                                if (tipo(i,j,k)==0) then !solido indietro
                                    f3(i,j,k)=0.
                                elseif (tipo(i,j,k+2)==0) then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                                end if !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                            end if !tipo

                        end if !rc
                    end do
                end do
            end do

        ! SMART/quick modificato upwind
        elseif (insc==2) then

            ! GIULIA UCS-CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,n2
                    do i=1,n1
                        !
                        if (myid==0) then
                            f3(i,j,0)=rc(i,j,0 )*r(i,j,0   )-cgra3(i,j,0 )
                        end if

                        if (myid==nproc-1) then
                            f3(i,j,n3)=rc(i,j,n3)*r(i,j,n3+1)-cgra3(i,j,n3)
                        end if
                    !
                    end do
                end do
            !
            end do

            if (myid==0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid==nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            end if


            do k=kparastal,kparaendl
                do j=1,n2
                    do i=1,n1
                        !
                        f3(i,j,k)=-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid==0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid==nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            end if

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if (myid==0) then
                do j=1,n2
                    do i=1,n1
                        k=1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,2)
                            rindietro1=r(i,j,1)
                            rindietro2=kp*(2.*r(i,j,0)-r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti=r(i,j,1)
                            rindietro1=r(i,j,2)
                            rindietro2=r(i,j,3)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )


                    end do
                end do

            end if !myid==0

            if (myid==nproc-1) then
                do j=1,n2
                    do i=1,n1

                        k=n3-1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,n3  )
                            rindietro1=r(i,j,n3-1)
                            rindietro2=r(i,j,n3-2)
                        else
                            ravanti=r(i,j,n3-1)
                            rindietro1=r(i,j,n3  )
                            rindietro2=kp*(2.*r(i,j,n3+1)-r(i,j,n3)) &
                                +(1.-kp)*r(i,j,n3+1)
                        end if


                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*( 3.*ravanti +6.*rindietro1-rindietro2 )

                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                            else  !(rc(i,j,k)<=0.)

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2)) !my index +1


                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j,k+1)==0) then !solido davanti

                                        f3(i,j,k)=0.

                                    elseif (tipo(i,j,k-1)==0) then !solido dietro

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                                    end if !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k-1)==1) then !la cella prima è un ib
                                    !       if (tipo(i,j,k-2)==0) then !perché considero anche questa?
                                    !         rindietro2=r(i,j,k)
                                    !       end if
                                    !         end if !fine k-1 IB


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))


                                end if !tip02

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j,k+1)==1) then !ib
                                    if (tipo(i,j,k)==0) then !solido avanti
                                        f3(i,j,k)=0.
                                    elseif (tipo(i,j,k+2)==0) then !solido indietro
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)
                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))


                                    end if !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k+1)==1) then
                                    !       if (tipo(i,j,k+2)==0) then
                                    !       rindietro2=r(i,j,k+1)
                                    !       end if
                                    !         end if

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        !    end if !SMART
        !-----SMART /quick modificato upwind
        !  upwind
        elseif (insc==3) then

            ! GIULIA UCS-CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,n2
                    do i=1,n1
                        !
                        if (myid==0) then
                            f3(i,j,0)=rc(i,j,0)*r(i,j,0)-cgra3(i,j,0 )
                        end if

                        if (myid==nproc-1) then
                            f3(i,j,n3)=rc(i,j,n3)*r(i,j,n3+1)-cgra3(i,j,n3)
                        end if
                    !
                    end do
                end do
            !
            end do

            if (myid==0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid==nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            end if


            do k=kparastal,kparaendl
                do j=1,n2
                    do i=1,n1
                        !
                        f3(i,j,k)=-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid==0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid==nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            end if

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if (myid==0) then
                do j=1,n2
                    do i=1,n1
                        k=1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,2)
                            rindietro1=r(i,j,1)
                            rindietro2=kp*(2.*r(i,j,0)-r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti=r(i,j,1)
                            rindietro1=r(i,j,2)
                            rindietro2=r(i,j,3)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1


                    end do
                end do

            end if !myid==0

            if (myid==nproc-1) then
                do j=1,n2
                    do i=1,n1

                        k=n3-1
                        if (rc(i,j,k)>0.) then
                            ravanti=r(i,j,n3  )
                            rindietro1=r(i,j,n3-1)
                            rindietro2=r(i,j,n3-2)
                        else
                            ravanti=r(i,j,n3-1)
                            rindietro1=r(i,j,n3  )
                            rindietro2=kp*(2.*r(i,j,n3+1)-r(i,j,n3)) &
                                +(1.-kp)*r(i,j,n3+1)
                        end if


                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1

                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k)<=0.)

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)


                            end if !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,n2
                        do i=1,n1
                            if (rc(i,j,k)>0.) then
                                if (tipo(i,j,k)==1) then !ib cell
                                    if (tipo(i,j,k+1)==0) then !solido davanti

                                        f3(i,j,k)=0.


                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    end if !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k-1)==1) then !la cella prima è un ib
                                    !       if (tipo(i,j,k-2)==0) then !perché considero anche questa?
                                    !         rindietro2=r(i,j,k)
                                    !       end if
                                    !         end if !fine k-1 IB


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)


                                end if !tip02

                            else !rc(i,j,k)<=0.

                                if (tipo(i,j,k+1)==1) then !ib
                                    if (tipo(i,j,k)==0) then !solido avanti
                                        f3(i,j,k)=0.

                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)


                                    end if !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if (tipo(i,j,k+1)==1) then
                                    !       if (tipo(i,j,k+2)==0) then
                                    !       rindietro2=r(i,j,k+1)
                                    !       end if
                                    !         end if

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

                                end if !tipo

                            end if !rc

                        end do
                    end do
                end do

            end if !body yes/no
        end if !SMART
        !-----SMART /quick modificato upwind

        !----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G31*DU/D(CSI)
        !
        !     sides back and front
        ks=0
        kss=0
        do k=1,2*kp

            if (myid==(k-1)*(nproc-1)) then

                do j=1,n2
                    !
                    do i=1+ip,n1-ip
                        !
                        fg=.5*(r(i+1,j,kss)-r(i-1,j,kss))
                        f3(i,j,ks)=-f3(i,j,ks) &
                            +annit(i,j,kss)*g31(i,j,ks)*fg
                    !
                    end do
                    !
                    do ii=1,ip
                        !
                        ! check derivative on sides left and right
                        !
                        i=1
                        fg=.5*(-3.*r(i,j,kss)+4.*r(i+1,j,kss)-r(i+2,j,kss))
                        f3(i,j,ks)=-f3(i,j,ks) &
                            +annit(i,j,kss)*g31(i,j,ks)*fg
                        !
                        i=n1
                        fg=.5*(3.*r(n1,j,kss)-4.*r(n1-1,j,kss)+r(n1-2,j,kss))
                        f3(i,j,ks)=-f3(i,j,ks) &
                            +annit(i,j,kss)*g31(i,j,ks)*fg
                    !
                    end do
                !
                end do
            !
            end if

            ks=n3
            kss=n3+1
        !
        end do
        !
        !     into the field
        do k=kparastal,kparaendl
            do j=1,n2
                !
                do i=1+ip,n1-ip
                    !
                    r1=.5*(r(i-1,j,k)+r(i-1,j,k+1))
                    r2=.5*(r(i+1,j,k)+r(i+1,j,k+1))
                    fg=.5*(r2-r1)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*an*fg
                !
                end do
                !
                do ii=1,ip
                    !
                    ! check derivative on sides left and right
                    !
                    i=1
                    r0=.5*(r(i  ,j,k)+r(i  ,j,k+1))
                    r1=.5*(r(i+1,j,k)+r(i+1,j,k+1))
                    r2=.5*(r(i+2,j,k)+r(i+2,j,k+1))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*an*fg
                    !
                    i=n1
                    r0=.5*(r(n1  ,j,k)+r(n1  ,j,k+1))
                    r1=.5*(r(n1-1,j,k)+r(n1-1,j,k+1))
                    r2=.5*(r(n1-2,j,k)+r(n1-2,j,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*an*fg
                !
                end do
            !
            end do
        end do
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G32*DU/D(ZITA)
        !
        !     sides back and front
        !
        ks=0
        kss=0

        do k=1,2*kp

            if (myid==(k-1)*(nproc-1)) then

                do i=1,n1
                    !
                    do j=2,n2-1
                        !
                        fg=.5*(r(i,j+1,kss)-r(i,j-1,kss))
                        f3(i,j,ks)=f3(i,j,ks) &
                            +annit(i,j,kss)*g32(i,j,ks)*fg
                    !
                    end do
                    !
                    ! direction 2 is always not periodic
                        !
                        ! check derivative on sides back and front
                        !
                        j=1
                        fg=.5*(-3.*r(i,j,kss)+4.*r(i,j+1,kss)-r(i,j+2,kss))
                        f3(i,j,ks)=f3(i,j,ks) &
                            +annit(i,j,kss)*g32(i,j,ks)*fg
                        !
                        j=n2
                        fg=.5*(3.*r(i,n2,kss)-4.*r(i,n2-1,kss)+r(i,n2-2,kss))
                        f3(i,j,ks)=f3(i,j,ks) &
                            +annit(i,j,kss)*g32(i,j,ks)*fg
                    !
                !
                end do
            !
            end if

            ks=n3
            kss=n3+1
        !
        end do

        !
        !     into the field
        do k=kparastal,kparaendl
            do i=1,n1
                !
                do j=2,n2-1
                    !
                    r1=.5*(r(i,j-1,k)+r(i,j-1,k+1))
                    r2=.5*(r(i,j+1,k)+r(i,j+1,k+1))
                    fg=.5*(r2-r1)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
                !
                end do
                !
                ! direction 2 is always not periodic
                    !
                    ! check derivative on sides back and front
                    !
                    j=1
                    r0=.5*(r(i,j  ,k)+r(i,j  ,k+1))
                    r1=.5*(r(i,j+1,k)+r(i,j+1,k+1))
                    r2=.5*(r(i,j+2,k)+r(i,j+2,k+1))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
                    !
                    j=n2
                    r0=.5*(r(i,n2  ,k)+r(i,n2  ,k+1))
                    r1=.5*(r(i,n2-1,k)+r(i,n2-1,k+1))
                    r2=.5*(r(i,n2-2,k)+r(i,n2-2,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
                !
            !
            end do
        end do
        !

        do j=1,n2
            do i=1,n1
                do k=kparastal,kparaendl
                    !
                    if (bodyforce) then
                        if (tipo(i,j,k)==0)f3(i,j,k)=0.
                        if (k<n3) then
                            if (tipo(i,j,k+1)==0)f3(i,j,k)=0.
                        end if
                    end if
                end do
            end do
        end do
        !     pass f3 at the border between procs
        if (myid==0) then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if (myid==nproc-1) then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else
            leftpem=leftpe
            rightpem=rightpe
        end if

        if (rightpem/=MPI_PROC_NULL) then
            call MPI_SEND(f3(1,1,kparaend),n1*n2,MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
        end if
        if (leftpem/=MPI_PROC_NULL) then
            call MPI_RECV(f3(1,1,kparasta-1),n1*n2,MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
        end if

        !     integral for f3
        !        bulk_loc=0.
        !
        !        do k=kparasta,kparaend
        !            do j=1,jy
        !                do i=1,jx
        !                    bulk_loc=bulk_loc+f3(i,j,k)-f3(i,j,k-1)
        !                end do
        !            end do
        !        end do
        !
        !        !     make the value known to all procs
        !
        !        call MPI_ALLREDUCE(bulk_loc,bulk3,1,MPI_REAL_SD, &
        !            MPI_SUM, &
        !            MPI_COMM_WORLD,ierr)
        !
        !        !    now bulk is known by all procs
        !
        !        bulk=bulk+bulk3
        !
        return

    end subroutine flux3

end module flu_module
