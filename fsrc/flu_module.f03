module flu_module

    use mysending
    use myarrays_metri3
    use myarrays_velo3
    use myarrays_wallmodel, only: att_mod_par, eseguo34, u_t, utangente
    use myarrays_ibm, only: bodyforce
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    real :: bulk

contains

    subroutine flucn(r,visualizzo,coef_wall,ti,tipo,tau_wind)
        !***********************************************************************
        ! compute explicit diffusive term
        !
        ! NNI*G11*D(u,v,w)/D(csi)+
        ! NNI*G22*D(u,v,w)/D(eta)+
        ! NNI*G33*D(u,v,w)/D(zita)
        !
        !-----------------------------------------------------------------------

        use inflow_module, only: areola3,areola4

        implicit none

        !     array declaration
        integer ierr,status(MPI_STATUS_SIZE)
        integer m
        integer kparastal,kparaendl
        double precision bulk_loc,bulkn
        integer i,j,k,ii,jj,kk,ktime
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        integer visualizzo!,windyes
        !
        integer req1,req2
        integer istatus(MPI_STATUS_SIZE)

        ! Giulia modificavento:
        !      real vel_tau(0:n1+1,kparasta-1:kparaend+1)
        real, optional, intent(in out) :: tau_wind(0:,kparasta-1:)!(0:n1+1,kparasta-1:kparaend+1)
        real rhow
        ! Giulia modificavento:

        integer iwall,ieseguo
        integer,allocatable::prov(:,:,:)
        integer coef_wall
        real ti

        ! giulia aggiungo questi per eliminare i flussi tra i solidi
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! term nni*g11*d/d(csi)
        !
        do ii=1,ip
            !
            !     side left
            do k=kparasta,kparaend
                do j=1,jy
                    f1(0,j,k)=annit(0,j,k)*g11(0,j,k)*&
                        (-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
                end do
            end do
            !
            !     side right
            do k=kparasta,kparaend
                do j=1,jy
                    f1(jx,j,k)=annit(jx+1,j,k)*g11(jx,j,k)*&
                        &   (8.*r(jx+1,j,k)-9.*r(jx,j,k)+r(jx-1,j,k))/3.
                end do
            end do
        !
        enddo
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    f1(i,j,k)=.5*(annit(i,j,k)+annit(i+1,j,k))&
                        &          *g11(i,j,k)*(r(i+1,j,k)-r(i,j,k))
                !
                enddo
            enddo
        enddo

        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
                        if(i.lt.jx)then
                            if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
                        endif
                    end if
                enddo
            enddo
        enddo

        !
        !-----------------------------------------------------------------------
        ! term nni*g22*d/d(eta)
        !
        do jj=1,jp
            !
            !     side bottom
            !
            !     for direction 2, allthough coef_wall=1, the wall model is switched off
            if(eseguo34==0.and.coef_wall.eq.1)then
                allocate(prov(jx,2,kparasta:kparaend))
                prov=att_mod_par
                att_mod_par=0
            end if

            !     wall model off
            if(coef_wall.eq.0)then
                do k=kparasta,kparaend
                    do i=1,jx
                        f2(i,0,k)=annitV(i,0,k)*g22(i,0,k)*&
                            &              (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                    end do
                end do
            !     wall model on
            elseif(coef_wall.eq.1)then
                do k=kparasta,kparaend
                    do i=1,jx
                        do ii=1,1-att_mod_par(i,1,k)
                            f2(i,0,k)=annitV(i,0,k)*g22(i,0,k)*&
                                &              (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                        end do
                        do ii=1,att_mod_par(i,1,k)
                            f2(i,0,k)=u_t(i,1,k)*u_t(i,1,k)*areola3(i,k)*&
                                &                 ((2-eseguo34)*u(i,1,k)&
                                &              -(1-eseguo34)*w(i,1,k))/utangente(i,1,k)
                        end do

                    end do
                end do
            end if
	
            !         if(coef_wall ==1 .and. imoist==1 .and. eseguo34.ne.0)then
            !            call readtau(ti,r)
            !         end if
	
            !
            ! side upper
            !
            if(present(tau_wind))then

                ! Giulia modificavento:
                rhow = 1000.  ! SMELLS LIKE MAGIC
                do k=kparasta,kparaend
                    do i=1,jx
                        !          f2(i,jy,k)=abs(vel_tau(i,k))*vel_tau(i,k)*areola4(i,k)
                        f2(i,jy,k)=tau_wind(i,k)*areola4(i,k)/rhow
                    end do
                end do
            ! Giulia modificavento:

            else

                !     wall model off
                if(coef_wall.eq.0)then
                    do k=kparasta,kparaend
                        do i=1,jx
                            f2(i,jy,k)=annitV(i,jy+1,k)*g22(i,jy,k)*&
                                &              (8.*r(i,jy+1,k)-9.*r(i,jy,k)+r(i,jy-1,k))/3.
                        end do
                    end do

                !      wall model on
                elseif(coef_wall.eq.1)then
                    do k=kparasta,kparaend
                        do i=1,jx
                            do ii=1,1-att_mod_par(i,2,k)
                                f2(i,jy,k)=annitV(i,jy+1,k)*g22(i,jy,k)*&
                                    &              (8.*r(i,jy+1,k)-9.*r(i,jy,k)+r(i,jy-1,k))/3.
                            end do
                            do ii=1,att_mod_par(i,2,k)
                                f2(i,jy,k)=-u_t(i,2,k)*u_t(i,2,k)*areola4(i,k)*&
                                    &                 ((2-eseguo34)*u(i,jy,k)&
                                    &                 -(1-eseguo34)*w(i,jy,k))/utangente(i,2,k)
                            enddo
                        end do
                    end do
                end if

                ! wall model on but switched for direction2
                if(eseguo34==0.and.coef_wall.eq.1)then
                    att_mod_par=prov
                    deallocate(prov)
                end if

            end if  !end if windyes
        !
        enddo
        !
        !     into the field
        do k=kparasta,kparaend
            do i=1,jx
                do j=jp,jy-jp
                    !
                    f2(i,j,k)=.5*(annitV(i,j,k)+annitV(i,j+1,k))&
                        &              *g22(i,j,k)*(r(i,j+1,k)-r(i,j,k))
                !
                enddo
            enddo
        enddo
        !
        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx
                    !

                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
                        if(j.lt.jy)then
                            if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! term nni*g33d/d(zita)
        !
        do kk=1,kp
            !
            !     side back
            if (myid.eq.0) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,0)=annit(i,j,0)*g33(i,j,0)*&
                            &     (-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
                    end do
                end do
  
            endif
            !
            !     side front
            if (myid.eq.nproc-1) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,jz)=annit(i,j,jz+1)*g33(i,j,jz)*&
                            &   (8.*r(i,j,jz+1)-9.*r(i,j,jz)+r(i,j,jz-1))/3.
                    end do
                end do
 
            endif
        !
        enddo
        !
        !     into the field
        if (myid.eq.0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid.eq.nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else
            kparastal=kparasta
            kparaendl=kparaend
        endif

        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    !
                    f3(i,j,k)=.5*(annit(i,j,k)+annit(i,j,k+1))&
                        &           *g33(i,j,k)*(r(i,j,k+1)-r(i,j,k))
                !
                enddo
            enddo
        enddo

        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
                        if(k.lt.jz)then
                            if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo
        !
        ! pass f3 at the border to all procs
        !
        if(myid.eq.0)then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if(myid.eq.nproc-1)then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else
            leftpem=leftpe
            rightpem=rightpe
        endif

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(f3(1,1,kparaend),jx*jy,MPI_REAL_SD,&
                &                 rightpem ,tagrs,MPI_COMM_WORLD,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(f3(1,1,kparasta-1),jx*jy,MPI_REAL_SD,&
                &                 leftpem  ,taglr,MPI_COMM_WORLD,status,ierr)
        endif

        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif

        !     integral
        bulk_loc=0.
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    bulk_loc=bulk_loc+f3(i,j,k)-f3(i  ,j  ,k-1)+&
                        &                  f2(i,j,k)-f2(i  ,j-1,k  )+&
                        &                  f1(i,j,k)-f1(i-1,j  ,k  )
                end do
            end do
        end do

        ! make the value known to all procs

        call MPI_ALLREDUCE(bulk_loc,bulkn,1,MPI_DOUBLE_PRECISION,&
            &                   MPI_SUM,&
            &                   MPI_COMM_WORLD,ierr)

        bulk=bulk+bulkn

    end subroutine flucn

    subroutine flu_turbo()
        !***********************************************************************
        ! compute explicit terms for turbulence model like
        ! d/dcsi(Uf*annit), periodic version
        !
        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,kper
        !
        integer ierr
        integer kparaendp

        !-----------------------------------------------------------------------
        ! compute controvariant flux Uf (jordan notation)
        !
        !     side left and right, periodic
        do k=kparasta,kparaend
            do j=1,jy
                !
                cgra1(0 ,j,k)=csx(0 ,j,k)*.5*(gra1(jx,j,k)+gra1(1,j,k))+ &
                    csy(0 ,j,k)*.5*(gra2(jx,j,k)+gra2(1,j,k))+ &
                    csz(0 ,j,k)*.5*(gra3(jx,j,k)+gra3(1,j,k))- &
                    cgra1(0,j,k)
                cgra1(0 ,j,k)=(1-ip)*cgra1(0 ,j,k)

                cgra1(jx,j,k)=csx(jx,j,k)*.5*(gra1(jx,j,k)+gra1(1,j,k))+ &
                    csy(jx,j,k)*.5*(gra2(jx,j,k)+gra2(1,j,k))+ &
                    csz(jx,j,k)*.5*(gra3(jx,j,k)+gra3(1,j,k))- &
                    cgra1(jx,j,k)
                cgra1(jx,j,k)=(1-ip)*cgra1(jx,j,k)
            !
            end do
        end do
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx-1
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
            do i=1,jx

                cgra2(i,0,k)=etx(i,0,k)*.5*(gra1(i,jy,k)+gra1(i,1,k))+ &
                    ety(i,0,k)*.5*(gra2(i,jy,k)+gra2(i,1,k))+ &
                    etz(i,0,k)*.5*(gra3(i,jy,k)+gra3(i,1,k))- &
                    cgra2(i,0,k)
                cgra2(i,0,k)=(1-jp)*cgra2(i,0,k)

                cgra2(i,jy,k)=etx(i,jy,k)*.5*(gra1(i,jy,k)+gra1(i,1,k))+ &
                    ety(i,jy,k)*.5*(gra2(i,jy,k)+gra2(i,1,k))+ &
                    etz(i,jy,k)*.5*(gra3(i,jy,k)+gra3(i,1,k))- &
                    cgra2(i,jy,k)
                cgra2(i,jy,k)=(1-jp)*cgra2(i,jy,k)

            end do
        end do

        !     into the field
        do k=kparasta,kparaend
            do j=1,jy-1
                do i=1,jx
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
        do j=1,jy
            do i=1,jx
                !
                if (myid.eq.0) then

                    cgra3(i,j,0 )=ztx(i,j,0 ) &
                        *.5*(gra1(i,j,1)+gra1_appoggio(i,j,jz))+ &
                        zty(i,j,0 ) &
                        *.5*(gra2(i,j,1)+gra2_appoggio(i,j,jz))+ &
                        ztz(i,j,0 ) &
                        *.5*(gra3(i,j,1)+gra3_appoggio(i,j,jz))- &
                        cgra3(i,j,0)
                    cgra3(i,j,0 )=(1-kp)*cgra3(i,j,0 )

                else if (myid.eq.nproc-1) then

                    cgra3(i,j,jz)=ztx(i,j,jz) &
                        *.5*(gra1_appoggio(i,j,1)+gra1(i,j,jz))+ &
                        zty(i,j,jz) &
                        *.5*(gra2_appoggio(i,j,1)+gra2(i,j,jz))+ &
                        ztz(i,j,jz) &
                        *.5*(gra3_appoggio(i,j,1)+gra3(i,j,jz))- &
                        cgra3(i,j,jz)
                    cgra3(i,j,jz)=(1-kp)*cgra3(i,j,jz)

                endif
            !
            end do
        end do
        !
        !     into the field
        if(myid.eq.nproc-1)then
            kparaendp=kparaend-1
        else
            kparaendp=kparaend
        endif
        do k=kparasta,kparaendp
            do j=1,jy
                do i=1,jx
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

    subroutine flucnesp(r,visualizzo,tipo,tau_wind)
        !***********************************************************************
        ! compute explicit diffusive term
        !
        ! NNI*G11*D(u,v,w)/D(csi)+
        ! NNI*G22*D(u,v,w)/D(eta)+
        ! NNI*G33*D(u,v,w)/D(zita)
        !
        use inflow_module, only: areola3,areola4

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer ierr,status(MPI_STATUS_SIZE)
        integer m
        integer kparastal,kparaendl
        double precision bulk_loc,bulkn
        integer i,j,k,ii,jj,kk
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        integer visualizzo!,windyes
        integer req1,req2
        integer istatus(MPI_STATUS_SIZE)

        ! Giulia modificavento: !SANTIAGO FIXED
        !      real vel_tau(0:n1+1,kparasta-1:kparaend+1)
        real,optional,intent(inout) :: tau_wind(0:,kparasta-1:)!(0:n1+1,kparasta-1:kparaend+1)

        real rhow
        ! Giulia modificavento:

        ! giulia aggiungo questi per eliminare i flussi tra i solidi
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! term nni*g11*d/d(csi)
        !
        do ii=1,ip
            !
            !     side left
            do k=kparasta,kparaend
                do j=1,jy
                    f1(0,j,k)=f1(0,j,k)+annit(0,j,k)*g11(0,j,k)* &
                        (-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
                end do
            end do
            !
            !     side right
            do k=kparasta,kparaend
                do j=1,jy
                    f1(jx,j,k)=f1(jx,j,k)+annit(jx+1,j,k)*g11(jx,j,k)* &
                        (8.*r(jx+1,j,k)-9.*r(jx,j,k)+r(jx-1,j,k))/3.
                end do
            end do
        !
        enddo
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    f1(i,j,k)=f1(i,j,k) &
                        +.5*(annit(i,j,k)+annit(i+1,j,k))*g11(i,j,k) &
                        *(r(i+1,j,k)-r(i,j,k))
                !
                enddo
            enddo
        enddo


        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !

                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
                        if(i.lt.jx)then
                            if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo

        !
        !-----------------------------------------------------------------------
        ! term nni*g22*d/d(eta)
        !
        do jj=1,jp
            !
            !       side bottom
            do k=kparasta,kparaend
                do i=1,jx

                    f2(i,0,k)=f2(i,0,k)+annitV(i,0,k)*g22(i,0,k)* &
                        (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                end do
            end do
            !  c
            !  c     side upper
            if(present(tau_wind))then

                !  c Giulia modificavento:
                rhow = 1000.

                do k=kparasta,kparaend
                    do i=1,jx

                        !    c      f2(i,jy,k) = abs(vel_tau(i,k))*vel_tau(i,k)*areola4(i,k)
                        f2(i,jy,k)=tau_wind(i,k)*areola4(i,k)/rhow

                    end do
                end do
            !  c Giulia modificavento:

            else

                do k=kparasta,kparaend
                    do i=1,jx

                        f2(i,jy,k)=f2(i,jy,k) &
                            +annitV(i,jy+1,k)*g22(i,jy,k)* &
                            (8.*r(i,jy+1,k)-9.*r(i,jy,k)+r(i,jy-1,k))/3.
                    end do
                end do

            end if
        !
        enddo
        !
        !     into the field
        do k=kparasta,kparaend
            do i=1,jx
                do j=jp,jy-jp
                    !
                    f2(i,j,k)=f2(i,j,k) &
                        +.5*(annitV(i,j,k)+annitV(i,j+1,k)) &
                        *g22(i,j,k)*(r(i,j+1,k)-r(i,j,k))
                !
                enddo
            enddo
        enddo

        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
                        if(j.lt.jy)then
                            if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! term nni*g33d/d(zita)
        !
        do kk=1,kp
            !
            !     side back
            if (myid.eq.0) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,0)=f3(i,j,0)+annit(i,j,0)*g33(i,j,0)* &
                            (-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
                    end do
                end do

            endif
            !
            !     side front
            if (myid.eq.nproc-1) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,jz)=f3(i,j,jz) &
                            +annit(i,j,jz+1)*g33(i,j,jz)* &
                            (8.*r(i,j,jz+1)-9.*r(i,j,jz)+r(i,j,jz-1))/3.
                    end do
                end do

            endif
        !
        enddo
        !
        !     into the field
        !
        if (myid.eq.0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid.eq.nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else
            kparastal=kparasta
            kparaendl=kparaend
        endif

        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    !
                    f3(i,j,k)=f3(i,j,k) &
                        +.5*(annit(i,j,k)+annit(i,j,k+1))*g33(i,j,k) &
                        *(r(i,j,k+1)-r(i,j,k))
                !
                enddo
            enddo
        enddo

        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
                        if(k.lt.jz)then
                            if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo

        !
        ! pass f3 at the border between procs

        if(myid.eq.0)then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if(myid.eq.nproc-1)then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else
            leftpem=leftpe
            rightpem=rightpe
        endif


        if(rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(f3(1,1,kparaend),jx*jy,MPI_REAL_SD, &
                rightpem ,tagrs,MPI_COMM_WORLD,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(f3(1,1,kparasta-1),jx*jy,MPI_REAL_SD, &
                leftpem  ,taglr,MPI_COMM_WORLD,status,ierr)
        endif

        if(rightpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req1,istatus,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
        !      call MPI_WAIT(req2,istatus,ierr)
        endif

        !     integral
        bulk_loc=0.
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    bulk_loc=bulk_loc+f3(i,j,k)-f3(i  ,j  ,k-1)+ &
                        f2(i,j,k)-f2(i  ,j-1,k  )+ &
                        f1(i,j,k)-f1(i-1,j  ,k  )
                end do
            end do
        end do

        ! make bulk known to all procs

        call MPI_ALLREDUCE(bulk_loc,bulkn,1,MPI_DOUBLE_PRECISION, &
            MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        bulk=bulk+bulkn

        return
    end subroutine flucnesp

    subroutine flucrho(ti,r,akaptV,akapt,isc,tipo)
        !***********************************************************************
        ! compute diffusive explicit terms
        !
        ! k*g11*D(rho)/D(csi) +
        ! k*g22*D(rho)/D(eta) +
        ! k*g33*D(rho)/D(zita)
        !
        ! the subroutine is for scalar eq.
        !
        use mysettings, only: coef_wall, insc, visualizzo
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer ierr,status(MPI_STATUS_SIZE)
        integer m
        integer kparastal,kparaendl
        integer i,j,k,ii,jj,kk,isc

        real ti
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real akaptV(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)

        real f2_loc,f2_tot,f2_mean

        ! giulia aggiungo questi per eliminare i flussi tra i solidi
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! term nni*g11*d/d(csi)
        !
        do ii=1,ip
            !
            !hicco usare un solo ciclo!!!!

            !     side left
            do k=kparasta,kparaend
                do j=1,jy
                    f1(0,j,k)=akapt(isc,0,j,k)*g11(0,j,k)* &
                        (-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
                end do
            end do
            !
            !     side right
            do k=kparasta,kparaend
                do j=1,jy
                    f1(jx,j,k)=akapt(isc,jx+1,j,k)*g11(jx,j,k)* &
                        (8.*r(jx+1,j,k)-9.*r(jx,j,k)+r(jx-1,j,k))/3.
                end do
            end do
        !
        end do
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    f1(i,j,k)=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)* &
                        (r(i+1,j,k)-r(i,j,k))
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !

                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
                        if(i.lt.jx)then
                            if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
                        endif
                    endif
                end do
            end do
        end do


        !
        !-----------------------------------------------------------------------
        ! term nni*g22*d/d(eta)
        !
        do jj=1,jp
            !
            !      if(coef_wall==1 .and. imoist==1)then
            !         call readflux(ti,isc)
            !      else

            !     side bottom
            do k=kparasta,kparaend
                do i=1,jx
                    f2(i,0,k)=akaptV(isc,i,0,k)*g22(i,0,k)* &
                        (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                end do
            end do

            !    end if
            !
            !     side upper
            do k=kparasta,kparaend
                do i=1,jx
                    f2(i,jy,k)=akaptV(isc,i,jy+1,k)*g22(i,jy,k)* &
                        (8.*r(i,jy+1,k)-9.*r(i,jy,k)+r(i,jy-1,k))/3.
                end do
            end do

        !hicco      end if
        !
        end do
        !
        !     into the field
        do k=kparasta,kparaend
            do i=1,jx
                do j=jp,jy-jp
                    !
                    f2(i,j,k)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)* &
                        (r(i,j+1,k)-r(i,j,k))
                !
                end do
            end do
        end do

        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
                        if(j.lt.jy)then
                            if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
                        endif
                    endif
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
            if (myid.eq.0) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,0)=akapt(isc,i,j,0)*g33(i,j,0)* &
                            (-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
                    end do
                end do
            !
            !     side front
            else if (myid.eq.nproc-1) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,jz)=akapt(isc,i,j,jz+1)*g33(i,j,jz)* &
                            (8.*r(i,j,jz+1)-9.*r(i,j,jz)+r(i,j,jz-1))/3.
                    end do
                end do

            endif
        !
        end do
        !
        !     into the field
        if (myid.eq.0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid.eq.nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else if ((myid.ne.0).and.(myid.ne.nproc-1)) then
            kparastal=kparasta
            kparaendl=kparaend
        endif


        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    !
                    f3(i,j,k)=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))*g33(i,j,k)* &
                        (r(i,j,k+1)-r(i,j,k))
                !
                end do
            end do
        end do

        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
                        if(k.lt.jz)then
                            if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
                        endif
                    endif
                end do
            end do
        end do
        !
        ! make border values f1, f2, f3 known to all procs

        call MPI_SENDRECV(f3(1,1,kparaend),jx*jy, &
            MPI_REAL_SD,rightpe,tagrs, &
            f3(1,1,kparasta-1),jx*jy, &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)

        !-----------------------------------------------------------------------
        !   check for the fluxes
        !    if(imoist==1)then
        !        !     Compute the mean fluxes at bottom
        !        f2_loc = 0.
        !        f2_tot = 0.
        !        f2_mean= 0.
        !        do k=kparasta,kparaend
        !            do i=1,jx
        !                f2_loc = f2_loc + f2(i,0,k)
        !            end do
        !        end do
        !        !     sum on f2_loc
        !        call MPI_ALLREDUCE(f2_loc,f2_tot,1,MPI_REAL_SD, &
        !            MPI_SUM,MPI_COMM_WORLD,ierr)
        !        f2_mean = f2_tot/real(jx*jz)
        !
        !        if(myid==0)then
        !            write(*,*)'scalar ',isc,' f2 mean bot: ',f2_mean
        !        end if
        !
        !        !     Compute the mean fluxes at top
        !        f2_loc = 0.
        !        f2_tot = 0.
        !        f2_mean= 0.
        !        do k=kparasta,kparaend
        !            do i=1,jx
        !                f2_loc = f2_loc + f2(i,jy,k)
        !            end do
        !        end do
        !        !     sum on f2_loc
        !        call MPI_ALLREDUCE(f2_loc,f2_tot,1,MPI_REAL_SD, &
        !            MPI_SUM,MPI_COMM_WORLD,ierr)
        !        f2_mean = f2_tot/real(jx*jz)
        !
        !        if(myid==0)then
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
        use myarrays_density
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer ierr,status(MPI_STATUS_SIZE)
        integer m
        integer kparastal,kparaendl
        integer i,j,k,ii,jj,kk,isc
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)

        ! giulia aggiungo questi per eliminare i flussi tra i solidi
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        ! term nni*g11*d/d(csi)
        !
        do ii=1,ip
            !hicco perche' usa due cicli invece di uno???? ridurre
            ! analogamente per sotto
            !
            ! side left

            do k=kparasta,kparaend
                do j=1,jy
                    f1(0,j,k)=f1(0,j,k)+akapt(isc,0,j,k)*g11(0,j,k)* &
                        (-8.*r(0,j,k)+9.*r(1,j,k)-r(2,j,k))/3.
                end do
            end do
            !
            ! side right
            !
            do k=kparasta,kparaend
                do j=1,jy
                    f1(jx,j,k)=f1(jx,j,k) &
                        +akapt(isc,jx+1,j,k)*g11(jx,j,k)* &
                        (8.*r(jx+1,j,k)-9.*r(jx,j,k)+r(jx-1,j,k))/3.
                end do
            end do
        !
        enddo
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    f1(i,j,k)=f1(i,j,k) &
                        +.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)* &
                        (r(i+1,j,k)-r(i,j,k))
                !
                enddo
            enddo
        enddo

        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
                        if(i.lt.jx)then
                            if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo

        !
        !-----------------------------------------------------------------------
        ! term nni*g22*d/d(eta)
        !
        do jj=1,jp
            !
            ! side bottom
            !
            do k=kparasta,kparaend
                do i=1,jx
                    f2(i,0,k)=f2(i,0,k)+akaptV(isc,i,0,k)*g22(i,0,k)* &
                        (-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
                end do
            end do
            !
            ! side upper
            !
            do k=kparasta,kparaend
                do i=1,jx
                    f2(i,jy,k)=f2(i,jy,k) &
                        +akaptV(isc,i,jy+1,k)*g22(i,jy,k)* &
                        (8.*r(i,jy+1,k)-9.*r(i,jy,k)+r(i,jy-1,k))/3.
                end do
            end do
        !
        enddo
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do i=1,jx
                do j=jp,jy-jp
                    !
                    f2(i,j,k)=f2(i,j,k) &
                        +.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)* &
                        (r(i,j+1,k)-r(i,j,k))
                !
                enddo
            enddo
        enddo

        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
                        if(j.lt.jy)then
                            if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! term nni*g33d/d(zita)
        !
        do kk=1,kp
            !
            ! side back
            !
            if (myid.eq.0) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,0)=f3(i,j,0)+akapt(isc,i,j,0)*g33(i,j,0)* &
                            (-8.*r(i,j,0)+9.*r(i,j,1)-r(i,j,2))/3.
                    end do
                end do

            ! side front
            !
            else if (myid.eq.nproc-1) then

                do i=1,jx
                    do j=1,jy
                        f3(i,j,jz)=f3(i,j,jz) &
                            +akapt(isc,i,j,jz+1)*g33(i,j,jz)* &
                            (8.*r(i,j,jz+1)-9.*r(i,j,jz)+r(i,j,jz-1))/3.
                    end do
                end do

            endif
        !
        enddo
        !
        ! into the field
        !
        if (myid.eq.0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid.eq.nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else if ((myid.ne.0).and.(myid.ne.nproc-1)) then
            kparastal=kparasta
            kparaendl=kparaend
        endif


        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    !
                    f3(i,j,k)=f3(i,j,k) &
                        +.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))*g33(i,j,k)* &
                        (r(i,j,k+1)-r(i,j,k))
                !
                enddo
            enddo
        enddo

        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
                        if(k.lt.jz)then
                            if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo


        !
        ! make border values f1, f2, f3 known to closer procs
        !
        call MPI_SENDRECV(f3(1,1,kparaend),jx*jy, &
            MPI_REAL_SD,rightpe,tagrs, &
            f3(1,1,kparasta-1),jx*jy, &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)
        !
        return
    end subroutine flucrhoesp

    subroutine flud1(rc,cgra1,r,akapt,insc,isc,tipo)
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
        integer ierr,insc,isc
        integer m
        integer kparastak,kparaendk
        integer i,j,k,is,iss,jj,kk,ii,iii,jjj,kkk
        !
        real ak,r0,r1,r2,fg
        real rc(0:n1,n2,kparasta:kparaend) !n3)
        real  r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real cgra1(0:n1,n2,kparasta-1:kparaend+1)
        real ravanti,rindietro1,rindietro2
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !      integer tipo(0:n1+1,0:n2+1,0:n3+1) c7

        real phic(n1,n2,kparasta:kparaend)
        real den, num, myv
        logical found_solid
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
                do j=1,jy
                    !
                    f1(0,j,k) = rc(0 ,j,k)*r(0   ,j,k)-cgra1(0 ,j,k)
                    f1(jx,j,k)= rc(jx,j,k)*r(jx+1,j,k)-cgra1(jx,j,k)
                !
                enddo
            enddo

        enddo
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    f1(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
                !
                end do
            end do
        end do
        !
        ! quick
        !
        if(insc.eq.1)then
            !     sides 1 and 2
            do k=kparasta,kparaend
                do j=1,jy

                    i=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(2,j,k)
                        rindietro1 =        r(1,j,k)
                        rindietro2 = ip*(2.*r(0,j,k) - r(1,j,k)) &
                            +(1.-ip)*r(0,j,k)

                        if(tipo(i,j,k)==1)then !ib cell
                            if(tipo(i+1,j,k)==0)then !solido davanti
                                f1(i,j,k)=0.
                            elseif(tipo(i-1,j,k)==0)then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-ravanti +2.*rindietro1 - rindietro2)
                            endif !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo

                    else
                        ravanti    =    r(1,j,k)
                        rindietro1 =    r(2,j,k)
                        rindietro2 =    r(3,j,k)

                        if(tipo(i+1,j,k)==1)then
                            if(tipo(i,j,k)==0)then !solido indietro
                                f1(i,j,k)=0.
                            elseif(tipo(i+2,j,k)==0)then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else !solido di lato
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*( -ravanti +2.*rindietro1 - rindietro2)
                            endif !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo
                    end if

                    !         if(tipo(i,j,k)==2)then
                    !    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    !         end if

                    i=jx-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =    r(jx  ,j,k)
                        rindietro1 =    r(jx-1,j,k)
                        rindietro2 =    r(jx-2,j,k)

                        if(tipo(i,j,k)==1)then !ib cell
                            if(tipo(i+1,j,k)==0)then !solido davanti
                                f1(i,j,k)=0.
                            elseif(tipo(i-1,j,k)==0)then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-ravanti +2.*rindietro1 - rindietro2)
                            endif !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo
                    else
                        ravanti    =        r(jx-1,j,k)
                        rindietro1 =        r(jx  ,j,k)
                        rindietro2 = ip*(2.*r(jx+1,j,k) - r(jx,j,k)) &
                            +(1.-ip)*r(jx+1,j,k)

                        if(tipo(i+1,j,k)==1)then
                            if(tipo(i,j,k)==0)then !solido indietro
                                f1(i,j,k)=0.
                            elseif(tipo(i+2,j,k)==0)then !uso diff centrate
                                f1(i,j,k)=f1(i,j,k)
                            else !solido di lato
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*( -ravanti +2.*rindietro1 - rindietro2)
                            endif !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                *.125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo
                    end if

                !         if(tipo(i,j,k)==2)then
                !   f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                !     >              *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                !         end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=1,jy
                    do i=2,jx-2
                        !Giulia non comprende tutti i casi, da un lato è quick dall'altro dc
                        !      if(tipo(i,j,k)==2)then
                        !     if (rc(i,j,k).gt.0.) then
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                        !      else
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                        !      end if
                        !         end if
                        !
                        if (rc(i,j,k).gt.0.) then

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i+1,j,k)==0)then !solido davanti
                                    f1(i,j,k)=0.
                                elseif(tipo(i-1,j,k)==0)then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                            endif !tipo

                        else !rc(i,j,k).le.0.

                            if(tipo(i+1,j,k)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f1(i,j,k)=0.
                                elseif(tipo(i+2,j,k)==0)then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else !solido di lato
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                            endif !tipo

                        endif !rc
                    end do
                end do
            end do


        ! QUICK MODIFICATO /SMART
        elseif(insc.eq.2)then

            ! GIULIA UCS - CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,jy
                        !
                        f1(0,j,k) = rc(0 ,j,k)*r(0   ,j,k)-cgra1(0 ,j,k)
                        f1(jx,j,k)= rc(jx,j,k)*r(jx+1,j,k)-cgra1(jx,j,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,jy
                    do i=ip,jx-ip
                        !
                        f1(i,j,k)=-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,jy
                    i=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(2,j,k)
                        rindietro1 =        r(1,j,k)
                        rindietro2 = ip*(2.*r(0,j,k) - r(1,j,k)) &
                            +(1.-ip)*r(0,j,k)
                    else
                        ravanti    =    r(1,j,k)
                        rindietro1 =    r(2,j,k)
                        rindietro2 =    r(3,j,k)
                    end if
                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)
                    !if (den.lt.0.0000000001.and.den.ge.0) den=0.0000000001
                    !if (den.gt.-0.0000000001.and.den.lt.0) den=-0.000000001
                    !        if(den.eq.0)then
                    !        phic(i,j,k)= 1.
                    !   else
                    !   phic(i,j,k)= (num /den)
                    !        end if
                    !   !if(phic(i,j,k).gt.1.or.phic(i,j,k).lt.0)phic(i,j,k)=-0.1
                    !
                    !    if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then !1/6<phic<4/5
                    !       f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                    !     >           *.125*( 3*ravanti +6*rindietro1 - rindietro2 )
                    !   else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then !0<phic<1/6
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                    !   else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then !4/5<phic<1
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                    !   else
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1
                    !   endif!phic
                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(3*ravanti+6*rindietro1-rindietro2)


                    i=jx-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =    r(jx   ,j,k)
                        rindietro1 =    r(jx-1,j, k)
                        rindietro2 =    r(jx-2,j, k)
                    else
                        ravanti    =        r(jx-1,j,k)
                        rindietro1 =        r(jx  ,j,k)
                        rindietro2 = ip*(2.*r(jx+1,j,k) - r(jx,j,k)) &
                            +(1.-ip)*r(jx+1,j,k)
                    end if

                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)
                    !   !if (den.lt.0.0000000001.and.den.ge.0) den=0.0000000001
                    !   !if (den.gt.-0.0000000001.and.den.lt.0) den=-0.000000001
                    !        if(den.eq.0)then
                    !        phic(i,j,k)= 1.
                    !   else
                    !   phic(i,j,k)= (num /den)
                    !        end if
                    !   !if(phic(i,j,k).gt.1.or.phic(i,j,k).lt.0)phic(i,j,k)=-0.1

                    !   if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then
                    !       f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
                    !     >           *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )
                    !   else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                    !   else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                    !   else
                    !        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1
                    !   endif!phic

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(3*ravanti+6*rindietro1-rindietro2)
                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                            else  !(rc(i,j,k).le.0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                                    .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                            endif
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i+1,j,k)==0)then !solido davanti

                                        f1(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    elseif(tipo(i-1,j,k)==0)then

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                            *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                                    !       if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                                    !              rindietro2 =r(i,j,k)

                                    !       endif
                                    !         endif

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                endif !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i+1,j,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f1(i,j,k)=0.
                                    elseif(tipo(i+2,j,k)==0)then !solido avanti
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                    else !solido di lato
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                                            .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))

                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib



                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if(tipo(i+1,j,k)==1)then
                                    !       if(tipo(i+2,j,k)==0)then
                                    !             rindietro2=r(i+1,j,k)
                                    !       endif
                                    !       endif

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                                        .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        !    endif !SMART /quick modificato upwind

        elseif(insc.eq.3)then !upwind

            ! GIULIA UCS - CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,jy
                        !
                        f1(0,j,k) = rc(0 ,j,k)*r(0   ,j,k)-cgra1(0 ,j,k)
                        f1(jx,j,k)= rc(jx,j,k)*r(jx+1,j,k)-cgra1(jx,j,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,jy
                    do i=ip,jx-ip
                        !
                        f1(i,j,k)=-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,jy
                    i=1
                    if(rc(i,j,k).gt.0.)then
                        rindietro1 =        r(1,j,k)
                    else
                        rindietro1 =    r(2,j,k)

                    end if
                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1



                    i=jx-1
                    if(rc(i,j,k).gt.0.)then
                        rindietro1 =    r(jx-1,j, k)
                    else
                        rindietro1 =        r(jx  ,j,k)
                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k).le.0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                            endif
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i+1,j,k)==0)then !solido davanti

                                        f1(i,j,k)=0.

                                    !       elseif(tipo(i-1,j,k)==0)then
                                    !
                                    !               f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione o dietro

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                                    !       if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                                    !              rindietro2 =r(i,j,k)

                                    !       endif
                                    !         endif

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                                endif !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i+1,j,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f1(i,j,k)=0.
                                        !        elseif(tipo(i+2,j,k)==0)then !solido avanti
                                        !              f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                        !               else !solido di lato
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib



                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if(tipo(i+1,j,k)==1)then
                                    !       if(tipo(i+2,j,k)==0)then
                                    !             rindietro2=r(i+1,j,k)
                                    !       endif
                                    !       endif

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        endif !quick





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
                do j=1+jp,jy-jp
                    !
                    fg=.5*(r(iss,j+1,k)-r(iss,j-1,k))
                    !             f1(is,j,k)=-f1(is,j,k)
                    !     >                   +  akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    !
                    if(tipo(iss,j-1,k).eq.0.or.tipo(iss,j+1,k).eq.0)then
                        f1(is,j,k)=-f1(is,j,k)
                    !         elseif(tipo(iss,j+1,k).eq.0)then
                    !           f1(is,j,k)=-f1(is,j,k)
                    else
                        f1(is,j,k)=-f1(is,j,k) &
                            +  akapt(isc,iss,j,k)*g12(is,j,k)*fg

                    end if

                end do

                !
                do jj=1,jp
                    !
                    !        check derivative at wall bottom and upper
                    !
                    j=1
                    fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j+1,k)-r(iss,j+2,k))
                    !         f1(is,j,k)=-f1(is,j,k)
                    !     >                +akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    !
                    if(tipo(iss,j+1,k).eq.0.or.tipo(iss,j+2,k).eq.0)then
                        f1(is,j,k)=-f1(is,j,k)
                    !         elseif(tipo(iss,j+2,k).eq.0)then
                    !           f1(is,j,k)=-f1(is,j,k)
                    else
                        f1(is,j,k)=-f1(is,j,k) &
                            +akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    end if


                    j=jy
                    fg=.5*(3.*r(iss,jy,k)-4.*r(iss,jy-1,k)+r(iss,jy-2,k))
                    !         f1(is,j,k)=-f1(is,j,k)
                    !     >                +akapt(isc,iss,j,k)*g12(is,j,k)*fg

                    if(tipo(iss,j-1,k).eq.0.or.tipo(iss,j-2,k).eq.0)then
                        f1(is,j,k)=-f1(is,j,k)
                    !         elseif(tipo(iss,j-2,k).eq.0)then
                    !           f1(is,j,k)=-f1(is,j,k)
                    else
                        f1(is,j,k)=-f1(is,j,k) &
                            +akapt(isc,iss,j,k)*g12(is,j,k)*fg
                    end if
                !
                end do
            !
            enddo
            !
            is=jx
            iss=jx+1
        !
        enddo
        !
        ! inside the field
        !
        do k=kparasta,kparaend
            do i=ip,jx-ip
                !
                do j=1+jp,jy-jp
                    !
                    r1=.5*(r(i,j-1,k)+r(i+1,j-1,k))
                    r2=.5*(r(i,j+1,k)+r(i+1,j+1,k))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))


                    found_solid=.false.
                    do kk=k,k
                        do jj=j-1,j+1
                            do ii=i,i+1
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f1(i,j,k)=-f1(i,j,k)
                    else
                        f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    end if

                !
                end do
                !
                !     check derivative on sides bottom and upper
                !
                do jjj=1,jp
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
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f1(i,j,k)=-f1(i,j,k)
                    else
                        f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    end if

                    !
                    j=jy
                    r0=.5*(r(i,jy  ,k)+r(i+1,jy  ,k))
                    r1=.5*(r(i,jy-1,k)+r(i+1,jy-1,k))
                    r2=.5*(r(i,jy-2,k)+r(i+1,jy-2,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
                    !      f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    !
                    found_solid=.false.
                    do kk=k,k
                        do jj=j-2,j
                            do ii=i,i+1
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f1(i,j,k)=-f1(i,j,k)
                    else
                        f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
                    end if
                end do
            !
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G13*Drho/D(ZITA)
        !
        ! sides left and right
        !
        ! define the limit dependeing on periodicity in z
        !
        if (myid.eq.0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid.eq.nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else if ((myid.ne.0).and.(myid.ne.nproc-1)) then
            kparastak=kparasta
            kparaendk=kparaend
        endif
        !
        is=0
        iss=0
        do i=1,2*ip
            !
            do j=1,jy
                !
                do k=kparastak,kparaendk
                    !
                    fg=.5*(r(iss,j,k+1)-r(iss,j,k-1))
                    !       f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
                    !

                    if(tipo(iss,j,k-1).eq.0.or.tipo(iss,j,k+1).eq.0)then
                        f1(is,j,k)=f1(is,j,k)
                    !         elseif(tipo(iss,j,k+1).eq.0)then
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
                    if (myid.eq.0) then

                        k=1
                        fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j,k+1)-r(iss,j,k+2))
                        !      f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg

                        if(tipo(iss,j,k+1).eq.0.or.tipo(iss,j,k+2).eq.0)then
                            f1(is,j,k)=f1(is,j,k)
                        !         elseif(tipo(iss,j,k+2).eq.0)then
                        !           f1(is,j,k)=f1(is,j,k)
                        else
                            f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
                        end if

                    endif
                    !
                    if (myid.eq.nproc-1) then
                        k=jz
                        fg=.5*(3.*r(iss,j,jz)-4.*r(iss,j,jz-1)+r(iss,j,jz-2))
                        !      f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg

                        if(tipo(iss,j,k-1).eq.0.or.tipo(iss,j,k-2).eq.0)then
                            f1(is,j,k)=f1(is,j,k)
                        !         elseif(tipo(iss,j,k-2).eq.0)then
                        !           f1(is,j,k)=f1(is,j,k)
                        else
                            f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
                        end if
                    endif
                !
                end do
            !
            enddo
            !
            is=jx
            iss=jx+1
        !
        enddo
        !
        ! into the field
        !

        do j=1,jy
            do i=ip,jx-ip
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
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f1(i,j,k)=f1(i,j,k)
                    else
                        f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
                    end if




                end do
                !
                do kkk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid.eq.0) then

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
                                    if(tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if(found_solid)then
                            f1(i,j,k)=f1(i,j,k)
                        else
                            f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
                        end if



                    endif
                    !
                    if (myid.eq.nproc-1) then

                        k=jz
                        r0=.5*(r(i,j,jz  )+r(i+1,j,jz  ))
                        r1=.5*(r(i,j,jz-1)+r(i+1,j,jz-1))
                        r2=.5*(r(i,j,jz-2)+r(i+1,j,jz-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
                        !      f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg

                        found_solid=.false.
                        do kk=k-2,k
                            do jj=j,j
                                do ii=i,i+1
                                    if(tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if(found_solid)then
                            f1(i,j,k)=f1(i,j,k)
                        else
                            f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
                        end if

                    endif
                !
                end do
            !
            enddo
        enddo


        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !

                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
                        if(i.lt.jx)then
                            if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo
        !
        return
    end subroutine flud1

    subroutine flud2(rc,cgra2,r,akapt,insc,isc,tipo)
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
        integer ierr,insc,isc
        integer m
        integer kparastak,kparaendk
        integer i,j,k,js,jss,ii,kk,jj,iii,jjj,kkk

        real ak,r0,r1,r2,fg
        real rc(n1,0:n2,kparasta:kparaend) !n3)
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real cgra2(n1,0:n2,kparasta-1:kparaend+1)
        real ravanti,rindietro1,rindietro2
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !      integer tipo(0:n1+1,0:n2+1,0:n3+1) c7
        !
        real phic(n1,n2,kparasta:kparaend)
        real den, num, myv
        logical found_solid
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Vc*rho:
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell

        !     sides bottom and upper
        do j=1,jp
            !
            do k=kparasta,kparaend
                do i=1,jx
                    !
                    f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
                    f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
                !
                enddo
            enddo

        enddo
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do i=1,jx
                do j=jp,jy-jp
                    !
                    f2(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j+1,k))-cgra2(i,j,k)
                !
                end do
            end do
        end do
        !
        ! quick
        !
        if(insc.eq.1)then

            !     sides 3 and 4
            do k=kparasta,kparaend
                do i=1,jx

                    j=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(i,2,k)
                        rindietro1 =        r(i,1,k)
                        rindietro2 = jp*(2.*r(i,0,k) - r(i,1,k)) &
                            +(1.-jp)*r(i,0,k)

                        if(tipo(i,j,k)==1)then !ib cell
                            if(tipo(i,j+1,k)==0)then !solido davanti
                                f2(i,j,k)=0.
                            elseif(tipo(i,j-1,k)==0)then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2)

                            endif !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo
                    else
                        ravanti    = r(i,1,k)
                        rindietro1 = r(i,2,k)
                        rindietro2 = r(i,3,k)

                        if(tipo(i,j+1,k)==1)then
                            if(tipo(i,j,k)==0)then !solido indietro
                                f2(i,j,k)=0.
                            elseif(tipo(i,j+2,k)==0)then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)

                            else !solido di lato
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2)


                            endif !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo
                    end if

                    !         if(tipo(i,j,k)==2)then
                    !    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    !         end if

                    j=jy-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    = r(i,jy  ,k)
                        rindietro1 = r(i,jy-1,k)
                        rindietro2 = r(i,jy-2,k)

                        if(tipo(i,j,k)==1)then !ib cell
                            if(tipo(i,j+1,k)==0)then !solido davanti
                                f2(i,j,k)=0.
                            elseif(tipo(i,j-1,k)==0)then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)
                            else ! solido sta nell'altra direzione quick
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2)

                            endif !solido avanti/indietro/lato
                        else !i è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo
                    else
                        ravanti    =        r(i,jy-1,k)
                        rindietro1 =        r(i,jy  ,k)
                        rindietro2 = jp*(2.*r(i,jy+1,k) - r(i,jy,k)) &
                            +(1.-jp)*r(i,jy+1,k)

                        if(tipo(i,j+1,k)==1)then
                            if(tipo(i,j,k)==0)then !solido indietro
                                f2(i,j,k)=0.
                            elseif(tipo(i,j+2,k)==0)then !uso diff centrate
                                f2(i,j,k)=f2(i,j,k)

                            else !solido di lato
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2)


                            endif !solido avanti/indietro/lato
                        else !i+1 è fluido allora quick normale
                            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                .125*(-ravanti +2.*rindietro1 - rindietro2)
                        endif !tipo
                    end if
                !         if(tipo(i,j,k)==2)then
                !   f2(i,j,k)=f2(i,j,k)+rc(i,j,k)
                !     >              *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                !         end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=2,jy-2
                    do i=1,jx
                        !
                        !         if(tipo(i,j,k)==2)then
                        !      if (rc(i,j,k).gt.0.) then
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                        !      else
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                        !      end if
                        !      end if
                        !

                        if (rc(i,j,k).gt.0.) then

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i,j+1,k)==0)then !solido davanti
                                    f2(i,j,k)=0.
                                elseif(tipo(i,j-1,k)==0)then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                            endif !tipo

                        else !rc(i,j,k).le.0.
                            if(tipo(i,j+1,k)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f2(i,j,k)=0.
                                elseif(tipo(i,j+2,k)==0)then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)

                                else !solido di lato
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))


                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                            endif !tipo

                        endif !rc
                    !
                    end do
                end do
            end do

        ! smart/quick modificato
        elseif(insc.eq.2)then

            ! GIULIA UCS - CGRA
            !     side left and right
            do j=1,jp
                !
                do k=kparasta,kparaend
                    do i=1,jx
                        !
                        f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
                        f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=jp,jy-jp
                    do i=1,jx
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,jx
                    j=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(i,2,k)
                        rindietro1 =        r(i,1,k)
                        rindietro2 = jp*(2.*r(i,0,k) - r(i,1,k)) &
                            +(1.-jp)*r(i,0,k)
                    else
                        ravanti    =    r(i,1,k)
                        rindietro1 =    r(i,2,k)
                        rindietro2 =    r(i,3,k)
                    end if
                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)
                    !   !if (den.lt.0.0000000001.and.den.ge.0) den=0.0000000001
                    !   !if (den.gt.-0.0000000001.and.den.lt.0) den=-0.000000001
                    !        if(den.eq.0)then
                    !        phic(i,j,k)= 1.
                    !   else
                    !   phic(i,j,k)= (num /den)
                    !        end if
                    !   !if(phic(i,j,k).gt.1.or.phic(i,j,k).lt.0)phic(i,j,k)=-0.1

                    !    if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then !1/6<phic<4/5
                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )
                    !   else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then !0<phic<1/6
                    !        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                    !   else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then !4/5<phic<1
                    !        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                    !   else
                    !        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1
                    !   endif!phic

                    j=jy-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =    r(i,jy  , k)
                        rindietro1 =    r(i,jy-1, k)
                        rindietro2 =    r(i,jy-2, k)
                    else
                        ravanti    =        r(i,jy-1,k)
                        rindietro1 =        r(i,jy  ,k)
                        rindietro2 = jp*(2.*r(i,jy+1,k) - r(i,jy,k)) &
                            +(1.-jp)*r(i,jy+1,k)
                    end if

                    !   den=(ravanti-rindietro2)
                    !   num=(rindietro1-rindietro2)

                    !        !if (den.lt.0.0000000001.and.den.ge.0) den=0.0000000001
                    !   !if (den.gt.-0.0000000001.and.den.lt.0) den=-0.000000001
                    !        if(den.eq.0)then
                    !        phic(i,j,k)= 1.
                    !   else
                    !   phic(i,j,k)= (num /den)
                    !        end if
                    !   !if(phic(i,j,k).gt.1.or.phic(i,j,k).lt.0)phic(i,j,k)=-0.1

                    !   if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then
                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )
                !  else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then
                !       f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
                !  else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then
                !       f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
                !  else
                !       f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1
                !  endif!phic
                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                            else  !(rc(i,j,k).le.0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k)) !my index +1

                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j+1,k)==0)then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo
                                    elseif(tipo(i,j-1,k)==0)then

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione


                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                                endif   !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j+1,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f2(i,j,k)=0.
                                    elseif(tipo(i,j+2,k)==0)then !solido avanti
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                    else !solido di lato

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))
                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        !    endif !SMART/QUICK modificato upwind

        ! smart/quick modificato
        elseif(insc.eq.3)then

            ! GIULIA UCS - CGRA
            !     side left and right
            do j=1,jp
                !
                do k=kparasta,kparaend
                    do i=1,jx
                        !
                        f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
                        f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=jp,jy-jp
                    do i=1,jx
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,jx
                    j=1
                    if(rc(i,j,k).gt.0.)then

                        rindietro1 =        r(i,1,k)

                    else

                        rindietro1 =    r(i,2,k)

                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1


                    j=jy-1
                    if(rc(i,j,k).gt.0.)then

                        rindietro1 =    r(i,jy-1, k)

                    else

                        rindietro1 =        r(i,jy  ,k)

                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)


                            else  !(rc(i,j,k).le.0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j+1,k)==0)then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo
                                    !       elseif(tipo(i,j-1,k)==0)then

                                    !               f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione


                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                endif   !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j+1,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f2(i,j,k)=0.
                                    else
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        endif !quick
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G21*Drho/D(CSI)
        !
        !     sides bottom and upper
        js=0
        jss=0
        do j=1,2*jp
            !
            do k=kparasta,kparaend
                !
                do  i=1+ip,jx-ip
                    !
                    fg=.5*(r(i+1,jss,k)-r(i-1,jss,k))
                    !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
                    !

                    if(tipo(i-1,jss,k).eq.0.or.tipo(i+1,jss,k).eq.0)then
                        f2(i,js,k)=-f2(i,js,k)
                    !         elseif(tipo(i+1,jss,k).eq.0)then
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

                    if(tipo(i+1,jss,k).eq.0.or.tipo(i+2,jss,k).eq.0)then
                        f2(i,js,k)=-f2(i,js,k)
                    !         elseif(tipo(i+2,jss,k).eq.0)then
                    !           f2(i,js,k)=-f2(i,js,k)
                    else
                        f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
                    end if
                    !
                    i=jx
                    fg=.5*(3.*r(jx,jss,k)-4.*r(jx-1,jss,k)+r(jx-2,jss,k))
                    !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg

                    if(tipo(i-1,jss,k).eq.0.or.tipo(i-2,jss,k).eq.0)then
                        f2(i,js,k)=-f2(i,js,k)
                    !         elseif(tipo(i-2,jss,k).eq.0)then
                    !           f2(i,js,k)=-f2(i,js,k)
                    else
                        f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
                    end if
                !
                end do
            !
            enddo
            !
            js=jy
            jss=jy+1
        !
        enddo
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=jp,jy-jp
                !
                do i=1+ip,jx-ip
                    !
                    r1=.5*(r(i-1,j,k)+r(i-1,j+1,k))
                    r2=.5*(r(i+1,j,k)+r(i+1,j+1,k))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))


                    found_solid=.false.
                    do kk=k,k
                        do jj=j,j+1
                            do ii=i-1,i+1
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f2(i,j,k)=-f2(i,j,k)
                    else
                        f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    end if


                end do
                !
                do iii=1,ip
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
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f2(i,j,k)=-f2(i,j,k)
                    else
                        f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    end if
                    !
                    i=jx
                    r0=.5*(r(jx  ,j,k)+r(jx  ,j+1,k))
                    r1=.5*(r(jx-1,j,k)+r(jx-1,j+1,k))
                    r2=.5*(r(jx-2,j,k)+r(jx-2,j+1,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
                    !      f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    !
                    found_solid=.false.
                    do kk=k,k
                        do jj=j,j+1
                            do ii=i-2,i
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f2(i,j,k)=-f2(i,j,k)
                    else
                        f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
                    end if

                end do
            !
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G23*Drho/D(ZITA)
        !
        ! sides bottom and upper
        !
        if (myid.eq.0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid.eq.nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else if ((myid.gt.0).and.(myid.lt.nproc-1)) then
            kparastak=kparasta
            kparaendk=kparaend
        endif

        js=0
        jss=0
        do j=1,2*jp
            !
            do i=1,jx
                !
                do k=kparastak,kparaendk
                    !
                    fg=.5*(r(i,jss,k+1)-r(i,jss,k-1))
                    !         f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

                    if(tipo(i,jss,k-1).eq.0.or.tipo(i+1,jss,k+1).eq.0)then
                        f2(i,js,k)=f2(i,js,k)
                    !         elseif(tipo(i+1,jss,k+1).eq.0)then
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
                    if (myid.eq.0) then

                        k=1
                        fg=.5*(-3.*r(i,jss,k)+4.*r(i,jss,k+1)-r(i,jss,k+2))
                        !      f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

                        if(tipo(i,jss,k+1).eq.0.or.tipo(i+1,jss,k+2).eq.0)then
                            f2(i,js,k)=f2(i,js,k)
                        !         elseif(tipo(i+1,jss,k+2).eq.0)then
                        !           f2(i,js,k)=f2(i,js,k)
                        else
                            f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
                        end if

                    endif
                    !
                    if (myid.eq.nproc-1) then

                        k=jz
                        fg=.5*(3.*r(i,jss,jz)-4.*r(i,jss,jz-1)+r(i,jss,jz-2))
                        !      f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

                        if(tipo(i,jss,k-1).eq.0.or.tipo(i+1,jss,k-2).eq.0)then
                            f2(i,js,k)=f2(i,js,k)
                        !         elseif(tipo(i+1,jss,k-2).eq.0)then
                        !           f2(i,js,k)=f2(i,js,k)
                        else
                            f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
                        end if

                    endif
                !
                end do
            !
            enddo
            !
            js=jy
            jss=jy+1
        !
        enddo
        !
        ! into the field
        !
        do j=jp,jy-jp
            do i=1,jx
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
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f2(i,j,k)=f2(i,j,k)
                    else
                        f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
                    end if



                end do
                !
                do kkk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid.eq.0) then

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
                                    if(tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if(found_solid)then
                            f2(i,j,k)=f2(i,j,k)
                        else
                            f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
                        end if


                    endif
                    !
                    if (myid.eq.nproc-1) then

                        k=jz
                        r0=.5*(r(i,j,jz  )+r(i,j+1,jz  ))
                        r1=.5*(r(i,j,jz-1)+r(i,j+1,jz-1))
                        r2=.5*(r(i,j,jz-2)+r(i,j+1,jz-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
                        !      f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg

                        found_solid=.false.
                        do kk=k-2,k
                            do jj=j,j+1
                                do ii=i,i
                                    if(tipo(ii,jj,kk)==0)found_solid=.true.
                                end do
                            end do
                        end do

                        if(found_solid)then
                            f2(i,j,k)=f2(i,j,k)
                        else
                            f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
                        end if
                    endif
                !
                end do
            !
            enddo
        enddo
        !


        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx
                    !

                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
                        if(j.lt.jy)then
                            if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo

        return
    end subroutine flud2

    subroutine flud3(rc,cgra3,r,akapt,insc,isc,tipo)
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
        integer ierr,status(MPI_STATUS_SIZE)
        integer m,insc,isc
        integer kparastal,kparaendl
        integer kparastall,kparaendll
        integer i,j,k,ks,kss,ii,jj,kk,iii,jjj,kkk
        !
        real ak,r0,r1,r2,fg
        real    rc(n1,n2,kparasta-1:kparaend) !0:n3)
        real    r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real    akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real cgra3(n1,n2,kparasta-1:kparaend+1)
        real ravanti,rindietro1,rindietro2
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !      integer tipo(0:n1+1,0:n2+1,0:n3+1) c7

        real phic(n1,n2,kparasta:kparaend)
        real den, num, myv

        integer req1,req2
        integer istatus(MPI_STATUS_SIZE)
        logical found_solid
        !
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Wc*rho:
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell
        !
        ! sides back and front
        do k=1,kp

            do j=1,jy
                do i=1,jx
                    !
                    if (myid.eq.0) then
                        f3(i,j,0)  = rc(i,j,0 )*r(i,j,0   )-cgra3(i,j,0 )
                    endif

                    if (myid.eq.nproc-1) then
                        f3(i,j,jz) = rc(i,j,jz)*r(i,j,jz+1)-cgra3(i,j,jz)
                    endif
                !
                enddo
            enddo

        enddo
        !
        ! into the field
        !
        if (myid.eq.0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid.eq.nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else
            kparastal=kparasta
            kparaendl=kparaend
        endif

        do k=kparastal,kparaendl
            do j=1,jy
                do i=1,jx
                    !
                    f3(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
                !
                end do
            end do
        end do

        ! quick
        !
        if(insc.eq.1)then

            if (myid.eq.0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid.eq.nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            endif

            !     sides 5 and 6
            if(myid .eq. 0)then
                do j=1,jy
                    do i=1,jx

                        k=1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    =        r(i,j,2)
                            rindietro1 =        r(i,j,1)
                            rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i,j,k+1)==0)then !solido davanti
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k-1)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1 - rindietro2 )

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2 )
                            endif !tipo
                        else
                            ravanti    =    r(i,j,1)
                            rindietro1 =    r(i,j,2)
                            rindietro2 =    r(i,j,3)

                            if(tipo(i,j,k+1)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k+2)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1 - rindietro2)
                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2)
                            endif !tipo
                        end if

                    !         if(tipo(i,j,k)==2)then
                    !    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    !         end if
                    end do
                end do
            end if

            if(myid .eq. nproc-1)then
                do j=1,jy
                    do i=1,jx

                        k=jz-1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    = r(i,j,jz  )
                            rindietro1 = r(i,j,jz-1)
                            rindietro2 = r(i,j,jz-2)

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i,j,k+1)==0)then !solido davanti
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k-1)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1 - rindietro2 )

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2 )
                            endif !tipo
                        else
                            ravanti    =        r(i,j,jz-1)
                            rindietro1 =        r(i,j,jz  )
                            rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                                +(1.-kp)*r(i,j,jz+1)


                            if(tipo(i,j,k+1)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k+2)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-ravanti +2.*rindietro1 - rindietro2)
                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-ravanti +2.*rindietro1 - rindietro2)
                            endif !tipo

                        end if

                    !         if(tipo(i,j,k)==2)then
                    !    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)
                    !     >           *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    !         end if
                    end do
                end do
            end if

            !     into the field

            do k=kparastall,kparaendll
                do j=1,jy
                    do i=1,jx
                        !
                        !      if(tipo(i,j,k)==2)then
                        !      if (rc(i,j,k).gt.0.) then
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                        !      else
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                        !      end if
                        !      end if

                        if (rc(i,j,k).gt.0.) then

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i,j,k+1)==0)then !solido davanti
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k-1)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                            endif !tipo

                        else !rc(i,j,k).le.0.

                            if(tipo(i,j,k+1)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k+2)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                            endif !tipo

                        endif !rc
                    !
                    end do
                end do
            end do


        ! quick modificato
        elseif(insc.eq.2)then

            ! GIULIA UCS - CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,jy
                    do i=1,jx
                        !
                        if (myid.eq.0) then
                            f3(i,j,0)  = rc(i,j,0 )*r(i,j,0   )-cgra3(i,j,0 )
                        endif

                        if (myid.eq.nproc-1) then
                            f3(i,j,jz) = rc(i,j,jz)*r(i,j,jz+1)-cgra3(i,j,jz)
                        endif
                    !
                    enddo
                enddo
            !
            enddo

            if (myid.eq.0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid.eq.nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            endif


            do k=kparastal,kparaendl
                do j=1,jy
                    do i=1,jx
                        !
                        f3(i,j,k)=-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid.eq.0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid.eq.nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            endif

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if(myid .eq. 0)then
                do j=1,jy
                    do i=1,jx
                        k=1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    =        r(i,j,2)
                            rindietro1 =        r(i,j,1)
                            rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti    =    r(i,j,1)
                            rindietro1 =    r(i,j,2)
                            rindietro2 =    r(i,j,3)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*(3*ravanti+6*rindietro1-rindietro2)

                    end do
                end do

            end if !myid==0

            if(myid .eq. nproc-1)then
                do j=1,jy
                    do i=1,jx

                        k=jz-1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    = r(i,j,jz  )
                            rindietro1 = r(i,j,jz-1)
                            rindietro2 = r(i,j,jz-2)
                        else
                            ravanti    =        r(i,j,jz-1)
                            rindietro1 =        r(i,j,jz  )
                            rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                                +(1.-kp)*r(i,j,jz+1)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*(3*ravanti+6*rindietro1-rindietro2)
                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                    *.125*(3*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                            else  !(rc(i,j,k).le.0.)
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(3*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j,k+1)==0)then !solido davanti

                                        f3(i,j,k)=0.

                                    elseif(tipo(i,j,k-1)==0)then !solido dietro

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else !solido lati
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                            *.125*(3*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))
                                    endif !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k-1)==1)then !la cella prima è un ib
                                    !       if(tipo(i,j,k-2)==0)then !perché considero anche questa?
                                    !         rindietro2 =r(i,j,k)
                                    !       endif
                                    !         endif !fine k-1 IB

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                        *.125*(3*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))


                                endif !tip02

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j,k+1)==1)then !ib
                                    if(tipo(i,j,k)==0)then !solido avanti
                                        f3(i,j,k)=0.
                                    elseif(tipo(i,j,k+2)==0)then !solido indietro
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)
                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                            .125*(3*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

                                    endif !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k+1)==1)then
                                    !       if(tipo(i,j,k+2)==0)then
                                    !       rindietro2=r(i,j,k+1)
                                    !       endif
                                    !         endif


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(3*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))
                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        !    endif !SMART
        !-----SMART

        ! quick modificato
        elseif(insc.eq.3)then

            ! GIULIA UCS - CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,jy
                    do i=1,jx
                        !
                        if (myid.eq.0) then
                            f3(i,j,0)  = rc(i,j,0 )*r(i,j,0   )-cgra3(i,j,0 )
                        endif

                        if (myid.eq.nproc-1) then
                            f3(i,j,jz) = rc(i,j,jz)*r(i,j,jz+1)-cgra3(i,j,jz)
                        endif
                    !
                    enddo
                enddo
            !
            enddo

            if (myid.eq.0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid.eq.nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            endif


            do k=kparastal,kparaendl
                do j=1,jy
                    do i=1,jx
                        !
                        f3(i,j,k)=-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid.eq.0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid.eq.nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            endif

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if(myid .eq. 0)then
                do j=1,jy
                    do i=1,jx
                        k=1
                        if(rc(i,j,k).gt.0.)then

                            rindietro1 =        r(i,j,1)

                        else

                            rindietro1 =    r(i,j,2)

                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1

                    end do
                end do

            end if !myid==0

            if(myid .eq. nproc-1)then
                do j=1,jy
                    do i=1,jx

                        k=jz-1
                        if(rc(i,j,k).gt.0.)then

                            rindietro1 = r(i,j,jz-1)

                        else

                            rindietro1 =        r(i,j,jz  )

                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1
                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k).le.0.)
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j,k+1)==0)then !solido davanti

                                        f3(i,j,k)=0.


                                    else !solido lati
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    endif !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k-1)==1)then !la cella prima è un ib
                                    !       if(tipo(i,j,k-2)==0)then !perché considero anche questa?
                                    !         rindietro2 =r(i,j,k)
                                    !       endif
                                    !         endif !fine k-1 IB

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)


                                endif !tip02

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j,k+1)==1)then !ib
                                    if(tipo(i,j,k)==0)then !solido avanti
                                        f3(i,j,k)=0.

                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

                                    endif !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib


                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k+1)==1)then
                                    !       if(tipo(i,j,k+2)==0)then
                                    !       rindietro2=r(i,j,k+1)
                                    !       endif
                                    !         endif


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)
                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        endif !SMART
        !-----SMART


        !----------------------------------------------------------------------
        ! DIFFUSIVE TERM  NNI*G31*Drho/D(CSI)
        !
        !     sides back and front
        ks=0
        kss=0
        do k=1,2*kp

            if (myid.eq.(k-1)*(nproc-1)) then

                do j=1,jy
                    !
                    do i=1+ip,jx-ip
                        !
                        fg=.5*(r(i+1,j,kss)-r(i-1,j,kss))
                        !         f3(i,j,ks)=-f3(i,j,ks)
                        !     >                 +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        !

                        if(tipo(i-1,j,kss).eq.0)then
                            f3(i,j,ks)=-f3(i,j,ks)
                        elseif(tipo(i+1,j,kss).eq.0)then
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

                        if(tipo(i+1,j,kss).eq.0.or.tipo(i+2,j,kss).eq.0)then
                            f3(i,j,ks)=-f3(i,j,ks)
                        !         elseif(tipo(i+2,j,kss).eq.0)then
                        !           f3(i,j,ks)=-f3(i,j,ks)
                        else
                            f3(i,j,ks)=-f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        end if
                        !
                        i=jx
                        fg=.5*(3.*r(jx,j,kss)-4.*r(jx-1,j,kss)+r(jx-2,j,kss))
                        !      f3(i,j,ks)=-f3(i,j,ks)
                        !     >               +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        !
                        if(tipo(i-1,j,kss).eq.0.or.tipo(i-2,j,kss).eq.0)then
                            f3(i,j,ks)=-f3(i,j,ks)
                        !         elseif(tipo(i-2,j,kss).eq.0)then
                        !           f3(i,j,ks)=-f3(i,j,ks)
                        else
                            f3(i,j,ks)=-f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g31(i,j,ks)*fg
                        end if
                    end do
                !
                enddo
            !
            endif

            ks=jz
            kss=jz+1
        !
        enddo
        !
        ! into the field
        !
        do k=kparastal,kparaendl
            do j=1,jy
                !
                do i=1+ip,jx-ip
                    !
                    r1=.5*(r(i-1,j,k)+r(i-1,j,k+1))
                    r2=.5*(r(i+1,j,k)+r(i+1,j,k+1))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))


                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j,j
                            do ii=i-1,i+1
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f3(i,j,k)=-f3(i,j,k)
                    else
                        f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg
                    end if



                !
                end do
                !
                do iii=1,ip
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
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f3(i,j,k)=-f3(i,j,k)
                    else
                        f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg
                    end if
                    !
                    i=jx
                    r0=.5*(r(jx  ,j,k)+r(jx  ,j,k+1))
                    r1=.5*(r(jx-1,j,k)+r(jx-1,j,k+1))
                    r2=.5*(r(jx-2,j,k)+r(jx-2,j,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))
                    !      f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg

                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j,j
                            do ii=i-2,i
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f3(i,j,k)=-f3(i,j,k)
                    else
                        f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*ak*fg
                    end if
                !
                end do
            !
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G32*Drho/D(eta)   !ZITA)
        !
        !     sides back and front
        ks=0
        kss=0
        do k=1,2*kp

            if (myid.eq.(k-1)*(nproc-1)) then

                do i=1,jx
                    !
                    do j=1+jp,jy-jp
                        !
                        fg=.5*(r(i,j+1,kss)-r(i,j-1,kss))
                        !         f3(i,j,ks)=f3(i,j,ks)
                        !     >                 +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        !
                        if(tipo(i,j-1,kss).eq.0.or.tipo(i,j+1,kss).eq.0)then
                            f3(i,j,ks)=f3(i,j,ks)
                        !         elseif(tipo(i,j+1,kss).eq.0)then
                        !           f3(i,j,ks)=f3(i,j,ks)
                        else
                            f3(i,j,ks)=f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        end if
                    end do
                    !
                    do jj=1,jp
                        !
                        ! check derivative on sides back and front
                        !
                        j=1
                        fg=.5*(-3.*r(i,j,kss)+4.*r(i,j+1,kss)-r(i,j+2,kss))
                        !      f3(i,j,ks)=f3(i,j,ks)
                        !     >              +akapt(isc,i,j,kss)*g32(i,j,ks)*fg

                        if(tipo(i,j+1,kss).eq.0.or.tipo(i,j+2,kss).eq.0)then
                            f3(i,j,ks)=f3(i,j,ks)
                        !         elseif(tipo(i,j+2,kss).eq.0)then
                        !           f3(i,j,ks)=f3(i,j,ks)
                        else
                            f3(i,j,ks)=f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        end if
                        !
                        j=jy
                        fg=.5*(3.*r(i,jy,kss)-4.*r(i,jy-1,kss)+r(i,jy-2,kss))
                        !      f3(i,j,ks)=f3(i,j,ks)
                        !     >              +akapt(isc,i,j,kss)*g32(i,j,ks)*fg

                        if(tipo(i,j-1,kss).eq.0.or.tipo(i,j-2,kss).eq.0)then
                            f3(i,j,ks)=f3(i,j,ks)
                        !         elseif(tipo(i,j-2,kss).eq.0)then
                        !           f3(i,j,ks)=f3(i,j,ks)
                        else
                            f3(i,j,ks)=f3(i,j,ks) &
                                +akapt(isc,i,j,kss)*g32(i,j,ks)*fg
                        end if
                    !
                    end do
                !
                enddo
            !
            endif

            ks=jz
            kss=jz+1
        !
        enddo
        !
        ! into the field
        !

        do k=kparastal,kparaendl
            do i=1,jx
                !
                do j=1+jp,jy-jp
                    !
                    r1=.5*(r(i,j-1,k)+r(i,j-1,k+1))
                    r2=.5*(r(i,j+1,k)+r(i,j+1,k+1))
                    fg=.5*(r2-r1)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))


                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j-1,j+1
                            do ii=i,i
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f3(i,j,k)=f3(i,j,k)
                    else
                        f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    end if



                !
                end do
                !
                do jjj=1,jp
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
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f3(i,j,k)=f3(i,j,k)
                    else
                        f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    end if

                    !
                    j=jy
                    r0=.5*(r(i,jy  ,k)+r(i,jy  ,k+1))
                    r1=.5*(r(i,jy-1,k)+r(i,jy-1,k+1))
                    r2=.5*(r(i,jy-2,k)+r(i,jy-2,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j,k+1))
                    !      f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    !

                    found_solid=.false.
                    do kk=k,k+1
                        do jj=j-2,j
                            do ii=i,i
                                if(tipo(ii,jj,kk)==0)found_solid=.true.
                            end do
                        end do
                    end do

                    if(found_solid)then
                        f3(i,j,k)=f3(i,j,k)
                    else
                        f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*ak*fg
                    end if

                end do
            !
            enddo
        enddo


        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    !
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
                        if(k.lt.jz)then
                            if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo


        !     pass f3 at the border between procs
        if(myid.eq.0)then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if(myid.eq.nproc-1)then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else
            leftpem=leftpe
            rightpem=rightpe
        endif


        if(rightpem /= MPI_PROC_NULL) then
            call MPI_ISEND(f3(1,1,kparaend),jx*jy,MPI_REAL_SD, &
                rightpem ,tagrs,MPI_COMM_WORLD,req1,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_IRECV(f3(1,1,kparasta-1),jx*jy,MPI_REAL_SD, &
                leftpem  ,taglr,MPI_COMM_WORLD,req2,ierr)
        endif

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_WAIT(req1,istatus,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_WAIT(req2,istatus,ierr)
        endif

        return
    end subroutine flud3

    subroutine flux1(rc,cgra1,r,insc,tipo)
        !***********************************************************************
        ! compute explict flux on csi component:
        !
        ! convective  uc*(u,v,w)  with cenetered scheme or quick
        ! diffusive   nni*g12*d(u,v,w)/d(eta)
        ! diffusive   nni*g13*d(u,v,w)/d(zita)
        !
        use myarrays_LC

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer ierr,insc
        integer m
        integer kparastak,kparaendk
        double precision bulk_loc
        !
        integer i,j,k,jj,kk,is,iss
        real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddi,dsi
        real fg,r0,r1,r2,an
        real rc(0:n1,1:n2,kparasta  :kparaend)
        real cgra1(0:n1,n2,kparasta-1:kparaend+1)
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        real ravanti,rindietro1,rindietro2
        !      integer tipo2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
        !      integer tipo2(0:n1+1,0:n2+1,0:n3+1) c7

        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        real phic(n1,n2,kparasta:kparaend)
        real den, num, myv
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
                do j=1,jy
                    !
                    f1(0,j,k) = (rc(0 ,j,k)+ucs(0,j,k))*r(0   ,j,k)- &
                        cgra1(0 ,j,k)
                    f1(jx,j,k)= (rc(jx,j,k)+ucs(jx,j,k))*r(jx+1,j,k)- &
                        cgra1(jx,j,k)
                !
                enddo
            enddo

        enddo
        !
        !     into the field
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    f1(i,j,k)=(rc(i,j,k)+ucs(i,j,k)) &
                        *.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
                !
                end do
            end do
        end do
        !
        !
        !
        ! quick
        !
        if(insc.eq.1)then
            !     sides 1 and 2
            do k=kparasta,kparaend
                do j=1,jy

                    i=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(2,j,k)
                        rindietro1 =        r(1,j,k)
                        rindietro2 = ip*(2.*r(0,j,k) - r(1,j,k)) &
                            +(1.-ip)*r(0,j,k)
                    else
                        ravanti    =    r(1,j,k)
                        rindietro1 =    r(2,j,k)
                        rindietro2 =    r(3,j,k)
                    end if

                    if(tipo(i,j,k)==2)then
                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    end if

                    i=jx-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    = r(jx  ,j,k)
                        rindietro1 = r(jx-1,j,k)
                        rindietro2 = r(jx-2,j,k)
                    else
                        ravanti    =        r(jx-1,j,k)
                        rindietro1 =        r(jx  ,j,k)
                        rindietro2 = ip*(2.*r(jx+1,j,k) - r(jx,j,k)) &
                            +(1.-ip)*r(jx+1,j,k)
                    end if

                    if(tipo(i,j,k)==2)then
                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=1,jy
                    do i=2,jx-2
                        !Giulia non comprende tutti i casi, da un lato è quick dall'altro dc
                        !      if(tipo(i,j,k)==2)then
                        !     if (rc(i,j,k).gt.0.) then
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                        !      else
                        !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                        !      end if
                        !         end if
                        !
                        if (rc(i,j,k).gt.0.) then

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i+1,j,k)==0)then !solido davanti
                                    f1(i,j,k)=0.
                                elseif(tipo(i-1,j,k)==0)then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                            endif !tipo

                        else !rc(i,j,k).le.0.

                            if(tipo(i+1,j,k)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f1(i,j,k)=0.
                                elseif(tipo(i+2,j,k)==0)then !uso diff centrate
                                    f1(i,j,k)=f1(i,j,k)
                                else !solido di lato
                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                            endif !tipo

                        endif !rc

                    end do
                end do
            end do
        ! SMART/QUICK modificato upwind
        elseif(insc.eq.2)then

            ! GIULIA UCS - CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,jy
                        !
                        f1(0,j,k) = (rc(0 ,j,k)+ucs(0,j,k))*r(0   ,j,k)- &
                            cgra1(0 ,j,k)
                        f1(jx,j,k)= (rc(jx,j,k)+ucs(jx,j,k))*r(jx+1,j,k)- &
                            cgra1(jx,j,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,jy
                    do i=ip,jx-ip
                        !
                        f1(i,j,k)=ucs(i,j,k) &
                            *.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,jy
                    i=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(2,j,k)
                        rindietro1 =        r(1,j,k)
                        rindietro2 = ip*(2.*r(0,j,k) - r(1,j,k)) &
                            +(1.-ip)*r(0,j,k)
                    else
                        ravanti    =    r(1,j,k)
                        rindietro1 =    r(2,j,k)
                        rindietro2 =    r(3,j,k)
                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )


                    i=jx-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =    r(jx   ,j,k)
                        rindietro1 =    r(jx-1,j, k)
                        rindietro2 =    r(jx-2,j, k)
                    else
                        ravanti    =        r(jx-1,j,k)
                        rindietro1 =        r(jx  ,j,k)
                        rindietro2 = ip*(2.*r(jx+1,j,k) - r(jx,j,k)) &
                            +(1.-ip)*r(jx+1,j,k)
                    end if


                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))



                            else  !(rc(i,j,k).le.0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k)) !my index +1

                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo



                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i+1,j,k)==0)then !solido davanti

                                        f1(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    elseif(tipo(i-1,j,k)==0)then

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)


                                    else ! solido sta nell'altra direzione

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib
                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                                    !       if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                                    !              rindietro2 =r(i,j,k)

                                    !       endif
                                    !         endif

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                                endif !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i+1,j,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f1(i,j,k)=0.
                                    elseif(tipo(i+2,j,k)==0)then !solido avanti
                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                    else !solido di lato

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))

                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if(tipo(i+1,j,k)==1)then
                                    !       if(tipo(i+2,j,k)==0)then
                                    !             rindietro2=r(i+1,j,k)
                                    !       endif
                                    !       endif


                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        !    endif !SMART/quick modificato upwind
        ! SMART/QUICK modificato upwind
        elseif(insc.eq.3)then

            ! GIULIA UCS - CGRA
            !     side left and right
            do i=1,ip
                !
                do k=kparasta,kparaend
                    do j=1,jy
                        !
                        f1(0,j,k) = (rc(0 ,j,k)+ucs(0,j,k))*r(0   ,j,k)- &
                            cgra1(0 ,j,k)
                        f1(jx,j,k)= (rc(jx,j,k)+ucs(jx,j,k))*r(jx+1,j,k)- &
                            cgra1(jx,j,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=1,jy
                    do i=ip,jx-ip
                        !
                        f1(i,j,k)=ucs(i,j,k) &
                            *.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do j=1,jy
                    i=1
                    if(rc(i,j,k).gt.0.)then
                        rindietro1 =        r(1,j,k)

                    else
                        rindietro1 =    r(2,j,k)
                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1


                    i=jx-1
                    if(rc(i,j,k).gt.0.)then
                        rindietro1 =    r(jx-1,j, k)

                    else

                        rindietro1 =        r(jx  ,j,k)

                    end if

                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k).le.0.)

                                f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=1,jy
                        do i=2,jx-2
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i+1,j,k)==0)then !solido davanti

                                        f1(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    else ! solido sta nell'altra direzione

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)


                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib
                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                                    !       if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                                    !              rindietro2 =r(i,j,k)

                                    !       endif
                                    !         endif

                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                                endif !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i+1,j,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f1(i,j,k)=0.
                                    !       elseif(tipo(i+2,j,k)==0)then !solido avanti
                                    !             f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                    else !solido di lato

                                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if(tipo(i+1,j,k)==1)then
                                    !       if(tipo(i+2,j,k)==0)then
                                    !             rindietro2=r(i+1,j,k)
                                    !       endif
                                    !       endif


                                    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        endif !SMART/quick modificato upwind

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
                do j=1+jp,jy-jp
                    !
                    fg=.5*(r(iss,j+1,k)-r(iss,j-1,k))
                    f1(is,j,k)=-f1(is,j,k) &
                        +annit(iss,j,k)*g12(is,j,k)*fg
                !
                end do
                !
                do jj=1,jp
                    !
                    !        check derivative at wall bottom and upper
                    !
                    j=1
                    fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j+1,k)-r(iss,j+2,k))
                    f1(is,j,k)=-f1(is,j,k) &
                        +annit(iss,j,k)*g12(is,j,k)*fg
                    !
                    j=jy
                    fg=.5*(3.*r(iss,jy,k)-4.*r(iss,jy-1,k)+r(iss,jy-2,k))
                    f1(is,j,k)=-f1(is,j,k) &
                        +annit(iss,j,k)*g12(is,j,k)*fg
                !
                end do
            !
            enddo
            !
            is=jx
            iss=jx+1
        !
        enddo
        !
        !
        !     into the field
        do k=kparasta,kparaend
            do i=ip,jx-ip
                !
                do j=1+jp,jy-jp
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
                do jj=1,jp
                    !
                    j=1
                    r0=.5*(r(i,j  ,k)+r(i+1,j  ,k))
                    r1=.5*(r(i,j+1,k)+r(i+1,j+1,k))
                    r2=.5*(r(i,j+2,k)+r(i+1,j+2,k))
                    fg=.5*(-3.*r0+4.*r1-r2)
                    an=.5*(annit(i,j,k)+annit(i+1,j,k))
                    f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
                    !
                    j=jy
                    r0=.5*(r(i,jy  ,k)+r(i+1,jy  ,k))
                    r1=.5*(r(i,jy-1,k)+r(i+1,jy-1,k))
                    r2=.5*(r(i,jy-2,k)+r(i+1,jy-2,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i+1,j,k))
                    f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
                !
                end do
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
        if (myid.eq.0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid.eq.nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else
            kparastak=kparasta
            kparaendk=kparaend
        endif

        is=0
        iss=0

        do i=1,2*ip
            !
            do j=1,jy
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
                    if (myid.eq.0) then

                        k=1
                        fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j,k+1)-r(iss,j,k+2))
                        f1(is,j,k)=f1(is,j,k) &
                            +annit(iss,j,k)*g13(is,j,k)*fg

                    endif
                    !
                    if (myid.eq.nproc-1) then

                        k=jz
                        fg=.5*(3.*r(iss,j,jz)-4.*r(iss,j,jz-1)+r(iss,j,jz-2))
                        f1(is,j,k)=f1(is,j,k) &
                            +annit(iss,j,k)*g13(is,j,k)*fg

                    endif
                !
                end do
            !
            enddo
            !
            is=jx
            iss=jx+1
        !
        enddo


        !
        ! into the field
        !
        do j=1,jy
            do i=ip,jx-ip
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
                    if (myid.eq.0) then

                        k=1
                        r0=.5*(r(i,j,k  )+r(i+1,j,k  ))
                        r1=.5*(r(i,j,k+1)+r(i+1,j,k+1))
                        r2=.5*(r(i,j,k+2)+r(i+1,j,k+2))
                        fg=.5*(-3.*r0+4.*r1-r2)
                        an=.5*(annit(i,j,k)+annit(i+1,j,k))
                        f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg

                    endif
                    !
                    if (myid.eq.nproc-1) then

                        k=jz
                        r0=.5*(r(i,j,jz  )+r(i+1,j,jz  ))
                        r1=.5*(r(i,j,jz-1)+r(i+1,j,jz-1))
                        r2=.5*(r(i,j,jz-2)+r(i+1,j,jz-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        an=.5*(annit(i,j,k)+annit(i+1,j,k))
                        f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg

                    endif
                !
                end do

            end do

        end do


        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !

                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
                        if(i.lt.jx)then
                            if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo

        !
        !
        ! integral of convective plus diffusive off diagonal
        !
        bulk_loc=0.
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    bulk_loc=bulk_loc+f1(i,j,k)-f1(i-1,j,k)
                end do
            end do
        end do
        !
        ! make the value known to all procs to compute the total
        !
        call MPI_ALLREDUCE(bulk_loc,bulk,1,MPI_DOUBLE_PRECISION, &
            MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        ! now bulk is known to all procs
        !
        return
    end subroutine flux1

    subroutine flux2(rc,cgra2,r,insc,tipo)
        !***********************************************************************
        ! compute explicit flux on eta component
        !
        ! convective  vc*(u,v,w)   with central schema or quick
        ! diffusive   nni*g21*d(u,v,w)/d(csi)
        ! diffusive   nni*g23*d(u,v,w)/d(zita)
        !
        use myarrays_wallmodel, only: wfp3, wfp4

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer ierr,insc
        integer m
        integer kparastak,kparaendk
        double precision bulk_loc,bulk2

        integer i,j,k,ii,kk,js,jss
        real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddj,dsj
        real fg,r0,r1,r2,an
        real rc(n1,0:n2,kparasta:kparaend)
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)    !0:n3+1)
        real cgra2(n1,0:n2,kparasta-1:kparaend+1)
        real ravanti,rindietro1,rindietro2
        integer iwall,iwfp
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !
        real phic(n1,n2,kparasta:kparaend)
        real den, num, myv
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM Vc*u
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell

        !     sides bottom and upper
        do j=1,jp
            !
            do k=kparasta,kparaend
                do i=1,jx
                    !
                    f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
                    f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
                !
                enddo
            enddo

        enddo

        !     into the field
        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx
                    !
                    f2(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j+1,k))-cgra2(i,j,k)
                !
                end do
            end do
        end do
        !
        ! quick
        !
        if(insc.eq.1)then

            !     sides 3 and 4
            do k=kparasta,kparaend
                do i=1,jx

                    j=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(i,2,k)
                        rindietro1 =        r(i,1,k)
                        rindietro2 = jp*(2.*r(i,0,k) - r(i,1,k)) &
                            +(1.-jp)*r(i,0,k)
                    else
                        ravanti    = r(i,1,k)
                        rindietro1 = r(i,2,k)
                        rindietro2 = r(i,3,k)
                    end if

                    if(tipo(i,j,k)==2)then
                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    end if

                    j=jy-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    = r(i,jy  ,k)
                        rindietro1 = r(i,jy-1,k)
                        rindietro2 = r(i,jy-2,k)
                    else
                        ravanti    =        r(i,jy-1,k)
                        rindietro1 =        r(i,jy  ,k)
                        rindietro2 = jp*(2.*r(i,jy+1,k) - r(i,jy,k)) &
                            +(1.-jp)*r(i,jy+1,k)
                    end if

                    if(tipo(i,j,k)==2)then
                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                            *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                    end if
                end do
            end do

            !     into the field
            do k=kparasta,kparaend
                do j=2,jy-2
                    do i=1,jx
                        !
                        !         if(tipo(i,j,k)==2)then
                        !      if (rc(i,j,k).gt.0.) then
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                        !      else
                        !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                        !      end if
                        !      end if
                        !

                        if (rc(i,j,k).gt.0.) then

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i,j+1,k)==0)then !solido davanti
                                    f2(i,j,k)=0.
                                elseif(tipo(i,j-1,k)==0)then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                            endif !tipo

                        else !rc(i,j,k).le.0.
                            if(tipo(i,j+1,k)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f2(i,j,k)=0.
                                elseif(tipo(i,j+2,k)==0)then !uso diff centrate
                                    f2(i,j,k)=f2(i,j,k)

                                else !solido di lato
                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))


                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                            endif !tipo

                        endif !rc
                    end do
                end do
            end do

        ! quick upwind
        elseif(insc.eq.2)then

            ! GIULIA UCS - CGRA
            !     side left and right
            do j=1,jp
                !
                do k=kparasta,kparaend
                    do i=1,jx
                        !
                        f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
                        f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=jp,jy-jp
                    do i=1,jx
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,jx
                    j=1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =        r(i,2,k)
                        rindietro1 =        r(i,1,k)
                        rindietro2 = jp*(2.*r(i,0,k) - r(i,1,k)) &
                            +(1.-jp)*r(i,0,k)
                    else
                        ravanti    =    r(i,1,k)
                        rindietro1 =    r(i,2,k)
                        rindietro2 =    r(i,3,k)
                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

                    j=jy-1
                    if(rc(i,j,k).gt.0.)then
                        ravanti    =    r(i,jy  , k)
                        rindietro1 =    r(i,jy-1, k)
                        rindietro2 =    r(i,jy-2, k)
                    else
                        ravanti    =        r(i,jy-1,k)
                        rindietro1 =        r(i,jy  ,k)
                        rindietro2 = jp*(2.*r(i,jy+1,k) - r(i,jy,k)) &
                            +(1.-jp)*r(i,jy+1,k)
                    end if


                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                            else  !(rc(i,j,k).le.0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k)) !my index +1


                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo


                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j+1,k)==0)then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo
                                    elseif(tipo(i,j-1,k)==0)then

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else ! solido sta nell'altra direzione

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                                    !       if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                                    !              rindietro2 =r(i,j,k)

                                    !       endif
                                    !         endif


                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                                endif !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j+1,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f2(i,j,k)=0.
                                    elseif(tipo(i,j+2,k)==0)then !solido avanti
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                    else !solido di lato

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if(tipo(i+1,j,k)==1)then
                                    !       if(tipo(i+2,j,k)==0)then
                                    !             rindietro2=r(i+1,j,k)
                                    !       endif
                                    !       endif


                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        !    endif !SMART/   quick modificato upwind

        ! quick upwind
        elseif(insc.eq.3)then

            ! GIULIA UCS - CGRA
            !     side left and right
            do j=1,jp
                !
                do k=kparasta,kparaend
                    do i=1,jx
                        !
                        f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
                        f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
                    !
                    enddo
                enddo

            enddo
            !
            !     into the field
            do k=kparasta,kparaend
                do j=jp,jy-jp
                    do i=1,jx
                        !
                        f2(i,j,k)=-cgra2(i,j,k)
                    !
                    end do
                end do
            end do


            !     sides 1 and 2  sia per bodyforce 0 o 1
            do k=kparasta,kparaend
                do i=1,jx
                    j=1
                    if(rc(i,j,k).gt.0.)then

                        rindietro1 =        r(i,1,k)

                    else

                        rindietro1 =    r(i,2,k)

                    end if

                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

                    j=jy-1
                    if(rc(i,j,k).gt.0.)then

                        rindietro1 =    r(i,jy-1, k)

                    else

                        rindietro1 =        r(i,jy  ,k)

                    end if


                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

                end do
            end do


            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)


                            else  !(rc(i,j,k).le.0.)

                                f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)


                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparasta,kparaend
                    do j=2,jy-2
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j+1,k)==0)then !solido davanti

                                        f2(i,j,k)=0.
                                    ! dopo diventerà nullo

                                    else ! solido sta nell'altra direzione

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    endif !solido avanti/indietro/lato


                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                                    !       if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                                    !              rindietro2 =r(i,j,k)

                                    !       endif
                                    !         endif

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                                endif !tipo

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j+1,k)==1)then
                                    if(tipo(i,j,k)==0)then !solido indietro
                                        f2(i,j,k)=0.
                                    elseif(tipo(i,j+2,k)==0)then !solido avanti
                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                                    else !solido di lato

                                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                                    endif !solido avanti/indietro/lato

                                else !i non è IB procedo come se non ci fossero gli ib

                                    ! se i+1 è ib faccio la correzione per evitare la scia
                                    !           if(tipo(i+1,j,k)==1)then
                                    !       if(tipo(i+2,j,k)==0)then
                                    !             rindietro2=r(i+1,j,k)
                                    !       endif
                                    !       endif

                                    f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        endif !SMART/    quick modificato upwind
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G21*DU/D(CSI)
        !
        !     sides bottom and upper
        !
        js=0
        jss=0
        iwfp = wfp3

        do j=1,2*jp
            !
            do k=kparasta,kparaend
                !


                do  i=1+ip,jx-ip
                    fg=.5*(r(i+1,jss,k)-r(i-1,jss,k))
                    !         WALL MODEL : all in flucn
                    do iwall = 1,iwfp
                        f2(i,js,k)=-f2(i,js,k)
                    end do
                    !         NO WALL MODEL
                    do iwall = 1,1-iwfp
                        f2(i,js,k)=-f2(i,js,k) &
                            +annit(i,jss,k)*g21(i,js,k)*fg
                    end do
                end do

                !
                do ii=1,ip
                    !     check derivative on sides left and right
                    i=1
                    fg=.5*(-3.*r(i,jss,k)+4.*r(i+1,jss,k)-r(i+2,jss,k))
                    f2(i,js,k)=-f2(i,js,k) &
                        +annit(i,jss,k)*g21(i,js,k)*fg
                    i=jx
                    fg=.5*(3.*r(jx,jss,k)-4.*r(jx-1,jss,k)+r(jx-2,jss,k))
                    f2(i,js,k)=-f2(i,js,k) &
                        +annit(i,jss,k)*g21(i,js,k)*fg
                end do


            !
            enddo
            !
            js=jy
            jss=jy+1
            iwfp = wfp4
        !
        enddo
        !
        ! into the field
        !
        do k=kparasta,kparaend
            do j=jp,jy-jp
                !
                do i=1+ip,jx-ip
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
                    i=jx
                    r0=.5*(r(jx  ,j,k)+r(jx  ,j+1,k))
                    r1=.5*(r(jx-1,j,k)+r(jx-1,j+1,k))
                    r2=.5*(r(jx-2,j,k)+r(jx-2,j+1,k))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i,j+1,k))
                    f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*an*fg
                !
                end do
            !
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G23*DU/D(ZITA)
        !
        ! sides bottom and upper
        !
        ! define computation limit depending on periodicity in z
        !
        if (myid.eq.0) then
            kparastak=kparasta+kp
            kparaendk=kparaend
        else if (myid.eq.nproc-1) then
            kparastak=kparasta
            kparaendk=kparaend-kp
        else
            kparastak=kparasta
            kparaendk=kparaend
        endif

        js=0
        jss=0
        iwfp = wfp3
        do j=1,2*jp
            !
            do i=1,jx
                !
                do k=kparastak,kparaendk
                    !
                    fg=.5*(r(i,jss,k+1)-r(i,jss,k-1))

                    !           WALL MODEL: all in flucn
                    do iwall=1,iwfp
                        f2(i,js,k)=f2(i,js,k)
                    end do
                    !           NO WALL MODEL
                    do iwall=1,1-iwfp
                        f2(i,js,k)=f2(i,js,k) &
                            +annit(i,jss,k)*g23(i,js,k)*fg
                    end do
                !
                end do
                !

                do kk=1,kp
                    !
                    ! check derivative on sides back and front
                    !
                    if (myid.eq.0) then

                        k=1
                        fg=.5*(-3.*r(i,jss,k)+4.*r(i,jss,k+1)-r(i,jss,k+2))
                        f2(i,js,k)=f2(i,js,k) &
                            +annit(i,jss,k)*g23(i,js,k)*fg

                    endif
                    !
                    if (myid.eq.nproc-1) then

                        k=jz
                        fg=.5*(3.*r(i,jss,jz)-4.*r(i,jss,jz-1)+r(i,jss,jz-2))
                        f2(i,js,k)=f2(i,js,k) &
                            +annit(i,jss,k)*g23(i,js,k)*fg

                    endif
                !
                end do
            !
            enddo
            !
            js=jy
            jss=jy+1
            iwfp = wfp4
        !
        enddo
        !
        ! into the field
        !
        do j=jp,jy-jp
            do i=1,jx
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
                    if (myid.eq.0) then

                        k=1
                        r0=.5*(r(i,j,k  )+r(i,j+1,k  ))
                        r1=.5*(r(i,j,k+1)+r(i,j+1,k+1))
                        r2=.5*(r(i,j,k+2)+r(i,j+1,k+2))
                        fg=.5*(-3.*r0+4.*r1-r2)
                        an=.5*(annit(i,j,k)+annit(i,j+1,k))
                        f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg

                    endif
                    !
                    if (myid.eq.nproc-1) then

                        k=jz
                        r0=.5*(r(i,j,jz  )+r(i,j+1,jz  ))
                        r1=.5*(r(i,j,jz-1)+r(i,j+1,jz-1))
                        r2=.5*(r(i,j,jz-2)+r(i,j+1,jz-2))
                        fg=.5*(3.*r0-4.*r1+r2)
                        an=.5*(annit(i,j,k)+annit(i,j+1,k))
                        f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg

                    endif
                !
                end do
            !
            enddo
        enddo


        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx
                    !

                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
                        if(j.lt.jy)then
                            if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo


        !
        ! integral on f2
        bulk_loc=0.
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    bulk_loc=bulk_loc+f2(i,j,k)-f2(i,j-1,k)
                end do
            end do
        end do

        ! make the value known to all procs
        call MPI_ALLREDUCE(bulk_loc,bulk2,1,MPI_DOUBLE_PRECISION, &
            MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        ! now bulk is known by all procs

        bulk=bulk+bulk2

        return
    end subroutine flux2

    subroutine flux3(rc,cgra3,r,insc,tipo)
        !***********************************************************************
        ! compute explicit flux on component zita
        !
        ! convective  wc*(u,v,w)   with central schema or quick
        ! diffusive   nni*g31*d(u,v,w)/d(csi)
        ! diffusive   nni*g32*d(u,v,w)/d(zita)
        !
        use myarrays_LC
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer ierr,status(MPI_STATUS_SIZE)
        integer m,insc
        integer kparastal,kparaendl
        integer kparastall,kparaendll
        double precision bulk_loc,bulk3
        integer i,j,k,ii,jj,ks,kss
        real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddk,dsk
        real fg,r0,r1,r2,an
        real rc(n1,n2,kparasta-1:kparaend) !0:n3)
        real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real cgra3(n1,n2,kparasta-1:kparaend+1)
        real ravanti,rindietro1,rindietro2
        integer req1,req2
        integer istatus(MPI_STATUS_SIZE)
        !      integer tipo2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
        !      integer tipo2(0:n1+1,0:n2+1,0:n3+1) c7
        integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

        real phic(n1,n2,kparasta:kparaend)
        real den, num, myv
        !-----------------------------------------------------------------------
        ! CONVECTIVE TERM  Wc*u:
        ! implemented with central scheme or with quick depending on settings
        ! in points closer to wall it uses uf as ghost cell
        !
        !     sides back and front
        do k=1,kp

            do j=1,jy
                do i=1,jx
                    !
                    if (myid.eq.0) then
                        f3(i,j,0)  = (rc(i,j,0 )+wcs(i,j,0))*r(i,j,0   )- &
                            cgra3(i,j,0 )
                    endif

                    if (myid.eq.nproc-1) then
                        f3(i,j,jz) = (rc(i,j,jz)+wcs(i,j,jz))*r(i,j,jz+1)- &
                            cgra3(i,j,jz)
                    endif
                !
                enddo
            enddo
        !
        enddo
        !
        !     into the field
        !
        if (myid.eq.0) then
            kparastal=kp
            kparaendl=kparaend
        else if (myid.eq.nproc-1) then
            kparastal=kparasta
            kparaendl=kparaend-kp
        else
            kparastal=kparasta
            kparaendl=kparaend
        endif


        do k=kparastal,kparaendl
            do j=1,jy
                do i=1,jx
                    !
                    f3(i,j,k)=(rc(i,j,k)+wcs(i,j,k)) &
                        *.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
                !
                end do
            end do
        end do

        ! quick
        if(insc.eq.1)then

            if (myid.eq.0) then
                kparastall=kparasta +2
                kparaendll=kparaend
            else if (myid.eq.nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            endif

            !     sides 5 and 6

            if(myid .eq. 0)then
                do j=1,jy
                    do i=1,jx

                        k=1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    =        r(i,j,2)
                            rindietro1 =        r(i,j,1)
                            rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti    =    r(i,j,1)
                            rindietro1 =    r(i,j,2)
                            rindietro2 =    r(i,j,3)
                        end if

                        if(tipo(i,j,k)==2)then
                            f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                        end if

                    end do
                end do
            end if

            if(myid .eq. nproc-1)then
                do j=1,jy
                    do i=1,jx

                        k=jz-1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    = r(i,j,jz  )
                            rindietro1 = r(i,j,jz-1)
                            rindietro2 = r(i,j,jz-2)
                        else
                            ravanti    =        r(i,j,jz-1)
                            rindietro1 =        r(i,j,jz  )
                            rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                                +(1.-kp)*r(i,j,jz+1)
                        end if

                        if(tipo(i,j,k)==2)then
                            f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                *.125*( -ravanti +2.*rindietro1 - rindietro2 )
                        end if
                    end do
                end do
            end if

            !     into the field
            do k=kparastall,kparaendll
                do j=1,jy
                    do i=1,jx
                        !
                        !      if(tipo(i,j,k)==2)then
                        !      if (rc(i,j,k).gt.0.) then
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                        !      else
                        !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
                        !     > .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                        !      end if
                        !      end if

                        if (rc(i,j,k).gt.0.) then

                            if(tipo(i,j,k)==1)then !ib cell
                                if(tipo(i,j,k+1)==0)then !solido davanti
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k-1)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else ! solido sta nell'altra direzione quick
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))

                                endif !solido avanti/indietro/lato
                            else !i è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                            endif !tipo

                        else !rc(i,j,k).le.0.

                            if(tipo(i,j,k+1)==1)then
                                if(tipo(i,j,k)==0)then !solido indietro
                                    f3(i,j,k)=0.
                                elseif(tipo(i,j,k+2)==0)then !uso diff centrate
                                    f3(i,j,k)=f3(i,j,k)
                                else !solido di lato
                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                        .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                                endif !solido avanti/indietro/lato
                            else !i+1 è fluido allora quick normale
                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                                    .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                            endif !tipo

                        endif !rc
                    end do
                end do
            end do

        ! SMART/quick modificato upwind
        elseif(insc.eq.2)then

            ! GIULIA UCS - CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,jy
                    do i=1,jx
                        !
                        if (myid.eq.0) then
                            f3(i,j,0)  = (rc(i,j,0 )+wcs(i,j,0))*r(i,j,0   )- &
                                cgra3(i,j,0 )
                        endif

                        if (myid.eq.nproc-1) then
                            f3(i,j,jz) = (rc(i,j,jz)+wcs(i,j,jz))*r(i,j,jz+1)- &
                                cgra3(i,j,jz)
                        endif
                    !
                    enddo
                enddo
            !
            enddo

            if (myid.eq.0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid.eq.nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            endif


            do k=kparastal,kparaendl
                do j=1,jy
                    do i=1,jx
                        !
                        f3(i,j,k)=wcs(i,j,k) &
                            *.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid.eq.0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid.eq.nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            endif

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if(myid .eq. 0)then
                do j=1,jy
                    do i=1,jx
                        k=1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    =        r(i,j,2)
                            rindietro1 =        r(i,j,1)
                            rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti    =    r(i,j,1)
                            rindietro1 =    r(i,j,2)
                            rindietro2 =    r(i,j,3)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )


                    end do
                end do

            end if !myid==0

            if(myid .eq. nproc-1)then
                do j=1,jy
                    do i=1,jx

                        k=jz-1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    = r(i,j,jz  )
                            rindietro1 = r(i,j,jz-1)
                            rindietro2 = r(i,j,jz-2)
                        else
                            ravanti    =        r(i,j,jz-1)
                            rindietro1 =        r(i,j,jz  )
                            rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                                +(1.-kp)*r(i,j,jz+1)
                        end if


                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                            *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                            else  !(rc(i,j,k).le.0.)

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                    *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2)) !my index +1


                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j,k+1)==0)then !solido davanti

                                        f3(i,j,k)=0.

                                    elseif(tipo(i,j,k-1)==0)then !solido dietro

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)
                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                                    endif !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k-1)==1)then !la cella prima è un ib
                                    !       if(tipo(i,j,k-2)==0)then !perché considero anche questa?
                                    !         rindietro2 =r(i,j,k)
                                    !       endif
                                    !         endif !fine k-1 IB


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))


                                endif !tip02

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j,k+1)==1)then !ib
                                    if(tipo(i,j,k)==0)then !solido avanti
                                        f3(i,j,k)=0.
                                    elseif(tipo(i,j,k+2)==0)then !solido indietro
                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)
                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                            *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))


                                    endif !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k+1)==1)then
                                    !       if(tipo(i,j,k+2)==0)then
                                    !       rindietro2=r(i,j,k+1)
                                    !       endif
                                    !         endif

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                                        *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        !    endif !SMART
        !-----SMART /quick modificato upwind
        !  upwind
        elseif(insc.eq.3)then

            ! GIULIA UCS - CGRA
            !     side left and right
            !     sides back and front
            do k=1,kp

                do j=1,jy
                    do i=1,jx
                        !
                        if (myid.eq.0) then
                            f3(i,j,0)  = (rc(i,j,0 )+wcs(i,j,0))*r(i,j,0   )- &
                                cgra3(i,j,0 )
                        endif

                        if (myid.eq.nproc-1) then
                            f3(i,j,jz) = (rc(i,j,jz)+wcs(i,j,jz))*r(i,j,jz+1)- &
                                cgra3(i,j,jz)
                        endif
                    !
                    enddo
                enddo
            !
            enddo

            if (myid.eq.0) then
                kparastal=kp
                kparaendl=kparaend
            else if (myid.eq.nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            endif


            do k=kparastal,kparaendl
                do j=1,jy
                    do i=1,jx
                        !
                        f3(i,j,k)=wcs(i,j,k) &
                            *.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
                    !
                    end do
                end do
            end do

            if (myid.eq.0) then
                kparastall=2
                kparaendll=kparaend
            else if (myid.eq.nproc-1) then
                kparastall=kparasta
                kparaendll=kparaend-2
            else
                kparastall=kparasta
                kparaendll=kparaend
            endif

            !     sides 5 and 6 sia per bodyforce 0 o 1
            if(myid .eq. 0)then
                do j=1,jy
                    do i=1,jx
                        k=1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    =        r(i,j,2)
                            rindietro1 =        r(i,j,1)
                            rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                                +(1.-kp)*r(i,j,0)
                        else
                            ravanti    =    r(i,j,1)
                            rindietro1 =    r(i,j,2)
                            rindietro2 =    r(i,j,3)
                        end if

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1


                    end do
                end do

            end if !myid==0

            if(myid .eq. nproc-1)then
                do j=1,jy
                    do i=1,jx

                        k=jz-1
                        if(rc(i,j,k).gt.0.)then
                            ravanti    = r(i,j,jz  )
                            rindietro1 = r(i,j,jz-1)
                            rindietro2 = r(i,j,jz-2)
                        else
                            ravanti    =        r(i,j,jz-1)
                            rindietro1 =        r(i,j,jz  )
                            rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                                +(1.-kp)*r(i,j,jz+1)
                        end if


                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1

                    end do
                end do
            end if



            !     into the field without bodyforce
            if (.not.bodyforce) then
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                            else  !(rc(i,j,k).le.0.)

                                f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)


                            endif !rc</>0
                        end do
                    end do
                end do

            !
            else !bodyforce attivo
                do k=kparastall,kparaendll
                    do j=1,jy
                        do i=1,jx
                            if(rc(i,j,k).gt.0.)then
                                if(tipo(i,j,k)==1)then !ib cell
                                    if(tipo(i,j,k+1)==0)then !solido davanti

                                        f3(i,j,k)=0.


                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                                    endif !solido avanti/indietro/lati

                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k-1)==1)then !la cella prima è un ib
                                    !       if(tipo(i,j,k-2)==0)then !perché considero anche questa?
                                    !         rindietro2 =r(i,j,k)
                                    !       endif
                                    !         endif !fine k-1 IB


                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)


                                endif !tip02

                            else !rc(i,j,k).le.0.

                                if(tipo(i,j,k+1)==1)then !ib
                                    if(tipo(i,j,k)==0)then !solido avanti
                                        f3(i,j,k)=0.

                                    else !solido lati

                                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)


                                    endif !solido avanti indietro e lati


                                else !k non è IB procedo come se non ci fossero gli ib

                                    ! se k-1 è ib faccio la correzione per evitare la scia
                                    !         if(tipo(i,j,k+1)==1)then
                                    !       if(tipo(i,j,k+2)==0)then
                                    !       rindietro2=r(i,j,k+1)
                                    !       endif
                                    !         endif

                                    f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

                                endif !tipo

                            endif !rc

                        end do
                    end do
                end do

            endif !body yes/no
        endif !SMART
        !-----SMART /quick modificato upwind

        !----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G31*DU/D(CSI)
        !
        !     sides back and front
        ks=0
        kss=0
        do k=1,2*kp

            if (myid.eq.(k-1)*(nproc-1)) then

                do j=1,jy
                    !
                    do i=1+ip,jx-ip
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
                        i=jx
                        fg=.5*(3.*r(jx,j,kss)-4.*r(jx-1,j,kss)+r(jx-2,j,kss))
                        f3(i,j,ks)=-f3(i,j,ks) &
                            +annit(i,j,kss)*g31(i,j,ks)*fg
                    !
                    end do
                !
                enddo
            !
            endif

            ks=jz
            kss=jz+1
        !
        enddo
        !
        !     into the field
        do k=kparastal,kparaendl
            do j=1,jy
                !
                do i=1+ip,jx-ip
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
                    i=jx
                    r0=.5*(r(jx  ,j,k)+r(jx  ,j,k+1))
                    r1=.5*(r(jx-1,j,k)+r(jx-1,j,k+1))
                    r2=.5*(r(jx-2,j,k)+r(jx-2,j,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*an*fg
                !
                end do
            !
            enddo
        enddo
        !
        !-----------------------------------------------------------------------
        ! DIFFUSIVE TERM NNI*G32*DU/D(ZITA)
        !
        !     sides back and front
        !
        ks=0
        kss=0

        do k=1,2*kp

            if (myid.eq.(k-1)*(nproc-1)) then

                do i=1,jx
                    !
                    do j=1+jp,jy-jp
                        !
                        fg=.5*(r(i,j+1,kss)-r(i,j-1,kss))
                        f3(i,j,ks)=f3(i,j,ks) &
                            +annit(i,j,kss)*g32(i,j,ks)*fg
                    !
                    end do
                    !
                    do jj=1,jp
                        !
                        ! check derivative on sides back and front
                        !
                        j=1
                        fg=.5*(-3.*r(i,j,kss)+4.*r(i,j+1,kss)-r(i,j+2,kss))
                        f3(i,j,ks)=f3(i,j,ks) &
                            +annit(i,j,kss)*g32(i,j,ks)*fg
                        !
                        j=jy
                        fg=.5*(3.*r(i,jy,kss)-4.*r(i,jy-1,kss)+r(i,jy-2,kss))
                        f3(i,j,ks)=f3(i,j,ks) &
                            +annit(i,j,kss)*g32(i,j,ks)*fg
                    !
                    end do
                !
                enddo
            !
            endif

            ks=jz
            kss=jz+1
        !
        enddo

        !
        !     into the field
        do k=kparastal,kparaendl
            do i=1,jx
                !
                do j=1+jp,jy-jp
                    !
                    r1=.5*(r(i,j-1,k)+r(i,j-1,k+1))
                    r2=.5*(r(i,j+1,k)+r(i,j+1,k+1))
                    fg=.5*(r2-r1)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
                !
                end do
                !
                do jj=1,jp
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
                    j=jy
                    r0=.5*(r(i,jy  ,k)+r(i,jy  ,k+1))
                    r1=.5*(r(i,jy-1,k)+r(i,jy-1,k+1))
                    r2=.5*(r(i,jy-2,k)+r(i,jy-2,k+1))
                    fg=.5*(3.*r0-4.*r1+r2)
                    an=.5*(annit(i,j,k)+annit(i,j,k+1))
                    f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
                !
                end do
            !
            enddo
        enddo
        !

        do j=1,jy
            do i=1,jx
                do k=kparastal,kparaendl
                    !
                    if (bodyforce) then
                        if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
                        if(k.lt.jz)then
                            if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
                        endif
                    endif
                enddo
            enddo
        enddo
        !     pass f3 at the border between procs
        if(myid.eq.0)then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if(myid.eq.nproc-1)then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else
            leftpem=leftpe
            rightpem=rightpe
        endif

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_ISEND(f3(1,1,kparaend),jx*jy,MPI_REAL_SD, &
                rightpem ,tagrs,MPI_COMM_WORLD,req1,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_IRECV(f3(1,1,kparasta-1),jx*jy,MPI_REAL_SD, &
                leftpem  ,taglr,MPI_COMM_WORLD,req2,ierr)
        endif

        if(rightpem /= MPI_PROC_NULL) then
            call MPI_WAIT(req1,istatus,ierr)
        endif
        if(leftpem /= MPI_PROC_NULL) then
            call MPI_WAIT(req2,istatus,ierr)
        endif

        !     integral for f3
        bulk_loc=0.

        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    bulk_loc=bulk_loc+f3(i,j,k)-f3(i,j,k-1)
                end do
            end do
        end do

        !     make the value known to all procs

        call MPI_ALLREDUCE(bulk_loc,bulk3,1,MPI_DOUBLE_PRECISION, &
            MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        !    now bulk is known by all procs

        bulk=bulk+bulk3
        !
        return

    end subroutine flux3

end module flu_module
