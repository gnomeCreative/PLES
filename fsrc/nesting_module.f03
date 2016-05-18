!***********************************************************************
module nesting_module
    !***********************************************************************
    !      for nesting procedure
    !-----------------------------------------------------------------------
    use contour_module
    use orlansky_module
    use buffer_bodyforce_module, only: ti_pom_new,ti_pom_old,ti_pom_fin
    !
    use mysending, only: kparasta,kparaend,myid,nproc,MPI_REAL_SD
    use scala3
    !
    use mpi

    use iso_c_binding

    real,allocatable,dimension(:,:) :: uo1,vo1,wo1,un1,vn1,wn1
    real,allocatable,dimension(:,:) :: uo2,vo2,wo2,un2,vn2,wn2
    real,allocatable,dimension(:,:) :: uo5,vo5,wo5,un5,vn5,wn5
    real,allocatable,dimension(:,:) :: uo6,vo6,wo6,un6,vn6,wn6

    integer n_ti_pom
    logical,bind(C) :: termina
    integer ntke

contains

    subroutine nesting(ti)
        !***********************************************************************
        ! subroutine for nesting procedure
        use inflow_module, only: redistribuzione
        use mysettings, only: freesurface
        !

        implicit none

        !-----------------------------------------------------------------------
        ! array declaration
        integer tii
        real ti
        !-----------------------------------------------------------------------
        ! check that dt not larger than t_new_pom

        if(ti<=ti_pom_new)then

            call interpola_pareti_pom(ti)
            call contrin_lat

            if(.not.freesurface)then !if freesurface off
                if(myid.eq.0)then
                    write(*,*)'freesurface is off and entering redistribuzione.'
                end if
                call redistribuzione()
            elseif(freesurface)then !if freesurface is on
                if(myid.eq.0)then
                    write(*,*)'free surface is on. redistribuzione skipped.'
                end if
            end if !if freesurface on/off


        elseif(ti>ti_pom_new.and.ti<ti_pom_fin)then

            if(myid==0)write(*,*)'NESTING READ NEW PLANE'

            ti=ti_pom_new
            ti_pom_old=ti_pom_new

            if(infout1 /=0)then
                uo1=un1
                vo1=vn1
                wo1=wn1
                rhovo1=rhovn1
            end if

            if(infout2 /=0)then
                uo2=un2
                vo2=vn2
                wo2=wn2
                rhovo2=rhovn2
            end if

            if(infout5 /=0)then
                uo5=un5
                vo5=vn5
                wo5=wn5
                rhovo5=rhovn5
            end if

            if(infout6 /=0)then
                uo6=un6
                vo6=vn6
                wo6=wn6
                rhovo6=rhovn6
            end if

            call leggi_pareti_pom(ti)
            call interpola_pareti_pom(ti)
            call contrin_lat

            if(.not.freesurface)then !if freesurface off
                if(myid.eq.0)then
                    write(*,*)'freesurface is off and entering redistribuzione.'
                end if
                call redistribuzione()
            elseif(freesurface)then !if freesurface is on
                if(myid.eq.0)then
                    write(*,*)'free surface is on. redistribuzione skipped.'
                end if
            end if !if freesurface on/off


        end if
        !
        ! if ti greater than t_pom_fin the computation end
        if(ti>=ti_pom_fin)then
            ti=ti_pom_fin
            termina=.true.
        !cc      close(1313)
        end if

        return
    end subroutine nesting

    subroutine prepare_nesting(ti,unit1,unit2,unit3,unit4,unit5,unit6)

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        !      integer imoist
        integer unit1,unit2,unit3
        integer unit4,unit5,unit6
        integer i,j,k,n,isc,kpst,kpnd

        real ti,ntireal
        real var1,var2,var3,var4,var(nscal)

        character*1 commento
        !-----------------------------------------------------------------------
        ntke = 1
        !-----------------------------------------------------------------------
        !     array allocation and initialization for nesting
        if(infout1 /=0)then
            !     side   1
            allocate(uo1(0:jy+1,0:jz+1))
            allocate(vo1(0:jy+1,0:jz+1))
            allocate(wo1(0:jy+1,0:jz+1))
            allocate(rhovo1(nscal,0:jy+1,0:jz+1))
            allocate(ucp1(jy,kparasta:kparaend))

            uo1=0.
            vo1=0.
            wo1=0.
            rhovo1=0.
            ucp1 = 0.


            read(unit1,*)commento
            read(unit1,*)ntireal,ti_pom_fin
            n_ti_pom = int(ntireal)
            read(unit1,*)ti_pom_old

            allocate(tkepom1(0:jy+1,kparasta-1:kparaend+1,n_ti_pom)) !face 1
            tkepom1 = 0.

            do k=0,jz+1
                do j=0,jy+1
                    read(unit1,*)var1,var2,var3,var4,(var(n),n=1,nscal)

                    if(k.ge.kparasta.and.k.le.kparaend)then
                        uo1(j,k) = var1*index_out1(j,k)
                        vo1(j,k) = var2*index_out1(j,k)
                        wo1(j,k) = var3*index_out1(j,k)
                        tkepom1(j,k,ntke)=var4*index_out1(j,k)
                        do isc = 1,nscal
                            rhovo1(isc,j,k) = var(isc)*index_rho1(j,k)
                        end do
                    endif
                end do
            end do

        end if

        !-----------------------------------------------------------------------
        !     side   2
        if(infout2 /=0)then
            allocate(uo2(0:jy+1,0:jz+1))
            allocate(vo2(0:jy+1,0:jz+1))
            allocate(wo2(0:jy+1,0:jz+1))
            allocate(rhovo2(nscal,0:jy+1,0:jz+1))
            allocate(ucp2(jy,kparasta:kparaend))

            uo2=0.
            vo2=0.
            wo2=0.
            rhovo2=0.
            ucp2 = 0.


            read(unit2,*)commento
            read(unit2,*)ntireal,ti_pom_fin
            n_ti_pom = int(ntireal)
            read(unit2,*)ti_pom_old

            allocate(tkepom2(0:jy+1,kparasta-1:kparaend+1,n_ti_pom)) !face 2
            tkepom2 = 0.

            do k=0,jz+1
                do j=0,jy+1
                    read(unit2,*)var1,var2,var3,var4,(var(n),n=1,nscal)
                    if(k.ge.kparasta.and.k.le.kparaend)then
                        uo2(j,k) = var1*index_out2(j,k)
                        vo2(j,k) = var2*index_out2(j,k)
                        wo2(j,k) = var3*index_out2(j,k)
                        tkepom2(j,k,ntke)=var4*index_out2(j,k)
                        do isc = 1,nscal
                            rhovo2(isc,j,k) = var(isc)*index_rho2(j,k)
                        end do
                    endif
                enddo
            enddo

        end if
        !-----------------------------------------------------------------------
        !     side   3
        ! THIS WAS WORKING ONLY WITH imoist
        !        if(imoist == 1)then
        !            if(infout3 /=0)then
        !                allocate(    tauw3nest(n_ti_pom))
        !                allocate(scalarflux3nest(nscal,n_ti_pom))
        !
        !                tauw3nest = 0.
        !                scalarflux3nest = 0.
        !
        !                read(unit3,*)commento
        !                do i=1,n_ti_pom
        !                    read(unit3,*) ti_pom_old
        !                    read(unit3,*) tauw3nest(i),(scalarflux3nest(isc,i),isc=1,nscal)
        !                end do
        !                close(83)
        !
        !            end if
        !
        !
        !            !      side 4
        !            if(infout4 /=0)then
        !
        !            !      allocate(    Qvflux4nest(n_ti_pom))
        !            !      allocate(TKEflux4nest(n_ti_pom))
        !            !      allocate( vapflux4nest(n_ti_pom))
        !            !      allocate(scalarflux4nest(nscal,n_ti_pom))
        !
        !
        !            !      do i=1,n_ti_pom
        !            !  read(84,*)
        !            !       end do
        !            !      close(84)
        !
        !            end if
        !        end if

        !-----------------------------------------------------------------------
        !     side   5
        if(infout5 /=0)then
            allocate(uo5(0:jx+1,0:jy+1))
            allocate(vo5(0:jx+1,0:jy+1))
            allocate(wo5(0:jx+1,0:jy+1))
            allocate(rhovo5(nscal,0:jx+1,0:jy+1))
            allocate(wcp5(jx,jy))

            uo5=0.
            vo5=0.
            wo5=0.
            rhovo5=0.
            wcp5 = 0.


            read(unit5,*)commento
            read(unit5,*)ntireal,ti_pom_fin
            n_ti_pom = int(ntireal)
            read(unit5,*)ti_pom_old

            allocate(tkepom5(0:jx+1,0:jy+1,n_ti_pom))
            tkepom5 = 0.

            do j=0,jy+1
                do i=0,jx+1

                    read(unit5,*)var1,var2,var3,var4,(var(n),n=1,nscal)

                    uo5(i,j) = var1*index_out5(i,j)
                    vo5(i,j) = var2*index_out5(i,j)
                    wo5(i,j) = var3*index_out5(i,j)
                    tkepom5(i,j,ntke)=var4*index_out5(i,j)
                    do isc = 1,nscal
                        rhovo5(isc,i,j) = var(isc)*index_rho5(i,j)
                    end do
                end do
            end do

        end if
        !-----------------------------------------------------------------------
        !     side   6
        if(infout6 /=0)then
            allocate(uo6(0:jx+1,0:jy+1))
            allocate(vo6(0:jx+1,0:jy+1))
            allocate(wo6(0:jx+1,0:jy+1))
            allocate(rhovo6(nscal,0:jx+1,0:jy+1))
            allocate(wcp6(jx,jy))

            uo6=0.
            vo6=0.
            wo6=0.
            rhovo6=0.
            wcp6=0.


            read(unit6,*)commento
            read(unit6,*)ntireal,ti_pom_fin
            n_ti_pom = int(ntireal)
            read(unit6,*)ti_pom_old

            allocate(tkepom6(0:jx+1,0:jy+1,n_ti_pom))
            tkepom6 = 0.

            do j=0,jy+1
                do i=0,jx+1

                    read(unit6,*)var1,var2,var3,var4,(var(n),n=1,nscal)

                    uo6(i,j) = var1*index_out6(i,j)
                    vo6(i,j) = var2*index_out6(i,j)
                    wo6(i,j) = var3*index_out6(i,j)
                    tkepom6(i,j,ntke)=var4*index_out6(i,j)
                    do isc = 1,nscal
                        rhovo6(isc,i,j) = var(isc)*index_rho6(i,j)
                    end do
                end do
            end do

        end if
        !-----------------------------------------------------------------------
        !c      ti_pom_old = ti_pom_old * 3600. ! in seconds
        !c      ti_pom_fin = ti_pom_fin * 3600.

        if(myid.eq.0)then
            write(*,*)'----------------------'
            write(*,*)'read pareti.dat '
            write(*,*)'POM time number:',n_ti_pom
            write(*,*)'POM initial time:',ti_pom_old
            write(*,*)'POM final time:',ti_pom_fin
            write(*,*)'----------------------'
        end if
        !-----------------------------------------------------------------------
        !
        if(.not.potenziale)then
            !        CASE START WITH RESTART FILE
            !        allocation for nesting
            !        side 1
            if(infout1 /=0)then
                allocate(un1(0:jy+1,0:jz+1))
                allocate(vn1(0:jy+1,0:jz+1))
                allocate(wn1(0:jy+1,0:jz+1))
                allocate(rhovn1(nscal,0:jy+1,0:jz+1))
                un1 = 0.
                vn1 = 0.
                wn1 = 0.
                rhovn1 = 0.
            end if
            !        side 2
            if(infout2 /=0)then
                allocate(un2(0:jy+1,0:jz+1))
                allocate(vn2(0:jy+1,0:jz+1))
                allocate(wn2(0:jy+1,0:jz+1))
                allocate(rhovn2(nscal,0:jy+1,0:jz+1))
                un2 = 0.
                vn2 = 0.
                wn2 = 0.
                rhovn2 = 0.
            end if
            !        side 5
            if(infout5 /=0)then
                allocate(un5(0:jx+1,0:jy+1))
                allocate(vn5(0:jx+1,0:jy+1))
                allocate(wn5(0:jx+1,0:jy+1))
                allocate(rhovn5(nscal,0:jx+1,0:jy+1))
                un5 = 0.
                vn5 = 0.
                wn5 = 0.
                rhovn5 = 0.
            end if
            !        side 6
            if(infout6 /=0)then
                allocate(un6(0:jx+1,0:jy+1))
                allocate(vn6(0:jx+1,0:jy+1))
                allocate(wn6(0:jx+1,0:jy+1))
                allocate(rhovn6(nscal,0:jx+1,0:jy+1))
                un6 = 0.
                vn6 = 0.
                wn6 = 0.
                rhovn6 = 0.
            end if

            call restart(ti)
            call leggi_pareti_pom(ti)
            call interpola_pareti_pom(ti)
        else ! if start with potential flow
            ti = ti_pom_old

            u=0.;v=0.;w=0. ! set a zero velocity field

            !         inizialize density field
            call interp()

        !          if(myid==0)then
        !             kpst = kparasta
        !         kpnd = kparaend + deepr
        !          elseif(myid==nproc-1)then
        !             kpst = kparasta - deepl
        !         kpnd = kparaend
        !          else
        !             kpst = kparasta - deepl
        !         kpnd = kparaend + deepr
        !          end if

        !          do k=kpst,kpnd
        !          do j=1,jy
        !          do i=1,jx
        !            do isc=1,nscal
        !           rhov(isc,i,j,k)=rhovo1(isc,j,kparasta)
        !            end do
        !          end do
        !          end do
        !          end do

        !          do k=kpst,kpnd
        !          do i=1,jx
        !            do isc=1,nscal
        !          rhov(isc,i,0,k)=rhov(isc,i,1,k)
        !            end do
        !          end do
        !          end do

        end if
        !-----------------------------------------------------------------------
        !     set boundary conditions
        !     lateral values for velocity vector
        !
        !     side 1
        if(infout1 /= 0)then
            do k=kparasta,kparaend
                do j=1,jy

                    u(0,j,k)=up1(j,k)
                    v(0,j,k)=vp1(j,k)
                    w(0,j,k)=wp1(j,k)

                    do isc=1,nscal
                        rhov(isc,0   ,j,k)=rhovp1(isc,j,k)
                    end do
                end do
            end do
        end if

        !     side 1
        if(infout2 /= 0)then
            do k=kparasta,kparaend
                do j=1,jy

                    u(jx+1,j,k)=up2(j,k)
                    v(jx+1,j,k)=vp2(j,k)
                    w(jx+1,j,k)=wp2(j,k)

                    do isc=1,nscal
                        rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
                    end do
                end do
            end do
        end if

        !     sides 5
        if(infout5 /= 0)then
            if(myid==0)then
                do j=1,jy
                    do i=1,jx
                        u(i,j,0)=up5(i,j)
                        v(i,j,0)=vp5(i,j)
                        w(i,j,0)=wp5(i,j)
                        do isc=1,nscal
                            rhov(isc,i,j,0   )=rhovp5(isc,i,j)
                        end do
                    end do
                end do
            end if
        end if

        !     sides 6
        if(infout6 /= 0)then
            if(myid==nproc-1)then
                do j=1,jy
                    do i=1,jx
                        u(i,j,jz+1)=up6(i,j)
                        v(i,j,jz+1)=vp6(i,j)
                        w(i,j,jz+1)=wp6(i,j)
                        do isc=1,nscal
                            rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)
                        end do
                    end do
                end do
            end if
        end if
        !-----------------------------------------------------------------------

        return
    end subroutine prepare_nesting

    subroutine interp()
        !***********************************************************************
        use myarrays_metri3
        use myarrays_velo3
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,isc

        real xp1,xp2
        real zp1,zp2
        real xc,zc

        real dist1,dist2,dist_tot

        real x_p0(0:jx,0:jy)
        real y_p0(0:jx,0:jy)
        real z_p0(0:jx,0:jy)

        real x_pjz(0:jx,0:jy)
        real y_pjz(0:jx,0:jy)
        real z_pjz(0:jx,0:jy)

        real r_ew,r_sn

        real buffer(3*(jx+1)*(jy+1))
        integer icount,numcount
        integer ierr
        !-----------------------------------------------------------------------

        numcount = (jx+1)*(jy+1)

        !     plane comunication
        !     plane 0
        buffer = 0.
        if(myid==0)then

            k=0
            do j=0,jy
                do i=0,jx
                    x_p0(i,j)=x(i,j,k)
                    y_p0(i,j)=y(i,j,k)
                    z_p0(i,j)=z(i,j,k)
                end do
            end do

            icount = 1
            do j=0,jy
                do i=0,jx
                    buffer(           icount) = x_p0(i,j)
                    buffer(  numcount+icount) = y_p0(i,j)
                    buffer(2*numcount+icount) = z_p0(i,j)

                    icount = icount + 1
                end do
            end do
        end if

        call MPI_BCAST(buffer(1),3*numcount,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)

        icount = 1
        do j=0,jy
            do i=0,jx
                x_p0(i,j) = buffer(           icount)
                y_p0(i,j) = buffer(  numcount+icount)
                z_p0(i,j) = buffer(2*numcount+icount)

                icount = icount + 1
            end do
        end do

        !-----------------------------------------------------------------------
        !     plane comunication
        !     plane jz
        buffer = 0.
        if(myid==nproc-1)then

            k=jz
            do j=0,jy
                do i=0,jx
                    x_pjz(i,j)=x(i,j,k)
                    y_pjz(i,j)=y(i,j,k)
                    z_pjz(i,j)=z(i,j,k)
                end do
            end do

            icount = 1
            do j=0,jy
                do i=0,jx
                    buffer(           icount) = x_pjz(i,j)
                    buffer(  numcount+icount) = y_pjz(i,j)
                    buffer(2*numcount+icount) = z_pjz(i,j)

                    icount = icount + 1
                end do
            end do
        end if

        call MPI_BCAST(buffer(1),3*numcount,MPI_REAL_SD,nproc-1,MPI_COMM_WORLD,ierr)

        icount = 1
        do j=0,jy
            do i=0,jx
                x_pjz(i,j) = buffer(           icount)
                y_pjz(i,j) = buffer(  numcount+icount)
                z_pjz(i,j) = buffer(2*numcount+icount)

                icount = icount + 1
            end do
        end do
        !
        !-----------------------------------------------------------------------
        !
        if(myid==0)write(*,*)'ENTRO'
        !     interpolation
        do j=1,jy
            do k=kparasta,kparaend
                do i=1,jx
                    do isc= 1,nscal

                        xp1 = .5*(x(0 ,j,k)+x(0 ,j,k-1))
                        zp1 = .5*(z(0 ,j,k)+z(0 ,j,k-1))

                        xp2 = .5*(x(jx,j,k)+x(jx,j,k-1))
                        zp2 = .5*(z(jx,j,k)+z(jx,j,k-1))

                        xc = .25*(x(i,j,k)+x(i,j,k-1)+x(i-1,j,k-1)+x(i-1,j,k))
                        zc = .25*(z(i,j,k)+z(i,j,k-1)+z(i-1,j,k-1)+z(i-1,j,k))

                        dist1 = sqrt((xc-xp1)*(xc-xp1) + (zc-zp1)*(zc-zp1))
                        dist2 = sqrt((xc-xp2)*(xc-xp2) + (zc-zp2)*(zc-zp2))
                        dist_tot = dist1 + dist2

                        r_ew = (rhovo1(isc,j,k)*dist2 + rhovo2(isc,j,k)*dist1 )/ dist_tot

                        xp1 = .5*(x_p0(i,j)+x_p0(i-1,j))
                        zp1 = .5*(z_p0(i,j)+z_p0(i-1,j))

                        xp2 = .5*(x_pjz(i,j)+x_pjz(i-1,j))
                        zp2 = .5*(z_pjz(i,j)+z_pjz(i-1,j))

                        dist1 = sqrt((xc-xp1)*(xc-xp1) + (zc-zp1)*(zc-zp1))
                        dist2 = sqrt((xc-xp2)*(xc-xp2) + (zc-zp2)*(zc-zp2))
                        dist_tot = dist1 + dist2

                        r_sn = (rhovo5(isc,i,j)*dist2 + rhovo6(isc,i,j)*dist1 )/ dist_tot

                        rhov(isc,i,j,k) = 0.5*(r_ew+r_sn)

                    !          if(isc ==1 .and. i==10 .and. j==10 .and. myid==0)then
                    !            write(*,*)k,r_sn,r_ew,
                    !     >        rhovo5(isc,i,j),rhov(isc,i,j,k),rhovo6(isc,i,j)
                    !          end if
                    end do
                end do
            end do
        end do
        !
        !-----------------------------------------------------------------------
        !
        return
    end subroutine interp

    subroutine interpola_pareti_pom(ti)
        !***********************************************************************
        use myarrays_velo3

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,ii
        real ti,inv_timewindow
        real delta_l,delta_r
        !-----------------------------------------------------------------------
        !     define some constant values

        inv_timewindow = 1./(ti_pom_new-ti_pom_old)
        delta_l = ti - ti_pom_old
        delta_r = ti_pom_new - ti

        !     side 1
        if(infout1 /=0)then
            i=0
            do k=kparasta-1,kparaend+1 !0,jz+1
                do j=0,jy+1
                    up1(j,k)=(delta_r*uo1(j,k)+delta_l*un1(j,k))*inv_timewindow*index_out1(j,k)

                    vp1(j,k)=(delta_r*vo1(j,k)+delta_l*vn1(j,k))*inv_timewindow*index_out1(j,k)

                    wp1(j,k)=(delta_r*wo1(j,k)+delta_l*wn1(j,k))*inv_timewindow*index_out1(j,k)

                    tke1(j,k)=(delta_r*tkepom1(j,k,ntke-1)+delta_l*tkepom1(j,k,ntke))*inv_timewindow*index_out1(j,k)

                    do ii=1,nscal
                        rhovp1(ii,j,k)=(delta_r*rhovo1(ii,j,k)+delta_l*rhovn1(ii,j,k))*inv_timewindow*index_rho1(j,k)
                    end do

                end do
            end do
        end if

        !     side 2
        if(infout2 /=0)then
            i=jx+1
            do k=kparasta-1,kparaend+1 !0,jz+1
                do j=0,jy+1
                    up2(j,k)=(delta_r*uo2(j,k)+delta_l*un2(j,k))*inv_timewindow*index_out2(j,k)

                    vp2(j,k)=(delta_r*vo2(j,k)+delta_l*vn2(j,k))*inv_timewindow*index_out2(j,k)

                    wp2(j,k)=(delta_r*wo2(j,k)+delta_l*wn2(j,k))*inv_timewindow*index_out2(j,k)

                    tke2(j,k)=(delta_r*tkepom2(j,k,ntke-1)+delta_l*tkepom2(j,k,ntke))*inv_timewindow*index_out2(j,k)

                    do ii=1,nscal
                        rhovp2(ii,j,k)=(delta_r*rhovo2(ii,j,k)+delta_l*rhovn2(ii,j,k))*inv_timewindow*index_rho2(j,k)

                    end do

                end do
            end do
        end if

        !     side 5
        if(infout5 /=0)then
            if(myid==0)then
                k=0
                do j=0,jy+1
                    do i=0,jx+1
                        up5(i,j)=(delta_r*uo5(i,j)+delta_l*un5(i,j))*inv_timewindow*index_out5(i,j)

                        vp5(i,j)=(delta_r*vo5(i,j)+delta_l*vn5(i,j))*inv_timewindow*index_out5(i,j)

                        wp5(i,j)=(delta_r*wo5(i,j)+delta_l*wn5(i,j))*inv_timewindow*index_out5(i,j)

                        tke5(i,j)=(delta_r*tkepom5(i,j,ntke-1)+delta_l*tkepom5(i,j,ntke))*inv_timewindow*index_out5(i,j)

                        do ii=1,nscal
                            rhovp5(ii,i,j)=(delta_r*rhovo5(ii,i,j)+delta_l*rhovn5(ii,i,j))*inv_timewindow*index_rho5(i,j)
                        end do

                    end do
                end do
            end if
        end if

        !     side 6
        if(infout6 /=0)then
            if(myid==nproc-1)then
                k=jz+1
                do j=0,jy+1
                    do i=0,jx+1
                        up6(i,j)=(delta_r*uo6(i,j)+delta_l*un6(i,j))*inv_timewindow*index_out6(i,j)

                        vp6(i,j)=(delta_r*vo6(i,j)+delta_l*vn6(i,j))*inv_timewindow*index_out6(i,j)

                        wp6(i,j)=(delta_r*wo6(i,j)+delta_l*wn6(i,j))*inv_timewindow*index_out6(i,j)

                        tke6(i,j)=(delta_r*tkepom6(i,j,ntke-1)+delta_l*tkepom6(i,j,ntke))*inv_timewindow*index_out6(i,j)

                        do ii=1,nscal
                            rhovp6(ii,i,j)=(delta_r*rhovo6(ii,i,j)+delta_l*rhovn6(ii,i,j))*inv_timewindow*index_rho6(i,j)
                        end do

                    end do
                end do
            end if
        end if

        return

    end subroutine interpola_pareti_pom

    subroutine contrin_pot()
        !***********************************************************************
        ! compute controvariant fluxes from interpolated velocity field
        ! this subroutine is used for nesting
        !
        use myarrays_velo3
        use myarrays_metri3

        implicit none
        !

        !
        !-----------------------------------------------------------------------
        ! array declaration
        integer i,j,k
        real diver,divmax
        !
        !-----------------------------------------------------------------------
        !
        ! flux J-1*U
        !
        !     side 1
        if(infout1 /= 0)then
            do k=kparasta,kparaend !1,jz
                do j=1,jy

                    uc(0,j,k)=csx(0,j,k)*uo1(j,k)+csy(0,j,k)*vo1(j,k)+csz(0,j,k)*wo1(j,k)
                    ucp1(j,k)=uc(0,j,k)

                enddo
            enddo
        end if

        !     side 2
        if(infout2 /= 0)then
            do k=kparasta,kparaend !1,jz
                do j=1,jy

                    uc(jx,j,k)=csx(jx,j,k)*uo2(j,k)+csy(jx,j,k)*vo2(j,k)+csz(jx,j,k)*wo2(j,k)
                    ucp2(j,k)=uc(jx,j,k)

                enddo
            enddo
        end if

        !
        ! flux J-1*V
        vc(:,0,:)=0.
        vc(:,jy,:)=0.
        !
        !
        ! flux J-1*W
        !
        !     side 5
        if(infout5 /= 0)then
            if(myid.eq.0)then

                do j=1,jy
                    do i=1,jx
                        wc(i,j,0)=ztx(i,j,0)*uo5(i,j)+zty(i,j,0)*vo5(i,j)+ztz(i,j,0)*wo5(i,j)
                        wcp5(i,j)=wc(i,j,0)
                    enddo
                enddo

            end if
        end if

        !     side 6
        if(infout6 /= 0)then
            if(myid.eq.nproc-1)then

                do j=1,jy
                    do i=1,jx
                        wc(i,j,jz)=ztx(i,j,jz)*uo6(i,j)+zty(i,j,jz)*vo6(i,j)+ztz(i,j,jz)*wo6(i,j)
                        wcp6(i,j)=wc(i,j,jz)
                    enddo
                enddo

            end if
        end if
        !
        return
    end

    subroutine contrin_lat()
        !***********************************************************************
        ! compute controvariant fluxes from interpolated velocity field
        ! this subroutine is used for nesting
        !
        use myarrays_velo3
        use myarrays_metri3

        implicit none
        !
        !
        !-----------------------------------------------------------------------
        ! array declaration
        integer i,j,k
        real diver,divmax
        !
        !-----------------------------------------------------------------------
        !
        ! flux J-1*U
        !
        !     side 1
        if(infout1 /= 0)then
            do k=kparasta,kparaend !1,jz
                do j=1,jy

                    uc(0,j,k)=csx(0,j,k)*up1(j,k)+csy(0,j,k)*vp1(j,k)+csz(0,j,k)*wp1(j,k)
                    ucp1(j,k)=uc(0,j,k)

                enddo
            enddo
        end if
        !    side 2
        if(infout2 /= 0)then
            do k=kparasta,kparaend !1,jz
                do j=1,jy

                    uc(jx,j,k)=csx(jx,j,k)*up2(j,k)+csy(jx,j,k)*vp2(j,k)+csz(jx,j,k)*wp2(j,k)
                    ucp2(j,k)=uc(jx,j,k)

                enddo
            enddo
        end if
        !
        ! flux J-1*V
        vc(:,0,:)=0.
        vc(:,jy,:)=0.
        !
        !
        ! flux J-1*W
        !
        !     side 5
        if(infout5 /= 0)then
            if(myid.eq.0)then

                do j=1,jy
                    do i=1,jx
                        wc(i,j,0)=ztx(i,j,0)*up5(i,j)+zty(i,j,0)*vp5(i,j)+ztz(i,j,0)*wp5(i,j)
                        wcp5(i,j)=wc(i,j,0)
                    enddo
                enddo

            end if
        end if

        !     side 6
        if(infout6 /= 0)then
            if(myid.eq.nproc-1)then

                do j=1,jy
                    do i=1,jx
                        wc(i,j,jz)=ztx(i,j,jz)*up6(i,j)+zty(i,j,jz)*vp6(i,j)+ztz(i,j,jz)*wp6(i,j)
                        wcp6(i,j)=wc(i,j,jz)
                    enddo
                enddo

            end if
        end if
        !
        return
    end

    subroutine leggi_pareti_pom(ti)
        !***********************************************************************
        !     read input data from pom simulation, used in nesting procedure
        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,n,isc

        real ti
        real yyy

        real var1,var2,var3,var4,var(nscal)

        logical read_plane
        !-----------------------------------------------------------------------

        read_plane = .true.

        do while(read_plane)
            if(myid==0)write(*,*)'READ PLANE NESTING'

            ntke = ntke + 1

            ! read side 1
            if(infout1 /= 0)then

                read(  81,*)ti_pom_new

                do k=0,jz+1
                    do j=0,jy+1
                        read(  81,*)var1,var2,var3,var4,(var(n),n=1,nscal)

                        if(k.ge.kparasta.and.k.le.kparaend)then
                            un1(j,k) = var1*index_out1(j,k)
                            vn1(j,k) = var2*index_out1(j,k)
                            wn1(j,k) = var3*index_out1(j,k)
                            tkepom1(j,k,ntke)=var4*index_out1(j,k)
                            do isc = 1,nscal
                                rhovn1(isc,j,k) = var(isc)*index_rho1(j,k)
                            end do
                        endif
                    end do
                end do

            end if

            ! read side 2
            if(infout2 /= 0)then

                read(  82,*)ti_pom_new

                do k=0,jz+1
                    do j=0,jy+1
                        read(   82,*)var1,var2,var3,var4,(var(n),n=1,nscal)
                        if(k.ge.kparasta.and.k.le.kparaend)then
                            un2(j,k) = var1*index_out2(j,k)
                            vn2(j,k) = var2*index_out2(j,k)
                            wn2(j,k) = var3*index_out2(j,k)
                            tkepom2(j,k,ntke)=var4*index_out2(j,k)
                            do isc = 1,nscal
                                rhovn2(isc,j,k) = var(isc)*index_rho2(j,k)
                            end do
                        endif
                    enddo
                enddo

            end if

            ! read side 5
            if(infout5 /= 0)then
                read(  85,*)ti_pom_new

                do j=1,jy
                    do i=1,jx

                        read( 85,*)var1,var2,var3,var4,(var(n),n=1,nscal)

                        un5(i,j) = var1*index_out5(i,j)
                        vn5(i,j) = var2*index_out5(i,j)
                        wn5(i,j) = var3*index_out5(i,j)
                        tkepom5(i,j,ntke)=var4*index_out5(i,j)
                        do isc = 1,nscal
                            rhovn5(isc,i,j) = var(isc)*index_rho5(i,j)
                        end do
                    end do
                end do

            end if

            ! read side 6
            if(infout6 /= 0)then
                read(  86,*)ti_pom_new

                do j=1,jy
                    do i=1,jx

                        read( 86,*)var1,var2,var3,var4,(var(n),n=1,nscal)

                        un6(i,j) = var1*index_out6(i,j)
                        vn6(i,j) = var2*index_out6(i,j)
                        wn6(i,j) = var3*index_out6(i,j)
                        tkepom6(i,j,ntke)=var4*index_out6(i,j)
                        do isc = 1,nscal
                            rhovn6(isc,i,j) = var(isc)*index_rho6(i,j)
                        end do
                    end do
                end do

            end if

            !cc         ti_pom_new = ti_pom_new * 3600.


            if(ti > ti_pom_new)then
                if(infout1 /= 0)then
                    uo1=un1
                    vo1=vn1
                    wo1=wn1
                    rhovo1=rhovn1
                end if

                if(infout2 /= 0)then
                    uo2=un2
                    vo2=vn2
                    wo2=wn2
                    rhovo2=rhovn2
                end if

                if(infout5 /= 0)then
                    uo5=un5
                    vo5=vn5
                    wo5=wn5
                    rhovo5=rhovn5
                end if

                if(infout6 /= 0)then
                    uo6=un6
                    vo6=vn6
                    wo6=wn6
                    rhovo6=rhovn6
                end if

                ti_pom_old=ti_pom_new
            else

                read_plane = .false.

            end if


        end do !do while

8753    format(6e14.8)

        return
    end

!***********************************************************************
end module nesting_module
!***********************************************************************
