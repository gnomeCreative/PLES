module inflow_module

    use orlansky_module
    use scala3
    use mysending, only: kparasta,kparaend,rightpe,leftpe,nproc,myid
    !
    use mpi

    use contour_module

    implicit none

    real :: area_bagnata1,area_bagnata2
    real :: area_bagnata3,area_bagnata4
    real :: area_bagnata5,area_bagnata6

    real :: area1,area2,area3,area4,area5,area6

    real,allocatable :: areola1(:,:),areola2(:,:)
    real,allocatable :: areola3(:,:),areola4(:,:)
    real,allocatable :: areola5(:,:),areola6(:,:)
    real,allocatable :: gg221(:,:),gg331(:,:),gg231(:,:)
    real,allocatable :: gg112(:,:),gg222(:,:),gg122(:,:)
    real,allocatable :: gg113(:,:),gg133(:,:),gg333(:,:)

contains

    subroutine intialize_myarrays2()
        ! compute total area for each sides and cell

        implicit none

        allocate(areola1(0:n2+1,kparasta-1:kparaend+1))
        allocate(areola2(0:n2+1,kparasta-1:kparaend+1))
        allocate(areola3(0:n1+1,kparasta-1:kparaend+1))
        allocate(areola4(0:n1+1,kparasta-1:kparaend+1))
        allocate(areola5(0:n1+1,0:n2+1))
        allocate(areola6(0:n1+1,0:n2+1))
        !
        allocate(gg221(0:n2+1,kparasta-1:kparaend+1))
        allocate(gg331(0:n2+1,kparasta-1:kparaend+1))
        allocate(gg231(0:n2+1,kparasta-1:kparaend+1))
        allocate(gg112(0:n1+1,0:n2+1))
        allocate(gg222(0:n1+1,0:n2+1))
        allocate(gg122(0:n1+1,0:n2+1))
        allocate(gg113(0:n1+1,kparasta-1:kparaend+1))
        allocate(gg133(0:n1+1,kparasta-1:kparaend+1))
        allocate(gg333(0:n1+1,kparasta-1:kparaend+1))

    end subroutine intialize_myarrays2

    subroutine read_inflow(ktime)

        use mysettings, only:lett,niter

        implicit none

        integer,intent(in) :: ktime
        integer :: i,j,k,isc
        integer :: pianoturbo

        if (lett) then
            !     side 1
            if (infout1==0) then
                if (npn1 > 1 .and. ktime /= niter) then
                    read(10011,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 1'
                    end if
                    do k=1,jz
                        do j=1,jy
                            read(10011,100)up1(j,k),vp1(j,k),wp1(j,k),(rhovp1(isc,j,k),isc=1,nscal)

                            up1(j,k)=up1(j,k)*real(index_out1(j,k))
                            vp1(j,k)=vp1(j,k)*real(index_out1(j,k))
                            wp1(j,k)=wp1(j,k)*real(index_out1(j,k))
                            do isc=1,nscal
                                rhovp1(isc,j,k)=rhovp1(isc,j,k)*real(index_out1(j,k))
                            end do
                        end do
                    end do
                elseif (ktime == niter) then
                    close(10011)
                end if
            end if
            !     side 2
            if (infout2==0) then
                if (npn2 > 1 .and. ktime /= niter) then
                    read(10012,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 2'
                    end if
                    do k=1,jz
                        do j=1,jy
                            read(10012,100)up2(j,k),vp2(j,k),wp2(j,k),(rhovp2(isc,j,k),isc=1,nscal)

                            up2(j,k)=up2(j,k)*real(index_out2(j,k))
                            vp2(j,k)=vp2(j,k)*real(index_out2(j,k))
                            wp2(j,k)=wp2(j,k)*real(index_out2(j,k))
                            do isc=1,nscal
                                rhovp2(isc,j,k)=rhovp2(isc,j,k)*real(index_out2(j,k))
                            end do
                        end do
                    end do
                elseif (ktime == niter) then
                    close(10012)
                end if
            end if
            !     side 3
            if (infout3==0) then
                if (npn3 > 1 .and. ktime /= niter) then
                    read(10013,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 3'
                    end if
                    do k=1,jz
                        do i=1,jx
                            read(10013,100)up3(i,k),vp3(i,k),wp3(i,k),(rhovp3(isc,i,k),isc=1,nscal)

                            up3(i,k)=up3(i,k)*real(index_out3(i,k))
                            vp3(i,k)=vp3(i,k)*real(index_out3(i,k))
                            wp3(i,k)=wp3(i,k)*real(index_out3(i,k))
                            do isc=1,nscal
                                rhovp3(isc,i,k)=rhovp3(isc,i,k)*real(index_out3(i,k))
                            end do

                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10013)
                end if
            end if
            !     side 4
            if (infout4==0) then
                if (npn4 > 1 .and. ktime /= niter) then
                    read(10014,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 4'
                    end if
                    do k=1,jz
                        do i=1,jx
                            read(10014,100)up4(i,k),vp4(i,k),wp4(i,k),(rhovp4(isc,i,k),isc=1,nscal)

                            up4(i,k)=up4(i,k)*real(index_out4(i,k))
                            vp4(i,k)=vp4(i,k)*real(index_out4(i,k))
                            wp4(i,k)=wp4(i,k)*real(index_out4(i,k))
                            do isc=1,nscal
                                rhovp4(isc,i,k)=rhovp4(isc,i,k)*real(index_out4(i,k))
                            end do

                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10014)
                end if
            end if
            !     side 5
            if (infout5==0) then
                if (npn5 > 1 .and. ktime /= niter) then
                    read(10015,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 5'
                    end if
                    do j=1,jy
                        do i=1,jx
                            read(10015,100)up5(i,j),vp5(i,j),wp5(i,j),(rhovp5(isc,i,j),isc=1,nscal)

                            up5(i,j)=up5(i,j)*real(index_out5(i,j))
                            vp5(i,j)=vp5(i,j)*real(index_out5(i,j))
                            wp5(i,j)=wp5(i,j)*real(index_out5(i,j))
                            do isc=1,nscal
                                rhovp5(isc,i,j)=rhovp5(isc,i,j)*real(index_out5(i,j))
                            end do

                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10015)
                end if
            end if
            !     side 6
            if (infout6==0) then
                if (npn6 > 1 .and. ktime /= niter) then
                    read(10016,*)pianoturbo
                    if (myid==0) then
                        write(*,*)'read plane: ',pianoturbo,' side 6'
                    end if
                    do j=1,jy
                        do i=1,jx
                            read(10016,100)up6(i,j),vp6(i,j),wp6(i,j),(rhovp6(isc,i,j),isc=1,nscal)

                            up6(i,j)=up6(i,j)*real(index_out6(i,j))
                            vp6(i,j)=vp6(i,j)*real(index_out6(i,j))
                            wp6(i,j)=wp6(i,j)*real(index_out6(i,j))
                            do isc=1,nscal
                                rhovp6(isc,i,j)=rhovp6(isc,i,j)*real(index_out6(i,j))
                            end do
                        enddo
                    enddo
                elseif (ktime == niter) then
                    close(10016)
                end if
            end if

100         format(10e13.5)

        end if

    end subroutine read_inflow

    subroutine inflow()

        use myarrays_velo3

        implicit none

        !-----------------------------------------------------------------------
        ! array declaration
        character*20 filepn
        integer i,j,k,isc,npn,ierr
        integer l,numero_celle_bloccate
        integer nnx,nny,nnz


        real parziale1,parziale2
        real parziale3,parziale4
        real parziale5,parziale6

        integer cont1,cont2,cont3,cont4,cont5,cont6

        real disct1A(jy,jz)
        real disct1B,disct1,disct2,a1,a2,a3
        real xib,yib,zib,xv,yv,zv,fattore
        integer kk,pianoturbo,npiani

        real massa,massatot
      
        real val_tke1,val_tke2
        real val_tke5,val_tke6
        !-----------------------------------------------------------------------
        !     read the tke

        !     side 1
        if(ibodybuffer1)then
            open(121,file='tke1.dat',status='unknown')
            do k=1,jz
                do j=1,jy
                    read(121,*)val_tke1
                    if(k.ge. kparasta-1 .and. k.le. kparaend+1)then
                        tke1(j,k) = val_tke1
                    end if
                end do
            end do
            close(121)
        end if

        !     side 2
        if(ibodybuffer2)then
            open(121,file='tke2.dat',status='unknown')
            do k=1,jz
                do j=1,jy
                    read(122,*)val_tke2
                    if(k.ge. kparasta-1 .and. k.le. kparaend+1)then
                        tke2(j,k) = val_tke2
                    end if
                end do
            end do
            close(122)
        end if

        !     side 5
        if(myid==0)then
            if(ibodybuffer5)then
                open(125,file='tke5.dat',status='unknown')
                do j=1,jy
                    do i=1,jx
                        read(125,*)val_tke5
                        tke5(i,j) = val_tke5
                    end do
                end do
                close(125)
            end if
        end if

        !     side 6
        if(myid==nproc-1)then
            if(ibodybuffer6)then
                open(126,file='tke6.dat',status='unknown')
                do j=1,jy
                    do i=1,jx
                        read(126,*)val_tke6
                        tke6(i,j) = val_tke6
                    end do
                end do
                close(126)
            end if
        end if
      
        !-----------------------------------------------------------------------
        !    OUTFLOW
        !-----------------------------------------------------------------------
        ! in case of partial outflow
        ! if no outflow index = 0
        ! if    outflow index = 1
        !
        !-----------------------------------------------------------------------
        ! inflow

        !     planes at constant csi
        if(infout1.eq.0)then
	
            filepn='piano1.dat'
            open(10011,file=filepn,status='unknown')
            read(10011,*)nny,nnz !jy,jz
            read(10011,*)npn1
            read(10011,*)npn

            do k=1,jz
                do j=1,jy
                    read(10011,100)up1(j,k),&
                        vp1(j,k),&
                        wp1(j,k),&
                        (rhovp1(isc,j,k),isc=1,nscal)
     
                    up1(j,k)=up1(j,k)*real(index_out1(j,k))
                    vp1(j,k)=vp1(j,k)*real(index_out1(j,k))
                    wp1(j,k)=wp1(j,k)*real(index_out1(j,k))
                    do isc=1,nscal
                        rhovp1(isc,j,k)=rhovp1(isc,j,k)*real(index_out1(j,k))
                    end do
                end do
            end do
	

	
	
            if(myid.eq.0) then
                print*,myid,'inflow plane face 1'
            endif
        endif
100     format(10e13.5)

        if(infout2.eq.0)then
	
            filepn='piano2.dat'
            open(10012,file=filepn,status='unknown')
            read(10012,*)nny,nnz !jy,jz
            read(10012,*)npn2
            read(10012,*)npn

            do k=1,jz
                do j=1,jy
                    read(10012,101)up2(j,k),&
                        vp2(j,k),&
                        wp2(j,k),&
                        (rhovp2(isc,j,k),isc=1,nscal)
     
                    up2(j,k)=up2(j,k)*real(index_out2(j,k))
                    vp2(j,k)=vp2(j,k)*real(index_out2(j,k))
                    wp2(j,k)=wp2(j,k)*real(index_out2(j,k))
                    do isc=1,nscal
                        rhovp2(isc,j,k)=rhovp2(isc,j,k)*real(index_out2(j,k))
                    end do
                enddo
            enddo
            if(myid.eq.0) then
                print*,myid,'inflow plane face 2'
            endif
        endif
101     format(10e13.5)

        !     planes at constant eta
        if(infout3.eq.0)then

            filepn='piano3.dat'
            open(10013,file=filepn,status='unknown')
            read(10013,*)nnx,nnz !jx,jz
            read(10013,*)npn3
            read(10013,*)npn
            do k=1,jz
                do i=1,jx
                    read(10013,102)up3(i,k),&
                        vp3(i,k),&
                        wp3(i,k),&
                        (rhovp3(isc,i,k),isc=1,nscal)

                    up3(i,k)=up3(i,k)*real(index_out3(i,k))
                    vp3(i,k)=vp3(i,k)*real(index_out3(i,k))
                    wp3(i,k)=wp3(i,k)*real(index_out3(i,k))
                    do isc=1,nscal
                        rhovp3(isc,i,k)=rhovp3(isc,i,k)*real(index_out3(i,k))
                    end do
                enddo
            enddo
            if(myid.eq.0) then
                print*,myid,'inflow plane face 3'
            endif
        endif
102     format(10e13.5)
 
        if(infout4.eq.0)then

            filepn='piano4.dat'
            open(10014,file=filepn,status='unknown')
            read(10014,*)nnx,nnz !jx,jz
            read(10014,*)npn4
            read(10014,*)npn
	
            do k=1,jz
                do i=1,jx
                    read(10014,103)up4(i,k),&
                        vp4(i,k),&
                        wp4(i,k),&
                        (rhovp4(isc,i,k),isc=1,nscal)

                    up4(i,k)=up4(i,k)*real(index_out4(i,k))
                    vp4(i,k)=vp4(i,k)*real(index_out4(i,k))
                    wp4(i,k)=wp4(i,k)*real(index_out4(i,k))
                    do isc=1,nscal
                        rhovp4(isc,i,k)=rhovp4(isc,i,k)*real(index_out4(i,k))
                    end do
                enddo
            enddo
            if(myid.eq.0) then
                print*,myid,'inflow plane face 4'
            endif
        endif
103     format(10e13.5)

        !     planes at constant eta
        if(infout5.eq.0)then

            filepn='piano5.dat'
            open(10015,file=filepn,status='unknown')
            read(10015,*)nnx,nny !jx,jy
            read(10015,*)npn5
            read(10015,*)npn
	
            do j=1,jy
                do i=1,jx
                    read(10015,104)up5(i,j),&
                        vp5(i,j),&
                        wp5(i,j),&
                        (rhovp5(isc,i,j),isc=1,nscal)

                    up5(i,j)=up5(i,j)*real(index_out5(i,j))
                    vp5(i,j)=vp5(i,j)*real(index_out5(i,j))
                    wp5(i,j)=wp5(i,j)*real(index_out5(i,j))
                    do isc=1,nscal
                        rhovp5(isc,i,j)=rhovp5(isc,i,j)*real(index_out5(i,j))
                    end do
                    if (myid==nproc-1) then
                        write(*,*)i*j,up5(i,j),&
                            vp5(i,j),&
                            wp5(i,j),&
                            (rhovp5(isc,i,j),isc=1,nscal)
                    end if
                enddo
            enddo
            if(myid.eq.0) then
                print*,myid,' inflow plane face 5'
            endif
        endif
104     format(10e13.5)

        if(infout6.eq.0)then

            filepn='piano6.dat'
            open(10016,file=filepn,status='unknown')
            read(10016,*)nnx,nny !jx,jy
            read(10016,*)npn6
            read(10016,*)npn
	
            do j=1,jy
                do i=1,jx
                    read(10016,105)up6(i,j),&
                        vp6(i,j),&
                        wp6(i,j),&
                        (rhovp6(isc,i,j),isc=1,nscal)

                    up6(i,j)=up6(i,j)*real(index_out6(i,j))
                    vp6(i,j)=vp6(i,j)*real(index_out6(i,j))
                    wp6(i,j)=wp6(i,j)*real(index_out6(i,j))
                    do isc=1,nscal
                        rhovp6(isc,i,j)=rhovp6(isc,i,j)*real(index_out6(i,j))
                    end do
                enddo
            enddo
            if(myid.eq.0) then
                print*,myid,'inflow plane face 6'
            endif
        endif
105     format(10e13.5)
        !-----------------------------------------------------------------------
        return
    end subroutine inflow

    subroutine contra_infout(ktime)

        !***********************************************************************
        ! compute Uc* intermediate controvariant quantities in case of
        ! inflow/outflow and boundary condition for pressure cs1 etc.
        !
        use myarrays_metri3
        use myarrays_velo3
        use output_module, only: info_run_file
        !
        use period


        implicit none


        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,ii,jj,kk,iii,ktime
        real uinter,vinter,winter

        real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddi,dsi,ddj,dsj,ddk,dsk
        !
        integer ierr
        integer ncolperproc,m
        integer leftpem,rightpem
        integer kparastal,kparaendl
        integer status(MPI_STATUS_SIZE)

        real :: am_in,am_in_loc
        real :: am_out,am_out_loc
        real :: am_out_loc1,am_out_loc2,am_out_loc3
        real :: am_out_loc4,am_out_loc5,am_out_loc6
        real :: am_out1,am_out2,am_out3
        real :: am_out4,am_out5,am_out6
        real :: am_out_ripartisco

        real :: am1,am2,am1t,am2t
        real :: am3,am4,am3t,am4t
        real :: am5,am6,am5t,am6t

        real :: giactot,giactot_loc
        real :: adel,adel_loc
        real :: del_mas

        real :: coef_massa1,coef_massa2,coef_massa3
        real :: coef_massa4,coef_massa5,coef_massa6

        real :: c_massa1(n2,n3),c_massa2(n2,n3)
        real :: c_massa3(n1,n3),c_massa4(n1,n3)
        real :: c_massa5(n1,n2),c_massa6(n1,n2)

        real :: contatore

        real area_bagnata
        real :: somma,sommat,amm,massa_tot

        integer faccio

        real ucc,vcc,wcc
        real area_insisto,fattore_flusso
        !
        !---------------------------------------------------------------------
        !
        contatore=0.
        if(infout1.eq.1)then
            contatore=contatore+1.
        endif
        if(infout2.eq.1)then
            contatore=contatore+1.
        endif
        if(infout3.eq.1)then
            contatore=contatore+1.
        endif
        if(infout4.eq.1)then
            contatore=contatore+1.
        endif
        if(infout5.eq.1)then
            contatore=contatore+1.
        endif
        if(infout6.eq.1)then
            contatore=contatore+1.
        endif
        if(myid.eq.0)then
            write(*,*)'INFLOW/OUTFLOW'
            write(info_run_file,*)'contatore pareti orlansky',contatore
        endif

        !-----------------------------------------------------------------------
        ! COMPUTE INCOMING MASS
        !-----------------------------------------------------------------------
        am_in=0.
        am_in_loc=0.
        del_mas=0.

        !     sides 1 and 2 constant csi
        do ii=1,ip

            if(infout1.ne.1)then
                do k=kparasta,kparaend
                    do j=1,jy
                        am_in_loc=am_in_loc+uc(0,j,k)
                    end do
                end do
            endif

            if(infout2.ne.1)then
                do k=kparasta,kparaend
                    do j=1,jy
                        am_in_loc=am_in_loc-uc(jx,j,k)
                    end do
                end do
            endif

        end do
        !
        !     sides 3 and 4 constant eta
        do jj=1,jp

            if(infout3.ne.1)then
                do k=kparasta,kparaend
                    do i=1,jx
                        am_in_loc=am_in_loc+vc(i,0,k)
                    end do
                end do
            endif

            if(infout4.ne.1)then
                do k=kparasta,kparaend
                    do i=1,jx
                        am_in_loc=am_in_loc-vc(i,jy,k)
                    end do
                end do
            endif

        end do
        !
        !     side 5 and 6 constant zita
        do kk=1,kp

            if(myid.eq.0)then
                if(infout5.ne.1)then
                    do j=1,jy
                        do i=1,jx
                            am_in_loc=am_in_loc+wc(i,j,0)
                        end do
                    end do
                endif

            elseif(myid.eq.nproc-1)then

                if(infout6.ne.1)then
                    do j=1,jy
                        do i=1,jx
                            am_in_loc=am_in_loc-wc(i,j,jz)
                        end do
                    end do
                endif
            endif

        end do
        !
        ! make the incoming mass known to all procs
        call MPI_ALLREDUCE(am_in_loc,am_in,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        !-----------------------------------------------------------------------
        ! COMPUTE OUTFLOW MASS
        !-----------------------------------------------------------------------
        ! index initialization

        am_out_loc=0.
        am_out=0.

        am_out_loc1=0.
        am_out_loc2=0.
        am_out_loc3=0.
        am_out_loc4=0.
        am_out_loc5=0.
        am_out_loc6=0.

        am_out1=0.
        am_out2=0.
        am_out3=0.
        am_out4=0.
        am_out5=0.
        am_out6=0.

        am1=0.
        am1t=0.

        am2=0.
        am2t=0.

        am3=0.
        am3t=0.

        am4=0.
        am4t=0.

        am5=0.
        am5t=0.

        am6=0.
        am6t=0.
        !
        !     side 1 and 2 constant csi
        do ii=1,ip

            if(infout1.eq.1)then
                do k=kparasta,kparaend
                    do j=1,jy

                        uc(0,j,k)=du_dx1(j,k)*csx(0,j,k)+dv_dx1(j,k)*csy(0,j,k)+dw_dx1(j,k)*csz(0,j,k) !U^(n+1)  with orlansky

                        am_out_loc1=am_out_loc1-uc(0,j,k)   ! outflow
                        am1=am1+abs(uc(0,j,k))              ! module outflow
                    end do
                end do
            endif

            if(infout2.eq.1)then
                do k=kparasta,kparaend
                    do j=1,jy

                        uc(jx,j,k)=du_dx2(j,k)*csx(jx,j,k)+dv_dx2(j,k)*csy(jx,j,k)+dw_dx2(j,k)*csz(jx,j,k)

                        am_out_loc2=am_out_loc2+uc(jx,j,k)
                        am2=am2+abs(uc(jx,j,k))
                    end do
                end do
            endif

        end do
        !
        !     sides 3 and 4 constant eta
        do jj=1,jp

            if(infout3.eq.1)then
                do k=kparasta,kparaend
                    do i=1,jx

                        vc(i,0,k)=du_dy3(i,k)*etx(i,0,k)+dv_dy3(i,k)*ety(i,0,k)+dw_dy3(i,k)*etz(i,0,k)

                        am_out_loc3=am_out_loc3-vc(i,0,k)
                        am3=am3+abs(vc(i,0,k))
                    end do
                end do
            endif

            if(infout4.eq.1)then
                do k=kparasta,kparaend
                    do i=1,jx

                        vc(i,jy,k)=du_dy4(i,k)*etx(i,jy,k)+dv_dy4(i,k)*ety(i,jy,k)+dw_dy4(i,k)*etz(i,jy,k)

                        am_out_loc4=am_out_loc4+vc(i,jy,k)
                        am4=am4+abs(vc(i,jy,k))
                    end do
                end do
            endif

        end do
        !
        !     sides 5 and 6 constant zita
        do kk=1,kp

            if(myid.eq.0)then

                if(infout5.eq.1)then
                    do j=1,jy
                        do i=1,jx

                            wc(i,j,0)=du_dz5(i,j)*ztx(i,j,0)+dv_dz5(i,j)*zty(i,j,0)+dw_dz5(i,j)*ztz(i,j,0)

                            am_out_loc5=am_out_loc5-wc(i,j,0)
                            am5=am5+abs(wc(i,j,0))
                        end do
                    end do
                endif

            elseif(myid.eq.nproc-1)then

                if(infout6.eq.1)then
                    do j=1,jy
                        do i=1,jx

                            wc(i,j,jz)=du_dz6(i,j)*ztx(i,j,jz)+dv_dz6(i,j)*zty(i,j,jz)+dw_dz6(i,j)*ztz(i,j,jz)

                            am_out_loc6=am_out_loc6+wc(i,j,jz)
                            am6=am6+abs(wc(i,j,jz))
                        end do
                    end do
                endif

            endif

        end do

        am_out_loc=am_out_loc1+am_out_loc2+am_out_loc3+am_out_loc4+am_out_loc5+am_out_loc6

        !
        ! make outflow known to all procs
        call MPI_ALLREDUCE(am_out_loc,am_out,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        !-----------------------------------------------------------------------
        ! COMPUTE MASS FRACTION
        !-----------------------------------------------------------------------
        !
        ! mass fraction for each sides, for distribution on cell area
        !
        !
        ! mass for each face

        call MPI_ALLREDUCE(am_out_loc1,am_out1,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc2,am_out2,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc3,am_out3,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc4,am_out4,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc5,am_out5,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc6,am_out6,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! mass for each face in modulus

        call MPI_ALLREDUCE(am1,am1t,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am2,am2t,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am3,am3t,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am4,am4t,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am5,am5t,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am6,am6t,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        amm = am1t + am2t + am3t + am4t + am5t + am6t
        !
        ! compute area with outflow
        area_bagnata=0.
        if(infout1.eq.1)then
            area_bagnata = area_bagnata+area_bagnata1
        end if
        if(infout2.eq.1)then
            area_bagnata = area_bagnata+area_bagnata2
        end if
        if(infout3.eq.1)then
            area_bagnata = area_bagnata+area_bagnata3
        end if
        if(infout4.eq.1)then
            area_bagnata = area_bagnata+area_bagnata4
        end if
        if(infout5.eq.1)then
            area_bagnata = area_bagnata+area_bagnata5
        end if
        if(infout6.eq.1)then
            area_bagnata = area_bagnata+area_bagnata6
        end if

        if(myid.eq.0)then
            write(*,*)'area_bagnata totale',area_bagnata
        end if
        !
        !-----------------------------------------------------------------------
        ! MASS DEFECT
        !-----------------------------------------------------------------------
        !
        del_mas=am_out-am_in
        !
        !-----------------------------------------------------------------------
        ! OUTPUT
        !-----------------------------------------------------------------------
        !
        if(myid.eq.0)then
            write(*,*)'inflow',am_in
            write(*,*)'outflow',am_out
            write(*,*)'delta',del_mas
            write(info_run_file,*)'massa attraverso faccia 1:',am_out1
            write(info_run_file,*)'massa attraverso faccia 2:',am_out2
            write(info_run_file,*)'massa attraverso faccia 3:',am_out3
            write(info_run_file,*)'massa attraverso faccia 4:',am_out4
            write(info_run_file,*)'massa attraverso faccia 5:',am_out5
            write(info_run_file,*)'massa attraverso faccia 6:',am_out6
        endif
        !
        !-----------------------------------------------------------------------
        ! MASS CORRECTION ON CONTROVARIANT FLUX
        !-----------------------------------------------------------------------
        !
        somma = 0.

        massa_tot = abs(am_out1)+abs(am_out2)+abs(am_out3)+abs(am_out4)+abs(am_out5)+abs(am_out6)


        do ii=1,ip

            !     side 1 constant csi
            if(massa_tot .lt. 1.d-8)then
                area_insisto   =   area_bagnata
                fattore_flusso =   1.
            else
                area_insisto   =   area_bagnata1
                fattore_flusso =  (abs(am_out1)/massa_tot)
            end if


            if(infout1 == 1)then

                do k=kparasta,kparaend
                    do j=1,jy
                        !
                        coef_massa1=(areola1(j,k)/area_insisto)*fattore_flusso

                        somma=somma+coef_massa1*index_out1(j,k)
                        !
                        !          Uc^n+1
                        ucc=uc(0,j,k)+del_mas*coef_massa1*index_out1(j,k)

                        uc1_orl(j,k)=ucc
                        !          pressure b.c. = Uc* - Uc^n+1
                        cs1(j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)-ucc
                        if(potenziale)then
                            cs1(j,k)=-ucc
                        end if
                        !       Uc*
                        uc(0,j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)

                    enddo
                enddo

            else

                do k=kparasta,kparaend
                    do j=1,jy
                        !          Uc stored to compute the divg
                        uc1_orl(j,k)=uc(0,j,k)

                        cs1(j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)-uc(0,j,k)
                        !          pressure b.c. = Uc* - Uc^n+1
                        if(potenziale)then
                            cs1(j,k)= -uc(0,j,k)
                        end if
                        !       Uc*
                        uc(0,j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)

                    !
                    enddo
                enddo

            end if !infout1
            !     ..................................................................
            !     side 2 constant csi
            if(massa_tot .lt. 1.d-8)then
                area_insisto   =   area_bagnata
                fattore_flusso =   1.
            else
                area_insisto   =   area_bagnata2
                fattore_flusso =  (abs(am_out2)/massa_tot)
            end if

            if(infout2 == 1)then

                do k=kparasta,kparaend
                    do j=1,jy
                        !
                        coef_massa2=(areola2(j,k)/area_insisto)*fattore_flusso
                        !
                        somma=somma+coef_massa2*index_out2(j,k)

                        !          Uc^n+1
                        ucc=uc(jx,j,k)-del_mas*coef_massa2*index_out2(j,k)
                        !
                        uc2_orl(j,k)=ucc

                        !          pressure b.c. = Uc* - Uc^n+1
                        cs2(j,k)=u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)-ucc
                        if(potenziale)then
                            cs2(j,k)=-ucc
                        end if
                        !       Uc*
                        uc(jx,j,k) = u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)

                    enddo
                enddo

            else

                do k=kparasta,kparaend
                    do j=1,jy
                        !          Uc stored to compute the divg
                        uc2_orl(j,k)=uc(jx,j,k)
                        !          pressure b.c. = Uc* - Uc^n+1
                        cs2(j,k)=u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)-uc(jx,j,k)

                        if(potenziale)then
                            cs2(j,k)=-uc(jx,j,k)
                        end if
                        !       Uc*
                        uc(jx,j,k)=u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)

                    enddo
                enddo

            end if !infout2

        enddo   !end loop ii=1,ip
        !.......................................................................

        do jj=1,jp

            if(massa_tot .lt. 1.d-8)then
                area_insisto   =   area_bagnata
                fattore_flusso =   1.
            else
                area_insisto   =   area_bagnata3
                fattore_flusso =  (abs(am_out3)/massa_tot)
            end if

            !     side 3 constant eta
            if(infout3 == 1)then

                do k=kparasta,kparaend
                    do i=1,jx

                        coef_massa3=(areola3(i,k)/area_insisto)*fattore_flusso

                        somma=somma+coef_massa3*index_out3(i,k)
                        !          Vc^n+1
                        vcc=vc(i,0,k)+del_mas*coef_massa3*index_out3(i,k)

                        vc3_orl(i,k)=vcc
                        !          pressure b.c. = Vc* - Vc^n+1
                        cs3(i,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)-vcc

                        if(potenziale)then
                            cs3(i,k)=-vcc
                        end if
                        !       Vc*
                        vc(i,0,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)

                    enddo
                enddo

            else

                do k=kparasta,kparaend
                    do i=1,jx
                        !          Vc stored to compute the divg
                        vc3_orl(i,k)=vc(i,0,k)
                        !          pressure b.c. = Vc* - Vc^n+1
                        cs3(i,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)-vc(i,0,k)

                        if(potenziale)then
                            cs3(i,k)=-vc(i,0,k)
                        end if
                        !       Vc*
                        vc(i,0,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)

                    enddo
                enddo

            end if !infout3
            !     ..................................................................
            !     side 4 constant eta
            if(massa_tot .lt. 1.d-8)then
                area_insisto   =   area_bagnata
                fattore_flusso =   1.
            else
                area_insisto   =   area_bagnata4
                fattore_flusso =  (abs(am_out4)/massa_tot)
            end if

            if(infout4 == 1)then

                do k=kparasta,kparaend
                    do i=1,jx

                        coef_massa4=(areola4(i,k)/area_insisto)*fattore_flusso

                        somma=somma+coef_massa4*index_out4(i,k)
                        !          Vc^n+1
                        vcc=vc(i,jy,k)-del_mas*coef_massa4*index_out4(i,k)

                        vc4_orl(i,k)=vcc
                        !          pressure b.c. = Vc* - Vc^n+1
                        cs4(i,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)-vcc
                        if(potenziale)then
                            cs4(i,k)=-vcc
                        end if
                        !       Vc*
                        vc(i,jy,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)

                    enddo
                enddo

            else

                do k=kparasta,kparaend
                    do i=1,jx
                        !          Vc stored to compute the divg
                        vc4_orl(i,k)=vc(i,jy,k)
                        !          pressure b.c. = Vc* - Vc^n+1
                        cs4(i,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)-vc(i,jy,k)

                        if(potenziale)then
                            cs4(i,k)=-vc(i,jy,k)
                        end if
                        !       Vc*
                        vc(i,jy,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)

                    enddo
                enddo

            end if !infout4

        enddo    !end loop jj=1,jp
        !.....................................

        do kk=1,kp
            !     side 5 constant zita
            if(myid==0)then

                if(massa_tot .lt. 1.d-8)then
                    area_insisto   =   area_bagnata
                    fattore_flusso =   1.
                else
                    area_insisto   =   area_bagnata5
                    fattore_flusso =  (abs(am_out5)/massa_tot)
                end if

                if(infout5 == 1)then

                    do j=1,jy
                        do i=1,jx
                            !
                            coef_massa5=(areola5(i,j)/area_insisto)*fattore_flusso
                            !
                            somma=somma+coef_massa5*index_out5(i,j)
                            !          Wc^n+1
                            wcc=wc(i,j,0)+del_mas*coef_massa5*index_out5(i,j)

                            wc5_orl(i,j)=wcc
                            !          pressure b.c. = Wc* - Wc^n+1
                            cs5(i,j)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)-wcc
                            if(potenziale)then
                                cs5(i,j)=-wcc
                            end if
                            !    Wc*
                            wc(i,j,0)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)
                        !
                        end do
                    end do

                else

                    do j=1,jy
                        do i=1,jx
                            !          wc stored to compute the divg
                            wc5_orl(i,j)=wc(i,j,0)
                            !          pressure b.c. = Wc* - Wc^n+1
                            cs5(i,j)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)-wc(i,j,0)

                            if(potenziale)then
                                cs5(i,j)=-wc(i,j,0)
                            end if
                            !    Wc*
                            wc(i,j,0)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)
                        !
                        enddo
                    enddo

                end if !infout5

            end if !myid

            !     side 6 constant zita
            if(myid.eq.nproc-1)then

                if(massa_tot .lt. 1.d-8)then
                    area_insisto   =   area_bagnata
                    fattore_flusso =   1.
                else
                    area_insisto   =   area_bagnata6
                    fattore_flusso =  (abs(am_out6)/massa_tot)
                end if

                if(infout6 == 1)then

                    do j=1,jy
                        do i=1,jx

                            coef_massa6=(areola6(i,j)/area_insisto)*fattore_flusso
                            !
                            somma=somma+coef_massa6*index_out6(i,j)
                            !          Wc^n+1
                            wcc=wc(i,j,jz)-del_mas*coef_massa6*index_out6(i,j)

                            wc6_orl(i,j)=wcc
                            !          pressure b.c. = Wc* - Wc^n+1
                            cs6(i,j)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)-wc(i,j,jz)
                            if(potenziale)then
                                cs6(i,j)=-wcc
                            end if
                            !    Wc*
                            wc(i,j,jz)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)
                        end do
                    end do

                else

                    do j=1,jy
                        do i=1,jx
                            !          Wc stored to compute the divg
                            wc6_orl(i,j)=wc(i,j,jz)
                            !          pressure b.c. = Wc* - Wc^n+1
                            cs6(i,j)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)-wc(i,j,jz)

                            if(potenziale)then
                                cs6(i,j)=-wc(i,j,jz)
                            end if
                            !    Wc*
                            wc(i,j,jz)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)
                        !
                        enddo
                    enddo

                end if !infout6

            end if !myid
        !
        enddo    !end loop kk=1,kp
        !
        !
        !
        call MPI_ALLREDUCE(somma,sommat,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !
        !
        if(myid.eq.0)then
            write(*,*)myid,'sum coef',sommat
        end if

        !
        !-----------------------------------------------------------------------
        ! INSIDE THE FIELD
        !-----------------------------------------------------------------------
        !
        !     compute Uc*
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    !
                    !
                    uinter=.5*(u(i,j,k)+u(i+1,j,k))
                    vinter=.5*(v(i,j,k)+v(i+1,j,k))
                    winter=.5*(w(i,j,k)+w(i+1,j,k))
                    !
                    uc(i,j,k)=uinter*csx(i,j,k)+vinter*csy(i,j,k)+winter*csz(i,j,k)
                !
                end do
            end do
        end do

        !     compute Vc*
        do k=kparasta,kparaend
            do j=jp,jy-jp
                do i=1,jx

                    !
                    uinter=.5*(u(i,j,k)+u(i,j+1,k))
                    vinter=.5*(v(i,j,k)+v(i,j+1,k))
                    winter=.5*(w(i,j,k)+w(i,j+1,k))
                    !
                    vc(i,j,k)=uinter*etx(i,j,k)+vinter*ety(i,j,k)+winter*etz(i,j,k)
                !
                end do
            end do
        end do

        !     compute Wc*
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
                    !
                    uinter=.5*(u(i,j,k)+u(i,j,k+1))
                    vinter=.5*(v(i,j,k)+v(i,j,k+1))
                    winter=.5*(w(i,j,k)+w(i,j,k+1))
                    !
                    wc(i,j,k)=uinter*ztx(i,j,k)+vinter*zty(i,j,k)+winter*ztz(i,j,k)
                !
                enddo
            enddo

        enddo
        !
        !-----------------------------------------------------------------------
        ! SEND PLANE wc(k-1)
        !-----------------------------------------------------------------------
        ! send kparaend plane to closer proc in order to compute rhs in diver
        !
        if(myid.eq.0)then
            leftpem=MPI_PROC_NULL
            rightpem=rightpe
        else if(myid.eq.nproc-1)then
            leftpem=leftpe
            rightpem=MPI_PROC_NULL
        else if((myid.ne.0).and.(myid.ne.nproc-1))then
            leftpem=leftpe
            rightpem=rightpe
        endif


        call MPI_SENDRECV(wc(1,1,kparaend),jx*jy,MPI_REAL_SD,rightpem,51+myid,&
            wc(1,1,kparasta-1),jx*jy,MPI_REAL_SD,leftpem,50+myid,MPI_COMM_WORLD,status,ierr)
        !
        !-----------------------------------------------------------------------
        ! MASS CHECK
        !-----------------------------------------------------------------------
        ! check the mass balance after mass distribution

        am_out_loc=0.
        am_out=0.

        am_out_loc1=0.
        am_out_loc2=0.
        am_out_loc3=0.
        am_out_loc4=0.
        am_out_loc5=0.
        am_out_loc6=0.

        am_out1=0.
        am_out2=0.
        am_out3=0.
        am_out4=0.
        am_out5=0.
        am_out6=0.

        !     sides 1 and 2 constant csi
        do ii=1,ip

            if(infout1==1)then
                do k=kparasta,kparaend
                    do j=1,jy

                        am_out_loc1=am_out_loc1-uc1_orl(j,k)

                    end do
                end do
            endif

            if(infout2==1)then
                do k=kparasta,kparaend
                    do j=1,jy

                        am_out_loc2=am_out_loc2+uc2_orl(j,k)

                    end do
                end do
            endif

        end do
        !
        !     sides 3 and 4 constant eta
        do jj=1,jp

            if(infout3==1)then
                do k=kparasta,kparaend
                    do i=1,jx

                        am_out_loc3=am_out_loc3-vc3_orl(i,k)

                    end do
                end do
            endif

            if(infout4==1)then
                do k=kparasta,kparaend
                    do i=1,jx

                        am_out_loc4=am_out_loc4+vc4_orl(i,k)

                    end do
                end do
            endif

        end do
        !
        !     sides 5 and 6 constant zita
        do kk=1,kp

            if(myid==0)then

                if(infout5==1)then
                    do j=1,jy
                        do i=1,jx

                            am_out_loc5=am_out_loc5-wc5_orl(i,j)

                        end do
                    end do
                endif

            elseif(myid==nproc-1)then

                if(infout6==1)then
                    do j=1,jy
                        do i=1,jx

                            am_out_loc6=am_out_loc6+wc6_orl(i,j)

                        end do
                    end do
                endif

            endif

        end do

        am_out_loc=am_out_loc1+am_out_loc2+am_out_loc3+am_out_loc4+am_out_loc5+am_out_loc6

        !
        ! all procs know outflow
        call MPI_ALLREDUCE(am_out_loc,am_out,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !

        ! outflow from each face
        call MPI_ALLREDUCE(am_out_loc1,am_out1,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc2,am_out2,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc3,am_out3,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc4,am_out4,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc5,am_out5,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(am_out_loc6,am_out6,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(myid.eq.0)then
            !        write(*,*)'verifica massa entrante',am_in
            !        write(*,*)'verifica distribuzione massa uscente',am_out
            write(*,*)'check delta ',am_out-am_in
        endif


        return
    end subroutine contra_infout

    subroutine redistribuzione()

        ! mass distribution on boundary to obtain a flow field with
        ! zero divergence. Used for nesting
        !
        !use myarrays_nesting
        use myarrays_velo3
        use mysettings, only: bodyforce


        implicit none
        !
        !-----------------------------------------------------------------------
        !     variables declaration
        integer i,j,k
        integer ierr,status(mpi_status_size)
        real areatot
        real massa1,massa2,massa3,massa4,massa5,massa6
        real massatot
        real massa1tot,massa2tot
        real massa3tot,massa4tot
        real massa5tot,massa6tot
        real bilancio

        real ar1,ar2,ar3,ar4,ar5,ar6
        real artot
        real ar1tot,ar2tot
        real ar3tot,ar4tot
        real ar5tot,ar6tot
        !
        !-----------------------------------------------------------------------
        ! put to zero fluxes at faces 3 and 4
        !
        vc(:, 0,kparasta:kparaend)=0.
        vc(:,jy,kparasta:kparaend)=0.

        !-----------------------------------------------------------------------
        ! mass check
        !     faces 1 and 2
        ar1 = 0.
        massa1=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa1=massa1+uc(0,j,k)*index_out1(j,k)
                if(uc(0,j,k).lt.0.)ar1 = ar1 + areola1(j,k)*index_out1(j,k)
            end do
        end do

        ar2 = 0.
        massa2=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa2=massa2-uc(jx,j,k)*index_out2(j,k)
                if(uc(jx,j,k).gt.0.)ar2 = ar2 + areola2(j,k)*index_out2(j,k)
            end do
        end do

        !     solid walls faces 3 and 4
        massa3=0.
        massa4=0.

        !     faces 5 and 6

        ar5 = 0.
        massa5=0.
        if(myid.eq.0)then
            do j=1,jy
                do i=1,jx
                    massa5=massa5+wc(i,j,0)*index_out5(i,j)
                    if(wc(i,j,0).lt.0.)ar5 = ar5 + areola5(i,j)*index_out5(i,j)
                end do
            end do
        end if

        ar6 = 0.
        massa6=0.
        if(myid.eq.nproc-1)then
            do j=1,jy
                do i=1,jx
                    massa6=massa6-wc(i,j,jz)*index_out6(i,j)
                    if(wc(i,j,jz).gt.0.)ar6 = ar6 + areola6(i,j)*index_out6(i,j)
                end do
            end do
        end if

        !     global mass at the faces
        call MPI_ALLREDUCE(massa1,massa1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa2,massa2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa3,massa3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa4,massa4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa5,massa5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa6,massa6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !     global mass at the faces
        call MPI_ALLREDUCE(ar1,ar1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(ar2,ar2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(ar3,ar3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(ar4,ar4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(ar5,ar5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(ar6,ar6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(myid.eq.0)then
            bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot

            artot = ar1tot + ar2tot + ar3tot + ar4tot + ar5tot + ar6tot
            write(*,*)'mass balance before redistribution: ',bilancio
            write(*,*)massa1tot,massa2tot,massa3tot,massa4tot,massa5tot,massa6tot
        end if

        call MPI_BCAST(bilancio,1,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(artot,1,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)
        !-----------------------------------------------------------------------
        ! excluding ibm and inflow
        areatot = 0.

        if(.not.bodyforce)then
            if(infout1 /= 0 .and. massa1tot<0.)areatot = areatot + area1  !goes out
            if(infout2 /= 0 .and. massa2tot<0.)areatot = areatot + area2
            if(infout5 /= 0 .and. massa5tot<0.)areatot = areatot + area5
            if(infout6 /= 0 .and. massa6tot<0.)areatot = areatot + area6
        else
            if(infout1/=0 .and. massa1tot<0.)areatot=areatot+area_bagnata1
            if(infout2/=0 .and. massa2tot<0.)areatot=areatot+area_bagnata2
            if(infout5/=0 .and. massa5tot<0.)areatot=areatot+area_bagnata5
            if(infout6/=0 .and. massa6tot<0.)areatot=areatot+area_bagnata6
        end if



        !-----------------------------------------------------------------------
        ! fluxes adjustment with mass defeact

        !     face 1
        if(infout1 /= 0)then
            do k=kparasta,kparaend
                do j=1,jy
                    if(uc(0,j,k).lt.0.)then
                        uc(0,j,k)=uc(0,j,k)-bilancio*areola1(j,k)/artot
                        uc(0,j,k)=uc(0,j,k)*index_out1(j,k)
                        ucp1(j,k)=uc(0,j,k)
                    end if
                end do
            end do
        end if
        !     face 2
        if(infout2 /= 0)then
            do k=kparasta,kparaend
                do j=1,jy
                    if(uc(jx,j,k).gt.0.)then
                        uc(jx,j,k)=uc(jx,j,k)+bilancio*areola2(j,k)/artot
                        uc(jx,j,k)=uc(jx,j,k)*index_out2(j,k)
                        ucp2(j,k)=uc(jx,j,k)
                    end if
                end do
            end do
        end if
        !     face 5
        if(infout5 /= 0)then
            if(myid.eq.0)then
                do j=1,jy
                    do i=1,jx
                        if(wc(i,j,0).lt.0.)then
                            wc(i,j,0)=wc(i,j,0)-bilancio*areola5(i,j)/artot
                            wc(i,j,0)=wc(i,j,0)*index_out5(i,j)
                            wcp5(i,j)=wc(i,j,0)
                        end if
                    end do
                end do
            end if
        end if
        !     face 6
        if(infout6 /= 0)then
            if(myid.eq.nproc-1)then
                do j=1,jy
                    do i=1,jx
                        if(wc(i,j,jz).gt.0.)then
                            wc(i,j,jz)=wc(i,j,jz)+bilancio*areola6(i,j)/artot
                            wc(i,j,jz)=wc(i,j,jz)*index_out6(i,j)
                            wcp6(i,j)=wc(i,j,jz)
                        end if
                    end do
                end do
            end if
        end if
        !
        !-----------------------------------------------------------------------
        ! mass check after redistribution
        !     faces 1 and 2
        massa1=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa1=massa1+uc(0,j,k)
            end do
        end do

        massa2=0.
        do k=kparasta,kparaend
            do j=1,jy
                massa2=massa2-uc(jx,j,k)
            end do
        end do

        !     solid walls faces 3 and 4
        massa3=0.
        massa4=0.

        !     faces 5 and 6
        massa5=0.
        if(myid.eq.0)then
            do j=1,jy
                do i=1,jx
                    massa5=massa5+wc(i,j,0)
                end do
            end do
        end if

        massa6=0.
        if(myid.eq.nproc-1)then
            do j=1,jy
                do i=1,jx
                    massa6=massa6-wc(i,j,jz)
                end do
            end do
        end if

        !     global mass at the faces
        call MPI_ALLREDUCE(massa1,massa1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa2,massa2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa3,massa3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa4,massa4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa5,massa5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(massa6,massa6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(myid.eq.0)then
            bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot
            write(*,*)'mass balance after redistribution: ',bilancio
        end if

        ! chicco     NOTA andrebbero forse cambiate le velocita' con i
        ! chicco     coseni direttori in seguito alla ridistribuzione massa

        return
    end subroutine redistribuzione

    subroutine facce()
        !     find area of cells at the sides

        use myarrays_metri3
        use output_module, only: info_run_file

        implicit none

        !-----------------------------------------------------------------------
        !     variables declaration
        integer ierr
        integer i,j,k    !,jx,jy,jz,n1,n2,n3
        integer status(mpi_status_size)
        real area1_loc, area2_loc, area3_loc
        real area4_loc, area5_loc, area6_loc
        real dxdeta,dydeta,dzdeta,dxdzet,dydzet,dzdzet
        real dxdcsi,dydcsi,dzdcsi
        real ra, sump
        !
        !-----------------------------------------------------------------------
        ! area face 1
        i = 0
        area1     = 0.
        area1_loc = 0.

        do k=kparasta,kparaend
            do j=1,n2
                areola1(j,k) = 0.
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2

                dxdeta = .5*( x(i,j,k)  + x(i,j,k-1)   )- &
                    .5*( x(i,j-1,k)+ x(i,j-1,k-1) )

                dydeta = .5*( y(i,j,k)  + y(i,j,k-1)   )- &
                    .5*( y(i,j-1,k)+ y(i,j-1,k-1) )

                dzdeta = .5*( z(i,j,k)  + z(i,j,k-1)   )- &
                    .5*( z(i,j-1,k)+ z(i,j-1,k-1) )

                dxdzet = .5*( x(i,j,k) + x(i,j-1,k)    )- &
                    .5*( x(i,j,k-1)+x(i,j-1,k-1)  )

                dydzet = .5*( y(i,j,k) + y(i,j-1,k)    )- &
                    .5*( y(i,j,k-1)+y(i,j-1,k-1)  )

                dzdzet = .5*( z(i,j,k) + z(i,j-1,k)    )- &
                    .5*( z(i,j,k-1)+z(i,j-1,k-1)  )

                gg221(j,k) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
                gg331(j,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
                gg231(j,k) = dxdeta*dxdzet + dydeta*dydzet + dzdeta*dzdzet

                areola1(j,k)=sqrt(gg221(j,k)*gg331(j,k)-gg231(j,k)*gg231(j,k) )

            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                area1_loc = area1_loc + areola1(j,k)
            end do
        end do
        !
        ! areas sum, the sum is made by myid=0
        !
        if (myid.eq.0) then
            area1 = area1_loc
            sump  = 0.
            ra    = 0.
        end if

        do i=1,nproc-1

            if (myid.eq.0) then

                call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
                sump = sump + ra

            else if(myid.eq.i) then

                call mpi_send(area1_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)

            end if

        end do

        if (myid.eq.0) then
            area1 = area1 + sump
            write(*,*)'area face 1 =',area1
        end if
        !
        ! myid 0 send the computed total area to all procs
        ! mpi broadcast

        call mpi_bcast(area1,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

        write(info_run_file,*)'PE=',myid,'area1=',area1

        !
        !-----------------------------------------------------------------------
        ! area face 2

        i = n1
        area2     = 0.
        area2_loc = 0.


        do k=kparasta,kparaend
            do j=1,n2
                areola2(j,k) = 0.
            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2

                dxdeta = .5*( x(i,j,k)  + x(i,j,k-1)   )- &
                    .5*( x(i,j-1,k)+ x(i,j-1,k-1) )

                dydeta = .5*( y(i,j,k)  + y(i,j,k-1)   )- &
                    .5*( y(i,j-1,k)+ y(i,j-1,k-1) )

                dzdeta = .5*( z(i,j,k)  + z(i,j,k-1)   )- &
                    .5*( z(i,j-1,k)+ z(i,j-1,k-1) )

                dxdzet = .5*( x(i,j,k) + x(i,j-1,k)    )- &
                    .5*( x(i,j,k-1)+x(i,j-1,k-1)  )

                dydzet = .5*( y(i,j,k) + y(i,j-1,k)    )- &
                    .5*( y(i,j,k-1)+y(i,j-1,k-1)  )

                dzdzet = .5*( z(i,j,k) + z(i,j-1,k)    )- &
                    .5*( z(i,j,k-1)+z(i,j-1,k-1)  )

                gg221(j,k) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
                gg331(j,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
                gg231(j,k) = dxdeta*dxdzet + dydeta*dydzet + dzdeta*dzdzet

                areola2(j,k)=sqrt(gg221(j,k)*gg331(j,k)-gg231(j,k)*gg231(j,k))

            end do
        end do

        do k=kparasta,kparaend
            do j=1,n2
                area2_loc = area2_loc + areola2(j,k)
            end do
        end do


        ! sum to have total area

        if (myid.eq.0) then
            area2 = area2_loc
            sump  = 0.
            ra    = 0.
        end if

        do i=1,nproc-1

            if (myid.eq.0) then
                call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
                sump = sump + ra

            else if(myid.eq.i) then
                call mpi_send(area2_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)

            end if

        end do

        if (myid.eq.0) then
            area2 = area2 +sump
            write(*,*)'area face 2 =',area2
        end if

        call mpi_bcast(area2,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

        write(info_run_file,*)'PE=',myid,'area2=',area2

        !
        !-----------------------------------------------------------------------
        ! area face 5

        if (myid.eq.0) then

            k = 0
            area5     = 0.
            area5_loc = 0.

            do j=1,n2
                do i=1,n1
                    areola5(i,j) = 0.
                end do
            end do

            do j=1,n2
                do i=1,n1

                    dxdcsi = .5*( x(i,j,k)   + x(i,j-1,k)   )- &
                        .5*( x(i-1,j,k) + x(i-1,j-1,k) )
                    dydcsi = .5*( y(i,j,k)   + y(i,j-1,k)   )- &
                        .5*( y(i-1,j,k) + y(i-1,j-1,k) )
                    dzdcsi = .5*( z(i,j,k)   + z(i,j-1,k)   )- &
                        .5*( z(i-1,j,k) + z(i-1,j-1,k) )
                    dxdeta = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
                        .5*( x(i,j-1,k) + x(i-1,j-1,k) )
                    dydeta = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
                        .5*( y(i,j-1,k) + y(i-1,j-1,k) )
                    dzdeta = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
                        .5*( z(i,j-1,k) + z(i-1,j-1,k) )

                    gg112(i,j) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
                    gg222(i,j) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
                    gg122(i,j) = dxdcsi*dxdeta + dydcsi*dydeta + dzdcsi*dzdeta

                    areola5(i,j)=sqrt(gg112(i,j)*gg222(i,j)-gg122(i,j)*gg122(i,j))

                end do
            end do

            do j=1,n2
                do i=1,n1
                    area5_loc = area5_loc + areola5(i,j)
                end do
            end do

            area5 = area5_loc

            write(info_run_file,*)'area face 5 =',area5

        end if ! close if PE=0

        call mpi_bcast(area5,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

        write(info_run_file,*)'PE=',myid,'area5=',area5

        !-----------------------------------------------------------------------
        ! area face 6

        if (myid.eq.(nproc-1)) then

            k = n3
            area6     = 0.
            area6_loc = 0.

            do j=1,n2
                do i=1,n1
                    areola6(i,j) = 0.
                end do
            end do

            do j=1,n2
                do i=1,n1

                    dxdcsi = .5*( x(i,j,k)   + x(i,j-1,k)   )- &
                        .5*( x(i-1,j,k) + x(i-1,j-1,k) )
                    dydcsi = .5*( y(i,j,k)   + y(i,j-1,k)   )- &
                        .5*( y(i-1,j,k) + y(i-1,j-1,k) )
                    dzdcsi = .5*( z(i,j,k)   + z(i,j-1,k)   )- &
                        .5*( z(i-1,j,k) + z(i-1,j-1,k) )
                    dxdeta = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
                        .5*( x(i,j-1,k) + x(i-1,j-1,k) )
                    dydeta = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
                        .5*( y(i,j-1,k) + y(i-1,j-1,k) )
                    dzdeta = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
                        .5*( z(i,j-1,k) + z(i-1,j-1,k) )

                    gg112(i,j) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
                    gg222(i,j) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
                    gg122(i,j) = dxdcsi*dxdeta + dydcsi*dydeta + dzdcsi*dzdeta

                    areola6(i,j)=sqrt(gg112(i,j)*gg222(i,j)-gg122(i,j)*gg122(i,j))

                end do
            end do

            do j=1,n2
                do i=1,n1
                    area6_loc = area6_loc + areola6(i,j)
                end do
            end do

            area6 = area6_loc

            write(info_run_file,*)'area della faccia 6 =',area6

        end if ! close if PE=nproc-1

        call mpi_bcast(area6,1,MPI_REAL_SD,nproc-1,mpi_comm_world,ierr)

        write(info_run_file,*)'PE=',myid,'area6=',area6

        !
        !-----------------------------------------------------------------------
        ! area face 3

        j = 0
        area3     = 0.
        area3_loc = 0.

        do k=kparasta,kparaend
            do i=1,n1
                areola3(i,k) = 0.
            end do
        end do

        do k=kparasta,kparaend
            do i=1,n1

                dxdcsi = .5*( x(i,j,k)   + x(i,j,k-1)   )- &
                    .5*( x(i-1,j,k) + x(i-1,j,k-1) )
                dydcsi = .5*( y(i,j,k)   + y(i,j,k-1)   )- &
                    .5*( y(i-1,j,k) + y(i-1,j,k-1) )
                dzdcsi = .5*( z(i,j,k)   + z(i,j,k-1)   )- &
                    .5*( z(i-1,j,k) + z(i-1,j,k-1) )
                dxdzet = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
                    .5*( x(i,j,k-1) + x(i-1,j,k-1) )
                dydzet = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
                    .5*( y(i,j,k-1) + y(i-1,j,k-1) )
                dzdzet = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
                    .5*( z(i,j,k-1) + z(i-1,j,k-1) )

                gg113(i,k) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
                gg333(i,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
                gg133(i,k) = dxdcsi*dxdzet + dydcsi*dydzet + dzdcsi*dzdzet

                areola3(i,k)=sqrt(gg113(i,k)*gg333(i,k)-gg133(i,k)*gg133(i,k))

            end do
        end do

        do k=kparasta,kparaend
            do i=1,n1
                area3_loc = area3_loc + areola3(i,k)
            end do
        end do

        ! call to sum total area

        if (myid.eq.0) then
            area3 = area3_loc
            sump  = 0.
            ra    = 0.
        end if

        do i=1,nproc-1

            if (myid.eq.0) then
                call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
                sump = sump + ra

            else if(myid.eq.i) then
                call mpi_send(area3_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)

            end if

        end do

        if (myid.eq.0) then
            area3 = area3 + sump
            write(*,*)'area face 3 =',area3
        end if

        call mpi_bcast(area3,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

        write(info_run_file,*)'PE=',myid,'area3=',area3

        !
        !-----------------------------------------------------------------------
        ! area face 4

        j = n2
        area4     = 0.
        area4_loc = 0.

        do k=kparasta,kparaend
            do i=1,n1
                areola4(i,k) = 0.
            end do
        end do

        do k=kparasta,kparaend
            do i=1,n1

                dxdcsi = .5*( x(i,j,k)   + x(i,j,k-1)   )- &
                    .5*( x(i-1,j,k) + x(i-1,j,k-1) )
                dydcsi = .5*( y(i,j,k)   + y(i,j,k-1)   )- &
                    .5*( y(i-1,j,k) + y(i-1,j,k-1) )
                dzdcsi = .5*( z(i,j,k)   + z(i,j,k-1)   )- &
                    .5*( z(i-1,j,k) + z(i-1,j,k-1) )
                dxdzet = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
                    .5*( x(i,j,k-1) + x(i-1,j,k-1) )
                dydzet = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
                    .5*( y(i,j,k-1) + y(i-1,j,k-1) )
                dzdzet = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
                    .5*( z(i,j,k-1) + z(i-1,j,k-1) )

                gg113(i,k) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
                gg333(i,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
                gg133(i,k) = dxdcsi*dxdzet + dydcsi*dydzet + dzdcsi*dzdzet

                areola4(i,k)=sqrt(gg113(i,k)*gg333(i,k)-gg133(i,k)*gg133(i,k))

            end do
        end do

        do k=kparasta,kparaend
            do i=1,n1
                area4_loc = area4_loc + areola4(i,k)
            end do
        end do

        ! call to know the total area

        if (myid.eq.0) then
            area4 = area4_loc
            sump  = 0.
            ra    = 0.
        end if

        do i=1,nproc-1

            if (myid.eq.0) then
                call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
                sump = sump + ra

            else if (myid.eq.i) then
                call mpi_send(area4_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)
            end if

        end do

        if (myid.eq.0) then
            area4 = area4 + sump
            write(*,*)'area face 4 =',area4
        end if

        call mpi_bcast(area4,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

        write(info_run_file,*)'PE=',myid,'area4=',area4


    end subroutine facce

    subroutine aree_parziali()
        !***********************************************************************
        use myarrays_metri3
        use ibm_module
        use mysettings, only: lett,i_rest,bodyforce

        implicit none

        !-----------------------------------------------------------------------
        ! array declaration
        integer i,j,k,isc,ierr,l
        integer cont1,cont2,cont3,cont4,cont5,cont6

        integer inf1,inf2,inf3,inf4,inf5,inf6

        real parziale1,parziale2
        real parziale3,parziale4
        real parziale5,parziale6
        !
        !-----------------------------------------------------------------------
        !    OUTFLOW
        !-----------------------------------------------------------------------
        ! in case of partial outflow
        ! if no outflow index = 0
        ! if    outflow index = 1

        !-----------------------------------------------------------------------
        !     INDEXES
        inf1 = 0
        inf2 = 0
        inf3 = 0
        inf4 = 0
        inf5 = 0
        inf6 = 0

        !case inflow outflow
        if(lett .and. infout1==1)inf1=1
        if(lett .and. infout2==1)inf2=1
        if(lett .and. infout3==1)inf3=1
        if(lett .and. infout4==1)inf4=1
        if(lett .and. infout5==1)inf5=1
        if(lett .and. infout6==1)inf6=1

        !case nesting
        if(i_rest==3 .and. infout1/=0)inf1=1
        if(i_rest==3 .and. infout2/=0)inf2=1
        if(i_rest==3 .and. infout5/=0)inf5=1
        if(i_rest==3 .and. infout6/=0)inf6=1



        !-----------------------------------------------------------------------
        !     initialization

        area_bagnata1 = area1
        area_bagnata2 = area2
        area_bagnata3 = area3
        area_bagnata4 = area4
        area_bagnata5 = area5
        area_bagnata6 = area6

        parziale1 = 0.
        parziale2 = 0.
        parziale3 = 0.
        parziale4 = 0.
        parziale5 = 0.
        parziale6 = 0.

        do k=1,jz
            do j=1,jy
                index_out1(j,k)=1
                index_out2(j,k)=1

                index_rho1(j,k)=1
                index_rho2(j,k)=1
            end do
        end do

        do k=1,jz
            do i=1,jx
                index_out3(i,k)=1
                index_out4(i,k)=1

                index_rho3(i,k)=1
                index_rho4(i,k)=1
            end do
        end do

        do j=1,jy
            do i=1,jx
                index_out5(i,j)=1
                index_out6(i,j)=1

                index_rho5(i,j)=1
                index_rho6(i,j)=1
            end do
        end do

        cont1=0
        cont2=0
        cont3=0
        cont4=0
        cont5=0
        cont6=0

        !-----------------------------------------------------------------------
        !     IBM settings for outflow condition
        if( (lett .or. i_rest==3) .and. bodyforce)then
            ! stop orlansky on solid cells
            if (.not.particles) then
                open (22,file='Celle_Bloccate_Indici.inp',status='old')
                read (22,*) num_solide
            end if
            do l=1,num_solide
                if (particles) then
                    i=indici_celle_bloccate(l,1)
                    j=indici_celle_bloccate(l,2)
                    k=indici_celle_bloccate(l,3)
                else
                    read(22,*)i,j,k
                end if

                !         face 1
                if(i.eq.1)then
                    index_out1(j,k)=0
                    cont1=cont1+1
                end if
                !         face 2
                if(i.eq.jx)then
                    index_out2(j,k)=0
                    cont2=cont2+1
                end if
                !         face 3
                if(j.eq.1)then
                    index_out3(i,k)=0
                    cont3=cont3+1
                end if
                !         face 4
                if(j.eq.jy)then
                    index_out4(i,k)=0
                    cont4=cont4+1
                end if
                !         face 5
                if(k.eq.1)then
                    index_out5(i,j)=0
                    cont5=cont5+1
                end if
                !         face 6
                if(k.eq.jz)then
                    index_out6(i,j)=0
                    cont6=cont6+1
                end if
            end do
            if (.not.particles) then
                close(22)
            end if

            !      stop orlansky on ib cells
            if (.not.particles) then
                open (22,file='Celle_IB_indici.inp',status='old')
                read (22,*) num_ib
            end if
            do l=1,num_ib
                if (particles) then
                    i = indici_CELLE_IB(l,1)
                    j = indici_CELLE_IB(l,2)
                    k = indici_CELLE_IB(l,3)
                else
                    read(22,*)i,j,k
                end if
                !         face 1
                if(i.eq.1)then
                    index_out1(j,k)=1 !0
                    cont1=cont1+1
                end if
                !         face 2
                if(i.eq.jx)then
                    index_out2(j,k)=1 !0
                    cont2=cont2+1
                end if
                !         face 3
                if(j.eq.1)then
                    index_out3(i,k)=1 !0
                    cont3=cont3+1
                end if
                !         face 4
                if(j.eq.jy)then
                    index_out4(i,k)=1 !0
                    cont4=cont4+1
                end if
                !         face 5
                if(k.eq.1)then
                    index_out5(i,j)=1 !0
                    cont5=cont5+1
                end if
                !         face 6
                if(k.eq.jz)then
                    index_out6(i,j)=1 !0
                    cont6=cont6+1
                end if
            end do
            if (.not.particles) then
                close(22)
            end if
        end if

        !-----------------------------------------------------------------------
        ! compute partial face surface (it will be reduced if IBM is on)

        !     faces 1 and 2
        if(inf1==1)then
            parziale1=0.
            area_bagnata1=0.
            do k=kparasta,kparaend  !1,jz
                do j=1,jy
                    parziale1 = parziale1+areola1(j,k)*real(index_out1(j,k))
                end do
            end do

            call MPI_ALLREDUCE(parziale1,area_bagnata1,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end if
        !
        if(inf2==1)then
            parziale2=0.
            area_bagnata2=0.
            do k=kparasta,kparaend  !1,jz
                do j=1,jy
                    parziale2 = parziale2+areola2(j,k)*real(index_out2(j,k))
                end do
            end do

            call MPI_ALLREDUCE(parziale2,area_bagnata2,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end if
        !
        !     faces 3 and 4
        if(inf3==1)then
            parziale3=0.
            area_bagnata3=0.
            do k=kparasta,kparaend  !1,jz
                do i=1,jx
                    parziale3 = parziale3+areola3(i,k)*real(index_out3(i,k))
                end do
            end do

            call MPI_ALLREDUCE(parziale3,area_bagnata3,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end if
        !
        if(inf4==1)then
            parziale4=0.
            area_bagnata4=0.
            do k=kparasta,kparaend  !1,jz
                do i=1,jx
                    parziale4 = parziale4+areola4(i,k)*real(index_out4(i,k))
                end do
            end do

            call MPI_ALLREDUCE(parziale4,area_bagnata4,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end if
        !
        !     faces 5 and 6
        if(inf5==1)then
            if(myid.eq.0)then
                parziale5=0.
                area_bagnata5=0.
                do j=1,jy
                    do i=1,jx
                        parziale5 = parziale5+areola5(i,j)*real(index_out5(i,j))
                    end do
                end do
                area_bagnata5=parziale5
            end if
            call MPI_BCAST(area_bagnata5,1,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)
        end if

        if(inf6==1)then
            if(myid.eq.nproc-1)then
                parziale6=0.
                area_bagnata6=0.
                do j=1,jy
                    do i=1,jx
                        parziale6 = parziale6+areola6(i,j)*real(index_out6(i,j))
                    end do
                end do
                area_bagnata6=parziale6
            end if
            call MPI_BCAST(area_bagnata6,1,MPI_REAL_SD,nproc-1,MPI_COMM_WORLD,ierr)
        end if


        return
    end subroutine aree_parziali

end module inflow_module
