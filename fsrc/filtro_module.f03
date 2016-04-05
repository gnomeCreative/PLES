module filtro_module

    use mysending
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    integer,bind(C) :: ifiltro,nfiltro

    integer,bind(C) :: filtrou,filtrov,filtrow,filtrorho,filtrofi

    integer,bind(C) :: xstart,xend,ystart,yend,zstart,zend

contains

    subroutine initialize_filtro()

        use output_module, only: info_run_file

        if (ifiltro==1) then

            if (myid==0) then
                write(*,*)'PAY ATTENTION YOU ARE FILTERING THE FLOW FIELD'
                write(info_run_file,*)'FILTERING ON',ifiltro
                if (xend>n1 .or. yend>n2 .or.zend>n3 &
                    .or. xstart<1 .or. ystart<1 .or. zstart<1) then
                    write(*,*)'FILTERING AREA TOO LARGE'
                    write(*,*)'OUT OF BOUNDS'
                    write(*,*)'xend:',xend,',n1:',n1
                    write(*,*)'yend:',yend,',n2:',n2
                    write(*,*)'zend:',zend,',n3:',n3
                    write(*,*)'xend:',xstart,1
                    write(*,*)'yend:',ystart,1
                    write(*,*)'zend:',zstart,1

                    write(info_run_file,*)'FILTERING AREA TOO LARGE'
                    write(info_run_file,*)'OUT OF BOUNDS'
                    write(info_run_file,*)'xend:',xend,',n1:',n1
                    write(info_run_file,*)'yend:',yend,',n2:',n2
                    write(info_run_file,*)'zend:',zend,',n3:',n3
                    write(info_run_file,*)'xend:',xstart,1
                    write(info_run_file,*)'yend:',ystart,1
                    write(info_run_file,*)'zend:',zstart,1

                    stop

                end if
            end if

            if (zstart<=kparasta) then
                zstart=kparasta
            end if
            if (zend>=kparaend) then
                zend=kparaend
            endif

            if (zstart>kparaend .or. zend<kparasta) then
                zstart=kparaend
                zend=  kparaend -1
                write(*,*)'proc: ',myid,'no filtering',zstart,zend
                write(info_run_file,*)'proc: ',myid,'no filtering',zstart,zend
            else
                write(*,*) 'proc:',myid,'filter in k between',zstart,'and',zend
                write(info_run_file,*)'proc:',myid,'filter in k between',zstart,'and',zend
            end if

        else
            if (myid == 0) then
                write(*,*)'FILTERING OFF',ifiltro
                write(info_run_file,*)'FILTERING OFF',ifiltro
            end if
        end if

    end subroutine initialize_filtro

    subroutine filtering_procedure(ktime)

        use myarrays_velo3, only: u,v,w,fi,rhov

        implicit none

        integer,intent(in) :: ktime
        !
        real,allocatable :: rho(:,:,:)

        integer i,j,k,isc

        ! --------------------------------------------

        if (ifiltro == 1 .and. ktime==nfiltro*(ktime/nfiltro)) then

            if (filtrou==1) then
                call filtro(u)
                if (myid==0) then
                    write(*,*)'call filtering for u'
                end if
            end if

            if (filtrov==1) then
                call filtro(v)
                if (myid==0) then
                    write(*,*)'call filtering for v'
                end if
            end if


            if (filtrow==1) then
                call filtro(w)
                if (myid==0) then
                    write(*,*)'call filtering for w'
                end if
            end if


            if (filtrorho==1) then
                allocate(rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)) !0:n3+1))
                do isc=1,nscal
                    do k=kparasta-1,kparaend+1
                        do j=0,n2+1
                            do i=0,n1+1
                                rho(i,j,k)=rhov(isc,i,j,k)
                            end do
                        end do
                    end do

                    call filtro(rho)
                    if (myid==0) then
                        write(*,*)'call filtering for rho'
                    end if

                    do k=kparasta-1,kparaend+1
                        do j=0,n2+1
                            do i=0,n1+1
                                rhov(isc,i,j,k)=rho(i,j,k)
                            end do
                        end do
                    end do
                end do
                deallocate(rho)
            end if

            if (filtrofi==1) then
                call filtro(fi)
                if (myid==0) then
                    write(*,*)'call filtering for fi'
                end if
            end if

        end if

    end subroutine filtering_procedure

    !***********************************************************************
    subroutine filtro(r)
        !***********************************************************************
        ! applica il test filter su csi ed eta
        !
        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k
          
        integer istatus,ierr
        integer req1,req2,req3,req4

        real  r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        real rf(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
     
        real r_destra,r_sinistra
        real r_sotto,r_sopra
        real r_indietro,r_avanti

        integer xstart,xend,ystart,yend,zstart,zend

        integer status(MPI_STATUS_SIZE)
        integer plantype
        !-----------------------------------------------------------------------
        !     initialize
        !
        do k=kparasta-1,kparaend+1
            do j=0,jy+1
                do i=0,jx+1
                    rf(i,j,k)=0.
                enddo
            enddo
        enddo

        !     direction k

        ! if filter is applied first to k it is not necessary an intermediate sending
        ! for planes
        do k=zstart,zend !kparasta,kparaend
            do j=ystart,yend
                do i=xstart,xend
      
                    r_indietro=0.5*(r(i,j,k)+r(i,j,k-1))
                    r_avanti=  0.5*(r(i,j,k)+r(i,j,k+1))

                    if (k==1) then
                        r_indietro=r(i,j,0)
                    end if

                    if (k==n3) then
                        r_avanti=r(i,j,n3+1)
                    end if
                
                    rf(i,j,k)=(r_indietro+2.*r(i,j,k)+r_avanti)*0.25
    
                enddo
            enddo
        enddo

        do k=zstart,zend !kparasta,kparaend
            do j=ystart,yend
                do i=xstart,xend
                    r(i,j,k)=rf(i,j,k)
                enddo
            enddo
        enddo

        !     direction i
        do k=zstart,zend !kparasta,kparaend
            do j=ystart,yend
                do i=xstart,xend
      
                    r_sinistra=0.5*(r(i,j,k)+r(i-1,j,k))
                    r_destra=0.5*(r(i,j,k)+r(i+1,j,k))

                    if (i==1) then
                        r_sinistra=r(0,j,k)
                    end if

                    if (i==n1) then
                        r_destra=r(n1+1,j,k)
                    end if

                    rf(i,j,k)=(r_sinistra+2.*r(i,j,k)+r_destra)*0.25
          
                enddo
            enddo
        enddo

        do k=zstart,zend !kparasta,kparaend
            do j=ystart,yend
                do i=xstart,xend
                    r(i,j,k)=rf(i,j,k)
                enddo
            enddo
        enddo
      
        !     direction j
        do k=zstart,zend !kparasta,kparaend
            do j=ystart,yend
                do i=xstart,xend
      
      
                    r_sotto= 0.5*(r(i,j,k)+r(i,j-1,k))
                    r_sopra= 0.5*(r(i,j,k)+r(i,j+1,k))

                    if (j==1) then
                        r_sotto=r(i,0,k)
                    end if

                    if (j==n2) then
                        r_sopra=r(i,n2+1,k)
                    end if
      
                    rf(i,j,k)=(r_sotto+2.*r(i,j,k)+r_sopra)*0.25

                enddo
            enddo
        enddo

        do k=zstart,zend !kparasta,kparaend
            do j=ystart,yend
                do i=xstart,xend
                    r(i,j,k)=rf(i,j,k)
                enddo
            enddo
        enddo
         
        !-----------------------------------------------------------------------
        !    send planes at the border
        !-----------------------------------------------------------------------

        if (myid.ne.0) then
            call MPI_SSEND(r(0,0,kparasta),(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,tagls, &
                MPI_COMM_WORLD,ierr)
        !         call MPI_WAIT(req1,istatus,ierr)
        endif
      
        if (myid.ne.nproc-1) then
            call MPI_RECV(r(0,0,kparaend+1),(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrr, &
                MPI_COMM_WORLD,status,ierr)
        !         call MPI_WAIT(req2,istatus,ierr)
        endif

        if (myid.ne.nproc-1) then
            call MPI_SSEND(r(0,0,kparaend),(jx+2)*(jy+2), &
                MPI_REAL_SD,rightpem,tagrs, &
                MPI_COMM_WORLD,ierr)
        !         call MPI_WAIT(req3,istatus,ierr)
        endif
        if (myid.ne.0) then
            call MPI_RECV(r(0,0,kparasta-1),(jx+2)*(jy+2), &
                MPI_REAL_SD,leftpem,taglr, &
                MPI_COMM_WORLD,status,ierr)
        !         call MPI_WAIT(req4,istatus,ierr)
        endif
        !
        !     if periodic
        do k=1,1-kp

            call MPI_TYPE_VECTOR(jy+2,jx,jx+2,MPI_REAL_SD,plantype,ierr)
            call MPI_TYPE_COMMIT(plantype,ierr)


            if (myid==nproc-1) then
                call MPI_SENDRECV(r(1,0,jz),1, &
                    plantype,0,11, &
                    r(1,0,jz+1),1, &
                    plantype,0,12, &
                    MPI_COMM_WORLD,status,ierr)


            endif

            if (myid==0) then
                call MPI_SENDRECV(r(1,0,1),1, &
                    plantype,nproc-1,12, &
                    r(1,0,0),1, &
                    plantype,nproc-1,11, &
                    MPI_COMM_WORLD,status,ierr)

            endif

            call MPI_TYPE_FREE(plantype,ierr)

        enddo
        !hicco serve riapplicare periodicity per i etc...
        return
    end subroutine filtro

end module filtro_module
