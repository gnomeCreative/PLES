!***********************************************************************
subroutine passo_ibm()

    ! used to be:
    ! passo_ibm(correggo_delu)

    !***********************************************************************
    ! send IB at the border between the procs
    use myarrays_ibm
    use mysending
    use myarrays_velo3
    !
    use scala3
    use period
    use tipologia
    !
    use mpi


    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer i,j,k,l,isc
    integer iper,kper,jper
    integer ireq1,ireq2,ireq3,ireq4
    integer status(MPI_STATUS_SIZE),ierror
    !integer correggo_delu
      
    logical icheck
      
    !-----------------------------------------------------------------------
    icheck = .false.
    !-----------------------------------------------------------------------
    if(icheck)then
        write(900+myid,*)myid,num_left_snd,num_right_snd,num_left_rcv,num_right_rcv
    end if

    !     left send, right recive
    if(num_left_snd .ne.0)then
        allocate(sbuff_ibm( (3+nscal)*num_left_snd )) 
        sbuff_ibm = 0.
    end if
    if(num_right_rcv .ne.0)then
        allocate(rbuff_ibm( (3+nscal)*num_right_rcv))
        rbuff_ibm = 0.
    end if

    if(num_left_snd .ne.0)then
        if(myid.ne.0)then
            do l=1,num_left_snd
                i = stencil_left_snd(l,1)
                j = stencil_left_snd(l,2)
                k = stencil_left_snd(l,3)
	 	
                sbuff_ibm(l               ) = u(i,j,k)
                sbuff_ibm(l+  num_left_snd) = v(i,j,k)
                sbuff_ibm(l+2*num_left_snd) = w(i,j,k)
                do isc=1,nscal
                    sbuff_ibm(l+3*num_left_snd &
                        +(isc-1)*num_left_snd) &
                        = rhov(isc,i,j,k)
                end do
            end do

            call MPI_SSEND(sbuff_ibm(1),(3+nscal)*num_left_snd,MPI_REAL_SD, &
                leftpe,tagls,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (ireq1,status,ierror)
        end if
    end if
      
    if(num_right_rcv.ne.0)then
        if(myid.ne.nproc-1)then
            call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_right_rcv,MPI_REAL_SD, &
                rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
            !      call MPI_WAIT (ireq2,status,ierror)
  
            do l=1,num_right_rcv
                i = stencil_right_rcv(l,1)
                j = stencil_right_rcv(l,2)
                k = stencil_right_rcv(l,3)
	 	
                u(i,j,k) = rbuff_ibm(l                )
                v(i,j,k) = rbuff_ibm(l+  num_right_rcv)
                w(i,j,k) = rbuff_ibm(l+2*num_right_rcv)
                do isc=1,nscal
                    rhov(isc,i,j,k)= rbuff_ibm(l+3*num_right_rcv  &
                        +(isc-1)*num_right_rcv)
                end do
            end do
        end if
    end if

    ! now periodicity
    do kper = 1,1-kp
        if(num_left_snd .ne.0)then
            if(myid.eq.0)then
                do l=1,num_left_snd
                    i = stencil_left_snd(l,1)
                    j = stencil_left_snd(l,2)
                    k = stencil_left_snd(l,3)
	 	
                    sbuff_ibm(l               ) = u(i,j,k)
                    sbuff_ibm(l+  num_left_snd) = v(i,j,k)
                    sbuff_ibm(l+2*num_left_snd) = w(i,j,k)
                    do isc=1,nscal
                        sbuff_ibm(l+3*num_left_snd+(isc-1)*num_left_snd) &
                            = rhov(isc,i,j,k)
                    end do
                end do

                call MPI_SSEND(sbuff_ibm(1),(3+nscal)*num_left_snd, &
                    MPI_REAL_SD,nproc-1,tagls,MPI_COMM_WORLD,ierror)
            !          call MPI_WAIT (ireq1,status,ierror)
            end if
        end if

        if(num_right_rcv .ne.0)then	
            if(myid.eq.nproc-1)then
                call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_right_rcv, &
                    MPI_REAL_SD,0,tagrr,MPI_COMM_WORLD,status,ierror)
                !          call MPI_WAIT (ireq2,status,ierror)


                ! chicco AAA nno puo' andare bene qui con la modifica applicata!
                ! deve inserire questi dati non in u(i,j,k) ma in u_piano(i,j)
                ! oppure direttamente su u(i,j,0)
 
                do l=1,num_right_rcv
                    i = stencil_right_rcv(l,1)
                    j = stencil_right_rcv(l,2)
                    k = stencil_right_rcv(l,3)
	 	
                    u(i,j,k) = rbuff_ibm(l                )
                    v(i,j,k) = rbuff_ibm(l+  num_right_rcv)
                    w(i,j,k) = rbuff_ibm(l+2*num_right_rcv)
                    do isc=1,nscal
                        rhov(isc,i,j,k)  &
                            = rbuff_ibm(l+3*num_right_rcv+(isc-1)*num_right_rcv)
                    end do
                end do
            end if
        end if
    end do


    if(num_left_snd .ne.0)then
        deallocate(sbuff_ibm) 
    end if
    if(num_right_rcv .ne.0)then
        deallocate(rbuff_ibm)
    end if


    !-----------------------------------------------------------------------
    !     right send, left recive
    !
    if(num_right_snd .ne.0)then
        allocate(sbuff_ibm( (3+nscal)*num_right_snd ))
    end if
    if(num_left_rcv.ne.0)then
        allocate(rbuff_ibm( (3+nscal)*num_left_rcv))
    end if

    if(num_right_snd .ne.0)then
        if(myid.ne.nproc-1)then
            do l=1,num_right_snd
                i = stencil_right_snd(l,1)
                j = stencil_right_snd(l,2)
                k = stencil_right_snd(l,3)
	 	
                sbuff_ibm(l                ) = u(i,j,k)
                sbuff_ibm(l+  num_right_snd) = v(i,j,k)
                sbuff_ibm(l+2*num_right_snd) = w(i,j,k)
                do isc=1,nscal
                    sbuff_ibm(l+3*num_right_snd+(isc-1)*num_right_snd) &
                        =rhov(isc,i,j,k)
                end do
            end do

            call MPI_SSEND(sbuff_ibm(1),(3+nscal)*num_right_snd,MPI_REAL_SD, &
                rightpe,tagrs,MPI_COMM_WORLD,ierror)
        !      call MPI_WAIT (ireq3,status,ierror)
        end if
    end if
      
    if(num_left_rcv.ne.0)then
        if(myid.ne.0)then
            call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_left_rcv,MPI_REAL_SD, &
                leftpe,taglr,MPI_COMM_WORLD,status,ierror)
            !      call MPI_WAIT (ireq4,status,ierror)
  
            do l=1,num_left_rcv
                i = stencil_left_rcv(l,1)
                j = stencil_left_rcv(l,2)
                k = stencil_left_rcv(l,3)
	 	
                u(i,j,k) = rbuff_ibm(l               )
                v(i,j,k) = rbuff_ibm(l+  num_left_rcv)
                w(i,j,k) = rbuff_ibm(l+2*num_left_rcv)
                do isc=1,nscal
                    rhov(isc,i,j,k)  &
                        = rbuff_ibm(l+3*num_left_rcv+(isc-1)*num_left_rcv)
                end do
            end do
        end if
    end if
      
    !     now periodicity
    do kper = 1,1-kp
        if(num_right_snd .ne.0)then
            if(myid.eq.nproc-1)then
                do l=1,num_right_snd
                    i = stencil_right_snd(l,1)
                    j = stencil_right_snd(l,2)
                    k = stencil_right_snd(l,3)
	 	
                    sbuff_ibm(l                ) = u(i,j,k)
                    sbuff_ibm(l+  num_right_snd) = v(i,j,k)
                    sbuff_ibm(l+2*num_right_snd) = w(i,j,k)
                    do isc=1,nscal
                        sbuff_ibm(l+3*num_right_snd+(isc-1)*num_right_snd) &
                            =rhov(isc,i,j,k)
                    end do
                end do

                call MPI_SSEND(sbuff_ibm(1),(3+nscal)*num_right_snd, &
                    MPI_REAL_SD,0,tagrs,MPI_COMM_WORLD,ierror)
            !           call MPI_WAIT (ireq3,status,ierror)
            end if
        end if
	 
        if(num_left_rcv.ne.0)then
            if(myid.eq.0)then
                call MPI_RECV(rbuff_ibm(1),(3+nscal)*num_left_rcv, &
                    MPI_REAL_SD,nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
                !           call MPI_WAIT (ireq4,status,ierror)
  
                do l=1,num_left_rcv
                    i = stencil_left_rcv(l,1)
                    j = stencil_left_rcv(l,2)
                    k = stencil_left_rcv(l,3)
	 	
                    u(i,j,k) = rbuff_ibm(l               )
                    v(i,j,k) = rbuff_ibm(l+  num_left_rcv)
                    w(i,j,k) = rbuff_ibm(l+2*num_left_rcv)
                    do isc=1,nscal
                        rhov(isc,i,j,k)  &
                            = rbuff_ibm(l+3*num_left_rcv+(isc-1)*num_left_rcv)
                    end do
                end do
            end if
        end if
    end do

    if(num_right_snd .ne.0)then
        deallocate(sbuff_ibm)
    end if
    if(num_left_rcv.ne.0)then
        deallocate(rbuff_ibm)
    end if

    !-----------------------------------------------------------------------
    !     i direction, impose periodicity or wall
    !
    !     parabolic extrapolation y=ax2+bx+c
    if(ip==1)then
      
!        if(correggo_delu==1)then
!            do k=kparasta,kparaend
!                do j=1,jy
!                    !
!                    u(0,j,k)=1.875*u(1,j,k) &
!                        -1.25*u(2,j,k) &
!                        +.375*u(3,j,k)
!                    v(0,j,k)=1.875*v(1,j,k) &
!                        -1.25*v(2,j,k) &
!                        +.375*v(3,j,k)
!                    w(0,j,k)=1.875*w(1,j,k) &
!                        -1.25*w(2,j,k) &
!                        +.375*w(3,j,k)
!                    !
!                    u(jx+1,j,k)=.375*u(jx-2,j,k) &
!                        -1.25*u(jx-1,j,k) &
!                        +1.875*u(jx  ,j,k)
!                    !
!                    v(jx+1,j,k)=.375*v(jx-2,j,k) &
!                        -1.25*v(jx-1,j,k) &
!                        +1.875*v(jx  ,j,k)
!                    !
!                    w(jx+1,j,k)=.375*w(jx-2,j,k) &
!                        -1.25*w(jx-1,j,k) &
!                        +1.875*w(jx  ,j,k)
!                !
!                end do
!            end do
!        end if
    !
    else
        !
        !     periodicity
        !
        do k=kparasta,kparaend
            do j=1,jy
      
                u(0   ,j,k)=u(jx,j,k)
                v(0   ,j,k)=v(jx,j,k)
                w(0   ,j,k)=w(jx,j,k)
                u(jx+1,j,k)=u(1 ,j,k)
                v(jx+1,j,k)=v(1 ,j,k)
                w(jx+1,j,k)=w(1 ,j,k)
            !
            end do
        end do
      
    end if
    !
    !-----------------------------------------------------------------------------
    !     j direction, impose periodicity or wall
    !
    !     parabolic extrapolation y=ax2+bx+c
    if(jp==1)then
        !
!        if(correggo_delu==1)then
!            do k=kparasta,kparaend
!                do i=1,jx
!
!                    u(i,0,k)=1.875*u(i,1,k) &
!                        -1.25*u(i,2,k) &
!                        +.375*u(i,3,k)
!                    v(i,0,k)=1.875*v(i,1,k) &
!                        -1.25*v(i,2,k) &
!                        +.375*v(i,3,k)
!                    w(i,0,k)=1.875*w(i,1,k) &
!                        -1.25*w(i,2,k) &
!                        +.375*w(i,3,k)
!                    !
!                    u(i,jy+1,k)=.375*u(i,jy-2,k) &
!                        -1.25*u(i,jy-1,k) &
!                        +1.875*u(i,jy  ,k)
!                    !
!                    v(i,jy+1,k)=.375*v(i,jy-2,k) &
!                        -1.25*v(i,jy-1,k) &
!                        +1.875*v(i,jy  ,k)
!                    !
!                    w(i,jy+1,k)=.375*w(i,jy-2,k) &
!                        -1.25*w(i,jy-1,k) &
!                        +1.875*w(i,jy  ,k)
!                !
!                end do
!            end do
!        end if
      
    else
        !
        !     periodicity on eta
        do k=kparasta,kparaend
            do i=1,jx
      
                u(i,   0,k)=u(i,jy,k)
                v(i,   0,k)=v(i,jy,k)
                w(i,   0,k)=w(i,jy,k)
                u(i,jy+1,k)=u(i, 1,k)
                v(i,jy+1,k)=v(i, 1,k)
                w(i,jy+1,k)=w(i, 1,k)
            !
            end do
        end do
      
    end if
    !
    !-----------------------------------------------------------------------------
    !    k direction, impose periodicity or wall
    !
    ! extrapolation on front and back sides (5 and 6)
    if(kp==1)then
      
!        if(correggo_delu==1)then
!            if (myid.eq.0) then
!                do j=1,jy
!                    do i=1,jx
!                        !
!                        u(i,j,0)=1.875*u(i,j,1) &
!                            -1.25*u(i,j,2) &
!                            +.375*u(i,j,3)
!                        v(i,j,0)=1.875*v(i,j,1) &
!                            -1.25*v(i,j,2) &
!                            +.375*v(i,j,3)
!                        w(i,j,0)=1.875*w(i,j,1) &
!                            -1.25*w(i,j,2) &
!                            +.375*w(i,j,3)
!                    end do
!                end do
!            endif
!            !
!            if (myid.eq.nproc-1) then
!                do j=1,jy
!                    do i=1,jx
!
!                        u(i,j,jz+1)=.375*u(i,j,jz-2) &
!                            -1.25*u(i,j,jz-1) &
!                            +1.875*u(i,j,jz  )
!                        !
!                        v(i,j,jz+1)=.375*v(i,j,jz-2) &
!                            -1.25*v(i,j,jz-1) &
!                            +1.875*v(i,j,jz  )
!                        !
!                        w(i,j,jz+1)=.375*w(i,j,jz-2) &
!                            -1.25*w(i,j,jz-1) &
!                            +1.875*w(i,j,jz  )
!                    !
!                    enddo
!                enddo
!            end if
!        end if

    end if
      
    ! the periodicity in k has been already applied with the message passing
    !-----------------------------------------------------------------------------
    if(icheck)then
      
        k=kparasta
        do j=1,jy
            do i=1,jx
                write(4500+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
            end do
        end do

        k=kparasta-1
        do j=1,jy
            do i=1,jx
                write(4600+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
            end do
        end do

        k=kparaend
        do j=1,jy
            do i=1,jx
                write(4700+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
            end do
        end do

        k=kparaend+1
        do j=1,jy
            do i=1,jx
                write(4800+myid,*)u(i,j,k),v(i,j,k),w(i,j,k)
            end do
        end do

 
    end if
    !-----------------------------------------------------------------------------
    return
end
