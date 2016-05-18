!***********************************************************************
subroutine vel_up()
    !***********************************************************************
    ! compute the cartesian and controvariant velocity at step n+1
    !
    use myarrays_velo3
    use mysending
    use mysettings, only: freesurface
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer ierr,m
    integer status(MPI_STATUS_SIZE)
    integer kparastal,kparaendl
    integer plantype

    integer i,j,k,lll
      
    !     for immersed boundary
    !      integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

    !-----------------------------------------------------------------------
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    !-----------------------------------------------------------------------
    ! velocity  at step n+1 into the field
    !
    if (.not.freesurface) then !no free surface <<<<<<<<<<<<<<<
        if (.not.potenziale) then
      
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        !
                        u(i,j,k)=u(i,j,k)+gra1(i,j,k)
                        v(i,j,k)=v(i,j,k)+gra2(i,j,k)
                        w(i,j,k)=w(i,j,k)+gra3(i,j,k)
                    !
                    end do
                end do
            end do

        else if (potenziale) then

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        !
                        u(i,j,k)=gra1(i,j,k)
                        v(i,j,k)=gra2(i,j,k)
                        w(i,j,k)=gra3(i,j,k)
                    !
                    end do
                end do
            end do

        end if

    else if (freesurface) then !free surface on<<<<<<<<<<<<<<<<
        !
        do k=kparasta,kparaend
            do j=1,jy
                do i=1,jx
                    u(i,j,k)=u(i,j,k)+gra1(i,j,k)
                    v(i,j,k)=v(i,j,k)+gra2(i,j,k)
                    w(i,j,k)=w(i,j,k)+gra3(i,j,k)
                end do
            end do
        end do
    !
    end if !free surface <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    !
    ! periodic cell in csi (also outside the domain)
    !
    do i=1,1-ip

        do k=kparasta,kparaend
            do j=0,jy+1
                !
                u(0   ,j,k)=u(jx,j,k)
                v(0   ,j,k)=v(jx,j,k)
                w(0   ,j,k)=w(jx,j,k)
                u(jx+1,j,k)=u(1 ,j,k)
                v(jx+1,j,k)=v(1 ,j,k)
                w(jx+1,j,k)=w(1 ,j,k)
            !
            end do
        end do

    end do
    !
    !
    ! periodic cell in eta (also outside the domain)
    !
    do j=1,1-jp
        do i=ip,jx+1-ip

            if(myid.eq.0)then
                kparastal=kp
                kparaendl=kparaend
            else if (myid.eq.nproc-1) then
                kparastal=kparasta
                kparaendl=kparaend+1-kp
            else
                kparastal=kparasta
                kparaendl=kparaend
            endif

            do k=kparastal,kparaendl
                u(i,0   ,k)=u(i,jy,k)
                v(i,0   ,k)=v(i,jy,k)
                w(i,0   ,k)=w(i,jy,k)
                u(i,jy+1,k)=u(i,1 ,k)
                v(i,jy+1,k)=v(i,1 ,k)
                w(i,jy+1,k)=w(i,1 ,k)
            end do
        end do
    end do
    !
    ! periodic cell in zita (also outside the domain)
    !
    do k=1,1-kp

        call MPI_TYPE_VECTOR(jy+2,jx,jx+2,MPI_REAL_SD,plantype,ierr)
        call MPI_TYPE_COMMIT(plantype,ierr)



        if (myid.eq.nproc-1) then
            call MPI_SENDRECV(u(1,0,jz),1,plantype,0,11,u(1,0,jz+1),1,plantype,0,12,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(v(1,0,jz),1,plantype,0,13,v(1,0,jz+1),1,plantype,0,14,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(w(1,0,jz),1,plantype,0,15,w(1,0,jz+1),1,plantype,0,16,MPI_COMM_WORLD,status,ierr)

        endif

        if (myid.eq.0) then
            call MPI_SENDRECV(u(1,0,1),1,plantype,nproc-1,12,u(1,0,0),1,plantype,nproc-1,11,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(v(1,0,1),1,plantype,nproc-1,14,v(1,0,0),1,plantype,nproc-1,13,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(w(1,0,1),1,plantype,nproc-1,16,w(1,0,0),1,plantype,nproc-1,15,MPI_COMM_WORLD,status,ierr)

        endif

        call MPI_TYPE_FREE(plantype,ierr)

    end do
    !
    ! compute controvariant velocity: cycle 0,jx for Uc
    !                                 cycle 0,jz for Wc
    !
    if (.not.freesurface) then !no free surface <<<<<<<<<<<<<<<
        if (.not.potenziale) then
            do k=kparasta,kparaend
                do j=1,jy
                    do i=ip,jx-ip
                        uc(i,j,k)=uc(i,j,k)-dt*cgra1(i,j,k)
                    end do
                end do
            end do
            !
            do k=kparasta,kparaend
                do j=jp,jy-jp
                    do i=1,jx
                        vc(i,j,k)=vc(i,j,k)-dt*cgra2(i,j,k)
                    end do
                end do
            end do
      
        else if (potenziale) then
      
            do k=kparasta,kparaend
                do j=1,jy
                    do i=ip,jx-ip
                        uc(i,j,k)=-dt*cgra1(i,j,k)
                    end do
                end do
            end do
            !
            do k=kparasta,kparaend
                do j=jp,jy-jp
                    do i=1,jx
                        vc(i,j,k)=-dt*cgra2(i,j,k)
                    end do
                end do
            end do
        end if
    else if (freesurface) then !free surface on<<<<<<<<<<<<<<<<
        do k=kparasta,kparaend
            do j=1,jy
                do i=ip,jx-ip
                    uc(i,j,k)=uc(i,j,k)-dt*cgra1(i,j,k)
                end do
            end do
        end do
        !
        do k=kparasta,kparaend
            do j=0,jy  !including also the surface gradient
                do i=1,jx
                    vc(i,j,k)=vc(i,j,k)-dt*cgra2(i,j,k)
                end do
            end do
        end do
    end if !free surface <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      
    !
    ! I need to specify the values in k=0
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

    if (.not.freesurface) then !no free surface <<<<<<<<<<<<<<<
        if (.not.potenziale) then
            do k=kparastal,kparaendl
                do j=1,jy
                    do i=1,jx
                        wc(i,j,k)=wc(i,j,k)-dt*cgra3(i,j,k)
                    end do
                end do
            end do

        else if (potenziale) then

            do k=kparastal,kparaendl
                do j=1,jy
                    do i=1,jx
                        wc(i,j,k)=-dt*cgra3(i,j,k)
                    end do
                end do
            end do
        end if

    else if (freesurface) then !free surface on<<<<<<<<<<<<<<<<
        do k=kparastal,kparaendl
            do j=1,jy
                do i=1,jx
                    wc(i,j,k)=wc(i,j,k)-dt*cgra3(i,j,k)
                end do
            end do
        end do
    end if !free surface <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !
    ! I need Wc(k-1) to compute divmax
    ! so sending procedure for kparaend
    !
    if (myid.eq.0) then
        leftpem=MPI_PROC_NULL
        rightpem=rightpe
    else if (myid.eq.nproc-1) then
        leftpem=leftpe
        rightpem=MPI_PROC_NULL
    else if ((myid.ne.0).and.(myid.ne.nproc-1)) then
        leftpem=leftpe
        rightpem=rightpe
    end if

    call MPI_SENDRECV(wc(1,1,kparaend),jx*jy, &
        MPI_REAL_SD,rightpem,71+myid, &
        wc(1,1,kparasta-1),jx*jy, &
        MPI_REAL_SD,leftpem,70+myid, &
        MPI_COMM_WORLD,status,ierr)

    return
end
