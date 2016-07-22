subroutine vel_up()

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
    integer :: ierr,status(MPI_STATUS_SIZE)
    integer :: kparastal,kparaendl
    integer plantype
    integer i,j,k
    !-----------------------------------------------------------------------

    ! velocity  at step n+1 into the field
    !
    if (.not.freesurface) then !no free surface <<<<<<<<<<<<<<<
        if (.not.potenziale) then
      
            do k=kparasta,kparaend
                do j=1,n2
                    do i=1,n1
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
                do j=1,n2
                    do i=1,n1
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
            do j=1,n2
                do i=1,n1
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
    if (ip==0) then

        do k=kparasta,kparaend
            do j=0,n2+1
                !
                u(0,j,k)=u(n1,j,k)
                v(0,j,k)=v(n1,j,k)
                w(0,j,k)=w(n1,j,k)
                u(n1+1,j,k)=u(1,j,k)
                v(n1+1,j,k)=v(1,j,k)
                w(n1+1,j,k)=w(1,j,k)
            !
            end do
        end do

    end if

    ! periodic cell in zita (also outside the domain)
    if (kp==0) then

        if (1==0) then

        call MPI_TYPE_VECTOR(n2+2,n1,n1+2,MPI_REAL_SD,plantype,ierr)
        call MPI_TYPE_COMMIT(plantype,ierr)



        if (myid.eq.nproc-1) then
            call MPI_SENDRECV(u(1,0,n3),1,plantype,0,11,u(1,0,n3+1),1,plantype,0,12,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(v(1,0,n3),1,plantype,0,13,v(1,0,n3+1),1,plantype,0,14,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(w(1,0,n3),1,plantype,0,15,w(1,0,n3+1),1,plantype,0,16,MPI_COMM_WORLD,status,ierr)

        endif

        if (myid.eq.0) then
            call MPI_SENDRECV(u(1,0,1),1,plantype,nproc-1,12,u(1,0,0),1,plantype,nproc-1,11,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(v(1,0,1),1,plantype,nproc-1,14,v(1,0,0),1,plantype,nproc-1,13,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(w(1,0,1),1,plantype,nproc-1,16,w(1,0,0),1,plantype,nproc-1,15,MPI_COMM_WORLD,status,ierr)

        endif

        call MPI_TYPE_FREE(plantype,ierr)

        else



        end if



    end if
    !
    ! compute controvariant velocity: cycle 0,jx for Uc
    !                                 cycle 0,jz for Wc
    !
    if (.not.freesurface) then !no free surface <<<<<<<<<<<<<<<
        if (.not.potenziale) then
            do k=kparasta,kparaend
                do j=1,n2
                    do i=ip,n1-ip
                        uc(i,j,k)=uc(i,j,k)-dt*cgra1(i,j,k)
                    end do
                end do
            end do
            !
            do k=kparasta,kparaend
                do j=1,n2-1
                    do i=1,n1
                        vc(i,j,k)=vc(i,j,k)-dt*cgra2(i,j,k)
                    end do
                end do
            end do
      
        else if (potenziale) then
      
            do k=kparasta,kparaend
                do j=1,n2
                    do i=ip,n1-ip
                        uc(i,j,k)=-dt*cgra1(i,j,k)
                    end do
                end do
            end do
            !
            do k=kparasta,kparaend
                do j=1,n2-1
                    do i=1,n1
                        vc(i,j,k)=-dt*cgra2(i,j,k)
                    end do
                end do
            end do
        end if
    else if (freesurface) then !free surface on<<<<<<<<<<<<<<<<
        do k=kparasta,kparaend
            do j=1,n2
                do i=ip,n1-ip
                    uc(i,j,k)=uc(i,j,k)-dt*cgra1(i,j,k)
                end do
            end do
        end do
        !
        do k=kparasta,kparaend
            do j=0,n2  !including also the surface gradient
                do i=1,n1
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
                do j=1,n2
                    do i=1,n1
                        wc(i,j,k)=wc(i,j,k)-dt*cgra3(i,j,k)
                    end do
                end do
            end do

        else if (potenziale) then

            do k=kparastal,kparaendl
                do j=1,n2
                    do i=1,n1
                        wc(i,j,k)=-dt*cgra3(i,j,k)
                    end do
                end do
            end do
        end if

    else if (freesurface) then !free surface on<<<<<<<<<<<<<<<<
        do k=kparastal,kparaendl
            do j=1,n2
                do i=1,n1
                    wc(i,j,k)=wc(i,j,k)-dt*cgra3(i,j,k)
                end do
            end do
        end do
    end if !free surface <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! I need Wc(k-1) to compute divmax so sending procedure for kparaend
    call MPI_SENDRECV(wc(1,1,kparaend),n1*n2,MPI_REAL_SD,rightpem,71+myid, &
        wc(1,1,kparasta-1),n1*n2,MPI_REAL_SD,leftpem,70+myid,MPI_COMM_WORLD,status,ierr)

    return
end
