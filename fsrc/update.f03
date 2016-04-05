!***********************************************************************
subroutine update()
    !***********************************************************************
    ! update for intermediate velocity into the flow field
    ! at the walls parabolic extrapolation
    !
    use myarrays_velo3
    use mysending
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer ierr
    integer status(MPI_STATUS_SIZE)
    integer plantypee
    integer i,j,k,ii,jj,kk
    !-----------------------------------------------------------------------
    real,parameter :: pcoef1=1.875
    real,parameter :: pcoef2=-1.250
    real,parameter :: pcoef3=0.375
    !-----------------------------------------------------------------------

    !     into the field
    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
                u(i,j,k)=u(i,j,k)+delu(i,j,k)
                v(i,j,k)=v(i,j,k)+delv(i,j,k)
                w(i,j,k)=w(i,j,k)+delw(i,j,k)
            enddo
        enddo
    enddo
    !
    !
    ! extrapolation left and right sides (1 and 2)
    !
    do ii=1,ip
        !
        do k=kparasta,kparaend
            do j=1,jy
                !     extrapolation with parable y=ax2+bx+c
                !
                u(0,j,k)=1.875*u(1,j,k)-1.25*u(2,j,k)+.375*u(3,j,k)
                v(0,j,k)=1.875*v(1,j,k)-1.25*v(2,j,k)+.375*v(3,j,k)
                w(0,j,k)=1.875*w(1,j,k)-1.25*w(2,j,k)+.375*w(3,j,k)
                !
                u(jx+1,j,k)=.375*u(jx-2,j,k)-1.25*u(jx-1,j,k)+1.875*u(jx,j,k)
                v(jx+1,j,k)=.375*v(jx-2,j,k)-1.25*v(jx-1,j,k)+1.875*v(jx,j,k)
                w(jx+1,j,k)=.375*w(jx-2,j,k)-1.25*w(jx-1,j,k)+1.875*w(jx,j,k)
            !
            end do
        end do
    !
    end do
    !
    !     periodicity on csi
    do ii=1,1-ip
      
        do k=kparasta,kparaend
            do j=1,jy
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
    !-----------------------------------------------------------------------
    !     extrapolation on bottom and upper sides (3 and 4)
    do jj=1,jp
        !
        do k=kparasta,kparaend
            do i=1,jx
                !
                u(i,0,k)=1.875*u(i,1,k)-1.25*u(i,2,k)+.375*u(i,3,k)
                v(i,0,k)=1.875*v(i,1,k)-1.25*v(i,2,k)+.375*v(i,3,k)
                w(i,0,k)=1.875*w(i,1,k)-1.25*w(i,2,k)+.375*w(i,3,k)
                !
                u(i,jy+1,k)=.375*u(i,jy-2,k)-1.25*u(i,jy-1,k)+1.875*u(i,jy,k)
                v(i,jy+1,k)=.375*v(i,jy-2,k)-1.25*v(i,jy-1,k)+1.875*v(i,jy,k)
                w(i,jy+1,k)=.375*w(i,jy-2,k)-1.25*w(i,jy-1,k)+1.875*w(i,jy,k)
            !
            end do
        end do

    end do
    !
    !     periodicity on eta
    do jj=1,1-jp
        !
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
    !
    end do

    !-----------------------------------------------------------------------
    ! extrapolation on front and back sides (5 and 6)
    do kk=1,kp
        !
        if (myid.eq.0) then
            do j=1,jy
                do i=1,jx
                    u(i,j,0)=1.875*u(i,j,1)-1.25*u(i,j,2)+.375*u(i,j,3)
                    v(i,j,0)=1.875*v(i,j,1)-1.25*v(i,j,2)+.375*v(i,j,3)
                    w(i,j,0)=1.875*w(i,j,1)-1.25*w(i,j,2)+.375*w(i,j,3)
                enddo
            enddo
        endif
        !
        if (myid.eq.nproc-1) then
            do j=1,jy
                do i=1,jx
                    u(i,j,jz+1)=.375*u(i,j,jz-2)-1.25*u(i,j,jz-1)+1.875*u(i,j,jz)
                    v(i,j,jz+1)=.375*v(i,j,jz-2)-1.25*v(i,j,jz-1)+1.875*v(i,j,jz)
                    w(i,j,jz+1)=.375*w(i,j,jz-2)-1.25*w(i,j,jz-1)+1.875*w(i,j,jz)
                enddo
            enddo
        endif
    !
    enddo

    !
    !     periodicity on zita

    do kk=1,1-kp
        !
        call MPI_TYPE_VECTOR(jy,jx,jx+2,MPI_REAL_SD,plantypee,ierr)
        call MPI_TYPE_COMMIT(plantypee,ierr)

        if (myid.eq.nproc-1) then

            call MPI_SENDRECV(u(1,1,jz),1,plantypee,0,11,u(1,1,jz+1),1,plantypee,0,12,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(v(1,1,jz),1,plantypee,0,13,v(1,1,jz+1),1,plantypee,0,14,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(w(1,1,jz),1,plantypee,0,15,w(1,1,jz+1),1,plantypee,0,16,MPI_COMM_WORLD,status,ierr)

        else if (myid.eq.0) then

            call MPI_SENDRECV(u(1,1,1),1,plantypee,nproc-1,12,u(1,1,0),1,plantypee,nproc-1,11,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(v(1,1,1),1,plantypee,nproc-1,14,v(1,1,0),1,plantypee,nproc-1,13,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(w(1,1,1),1,plantypee,nproc-1,16,w(1,1,0),1,plantypee,nproc-1,15,MPI_COMM_WORLD,status,ierr)

        endif


        call MPI_TYPE_FREE(plantypee,ierr)

    enddo
    !
    ! every proc needs to know values in kparaend+1 to compute
    ! controvariant velocity in contra

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


    if(leftpem /= MPI_PROC_NULL) then
        call MPI_SSEND(u(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
        call MPI_SSEND(v(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
        call MPI_SSEND(w(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
    endif
    if(rightpem /= MPI_PROC_NULL) then
        call MPI_RECV(u(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(v(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(w(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
    endif

!    ! ALESSANDRO
!    if(leftpem /= MPI_PROC_NULL) then
!        call MPI_SSEND(u(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
!        call MPI_SSEND(v(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
!        call MPI_SSEND(w(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
!    endif
!    if(rightpem /= MPI_PROC_NULL) then
!        call MPI_RECV(u(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
!        call MPI_RECV(v(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
!        call MPI_RECV(w(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
!    endif


    return
end
