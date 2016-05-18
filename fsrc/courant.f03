!***********************************************************************
subroutine courant(ktime)
    !***********************************************************************
    ! compute the time step of the simulation or the max courant number
    ! depending on the simulation settings
    !
    use myarrays_metri3, only: annit, g11, g22, g33, giac
    use myarrays_velo3, only: uc, vc, wc
    use mysending, only: kparasta,kparaend,nproc,myid,MPI_REAL_SD
    use mysettings, only: ind_cou,cou,espl,i_rest
    !
    use scala3
    !
    use mpi

    implicit none

    integer,intent(in) :: ktime

    !-----------------------------------------------------------------------
    !     array declaration
    integer ierr,i,j,k,lll
    integer imax,jmax,kmax
      
    real auc,avc,awc
    real cou_int,dtint,epp,den
    real cou_loc,cou_int_loc
    real dt_loc,dt_loc2,dtint_loc
    real dtc,dtd,coef
    real cmax
    !-----------------------------------------------------------------------
    !     potential flow time step equal to 1.
    if(potenziale)then
        dt = 1.
    else
        !
        !-----------------------------------------------------------------------
        !     FIXED TIME STEP
        !     compute cou max with fixed dt
        if (ind_cou==0) then

            dt=dt_start

            cou_loc=0.

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        !
                        auc=max( (abs(uc(i,j,k))) , (abs(uc(i-1,j,k))) )
                        avc=max( (abs(vc(i,j,k))) , (abs(vc(i,j-1,k))) )
                        awc=max( (abs(wc(i,j,k))) , (abs(wc(i,j,k-1))) )
                        cou_int_loc=(auc+avc+awc)*dt/giac(i,j,k)
                        cou_loc=max(cou_loc,cou_int_loc)
                    !
                    end do
                end do
            end do

            call MPI_ALLREDUCE(cou_loc,cou,1,MPI_REAL_SD,MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (myid==0) then
                write(*,*)'test courant max = ',cou
                write(*,*)'dt',dt
            end if
      


            ! put a control, if cou greater than a value, move from
            ! fixed dt to fixed courant

!            if (cou>0.75) then
!                cou = 0.75
!                if(ktime==1.and.i_rest==0)then
!                    dt_loc=0.01
!                    do k=kparasta,kparaend
!                        do j=1,jy
!                            do i=1,jx
!                                auc=abs(uc(i,j,k))
!                                avc=abs(vc(i,j,k))
!                                awc=abs(wc(i,j,k))
!                                den=auc+avc+awc
!                                epp=0.000001
!                                den=max(epp,den)
!                                dtint_loc=cou*giac(i,j,k)/den
!                                dt_loc=min(dt_loc,dtint_loc)
!                            end do
!                        end do
!                    end do
!                else
!                    auc=abs(uc(1,1,kparasta))
!                    avc=abs(vc(1,1,kparasta))
!                    awc=abs(wc(1,1,kparasta))
!                    den=auc+avc+awc
!                    epp=0.000001
!                    den=max(epp,den)
!                    dt_loc=cou*giac(1,1,kparasta)/den
!                    do k=kparasta,kparaend
!                        do j=1,jy
!                            do i=1,jx
!                                auc=abs(uc(i,j,k))
!                                avc=abs(vc(i,j,k))
!                                awc=abs(wc(i,j,k))
!                                den=auc+avc+awc
!                                epp=0.000001
!                                den=max(epp,den)
!                                dtint_loc=cou*giac(i,j,k)/den
!                                dt_loc=min(dt_loc,dtint_loc)
!                            end do
!                        end do
!                    end do
!
!                end if
!
!                ! I need to find the min dt
!                call MPI_ALLREDUCE(dt_loc,dtc,1,MPI_REAL_SD,MPI_MIN, &
!                    MPI_COMM_WORLD,ierr)
!                if(myid==0)then
!                    write(*,*)'PAY ATTENTION dt MOVE TO NOT CONSTANT '
!                    write(*,*)'dt convective = ',dtc
!                endif
!
!                dt=dtc
!
!                if(espl==1)then
!                    !          0.4<0.5 requested diffusive condition
!                    coef=0.4
!                else
!                    coef=10.
!                end if
!
!                ! now the diffusive time step
!                do k=kparasta,kparaend
!                    do j=1,jy
!                        do i=1,jx
!                            epp=0.000001
!                            den=( (g11(i,j,k)+g22(i,j,k)+g33(i,j,k))* &
!                                annit(i,j,k) )/giac(i,j,k)
!                            den=max(epp,den)
!                            dt_loc2=coef/den
!                            dt=min(dt,dt_loc2)
!                        end do
!                    end do
!                end do
!                !        the min value between all procs
!                dtd=0.
!                call MPI_ALLREDUCE(dt,dtd,1,MPI_REAL_SD,MPI_MIN, &
!                    MPI_COMM_WORLD,ierr)
!                dt=dtd
!                if(myid==0)then
!                    write(*,*)'dt diffusive = ',dt
!                    write(*,*)'time step dt = ',dt
!                endif
!
!            end if
        !
        !-----------------------------------------------------------------------
        !     CONSTANT COURANT NUMBER
        !     compute dt with a fixed couraant number
        else if (ind_cou==1) then
      
            if(ktime==1.and.i_rest==0)then
                dt_loc=0.01

                do k=kparasta,kparaend
                    do j=1,jy
                        do i=1,jx
                            auc=abs(uc(i,j,k))
                            avc=abs(vc(i,j,k))
                            awc=abs(wc(i,j,k))
                            den=auc+avc+awc
                            epp=0.000001
                            den=max(epp,den)
                            dtint_loc=cou*giac(i,j,k)/den
                            dt_loc=min(dt_loc,dtint_loc)
                        end do
                    end do
                end do
            else

                auc=abs(uc(1,1,kparasta))
                avc=abs(vc(1,1,kparasta))
                awc=abs(wc(1,1,kparasta))
                den=auc+avc+awc
                epp=0.000001
                den=max(epp,den)
                dt_loc=cou*giac(1,1,kparasta)/den

                do k=kparasta,kparaend
                    do j=1,jy
                        do i=1,jx
                            auc=abs(uc(i,j,k))
                            avc=abs(vc(i,j,k))
                            awc=abs(wc(i,j,k))
                            den=auc+avc+awc
                            epp=0.000001
                            den=max(epp,den)
                            dtint_loc=cou*giac(i,j,k)/den
                            dt_loc=min(dt_loc,dtint_loc)
                        end do
                    end do
                end do
            end if

            !
            ! I need to find the min dt
            call MPI_ALLREDUCE(dt_loc,dtc,1,MPI_REAL_SD,MPI_MIN,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if(myid==0)then
                write(*,*)'dt convective = ',dtc
            endif

            dt=dtc

            if(espl==1)then
                !       0.4<0.5 requested diffusive condition
                coef=0.4
            else
                coef=10.
            end if

            !     now a check on the diffusive time step
            !     (I don't want a too large difference between diff and conv
            !     to avoid errors, although the implicit is stable)
            do k=kparasta,kparaend
                do j=1,jy
                    do i=1,jx
                        !
                        epp=0.000001

                        den=( (g11(i,j,k)+g22(i,j,k)+g33(i,j,k))* &
                            annit(i,j,k) )/giac(i,j,k)

                        den=max(epp,den)

                        dt_loc2=coef/den

                        dt=min(dt,dt_loc2)

                    end do
                end do
            end do

            !    the min between the procs
            dtd=0.
            call MPI_ALLREDUCE(dt,dtd,1,MPI_REAL_SD,MPI_MIN,MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            dt=dtc

            if(myid==0)then
                write(*,*)'dt diffusive = ',dtd
                write(*,*)'time step dt = ',dt
            endif
        !-----------------------------------------------------------------------
        !
        end if
      
    end if ! potential flow

    return
end
