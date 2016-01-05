!***********************************************************************
subroutine restart(ti,dbbx)
    !***********************************************************************
    !
    ! read restart file old_res_form
    !
    use myarrays_velo3
    use mysending
    use scala3

    implicit none
    !-----------------------------------------------------------------------
    ! variables declaration
    integer i,j,k,isc
    real ti
    real bbx,dbbx
    integer nscali
    integer kpsta,kpsta_deep,kpend_deep
    !
    real val_u,val_v,val_w,val_fi
    real value
    real, allocatable :: val_rho(:)
    !-----------------------------------------------------------------------

    kpsta_deep = kparasta - deepl
    kpend_deep = kparaend + deepr
    write(11,*)myid,'DEEP',deepl,deepr

    kpsta = kparasta
    if(myid.eq.0)kpsta=kparasta-1


    open(10,file='old_res_form',status='old')

    read(10,*)nscali
    read(10,*)ti
    read(10,*)bbx,dbbx !bbx will be taken from Agenerale.in

    allocate(val_rho(nscal))

    !     velocity field
    do k=0,jz+1
        do j=0,jy+1
            do i=0,jx+1
                !      read(10,1000)u(i,j,k),v(i,j,k),w(i,j,k),fi(i,j,k),
                !     >             (rhov(isc,i,j,k),isc=1,nscal)
                read(10,1000)val_u,val_v,val_w,val_fi, &
                    (val_rho(isc),isc=1,nscal)
      
                if(k.ge. kpsta_deep .and. k.le.kpend_deep)then
                    u(i,j,k)   = val_u
                    v(i,j,k)   = val_v
                    w(i,j,k)   = val_w
                    fi(i,j,k)  = val_fi
                    do isc=1,nscal
                        rhov(isc,i,j,k) = val_rho(isc)
                    end do
                end if
            end do
        end do
    end do
    !
    !     controvariant flux
    do k=1,jz
        do j=1,jy
            do i=0,jx
                read(10,1000)value
                if(k.ge.kparasta .and. k.le.kparaend)then
                    uc(i,j,k) = value
                end if
            end do
        end do
    end do

    do k=1,jz
        do j=0,jy
            do i=1,jx
                read(10,1000)value
                if(k.ge.kparasta .and. k.le.kparaend)then
                    vc(i,j,k) = value
                end if
            end do
        end do
    end do
      
    do k=0,jz
        do j=1,jy
            do i=1,jx
                read(10,1000)value
                if(k.ge.kparasta-1 .and. k.le.kparaend)then
                    wc(i,j,k) = value
                end if
      
            end do
        end do
    end do
    close(10)
1000 format(10e18.10)
    !
    return
end
