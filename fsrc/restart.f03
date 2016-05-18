!***********************************************************************
subroutine restart(restart_file,ti)
    !***********************************************************************
    !
    ! read restart file old_res_form
    !
    use myarrays_velo3
    use mysending
    use scala3
    use output_module, only: info_run_file
    use mysettings, only: bbx

    implicit none
    !-----------------------------------------------------------------------
    ! variables declaration
    real,intent(out) :: ti
    character(len=500),intent(in) :: restart_file

    integer i,j,k,isc
    integer nscali ! dump variable
    integer kpsta,kpsta_deep,kpend_deep
    integer,parameter :: oldresfile_id=10
    !
    real val_u,val_v,val_w,val_fi
    real value,dumpvalue
    real, allocatable :: val_rho(:)
    !-----------------------------------------------------------------------

    kpsta_deep=kparasta-deepl
    kpend_deep=kparaend+deepr

    write(info_run_file,*)myid,'DEEP',deepl,deepr

    kpsta = kparasta
    if(myid.eq.0) then
        kpsta=kparasta-1
    end if

    open(oldresfile_id,file=trim(restart_file),status='old')

    read(oldresfile_id,*) nscali
    read(oldresfile_id,*) ti
    read(oldresfile_id,*) bbx,dumpvalue !bbx will be taken from Agenerale.in

    allocate(val_rho(nscal))

    !     velocity field
    do k=0,jz+1
        do j=0,jy+1
            do i=0,jx+1
                !      read(10,*)u(i,j,k),v(i,j,k),w(i,j,k),fi(i,j,k),
                !     >             (rhov(isc,i,j,k),isc=1,nscal)
                read(oldresfile_id,*)val_u,val_v,val_w,val_fi, &
                    (val_rho(isc),isc=1,nscal)
      
                if(k.ge.kpsta_deep .and. k.le.kpend_deep)then
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
                read(oldresfile_id,*)value
                if(k.ge.kparasta .and. k.le.kparaend)then
                    uc(i,j,k) = value
                end if
            end do
        end do
    end do

    do k=1,jz
        do j=0,jy
            do i=1,jx
                read(oldresfile_id,*)value
                if(k.ge.kparasta .and. k.le.kparaend)then
                    vc(i,j,k) = value
                end if
            end do
        end do
    end do
      
    do k=0,jz
        do j=1,jy
            do i=1,jx
                read(oldresfile_id,*)value
                if(k.ge.kparasta-1 .and. k.le.kparaend)then
                    wc(i,j,k) = value
                end if
      
            end do
        end do
    end do
    close(oldresfile_id)

    return
end
