subroutine restart(restart_file,ti)
    !
    ! read restart file old_res_form
    !
    use myarrays_velo3, only: fi,rhov,u,uc,v,vc,w,wc
    use scala3, only: n1,n2,n3,nscal,deepl,deepr,kparaend,kparasta

    implicit none
    !-----------------------------------------------------------------------
    ! variables declaration
    real,intent(out) :: ti
    character(len=500),intent(in) :: restart_file

    integer i,j,k,isc
    integer nscali ! dump variable
    integer kpsta_deep,kpend_deep
    integer,parameter :: oldresfile_id=10
    !
    real val_u,val_v,val_w,val_fi
    real value,dumpvalue
    real, allocatable :: val_rho(:)
    !-----------------------------------------------------------------------

    kpsta_deep=kparasta-deepl
    kpend_deep=kparaend+deepr

    open(oldresfile_id,file=trim(restart_file),status='old')

    read(oldresfile_id,*) nscali
    read(oldresfile_id,*) ti
    read(oldresfile_id,*) dumpvalue,dumpvalue !bbx will be taken from Agenerale.in

    allocate(val_rho(nscal))

    !     velocity field
    do k=0,n3+1
        do j=0,n2+1
            do i=0,n1+1

                read(oldresfile_id,*) val_u,val_v,val_w,val_fi,(val_rho(isc),isc=1,nscal)
      
                if (k>=kpsta_deep .and. k<=kpend_deep) then
                    u(i,j,k)=val_u
                    v(i,j,k)=val_v
                    w(i,j,k)=val_w
                    fi(i,j,k)=val_fi
                    do isc=1,nscal
                        rhov(isc,i,j,k)=val_rho(isc)
                    end do
                end if
            end do
        end do
    end do

    !  controvariant flux
    do k=1,n3
        do j=1,n2
            do i=0,n1
                read(oldresfile_id,*) value
                if (k>=kparasta .and. k<=kparaend) then
                    uc(i,j,k)=value
                end if
            end do
        end do
    end do

    do k=1,n3
        do j=0,n2
            do i=1,n1
                read(oldresfile_id,*) value
                if (k>=kparasta .and. k<=kparaend) then
                    vc(i,j,k)=value
                end if
            end do
        end do
    end do
      
    do k=0,n3
        do j=1,n2
            do i=1,n1
                read(oldresfile_id,*) value
                if (k>=kparasta-1 .and. k<=kparaend) then
                    wc(i,j,k)=value
                end if
      
            end do
        end do
    end do
    close(oldresfile_id)

    return
end
