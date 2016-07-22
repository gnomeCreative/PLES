module orlansky_module

    use,intrinsic :: iso_c_binding

    use myarrays_velo3, only: delu,delv,delw,rhov,u,uc,v,w
    use mysending
    !
    use scala3
    use period
    !
    use mpi

    implicit none

    ! quantita' necessarie per subroutine orlansky ed inflow

    integer,bind(C) :: infout1,infout2,infout3
    integer,bind(C) :: infout4,infout5,infout6

    real,allocatable :: du_dx1(:,:),dv_dx1(:,:),dw_dx1(:,:)
    real,allocatable :: du_dx2(:,:),dv_dx2(:,:),dw_dx2(:,:)

    real,allocatable :: du_dy3(:,:),dv_dy3(:,:),dw_dy3(:,:)
    real,allocatable :: du_dy4(:,:),dv_dy4(:,:),dw_dy4(:,:)

    real,allocatable :: du_dz5(:,:),dv_dz5(:,:),dw_dz5(:,:)
    real,allocatable :: du_dz6(:,:),dv_dz6(:,:),dw_dz6(:,:)

    real,allocatable :: drho_dx1(:,:,:),drho_dx2(:,:,:)
    real,allocatable :: drho_dy3(:,:,:),drho_dy4(:,:,:)
    real,allocatable :: drho_dz5(:,:,:),drho_dz6(:,:,:)

    integer,allocatable :: index_out1(:,:),index_out2(:,:)
    integer,allocatable :: index_out3(:,:),index_out4(:,:)
    integer,allocatable :: index_out5(:,:),index_out6(:,:)

    integer,allocatable :: index_rho1(:,:),index_rho2(:,:)
    integer,allocatable :: index_rho3(:,:),index_rho4(:,:)
    integer,allocatable :: index_rho5(:,:),index_rho6(:,:)
      
    integer :: npn1,npn2,npn3,npn4,npn5,npn6

contains

    subroutine initialize_Orlansky()

        implicit none

        allocate(du_dx1(n2,n3),dv_dx1(n2,n3),dw_dx1(n2,n3))
        allocate(du_dx2(n2,n3),dv_dx2(n2,n3),dw_dx2(n2,n3))

        allocate(du_dy3(n1,n3),dv_dy3(n1,n3),dw_dy3(n1,n3))
        allocate(du_dy4(n1,n3),dv_dy4(n1,n3),dw_dy4(n1,n3))

        allocate(du_dz5(n1,n2),dv_dz5(n1,n2),dw_dz5(n1,n2))
        allocate(du_dz6(n1,n2),dv_dz6(n1,n2),dw_dz6(n1,n2))

        allocate(drho_dx1(nscal,n2,n3),drho_dx2(nscal,n2,n3))
        allocate(drho_dy3(nscal,n1,n3),drho_dy4(nscal,n1,n3))
        allocate(drho_dz5(nscal,n1,n2),drho_dz6(nscal,n1,n2))

        allocate(index_out1(0:n2+1,0:n3+1),index_out2(0:n2+1,0:n3+1))
        allocate(index_out3(0:n1+1,0:n3+1),index_out4(0:n1+1,0:n3+1))
        allocate(index_out5(0:n1+1,0:n2+1),index_out6(0:n1+1,0:n2+1))

        allocate(index_rho1(0:n2+1,0:n3+1),index_rho2(0:n2+1,0:n3+1))
        allocate(index_rho3(0:n1+1,0:n3+1),index_rho4(0:n1+1,0:n3+1))
        allocate(index_rho5(0:n1+1,0:n2+1),index_rho6(0:n1+1,0:n2+1))

        index_out1=1
        index_out2=1

        index_rho1=1
        index_rho2=1

        index_out3=1
        index_out4=1

        index_rho3=1
        index_rho4=1

        index_out5=1
        index_out6=1

        index_rho5=1
        index_rho6=1

    end subroutine initialize_Orlansky

    subroutine orlansky_generale(giac)
        ! compute boundary conditions with orlansky
        !
        implicit none

        !-----------------------------------------------------------------------
        real,intent(in) :: giac(1:n1,1:n2,kparasta:kparaend)
        !-----------------------------------------------------------------------
        integer i,j,k,ii,iii,isc
        integer ierr
        integer,parameter :: i_ob = 1
        real c_orl ,amass,amass_loc
        real vol,vol_loc,c_orlansky(n2)
        real deltax
        !-----------------------------------------------------------------------

        ! initialization
        do k=1,n3
            do j=1,n2
                du_dx1(j,k)=0.
                dv_dx1(j,k)=0.
                dw_dx1(j,k)=0.

                du_dx2(j,k)=0.
                dv_dx2(j,k)=0.
                dw_dx2(j,k)=0.
                do isc=1,nscal
                    drho_dx1(isc,j,k)=0.
                    drho_dx2(isc,j,k)=0.
                end do
            end do
        end do

        do k=1,n3
            do i=1,n1
                du_dy3(i,k)=0.
                dv_dy3(i,k)=0.
                dw_dy3(i,k)=0.

                du_dy4(i,k)=0.
                dv_dy4(i,k)=0.
                dw_dy4(i,k)=0.
                do isc=1,nscal
                    drho_dy3(isc,i,k)=0.
                    drho_dy4(isc,i,k)=0.
                end do
            end do
        end do

        do j=1,n2
            do i=1,n1
                du_dz5(i,j)=0.
                dv_dz5(i,j)=0.
                dw_dz5(i,j)=0.

                du_dz6(i,j)=0.
                dv_dz6(i,j)=0.
                dw_dz6(i,j)=0.
                do isc=1,nscal
                    drho_dz5(isc,i,j)=0.
                    drho_dz6(isc,i,j)=0.
                end do
            end do
        end do


        !-----------------------------------------------------------------------
        !     ORLANSKY
        !-----------------------------------------------------------------------
        if(i_ob == 1)then


            ! side 1 constant csi
            do iii=1,ip
                if(infout1.eq.1)then

                    do k=kparasta,kparaend
                        do j=1,n2

                            du_dx1(j,k)=u(1,j,k)+delu(1,j,k)
                            dv_dx1(j,k)=v(1,j,k)+delv(1,j,k)
                            dw_dx1(j,k)=w(1,j,k)+delw(1,j,k)

                            du_dx1(j,k)=index_out1(j,k)*du_dx1(j,k)
                            dv_dx1(j,k)=index_out1(j,k)*dv_dx1(j,k)
                            dw_dx1(j,k)=index_out1(j,k)*dw_dx1(j,k)
                            do isc=1,nscal
                                drho_dx1(isc,j,k)=rhov(isc,1,j,k)
                                drho_dx1(isc,j,k)=index_out1(j,k)*drho_dx1(isc,j,k)
                            enddo
                        enddo
                    enddo

                endif
                !
                !.......................................................................
                ! side 2 constant csi
                if(infout2.eq.1)then

                    amass_loc = 0.
                    amass = 0.
                    vol_loc = 0.
                    vol = 0.

                    do j=1,n2
                        do k=kparasta,kparaend
                            amass_loc = amass_loc + uc(n1,j,k)
                            vol_loc = vol_loc + giac(n1,j,k)
                        end do
                    end do

                    call MPI_ALLREDUCE(amass_loc,amass,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
                    call MPI_ALLREDUCE(vol_loc,vol,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
                    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

                    !      c_orl = amass/( (z(1,1,jz)-z(1,1,0)) * (y(1,jy,1)-y(1,0,1))  )
                    !      if(myid.eq.0)write(*,*)'C_ORL',c_orl

                    do k=kparasta,kparaend
                        do j=1,n2

                            !      c_orl=c_orlansky(j) !uc(jx,j,k)/giac(jx,j,k)
                            c_orl=uc(n1,j,k)/giac(n1,j,k)

                            if(c_orl.gt.0) then
                                du_dx2(j,k)=u(n1+1,j,k)-2.*dt*c_orl*(u(n1+1,j,k)-u(n1,j,k)) ! /deltax
                                dv_dx2(j,k)=v(n1,j,k)    !v(jx+1,j,k)-2.*dt*c_orl*(v(jx+1,j,k)-v(jx,j,k)) ! /deltax
                                dw_dx2(j,k)=w(n1,j,k)    !w(jx+1,j,k)-2.*dt*c_orl*(w(jx+1,j,k)-w(jx,j,k)) ! /deltax
                            ! drho_dx2(j,k)=rho(jx+1,j,k)  ! -2.*dt*c_orl*(rho(jx+1,j,k)-rho(jx,j,k))
                            else
                                du_dx2(j,k)=u(n1+1,j,k)-2.*dt*c_orl*(u(n1,j,k)-u(n1+1,j,k)) ! /deltax
                                dv_dx2(j,k)=v(n1+1,j,k)-2.*dt*c_orl*(v(n1,j,k)-v(n1+1,j,k)) ! /deltax
                                dw_dx2(j,k)=w(n1+1,j,k)-2.*dt*c_orl*(w(n1,j,k)-w(n1+1,j,k)) ! /deltax
                            ! drho_dx2(j,k)=rho(jx+1,j,k) ! -2.*dt*c_orl*(rho(jx,j,k)-rho(jx+1,j,k))
                            end if

                            du_dx2(j,k)=index_out2(j,k)*du_dx2(j,k)
                            dv_dx2(j,k)=index_out2(j,k)*dv_dx2(j,k)
                            dw_dx2(j,k)=index_out2(j,k)*dw_dx2(j,k)


                            do isc=1,nscal
                                drho_dx2(isc,j,k)=rhov(isc,n1,j,k)
                                drho_dx2(isc,j,k)=index_out2(j,k)*drho_dx2(isc,j,k)
                            end do
                        enddo
                    enddo

                endif
            enddo !fine loop ii=1,ip
            !-----------------------------------------------------------------------
            ! side 3 constant eta
            ! always not periodic
                if(infout3.eq.1)then

                    do k=kparasta,kparaend
                        do i=1,n1

                            du_dy3(i,k)=u(i,1,k)+delu(i,1,k)
                            dv_dy3(i,k)=v(i,1,k)+delv(i,1,k)
                            dw_dy3(i,k)=w(i,1,k)+delw(i,1,k)

                            du_dy3(i,k)=index_out3(i,k)*du_dy3(i,k)
                            dv_dy3(i,k)=index_out3(i,k)*dv_dy3(i,k)
                            dw_dy3(i,k)=index_out3(i,k)*dw_dy3(i,k)

                            do isc=1,nscal
                                drho_dy3(isc,i,k)=rhov(isc,i,1,k)
                                drho_dy3(isc,i,k)=index_out3(i,k)*drho_dy3(isc,i,k)
                            end do
                        enddo
                    enddo

                endif
                !.......................................................................
                ! side 4 constant eta
                if(infout4.eq.1)then

                    do k=kparasta,kparaend
                        do i=1,n1

                            du_dy4(i,k)=u(i,n2,k)+delu(i,n2,k)
                            dv_dy4(i,k)=v(i,n2,k)+delv(i,n2,k)
                            dw_dy4(i,k)=w(i,n2,k)+delw(i,n2,k)


                            du_dy4(i,k)=index_out4(i,k)*du_dy4(i,k)
                            dv_dy4(i,k)=index_out4(i,k)*dv_dy4(i,k)
                            dw_dy4(i,k)=index_out4(i,k)*dw_dy4(i,k)

                            do isc=1,nscal
                                drho_dy4(isc,i,k)=rhov(isc,i,n2,k)
                                drho_dy4(isc,i,k)=index_out4(i,k)*drho_dy4(isc,i,k)
                            end do
                        enddo
                    enddo

                endif
            !-----------------------------------------------------------------------
            ! side 5 constant zita
            do iii=1,kp

                if(myid.eq.0)then
                    if(infout5.eq.1)then

                        do j=1,n2
                            do i=1,n1

                                du_dz5(i,j)=u(i,j,1)+delu(i,j,1)
                                dv_dz5(i,j)=v(i,j,1)+delv(i,j,1)
                                dw_dz5(i,j)=w(i,j,1)+delw(i,j,1)


                                du_dz5(i,j)=index_out5(i,j)*du_dz5(i,j)
                                dv_dz5(i,j)=index_out5(i,j)*dv_dz5(i,j)
                                dw_dz5(i,j)=index_out5(i,j)*dw_dz5(i,j)

                                do isc=1,nscal
                                    drho_dz5(isc,i,j)=rhov(isc,i,j,1)
                                    drho_dz5(isc,i,j)=index_out5(i,j)*drho_dz5(isc,i,j)
                                end do
                            enddo
                        enddo

                    endif
                endif
                !.......................................................................
                ! side 6 constant zita
                if(myid.eq.nproc-1)then
                    if(infout6.eq.1)then

                        do j=1,n2
                            do i=1,n1

                                du_dz6(i,j)=u(i,j,n3)+delu(i,j,n3)
                                dv_dz6(i,j)=v(i,j,n3)+delv(i,j,n3)
                                dw_dz6(i,j)=w(i,j,n3)+delw(i,j,n3)


                                du_dz6(i,j)=index_out6(i,j)*du_dz6(i,j)
                                dv_dz6(i,j)=index_out6(i,j)*dv_dz6(i,j)
                                dw_dz6(i,j)=index_out6(i,j)*dw_dz6(i,j)

                                do isc=1,nscal
                                    drho_dz6(isc,i,j)=rhov(isc,i,j,n3)
                                    drho_dz6(isc,i,j)=index_out6(i,j)*drho_dz6(isc,i,j)
                                end do
                            enddo
                        enddo

                    endif
                endif

            enddo !end loop kk=1,kp


        end if ! i_ob = 1





        !-----------------------------------------------------------------------
        !     FREE SLIP
        !-----------------------------------------------------------------------
        if(i_ob == 0)then

            ! side 1 constant csi
            do iii=1,ip
                if(infout1==1)then

                    do k=kparasta,kparaend
                        do j=1,n2

                            du_dx1(j,k)=u(1,j,k)+delu(1,j,k)
                            dv_dx1(j,k)=v(1,j,k)+delv(1,j,k)
                            dw_dx1(j,k)=w(1,j,k)+delw(1,j,k)

                            du_dx1(j,k)=index_out1(j,k)*du_dx1(j,k)
                            dv_dx1(j,k)=index_out1(j,k)*dv_dx1(j,k)
                            dw_dx1(j,k)=index_out1(j,k)*dw_dx1(j,k)
                            do isc=1,nscal
                                drho_dx1(isc,j,k)=rhov(isc,1,j,k)
                                drho_dx1(isc,j,k)=index_out1(j,k)*drho_dx1(isc,j,k)
                            enddo
                        enddo
                    enddo

                endif
                !
                !.......................................................................
                ! side 2 constant csi
                if(infout2==1)then

                    do k=kparasta,kparaend
                        do j=1,n2

                            du_dx2(j,k)=u(n1,j,k)+delu(n1,j,k)
                            dv_dx2(j,k)=v(n1,j,k)+delv(n1,j,k)
                            dw_dx2(j,k)=w(n1,j,k)+delw(n1,j,k)


                            du_dx2(j,k)=index_out2(j,k)*du_dx2(j,k)
                            dv_dx2(j,k)=index_out2(j,k)*dv_dx2(j,k)
                            dw_dx2(j,k)=index_out2(j,k)*dw_dx2(j,k)

                            do isc=1,nscal
                                drho_dx2(isc,j,k)=rhov(isc,n1,j,k)
                                drho_dx2(isc,j,k)=index_out2(j,k)*drho_dx2(isc,j,k)
                            end do

                        enddo
                    enddo

                endif
            enddo !fine loop ii=1,ip
            !-----------------------------------------------------------------------
            ! side 3 constant eta
            ! alwayws not periodic
                if(infout3==1)then

                    do k=kparasta,kparaend
                        do i=1,n1

                            du_dy3(i,k)=u(i,1,k)+delu(i,1,k)
                            dv_dy3(i,k)=v(i,1,k)+delv(i,1,k)
                            dw_dy3(i,k)=w(i,1,k)+delw(i,1,k)

                            du_dy3(i,k)=index_out3(i,k)*du_dy3(i,k)
                            dv_dy3(i,k)=index_out3(i,k)*dv_dy3(i,k)
                            dw_dy3(i,k)=index_out3(i,k)*dw_dy3(i,k)

                            do isc=1,nscal
                                drho_dy3(isc,i,k)=rhov(isc,i,1,k)
                                drho_dy3(isc,i,k)=index_out3(i,k)*drho_dy3(isc,i,k)
                            end do
                        enddo
                    enddo

                endif
                !.......................................................................
                ! side 4 constant eta
                if(infout4==1)then

                    do k=kparasta,kparaend
                        do i=1,n1

                            du_dy4(i,k)=u(i,n2,k)+delu(i,n2,k)
                            dv_dy4(i,k)=v(i,n2,k)+delv(i,n2,k)
                            dw_dy4(i,k)=w(i,n2,k)+delw(i,n2,k)


                            du_dy4(i,k)=index_out4(i,k)*du_dy4(i,k)
                            dv_dy4(i,k)=index_out4(i,k)*dv_dy4(i,k)
                            dw_dy4(i,k)=index_out4(i,k)*dw_dy4(i,k)

                            do isc=1,nscal
                                drho_dy4(isc,i,k)=rhov(isc,i,n2,k)
                                drho_dy4(isc,i,k)=index_out4(i,k)*drho_dy4(isc,i,k)
                            end do
                        enddo
                    enddo

                endif
            !-----------------------------------------------------------------------
            ! side 5 constant zita
            do iii=1,kp

                if(myid==0)then
                    if(infout5==1)then

                        do j=1,n2
                            do i=1,n1

                                du_dz5(i,j)=u(i,j,1)+delu(i,j,1)
                                dv_dz5(i,j)=v(i,j,1)+delv(i,j,1)
                                dw_dz5(i,j)=w(i,j,1)+delw(i,j,1)


                                du_dz5(i,j)=index_out5(i,j)*du_dz5(i,j)
                                dv_dz5(i,j)=index_out5(i,j)*dv_dz5(i,j)
                                dw_dz5(i,j)=index_out5(i,j)*dw_dz5(i,j)

                                do isc=1,nscal
                                    drho_dz5(isc,i,j)=rhov(isc,i,j,1)
                                    drho_dz5(isc,i,j)=index_out5(i,j)*drho_dz5(isc,i,j)
                                end do
                            enddo
                        enddo

                    endif
                endif
                !.......................................................................
                ! side 6 constant zita
                if(myid==nproc-1)then
                    if(infout6==1)then

                        do j=1,n2
                            do i=1,n1

                                du_dz6(i,j)=u(i,j,n3)+delu(i,j,n3)
                                dv_dz6(i,j)=v(i,j,n3)+delv(i,j,n3)
                                dw_dz6(i,j)=w(i,j,n3)+delw(i,j,n3)


                                du_dz6(i,j)=index_out6(i,j)*du_dz6(i,j)
                                dv_dz6(i,j)=index_out6(i,j)*dv_dz6(i,j)
                                dw_dz6(i,j)=index_out6(i,j)*dw_dz6(i,j)

                                do isc=1,nscal
                                    drho_dz6(isc,i,j)=rhov(isc,i,j,n3)
                                    drho_dz6(isc,i,j)=index_out6(i,j)*drho_dz6(isc,i,j)
                                end do
                            enddo
                        enddo

                    endif
                endif

            enddo !end loop kk=1,kp



        end if !i_ob == 0

        return

    end subroutine orlansky_generale

end module orlansky_module

