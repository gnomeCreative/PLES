module tridag_module

    use mysending
    use scala3
    use period

    use mpi

    implicit none

    private

    ! for transposed tridiag and approximate factorization
    real,allocatable :: g33_tr(:,:,:), giac_tr(:,:,:)
    real,allocatable :: aaa(:,:),rh(:)
    real,allocatable :: aa(:),bb(:),cc(:)

    public :: initialize_tridiag,factorization,tridag

contains

    subroutine initialize_tridiag()

        implicit none

        !-----------------------------------------------------------------------
        ! no idea about what these are used for...
        allocate(aaa(3,n1+n2+n3),rh(n1+n2+n3))
        allocate(aa(n1+n2+n3),bb(n1+n2+n3),cc(n1+n2+n3))

        !-----------------------------------------------------------------------
        ! if semimplict allocation for tridiag of g33 and giac

        iparasta=(myid* int(n1/nproc)  +1)
        iparaend=((myid+1)* int(n1/nproc))

        allocate(giac_tr(n3,n2,iparasta:iparaend))
        giac_tr(:,:,:) = 0.0
        allocate(g33_tr(0:n3,n2,iparasta:iparaend))
        g33_tr(:,:,:) = 0.0

        call set_transpose_implicit(g33_tr,giac_tr)

    end subroutine initialize_tridiag

    subroutine factorization(scalar,ktime,delx,isc)
        !  approximate factorization

        use mysettings, only: pran
        use myarrays_metri3, only: annit
        use myarrays_velo3, only: akapt,akapt_piano

        implicit none

        integer,intent(in) :: ktime
        logical,intent(in) :: scalar
        integer,intent(in) :: isc
        real,intent(inout) :: delx(0:n1+1,0:n2+1,kparasta-1:kparaend+1)

        integer :: i,j,k,ii,jj

        !  upload csi
        do k=kparasta,kparaend
            do j=1,jy
                ! coefficent construction
                if (scalar) then
                    call coed1(j,k,delx,aaa,rh,isc)
                else
                    call coef1_par(j,k,delx,aaa,rh)
                end if
                do ii=1,jx
                    aa(ii)=aaa(1,ii)
                    bb(ii)=aaa(2,ii)
                    cc(ii)=aaa(3,ii)
                end do
                do ii=1,1-ip
                    call triper(aa,bb,cc,rh,jx-1)
                end do
                do ii=1,ip
                    call tridag(aa,bb,cc,rh,jx)
                end do
                ! put out in delu
                do i=1,jx
                    delx(i,j,k)=rh(i)
                end do
            end do
        end do

        !  upload eta
        do k=kparasta,kparaend
            do i=1,jx
                !coefficent construction upload eta
                if (scalar) then
                    call coed2(i,k,delx,aaa,rh,isc)
                else
                    call coef2_par(i,k,delx,aaa,rh)
                end if
                do ii=1,jy
                    aa(ii)=aaa(1,ii)
                    bb(ii)=aaa(2,ii)
                    cc(ii)=aaa(3,ii)
                end do
                do jj=1,1-jp
                    call triper(aa,bb,cc,rh,jy-1)
                end do
                do jj=1,jp
                    call tridag(aa,bb,cc,rh,jy)
                end do
                do j=1,jy
                    delx(i,j,k)=rh(j)          ! put out in delu
                end do
            end do
        end do

        if (scalar) then
            call tridiag_trasp_para_rho(akapt,g33_tr,giac_tr,delx,akapt_piano,pran,isc)
        else
            call tridiag_trasp_para(annit,g33_tr,giac_tr,delx,ktime)!annit_piano,
        end if

        return

    end subroutine factorization

    subroutine set_transpose_implicit(g33_tr,giac_tr)

        use myarrays_metri3, only: annit,g11,g22,g33,giac

        implicit none

        !-----------------------------------------------------------------------
        !     array declaration
        real giac_tr(n3,n2,iparasta:iparaend)
        real g33_tr(0:n3,n2,iparasta:iparaend)
        real, allocatable :: buffer(:),buffer_tot(:)
        real, allocatable :: buffer_plane(:),buffer_plane_loc(:)

        integer count,count_plane
        integer iparasta_tmp,iparaend_tmp
        integer n,m
        integer i,j,k
        integer ierr,status(MPI_STATUS_SIZE)
        !
        !-----------------------------------------------------------------------
        !     TRANSPOSE GIAC
        !-----------------------------------------------------------------------
        do n=0,nproc-1

            iparasta_tmp= (n* int(n1/nproc)  +1)
            iparaend_tmp= ((n+1)* int(n1/nproc))

            count = jy*(iparaend_tmp-iparasta_tmp+1)*(kparaend-kparasta+1)

            if(n==0)then
                allocate(buffer(count))
                buffer = 0.
                allocate(buffer_tot(nproc*count))
                buffer_tot=0.
            else
                buffer     = 0.
                buffer_tot = 0.
            end if

            m = 0
            do k=kparasta,kparaend
                do j=1,jy
                    do i=iparasta_tmp,iparaend_tmp
                        m = m + 1
                        buffer(m) = giac(i,j,k)
                    end do
                end do
            end do

            call MPI_GATHER(buffer(1)    ,count,MPI_REAL_SD,  &
                buffer_tot(1),count,MPI_REAL_SD,  &
                n,MPI_COMM_WORLD,ierr)

            if(myid == n)then

                m = 0
                do k=1,jz
                    do j=1,jy
                        do i=iparasta,iparaend
                            m = m +1
                            giac_tr(k,j,i)=buffer_tot(m)
                        end do
                    end do
                end do

            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        end do

        deallocate(buffer,buffer_tot)
        !-----------------------------------------------------------------------
        !     TRANSPOSE G33
        !-----------------------------------------------------------------------
        do n=0,nproc-1

            iparasta_tmp= (n* int(n1/nproc)  +1)
            iparaend_tmp= ((n+1)* int(n1/nproc))

            count = jy*(iparaend_tmp-iparasta_tmp+1)*(kparaend-kparasta+1)

            if(n==0)then
                allocate(buffer(count))
                buffer = 0.
                allocate(buffer_tot(nproc*count))
                buffer_tot=0.
            else
                buffer     = 0.
                buffer_tot = 0.
            end if

            m = 0
            do k=kparasta,kparaend
                do j=1,jy
                    do i=iparasta_tmp,iparaend_tmp
                        m = m + 1
                        buffer(m) = g33(i,j,k)
                    end do
                end do
            end do

            call MPI_GATHER(buffer(1)    ,count,MPI_REAL_SD,  &
                buffer_tot(1),count,MPI_REAL_SD,  &
                n,MPI_COMM_WORLD,ierr)

            if(myid == n)then

                m = 0
                do k=1,jz
                    do j=1,jy
                        do i=iparasta,iparaend
                            m = m +1
                            g33_tr(k,j,i)=buffer_tot(m)
                        end do
                    end do
                end do

            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        end do

        deallocate(buffer,buffer_tot)

        !.......................................................................
        !     proc 0 send the data
        !      count_plane = jy*jx

        !      allocate(buffer_plane(count_plane))
        !      allocate(buffer_plane_loc(count_plane/nproc))
        !      buffer_plane = 0
        !      buffer_plane_loc = 0

        !      if(myid==0)then
        !      m=0
        !      do i=1,jx
        !      do j=1,jy
        !         m = m + 1
        !         buffer_plane(m) = g33(i,j,0)
        !      end do
        !      end do
        !      end if

        !      call MPI_SCATTER(buffer_plane(1),count_plane/nproc,MPI_REAL_SD,
        !     >        buffer_plane_loc(1),count_plane/nproc,MPI_REAL_SD,
        !     >                 0,MPI_COMM_WORLD,ierr)

        !      m=0
        !      do i=iparasta,iparaend
        !      do j=1,jy
        !         m = m + 1
        !         g33_tr(0,j,i)=buffer_plane_loc(m)
        !      end do
        !      end do

        !      deallocate(buffer_plane,buffer_plane_loc)
        !.......................................................................
        !     proc 0 send the data
        !      g33_tr = 0.

        if(myid==0)then
            do j=1,jy
                do i=1,iparasta,iparaend !jx
                    m = m + 1
                    g33_tr(0,j,i) = g33(i,j,0)
                end do
            end do
        end if

        count_plane = jy*(iparaend-iparasta+1)

        allocate(buffer_plane(count_plane))
        allocate(buffer_plane_loc(count_plane))
        buffer_plane = 0
        buffer_plane_loc = 0

        do n=1,nproc-1
            iparasta_tmp= (n* int(n1/nproc)  +1)
            iparaend_tmp= ((n+1)* int(n1/nproc))

            buffer_plane = 0
            buffer_plane_loc = 0

            if(myid==0)then
                m=0
                do j=1,jy
                    do i=1,iparasta_tmp,iparaend_tmp !jx
                        m = m + 1
                        buffer_plane(m) = g33(i,j,0)
                    end do
                end do
            end if


            if(myid == 0)then
                call MPI_SEND(buffer_plane(1),count_plane,MPI_REAL_SD, &
                    n,1001,MPI_COMM_WORLD,ierr)

            elseif(myid==n)then
                call MPI_RECV(buffer_plane_loc(1),count_plane,MPI_REAL_SD, &
                    0,1001,MPI_COMM_WORLD,status,ierr)
            endif


            call MPI_BARRIER(MPI_COMM_WORLD,ierr)


            if(myid==n)then
                m=0
                do j=1,jy
                    do i=iparasta_tmp,iparaend_tmp
                        m = m + 1
                        g33_tr(0,j,i)=buffer_plane_loc(m)
                    end do
                end do
            end if

        end do

        deallocate(buffer_plane,buffer_plane_loc)

        return
    end subroutine set_transpose_implicit

    subroutine tridag(aa,bb,cc,ff,n)

        ! solution for a tridiagonal system

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer n1,n2,k,k1,n,n1n
        real aa(*),bb(*),cc(*),ff(*)
        !-----------------------------------------------------------------------
        n1=1
        bb(n1)=1./bb(n1)
        aa(n1)=ff(n1)*bb(n1)
        n2=n1+1
        n1n=n1+n
        do k=n2,n
            k1=k-1
            cc(k1)=cc(k1)*bb(k1)
            bb(k)=bb(k)-aa(k)*cc(k1)
            bb(k)=1./bb(k)
            aa(k)=(ff(k)-aa(k)*aa(k1))*bb(k)
        end do
        !
        ! back substitution
        !
        ff(n)=aa(n)
        do k1=n2,n
            k=n1n-k1
            ff(k)=aa(k)-cc(k)*ff(k+1)
        end do
        return
    end subroutine tridag

    subroutine tridiag_trasp_para_rho(akapt,g33_tr,giac_tr,del,akapt_piano,pran,isc)
        !*************************************************************************
        ! transpose operation for matrix.
        ! each plane XZ it is seen as two dimensional matrix and it is transposed,
        ! in this way the solution of the tridiagonal matrix is made in parallel in
        ! i instead of k, then it will be necesseary to transpose again the result
        ! (it is necessary to transpose the matrix used in "coed3", so g33, akapt,
        ! giac)

        implicit none
        !

        !-----------------------------------------------------------------------
        !     array declaration
        real pran(nscal)
        real akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        real akapt_piano(nscal,0:n1+1,0:n2+1,1)

        real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
        real aa(n1+n2+n3),bb(n1+n2+n3),cc(n1+n2+n3)
        real aaa(3,n1+n2+n3),rh(n1+n2+n3)

        real akapt_tr(0:n3+1,0:n2+1,iparasta-1:iparaend+1)
        real g33_tr(0:n3,n2,iparasta:iparaend) !n1)
        real giac_tr(n3,n2,iparasta:iparaend) !n1)

        !     this is not really a transposed matrix at the first definition
        !    it is a block transposed, then with ALLTOALL there is the final transposed
        real del_tr(0:n3+1,0:n2+1,iparasta:iparaend)

        integer i,j,k,it,jt,kt,m,ii,isc
        integer ncoljx,ncoljz,n
        integer lll

        real,allocatable :: del_trcols(:),del_trcolr(:)
        real,allocatable :: del_cols(:),del_colr(:)

        integer ierr
        integer nn,icheck
        integer iparasta_tmp,iparaend_tmp,count,count_plane
        !-----------------------------------------------------------------------
        real, allocatable :: buffer(:),buffer_tot(:)
        real, allocatable :: buffer_plane(:),buffer_plane_loc(:)
        integer ifaccio
        !-----------------------------------------------------------------------
        ! transposed for the matrix


        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !     AKAPT TRANSPOSE
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------

        akapt_tr = 1./re/pran(isc)

        do n=0,nproc-1
            iparasta_tmp= (n* int(n1/nproc)  +1)
            iparaend_tmp= ((n+1)* int(n1/nproc))

            count=jy*(iparaend_tmp-iparasta_tmp+1+2)*(kparaend-kparasta+1)
            count_plane = jy*(iparaend_tmp-iparasta_tmp+1)

            if(n==0)then
                allocate(buffer(count))
                buffer = 0.
                allocate(buffer_tot(nproc*count))
                buffer_tot = 0.
            else
                buffer = 0.
                buffer_tot = 0.
            end if

            m = 0
            do k=kparasta,kparaend
                do j=1,jy
                    do i=iparasta_tmp-1,iparaend_tmp+1
                        m = m + 1
                        buffer(m) = akapt(isc,i,j,k)
                    end do
                end do
            end do

            call MPI_GATHER(buffer(1)    ,count,MPI_REAL_SD,  &
                buffer_tot(1),count,MPI_REAL_SD,  &
                n,MPI_COMM_WORLD,ierr)

            if(myid==n)then
                m = 0
                do k=1,jz
                    do j=1,jy
                        do i=iparasta-1,iparaend+1
                            m = m +1
                            akapt_tr(k,j,i) = buffer_tot(m)
                        end do
                    end do
                end do
            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do

        !..................................................................

        !     proc 0 send the data
        count_plane = jy*jx

        allocate(buffer_plane(count_plane))
        allocate(buffer_plane_loc(count_plane/nproc))
        buffer_plane = 0
        buffer_plane_loc = 0

        if(myid==0)then
            m=0
            do i=1,jx
                do j=1,jy
                    m = m + 1
                    buffer_plane(m) = akapt(isc,i,j,0)
                end do
            end do
        end if

        call MPI_SCATTER(buffer_plane(1),count_plane/nproc,MPI_REAL_SD, &
            buffer_plane_loc(1),count_plane/nproc,MPI_REAL_SD, &
            0,MPI_COMM_WORLD,ierr)

        m=0
        do i=iparasta,iparaend
            do j=1,jy
                m = m + 1
                akapt_tr(0,j,i)=buffer_plane_loc(m)
            end do
        end do

        deallocate(buffer_plane,buffer_plane_loc)

        !     nproc-1 send the data
        count_plane = jy*jx

        allocate(buffer_plane(count_plane))
        allocate(buffer_plane_loc(count_plane/nproc))
        buffer_plane = 0
        buffer_plane_loc = 0

        if(myid==nproc-1)then
            m=0
            do i=1,jx
                do j=1,jy
                    m = m + 1
                    buffer_plane(m) = akapt(isc,i,j,jz+1)
                end do
            end do
        end if

        call MPI_SCATTER(buffer_plane(1),count_plane/nproc,MPI_REAL_SD, &
            buffer_plane_loc(1),count_plane/nproc,MPI_REAL_SD, &
            nproc-1,MPI_COMM_WORLD,ierr)

        m=0
        do i=iparasta,iparaend
            do j=1,jy
                m = m + 1
                akapt_tr(jz+1,j,i)=buffer_plane_loc(m)
            end do
        end do



        deallocate(buffer_plane,buffer_plane_loc)
        deallocate(buffer,buffer_tot)
        !
        !----------------------------------------------------------------------
        !     DEL TRANSPOSE
        !----------------------------------------------------------------------
        !
        del_tr = 0.

        do n=0,nproc-1
            iparasta_tmp= (n* int(n1/nproc)  +1)
            iparaend_tmp= ((n+1)* int(n1/nproc))

            count = jy*(iparaend_tmp-iparasta_tmp+1)*(kparaend-kparasta+1)
            count_plane = jy*(iparaend_tmp-iparasta_tmp+1)

            if(n==0)then
                allocate(buffer(count))
                buffer = 0.
                allocate(buffer_tot(nproc*count))
                buffer_tot = 0.
            else
                buffer = 0.
                buffer_tot = 0.
            end if

            m = 0
            do k=kparasta,kparaend
                do j=1,jy
                    do i=iparasta_tmp,iparaend_tmp
                        m = m + 1
                        buffer(m) = del(i,j,k)
                    end do
                end do
            end do

            call MPI_GATHER(buffer(1)    ,count,MPI_REAL_SD,  &
                buffer_tot(1),count,MPI_REAL_SD,  &
                n,MPI_COMM_WORLD,ierr)

            if(myid==n)then
                m = 0
                do k=1,jz
                    do j=1,jy
                        do i=iparasta,iparaend
                            m = m +1

                            del_tr(k,j,i) = buffer_tot(m)
                        end do
                    end do
                end do
            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do

        deallocate(buffer,buffer_tot)

        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        ncoljx=jx/nproc
        ncoljz=jz/nproc
        ! now del_tr is ready to be accepted in the tridiag system resolution,
        ! del_tr is subdivided according to 1st coordinate
        !
        !
        ! 5 - parallel computation for the tridiagonal system
        !
        !
        !     take up zita ---> in reality is csi
        !
        do kt=1+myid*ncoljx,(myid+1)*ncoljx
            do jt=1,jy
                !
                call coed3_tr(jt,kt,akapt_tr,del_tr, &
                    aaa,rh,myid,g33_tr,giac_tr,iparasta,iparaend)
                !
                do ii=1,jz
                    aa(ii)=aaa(1,ii)
                    bb(ii)=aaa(2,ii)
                    cc(ii)=aaa(3,ii)
                end do
                do ii=1,1-kp
                    call triper(aa,bb,cc,rh,jz-1)
                end do
                do ii=1,kp
                    call tridag(aa,bb,cc,rh,jz)
                end do
                do it=1,jz
                    del_tr(it,jt,kt)=rh(it)          ! put output in del_tr
                end do

            end do
        end do


        !
        ! 6 - now transpose again the computed quantities
        !
        allocate (del_cols(jx*jy*ncoljz))
        allocate (del_colr(jx*jy*ncoljz))

        do n=1,nproc

            do kt=1+myid*ncoljx,(myid+1)*ncoljx
                do jt=1,jy
                    do it=1+(n-1)*ncoljz,n*ncoljz

                        k=kt+(n-1-myid)*ncoljx
                        j=jt
                        i=it-(n-1-myid)*ncoljz
                        del(k,j,i)=del_tr(it,jt,kt)

                    enddo
                enddo
            enddo
            !
            ! 7 - local allocation for re-transposed elements in vector
            !

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1+(n-1)*ncoljx,n*ncoljx

                        m=i-(n-1)*ncoljx+ncoljx*(j-1+jy*(k-kparasta))
                        m=m+(n-1)*jy*ncoljx*ncoljz
                        del_cols(m)=del(i,j,k)

                    enddo
                enddo
            enddo

        enddo

        !
        ! 8 - call ALLTOALL for exchange between procs
        !
        call MPI_ALLTOALL(del_cols,ncoljx*jy*ncoljz,MPI_REAL_SD, &
            del_colr,ncoljx*jy*ncoljz,MPI_REAL_SD, &
            MPI_COMM_WORLD,ierr)
        !
        ! 9 - re construct matrix exchanged between procs
        !
        do n=1,nproc

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1+(n-1)*ncoljx,n*ncoljx

                        m=i-(n-1)*ncoljx+ncoljx*(j-1+jy*(k-kparasta))
                        m=m+(n-1)*jy*ncoljx*ncoljz
                        del(i,j,k)=del_colr(m)


                    enddo
                enddo
            enddo

        enddo


        deallocate(del_cols,del_colr)
        !
        ! now del is ready for the computation of rho-u-v-w
        !
        return
    end subroutine tridiag_trasp_para_rho

    subroutine tridiag_trasp_para(akapt,g33_tr,giac_tr,del,ktime)
        !***********************************************************************
        ! transpose operation for matrix.
        ! each plane XZ it is seen as two dimensional matrix and it is transposed,
        ! in this way the solution of the tridiagonal matrix is made in parallel in
        ! i instead of k, then it will be necesseary to transpose again the result
        ! (it is necessary to transpose the matrix used in "coed3", so g33, akapt,
        ! giac)
        !
        implicit none


        !-----------------------------------------------------------------------
        !     array declaration
        real akapt(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
        !      real annit_piano(0:n1+1,0:n2+1,1)

        real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
        real aa(n1+n2+n3),bb(n1+n2+n3),cc(n1+n2+n3)
        real aaa(3,n1+n2+n3),rh(n1+n2+n3)


        real akapt_tr(0:n3+1,0:n2+1,iparasta-1:iparaend+1)
        real g33_tr(0:n3,n2,iparasta:iparaend)
        real giac_tr(n3,n2,iparasta:iparaend) ! n1)

        !     this is not really a transposed matrix at the first definition
        !    it is a block transposed, then with ALLTOALL there is the final transposed

        real del_tr(0:n3+1,0:n2+1,iparasta:iparaend)


        integer i,j,k,it,jt,kt,m,ii,iii
        integer ncoljx,ncoljz,n
        integer lll

        real,allocatable :: del_trcols(:),del_trcolr(:)
        real,allocatable :: del_cols(:),del_colr(:)

        integer ierr

        integer nn,ktime
        integer iparasta_tmp,iparaend_tmp,count,count_plane
        !-----------------------------------------------------------------------
        real, allocatable :: buffer(:),buffer_tot(:)
        real, allocatable :: buffer_plane(:),buffer_plane_loc(:)
        integer ifaccio
        !-----------------------------------------------------------------------
        !
        !

        !-----------------------------------------------------------------------
        ! 0 - transposed for the matrix
        !-----------------------------------------------------------------------


        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !     AKAPT TRANSPOSE
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------

        akapt_tr = 1/re

        do n=0,nproc-1
            iparasta_tmp= (n* int(n1/nproc)  +1)
            iparaend_tmp= ((n+1)* int(n1/nproc))

            count=jy*(iparaend_tmp-iparasta_tmp+1+2)*(kparaend-kparasta+1)
            count_plane = jy*(iparaend_tmp-iparasta_tmp+1)

            if(n==0)then
                allocate(buffer(count))
                buffer = 0.
                allocate(buffer_tot(nproc*count))
                buffer_tot = 0.
            else
                buffer = 0.
                buffer_tot = 0.
            end if

            m = 0
            do k=kparasta,kparaend
                do j=1,jy
                    do i=iparasta_tmp-1,iparaend_tmp+1
                        m = m + 1
                        buffer(m) = akapt(i,j,k)
                    end do
                end do
            end do

            call MPI_GATHER(buffer(1)    ,count,MPI_REAL_SD,  &
                buffer_tot(1),count,MPI_REAL_SD,  &
                n,MPI_COMM_WORLD,ierr)

            if(myid==n)then
                m = 0
                do k=1,jz
                    do j=1,jy
                        do i=iparasta-1,iparaend+1
                            m = m +1
                            akapt_tr(k,j,i) = buffer_tot(m)
                        end do
                    end do
                end do
            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do

        !..................................................................

        !     proc 0 send the data
        count_plane = jy*jx

        allocate(buffer_plane(count_plane))
        allocate(buffer_plane_loc(count_plane/nproc))
        buffer_plane = 0
        buffer_plane_loc = 0

        if(myid==0)then
            m=0
            do i=1,jx
                do j=1,jy
                    m = m + 1
                    buffer_plane(m) = akapt(i,j,0)
                end do
            end do
        end if

        call MPI_SCATTER(buffer_plane(1),count_plane/nproc,MPI_REAL_SD, &
            buffer_plane_loc(1),count_plane/nproc,MPI_REAL_SD, &
            0,MPI_COMM_WORLD,ierr)

        m=0
        do i=iparasta,iparaend
            do j=1,jy
                m = m + 1
                akapt_tr(0,j,i)=buffer_plane_loc(m)
            end do
        end do

        deallocate(buffer_plane,buffer_plane_loc)

        !     nproc-1 send the data
        count_plane = jy*jx

        allocate(buffer_plane(count_plane))
        allocate(buffer_plane_loc(count_plane/nproc))
        buffer_plane = 0
        buffer_plane_loc = 0

        if(myid==nproc-1)then
            m=0
            do i=1,jx
                do j=1,jy
                    m = m + 1
                    buffer_plane(m) = akapt(i,j,jz+1)
                end do
            end do
        end if

        call MPI_SCATTER(buffer_plane(1),count_plane/nproc,MPI_REAL_SD, &
            buffer_plane_loc(1),count_plane/nproc,MPI_REAL_SD, &
            nproc-1,MPI_COMM_WORLD,ierr)

        m=0
        do i=iparasta,iparaend
            do j=1,jy
                m = m + 1
                akapt_tr(jz+1,j,i)=buffer_plane_loc(m)
            end do
        end do



        deallocate(buffer_plane,buffer_plane_loc)
        deallocate(buffer,buffer_tot)
        !
        !----------------------------------------------------------------------
        !     DEL TRANSPOSE
        !----------------------------------------------------------------------
        !
        del_tr = 0.

        do n=0,nproc-1
            iparasta_tmp= (n* int(n1/nproc)  +1)
            iparaend_tmp= ((n+1)* int(n1/nproc))

            count = jy*(iparaend_tmp-iparasta_tmp+1)*(kparaend-kparasta+1)
            count_plane = jy*(iparaend_tmp-iparasta_tmp+1)

            if(n==0)then
                allocate(buffer(count))
                buffer = 0.
                allocate(buffer_tot(nproc*count))
                buffer_tot = 0.
            else
                buffer = 0.
                buffer_tot = 0.
            end if

            m = 0
            do k=kparasta,kparaend
                do j=1,jy
                    do i=iparasta_tmp,iparaend_tmp
                        m = m + 1
                        buffer(m) = del(i,j,k)
                    end do
                end do
            end do

            call MPI_GATHER(buffer(1)    ,count,MPI_REAL_SD,  &
                buffer_tot(1),count,MPI_REAL_SD,  &
                n,MPI_COMM_WORLD,ierr)

            if(myid==n)then
                m = 0
                do k=1,jz
                    do j=1,jy
                        do i=iparasta,iparaend
                            m = m +1

                            del_tr(k,j,i) = buffer_tot(m)
                        end do
                    end do
                end do
            end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do

        deallocate(buffer,buffer_tot)



        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !
        ncoljx=jx/nproc
        ncoljz=jz/nproc

        ! now del_tr is ready to be accepted in the tridiag system resolution,
        ! del_tr is subdivided according to 1st coordinate
        !
        !
        ! 5 - parallel computation for the tridiagonal system
        !
        !
        !     take up zita ---> in reality is csi
        !
        do kt=1+myid*ncoljx,(myid+1)*ncoljx

            do jt=1,jy
                !
                call coed3_tr(jt,kt,akapt_tr,del_tr, &
                    aaa,rh,myid,g33_tr,giac_tr,iparasta,iparaend)
                !
                do ii=1,jz
                    aa(ii)=aaa(1,ii)
                    bb(ii)=aaa(2,ii)
                    cc(ii)=aaa(3,ii)
                end do
                do ii=1,1-kp
                    call triper(aa,bb,cc,rh,jz-1)
                end do
                do ii=1,kp
                    call tridag(aa,bb,cc,rh,jz)
                end do
                do it=1,jz
                    del_tr(it,jt,kt)=rh(it)      ! put output in del_tr
                end do

            end do
        end do


        !
        ! 6 - now transpose again the computed quantities
        !
        allocate (del_cols(jx*jy*ncoljz))
        allocate (del_colr(jx*jy*ncoljz))

        do n=1,nproc

            do kt=1+myid*ncoljx,(myid+1)*ncoljx
                do jt=1,jy
                    do it=1+(n-1)*ncoljz,n*ncoljz

                        k=kt+(n-1-myid)*ncoljx
                        j=jt
                        i=it-(n-1-myid)*ncoljz

                        del(k,j,i)=del_tr(it,jt,kt)
                    enddo
                enddo
            enddo
            !
            ! 7 - local allocation for re-transposed elements in vector
            !

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1+(n-1)*ncoljx,n*ncoljx

                        m=i-(n-1)*ncoljx+ncoljx*(j-1+jy*(k-kparasta))
                        m=m+(n-1)*jy*ncoljx*ncoljz
                        del_cols(m)=del(i,j,k)

                    enddo
                enddo
            enddo

        enddo

        !
        ! 8 - call ALLTOALL for exchange between procs
        !
        call MPI_ALLTOALL(del_cols,ncoljx*jy*ncoljz,MPI_REAL_SD, &
            del_colr,ncoljx*jy*ncoljz,MPI_REAL_SD, &
            MPI_COMM_WORLD,ierr)
        !
        ! 9 - re construct matrix exchanged between procs
        !
        do n=1,nproc

            do k=kparasta,kparaend
                do j=1,jy
                    do i=1+(n-1)*ncoljx,n*ncoljx

                        m=i-(n-1)*ncoljx+ncoljx*(j-1+jy*(k-kparasta))
                        m=m+(n-1)*jy*ncoljx*ncoljz
                        del(i,j,k)=del_colr(m)

                    enddo
                enddo
            enddo

        enddo

        deallocate(del_cols,del_colr)
        !
        ! now del is ready for the computation of rho-u-v-w
        !
        return
    end subroutine tridiag_trasp_para

    subroutine triper(aa,bb,cc,ff,n)
        !***********************************************************************
        ! subroutine from book ...
        !
        ! aa,bb,cc coefficent vectors: aa(i)=a(1,i)
        !                              bb(i)=a(2,i)
        ! the periodic term up/right is in aa(N1)
        ! the periodic term bottom/left is in cc(N+1)
        ! ff is the vector with the known terms, then transformed
        ! in the solution vector.
        ! aa,bb,cc,ff,gam2 are vectors with dimension N+2 to be specified in the
        ! calling program
        !
        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer k,k1,k2,n,np1,np2,n1n
        real zaa
        real aa(*),bb(*),cc(*),ff(*)
        real gam2(n1+n2+n3)
        !-----------------------------------------------------------------------
        np1=1
        bb(np1)=1./bb(np1)
        gam2(np1)=-aa(np1)*bb(np1)
        aa(np1)=ff(np1)*bb(np1)
        np2=np1+1
        n1n=np1+n
        do k=np2,n
            k1=k-1
            cc(k1)=cc(k1)*bb(k1)
            bb(k)=bb(k)-aa(k)*cc(k1)
            bb(k)=1./bb(k)
            gam2(k)=-aa(k)*gam2(k1)*bb(k)
            aa(k)=(ff(k)-aa(k)*aa(k1))*bb(k)
        end do
        gam2(n)=gam2(n)-cc(n)*bb(n)
        !
        ! back substitution
        !
        ff(n)=aa(n)
        bb(n)=gam2(n)
        do k1=np2,n
            k=n1n-k1
            k2=k+1
            ff(k)=aa(k)-cc(k)*ff(k2)
            bb(k)=gam2(k)-cc(k)*bb(k2)
        end do
        !
        k1=n+1
        zaa=ff(k1)-cc(k1)*ff(np1)-aa(k1)*ff(n)
        zaa=zaa/(bb(k1)+aa(k1)*bb(n)+cc(k1)*bb(np1))
        ff(k1)=zaa
        do k=np1,n
            ff(k)=ff(k)+bb(k)*zaa
        end do    !
        ff(n+2)=ff(np1)
        return
    end subroutine triper

    subroutine coed1(j,k,del,aa,rh,isc)
        !***********************************************************************
        ! compute the coefficents for band tridiagonal matrix for the solution
        ! along csi of the scalar eq.
        !
        use myarrays_metri3, only: annit,g11,g22,g33,giac
        use myarrays_velo3, only: akapt,akapt_piano,rhs

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,ii,isc
        real aa(3,*),rh(*)
        real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !0:n3+1)
        !-----------------------------------------------------------------------
        !
        do ii=1,ip
            !     computation on left side
            i=1
            !
            aa(1,i)=0.
            !
            aa(2,i)=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)+ &
                3.*akapt(isc,i-1,j,k)*g11(i-1,j,k)
            aa(2,i)=1.+aa(2,i)*dt/giac(i,j,k)/2.
            !
            aa(3,i)=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)+ &
                akapt(isc,i-1,j,k)*g11(i-1,j,k)/3.
            aa(3,i)=-dt*aa(3,i)/giac(i,j,k)/2.
            !
            !
            rh(i)=rhs(i,j,k)
            !
            !     computation on right side
            i=jx
            !
            aa(1,i)=akapt(isc,i+1,j,k)*g11(i,j,k)/3.+ &
                .5*(akapt(isc,i,j,k)+akapt(isc,i-1,j,k))*g11(i-1,j,k)
            aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
            !
            aa(2,i)=akapt(isc,i+1,j,k)*g11(i,j,k)*3.+ &
                .5*(akapt(isc,i,j,k)+akapt(isc,i-1,j,k))*g11(i-1,j,k)
            aa(2,i)=1.+dt*aa(2,i)/giac(i,j,k)/2.
            !
            aa(3,i)=0.
            !
            rh(i)=rhs(i,j,k)
        !
        enddo
        !
        !     computation into the field
        do i=1+ip,jx-ip
            !
            aa(1,i)=.5*(akapt(isc,i,j,k)+akapt(isc,i-1,j,k))*g11(i-1,j,k)
            aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
            !
            aa(3,i)=.5*(akapt(isc,i+1,j,k)+akapt(isc,i,j,k))*g11(i,j,k)
            aa(3,i)=-dt*aa(3,i)/giac(i,j,k)/2.
            !
            aa(2,i)=1.-aa(1,i)-aa(3,i)
            !
            rh(i)=rhs(i,j,k)
        !
        enddo
        !
        return
    end subroutine coed1

    subroutine coed2(i,k,del,aa,rh,isc)
        !***********************************************************************
        ! compute the coefficents for band tridiagonal matrix for the solution
        ! along eta of the scalar eq.
        !
        use myarrays_metri3, only: annit,g11,g22,g33,giac
        use myarrays_velo3, only: akapt,akaptV

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,jj,isc
        real aa(3,*),rh(*)
        real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !0:n3+1)
        !-----------------------------------------------------------------------
        !
        !
        do jj=1,jp
            !     computation on bottom side
            j=1
            !
            aa(1,j)=0.
            !
            aa(2,j)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)+ &
                3.*akaptV(isc,i,j-1,k)*g22(i,j-1,k)
            aa(2,j)=1.+aa(2,j)*dt/giac(i,j,k)/2.
            !
            aa(3,j)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)+ &
                akaptV(isc,i,j-1,k)*g22(i,j-1,k)/3.
            aa(3,j)=-dt*aa(3,j)/giac(i,j,k)/2.
            !
            rh(j)=del(i,j,k)
            !
            !     computation on upper side
            j=jy
            !
            aa(1,j)=akaptV(isc,i,j+1,k)*g22(i,j,k)/3.+ &
                .5*(akaptV(isc,i,j,k)+akaptV(isc,i,j-1,k))*g22(i,j-1,k)
            aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
            !
            aa(2,j)=akaptV(isc,i,j+1,k)*g22(i,j,k)*3.+ &
                .5*(akaptV(isc,i,j,k)+akaptV(isc,i,j-1,k))*g22(i,j-1,k)
            aa(2,j)=1.+dt*aa(2,j)/giac(i,j,k)/2.
            !
            aa(3,j)=0.
            rh(j)=del(i,j,k)
        !
        enddo
        !
        !     computation into the field
        do j=1+jp,jy-jp
            !
            aa(1,j)=.5*(akaptV(isc,i,j,k) &
                +akaptV(isc,i,j-1,k))*g22(i,j-1,k)
            aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
            !
            aa(3,j)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)
            aa(3,j)=-dt*aa(3,j)/giac(i,j,k)/2.
            !
            aa(2,j)=1.-aa(1,j)-aa(3,j)
            !
            rh(j)=del(i,j,k)
        !
        enddo
        !
        return
    end subroutine coed2

    subroutine coed3_tr(j,k,aka_tr,del,aa,rh,myid,g33_tr,giac_tr,iparasta,iparaend)
        !***********************************************************************
        ! compute the coefficent for the banded tridiagonal matrix
        ! for solution in zita

        implicit none
        !
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,kk,myid
        integer iparasta,iparaend
        real aa(3,*),rh(*)
        real aka_tr(0:n3+1,0:n2+1,iparasta-1:iparaend+1)
        real del(0:n3+1,0:n2+1,iparasta:iparaend)
        real g33_tr(0:n3,n2,iparasta:iparaend) !n1)
        real giac_tr(n3,n2,iparasta:iparaend) !n1)
        !-----------------------------------------------------------------------
        !
        !     inside the field
        !
        do i=1,jz
            !
            aa(1,i)=.5*(aka_tr(i,j,k)+aka_tr(i-1,j,k))*g33_tr(i-1,j,k)
            aa(1,i)=-dt*aa(1,i)/giac_tr(i,j,k)/2.
            !
            aa(3,i)=.5*(aka_tr(i,j,k)+aka_tr(i+1,j,k))*g33_tr(i,j,k)
            aa(3,i)=-dt*aa(3,i)/giac_tr(i,j,k)/2.
            !
            aa(2,i)=1.-aa(1,i)-aa(3,i)
            !
            rh(i)=del(i,j,k)
        !
        enddo
        !
        return
    end subroutine coed3_tr

    subroutine coef1_par(j,k,del,aa,rh)
        !***********************************************************************
        ! compute the coefficent of the band tridiagonal matrix for solution
        ! along csi for equation u,v,  w
        !
        use myarrays_metri3, only: annit,g11,g22,g33,giac
        use myarrays_velo3, only: akapt,akapt_piano,rhs

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,ii
        real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1),aa(3,*),rh(*)
        !-----------------------------------------------------------------------
        !
        do ii=1,ip
            !     side left
            i=1
            !
            aa(1,i)=0.
            !
            aa(2,i)=.5*(annit(i,j,k)+annit(i+1,j,k))*g11(i,j,k)+ &
                3.*annit(i-1,j,k)*g11(i-1,j,k)
            aa(2,i)=1.+aa(2,i)*dt/giac(i,j,k)/2.
            !
            aa(3,i)=.5*(annit(i,j,k)+annit(i+1,j,k))*g11(i,j,k)+ &
                annit(i-1,j,k)*g11(i-1,j,k)/3.
            aa(3,i)=-dt*aa(3,i)/giac(i,j,k)/2.
            !
            rh(i)=rhs(i,j,k)+ &
                annit(i-1,j,k)*g11(i-1,j,k)* &
                del(i-1,j,k)*8.*dt/giac(i,j,k)/6.
            !
            !     side right
            i=jx
            !
            aa(1,i)=annit(i+1,j,k)*g11(i,j,k)/3.+ &
                .5*(annit(i,j,k)+annit(i-1,j,k))*g11(i-1,j,k)
            aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
            !
            aa(2,i)=annit(i+1,j,k)*g11(i,j,k)*3.+ &
                .5*(annit(i,j,k)+annit(i-1,j,k))*g11(i-1,j,k)
            aa(2,i)=1.+dt*aa(2,i)/giac(i,j,k)/2.
            !
            aa(3,i)=0.
            !
            rh(i)=rhs(i,j,k)+ &
                annit(i+1,j,k)*g11(i,j,k)* &
                del(i+1,j,k)*8.*dt/giac(i,j,k)/6.
        !
        enddo
        !
        !     into the field
        do i=1+ip,jx-ip
            !
            aa(1,i)=.5*(annit(i,j,k)+annit(i-1,j,k))*g11(i-1,j,k)
            aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
            !
            aa(3,i)=.5*(annit(i+1,j,k)+annit(i,j,k))*g11(i,j,k)
            aa(3,i)=-dt*aa(3,i)/giac(i,j,k)/2.
            !
            aa(2,i)=1.-aa(1,i)-aa(3,i)
            !
            rh(i)=rhs(i,j,k)
        !
        enddo
        !
        return
    end subroutine coef1_par

    subroutine coef2_par(i,k,del,aa,rh)
        !***********************************************************************
        ! compute the coefficent of the band tridiagonal matrix for solution
        ! along eta for equations in u,v,w
        !
        use myarrays_metri3, only: annitV,g11,g22,g33,giac
        use myarrays_velo3, only: akapt,akapt_piano,rhs

        implicit none
        !-----------------------------------------------------------------------
        !     array declaration
        integer i,j,k,jj
        real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1),aa(3,*),rh(*)
        !-----------------------------------------------------------------------
        !
        do jj=1,jp
            !
            !     side bottom
            j=1
            !
            aa(1,j)=0.
            !
            aa(2,j)=.5*(annitV(i,j,k)+annitV(i,j+1,k))*g22(i,j,k) &
                + 3.*annitV(i,j-1,k)*g22(i,j-1,k)
            aa(2,j)=1.+aa(2,j)*dt/giac(i,j,k)/2.
            !
            aa(3,j)=.5*(annitV(i,j,k)+annitV(i,j+1,k))*g22(i,j,k) &
                + annitV(i,j-1,k)*g22(i,j-1,k)/3.
            aa(3,j)=-dt*aa(3,j)/giac(i,j,k)/2.
            !
            rh(j)=del(i,j,k)+ &
                annitV(i,j-1,k)*g22(i,j-1,k)* &
                del(i,j-1,k)*8.*dt/giac(i,j,k)/6.
            !
            !     side upper
            j=jy
            !
            aa(1,j)=annitV(i,j+1,k)*g22(i,j,k)/3.+ &
                .5*(annitV(i,j,k)+annitV(i,j-1,k))*g22(i,j-1,k)
            aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
            !
            aa(2,j)=annitV(i,j+1,k)*g22(i,j,k)*3.+ &
                .5*(annitV(i,j,k)+annitV(i,j-1,k))*g22(i,j-1,k)
            aa(2,j)=1.+dt*aa(2,j)/giac(i,j,k)/2.
            !
            aa(3,j)=0.
            !
            rh(j)=del(i,j,k)+ &
                annitV(i,j+1,k)*g22(i,j,k) &
                *del(i,j+1,k)*8.*dt/giac(i,j,k)/6.
        !
        enddo
        !
        !     into the field
        do j=1+jp,jy-jp
            !
            aa(1,j)=.5*(annitV(i,j,k) &
                +annitV(i,j-1,k))*g22(i,j-1,k)
            aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
            !
            aa(3,j)=.5*(annitV(i,j,k) &
                +annitV(i,j+1,k))*g22(i,j,k)
            aa(3,j)=-dt*aa(3,j)/giac(i,j,k)/2.
            !
            aa(2,j)=1.-aa(1,j)-aa(3,j)
            !
            rh(j)=del(i,j,k)
        !
        enddo
        !
        return
    end subroutine coef2_par

end module tridag_module

