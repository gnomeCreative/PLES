!*************************************************************************
subroutine tridiag_trasp_para_rho(akapt,g33_tr,giac_tr,del, &
    akapt_piano,pran,isc)
    !*************************************************************************
    ! transpose operation for matrix.
    ! each plane XZ it is seen as two dimensional matrix and it is transposed,
    ! in this way the solution of the tridiagonal matrix is made in parallel in
    ! i instead of k, then it will be necesseary to transpose again the result
    ! (it is necessary to transpose the matrix used in "coed3", so g33, akapt,
    ! giac)

    use mysending
    !
    use scala3
    use period
    use tipologia
    !
    use mpi

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
end

