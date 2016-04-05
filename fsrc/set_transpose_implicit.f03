!***********************************************************************
subroutine set_transpose_implicit(g33_tr,giac_tr)
    !***********************************************************************
    use mysending
    use myarrays_metri3
    !
    use scala3
    !
    use mpi

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
    !     >    	   buffer_plane_loc(1),count_plane/nproc,MPI_REAL_SD,
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
            call MPI_SSEND(buffer_plane(1),count_plane,MPI_REAL_SD, &
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
end
