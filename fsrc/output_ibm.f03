subroutine output_ibm (MN,MP,num_ib,num_solide,nproc,myid, &
    indici_celle_bloccate,indici_CELLE_IB)

    use mpi

    implicit none

    !------------------------------------------------------------------
    ! DICHIARAZIONI TIPI
    integer MN,MP,i
    integer num_ib,num_solide,nproc,myid,ierr

    integer num_ib_totale,num_solide_totale

    integer indici_CELLE_IB(MN,6)
    integer indici_celle_bloccate(MP,3)

    integer idisp(0:nproc-1)

    integer,allocatable :: send_solide1(:)
    integer,allocatable :: send_solide2(:)
    integer,allocatable :: send_solide3(:)
    integer,allocatable :: send_ib1(:),send_ib2(:),send_ib3(:)
    integer,allocatable :: send_ib4(:),send_ib5(:),send_ib6(:)

    integer,allocatable :: recv_solide1(:)
    integer,allocatable :: recv_solide2(:)
    integer,allocatable :: recv_solide3(:)
    integer,allocatable :: recv_ib1(:),recv_ib2(:),recv_ib3(:)
    integer,allocatable :: recv_ib4(:),recv_ib5(:),recv_ib6(:)

    integer,allocatable :: num_solide_array(:),num_ib_array(:)

    !------------------------------------------------------------------
    !     rendo noto il numero di celle solide e ib che ciascun proc ha
    allocate(num_solide_array(0:nproc-1))
    allocate(num_ib_array(0:nproc-1))

    allocate(send_solide1(1:num_solide))
    allocate(send_solide2(1:num_solide))
    allocate(send_solide3(1:num_solide))

    allocate(send_ib1(1:num_ib))
    allocate(send_ib2(1:num_ib))
    allocate(send_ib3(1:num_ib))
    allocate(send_ib4(1:num_ib))
    allocate(send_ib5(1:num_ib))
    allocate(send_ib6(1:num_ib))

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call MPI_ALLGATHER(num_solide,1,MPI_INTEGER, &
        num_solide_array,1,MPI_INTEGER, &
        MPI_COMM_WORLD,ierr)

    call MPI_ALLGATHER(num_ib,1,MPI_INTEGER, &
        num_ib_array,1,MPI_INTEGER, &
        MPI_COMM_WORLD,ierr)

    num_solide_totale=0
    num_ib_totale=0

    do i=0,nproc-1
        num_solide_totale=num_solide_totale+num_solide_array(i)
        num_ib_totale=num_ib_totale+num_ib_array(i)
    end do

    if(myid.eq.0)then
        write(*,*)'scrittura indici solide ib per postprocessing'
        write(*,*)'numero celle solide e ib totali: ', &
            num_solide_totale,' ',num_ib_totale
    end if

    allocate(recv_solide1(1:num_solide_totale))
    allocate(recv_solide2(1:num_solide_totale))
    allocate(recv_solide3(1:num_solide_totale))

    allocate(recv_ib1(1:num_ib_totale))
    allocate(recv_ib2(1:num_ib_totale))
    allocate(recv_ib3(1:num_ib_totale))
    allocate(recv_ib4(1:num_ib_totale))
    allocate(recv_ib5(1:num_ib_totale))
    allocate(recv_ib6(1:num_ib_totale))

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(myid.eq.0)then
        write(*,*)'controllo 1'
    end if
    !...................................................................
    !       passo indici celle solide tra i processori
    idisp(0)=0
    do i=1,nproc-1
        idisp(i)=idisp(i-1)+num_solide_array(i)
    end do

    do i=1,num_solide
        send_solide1(i)=indici_celle_bloccate(i,1)
        send_solide2(i)=indici_celle_bloccate(i,2)
        send_solide3(i)=indici_celle_bloccate(i,3)
    end do

    call MPI_GATHERV(send_solide1(1),num_solide,MPI_INTEGER, &
        recv_solide1,num_solide_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(send_solide2(1),num_solide,MPI_INTEGER, &
        recv_solide2,num_solide_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(send_solide3(1),num_solide,MPI_INTEGER, &
        recv_solide3,num_solide_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)

    if(myid.eq.0)then
        open(1000,file='output_indici_solide',status='unknown')
        write(1000,*)num_solide_totale
        do i=1,num_solide_totale
            write(1000,*)recv_solide1(i),recv_solide2(i),recv_solide3(i)
        end do
        close(1000)
    end if

    if(myid.eq.0)then
        write(*,*)'controllo 2'
    end if
    !.............................................................................
    !       passo indici celle ib tra i processori
    idisp=0
    idisp(0)=0
    do i=1,nproc-1
        idisp(i)=idisp(i-1)+num_ib_array(i)
    end do

    do i=1,num_ib
        send_ib1(i)=indici_CELLE_IB(i,1)
        send_ib2(i)=indici_CELLE_IB(i,2)
        send_ib3(i)=indici_CELLE_IB(i,3)
        send_ib4(i)=indici_CELLE_IB(i,4)
        send_ib5(i)=indici_CELLE_IB(i,5)
        send_ib6(i)=indici_CELLE_IB(i,6)
    end do

    if(myid.eq.0)then
        write(*,*)'controllo 3'
    end if


    call MPI_GATHERV(send_ib1(1),num_ib,MPI_INTEGER, &
        recv_ib1,num_ib_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(send_ib2(1),num_ib,MPI_INTEGER, &
        recv_ib2,num_ib_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(send_ib3(1),num_ib,MPI_INTEGER, &
        recv_ib3,num_ib_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(send_ib4(1),num_ib,MPI_INTEGER, &
        recv_ib4,num_ib_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(send_ib5(1),num_ib,MPI_INTEGER, &
        recv_ib5,num_ib_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(send_ib6(1),num_ib,MPI_INTEGER, &
        recv_ib6,num_ib_array,idisp,MPI_INTEGER, &
        0,MPI_COMM_WORLD,ierr)
    if(myid.eq.0)then
        open(1001,file='output_indici_ib',status='unknown')
        write(1001,*)num_ib_totale
        do i=1,num_ib_totale
            write(1001,*)recv_ib1(i),recv_ib2(i),recv_ib3(i), &
                recv_ib4(i),recv_ib5(i),recv_ib6(i)
        end do
        close(1001)
    end if

    if(myid.eq.0)then
        write(*,*)'controllo 4'
    end if

    deallocate(send_solide1,send_solide2,send_solide3)
    deallocate(recv_solide1,recv_solide2,recv_solide3)
    deallocate(send_ib1,send_ib2,send_ib3,send_ib4,send_ib5,send_ib6)
    deallocate(recv_ib1,recv_ib2,recv_ib3,recv_ib4,recv_ib5,recv_ib6)
    deallocate(num_solide_array,num_ib_array)

    if(myid.eq.0)then
        write(*,*)'controllo 5'
    end if


    return
end
