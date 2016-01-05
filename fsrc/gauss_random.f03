!************************************************************
subroutine gauss_random(ktime,n1,n2,n3,myid,nproc,kparasta,kparaend)
    !************************************************************
    !
    use myarrays_WB
    !
    use mpi

    implicit none
      

    real, allocatable :: eps(:,:)
    real, allocatable :: stepg(:)
    real, allocatable :: A(:)
    integer, allocatable :: countf(:)

    integer count, i, j, k, ktime
    integer kk, ii, myid, nproc,n1,n2,n3
    integer kparasta, kparaend
    integer nreal, nstep, npts, idly
    integer iverifico,ierr

    real smean,deviaz,sum, percentual, passo,sumnew
    real massimo, minimo, ent, cortim
    real std,mean,cgauss,dt
    real varianza      ! GORAN
    real time1,time2,time3,time4,time5,timeB,timeA
    !
    !      nreal = 1000        !number of realizations=1000
    !      nreal = n1*n3
    nreal = n1*(kparaend-kparasta+1) !n3
    nstep = 1          !max delay in corr. func=10
    dt    = .5          !time step size=.5
    cortim= 5.          !corr. time in the same units as DT=5
    iverifico = 1
 
    allocate(eps(nreal,-1:nstep*2))
    allocate(A(1:nreal*(nstep+1)*2))
    allocate(stepg(1:100))
    allocate(countf(1:100))
    !
    !--------------------------------------------------------------
    !**************************************************************
    !--------------------------------------------------------------
    ! initialize
    !      call cpu_time(time1)
    call cgaus0(dt,cortim)
    !      call cpu_time(time2)
    !      write(*,*)'cgaus0',time2-time1
    ! store several series of Gaussian noise values in array EPS.
    !.......................................................................
    !      call cpu_time(timeA)
      
    count   = 1
    sum     = 0.
    do i=1,nreal             !realizations
        do j=0,nstep*2          !time delays
            eps(i,j)= cgauss()
            sum     = sum + eps(i,j)
        end do
    end do
    sum = sum/float(nreal*(nstep*2+1))
    !      write(*,*)myid,'sum',sum
    !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
    sumnew = 0.
    do i=1,nreal             !realizations
        do j=0,nstep*2          !time delays
            eps(i,j)= eps(i,j) - sum
            sumnew = sumnew + eps(i,j)
        end do
    end do
    !      sumnew = sumnew/float(nreal*(nstep*2+1))
    !      write(*,*)myid,'sum new',sumnew
    !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !
    !.......................................................................
    !
    !      call cpu_time(timeA)
    deviaz  = 0.
    do i=1,nreal
        do j=0,nstep*2
            deviaz = deviaz + (eps(i,j)-sum)*(eps(i,j)-sum)
        end do
    end do
    deviaz = deviaz/float(nreal*(nstep*2+1))
    deviaz = sqrt(deviaz)
    !      write(*,*)myid,'deviazione',deviaz
    !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !      call cpu_time(timeB)
    !      write(*,*)'loop 2',timeB-timeA
    !
    !.......................................................................
    !
    !      call cpu_time(timeA)
    !      count = 0
    !      do i=1,nreal
    !      do j=0,nstep*2
    !          if (abs(eps(i,j)).lt.deviaz) then
    !           count = count +1
    !          end if
    !      end do
    !      end do
    !      percentual = float(count)/float(nreal*nstep*2)
    !      call cpu_time(timeB)
    !      write(*,*)'loop 3',timeB-timeA
    !
    !.......................................................................
    !
    !      call cpu_time(time3)
    !      write(*,*)'store',time3-time2
	
    !
    !--------------------------------------------------------------
    !**************************************************************
    !--------------------------------------------------------------
    ! qui 'normalizzo' i dati in modo che rispondano alle
    ! richieste di GORAN: G(0,varianza)
    ! la deviazione standard richiesta vale "varianza" e media 0 (Noh03)
    !    GORAN

    varianza = 0.2        ! GORAN : immettere la varianza desiderata
    varianza = varianza/deviaz
    do i=1,nreal
        do j=0,nstep*2
            eps(i,j)= eps(i,j)*varianza
        end do
    end do

    ! GORAN : verifichiamo varianza e media...

    if(iverifico==1)then
        sum = 0.
        do i=1,nreal
            do j=0,nstep*2
                sum = sum + eps(i,j)
            end do
        end do
        sum = sum/float(nreal*(nstep*2+1))
        !     	 write(*,*)myid,'sum',sum
        !     	 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        deviaz = 0.
        do i=1,nreal
            do j=0,nstep*2
                deviaz = deviaz + (eps(i,j)-sum)*(eps(i,j)-sum)
            end do
        end do
        deviaz = deviaz/float(nreal*(nstep*2+1))
        deviaz = sqrt(deviaz)
        !     	 write(*,*)myid,'deviaz',deviaz
        !     	 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     	 
        if(myid==0)then
            write(*,*)'Media vento= ',sum,'Deviazione= ',deviaz
        end if
    end if
    !
    !      call cpu_time(time4)
    !      write(*,*)'normalize',time4-time3
    !------------------------------------------------------
    !------------------------------------------------------
    !
    nreal = n1*(kparaend-kparasta+1)

    old_j = old_j + 1 + myid
    j = old_j
    if (j>=nstep) then
        old_j = 1
        j = old_j
    end if

    do i=1,nreal !(n1*n3)
        A(i) = eps(i,j)
    end do

    !
    ! permutazioni: assegno i valori ad Fx(i,k)
    !
    count = 0
    do k=kparasta,kparaend
        do i=1,n1
            count = count +1
            Fx(i,k) = A(count)
        end do
    end do
    !
    ! ora assegno i valori ad Fz(i,k)
    !
    j = old_j+1

    do i=1,nreal !(n1*n3)
        A(i) = eps(i,j)
    end do

    count = 0
    do k=kparasta,kparaend
        do i=1,n1
            count = count +1
            Fz(i,k) = A(count)
        end do
    end do
    !
    !
    !--------------------------------------------------------
    !********************************************************
    !--------------------------------------------------------
    !
    deallocate(eps)
    deallocate(countf)
    deallocate(stepg)
    deallocate(A)

!
end subroutine gauss_random
!
!-----------------------------------------------------------------
!*****************************************************************
!-----------------------------------------------------------------
! initialize the RNG's
! and set the color of gaussian noise
! DT is the time step used in whatever process the colored Gaussian noise
! is used.
! CORTIM is correlation time in the same units as time step DT.
! WHITE=.true. means generate white gaussian noise which happens when
! CORTIM=0. This flag is used in CGAUSS.
! Here we use the flat distribution RAN1 also taken from Numerical Recipe
! but any other good flat distribution random number generator will do.

subroutine cgaus0(dt,cortim)

    implicit none

    real cape,dt,l1me2,cgauss
    real cortim,x
    logical white
    common /color/l1me2,cape,white

    if(cortim==0.)then
        white=.true.
        l1me2=-2.000			   !white noise
        cape=0.0
    else
        white=.false.
        cape=exp(-dt/real(cortim))
        !parameter needed in CGAUSS
        l1me2=-(real(1.)-cape*cape)*real(2./cortim)
    endif
    x=cgauss()	    !initialize CGAUSS value

    return
end
!
!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------
!
real function cgauss()

    implicit none
    !
    integer iset
    integer n, iseed
    logical white
    real  fac,gset,rsq,v1,v2,l1me2,h,cape,rinv,prev
    common /color/l1me2,cape,white
      
    save iset,gset,prev
    data iset/0/
    data prev/0.0d0/
    !
    iseed = 12345.6

    rinv =  0.000014874754566549652 ! 2./134456.
    if (iset==0) then
      
        !      n  = abs(iseed)
        n  = mod(8121*n+28411,134456)
        !      v1 = 2.*real(n)/134456.
        v1 = real(n)*rinv

        call random_number(v1)
        v1=2.*v1-1

        call random_number(v2)
        v2=2.*v2-1
        !
        !        rsq=v1**2.+v2**2.
        rsq=v1*v1+v2*v2
     
        do while (rsq>=1. .or. rsq==0.) ! Alessandro there was a horrible goto here

            call random_number(v1)
            v1=2.*v1-1

            call random_number(v2)
            v2=2.*v2-1
            !
            !        rsq=v1**2.+v2**2.
            rsq=v1*v1+v2*v2

        end do

        fac=sqrt(l1me2*log(rsq)/rsq)
        gset=v1*fac
        h=v2*fac
        iset=1
    else
        h=gset
        iset=0
    endif
      
    if(white)then  !please note that the time step vs its sqrt
        cgauss=h      !in integration is previously set in PARAM
    else
        cgauss=prev*cape+h
        prev=cgauss
    end if
      
    return
end
