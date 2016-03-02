!***********************************************************************
subroutine sett(myid,nproc,n,ksta,kend,kfsta,ix,iy,iz,nx,ny,nz, &
    i_dx,i_sn,i_sp,i_st,i_av,i_in, &
    r11,r12,r13,r21,r22,r23,r31,r32,r33,aa,bb,cc)
    !***********************************************************************
    !     set index
    !
    use period
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     variables declaration
    integer ierr
    integer ksta,kend,kfsta
    integer i,j,k,ij
    integer ix,iy,iz
    integer nx,ny,nz
      
    integer i_dx(nx,ny,ksta:kend)
    integer i_sn(nx,ny,ksta:kend)
    integer i_sp(nx,ny,ksta:kend)
    integer i_st(nx,ny,ksta:kend)
    integer i_av(nx,ny,ksta:kend)
    integer i_in(nx,ny,ksta:kend)

    real  aa(nx,0:ny+1,ksta:kend)
    real  bb(nx,0:ny+1,ksta:kend)
    real  cc(nx,0:ny+1,ksta:kend)
      
    real r11(0:nx,ny,ksta:kend)
    real r12(0:nx,ny,ksta:kend)
    real r13(0:nx,ny,ksta:kend)
      
    real r21(nx,0:ny,ksta:kend)
    real r22(nx,0:ny,ksta:kend)
    real r23(nx,0:ny,ksta:kend)
      
    real r31(nx,ny,kfsta:kend)
    real r32(nx,ny,kfsta:kend)
    real r33(nx,ny,kfsta:kend)
      
    integer ipot
    real pot,ppot1,ppot2
    integer n,myid,nproc
      
    !-----------------------------------------------------------------------
    !
    ! in all the field
    !
    do k=ksta,kend
        do j=1,ny
            do i=1,nx
                !
                i_dx(i,j,k) =1
                i_sn(i,j,k) =1
                i_sp(i,j,k) =1
                i_st(i,j,k) =1
                i_av(i,j,k) =1
                i_in(i,j,k) =1
            !
            end do
        end do
    end do
    !
    ! index sides sn and dx (1 and 2), general periodicity
    !
    do j=1,ny
        do k=ksta,kend
            !
            i_sn(1 ,j,k) =1-ip
            i_dx(ix,j,k) =1-ip
                !
        end do
    end do
    !
    ! index sides bottom and upper (3 and 4), general periodicity
    !
    do i=1,nx
        do k=ksta,kend
            !
            i_st(i,1 ,k) =1-jp
            i_sp(i,iy,k) =1-jp
                    !
        end do
    end do
    !
    ! index sides back and front (6 and 5), general periodicity

    !hicco MA SIAMO SICUI CHE SIANO COSI' INDIETRO E AVANTI?????
    !hicco DOVREBBE ESSERE IL CONTRARIO !!!!!!!

    ! indici parete  indietro e avanti: (periodicita generalizzata)
    !
    if(myid.eq.0)then
        do j=1,ny
            do i=1,nx
                i_in(i,j,ksta) =1-kp
            end do
        end do
    elseif(myid.eq.nproc-1)then
        write(*,*)'IZ,NZ',iz,nz,ksta,kend
        do j=1,ny
            do i=1,nx
                i_av(i,j,kend) =1-kp
            end do
        end do
    end if

    !-----------------------------------------------------------------------
    ! for multigrid in case of line sor in j direction (implict solution)
    ! build fixed coefficent for the tridiagonal matrix
    !
    ipot=2**(n-1)
    pot=float(ipot)
    ppot1=1/pot
    ppot2=1/pot/pot
            
    do k=ksta,kend !1,nz
        do i=1,nx

            do j=1,ny
                aa(i,j,k)=ppot2*r22(i,j-1,k)
                bb(i,j,k)=-ppot2*(r22(i  ,j,k  )+r22(i,j-1,k) &
                    +r11(i  ,j,k  )*i_dx(i,j,k) &
                    +r11(i-1,j,k  )*i_sn(i,j,k) &
                    +r33(i  ,j,k  )*i_av(i,j,k) &
                    +r33(i  ,j,k-1)*i_in(i,j,k))
                cc(i,j,k)=ppot2*r22(i,j,k)
            end do

            !     side bottom (3) aa=0
            aa(i,0,k)= 0.
            bb(i,0,k)= r22(i,0,k)*ppot2
            cc(i,0,k)=-r22(i,0,k)*ppot2
     
            !     side upper (4) cc=0
            aa(i,ny+1,k)=-r22(i,ny,k)*ppot2
            bb(i,ny+1,k)= r22(i,ny,k)*ppot2
            cc(i,ny+1,k)= 0.

        end do
    end do

    return
end
