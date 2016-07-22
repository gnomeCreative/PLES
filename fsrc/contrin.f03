subroutine contrin

   ! controvariant fluxes compuatation from interpolation data (dns)
   ! generalized periodicity
   !
   use myarrays_velo3
   use myarrays_metri3
   use mysending
   !
   use scala3
   use period
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   integer i,j,k
   integer ierr
   real diver,divmax,divmax_loc
   integer kpsta,kpend
   !-----------------------------------------------------------------------
   !
   ! flux J-1*U
   !
   do k=kparasta,kparaend !1,jz
      do j=1,n2
         !        face 1
         uc(0,j,k)=csx(0,j,k)*u(0,j,k) &
            +csy(0,j,k)*v(0,j,k) &
            +csz(0,j,k)*w(0,j,k)
         !        face 2
         uc(n1,j,k)=csx(n1,j,k)*u(n1+1,j,k) &
            +csy(n1,j,k)*v(n1+1,j,k) &
            +csz(n1,j,k)*w(n1+1,j,k)
         do i=ip,n1-ip
            !           inside the domain and face 1 and 2 if periodic
            uc(i,j,k)= &
               csx(i,j,k)*.5*(u(i,j,k)+u(i+1,j,k)) &
               +csy(i,j,k)*.5*(v(i,j,k)+v(i+1,j,k)) &
               +csz(i,j,k)*.5*(w(i,j,k)+w(i+1,j,k))
         enddo
      enddo
   enddo
   !
   ! flusso J-1*V
   !
   do k=kparasta,kparaend !1,jz
      do i=1,n1
         !        face 3
         vc(i,0,k)=etx(i,0,k)*u(i,0,k) &
            +ety(i,0,k)*v(i,0,k) &
            +etz(i,0,k)*w(i,0,k)
         !        face 4
         vc(i,n2,k)=etx(i,n2,k)*u(i,n2+1,k) &
            +ety(i,n2,k)*v(i,n2+1,k) &
            +etz(i,n2,k)*w(i,n2+1,k)
         do j=1,n2-1
            !           inside the domain and face 3 and 4 if periodic
            vc(i,j,k)= &
               etx(i,j,k)*.5*(u(i,j,k)+u(i,j+1,k)) &
               +ety(i,j,k)*.5*(v(i,j,k)+v(i,j+1,k)) &
               +etz(i,j,k)*.5*(w(i,j,k)+w(i,j+1,k))
         enddo
      enddo
   enddo
   !
   ! flux J-1*W
   !
   kpsta = kparasta-1
   kpend = kparaend
   if(myid.eq.0)then
      kpsta = kpsta + kp
      do j=1,n2
         do i=1,n1
            !          face 5
            wc(i,j,0)=ztx(i,j,0)*u(i,j,0) &
               +zty(i,j,0)*v(i,j,0) &
               +ztz(i,j,0)*w(i,j,0)
            do k=kpsta,kpend  !kp,jz-kp
               !
               wc(i,j,k)= &
                  ztx(i,j,k)*.5*(u(i,j,k)+u(i,j,k+1)) &
                  +zty(i,j,k)*.5*(v(i,j,k)+v(i,j,k+1)) &
                  +ztz(i,j,k)*.5*(w(i,j,k)+w(i,j,k+1))
            !
            enddo
         enddo
      enddo
   elseif(myid.eq. nproc-1)then
      kpend = kparaend-kp
      do j=1,n2
         do i=1,n1
            !          face 6
            wc(i,j,n3)=ztx(i,j,n3)*u(i,j,n3+1) &
               +zty(i,j,n3)*v(i,j,n3+1) &
               +ztz(i,j,n3)*w(i,j,n3+1)
            do k=kpsta,kpend !kp,jz-kp
               !             inside the domain and at face 6 if periodic
               wc(i,j,k)= &
                  ztx(i,j,k)*.5*(u(i,j,k)+u(i,j,k+1)) &
                  +zty(i,j,k)*.5*(v(i,j,k)+v(i,j,k+1)) &
                  +ztz(i,j,k)*.5*(w(i,j,k)+w(i,j,k+1))
            enddo
         enddo
      enddo
   else
      do j=1,n2
         do i=1,n1
            do k=kpsta,kpend !kp,jz-kp
               !          inside the domain
               wc(i,j,k)= &
                  ztx(i,j,k)*.5*(u(i,j,k)+u(i,j,k+1)) &
                  +zty(i,j,k)*.5*(v(i,j,k)+v(i,j,k+1)) &
                  +ztz(i,j,k)*.5*(w(i,j,k)+w(i,j,k+1))
            enddo
         enddo
      enddo
   end if
   !
   ! compute initial divergence
   !
   divmax=0.
   !
   !     AAA needs isend for wc(kparasta-1)
   do k=kparasta,kparaend !1,jz
      do j=1,n2
         do i=1,n1
            diver=uc(i,j,k)-uc(i-1,j,k) &
               +vc(i,j,k)-vc(i,j-1,k) &
               +wc(i,j,k)-wc(i,j,k-1)
            !
            diver=abs(diver)
            divmax=max(divmax,diver)
         !
         enddo
      enddo
   enddo
   !
   call MPI_REDUCE(divmax_loc,divmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     
   write(*,*)'divergence dns-interpolato:',divmax
   !
   return
end
