!***********************************************************************
subroutine contrin
   !***********************************************************************
   ! controvariant fluxes compuatation from interpolation data (dns)
   ! generalized periodicity
   !
   use myarrays_velo3
   use myarrays_metri3
   use mysending
   !
   use scala3
   use period
   use tipologia
   !
   use mpi

   implicit none
   !


   !
   !-----------------------------------------------------------------------
   ! array declaration

   integer i,j,k
   integer ierr
   real diver,divmax,divmax_loc
   integer kpsta,kpend
   !
   !-----------------------------------------------------------------------
   !
   ! flux J-1*U
   !
   do k=kparasta,kparaend !1,jz
      do j=1,jy
         !        face 1
         uc(0,j,k)=csx(0,j,k)*u(0,j,k) &
            +csy(0,j,k)*v(0,j,k) &
            +csz(0,j,k)*w(0,j,k)
         !        face 2
         uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k) &
            +csy(jx,j,k)*v(jx+1,j,k) &
            +csz(jx,j,k)*w(jx+1,j,k)
         do i=ip,jx-ip
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
      do i=1,jx
         !        face 3
         vc(i,0,k)=etx(i,0,k)*u(i,0,k) &
            +ety(i,0,k)*v(i,0,k) &
            +etz(i,0,k)*w(i,0,k)
         !        face 4
         vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k) &
            +ety(i,jy,k)*v(i,jy+1,k) &
            +etz(i,jy,k)*w(i,jy+1,k)
         do j=jp,jy-jp
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
      do j=1,jy
         do i=1,jx
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
      do j=1,jy
         do i=1,jx
            !          face 6
            wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1) &
               +zty(i,j,jz)*v(i,j,jz+1) &
               +ztz(i,j,jz)*w(i,j,jz+1)
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
      do j=1,jy
         do i=1,jx
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
      do j=1,jy
         do i=1,jx
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
   call MPI_REDUCE(divmax_loc,divmax,1,MPI_REAL_SD,MPI_MAX,0, &
      MPI_COMM_WORLD,ierr)
     
   write(*,*)'divergence dns-interpolato:',divmax
   !
   return
end
