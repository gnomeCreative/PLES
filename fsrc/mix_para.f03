!***********************************************************************
subroutine mix_para(inmod,iq,ktime,i_print,lagr)
   !***********************************************************************
   ! compute scale similar part for momentum eq.
   !
   use filter_module
   use turbo_module
   use myarrays_velo3
   use myarrays_metri3
   use mysending
   !
   use scala3
   use subgrid
   use period
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr
   integer kparastam,kparaendm
   integer status(MPI_STATUS_SIZE)
   !
   real sbuff((n1+2)*(n2+2)*40)
   real rbuff((n1+2)*(n2+2)*40)
   !
   integer i,j,k,kparastal,iq,inmod,m,lll
   integer kparastan
   integer ii,jj,kk
   integer,parameter :: debugg = 0
   integer ktime,i_print,lagr
   real somma
   !
   real,allocatable :: sus_loc11(:),sus_loc12(:),sus_loc13(:)
   real,allocatable :: sus_loc21(:),sus_loc22(:),sus_loc23(:)
   real,allocatable :: sus_loc31(:),sus_loc32(:),sus_loc33(:)
   !-----------------------------------------------------------------------


   do k=kparasta,kparaend
      do j=1,jy
         do i=0,jx
            cgra1(i,j,k)=0.
         enddo
      enddo
   enddo
   !
   do k=kparasta,kparaend
      do j=0,jy
         do i=1,jx
            cgra2(i,j,k)=0.
         enddo
      enddo
   enddo
   !
   if(myid.eq.0)then
      kparastal=0
   else
      kparastal=kparasta
   endif

   do k=kparastal,kparaend
      do j=1,jy
         do i=1,jx
            cgra3(i,j,k)=0.
         enddo
      enddo
   enddo
   !
   !-----------------------------------------------------------------------
   ! compute controvariant term for velocity and for mixed product only
   ! for scale similar
   !
   if(inmod.eq.1)then

      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx

               apcsx(i,j,k)=0.5*(csx(i,j,k)+csx(i-1,j,k))
               apcsy(i,j,k)=0.5*(csy(i,j,k)+csy(i-1,j,k))
               apcsz(i,j,k)=0.5*(csz(i,j,k)+csz(i-1,j,k))
               !
               apetx(i,j,k)=0.5*(etx(i,j,k)+etx(i,j-1,k))
               apety(i,j,k)=0.5*(ety(i,j,k)+ety(i,j-1,k))
               apetz(i,j,k)=0.5*(etz(i,j,k)+etz(i,j-1,k))
               !
               apztx(i,j,k)=0.5*(ztx(i,j,k)+ztx(i,j,k-1))
               apzty(i,j,k)=0.5*(zty(i,j,k)+zty(i,j,k-1))
               apztz(i,j,k)=0.5*(ztz(i,j,k)+ztz(i,j,k-1))
               !
               ! compute controvariant component at the centroid
               !
               uco(i,j,k)= apcsx(i,j,k)*u(i,j,k)+ &
                  apcsy(i,j,k)*v(i,j,k)+ &
                  apcsz(i,j,k)*w(i,j,k)
               !
               vco(i,j,k)= apetx(i,j,k)*u(i,j,k)+ &
                  apety(i,j,k)*v(i,j,k)+ &
                  apetz(i,j,k)*w(i,j,k)
               !
               wco(i,j,k)= apztx(i,j,k)*u(i,j,k)+ &
                  apzty(i,j,k)*v(i,j,k)+ &
                  apztz(i,j,k)*w(i,j,k)
               !
               ! compute mixed product
               !      uuco(i,j,k)=uco(i,j,k)*u(i,j,k)
               vuco(i,j,k)=uco(i,j,k)*v(i,j,k)
               wuco(i,j,k)=uco(i,j,k)*w(i,j,k)
               !
               uvco(i,j,k)=vco(i,j,k)*u(i,j,k)
               vvco(i,j,k)=vco(i,j,k)*v(i,j,k)
               wvco(i,j,k)=vco(i,j,k)*w(i,j,k)
               !
               uwco(i,j,k)=wco(i,j,k)*u(i,j,k)
               vwco(i,j,k)=wco(i,j,k)*v(i,j,k)
               wwco(i,j,k)=wco(i,j,k)*w(i,j,k)

            enddo
         enddo
      enddo
      !
      ! generalized periodicity
      !
      ! extrapolation for velocity on sides i=0 and i=jx+1
      !
      do k=kparasta,kparaend
         do j=1,jy

            uco(0,j,k)=(1-ip)*uco(jx,j,k) + ip*(2.*uco(1,j,k)-uco(2,j,k))
            vco(0,j,k)=(1-ip)*vco(jx,j,k) + ip*(2.*vco(1,j,k)-vco(2,j,k))
            wco(0,j,k)=(1-ip)*wco(jx,j,k) + ip*(2.*wco(1,j,k)-wco(2,j,k))
            !
            uuco(0,j,k)=(1-ip)*uuco(jx,j,k) + ip*(2.*uuco(1,j,k)-uuco(2,j,k))
            uvco(0,j,k)=(1-ip)*uvco(jx,j,k) + ip*(2.*uvco(1,j,k)-uvco(2,j,k))
            uwco(0,j,k)=(1-ip)*uwco(jx,j,k) + ip*(2.*uwco(1,j,k)-uwco(2,j,k))
            vuco(0,j,k)=(1-ip)*vuco(jx,j,k) + ip*(2.*vuco(1,j,k)-vuco(2,j,k))
            vvco(0,j,k)=(1-ip)*vvco(jx,j,k) + ip*(2.*vvco(1,j,k)-vvco(2,j,k))
            vwco(0,j,k)=(1-ip)*vwco(jx,j,k) + ip*(2.*vwco(1,j,k)-vwco(2,j,k))
            wuco(0,j,k)=(1-ip)*wuco(jx,j,k) + ip*(2.*wuco(1,j,k)-wuco(2,j,k))
            wvco(0,j,k)=(1-ip)*wvco(jx,j,k) + ip*(2.*wvco(1,j,k)-wvco(2,j,k))
            wwco(0,j,k)=(1-ip)*wwco(jx,j,k) + ip*(2.*wwco(1,j,k)-wwco(2,j,k))
            !
            apcsx(0,j,k)=(1-ip)*apcsx(jx,j,k) + ip*(2.*apcsx(1,j,k)-apcsx(2,j,k))
            apcsy(0,j,k)=(1-ip)*apcsy(jx,j,k) + ip*(2.*apcsy(1,j,k)-apcsy(2,j,k))
            apcsz(0,j,k)=(1-ip)*apcsz(jx,j,k) + ip*(2.*apcsz(1,j,k)-apcsz(2,j,k))
            apetx(0,j,k)=(1-ip)*apetx(jx,j,k) + ip*(2.*apetx(1,j,k)-apetx(2,j,k))
            apety(0,j,k)=(1-ip)*apety(jx,j,k) + ip*(2.*apety(1,j,k)-apety(2,j,k))
            apetz(0,j,k)=(1-ip)*apetz(jx,j,k) + ip*(2.*apetz(1,j,k)-apetz(2,j,k))
            apztx(0,j,k)=(1-ip)*apztx(jx,j,k) + ip*(2.*apztx(1,j,k)-apztx(2,j,k))
            apzty(0,j,k)=(1-ip)*apzty(jx,j,k) + ip*(2.*apzty(1,j,k)-apzty(2,j,k))
            apztz(0,j,k)=(1-ip)*apztz(jx,j,k) + ip*(2.*apztz(1,j,k)-apztz(2,j,k))
            !
            !
            apcsx(jx+1,j,k)=(1-ip)*apcsx(1,j,k) + ip*(2.*apcsx(jx,j,k)-apcsx(jx-1,j,k))
            apcsy(jx+1,j,k)=(1-ip)*apcsy(1,j,k) + ip*(2.*apcsy(jx,j,k)-apcsy(jx-1,j,k))
            apcsz(jx+1,j,k)=(1-ip)*apcsz(1,j,k) + ip*(2.*apcsz(jx,j,k)-apcsz(jx-1,j,k))
            apetx(jx+1,j,k)=(1-ip)*apetx(1,j,k) + ip*(2.*apetx(jx,j,k)-apetx(jx-1,j,k))
            apety(jx+1,j,k)=(1-ip)*apety(1,j,k) + ip*(2.*apety(jx,j,k)-apety(jx-1,j,k))
            apetz(jx+1,j,k)=(1-ip)*apetz(1,j,k) + ip*(2.*apetz(jx,j,k)-apetz(jx-1,j,k))
            apztx(jx+1,j,k)=(1-ip)*apztx(1,j,k) + ip*(2.*apztx(jx,j,k)-apztx(jx-1,j,k))
            apzty(jx+1,j,k)=(1-ip)*apzty(1,j,k) + ip*(2.*apzty(jx,j,k)-apzty(jx-1,j,k))
            apztz(jx+1,j,k)=(1-ip)*apztz(1,j,k) + ip*(2.*apztz(jx,j,k)-apztz(jx-1,j,k))
            !
            uco(jx+1,j,k)=(1-ip)*uco(1,j,k) + ip*(2.*uco(jx,j,k)-uco(jx-1,j,k))
            vco(jx+1,j,k)=(1-ip)*vco(1,j,k) + ip*(2.*vco(jx,j,k)-vco(jx-1,j,k))
            wco(jx+1,j,k)=(1-ip)*wco(1,j,k) + ip*(2.*wco(jx,j,k)-wco(jx-1,j,k))
            !
            uuco(jx+1,j,k)=(1-ip)*uuco(1,j,k) + ip*(2.*uuco(jx,j,k)-uuco(jx-1,j,k))
            uvco(jx+1,j,k)=(1-ip)*uvco(1,j,k) + ip*(2.*uvco(jx,j,k)-uvco(jx-1,j,k))
            uwco(jx+1,j,k)=(1-ip)*uwco(1,j,k) + ip*(2.*uwco(jx,j,k)-uwco(jx-1,j,k))
            vuco(jx+1,j,k)=(1-ip)*vuco(1,j,k) + ip*(2.*vuco(jx,j,k)-vuco(jx-1,j,k))
            vvco(jx+1,j,k)=(1-ip)*vvco(1,j,k) + ip*(2.*vvco(jx,j,k)-vvco(jx-1,j,k))
            vwco(jx+1,j,k)=(1-ip)*vwco(1,j,k) + ip*(2.*vwco(jx,j,k)-vwco(jx-1,j,k))
            wuco(jx+1,j,k)=(1-ip)*wuco(1,j,k) + ip*(2.*wuco(jx,j,k)-wuco(jx-1,j,k))
            wvco(jx+1,j,k)=(1-ip)*wvco(1,j,k) + ip*(2.*wvco(jx,j,k)-wvco(jx-1,j,k))
            wwco(jx+1,j,k)=(1-ip)*wwco(1,j,k) + ip*(2.*wwco(jx,j,k)-wwco(jx-1,j,k))

         end do
      end do
      !
      ! extrapolation for velocity on sides j=0 and j=jy+1
      !
      do k=kparasta,kparaend
         do i=1,jx

            uco(i,0,k)=(1-jp)*uco(i,jy,k) + jp*(2.*uco(i,1,k)-uco(i,2,k))
            vco(i,0,k)=(1-jp)*vco(i,jy,k) + jp*(2.*vco(i,1,k)-vco(i,2,k))
            wco(i,0,k)=(1-jp)*wco(i,jy,k) + jp*(2.*wco(i,1,k)-wco(i,2,k))
            !
            uuco(i,0,k)=(1-jp)*uuco(i,jy,k) + jp*(2.*uuco(i,1,k)-uuco(i,2,k))
            uvco(i,0,k)=(1-jp)*uvco(i,jy,k) + jp*(2.*uvco(i,1,k)-uvco(i,2,k))
            uwco(i,0,k)=(1-jp)*uwco(i,jy,k) + jp*(2.*uwco(i,1,k)-uwco(i,2,k))
            vuco(i,0,k)=(1-jp)*vuco(i,jy,k) + jp*(2.*vuco(i,1,k)-vuco(i,2,k))
            vvco(i,0,k)=(1-jp)*vvco(i,jy,k) + jp*(2.*vvco(i,1,k)-vvco(i,2,k))
            vwco(i,0,k)=(1-jp)*vwco(i,jy,k) + jp*(2.*vwco(i,1,k)-vwco(i,2,k))
            wuco(i,0,k)=(1-jp)*wuco(i,jy,k) + jp*(2.*wuco(i,1,k)-wuco(i,2,k))
            wvco(i,0,k)=(1-jp)*wvco(i,jy,k) + jp*(2.*wvco(i,1,k)-wvco(i,2,k))
            wwco(i,0,k)=(1-jp)*wwco(i,jy,k) + jp*(2.*wwco(i,1,k)-wwco(i,2,k))
            !
            apcsx(i,0,k)=(1-jp)*apcsx(i,jy,k) + jp*(2.*apcsx(i,1,k)-apcsx(i,2,k))
            apcsy(i,0,k)=(1-jp)*apcsy(i,jy,k) + jp*(2.*apcsy(i,1,k)-apcsy(i,2,k))
            apcsz(i,0,k)=(1-jp)*apcsz(i,jy,k) + jp*(2.*apcsz(i,1,k)-apcsz(i,2,k))
            apetx(i,0,k)=(1-jp)*apetx(i,jy,k) + jp*(2.*apetx(i,1,k)-apetx(i,2,k))
            apety(i,0,k)=(1-jp)*apety(i,jy,k) + jp*(2.*apety(i,1,k)-apety(i,2,k))
            apetz(i,0,k)=(1-jp)*apetz(i,jy,k) + jp*(2.*apetz(i,1,k)-apetz(i,2,k))
            apztx(i,0,k)=(1-jp)*apztx(i,jy,k) + jp*(2.*apztx(i,1,k)-apztx(i,2,k))
            apzty(i,0,k)=(1-jp)*apzty(i,jy,k) + jp*(2.*apzty(i,1,k)-apzty(i,2,k))
            apztz(i,0,k)=(1-jp)*apztz(i,jy,k) + jp*(2.*apztz(i,1,k)-apztz(i,2,k))
            !
            !
            apcsx(i,jy+1,k)=(1-jp)*apcsx(i,1,k) + jp*(2.*apcsx(i,jy,k)-apcsx(i,jy-1,k))
            apcsy(i,jy+1,k)=(1-jp)*apcsy(i,1,k) + jp*(2.*apcsy(i,jy,k)-apcsy(i,jy-1,k))
            apcsz(i,jy+1,k)=(1-jp)*apcsz(i,1,k) + jp*(2.*apcsz(i,jy,k)-apcsz(i,jy-1,k))
            apetx(i,jy+1,k)=(1-jp)*apetx(i,1,k) + jp*(2.*apetz(i,jy,k)-apetx(i,jy-1,k))
            apety(i,jy+1,k)=(1-jp)*apety(i,1,k) + jp*(2.*apety(i,jy,k)-apety(i,jy-1,k))
            apetz(i,jy+1,k)=(1-jp)*apetz(i,1,k) + jp*(2.*apetz(i,jy,k)-apetz(i,jy-1,k))
            apztx(i,jy+1,k)=(1-jp)*apztx(i,1,k) + jp*(2.*apztx(i,jy,k)-apztx(i,jy-1,k))
            apzty(i,jy+1,k)=(1-jp)*apzty(i,1,k) + jp*(2.*apzty(i,jy,k)-apzty(i,jy-1,k))
            apztz(i,jy+1,k)=(1-jp)*apztz(i,1,k) + jp*(2.*apztz(i,jy,k)-apztz(i,jy-1,k))
            !
            uco(i,jy+1,k)=(1-jp)*uco(i,1,k) + jp*(2.*uco(i,jy,k)-uco(i,jy-1,k))
            vco(i,jy+1,k)=(1-jp)*vco(i,1,k) + jp*(2.*vco(i,jy,k)-vco(i,jy-1,k))
            wco(i,jy+1,k)=(1-jp)*wco(i,1,k) + jp*(2.*wco(i,jy,k)-wco(i,jy-1,k))
            !
            uuco(i,jy+1,k)=(1-jp)*uuco(i,1,k) + jp*(2.*uuco(i,jy,k)-uuco(i,jy-1,k))
            uvco(i,jy+1,k)=(1-jp)*uvco(i,1,k) + jp*(2.*uvco(i,jy,k)-uvco(i,jy-1,k))
            uwco(i,jy+1,k)=(1-jp)*uwco(i,1,k) + jp*(2.*uwco(i,jy,k)-uwco(i,jy-1,k))
            vuco(i,jy+1,k)=(1-jp)*vuco(i,1,k) + jp*(2.*vuco(i,jy,k)-vuco(i,jy-1,k))
            vvco(i,jy+1,k)=(1-jp)*vvco(i,1,k) + jp*(2.*vvco(i,jy,k)-vvco(i,jy-1,k))
            vwco(i,jy+1,k)=(1-jp)*vwco(i,1,k) + jp*(2.*vwco(i,jy,k)-vwco(i,jy-1,k))
            wuco(i,jy+1,k)=(1-jp)*wuco(i,1,k) + jp*(2.*wuco(i,jy,k)-wuco(i,jy-1,k))
            wvco(i,jy+1,k)=(1-jp)*wvco(i,1,k) + jp*(2.*wvco(i,jy,k)-wvco(i,jy-1,k))
            wwco(i,jy+1,k)=(1-jp)*wwco(i,1,k) + jp*(2.*wwco(i,jy,k)-wwco(i,jy-1,k))

         end do
      end do
      !
      ! periodicity for sides front and back
      !
      ! extrapolation for velocity on sides k=0 k=jz+1
      !
      ! I make a sending buffer to contain all the jx*jy planes for
      ! using only a sending procedure for the transfer between P0 and Pn-1

      do m=1,40*(jx+2)*(jy+2)
         sbuff(m)=0.
         rbuff(m)=0.
      enddo

      if (myid.eq.nproc-1) then

         call buffer1g(uco,1,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(vco,2,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(wco,3,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(uuco,4,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(uvco,5,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(uwco,6,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(vuco,7,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(vvco,8,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(vwco,9,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(wuco,10,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(wvco,11,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(wwco,12,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsx,13,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsy,14,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsz,15,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apetx,16,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apety,17,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apetz,18,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apztx,19,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apzty,20,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apztz,21,jz,myid,nproc,kparasta,kparaend)

      else if (myid.eq.0) then

         call buffer1g(uco,1,1,myid,nproc,kparasta,kparaend)
         call buffer1g(vco,2,1,myid,nproc,kparasta,kparaend)
         call buffer1g(wco,3,1,myid,nproc,kparasta,kparaend)
         call buffer1g(uuco,4,1,myid,nproc,kparasta,kparaend)
         call buffer1g(uvco,5,1,myid,nproc,kparasta,kparaend)
         call buffer1g(uwco,6,1,myid,nproc,kparasta,kparaend)
         call buffer1g(vuco,7,1,myid,nproc,kparasta,kparaend)
         call buffer1g(vvco,8,1,myid,nproc,kparasta,kparaend)
         call buffer1g(vwco,9,1,myid,nproc,kparasta,kparaend)
         call buffer1g(wuco,10,1,myid,nproc,kparasta,kparaend)
         call buffer1g(wvco,11,1,myid,nproc,kparasta,kparaend)
         call buffer1g(wwco,12,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsx,13,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsy,14,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsz,15,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apetx,16,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apety,17,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apetz,18,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apztx,19,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apzty,20,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apztz,21,1,myid,nproc,kparasta,kparaend)

      endif

      ! now the exchange so that P0 knows k=jz and Pn-1 knows k=1

      if (myid.eq.nproc-1) then

         call MPI_SENDRECV(sbuff1(1),21*(jx+2)*(jy+2), &
            MPI_REAL_SD,0,9101, &
            rbuff1(1),21*(jx+2)*(jy+2), &
            MPI_REAL_SD,0,8101, &
            MPI_COMM_WORLD,status,ierr)

      else if (myid.eq.0) then

         call MPI_SENDRECV(sbuff1(1),21*(jx+2)*(jy+2),MPI_REAL_SD, &
            nproc-1,8101, &
            rbuff1(1),21*(jx+2)*(jy+2),MPI_REAL_SD, &
            nproc-1,9101, &
            MPI_COMM_WORLD,status,ierr)
      
      endif

      if (myid.eq.0) then

         call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano4,4,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano5,5,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano6,6,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano7,7,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano8,8,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano9,9,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano10,10,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano11,11,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano12,12,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano13,13,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano14,14,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano15,15,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano16,16,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano17,17,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano18,18,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano19,19,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano20,20,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano21,21,myid,nproc,kparasta,kparaend)


      else if (myid.eq.nproc-1) then


         call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano4,4,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano5,5,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano6,6,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano7,7,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano8,8,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano9,9,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano10,10,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano11,11,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano12,12,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano13,13,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano14,14,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano15,15,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano16,16,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano17,17,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano18,18,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano19,19,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano20,20,myid,nproc,kparasta,kparaend)
         call buffer2gg(piano21,21,myid,nproc,kparasta,kparaend)



      endif

      ! now P0 knows the plane k=jz
      ! and Pn-1 knows k=1

      if(myid.eq.0)then

         do j=1,jy
            do i=1,jx
      

               uco(i,j,0)=(1-kp)*piano1(i,j) + kp*(2.*uco(i,j,1)-uco(i,j,2))
               vco(i,j,0)=(1-kp)*piano2(i,j) + kp*(2.*vco(i,j,1)-vco(i,j,2))
               wco(i,j,0)=(1-kp)*piano3(i,j) + kp*(2.*wco(i,j,1)-wco(i,j,2))
               !
               !
               uuco(i,j,0)=(1-kp)*piano4(i,j) + kp*(2.*uuco(i,j,1)-uuco(i,j,2))
               uvco(i,j,0)=(1-kp)*piano5(i,j) + kp*(2.*uvco(i,j,1)-uvco(i,j,2))
               uwco(i,j,0)=(1-kp)*piano6(i,j) + kp*(2.*uwco(i,j,1)-uwco(i,j,2))
               vuco(i,j,0)=(1-kp)*piano7(i,j) + kp*(2.*vuco(i,j,1)-vuco(i,j,2))
               vvco(i,j,0)=(1-kp)*piano8(i,j) + kp*(2.*vvco(i,j,1)-vvco(i,j,2))
               vwco(i,j,0)=(1-kp)*piano9(i,j) + kp*(2.*vwco(i,j,1)-vwco(i,j,2))
               wuco(i,j,0)=(1-kp)*piano10(i,j) + kp*(2.*wuco(i,j,1)-wuco(i,j,2))
               wvco(i,j,0)=(1-kp)*piano11(i,j) + kp*(2.*wvco(i,j,1)-wvco(i,j,2))
               wwco(i,j,0)=(1-kp)*piano12(i,j) + kp*(2.*wwco(i,j,1)-wwco(i,j,2))
               !
               apcsx(i,j,0)=(1-kp)*piano13(i,j) + kp*(2.*apcsx(i,j,1)-apcsx(i,j,2))
               apcsy(i,j,0)=(1-kp)*piano14(i,j) + kp*(2.*apcsy(i,j,1)-apcsy(i,j,2))
               apcsz(i,j,0)=(1-kp)*piano15(i,j) + kp*(2.*apcsz(i,j,1)-apcsz(i,j,2))
               apetx(i,j,0)=(1-kp)*piano16(i,j) + kp*(2.*apetx(i,j,1)-apetx(i,j,2))
               apety(i,j,0)=(1-kp)*piano17(i,j) + kp*(2.*apety(i,j,1)-apety(i,j,2))
               apetz(i,j,0)=(1-kp)*piano18(i,j) + kp*(2.*apetz(i,j,1)-apetz(i,j,2))
               apztx(i,j,0)=(1-kp)*piano19(i,j) + kp*(2.*apztx(i,j,1)-apztx(i,j,2))
               apzty(i,j,0)=(1-kp)*piano20(i,j) + kp*(2.*apzty(i,j,1)-apzty(i,j,2))
               apztz(i,j,0)=(1-kp)*piano21(i,j) + kp*(2.*apztz(i,j,1)-apztz(i,j,2))
            !
            end do
         end do

      endif
      !
      !
      if(myid.eq.nproc-1)then

         do j=1,jy
            do i=1,jx

               uco(i,j,jz+1)=(1-kp)*piano1(i,j) + kp*(2.*uco(i,j,jz)-uco(i,j,jz-1))
               vco(i,j,jz+1)=(1-kp)*piano2(i,j) + kp*(2.*vco(i,j,jz)-vco(i,j,jz-1))
               wco(i,j,jz+1)=(1-kp)*piano3(i,j) + kp*(2.*wco(i,j,jz)-wco(i,j,jz-1))

               uuco(i,j,jz+1)=(1-kp)*piano4(i,j) + kp*(2.*uuco(i,j,jz)-uuco(i,j,jz-1))
               uvco(i,j,jz+1)=(1-kp)*piano5(i,j) + kp*(2.*uvco(i,j,jz)-uvco(i,j,jz-1))
               uwco(i,j,jz+1)=(1-kp)*piano6(i,j) + kp*(2.*uwco(i,j,jz)-uwco(i,j,jz-1))
               vuco(i,j,jz+1)=(1-kp)*piano7(i,j) + kp*(2.*vuco(i,j,jz)-vuco(i,j,jz-1))
               vvco(i,j,jz+1)=(1-kp)*piano8(i,j) + kp*(2.*vvco(i,j,jz)-vvco(i,j,jz-1))
               vwco(i,j,jz+1)=(1-kp)*piano9(i,j) + kp*(2.*vwco(i,j,jz)-vwco(i,j,jz-1))
               wuco(i,j,jz+1)=(1-kp)*piano10(i,j) + kp*(2.*wuco(i,j,jz)-wuco(i,j,jz-1))
               wvco(i,j,jz+1)=(1-kp)*piano11(i,j) + kp*(2.*wvco(i,j,jz)-wvco(i,j,jz-1))
               wwco(i,j,jz+1)=(1-kp)*piano12(i,j) + kp*(2.*wwco(i,j,jz)-wwco(i,j,jz-1))



               apcsx(i,j,jz+1)=(1-kp)*piano13(i,j) + kp*(2.*apcsx(i,j,jz)-apcsx(i,j,jz-1))
               apcsy(i,j,jz+1)=(1-kp)*piano14(i,j) + kp*(2.*apcsy(i,j,jz)-apcsy(i,j,jz-1))
               apcsz(i,j,jz+1)=(1-kp)*piano15(i,j) + kp*(2.*apcsz(i,j,jz)-apcsz(i,j,jz-1))
               apetx(i,j,jz+1)=(1-kp)*piano16(i,j) + kp*(2.*apetx(i,j,jz)-apetx(i,j,jz-1))
               apety(i,j,jz+1)=(1-kp)*piano17(i,j) + kp*(2.*apety(i,j,jz)-apety(i,j,jz-1))
               apetz(i,j,jz+1)=(1-kp)*piano18(i,j) + kp*(2.*apetz(i,j,jz)-apetz(i,j,jz-1))
               apztx(i,j,jz+1)=(1-kp)*piano19(i,j) + kp*(2.*apztx(i,j,jz)-apztx(i,j,jz-1))
               apzty(i,j,jz+1)=(1-kp)*piano20(i,j) + kp*(2.*apzty(i,j,jz)-apzty(i,j,jz-1))
               apztz(i,j,jz+1)=(1-kp)*piano21(i,j) + kp*(2.*apztz(i,j,jz)-apztz(i,j,jz-1))


            end do
         end do

      endif


      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(7030+myid,*)uco(i,j,k),vco(i,j,k),wco(i,j,k)
                  write(7035+myid,*)uuco(i,j,k),uvco(i,j,k),uwco(i,j,k)
                  write(7040+myid,*)vuco(i,j,k),vvco(i,j,k),vwco(i,j,k)
                  write(7045+myid,*)apcsx(i,j,k),apcsy(i,j,k),apcsz(i,j,k)
                  write(7050+myid,*)apetx(i,j,k),apety(i,j,k),apetz(i,j,k)
                  write(7055+myid,*)apztx(i,j,k),apzty(i,j,k),apztz(i,j,k)
               end do
            end do
         end do
      end if

      ! bar filtering for scale similar

      call filterb_csieta(uco,m21,kparasta,kparaend,myid,nproc)
      call filterb_csieta(vco,m22,kparasta,kparaend,myid,nproc)
      call filterb_csieta(wco,m23,kparasta,kparaend,myid,nproc)
      !
      call filterb_csieta(uuco,l11,kparasta,kparaend,myid,nproc)
      call filterb_csieta(vuco,l12,kparasta,kparaend,myid,nproc)
      call filterb_csieta(wuco,l13,kparasta,kparaend,myid,nproc)
      call filterb_csieta(uvco,l21,kparasta,kparaend,myid,nproc)
      call filterb_csieta(vvco,l22,kparasta,kparaend,myid,nproc)
      call filterb_csieta(wvco,l23,kparasta,kparaend,myid,nproc)
      call filterb_csieta(uwco,l31,kparasta,kparaend,myid,nproc)
      call filterb_csieta(vwco,l32,kparasta,kparaend,myid,nproc)
      call filterb_csieta(wwco,l33,kparasta,kparaend,myid,nproc)
      !
      call filterb_csieta(apcsx,m11m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apcsy,m12m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apcsz,m13m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apetx,m21m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apety,m22m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apetz,m23m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apztx,m31m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apzty,m32m,kparasta,kparaend,myid,nproc)
      call filterb_csieta(apztz,m33m,kparasta,kparaend,myid,nproc)



      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(7065+myid,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                  !v      write(7070+myid,*)m31(i,j,k),m32(i,j,k),m33(i,j,k)
                  write(7075+myid,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                  write(7080+myid,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                  write(7085+myid,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                  write(7090+myid,*)m11m(i,j,k),m12m(i,j,k),m13m(i,j,k)
                  write(7095+myid,*)m21m(i,j,k),m22m(i,j,k),m23m(i,j,k)
                  write(7100+myid,*)m31m(i,j,k),m32m(i,j,k),m33m(i,j,k)
               end do
            end do
         end do
      end if

      ! now I need to exchange "ghost" planes to impose periodicity
      ! to conclude the filtering in zita (k+1 and k-1)
      !
      ! first kparasta
      !
      do m=1,40*(jx+2)*(jy+2)
         sbuff(m)=0.
         rbuff(m)=0.
      enddo

      call buffer1g(m21,1,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m22,2,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m23,3,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l11,4,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l12,5,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l13,6,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l21,7,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l22,8,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l23,9,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l31,10,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l32,11,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l33,12,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m11m,13,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m12m,14,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m13m,15,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m21m,16,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m22m,17,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m23m,18,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m31m,19,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m32m,20,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(m33m,21,kparasta,myid,nproc,kparasta,kparaend)
      !
      ! if I use pe and not pem I have implicit periodicity in k
      ! and the values in k=0 and k=jz+1 are defined in filter06
      !
      if (kp.eq.0) then

         call MPI_SENDRECV(sbuff1(1),21*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,tagls, &
            rbuff1(1),21*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrr, &
            MPI_COMM_WORLD,status,ierr)

      else if (kp.eq.1) then

         if(leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(sbuff1(1),21*(jx+2)*(jy+2), &
               MPI_REAL_SD,leftpem,tagls, &
               MPI_COMM_WORLD,ierr)
         endif
         if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),21*(jx+2)*(jy+2), &
               MPI_REAL_SD,rightpem,tagrr, &
               MPI_COMM_WORLD,status,ierr)
         endif

      endif
      !
      call buffer2(rbuff1,m21,1,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m22,2,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m23,3,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l11,4,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l12,5,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l13,6,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l21,7,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l22,8,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l23,9,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l31,10,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l32,11,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l33,12,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m11m,13,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m12m,14,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m13m,15,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m21m,16,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m22m,17,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m23m,18,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m31m,19,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m32m,20,kparaend+1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m33m,21,kparaend+1,myid,nproc,kparasta,kparaend)
      !
      ! now kparaend
      !
      do m=1,40*(jx+2)*(jy+2)
         sbuff(m)=0.
         rbuff(m)=0.
      enddo

      call buffer1g(m21,1,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m22,2,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m23,3,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l11,4,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l12,5,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l13,6,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l21,7,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l22,8,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l23,9,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l31,10,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l32,11,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l33,12,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m11m,13,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m12m,14,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m13m,15,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m21m,16,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m22m,17,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m23m,18,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m31m,19,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m32m,20,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m33m,21,kparaend,myid,nproc,kparasta,kparaend)

      if (kp.eq.0) then

         call MPI_SENDRECV(sbuff1(1),21*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrs, &
            rbuff1(1),21*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)

      else if (kp.eq.1) then

         if(rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(sbuff1(1),21*(jx+2)*(jy+2), &
               MPI_REAL_SD,rightpem,tagrs, &
               MPI_COMM_WORLD,ierr)
         endif
         if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),21*(jx+2)*(jy+2), &
               MPI_REAL_SD,leftpem,taglr, &
               MPI_COMM_WORLD,status,ierr)
         endif

      endif

      call buffer2(rbuff1,m21,1,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m22,2,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m23,3,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l11,4,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l12,5,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l13,6,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l21,7,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l22,8,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l23,9,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l31,10,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l32,11,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,l33,12,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m11m,13,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m12m,14,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m13m,15,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m21m,16,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m22m,17,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m23m,18,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m31m,19,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m32m,20,kparasta-1,myid,nproc,kparasta,kparaend)
      call buffer2(rbuff1,m33m,21,kparasta-1,myid,nproc,kparasta,kparaend)



      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(7110+myid,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                  !      write(7115+myid,*)m31(i,j,k),m32(i,j,k),m33(i,j,k)
                  write(7120+myid,*)l11(i,j,k),l12(i,j,k),l13(i,j,k)
                  write(7125+myid,*)l21(i,j,k),l22(i,j,k),l23(i,j,k)
                  write(7130+myid,*)l31(i,j,k),l32(i,j,k),l33(i,j,k)
                  write(7135+myid,*)m11m(i,j,k),m12m(i,j,k),m13m(i,j,k)
                  write(7140+myid,*)m21m(i,j,k),m22m(i,j,k),m23m(i,j,k)
                  write(7145+myid,*)m31m(i,j,k),m32m(i,j,k),m33m(i,j,k)
               end do
            end do
         end do
      end if
      !
      ! complete the filtering in zita
      !
      call filterb_zita(m21,ucof,kparasta,kparaend,myid,nproc)
      call filterb_zita(m22,vcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(m23,wcof,kparasta,kparaend,myid,nproc)
      !
      call filterb_zita(l11,uucof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l12,vucof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l13,wucof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l21,uvcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l22,vvcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l23,wvcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l31,uwcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l32,vwcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l33,wwcof,kparasta,kparaend,myid,nproc)
      !
      call filterb_zita(m11m,lmf11,kparasta,kparaend,myid,nproc)
      call filterb_zita(m12m,lmf12,kparasta,kparaend,myid,nproc)
      call filterb_zita(m13m,lmf13,kparasta,kparaend,myid,nproc)
      call filterb_zita(m21m,lmf21,kparasta,kparaend,myid,nproc)
      call filterb_zita(m22m,lmf22,kparasta,kparaend,myid,nproc)
      call filterb_zita(m23m,lmf23,kparasta,kparaend,myid,nproc)
      call filterb_zita(m31m,lmf31,kparasta,kparaend,myid,nproc)
      call filterb_zita(m32m,lmf32,kparasta,kparaend,myid,nproc)
      call filterb_zita(m33m,lmf33,kparasta,kparaend,myid,nproc)
      !
      ! cartesian component from the controvariant
      !
      do k=kparasta,kparaend
         do j=1,n2
            do i=1,n1
               pc1(i,j,k)  = ucof(i,j,k)
               pc2(i,j,k)  = vcof(i,j,k)
               pc3(i,j,k)  = wcof(i,j,k)
               ap11(i,j,k) = lmf11(i,j,k)
               ap12(i,j,k) = lmf12(i,j,k)
               ap13(i,j,k) = lmf13(i,j,k)
               ap21(i,j,k) = lmf21(i,j,k)
               ap22(i,j,k) = lmf22(i,j,k)
               ap23(i,j,k) = lmf23(i,j,k)
               ap31(i,j,k) = lmf31(i,j,k)
               ap32(i,j,k) = lmf32(i,j,k)
               ap33(i,j,k) = lmf33(i,j,k)
            end do
         end do
      end do

      
      call inverse_para2(myid,nproc,kparasta,kparaend)
      
      do k=kparasta,kparaend
         do j=1,n2
            do i=1,n1
               uf(i,j,k) = pp1(i,j,k)
               vf(i,j,k) = pp2(i,j,k)
               wf(i,j,k) = pp3(i,j,k)
            end do
         end do
      end do


      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(7150+myid,*)uf(i,j,k),vf(i,j,k),wf(i,j,k)
                  write(7155+myid,*)uucof(i,j,k),uvcof(i,j,k),uwcof(i,j,k)
                  write(7160+myid,*)vucof(i,j,k),vvcof(i,j,k),vwcof(i,j,k)
                  write(7165+myid,*)wucof(i,j,k),wvcof(i,j,k),wwcof(i,j,k)
               end do
            end do
         end do
      end if

      ! Leonard term Aik
      !
      if(myid.eq.0)then
         kparastam=kparasta+kp
         kparaendm=kparaend
      else if (myid.eq.nproc-1) then
         kparastam=kparasta
         kparaendm=kparaend-kp
      else
         kparastam=kparasta
         kparaendm=kparaend
      endif   
    
           
      do k=kparastam,kparaendm
         do j=1+jp,jy-jp
            do i=1+ip,jx-ip
               ass11(i,j,k)=uucof(i,j,k)-uf(i,j,k)*ucof(i,j,k)       !uU
               ass12(i,j,k)=vucof(i,j,k)-vf(i,j,k)*ucof(i,j,k)       !vU
               ass13(i,j,k)=wucof(i,j,k)-wf(i,j,k)*ucof(i,j,k)       !wU
               ass21(i,j,k)=uvcof(i,j,k)-uf(i,j,k)*vcof(i,j,k)       !uV
               ass22(i,j,k)=vvcof(i,j,k)-vf(i,j,k)*vcof(i,j,k)       !vV
               ass23(i,j,k)=wvcof(i,j,k)-wf(i,j,k)*vcof(i,j,k)       !wV
               ass31(i,j,k)=uwcof(i,j,k)-uf(i,j,k)*wcof(i,j,k)       !uW
               ass32(i,j,k)=vwcof(i,j,k)-vf(i,j,k)*wcof(i,j,k)       !vW
               ass33(i,j,k)=wwcof(i,j,k)-wf(i,j,k)*wcof(i,j,k)       !wW
            enddo
         enddo
      enddo
      !
      ! impose the value i=jx equal to jx-1 (if not periodic in i)
      !
      ! chicco perche' no per i=1 ???
      do ii=1,ip
      
         do k=kparasta,kparaend
            do j=1,jy
               ass11(jx,j,k)=ass11(jx-1,j,k)
               ass12(jx,j,k)=ass12(jx-1,j,k)
               ass13(jx,j,k)=ass13(jx-1,j,k)
               ass21(jx,j,k)=ass21(jx-1,j,k)
               ass22(jx,j,k)=ass22(jx-1,j,k)
               ass23(jx,j,k)=ass23(jx-1,j,k)
               ass31(jx,j,k)=ass31(jx-1,j,k)
               ass32(jx,j,k)=ass32(jx-1,j,k)
               ass33(jx,j,k)=ass33(jx-1,j,k)
            enddo
         enddo

      enddo
      
      ! impose the value j=jy equal to jj-1 (if not periodic in j)
      !
      do jj=1,jp
      
         do k=kparasta,kparaend
            do i=1,jx
               ass11(i,jy,k)=ass11(i,jy-1,k)
               ass12(i,jy,k)=ass12(i,jy-1,k)
               ass13(i,jy,k)=ass13(i,jy-1,k)
               ass21(i,jy,k)=ass21(i,jy-1,k)
               ass22(i,jy,k)=ass22(i,jy-1,k)
               ass23(i,jy,k)=ass23(i,jy-1,k)
               ass31(i,jy,k)=ass31(i,jy-1,k)
               ass32(i,jy,k)=ass32(i,jy-1,k)
               ass33(i,jy,k)=ass33(i,jy-1,k)
            enddo
         enddo
      
      enddo
      ! impose the value k=jz equal to jz-1 (if not periodic in k)
      !
      do kk=1,kp
      
         if(myid.eq.nproc-1)then
      
            do j=1,jy
               do i=1,jx
                  ass11(i,j,jz)=ass11(i,j,jz-1)
                  ass12(i,j,jz)=ass12(i,j,jz-1)
                  ass13(i,j,jz)=ass13(i,j,jz-1)
                  ass21(i,j,jz)=ass21(i,j,jz-1)
                  ass22(i,j,jz)=ass22(i,j,jz-1)
                  ass23(i,j,jz)=ass23(i,j,jz-1)
                  ass31(i,j,jz)=ass31(i,j,jz-1)
                  ass32(i,j,jz)=ass32(i,j,jz-1)
                  ass33(i,j,jz)=ass33(i,j,jz-1)
               enddo
            enddo
      
         endif
      
      enddo
      
      !
      ! periodicity Aik term
      !
      call periodic(ass11,ass12,ass13,myid,nproc,kparasta,kparaend)
      call periodic(ass21,ass22,ass23,myid,nproc,kparasta,kparaend)
      call periodic(ass31,ass32,ass33,myid,nproc,kparasta,kparaend)
      !
      ! periodicity in zita
      !
      if(kp.eq.0)then
         !
         ! here I need to exchange planes k=0 and k=jz+1
         ! so that P0 knows k=jz and Pn-1 knows k=1
         !
         do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
         enddo

         if (myid.eq.nproc-1) then

            call buffer1g(ass11,1,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass12,2,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass13,3,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass21,4,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass22,5,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass23,6,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass31,7,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass32,8,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass33,9,jz,myid,nproc,kparasta,kparaend)

         else if (myid.eq.0) then

            call buffer1g(ass11,1,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass12,2,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass13,3,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass21,4,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass22,5,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass23,6,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass31,7,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass32,8,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass33,9,1,myid,nproc,kparasta,kparaend)

         endif

         ! now the exchange so that P0 will know k=jz and Pn-1 k=1

         if (myid.eq.nproc-1) then

            call MPI_SENDRECV(sbuff1(1),9*(jx+2)*(jy+2),MPI_REAL_SD,0,951, &
               rbuff1(1),9*(jx+2)*(jy+2),MPI_REAL_SD,0,851, &
               MPI_COMM_WORLD,status,ierr)

         else if (myid.eq.0) then

            call MPI_SENDRECV( &
               sbuff1(1),9*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,851, &
               rbuff1(1),9*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,951, &
               MPI_COMM_WORLD,status,ierr)

         endif

         if (myid.eq.0) then

            call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano4,4,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano5,5,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano6,6,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano7,7,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano8,8,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano9,9,myid,nproc,kparasta,kparaend)


         else if (myid.eq.nproc-1) then

            call buffer2gg(piano1,1,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano2,2,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano3,3,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano4,4,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano5,5,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano6,6,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano7,7,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano8,8,myid,nproc,kparasta,kparaend)
            call buffer2gg(piano9,9,myid,nproc,kparasta,kparaend)


         endif

         ! now P0 knows the plane k=jz and Pn-1 knows k=1

         if(myid.eq.0)then
            do j=1,jy
               do i=1,jx

                  ass11(i,j,0)=piano1(i,j)
                  ass12(i,j,0)=piano2(i,j)
                  ass13(i,j,0)=piano3(i,j)
                  ass21(i,j,0)=piano4(i,j)
                  ass22(i,j,0)=piano5(i,j)
                  ass23(i,j,0)=piano6(i,j)
                  ass31(i,j,0)=piano7(i,j)
                  ass32(i,j,0)=piano8(i,j)
                  ass33(i,j,0)=piano9(i,j)

               enddo
            enddo
         endif
         !
         if(myid.eq.nproc-1)then
            do j=1,jy
               do i=1,jx

                  ass11(i,j,jz+1)=piano1(i,j)
                  ass12(i,j,jz+1)=piano2(i,j)
                  ass13(i,j,jz+1)=piano3(i,j)
                  ass21(i,j,jz+1)=piano4(i,j)
                  ass22(i,j,jz+1)=piano5(i,j)
                  ass23(i,j,jz+1)=piano6(i,j)
                  ass31(i,j,jz+1)=piano7(i,j)
                  ass32(i,j,jz+1)=piano8(i,j)
                  ass33(i,j,jz+1)=piano9(i,j)

               enddo
            enddo
         endif
      !
      endif
      !
      ! finally the computation of cgra123 and I need to pass the
      ! plane kparaend+1 of ass31, ass32 and ass33 to compute cgra3
      !
      if(iq.eq.1)then

         if(leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(ass31(1,1,kparasta),(jx+2)*(jy+2), &
               MPI_REAL_SD,leftpem,tagls, &
               MPI_COMM_WORLD,ierr)
         endif
         if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(ass31(1,1,kparaend+1),(jx+2)*(jy+2), &
               MPI_REAL_SD,rightpem,tagrr, &
               MPI_COMM_WORLD,status,ierr)
         endif

      else if(iq.eq.2)then

         if(leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(ass32(1,1,kparasta),(jx+2)*(jy+2), &
               MPI_REAL_SD,leftpem,tagls, &
               MPI_COMM_WORLD,ierr)
         endif
         if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(ass32(1,1,kparaend+1),(jx+2)*(jy+2), &
               MPI_REAL_SD,rightpem,tagrr, &
               MPI_COMM_WORLD,status,ierr)
         endif

      else if(iq.eq.3)then
 
         if(leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(ass33(1,1,kparasta),(jx+2)*(jy+2), &
               MPI_REAL_SD,leftpem,tagls, &
               MPI_COMM_WORLD,ierr)
         endif
         if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(ass33(1,1,kparaend+1),(jx+2)*(jy+2), &
               MPI_REAL_SD,rightpem,tagrr, &
               MPI_COMM_WORLD,status,ierr)
         endif

      endif
      !-----------------------------------------------------------------------
      ! compute sgs stress for scale similar part of momentum
      ! accumulation to compute sgs stress

      stampa: if(ktime.eq. i_print*(ktime/i_print) .and. iq.eq.1 .and. lagr.eq.0) then

         somma=float(jx)*float(jz)

         allocate(fil11(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil12(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil13(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil21(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil22(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil23(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil31(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil32(0:n1+1,0:n2+1,0:n3+1))
         allocate(fil33(0:n1+1,0:n2+1,0:n3+1))
         allocate(sus_loc11(n2))
         allocate(sus_loc12(n2))
         allocate(sus_loc13(n2))
         allocate(sus_loc21(n2))
         allocate(sus_loc22(n2))
         allocate(sus_loc23(n2))
         allocate(sus_loc31(n2))
         allocate(sus_loc32(n2))
         allocate(sus_loc33(n2))
         !
         ! compute sgs (scale similar) cartesian stress momentum

         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  pc1(i,j,k)  = ass11(i,j,k)
                  pc2(i,j,k)  = ass12(i,j,k)
                  pc3(i,j,k)  = ass13(i,j,k)
                  ap11(i,j,k) = apcsx(i,j,k)
                  ap12(i,j,k) = apcsy(i,j,k)
                  ap13(i,j,k) = apcsz(i,j,k)
                  ap21(i,j,k) = apetx(i,j,k)
                  ap22(i,j,k) = apety(i,j,k)
                  ap23(i,j,k) = apetz(i,j,k)
                  ap31(i,j,k) = apztx(i,j,k)
                  ap32(i,j,k) = apzty(i,j,k)
                  ap33(i,j,k) = apztz(i,j,k)
               end do
            end do
         end do

         call inverse_para2(myid,nproc,kparasta,kparaend)
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  fil11(i,j,k) = pp1(i,j,k)
                  fil12(i,j,k) = pp2(i,j,k)
                  fil13(i,j,k) = pp3(i,j,k)
               end do
            end do
         end do

         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  pc1(i,j,k)  = ass21(i,j,k)
                  pc2(i,j,k)  = ass22(i,j,k)
                  pc3(i,j,k)  = ass23(i,j,k)
                  ap11(i,j,k) = apcsx(i,j,k)
                  ap12(i,j,k) = apcsy(i,j,k)
                  ap13(i,j,k) = apcsz(i,j,k)
                  ap21(i,j,k) = apetx(i,j,k)
                  ap22(i,j,k) = apety(i,j,k)
                  ap23(i,j,k) = apetz(i,j,k)
                  ap31(i,j,k) = apztx(i,j,k)
                  ap32(i,j,k) = apzty(i,j,k)
                  ap33(i,j,k) = apztz(i,j,k)
               end do
            end do
         end do

         call inverse_para2(myid,nproc,kparasta,kparaend)
  
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  fil21(i,j,k) = pp1(i,j,k)
                  fil22(i,j,k) = pp2(i,j,k)
                  fil23(i,j,k) = pp3(i,j,k)
               end do
            end do
         end do

         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  pc1(i,j,k)  = ass31(i,j,k)
                  pc2(i,j,k)  = ass32(i,j,k)
                  pc3(i,j,k)  = ass33(i,j,k)
                  ap11(i,j,k) = apcsx(i,j,k)
                  ap12(i,j,k) = apcsy(i,j,k)
                  ap13(i,j,k) = apcsz(i,j,k)
                  ap21(i,j,k) = apetx(i,j,k)
                  ap22(i,j,k) = apety(i,j,k)
                  ap23(i,j,k) = apetz(i,j,k)
                  ap31(i,j,k) = apztx(i,j,k)
                  ap32(i,j,k) = apzty(i,j,k)
                  ap33(i,j,k) = apztz(i,j,k)
               end do
            end do
         end do

         call inverse_para2(myid,nproc,kparasta,kparaend)
  
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  fil31(i,j,k) = pp1(i,j,k)
                  fil32(i,j,k) = pp2(i,j,k)
                  fil33(i,j,k) = pp3(i,j,k)
               end do
            end do
         end do

         ! average on homogeneus plane
         !
         do j=1,jy

            sus_loc11(j)=0.
            sus_loc12(j)=0.
            sus_loc13(j)=0.
            sus_loc22(j)=0.
            sus_loc23(j)=0.
            sus_loc33(j)=0.

            do k=kparasta,kparaend
               do i=1,jx
                  sus_loc11(j)=sus_loc11(j)+fil11(i,j,k)
                  sus_loc22(j)=sus_loc22(j)+fil22(i,j,k)
                  sus_loc33(j)=sus_loc33(j)+fil33(i,j,k)
                  sus_loc12(j)=sus_loc12(j)+.5*(fil12(i,j,k)+fil21(i,j,k))
                  sus_loc13(j)=sus_loc13(j)+.5*(fil13(i,j,k)+fil31(i,j,k))
                  sus_loc23(j)=sus_loc23(j)+.5*(fil23(i,j,k)+fil32(i,j,k))
               enddo
            enddo
         enddo

         ! now the sum of local sum and send to all procs for each plane
         !
         do m=1,40*(jx+2)*(jy+2)
            sbuff(m)=0.
            rbuff(m)=0.
         enddo

         call buffvect1d(sbuff,sus_loc11,1)
         call buffvect1d(sbuff,sus_loc12,2)
         call buffvect1d(sbuff,sus_loc13,3)
         call buffvect1d(sbuff,sus_loc22,4)
         call buffvect1d(sbuff,sus_loc23,5)
         call buffvect1d(sbuff,sus_loc33,6)

         call MPI_ALLREDUCE(sbuff(1),rbuff(1),6*jy,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

         call buffvect2d(rbuff,sus11,1)
         call buffvect2d(rbuff,sus12,2)
         call buffvect2d(rbuff,sus13,3)
         call buffvect2d(rbuff,sus22,4)
         call buffvect2d(rbuff,sus23,5)
         call buffvect2d(rbuff,sus33,6)
         !
         ! finally the average value, therefore the sum divided by the number of points
         ! on j-plane (known by all procs)
         !
         do j=1,jy
            sus11(j)=sus11(j)/somma
            sus22(j)=sus22(j)/somma
            sus33(j)=sus33(j)/somma
            sus12(j)=sus12(j)/somma
            sus13(j)=sus13(j)/somma
            sus23(j)=sus23(j)/somma
         enddo

         deallocate(fil11)
         deallocate(fil12)
         deallocate(fil13)
         deallocate(fil21)
         deallocate(fil22)
         deallocate(fil23)
         deallocate(fil31)
         deallocate(fil32)
         deallocate(fil33)
         deallocate(sus_loc11)
         deallocate(sus_loc12)
         deallocate(sus_loc13)
         deallocate(sus_loc21)
         deallocate(sus_loc22)
         deallocate(sus_loc23)
         deallocate(sus_loc31)
         deallocate(sus_loc32)
         deallocate(sus_loc33)
      end if stampa

      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(7170+myid,*)ass11(i,j,k),ass12(i,j,k),ass13(i,j,k)
                  write(7175+myid,*)ass21(i,j,k),ass22(i,j,k),ass23(i,j,k)
                  write(7180+myid,*)ass31(i,j,k),ass32(i,j,k),ass33(i,j,k)
               end do
            end do
         end do
      end if

      if(myid.eq.0)then
         kparastan=2*kp
      else
         kparastan=kparasta
      endif
      
      if (iq.eq.1) then
         !
         do k=kparastam,kparaendm
            ! chicco this loop better if insert in if, also for the next loops
            do j=1+jp,jy-jp
               do i=2*ip,jx-ip
                  cgra1(i,j,k)=.5*(ass11(i,j,k)+ass11(i+1,j,k))
               enddo
            enddo
         enddo
         !
         do k=kparastam,kparaendm
            do j=2*jp,jy-jp
               do i=1+ip,jx-ip
                  cgra2(i,j,k)=.5*(ass21(i,j,k)+ass21(i,j+1,k))
               enddo
            enddo
         enddo
         !
         do k=kparastan,kparaendm
            do j=1+jp,jy-jp
               do i=1+ip,jx-ip
                  cgra3(i,j,k)=.5*(ass31(i,j,k)+ass31(i,j,k+1))
               enddo
            enddo
         enddo
      !
      else if (iq.eq.2) then
         !
         do k=kparastam,kparaendm
            do j=1+jp,jy-jp
               do i=2*ip,jx-ip
                  cgra1(i,j,k)=.5*(ass12(i,j,k)+ass12(i+1,j,k))
               enddo
            enddo
         enddo
         !
         do k=kparastam,kparaendm
            do j=2*jp,jy-jp
               do i=1+ip,jx-ip
                  cgra2(i,j,k)=.5*(ass22(i,j,k)+ass22(i,j+1,k))
               enddo
            enddo
         enddo
         !
         do k=kparastan,kparaendm
            do j=1+jp,jy-jp
               do i=1+ip,jx-ip
                  cgra3(i,j,k)=.5*(ass32(i,j,k)+ass32(i,j,k+1))
               enddo
            enddo
         enddo

      !
      else if (iq.eq.3) then
         !
         do k=kparastam,kparaendm
            do j=1+jp,jy-jp
               do i=2*ip,jx-ip
                  cgra1(i,j,k)=.5*(ass13(i,j,k)+ass13(i+1,j,k))
               enddo
            enddo
         enddo
         !
         do k=kparastam,kparaendm
            do j=2*jp,jy-jp
               do i=1+ip,jx-ip
                  cgra2(i,j,k)=.5*(ass23(i,j,k)+ass23(i,j+1,k))
               enddo
            enddo
         enddo
         !
         do k=kparastan,kparaendm
            do j=1+jp,jy-jp
               do i=1+ip,jx-ip
                  cgra3(i,j,k)=.5*(ass33(i,j,k)+ass33(i,j,k+1))
               enddo
            enddo
         enddo

      endif
   !
   endif
   !
   return
end subroutine mix_para
