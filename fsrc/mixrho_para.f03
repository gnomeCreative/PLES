!***********************************************************************
subroutine mixrho_para(inmod,rho)
   !***********************************************************************
   ! compute scale similar part for scalar equation
   !
   use filter_module
   use mysending
   use turbo_module
   use myarrays_velo3
   use myarrays_metri3
   !
   use scala3
   use period
   !
   use mpi


   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr
   integer kparastam,kparaendm
   integer kparastan
   integer status(MPI_STATUS_SIZE)
   integer req1,req2
   integer istatus(MPI_STATUS_SIZE)
   !
   real sbuff((n1+2)*(n2+2)*40)
   real rbuff((n1+2)*(n2+2)*40)
   real rho(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   !
   integer i,j,k,kparastal,inmod,m,lll
   integer ii,jj,kk
   integer debugg
   real aden
   !-----------------------------------------------------------------------

   debugg = 0
   !
   aden=float(jx)*float(jz)

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
   !-------------------------------------------------
   ! compute controvariant term for product rho*velocity
   ! only if scale similar
   !
   if (inmod.eq.1) then

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
               ! compute rhoU rhoV rhoW
               ! and rho, to filter
               !
               uuco(i,j,k)=uco(i,j,k)*rho(i,j,k)
               vvco(i,j,k)=vco(i,j,k)*rho(i,j,k)
               wwco(i,j,k)=wco(i,j,k)*rho(i,j,k)

               rhofl(i,j,k)=rho(i,j,k)

            enddo
         enddo
      enddo
      !
      ! generalized periodicity
      !
      ! extrapolation for velocity on side i=0 and i=jx+1
      !
      do k=kparasta,kparaend
         do j=1,jy

            uco(0,j,k)=(1-ip)*uco(jx,j,k) + &
               ip*(2.*uco(1,j,k)-uco(2,j,k))
            vco(0,j,k)=(1-ip)*vco(jx,j,k) + &
               ip*(2.*vco(1,j,k)-vco(2,j,k))
            wco(0,j,k)=(1-ip)*wco(jx,j,k) + &
               ip*(2.*wco(1,j,k)-wco(2,j,k))
            !
            uuco(0,j,k)=(1-ip)*uuco(jx,j,k) + &
               ip*(2.*uuco(1,j,k)-uuco(2,j,k))
            vvco(0,j,k)=(1-ip)*vvco(jx,j,k) + &
               ip*(2.*vvco(1,j,k)-vvco(2,j,k))
            wwco(0,j,k)=(1-ip)*wwco(jx,j,k) + &
               ip*(2.*wwco(1,j,k)-wwco(2,j,k))
            rhofl(0,j,k)=(1-ip)*rhofl(jx,j,k) + &
               ip*(2.*rhofl(1,j,k)-rhofl(2,j,k))
            !
            apcsx(0,j,k)=(1-ip)*apcsx(jx,j,k) + &
               ip*(2.*apcsx(1,j,k)-apcsx(2,j,k))
            apcsy(0,j,k)=(1-ip)*apcsy(jx,j,k) + &
               ip*(2.*apcsy(1,j,k)-apcsy(2,j,k))
            apcsz(0,j,k)=(1-ip)*apcsz(jx,j,k) + &
               ip*(2.*apcsz(1,j,k)-apcsz(2,j,k))
            apetx(0,j,k)=(1-ip)*apetx(jx,j,k) + &
               ip*(2.*apetx(1,j,k)-apetx(2,j,k))
            apety(0,j,k)=(1-ip)*apety(jx,j,k) + &
               ip*(2.*apety(1,j,k)-apety(2,j,k))
            apetz(0,j,k)=(1-ip)*apetz(jx,j,k) + &
               ip*(2.*apetz(1,j,k)-apetz(2,j,k))
            apztx(0,j,k)=(1-ip)*apztx(jx,j,k) + &
               ip*(2.*apztx(1,j,k)-apztx(2,j,k))
            apzty(0,j,k)=(1-ip)*apzty(jx,j,k) + &
               ip*(2.*apzty(1,j,k)-apzty(2,j,k))
            apztz(0,j,k)=(1-ip)*apztz(jx,j,k) + &
               ip*(2.*apztz(1,j,k)-apztz(2,j,k))
            !
            !
            apcsx(jx+1,j,k)=(1-ip)*apcsx(1,j,k) + &
               ip*(2.*apcsx(jx,j,k)-apcsx(jx-1,j,k))
            apcsy(jx+1,j,k)=(1-ip)*apcsy(1,j,k) + &
               ip*(2.*apcsy(jx,j,k)-apcsy(jx-1,j,k))
            apcsz(jx+1,j,k)=(1-ip)*apcsz(1,j,k) + &
               ip*(2.*apcsz(jx,j,k)-apcsz(jx-1,j,k))
            apetx(jx+1,j,k)=(1-ip)*apetx(1,j,k) + &
               ip*(2.*apetx(jx,j,k)-apetx(jx-1,j,k))
            apety(jx+1,j,k)=(1-ip)*apety(1,j,k) + &
               ip*(2.*apety(jx,j,k)-apety(jx-1,j,k))
            apetz(jx+1,j,k)=(1-ip)*apetz(1,j,k) + &
               ip*(2.*apetz(jx,j,k)-apetz(jx-1,j,k))
            apztx(jx+1,j,k)=(1-ip)*apztx(1,j,k) + &
               ip*(2.*apztx(jx,j,k)-apztx(jx-1,j,k))
            apzty(jx+1,j,k)=(1-ip)*apzty(1,j,k) + &
               ip*(2.*apzty(jx,j,k)-apzty(jx-1,j,k))
            apztz(jx+1,j,k)=(1-ip)*apztz(1,j,k) + &
               ip*(2.*apztz(jx,j,k)-apztz(jx-1,j,k))
            !
            uco(jx+1,j,k)=(1-ip)*uco(1,j,k) + &
               ip*(2.*uco(jx,j,k)-uco(jx-1,j,k))
            vco(jx+1,j,k)=(1-ip)*vco(1,j,k) + &
               ip*(2.*vco(jx,j,k)-vco(jx-1,j,k))
            wco(jx+1,j,k)=(1-ip)*wco(1,j,k) + &
               ip*(2.*wco(jx,j,k)-wco(jx-1,j,k))
            !
            uuco(jx+1,j,k)=(1-ip)*uuco(1,j,k) + &
               ip*(2.*uuco(jx,j,k)-uuco(jx-1,j,k))
            vvco(jx+1,j,k)=(1-ip)*vvco(1,j,k) + &
               ip*(2.*vvco(jx,j,k)-vvco(jx-1,j,k))
            wwco(jx+1,j,k)=(1-ip)*wwco(1,j,k) + &
               ip*(2.*wwco(jx,j,k)-wwco(jx-1,j,k))
            rhofl(jx+1,j,k)=(1-ip)*rhofl(1,j,k) + &
               ip*(2.*rhofl(jx,j,k)-rhofl(jx-1,j,k))

         enddo
      enddo
      !
      ! extrapolation for velocity on side j=0 j=jy+1
      !
      do k=kparasta,kparaend
         do i=1,jx

            uco(i,0,k)=(1-jp)*uco(i,jy,k) + &
               jp*(2.*uco(i,1,k)-uco(i,2,k))
            vco(i,0,k)=(1-jp)*vco(i,jy,k) + &
               jp*(2.*vco(i,1,k)-vco(i,2,k))
            wco(i,0,k)=(1-jp)*wco(i,jy,k) + &
               jp*(2.*wco(i,1,k)-wco(i,2,k))
            !
            uuco(i,0,k)=(1-jp)*uuco(i,jy,k) + &
               jp*(2.*uuco(i,1,k)-uuco(i,2,k))
            vvco(i,0,k)=(1-jp)*vvco(i,jy,k) + &
               jp*(2.*vvco(i,1,k)-vvco(i,2,k))
            wwco(i,0,k)=(1-jp)*wwco(i,jy,k) + &
               jp*(2.*wwco(i,1,k)-wwco(i,2,k))
            rhofl(i,0,k)=(1-jp)*rhofl(i,jy,k) + &
               jp*(2.*rhofl(i,1,k)-rhofl(i,2,k))
            !
            apcsx(i,0,k)=(1-jp)*apcsx(i,jy,k) + &
               jp*(2.*apcsx(i,1,k)-apcsx(i,2,k))
            apcsy(i,0,k)=(1-jp)*apcsy(i,jy,k) + &
               jp*(2.*apcsy(i,1,k)-apcsy(i,2,k))
            apcsz(i,0,k)=(1-jp)*apcsz(i,jy,k) + &
               jp*(2.*apcsz(i,1,k)-apcsz(i,2,k))
            apetx(i,0,k)=(1-jp)*apetx(i,jy,k) + &
               jp*(2.*apetx(i,1,k)-apetx(i,2,k))
            apety(i,0,k)=(1-jp)*apety(i,jy,k) + &
               jp*(2.*apety(i,1,k)-apety(i,2,k))
            apetz(i,0,k)=(1-jp)*apetz(i,jy,k) + &
               jp*(2.*apetz(i,1,k)-apetz(i,2,k))
            apztx(i,0,k)=(1-jp)*apztx(i,jy,k) + &
               jp*(2.*apztx(i,1,k)-apztx(i,2,k))
            apzty(i,0,k)=(1-jp)*apzty(i,jy,k) + &
               jp*(2.*apzty(i,1,k)-apzty(i,2,k))
            apztz(i,0,k)=(1-jp)*apztz(i,jy,k) + &
               jp*(2.*apztz(i,1,k)-apztz(i,2,k))
            !
            !
            apcsx(i,jy+1,k)=(1-jp)*apcsx(i,1,k) + &
               jp*(2.*apcsx(i,jy,k)-apcsx(i,jy-1,k))
            apcsy(i,jy+1,k)=(1-jp)*apcsy(i,1,k) + &
               jp*(2.*apcsy(i,jy,k)-apcsy(i,jy-1,k))
            apcsz(i,jy+1,k)=(1-jp)*apcsz(i,1,k) + &
               jp*(2.*apcsz(i,jy,k)-apcsz(i,jy-1,k))
            apetx(i,jy+1,k)=(1-jp)*apetx(i,1,k) + &
               jp*(2.*apetz(i,jy,k)-apetx(i,jy-1,k))
            apety(i,jy+1,k)=(1-jp)*apety(i,1,k) + &
               jp*(2.*apety(i,jy,k)-apety(i,jy-1,k))
            apetz(i,jy+1,k)=(1-jp)*apetz(i,1,k) + &
               jp*(2.*apetz(i,jy,k)-apetz(i,jy-1,k))
            apztx(i,jy+1,k)=(1-jp)*apztx(i,1,k) + &
               jp*(2.*apztx(i,jy,k)-apztx(i,jy-1,k))
            apzty(i,jy+1,k)=(1-jp)*apzty(i,1,k) + &
               jp*(2.*apzty(i,jy,k)-apzty(i,jy-1,k))
            apztz(i,jy+1,k)=(1-jp)*apztz(i,1,k) + &
               jp*(2.*apztz(i,jy,k)-apztz(i,jy-1,k))
            !
            uco(i,jy+1,k)=(1-jp)*uco(i,1,k) + &
               jp*(2.*uco(i,jy,k)-uco(i,jy-1,k))
            vco(i,jy+1,k)=(1-jp)*vco(i,1,k) + &
               jp*(2.*vco(i,jy,k)-vco(i,jy-1,k))
            wco(i,jy+1,k)=(1-jp)*wco(i,1,k) + &
               jp*(2.*wco(i,jy,k)-wco(i,jy-1,k))
            !
            uuco(i,jy+1,k)=(1-jp)*uuco(i,1,k) + &
               jp*(2.*uuco(i,jy,k)-uuco(i,jy-1,k))
            vvco(i,jy+1,k)=(1-jp)*vvco(i,1,k) + &
               jp*(2.*vvco(i,jy,k)-vvco(i,jy-1,k))
            wwco(i,jy+1,k)=(1-jp)*wwco(i,1,k) + &
               jp*(2.*wwco(i,jy,k)-wwco(i,jy-1,k))
            rhofl(i,jy+1,k)=(1-jp)*rhofl(i,1,k) + &
               jp*(2.*rhofl(i,jy,k)-rhofl(i,jy-1,k))

         enddo
      enddo
      !
      ! periodicity on sides front and back
      !
      ! extrapolation for velocity on sides k=0 and k=jz+1
      !
      ! I use a sending buffer to contain all the plane jx*jy
      ! (16) for the exchange between P0 and Pn-1

      do m=1,40*(jx+2)*(jy+2)
         sbuff(m)=0.
         rbuff(m)=0.
      enddo


      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(6000+myid,*)uco(i,j,k),vco(i,j,k),wco(i,j,k)
                  write(6005+myid,*)uuco(i,j,k),vvco(i,j,k),wwco(i,j,k)
                  write(6010+myid,*)rhofl(i,j,k)
                  write(6015+myid,*)apcsx(i,j,k),apcsy(i,j,k),apcsz(i,j,k)
                  write(6020+myid,*)apetx(i,j,k),apety(i,j,k),apetz(i,j,k)
                  write(6015+myid,*)apztx(i,j,k),apzty(i,j,k),apztz(i,j,k)
               end do
            end do
         end do
      end if


      if (myid.eq.nproc-1) then

         call buffer1g(uco,1,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(vco,2,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(wco,3,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(uuco,4,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(vvco,5,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(wwco,6,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(rhofl,7,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsx,8,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsy,9,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsz,10,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apetx,11,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apety,12,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apetz,13,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apztx,14,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apzty,15,jz,myid,nproc,kparasta,kparaend)
         call buffer1g(apztz,16,jz,myid,nproc,kparasta,kparaend)

      else if (myid.eq.0) then

         call buffer1g(uco,1,1,myid,nproc,kparasta,kparaend)
         call buffer1g(vco,2,1,myid,nproc,kparasta,kparaend)
         call buffer1g(wco,3,1,myid,nproc,kparasta,kparaend)
         call buffer1g(uuco,4,1,myid,nproc,kparasta,kparaend)
         call buffer1g(vvco,5,1,myid,nproc,kparasta,kparaend)
         call buffer1g(wwco,6,1,myid,nproc,kparasta,kparaend)
         call buffer1g(rhofl,7,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsx,8,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsy,9,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apcsz,10,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apetx,11,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apety,12,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apetz,13,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apztx,14,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apzty,15,1,myid,nproc,kparasta,kparaend)
         call buffer1g(apztz,16,1,myid,nproc,kparasta,kparaend)

      endif
      !
      ! now exchange so that P0 knows k=jz and Pn-1 knows k=1
      !
      if (myid.eq.nproc-1) then

         call MPI_SENDRECV(sbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,0,9101, &
            rbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,0,8101, &
            MPI_COMM_WORLD,status,ierr)

      else if (myid.eq.0) then

         call MPI_SENDRECV(sbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,nproc-1,8101, &
            rbuff1(1),16*(jx+2)*(jy+2), &
            MPI_REAL_SD,nproc-1,9101, &
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
         call buffer2gg(piano10,10,myid,nproc, &
            kparasta,kparaend)
         call buffer2gg(piano11,11,myid,nproc, &
            kparasta,kparaend)
         call buffer2gg(piano12,12,myid,nproc, &
            kparasta,kparaend)
         call buffer2gg(piano13,13,myid,nproc, &
            kparasta,kparaend)
         call buffer2gg(piano14,14,myid,nproc, &
            kparasta,kparaend)
         call buffer2gg(piano15,15,myid,nproc, &
            kparasta,kparaend)
         call buffer2gg(piano16,16,myid,nproc, &
            kparasta,kparaend)

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

      endif

      ! now P0 knows plane k=jz
      ! and Pn-1 knows plane k=1

      if(myid.eq.0)then

         do j=1,jy
            do i=1,jx
               !
               uco(i,j,0)=(1-kp)*piano1(i,j) + &
                  kp*(2.*uco(i,j,1)-uco(i,j,2))
               vco(i,j,0)=(1-kp)*piano2(i,j) + &
                  kp*(2.*vco(i,j,1)-vco(i,j,2))
               wco(i,j,0)=(1-kp)*piano3(i,j) + &
                  kp*(2.*wco(i,j,1)-wco(i,j,2))
               !
               uuco(i,j,0)=(1-kp)*piano4(i,j) + &
                  kp*(2.*uuco(i,j,1)-uuco(i,j,2))
               vvco(i,j,0)=(1-kp)*piano5(i,j) + &
                  kp*(2.*vvco(i,j,1)-vvco(i,j,2))
               wwco(i,j,0)=(1-kp)*piano6(i,j) + &
                  kp*(2.*wwco(i,j,1)-wwco(i,j,2))
               rhofl(i,j,0)=(1-kp)*piano7(i,j) + &
                  kp*(2.*rhofl(i,j,1)-rhofl(i,j,2))
               !
               apcsx(i,j,0)=(1-kp)*piano8(i,j) + &
                  kp*(2.*apcsx(i,j,1)-apcsx(i,j,2))
               apcsy(i,j,0)=(1-kp)*piano9(i,j) + &
                  kp*(2.*apcsy(i,j,1)-apcsy(i,j,2))
               apcsz(i,j,0)=(1-kp)*piano19(i,j) + &
                  kp*(2.*apcsz(i,j,1)-apcsz(i,j,2))
               apetx(i,j,0)=(1-kp)*piano11(i,j) + &
                  kp*(2.*apetx(i,j,1)-apetx(i,j,2))
               apety(i,j,0)=(1-kp)*piano12(i,j) + &
                  kp*(2.*apety(i,j,1)-apety(i,j,2))
               apetz(i,j,0)=(1-kp)*piano13(i,j) + &
                  kp*(2.*apetz(i,j,1)-apetz(i,j,2))
               apztx(i,j,0)=(1-kp)*piano14(i,j) + &
                  kp*(2.*apztx(i,j,1)-apztx(i,j,2))
               apzty(i,j,0)=(1-kp)*piano15(i,j) + &
                  kp*(2.*apzty(i,j,1)-apzty(i,j,2))
               apztz(i,j,0)=(1-kp)*piano16(i,j) + &
                  kp*(2.*apztz(i,j,1)-apztz(i,j,2))
            !
            end do
         end do

      endif
      !
      !
      if(myid.eq.nproc-1)then

         do j=1,jy
            do i=1,jx
               !
               uco(i,j,jz+1)=(1-kp)*piano1(i,j) + &
                  kp*(2.*uco(i,j,jz)-uco(i,j,jz-1))
               vco(i,j,jz+1)=(1-kp)*piano2(i,j) + &
                  kp*(2.*vco(i,j,jz)-vco(i,j,jz-1))
               wco(i,j,jz+1)=(1-kp)*piano3(i,j) + &
                  kp*(2.*wco(i,j,jz)-wco(i,j,jz-1))
               !
               uuco(i,j,jz+1)=(1-kp)*piano4(i,j) + &
                  kp*(2.*uuco(i,j,jz)-uuco(i,j,jz-1))
               vvco(i,j,jz+1)=(1-kp)*piano5(i,j) + &
                  kp*(2.*vvco(i,j,jz)-vvco(i,j,jz-1))
               wwco(i,j,jz+1)=(1-kp)*piano6(i,j) + &
                  kp*(2.*wwco(i,j,jz)-wwco(i,j,jz-1))
               rhofl(i,j,jz+1)=(1-kp)*piano7(i,j) + &
                  kp*(2.*rhofl(i,j,jz)-rhofl(i,j,jz-1))

               apcsx(i,j,jz+1)=(1-kp)*piano8(i,j) + &
                  kp*(2.*apcsx(i,j,jz)-apcsx(i,j,jz-1))
               apcsy(i,j,jz+1)=(1-kp)*piano9(i,j) + &
                  kp*(2.*apcsy(i,j,jz)-apcsy(i,j,jz-1))
               apcsz(i,j,jz+1)=(1-kp)*piano10(i,j) + &
                  kp*(2.*apcsz(i,j,jz)-apcsz(i,j,jz-1))
               apetx(i,j,jz+1)=(1-kp)*piano11(i,j) + &
                  kp*(2.*apetx(i,j,jz)-apetx(i,j,jz-1))
               apety(i,j,jz+1)=(1-kp)*piano12(i,j) + &
                  kp*(2.*apety(i,j,jz)-apety(i,j,jz-1))
               apetz(i,j,jz+1)=(1-kp)*piano13(i,j) + &
                  kp*(2.*apetz(i,j,jz)-apetz(i,j,jz-1))
               apztx(i,j,jz+1)=(1-kp)*piano14(i,j) + &
                  kp*(2.*apztx(i,j,jz)-apztx(i,j,jz-1))
               apzty(i,j,jz+1)=(1-kp)*piano15(i,j) + &
                  kp*(2.*apzty(i,j,jz)-apzty(i,j,jz-1))
               apztz(i,j,jz+1)=(1-kp)*piano16(i,j) + &
                  kp*(2.*apztz(i,j,jz)-apztz(i,j,jz-1))
            enddo
         enddo

      endif 


      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(6020+myid,*)uco(i,j,k),vco(i,j,k),wco(i,j,k)
                  write(6025+myid,*)uuco(i,j,k),vvco(i,j,k),wwco(i,j,k)
                  write(6030+myid,*)rhofl(i,j,k)
                  write(6035+myid,*)apcsx(i,j,k),apcsy(i,j,k),apcsz(i,j,k)
                  write(6040+myid,*)apetx(i,j,k),apety(i,j,k),apetz(i,j,k)
                  write(6045+myid,*)apztx(i,j,k),apzty(i,j,k),apztz(i,j,k)
               end do
            end do
         end do
      end if
      !
      ! bar filtering for the scale similar part
      !
      call filterb_csieta(uco,m21,kparasta,kparaend,myid,nproc)
      call filterb_csieta(vco,m22,kparasta,kparaend,myid,nproc)
      call filterb_csieta(wco,m23,kparasta,kparaend,myid,nproc)
      !
      call filterb_csieta(rhofl,m33,kparasta,kparaend,myid,nproc)
      !
      call filterb_csieta(uuco,l11,kparasta,kparaend,myid,nproc)
      call filterb_csieta(vvco,l22,kparasta,kparaend,myid,nproc)
      call filterb_csieta(wwco,l33,kparasta,kparaend,myid,nproc)


      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(6050+myid,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                  write(6055+myid,*)m33(i,j,k)
                  write(6060+myid,*)l11(i,j,k),l22(i,j,k),l33(i,j,k)
               end do
            end do
         end do
      end if
      !
      ! now I need to exchange the "ghost" plane to impose periodicity
      ! and to complete the filtering in zita (k+1 and k-1)
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
      call buffer1g(m33,4,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l11,5,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l22,6,kparasta,myid,nproc,kparasta,kparaend)
      call buffer1g(l33,7,kparasta,myid,nproc,kparasta,kparaend)
      !
      ! if I use pe and not pem I have implicitly the periodicity on k
      ! and the values in k=0 and k=jz+1 are defined in filter06
      !
      if (kp.eq.0) then

         call MPI_SENDRECV(sbuff1(1),7*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,tagls, &
            rbuff1(1),7*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrr, &
            MPI_COMM_WORLD,status,ierr)

      else if (kp.eq.1) then

         if(leftpem /= MPI_PROC_NULL) then
            call MPI_SSEND(sbuff1(1),7*(jx+2)*(jy+2), &
               MPI_REAL_SD,leftpem,tagls, &
               MPI_COMM_WORLD,ierr)
         endif
         if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),7*(jx+2)*(jy+2), &
               MPI_REAL_SD,rightpem,tagrr, &
               MPI_COMM_WORLD,status,ierr)
         endif

         if(leftpem /= MPI_PROC_NULL) then
         !      call MPI_WAIT(req1,istatus,ierr)
         endif
         if(rightpem /= MPI_PROC_NULL) then
         !      call MPI_WAIT(req2,istatus,ierr)
         endif


      endif
      !
      call buffer2(rbuff1,m21,1,kparaend+1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,m22,2,kparaend+1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,m23,3,kparaend+1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,m33,4,kparaend+1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,l11,5,kparaend+1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,l22,6,kparaend+1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,l33,7,kparaend+1,myid,nproc, &
         kparasta,kparaend)
      !
      ! now kparaend
      !
      do m=1,40*(jx+2)*(jy+2)
         sbuff(m)=0.
         rbuff(m)=0.
      enddo
      !
      call buffer1g(m21,1,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m22,2,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m23,3,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(m33,4,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l11,5,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l22,6,kparaend,myid,nproc,kparasta,kparaend)
      call buffer1g(l33,7,kparaend,myid,nproc,kparasta,kparaend)

      if (kp.eq.0) then

         call MPI_SENDRECV(sbuff1(1),7*(jx+2)*(jy+2), &
            MPI_REAL_SD,rightpe,tagrs, &
            rbuff1(1),7*(jx+2)*(jy+2), &
            MPI_REAL_SD,leftpe,taglr, &
            MPI_COMM_WORLD,status,ierr)

      else if (kp.eq.1) then

         if(rightpem /= MPI_PROC_NULL) then
            call MPI_SSEND(sbuff1(1),7*(jx+2)*(jy+2), &
               MPI_REAL_SD,rightpem,tagrs, &
               MPI_COMM_WORLD,ierr)
         endif
         if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),7*(jx+2)*(jy+2), &
               MPI_REAL_SD,leftpem,taglr, &
               MPI_COMM_WORLD,status,ierr)
         endif

         if(rightpem /= MPI_PROC_NULL) then
         !      call MPI_WAIT(req1,istatus,ierr)
         endif
         if(leftpem /= MPI_PROC_NULL) then
         !      call MPI_WAIT(req2,istatus,ierr)
         endif


      endif

      call buffer2(rbuff1,m21,1,kparasta-1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,m22,2,kparasta-1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,m23,3,kparasta-1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,m33,4,kparasta-1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,l11,5,kparasta-1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,l22,6,kparasta-1,myid,nproc, &
         kparasta,kparaend)
      call buffer2(rbuff1,l33,7,kparasta-1,myid,nproc, &
         kparasta,kparaend)


      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(6065+myid,*)m21(i,j,k),m22(i,j,k),m23(i,j,k)
                  write(6070+myid,*)m33(i,j,k)
                  write(6075+myid,*)l11(i,j,k),l22(i,j,k),l33(i,j,k)
               end do
            end do
         end do
      end if
      !
      ! complete filtering in zita
      !
      call filterb_zita(m21,ucof,kparasta,kparaend,myid,nproc)
      call filterb_zita(m22,vcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(m23,wcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(m33,rhof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l11,uucof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l22,vvcof,kparasta,kparaend,myid,nproc)
      call filterb_zita(l33,wwcof,kparasta,kparaend,myid,nproc)



      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(6080+myid,*)ucof(i,j,k),vcof(i,j,k),wcof(i,j,k)
                  write(6085+myid,*)rhof(i,j,k)
                  write(6090+myid,*)uucof(i,j,k),vvcof(i,j,k),wwcof(i,j,k)
               end do
            end do
         end do
      end if
      !
      ! compute Leonard term Arhoik
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
               ass11(i,j,k)=uucof(i,j,k)-rhof(i,j,k)*ucof(i,j,k)   !rhoU
               ass22(i,j,k)=vvcof(i,j,k)-rhof(i,j,k)*vcof(i,j,k)   !rhoV
               ass33(i,j,k)=wwcof(i,j,k)-rhof(i,j,k)*wcof(i,j,k)   !rhoW
            enddo
         enddo
      enddo
      !
      ! impose value i=jx equal to jx-1
      !hicco and i=0 uguale a i=1 ????
      do ii=1,ip

         do k=kparasta,kparaend
            do j=1,jy
               ass11(jx,j,k)=ass11(jx-1,j,k)
               ass22(jx,j,k)=ass22(jx-1,j,k)
               ass33(jx,j,k)=ass33(jx-1,j,k)

            enddo
         enddo

      enddo

      ! impose value j=jy equal to jy-1
      !hicco and j=0 equal to j=1 ????
      do jj=1,jp

         do k=kparasta,kparaend
            do i=1,jx
               ass11(i,jy,k)=ass11(i,jy-1,k)
               ass22(i,jy,k)=ass22(i,jy-1,k)
               ass33(i,jy,k)=ass33(i,jy-1,k)

            enddo
         enddo

      enddo

      !
      ! impose value k=jz equal to jz-1
      do kk=1,kp

         if(myid.eq.nproc-1)then
            do j=1,jy
               do i=1,jx
                  ass11(i,j,jz)=ass11(i,j,jz-1)
                  ass22(i,j,jz)=ass22(i,j,jz-1)
                  ass33(i,j,jz)=ass33(i,j,jz-1)
               enddo
            enddo
         endif

      enddo
      !
      ! periodicity for term Arhoik
      !
      call periodic(ass11,ass22,ass33,myid,nproc, &
         kparasta,kparaend)

      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(6100+myid,*)ass11(i,j,k),ass22(i,j,k),ass33(i,j,k)
               end do
            end do
         end do
      end if
      !
      ! periodicity on zita
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
            call buffer1g(ass22,2,jz,myid,nproc,kparasta,kparaend)
            call buffer1g(ass33,3,jz,myid,nproc,kparasta,kparaend)

         else if (myid.eq.0) then

            call buffer1g(ass11,1,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass22,2,1,myid,nproc,kparasta,kparaend)
            call buffer1g(ass33,3,1,myid,nproc,kparasta,kparaend)
         endif
         !
         ! now the exchange so taht P0 knows k=jz and Pn-1 knows k=1
         !
         if (myid.eq.nproc-1) then

            call MPI_SENDRECV(sbuff1(1),3*(jx+2)*(jy+2),MPI_REAL_SD,0,901, &
               rbuff1(1),3*(jx+2)*(jy+2),MPI_REAL_SD,0,801, &
               MPI_COMM_WORLD,status,ierr)

         else if (myid.eq.0) then

            call MPI_SENDRECV( &
               sbuff1(1),3*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,801, &
               rbuff1(1),3*(jx+2)*(jy+2),MPI_REAL_SD,nproc-1,901, &
               MPI_COMM_WORLD,status,ierr)

         endif

         if (myid.eq.0) then


            call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
            call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
            call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)

         else if (myid.eq.nproc-1) then

            call buffer2g(rbuff1,piano1,1,myid,nproc,kparasta,kparaend)
            call buffer2g(rbuff1,piano2,2,myid,nproc,kparasta,kparaend)
            call buffer2g(rbuff1,piano3,3,myid,nproc,kparasta,kparaend)

         endif
         !
         ! now P0 knows plane k=jz and Pn-1 knows plane k=1
         !
         if(myid.eq.0)then
            do j=1,jy
               do i=1,jx

                  ass11(i,j,0)=piano1(i,j)
                  ass22(i,j,0)=piano2(i,j)
                  ass33(i,j,0)=piano3(i,j)
               enddo
            enddo
         endif
         !
         if(myid.eq.nproc-1)then
            do j=1,jy
               do i=1,jx

                  ass11(i,j,jz+1)=piano1(i,j)
                  ass22(i,j,jz+1)=piano2(i,j)
                  ass33(i,j,jz+1)=piano3(i,j)

               enddo
            enddo
         endif
      !
      endif
      !
      ! finally computation for cgra123
      ! I need to comunicate the plane kparaend+1 of ass33 to compute cgra3
      !

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

      if(leftpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req1,istatus,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req2,istatus,ierr)
      endif

      !********************************************************
      ! calcola sgs stress parte scale similar
      !v      do k=kparasta,kparaend
      !v    do j=1,n2
      !v    do i=1,n1
      !v    pc1(i,j,k)  = sgs31(i,j,k)
      !v    pc2(i,j,k)  = sgs32(i,j,k)
      !v    pc3(i,j,k)  = sgs33(i,j,k)
      !v    ap11(i,j,k) = apcsx(i,j,k)
      !v    ap12(i,j,k) = apcsy(i,j,k)
      !v    ap13(i,j,k) = apcsz(i,j,k)
      !v    ap21(i,j,k) = apetx(i,j,k)
      !v    ap22(i,j,k) = apety(i,j,k)
      !v    ap23(i,j,k) = apetz(i,j,k)
      !v    ap31(i,j,k) = apztx(i,j,k)
      !v    ap32(i,j,k) = apzty(i,j,k)
      !v    ap33(i,j,k) = apztz(i,j,k)
      !v    end do
      !v    end do
      !v    end do

      !v     call inverse_para2(fil31,fil32,fil33,sgs31,sgs32,sgs33,
      !v    > apcsx,apcsy,apcsz,apetx,apety,apetz,apztx,apzty,apztz,
      !v    > kparasta,kparaend)
      !v     call inverse_para2(myid,nproc,kparasta,kparaend)
      !v
      !v    do k=kparasta,kparaend
      !v    do j=1,n2
      !v    do i=1,n1
      !v    fil31(i,j,k) = pp1(i,j,k)
      !v    fil32(i,j,k) = pp2(i,j,k)
      !v    fil33(i,j,k) = pp3(i,j,k)
      !v    end do
      !v    end do
      !v    end do

      !vmedia sui piani di omogeneita'
      !

      !vcc      if(ktime.eq.i_print*(ktime/i_print)) then
      !vcc      if(ti.gt.tp)then

      !v    do j=1,jy

      !v    sus_loc11(j)=0.
      !v    sus_loc12(j)=0.
      !v    sus_loc13(j)=0.
      !v    sus_loc22(j)=0.
      !v    sus_loc23(j)=0.
      !v    sus_loc33(j)=0.

      !v    do k=kparasta,kparaend
      !v    do i=1,jx
      !v    sus_loc11(j)=sus_loc11(j)+fil11(i,j,k)
      !v    sus_loc22(j)=sus_loc22(j)+fil22(i,j,k)
      !v    sus_loc33(j)=sus_loc33(j)+fil33(i,j,k)
      !v    sus_loc12(j)=sus_loc12(j)+.5*(fil12(i,j,k)+fil21(i,j,k))
      !v    sus_loc13(j)=sus_loc13(j)+.5*(fil13(i,j,k)+fil31(i,j,k))
      !v    sus_loc23(j)=sus_loc23(j)+.5*(fil23(i,j,k)+fil32(i,j,k))
      !v    enddo
      !v    enddo

      !v    enddo
      !vora occorre sommare le somme locali e passare
      !vil risultato a tutti i PEs (per ciascun piano)
      !
      !v    do m=1,40*(jx+2)*(jy+2)
      !v     sbuffd(m)=0.
      !v     rbuffd(m)=0.
      !v    enddo

      !v    call buffvect1d(sbuffd,sus_loc11,1)
      !v    call buffvect1d(sbuffd,sus_loc12,2)
      !v    call buffvect1d(sbuffd,sus_loc13,3)
      !v    call buffvect1d(sbuffd,sus_loc22,4)
      !v    call buffvect1d(sbuffd,sus_loc23,5)
      !v    call buffvect1d(sbuffd,sus_loc33,6)

      !v    call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),6*jy,MPI_DOUBLE_PRECISION,
      !v   >                   MPI_SUM,MPI_COMM_WORLD,ierr)

      !v    call buffvect2d(rbuffd,sus11,1)
      !v    call buffvect2d(rbuffd,sus12,2)
      !v    call buffvect2d(rbuffd,sus13,3)
      !v    call buffvect2d(rbuffd,sus22,4)
      !v    call buffvect2d(rbuffd,sus23,5)
      !v    call buffvect2d(rbuffd,sus33,6)
      !
      !ve infine il valore medio, ovvero la somma diviso il
      !vnumero di punti sul piano j-esimo (noto a tutti i PEs)
      !
      !v    do j=1,jy
      !v    sus11(j)=cb(j)*sus11(j)/somma
      !v    sus22(j)=cb(j)*sus22(j)/somma
      !v    sus33(j)=cb(j)*sus33(j)/somma
      !v    sus12(j)=cb(j)*sus12(j)/somma
      !v    sus13(j)=cb(j)*sus13(j)/somma
      !v    sus23(j)=cb(j)*sus23(j)/somma
      !v    enddo

      !vcc      end if

      if (debugg.eq.1) then
         do k=kparasta,kparaend
            do j=1,n2
               do i=1,n1
                  write(6105+myid,*)ass11(i,j,k),ass22(i,j,k),ass33(i,j,k)
               end do
            end do
         end do
      end if

      if(myid.eq.0)then
         kparastan=2*kp
      else
         kparastan=kparasta
      endif


      do k=kparastam,kparaendm
         do j=1+jp,jy-jp
            do i=2*ip,jx-ip
               cgra1(i,j,k)=.5*(ass11(i,j,k)+ass11(i+1,j,k))
            enddo
         enddo
      enddo

      do k=kparastam,kparaendm
         do j=2*jp,jy-jp
            do i=1+ip,jx-ip
               cgra2(i,j,k)=.5*(ass22(i,j,k)+ass22(i,j+1,k))
            enddo
         enddo
      enddo

      do k=kparastan,kparaendm
         do j=1+jp,jy-jp
            do i=1+ip,jx-ip
               cgra3(i,j,k)=.5*(ass33(i,j,k)+ass33(i,j,k+1))
            enddo
         enddo
      enddo

   endif
   return
end
