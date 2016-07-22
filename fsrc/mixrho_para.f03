!***********************************************************************
subroutine mixrho_para(rho)
   !***********************************************************************
   ! compute scale similar part for scalar equation
   !
   use filter_module
   use mysending
   use turbo_module
   use myarrays_velo3
   use myarrays_metri3
   use mysettings, only: inmodrho
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
   integer i,j,k,kparastal,m
   integer ii,jj,kk
   integer debugg
   real aden
   !-----------------------------------------------------------------------

   debugg = 0
   !
   aden=float(n1)*float(n3)

   do k=kparasta,kparaend
      do j=1,n2
         do i=0,n1
            cgra1(i,j,k)=0.
         enddo
      enddo
   enddo
   !
   do k=kparasta,kparaend
      do j=0,n2
         do i=1,n1
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
      do j=1,n2
         do i=1,n1
            cgra3(i,j,k)=0.
         enddo
      enddo
   enddo
   !
   !-------------------------------------------------
   ! compute controvariant term for product rho*velocity
   ! only if scale similar
   !
   if (inmodrho) then

      do k=kparasta,kparaend
         do j=1,n2
            do i=1,n1

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
               uco(i,j,k)= apcsx(i,j,k)*u(i,j,k)+apcsy(i,j,k)*v(i,j,k)+apcsz(i,j,k)*w(i,j,k)
               !
               vco(i,j,k)= apetx(i,j,k)*u(i,j,k)+apety(i,j,k)*v(i,j,k)+apetz(i,j,k)*w(i,j,k)
               !
               wco(i,j,k)= apztx(i,j,k)*u(i,j,k)+apzty(i,j,k)*v(i,j,k)+apztz(i,j,k)*w(i,j,k)
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
         do j=1,n2

            uco(0,j,k)=(1-ip)*uco(n1,j,k) +ip*(2.*uco(1,j,k)-uco(2,j,k))
            vco(0,j,k)=(1-ip)*vco(n1,j,k) +ip*(2.*vco(1,j,k)-vco(2,j,k))
            wco(0,j,k)=(1-ip)*wco(n1,j,k) +ip*(2.*wco(1,j,k)-wco(2,j,k))
            !
            uuco(0,j,k)=(1-ip)*uuco(n1,j,k) +ip*(2.*uuco(1,j,k)-uuco(2,j,k))
            vvco(0,j,k)=(1-ip)*vvco(n1,j,k) +ip*(2.*vvco(1,j,k)-vvco(2,j,k))
            wwco(0,j,k)=(1-ip)*wwco(n1,j,k) +ip*(2.*wwco(1,j,k)-wwco(2,j,k))
            rhofl(0,j,k)=(1-ip)*rhofl(n1,j,k) +ip*(2.*rhofl(1,j,k)-rhofl(2,j,k))
            !
            apcsx(0,j,k)=(1-ip)*apcsx(n1,j,k) +ip*(2.*apcsx(1,j,k)-apcsx(2,j,k))
            apcsy(0,j,k)=(1-ip)*apcsy(n1,j,k) +ip*(2.*apcsy(1,j,k)-apcsy(2,j,k))
            apcsz(0,j,k)=(1-ip)*apcsz(n1,j,k) +ip*(2.*apcsz(1,j,k)-apcsz(2,j,k))
            apetx(0,j,k)=(1-ip)*apetx(n1,j,k) +ip*(2.*apetx(1,j,k)-apetx(2,j,k))
            apety(0,j,k)=(1-ip)*apety(n1,j,k) +ip*(2.*apety(1,j,k)-apety(2,j,k))
            apetz(0,j,k)=(1-ip)*apetz(n1,j,k) +ip*(2.*apetz(1,j,k)-apetz(2,j,k))
            apztx(0,j,k)=(1-ip)*apztx(n1,j,k) +ip*(2.*apztx(1,j,k)-apztx(2,j,k))
            apzty(0,j,k)=(1-ip)*apzty(n1,j,k) +ip*(2.*apzty(1,j,k)-apzty(2,j,k))
            apztz(0,j,k)=(1-ip)*apztz(n1,j,k) +ip*(2.*apztz(1,j,k)-apztz(2,j,k))
            !
            !
            apcsx(n1+1,j,k)=(1-ip)*apcsx(1,j,k) +ip*(2.*apcsx(n1,j,k)-apcsx(n1-1,j,k))
            apcsy(n1+1,j,k)=(1-ip)*apcsy(1,j,k) +ip*(2.*apcsy(n1,j,k)-apcsy(n1-1,j,k))
            apcsz(n1+1,j,k)=(1-ip)*apcsz(1,j,k) +ip*(2.*apcsz(n1,j,k)-apcsz(n1-1,j,k))
            apetx(n1+1,j,k)=(1-ip)*apetx(1,j,k) +ip*(2.*apetx(n1,j,k)-apetx(n1-1,j,k))
            apety(n1+1,j,k)=(1-ip)*apety(1,j,k) +ip*(2.*apety(n1,j,k)-apety(n1-1,j,k))
            apetz(n1+1,j,k)=(1-ip)*apetz(1,j,k) +ip*(2.*apetz(n1,j,k)-apetz(n1-1,j,k))
            apztx(n1+1,j,k)=(1-ip)*apztx(1,j,k) +ip*(2.*apztx(n1,j,k)-apztx(n1-1,j,k))
            apzty(n1+1,j,k)=(1-ip)*apzty(1,j,k) +ip*(2.*apzty(n1,j,k)-apzty(n1-1,j,k))
            apztz(n1+1,j,k)=(1-ip)*apztz(1,j,k) +ip*(2.*apztz(n1,j,k)-apztz(n1-1,j,k))
            !
            uco(n1+1,j,k)=(1-ip)*uco(1,j,k) +ip*(2.*uco(n1,j,k)-uco(n1-1,j,k))
            vco(n1+1,j,k)=(1-ip)*vco(1,j,k) +ip*(2.*vco(n1,j,k)-vco(n1-1,j,k))
            wco(n1+1,j,k)=(1-ip)*wco(1,j,k) +ip*(2.*wco(n1,j,k)-wco(n1-1,j,k))
            !
            uuco(n1+1,j,k)=(1-ip)*uuco(1,j,k) +ip*(2.*uuco(n1,j,k)-uuco(n1-1,j,k))
            vvco(n1+1,j,k)=(1-ip)*vvco(1,j,k) +ip*(2.*vvco(n1,j,k)-vvco(n1-1,j,k))
            wwco(n1+1,j,k)=(1-ip)*wwco(1,j,k) +ip*(2.*wwco(n1,j,k)-wwco(n1-1,j,k))
            rhofl(n1+1,j,k)=(1-ip)*rhofl(1,j,k) +ip*(2.*rhofl(n1,j,k)-rhofl(n1-1,j,k))

         enddo
      enddo
      !
      ! extrapolation for velocity on side j=0 j=jy+1
      !
      do k=kparasta,kparaend
         do i=1,n1

            uco(i,0,k)=2.*uco(i,1,k)-uco(i,2,k)
            vco(i,0,k)=2.*vco(i,1,k)-vco(i,2,k)
            wco(i,0,k)=2.*wco(i,1,k)-wco(i,2,k)
            !
            uuco(i,0,k)=2.*uuco(i,1,k)-uuco(i,2,k)
            vvco(i,0,k)=2.*vvco(i,1,k)-vvco(i,2,k)
            wwco(i,0,k)=2.*wwco(i,1,k)-wwco(i,2,k)
            rhofl(i,0,k)=2.*rhofl(i,1,k)-rhofl(i,2,k)
            !
            apcsx(i,0,k)=2.*apcsx(i,1,k)-apcsx(i,2,k)
            apcsy(i,0,k)=2.*apcsy(i,1,k)-apcsy(i,2,k)
            apcsz(i,0,k)=2.*apcsz(i,1,k)-apcsz(i,2,k)
            apetx(i,0,k)=2.*apetx(i,1,k)-apetx(i,2,k)
            apety(i,0,k)=2.*apety(i,1,k)-apety(i,2,k)
            apetz(i,0,k)=2.*apetz(i,1,k)-apetz(i,2,k)
            apztx(i,0,k)=2.*apztx(i,1,k)-apztx(i,2,k)
            apzty(i,0,k)=2.*apzty(i,1,k)-apzty(i,2,k)
            apztz(i,0,k)=2.*apztz(i,1,k)-apztz(i,2,k)
            !
            apcsx(i,n2+1,k)=2.*apcsx(i,n2,k)-apcsx(i,n2-1,k)
            apcsy(i,n2+1,k)=2.*apcsy(i,n2,k)-apcsy(i,n2-1,k)
            apcsz(i,n2+1,k)=2.*apcsz(i,n2,k)-apcsz(i,n2-1,k)
            apetx(i,n2+1,k)=2.*apetz(i,n2,k)-apetx(i,n2-1,k)
            apety(i,n2+1,k)=2.*apety(i,n2,k)-apety(i,n2-1,k)
            apetz(i,n2+1,k)=2.*apetz(i,n2,k)-apetz(i,n2-1,k)
            apztx(i,n2+1,k)=2.*apztx(i,n2,k)-apztx(i,n2-1,k)
            apzty(i,n2+1,k)=2.*apzty(i,n2,k)-apzty(i,n2-1,k)
            apztz(i,n2+1,k)=2.*apztz(i,n2,k)-apztz(i,n2-1,k)
            !
            uco(i,n2+1,k)=2.*uco(i,n2,k)-uco(i,n2-1,k)
            vco(i,n2+1,k)=2.*vco(i,n2,k)-vco(i,n2-1,k)
            wco(i,n2+1,k)=2.*wco(i,n2,k)-wco(i,n2-1,k)
            !
            uuco(i,n2+1,k)=2.*uuco(i,n2,k)-uuco(i,n2-1,k)
            vvco(i,n2+1,k)=2.*vvco(i,n2,k)-vvco(i,n2-1,k)
            wwco(i,n2+1,k)=2.*wwco(i,n2,k)-wwco(i,n2-1,k)
            rhofl(i,n2+1,k)=2.*rhofl(i,n2,k)-rhofl(i,n2-1,k)

         enddo
      enddo
      !
      ! periodicity on sides front and back
      !
      ! extrapolation for velocity on sides k=0 and k=jz+1
      !
      ! I use a sending buffer to contain all the plane jx*jy
      ! (16) for the exchange between P0 and Pn-1

      do m=1,40*(n1+2)*(n2+2)
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

         call buffer1g(uco,1,n3)
         call buffer1g(vco,2,n3)
         call buffer1g(wco,3,n3)
         call buffer1g(uuco,4,n3)
         call buffer1g(vvco,5,n3)
         call buffer1g(wwco,6,n3)
         call buffer1g(rhofl,7,n3)
         call buffer1g(apcsx,8,n3)
         call buffer1g(apcsy,9,n3)
         call buffer1g(apcsz,10,n3)
         call buffer1g(apetx,11,n3)
         call buffer1g(apety,12,n3)
         call buffer1g(apetz,13,n3)
         call buffer1g(apztx,14,n3)
         call buffer1g(apzty,15,n3)
         call buffer1g(apztz,16,n3)

      else if (myid.eq.0) then

         call buffer1g(uco,1,1)
         call buffer1g(vco,2,1)
         call buffer1g(wco,3,1)
         call buffer1g(uuco,4,1)
         call buffer1g(vvco,5,1)
         call buffer1g(wwco,6,1)
         call buffer1g(rhofl,7,1)
         call buffer1g(apcsx,8,1)
         call buffer1g(apcsy,9,1)
         call buffer1g(apcsz,10,1)
         call buffer1g(apetx,11,1)
         call buffer1g(apety,12,1)
         call buffer1g(apetz,13,1)
         call buffer1g(apztx,14,1)
         call buffer1g(apzty,15,1)
         call buffer1g(apztz,16,1)

      endif
      !
      ! now exchange so that P0 knows k=jz and Pn-1 knows k=1
      !
      if (myid.eq.nproc-1) then

         call MPI_SENDRECV(sbuff1(1),16*(n1+2)*(n2+2),MPI_REAL_SD,0,9101, &
            rbuff1(1),16*(n1+2)*(n2+2),MPI_REAL_SD,0,8101,MPI_COMM_WORLD,status,ierr)

      else if (myid.eq.0) then

         call MPI_SENDRECV(sbuff1(1),16*(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,8101, &
            rbuff1(1),16*(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,9101,MPI_COMM_WORLD,status,ierr)
      
      endif

      if (myid.eq.0) then

         call buffer2gg(piano1,1)
         call buffer2gg(piano2,2)
         call buffer2gg(piano3,3)
         call buffer2gg(piano4,4)
         call buffer2gg(piano5,5)
         call buffer2gg(piano6,6)
         call buffer2gg(piano7,7)
         call buffer2gg(piano8,8)
         call buffer2gg(piano9,9)
         call buffer2gg(piano10,10)
         call buffer2gg(piano11,11)
         call buffer2gg(piano12,12)
         call buffer2gg(piano13,13)
         call buffer2gg(piano14,14)
         call buffer2gg(piano15,15)
         call buffer2gg(piano16,16)

      else if (myid.eq.nproc-1) then

         call buffer2gg(piano1,1)
         call buffer2gg(piano2,2)
         call buffer2gg(piano3,3)
         call buffer2gg(piano4,4)
         call buffer2gg(piano5,5)
         call buffer2gg(piano6,6)
         call buffer2gg(piano7,7)
         call buffer2gg(piano8,8)
         call buffer2gg(piano9,9)
         call buffer2gg(piano10,10)
         call buffer2gg(piano11,11)
         call buffer2gg(piano12,12)
         call buffer2gg(piano13,13)
         call buffer2gg(piano14,14)
         call buffer2gg(piano15,15)
         call buffer2gg(piano16,16)

      endif

      ! now P0 knows plane k=jz
      ! and Pn-1 knows plane k=1

      if(myid.eq.0)then

         do j=1,n2
            do i=1,n1
               !
               uco(i,j,0)=(1-kp)*piano1(i,j) +kp*(2.*uco(i,j,1)-uco(i,j,2))
               vco(i,j,0)=(1-kp)*piano2(i,j) +kp*(2.*vco(i,j,1)-vco(i,j,2))
               wco(i,j,0)=(1-kp)*piano3(i,j) +kp*(2.*wco(i,j,1)-wco(i,j,2))
               !
               uuco(i,j,0)=(1-kp)*piano4(i,j) +kp*(2.*uuco(i,j,1)-uuco(i,j,2))
               vvco(i,j,0)=(1-kp)*piano5(i,j) +kp*(2.*vvco(i,j,1)-vvco(i,j,2))
               wwco(i,j,0)=(1-kp)*piano6(i,j) +kp*(2.*wwco(i,j,1)-wwco(i,j,2))
               rhofl(i,j,0)=(1-kp)*piano7(i,j) +kp*(2.*rhofl(i,j,1)-rhofl(i,j,2))
               !
               apcsx(i,j,0)=(1-kp)*piano8(i,j) +kp*(2.*apcsx(i,j,1)-apcsx(i,j,2))
               apcsy(i,j,0)=(1-kp)*piano9(i,j) +kp*(2.*apcsy(i,j,1)-apcsy(i,j,2))
               apcsz(i,j,0)=(1-kp)*piano19(i,j) +kp*(2.*apcsz(i,j,1)-apcsz(i,j,2))
               apetx(i,j,0)=(1-kp)*piano11(i,j) +kp*(2.*apetx(i,j,1)-apetx(i,j,2))
               apety(i,j,0)=(1-kp)*piano12(i,j) +kp*(2.*apety(i,j,1)-apety(i,j,2))
               apetz(i,j,0)=(1-kp)*piano13(i,j) +kp*(2.*apetz(i,j,1)-apetz(i,j,2))
               apztx(i,j,0)=(1-kp)*piano14(i,j) +kp*(2.*apztx(i,j,1)-apztx(i,j,2))
               apzty(i,j,0)=(1-kp)*piano15(i,j) +kp*(2.*apzty(i,j,1)-apzty(i,j,2))
               apztz(i,j,0)=(1-kp)*piano16(i,j) +kp*(2.*apztz(i,j,1)-apztz(i,j,2))
            !
            end do
         end do

      endif
      !
      !
      if(myid.eq.nproc-1)then

         do j=1,n2
            do i=1,n1
               !
               uco(i,j,n3+1)=(1-kp)*piano1(i,j) +kp*(2.*uco(i,j,n3)-uco(i,j,n3-1))
               vco(i,j,n3+1)=(1-kp)*piano2(i,j) +kp*(2.*vco(i,j,n3)-vco(i,j,n3-1))
               wco(i,j,n3+1)=(1-kp)*piano3(i,j) +kp*(2.*wco(i,j,n3)-wco(i,j,n3-1))
               !
               uuco(i,j,n3+1)=(1-kp)*piano4(i,j) +kp*(2.*uuco(i,j,n3)-uuco(i,j,n3-1))
               vvco(i,j,n3+1)=(1-kp)*piano5(i,j) +kp*(2.*vvco(i,j,n3)-vvco(i,j,n3-1))
               wwco(i,j,n3+1)=(1-kp)*piano6(i,j) +kp*(2.*wwco(i,j,n3)-wwco(i,j,n3-1))
               rhofl(i,j,n3+1)=(1-kp)*piano7(i,j) +kp*(2.*rhofl(i,j,n3)-rhofl(i,j,n3-1))

               apcsx(i,j,n3+1)=(1-kp)*piano8(i,j) +kp*(2.*apcsx(i,j,n3)-apcsx(i,j,n3-1))
               apcsy(i,j,n3+1)=(1-kp)*piano9(i,j) +kp*(2.*apcsy(i,j,n3)-apcsy(i,j,n3-1))
               apcsz(i,j,n3+1)=(1-kp)*piano10(i,j) +kp*(2.*apcsz(i,j,n3)-apcsz(i,j,n3-1))
               apetx(i,j,n3+1)=(1-kp)*piano11(i,j) +kp*(2.*apetx(i,j,n3)-apetx(i,j,n3-1))
               apety(i,j,n3+1)=(1-kp)*piano12(i,j) +kp*(2.*apety(i,j,n3)-apety(i,j,n3-1))
               apetz(i,j,n3+1)=(1-kp)*piano13(i,j) +kp*(2.*apetz(i,j,n3)-apetz(i,j,n3-1))
               apztx(i,j,n3+1)=(1-kp)*piano14(i,j) +kp*(2.*apztx(i,j,n3)-apztx(i,j,n3-1))
               apzty(i,j,n3+1)=(1-kp)*piano15(i,j) +kp*(2.*apzty(i,j,n3)-apzty(i,j,n3-1))
               apztz(i,j,n3+1)=(1-kp)*piano16(i,j) +kp*(2.*apztz(i,j,n3)-apztz(i,j,n3-1))
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
      call filterb_csieta(uco,m21)
      call filterb_csieta(vco,m22)
      call filterb_csieta(wco,m23)
      !
      call filterb_csieta(rhofl,m33)
      !
      call filterb_csieta(uuco,l11)
      call filterb_csieta(vvco,l22)
      call filterb_csieta(wwco,l33)


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
      do m=1,40*(n1+2)*(n2+2)
         sbuff(m)=0.
         rbuff(m)=0.
      enddo

      call buffer1g(m21,1,kparasta)
      call buffer1g(m22,2,kparasta)
      call buffer1g(m23,3,kparasta)
      call buffer1g(m33,4,kparasta)
      call buffer1g(l11,5,kparasta)
      call buffer1g(l22,6,kparasta)
      call buffer1g(l33,7,kparasta)
      !
      ! if I use pe and not pem I have implicitly the periodicity on k
      ! and the values in k=0 and k=jz+1 are defined in filter06
      !
      if (kp.eq.0) then

         call MPI_SENDRECV(sbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,leftpe,tagls, &
            rbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,rightpe,tagrr,MPI_COMM_WORLD,status,ierr)

      else if (kp.eq.1) then

         if(leftpem /= MPI_PROC_NULL) then
            call MPI_SEND(sbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
         endif
         if(rightpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
         endif

      endif
      !
      call buffer2(rbuff1,m21,1,kparaend+1)
      call buffer2(rbuff1,m22,2,kparaend+1)
      call buffer2(rbuff1,m23,3,kparaend+1)
      call buffer2(rbuff1,m33,4,kparaend+1)
      call buffer2(rbuff1,l11,5,kparaend+1)
      call buffer2(rbuff1,l22,6,kparaend+1)
      call buffer2(rbuff1,l33,7,kparaend+1)
      !
      ! now kparaend
      !
      do m=1,40*(n1+2)*(n2+2)
         sbuff(m)=0.
         rbuff(m)=0.
      enddo
      !
      call buffer1g(m21,1,kparaend)
      call buffer1g(m22,2,kparaend)
      call buffer1g(m23,3,kparaend)
      call buffer1g(m33,4,kparaend)
      call buffer1g(l11,5,kparaend)
      call buffer1g(l22,6,kparaend)
      call buffer1g(l33,7,kparaend)

      if (kp.eq.0) then

         call MPI_SENDRECV(sbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,rightpe,tagrs, &
            rbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,leftpe,taglr,MPI_COMM_WORLD,status,ierr)

      else if (kp.eq.1) then

         if(rightpem /= MPI_PROC_NULL) then
            call MPI_SEND(sbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,rightpem,tagrs,MPI_COMM_WORLD,ierr)
         endif
         if(leftpem /= MPI_PROC_NULL) then
            call MPI_RECV(rbuff1(1),7*(n1+2)*(n2+2),MPI_REAL_SD,leftpem,taglr,MPI_COMM_WORLD,status,ierr)
         endif


      endif

      call buffer2(rbuff1,m21,1,kparasta-1)
      call buffer2(rbuff1,m22,2,kparasta-1)
      call buffer2(rbuff1,m23,3,kparasta-1)
      call buffer2(rbuff1,m33,4,kparasta-1)
      call buffer2(rbuff1,l11,5,kparasta-1)
      call buffer2(rbuff1,l22,6,kparasta-1)
      call buffer2(rbuff1,l33,7,kparasta-1)


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
      call filterb_zita(m21,ucof)
      call filterb_zita(m22,vcof)
      call filterb_zita(m23,wcof)
      call filterb_zita(m33,rhof)
      call filterb_zita(l11,uucof)
      call filterb_zita(l22,vvcof)
      call filterb_zita(l33,wwcof)



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
         do j=2,n2-1
            do i=1+ip,n1-ip
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
            do j=1,n2
               ass11(n1,j,k)=ass11(n1-1,j,k)
               ass22(n1,j,k)=ass22(n1-1,j,k)
               ass33(n1,j,k)=ass33(n1-1,j,k)

            enddo
         enddo

      enddo

      ! impose value j=jy equal to jy-1
      !hicco and j=0 equal to j=1 ????
      ! direction j is always not periodic

      do k=kparasta,kparaend
          do i=1,n1
              ass11(i,n2,k)=ass11(i,n2-1,k)
              ass22(i,n2,k)=ass22(i,n2-1,k)
              ass33(i,n2,k)=ass33(i,n2-1,k)

          enddo
      enddo

      !
      ! impose value k=jz equal to jz-1
      do kk=1,kp

         if(myid.eq.nproc-1)then
            do j=1,n2
               do i=1,n1
                  ass11(i,j,n3)=ass11(i,j,n3-1)
                  ass22(i,j,n3)=ass22(i,j,n3-1)
                  ass33(i,j,n3)=ass33(i,j,n3-1)
               enddo
            enddo
         endif

      enddo
      !
      ! periodicity for term Arhoik
      !
      call periodic(ass11,ass22,ass33)

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
         do m=1,40*(n1+2)*(n2+2)
            sbuff(m)=0.
            rbuff(m)=0.
         enddo

         if (myid.eq.nproc-1) then

            call buffer1g(ass11,1,n3)
            call buffer1g(ass22,2,n3)
            call buffer1g(ass33,3,n3)

         else if (myid.eq.0) then

            call buffer1g(ass11,1,1)
            call buffer1g(ass22,2,1)
            call buffer1g(ass33,3,1)
         endif
         !
         ! now the exchange so taht P0 knows k=jz and Pn-1 knows k=1
         !
         if (myid.eq.nproc-1) then

            call MPI_SENDRECV(sbuff1(1),3*(n1+2)*(n2+2),MPI_REAL_SD,0,901, &
               rbuff1(1),3*(n1+2)*(n2+2),MPI_REAL_SD,0,801,MPI_COMM_WORLD,status,ierr)

         else if (myid.eq.0) then

            call MPI_SENDRECV(sbuff1(1),3*(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,801, &
               rbuff1(1),3*(n1+2)*(n2+2),MPI_REAL_SD,nproc-1,901,MPI_COMM_WORLD,status,ierr)

         endif

         if (myid.eq.0) then


            call buffer2g(rbuff1,piano1,1)
            call buffer2g(rbuff1,piano2,2)
            call buffer2g(rbuff1,piano3,3)

         else if (myid.eq.nproc-1) then

            call buffer2g(rbuff1,piano1,1)
            call buffer2g(rbuff1,piano2,2)
            call buffer2g(rbuff1,piano3,3)

         endif
         !
         ! now P0 knows plane k=jz and Pn-1 knows plane k=1
         !
         if(myid.eq.0)then
            do j=1,n2
               do i=1,n1

                  ass11(i,j,0)=piano1(i,j)
                  ass22(i,j,0)=piano2(i,j)
                  ass33(i,j,0)=piano3(i,j)
               enddo
            enddo
         endif
         !
         if(myid.eq.nproc-1)then
            do j=1,n2
               do i=1,n1

                  ass11(i,j,n3+1)=piano1(i,j)
                  ass22(i,j,n3+1)=piano2(i,j)
                  ass33(i,j,n3+1)=piano3(i,j)

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
         call MPI_SEND(ass33(1,1,kparasta),(n1+2)*(n2+2),MPI_REAL_SD,leftpem,tagls,MPI_COMM_WORLD,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
         call MPI_RECV(ass33(1,1,kparaend+1),(n1+2)*(n2+2),MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
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

      !v    call MPI_ALLREDUCE(sbuffd(1),rbuffd(1),6*jy,MPI_REAL_SD,
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
         do j=2,n2-1
            do i=2*ip,n1-ip
               cgra1(i,j,k)=.5*(ass11(i,j,k)+ass11(i+1,j,k))
            enddo
         enddo
      enddo

      do k=kparastam,kparaendm
         do j=2,n2-1
            do i=1+ip,n1-ip
               cgra2(i,j,k)=.5*(ass22(i,j,k)+ass22(i,j+1,k))
            enddo
         enddo
      enddo

      do k=kparastan,kparaendm
         do j=2,n2-1
            do i=1+ip,n1-ip
               cgra3(i,j,k)=.5*(ass33(i,j,k)+ass33(i,j,k+1))
            enddo
         enddo
      enddo

   endif
   return
end
