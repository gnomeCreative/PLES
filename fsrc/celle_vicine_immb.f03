!***********************************************************************
subroutine celle_vicine_immb (i,j,k, &
   usn,udx,usp,ust,uav,uin, &
   vsn,vdx,vsp,vst,vav,vin, &
   wsn,wdx,wsp,wst,wav,win, &
   rsn,rdx,rsp,rst,rav,rin, &
   asn,adx,asp,ast,aav,ain,      &
   myid,nproc,ktime,kiter,l)
   !***********************************************************************
   ! find velocity for the stencil around V, needed by incrementi
   !
   !-----------------------------------------------------------------------
   use myarrays_velo3
   use myarrays_metri3
   use scala3
   use period

   implicit none

   !
   !-----------------------------------------------------------------------
   !     array declaration

   integer i,j,k  !index V i,j,k
   integer isc
   integer indice_per,l
   integer myid,nproc
   integer ktime,kiter

   real usn,udx,usp,ust,uav,uin
   real vsn,vdx,vsp,vst,vav,vin
   real wsn,wdx,wsp,wst,wav,win
   real rsn(nscal),rdx(nscal),rsp(nscal),rst(nscal),rav(nscal)
   real rin(nscal)
   real asn,adx,asp,ast,aav,ain

   !-----------------------------------------------------------------------
   !
   if(ktime.eq.1.and.myid.eq.0.and.kiter.eq.1.and.l.eq.1)then
      write(*,*)'CHECK SUBROUTINE CELLE VICINE'
   end if

   !-----------------------------------------------------------------------
   ! inside the flow field
   !
   ! periodicity correction at the end

   ! direction i left, right
   usn = .5*(u(i,j,k)+u(i-1,j,k))
   vsn = .5*(v(i,j,k)+v(i-1,j,k))
   wsn = .5*(w(i,j,k)+w(i-1,j,k))
   asn = .5*(annit(i,j,k)+annit(i-1,j,k))

   udx = .5*(u(i+1,j,k)+u(i,j,k))
   vdx = .5*(v(i+1,j,k)+v(i,j,k))
   wdx = .5*(w(i+1,j,k)+w(i,j,k))
   adx = .5*(annit(i+1,j,k)+annit(i,j,k))

   do isc=1,nscal
      rsn(isc) = .5*(rhov(isc,i  ,j,k)+rhov(isc,i-1,j,k))
      rdx(isc) = .5*(rhov(isc,i+1,j,k)+rhov(isc,i  ,j,k))
   end do

   ! direction j bottom, upper
   ust = .5*(u(i,j,k)+u(i,j-1,k))
   vst = .5*(v(i,j,k)+v(i,j-1,k))
   wst = .5*(w(i,j,k)+w(i,j-1,k))
   ast = .5*(annit(i,j,k)+annit(i,j-1,k))

   usp = .5*(u(i,j+1,k)+u(i,j,k))
   vsp = .5*(v(i,j+1,k)+v(i,j,k))
   wsp = .5*(w(i,j+1,k)+w(i,j,k))
   asp = .5*(annit(i,j+1,k)+annit(i,j,k))

   do isc=1,nscal
      rst(isc) = .5*(rhov(isc,i,j  ,k)+rhov(isc,i,j-1,k))
      rsp(isc) = .5*(rhov(isc,i,j+1,k)+rhov(isc,i,j  ,k))
   end do

   ! direction k back, front
   !hicco mettere in modo consistente col resto del codice !!!!!
   uin = .5*(u(i,j,k)+u(i,j,k-1))
   vin = .5*(v(i,j,k)+v(i,j,k-1))
   win = .5*(w(i,j,k)+w(i,j,k-1))
   ain = .5*(annit(i,j,k)+annit(i,j,k-1))

   uav = .5*(u(i,j,k+1)+u(i,j,k))
   vav = .5*(v(i,j,k+1)+v(i,j,k))
   wav = .5*(w(i,j,k+1)+w(i,j,k))
   aav = .5*(annit(i,j,k+1)+annit(i,j,k))

   do isc=1,nscal
      rin(isc) = .5*(rhov(isc,i,j,k  )+rhov(isc,i,j,k-1))
      rav(isc) = .5*(rhov(isc,i,j,k+1)+rhov(isc,i,j,k  ))
   end do
   !
   !-----------------------------------------------------------------------
   ! correction for periodicity
   !
   ! direction i
   !
   !     periodic
   do indice_per=1,1-ip
      if(i .eq. 1 )then
         usn = .5*(u(jx,j,k)+u(i,j,k))
         vsn = .5*(v(jx,j,k)+v(i,j,k))
         wsn = .5*(w(jx,j,k)+w(i,j,k))
         do isc=1,nscal
            rsn(isc) = .5*(rhov(isc,jx,j,k)+rhov(isc,i,j,k))
         end do
         asn = .5*(annit(jx,j,k)+annit(i,j,k))
      end if
      if(i .eq. jx)then
         udx = .5*(u(1,j,k)+u(i,j,k))
         vdx = .5*(v(1,j,k)+v(i,j,k))
         wdx = .5*(w(1,j,k)+w(i,j,k))
         do isc=1,nscal
            rdx(isc) = .5*(rhov(isc,1,j,k)+rhov(isc,i,j,k))
         end do
         adx = .5*(annit(1,j,k)+annit(i,j,k))
      end if
   end do

   !     not periodic
   do indice_per=1,ip
      if(i .eq. 1 )then
         usn = u(0,j,k)
         vsn = v(0,j,k)
         wsn = w(0,j,k)
         do isc=1,nscal
            rsn(isc) = rhov(isc,0,j,k)
         end do
         asn = annit(0,j,k)
      end if
      if(i .eq. jx)then
         udx = u(jx+1,j,k)
         vdx = v(jx+1,j,k)
         wdx = w(jx+1,j,k)
         do isc=1,nscal
            rdx(isc) = rhov(isc,jx+1,j,k)
         end do
         adx = annit(jx+1,j,k)
      end if
   end do
   !
   !-----------------------------------------------------------------------
   !
   ! direction j
   !
   !     periodic
   do indice_per=1,1-jp
      if(j .eq. 1 )then
         ust = .5*(u(i,jy,k)+u(i,j,k))
         vst = .5*(v(i,jy,k)+v(i,j,k))
         wst = .5*(w(i,jy,k)+w(i,j,k))
         do isc=1,nscal
            rst(isc) = .5*(rhov(isc,i,jy,k)+rhov(isc,i,j,k))
         end do
         ast = .5*(annit(i,jy,k)+annit(i,j,k))
      end if
      if(j .eq. jy)then
         usp = .5*(u(i,1,k)+u(i,j,k))
         vsp = .5*(v(i,1,k)+v(i,j,k))
         wsp = .5*(w(i,1,k)+w(i,j,k))
         do isc=1,nscal
            rsp(isc) = .5*(rhov(isc,i,1,k)+rhov(isc,i,j,k))
         end do
         asp = .5*(annit(i,1,k)+annit(i,j,k))
      end if
   end do
   !
   !     not periodic
   do indice_per=1,jp
      if(j .eq. 1 )then
         ust = u(i,0,k)
         vst = v(i,0,k)
         wst = w(i,0,k)
         do isc=1,nscal
            rst(isc) = rhov(isc,i,1,k) ! zero derivative   rho(i,0,k)
         end do
         ast = annit(i,1,k)
      end if
      if(j .eq. jy)then
         usp = u(i,jy+1,k)
         vsp = v(i,jy+1,k)
         wsp = w(i,jy+1,k)
         do isc=1,nscal
            rsp(isc) = rhov(isc,i,jy,k)  ! zero derivative  rho(i,jy+1,k)
         end do
         asp = annit(i,jy,k)
      end if
   end do

   !---------------------------------------------------------------
   ! direction k
   !
   ! periodic
   !
   ! with passo_piani I know the velocity at the border
   ! probably the next part of the code is not necessary
   !
   do indice_per=1,1-kp
      if(myid.eq.0)then
         if(k .eq. 1 )then
            uin = .5*(u(i,j,0)+u(i,j,k))
            vin = .5*(v(i,j,0)+v(i,j,k))
            win = .5*(w(i,j,0)+w(i,j,k))
            do isc=1,nscal
               rin(isc) = .5*(rhov(isc,i,j,0)+rhov(isc,i,j,k))
            end do
            ain = .5*(annit(i,j,0)+annit(i,j,k))
         end if
      end if
      if(myid.eq.nproc-1)then
         if(k .eq. jz)then
            uav = .5*(u(i,j,jz+1)+u(i,j,k))
            vav = .5*(v(i,j,jz+1)+v(i,j,k))
            wav = .5*(w(i,j,jz+1)+w(i,j,k))
            do isc=1,nscal
               rav(isc) = .5*(rhov(isc,i,j,jz+1) &
                  +rhov(isc,i,j,k   ))
            end do
            aav = .5*(annit(i,j,jz+1)+annit(i,j,k))
         end if
      end if
   end do
   !
   !     not periodic
   do indice_per=1,kp
      if(myid.eq.0)then
         if(k .eq. 1 )then
            uin = u(i,j,0)
            vin = v(i,j,0)
            win = w(i,j,0)
            do isc=1,nscal
               rin(isc) = rhov(isc,i,j,0)
            end do
            ain = annit(i,j,0)
         end if
      end if
      if(myid.eq.nproc-1)then
         if(k .eq. jz)then
            uav = u(i,j,jz+1)
            vav = v(i,j,jz+1)
            wav = w(i,j,jz+1)
            do isc=1,nscal
               rav(isc) = rhov(isc,i,j,jz+1)
            end do
            aav = annit(i,j,jz+1)
         end if
      end if
   end do
   !-----------------------------------------------------------------------

   return
end
