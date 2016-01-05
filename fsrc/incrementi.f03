!***********************************************************************
subroutine incrementi (i,j,k, &
   usn,udx,usp,ust,uav,uin, &
   vsn,vdx,vsp,vst,vav,vin, &
   wsn,wdx,wsp,wst,wav,win, &
   rsn,rdx,rsp,rst,rav,rin, &
   asn,adx,asp,ast,aav,ain, &
   dudx,dudy,dudz, &
   dvdx,dvdy,dvdz, &
   dwdx,dwdy,dwdz, &
   drdx,drdy,drdz, &
   dadx,dady,dadz,      &
   l,myid,kiter,ktime,visualizzo)
   !***********************************************************************

   use myarrays_metri3
   use scala3
   use period

   implicit none
   !--------------------------------------------------------------
   !     array declaration
   integer i,j,k,l,visualizzo,isc
   integer indice_per
   integer myid,kiter,ktime
   !     for metric
   real xdx,xsn,xsop,xsot,xav,xdt,x2,x1,x0
   real xcsi,xeta,xzet
   real ydx,ysn,ysop,ysot,yav,ydt,y2,y1,y0
   real ycsi,yeta,yzet
   real zdx,zsn,zsop,zsot,zav,zdt,z2,z1,z0
   real zcsi,zeta,zzet

   real apcsx,apcsy,apcsz
   real apetx,apety,apetz
   real apztx,apzty,apztz

   real giaco
   !     increments
   real usn,udx,usp,ust,uav,uin
   real vsn,vdx,vsp,vst,vav,vin
   real wsn,wdx,wsp,wst,wav,win
   real rsn(nscal),rdx(nscal),rsp(nscal),rst(nscal),rav(nscal)
   real rin(nscal)
   real asn,adx,asp,ast,aav,ain

   real dudx,dudy,dudz
   real dvdx,dvdy,dvdz
   real dwdx,dwdy,dwdz
   real drdx(nscal),drdy(nscal),drdz(nscal)
   real dadx,dady,dadz
   !
   integer ptx,pty,ptz
   !-------------------------------------------------------------

   ptx=jx/2
   pty=3
   ptz=4

   !--------------------------------------------------------------
   ! NOTE:
   ! ijk is the V point
   !--------------------------------------------------------------

   if(ktime.eq.1.and.myid.eq.0.and.kiter.eq.1 .and.l.eq.1)then
      write(*,*)'CHECK SUBROUTINE INCREMENTI'
   end if

   !--------------------------------------------------------------
   ! initialize

   dudx=0.
   dudy=0.
   dudz=0.

   dvdx=0.
   dvdy=0.
   dvdz=0.

   dwdx=0.
   dwdy=0.
   dwdz=0.
   do isc=1,nscal
      drdx(isc)=0.
      drdy(isc)=0.
      drdz(isc)=0.
   end do

   dadx=0.
   dady=0.
   dadz=0.
   !-----------------------------------------------------------------------
   ! COMPUTE METRIC
   !
   ! compute terms xcsi at the centroid
   !
   xdx=.25*(x(i  ,j  ,k)+x(i  ,j  ,k-1) &
      +x(i  ,j-1,k)+x(i  ,j-1,k-1))
   xsn=.25*(x(i-1,j  ,k)+x(i-1,j  ,k-1) &
      +x(i-1,j-1,k)+x(i-1,j-1,k-1))
   !
   ydx=.25*(y(i  ,j  ,k)+y(i  ,j  ,k-1) &
      +y(i  ,j-1,k)+y(i  ,j-1,k-1))
   ysn=.25*(y(i-1,j  ,k)+y(i-1,j  ,k-1) &
      +y(i-1,j-1,k)+y(i-1,j-1,k-1))
   !
   zdx=.25*(z(i  ,j  ,k)+z(i  ,j  ,k-1) &
      +z(i  ,j-1,k)+z(i  ,j-1,k-1))
   zsn=.25*(z(i-1,j  ,k)+z(i-1,j  ,k-1) &
      +z(i-1,j-1,k)+z(i-1,j-1,k-1))
   !
   xcsi = xdx-xsn
   ycsi = ydx-ysn
   zcsi = zdx-zsn
   !
   !
   xsop=.25*(x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
      +x(i-1,j  ,k-1)+x(i-1,j  ,k  ))
   xsot=.25*(x(i  ,j-1,k  )+x(i  ,j-1,k-1) &
      +x(i-1,j-1,k-1)+x(i-1,j-1,k  ))
   !
   ysop=.25*(y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
      +y(i-1,j  ,k-1)+y(i-1,j  ,k  ))
   ysot=.25*(y(i  ,j-1,k  )+y(i  ,j-1,k-1) &
      +y(i-1,j-1,k-1)+y(i-1,j-1,k  ))
   !
   zsop=.25*(z(i  ,j  ,k  )+z(i  ,j  ,k-1) &
      +z(i-1,j  ,k-1)+z(i-1,j  ,k  ))
   zsot=.25*(z(i  ,j-1,k  )+z(i  ,j-1,k-1) &
      +z(i-1,j-1,k-1)+z(i-1,j-1,k  ))
   !
   xeta = xsop-xsot
   yeta = ysop-ysot
   zeta = zsop-zsot
   !
   !
   xav=.25*(x(i  ,j  ,k  )+x(i  ,j-1,k  ) &
      +x(i-1,j-1,k  )+x(i-1,j  ,k  ))
   xdt=.25*(x(i  ,j  ,k-1)+x(i  ,j-1,k-1) &
      +x(i-1,j-1,k-1)+x(i-1,j  ,k-1))
   !
   yav=.25*(y(i  ,j  ,k  )+y(i  ,j-1,k  ) &
      +y(i-1,j-1,k  )+y(i-1,j  ,k  ))
   ydt=.25*(y(i  ,j  ,k-1)+y(i  ,j-1,k-1) &
      +y(i-1,j-1,k-1)+y(i-1,j  ,k-1))
   !
   zav=.25*(z(i  ,j  ,k  )+z(i  ,j-1,k  ) &
      +z(i-1,j-1,k  )+z(i-1,j  ,k  ))
   zdt=.25*(z(i  ,j  ,k-1)+z(i  ,j-1,k-1) &
      +z(i-1,j-1,k-1)+z(i-1,j  ,k-1))
   !
   xzet = xav-xdt
   yzet = yav-ydt
   zzet = zav-zdt
   !
   !...............................................................
   ! compute the term J-1*csx etc
   !
   apcsx=(yeta*zzet-yzet*zeta)
   apcsy=(xzet*zeta-xeta*zzet)
   apcsz=(xeta*yzet-xzet*yeta)
   !
   apetx=(yzet*zcsi-ycsi*zzet)
   apety=(xcsi*zzet-xzet*zcsi)
   apetz=(xzet*ycsi-xcsi*yzet)
   !
   apztx=(ycsi*zeta-yeta*zcsi)
   apzty=(xeta*zcsi-xcsi*zeta)
   apztz=(xcsi*yeta-xeta*ycsi)
   !...............................................................
   ! giacobian determinant
   !
   giaco=xcsi*(yeta*zzet-yzet*zeta)- &
      xeta*(ycsi*zzet-yzet*zcsi)+ &
      xzet*(ycsi*zeta-yeta*zcsi)
   !
   !--------------------------------------------------------------
   ! now compute the increments

   dudx=((udx-usn)*apcsx+ &
      (usp-ust)*apetx+ &
      (uav-uin)*apztx)/giaco

   dudy=((udx-usn)*apcsy+ &
      (usp-ust)*apety+ &
      (uav-uin)*apzty)/giaco

   dudz=((udx-usn)*apcsz+ &
      (usp-ust)*apetz+ &
      (uav-uin)*apztz)/giaco



   dvdx=((vdx-vsn)*apcsx+ &
      (vsp-vst)*apetx+ &
      (vav-vin)*apztx)/giaco

   dvdy=((vdx-vsn)*apcsy+ &
      (vsp-vst)*apety+ &
      (vav-vin)*apzty)/giaco

   dvdz=((vdx-vsn)*apcsz+ &
      (vsp-vst)*apetz+ &
      (vav-vin)*apztz)/giaco



   dwdx=((wdx-wsn)*apcsx+ &
      (wsp-wst)*apetx+ &
      (wav-win)*apztx)/giaco

   dwdy=((wdx-wsn)*apcsy+ &
      (wsp-wst)*apety+ &
      (wav-win)*apzty)/giaco

   dwdz=((wdx-wsn)*apcsz+ &
      (wsp-wst)*apetz+ &
      (wav-win)*apztz)/giaco


   do isc=1,nscal
      drdx(isc)=((rdx(isc)-rsn(isc))*apcsx+ &
         (rsp(isc)-rst(isc))*apetx+ &
         (rav(isc)-rin(isc))*apztx)/giaco

      drdy(isc)=((rdx(isc)-rsn(isc))*apcsy+ &
         (rsp(isc)-rst(isc))*apety+ &
         (rav(isc)-rin(isc))*apzty)/giaco

      drdz(isc)=((rdx(isc)-rsn(isc))*apcsz+ &
         (rsp(isc)-rst(isc))*apetz+ &
         (rav(isc)-rin(isc))*apztz)/giaco
   end do

   !hicco togliere il lavoro su annit non serve!
   dadx=((adx-asn)*apcsx+ &
      (asp-ast)*apetx+ &
      (aav-ain)*apztx)/giaco

   dady=((adx-asn)*apcsy+ &
      (asp-ast)*apety+ &
      (aav-ain)*apzty)/giaco

   dadz=((adx-asn)*apcsz+ &
      (asp-ast)*apetz+ &
      (aav-ain)*apztz)/giaco

   !--------------------------------------------------------------
   return
end
