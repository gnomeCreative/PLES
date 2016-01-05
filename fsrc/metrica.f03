!***********************************************************************
subroutine metrica
   !***********************************************************************
   !
   ! build metric term 3D
   ! for solution Zang, Street, Koseff
   !
   !-----------------------------------------------------------------------
   !
   !    J(-1)csx -----> csx(i,j,k) defined on plane at constant csi (i constant)
   !    J(-1)csy -----> csy(i,j,k)
   !    J(-1)csz -----> csz(i,j,k)
   !
   !    J(-1)etx -----> etx(i,j,k) defined on plane at constant eta (j constant)
   !    J(-1)ety -----> ety(i,j,k)
   !    J(-1)etz -----> etz(i,j,k)
   !
   !    J(-1)ztx -----> ztx(i,j,k) defined on plane at constant zita (k constant)
   !    J(-1)zty -----> zty(i,j,k)
   !    J(-1)ztz -----> ztz(i,j,k)
   !
   !-----------------------------------------------------------------------
   !
   ! controvariant metric tensor *J^-1
   !
   ! J(-1)*g11  : g11(i,j,k) defined on plane at constant csi
   ! J(-1)*g12  : g12(i,j,k)
   ! J(-1)*g13  : g13(i,j,k)
   !
   ! J(-1)*g21  : g21(i,j,k) defined on plane at constant eta
   ! J(-1)*g22  : g22(i,j,k)
   ! J(-1)*g23  : g23(i,j,k)
   !
   ! J(-1)*g31  : g31(i,j,k) defined on plane at constant zeta
   ! J(-1)*g32  : g32(i,j,k)
   ! J(-1)*g33  : g33(i,j,k)
   !
   ! a_giac working giacobian: is centered on face, it does not coincide with giac(i,j,k)
   !-----------------------------------------------------------------------
   !CORIOLIS oil
   !
   ! J(-1)*g_co11  : g_co11(i,k)=0 defined on plane at constant csi
   ! J(-1)*g_co12  : g_co12(i,k)=0
   ! J(-1)*g_co13  : g_co13(i,k)
   !
   ! J(-1)*g_co31  : g_co31(i,k) defined on plane at constant zeta
   ! J(-1)*g_co32  : g_co32(i,k)=0
   ! J(-1)*g_co33  : g_co33(i,k)=0
   !
   !
   !-----------------------------------------------------------------------
   !
   use myarrays_metri3
   use mysending
   !
   use scala3
   use period

   implicit none
   !
   !-----------------------------------------------------------------------
   !     variables declaration
   integer i,j,k
   real a_csx,a_csy,a_csz
   real a_etx,a_ety,a_etz
   real a_ztx,a_zty,a_ztz
   real a_giac
   real xcsi,ycsi,zcsi
   real xeta,yeta,zeta
   real xzet,yzet,zzet
   real xsop,xsot,ysop,ysot,zsop,zsot
   real xav,xdt,yav,ydt,zav,zdt
   real xsn,xdx,ysn,ydx,zsn,zdx
   real x0,x1,x2
   real y0,y1,y2
   real z0,z1,z2
   real giac_max,giac_min

   integer kinizio,kfine

   !-----------------------------------------------------------------------
   !
   ! planes at csi constant
   !
   do i=0,jx
      do j=1,jy
         do k=kparasta,kparaend !kinizio,kfine !1,jz
            !
            if (i.eq.0.and.ip.eq.1) then
               !
               x2=.25*(x(i+2,j,k)+x(i+2,j,k-1)+x(i+2,j-1,k)+x(i+2,j-1,k-1))
               x1=.25*(x(i+1,j,k)+x(i+1,j,k-1)+x(i+1,j-1,k)+x(i+1,j-1,k-1))
               x0=.25*(x(i  ,j,k)+x(i  ,j,k-1)+x(i  ,j-1,k)+x(i  ,j-1,k-1))
               !
               y2=.25*(y(i+2,j,k)+y(i+2,j,k-1)+y(i+2,j-1,k)+y(i+2,j-1,k-1))
               y1=.25*(y(i+1,j,k)+y(i+1,j,k-1)+y(i+1,j-1,k)+y(i+1,j-1,k-1))
               y0=.25*(y(i  ,j,k)+y(i  ,j,k-1)+y(i  ,j-1,k)+y(i  ,j-1,k-1))
               !
               z2=.25*(z(i+2,j,k)+z(i+2,j,k-1)+z(i+2,j-1,k)+z(i+2,j-1,k-1))
               z1=.25*(z(i+1,j,k)+z(i+1,j,k-1)+z(i+1,j-1,k)+z(i+1,j-1,k-1))
               z0=.25*(z(i  ,j,k)+z(i  ,j,k-1)+z(i  ,j-1,k)+z(i  ,j-1,k-1))
               !
               xcsi=.5 * ( -3.*x0 + 4.*x1 - x2 )
               ycsi=.5 * ( -3.*y0 + 4.*y1 - y2 )
               zcsi=.5 * ( -3.*z0 + 4.*z1 - z2 )
            !
            else if (i.eq.jx.and.ip.eq.1) then
               !
               x2=.25*(x(i-2,j,k)+x(i-2,j,k-1)+x(i-2,j-1,k)+x(i-2,j-1,k-1))
               x1=.25*(x(i-1,j,k)+x(i-1,j,k-1)+x(i-1,j-1,k)+x(i-1,j-1,k-1))
               x0=.25*(x(i  ,j,k)+x(i  ,j,k-1)+x(i  ,j-1,k)+x(i  ,j-1,k-1))
               !
               y2=.25*(y(i-2,j,k)+y(i-2,j,k-1)+y(i-2,j-1,k)+y(i-2,j-1,k-1))
               y1=.25*(y(i-1,j,k)+y(i-1,j,k-1)+y(i-1,j-1,k)+y(i-1,j-1,k-1))
               y0=.25*(y(i  ,j,k)+y(i  ,j,k-1)+y(i  ,j-1,k)+y(i  ,j-1,k-1))
               !
               z2=.25*(z(i-2,j,k)+z(i-2,j,k-1)+z(i-2,j-1,k)+z(i-2,j-1,k-1))
               z1=.25*(z(i-1,j,k)+z(i-1,j,k-1)+z(i-1,j-1,k)+z(i-1,j-1,k-1))
               z0=.25*(z(i  ,j,k)+z(i  ,j,k-1)+z(i  ,j-1,k)+z(i  ,j-1,k-1))
               !
               xcsi=.5 * ( 3.*x0 - 4.*x1 + x2 )
               ycsi=.5 * ( 3.*y0 - 4.*y1 + y2 )
               zcsi=.5 * ( 3.*z0 - 4.*z1 + z2 )
            !
            else
               !
               xdx=.25*(x(i+1,j,k)+x(i+1,j,k-1)+x(i+1,j-1,k)+x(i+1,j-1,k-1))
               xsn=.25*(x(i-1,j,k)+x(i-1,j,k-1)+x(i-1,j-1,k)+x(i-1,j-1,k-1))
               !
               ydx=.25*(y(i+1,j,k)+y(i+1,j,k-1)+y(i+1,j-1,k)+y(i+1,j-1,k-1))
               ysn=.25*(y(i-1,j,k)+y(i-1,j,k-1)+y(i-1,j-1,k)+y(i-1,j-1,k-1))
               !
               zdx=.25*(z(i+1,j,k)+z(i+1,j,k-1)+z(i+1,j-1,k)+z(i+1,j-1,k-1))
               zsn=.25*(z(i-1,j,k)+z(i-1,j,k-1)+z(i-1,j-1,k)+z(i-1,j-1,k-1))
               !
               xcsi=.5*(xdx-xsn)
               ycsi=.5*(ydx-ysn)
               zcsi=.5*(zdx-zsn)
            !
            end if
            !
            xeta=.5*(x(i,j  ,k) + x(i,j  ,k-1))-.5*(x(i,j-1,k) + x(i,j-1,k-1))
            !
            yeta=.5*(y(i,j  ,k) + y(i,j  ,k-1))-.5*(y(i,j-1,k) + y(i,j-1,k-1))
            !
            zeta=.5*(z(i,j  ,k) + z(i,j  ,k-1))-.5*(z(i,j-1,k) + z(i,j-1,k-1))
            !
            xzet=.5*(x(i,j,k  ) + x(i,j-1,k  ))-.5*(x(i,j,k-1) + x(i,j-1,k-1))
            !
            yzet=.5*(y(i,j,k  ) + y(i,j-1,k  ))-.5*(y(i,j,k-1) + y(i,j-1,k-1))
            !
            zzet=.5*(z(i,j,k  ) + z(i,j-1,k  ))-.5*(z(i,j,k-1) + z(i,j-1,k-1))
            !
            a_giac=     xcsi*(yeta*zzet-yzet*zeta)- &
               xeta*(ycsi*zzet-yzet*zcsi)+ &
               xzet*(ycsi*zeta-yeta*zcsi)
            !
            csx(i,j,k) = yeta*zzet - yzet*zeta
            csy(i,j,k) = xzet*zeta - xeta*zzet
            csz(i,j,k) = xeta*yzet - xzet*yeta
            !
            a_etx = yzet*zcsi - ycsi*zzet
            a_ety = xcsi*zzet - xzet*zcsi
            a_etz = xzet*ycsi - xcsi*yzet
            !
            a_ztx=ycsi*zeta-yeta*zcsi
            a_zty=xeta*zcsi-xcsi*zeta
            a_ztz=xcsi*yeta-xeta*ycsi
            !
            g11(i,j,k)=(csx(i,j,k)**2 &
               +csy(i,j,k)**2 &
               +csz(i,j,k)**2)/a_giac
            !
            g12(i,j,k)=(csx(i,j,k)*a_etx &
               +csy(i,j,k)*a_ety &
               +csz(i,j,k)*a_etz)/a_giac
            !
            g13(i,j,k)=(csx(i,j,k)*a_ztx &
               +csy(i,j,k)*a_zty &
               +csz(i,j,k)*a_ztz)/a_giac
            !

            ! Coriolis oil
            if (j.eq.jy) then

               g_co11(i,k)=0
               !
               g_co12(i,k)=0.
               !
               g_co13(i,k)=(-csx(i,j,k)*a_ztz &
                  +csz(i,j,k)*a_ztz)/a_giac

            end if
      
         end do
      end do
   end do
   !
   !-----------------------------------------------------------------------
   !
   ! planes at eta constant
   !
   do i=1,jx
      do j=0,jy
         do k=kparasta,kparaend !kinizio,kfine  !1,jz
            !
            if      (j.eq.0.and.jp.eq.1)  then
               !
               x2=.25*(x(i,j+2,k)+x(i,j+2,k-1)+x(i-1,j+2,k-1)+x(i-1,j+2,k))
               x1=.25*(x(i,j+1,k)+x(i,j+1,k-1)+x(i-1,j+1,k-1)+x(i-1,j+1,k))
               x0=.25*(x(i,j  ,k)+x(i,j  ,k-1)+x(i-1,j  ,k-1)+x(i-1,j  ,k))
               !
               y2=.25*(y(i,j+2,k)+y(i,j+2,k-1)+y(i-1,j+2,k-1)+y(i-1,j+2,k))
               y1=.25*(y(i,j+1,k)+y(i,j+1,k-1)+y(i-1,j+1,k-1)+y(i-1,j+1,k))
               y0=.25*(y(i,j  ,k)+y(i,j  ,k-1)+y(i-1,j  ,k-1)+y(i-1,j  ,k))
               !
               z2=.25*(z(i,j+2,k)+z(i,j+2,k-1)+z(i-1,j+2,k-1)+z(i-1,j+2,k))
               z1=.25*(z(i,j+1,k)+z(i,j+1,k-1)+z(i-1,j+1,k-1)+z(i-1,j+1,k))
               z0=.25*(z(i,j  ,k)+z(i,j  ,k-1)+z(i-1,j  ,k-1)+z(i-1,j  ,k))
               !
               xeta=.5 * ( -3.*x0 + 4.*x1 - x2 )
               yeta=.5 * ( -3.*y0 + 4.*y1 - y2 )
               zeta=.5 * ( -3.*z0 + 4.*z1 - z2 )
            !
            else if (j.eq.jy.and.jp.eq.1) then
               !
               x2=.25*(x(i,j-2,k)+x(i,j-2,k-1)+x(i-1,j-2,k-1)+x(i-1,j-2,k))
               x1=.25*(x(i,j-1,k)+x(i,j-1,k-1)+x(i-1,j-1,k-1)+x(i-1,j-1,k))
               x0=.25*(x(i,j  ,k)+x(i,j  ,k-1)+x(i-1,j  ,k-1)+x(i-1,j  ,k))
               !
               y2=.25*(y(i,j-2,k)+y(i,j-2,k-1)+y(i-1,j-2,k-1)+y(i-1,j-2,k))
               y1=.25*(y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k))
               y0=.25*(y(i,j  ,k)+y(i,j  ,k-1)+y(i-1,j  ,k-1)+y(i-1,j  ,k))
               !
               z2=.25*(z(i,j-2,k)+z(i,j-2,k-1)+z(i-1,j-2,k-1)+z(i-1,j-2,k))
               z1=.25*(z(i,j-1,k)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j-1,k))
               z0=.25*(z(i,j  ,k)+z(i,j  ,k-1)+z(i-1,j  ,k-1)+z(i-1,j  ,k))
               !
               xeta=.5 * ( 3.*x0 - 4.*x1 + x2 )
               yeta=.5 * ( 3.*y0 - 4.*y1 + y2 )
               zeta=.5 * ( 3.*z0 - 4.*z1 + z2 )
            !
            else
               !
               xsop=.25*(x(i,j+1,k)+x(i,j+1,k-1)+x(i-1,j+1,k-1)+x(i-1,j+1,k))
               xsot=.25*(x(i,j-1,k)+x(i,j-1,k-1)+x(i-1,j-1,k-1)+x(i-1,j-1,k))
               !
               ysop=.25*(y(i,j+1,k)+y(i,j+1,k-1)+y(i-1,j+1,k-1)+y(i-1,j+1,k))
               ysot=.25*(y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k))
               !
               zsop=.25*(z(i,j+1,k)+z(i,j+1,k-1)+z(i-1,j+1,k-1)+z(i-1,j+1,k))
               zsot=.25*(z(i,j-1,k)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j-1,k))
               !
               xeta=.5*(xsop-xsot)
               yeta=.5*(ysop-ysot)
               zeta=.5*(zsop-zsot)
            !
            end if
            !
            yzet=.5*(y(i,j,k  ) + y(i-1,j,k  ))-.5*(y(i,j,k-1) + y(i-1,j,k-1))
            !
            zcsi=.5*(z(i  ,j,k) + z(i  ,j,k-1))-.5*(z(i-1,j,k) + z(i-1,j,k-1))
            !
            ycsi=.5*(y(i  ,j,k) + y(i  ,j,k-1))-.5*(y(i-1,j,k) + y(i-1,j,k-1))
            !
            zzet=.5*(z(i,j,k  ) + z(i-1,j,k  ))-.5*(z(i,j,k-1) + z(i-1,j,k-1))
            !
            xcsi=.5*(x(i  ,j,k) + x(i  ,j,k-1))-.5*(x(i-1,j,k) + x(i-1,j,k-1))
            !
            xzet=.5*(x(i,j,k  ) + x(i-1,j,k  ))-.5*(x(i,j,k-1) + x(i-1,j,k-1))
            !
            etx(i,j,k) = yzet*zcsi - ycsi*zzet
            ety(i,j,k) = xcsi*zzet - xzet*zcsi
            etz(i,j,k) = xzet*ycsi - xcsi*yzet
            !
            a_csx = yeta*zzet - yzet*zeta
            a_csy = xzet*zeta - xeta*zzet
            a_csz = xeta*yzet - xzet*yeta
            !
            a_ztx=ycsi*zeta-yeta*zcsi
            a_zty=xeta*zcsi-xcsi*zeta
            a_ztz=xcsi*yeta-xeta*ycsi
            !
            a_giac=xcsi*(yeta*zzet-yzet*zeta)-xeta*(ycsi*zzet-yzet*zcsi)+xzet*(ycsi*zeta-yeta*zcsi)
            !
            g21(i,j,k)=(etx(i,j,k)*a_csx+ety(i,j,k)*a_csy+etz(i,j,k)*a_csz)/a_giac
            !
            g22(i,j,k)=(etx(i,j,k)**2+ety(i,j,k)**2+etz(i,j,k)**2)/a_giac
            !
            g23(i,j,k)=(etx(i,j,k)*a_ztx+ety(i,j,k)*a_zty+etz(i,j,k)*a_ztz)/a_giac
                  !
         end do
      end do
   end do
   !-----------------------------------------------------------------------
   ! temporary index for k direction
   !
   if(myid.eq.0)then
      kinizio = kparasta-1
      kfine   = kparaend
   elseif(myid.eq.nproc-1)then
      kinizio = kparasta-1
      kfine   = kparaend
   else
      kinizio = kparasta-1
      kfine   = kparaend
   end if
   !
   !-----------------------------------------------------------------------
   !
   ! planes at zita constant
   !
   ! myid=0 need to start with k=0
   !

   do i=1,jx
      do j=1,jy
         do k=kparasta-1,kparaend ! needs z at kparaend+1
            !
            if (k.eq.0.and.kp.eq.1) then   !only myid=0
               !
               x2=.25*(x(i,j,k+2)+x(i,j-1,k+2)+x(i-1,j-1,k+2)+x(i-1,j,k+2))
               x1=.25*(x(i,j,k+1)+x(i,j-1,k+1)+x(i-1,j-1,k+1)+x(i-1,j,k+1))
               x0=.25*(x(i,j,k  )+x(i,j-1,k  )+x(i-1,j-1,k  )+x(i-1,j,k  ))
               !
               y2=.25*(y(i,j,k+2)+y(i,j-1,k+2)+y(i-1,j-1,k+2)+y(i-1,j,k+2))
               y1=.25*(y(i,j,k+1)+y(i,j-1,k+1)+y(i-1,j-1,k+1)+y(i-1,j,k+1))
               y0=.25*(y(i,j,k  )+y(i,j-1,k  )+y(i-1,j-1,k  )+y(i-1,j,k  ))
               !
               z2=.25*(z(i,j,k+2)+z(i,j-1,k+2)+z(i-1,j-1,k+2)+z(i-1,j,k+2))
               z1=.25*(z(i,j,k+1)+z(i,j-1,k+1)+z(i-1,j-1,k+1)+z(i-1,j,k+1))
               z0=.25*(z(i,j,k  )+z(i,j-1,k  )+z(i-1,j-1,k  )+z(i-1,j,k  ))
               !
               xzet=.5 * ( -3.*x0 + 4.*x1 - x2 )
               yzet=.5 * ( -3.*y0 + 4.*y1 - y2 )
               zzet=.5 * ( -3.*z0 + 4.*z1 - z2 )
            !
            else if (k.eq.jz.and.kp.eq.1) then !only myid=nproc-1
               !
               x2=.25*(x(i,j,k-2)+x(i,j-1,k-2)+x(i-1,j-1,k-2)+x(i-1,j,k-2))
               x1=.25*(x(i,j,k-1)+x(i,j-1,k-1)+x(i-1,j-1,k-1)+x(i-1,j,k-1))
               x0=.25*(x(i,j,k  )+x(i,j-1,k  )+x(i-1,j-1,k  )+x(i-1,j,k  ))
               !
               y2=.25*(y(i,j,k-2)+y(i,j-1,k-2)+y(i-1,j-1,k-2)+y(i-1,j,k-2))
               y1=.25*(y(i,j,k-1)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j,k-1))
               y0=.25*(y(i,j,k  )+y(i,j-1,k  )+y(i-1,j-1,k  )+y(i-1,j,k  ))
               !
               z2=.25*(z(i,j,k-2)+z(i,j-1,k-2)+z(i-1,j-1,k-2)+z(i-1,j,k-2))
               z1=.25*(z(i,j,k-1)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j,k-1))
               z0=.25*(z(i,j,k  )+z(i,j-1,k  )+z(i-1,j-1,k  )+z(i-1,j,k  ))
               !
               xzet=.5 * ( 3.*x0 - 4.*x1 + x2 )
               yzet=.5 * ( 3.*y0 - 4.*y1 + y2 )
               zzet=.5 * ( 3.*z0 - 4.*z1 + z2 )
            !
            else
               !
               xav=.25*(x(i,j,k+1)+x(i,j-1,k+1)+x(i-1,j-1,k+1)+x(i-1,j,k+1))
               xdt=.25*(x(i,j,k-1)+x(i,j-1,k-1)+x(i-1,j-1,k-1)+x(i-1,j,k-1))
               !
               yav=.25*(y(i,j,k+1)+y(i,j-1,k+1)+y(i-1,j-1,k+1)+y(i-1,j,k+1))
               ydt=.25*(y(i,j,k-1)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j,k-1))
               !
               zav=.25*(z(i,j,k+1)+z(i,j-1,k+1)+z(i-1,j-1,k+1)+z(i-1,j,k+1))
               zdt=.25*(z(i,j,k-1)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j,k-1))
               !
               xzet=.5*(xav-xdt)
               yzet=.5*(yav-ydt)
               zzet=.5*(zav-zdt)
            !
            end if
            !
            ycsi=.5*(y(i  ,j,k) + y(i  ,j-1,k))-.5*(y(i-1,j,k) + y(i-1,j-1,k))
            !
            zeta=.5*(z(i,j  ,k) + z(i-1,j  ,k))-.5*(z(i,j-1,k) + z(i-1,j-1,k))
            !
            yeta=.5*(y(i,j  ,k) + y(i-1,j  ,k))-.5*(y(i,j-1,k) + y(i-1,j-1,k))
            !
            zcsi=.5*(z(i  ,j,k) + z(i  ,j-1,k))-.5*(z(i-1,j,k) + z(i-1,j-1,k))
            !
            xcsi=.5*(x(i  ,j,k) + x(i  ,j-1,k))-.5*(x(i-1,j,k) + x(i-1,j-1,k))
            !
            xeta=.5*(x(i,j  ,k) + x(i-1,j  ,k))-.5*(x(i,j-1,k) + x(i-1,j-1,k))
            !
            a_csx = yeta*zzet - yzet*zeta
            a_csy = xzet*zeta - xeta*zzet
            a_csz = xeta*yzet - xzet*yeta
            !
            a_etx = yzet*zcsi - ycsi*zzet
            a_ety = xcsi*zzet - xzet*zcsi
            a_etz = xzet*ycsi - xcsi*yzet
            !
            ztx(i,j,k)=ycsi*zeta-yeta*zcsi
            zty(i,j,k)=xeta*zcsi-xcsi*zeta
            ztz(i,j,k)=xcsi*yeta-xeta*ycsi
            !
            a_giac=xcsi*(yeta*zzet-yzet*zeta)-xeta*(ycsi*zzet-yzet*zcsi)+xzet*(ycsi*zeta-yeta*zcsi)
            !
            g31(i,j,k)=(ztx(i,j,k)*a_csx+zty(i,j,k)*a_csy+ztz(i,j,k)*a_csz)/a_giac
            !
            g32(i,j,k)=(ztx(i,j,k)*a_etx+zty(i,j,k)*a_ety+ztz(i,j,k)*a_etz)/a_giac
            !
            g33(i,j,k)=(ztx(i,j,k)**2+zty(i,j,k)**2+ztz(i,j,k)**2)/a_giac
            !
            ! Coriolis oil
            if(j.eq.jy)then

               g_co31(i,k)=(-ztx(i,j,k)*a_csz+ztz(i,j,k)*a_csx)/a_giac
               !
               g_co32(i,k)=0.
               !
               g_co33(i,k)=0.

            end if


         end do
      end do
   end do
   !
   !-----------------------------------------------------------------------
   !
   ! giacobian J(-1) : giac(i,j,k) at the cell centroid
   !
   !-----------------------------------------------------------------------
   ! compute giacobian at the centroid

   !     check on grid quality
   giac_min=10000000000.
   giac_max=0.
   !
   do i=1,jx
      do j=1,jy
         do k=kparasta,kparaend  !1,jz
            !
            xcsi=.25*(x(i  ,j,k)+x(i  ,j,k-1) &
               +x(i  ,j-1,k-1)+x(i  ,j-1,k))- &
               .25*(x(i-1,j,k)+x(i-1,j,k-1) &
               +x(i-1,j-1,k-1)+x(i-1,j-1,k))
            !
            ycsi=.25*(y(i  ,j,k)+y(i  ,j,k-1) &
               +y(i  ,j-1,k-1)+y(i  ,j-1,k))- &
               .25*(y(i-1,j,k)+y(i-1,j,k-1) &
               +y(i-1,j-1,k-1)+y(i-1,j-1,k))
            !
            zcsi=.25*(z(i  ,j,k)+z(i  ,j,k-1) &
               +z(i  ,j-1,k-1)+z(i  ,j-1,k))- &
               .25*(z(i-1,j,k)+z(i-1,j,k-1) &
               +z(i-1,j-1,k-1)+z(i-1,j-1,k))
            !
            xeta=.25*(x(i,j  ,k)+x(i,j  ,k-1) &
               +x(i-1,j  ,k-1)+x(i-1,j  ,k))- &
               .25*(x(i,j-1,k)+x(i,j-1,k-1) &
               +x(i-1,j-1,k-1)+x(i-1,j-1,k))
            !
            yeta=.25*(y(i,j  ,k)+y(i,j  ,k-1) &
               +y(i-1,j  ,k-1)+y(i-1,j  ,k))- &
               .25*(y(i,j-1,k)+y(i,j-1,k-1) &
               +y(i-1,j-1,k-1)+y(i-1,j-1,k))
            !
            zeta=.25*(z(i,j  ,k)+z(i,j  ,k-1) &
               +z(i-1,j  ,k-1)+z(i-1,j  ,k))- &
               .25*(z(i,j-1,k)+z(i,j-1,k-1) &
               +z(i-1,j-1,k-1)+z(i-1,j-1,k))
            !
            xzet=.25*(x(i,j,k  )+x(i,j-1,k  ) &
               +x(i-1,j-1,k  )+x(i-1,j,k  ))- &
               .25*(x(i,j,k-1)+x(i,j-1,k-1) &
               +x(i-1,j-1,k-1)+x(i-1,j,k-1))
            !
            yzet=.25*(y(i,j,k  )+y(i,j-1,k  ) &
               +y(i-1,j-1,k  )+y(i-1,j,k  ))- &
               .25*(y(i,j,k-1)+y(i,j-1,k-1) &
               +y(i-1,j-1,k-1)+y(i-1,j,k-1))
            !
            zzet=.25*(z(i,j,k  )+z(i,j-1,k  ) &
               +z(i-1,j-1,k  )+z(i-1,j,k  ))- &
               .25*(z(i,j,k-1)+z(i,j-1,k-1) &
               +z(i-1,j-1,k-1)+z(i-1,j,k-1))
            !
            giac(i,j,k)=xcsi*(yeta*zzet-yzet*zeta)- &
               xeta*(ycsi*zzet-yzet*zcsi)+ &
               xzet*(ycsi*zeta-yeta*zcsi)
            !
            giac_min=min(giac_min,giac(i,j,k))
            giac_max=max(giac_max,giac(i,j,k))

            if(giac(i,j,k)<1.e-10) then
               write(*,*) 'METRICA: grid not good cell very small'
            end if
      
         end do
      end do
   end do

   if(myid.eq.0) then
      write(*,*)'giac min, giac max',giac_min,giac_max
   end if
   !
   return
end
