!***********************************************************************
subroutine gprima(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet, &
   t11,t22,t33)
   !***********************************************************************
   !
   ! compute the first line elements of controvariant metric tensor on
   ! higher grid level for multigrid
   !
   implicit none
   !
   !-----------------------------------------------------------------------
   !     variables declaration
   real giac,t11,t22,t33
   real xcsi,ycsi,zcsi
   real xeta,yeta,zeta
   real xzet,yzet,zzet
   real csx,csy,csz
   real etx,ety,etz
   real ztx,zty,ztz
   !-----------------------------------------------------------------------
   !
   giac=     xcsi*(yeta*zzet-yzet*zeta)- &
      xeta*(ycsi*zzet-yzet*zcsi)+ &
      xzet*(ycsi*zeta-yeta*zcsi)
   !
   csx = yeta*zzet - yzet*zeta
   csy = xzet*zeta - xeta*zzet
   csz = xeta*yzet - xzet*yeta
   !
   etx = yzet*zcsi - ycsi*zzet
   ety = xcsi*zzet - xzet*zcsi
   etz = xzet*ycsi - xcsi*yzet
   !
   ztx=ycsi*zeta-yeta*zcsi
   zty=xeta*zcsi-xcsi*zeta
   ztz=xcsi*yeta-xeta*ycsi
   !
   t11=(csx**2+csy**2+csz**2)/giac
   !
   t22=(csx*etx+csy*ety+csz*etz)/giac
   !
   t33=(csx*ztx+csy*zty+csz*ztz)/giac
   !
   return
end
