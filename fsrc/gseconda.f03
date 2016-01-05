!***********************************************************************
subroutine gseconda(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet, &
   t11,t22,t33)
   !***********************************************************************
   !
   ! compute the second line elements of controvariant metric tensor on
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
   etx = yzet*zcsi - ycsi*zzet
   ety = xcsi*zzet - xzet*zcsi
   etz = xzet*ycsi - xcsi*yzet
   !
   csx = yeta*zzet - yzet*zeta
   csy = xzet*zeta - xeta*zzet
   csz = xeta*yzet - xzet*yeta
   !
   ztx=ycsi*zeta-yeta*zcsi
   zty=xeta*zcsi-xcsi*zeta
   ztz=xcsi*yeta-xeta*ycsi
   !
   giac=xcsi*(yeta*zzet-yzet*zeta)- &
      xeta*(ycsi*zzet-yzet*zcsi)+ &
      xzet*(ycsi*zeta-yeta*zcsi)
   !
   t11=(etx*csx+ety*csy+etz*csz)/giac
   !
   t22=(etx**2+ety**2+etz**2)/giac
   !
   t33=(etx*ztx+ety*zty+etz*ztz)/giac
   !
   return
end
