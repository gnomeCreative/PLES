!***********************************************************************
subroutine gterza(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet, &
   t11,t22,t33)
   !***********************************************************************
   !
   ! compute the elements for the third line of the controvariant
   ! metric tensor for upper grid for multigrid cycle
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
   csx = yeta*zzet - yzet*zeta
   csy = xzet*zeta - xeta*zzet
   csz = xeta*yzet - xzet*yeta
   !
   etx = yzet*zcsi - ycsi*zzet
   ety = xcsi*zzet - xzet*zcsi
   etz = xzet*ycsi - xcsi*yzet
   !
   ztx=  ycsi*zeta - yeta*zcsi
   zty=  xeta*zcsi - xcsi*zeta
   ztz=  xcsi*yeta - xeta*ycsi
   !
   giac=     xcsi*(yeta*zzet-yzet*zeta)- &
      xeta*(ycsi*zzet-yzet*zcsi)+ &
      xzet*(ycsi*zeta-yeta*zcsi)
   !
   t11=(ztx*csx+zty*csy+ztz*csz)/giac
   !
   t22=(ztx*etx+zty*ety+ztz*etz)/giac
   !
   t33=(ztx**2+zty**2+ztz**2)/giac
   !
   return
end
