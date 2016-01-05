module velpar

   use scala3, only: nscal,n1,n2,n3

   !
   !  common per parametri relativi alle velocità parete
   !

   real   usn(n2,n3)
   real   vsn(n2,n3)
   real   wsn(n2,n3)
   real rhosn(nscal,n2,n3)
   real   udx(n2,n3)
   real   vdx(n2,n3)
   real   wdx(n2,n3)
   real rhodx(nscal,n2,n3)
   !
   real   usp(n1,n3)
   real   vsp(n1,n3)
   real   wsp(n1,n3)
   real rhosp(nscal,n1,n3)
   real   ust(n1,n3)
   real   vst(n1,n3)
   real   wst(n1,n3)
   real rhost(nscal,n1,n3)
   !
   real   uav(n1,n2)
   real   vav(n1,n2)
   real   wav(n1,n2)
   real rhoav(nscal,n1,n2)
   real   uin(n1,n2)
   real   vin(n1,n2)
   real   win(n1,n2)
   real rhoin(nscal,n1,n2)
   !
   common/velparcommon/usn,vsn,wsn,rhosn,udx,vdx,wdx,rhodx, &
      ust,vst,wst,rhost,usp,vsp,wsp,rhosp, &
      uav,vav,wav,rhoav,uin,vin,win,rhoin

end module velpar
