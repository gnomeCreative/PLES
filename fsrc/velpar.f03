module velpar

   !
   !  common per parametri relativi alle velocità parete
   !

   real,allocatable :: usn(:,:),vsn(:,:),wsn(:,:)
   real,allocatable :: udx(:,:),vdx(:,:),wdx(:,:)
   real,allocatable :: rhosn(:,:,:),rhodx(:,:,:)
   !
   real,allocatable :: usp(:,:),vsp(:,:),wsp(:,:)
   real,allocatable :: ust(:,:),vst(:,:),wst(:,:)
   real,allocatable :: rhosp(:,:,:),rhost(:,:,:)
   !
   real,allocatable :: uav(:,:),vav(:,:),wav(:,:)
   real,allocatable :: uin(:,:),vin(:,:),win(:,:)
   real,allocatable :: rhoav(:,:,:),rhoin(:,:,:)
   !
   !common/velparcommon/usn,vsn,wsn,rhosn,udx,vdx,wdx,rhodx, &
   !   ust,vst,wst,rhost,usp,vsp,wsp,rhosp, &
   !   uav,vav,wav,rhoav,uin,vin,win,rhoin

end module velpar
