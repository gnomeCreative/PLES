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

contains

    subroutine initialize_velpar()

        use scala3, only: n1,n2,n3,nscal

        implicit none

        allocate(usn(n2,n3),vsn(n2,n3),wsn(n2,n3))
        allocate(udx(n2,n3),vdx(n2,n3),wdx(n2,n3))
        allocate(rhosn(nscal,n2,n3),rhodx(nscal,n2,n3))
        !
        allocate(usp(n1,n3),vsp(n1,n3),wsp(n1,n3))
        allocate(ust(n1,n3),vst(n1,n3),wst(n1,n3))
        allocate(rhosp(nscal,n1,n3),rhost(nscal,n1,n3))
        !
        allocate(uav(n1,n2),vav(n1,n2),wav(n1,n2))
        allocate(uin(n1,n2),vin(n1,n2),win(n1,n2))
        allocate(rhoav(nscal,n1,n2),rhoin(nscal,n1,n2))

    end subroutine initialize_velpar

end module velpar
