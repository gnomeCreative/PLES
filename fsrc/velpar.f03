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

        use scala3, only: n1,n2,nscal,kparasta,kparaend

        implicit none

        allocate(usn(1:n2,kparasta:kparaend),vsn(1:n2,kparasta:kparaend),wsn(1:n2,kparasta:kparaend))
        allocate(udx(1:n2,kparasta:kparaend),vdx(1:n2,kparasta:kparaend),wdx(1:n2,kparasta:kparaend))
        allocate(rhosn(nscal,1:n2,kparasta:kparaend),rhodx(nscal,1:n2,kparasta:kparaend))
        !
        allocate(usp(1:n1,kparasta:kparaend),vsp(1:n1,kparasta:kparaend),wsp(1:n1,kparasta:kparaend))
        allocate(ust(1:n1,kparasta:kparaend),vst(1:n1,kparasta:kparaend),wst(1:n1,kparasta:kparaend))
        allocate(rhosp(nscal,1:n1,kparasta:kparaend),rhost(nscal,1:n1,kparasta:kparaend))
        !
        allocate(uav(n1,n2),vav(n1,n2),wav(n1,n2))
        allocate(uin(n1,n2),vin(n1,n2),win(n1,n2))
        allocate(rhoav(nscal,n1,n2),rhoin(nscal,n1,n2))

    end subroutine initialize_velpar

end module velpar
