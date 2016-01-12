!***********************************************************************
module myarrays_wallmodel
    !***********************************************************************
    !      for wall modeling
    !-----------------------------------------------------------------------
    use iso_c_binding

    real,allocatable :: tauw3x(:,:),tauw4x(:,:)
    real,allocatable :: tauw3z(:,:),tauw4z(:,:)
    real,allocatable :: punto_wfp3(:,:,:,:)
    real,allocatable :: punto_wfp4(:,:,:,:)
    real,allocatable :: u_t(:,:,:)

    integer,allocatable ::att_mod_par(:,:,:)
    real,allocatable :: utangente(:,:,:)

    integer,bind(C) :: wfp1,wfp2,wfp3,wfp4,wfp5,wfp6
    integer :: eseguo34
    integer,bind(C) :: att_wm_sgs,rough
    real,bind(C) :: z0

!***********************************************************************
end module myarrays_wallmodel
!***********************************************************************
