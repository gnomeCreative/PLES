!***********************************************************************
       module myarrays_wallmodel
!***********************************************************************
!      for wall modeling
!-----------------------------------------------------------------------
       real,allocatable :: tauw3x(:,:),tauw4x(:,:)
       real,allocatable :: tauw3z(:,:),tauw4z(:,:)
       real,allocatable :: punto_wfp3(:,:,:,:)
       real,allocatable :: punto_wfp4(:,:,:,:)  
       real,allocatable :: u_t(:,:,:)

       integer,allocatable ::att_mod_par(:,:,:)  
       real,allocatable :: utangente(:,:,:)

       integer wfp1,wfp2,wfp3,wfp4,wfp5,wfp6
       integer eseguo34,att_wm_sgs,rough
       real z0

!***********************************************************************
       end module myarrays_wallmodel
!***********************************************************************
