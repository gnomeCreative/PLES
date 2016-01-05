!***********************************************************************
      module myarrays_buffer_bodyforce
!***********************************************************************
      real, allocatable :: tkepom1(:,:,:),tkepom2(:,:,:)
      real, allocatable :: tkepom5(:,:,:),tkepom6(:,:,:)

      real, allocatable :: tke1(:,:),tke2(:,:)
      real, allocatable :: tke5(:,:),tke6(:,:)

!     face 1
      real, allocatable :: noise_f1_xtime(:,:),old_noise_f1_xtime(:)
      real, allocatable :: noise_f1_ytime(:,:),old_noise_f1_ytime(:)
      real, allocatable :: noise_f1_ztime(:,:),old_noise_f1_ztime(:)      
      real, allocatable :: noise_f1_xspace(:,:)
      real, allocatable :: noise_f1_yspace(:,:)
      real, allocatable :: noise_f1_zspace(:,:)	
      real, allocatable :: bcsi_f1(:,:,:)
      real, allocatable :: beta_f1(:,:,:)
      real, allocatable :: bzet_f1(:,:,:)	       

!     face 2
      real, allocatable :: noise_f2_xtime(:,:),old_noise_f2_xtime(:)
      real, allocatable :: noise_f2_ytime(:,:),old_noise_f2_ytime(:)
      real, allocatable :: noise_f2_ztime(:,:),old_noise_f2_ztime(:)      
      real, allocatable :: noise_f2_xspace(:,:)
      real, allocatable :: noise_f2_yspace(:,:)
      real, allocatable :: noise_f2_zspace(:,:)	 
      real, allocatable :: bcsi_f2(:,:,:)
      real, allocatable :: beta_f2(:,:,:)
      real, allocatable :: bzet_f2(:,:,:)
      
!     face 3
      real, allocatable :: noise_f5_xtime(:,:),old_noise_f5_xtime(:)
      real, allocatable :: noise_f5_ytime(:,:),old_noise_f5_ytime(:)
      real, allocatable :: noise_f5_ztime(:,:),old_noise_f5_ztime(:)      
      real, allocatable :: noise_f5_xspace(:,:)
      real, allocatable :: noise_f5_yspace(:,:)
      real, allocatable :: noise_f5_zspace(:,:)	 
      real, allocatable :: bcsi_f5(:,:,:)
      real, allocatable :: beta_f5(:,:,:)
      real, allocatable :: bzet_f5(:,:,:)
      
!     face 4
      real, allocatable :: noise_f6_xtime(:,:),old_noise_f6_xtime(:)
      real, allocatable :: noise_f6_ytime(:,:),old_noise_f6_ytime(:)
      real, allocatable :: noise_f6_ztime(:,:),old_noise_f6_ztime(:)	  
      real, allocatable :: noise_f6_xspace(:,:)
      real, allocatable :: noise_f6_yspace(:,:)
      real, allocatable :: noise_f6_zspace(:,:)	 
      real, allocatable :: bcsi_f6(:,:,:)
      real, allocatable :: beta_f6(:,:,:)
      real, allocatable :: bzet_f6(:,:,:)
      
      real corr_factor
      
      integer ispon,kspon
      integer ntpom

!***********************************************************************
      end module myarrays_buffer_bodyforce
!***********************************************************************
