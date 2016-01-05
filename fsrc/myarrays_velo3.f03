!***********************************************************************
      module myarrays_velo3
!***********************************************************************
! contains array like velocity etc
      implicit none
!-----------------------------------------------------------------------
      real,allocatable :: u(:,:,:),v(:,:,:),w(:,:,:),rhov(:,:,:,:)
      real,allocatable :: uc(:,:,:),vc(:,:,:),wc(:,:,:)
      real,allocatable :: rhs(:,:,:),fi(:,:,:),next_prs(:,:)

      real,allocatable ::  u_piano(:,:,:),v_piano(:,:,:),w_piano(:,:,:)    
      real,allocatable ::  rhov_piano(:,:,:,:)
      real,allocatable ::  uc_piano(:,:,:)
      real,allocatable ::  vc_piano(:,:,:)   
      real,allocatable ::  wc_piano(:,:,:)

      real,allocatable ::  cgra1(:,:,:),cgra2(:,:,:),cgra3(:,:,:)
      real,allocatable ::  gra1(:,:,:), gra2(:,:,:), gra3(:,:,:)

      real,allocatable ::  gra1_appoggio(:,:,:)
      real,allocatable ::  gra2_appoggio(:,:,:)
      real,allocatable ::  gra3_appoggio(:,:,:)      

      real,allocatable :: delu(:,:,:),delv(:,:,:),delw(:,:,:)
	
! arrays for inflow
      real,allocatable :: up1(:,:),vp1(:,:),wp1(:,:),rhovp1(:,:,:)
      real,allocatable :: up2(:,:),vp2(:,:),wp2(:,:),rhovp2(:,:,:)
      real,allocatable :: up3(:,:),vp3(:,:),wp3(:,:),rhovp3(:,:,:)
      real,allocatable :: up4(:,:),vp4(:,:),wp4(:,:),rhovp4(:,:,:)
      real,allocatable :: up5(:,:),vp5(:,:),wp5(:,:),rhovp5(:,:,:)
      real,allocatable :: up6(:,:),vp6(:,:),wp6(:,:),rhovp6(:,:,:)
	
! array for pressure
      real,allocatable :: cs1(:,:),cs2(:,:)
      real,allocatable :: cs3(:,:),cs4(:,:)
      real,allocatable :: cs5(:,:),cs6(:,:)   

      real,allocatable :: uc1_orl(:,:),uc2_orl(:,:)
      real,allocatable :: vc3_orl(:,:),vc4_orl(:,:)
      real,allocatable :: wc5_orl(:,:),wc6_orl(:,:)   
		
!***********************************************************************
      end module myarrays_velo3
!***********************************************************************





