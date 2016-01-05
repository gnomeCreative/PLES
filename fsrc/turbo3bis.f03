!***********************************************************************
      module turbo3bis
!***********************************************************************
!  array for turbulence model
!-----------------------------------------------------------------------
      real,allocatable :: smod(:,:,:)
      real,allocatable :: smodV(:,:,:)
      real,allocatable :: smodH(:,:,:)

      real,allocatable :: s11(:,:,:)
      real,allocatable :: s12(:,:,:)
      real,allocatable :: s13(:,:,:)
      real,allocatable :: s21(:,:,:)
      real,allocatable :: s22(:,:,:)
      real,allocatable :: s23(:,:,:)
      real,allocatable :: s31(:,:,:)
      real,allocatable :: s32(:,:,:)
      real,allocatable :: s33(:,:,:)
!
      real,allocatable :: uco(:,:,:)
      real,allocatable :: vco(:,:,:)
      real,allocatable :: wco(:,:,:)
      real,allocatable :: uuco(:,:,:)
      real,allocatable :: uvco(:,:,:)
      real,allocatable :: uwco(:,:,:)
      real,allocatable :: vuco(:,:,:)
      real,allocatable :: vvco(:,:,:)
      real,allocatable :: vwco(:,:,:)
      real,allocatable :: wuco(:,:,:)
      real,allocatable :: wvco(:,:,:)
      real,allocatable :: wwco(:,:,:)
!
      real,allocatable :: rhofl(:,:,:)
      real,allocatable :: rho11(:,:,:,:)
      real,allocatable :: rho22(:,:,:,:)
      real,allocatable :: rho33(:,:,:,:)
      real,allocatable :: rhouco(:,:,:)
      real,allocatable :: rhovco(:,:,:)
      real,allocatable :: rhowco(:,:,:)
!
      real,allocatable :: apcsx(:,:,:)
      real,allocatable :: apcsy(:,:,:)
      real,allocatable :: apcsz(:,:,:)
      real,allocatable :: apetx(:,:,:)
      real,allocatable :: apety(:,:,:)
      real,allocatable :: apetz(:,:,:)
      real,allocatable :: apztx(:,:,:)
      real,allocatable :: apzty(:,:,:)
      real,allocatable :: apztz(:,:,:)
!
      real,allocatable :: ap11(:,:,:)
      real,allocatable :: ap12(:,:,:)
      real,allocatable :: ap13(:,:,:)
      real,allocatable :: ap21(:,:,:)
      real,allocatable :: ap22(:,:,:)
      real,allocatable :: ap23(:,:,:)
      real,allocatable :: ap31(:,:,:)
      real,allocatable :: ap32(:,:,:)
      real,allocatable :: ap33(:,:,:)
      real,allocatable :: pp0(:,:,:)
      real,allocatable :: pp1(:,:,:)
      real,allocatable :: pp2(:,:,:)
      real,allocatable :: pp3(:,:,:)
      real,allocatable :: pc1(:,:,:)
      real,allocatable :: pc2(:,:,:)
      real,allocatable :: pc3(:,:,:)	  
      real,allocatable :: p0a(:,:,:)
      real,allocatable :: p0b(:,:,:)

      real,allocatable :: sbuff1(:)
      real,allocatable :: rbuff1(:)
!
!***********************************************************************
      end module turbo3bis 
!***********************************************************************
