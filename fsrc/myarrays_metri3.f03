!***********************************************************************
      module myarrays_metri3
!***********************************************************************
! coordinates x,y,z
! 19 metric terms, 6 for each face plus the giacobian
! fluxes
! eddy viscosity
! NOTE better to move annit and fluxes in an other module
      implicit none
!-----------------------------------------------------------------------
      real,allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
      real,allocatable :: csx(:,:,:),csy(:,:,:),csz(:,:,:)
      real,allocatable :: etx(:,:,:),ety(:,:,:),etz(:,:,:)
      real,allocatable :: ztx(:,:,:),zty(:,:,:),ztz(:,:,:)

      real,allocatable :: g11(:,:,:),g12(:,:,:),g13(:,:,:)
      real,allocatable :: g21(:,:,:),g22(:,:,:),g23(:,:,:)
      real,allocatable :: g31(:,:,:),g32(:,:,:),g33(:,:,:)
      
      real,allocatable :: giac(:,:,:)                      

      real,allocatable :: f1(:,:,:),f2(:,:,:),f3(:,:,:)
      real,allocatable :: annit(:,:,:),annitV(:,:,:)
      real,allocatable :: annit_piano(:,:,:),annitV_piano(:,:,:)     

!CORIOLIS OIL
      real,allocatable :: g_co11(:,:),g_co12(:,:),g_co13(:,:)
      real,allocatable :: g_co31(:,:),g_co32(:,:),g_co33(:,:)

                       
!***********************************************************************
      end module myarrays_metri3
!***********************************************************************
