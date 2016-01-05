!***********************************************************************
      module myarrays_moisture
!***********************************************************************
!      for moisture
!-----------------------------------------------------------------------
      integer imoist	  
      real Tref
      real Qref 
      real betaT
      real betaQ
      real Ma
      real Mv
      
      real Lv
      real Gdry
      real Rd
      real cpd
      
      real, allocatable :: tpotm(:), qm(:)
      real, allocatable :: tauw3nest(:)
!      real, allocatable :: heatflux3nest(:)
!      real, allocatable :: vapflux3nest(:)
      real, allocatable :: scalarflux3nest(:,:)    
!***********************************************************************
      end module myarrays_moisture
!***********************************************************************
