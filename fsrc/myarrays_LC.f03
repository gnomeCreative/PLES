!***********************************************************************
      module myarrays_LC
!***********************************************************************
!     for langmuir circulation
!-----------------------------------------------------------------------
      real :: betaw
      real :: lamb
      real :: h_0
      
      real, allocatable :: u_drift(:,:,:)
      real, allocatable :: w_drift(:,:,:)
      real, allocatable :: vortx(:,:,:)
      real ,allocatable :: vorty(:,:,:)
      real, allocatable :: vortz(:,:,:)  
   
      real, allocatable :: ucs(:,:,:)
      real, allocatable :: wcs(:,:,:)

!***********************************************************************
      end module myarrays_LC
!***********************************************************************
