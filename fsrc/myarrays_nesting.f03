!***********************************************************************
       module myarrays_nesting
!***********************************************************************
!      for nesting procedure
!-----------------------------------------------------------------------
       use iso_c_binding

       real,allocatable,dimension(:,:) :: uo1,vo1,wo1,un1,vn1,wn1
       real,allocatable,dimension(:,:) :: uo2,vo2,wo2,un2,vn2,wn2
       real,allocatable,dimension(:,:) :: uo5,vo5,wo5,un5,vn5,wn5
       real,allocatable,dimension(:,:) :: uo6,vo6,wo6,un6,vn6,wn6
       real,allocatable,dimension(:,:,:)::rhovo1,rhovo2,rhovo5,rhovo6
       real,allocatable,dimension(:,:,:)::rhovn1,rhovn2,rhovn5,rhovn6
       
!       real,allocatable,dimension(:,:) :: up1,vp1,wp1
!       real,allocatable,dimension(:,:) :: up2,vp2,wp2
!       real,allocatable,dimension(:,:) :: up5,vp5,wp5
!       real,allocatable,dimension(:,:) :: up6,vp6,wp6

       real,allocatable,dimension(:,:) :: ucp1
       real,allocatable,dimension(:,:) :: ucp2
       real,allocatable,dimension(:,:) :: wcp5
       real,allocatable,dimension(:,:) :: wcp6

!       real,allocatable,dimension(:,:,:)::rhovp1,rhovp2,rhovp5,rhovp6
     
       
       real ti_pom_old,ti_pom_new,ti_pom_fin
       integer n_ti_pom       
       logical,bind(C) :: termina
       integer ntke
!
!***********************************************************************
       end module myarrays_nesting
!***********************************************************************
