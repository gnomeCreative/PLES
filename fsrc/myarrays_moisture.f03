!***********************************************************************
module myarrays_moisture
    !***********************************************************************
    !      for moisture
    !-----------------------------------------------------------------------
    use,intrinsic :: iso_c_binding

    integer,bind(C) :: imoist
    real,bind(C) :: Tref
    real,bind(C) :: Qref
    real,bind(C) :: betaT
    real,bind(C) :: betaQ
    real,bind(C) :: Ma
    real,bind(C) :: Mv
      
    real,bind(C) :: Lv
    real,bind(C) :: Gdry
    real,bind(C) :: Rd
    real,bind(C) :: cpd
      
    real, allocatable :: tpotm(:), qm(:)
    real, allocatable :: tauw3nest(:)
    !      real, allocatable :: heatflux3nest(:)
    !      real, allocatable :: vapflux3nest(:)
    real, allocatable :: scalarflux3nest(:,:)
!***********************************************************************
end module myarrays_moisture
!***********************************************************************
