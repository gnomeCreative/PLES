!***********************************************************************
module myarrays_WB
    !***********************************************************************
    !     for wavebreaking
    !-----------------------------------------------------------------------
    use iso_c_binding

    integer :: old_j
    real,bind(C)    :: l_0
    real,bind(C)    :: alpha
    real,bind(C)    :: c10
    !
    real a_wind
    real, allocatable :: u_wind(:,:)
    real, allocatable :: w_wind(:,:)
    ! Giulia modificavento: variabili nuove
    real, allocatable :: tauu_att(:,:)
    real, allocatable :: tauw_att(:,:)
    ! Giulia modificavento:
    real, allocatable :: u_att(:,:)
    real, allocatable :: w_att(:,:)
    real, allocatable :: v_att(:,:)
    real, allocatable :: Fx(:,:)
    real, allocatable :: Fz(:,:)
    logical :: cf_exists
    real, allocatable :: cf(:,:)
    real :: varcf
!***********************************************************************
end module myarrays_WB
!***********************************************************************

