!***********************************************************************
module myarrays_ibm
    !***********************************************************************
    !     array for immersed boundary method
    !-----------------------------------------------------------------------
    !

    use,intrinsic :: iso_c_binding

    integer(kind=c_int),bind(C) :: num_iter
    logical(kind=c_bool),bind(C) :: particles,bodyforce,bodypressure
    !MN,MP matrix dimension to store mesh: points and triangle
    integer :: MN,MP
    integer :: numero_celle_IB,numero_celle_bloccate,numero_celle_bloccate_real
    integer :: num_ib,num_solide

    integer,allocatable :: indici_CELLE_IB(:,:)
    integer,allocatable :: indici_celle_bloccate(:,:)

    real,allocatable :: distanze_CELLE_IB(:,:)
    real,allocatable :: dist_pp_ib(:)
    real,allocatable :: dist_ib_parete(:)
    real,allocatable :: ustar(:),pressure_ib(:,:),shear_ib(:,:)
    integer,allocatable :: caso_ib(:)
    real,allocatable :: proiezioni(:,:)
      
    !     array for rotation with eulerian angles c
    real, allocatable :: rot(:,:,:) ! PROBLEMATIC
    real, allocatable :: rot_inverse(:,:,:) ! PROBLEMATIC
      
    !     array for trilinear
    real, allocatable :: tricoef(:,:)
    integer, allocatable :: trind(:,:,:)
      
    !     buffer
    real, allocatable :: sbuff_ibm(:)
    real, allocatable :: rbuff_ibm(:)
    !
    !     to send the ib stencil
    integer :: num_left_snd,num_right_snd
    integer :: num_left_rcv,num_right_rcv
    integer, allocatable :: stencil_left_snd(:,:)
    integer, allocatable :: stencil_left_rcv(:,:)
    integer, allocatable :: stencil_right_snd(:,:)
    integer, allocatable :: stencil_right_rcv(:,:)
    integer, allocatable :: tipo_spedito(:,:,:)

    !     to send the solid cells
    integer :: numsolid_left_snd,numsolid_right_snd
    integer :: numsolid_left_rcv,numsolid_right_rcv
    integer, allocatable :: solid_left_snd(:,:)
    integer, allocatable :: solid_left_rcv(:,:)
    integer, allocatable :: solid_right_snd(:,:)
    integer, allocatable :: solid_right_rcv(:,:)
    ! aggiunta giulia per particelle tar
    !     initial scalar value at IB cells
    real, allocatable :: scalar_bottom_ib(:,:)

    real, allocatable :: r_solid(:,:)

    integer :: ipressione_ibm

    ! ----------------------------------------------------
    ! For particles

    ! Identifier for the solid object the IB belongs to
    integer,allocatable :: solidIndex(:,:,:,:)
    ! Identifier for the solid object the IB belongs to
    integer,allocatable :: solidIndexSize(:,:,:)
    ! IP coordinates
    real,allocatable :: nodo_vicino_x_array(:,:,:), nodo_vicino_y_array(:,:,:), nodo_vicino_z_array(:,:,:)
    ! move this to the subroutine in ricerca
    real,allocatable :: normalVectorX(:,:,:), normalVectorY(:,:,:), normalVectorZ(:,:,:)
    real,allocatable :: surfaceDistance(:,:,:)


!***********************************************************************
end module myarrays_ibm
!***********************************************************************
