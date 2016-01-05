!***********************************************************************
      module myarrays_ibm
!***********************************************************************
!     array for immersed boundary method
!-----------------------------------------------------------------------
!
      integer bodyforce,num_iter
      integer MN,MP
      integer numero_celle_IB,numero_celle_bloccate
      integer num_ib,num_solide 

      integer,allocatable :: indici_CELLE_IB(:,:)
      integer,allocatable :: indici_celle_bloccate(:,:)

      real,allocatable :: distanze_CELLE_IB(:,:)
      real,allocatable :: dist_pp_ib(:)               
      real,allocatable :: dist_ib_parete(:)
      real,allocatable :: dist_nulle(:)
      real,allocatable :: ustar(:)       
      real,allocatable :: proiezioni(:,:)
      
!     array for rotation with eulerian angles c
      real, allocatable :: rot(:,:,:)
      real, allocatable :: rot_inverse(:,:,:)
      integer, allocatable :: id_ib(:,:,:)
      integer, allocatable :: id_ib_solide(:,:,:)
      
!     array for trilinear
      real, allocatable :: tricoef(:,:)
      integer, allocatable :: trind(:,:,:)
      
!
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


      integer ipressione_ibm
!***********************************************************************
      end module myarrays_ibm
!***********************************************************************
