!***********************************************************************
module myarrays_density
    !***********************************************************************
    !     eddy diffusivity
    !-----------------------------------------------------------------------
    !
    use,intrinsic ::  iso_c_binding

    real,allocatable :: akapt(:,:,:,:),akaptV(:,:,:,:)
    real,allocatable :: akapt_piano(:,:,:,:),akaptV_piano(:,:,:,:)
    real,pointer :: pran(:),prsc(:)
    type(c_ptr),bind(C) :: c_pran,c_prsc
    integer,bind(c) :: re_analogy
!
!***********************************************************************
end module myarrays_density
!***********************************************************************
