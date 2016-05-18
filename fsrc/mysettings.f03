!***********************************************************************
module mysettings
    !***********************************************************************
    !     contains settings for the simulation not already define in other
    !     module or include file.h
    !
    use iso_c_binding

    implicit none

    logical(kind=c_bool),bind(C) :: particles,bodyforce !,bodypressure

    real(kind=c_double),pointer :: pran(:),prsc(:)
    logical(kind=c_bool),bind(C) :: re_analogy

    real(kind=c_double),bind(C) :: alx,aly,alz

    integer(kind=c_int),bind(C) :: niter,i_rest
    logical(kind=c_bool),bind(C) :: inf
    real(kind=c_double),bind(C) :: rich

    ! time step
    integer(kind=c_int),bind(C) :: ind_cou,espl
    real(kind=c_double),bind(C) :: cou

    ! numerical scheme and equa
    integer(kind=c_int),bind(C) :: insc
    logical(kind=c_bool),bind(C) :: attiva_scal

    ! pressure (poisson eq.)
    real(kind=c_double),bind(C) :: bbx,bby,bbz,eps
    integer(kind=c_int),bind(C) :: ficycle
    logical(kind=c_bool),bind(C) :: freesurface

    ! boundary condition
    integer(kind=c_int),bind(C) :: coef_wall
    logical(kind=c_bool),bind(C) :: lett,ibb

    ! turbulence model
    integer(kind=c_int),bind(C) :: nsgs
    logical(kind=c_bool),bind(C) :: inmod,inmodrho,isotropo
    real(kind=c_double),bind(C) :: cost,costH,costV
     
    ! forcing
    integer(kind=c_int),bind(C) :: windyes
    real(kind=c_double),bind(C) :: latitude

end module mysettings
!***********************************************************************
