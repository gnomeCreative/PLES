!***********************************************************************
module mysettings
    !***********************************************************************
    !     contains settings for the simulation not already define in other
    !     module or include file.h
    !
    use iso_c_binding

    implicit none

    real(kind=c_double),bind(C) :: alx,aly,alz

    integer(kind=c_int),bind(C) :: niter,i_rest,inf
    real(kind=c_double),bind(C) :: rich
    real(kind=c_double),bind(C) :: latitude

    ! time step
    integer(kind=c_int),bind(C) :: ind_cou,espl
    real(kind=c_double),bind(C) :: cou

    ! numerical scheme and equa
    integer(kind=c_int),bind(C) :: attiva_scal,insc

    ! pressure (poisson eq.)
    real(kind=c_double),bind(C) :: bbx,bby,bbz,eps
    integer(kind=c_int),bind(C) :: ficycle,freesurface

    ! boundary condition
    integer(kind=c_int),bind(C) :: lett,ibb,coef_wall,integrale

    ! turbulence model
    integer(kind=c_int),bind(C) :: nsgs,inmod,inmodrho,isotropo
    real(kind=c_double),bind(C) :: cost,costH,costV
     
    ! forcing
    integer(kind=c_int),bind(C) :: windyes,wavebk,langyes

    ! other
    integer(kind=c_int),bind(C) :: visualizzo,lagr,i_sta

    ! movie
    integer(kind=c_int),bind(C) :: imovieprint,ktime_movie
    real(kind=c_double),bind(C) :: dt_movie,dt_delay
    integer(kind=c_int),bind(C) :: i_movie,j_movie,k_movie

end module mysettings
!***********************************************************************
