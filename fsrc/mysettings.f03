!***********************************************************************
module mysettings
    !***********************************************************************
    !     contains settings for the simulation not already define in other
    !     module or include file.h
    !
    use iso_c_binding

    implicit none
    !
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !     AGENERALE
    !
    integer,bind(C) :: niter,i_rest,inf
    !     i_respa         startpa.h
    !     i_print         print.h
    !      integer i_printfile
    !      integer iformat_newres
    !      integer iformat_grid
    !      character*12 string_newres_format
    !      character*12 string_grid_format

    !     characteristics number
    !     re                scala3
    !      real pran
    real,bind(C) :: rich

    !     time step
    integer,bind(C) :: ind_cou,espl
    real,bind(C) :: cou
    !      real dt     scala3

    !     numerical scheme and equa
    integer,bind(C) :: attiva_scal,insc
    !      integer potenziale    scala3

    !     pressure (poisson eq.)
    real,bind(C) :: bbx,bby,bbz,eps
    integer,bind(C) :: ficycle,bodypressure,ipress_cart,freesurface

    !     boundary condition, see a
    integer,bind(C) :: lett,ibb,coef_wall,integrale
    !     bodyforce	          myarrays_bodyforce
    !     num_iter	          myarrays_bodyforce
    !     rough 	          myarrays_wallmodel
    !     Z0                  myarrays_wallmodel
    !     att_wm_sgs	  myarrays_wallmodel

    !     turbulence model
    integer,bind(C) :: nsgs,inmod,inmodrho,isotropo
    !     cost	           turbo2_data
    !     costH	           turbo2_data
    !     costV	           turbo2_data
     
    !     forcing
    integer,bind(C) :: indm,windyes,wavebk,langyes
    !     alpha      myarrays_WB
    !     c10	 myarrays_WB
    !     l_0	 myarrays_WB
    !     lamb       myarrays_LC
    !     h_0        myarrays_LC
    !     A1         myarrays_cor
    !     A2	 myarrays_cor
    !     A3	 myarrays_cor
    !     U0	 myarrays_cor
    !     V0	 myarrays_cor
    !     W0	 myarrays_cor
    !     omega1	 myarrays_cor
    !     omega2	 myarrays_cor
    !     omega3	 myarrays_cor
    !     omegaM2	 myarrays_cor

    !     other
    integer,bind(C) :: visualizzo,lagr,i_sta

    !     movie
    integer,bind(C) :: imovieprint,ktime_movie
    real,bind(C) :: dt_movie,dt_delay
    integer,bind(C) :: i_movie,j_movie,k_movie

    !-----------------------------------------------------------------------
    !     ABOUNDARY
    !     read the boundary condition

    !     periodicity in i,j,k
    !     ip        period.h
    !     jp        period.h
    !     kp        period.h
      
    
    !     inflow and outflow conditions
    !      infout1         orl.h
    !      infout2         orl.h
    !      infout3         orl.h
    !      infout4         orl.h
    !      infout5         orl.h
    !      infout6         orl.h
      

    !     wall functions condition
    !     wfp1    myarrays_wallmodel
    !     wfp2    myarrays_wallmodel
    !     wfp3    myarrays_wallmodel
    !     wfp4    myarrays_wallmodel
    !     wfp5    myarrays_wallmodel
    !     wfp6    myarrays_wallmodel

    !-----------------------------------------------------------------------
    !     APIANISONDE
    !     read index for tracer

    integer,bind(C) :: npiani,nsonde
    integer,pointer :: sonde(:,:),piani(:)
    integer,pointer :: sondeindexi(:),sondeindexj(:),sondeindexk(:)
    type(c_ptr),bind(C) :: c_sondeindexi,c_sondeindexj,c_sondeindexk,c_piani

    !-----------------------------------------------------------------------
    !     AFILTRAGGIO
    !     read area for filtering
            
    integer,bind(C) :: ifiltro,nfiltro
      
    integer,bind(C) :: filtrou,filtrov,filtrow,filtrorho,filtrofi

    integer,bind(C) :: xstart,xend,ystart,yend,zstart,zend
           

!***********************************************************************
end module mysettings
!***********************************************************************
