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
    integer,bind(C) :: niter
    integer i_rest
    !     i_respa         startpa.h
    !     i_print         print.h
    !      integer i_printfile
    !      integer iformat_newres
    !      integer iformat_grid
    !      character*12 string_newres_format
    !      character*12 string_grid_format
    integer inf

    !     characteristics number
    !     re                scala3
    !      real pran
    real rich

    !     time step
    integer ind_cou
    real cou
    !      real dt     scala3
    integer espl

    !     numerical scheme and equa
    integer attiva_scal
    !      integer potenziale    scala3
    integer insc

    !     pressure (poisson eq.)
    real bbx
    real bby
    real bbz
    real eps
    integer ficycle
    integer bodypressure
    integer ipress_cart
    integer freesurface

    !     boundary condition, see a
    integer lett
    integer ibb
    !     bodyforce	          myarrays_bodyforce
    !     num_iter	          myarrays_bodyforce
    integer coef_wall
    integer integrale
    !     rough 	          myarrays_wallmodel
    !     Z0                  myarrays_wallmodel
    !     att_wm_sgs	  myarrays_wallmodel

    !     turbulence model
    integer nsgs
    integer inmod
    integer inmodrho
    integer isotropo
    !     cost	           turbo2_data
    !     costH	           turbo2_data
    !     costV	           turbo2_data
     
    !     forcing
    integer indm
    integer windyes
    integer wavebk
    !     alpha      myarrays_WB
    !     c10	 myarrays_WB
    !     l_0	 myarrays_WB
    integer langyes
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
    integer visualizzo
    integer lagr
    integer i_sta

    !     movie
    integer imovieprint
    integer ktime_movie
    real dt_movie
    real dt_delay
    integer i_movie
    integer j_movie
    integer k_movie

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

    integer npiani
    integer nsonde
    integer, allocatable :: sonde(:,:),piani(:)

    !-----------------------------------------------------------------------
    !     AFILTRAGGIO
    !     read area for filtering
            
    integer ifiltro
    integer nfiltro
      
    integer filtrou
    integer filtrov
    integer filtrow
    integer filtrorho
    integer filtrofi

    integer xstart
    integer xend
    integer ystart
    integer yend
    integer zstart
    integer zend
           

!***********************************************************************
end module mysettings
!***********************************************************************
