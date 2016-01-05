!***********************************************************************
subroutine read_simulation_setting
    !***********************************************************************
    ! read the settings for the simulation
    use mysending, only: kparaend, kparasta, myid
    use mysettings, only: piani, sonde, attiva_scal, bbx, bby, bbz, bodypressure, coef_wall, &
        imovieprint, dt_delay, dt_movie, i_movie, i_sta, ibb, ifiltro, &
        ind_cou, indm, inmod, inmodrho, insc, integrale, inf, espl, ficycle, &
        filtrofi, filtrorho, filtrou, filtrov, filtrow, freesurface, &
        cou, eps, i_rest, ipress_cart, isotropo, j_movie, &
        k_movie, ktime_movie, lagr, langyes, lett, nfiltro, niter, &
        npiani, nsgs, nsonde, rich, visualizzo, &
        wavebk, windyes, xend, xstart, yend, ystart, zend, zstart

    use multigrid_module, only: omega, nlevmultimax, islor, jpos

    use mysettings_boundary, only: ibodybuffer1, ibodybuffer2, ibodybuffer3, &
        ibodybuffer4, ibodybuffer5, ibodybuffer6, &
        iboun1, iboun2, iboun3, iboun4, iboun5, iboun6

    use myarrays_LC, only: h_0, lamb
    use myarrays_WB, only: alpha, c10, l_0
    use myarrays_cor, only: A1, A2, A3, LATITUDE, omegaM2, U0, V0, W0
    use myarrays_wallmodel, only: att_wm_sgs, rough, Z0, wfp1, wfp2, wfp3, wfp4, wfp5, wfp6
    use turbo2_data, only: cost, costH, costV
    use myarrays_ibm, only: bodyforce, num_iter
    use myarrays_buffer_bodyforce, only: corr_factor, ispon, kspon
    use myarrays_density, only: pran, prsc,re_analogy
    use myarrays_moisture, only: betaQ, betaT, cpd, Gdry, imoist, Lv, Ma, Mv, Qref, Rd, Tref
    use scala3, only: nscal,n1,n2,n3,dt,potenziale,re
    use period, only: ip,jp,kp
    use print, only: i_print
    use orl, only: infout1, infout2, infout3, infout4, infout5, infout6
    !
    use output_module, only: string_newres_format, i_cumulative, i_printfile, i_time, &
        ifolder, iformat_grid, iformat_newres, i_paraview, &
        print_iter_or_time, string_grid_format, string_newres_format


    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer ipiani,isonde,i_oil_species,oil_x0,oil_z0
    integer i,ichar,iic
    integer count_line,count_line_check,ierr
    character*1 commento
    !-----------------------------------------------------------------------
    !
    !     READ THE FILE AGENERALE.IN
    !
    !     general settings for the simulation
    !     value / variable / comment
    open(13,file='Agenerale.in',status='old')
    !     iteration / restart / print
    do i=1,5
        read(13,*)commento  !comments
    end do
    read(13,*) niter         !iteration number
    read(13,*) i_rest        !restart (0:new copmutation; 1:restart; 2: from interpolation; 3: nesting)
    !      read(13,*) i_respa       !restart for particles (0: start from the begin; 1:restart)
    read(13,*) print_iter_or_time
    read(13,*) i_print       !print the flow field every i_print iteration
    read(13,*) i_time
    read(13,*) i_printfile
    read(13,*) i_paraview
    read(13,*) i_cumulative
    read(13,*) iformat_newres
    read(13,*) iformat_grid
    read(13,*) ifolder
    read(13,*) string_newres_format

    string_newres_format = adjustl(string_newres_format)
    do ichar = 1,len(string_newres_format)
        if(string_newres_format(ichar:ichar)==';')then
            do iic = ichar,len(string_newres_format)
                string_newres_format(iic:iic) = ' '
            end do
        end if
    end do
      
      
      
      
    read(13,*) string_grid_format
    read(13,*) inf           !print flow plane to generate inflow condition (0:no ; 1: yes)

    !     characteristics number
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) re            !Reynolds (1/viscosity)
    !      read(13,*) pran          !Prandtl
    read(13,*) rich          !Richardson ( zero = passive scalar )

    !     time step
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) ind_cou       !(1: constant courant; 0: constant delta t)
    read(13,*) cou           !courant number  (only if ind_cou = 1)
    read(13,*) dt            !time step (only if ind_cou = 0)
    read(13,*) espl          !time advancment (0:semimplicit AB+CN ; 1:explicit AB)

    !     numerical scheme and equation switch on/off
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) attiva_scal   !(0: skip scalar eq.; 1:solve scalar eq.)
    read(13,*) potenziale    !potential flow (0:no ; 1: yes)
    read(13,*) insc          !spatial derivative scheme for convective term (0: central difference; 1: QUICK)

    !     pressure (poisson eq.)
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) omega        !relaxation factor
    read(13,*) jpos         ! poisson solver (0:sor; 1: sor/slor+multigrid)
    read(13,*) islor         !only for jpos=1, line sor in eta or sor (0:sor; 1:slor)
    read(13,*) nlevmultimax  !multigrid levell from 1 to 4
    read(13,*) bbx           !impose pressure gradient in csi
    read(13,*) bby           !impose pressure gradient in eta
    read(13,*) bbz           !impose pressure gradient in zita
    read(13,*) eps           !residual error for poisson (out from cycle)
    read(13,*) ficycle       !number of max cycle for pressure
    read(13,*) bodypressure  !if bodyforce=1 pressure correction on IBM (0:off; 1:on), only for 1 multigrid level with sor
    read(13,*) ipress_cart
    read(13,*) freesurface   ! free surface at top surface (side4)

    !     boundary condition, see also Aboundary.in
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) lett          !open boundary (0: off; 1: open boundary: 2: open boundary+IBM)
    read(13,*) ibb           !bodyforce on buffe (0:off, 1:on)
    read(13,*) bodyforce     !immersed boundary (0: off; 1: on)
    read(13,*) num_iter      !only if bodyforce=1, number of cycle for iterative cycle on taylor serie
    read(13,*) coef_wall     !wall function (0: off, 1: on), see also Aboundary.in
    read(13,*) integrale     !only if coef_wall=1, option for Werner Wengle (1: integral; 0: point)
    read(13,*) rough         !only if coef_wall=1, smooth or rough surface (0:smooth; 1: roughness)
    read(13,*) Z0            !only if coef_wall=1, roughness height
    read(13,*) att_wm_sgs    !only if coef_wall=1, correction on nu_t for wall function

    !     turbulence model
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) nsgs          !Smagorinsky sgs model (0:dns; 1: standard smagorinsky; 2: dynamic smagorinsky)
    read(13,*) inmod         !scale similar (0:off; 1:on), means mixed model
    read(13,*) inmodrho      !scale similar for scalars (0:off; 1:on)
    read(13,*) isotropo      !anisotropic/isotropic grid (0: two-eddy-viscosity, 1: one-eddy-viscosity)
    read(13,*) cost          !only if nsgs=1, isotropic model constant
    read(13,*) costH         !only if nsgs=1, anisotropic model constant for horizontal direction (1 and 3)
    read(13,*) costV         !only if nsgs=1, anisotropic model constant for vertical direction( 2 )
     
    !     forcing
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) indm          !      (1: constant mass, 0: constant force)
    read(13,*) windyes       !      wind on side 4 (0:off; 1:on)
    read(13,*) wavebk        !      wavebreaking (0:off; 1:on)
    read(13,*) alpha         !
    read(13,*) c10           !
    read(13,*) l_0           !
    read(13,*) langyes       !
    read(13,*) lamb          !
    read(13,*) h_0           !
    read(13,*) A1            !
    read(13,*) A2            !
    read(13,*) A3            !
    read(13,*) U0            !
    read(13,*) V0            !
    read(13,*) W0            !
    read(13,*) LATITUDE
    read(13,*) omegaM2       !      tide

    !     other
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*) visualizzo    !     print on window some info about the computation
    read(13,*) lagr          !     homogeneous output in xz (0: write ,1 no)
    read(13,*) i_sta         !     NOT USED

    !     movie
    do i=1,3
        read(13,*)commento  !comments
    end do
    read(13,*)imovieprint
    read(13,*)ktime_movie
    read(13,*)dt_movie
    read(13,*)dt_delay
    read(13,*)i_movie
    read(13,*)j_movie
    read(13,*)k_movie

    close(13)
      
      
    !     check
    if(imovieprint==1 .and. i_movie .gt. n1)then
        if(myid==0)then
            write(*,*)'check Agenerale.in'
            write(*,*)'wrong index i_movie'
        end if
        stop
    end if

    if(imovieprint==1 .and. j_movie .gt. n2)then
        if(myid==0)then
            write(*,*)'check Agenerale.in'
            write(*,*)'wrong index j_movie'
        end if
        stop
    end if

    if(imovieprint==1 .and. k_movie .gt. n3)then
        if(myid==0)then
            write(*,*)'check Agenerale.in'
            write(*,*)'wrong index k_movie'
        end if
        stop
    end if
      
    !-----------------------------------------------------------------------
 
    !     read the boundary condition
    open(15,file='Aboundary.in',status='old')
      
    do i=1,5
        read(15,*)commento  !comments
    end do
    !     periodicity in i,j,k
    read(15,*)ip
    read(15,*)jp
    read(15,*)kp

    !     inflow and outflow conditions
    do i=1,3
        read(15,*)commento  !comments
    end do
    read(15,*) infout1  ! for each side  0:inflow ; 1:outflow ; 2:wall
    read(15,*) infout2
    read(15,*) infout3
    read(15,*) infout4
    read(15,*) infout5
    read(15,*) infout6

    !     b.c. at the wall condition
    do i=1,3
        read(15,*)commento  !comments
    end do
    read(15,*) iboun1
    read(15,*) iboun2
    read(15,*) iboun3
    read(15,*) iboun4
    read(15,*) iboun5
    read(15,*) iboun6


    !     wall functions condition
    do i=1,3
        read(15,*)commento  !comments
    end do
    read(15,*) wfp1
    read(15,*) wfp2
    read(15,*) wfp3
    read(15,*) wfp4
    read(15,*) wfp5
    read(15,*) wfp6

    !     disturbance on buffer
    do i=1,3
        read(15,*)commento  !comments
    end do
    read(15,*) ibodybuffer1
    read(15,*) ibodybuffer2
    read(15,*) ibodybuffer3
    read(15,*) ibodybuffer4
    read(15,*) ibodybuffer5
    read(15,*) ibodybuffer6
      
      
    if(ibb == 0)then
        ibodybuffer1 = 0
        ibodybuffer2 = 0
        ibodybuffer3 = 0
        ibodybuffer4 = 0
        ibodybuffer5 = 0
        ibodybuffer6 = 0
    end if
            
    !     read sponge dimension for nesting
    do i=1,3
        read(15,*)commento  !comments
    end do
    read(15,*)ispon    !in myarrays_buffer_bodyforce
    read(15,*)kspon
    read(15,*)corr_factor

    close(15)

    !-----------------------------------------------------------------------
    !     Ascalare.in
    allocate(pran(nscal))
    allocate(prsc(nscal))
     
     
    !     check the file Ascalare.in
    open(15,file='Ascalare.in',status='old',iostat=ierr)
    count_line=0
    do while(ierr==0)
        read(15,*,iostat=ierr)commento
        if(commento=='e')ierr = 1
        count_line = count_line + 1
    end do
    count_line = count_line - 1
    !      if(myid==0)write(*,*)'COUNT LINE',count_line
    close(15)
    count_line_check = 6 + nscal + nscal
    if(count_line < count_line_check)then
        if(myid==0)then
            write(*,*)'CHECK Ascalare.in too few values'
            write(*,*)'with respect to the number of scalar'
            stop
        end if
    elseif(count_line > count_line_check)then
        if(myid==0)then
            write(*,*)'CHECK Ascalare.in too many values'
            write(*,*)'with respect to the number of scalar'
            stop
        end if
    end if
      
      
    open(15,file='Ascalare.in',status='old')
    do i=1,5
        read(15,*)commento  !comments
    end do
    read(15,*)re_analogy
    do i=1,nscal
        read(15,*)pran(i) !prandtl number
    end do
    do i=1,nscal
        read(15,*)prsc(i) !prandtl sgs
    end do
    read(15,*)commento  !fine scalari
	       
    do i=1,3
        read(15,*)commento  !comments
    end do
    read(15,*)imoist
      
    do i=1,3
        read(15,*)commento  !comments
    end do
    read(15,*)Tref
    read(15,*)Qref
    read(15,*)betaT
    read(15,*)betaQ
    read(15,*)Ma
    read(15,*)Mv
    read(15,*)Lv
    read(15,*)Gdry
    read(15,*)Rd
    read(15,*)cpd

    close(15)
    !-----------------------------------------------------------------------
    !     read index for tracer
    open(16,file='Apianisonde.in',status='old')
    do i=1,5
        read(16,*)commento  !comments
    end do
    !     inflow planes
    read(16,*)npiani
    do i=1,3
        read(16,*)commento  !comments
    end do
      
    if(inf .ne.0)then
        allocate(piani(npiani))
             
        do ipiani=1,npiani
            read(16,*)piani(ipiani)
            if(myid.eq.0 .and. piani(ipiani).gt.n1+1)then
                write(*,*)'piano',ipiani,'out of bounds i'
                stop
            end if
        end do
       
    end if
      
    do ipiani=1,100
        read(16,*)commento
        if(commento .eq. 'F')exit
    end do
    if(myid.eq.0)write(*,*)'plane reading end',commento

    !     tracers
    do i=1,3
        read(16,*)commento  !comments
    end do

    read(16,*)nsonde
      
    allocate(sonde(3,nsonde))
      
    do isonde=1,nsonde
      
        do i=1,3
            read(16,*)commento  !comments
        end do
	 
        read(16,*)sonde(1,isonde)  !index i for tracers
        if(myid.eq.0 .and. sonde(1,isonde).gt.n1+1)then
            write(*,*)'sonda',isonde,'out of bounds for index i'
            stop
        end if
	 
        read(16,*)sonde(2,isonde)  !index j for tracers
        if(myid.eq.0 .and. sonde(2,isonde).gt.n2+1)then
            write(*,*)'sonda',isonde,'out of bounds for index j'
            stop
        end if
	 
        read(16,*)sonde(3,isonde)  !index k for tracers
        if(myid.eq.0 .and. sonde(3,isonde).gt.n3+1)then
            write(*,*)'sonda',isonde,'out of bounds for index k'
            stop
        end if
	       
    end do

    do isonde=1,100
        read(16,*)commento
        if(commento .eq. 'F')exit
    end do
    if(myid.eq.0)write(*,*)'reading tracers end',commento
      
    close(16)

    !-----------------------------------------------------------------------
    !     read area for filtering
    open(14,file='Afiltraggio.in',status='old')
    read(14,*) ifiltro
    read(14,*) nfiltro
      
    read(14,*) filtrou
    read(14,*) filtrov
    read(14,*) filtrow
    read(14,*) filtrorho
    read(14,*) filtrofi

    read(14,*) xstart
    read(14,*) xend
    read(14,*) ystart
    read(14,*) yend
    read(14,*) zstart
    read(14,*) zend
           
    close(14)
      
    if(ifiltro .eq. 1)then
        write(*,*)'PAY ATTENTION YOU ARE FILTERING THE FLOW FIELD'
        write(11,*)'FILTERING ON',ifiltro

      
        if(myid .eq. 0)then
            if(xend.gt.n1 .or. yend.gt.n2 .or.zend.gt.n3 &
                .or. xstart.lt.1 .or. ystart.lt.1 .or. zstart.lt.1)then
                write(*,*)'FILTERING AREA TOO LARGE'
                write(*,*)'OUT OF BOUNDS'
                write(*,*)'xend:',xend,',n1:',n1
                write(*,*)'yend:',yend,',n2:',n2
                write(*,*)'zend:',zend,',n3:',n3
                write(*,*)'xend:',xstart,1
                write(*,*)'yend:',ystart,1
                write(*,*)'zend:',zstart,1

                write(11,*)'FILTERING AREA TOO LARGE'
                write(*,*)'OUT OF BOUNDS'
                write(11,*)'xend:',xend,',n1:',n1
                write(11,*)'yend:',yend,',n2:',n2
                write(11,*)'zend:',zend,',n3:',n3
                write(11,*)'xend:',xstart,1
                write(11,*)'yend:',ystart,1
                write(11,*)'zend:',zstart,1
	  
                stop
	  
            end if
        end if

        if(zstart .le. kparasta)then
            zstart=kparasta
        end if
        if(zend .ge. kparaend)then
            zend = kparaend
        endif

        if(zstart .gt. kparaend .or. zend .lt. kparasta)then
            zstart = kparaend
            zend =   kparaend -1
            write(*,*)'proc: ',myid,'no filtering',zstart,zend
            write(11,*)'proc: ',myid,'no filtering',zstart,zend
        else
            write(*,*) 'proc:',myid,'filter in k between',zstart,'and',zend
            write(11,*)'proc:',myid,'filter in k between',zstart,'and',zend
        end if

    else
        write(11,*)'FILTERING OFF',ifiltro
    end if
    !
    !-----------------------------------------------------------------------
     
    return
end
