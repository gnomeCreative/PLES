module wallmodel_module

    !      for wall modeling

    use geometricRoutines
    use iso_c_binding

    private

    !real,allocatable,public :: punto_wfp3(:,:,:,:),punto_wfp4(:,:,:,:)
    real,allocatable,public :: wf_points(:,:,:,:),wf_projections(:,:,:,:),wf_distance(:,:,:)
    ! new version of same vectors
    real,allocatable,public :: ip_wall3(:,:,:),ib_wall3(:,:,:)
    real,allocatable,public :: ip_wall4(:,:,:),ib_wall4(:,:,:)
    real,allocatable,public :: u_t(:,:,:)

    logical,allocatable,public ::att_mod_par(:,:,:)
    real,allocatable,public :: utangente(:,:,:)

    logical(kind=c_bool),public,bind(C) :: wfp1,wfp2,wfp3,wfp4,wfp5,wfp6
    integer,public :: eseguo34
    logical(kind=c_bool),public,bind(C) :: att_wm_sgs
    logical(kind=c_bool),public,bind(C) :: rough
    real(kind=c_double),public,bind(C) :: z0

    ! costants for the law of the wall
    real,parameter :: vonKarman=0.41
    ! 1/k with k being thevon karman constant
    real,parameter :: kdynamic=1/vonKarman
    real,parameter :: coef_rough_const=5.1
    real,parameter :: viscous_threshold=5.0
    real,parameter :: log_threshold=30.0

    ! constants for Werner-Wengle approximated profile
    ! only const3 is actually used
    real,parameter :: wwa=8.3
    real,parameter :: wwb=1.0/7.0
    !real,parameter :: const1=0.5*(1.0-wwb)*wwa**((1.0+wwb)/(1.0-wwb))
    !real,parameter :: const2=(1.0+wwb)/wwa
    real,parameter :: const3=wwa**(2.0/(1.0-wwb))
    !real,parameter :: const4=2.0/(1.0+wwb)
    real,parameter :: const5=1/(1+wwb)
    real,parameter :: ww_viscous_threshold=11.81

    public :: initialize_wallmodel,ustar_lawofwall,speed_lawofwall
    public :: ustar_linprofile,ustar_logprofile
    public :: ustar_wernerwengle,speed_wernerwengle
    public :: correggi_walls


contains

    subroutine initialize_wallmodel(coef_wall)

        use scala3, only: jx
        use mysending, only: kparasta,kparaend

        implicit none

        integer,intent(in) :: coef_wall

        if (coef_wall == 0) then
            wfp1 = .false.
            wfp2 = .false.
            wfp3 = .false.
            wfp4 = .false.
            wfp5 = .false.
            wfp6 = .false.
            att_wm_sgs=.false.
        else if (coef_wall==1) then
            allocate(att_mod_par(1:jx,2,kparasta:kparaend))
            allocate(u_t(1:jx,2,kparasta:kparaend))
            allocate(utangente(1:jx,2,kparasta:kparaend))
            !allocate(punto_wfp3(3,3,jx,kparasta:kparaend))
            !allocate(punto_wfp4(3,3,jx,kparasta:kparaend))
            allocate(wf_points(3,jx,kparasta:kparaend,2))
            allocate(wf_projections(3,jx,kparasta:kparaend,2))
            allocate(wf_distance(jx,kparasta:kparaend,2))

            att_mod_par(:,:,:)=.false.
            utangente(:,:,:)=1.0
            u_t(:,:,:)=1.0

            if (wfp3) then
                att_mod_par(:,1,:)=.true. ! first index=1 : face 3
            end if

            if (wfp4) then
                att_mod_par(:,2,:)=.true. ! second index=2 : face 4
            end if

        end if

  end subroutine initialize_wallmodel


    subroutine correggi_walls(ktime,niter,tipo,i_rest)

        use myarrays_velo3, only: u,v,w
        use myarrays_metri3, only: x,y,z,xcd,ycd,zcd
        !
        use mysending
        use scala3, only: jx,jy,jz,n1,n2,re
        !
        use mpi

        implicit none

        !-----------------------------------------------------------------------
        integer,intent(in) :: ktime,niter,i_rest
        integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
        !-----------------------------------------------------------------------
        real,dimension(3) :: p1,p2,p3
        real,dimension(3) :: speed
        real,dimension(3) :: wf_point,wf_projection
        !real :: wf_distance
        real :: a,b,c,d
        real :: distanza,alfa
        integer :: i,j0,j1,k,l,ierr
        integer :: start_loop,end_loop
        integer,dimension(2) :: wallpoint,centroid
        !-----------------------------------------------------------------------
        att_mod_par(:,:,:)=.true.

        ! define how to execute the loop
        ! being that the array structure is weird, this is a trick to execute the correct indices
        if (wfp3) then
            start_loop=1
            wallpoint(1)=0
            centroid(1)=1
        else
            start_loop=2
        end if
        if (wfp4) then
            end_loop=2
            wallpoint(2)=jy
            centroid(2)=jy
        else
            end_loop=1
        end if

        ! first value for ustar
        if (ktime==1) then
            u_t(:,:,:)=1
        end if

        if (ktime==1) then
            ! chicco questa va bene per cartesiano!!! per curvilineo dovrebbe trovare
            ! la coordinata sulla faccia con la normale al piano !!!!!!
            do l=start_loop,end_loop

                j0=wallpoint(l)
                j1=centroid(l)

                do k=kparasta,kparaend
                    do i=1,jx

                        ! point 1, cell centroid
                        wf_point(1)=xcd(i,j1,k)
                        wf_point(2)=ycd(i,j1,k)
                        wf_point(3)=zcd(i,j1,k)

                        ! triangle at the wall
                        p1=(/x(i,j0,k),y(i,j0,k),z(i,j0,k)/)
                        p2=(/x(i-1,j0,k),y(i-1,j0,k),z(i-1,j0,k)/)
                        p3=(/x(i,j0,k-1),y(i,j0,k-1),z(i,j0,k-1)/)

                        call plane_3points(p1,p2,p3,a,b,c,d)

                        call find_plane_projection(a,b,c,d,wf_point,wf_projection)

                        ! distance between point and projection
                        wf_distance(i,k,l)=norm2(wf_point(:)-wf_projection(:))

                        wf_points(:,i,k,l)=wf_point(:)
                        wf_projections(:,i,k,l)=wf_projection(:)

                    end do
                end do
            end do
        end if


        if (ktime==1.and.i_rest==0) then

            u_t(:,:,:)=1.0
            utangente(:,:,:)=1.0
            att_mod_par(:,:,:)=.false.

        else

            att_mod_par(:,:,:)=.true.

            do l=start_loop,end_loop

                j1=centroid(l)

                do k=kparasta,kparaend
                    do i=1,jx

                        ! solo su celle fluide
                        if (tipo(i,j1,k)==2) then

                            p1(:)=wf_points(:,i,k,l)
                            p2(:)=wf_projections(:,i,k,l)

                            speed=(/u(i,j1,k),v(i,j1,k),w(i,j1,k)/)

                            ! this is a position vector
                            p3(:)=p1(:)+speed(:)

                            call line_angle(p1,p2,p1,p3,alfa)

                            distanza=wf_distance(i,k,l)

                            utangente(i,l,k)=norm2(speed)*sin(alfa)

                            if (utangente(i,l,k)>0.000001) then

                                u_t(i,l,k)=ustar_lawofwall(distanza,utangente(i,l,k))

                                    ! switch off the wall function if y+<11 or if it is not a fluid node
                                    ! in case of ibm
                                !if (u_t(i,l,k)*distanza*re<=11) then
                                !    att_mod_par(i,l,k)=.false.
                                !end if
                            else
                                att_mod_par(i,l,k)=.false.
                            end if
                        else
                            !punto solido
                            u_t(i,l,k)=1.0
                            utangente(i,l,k)=1.0
                            att_mod_par(i,l,k)=.false.
                        end if

                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------
        !             FACE 4
        !-----------------------------------------------------------------------
!
!        if (wfp4.and.ktime==1) then
!            do k=kparasta,kparaend
!                do i=1,jx
!
!                    wf_point(1)=xcd(i,jy,k)
!                    wf_point(2)=ycd(i,jy,k)
!                    wf_point(3)=zcd(i,jy,k)
!
!                    ! triangle at the wall
!                    p1=(/x(i,jy,k),y(i,jy,k),z(i,jy,k)/)
!                    p2=(/x(i-1,jy,k),y(i-1,jy,k),z(i-1,jy,k)/)
!                    p3=(/x(i,jy,k-1),y(i,jy,k-1),z(i,jy,k-1)/)
!
!                    call plane_3points(p1,p2,p3,a,b,c,d)
!
!                    call find_plane_projection(a,b,c,d,wf_point,wf_projection)
!
!                    ! distance between point and projection
!                    wf_distance(i,k,2)=norm2(wf_point(:)-wf_projection(:))
!
!                    wf_points(:,i,k,2)=wf_point(:)
!                   wf_projections(:,i,k,2)=wf_projection(:)
!
!                end do
!            end do
!        end if
!
!
!        if (ktime==1 .and. i_rest==0) then
!            u_t(:,2,:)=1.0
!            utangente(:,2,:)=1.0
!            att_mod_par(:,2,:)=.false.
!        else
!            if (wfp4) then
!                do k=kparasta,kparaend
!                    do i=1,jx
!
!                        if (tipo(i,jy,k)==2) then
!
!                            p1(:)=wf_points(:,i,k,2)
!                            p2(:)=wf_projections(:,i,k,2)
!
!                            speed=(/ u(i,jy,k),v(i,jy,k),w(i,jy,k) /)
!
!                            ! this is a position vector
!                            p3(:)=p1(:)+speed(:)
!
!                            call line_angle(p1,p2,p1,p3,alfa)
!
!                            distanza=wf_distance(i,k,2)
!
!                            utangente(i,2,k)=norm2(speed)*sin(alfa)
!
!                            if (abs(utangente(i,2,k))>0.000001) then
!                                u_t(i,2,k)=ustar_logprofile(distanza,utangente(i,2,k))
!
!                                    ! switch off the wall function if y+<11 or if it is not a fluid node
!                                    ! in case of ibm
!                                if (u_t(i,2,k)*distanza*re <= 11. .or.tipo(i,jy,k).ne.2) then
!                                    att_mod_par(i,2,k)=.false.
!                                end if
!                            else
!                                att_mod_par(i,2,k)=.false.
!                            end if
!
!                        else !punto solido
!                            u_t(i,2,k)=1.0
!                            utangente(i,2,k)=1.0
!                            att_mod_par(i,2,k)=.false.
!                        end if
!                    end do
!                end do
!            else
!                do k=kparasta,kparaend
!                    do i=1,jx
!                        att_mod_par(i,2,k)=.false.
!                    end do
!                end do
!            end if
!
!        end if

        return

    end subroutine correggi_walls


    function ustar_linprofile(distance,speed) result (ustar)

        use scala3, only: re

        ! simple linear profile

        implicit none

        ! output: shear velocity
        real :: ustar
        !---------------------------------------------------------------
        real,intent(in) :: distance
        real,intent(in) :: speed
        !---------------------------------------------------------------
        real :: visco
        !---------------------------------------------------------------

        visco=1.0/re

        ustar=sqrt(speed*visco/distance)

    end function ustar_linprofile

    function ustar_logprofile(distance,speed,ustar_old) result (ustar)

        use scala3, only: re

        implicit none

        ! output is shear velocity of log profile
        real :: ustar
        ! known point of the profile
        real,intent(in) :: distance,speed
        ! value at previous time step (if known)
        real,optional,intent(in) :: ustar_old

        ! ----------------------------------------------------------------------
        ! stuff for Newton-Raphson
        real :: f,fprime
        real :: errore_star
        real :: distance_plus,speed_plus
        integer :: contatore_star
        integer,parameter :: max_iter=100
        real,parameter :: toll=1.e-6
        real,parameter :: init_error=1.e-1
        real,parameter :: ustar_toll=1.e-4


        ! if the wall is smooth, an internal subcycle begins
        contatore_star=0
        errore_star=init_error

        ! first value from input, or Werner-Wengle
        if (present(ustar_old)) then
            ustar=ustar_old
        else
            ustar=ustar_wernerwengle(distance,speed)
        end if

        ! -------------------------------------------------------
        ! Newton - Raphson iterative procedure
        do while (errore_star>toll .and. contatore_star<max_iter) !errore_star>1.d-6
            !write(*,*) 'cont=',contatore_star

            !ustar_log=abs(ustar_log)

            !if (ustar_new<1.e-8) then !ustar(l).lt.1.d-8
            !    caso=0
            !    errore_star=1.e-10 !errore_star = 1.d-10
            !    exit
            !end if

            distance_plus=distance*ustar*re
            speed_plus=speed_plus_logprofile(distance_plus)
            ! since transition is not implemented, we use a simplified version
            ! logarithm (so that it's computed only once)
            f=ustar*speed_plus-speed

            ! fprime is the derivative of f (f')
            fprime=speed_plus+kdynamic!/re*distance

            ustar=ustar-f/fprime

            if (ustar>ustar_toll) then
                errore_star=abs(f/sqrt(ustar))
            end if


            ! older version -------
            ! argomentolog = 1.0 + ustar(l)*rougheight*re
            ! coef_rough = coef_rough - real(transition)*kdynamic*log(argomentolog)
            ! f=ustar(l)*(kdynamic*log(abs(distanza*ustar(l)*re))+coef_rough)-vtan
            ! ! fprime is the derivative of f'
            ! fprime = kdynamic*(log(abs(distanza*ustar(l)*re)))+kdynamic+coef_rough &
            !   -real(transition)*ustar(l)*kdynamic*rougheight/(1./re + rougheight*ustar(l))
            ! ustar(l) = ustar(l) - f/fprime
            ! if (ustar(l)>1.E-8) then
            !    errore_star = abs(f/sqrt(ustar(l)))
            ! end if
            ! coef_rough = coef_rough_const

            contatore_star=contatore_star+1

        end do   !end for do while#

    end function ustar_logprofile

    function ustar_lawofwall(distance,speed,ustar_old) result (ustar)

        implicit none

        ! output is shear velocity of log profile
        real :: ustar
        ! known point of the profile
        real,intent(in) :: distance,speed
        ! value at previous time step
        real,optional,intent(in) :: ustar_old

        real :: ustar_linear,ustar_log
        real :: right_hand_side
        ! ----------------------------------------------------------------------

        if (rough) then

            ! if the wall is rough, compute directly
            right_hand_side=kdynamic*log(distance/z0)
            ustar=abs(speed/right_hand_side)

        else

            if (present(ustar_old)) then
                ustar_log=ustar_logprofile(distance,speed,ustar_old)
            else
                ustar_log=ustar_logprofile(distance,speed)
            end if

            ustar_linear=ustar_linprofile(distance,speed)

            ustar=max(ustar_log,ustar_linear)

        end if


    end function ustar_lawofwall

    function ustar_wernerwengle(distance,speed) result (ustar)

        !-----------------------------------------------------------------------
        ! compute wall shear stress with Werner-Wengle power law (1/7)
        ! for reference, see:
        ! http://aerojet.engr.ucdavis.edu/fluenthelp/html/ug/node516.htm
        ! http://www.tfd.chalmers.se/~lada/postscript_files/vlad_thesis.pdf
        ! but attention, read note at end of subroutine
        !---------------------------------------------------------------
        ! wall function with Werner-Wengle:
        ! if y<11.8       u+ = y+              (linear)
        !    y>11.8       u+ = 8.3*y+**1/7     (exponential)
        !---------------------------------------------------------------

        use scala3, only: re

        implicit none

        ! output: shear velocity
        real :: ustar
        !---------------------------------------------------------------
        real,intent(in) :: distance
        real,intent(in) :: speed

        real :: visco
        real :: u_limit

        visco=1.0/re

        ! threshold velocity for viscous sublayer
        u_limit=visco/distance*const3

        if (speed<u_limit) then
            ! linear profile
            ustar=ustar_linprofile(distance,speed)
        else
            ! exponential profile
            ustar=(speed/wwa*(visco/distance)**wwb)**const5
        end if

        ! this version is the one contained in the references
        ! it works for difite volumes, where the velocity is the average velocity in the cells
        ! in finite differences, velocity is in the centerpoint, so no need to follow this complicated procedure
        ! threshold velocity for viscous sublayer
        !delta_y=2.0*distance
        !ratio=visco/delta_y
        !u_limit=visco/(2.0*delta_y)*const3
        !if (speed<u_limit) then
        !    ! linear profile
        !    tau_w=2.0*visco*speed/delta_y
        !else
        !    ! exponential profile
        !    tau_w=(const1*ratio**(1.0+wwb)+const2*speed*ratio**wwb)**const4
        !end if


    end function ustar_wernerwengle

    function speed_wernerwengle(distance,ustar) result (speed)

        ! computes velocity at a distance from the wall using the Werner-Wengle
        ! for references, see function ustar_wernerwengle

        !-----------------------------------------------------------------------
        ! compute wall shear stress with Werner-Wengle power law (1/7)
        ! for reference, see:
        ! http://aerojet.engr.ucdavis.edu/fluenthelp/html/ug/node516.htm
        ! http://www.tfd.chalmers.se/~lada/postscript_files/vlad_thesis.pdf
        !
        !---------------------------------------------------------------

        use scala3, only: re

        implicit none

        real :: speed
        ! ------------------------------------------------------
        real,intent(in) :: distance,ustar
        ! ----------------------------------------------------

        real :: y_plus

        ! compute y+ of point
        y_plus=distance*ustar*re

        ! check position in the boundary layer
        if (y_plus<=ww_viscous_threshold) then
            ! viscous layer
            speed=y_plus*ustar
        else
            ! exponential approximation
            speed=ustar*wwa*y_plus**wwb
        end if

    end function speed_wernerwengle

    function speed_lawofwall(distance,ustar) result (speed)

        ! computes velocity in wall units for the whole boundary layer
        ! y+ < 5 -> viscous sublayer (linear profile)
        ! 5 < y+ < 30 -> buffer layer (smooth)
        ! y+ > 30 -> log layer (logarithmic profile)

        use scala3, only: re

        implicit none

        real :: speed
        ! ------------------------------------------------------
        real,intent(in) :: distance,ustar
        ! ----------------------------------------------------
        real :: distance_plus,speed_plus

        ! compute y+ of point
        distance_plus=distance*ustar*re

        ! check position in the limit layer
        if (distance_plus<=viscous_threshold) then
            ! viscous layer
            speed_plus=distance_plus
            !speed_pluswrite(*,*) 'viscous',speed_plus
        else if (distance_plus>=log_threshold) then
            ! log layer
            speed_plus=speed_plus_logprofile(distance_plus)
            !write(*,*) 'log',speed_plus
            ! old:
            !vtan_ib = (kdynamic*log(yib_plus)+coef_rough)*ustar(l)
        else
            !buffer layer
            speed_plus=smooth_ibm(distance_plus)
            !write(*,*) 'buffer',speed_plus
            ! old:
            !call smooth_ibm(yib_plus,vtan_ib,kdynamic,coef_rough)
            !speed=speed
        end if

        speed=speed_plus*ustar

    end function speed_lawofwall

    function speed_plus_logprofile(distance_plus) result (speed_plus)

        ! speed accordin gto logarithm law, expressed in wall units

        use scala3, only: re

        implicit none

        real :: speed_plus
        ! ------------------------------------------------------
        real,intent(in) :: distance_plus

        speed_plus=kdynamic*log(distance_plus)+coef_rough_const

    end function speed_plus_logprofile

    function smooth_ibm(distance_plus) result(speed_plus)

        ! computes velocity in wall units for buffer layer
        ! 5 < y+ < 30

        implicit none

        !-----------------------------------------------------------------------
        real :: speed_plus
        !-----------------------------------------------------------------------
        real,intent(in) :: distance_plus
        !-----------------------------------------------------------------------
        integer,parameter :: n=2
        integer ibcbeg, ibcend
        real :: ybcbeg,ybcend
        real,dimension(2) :: t(2)
        real,dimension(2) :: ubuf(2),ypp(2)
        real :: tval, yval, ypval, yppval
        !-----------------------------------------------------------------------

        t(1) = viscous_threshold
        t(2) = log_threshold

        ubuf(1) = viscous_threshold
        ubuf(2) = speed_plus_logprofile(log_threshold)

        ibcbeg = 1
        ybcbeg = 1. !this is the derivative on d u+/d y+

        ibcend = 1
        ybcend = kdynamic * (1./distance_plus) ! du+/dy+ for log law

        !-----------------------------------------------------------------------

        call spline_cubic_set(n,t,ubuf,ibcbeg,ybcbeg,ibcend,ybcend,ypp)

        tval = distance_plus

        call spline_cubic_val(n,t,ubuf,ypp,tval,yval,ypval,yppval)

        speed_plus = yval

        return
    end function smooth_ibm

end module wallmodel_module
