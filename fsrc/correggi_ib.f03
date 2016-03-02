!***********************************************************************
subroutine correggi_ib(ktime,tipo)


    !***********************************************************************
    !
    ! Roman et al. 2008 Computer and Fluids
    ! Roman et al. 2009 Physics of Fluids
    !
    ! subroutine to correct velocity with Immersed Boundary Method (IBM)
    ! correction to zero for velocity and pressure inside the body
    !
    ! Interpolation on PP point with Taylor series
    !
    ! IB velocity reconstruction with linear or log profile
    !-----------------------------------------------------------------------
    ! used to be:
    ! correggi_ib(ktime,coef_wall,visualizzo,integrale,correggo_rho,correggo_delu,delrhov,tipo,z0,ti,nfinale)

    use myarrays_WB
    use myarrays_ibm
    use mysending
    use myarrays_velo3
    use myarrays_metri3
    !
    use scala3
    use period
    use tipologia
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    ! input stuff
    integer,intent(in) :: ktime
    integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

    !-----------------------------------------------------------------------
    ! local variables
    integer :: kiter
    integer :: i0,j0,k0
    integer :: i,j,k,l
    integer :: coincide_check !,caso,contatore_star
    logical :: node_is_ok
                               
    ! check for errors
    integer :: ierr
    real :: errore,errore_max,errore_max_loc
    real :: variazione1,variazione2,variazione3

    ! velocity at ib points
    real :: uib,vib,wib
    real,allocatable :: u_ib(:),v_ib(:),w_ib(:)

    ! for output
    real :: tautot,tautot_loc
     
    integer :: solidIndexHere,p,index_here
    real :: ib_counter
    real :: x_force_here,y_force_here,z_force_here

    !-----------------------------------------------------------------------
    ! check call subroutine

    if (myid==0) then
        write(*,*)'-----------------------------------------------'
        write(*,*)'       IMMERSED BOUNDARY'

        if (ktime==1) then
            write(*,*)'IB proc0',num_ib,'su',numero_celle_IB
            write(*,*)'solide proc0',num_solide,'on',numero_celle_bloccate
        end if
    end if

    allocate(u_ib(num_ib),v_ib(num_ib),w_ib(num_ib))

    !-----------------------------------------------------------------------
    !     solid border
    do j=1,jy
        do i=1,jx
            k = kparasta -1
            if(tipo(i,j,k)==0)then
                u(i,j,k)=0.
                v(i,j,k)=0.
                w(i,j,k)=0.
            end if
            k = kparaend +1
            if(tipo(i,j,k)==0)then
                u(i,j,k)=0.
                v(i,j,k)=0.
                w(i,j,k)=0.
            end if
        end do
    end do

    !-----------------------------------------------------------------------
    ! correction to solid for IB without V
    do l=1,num_ib

        ! index ib
        i0=indici_CELLE_IB(l,1)
        j0=indici_CELLE_IB(l,2)
        k0=indici_CELLE_IB(l,3)

        ! index V
        i=indici_CELLE_IB(l,4)
        j=indici_CELLE_IB(l,5)
        k=indici_CELLE_IB(l,6)

        coincide_check=(i0-i)*(i0-i)+(j0-j)*(j0-j)+(k0-k)*(k0-k)

        if (coincide_check==0) then
            ! node_is_ok==.false.
            u(i0,j0,k0)=0.0
            v(i0,j0,k0)=0.0
            w(i0,j0,k0)=0.0
        end if

    end do

    !-----------------------------------------------------------------------
    ! CORRECTION ON IB
    !
    ! initialize
    do l=1,num_ib

        i0=indici_CELLE_IB(l,1)
        j0=indici_CELLE_IB(l,2)
        k0=indici_CELLE_IB(l,3)
           
        u_ib(l)=u(i0,j0,k0)
        v_ib(l)=v(i0,j0,k0)
        w_ib(l)=w(i0,j0,k0)

    end do

    !-----------------------------------------------------------------------
    ! iterative procedure, necessary if the stencil has IB nodes
    !
    ! chicco mettere errore adimensionale !!!!
    errore_max=1.1e-7 !-5 !1.1d-7
    kiter=1

    ! do while(errore_max.gt.1.0d-7.and.kiter.le.num_iter)
    do while (kiter<num_iter)

        errore_max=1.0e-7 !-5 !1.0d-7
        errore_max_loc=errore_max
        errore=errore_max

        ! exchange planes between procs:
        ! - to compute derivative for Taylor
        ! - to coorect IB values in the iterative procedure
        do l=1,num_ib

            i0 = indici_CELLE_IB(l,1)
            j0 = indici_CELLE_IB(l,2)
            k0 = indici_CELLE_IB(l,3)

            u(i0,j0,k0) = u_ib(l)
            v(i0,j0,k0) = v_ib(l)
            w(i0,j0,k0) = w_ib(l)

        end do

        call passo_ibm()

        do l=1,num_ib

            ! index IB
            i0 = indici_CELLE_IB(l,1)
            j0 = indici_CELLE_IB(l,2)
            k0 = indici_CELLE_IB(l,3)

            ! index V
            i = indici_CELLE_IB(l,4)
            j = indici_CELLE_IB(l,5)
            k = indici_CELLE_IB(l,6)

            ! -------------------------------------------------------
            ! if IB coincide with V, IB node off
            coincide_check = (i0-i)*(i0-i)+(j0-j)*(j0-j)+(k0-k)*(k0-k)

            if (coincide_check/=0) then
                node_is_ok=.true.
            else
                node_is_ok=.false.
            end if

            ! giulia aggiungo questo controllo
            !            if (.not.node_is_ok .and. ktime==1 .and. myid==0) then
            !                write(*,*) 'IB=V in',i0,j0,k0,coincide_check
            !            end if


            if (node_is_ok) then !only for IB /= V

                errore=0.0

                call compute_u_ib(kiter,ktime,uib,vib,wib,l,i0,j0,k0)

                ! chicco mettere errore adimensionale !!!!
                ! error at ib
                variazione1=uib-u(i0,j0,k0)
                variazione2=vib-v(i0,j0,k0)
                variazione3=wib-w(i0,j0,k0)

                errore=max(abs(variazione1),abs(variazione2),abs(variazione3))

                ! update velocity at IB
                u_ib(l)=variazione1+u(i0,j0,k0)
                v_ib(l)=variazione2+v(i0,j0,k0)
                w_ib(l)=variazione3+w(i0,j0,k0)

                ! update error
                if (errore>=errore_max) then
                    errore_max=errore
                    errore_max_loc=errore
                end if

            end if

        end do !fine loop su celle ib

        ! ----------------------------------------------------------------

        call MPI_ALLREDUCE(errore_max_loc,errore_max,1,MPI_REAL_SD,MPI_MAX,MPI_COMM_WORLD,ierr)

        if(myid.eq.0)then
            write(*,*)'***',kiter,'errore ciclo IBM',errore_max,'***'
        end if

        kiter=kiter+1

    end do ! end loop correzione

    !-----------------------------------------------------------------------
    ! correction for solid cell
    ! giulia aggiungo aggiungo come in 3       if(ktime==1)then .....
    !    if(ktime==1)then
    !        allocate(r_solid(nscal,num_solide))
    !        do l=1,num_solide
    !            i=indici_celle_bloccate(l,1)
    !            j=indici_celle_bloccate(l,2)
    !            k=indici_celle_bloccate(l,3)
    !            do isc=1,nscal
    !                r_solid(isc,l)=rhov(isc,i,j,k)
    !            end do
    !        end do
    !    end if

    do l=1,num_solide

        i=indici_celle_bloccate(l,1)
        j=indici_celle_bloccate(l,2)
        k=indici_celle_bloccate(l,3)

        u(i,j,k)=0.0
        v(i,j,k)=0.0
        w(i,j,k)=0.0
        ! giulia
        if(i.eq.1)then
            u(i-1,j,k)=0.
            v(i-1,j,k)=0.
            w(i-1,j,k)=0.
        end if
        if(i.eq.jx)then
            u(i+1,j,k)=0.
            v(i+1,j,k)=0.
            w(i+1,j,k)=0.
        end if
        if(j.eq.1)then
            u(i,j-1,k)=0.
            v(i,j-1,k)=0.
            w(i,j-1,k)=0.
        end if
        if(j.eq.jy)then
            u(i,j+1,k)=0.
            v(i,j+1,k)=0.
            w(i,j+1,k)=0.
        end if
        if(myid.eq.0.and.k.eq.1)then
            u(i,j,k-1)=0.
            v(i,j,k-1)=0.
            w(i,j,k-1)=0.
        end if
        if(myid.eq.(nproc-1).and.k.eq.jz)then
            u(i,j,k+1)=0.
            v(i,j,k+1)=0.
            w(i,j,k+1)=0.
        end if

    !        fi(i,j,k)=0.
    end do
            
    !-----------------------------------------------------------------------
    !     solid border
    do j=1,jy
        do i=1,jx
            k = kparasta -1
            if(tipo(i,j,k)==0)then
                u(i,j,k)=0.
                v(i,j,k)=0.
                w(i,j,k)=0.
            end if
            k = kparaend +1
            if(tipo(i,j,k)==0)then
                u(i,j,k)=0.
                v(i,j,k)=0.
                w(i,j,k)=0.
            end if
        end do
    end do
    !
    !-----------------------------------------------------------------------
    !     pass data
    call passo_ibm()

    !-----------------------------------------------------------------------
    !     COMPUTE PRESSURE ON IMMERSED BODY
    if (ipressione_ibm==1) then

        ! WHAT ABOUT THIS??? :

    !        do l=1,num_ib
    !
    !            !        index IB
    !            i0 = indici_CELLE_IB(l,1)
    !            j0 = indici_CELLE_IB(l,2)
    !            k0 = indici_CELLE_IB(l,3)
    !
    !            !        index V
    !            i = indici_CELLE_IB(l,4)
    !            j = indici_CELLE_IB(l,5)
    !            k = indici_CELLE_IB(l,6)
    !
    !            fip= tricoef(l,1)*fi(trind(l,1,1),trind(l,1,2),trind(l,1,3)) &
    !                +tricoef(l,2)*fi(trind(l,2,1),trind(l,2,2),trind(l,2,3)) &
    !                +tricoef(l,3)*fi(trind(l,3,1),trind(l,3,2),trind(l,3,3)) &
    !                +tricoef(l,4)*fi(trind(l,4,1),trind(l,4,2),trind(l,4,3))
    !
    !            fib = fi(i0,j0,k0)
    !
    !            fi_pro =  ((dist_ib_parete(l)+dist_pp_ib(l))*fib &
    !                - fip*dist_ib_parete(l)) / dist_pp_ib(l)
    !
    !        end do



       !end if
      
       !if(ktime==nfinale)close(5000+myid)
       !-----------------------------------------------------------------------
       !     print the projection point for pressure on ibm ?? ALE: don't understand
    !      if(ktime==1)then
    !
    !         write(char_myid,'(i6)')myid
    !         char_myid = adjustl(char_myid)
    !
    !         filename='fi_punti'//char_myid(1:len_trim(char_myid))//'.dat'
    !         open(5300+myid,file=filename,status='unknown')
    !
    !         do l=1,num_ib
    !            write(5300+myid,'(3e18.10)')proiezioni(l,1), &
    !               proiezioni(l,2), &
    !               proiezioni(l,3)
    !         end do
    !
    !         close(5300+myid)
    !      end if

    end if
    !-----------------------------------------------------------------------
    !     TAU TOTAL
    tautot_loc=0.
    do l=1,num_ib
        tautot_loc = tautot_loc + ustar(l)*ustar(l)
    end do

    call MPI_ALLREDUCE(tautot_loc,tautot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

    tautot = tautot/(2.*real(jx*jz))

    if(myid.eq.0)then
        write(*,*)'IBM stress tot', tautot
        write(*,*)'IBM u star',sqrt(tautot)
        write(*,*)'-----------------------------------------------'
    end if


    ! Force on particles
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !write(*,*) 'Proc = ',myid,' num_ib=',num_ib
    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    do p=1,totParticles

        x_force_here=0.0
        y_force_here=0.0
        z_force_here=0.0

        index_here=sphereIndex(p)

        ib_counter=0.0

        do l=1,num_ib

            i0=indici_CELLE_IB(l,1)
            j0=indici_CELLE_IB(l,2)
            k0=indici_CELLE_IB(l,3)

            solidIndexHere=solidIndex(i0,j0,k0,1)

            if (solidIndexHere==index_here) then

                x_force_here=y_force_here+shear_ib(l,1)
                y_force_here=y_force_here+shear_ib(l,2)
                z_force_here=z_force_here+shear_ib(l,3)
                ib_counter=ib_counter+1.0

            end if

        end do

        if (ib_counter>0.001) then
            x_force_here=x_force_here/ib_counter
            y_force_here=y_force_here/ib_counter
            z_force_here=z_force_here/ib_counter
        end if

        call MPI_REDUCE(x_force_here,x_force_sphere(p),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(y_force_here,y_force_sphere(p),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(z_force_here,z_force_sphere(p),1,MPI_REAL_SD,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    end do


    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (myid==0) then
        do p=1,totParticles
            x_force_sphere(p)=x_force_sphere(p)*surface_sphere(p)
            y_force_sphere(p)=y_force_sphere(p)*surface_sphere(p)
            z_force_sphere(p)=z_force_sphere(p)*surface_sphere(p)
            write(*,*) '--> Force = (',x_force_sphere(p),',',y_force_sphere(p),',',z_force_sphere(p),') surf = ',surface_sphere(p)
        end do
    end if

    !-----------------------------------------------------------------------
    deallocate(u_ib)
    deallocate(v_ib)
    deallocate(w_ib)
    !-----------------------------------------------------------------------

    return
end

subroutine compute_u_ib(kiter,ktime,uib,vib,wib,l,i0,j0,k0)

    use myarrays_ibm
    use myarrays_velo3
    use myarrays_metri3

    use mysettings, only: coef_wall
    use scala3, only: re
    use myarrays_wallmodel, only: z0

    implicit none

    integer,intent(in) :: l,i0,j0,k0
    integer,intent(in) :: kiter,ktime
    real,intent(inout) :: uib,vib,wib

    !-----------------------------------------------------------------------

    ! local parameters
    real, parameter :: kdynamic=2.44        ! 1/k with k von karman constant
    integer, parameter :: modelloparete = 2   ! 1 wernerwengle, 2 log
    real, parameter :: rough = 0
    integer, parameter :: transition = 0
    real,parameter :: coef_rough_const = 5.1

    !-----------------------------------------------------------------------

    real :: rougheight
    real :: distanza,f,fprime,fattore_distanza
    real :: uip,vip,wip,fip,fib,fi_pro
    integer :: caso,contatore_star
    real :: errore_star
    real :: pro1,pro2,pro3
    real :: local_shear1,local_shear2
    real :: right_hand_side
    real :: u1,u2,u3
    real :: ypp,yib_plus,u_tau_linear
    real :: vtangente1_pp,vtangente2_pp,vnormale_pp
    real :: vtangente1_ib,vtangente2_ib,vnormale_ib
    real :: vtan,vtan_ib
    real :: xib,yib,zib

    ! ---------------------------------------------------------------------

    rougheight = z0

    ! distanza IP-PP
    distanza = dist_ib_parete(l)+dist_pp_ib(l)
    fattore_distanza =  dist_ib_parete(l)/distanza

    ! velocity at interpolation points (V?)
    uip=tricoef(l,1)*u(trind(l,1,1),trind(l,1,2),trind(l,1,3)) &
       +tricoef(l,2)*u(trind(l,2,1),trind(l,2,2),trind(l,2,3)) &
       +tricoef(l,3)*u(trind(l,3,1),trind(l,3,2),trind(l,3,3)) &
       +tricoef(l,4)*u(trind(l,4,1),trind(l,4,2),trind(l,4,3))
    vip=tricoef(l,1)*v(trind(l,1,1),trind(l,1,2),trind(l,1,3)) &
       +tricoef(l,2)*v(trind(l,2,1),trind(l,2,2),trind(l,2,3)) &
       +tricoef(l,3)*v(trind(l,3,1),trind(l,3,2),trind(l,3,3)) &
       +tricoef(l,4)*v(trind(l,4,1),trind(l,4,2),trind(l,4,3))
    wip=tricoef(l,1)*w(trind(l,1,1),trind(l,1,2),trind(l,1,3)) &
       +tricoef(l,2)*w(trind(l,2,1),trind(l,2,2),trind(l,2,3)) &
       +tricoef(l,3)*w(trind(l,3,1),trind(l,3,2),trind(l,3,3)) &
       +tricoef(l,4)*w(trind(l,4,1),trind(l,4,2),trind(l,4,3))

    ! velocity at IB points
    uib=u(i0,j0,k0)
    vib=v(i0,j0,k0)
    wib=w(i0,j0,k0)

    ! location of IB points
    xib=xcd(i0,j0,k0)
    yib=ycd(i0,j0,k0)
    zib=zcd(i0,j0,k0)

    pro1=proiezioni(l,1)
    pro2=proiezioni(l,2)
    pro3=proiezioni(l,3)

    ! velocity at PP
    u1=uip
    u2=wip
    u3=vip

    ! velocity at PP in local reference frame
    vtangente1_pp=u1*rot(l,1,1)+u2*rot(l,1,2)+u3*rot(l,1,3)
    vtangente2_pp=u1*rot(l,2,1)+u2*rot(l,2,2)+u3*rot(l,2,3)
    vnormale_pp  =u1*rot(l,3,1)+u2*rot(l,3,2)+u3*rot(l,3,3)

    ! tangential speed at PP
    vtan=sqrt(vtangente1_pp**2.0+vtangente2_pp**2.0)

    ! determine case
    ! (0 = no vel; 1 = linear profile; 2 = full law)
    if (coef_wall<1) then
        caso=1 ! linear profile
    else if ((abs(xib-pro1)<1.e-9 .and. abs(yib-pro2)<1.e-9 .and. abs(zib-pro3)<1.e-9) .or. &
        (abs(uip)<1.e-9 .and. abs(vip)<1.e-9 .and. abs(wip)<1.e-9 )) then !1.d-9
        caso=0 ! zero velocity
    else
        caso=2
    end if

    !write(*,*) 'Case = ',caso
    !-----------------------------------------------------------------------
    ! LOG PROFILE
    if (caso==2) then

        ! find the u_tau with linear profile
        ypp=distanza
        u_tau_linear=sqrt(vtan/(ypp*re))

        ! find the u_tau with log profile
        !if (modelloparete==1) then
        !call wernerwengle(l,vtan)
        !else if (modelloparete==2) then

        ! first value for ustar
        if (ktime==1 .and. kiter==1) then
            !do ntime=ktime,1
            !do niter=kiter,1
            call wernerwengle(distanza,ustar(l),vtan)
            !ustar(l)=0.0
            !end do
        !end do
        end if

        ! if the wall is smooth, an internal subcycle begins
        if (rough==0) then
            contatore_star=0
            errore_star=1.1e-6 !errore_star=1.1d-6

            ! -------------------------------------------------------
            ! Newton - Raphson iterative procedure
            do while (errore_star>1.e-6 .and. contatore_star<50) !errore_star>1.d-6

                ustar(l)=abs(ustar(l))

                if(ustar(l)<1.e-8)then !ustar(l).lt.1.d-8
                    caso=0
                    errore_star=1.e-10 !errore_star = 1.d-10
                    exit
                end if

                ! since transition is not implemented, we use a simplified version
                ! logarithm (so that it's computed only once)

                f=ustar(l)*(kdynamic*log(abs(distanza*ustar(l)*re))+coef_rough_const)-vtan

                ! fprime is the derivative of f'
                fprime = kdynamic*(log(abs(distanza*ustar(l)*re)))+kdynamic+coef_rough_const

                ustar(l) = ustar(l) - f/fprime

                if(ustar(l)>1.E-8)then
                    errore_star = abs(f/sqrt(ustar(l)))
                end if


                ! older version -------
                ! argomentolog = 1.0 + ustar(l)*rougheight*re
                ! coef_rough = coef_rough - real(transition)*kdynamic*log(argomentolog)
                ! f=ustar(l)*(kdynamic*log(abs(distanza*ustar(l)*re))+coef_rough)-vtan
                ! ! fprime is the derivative of f'
                ! fprime = kdynamic*(log(abs(distanza*ustar(l)*re)))+kdynamic+coef_rough &
                !   -real(transition)*ustar(l)*kdynamic*rougheight/(1./re + rougheight*ustar(l))
                ! ustar(l) = ustar(l) - f/fprime
                ! if(ustar(l)>1.E-8)then
                !    errore_star = abs(f/sqrt(ustar(l)))
                ! end if
                ! coef_rough = coef_rough_const

                contatore_star = contatore_star +1

            end do   !end for do while#

        else if (rough==1) then
            right_hand_side=kdynamic*log(distanza/rougheight)
            ustar(l)=abs(vtan/right_hand_side)
        end if

        !end if

        ! choose the ustar
        !ustar(l) = max(ustar(l),u_tau_linear)
        ustar(l) = u_tau_linear

        if (ustar(l)>u_tau_linear) then
            write(*,*) 'U* Log'
        else
            write(*,*) 'U* Linear'
        end if

        ! I know ustar which is fix for the normal

        ! compute y+ of IB
        yib_plus = dist_ib_parete(l)*ustar(l)*re

        ! viscous layer
        if (yib_plus<=5.0) then
            vtan_ib=yib_plus*ustar(l)
            write(*,*) 'y+ Linear',yib_plus
        end if

        ! log layer
        if(yib_plus>=30.0)then
            vtan_ib = (kdynamic*log(yib_plus)+coef_rough_const)*ustar(l)
            write(*,*) 'y+ Log',yib_plus
            ! old:
            !vtan_ib = (kdynamic*log(yib_plus)+coef_rough)*ustar(l)
        end if

        ! buffer layer
        if(yib_plus > 5.0 .and. yib_plus < 30.0)then
            call smooth_ibm(yib_plus,vtan_ib,kdynamic,coef_rough_const)
            write(*,*) 'y+ Trans',yib_plus
            ! old:
            !call smooth_ibm(yib_plus,vtan_ib,kdynamic,coef_rough)
            vtan_ib = vtan_ib * ustar(l)
        end if

    end if !caso2


    !-----------------------------------------------------------------------
    ! update velocity at IB

    !write(*,*) 'Caso = ',caso
    ! zero velocity
    if (caso==0) then

        uib=0.0
        vib=0.0
        wib=0.0

    else if (caso==1) then

        ! linear profile
        uib=fattore_distanza*uip
        vib=fattore_distanza*vip
        wib=fattore_distanza*wip

        ustar(l)=sqrt( abs(uip) / (distanza*re) )

    else if (caso==2) then
        ! log profile

        ! normal and tangential component at IB in the local frame of reference
        vtangente1_ib =  vtangente1_pp * (vtan_ib/vtan)
        vtangente2_ib =  vtangente2_pp * (vtan_ib/vtan)
        vnormale_ib   =  vnormale_pp   * (vtan_ib/vtan) !  fattore_distanza**2.

        ! construct the cartesian component in the general frame of reference
        u1  = vtangente1_ib*rot_inverse(l,1,1)  &
            + vtangente2_ib*rot_inverse(l,1,2)  &
            + vnormale_ib  *rot_inverse(l,1,3)

        u2  = vtangente1_ib*rot_inverse(l,2,1)  &
            + vtangente2_ib*rot_inverse(l,2,2)  &
            + vnormale_ib  *rot_inverse(l,2,3)

        u3  = vtangente1_ib*rot_inverse(l,3,1)  &
            + vtangente2_ib*rot_inverse(l,3,2)  &
            + vnormale_ib  *rot_inverse(l,3,3)

        uib = u1
        vib = u3
        wib = u2

    end if

    ! compute components of shear stress
    local_shear1 =  vtangente1_pp * (ustar(l)*ustar(l)/vtan)
    local_shear2=  vtangente2_pp * (ustar(l)*ustar(l)/vtan)

    ! and pressure by interpolation
    fip= tricoef(l,1)*fi(trind(l,1,1),trind(l,1,2),trind(l,1,3)) &
        +tricoef(l,2)*fi(trind(l,2,1),trind(l,2,2),trind(l,2,3)) &
        +tricoef(l,3)*fi(trind(l,3,1),trind(l,3,2),trind(l,3,3)) &
        +tricoef(l,4)*fi(trind(l,4,1),trind(l,4,2),trind(l,4,3))

    fib = fi(i0,j0,k0)

    fi_pro = ((dist_ib_parete(l)+dist_pp_ib(l))*fib - fip*dist_ib_parete(l)) /dist_pp_ib(l)


    shear_ib(l,1) = local_shear1*rot_inverse(l,1,1)  &
        + local_shear2*rot_inverse(l,1,2) &
        + fi_pro*rot_inverse(l,1,3)

    shear_ib(l,2)= local_shear1*rot_inverse(l,2,1)  &
        + local_shear2*rot_inverse(l,2,2)  &
        + fi_pro*rot_inverse(l,2,3)

    shear_ib(l,3)= local_shear1*rot_inverse(l,3,1)  &
        + local_shear2*rot_inverse(l,3,2)  &
        + fi_pro*rot_inverse(l,3,3)



end subroutine compute_u_ib


!***********************************************************************
subroutine smooth_ibm(yib_plus,vtan_ib,kdynamic,coef_rough)
    !***********************************************************************

    implicit none

    !-----------------------------------------------------------------------

    integer,parameter :: n=2
    integer ibcbeg, ibcend
      
    real yib_plus,vtan_ib
    real kdynamic
    real coef_rough
      
    real ybcbeg,ybcend
    real t(2)
    real ubuf(2),ypp(2)
    real tval, yval, ypval, yppval
  
    !-----------------------------------------------------------------------
 
    t(1) = 5.
    t(2) = 30.

    ubuf(1) = 5.
    ubuf(2) = kdynamic*log(30.)+coef_rough

    ibcbeg = 1
    ybcbeg = 1. !this is the derivative on d u+/d y+
      
    ibcend = 1
    ybcend = kdynamic * (1./yib_plus) ! du+/dy+ for log law
      
    !-----------------------------------------------------------------------

    call spline_cubic_set(n,t,ubuf,ibcbeg,ybcbeg,ibcend,ybcend,ypp)
 
    tval = yib_plus

    call spline_cubic_val(n,t,ubuf,ypp,tval,yval,ypval,yppval)

    vtan_ib = yval

    return
end

!***********************************************************************
subroutine spline_cubic_set(n,t,y,ibcbeg,ybcbeg,ibcend,ybcend,ypp)
    !***********************************************************************
    !
    !c SPLINE_CUBIC_SET computes the second derivatives of a cubic spline.
    !
    !
    !  Discussion:
    !
    !    For data interpolation, the user must call SPLINE_CUBIC_SET to
    !    determine the second derivative data, passing in the data to be
    !    interpolated, and the desired boundary conditions.
    !
    !    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
    !    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
    !    evaluate the spline at any point.
    !
    !    The cubic spline is a piecewise cubic polynomial.  The intervals
    !    are determined by the "knots" or abscissas of the data to be
    !    interpolated.  The cubic spline has continous first and second
    !    derivatives over the entire interval of interpolation.
    !
    !    For any point T in the interval T(IVAL), T(IVAL+1), the form of
    !    the spline is
    !
    !      SPL(T) = A(IVAL)
    !             + B(IVAL) * ( T - T(IVAL) )
    !             + C(IVAL) * ( T - T(IVAL) )**2
    !             + D(IVAL) * ( T - T(IVAL) )**3
    !
    !    If we assume that we know the values Y(*) and YPP(*), which represent
    !    the values and second derivatives of the spline at each knot, then
    !    the coefficients can be computed as:
    !
    !      A(IVAL) = Y(IVAL)
    !      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
    !        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
    !      C(IVAL) = YPP(IVAL) / 2
    !      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
    !
    !    Since the first derivative of the spline is
    !
    !      SPL'(T) =     B(IVAL)
    !              + 2 * C(IVAL) * ( T - T(IVAL) )
    !              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
    !
    !    the requirement that the first derivative be continuous at interior
    !    knot I results in a total of N-2 equations, of the form:
    !
    !      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
    !      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
    !
    !    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
    !
    !      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
    !      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
    !      + YPP(IVAL-1) * H(IVAL-1)
    !      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
    !      =
    !      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
    !      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
    !
    !    or
    !
    !      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
    !      + YPP(IVAL) * H(IVAL)
    !      =
    !      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
    !      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
    !
    !    Boundary conditions must be applied at the first and last knots.
    !    The resulting tridiagonal system can be solved for the YPP values.
    !
    !  Modified:
    !
    !    20 November 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, the number of data points; N must be at least 2.
    !
    !    Input, real T(N), the points where data is specified.
    !    The values should be distinct, and increasing.
    !
    !    Input, real Y(N), the data values to be interpolated.
    !
    !    Input, integer IBCBEG, the left boundary condition flag:
    !
    !      0: the spline should be a quadratic over the first interval;
    !      1: the first derivative at the left endpoint should be YBCBEG;
    !      2: the second derivative at the left endpoint should be YBCBEG.
    !
    !    Input, real YBCBEG, the left boundary value, if needed.
    !
    !    Input, integer IBCEND, the right boundary condition flag:
    !
    !      0: the spline should be a quadratic over the last interval;
    !      1: the first derivative at the right endpoint should be YBCEND;
    !      2: the second derivative at the right endpoint should be YBCEND.
    !
    !    Input, real YBCEND, the right boundary value, if needed.
    !
    !    Output, real YPP(N), the second derivatives of the cubic spline.
    !
    implicit none
    !
    integer n
    !
    real diag(n)
    integer i
    integer ibcbeg
    integer ibcend
    real sub(2:n)
    real sup(1:n-1)
    real t(n)
    real y(n)
    real ybcbeg
    real ybcend
    real ypp(n)
    !
    !  Check.
    !
    if ( n <= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal errorc'
        write ( *, '(a)' ) '  The number of knots must be at least 2.'
        write ( *, '(a,i6)' ) '  The input value of N = ', n
        stop
    end if

    do i = 1, n-1
        if ( t(i) >= t(i+1) ) then
            write(*,'(a)')' '
            write(*,'(a)')'SPLINE_CUBIC_SET - Fatal errorc'
            write(*,'(a)')'  The knots must be strictly increasing, but'
            write(*,'(a,i6,a,g14.6)' ) '  T(',  i,') = ', t(i)
            write(*,'(a,i6,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
            stop
        end if
    end do
    !
    !  Set the first equation.
    !
    if ( ibcbeg == 0 ) then
        ypp(1) = 0.0E+00
        diag(1) = 1.0E+00
        sup(1) = -1.0E+00
    else if ( ibcbeg == 1 ) then
        ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
        diag(1) = ( t(2) - t(1) ) / 3.0E+00
        sup(1) = ( t(2) - t(1) ) / 6.0E+00
    else if ( ibcbeg == 2 ) then
        ypp(1) = ybcbeg
        diag(1) = 1.0E+00
        sup(1) = 0.0E+00
    else
        write(*,'(a)') ' '
        write(*,'(a)') 'SPLINE_CUBIC_SET - Fatal errorc'
        write(*,'(a)') '  The boundary flag IBCBEG must be 0, 1 or 2.'
        write(*,'(a,i6)' ) '  The input value is IBCBEG = ', ibcbeg
        stop
    end if
    !
    !  Set the intermediate equations.
    !
    do i = 2, n-1
        ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) )  &
            - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
        sub(i) = ( t(i) - t(i-1) ) / 6.0E+00
        diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00
        sup(i) = ( t(i+1) - t(i) ) / 6.0E+00
    end do
    !
    !  Set the last equation.
    !
    if ( ibcend == 0 ) then
        ypp(n) = 0.0E+00
        sub(n) = -1.0E+00
        diag(n) = 1.0E+00
    else if ( ibcend == 1 ) then
        ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
        sub(n) = ( t(n) - t(n-1) ) / 6.0E+00
        diag(n) = ( t(n) - t(n-1) ) / 3.0E+00
    else if ( ibcend == 2 ) then
        ypp(n) = ybcend
        sub(n) = 0.0E+00
        diag(n) = 1.0E+00
    else
        write (*,'(a)') ' '
        write (*,'(a)') 'SPLINE_CUBIC_SET - Fatal errorc'
        write (*,'(a)') '  The boundary flag IBCEND must be 0, 1 or 2.'
        write (*,'(a,i6)') '  The input value is IBCEND = ', ibcend
        stop
    end if
    !
    !  Special case:
    !    N = 2, IBCBEG = IBCEND = 0.
    !
    if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

        ypp(1) = 0.0E+00
        ypp(2) = 0.0E+00
    !
    !  Solve the linear system.
    !
    else

        call s3_fs ( sub, diag, sup, n, ypp, ypp )

    end if

    return
end

!***********************************************************************
subroutine spline_cubic_val(n,t,y,ypp,tval,yval,ypval,yppval)
    !***********************************************************************
    !
    !c SPLINE_CUBIC_VAL evaluates a cubic spline at a specific point.
    !
    !
    !  Discussion:
    !
    !    SPLINE_CUBIC_SET must have already been called to define the
    !    values of YPP.
    !
    !    For any point T in the interval T(IVAL), T(IVAL+1), the form of
    !    the spline is
    !
    !      SPL(T) = A
    !             + B * ( T - T(IVAL) )
    !             + C * ( T - T(IVAL) )**2
    !             + D * ( T - T(IVAL) )**3
    !
    !    Here:
    !      A = Y(IVAL)
    !      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
    !        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
    !      C = YPP(IVAL) / 2
    !      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
    !
    !  Modified:
    !
    !    20 November 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, the number of data values.
    !
    !    Input, real T(N), the knot values.
    !
    !    Input, real Y(N), the data values at the knots.
    !
    !    Input, real YPP(N), the second derivatives of the spline at the knots.
    !
    !    Input, real TVAL, a point, typically between T(1) and T(N), at
    !    which the spline is to be evalulated.  If TVAL lies outside
    !    this range, extrapolation is used.
    !
    !    Output, real YVAL, YPVAL, YPPVAL, the value of the spline, and
    !    its first two derivatives at TVAL.
    !
    implicit none
    !
    integer n
    !
    real dt
    real h
    integer left
    integer right
    real t(n)
    real tval
    real y(n)
    real ypp(n)
    real yppval
    real ypval
    real yval
    !
    !  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
    !  Values below T(1) or above T(N) use extrapolation.
    !
    call rvec_bracket ( n, t, tval, left, right )
    !
    !  Evaluate the polynomial.
    !
    dt = tval - t(left)
    h = t(right) - t(left)

    yval = y(left)   &
        + dt * ( ( y(right) - y(left) ) / h  &
        - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h  &
        + dt * ( 0.5E+00 * ypp(left)  &
        + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0E+00 * h ) ) ) )

    ypval = ( y(right) - y(left) ) / h  &
        - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h  &
        + dt * ( ypp(left)  &
        + dt * ( 0.5E+00 * ( ypp(right) - ypp(left) ) / h ) )

    yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

    return
end

!***********************************************************************
subroutine s3_fs ( a1, a2, a3, n, b, x )
    !***********************************************************************
    !
    !c S3_FS factors and solves a tridiagonal linear system.
    !
    !
    !  Note:
    !
    !    This algorithm requires that each diagonal entry be nonzero.
    !
    !  Modified:
    !
    !    05 December 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, real A1(2:N), A2(1:N), A3(1:N-1).
    !    On input, the nonzero diagonals of the linear system.
    !    On output, the data in these vectors has been overwritten
    !    by factorization information.
    !
    !    Input, integer N, the order of the linear system.
    !
    !    Input/output, real B(N).
    !    On input, B contains the right hand side of the linear system.
    !    On output, B has been overwritten by factorization information.
    !
    !    Output, real X(N), the solution of the linear system.
    !
    implicit none
    !
    integer n
    !
    real a1(2:n)
    real a2(1:n)
    real a3(1:n-1)
    real b(n)
    integer i
    real x(n)
    real xmult
    !
    !  The diagonal entries can't be zero.
    !
    do i = 1, n
        if ( a2(i) == 0.0E+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'S3_FS - Fatal errorc'
            write ( *, '(a,i6,a)' ) '  A2(', i, ') = 0.'
            return
        end if
    end do

    do i = 2, n-1

        xmult = a1(i) / a2(i-1)
        a2(i) = a2(i) - xmult * a3(i-1)

        b(i) = b(i) - xmult * b(i-1)

    end do

    xmult = a1(n) / a2(n-1)
    a2(n) = a2(n) - xmult * a3(n-1)

    x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
    do i = n-1, 1, -1
        x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
    end do

    return
end
      
!***********************************************************************
subroutine rvec_bracket ( n, x, xval, left, right )
    !***********************************************************************
    !
    !c RVEC_BRACKET searches a sorted array for successive brackets of a value.
    !
    !
    !  Discussion:
    !
    !    If the values in the vector are thought of as defining intervals
    !    on the real line, then this routine searches for the interval
    !    nearest to or containing the given value.
    !
    !  Modified:
    !
    !    06 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, length of input array.
    !
    !    Input, real X(N), an array sorted into ascending order.
    !
    !    Input, real XVAL, a value to be bracketed.
    !
    !    Output, integer LEFT, RIGHT, the results of the search.
    !    Either:
    !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
    !      XVAL > X(N), when LEFT = N-1, RIGHT = N;
    !    or
    !      X(LEFT) <= XVAL <= X(RIGHT).
    !
    implicit none
    !
    integer n
    !
    integer i
    integer left
    integer right
    real x(n)
    real xval
    !
    do i = 2, n - 1

        if ( xval < x(i) ) then
            left = i - 1
            right = i
            return
        end if

    end do

    left = n - 1
    right = n

    return
end




















