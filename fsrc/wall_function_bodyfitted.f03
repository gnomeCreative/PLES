!***********************************************************************
subroutine wall_function_bodyfitted(ktime,niter,tipo,i_rest)
    !***********************************************************************
    use myarrays_velo3
    use myarrays_metri3
    use mysending
    use myarrays_wallmodel
    !
    use scala3
    use tipologia
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     variables declaration
    real x1,y1,z1
    real x2,y2,z2
    real x3,y3,z3
    real a,b,c,d
      
    real distanza,alfa

    real u_t_sotto,u_t_sopra
    real tau_sotto,tau_sopra
      
    integer i,j,k,kk
    integer ktime,niter,ierr,iwall
    integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
    integer i_rest
    !-----------------------------------------------------------------------
    att_mod_par = 1
    !-----------------------------------------------------------------------
    !             FACE 3
    !-----------------------------------------------------------------------
    if (ktime==1) then

        !     first value for ustar

        ! chicco TROVARE UN CRITERIO PIU' GENERALE PER U_T????
        u_t=1


        ! chicco questa va bene per cartesiano!!! per curvilineo dovrebbe trovare
        !      la coordinata sulla faccia con la normale al piano !!!!!!
     
        do iwall = 1,wfp3
            j=1
            do k=kparasta,kparaend
                do i=1,jx
	 
                    !           POINT 1  cell centroid
                    punto_wfp3(1,1,i,k) = .125*( &
                        x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
                        +x(i  ,j-1,k-1)+x(i  ,j-1,k  ) &
                        +x(i-1,j  ,k  )+x(i-1,j  ,k-1) &
                        +x(i-1,j-1,k-1)+x(i-1,j-1,k  ))

                    punto_wfp3(1,2,i,k) = .125*( &
                        y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
                        +y(i  ,j-1,k-1)+y(i  ,j-1,k  ) &
                        +y(i-1,j  ,k  )+y(i-1,j  ,k-1) &
                        +y(i-1,j-1,k-1)+y(i-1,j-1,k  ))

                    punto_wfp3(1,3,i,k) = .125*( &
                        z(i  ,j  ,k  )+z(i  ,j  ,k-1) &
                        +z(i  ,j-1,k-1)+z(i  ,j-1,k  ) &
                        +z(i-1,j  ,k  )+z(i-1,j  ,k-1) &
                        +z(i-1,j-1,k-1)+z(i-1,j-1,k  ))


                    !          triangle at the wall
                    x1 = x(i,0,k)
                    y1 = y(i,0,k)
                    z1 = z(i,0,k)
	   
                    x2 = x(i-1,0,k)
                    y2 = y(i-1,0,k)
                    z2 = z(i-1,0,k)

                    x3 = x(i,0,k-1)
                    y3 = y(i,0,k-1)
                    z3 = z(i,0,k-1)

                    call piano_per3punti ( x1, y1, z1,  &
                        x2, y2, z2,  &
                        x3, y3, z3, &
                        a, b, c, d )

                    call punto_proiezione_piano(a,b,c,d,  &
                        punto_wfp3(1,1,i,k),punto_wfp3(1,2,i,k),punto_wfp3(1,3,i,k), &
                        punto_wfp3(2,1,i,k),punto_wfp3(2,2,i,k),punto_wfp3(2,3,i,k))

                    !          check if the point is in the triangle


                    !           POINT 2
                    !            punto_wfp3(2,1,i,k) =
                    !     >                    .25*(x(i  ,0  ,k  )+x(i  ,0  ,k-1)
                    !     >                        +x(i-1,0  ,k  )+x(i-1,0  ,k-1))
            
                    !	    punto_wfp3(2,2,i,k) =
                    !     >	                  .25*(y(i  ,0  ,k  )+y(i  ,0  ,k-1)
                    !     >                        +y(i-1,0  ,k  )+y(i-1,0  ,k-1))
            
                    !	    punto_wfp3(2,3,i,k) =
                    !     >                    .25*(z(i  ,0  ,k  )+z(i  ,0  ,k-1)
                    !     >                        +z(i-1,0  ,k  )+z(i-1,0  ,k-1))

                    !           DISTANCE BETWEEN POINT 1 AND 2
                    punto_wfp3(3,1,i,k) = sqrt( &
                        (punto_wfp3(1,1,i,k)-punto_wfp3(2,1,i,k))**2.     &
                        +(punto_wfp3(1,2,i,k)-punto_wfp3(2,2,i,k))**2.  &
                        +(punto_wfp3(1,3,i,k)-punto_wfp3(2,3,i,k))**2.)
	    
                end do
            end do
        end do
      
    end if


    if (ktime==1 .and. i_rest==0)then
        u_t = 1.
        utangente = 1.
        att_mod_par = 0
    else
        if (wfp3==1) then

            j = 1
            do k=kparasta,kparaend
                do i=1,jx
                    !giulia solo su celle fluide
                    if(tipo(i,j,k).eq.2)then
                        x1 = punto_wfp3(1,1,i,k)
                        y1 = punto_wfp3(1,2,i,k)
                        z1 = punto_wfp3(1,3,i,k)

                        x2 = punto_wfp3(2,1,i,k)
                        y2 = punto_wfp3(2,2,i,k)
                        z2 = punto_wfp3(2,3,i,k)

                        !        this is a position vector
                        x3 = x1 + u(i,j,k)
                        y3 = y1 + v(i,j,k)
                        z3 = z1 + w(i,j,k)
     
                        call angolo(x1,y1,z1,x2,y2,z2,x1,y1,z1,x3,y3,z3,alfa)
           
                        distanza = punto_wfp3(3,1,i,k)
             
                        utangente(i,1,k)= sqrt( u(i,j,k)*u(i,j,k) &
                            +v(i,j,k)*v(i,j,k)      &
                            +w(i,j,k)*w(i,j,k) ) &
                            *sin(alfa)

                        if(utangente(i,1,k)>0.000001)then
                            call compute_ustar(distanza,re,utangente(i,1,k),u_t(i,1,k),rough,z0)

                            !    switch off the wall function if y+<11 or if it is not a fluid node
                            !    in case of ibm
                            if(u_t(i,1,k)*distanza*re.le.11. .or.tipo(i,j,k).ne.2)then
                                att_mod_par(i,1,k)=0
                            end if
                        else
                            att_mod_par(i,1,k)=0
                        end if
                    else !punto solido
                        ! ALE: I add (i,1,k) to the variables
                        u_t(i,1,k) = 1.
                        utangente(i,1,k) = 1.
                        att_mod_par(i,1,k) = 0
                    end if
                end do
            end do
        end if  ! wfp3


        if (wfp3==1) then
            do k=kparasta,kparaend
                do i=1,jx
                    att_mod_par(i,1,k)=0
                end do
            end do
        end if
      
    end if
    !-----------------------------------------------------------------------
    !             FACE 4
    !-----------------------------------------------------------------------

    if (ktime==1) then

     
        if (wfp4==1) then
            j=jy
            do k=kparasta,kparaend
                do i=1,jx
	 
                    !           POINT 1
                    punto_wfp4(1,1,i,k) = .125*( &
                        x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
                        +x(i  ,j-1,k-1)+x(i  ,j-1,k  ) &
                        +x(i-1,j  ,k  )+x(i-1,j  ,k-1) &
                        +x(i-1,j-1,k-1)+x(i-1,j-1,k  ))

                    punto_wfp4(1,2,i,k) = .125*( &
                        y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
                        +y(i  ,j-1,k-1)+y(i  ,j-1,k  ) &
                        +y(i-1,j  ,k  )+y(i-1,j  ,k-1) &
                        +y(i-1,j-1,k-1)+y(i-1,j-1,k  ))

                    punto_wfp4(1,3,i,k) = .125*( &
                        z(i  ,j  ,k  )+z(i  ,j  ,k-1) &
                        +z(i  ,j-1,k-1)+z(i  ,j-1,k  ) &
                        +z(i-1,j  ,k  )+z(i-1,j  ,k-1) &
                        +z(i-1,j-1,k-1)+z(i-1,j-1,k  ))
                    !          triangle at the wall
                    x1 = x(i,jy,k)
                    y1 = y(i,jy,k)
                    z1 = z(i,jy,k)
	   
                    x2 = x(i-1,jy,k)
                    y2 = y(i-1,jy,k)
                    z2 = z(i-1,jy,k)

                    x3 = x(i,jy,k-1)
                    y3 = y(i,jy,k-1)
                    z3 = z(i,jy,k-1)

                    call piano_per3punti ( x1, y1, z1,  &
                        x2, y2, z2,  &
                        x3, y3, z3, &
                        a, b, c, d )

                    call punto_proiezione_piano(a,b,c,d,  &
                        punto_wfp4(1,1,i,k),punto_wfp4(1,2,i,k),punto_wfp4(1,3,i,k), &
                        punto_wfp4(2,1,i,k),punto_wfp4(2,2,i,k),punto_wfp4(2,3,i,k))

                    !           POINT 2
                    !            punto_wfp4(2,1,i,k) =
                    !     >                       .25*(x(i  ,jy,k  )+x(i  ,jy,k-1)
                    !     >                           +x(i-1,jy,k  )+x(i-1,jy,k-1))
            
                    !	    punto_wfp4(2,2,i,k) =
                    !     >	                     .25*(y(i  ,jy,k  )+y(i  ,jy,k-1)
                    !     >                           +y(i-1,jy,k  )+y(i-1,jy,k-1))
            
                    !	    punto_wfp4(2,3,i,k) =
                    !     >                       .25*(z(i  ,jy,k  )+z(i  ,jy,k-1)
                    !     >                           +z(i-1,jy,k  )+z(i-1,jy,k-1))

                    !           DISTANCE BETWEEN POINTS 1 AND 2
                    punto_wfp4(3,1,i,k) = sqrt( &
                        (punto_wfp4(1,1,i,k)-punto_wfp4(2,1,i,k))**2.     &
                        +(punto_wfp4(1,2,i,k)-punto_wfp4(2,2,i,k))**2.  &
                        +(punto_wfp4(1,3,i,k)-punto_wfp4(2,3,i,k))**2.)
	    
                end do
            end do
        end if !wpf4
    end if


    if(ktime .eq. 1 .and. i_rest==0)then
        u_t = 1.
        utangente = 1.
        att_mod_par = 0
    else
        if (wfp4==1) then
            j = jy
            do k=kparasta,kparaend
                do i=1,jx

                    if(tipo(i,j,k).eq.2)then
	     	      	           
                        x1 = punto_wfp4(1,1,i,k)
                        y1 = punto_wfp4(1,2,i,k)
                        z1 = punto_wfp4(1,3,i,k)

                        x2 = punto_wfp4(2,1,i,k)
                        y2 = punto_wfp4(2,2,i,k)
                        z2 = punto_wfp4(2,3,i,k)

                        !        this is a position vector
                        x3 = x1 + u(i,j,k)
                        y3 = y1 + v(i,j,k)
                        z3 = z1 + w(i,j,k)
      
                        call angolo(x1,y1,z1,x2,y2,z2,x1,y1,z1,x3,y3,z3,alfa)
            
                        distanza = punto_wfp4(3,1,i,k)
            
                        utangente(i,2,k)= sqrt( u(i,j,k)*u(i,j,k) &
                            +v(i,j,k)*v(i,j,k)      &
                            +w(i,j,k)*w(i,j,k) ) &
                            *sin(alfa)
                        !         if(i==10 .and. k==10)write(*,*)'AAA',utangente(i,2,k)

                        if(abs(utangente(i,2,k))>0.000001)then
                            call compute_ustar(distanza,re,utangente(i,2,k),u_t(i,2,k), &
                                rough,z0)

                            !    switch off the wall function if y+<11 or if it is not a fluid node
                            !    in case of ibm
                            if(u_t(i,2,k)*distanza*re .le. 11. .or.tipo(i,j,k).ne.2)then
                                att_mod_par(i,2,k)=0
                            end if
                        else
                            att_mod_par(i,2,k)=0
                        end if
	 
	 
                    !         if(i==10 .and. k==10)write(*,*)'AAA',u_t(i,2,k),
                    !     >	                               att_mod_par(i,2,k)
	 
	 
                    else !punto solido
                        u_t = 1.
                        utangente = 1.
                        att_mod_par = 0
                    end if
                end do
            end do
        end if  ! wfp4
      
      
        if (wfp4==0) then
            do k=kparasta,kparaend
                do i=1,jx
                    att_mod_par(i,2,k)=0
                end do
            end do
        end if
           
    end if
    !-----------------------------------------------------------------------
    !     COMUNICATION BETWEEN PROCS
    !-----------------------------------------------------------------------

    u_t_sotto = 0.
    u_t_sopra = 0.
    do k=kparasta,kparaend
        do i=1,jx
            do kk=1,att_mod_par(i,1,k)
                u_t_sotto = u_t_sotto +u_t(i,1,k)*u_t(i,1,k)
            enddo
            do kk=1,att_mod_par(i,2,k)
                u_t_sopra = u_t_sopra +u_t(i,2,k)*u_t(i,2,k)
            end do

        end do
    end do


    call MPI_ALLREDUCE(u_t_sotto,tau_sotto,1,MPI_REAL_SD, &
        MPI_SUM,MPI_COMM_WORLD,ierr)
    
    call MPI_ALLREDUCE(u_t_sopra,tau_sopra,1,MPI_REAL_SD, &
        MPI_SUM,MPI_COMM_WORLD,ierr)

    if(myid.eq.0)then
        tau_sotto = tau_sotto/real(jx)/real(jz)
        tau_sopra = tau_sopra/real(jx)/real(jz)
        write(*,*)'tau',tau_sotto,tau_sopra
    end if
              
    return
end

!***********************************************************************
subroutine compute_ustar(distanza,reynol,umedia,u_t,rough,rougheight)

    implicit none

    !-----------------------------------------------------------------------
    !     variables declaration

    real,intent(in) :: rougheight,reynol,umedia,distanza
    integer,intent(in) :: rough
    real,intent(out) :: u_t
    !
    real errore_star,argomentolog,f,fprime,tempv
    integer contatore_star


    !-----------------------------------------------------------------------
    ! parameters law of the wall
    real,parameter :: vonKarman=0.41 ! k
    real,parameter :: coef_rough=5.1
    real,parameter :: kdynamic=1/vonKarman   ! 1/k with k von karman constant
    real,parameter :: transition=0.0

    ! parameters iterative procedure
    integer,parameter :: max_iterations=100
    real,parameter :: convergence_error=1.e-4
    real,parameter :: min_ut=1.0e-8

    ! coef_rough=5.1

    !-----------------------------------------------------------------------
    !     CASE 1: smooth surface

    if (rough==0) then
        !     Newton - Rapshon iterative procedure
        contatore_star = 0
        errore_star = 1000.
        do while (errore_star>convergence_error .and. contatore_star<max_iterations)
            tempv=u_t
            u_t = abs(u_t)

            ! argomentolog = 1.0+ u_t*rougheight*reynol ! ALE

            ! coef_rough = coef_rough  & ! ALE
                !- real(transition)*kdynamic*log(argomentolog)! ALE

            f=u_t*(kdynamic*log(abs(distanza*u_t*reynol))+ coef_rough) - umedia

            ! fprime is f derivative
            fprime = kdynamic*log(abs(distanza*u_t*reynol))+kdynamic+coef_rough !& ! ALE
            ! -real(transition)*u_t*kdynamic*rougheight/(1./reynol + rougheight*u_t) ! ALE

            u_t = u_t - f/fprime

            if (u_t>min_ut) then
                errore_star=abs(u_t-tempv)
            !   errore_star = abs(f/sqrt(u_t))
            end if

            contatore_star = contatore_star +1
            ! coef_rough = 5.1

        end do   !end loop do while
    !end do !smooth surface

    !-----------------------------------------------------------------------
    !     CASE 2: rough surface

    elseif (rough==1) then
        u_t=vonKarman*umedia/log(distanza/rougheight)
    end if


    return
end


!***********************************************************************
subroutine angolo (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,angle_deg)
    !***********************************************************************
    ! THIS SUBROUTINE IS TAKEN FROM LIBRARY GEOMETRY PACK
    !
    !c LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
    !
    !  Formula:
    !
    !    The explicit form of a line in 3D is:
    !
    !      (X1,Y1,Z1), (X2,Y2,Z2).
    !
    !  Modified:
    !
    !    24 January 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real X1, Y1, Z1, X2, Y2, Z2, two distince points on the first line.
    !
    !    Input, real X3, Y3, Z3, X4, Y4, Z4, two distinct points on the second line.
    !
    !    Output, real ANGLE, the angle in radians between the two lines.
    !    ANGLE is computed using the ACOS function, and so lies between 0 and PI.
    !    But if one of the lines is degenerate, ANGLE is returned as -1.0.
    !
    implicit none
    !
    real,intent(in) :: x1,x2,x3,x4
    real,intent(in) :: y1,y2,y3,y4
    real,intent(in) :: z1,z2,z3,z4
    real,intent(out) :: angle_deg


    real :: angle,arc_cosine
    real :: ctheta,enorm0_3d,pdotq,pnorm,qnorm
    real,parameter :: tol=0.000001
    real,parameter :: pi=acos(-1.0)

    pnorm = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )
    qnorm = sqrt ( ( x3 - x4 )**2 + ( y3 - y4 )**2 + ( z3 - z4 )**2 )

    pdotq =    ( x2 - x1 ) * ( x4 - x3 ) &
        + ( y2 - y1 ) * ( y4 - y3 ) &
        + ( z2 - z1 ) * ( z4 - z3 )

    if (pnorm==0.0 .or. qnorm==0.0) then
        write (*,*) ' '
        write (*,*) 'LINES_EXP_ANGLE_3D - Fatal error!'
        write (*,*) '  One of the lines is degenerate!'
        angle = - 1.0
    else

        !      write(*,*)'pi: ',pi
        ctheta = pdotq / ( pnorm * qnorm )

        if(abs(ctheta).ge.1.)then
            write(1000,*)'TROVATO ANGOLO DEGENERE',ctheta
            angle_deg = 0.0
        else
            angle =abs(acos(ctheta)) ! arc_cosine ( ctheta )
            !         write(*,*)'angolo',angle
            angle_deg=angle*180./pi
            angle_deg=angle
        !         write(*,*)'angolo deg: ',angle_deg
        end if
    end if


    return
end
