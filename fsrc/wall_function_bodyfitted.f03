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
      
    real distanza,alfa,cor

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
    if(ktime .eq. 1)then

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


    if(ktime .eq. 1 .and. i_rest==0)then
        u_t = 1.
        utangente = 1.
        att_mod_par = 0
    else
        do iwall = 1,wfp3
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
     
                        call angolo(x1,y1,z1,x2,y2,z2,x1,y1,z1,x3,y3,z3,alfa,cor)
           
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
                        u_t = 1.
                        utangente = 1.
                        att_mod_par = 0
                    end if
                end do
            end do
        end do  ! iwall


        do iwall = 1,1-wfp3
            do k=kparasta,kparaend
                do i=1,jx
                    att_mod_par(i,1,k)=0
                end do
            end do
        end do
      
    end if
    !-----------------------------------------------------------------------
    !             FACE 4
    !-----------------------------------------------------------------------

    if(ktime .eq. 1)then

     
        do iwall = 1,wfp4
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
        end do
    end if


    if(ktime .eq. 1 .and. i_rest==0)then
        u_t = 1.
        utangente = 1.
        att_mod_par = 0
    else
        do iwall = 1,wfp4
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
      
                        call angolo(x1,y1,z1,x2,y2,z2,x1,y1,z1,x3,y3,z3,alfa,cor)
            
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
        end do  ! iwall
      
      
        do iwall = 1,1-wfp4
            do k=kparasta,kparaend
                do i=1,jx
                    att_mod_par(i,2,k)=0
                end do
            end do
        end do
           
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

!-----------------------------------------------------------------------
! SUBROUTINE
!-----------------------------------------------------------------------

!***********************************************************************
subroutine piano_per3punti ( x1, y1, z1,  &
    x2, y2, z2,  &
    x3, y3, z3,  &
    a, b, c, d )
    !-----------------------------------------------------------------------
    ! dati tre punti di coordinate (x1,y1,z1),(x2,y2,z2),(x3,y3,z3) trovo
    ! il piano che vi passa di eq: ax+by+cz+d=0
    !-----------------------------------------------------------------------
    implicit none
    !
    real a,b,c,d
    real x1,y1,z1
    real x2,y2,z2
    real x3,y3,z3

    a = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
    b = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
    c = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )
    d = - x2 * a - y2 * b - z2 * c

    return
end


!***********************************************************************
subroutine punto_proiezione_piano(a,b,c,d,x,y,z,xn,yn,zn)
    !***********************************************************************
    !------------------------------------------------------------
    ! dato un piano di equazione ax+by+cz+d=0 ed un punto (x,y,z)
    ! trovo il punto appartenente al piano (xn,yn,zn) piu' vicino
    ! al punto sopra
    !------------------------------------------------------------
    ! COME PROCEDO:
    ! normale N al piano (a,b,c)
    ! la linea definita da (xn-x)/a = (yn-y)/b = (zn-z)/c = t
    ! passa per il punto (x,y,z) ed e' parallela a N
    !
    ! risolvendo per il punto (xn,yn,zn) ottengo:
    !
    !  xn = a * t + x
    !  yn = b * t + y
    !  zn = c * t + z
    !
    ! pongo questi valori nell'equazione del piano e risolvo per t
    !
    ! a*(a * t + x) + b*(b * t + y) + c*(c * t + z) + d = 0
    !
    ! t=(-a*x -b*y -c*z -d) / (a*a + b*b + c*c)
    !------------------------------------------------------------

    implicit none

    real a,b,c,d
    real x,y,z
    real xn,yn,zn
    real t

    if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_IMP_POINT_NEAR_3D - Fatal error!'
        write ( *, '(a)' ) '  A = B = C = 0.'
        stop
    else

        t = -( a * x + b * y + c * z + d )/( a * a + b * b + c * c )

        xn = x + a * t
        yn = y + b * t
        zn = z + c * t

    end if
    return
end


!      subroutine trovopunti













