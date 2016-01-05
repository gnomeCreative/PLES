!***********************************************************************
subroutine readflux(ti,isc)
    !***********************************************************************
    ! read flux BC at surface and top layer, M-O similarity
    use mysending
    use myarrays_metri3
    use myarrays_velo3
    use myarrays2
    use myarrays_wallmodel
    use myarrays_WB
      
    use myarrays_moisture
    use myarrays_nesting
    use scala3

    implicit none
    !
    !------------------------------------------------------------------------
    !     array declaration
    integer isc
    integer i,j,k

    real ti

    ! jk: per calcolare theta_e
    real theta_e !temperatura potenziale equivalente
    real rt      !rapporto di mescolanza dell'acqua
    real e,es		! pressione di vapore, pressione di vapore saturo
    real ustar,thetastar
      
    real scaflux3,scaflux4
    !------------------------------------------------------------------------
    !     side bottom
    !
    !     scalar flux from wrf interpolated in time
    scaflux3 =       &
        ((ti_pom_new - ti)*scalarflux3nest(isc,ntke-1) &
        +(ti - ti_pom_old)*scalarflux3nest(isc,ntke  ))  &
        /(ti_pom_new-ti_pom_old)
           	 
    do k=kparasta,kparaend
        do i=1,jx
      
            f2(i,0,k)= scaflux3*(0.1*Fx(i,k) + 1.)*areola3(i,k)
      
        end do
    end do
      
    if(myid==0)write(*,*)'scalare:',isc, 'flux: ',scaflux3,f2(1,0,1)
    !
    !------------------------------------------------------------------------
    !     side top
    !
    !     scalar flux from wrf interpolated in time
    !      scaflux4 =
    !     >  ((ti_pom_new - ti)*scalarflux4nest(isc,ntke-1)
    !     >  +(ti - ti_pom_old)*scalarflux4nest(isc,ntke  ))
    !     >  /(ti_pom_new-ti_pom_old)
      
    !      do k=kparasta,kparaend
    !      do i=1,jx

    !         f2(i,jy,k) = scaflux4*areola4(i,k)*(0.1*Fx(i,k)+ 1.)

    !      end do
    !      end do
      
    !      if(myid==0)write(*,*)'scalare:',isc, 'flux: ',scaflux3,f2(1,0,1)


    return
end
