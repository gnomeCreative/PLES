!***********************************************************************
subroutine readtau(ti,r)
    !***********************************************************************
    ! read flux tau_wall at surface and top layer, M-O similarity
    use mysending
    use myarrays_metri3
    use myarrays_velo3
    use myarrays2           !for areola
    use myarrays_wallmodel  !for att_mod_par
    use myarrays_WB         !for Fx
      
    use myarrays_moisture
    use myarrays_nesting
    use scala3

    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer i,j,k,ii,jj,kk
    real ti,taulin
    real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
    !-----------------------------------------------------------------------

    !     tau wall from wrf interpolated in time
    taulin = &
        ((ti-ti_pom_new)/(ti_pom_old-ti_pom_new))*tauw3nest(ntke-1) &
        -((ti-ti_pom_old)/(ti_pom_old-ti_pom_new))*tauw3nest(ntke)

    if(myid==0)then
        if(att_mod_par(1,1,1)==0)then
            write(*,*) 'MO bot spento '
        elseif(att_mod_par(1,1,1)==1)then
            write(*,*) 'MO bot attivo '
        end if
    end if
    !hicco NOTA se uso u che e' turbolenta forse non serve Fx
    do k=kparasta,kparaend
        do i=1,jx
            do ii=1,1-att_mod_par(i,1,k)
                f2(i,0,k)=annitV(i,0,k)*g22(i,0,k) &
                    *(-8.*r(i,0,k)+9.*r(i,1,k)-r(i,2,k))/3.
            end do
            do ii=1,att_mod_par(i,1,k)
                f2(i,0,k)= taulin		&		  !wall stress &
                    *(r(i,1,k)/abs(r(i,1,k)))    &              !sign &
                    *(r(i,1,k)**2/(u(i,1,k)**2.+w(i,1,k)**2.)) & !fraction on the total stress &
                    *areola3(i,k)	&			  !for the integral &
                    *(0.1*Fx(i,k) + 1.)			  !random disturbance
            end do
        end do
    end do
      
    if(myid==0) write(*,*) 'tau_w ', taulin, ' f2 ', f2(1,0,1)
    !-----------------------------------------------------------------------
    return
end
      
      
      
      
      
