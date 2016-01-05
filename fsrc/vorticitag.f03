!***********************************************************************
subroutine vorticitag(myid,nproc,kparasta,kparaend)
    !***********************************************************************
    ! compute vorticity (3D vector) on curvilinear grid
    ! for each point
  
    use myarrays_metri3
    use myarrays_velo3
    use myarrays_LC
    !
    use scala3
    !
    use mpi

    implicit none
    !-----------------------------------------------------------------------
    !     array declaration
    integer i,j,k,ii,jj,kk
    integer myid,nproc,kparasta,kparaend

    real apxcsi,apycsi,apzcsi
    real apxeta,apyeta,apzeta
    real apxzet,apyzet,apzzet
    real apgiac
    real apcsx_g1,apcsy_g1,apcsz_g1
    real apetx_g1,apety_g1,apetz_g1
    real apztx_g1,apzty_g1,apztz_g1
    real dudcsi,dvdcsi,dwdcsi
    real dudeta,dvdeta,dwdeta
    real dudzita,dvdzita,dwdzita

    real dudx,dudy,dudz
    real dvdx,dvdy,dvdz
    real dwdx,dwdy,dwdz
    !
    !-----------------------------------------------------------------------
    !
    do k=kparasta,kparaend
        do j=1,n2
            do i=1,n1

                apxcsi=.25*(x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
                    +x(i  ,j-1,k-1)+x(i  ,j-1,k  ))-   &
                    .25*(x(i-1,j  ,k  )+x(i-1,j  ,k-1) &
                    +x(i-1,j-1,k-1)+x(i-1,j-1,k  ))
       
                apycsi=.25*(y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
                    +y(i  ,j-1,k-1)+y(i  ,j-1,k  ))-   &
                    .25*(y(i-1,j  ,k  )+y(i-1,j  ,k-1) &
                    +y(i-1,j-1,k-1)+y(i-1,j-1,k  ))
         
                apzcsi=.25*(z(i  ,j  ,k  )+z(i  ,j  ,k-1) &
                    +z(i  ,j-1,k-1)+z(i  ,j-1,k  ))-   &
                    .25*(z(i-1,j  ,k  )+z(i-1,j  ,k-1) &
                    +z(i-1,j-1,k-1)+z(i-1,j-1,k  ))
        
                apxeta=.25*(x(i  ,j  ,k  )+x(i  ,j  ,k-1) &
                    +x(i-1,j  ,k-1)+x(i-1,j  ,k  ))-   &
                    .25*(x(i  ,j-1,k  )+x(i  ,j-1,k-1) &
                    +x(i-1,j-1,k-1)+x(i-1,j-1,k  ))
        
                apyeta=.25*(y(i  ,j  ,k  )+y(i  ,j  ,k-1) &
                    +y(i-1,j  ,k-1)+y(i-1,j  ,k  ))-   &
                    .25*(y(i  ,j-1,k  )+y(i  ,j-1,k-1) &
                    +y(i-1,j-1,k-1)+y(i-1,j-1,k  ))
        
                apzeta=.25*(z(i  ,j  ,k  )+z(i  ,j  ,k-1) &
                    +z(i-1,j  ,k-1)+z(i-1,j  ,k  ))-   &
                    .25*(z(i  ,j-1,k  )+z(i  ,j-1,k-1) &
                    +z(i-1,j-1,k-1)+z(i-1,j-1,k  ))
        
                apxzet=.25*(x(i  ,j  ,k  )+x(i  ,j-1,k  ) &
                    +x(i-1,j-1,k  )+x(i-1,j  ,k  ))-   &
                    .25*(x(i  ,j  ,k-1)+x(i  ,j-1,k-1) &
                    +x(i-1,j-1,k-1)+x(i-1,j  ,k-1))
         
                apyzet=.25*(y(i  ,j  ,k  )+y(i  ,j-1,k  ) &
                    +y(i-1,j-1,k  )+y(i-1,j  ,k  ))-   &
                    .25*(y(i  ,j  ,k-1)+y(i  ,j-1,k-1) &
                    +y(i-1,j-1,k-1)+y(i-1,j  ,k-1))
         
                apzzet=.25*(z(i  ,j  ,k  )+z(i  ,j-1,k  ) &
                    +z(i-1,j-1,k  )+z(i-1,j  ,k  ))-   &
                    .25*(z(i  ,j  ,k-1)+z(i  ,j-1,k-1) &
                    +z(i-1,j-1,k-1)+z(i-1,j  ,k-1))

                apgiac=apxcsi*(apyeta*apzzet-apyzet*apzeta)-        &
                    apxeta*(apycsi*apzzet-apyzet*apzcsi)+   &
                    apxzet*(apycsi*apzeta-apyeta*apzcsi)

                apcsx_g1 = (apyeta*apzzet - apyzet*apzeta)/apgiac
                apcsy_g1 = (apxzet*apzeta - apxeta*apzzet)/apgiac
                apcsz_g1 = (apxeta*apyzet - apxzet*apyeta)/apgiac

                apetx_g1 = (apyzet*apzcsi - apycsi*apzzet)/apgiac
                apety_g1 = (apxcsi*apzzet - apxzet*apzcsi)/apgiac
                apetz_g1 = (apxzet*apycsi - apxcsi*apyzet)/apgiac

                apztx_g1 = (apycsi*apzeta - apyeta*apzcsi)/apgiac
                apzty_g1 = (apxeta*apzcsi - apxcsi*apzeta)/apgiac
                apztz_g1 = (apxcsi*apyeta - apxeta*apycsi)/apgiac

                dudcsi = .5*( u(i+1,j  ,k  ) - u(i-1,j  ,k  ) )
                dudeta = .5*( u(i  ,j+1,k  ) - u(i  ,j-1,k  ) )
                dudzita= .5*( u(i  ,j  ,k+1) - u(i  ,j  ,k-1) )

                dvdcsi = .5*( v(i+1,j  ,k  ) - v(i-1,j  ,k  ) )
                dvdeta = .5*( v(i  ,j+1,k  ) - v(i  ,j-1,k  ) )
                dvdzita= .5*( v(i  ,j  ,k+1) - v(i  ,j  ,k-1) )

                dwdcsi = .5*( w(i+1,j  ,k  ) - w(i-1,j  ,k  ) )
                dwdeta = .5*( w(i  ,j+1,k  ) - w(i  ,j-1,k  ) )
                dwdzita= .5*( w(i  ,j  ,k+1) - w(i  ,j  ,k-1) )

                dudx = dudcsi*apcsx_g1 + dudeta*apetx_g1 + dudzita*apztx_g1
                dudy = dudcsi*apcsy_g1 + dudeta*apety_g1 + dudzita*apzty_g1
                dudz = dudcsi*apcsz_g1 + dudeta*apetz_g1 + dudzita*apztz_g1
                dvdx = dvdcsi*apcsx_g1 + dvdeta*apetx_g1 + dvdzita*apztx_g1
                dvdy = dvdcsi*apcsy_g1 + dvdeta*apety_g1 + dvdzita*apzty_g1
                dvdz = dvdcsi*apcsz_g1 + dvdeta*apetz_g1 + dvdzita*apztz_g1
                dwdx = dwdcsi*apcsx_g1 + dwdeta*apetx_g1 + dwdzita*apztx_g1
                dwdy = dwdcsi*apcsy_g1 + dwdeta*apety_g1 + dwdzita*apzty_g1
                dwdz = dwdcsi*apcsz_g1 + dwdeta*apetz_g1 + dwdzita*apztz_g1
                !
                vortx(i,j,k) = -dwdy+dvdz
                vorty(i,j,k) = -dudz+dwdx
                vortz(i,j,k) = -dvdx+dudy

            end do
        end do
    end do

    return
end

