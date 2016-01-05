!----------------------------------------------------------------------
subroutine rhostress(rhostre,pran,x,y,z,rhov)
    !----------------------------------------------------------------------
    ! calcola il gradiente di densità a parete sopra e sotto
    !
    !
    use scala3

    implicit none
    !
    integer isc, i,j,k
    !
    real rhoav(nscal)
    real rhoind(nscal)
    real rhostre(nscal),drhody
    real ak,pran
    real csx,csy,csz
    real etx,ety,etz
    real ztx,zty,ztz
    real giac
    real xcsi,ycsi,zcsi
    real xeta,yeta,zeta
    real xzet,yzet,zzet
    real x0,x1,x2
    real y0,y1,y2
    real z0,z1,z2
    real x(-8:n1+8,-8:n2+8,-8:n3+8)
    real y(-8:n1+8,-8:n2+8,-8:n3+8)
    real z(-8:n1+8,-8:n2+8,-8:n3+8)
    real rhov(nscal,0:n1+1,0:n2+1,0:n3+1)
    !

    ak=1./re/pran
    write(*,*)1./ak,re,pran
    !
    ! parete sotto
    !
    rhoav=0.
    j=0
    do i=1,jx
        do k=1,jz
            !
            ! metrica
            !
            x2=.25*(x(i,j+2,k)+x(i,j+2,k-1)+x(i-1,j+2,k-1)+x(i-1,j+2,k))
            x1=.25*(x(i,j+1,k)+x(i,j+1,k-1)+x(i-1,j+1,k-1)+x(i-1,j+1,k))
            x0=.25*(x(i,j  ,k)+x(i,j  ,k-1)+x(i-1,j  ,k-1)+x(i-1,j  ,k))
            !
            y2=.25*(y(i,j+2,k)+y(i,j+2,k-1)+y(i-1,j+2,k-1)+y(i-1,j+2,k))
            y1=.25*(y(i,j+1,k)+y(i,j+1,k-1)+y(i-1,j+1,k-1)+y(i-1,j+1,k))
            y0=.25*(y(i,j  ,k)+y(i,j  ,k-1)+y(i-1,j  ,k-1)+y(i-1,j  ,k))
            !
            z2=.25*(z(i,j+2,k)+z(i,j+2,k-1)+z(i-1,j+2,k-1)+z(i-1,j+2,k))
            z1=.25*(z(i,j+1,k)+z(i,j+1,k-1)+z(i-1,j+1,k-1)+z(i-1,j+1,k))
            z0=.25*(z(i,j  ,k)+z(i,j  ,k-1)+z(i-1,j  ,k-1)+z(i-1,j  ,k))
            !
            xeta=.5 * ( -3.*x0 + 4.*x1 - x2 )
            yeta=.5 * ( -3.*y0 + 4.*y1 - y2 )
            zeta=.5 * ( -3.*z0 + 4.*z1 - z2 )
            !
            yzet=.5*(y(i,j,k  ) + y(i-1,j,k  ))- &
                .5*(y(i,j,k-1) + y(i-1,j,k-1))
            !
            zcsi=.5*(z(i  ,j,k) + z(i  ,j,k-1))- &
                .5*(z(i-1,j,k) + z(i-1,j,k-1))
            !
            ycsi=.5*(y(i  ,j,k) + y(i  ,j,k-1))- &
                .5*(y(i-1,j,k) + y(i-1,j,k-1))
            !
            zzet=.5*(z(i,j,k  ) + z(i-1,j,k  ))- &
                .5*(z(i,j,k-1) + z(i-1,j,k-1))
            !
            xcsi=.5*(x(i  ,j,k) + x(i  ,j,k-1))- &
                .5*(x(i-1,j,k) + x(i-1,j,k-1))
            !
            xzet=.5*(x(i,j,k  ) + x(i-1,j,k  ))- &
                .5*(x(i,j,k-1) + x(i-1,j,k-1))
            !
            ! termini metrici controvarianti e Jacobiano
            !
            giac= xcsi*(yeta*zzet-yzet*zeta)- &
                xeta*(ycsi*zzet-yzet*zcsi)+ &
                xzet*(ycsi*zeta-yeta*zcsi)
            !
            etx = (yzet*zcsi - ycsi*zzet)/giac
            ety = (xcsi*zzet - xzet*zcsi)/giac
            etz = (xzet*ycsi - xcsi*yzet)/giac
            !
            csx = (yeta*zzet - yzet*zeta)/giac
            csy = (xzet*zeta - xeta*zzet)/giac
            csz = (xeta*yzet - xzet*yeta)/giac
            !
            ztx = (ycsi*zeta-yeta*zcsi)/giac
            zty = (xeta*zcsi-xcsi*zeta)/giac
            ztz = (xcsi*yeta-xeta*ycsi)/giac
            !
            ! calcolo delle derivate
            !
            do isc=1,nscal
                drhody=.5*(rhov(isc,i+1,j,k)-rhov(isc,i-1,j,k))       *csy+ &
                    (-8*rhov(isc,i,j,k)+9.*rhov(isc,i,j+1,k)-rhov(isc,i,j+2,k))/3. &
                    *ety+ &
                    .5*(rhov(isc,i,j,k+1)-rhov(isc,i,j,k-1))                 *zty

                !ccc        print*,i,k,drhody,'--',csy,ety,zty

                rhoav(isc)=rhoav(isc)+ak*drhody
            end do
        end do
    end do
       
    do isc=1,nscal
        rhoav(isc)=rhoav(isc)/float(jx)/float(jz)
    end do
    !
    ! parete sopra
    !
    j=jy
    rhoind=0.
    do i=1,jx
        do k=1,jz
            !
            ! metrica
            !
            x2=.25*(x(i,j-2,k)+x(i,j-2,k-1)+x(i-1,j-2,k-1)+x(i-1,j-2,k))
            x1=.25*(x(i,j-1,k)+x(i,j-1,k-1)+x(i-1,j-1,k-1)+x(i-1,j-1,k))
            x0=.25*(x(i,j  ,k)+x(i,j  ,k-1)+x(i-1,j  ,k-1)+x(i-1,j  ,k))
            !
            y2=.25*(y(i,j-2,k)+y(i,j-2,k-1)+y(i-1,j-2,k-1)+y(i-1,j-2,k))
            y1=.25*(y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k))
            y0=.25*(y(i,j  ,k)+y(i,j  ,k-1)+y(i-1,j  ,k-1)+y(i-1,j  ,k))
            !
            z2=.25*(z(i,j-2,k)+z(i,j-2,k-1)+z(i-1,j-2,k-1)+z(i-1,j-2,k))
            z1=.25*(z(i,j-1,k)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j-1,k))
            z0=.25*(z(i,j  ,k)+z(i,j  ,k-1)+z(i-1,j  ,k-1)+z(i-1,j  ,k))
            !
            xeta=.5 * ( 3.*x0 - 4.*x1 + x2 )
            yeta=.5 * ( 3.*y0 - 4.*y1 + y2 )
            zeta=.5 * ( 3.*z0 - 4.*z1 + z2 )
            !
            yzet=.5*(y(i,j,k  ) + y(i-1,j,k  ))- &
                .5*(y(i,j,k-1) + y(i-1,j,k-1))
            !
            zcsi=.5*(z(i  ,j,k) + z(i  ,j,k-1))- &
                .5*(z(i-1,j,k) + z(i-1,j,k-1))
            !
            ycsi=.5*(y(i  ,j,k) + y(i  ,j,k-1))- &
                .5*(y(i-1,j,k) + y(i-1,j,k-1))
            !
            zzet=.5*(z(i,j,k  ) + z(i-1,j,k  ))- &
                .5*(z(i,j,k-1) + z(i-1,j,k-1))
            !
            xcsi=.5*(x(i  ,j,k) + x(i  ,j,k-1))- &
                .5*(x(i-1,j,k) + x(i-1,j,k-1))
            !
            xzet=.5*(x(i,j,k  ) + x(i-1,j,k  ))- &
                .5*(x(i,j,k-1) + x(i-1,j,k-1))
            !
            ! termini metrici controvarianti e Jacobiano
            !
            giac= xcsi*(yeta*zzet-yzet*zeta)- &
                xeta*(ycsi*zzet-yzet*zcsi)+ &
                xzet*(ycsi*zeta-yeta*zcsi)
            !
            etx = (yzet*zcsi - ycsi*zzet)/giac
            ety = (xcsi*zzet - xzet*zcsi)/giac
            etz = (xzet*ycsi - xcsi*yzet)/giac
            !
            csx = (yeta*zzet - yzet*zeta)/giac
            csy = (xzet*zeta - xeta*zzet)/giac
            csz = (xeta*yzet - xzet*yeta)/giac
            !
            ztx = (ycsi*zeta-yeta*zcsi)/giac
            zty = (xeta*zcsi-xcsi*zeta)/giac
            ztz = (xcsi*yeta-xeta*ycsi)/giac
            !
            ! calcolo delle derivate
            !
            do isc=1,nscal
                drhody=.5*(rhov(isc,i+1,j+1,k)-rhov(isc,i-1,j+1,k))*csy+ &
                    (8.*rhov(isc,i,j+1,k)-9.*rhov(isc,i,j,k)+rhov(isc,i,j-1,k))/3.*ety &
                    +5*(rhov(isc,i,j+1,k+1)-rhov(isc,i,j+1,k-1))           *zty

                rhoind(isc)=rhoind(isc)+ak*drhody
            end do
        end do
    end do

    do isc=1,nscal
        rhoind(isc)=rhoind(isc)/float(jx)/float(jz)
        rhostre(isc)=abs(rhoav(isc)+rhoind(isc))/2.
    end do
    !
    return
end
