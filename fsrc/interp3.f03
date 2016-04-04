!***********************************************************************
subroutine interp3(xx1,yy1,zz1,vvf,kparasta,kparaend,alalm,alamm,alalmrho,alammrho)
   !***********************************************************************
   ! interpolazione delle funzioni num e den per modello
   ! lagrangiano in spazio computazionale
   !
   use myarrays_metri3
   use scala3

   implicit none
   !-----------------------------------------------------------------------
   !     variable declaration
   integer i,j,k,inp,jnp,knp,ist,jst,kst,ii
   integer kparasta,kparaend

   double precision alalm(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   double precision alamm(0:n1+1,0:n2+1,kparasta-1:kparaend+1)

   double precision alalmrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   double precision alammrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   !
   real dx1,dx2,dy1,dy2,dz1,dz2
   real epsi,ax,ay,az
   real vv1(4),vv2(4),vv3(4),vv4(4),vv5(4),vv6(4),vv7(4),vv8(4)
   real vv13(4),vv24(4),vv57(4),vv68(4)
   real vv1324(4),vv5768(4)
   real vvf(n1,n2,kparasta:kparaend,4)
   real xx1(n1,n2,kparasta:kparaend)
   real yy1(n1,n2,kparasta:kparaend)
   real zz1(n1,n2,kparasta:kparaend)
   !-----------------------------------------------------------------------
   !
   !     start cycling on plane
   !
   epsi=1.e-12
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            !
            ax=abs(xx1(i,j,k))
            ist=int((xx1(i,j,k)+ax)/(2.*ax+epsi)+.5)
            inp=i-1+ist
            !
            ay=abs(yy1(i,j,k))
            jst=int((yy1(i,j,k)+ay)/(2.*ay+epsi)+.5)
            jnp=j-1+jst
            !
            az=abs(zz1(i,j,k))
            kst=int((zz1(i,j,k)+az)/(2.*az+epsi)+.5)
            knp=k-1+kst
            !
            vv1(1)=alalm(inp,jnp,knp)
            vv1(2)=alamm(inp,jnp,knp)
            vv1(3)=alalmrho(inp,jnp,knp)
            vv1(4)=alammrho(inp,jnp,knp)
            !
            vv2(1)=alalm(inp+1,jnp,knp)
            vv2(2)=alamm(inp+1,jnp,knp)
            vv2(3)=alalmrho(inp+1,jnp,knp)
            vv2(4)=alammrho(inp+1,jnp,knp)

            !
            vv3(1)=alalm(inp,jnp,knp+1)
            vv3(2)=alamm(inp,jnp,knp+1)
            vv3(3)=alalmrho(inp,jnp,knp+1)
            vv3(4)=alammrho(inp,jnp,knp+1)

            !
            vv4(1)=alalm(inp+1,jnp,knp+1)
            vv4(2)=alamm(inp+1,jnp,knp+1)
            vv4(3)=alalmrho(inp+1,jnp,knp+1)
            vv4(4)=alammrho(inp+1,jnp,knp+1)
            !
            vv5(1)=alalm(inp,jnp+1,knp)
            vv5(2)=alamm(inp,jnp+1,knp)
            vv5(3)=alalmrho(inp,jnp+1,knp)
            vv5(4)=alammrho(inp,jnp+1,knp)
            !
            vv6(1)=alalm(inp+1,jnp+1,knp)
            vv6(2)=alamm(inp+1,jnp+1,knp)
            vv6(3)=alalmrho(inp+1,jnp+1,knp)
            vv6(4)=alammrho(inp+1,jnp+1,knp)
            !
            vv7(1)=alalm(inp,jnp+1,knp+1)
            vv7(2)=alamm(inp,jnp+1,knp+1)
            vv7(3)=alalmrho(inp,jnp+1,knp+1)
            vv7(4)=alammrho(inp,jnp+1,knp+1)
            !
            vv8(1)=alalm(inp+1,jnp+1,knp+1)
            vv8(2)=alamm(inp+1,jnp+1,knp+1)
            vv8(3)=alalmrho(inp+1,jnp+1,knp+1)
            vv8(4)=alammrho(inp+1,jnp+1,knp+1)
            !
            ! trova distanze in spazio computazionale
            !
            !      if(xx1(i,j,k).gt.0.) then
            !      dx1=xx1(i,j,k)
            !      dx2=1-xx1(i,j,k)
            !      else
            !      dx1=1.+xx1(i,j,k)
            !      dx2=-xx1(i,j,k)
            !      end if
            !
            dx1=    xx1(i,j,k) *ist + (1.+xx1(i,j,k))*(1-ist)
            dx2=(1.-xx1(i,j,k))*ist + (  -xx1(i,j,k))*(1-ist)
            !
            !      if(yy1(i,j,k).gt.0.) then
            !      dy1=yy1(i,j,k)
            !      dy2=1-yy1(i,j,k)
            !      else
            !      dy1=1.+yy1(i,j,k)
            !      dy2=-yy1(i,j,k)
            !      end if
            !
            dy1=    yy1(i,j,k) *jst + (1.+yy1(i,j,k))*(1-jst)
            dy2=(1.-yy1(i,j,k))*jst + (  -yy1(i,j,k))*(1-jst)
            !
            !      if(zz1(i,j,k).gt.0.) then
            !      dz1=zz1(i,j,k)
            !      dz2=1.-zz1(i,j,k)
            !      else
            !      dz1=1.+zz1(i,j,k)
            !      dz2=-zz1(i,j,k)
            !      end if
            !
            dz1=    zz1(i,j,k) *kst + (1.+zz1(i,j,k))*(1-kst)
            dz2=(1.-zz1(i,j,k))*kst + (  -zz1(i,j,k))*(1-kst)
            !
            ! parte interpolazione
            ! interseca il volume sul piano z=zp
            !

            ! interpola valori delle funzioni sui punti
            !
            do ii=1,4
               vv13(ii)=vv1(ii)*dz2+vv3(ii)*dz1
               vv24(ii)=vv2(ii)*dz2+vv4(ii)*dz1
               vv57(ii)=vv5(ii)*dz2+vv7(ii)*dz1
               vv68(ii)=vv6(ii)*dz2+vv8(ii)*dz1
           
            end do
                                                          
            !
            ! adesso interseco il piano con linea x=xp
            !
            do ii=1,4
               vv5768(ii)=vv57(ii)*dx2+vv68(ii)*dx1
               vv1324(ii)=vv13(ii)*dx2+vv24(ii)*dx1
            end do
            !
            ! infine interpola sul punto
            !
            do ii=1,4
               vvf(i,j,k,ii)=vv1324(ii)*dy2+vv5768(ii)*dy1
            end do
         !
         enddo
      enddo
   enddo
   !
   return
end
