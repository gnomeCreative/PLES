!************************************************************************
subroutine contour(kparasta,kparaend,nproc,myid)
   !***********************************************************************
   ! set boundary condition for cartesian velocity and controvariant flux
   !
   use mysettings
   use mysettings_boundary
   use myarrays_velo3
   use myarrays_metri3
   use myarrays_WB
   use myarrays_nesting
   !
   use scala3
   use period
   use orl
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,ii,jj,kk,isc
   integer ierr,myid,nproc
   integer kparasta,kparaend
   integer kparastal,kparaendl
   !      integer windyes,ktime,i_rest

   real theta(jx)
   real upp,delta,yc

   real delmassa1
   real delmassa2
   real delmassa3
   real delmassa4
   real delmassa5
   real delmassa6
      
   real segno_parete,segno_pre
   real segno_u,segno_v,segno_w
   real u_int,v_int,w_int
   real ucp_old,wcp_old
   !
   !-----------------------------------------------------------------------
   if(i_rest==3)then
      if(infout1 /= 0)iboun1 = 2
      if(infout2 /= 0)iboun2 = 2
      if(infout5 /= 0)iboun5 = 2
      if(infout6 /= 0)iboun6 = 2
   end if
   !-----------------------------------------------------------------------
   ! (for periodicity see vel_up)
   !chicco portare right and left dentro un blocco do unico
   !     left side
   do ii=1,ip
      !
      if(infout1.eq.0)then        !inflow
         do k=kparasta,kparaend
            do j=1,jy

               u(0,j,k)=up1(j,k)
               v(0,j,k)=vp1(j,k)
               w(0,j,k)=wp1(j,k)
               do isc=1,nscal
                  rhov(isc,0,j,k)=rhovp1(isc,j,k)
               end do
               uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)


            enddo
         enddo

      elseif(infout1.eq.1)then   !outflow
         do k=kparasta,kparaend
            do j=1,jy

               uc(0,j,k) = uc1_orl(j,k)

               delmassa1= uc(0,j,k)-(du_dx1(j,k)*csx(0,j,k)+dv_dx1(j,k)*csy(0,j,k)+dw_dx1(j,k)*csz(0,j,k))

               u(0,j,k)= du_dx1(j,k) + delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
               v(0,j,k)= dv_dx1(j,k) + delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
               w(0,j,k)= dw_dx1(j,k) + delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
	
               do isc=1,nscal
                  rhov(isc,0,j,k)=rhov(isc,1,j,k)
               end do
	
            enddo
         enddo

      elseif(infout1.eq.2)then      ! wall
         if(iboun1==0)then
            do k=kparasta,kparaend
               do j=1,jy
                  u(0,j,k)=0.
                  v(0,j,k)=0.
                  w(0,j,k)=0.
                  do isc=1,nscal
                     rhov(isc,0,j,k)=rhov(isc,1,j,k)
                  end do
                  uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
               enddo
            enddo
         elseif(iboun1==1)then
            do k=kparasta,kparaend
               do j=1,jy
                  u(0,j,k)=0.
                  v(0,j,k)=v(1,j,k)
                  w(0,j,k)=w(1,j,k)
                  do isc=1,nscal
                     rhov(isc,0,j,k)=rhov(isc,1,j,k)
                  end do
                  uc(0,j,k)=csx(0,j,k)*u(0,j,k)+csy(0,j,k)*v(0,j,k)+csz(0,j,k)*w(0,j,k)
               enddo
            enddo
         elseif(iboun1==2)then
            do k=kparasta,kparaend
               do j=1,jy



                  if( ucp1(j,k) .gt. 0.)then

                     ucp_old = csx(0,j,k)*up1(j,k) + csy(0,j,k)*vp1(j,k) + csz(0,j,k)*wp1(j,k)

                     delmassa1 = ucp1(j,k) - ucp_old

                     u(0,j,k)= up1(j,k) + delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                     !             v(0,j,k)=v(1,j,k)
                     v(0,j,k)= vp1(j,k) + delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
     
                     w(0,j,k)= wp1(j,k) + delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                  else

                     ucp_old = csx(0,j,k)*u(1,j,k)+ csy(0,j,k)*v(1,j,k)+ csz(0,j,k)*w(1,j,k)

                     delmassa1 = ucp1(j,k) - ucp_old


                     u(0,j,k)= u(1,j,k) + delmassa1*csx(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)

                     !             v(0,j,k)=v(1,j,k)
                     v(0,j,k)= v(1,j,k) + delmassa1*csy(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
     
                     w(0,j,k)= w(1,j,k) + delmassa1*csz(0,j,k)/(csx(0,j,k)**2+csy(0,j,k)**2+csz(0,j,k)**2)
                  end if


                  uc(0,j,k)=ucp1(j,k)

                  do isc=1,nscal
                     if(potenziale==0)then
                        rhov(isc,0,j,k)=rhovp1(isc,j,k)
                        if(ucp1(j,k).lt.0.) rhov(isc,0,j,k)=rhov(isc,1,j,k)

                     else
                        rhov(isc,0,j,k)=rhovo1(isc,j,k)
                     end if
                  end do
	     
               enddo
            enddo
         end if
      end if

      !     right side
      if(infout2.eq.0)then        !inflow
         do k=kparasta,kparaend
            do j=1,jy

               u(jx+1,j,k)=up2(j,k)
               v(jx+1,j,k)=vp2(j,k)
               w(jx+1,j,k)=wp2(j,k)
               do isc=1,nscal
                  rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
               end do
               uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
            enddo
         enddo
      elseif(infout2.eq.1)then   !outflow
         do k=kparasta,kparaend
            do j=1,jy
	
               uc(jx,j,k) = uc2_orl(j,k)
		
               delmassa2= uc(jx,j,k)-(du_dx2(j,k)*csx(jx,j,k)+dv_dx2(j,k)*csy(jx,j,k)+dw_dx2(j,k)*csz(jx,j,k))

               u(jx+1,j,k)= du_dx2(j,k) + delmassa2*csx(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
               v(jx+1,j,k)= dv_dx2(j,k) + delmassa2*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
               w(jx+1,j,k)= dw_dx2(j,k) + delmassa2*csz(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
               do isc=1,nscal
                  rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
               end do
	

            enddo
         enddo

      elseif(infout2.eq.2)then      !wall
         if(iboun2==0)then
            do k=kparasta,kparaend
               do j=1,jy
                  u(jx+1,j,k)=0.
                  v(jx+1,j,k)=0.
                  w(jx+1,j,k)=0.
                  do isc=1,nscal
                     rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                  end do
                  uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
               enddo
            enddo
         elseif(iboun2==1)then
            do k=kparasta,kparaend
               do j=1,jy
                  u(jx+1,j,k)=0.
                  v(jx+1,j,k)=v(jx,j,k)
                  w(jx+1,j,k)=w(jx,j,k)
                  do isc=1,nscal
                     rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                  end do
                  uc(jx,j,k)=csx(jx,j,k)*u(jx+1,j,k)+csy(jx,j,k)*v(jx+1,j,k)+csz(jx,j,k)*w(jx+1,j,k)
               enddo
            enddo
         elseif(iboun2==2)then
            do k=kparasta,kparaend
               do j=1,jy

                  if(ucp2(j,k) .lt. 0.)then
            
                     ucp_old = csx(jx,j,k)*up2(j,k)+ csy(jx,j,k)*vp2(j,k)+ csz(jx,j,k)*wp2(j,k)

                     delmassa2 = ucp2(j,k) - ucp_old

                     u(jx+1,j,k)= up2(j,k) + delmassa2*csx(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                     !             v(jx+1,j,k)=v(1,j,k)
                     v(jx+1,j,k)= vp2(j,k) + delmassa2*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
     
                     w(jx+1,j,k)= wp2(j,k) + delmassa2*csz(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
	     
                     uc(jx,j,k)=ucp2(j,k)
	     
                  else
	     
                     ucp_old = csx(jx,j,k)*u(jx,j,k)+ csy(jx,j,k)*v(jx,j,k) + csz(jx,j,k)*w(jx,j,k)

                     delmassa2 = ucp2(j,k) - ucp_old

                     u(jx+1,j,k)= u(jx,j,k) + delmassa2*csx(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)

                     !             v(jx+1,j,k)=v(1,j,k)
                     v(jx+1,j,k)= v(jx,j,k) + delmassa2*csy(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
     
                     w(jx+1,j,k)= w(jx,j,k) + delmassa2*csz(jx,j,k)/(csx(jx,j,k)**2+csy(jx,j,k)**2+csz(jx,j,k)**2)
	     
                     uc(jx,j,k)=ucp2(j,k)
	     
                  end if

                  do isc=1,nscal
                     if(potenziale==0)then
                        rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
                        if( ucp2(j,k) .gt. 0.)then
                           rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
                        end if

                     else
                        rhov(isc,jx+1,j,k)=rhovo2(isc,j,k)
                     end if
                  end do

               enddo
            enddo
         end if
      end if

   enddo   !end loop  ii=1,ip
   !-----------------------------------------------------------------------
   !
   !     bottom side
   do jj=1,jp

      if(infout3.eq.0)then             !inflow
         do k=kparasta,kparaend
            do i=1,jx

               u(i,0,k)=up3(i,k)
               v(i,0,k)=vp3(i,k)
               w(i,0,k)=wp3(i,k)
               do isc=1,nscal
                  rhov(isc,i,0,k)=rhovp3(isc,i,k)
               end do
               vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
            enddo
         enddo
      elseif(infout3.eq.1)then        !outflow
         do k=kparasta,kparaend
            do i=1,jx

               vc(i,0,k) = vc3_orl(i,k)

               delmassa3= vc(i,0,k)-(du_dy3(i,k)*etx(i,0,k)+dv_dy3(i,k)*ety(i,0,k)+dw_dy3(i,k)*etz(i,0,k))

               u(i,0,k)=  du_dy3(i,k)+delmassa3*etx(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)
               v(i,0,k)=  dv_dy3(i,k)+delmassa3*ety(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)
               w(i,0,k)=  dw_dy3(i,k)+delmassa3*etz(i,0,k)/(etx(i,0,k)**2+ety(i,0,k)**2+etz(i,0,k)**2)

               do isc=1,nscal
                  rhov(isc,i,0,k)=rhov(isc,i,1,k)
               end do


            enddo
         enddo

      elseif(infout3.eq.2)then           !wall
         if(iboun3==0)then
            do k=kparasta,kparaend
               do i=1,jx
                  u(i,0,k)=0.
                  v(i,0,k)=0.
                  w(i,0,k)=0.
                  do isc=1,nscal
                     rhov(isc,i,0,k)=rhov(isc,i,1,k)
                  end do
                  vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
               enddo
            enddo
         elseif(iboun3==1)then
            do k=kparasta,kparaend
               do i=1,jx
                  u(i,0,k)=u(i,1,k)
                  v(i,0,k)=0.
                  w(i,0,k)=w(i,1,k)
                  do isc=1,nscal
                     rhov(isc,i,0,k)=rhov(isc,i,1,k)
                  end do
                  vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
               enddo
            enddo
         elseif(iboun3==2)then
            do k=kparasta,kparaend
               do i=1,jx
                  u(i,0,k)=0.
                  v(i,0,k)=0.
                  w(i,0,k)=0.
                  do isc=1,nscal
                     rhov(isc,i,0,k)=rhov(isc,i,1,k)
                  end do
                  vc(i,0,k)=etx(i,0,k)*u(i,0,k)+ety(i,0,k)*v(i,0,k)+etz(i,0,k)*w(i,0,k)
               enddo
            enddo
         end if
      end if
      !
      !     upper side
      if(freesurface == 0)then !no freesurface<<<<<<<<<<<<<<<<<<<<<<<<<<

         if(infout4.eq.0)then             !inflow
            do k=kparasta,kparaend
               do i=1,jx

                  u(i,jy+1,k)=up4(i,k)
                  v(i,jy+1,k)=vp4(i,k)
                  w(i,jy+1,k)=wp4(i,k)
                  do isc=1,nscal
                     rhov(isc,i,jy+1,k)=rhovp4(isc,i,k)
                  end do
                  vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
               enddo
            enddo
         elseif(infout4.eq.1)then        !outflow
            do k=kparasta,kparaend
               do i=1,jx

                  vc(i,jy,k) = vc4_orl(i,k)
			
                  delmassa4= vc(i,jy,k)-(du_dy4(i,k)*etx(i,jy,k)+dv_dy4(i,k)*ety(i,jy,k)+dw_dy4(i,k)*etz(i,jy,k))

                  u(i,jy+1,k)= du_dy4(i,k)+delmassa4*etx(i,jy,k)/(etx(i,jy,k)**2+ety(i,jy,k)**2+etz(i,jy,k)**2)
                  v(i,jy+1,k)= dv_dy4(i,k)+delmassa4*ety(i,jy,k)/(etx(i,jy,k)**2+ety(i,jy,k)**2+etz(i,jy,k)**2)
                  w(i,jy+1,k)= dw_dy4(i,k)+delmassa4*etz(i,jy,k)/(etx(i,jy,k)**2+ety(i,jy,k)**2+etz(i,jy,k)**2)

                  do isc=1,nscal
                     rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                  end do
	
	

               enddo
            enddo

         ! windyes=1 option, if at the upper side there is a wind stress

         elseif(infout4.eq.2 .and. windyes .eq. 0)then ! wall without wind
            if(iboun4==0)then
               do k=kparasta,kparaend
                  do i=1,jx
                     u(i,jy+1,k)=0.
                     v(i,jy+1,k)=0.
                     w(i,jy+1,k)=0.
                     do isc=1,nscal
                        rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                     end do
                     vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k) +ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                  enddo
               enddo
            elseif(iboun4==1)then
               do k=kparasta,kparaend
                  do i=1,jx
                     u(i,jy+1,k)=u(i,jy,k)
                     v(i,jy+1,k)=0.
                     w(i,jy+1,k)=w(i,jy,k)
                     do isc=1,nscal
                        rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                     end do
                     vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                     if((i.eq.4).and.(k.eq.17))then
                        write(*,*)'VC_contour:',vc(i,32,k)
                     endif
                  enddo
               enddo
            elseif(iboun4==2)then
               do k=kparasta,kparaend
                  do i=1,jx
                     u(i,jy+1,k)=0.
                     v(i,jy+1,k)=0.
                     w(i,jy+1,k)=0.
                     do isc=1,nscal
                        rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                     end do
                     vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
                  enddo
               enddo
            end if

         elseif(infout4.eq.2 .and. windyes .eq. 1)then ! wall with wind extrapolation from the interior
    
            do k=kparasta,kparaend
               do i=1,jx
                  u(i,jy+1,k) = 1.5*u(i,jy,k)-.5*u(i,jy-1,k)
                  v(i,jy+1,k)= 0.0
                  w(i,jy+1,k) = 1.5*w(i,jy,k)-.5*w(i,jy-1,k)
                  do isc=1,nscal
                     rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                  end do
                  vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)+ety(i,jy,k)*v(i,jy+1,k)+etz(i,jy,k)*w(i,jy+1,k)
               enddo
            enddo
         end if

      elseif(freesurface == 1)then !free surface ON.<<<<<<<<<
         if(windyes .eq. 0)then
            do k=kparasta,kparaend
               do i=1,jx
                  u(i,jy+1,k)=u(i,jy,k)
                  v(i,jy+1,k)=v(i,jy,k)
                  w(i,jy+1,k)=w(i,jy,k)
                  do isc=1,nscal
                     rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                  end do
               !           vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)
               !     >               +ety(i,jy,k)*v(i,jy+1,k)
               !     >               +etz(i,jy,k)*w(i,jy+1,k)
               enddo
            enddo
         elseif(windyes .eq. 1)then
            do k=kparasta,kparaend
               do i=1,jx
                  u(i,jy+1,k) = 1.5*u(i,jy,k)-.5*u(i,jy-1,k)
                  v(i,jy+1,k)= v(i,jy,k)
                  w(i,jy+1,k) = 1.5*w(i,jy,k)-.5*w(i,jy-1,k)
                  do isc=1,nscal
                     rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                  end do
               !           vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k)
               !     >               +ety(i,jy,k)*v(i,jy+1,k)
               !     >               +etz(i,jy,k)*w(i,jy+1,k)
               enddo
            enddo
         end if
      end if !free surface<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   enddo !end loop  jj=1,jp
   !-----------------------------------------------------------------------
   !
   !
   do kk=1,kp
      !
      !     back side
      if (myid.eq.0) then
         if(infout5.eq.0)then            !inflow
            do j=1,jy
               do i=1,jx

                  u(i,j,0)=up5(i,j)
                  v(i,j,0)=vp5(i,j)
                  w(i,j,0)=wp5(i,j)
                  do isc=1,nscal
                     rhov(isc,i,j,0)=rhovp5(isc,i,j)
                  end do
                  wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)

               enddo
            enddo
         elseif(infout5.eq.1)then            !outflow
            do j=1,jy
               do i=1,jx

                  wc(i,j,0) = wc5_orl(i,j)

                  delmassa5= wc(i,j,0)-(du_dz5(i,j)*ztx(i,j,0)+dv_dz5(i,j)*zty(i,j,0)+dw_dz5(i,j)*ztz(i,j,0))

                  u(i,j,0)= du_dz5(i,j)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
                  v(i,j,0)= dv_dz5(i,j)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
                  w(i,j,0)= dw_dz5(i,j)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                  do isc=1,nscal
                     rhov(isc,i,j,0)=rhov(isc,i,j,1)
                  end do
	

               enddo
            enddo


         elseif(infout5.eq.2)then            !wall
            if(iboun5==0)then
               do j=1,jy
                  do i=1,jx
                     u(i,j,0)=0.
                     v(i,j,0)=0.
                     w(i,j,0)=0.
                     do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                     end do
                     wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                  enddo
               enddo
            elseif(iboun5==1)then
               do j=1,jy
                  do i=1,jx
                     u(i,j,0)=u(i,j,1)
                     v(i,j,0)=v(i,j,1)
                     w(i,j,0)=0.
                     do isc=1,nscal
                        rhov(isc,i,j,0)=rhov(isc,i,j,1)
                     end do
                     wc(i,j,0)=ztx(i,j,0)*u(i,j,0)+zty(i,j,0)*v(i,j,0)+ztz(i,j,0)*w(i,j,0)
                  enddo
               enddo
            elseif(iboun5==2)then
               do j=1,jy
                  do i=1,jx

                     if( wcp5(i,j) .gt. 0.)then

                        wcp_old =ztx(i,j,0)*up5(i,j)+zty(i,j,0)*vp5(i,j)+ztz(i,j,0)*wp5(i,j)

                        delmassa5 = wcp5(i,j) - wcp_old

                        u(i,j,0)= up5(i,j)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                        !             v(i,j,0)=v(i,j,1)
                        v(i,j,0)= vp5(i,j)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                        w(i,j,0)= wp5(i,j)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
	        
                        wc(i,j,0)=wcp5(i,j)
	     
                     else
	     
                        wcp_old =ztx(i,j,0)*u(i,j,1)+zty(i,j,0)*v(i,j,1)+ztz(i,j,0)*w(i,j,1)

                        delmassa5 = wcp5(i,j) - wcp_old

                        u(i,j,0)= u(i,j,1)+delmassa5*ztx(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                        !             v(i,j,0)=v(i,j,1)
                        v(i,j,0)= v(i,j,1)+delmassa5*zty(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)

                        w(i,j,0)= w(i,j,1)+delmassa5*ztz(i,j,0)/(ztx(i,j,0)**2+zty(i,j,0)**2+ztz(i,j,0)**2)
	        
                        wc(i,j,0)=wcp5(i,j)
	     	     
                     end if
	     
                     do isc=1,nscal
                        if(potenziale==0)then
                           rhov(isc,i,j,0)=rhovp5(isc,i,j)
                           if(wcp5(i,j).lt.0.)then
                              rhov(isc,i,j,0)=rhov(isc,i,j,1)
                           end if

                        else
                           rhov(isc,i,j,0)=rhovo5(isc,i,j)
                        end if
                     end do

                  enddo
               enddo
            end if
         endif
      endif
      !
      !     front side
      if (myid.eq.nproc-1) then
         if(infout6.eq.0)then            !inflow
            do j=1,jy
               do i=1,jx

                  u(i,j,jz+1)=up6(i,j)
                  v(i,j,jz+1)=vp6(i,j)
                  w(i,j,jz+1)=wp6(i,j)
                  do isc=1,nscal
                     rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)
                  end do
                  wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)

               enddo
            enddo
         elseif(infout6.eq.1)then            !outflow
            do j=1,jy
               do i=1,jx

                  wc(i,j,jz) = wc6_orl(i,j)

                  delmassa6= wc(i,j,jz)-(du_dz6(i,j)*ztx(i,j,jz) +dv_dz6(i,j)*zty(i,j,jz)+dw_dz6(i,j)*ztz(i,j,jz))

                  u(i,j,jz+1)= du_dz6(i,j)+delmassa6*ztx(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)
                  v(i,j,jz+1)= dv_dz6(i,j)+delmassa6*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)
                  w(i,j,jz+1)= dw_dz6(i,j)+delmassa6*ztz(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                  do isc=1,nscal
                     rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                  end do

               enddo
            enddo


         elseif(infout6.eq.2)then            !wall
            if(iboun6==0)then
               do j=1,jy
                  do i=1,jx
                     u(i,j,jz+1)=0.
                     v(i,j,jz+1)=0.
                     w(i,j,jz+1)=0.
                     do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                     end do
                     wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                  enddo
               enddo
            elseif(iboun6==1)then
               do j=1,jy
                  do i=1,jx
                     u(i,j,jz+1)=u(i,j,jz)
                     v(i,j,jz+1)=v(i,j,jz)
                     w(i,j,jz+1)=0.



                     do isc=1,nscal
                        rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                     end do
                     wc(i,j,jz)=ztx(i,j,jz)*u(i,j,jz+1)+zty(i,j,jz)*v(i,j,jz+1)+ztz(i,j,jz)*w(i,j,jz+1)
                  enddo
               enddo
            elseif(iboun6==2)then
               do j=1,jy
                  do i=1,jx

                     if(wcp6(i,j) .lt. 0.)then
                        wcp_old =ztx(i,j,jz)*up6(i,j)+zty(i,j,jz)*vp6(i,j)+ztz(i,j,jz)*wp6(i,j)

                        delmassa6 = wcp6(i,j) - wcp_old

                        u(i,j,jz+1)= up6(i,j)+delmassa6*ztx(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                        !             v(i,j,jz+1)=v(i,j,jz)
                        v(i,j,jz+1)= vp6(i,j)+delmassa6*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                        w(i,j,jz+1)= wp6(i,j)+delmassa6*ztz(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                        wc(i,j,jz)=wcp6(i,j)
	     
                     else

                        wcp_old =ztx(i,j,jz)*u(i,j,jz)+zty(i,j,jz)*v(i,j,jz)+ztz(i,j,jz)*w(i,j,jz)

                        delmassa6 = wcp6(i,j) - wcp_old

                        u(i,j,jz+1)= u(i,j,jz)+delmassa6*ztx(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                        !             v(i,j,jz+1)=v(i,j,jz)
                        v(i,j,jz+1)= v(i,j,jz)+delmassa6*zty(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                        w(i,j,jz+1)= w(i,j,jz)+delmassa6*ztz(i,j,jz)/(ztx(i,j,jz)**2+zty(i,j,jz)**2+ztz(i,j,jz)**2)

                        wc(i,j,jz)=wcp6(i,j)
	     	     
                     end if
	     
                     do isc=1,nscal
                        if(potenziale==0)then
                           rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)

                           if(wcp6(i,j).gt.0.)then
                              rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                           end if
                        else
                           rhov(isc,i,j,jz+1)= rhovo6(isc,i,j)
                        end if
                     end do

                  enddo
               enddo
            end if
         endif
      endif

   enddo !end loop  kk=1,kp
   !-----------------------------------------------------------------------
   !     set boundary condition at the corner

   if(myid.eq.0)then
      kparastal=0
      kparaendl=kparaend
   else if (myid.eq.nproc-1) then
      kparastal=kparasta
      kparaendl=kparaend+1
   else
      kparastal=kparasta
      kparaendl=kparaend
   endif

   ! FIRST: corner with two sides like orlansky (infout=1)

   !     corner close to side 5
   if (myid .eq. 0) then
      
      !     corner 5/1
      if (infout5 .eq. 1 .and. infout1 .eq. 1 ) then
         do j=0,jy+1
            u(0,j,0)=u(1,j,1)
            v(0,j,0)=v(1,j,1)
            w(0,j,0)=w(1,j,1)
            do isc=1,nscal
               rhov(isc,0,j,0)=rhov(isc,1,j,1)
            end do
         end do
      end if  !corner 5/1
      
      !     corner 5/2
      if (infout5 .eq. 1 .and. infout2 .eq. 1 ) then
         do j=0,jy+1
            u(jx+1,j,0)=u(jx,j,1)
            v(jx+1,j,0)=v(jx,j,1)
            w(jx+1,j,0)=w(jx,j,1)
            do isc=1,nscal
               rhov(isc,jx+1,j,0)=rhov(isc,jx,j,1)
            end do
         end do
      end if  !corner 5/2
      
      !     corner 5/3
      if (infout5 .eq. 1 .and. infout3 .eq. 1 ) then
         do i=0,jx+1
            u(i,0,0)=u(i,1,1)
            v(i,0,0)=v(i,1,1)
            w(i,0,0)=w(i,1,1)
            do isc=1,nscal
               rhov(isc,i,0,0)=rhov(isc,i,1,1)
            end do
         end do
      end if  !corner 5/3
      
      !     corner 5/4
      if (infout5 .eq. 1 .and. infout4 .eq. 1 ) then
         do i=0,jx+1
            u(i,jy+1,0)=u(i,jy,1)
            v(i,jy+1,0)=v(i,jy,1)
            w(i,jy+1,0)=w(i,jy,1)
            do isc=1,nscal
               rhov(isc,i,jy+1,0)=rhov(isc,i,jy,1)
            end do
         end do
      end if  !corner 5/4
      
   end if !proc 0
   !-----------------------------------------------------------------------
   !     corner close to side 6
   if (myid .eq. nproc-1) then
      !     corner 6/1
      if (infout6 .eq. 1 .and. infout1 .eq. 1 ) then
         do j=0,jy+1
            u(0,j,jz+1)=u(1,j,jz)
            v(0,j,jz+1)=v(1,j,jz)
            w(0,j,jz+1)=w(1,j,jz)
            do isc=1,nscal
               rhov(isc,0,j,jz+1)=rhov(isc,1,j,jz)
            end do
         end do
      end if  !corner 6/1
      
      !     corner 6/2
      if (infout6 .eq. 1 .and. infout2 .eq. 1 ) then
         do j=0,jy+1
            u(jx+1,j,jz+1)=u(jx,j,jz)
            v(jx+1,j,jz+1)=v(jx,j,jz)
            w(jx+1,j,jz+1)=w(jx,j,jz)
            do isc=1,nscal
               rhov(isc,jx+1,j,jz+1)=rhov(isc,jx,j,jz)
            end do
         end do
      end if  !corner 6/2
      
      !     corner 6/3
      if (infout6 .eq. 1 .and. infout3 .eq. 1 ) then
         do i=0,jx+1
            u(i,0,jz+1)=u(i,1,jz)
            v(i,0,jz+1)=v(i,1,jz)
            w(i,0,jz+1)=w(i,1,jz)
            do isc=1,nscal
               rhov(isc,i,0,jz+1)=rhov(isc,i,1,jz)
            end do
         end do
      end if  !corner 6/3
      
      !     corner 6/4
      if (infout6 .eq. 1 .and. infout4 .eq. 1 ) then
         do i=0,jx+1
            u(i,jy+1,jz+1)=u(i,jy,jz)
            v(i,jy+1,jz+1)=v(i,jy,jz)
            w(i,jy+1,jz+1)=w(i,jy,jz)
            do isc=1,nscal
               rhov(isc,i,jy+1,jz+1)=rhov(isc,i,jy,jz)
            end do
         end do
      end if  !corner 6/4

   end if ! nproc-1

   !-----------------------------------------------------------------------
   ! corner 1/3 e 1/4

   !     corner 1/3
   if (infout1 .eq. 1 .and. infout3 .eq. 1 ) then
      do k=kparastal,kparaendl
         u(0,0,k)=u(1,1,k)
         v(0,0,k)=v(1,1,k)
         w(0,0,k)=w(1,1,k)
         do isc=1,nscal
            rhov(isc,0,0,k)=rhov(isc,1,1,k)
         end do
      end do
   end if  !corner 1/3

   !     corner 1/4
   if (infout1 .eq. 1 .and. infout4 .eq. 1 ) then
      do k=kparastal,kparaendl
         u(0,jy+1,k)=u(1,jy,k)
         v(0,jy+1,k)=v(1,jy,k)
         w(0,jy+1,k)=w(1,jy,k)
         do isc=1,nscal
            rhov(isc,0,jy+1,k)=rhov(isc,1,jy,k)
         end do
      end do
   end if  !corner 1/4

   !-----------------------------------------------------------------------
   ! corner 2/3 e 2/4
   !     corner 2/3
   if (infout2 .eq. 1 .and. infout3 .eq. 1 ) then
      do k=kparastal,kparaendl
         u(jx+1,0,k)=u(jx,1,k)
         v(jx+1,0,k)=v(jx,1,k)
         w(jx+1,0,k)=w(jx,1,k)
         do isc=1,nscal
            rhov(isc,jx+1,0,k)=rhov(isc,jx,1,k)
         end do
      end do
   end if  !corner 2/3
      
   !     corner 2/4
   if (infout2 .eq. 1 .and. infout4 .eq. 1 ) then
      do k=kparastal,kparaendl
         u(jx+1,jy+1,k)=u(jx,jy,k)
         v(jx+1,jy+1,k)=v(jx,jy,k)
         w(jx+1,jy+1,k)=w(jx,jy,k)
         do isc=1,nscal
            rhov(isc,jx+1,jy+1,k)=rhov(isc,jx,jy,k)
         end do
      end do
   end if  !corner 2/4
   !-----------------------------------------------------------------------
   ! SECOND: case sides with inlow (infout=0)

   !     side 1
   if (infout1 .eq. 0) then
      i=0
      do k=kparastal,kparaendl
         !     corner 1/3
         u(i,0,k)=u(i,1,k)
         v(i,0,k)=v(i,1,k)
         w(i,0,k)=w(i,1,k)
         do isc=1,nscal
            rhov(isc,i,0,k)=rhov(isc,i,1,k)
         end do
      
         !     corner 1/4
         u(i,jy+1,k)=u(i,jy,k)
         v(i,jy+1,k)=v(i,jy,k)
         w(i,jy+1,k)=w(i,jy,k)
         do isc=1,nscal
            rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
         end do
      end do

      do j=0,jy+1
         !     corner 1/5
         if(myid.eq.0)then
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !     corner 1/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            end do
         end if
      end do

   end if !side 1
      
   !-----------------------------------------------------------------------
   !     side 2
   if (infout2 .eq. 0) then
      i=jx+1
      do k=kparastal,kparaendl
         !     corner 2/3
         u(i,0,k)=u(i,1,k)
         v(i,0,k)=v(i,1,k)
         w(i,0,k)=w(i,1,k)
         do isc=1,1,nscal
            rhov(isc,i,0,k)=rhov(isc,i,1,k)
         end do

         !     corner 2/4
         u(i,jy+1,k)=u(i,jy,k)
         v(i,jy+1,k)=v(i,jy,k)
         w(i,jy+1,k)=w(i,jy,k)
         do isc=1,nscal
            rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
         end do
      end do

      do j=0,jy+1
         if(myid.eq.0)then
            !     corner 2/5
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !     corner 2/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            end do
         end if
      end do

   end if !side 2
      
   !-----------------------------------------------------------------------
   !     side 3
   if (infout3 .eq. 0) then
      j=0
      do k=kparastal,kparaendl
         !     corner 3/1
         u(0,j,k)=u(1,j,k)
         v(0,j,k)=v(1,j,k)
         w(0,j,k)=w(1,j,k)
         do isc=1,nscal
            rhov(isc,0,j,k)=rhov(isc,1,j,k)
         end do

         !     corner 3/2
         u(jx+1,j,k)=u(jx,j,k)
         v(jx+1,j,k)=v(jx,j,k)
         w(jx+1,j,k)=w(jx,j,k)
         do isc=1,nscal
            rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
         end do
      end do

      do i=0,jx+1
         if(myid.eq.0)then
            !     corner 3/5
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !     corner 3/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            end do
         end if
      end do

   end if !side 3
      
   !-----------------------------------------------------------------------
   !     side 4
   if (infout4 .eq. 0) then
      j=jy+1
      do k=kparastal,kparaendl
         !     corner 4/1
         u(0,j,k)=u(1,j,k)
         v(0,j,k)=v(1,j,k)
         w(0,j,k)=w(1,j,k)
         do isc=1,nscal
            rhov(isc,0,j,k)=rhov(isc,1,j,k)
         end do

         !     corner 4/2
         u(jx+1,j,k)=u(jx,j,k)
         v(jx+1,j,k)=v(jx,j,k)
         w(jx+1,j,k)=w(jx,j,k)
         do isc=1,nscal
            rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
         end do
      end do

      do i=0,jx+1
         if(myid.eq.0)then
            !     corner 4/5
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !     corner 4/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            end do
         end if
      end do

   end if !side 4
      
   !-----------------------------------------------------------------------
   !     side 5
   if (infout5 .eq. 0) then
      if(myid.eq.0)then
         k=0
         do j=0,jy+1
            !     corner 5/1
            u(0,j,k)=u(1,j,k)
            v(0,j,k)=v(1,j,k)
            w(0,j,k)=w(1,j,k)
            do isc=1,nscal
               rhov(isc,0,j,k)=rhov(isc,1,j,k)
            end do

            !     corner 5/2
            u(jx+1,j,k)=u(jx,j,k)
            v(jx+1,j,k)=v(jx,j,k)
            w(jx+1,j,k)=w(jx,j,k)
            do isc=1,nscal
               rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
            end do
         end do

         do i=0,jx+1
            !     corner 5/3
            u(i,0,k)=u(i,1,k)
            v(i,0,k)=v(i,1,k)
            w(i,0,k)=w(i,1,k)
            do isc=1,nscal
               rhov(isc,i,0,k)=rhov(isc,i,1,k)
            end do
            !     corner 5/4
            u(i,jy+1,k)=u(i,jy,k)
            v(i,jy+1,k)=v(i,jy,k)
            w(i,jy+1,k)=w(i,jy,k)
            do isc=1,nscal
               rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
            end do
         end do
      end if
   end if !side 5
      
   !-----------------------------------------------------------------------
   !     side 6
   if (infout6 .eq. 0) then
      if(myid.eq.nproc-1)then
         k=jz+1
         do j=0,jy+1
            !     corner 6/1
            u(0,j,k)=u(1,j,k)
            v(0,j,k)=v(1,j,k)
            w(0,j,k)=w(1,j,k)
            do isc=1,nscal
               rhov(isc,0,j,k)=rhov(isc,1,j,k)
            end do
            !     corner 6/2
            u(jx+1,j,k)=u(jx,j,k)
            v(jx+1,j,k)=v(jx,j,k)
            w(jx+1,j,k)=w(jx,j,k)
            do isc=1,nscal
               rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
            end do
         end do

         do i=0,jx+1
            if(myid.eq.0)then
               !     corner 6/3
               u(i,0,k)=u(i,1,k)
               v(i,0,k)=v(i,1,k)
               w(i,0,k)=w(i,1,k)
               do isc=1,nscal
                  rhov(isc,i,0,k)=rhov(isc,i,1,k)
               end do
            elseif(myid.eq.nproc-1)then
               !     corner 6/4
               u(i,jy+1,k)=u(i,jy,k)
               v(i,jy+1,k)=v(i,jy,k)
               w(i,jy+1,k)=w(i,jy,k)
               do isc=1,nscal
                  rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
               end do
            end if
         end do
      end if
   end if !side 6
   !-----------------------------------------------------------------------
   ! THIRD: case wall side (infout=2)
   !
   !     side 1
   if (infout1 .eq. 2) then
      i=0
      do k=kparastal,kparaendl
         !     corner 1/3
         u(i,0,k)=u(i,1,k)
         v(i,0,k)=v(i,1,k)
         w(i,0,k)=w(i,1,k)
         do isc=1,nscal
            rhov(isc,i,0,k)=rhov(isc,i,1,k)
         end do
         !     corner 1/4
         u(i,jy+1,k)=u(i,jy,k)
         v(i,jy+1,k)=v(i,jy,k)
         w(i,jy+1,k)=w(i,jy,k)
         do isc=1,nscal
            rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
         end do
      end do

      do j=0,jy+1
         if(myid.eq.0)then
            !       corner 1/5
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !       corner 1/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            enddo
         end if
      end do

   end if !side 1

   !-----------------------------------------------------------------------
   !     side 2
   if (infout2 .eq. 2) then
      i=jx+1
      do k=kparastal,kparaendl
         !       corner 2/3
         u(i,0,k)=u(i,1,k)
         v(i,0,k)=v(i,1,k)
         w(i,0,k)=w(i,1,k)
         do isc=1,nscal
            rhov(isc,i,0,k)=rhov(isc,i,1,k)
         end do

         !       corner 2/4
         u(i,jy+1,k)=u(i,jy,k)
         v(i,jy+1,k)=v(i,jy,k)
         w(i,jy+1,k)=w(i,jy,k)
         do isc=1,nscal
            rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
         end do
      end do

      do j=0,jy+1
         if(myid.eq.0)then
            !       corner 2/5
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !       corner 2/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            end do
         end if
      end do

   end if !side 2
      
   !-----------------------------------------------------------------------
   !     side 5
   if (infout5 .eq. 2) then
      if(myid.eq.0)then
         k=0
         do j=0,jy+1
            !       corner 5/1
            u(0,j,k)=u(1,j,k)
            v(0,j,k)=v(1,j,k)
            w(0,j,k)=w(1,j,k)
            do isc=1,nscal
               rhov(isc,0,j,k)=rhov(isc,1,j,k)
            end do

            !       corner 5/2
            u(jx+1,j,k)=u(jx,j,k)
            v(jx+1,j,k)=v(jx,j,k)
            w(jx+1,j,k)=w(jx,j,k)
            do isc=1,nscal
               rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
            end do
         end do

         do i=0,jx+1
            !       corner 5/3
            u(i,0,k)=u(i,1,k)
            v(i,0,k)=v(i,1,k)
            w(i,0,k)=w(i,1,k)
            do isc=1,nscal
               rhov(isc,i,0,k)=rhov(isc,i,1,k)
            end do

            !       corner 5/4
            u(i,jy+1,k)=u(i,jy,k)
            v(i,jy+1,k)=v(i,jy,k)
            w(i,jy+1,k)=w(i,jy,k)
            do isc=1,nscal
               rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
            end do
         end do
      end if
   end if !side 5

   !-----------------------------------------------------------------------
   !     side 6
   if (infout6 .eq. 2) then
      if(myid.eq.nproc-1)then
         k=jz+1
         do j=0,jy+1
            !       corner 6/1
            u(0,j,k)=u(1,j,k)
            v(0,j,k)=v(1,j,k)
            w(0,j,k)=w(1,j,k)
            do isc=1,nscal
               rhov(isc,0,j,k)=rhov(isc,1,j,k)
            end do

            !       corner 6/2
            u(jx+1,j,k)=u(jx,j,k)
            v(jx+1,j,k)=v(jx,j,k)
            w(jx+1,j,k)=w(jx,j,k)
            do isc=1,nscal
               rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
            end do
         end do
         do i=0,jx+1
            !       corner 6/3
            u(i,0,k)=u(i,1,k)
            v(i,0,k)=v(i,1,k)
            w(i,0,k)=w(i,1,k)
            do isc=1,nscal
               rhov(isc,i,0,k)=rhov(isc,i,1,k)
            end do
            !       corner 6/4
            u(i,jy+1,k)=u(i,jy,k)
            v(i,jy+1,k)=v(i,jy,k)
            w(i,jy+1,k)=w(i,jy,k)
            do isc=1,nscal
               rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
            end do
         end do
      end if
   end if !side 6


   !-----------------------------------------------------------------------
   !     side 3
   if (infout3 .eq. 2) then
      j=0
      do k=kparastal,kparaendl
         !       corner 3/1
         u(0,j,k)=u(1,j,k)
         v(0,j,k)=v(1,j,k)
         w(0,j,k)=w(1,j,k)
         do isc=1,nscal
            rhov(isc,0,j,k)=rhov(isc,1,j,k)
         enddo

         !       corner 3/2
         u(jx+1,j,k)=u(jx,j,k)
         v(jx+1,j,k)=v(jx,j,k)
         w(jx+1,j,k)=w(jx,j,k)
         do isc=1,nscal
            rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
         end do
      end do

      do i=0,jx+1
         if(myid.eq.0)then
            !       corner 3/5
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !       corner 3/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            end do
         end if
      end do

   end if !side 3

   !-----------------------------------------------------------------------
   !     side 4
   if (infout4 .eq. 2) then
      j=jy+1
      do k=kparastal,kparaendl
         !       corner 4/1
         u(0,j,k)=u(1,j,k)
         v(0,j,k)=v(1,j,k)
         w(0,j,k)=w(1,j,k)
         do isc=1,nscal
            rhov(isc,0,j,k)=rhov(isc,1,j,k)
         end do

         !       corner 4/2
         u(jx+1,j,k)=u(jx,j,k)
         v(jx+1,j,k)=v(jx,j,k)
         w(jx+1,j,k)=w(jx,j,k)
         do isc=1,nscal
            rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
         end do
      end do

      do i=0,jx+1
         if(myid.eq.0)then
            !       corner 4/5
            u(i,j,0)=u(i,j,1)
            v(i,j,0)=v(i,j,1)
            w(i,j,0)=w(i,j,1)
            do isc=1,nscal
               rhov(isc,i,j,0)=rhov(isc,i,j,1)
            end do
         elseif(myid.eq.nproc-1)then
            !       corner 4/6
            u(i,j,jz+1)=u(i,j,jz)
            v(i,j,jz+1)=v(i,j,jz)
            w(i,j,jz+1)=w(i,j,jz)
            do isc=1,nscal
               rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
            end do
         end if
      end do

   end if !side 4


   return
end
