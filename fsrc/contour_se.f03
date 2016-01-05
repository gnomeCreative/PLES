!************************************************************************
subroutine contour_se
   !***********************************************************************
   ! boundary condition on cartesian velocities and scalars and
   ! for controvariant fluxes
   !
   use mysettings_boundary
   use myarrays_metri3
   use myarrays_velo3
   use mysending
   use scala3
   use period
   use orl

   implicit none
   !
   !-----------------------------------------------------------------------
   ! arrays declaration
   !
   integer i,j,k
   integer ii,jj,kk
   integer isc
   !
   !-----------------------------------------------------------------------
   !     face 1 'sinistra' (for periodicity see vel_up)
   do ii=1,ip
      if(infout1.eq.0)then  !inflow
         do k=kparasta,kparaend !1,jz
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
         if(myid.eq.0)write(*,*)'inflow face 1'
      elseif(infout1.eq.1)then  !outflow
         do k=kparasta,kparaend
            do j=1,jy
               u(0,j,k)=u(0,j,k)
               v(0,j,k)=v(0,j,k)
               w(0,j,k)=w(0,j,k)
               do isc=1,nscal
                  rhov(isc,0,j,k)=rhov(isc,1,j,k)
               end do
            enddo
         enddo
         if(myid.eq.0)write(*,*)'outflow face 1'
      elseif(infout1.eq.2)then  !wall
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
	
         end if
         if(myid.eq.0)write(*,*)'solid wall face 1'
      end if

      !     face 2 'destra' (for periodicity see vel_up)
      if(infout2.eq.0)then  !inflow
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
         if(myid.eq.0)write(*,*)'inflow face 2'
      elseif(infout2.eq.1)then  !outflow
         do k=kparasta,kparaend
            do j=1,jy
               u(jx+1,j,k)=u(jx+1,j,k)
               v(jx+1,j,k)=v(jx+1,j,k)
               w(jx+1,j,k)=w(jx+1,j,k)
               do isc=1,nscal
                  rhov(isc,jx+1,j,k)=rhov(isc,jx,j,k)
               end do
            enddo
         enddo
         if(myid.eq.0)write(*,*)'outflow face 2'
      elseif(infout2.eq.2)then  !wall
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
	
         end if
         if(myid.eq.0)write(*,*)'solid wall face 2'
      end if
   enddo   !end loop  ii=1,ip
   !
   !-----------------------------------------------------------------------
   !     face 3 'sotto' (for periodicity see vel_up)
   !
   do jj=1,jp
      if(infout3.eq.0)then  !inflow
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
         if(myid.eq.0)write(*,*)'inflow face 3'
      elseif(infout3.eq.1)then  !outflow
         do k=kparasta,kparaend
            do i=1,jx
               u(i,0,k)=u(i,0,k)
               v(i,0,k)=v(i,0,k)
               w(i,0,k)=w(i,0,k)
               do isc=1,nscal
                  rhov(isc,i,0,k)=rhov(isc,i,1,k)
               end do
            enddo
         enddo
         if(myid.eq.0)write(*,*)'outflow face 3'
      elseif(infout3.eq.2)then  !wall
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
	
         end if
         if(myid.eq.0)write(*,*)'solid wall face 3'
      end if

      !     face 4 'sopra' (for periodicity see vel_up)
      if(infout4.eq.0)then  !inflow
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
         if(myid.eq.0)write(*,*)'inflow face 4'
      elseif(infout4.eq.1)then  !outflow
         do k=kparasta,kparaend
            do i=1,jx
               u(i,jy+1,k)=u(i,jy+1,k)
               v(i,jy+1,k)=v(i,jy+1,k)
               w(i,jy+1,k)=w(i,jy+1,k)
               do isc=1,nscal
                  rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
               end do
            enddo
         enddo
         if(myid.eq.0)write(*,*)'outflow face 4'
      elseif(infout4.eq.2)then  !wall
         if(iboun4==0)then
            do k=kparasta,kparaend
               do i=1,jx
                  u(i,jy+1,k)=1.
                  v(i,jy+1,k)=0.
                  w(i,jy+1,k)=1.
                  do isc=1,nscal
                     rhov(isc,i,jy+1,k)=rhov(isc,i,jy,k)
                  end do
                  vc(i,jy,k)=etx(i,jy,k)*u(i,jy+1,k) +ety(i,jy,k)*v(i,jy+1,k) +etz(i,jy,k)*w(i,jy+1,k)
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
               enddo
            enddo
         elseif(iboun4==2)then

         end if
         if(myid.eq.0)write(*,*)'solid wall face 4'
      end if
   enddo !end loop  jj=1,jp
   !
   !-----------------------------------------------------------------------
   !     face 5 'indietro' (for periodicity see vel_up)
   !
   do kk=1,kp
      if (myid.eq.0) then
         if(infout5.eq.0)then  !inflow
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
            write(*,*)'inflow face 5'
         elseif(infout5.eq.1)then  !outflow
            do j=1,jy
               do i=1,jx
                  u(i,j,0)=u(i,j,0)
                  v(i,j,0)=v(i,j,0)
                  w(i,j,0)=w(i,j,0)
                  do isc=1,nscal
                     rhov(isc,i,j,0)=rhov(isc,i,j,1)
                  end do
               enddo
            enddo
            write(*,*)'outflow face 5'
         elseif(infout5.eq.2)then  !wall
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

            end if
            write(*,*)'solid wall face 5'
         endif
      endif ! myid=0
      !................................................
      !     face 6 'avanti' (for periodicity see vel_up)
      !
      if (myid.eq.nproc-1) then
         if(infout6.eq.0)then  !inflow
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
            write(*,*)'inflow face 6'
         elseif(infout6.eq.1)then  !outflow
            do j=1,jy
               do i=1,jx
                  u(i,j,jz+1)=u(i,j,jz+1)
                  v(i,j,jz+1)=v(i,j,jz+1)
                  w(i,j,jz+1)=w(i,j,jz+1)
                  do isc=1,nscal
                     rhov(isc,i,j,jz+1)=rhov(isc,i,j,jz)
                  end do
               enddo
            enddo
            write(*,*)'outflow face 6'
         elseif(infout6.eq.2)then  !wall
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
	
            end if
            write(*,*)'solid wall face 6'
         endif
      endif !myid=nproc-1
   enddo !end loop  kk=1,kp
   !-----------------------------------------------------------------------

   return
end
