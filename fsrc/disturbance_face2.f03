!***********************************************************************
subroutine disturbance_face2(ti,ktime,cx,cy,cz,nti)
   !***********************************************************************
   !     generate disturbance for face 2 for nesting procedure
   use mysending
   use myarrays_buffer_bodyforce
   !      use myarrays_nesting
   use myarrays_velo3
   !
   use scala3
   use orl

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer ktime,ntime
   integer i,j,k,n
   integer nti
   !      integer ispon
   integer istampo
   integer cont

      
   integer nreal,nstep
      
   real corr_factor_sp
   real cx,cy,cz
      
   real uprime,vprime,wprime
   
   real coef_x(n1+1-ispon:n1,n2,kparasta:kparaend)
   real coef_y(n1+1-ispon:n1,n2,kparasta:kparaend)
   real coef_z(n1+1-ispon:n1,n2,kparasta:kparaend)
       
   real ti,dt1,dt2,dt_tot
   real tke,tke_new

   integer cont_mean
   real meanx,meany,meanz
   real dev_x,dev_y,dev_z
   real var_x,var_y,var_z
   real velrms
   real segno_u,segno_v,segno_w
   real v_scale
      
   real a1,a2,a3,a_den
   real inv_dt
      
   real ke,ke_ref
   !-----------------------------------------------------------------------
   inv_dt = 1./dt
   !-----------------------------------------------------------------------
   ! allocation and initialization
   if(ktime .eq. 1)then
      
      corr_factor_sp = 0.
      
      !       face 1
      !       in x in space and time
      nreal = jy*jz/nproc
      nstep = nti
      allocate(noise_f2_xtime(nreal,-1:nstep*2))
      allocate(old_noise_f2_xtime(nreal))
      noise_f2_xtime      = 0.
      old_noise_f2_xtime  = 0.

      nreal = jy*jz/nproc
      nstep = ispon
      allocate(noise_f2_xspace(nreal,-1:nstep*2))
      noise_f2_xspace=0.
      istampo = 0
      call genero_random(corr_factor_sp,nreal,nstep,noise_f2_xspace,istampo,myid)
     
      !       in y in space and time
      nreal = ispon*jz/nproc
      nstep = nti
      allocate(noise_f2_ytime(nreal,-1:nstep*2))
      allocate(old_noise_f2_ytime(nreal))
      noise_f2_ytime     = 0.
      old_noise_f2_ytime = 0.
	 
      nreal = ispon*jz/nproc
      nstep = jy
      allocate(noise_f2_yspace(nreal,-1:nstep*2))
      noise_f2_yspace=0.
      istampo = 0
      call genero_random(corr_factor_sp,nreal,nstep,noise_f2_yspace,istampo,myid)

      !       in z in space and time
      nreal = jy*ispon
      nstep = nti
      allocate(noise_f2_ztime(nreal,-1:nstep*2))
      allocate(old_noise_f2_ztime(nreal))
      noise_f2_ztime     = 0.
      old_noise_f2_ztime = 0.
	
      nreal = jy*ispon
      nstep = jz/nproc
      allocate(noise_f2_zspace(nreal,-1:nstep*2))
      noise_f2_zspace=0.
      istampo = 0
      call genero_random(corr_factor_sp,nreal,nstep,noise_f2_zspace,istampo,myid)
   end if
   !-----------------------------------------------------------------------
   !     to have a colored noise connected to the previous number generation
   !
   if(ktime .gt.1 .and. ktime .eq. nti*(ktime/nti))then
      ntime = nti
      if(myid.eq.0)write(*,*)'TIME',ktime,ntime
      cont = 0
      do k=kparasta,kparaend
         do j=1,jy
            cont = cont + 1
            old_noise_f2_xtime(cont) = noise_f2_xtime(cont,ntime)
         end do
      end do

      cont = 0
      do k=kparasta,kparaend
         do i=1,ispon
            cont = cont + 1
            old_noise_f2_ytime(cont) = noise_f2_ytime(cont,ntime)
         end do
      end do

      cont = 0
      do j=1,jy
         do i=1,ispon
            cont = cont + 1
            old_noise_f2_ztime(cont) = noise_f2_ztime(cont,ntime)
         end do
      end do
   end if


   !-----------------------------------------------------------------------
   !     DISTURBANCE GENERATION IN TIME
   if(ktime .eq. 1 .or. ktime.eq.nti*(ktime/nti))then

      !      corr_factor = 100.
      
      !     in x in time
      nreal = jy*jz/nproc
      nstep = nti
      noise_f2_xtime=0.
      call genero_random(corr_factor,nreal,nstep,noise_f2_xtime,istampo,myid)
      
      !     in y in time
      nreal = ispon*jz/nproc
      nstep = nti
      noise_f2_ytime=0.
      istampo = 0
      call genero_random(corr_factor,nreal,nstep,noise_f2_ytime,istampo,myid)
   
      !     in z in time
      nreal = jy*ispon
      nstep = nti
      noise_f2_ztime=0.
      istampo = 0
      call genero_random(corr_factor,nreal,nstep,noise_f2_ztime,istampo,myid)

   end if  !if ktime

   !-----------------------------------------------------------------------
   !     shift for disturbance
   if(ktime .gt. 1 .and. ktime .eq. nti*(ktime/nti))then
      do ntime = 1,nti*2
         cont = 0
         do k=kparasta,kparaend
            do j=1,jy
               cont = cont + 1
               noise_f2_xtime(cont,ntime)=noise_f2_xtime(cont,ntime)+(old_noise_f2_xtime(cont)-noise_f2_xtime(cont,1))
            end do
         end do

         cont = 0
         do k=kparasta,kparaend
            do i=1,ispon
               cont = cont + 1
               noise_f2_ytime(cont,ntime)=noise_f2_ytime(cont,ntime)+(old_noise_f2_ytime(cont)-noise_f2_ytime(cont,1))
            end do
         end do

         cont = 0
         do j=1,jy
            do i=1,ispon
               cont = cont + 1
               noise_f2_ztime(cont,ntime)=noise_f2_ztime(cont,ntime)+(old_noise_f2_ztime(cont)-noise_f2_ztime(cont,1))
            end do
         end do
      end do
   end if
   !-----------------------------------------------------------------------
   !     noise in space and time in the sponge region

   ntime = ktime -(ktime/nti)*nti
   call creo_xyz_f2(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)
     
     
   meanx = 0.
   meany = 0.
   meanz = 0.
   cont_mean = 0
   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx
            cont_mean = cont_mean + 1
            meanx = meanx + coef_x(i,j,k)
            meany = meany + coef_y(i,j,k)
            meanz = meanz + coef_z(i,j,k)
         end do
      end do
   end do
   meanx = meanx/real(cont_mean)
   meany = meany/real(cont_mean)
   meanz = meanz/real(cont_mean)
   !      if(myid ==0)write(*,*)'MEAN:',meanx,meany,meanz
   dev_x = 0.
   dev_y = 0.
   dev_z = 0.
   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx
            coef_x(i,j,k) = coef_x(i,j,k)-meanx
            coef_y(i,j,k) = coef_y(i,j,k)-meany
            coef_z(i,j,k) = coef_z(i,j,k)-meanz
	 
            dev_x = dev_x + coef_x(i,j,k)*coef_x(i,j,k)
            dev_y = dev_y + coef_y(i,j,k)*coef_y(i,j,k)
            dev_z = dev_z + coef_z(i,j,k)*coef_z(i,j,k)
         end do
      end do
   end do
  
   dev_x = sqrt(dev_x/real(cont_mean))
   dev_y = sqrt(dev_y/real(cont_mean))
   dev_z = sqrt(dev_z/real(cont_mean))
   var_x = 0.05
   var_y = 0.05
   var_z = 0.05
   !      if(myid ==0)write(*,*)'DEV:',dev_x,dev_y,dev_z
   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx
            coef_x(i,j,k) = coef_x(i,j,k)*var_x/dev_x
            coef_y(i,j,k) = coef_y(i,j,k)*var_y/dev_y
            coef_z(i,j,k) = coef_z(i,j,k)*var_z/dev_z
         end do
      end do
   end do
   !-----------------------------------------------------------------------
   !     for interpolation in time of tke
   !      dt1 = ti - ti_pom_old
   !      dt2 = ti_pom_new - ti
   !      dt_tot = 1./(dt1+dt2)

   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx  !1,ispon
            if(index_out2(j,k) .gt. 0.01)then

               !fluctuation
               uprime = u(i,j,k)-up2(j,k)
               vprime = v(i,j,k)-vp2(j,k)
               wprime = w(i,j,k)-wp2(j,k)

               ! coef to scale tke for each direction a1+a2+a3 = 1
               a_den=1./( up2(j,k)*up2(j,k)+vp2(j,k)*vp2(j,k)+wp2(j,k)*wp2(j,k) )
     	 
               a1 = up2(j,k)*up2(j,k) * a_den
               a2 = vp2(j,k)*vp2(j,k) * a_den
               a3 = wp2(j,k)*wp2(j,k) * a_den
	 
                        !kinetic energy
               ke_ref =   up2(j,k)*up2(j,k)+vp2(j,k)*vp2(j,k)+wp2(j,k)*wp2(j,k)

               ke =   u(i,j,k)*u(i,j,k) +v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
     
               if(ke .gt. 1.3*ke_ref)cycle
	 
	 	 
               ! left hand side ( TKE(n+1)-TKE(n) ) / Delta T
               tke = tke2(j,k)

               tke_new = .5*(  uprime*uprime + vprime*vprime + wprime*wprime )
	 	 
               tke = abs((tke_new - tke)*inv_dt)
     
               ! fluctuation velocity scale
               velrms = sqrt(2.*tke_new)
     
               !         v_scale = up2(j,k)*up2(j,k) + vp2(j,k)*vp2(j,k) + wp2(j,k)*wp2(j,k)
     
               !         v_scale = 0.01*sqrt(v_scale)
     
               !         velrms = v_scale
     
               segno_u = sign(1.,up2(j,k))
               segno_v = sign(1.,vp2(j,k))
               segno_w = sign(1.,wp2(j,k))
	 
               segno_u = sign(1.,uprime)
               segno_v = sign(1.,vprime)
               segno_w = sign(1.,wprime)
	 	 

               if(abs(uprime) .lt. 0.1*velrms )then
                  uprime = 0.1*segno_u*velrms
               end if
               if(abs(vprime) .lt. 0.1*velrms )then
                  vprime = 0.1*segno_v*velrms
               end if
               if(abs(wprime) .lt. 0.1*velrms )then
                  wprime = 0.1*segno_w*velrms
               end if


               ! bodyforce
               bcsi_f2(i,j,k)=cx*(a1*tke/abs(uprime))*coef_x(i,j,k)

               beta_f2(i,j,k)=cy*(a2*tke/abs(vprime))*coef_y(i,j,k)

               bzet_f2(i,j,k)=cz*(a3*tke/abs(wprime))*coef_z(i,j,k)
	 

            end if
         end do
      end do
   end do
9200 format(21e18.10)

   return
end
      
!***********************************************************************
subroutine creo_xyz_f2(n1,n2,n3,kparasta,kparaend,coef_x,coef_y,coef_z,ntime,myid)
   !***********************************************************************
   use myarrays_buffer_bodyforce
      
   implicit none
   !-----------------------------------------------------------------------
   integer i,j,k,cont,icont,ntime
   integer n1,n2,n3,kparasta,kparaend
   integer myid
      
   real coef_x(n1+1-ispon:n1,n2,kparasta:kparaend)
   real coef_y(n1+1-ispon:n1,n2,kparasta:kparaend)
   real coef_z(n1+1-ispon:n1,n2,kparasta:kparaend)
   !-----------------------------------------------------------------------
   !     in x
   icont = 0
   do i=n1+1-ispon,n1 !1,ispon
      cont = 0
      icont = icont + 1
      do k=kparasta,kparaend
         do j=1,n2
            cont = cont + 1
            coef_x(i,j,k) = noise_f2_xtime(cont,ntime)+noise_f2_xspace(cont,icont)
         end do
      end do
   end do
                 
   !     in y
   do j=1,n2
      cont = 0
      do i=n1+1-ispon,n1 !1,ispon
         do k=kparasta,kparaend
            cont = cont + 1
            coef_y(i,j,k) = noise_f2_ytime(cont,ntime)+noise_f2_yspace(cont,j)
         end do
      end do
   end do
      
   !     in z
   do k=kparasta,kparaend
      cont = 0
      do i=n1+1-ispon,n1 !1,ispon
         do j=1,n2
            cont = cont + 1
            coef_z(i,j,k) = noise_f2_ztime(cont,ntime)+noise_f2_zspace(cont,k-kparasta+1)
         end do
      end do
   end do

   return
end
