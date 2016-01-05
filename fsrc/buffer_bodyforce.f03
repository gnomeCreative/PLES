!***********************************************************************
subroutine buffer_bodyforce(ti,ktime,bcsi,beta,bzet)
   !***********************************************************************
   ! generates a coloured noise to start turbulence from tke information
   ! tke comes from a RANS model, this sub is used for nesting procedure

   use mysending
   use mysettings_boundary
   use myarrays_buffer_bodyforce
   use myarrays_nesting
   use myarrays_metri3
   use myarrays_velo3
   !
   use scala3
   use tipologia
   use orl
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer ktime
   integer nti,ntime
      
   integer i,j,k,n
   integer cont
   !      integer ispon,kspon
   integer istampo

   real cx,cy,cz
   real L_d
   real bcsi(n1,n2,kparasta:kparaend)
   real beta(n1,n2,kparasta:kparaend)
   real bzet(n1,n2,kparasta:kparaend)


   real apply_dist1(n1,kparasta:kparaend)
   real apply_dist2(n1,kparasta:kparaend)
   real apply_dist5(n1,n2)
   real apply_dist6(n1,n2)
   real mass_sign
      
   real, allocatable :: var_prov(:,:,:),var_prov2(:,:,:)
   real, allocatable :: local_dump(:,:,:)
      
   real v1,v2,ti
   !
   !-----------------------------------------------------------------------
   ! some parameter
   !      n_ti_pom = 1
   !
   !-----------------------------------------------------------------------
   !for tkepom
   if(ktime .eq.1)then
      ntpom = 1
   end if
   if(myid==0)write(*,*)'BUFFER TURBULENCE on'
   if(ibodybuffer1==1)write(*,*)'side 1 on'
   if(ibodybuffer2==1)write(*,*)'side 2 on'
   if(ibodybuffer5==1)write(*,*)'side 5 on'
   if(ibodybuffer6==1)write(*,*)'side 6 on'
      
      
      
   if(ti .gt. ti_pom_new)ntpom = ntpom + 1
   !-----------------------------------------------------------------------

   !     scaling factor for forcing term
   cx = 0.1
   cy = 0.1
   cz = 0.1

   !     time of colored disturbance
   nti = 100
      
   !     check sponge region dimension in i and k directions
   if(kspon .gt. jz/nproc .and. myid.eq.0)then
      write(*,*)'NESTING PROBLEM'
      write(*,*)'program not suited to have a sponge region'
      write(*,*)'larger than the processor z allocation'
      stop
   end if
   !
   !-----------------------------------------------------------------------
   !     allocation and initialization
   if(ktime .eq. 1)then
      !        allocate(tkepom1(jy,kparasta:kparaend,n_ti_pom)) !face 1
      !        allocate(tkepom2(jy,kparasta:kparaend,n_ti_pom)) !face 2
      !        allocate(tkepom5(jx,jy               ,n_ti_pom)) !face 3
      !        allocate(tkepom6(jx,jy               ,n_ti_pom)) !face 4
      !	tkepom1 = 0.
      !	tkepom2 = 0.
      !	tkepom5 = 0.
      !	tkepom6 = 0.

      allocate(bcsi_f1(ispon,jy,kparasta:kparaend))
      allocate(beta_f1(ispon,jy,kparasta:kparaend))
      allocate(bzet_f1(ispon,jy,kparasta:kparaend))
      bcsi_f1 = 0.
      beta_f1 = 0.
      bzet_f1 = 0.
			
      allocate(bcsi_f2(jx-ispon+1:jx,jy,kparasta:kparaend))
      allocate(beta_f2(jx-ispon+1:jx,jy,kparasta:kparaend))
      allocate(bzet_f2(jx-ispon+1:jx,jy,kparasta:kparaend))
      bcsi_f2 = 0.
      beta_f2 = 0.
      bzet_f2 = 0.

      allocate(bcsi_f5(jx,jy,kspon))
      allocate(beta_f5(jx,jy,kspon))
      allocate(bzet_f5(jx,jy,kspon))
      bcsi_f5 = 0.
      beta_f5 = 0.
      bzet_f5 = 0.
	
      allocate(bcsi_f6(jx,jy,jz-kspon+1:jz))
      allocate(beta_f6(jx,jy,jz-kspon+1:jz))
      allocate(bzet_f6(jx,jy,jz-kspon+1:jz))
      bcsi_f6 = 0.
      beta_f6 = 0.
      bzet_f6 = 0.
   end if
   !
   !-----------------------------------------------------------------------
   !     DSITURBANCE
   !     on face 1
   !      if(infout1 /= 0 .or. ibodybuffer1==1)then
   if( ibodybuffer1==1)then
      call disturbance_face1(ti,ktime,cx,cy,cz,nti)
   end if
      
   !     on face 2
   if(ibodybuffer2==1)then
      call disturbance_face2(ti,ktime,cx,cy,cz,nti)
   end if
      
   !     on face 5
   if(ibodybuffer5==1)then
      if(myid==0)then
         call disturbance_face5(ti,ktime,cx,cy,cz,nti)
      end if
   end if
      
   !     on face 6
   if(ibodybuffer6==1)then
      if(myid==nproc-1)then
         call disturbance_face6(ti,ktime,cx,cy,cz,nti)
      end if
   end if
   !
   !-----------------------------------------------------------------------
   !     SPOONGE DUMPING

   !     face 1
   istampo = 1
   allocate(local_dump(ispon,n2,kparasta:kparaend))
   L_d = (x(ispon,1,kparasta)-x(0,1,kparasta))
   !      if(myid.eq.0)write(*,*)'L_d face 1 ',L_d
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,ispon
            local_dump(i,j,k) = x(i,j,k) - x(0,j,k)
            local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
         end do
      end do
   end do
   !quando faccia la f2 dovrebbe bastare invertire la local_dump?
   call sponge_dumping(ispon,n2,n3/nproc,local_dump,bcsi_f1,beta_f1,bzet_f1,istampo)
   deallocate(local_dump)



   !     face 2
   istampo = 1
   allocate(local_dump(n1+1-ispon:n1,n2,kparasta:kparaend))
   L_d = (x(jx,1,kparasta)-x(jx-ispon,1,kparasta))
   !      if(myid.eq.0)write(*,*)'L_d face 2 ',L_d
   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx  !1,ispon
            local_dump(i,j,k) = x(jx,j,k) - x(i-1,j,k)
            local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
         end do
      end do
   end do
   !quando faccia la f2 dovrebbe bastare invertire la local_dump?
   call sponge_dumping(ispon,n2,n3/nproc,local_dump,bcsi_f2,beta_f2,bzet_f2,istampo)
   deallocate(local_dump)

   !     face 5
   if(myid.eq.0)then
      istampo = 1
      allocate(local_dump(n1,n2,kspon))
      L_d = (z(1,1,kspon)-z(1,1,0))
      !      if(myid.eq.0)write(*,*)'L_d face 5 ',L_d
      do k=1,kspon !kparasta,kparaend
         do j=1,jy
            do i=1,jx  !ispon
               local_dump(i,j,k) = z(i,j,k) - z(i,j,0)
               local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
            !	 if(i.eq.1 .and. j.eq.10)then
            !	    write(452,*)k,local_dump(i,j,k)
            !	 end if
            end do
         end do
      end do

      call sponge_dumping(n1,n2,kspon,local_dump,bcsi_f5,beta_f5,bzet_f5,istampo)
      deallocate(local_dump)
   end if



   !     face 6
   if(myid.eq.nproc-1)then
      istampo = 1
      allocate(local_dump(n1,n2,n3+1-kspon:n3))
      L_d = (z(1,1,jz)-z(1,1,jz-kspon))
      !      if(myid.eq.nproc-1)write(*,*)'L_d face 6 ',L_d
      do k=jz+1-kspon,jz !1,kspon !kparasta,kparaend
         do j=1,jy
            do i=1,jx  !ispon
               local_dump(i,j,k) = z(i,j,jz) - z(i,j,k-1)
               local_dump(i,j,k) = 1. - local_dump(i,j,k) / L_d
            !	 if(i.eq.1 .and. j.eq.10)then
            !	    write(453,*)k,local_dump(i,j,k)
            !	 end if
	 
            end do
         end do
      end do

      call sponge_dumping(n1,n2,kspon,local_dump,bcsi_f6,beta_f6,bzet_f6,istampo)
      deallocate(local_dump)
   end if


   !----------------------------------------------------------------
   allocate( var_prov(jx+1-ispon:jx,jy,kparasta:kparaend))
   allocate(var_prov2(jx+1-ispon:jx,jy,kparasta:kparaend))
   allocate(local_dump(n1+1-ispon:n1,n2,kparasta:kparaend))
      
   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx  !1,ispon
            var_prov(i,j,k) = annit(i,j,k)
            var_prov2(i,j,k) = annitV(i,j,k)
         end do
      end do
   end do

   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx  !1,ispon
            !         local_dump(i,j,k) = x(jx,j,k) - x(i-1,j,k)
            local_dump(i,j,k) = (x(jx,j,k) - x(i-1,j,k) ) -  L_d
            local_dump(i,j,k) = 1. - local_dump(i,j,k) * L_d
         end do
      end do
   end do
      
   call sponge_dumping_visco(ispon,n2,n3/nproc,local_dump,var_prov,var_prov,var_prov2,istampo)

   do k=kparasta,kparaend
      do j=1,jy
         do i=jx+1-ispon,jx  !1,ispon
         !          annit(i,j,k) =  var_prov(i,j,k)
         !         annitV(i,j,k) = var_prov2(i,j,k)
         end do
      end do
   end do
   deallocate(var_prov)
   deallocate(var_prov2)
      
   deallocate(local_dump)

     
   !-----------------------------------------------------------------------
   !     CORNER CORRECTION BETWEEN SPONGE REGION IN i AND k
   !-----------------------------------------------------------------------
   !     CHECK IF DISTURBANCE MUST BE APPLIED
   apply_dist1 = 0.
   apply_dist2 = 0.
   do k=kparasta,kparaend
      do j=1,jy
                  !side 1
         ! uc>0 inflow; uc<0 outflow
         !inflow  mass_sign = 1
         !outflow mass_sign =-1
         mass_sign = INT(sign(1.,uc(0,j,k)))

         !inflow  apply_dist = 1.
         !outflow apply_dist = 0.
         apply_dist1(j,k) = 0.5*(1+mass_sign)
         !--------------------------------------

                  !side 2
         ! uc<0 inflow; uc>0 outflow
         !inflow  mass_sign = 1
         !outflow mass_sign =-1
         mass_sign = - INT(sign(1.,uc(jx,j,k)))

         !inflow  apply_dist = 1.
         !outflow apply_dist = 0.
         apply_dist2(j,k) = 0.5*(1+mass_sign)
      end do
   end do

   if(myid==0)then
      apply_dist5 = 0.
      do j=1,jy
         do i=1,jx
                     !side 5
            ! wc>0 inflow; wc<0 outflow
            !inflow  mass_sign = 1
            !outflow mass_sign =-1
            mass_sign = INT(sign(1.,wc(i,j,0)))
	  
            !inflow  apply_dist = 1.
            !outflow apply_dist = 0.
            apply_dist5(i,j) = 0.5*(1+mass_sign)
         end do
      end do
   end if

   if(myid==nproc-1)then
      apply_dist6 = 0.
      do j=1,jy
         do i=1,jx
                     !side 6
            ! wc<0 inflow; wc>0 outflow
            !inflow  mass_sign = 1
            !outflow mass_sign =-1
            mass_sign = - INT(sign(1.,wc(i,j,jz)))
	  
            !inflow  apply_dist = 1.
            !outflow apply_dist = 0.
            apply_dist6(i,j) = 0.5*(1+mass_sign)
         end do
      end do
   end if
   !
   !-----------------------------------------------------------------------
   !     GIVE THE VALUE TO THE BODYFORCE

   !     sponge in front of face 1
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,ispon
      
            bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f1(i,j,k)*apply_dist1(j,k)
            beta(i,j,k) = beta(i,j,k) + beta_f1(i,j,k)*apply_dist1(j,k)
            bzet(i,j,k) = bzet(i,j,k) + bzet_f1(i,j,k)*apply_dist1(j,k)

         end do
      end do
   end do
      
   !     sponge in front of face 2
   do k=kparasta,kparaend
      do j=1,jy
         do i=jx-ispon+1,jx
      
            bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f2(i,j,k)*apply_dist2(j,k)
            beta(i,j,k) = beta(i,j,k) + beta_f2(i,j,k)*apply_dist2(j,k)
            bzet(i,j,k) = bzet(i,j,k) + bzet_f2(i,j,k)*apply_dist2(j,k)

         end do
      end do
   end do

   !     sponge in front of face 5
   if(myid.eq.0)then
      do k=1,kspon
         do j=1,jy
            do i=1,jx
      
               bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f5(i,j,k)*apply_dist5(i,j)
               beta(i,j,k) = beta(i,j,k) + beta_f5(i,j,k)*apply_dist5(i,j)
               bzet(i,j,k) = bzet(i,j,k) + bzet_f5(i,j,k)*apply_dist5(i,j)

            end do
         end do
      end do
   end if
      
   !     sponge in front of face 6
   if(myid.eq.nproc-1)then
      do k=jz-kspon+1,jz
         do j=1,jy
            do i=1,jx

               bcsi(i,j,k) = bcsi(i,j,k) + bcsi_f6(i,j,k)*apply_dist6(i,j)
               beta(i,j,k) = beta(i,j,k) + beta_f6(i,j,k)*apply_dist6(i,j)
               bzet(i,j,k) = bzet(i,j,k) + bzet_f6(i,j,k)*apply_dist6(i,j)

            end do
         end do
      end do
   end if
   !-----------------------------------------------------------------------
   return
end


!***********************************************************************
subroutine sponge_dumping(nx,ny,nz,local_dump,var_csi,var_eta,var_zet,istampo)
   !***********************************************************************
   !      I use the form:
   !            (u^(n+1)-u^(n))/dt = - sigma(u^(n)-u_desired)
   !      this means that sigma at i=end_domain must be of order 1/dt
   !      I use an exponential form for sigma:
   !            sigma = e^beta*g(x)
   !      with g(x)= L - (Xend - X) with L the length of the sponge region
   !      so at the domain end g(x) = L
   !      beta is chosen so that sigma has order 1/dt

   !       L_sponge = x(ispon)-x(0)

   !       angbeta = (1/L_sponge)*log(1/dt)
   integer i,j,k
   integer nx,ny,nz
   integer istampo

   real local_dump(nx,ny,nz)

   real var_csi(nx,ny,nz)
   real var_eta(nx,ny,nz)
   real var_zet(nx,ny,nz)

   real var_csi_damp(nx,ny,nz)
   real var_eta_damp(nx,ny,nz)
   real var_zet_damp(nx,ny,nz)

   real esponente,sigma,Acoef

   Acoef = 1.
       
   do k=1,nz
      do j=1,ny
         do i=1,nx
            esponente = -.5*( 3.5 * local_dump(i,j,k) )**2.
            sigma = Acoef * 2.7183**esponente
            !	if(j.eq.1 .and. k.eq.1 .and. istampo.eq.1)then
            !	  write(600,*)i,esponente,sigma
            !	end if
            var_csi_damp(i,j,k) = var_csi(i,j,k)-sigma*(var_csi(i,j,k)-0.)
            var_eta_damp(i,j,k) = var_eta(i,j,k)-sigma*(var_eta(i,j,k)-0.)
            var_zet_damp(i,j,k) = var_zet(i,j,k)-sigma*(var_zet(i,j,k)-0.)

         end do
      end do
   end do
       
   do k=1,nz
      do j=1,ny
         do i=1,nx
            var_csi(i,j,k) = var_csi_damp(i,j,k)
            var_eta(i,j,k) = var_eta_damp(i,j,k)
            var_zet(i,j,k) = var_zet_damp(i,j,k)
         end do
      end do
   end do
      
   return
end

!***********************************************************************
subroutine sponge_dumping_visco(nx,ny,nz,local_dump,var_csi,var_eta,var_zet,istampo)
   !***********************************************************************
   !      I use the form:
   !            (u^(n+1)-u^(n))/dt = - sigma(u^(n)-u_desired)
   !      this means that sigma at i=end_domain must be of order 1/dt
   !      I use an exponential form for sigma:
   !            sigma = e^beta*g(x)
   !      with g(x)= L - (Xend - X) with L the length of the sponge region
   !      so at the domain end g(x) = L
   !      beta is chosen so that sigma has order 1/dt

   !       L_sponge = x(ispon)-x(0)

   !       angbeta = (1/L_sponge)*log(1/dt)
   integer i,j,k
   integer nx,ny,nz
   integer istampo

   real local_dump(nx,ny,nz)

   real var_csi(nx,ny,nz)
   real var_eta(nx,ny,nz)
   real var_zet(nx,ny,nz)

   real var_csi_damp(nx,ny,nz)
   real var_eta_damp(nx,ny,nz)
   real var_zet_damp(nx,ny,nz)

   real esponente,sigma,Acoef

   Acoef = 1.
       
   do k=1,nz
      do j=1,ny
         do i=1,nx
            esponente = -.5*( 3.5 * local_dump(i,j,k) )**2.
            sigma = Acoef * 2.7183**esponente
            !	if(j.eq.1 .and. k.eq.1 .and. istampo.eq.1)then
            !	  write(600,*)i,esponente,sigma
            !	end if



            var_csi_damp(i,j,k) = var_csi(i,j,k)-sigma*(var_csi(i,j,k)-4.*var_csi(i,j,k))
            var_eta_damp(i,j,k) = var_eta(i,j,k)-sigma*(var_eta(i,j,k)-4.*var_eta(i,j,k))
            var_zet_damp(i,j,k) = var_zet(i,j,k)-sigma*(var_zet(i,j,k)-4.*var_zet(i,j,k))



         end do
      end do
   end do
       
   do k=1,nz
      do j=1,ny
         do i=1,nx
            var_csi(i,j,k) = var_csi_damp(i,j,k)
            var_eta(i,j,k) = var_eta_damp(i,j,k)
            var_zet(i,j,k) = var_zet_damp(i,j,k)
         end do
      end do
   end do
      
   return
end


