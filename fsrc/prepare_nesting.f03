!***********************************************************************
subroutine prepare_nesting(ti,dbbx,unit1,unit2,unit3,unit4,unit5,unit6)
   !***********************************************************************
   use mysending
   use myarrays_nesting
   use myarrays_buffer_bodyforce
   use myarrays_velo3
   use myarrays_moisture
   !
   use scala3
   use orl

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   !      integer imoist
   integer unit1,unit2,unit3
   integer unit4,unit5,unit6
   integer i,j,k,n,isc,kpst,kpnd
      
   real ti,dbbx,ntireal
   real var1,var2,var3,var4,var(nscal)
            
   character*1 commento
   !-----------------------------------------------------------------------
   ntke = 1
   !-----------------------------------------------------------------------
   !     array allocation and initialization for nesting
   if(infout1 /=0)then
      !     side   1
      allocate(uo1(0:jy+1,0:jz+1))
      allocate(vo1(0:jy+1,0:jz+1))
      allocate(wo1(0:jy+1,0:jz+1))
      allocate(rhovo1(nscal,0:jy+1,0:jz+1))
      allocate(ucp1(jy,kparasta:kparaend))
              
      uo1=0.
      vo1=0.
      wo1=0.
      rhovo1=0.
      ucp1 = 0.
  

      read(unit1,*)commento
      read(unit1,*)ntireal,ti_pom_fin
      n_ti_pom = int(ntireal)
      read(unit1,*)ti_pom_old

      allocate(tkepom1(0:jy+1,kparasta-1:kparaend+1,n_ti_pom)) !face 1
      tkepom1 = 0.

      do k=0,jz+1
         do j=0,jy+1
            read(unit1,*)var1,var2,var3,var4,(var(n),n=1,nscal)
         
            if(k.ge.kparasta.and.k.le.kparaend)then
               uo1(j,k) = var1*index_out1(j,k)
               vo1(j,k) = var2*index_out1(j,k)
               wo1(j,k) = var3*index_out1(j,k)
               tkepom1(j,k,ntke)=var4*index_out1(j,k)
               do isc = 1,nscal
                  rhovo1(isc,j,k) = var(isc)*index_rho1(j,k)
               end do
            endif
         end do
      end do

   end if
       
   !-----------------------------------------------------------------------
   !     side   2
   if(infout2 /=0)then
      allocate(uo2(0:jy+1,0:jz+1))
      allocate(vo2(0:jy+1,0:jz+1))
      allocate(wo2(0:jy+1,0:jz+1))
      allocate(rhovo2(nscal,0:jy+1,0:jz+1))
      allocate(ucp2(jy,kparasta:kparaend))
              
      uo2=0.
      vo2=0.
      wo2=0.
      rhovo2=0.
      ucp2 = 0.
    

      read(unit2,*)commento
      read(unit2,*)ntireal,ti_pom_fin
      n_ti_pom = int(ntireal)
      read(unit2,*)ti_pom_old

      allocate(tkepom2(0:jy+1,kparasta-1:kparaend+1,n_ti_pom)) !face 2
      tkepom2 = 0.

      do k=0,jz+1
         do j=0,jy+1
            read(unit2,*)var1,var2,var3,var4,(var(n),n=1,nscal)
            if(k.ge.kparasta.and.k.le.kparaend)then
               uo2(j,k) = var1*index_out2(j,k)
               vo2(j,k) = var2*index_out2(j,k)
               wo2(j,k) = var3*index_out2(j,k)
               tkepom2(j,k,ntke)=var4*index_out2(j,k)
               do isc = 1,nscal
                  rhovo2(isc,j,k) = var(isc)*index_rho2(j,k)
               end do
            endif
         enddo
      enddo

   end if
   !-----------------------------------------------------------------------
   !     side   3
   if(imoist == 1)then
      if(infout3 /=0)then
         allocate(    tauw3nest(n_ti_pom))
         allocate(scalarflux3nest(nscal,n_ti_pom))
	 	 	 
         tauw3nest = 0.
         scalarflux3nest = 0.

         read(unit3,*)commento
         do i=1,n_ti_pom
            read(unit3,*) ti_pom_old
            read(unit3,*) tauw3nest(i),(scalarflux3nest(isc,i),isc=1,nscal)
         end do
         close(83)
       
      end if


      !      side 4
      if(infout4 /=0)then

      !      allocate(    Qvflux4nest(n_ti_pom))
      !      allocate(TKEflux4nest(n_ti_pom))
      !      allocate( vapflux4nest(n_ti_pom))
      !      allocate(scalarflux4nest(nscal,n_ti_pom))
	      
		      
      !      do i=1,n_ti_pom
      !	 read(84,*)
      !       end do
      !      close(84)

      end if
   end if
        
   !-----------------------------------------------------------------------
   !     side   5
   if(infout5 /=0)then
      allocate(uo5(0:jx+1,0:jy+1))
      allocate(vo5(0:jx+1,0:jy+1))
      allocate(wo5(0:jx+1,0:jy+1))
      allocate(rhovo5(nscal,0:jx+1,0:jy+1))
      allocate(wcp5(jx,jy))

      uo5=0.
      vo5=0.
      wo5=0.
      rhovo5=0.
      wcp5 = 0.


      read(unit5,*)commento
      read(unit5,*)ntireal,ti_pom_fin
      n_ti_pom = int(ntireal)
      read(unit5,*)ti_pom_old

      allocate(tkepom5(0:jx+1,0:jy+1,n_ti_pom))
      tkepom5 = 0.

      do j=0,jy+1
         do i=0,jx+1

            read(unit5,*)var1,var2,var3,var4,(var(n),n=1,nscal)

            uo5(i,j) = var1*index_out5(i,j)
            vo5(i,j) = var2*index_out5(i,j)
            wo5(i,j) = var3*index_out5(i,j)
            tkepom5(i,j,ntke)=var4*index_out5(i,j)
            do isc = 1,nscal
               rhovo5(isc,i,j) = var(isc)*index_rho5(i,j)
            end do
         end do
      end do
	   
   end if
   !-----------------------------------------------------------------------
   !     side   6
   if(infout6 /=0)then
      allocate(uo6(0:jx+1,0:jy+1))
      allocate(vo6(0:jx+1,0:jy+1))
      allocate(wo6(0:jx+1,0:jy+1))
      allocate(rhovo6(nscal,0:jx+1,0:jy+1))
      allocate(wcp6(jx,jy))

      uo6=0.
      vo6=0.
      wo6=0.
      rhovo6=0.
      wcp6=0.


      read(unit6,*)commento
      read(unit6,*)ntireal,ti_pom_fin
      n_ti_pom = int(ntireal)
      read(unit6,*)ti_pom_old

      allocate(tkepom6(0:jx+1,0:jy+1,n_ti_pom))
      tkepom6 = 0.

      do j=0,jy+1
         do i=0,jx+1

            read(unit6,*)var1,var2,var3,var4,(var(n),n=1,nscal)

            uo6(i,j) = var1*index_out6(i,j)
            vo6(i,j) = var2*index_out6(i,j)
            wo6(i,j) = var3*index_out6(i,j)
            tkepom6(i,j,ntke)=var4*index_out6(i,j)
            do isc = 1,nscal
               rhovo6(isc,i,j) = var(isc)*index_rho6(i,j)
            end do
         end do
      end do

   end if
   !-----------------------------------------------------------------------
   !c      ti_pom_old = ti_pom_old * 3600. ! in seconds
   !c      ti_pom_fin = ti_pom_fin * 3600.

   if(myid.eq.0)then
      write(*,*)'----------------------'
      write(*,*)'read pareti.dat '
      write(*,*)'POM time number:',n_ti_pom
      write(*,*)'POM initial time:',ti_pom_old
      write(*,*)'POM final time:',ti_pom_fin
      write(*,*)'----------------------'
   end if
   !-----------------------------------------------------------------------
   !
   if(potenziale==0)then
      !        CASE START WITH RESTART FILE
      !        allocation for nesting
      !        side 1
      if(infout1 /=0)then
         allocate(un1(0:jy+1,0:jz+1))
         allocate(vn1(0:jy+1,0:jz+1))
         allocate(wn1(0:jy+1,0:jz+1))
         allocate(rhovn1(nscal,0:jy+1,0:jz+1))
         un1 = 0.
         vn1 = 0.
         wn1 = 0.
         rhovn1 = 0.
      end if
      !        side 2
      if(infout2 /=0)then
         allocate(un2(0:jy+1,0:jz+1))
         allocate(vn2(0:jy+1,0:jz+1))
         allocate(wn2(0:jy+1,0:jz+1))
         allocate(rhovn2(nscal,0:jy+1,0:jz+1))
         un2 = 0.
         vn2 = 0.
         wn2 = 0.
         rhovn2 = 0.
      end if
      !        side 5
      if(infout5 /=0)then
         allocate(un5(0:jx+1,0:jy+1))
         allocate(vn5(0:jx+1,0:jy+1))
         allocate(wn5(0:jx+1,0:jy+1))
         allocate(rhovn5(nscal,0:jx+1,0:jy+1))
         un5 = 0.
         vn5 = 0.
         wn5 = 0.
         rhovn5 = 0.
      end if
      !        side 6
      if(infout6 /=0)then
         allocate(un6(0:jx+1,0:jy+1))
         allocate(vn6(0:jx+1,0:jy+1))
         allocate(wn6(0:jx+1,0:jy+1))
         allocate(rhovn6(nscal,0:jx+1,0:jy+1))
         un6 = 0.
         vn6 = 0.
         wn6 = 0.
         rhovn6 = 0.
      end if

      call restart(ti,dbbx)
      call leggi_pareti_pom(ti)
      call interpola_pareti_pom(ti)
   else ! if start with potential flow
      ti = ti_pom_old
      
      u=0.;v=0.;w=0. ! set a zero velocity field

      !         inizialize density field
      call interp

   !          if(myid==0)then
   !             kpst = kparasta
   !	     kpnd = kparaend + deepr
   !          elseif(myid==nproc-1)then
   !             kpst = kparasta - deepl
   !	     kpnd = kparaend
   !          else
   !             kpst = kparasta - deepl
   !	     kpnd = kparaend + deepr
   !          end if
      
   !          do k=kpst,kpnd
   !          do j=1,jy
   !          do i=1,jx
   !            do isc=1,nscal
   !	       rhov(isc,i,j,k)=rhovo1(isc,j,kparasta)
   !            end do
   !          end do
   !          end do
   !          end do
      
   !          do k=kpst,kpnd
   !          do i=1,jx
   !            do isc=1,nscal
   !	      rhov(isc,i,0,k)=rhov(isc,i,1,k)
   !            end do
   !          end do
   !          end do
	  	    
   end if
   !-----------------------------------------------------------------------
   !     set boundary conditions
   !     lateral values for velocity vector
   !
   !     side 1
   if(infout1 /= 0)then
      do k=kparasta,kparaend
         do j=1,jy
      
            u(0,j,k)=up1(j,k)
            v(0,j,k)=vp1(j,k)
            w(0,j,k)=wp1(j,k)
         
            do isc=1,nscal
               rhov(isc,0   ,j,k)=rhovp1(isc,j,k)
            end do
         end do
      end do
   end if
          
   !     side 1
   if(infout2 /= 0)then
      do k=kparasta,kparaend
         do j=1,jy

            u(jx+1,j,k)=up2(j,k)
            v(jx+1,j,k)=vp2(j,k)
            w(jx+1,j,k)=wp2(j,k)
         
            do isc=1,nscal
               rhov(isc,jx+1,j,k)=rhovp2(isc,j,k)
            end do
         end do
      end do
   end if

   !     sides 5
   if(infout5 /= 0)then
      if(myid==0)then
         do j=1,jy
            do i=1,jx
               u(i,j,0)=up5(i,j)
               v(i,j,0)=vp5(i,j)
               w(i,j,0)=wp5(i,j)
               do isc=1,nscal
                  rhov(isc,i,j,0   )=rhovp5(isc,i,j)
               end do
            end do
         end do
      end if
   end if
      
   !     sides 6
   if(infout6 /= 0)then
      if(myid==nproc-1)then
         do j=1,jy
            do i=1,jx
               u(i,j,jz+1)=up6(i,j)
               v(i,j,jz+1)=vp6(i,j)
               w(i,j,jz+1)=wp6(i,j)
               do isc=1,nscal
                  rhov(isc,i,j,jz+1)=rhovp6(isc,i,j)
               end do
            end do
         end do
      end if
   end if
   !-----------------------------------------------------------------------

   return
end


!***********************************************************************
subroutine interp
   !***********************************************************************
   use mysending
   use myarrays_metri3
   use myarrays_velo3
   use myarrays_nesting
   !
   use scala3
   use tipologia
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,isc
      
   real xp1,xp2
   real zp1,zp2
   real xc,zc
      
   real dist1,dist2,dist_tot
      
   real x_p0(0:jx,0:jy)
   real y_p0(0:jx,0:jy)
   real z_p0(0:jx,0:jy)

   real x_pjz(0:jx,0:jy)
   real y_pjz(0:jx,0:jy)
   real z_pjz(0:jx,0:jy)
      
   real r_ew,r_sn

   real buffer(3*(jx+1)*(jy+1))
   integer icount,numcount
   integer ierr
   !-----------------------------------------------------------------------

   numcount = (jx+1)*(jy+1)

   !     plane comunication
   !     plane 0
   buffer = 0.
   if(myid==0)then
      
      k=0
      do j=0,jy
         do i=0,jx
            x_p0(i,j)=x(i,j,k)
            y_p0(i,j)=y(i,j,k)
            z_p0(i,j)=z(i,j,k)
         end do
      end do
	     
      icount = 1
      do j=0,jy
         do i=0,jx
            buffer(           icount) = x_p0(i,j)
            buffer(  numcount+icount) = y_p0(i,j)
            buffer(2*numcount+icount) = z_p0(i,j)
	   
            icount = icount + 1
         end do
      end do
   end if
     
   call MPI_BCAST(buffer(1),3*numcount,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)

   icount = 1
   do j=0,jy
      do i=0,jx
         x_p0(i,j) = buffer(	       icount)
         y_p0(i,j) = buffer(  numcount+icount)
         z_p0(i,j) = buffer(2*numcount+icount)
         
         icount = icount + 1
      end do
   end do
      
   !-----------------------------------------------------------------------
   !     plane comunication
   !     plane jz
   buffer = 0.
   if(myid==nproc-1)then
      
      k=jz
      do j=0,jy
         do i=0,jx
            x_pjz(i,j)=x(i,j,k)
            y_pjz(i,j)=y(i,j,k)
            z_pjz(i,j)=z(i,j,k)
         end do
      end do
	      
      icount = 1
      do j=0,jy
         do i=0,jx
            buffer(           icount) = x_pjz(i,j)
            buffer(  numcount+icount) = y_pjz(i,j)
            buffer(2*numcount+icount) = z_pjz(i,j)
	   
            icount = icount + 1
         end do
      end do
   end if
      
   call MPI_BCAST(buffer(1),3*numcount,MPI_REAL_SD,nproc-1,MPI_COMM_WORLD,ierr)
 
   icount = 1
   do j=0,jy
      do i=0,jx
         x_pjz(i,j) = buffer(	        icount)
         y_pjz(i,j) = buffer(  numcount+icount)
         z_pjz(i,j) = buffer(2*numcount+icount)
         
         icount = icount + 1
      end do
   end do
   !
   !-----------------------------------------------------------------------
   !
   if(myid==0)write(*,*)'ENTRO'
   !     interpolation
   do j=1,jy
      do k=kparasta,kparaend
         do i=1,jx
            do isc= 1,nscal
           
               xp1 = .5*(x(0 ,j,k)+x(0 ,j,k-1))
               zp1 = .5*(z(0 ,j,k)+z(0 ,j,k-1))
      
               xp2 = .5*(x(jx,j,k)+x(jx,j,k-1))
               zp2 = .5*(z(jx,j,k)+z(jx,j,k-1))
      
               xc = .25*(x(i,j,k)+x(i,j,k-1)+x(i-1,j,k-1)+x(i-1,j,k))
               zc = .25*(z(i,j,k)+z(i,j,k-1)+z(i-1,j,k-1)+z(i-1,j,k))
      
               dist1 = sqrt((xc-xp1)*(xc-xp1) + (zc-zp1)*(zc-zp1))
               dist2 = sqrt((xc-xp2)*(xc-xp2) + (zc-zp2)*(zc-zp2))
               dist_tot = dist1 + dist2
      
               r_ew = (rhovo1(isc,j,k)*dist2 + rhovo2(isc,j,k)*dist1 )/ dist_tot
                 
               xp1 = .5*(x_p0(i,j)+x_p0(i-1,j))
               zp1 = .5*(z_p0(i,j)+z_p0(i-1,j))
      
               xp2 = .5*(x_pjz(i,j)+x_pjz(i-1,j))
               zp2 = .5*(z_pjz(i,j)+z_pjz(i-1,j))
            
               dist1 = sqrt((xc-xp1)*(xc-xp1) + (zc-zp1)*(zc-zp1))
               dist2 = sqrt((xc-xp2)*(xc-xp2) + (zc-zp2)*(zc-zp2))
               dist_tot = dist1 + dist2
      
               r_sn = (rhovo5(isc,i,j)*dist2 + rhovo6(isc,i,j)*dist1 )/ dist_tot
 
               rhov(isc,i,j,k) = 0.5*(r_ew+r_sn)
	       
            !	       if(isc ==1 .and. i==10 .and. j==10 .and. myid==0)then
            !	         write(*,*)k,r_sn,r_ew,
            !     >		   rhovo5(isc,i,j),rhov(isc,i,j,k),rhovo6(isc,i,j)
            !	       end if
            end do
         end do
      end do
   end do
   !
   !-----------------------------------------------------------------------
   !
   return
end
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
