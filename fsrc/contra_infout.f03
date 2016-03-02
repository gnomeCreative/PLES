!***********************************************************************
subroutine contra_infout(kparasta,kparaend,rightpe,leftpe,nproc,myid,&
   area1,area2,area3,area4,area5,area6,ktime)
   !***********************************************************************
   ! compute Uc* intermediate controvariant quantities in case of
   ! inflow/outflow and boundary condition for pressure cs1 etc.
   !
   use myarrays_metri3
   use myarrays_velo3
   use myarrays2
   use output_module, only: info_run_file
   !
   use scala3
   use period
   use tipologia
   use orl
   !
   use mpi

   implicit none


   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,ii,jj,kk,iii,ktime
   real uinter,vinter,winter

   real area1,area2,area3,area4,area5,area6

   real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddi,dsi,ddj,dsj,ddk,dsk
   !
   integer ierr,myid,nproc
   integer ncolperproc,kparasta,kparaend,m
   integer kparastal,kparaendl
   integer status(MPI_STATUS_SIZE)
   integer rightpe,leftpe,rightpem,leftpem

   double precision am_in,am_in_loc
   double precision am_out,am_out_loc
   double precision am_out_loc1,am_out_loc2,am_out_loc3
   double precision am_out_loc4,am_out_loc5,am_out_loc6
   double precision am_out1,am_out2,am_out3
   double precision am_out4,am_out5,am_out6
   double precision am_out_ripartisco
      
   double precision am1,am2,am1t,am2t
   double precision am3,am4,am3t,am4t
   double precision am5,am6,am5t,am6t

   double precision giactot,giactot_loc
   double precision adel,adel_loc
   double precision del_mas

   double precision coef_massa1,coef_massa2,coef_massa3
   double precision coef_massa4,coef_massa5,coef_massa6

   double precision c_massa1(n2,n3),c_massa2(n2,n3)
   double precision c_massa3(n1,n3),c_massa4(n1,n3)
   double precision c_massa5(n1,n2),c_massa6(n1,n2)
      
   double precision contatore
      
   real area_bagnata
   double precision somma,sommat,amm,massa_tot
      
   integer faccio
      
   real ucc,vcc,wcc
   real area_insisto,fattore_flusso
   !
   !---------------------------------------------------------------------
   !
   contatore=0.
   if(infout1.eq.1)then
      contatore=contatore+1.
   endif
   if(infout2.eq.1)then
      contatore=contatore+1.
   endif
   if(infout3.eq.1)then
      contatore=contatore+1.
   endif
   if(infout4.eq.1)then
      contatore=contatore+1.
   endif
   if(infout5.eq.1)then
      contatore=contatore+1.
   endif
   if(infout6.eq.1)then
      contatore=contatore+1.
   endif
   if(myid.eq.0)then
      write(*,*)'INFLOW/OUTFLOW'
      write(info_run_file,*)'contatore pareti orlansky',contatore
   endif

   !-----------------------------------------------------------------------
   ! COMPUTE INCOMING MASS
   !-----------------------------------------------------------------------
   am_in=0.
   am_in_loc=0.
   del_mas=0.

   !     sides 1 and 2 constant csi
   do ii=1,ip

      if(infout1.ne.1)then
         do k=kparasta,kparaend
            do j=1,jy
               am_in_loc=am_in_loc+uc(0,j,k)
            end do
         end do
      endif

      if(infout2.ne.1)then
         do k=kparasta,kparaend
            do j=1,jy
               am_in_loc=am_in_loc-uc(jx,j,k)
            end do
         end do
      endif

   end do
   !
   !     sides 3 and 4 constant eta
   do jj=1,jp

      if(infout3.ne.1)then
         do k=kparasta,kparaend
            do i=1,jx
               am_in_loc=am_in_loc+vc(i,0,k)
            end do
         end do
      endif

      if(infout4.ne.1)then
         do k=kparasta,kparaend
            do i=1,jx
               am_in_loc=am_in_loc-vc(i,jy,k)
            end do
         end do
      endif

   end do
   !
   !     side 5 and 6 constant zita
   do kk=1,kp

      if(myid.eq.0)then
         if(infout5.ne.1)then
            do j=1,jy
               do i=1,jx
                  am_in_loc=am_in_loc+wc(i,j,0)
               end do
            end do
         endif

      elseif(myid.eq.nproc-1)then

         if(infout6.ne.1)then
            do j=1,jy
               do i=1,jx
                  am_in_loc=am_in_loc-wc(i,j,jz)
               end do
            end do
         endif
      endif

   end do
   !
   ! make the incoming mass known to all procs
   call MPI_ALLREDUCE(am_in_loc,am_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   !
   !-----------------------------------------------------------------------
   ! COMPUTE OUTFLOW MASS
   !-----------------------------------------------------------------------
   ! index initialization

   am_out_loc=0.
   am_out=0.

   am_out_loc1=0.
   am_out_loc2=0.
   am_out_loc3=0.
   am_out_loc4=0.
   am_out_loc5=0.
   am_out_loc6=0.

   am_out1=0.
   am_out2=0.
   am_out3=0.
   am_out4=0.
   am_out5=0.
   am_out6=0.
      
   am1=0.
   am1t=0.
            
   am2=0.
   am2t=0.

   am3=0.
   am3t=0.

   am4=0.
   am4t=0.

   am5=0.
   am5t=0.

   am6=0.
   am6t=0.
   !
   !     side 1 and 2 constant csi
   do ii=1,ip

      if(infout1.eq.1)then
         do k=kparasta,kparaend
            do j=1,jy

               uc(0,j,k)=du_dx1(j,k)*csx(0,j,k)+dv_dx1(j,k)*csy(0,j,k)+dw_dx1(j,k)*csz(0,j,k) !U^(n+1)  with orlansky

               am_out_loc1=am_out_loc1-uc(0,j,k)   ! outflow
               am1=am1+abs(uc(0,j,k))              ! module outflow
            end do
         end do
      endif

      if(infout2.eq.1)then
         do k=kparasta,kparaend
            do j=1,jy

               uc(jx,j,k)=du_dx2(j,k)*csx(jx,j,k)+dv_dx2(j,k)*csy(jx,j,k)+dw_dx2(j,k)*csz(jx,j,k)

               am_out_loc2=am_out_loc2+uc(jx,j,k)
               am2=am2+abs(uc(jx,j,k))
            end do
         end do
      endif

   end do
   !
   !     sides 3 and 4 constant eta
   do jj=1,jp

      if(infout3.eq.1)then
         do k=kparasta,kparaend
            do i=1,jx

               vc(i,0,k)=du_dy3(i,k)*etx(i,0,k)+dv_dy3(i,k)*ety(i,0,k)+dw_dy3(i,k)*etz(i,0,k)

               am_out_loc3=am_out_loc3-vc(i,0,k)
               am3=am3+abs(vc(i,0,k))
            end do
         end do
      endif

      if(infout4.eq.1)then
         do k=kparasta,kparaend
            do i=1,jx

               vc(i,jy,k)=du_dy4(i,k)*etx(i,jy,k)+dv_dy4(i,k)*ety(i,jy,k)+dw_dy4(i,k)*etz(i,jy,k)

               am_out_loc4=am_out_loc4+vc(i,jy,k)
               am4=am4+abs(vc(i,jy,k))
            end do
         end do
      endif

   end do
   !
   !     sides 5 and 6 constant zita
   do kk=1,kp

      if(myid.eq.0)then

         if(infout5.eq.1)then
            do j=1,jy
               do i=1,jx

                  wc(i,j,0)=du_dz5(i,j)*ztx(i,j,0)+dv_dz5(i,j)*zty(i,j,0)+dw_dz5(i,j)*ztz(i,j,0)

                  am_out_loc5=am_out_loc5-wc(i,j,0)
                  am5=am5+abs(wc(i,j,0))
               end do
            end do
         endif

      elseif(myid.eq.nproc-1)then

         if(infout6.eq.1)then
            do j=1,jy
               do i=1,jx

                  wc(i,j,jz)=du_dz6(i,j)*ztx(i,j,jz)+dv_dz6(i,j)*zty(i,j,jz)+dw_dz6(i,j)*ztz(i,j,jz)

                  am_out_loc6=am_out_loc6+wc(i,j,jz)
                  am6=am6+abs(wc(i,j,jz))
               end do
            end do
         endif

      endif

   end do

   am_out_loc=am_out_loc1+am_out_loc2+am_out_loc3+am_out_loc4+am_out_loc5+am_out_loc6

   !
   ! make outflow known to all procs
   call MPI_ALLREDUCE(am_out_loc,am_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   !
   !-----------------------------------------------------------------------
   ! COMPUTE MASS FRACTION
   !-----------------------------------------------------------------------
   !
   ! mass fraction for each sides, for distribution on cell area
   !
   !
   ! mass for each face

   call MPI_ALLREDUCE(am_out_loc1,am_out1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc2,am_out2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc3,am_out3,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc4,am_out4,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc5,am_out5,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc6,am_out6,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   ! mass for each face in modulus
   
   call MPI_ALLREDUCE(am1,am1t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am2,am2t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am3,am3t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am4,am4t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am5,am5t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am6,am6t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     
   amm = am1t + am2t + am3t + am4t + am5t + am6t
   !
   ! compute area with outflow
   area_bagnata=0.
   if(infout1.eq.1)then
      area_bagnata = area_bagnata+area_bagnata1
   end if
   if(infout2.eq.1)then
      area_bagnata = area_bagnata+area_bagnata2
   end if
   if(infout3.eq.1)then
      area_bagnata = area_bagnata+area_bagnata3
   end if
   if(infout4.eq.1)then
      area_bagnata = area_bagnata+area_bagnata4
   end if
   if(infout5.eq.1)then
      area_bagnata = area_bagnata+area_bagnata5
   end if
   if(infout6.eq.1)then
      area_bagnata = area_bagnata+area_bagnata6
   end if
      
   if(myid.eq.0)then
      write(*,*)'area_bagnata totale',area_bagnata
   end if
   !
   !-----------------------------------------------------------------------
   ! MASS DEFECT
   !-----------------------------------------------------------------------
   !
   del_mas=am_out-am_in
   !
   !-----------------------------------------------------------------------
   ! OUTPUT
   !-----------------------------------------------------------------------
   !
   if(myid.eq.0)then
      write(*,*)'inflow',am_in
      write(*,*)'outflow',am_out
      write(*,*)'delta',del_mas
      write(info_run_file,*)'massa attraverso faccia 1:',am_out1
      write(info_run_file,*)'massa attraverso faccia 2:',am_out2
      write(info_run_file,*)'massa attraverso faccia 3:',am_out3
      write(info_run_file,*)'massa attraverso faccia 4:',am_out4
      write(info_run_file,*)'massa attraverso faccia 5:',am_out5
      write(info_run_file,*)'massa attraverso faccia 6:',am_out6
   endif
   !
   !-----------------------------------------------------------------------
   ! MASS CORRECTION ON CONTROVARIANT FLUX
   !-----------------------------------------------------------------------
   !
   somma = 0.

   massa_tot = abs(am_out1)+abs(am_out2)+abs(am_out3)+abs(am_out4)+abs(am_out5)+abs(am_out6)
          

   do ii=1,ip
      
      !     side 1 constant csi
      if(massa_tot .lt. 1.d-8)then
         area_insisto   =   area_bagnata
         fattore_flusso =   1.
      else     
         area_insisto   =   area_bagnata1
         fattore_flusso =  (abs(am_out1)/massa_tot)
      end if


      if(infout1 == 1)then
      
         do k=kparasta,kparaend
            do j=1,jy
               !
               coef_massa1=(areola1(j,k)/area_insisto)*fattore_flusso
		
               somma=somma+coef_massa1*index_out1(j,k)
               !
               !          Uc^n+1
               ucc=uc(0,j,k)+del_mas*coef_massa1*index_out1(j,k)

               uc1_orl(j,k)=ucc
               !          pressure b.c. = Uc* - Uc^n+1
               cs1(j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)-ucc
               if(potenziale==1)then
                  cs1(j,k)=-ucc
               end if
               !	   Uc*
               uc(0,j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)
	   
            enddo
         enddo
      
      else
      
         do k=kparasta,kparaend
            do j=1,jy
               !          Uc stored to compute the divg
               uc1_orl(j,k)=uc(0,j,k)
	
               cs1(j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)-uc(0,j,k)
               !          pressure b.c. = Uc* - Uc^n+1
               if(potenziale==1)then
                  cs1(j,k)= -uc(0,j,k)
               end if
               !	   Uc*
               uc(0,j,k)=u(0,j,k)*csx(0,j,k)+v(0,j,k)*csy(0,j,k)+w(0,j,k)*csz(0,j,k)
        
            !
            enddo
         enddo
      
      end if !infout1
      !     ..................................................................
      !     side 2 constant csi
      if(massa_tot .lt. 1.d-8)then
         area_insisto   =   area_bagnata
         fattore_flusso =   1.
      else     
         area_insisto   =   area_bagnata2
         fattore_flusso =  (abs(am_out2)/massa_tot)
      end if

      if(infout2 == 1)then
      
         do k=kparasta,kparaend
            do j=1,jy
               !
               coef_massa2=(areola2(j,k)/area_insisto)*fattore_flusso
               !
               somma=somma+coef_massa2*index_out2(j,k)

               !          Uc^n+1
               ucc=uc(jx,j,k)-del_mas*coef_massa2*index_out2(j,k)
               !
               uc2_orl(j,k)=ucc
	   
               !          pressure b.c. = Uc* - Uc^n+1
               cs2(j,k)=u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)-ucc
               if(potenziale==1)then
                  cs2(j,k)=-ucc
               end if
               !	   Uc*
               uc(jx,j,k) = u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)
	   	   
            enddo
         enddo
       
      else 
	
         do k=kparasta,kparaend
            do j=1,jy
               !          Uc stored to compute the divg
               uc2_orl(j,k)=uc(jx,j,k)
               !          pressure b.c. = Uc* - Uc^n+1
               cs2(j,k)=u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)-uc(jx,j,k)

               if(potenziale==1)then
                  cs2(j,k)=-uc(jx,j,k)
               end if
               !	   Uc*
               uc(jx,j,k)=u(jx+1,j,k)*csx(jx,j,k)+v(jx+1,j,k)*csy(jx,j,k)+w(jx+1,j,k)*csz(jx,j,k)

            enddo
         enddo
      
      end if !infout2

   enddo   !end loop ii=1,ip
   !.......................................................................
      
   do jj=1,jp

      if(massa_tot .lt. 1.d-8)then
         area_insisto   =   area_bagnata
         fattore_flusso =   1.
      else     
         area_insisto   =   area_bagnata3
         fattore_flusso =  (abs(am_out3)/massa_tot)
      end if

      !     side 3 constant eta
      if(infout3 == 1)then 
      
         do k=kparasta,kparaend
            do i=1,jx
	       
               coef_massa3=(areola3(i,k)/area_insisto)*fattore_flusso
	       
               somma=somma+coef_massa3*index_out3(i,k)
               !          Vc^n+1
               vcc=vc(i,0,k)+del_mas*coef_massa3*index_out3(i,k)

               vc3_orl(i,k)=vcc
               !          pressure b.c. = Vc* - Vc^n+1
               cs3(i,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)-vcc

               if(potenziale==1)then
                  cs3(i,k)=-vcc
               end if
               !	   Vc*
               vc(i,0,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)

            enddo
         enddo
      
      else
      	   	     
         do k=kparasta,kparaend
            do i=1,jx
               !          Vc stored to compute the divg
               vc3_orl(i,k)=vc(i,0,k)
               !          pressure b.c. = Vc* - Vc^n+1
               cs3(i,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)-vc(i,0,k)

               if(potenziale==1)then
                  cs3(i,k)=-vc(i,0,k)
               end if
               !	   Vc*
               vc(i,0,k)=u(i,0,k)*etx(i,0,k)+v(i,0,k)*ety(i,0,k)+w(i,0,k)*etz(i,0,k)
		     
            enddo
         enddo
      
      end if !infout3
      !     ..................................................................
      !     side 4 constant eta
      if(massa_tot .lt. 1.d-8)then
         area_insisto   =   area_bagnata
         fattore_flusso =   1.
      else     
         area_insisto   =   area_bagnata4
         fattore_flusso =  (abs(am_out4)/massa_tot)
      end if
      
      if(infout4 == 1)then
      
         do k=kparasta,kparaend
            do i=1,jx
	
               coef_massa4=(areola4(i,k)/area_insisto)*fattore_flusso
	
               somma=somma+coef_massa4*index_out4(i,k)
               !          Vc^n+1
               vcc=vc(i,jy,k)-del_mas*coef_massa4*index_out4(i,k)

               vc4_orl(i,k)=vcc
               !          pressure b.c. = Vc* - Vc^n+1
               cs4(i,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)-vcc
               if(potenziale==1)then
                  cs4(i,k)=-vcc
               end if
               !	   Vc*
               vc(i,jy,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)

            enddo
         enddo
      
      else
      
         do k=kparasta,kparaend
            do i=1,jx
               !          Vc stored to compute the divg
               vc4_orl(i,k)=vc(i,jy,k)
               !          pressure b.c. = Vc* - Vc^n+1
               cs4(i,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)-vc(i,jy,k)

               if(potenziale==1)then
                  cs4(i,k)=-vc(i,jy,k)
               end if
               !	   Vc*
               vc(i,jy,k)=u(i,jy+1,k)*etx(i,jy,k)+v(i,jy+1,k)*ety(i,jy,k)+w(i,jy+1,k)*etz(i,jy,k)
          
            enddo
         enddo
      
      end if !infout4

   enddo    !end loop jj=1,jp
   !.....................................
      
   do kk=1,kp
      !     side 5 constant zita
      if(myid==0)then

         if(massa_tot .lt. 1.d-8)then
            area_insisto   =   area_bagnata
            fattore_flusso =   1.
         else
            area_insisto   =   area_bagnata5
            fattore_flusso =  (abs(am_out5)/massa_tot)
         end if

         if(infout5 == 1)then

            do j=1,jy
               do i=1,jx
                  !
                  coef_massa5=(areola5(i,j)/area_insisto)*fattore_flusso
                  !
                  somma=somma+coef_massa5*index_out5(i,j)
                  !          Wc^n+1
                  wcc=wc(i,j,0)+del_mas*coef_massa5*index_out5(i,j)

                  wc5_orl(i,j)=wcc
                  !          pressure b.c. = Wc* - Wc^n+1
                  cs5(i,j)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)-wcc
                  if(potenziale==1)then
                     cs5(i,j)=-wcc
                  end if
                  !	   Wc*
                  wc(i,j,0)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)
               !
               end do
            end do

         else

            do j=1,jy
               do i=1,jx
                  !          wc stored to compute the divg
                  wc5_orl(i,j)=wc(i,j,0)
                  !          pressure b.c. = Wc* - Wc^n+1
                  cs5(i,j)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)-wc(i,j,0)

                  if(potenziale==1)then
                     cs5(i,j)=-wc(i,j,0)
                  end if
                  !	   Wc*
                  wc(i,j,0)=u(i,j,0)*ztx(i,j,0)+v(i,j,0)*zty(i,j,0)+w(i,j,0)*ztz(i,j,0)
               !
               enddo
            enddo

         end if !infout5

      end if !myid

      !     side 6 constant zita
      if(myid.eq.nproc-1)then

         if(massa_tot .lt. 1.d-8)then
            area_insisto   =   area_bagnata
            fattore_flusso =   1.
         else
            area_insisto   =   area_bagnata6
            fattore_flusso =  (abs(am_out6)/massa_tot)
         end if

         if(infout6 == 1)then

            do j=1,jy
               do i=1,jx
	 
                  coef_massa6=(areola6(i,j)/area_insisto)*fattore_flusso
                  !
                  somma=somma+coef_massa6*index_out6(i,j)
                  !          Wc^n+1
                  wcc=wc(i,j,jz)-del_mas*coef_massa6*index_out6(i,j)

                  wc6_orl(i,j)=wcc
                  !          pressure b.c. = Wc* - Wc^n+1
                  cs6(i,j)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)-wc(i,j,jz)
                  if(potenziale==1)then
                     cs6(i,j)=-wcc
                  end if
                  !	   Wc*
                  wc(i,j,jz)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)
               end do
            end do
      
         else
     
            do j=1,jy
               do i=1,jx
                  !          Wc stored to compute the divg
                  wc6_orl(i,j)=wc(i,j,jz)
                  !          pressure b.c. = Wc* - Wc^n+1
                  cs6(i,j)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)-wc(i,j,jz)

                  if(potenziale==1)then
                     cs6(i,j)=-wc(i,j,jz)
                  end if
                  !	   Wc*
                  wc(i,j,jz)=u(i,j,jz+1)*ztx(i,j,jz)+v(i,j,jz+1)*zty(i,j,jz)+w(i,j,jz+1)*ztz(i,j,jz)
               !
               enddo
            enddo
      
         end if !infout6
      
      end if !myid
   !
   enddo    !end loop kk=1,kp
   !
   !
   !
   call MPI_ALLREDUCE(somma,sommat,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   !
   !
   if(myid.eq.0)then
      write(*,*)myid,'sum coef',sommat
   end if
    
   !
   !-----------------------------------------------------------------------
   ! INSIDE THE FIELD
   !-----------------------------------------------------------------------
   !
   !     compute Uc*
   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
            !
            !
            uinter=.5*(u(i,j,k)+u(i+1,j,k))
            vinter=.5*(v(i,j,k)+v(i+1,j,k))
            winter=.5*(w(i,j,k)+w(i+1,j,k))
            !
            uc(i,j,k)=uinter*csx(i,j,k)+vinter*csy(i,j,k)+winter*csz(i,j,k)
         !
         end do
      end do
   end do

   !     compute Vc*
   do k=kparasta,kparaend
      do j=jp,jy-jp
         do i=1,jx

            !
            uinter=.5*(u(i,j,k)+u(i,j+1,k))
            vinter=.5*(v(i,j,k)+v(i,j+1,k))
            winter=.5*(w(i,j,k)+w(i,j+1,k))
            !
            vc(i,j,k)=uinter*etx(i,j,k)+vinter*ety(i,j,k)+winter*etz(i,j,k)
         !
         end do
      end do
   end do

   !     compute Wc*
   if (myid.eq.0) then
      kparastal=kp
      kparaendl=kparaend
   else if (myid.eq.nproc-1) then
      kparastal=kparasta
      kparaendl=kparaend-kp
   else
      kparastal=kparasta
      kparaendl=kparaend
   endif

   do k=kparastal,kparaendl

      do j=1,jy
         do i=1,jx
            !
            !
            uinter=.5*(u(i,j,k)+u(i,j,k+1))
            vinter=.5*(v(i,j,k)+v(i,j,k+1))
            winter=.5*(w(i,j,k)+w(i,j,k+1))
            !
            wc(i,j,k)=uinter*ztx(i,j,k)+vinter*zty(i,j,k)+winter*ztz(i,j,k)
         !
         enddo
      enddo

   enddo
   !
   !-----------------------------------------------------------------------
   ! SEND PLANE wc(k-1)
   !-----------------------------------------------------------------------
   ! send kparaend plane to closer proc in order to compute rhs in diver
   !
   if(myid.eq.0)then
      leftpem=MPI_PROC_NULL
      rightpem=rightpe
   else if(myid.eq.nproc-1)then
      leftpem=leftpe
      rightpem=MPI_PROC_NULL
   else if((myid.ne.0).and.(myid.ne.nproc-1))then
      leftpem=leftpe
      rightpem=rightpe
   endif


   call MPI_SENDRECV(wc(1,1,kparaend),jx*jy,MPI_REAL_SD,rightpem,51+myid,&
      wc(1,1,kparasta-1),jx*jy,MPI_REAL_SD,leftpem,50+myid,MPI_COMM_WORLD,status,ierr)
   !
   !-----------------------------------------------------------------------
   ! MASS CHECK
   !-----------------------------------------------------------------------
   ! check the mass balance after mass distribution

   am_out_loc=0.
   am_out=0.

   am_out_loc1=0.
   am_out_loc2=0.
   am_out_loc3=0.
   am_out_loc4=0.
   am_out_loc5=0.
   am_out_loc6=0.

   am_out1=0.
   am_out2=0.
   am_out3=0.
   am_out4=0.
   am_out5=0.
   am_out6=0.

   !     sides 1 and 2 constant csi
   do ii=1,ip

      if(infout1==1)then
         do k=kparasta,kparaend
            do j=1,jy

               am_out_loc1=am_out_loc1-uc1_orl(j,k)

            end do
         end do
      endif

      if(infout2==1)then
         do k=kparasta,kparaend
            do j=1,jy

               am_out_loc2=am_out_loc2+uc2_orl(j,k)

            end do
         end do
      endif

   end do
   !
   !     sides 3 and 4 constant eta
   do jj=1,jp

      if(infout3==1)then
         do k=kparasta,kparaend
            do i=1,jx

               am_out_loc3=am_out_loc3-vc3_orl(i,k)

            end do
         end do
      endif

      if(infout4==1)then
         do k=kparasta,kparaend
            do i=1,jx

               am_out_loc4=am_out_loc4+vc4_orl(i,k)

            end do
         end do
      endif

   end do
   !
   !     sides 5 and 6 constant zita
   do kk=1,kp

      if(myid==0)then

         if(infout5==1)then
            do j=1,jy
               do i=1,jx

                  am_out_loc5=am_out_loc5-wc5_orl(i,j)

               end do
            end do
         endif

      elseif(myid==nproc-1)then

         if(infout6==1)then
            do j=1,jy
               do i=1,jx

                  am_out_loc6=am_out_loc6+wc6_orl(i,j)

               end do
            end do
         endif

      endif

   end do

   am_out_loc=am_out_loc1+am_out_loc2+am_out_loc3+am_out_loc4+am_out_loc5+am_out_loc6

   !
   ! all procs know outflow
   call MPI_ALLREDUCE(am_out_loc,am_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   !

   ! outflow from each face
   call MPI_ALLREDUCE(am_out_loc1,am_out1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc2,am_out2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc3,am_out3,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc4,am_out4,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc5,am_out5,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   call MPI_ALLREDUCE(am_out_loc6,am_out6,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

   if(myid.eq.0)then
      !        write(*,*)'verifica massa entrante',am_in
      !        write(*,*)'verifica distribuzione massa uscente',am_out
      write(*,*)'check delta ',am_out-am_in
   endif


   return
end
