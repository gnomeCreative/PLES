!***********************************************************************
subroutine interpola_pareti_pom(ti)
   !***********************************************************************
   use myarrays_velo3
   use myarrays_nesting
   use myarrays_buffer_bodyforce
   use mysending
   !
   use scala3
   use orl

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,ii
   real ti,inv_timewindow
   real delta_l,delta_r
   !-----------------------------------------------------------------------
   !     define some constant values

   inv_timewindow = 1./(ti_pom_new-ti_pom_old)
   delta_l = ti - ti_pom_old
   delta_r = ti_pom_new - ti
      
   !     side 1
   if(infout1 /=0)then
      i=0
      do k=kparasta-1,kparaend+1 !0,jz+1
         do j=0,jy+1
            up1(j,k)=(delta_r*uo1(j,k)+delta_l*un1(j,k))*inv_timewindow*index_out1(j,k)
     
            vp1(j,k)=(delta_r*vo1(j,k)+delta_l*vn1(j,k))*inv_timewindow*index_out1(j,k)
     
            wp1(j,k)=(delta_r*wo1(j,k)+delta_l*wn1(j,k))*inv_timewindow*index_out1(j,k)
     
            tke1(j,k)=(delta_r*tkepom1(j,k,ntke-1)+delta_l*tkepom1(j,k,ntke))*inv_timewindow*index_out1(j,k)
     
            do ii=1,nscal
               rhovp1(ii,j,k)=(delta_r*rhovo1(ii,j,k)+delta_l*rhovn1(ii,j,k))*inv_timewindow*index_rho1(j,k)
            end do
       
         end do
      end do
   end if
             
   !     side 2
   if(infout2 /=0)then
      i=jx+1
      do k=kparasta-1,kparaend+1 !0,jz+1
         do j=0,jy+1
            up2(j,k)=(delta_r*uo2(j,k)+delta_l*un2(j,k))*inv_timewindow*index_out2(j,k)
          
            vp2(j,k)=(delta_r*vo2(j,k)+delta_l*vn2(j,k))*inv_timewindow*index_out2(j,k)
          
            wp2(j,k)=(delta_r*wo2(j,k)+delta_l*wn2(j,k))*inv_timewindow*index_out2(j,k)

            tke2(j,k)=(delta_r*tkepom2(j,k,ntke-1)+delta_l*tkepom2(j,k,ntke))*inv_timewindow*index_out2(j,k)

            do ii=1,nscal
               rhovp2(ii,j,k)=(delta_r*rhovo2(ii,j,k)+delta_l*rhovn2(ii,j,k))*inv_timewindow*index_rho2(j,k)
          
            end do
      
         end do
      end do
   end if
            
   !     side 5
   if(infout5 /=0)then
      if(myid==0)then
         k=0
         do j=0,jy+1
            do i=0,jx+1
               up5(i,j)=(delta_r*uo5(i,j)+delta_l*un5(i,j))*inv_timewindow*index_out5(i,j)
          
               vp5(i,j)=(delta_r*vo5(i,j)+delta_l*vn5(i,j))*inv_timewindow*index_out5(i,j)
          
               wp5(i,j)=(delta_r*wo5(i,j)+delta_l*wn5(i,j))*inv_timewindow*index_out5(i,j)
     
               tke5(i,j)=(delta_r*tkepom5(i,j,ntke-1)+delta_l*tkepom5(i,j,ntke))*inv_timewindow*index_out5(i,j)
     
               do ii=1,nscal
                  rhovp5(ii,i,j)=(delta_r*rhovo5(ii,i,j)+delta_l*rhovn5(ii,i,j))*inv_timewindow*index_rho5(i,j)
               end do
      
            end do
         end do
      end if
   end if
      
   !     side 6
   if(infout6 /=0)then
      if(myid==nproc-1)then
         k=jz+1
         do j=0,jy+1
            do i=0,jx+1
               up6(i,j)=(delta_r*uo6(i,j)+delta_l*un6(i,j))*inv_timewindow*index_out6(i,j)
          
               vp6(i,j)=(delta_r*vo6(i,j)+delta_l*vn6(i,j))*inv_timewindow*index_out6(i,j)
          
               wp6(i,j)=(delta_r*wo6(i,j)+delta_l*wn6(i,j))*inv_timewindow*index_out6(i,j)
          
               tke6(i,j)=(delta_r*tkepom6(i,j,ntke-1)+delta_l*tkepom6(i,j,ntke))*inv_timewindow*index_out6(i,j)
          
               do ii=1,nscal
                  rhovp6(ii,i,j)=(delta_r*rhovo6(ii,i,j)+delta_l*rhovn6(ii,i,j))*inv_timewindow*index_rho6(i,j)
               end do
      
            end do
         end do
      end if
   end if

   return
end
