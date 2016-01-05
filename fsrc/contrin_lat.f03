!***********************************************************************
subroutine contrin_lat
   !***********************************************************************
   ! compute controvariant fluxes from interpolated velocity field
   ! this subroutine is used for nesting
   !
   use myarrays_velo3
   use myarrays_metri3
   use myarrays_nesting
   use mysending
   !
   use scala3
   use period
   use orl

   implicit none
   !
   !
   !-----------------------------------------------------------------------
   ! array declaration
   integer i,j,k
   real diver,divmax
   !
   !-----------------------------------------------------------------------
   !
   ! flux J-1*U
   !
   !     side 1
   if(infout1 /= 0)then
      do k=kparasta,kparaend !1,jz
         do j=1,jy

            uc(0,j,k)=csx(0,j,k)*up1(j,k)+csy(0,j,k)*vp1(j,k)+csz(0,j,k)*wp1(j,k)
            ucp1(j,k)=uc(0,j,k)
    	 
         enddo
      enddo
   end if
   !    side 2
   if(infout2 /= 0)then
      do k=kparasta,kparaend !1,jz
         do j=1,jy
    	 
            uc(jx,j,k)=csx(jx,j,k)*up2(j,k)+csy(jx,j,k)*vp2(j,k)+csz(jx,j,k)*wp2(j,k)
            ucp2(j,k)=uc(jx,j,k)

         enddo
      enddo
   end if
   !
   ! flux J-1*V
   vc(:,0,:)=0.
   vc(:,jy,:)=0.
   !
   !
   ! flux J-1*W
   !
   !     side 5
   if(infout5 /= 0)then
      if(myid.eq.0)then
      
         do j=1,jy
            do i=1,jx
               wc(i,j,0)=ztx(i,j,0)*up5(i,j)+zty(i,j,0)*vp5(i,j)+ztz(i,j,0)*wp5(i,j)
               wcp5(i,j)=wc(i,j,0)
            enddo
         enddo
      
      end if
   end if

   !     side 6
   if(infout6 /= 0)then
      if(myid.eq.nproc-1)then

         do j=1,jy
            do i=1,jx
               wc(i,j,jz)=ztx(i,j,jz)*up6(i,j)+zty(i,j,jz)*vp6(i,j)+ztz(i,j,jz)*wp6(i,j)
               wcp6(i,j)=wc(i,j,jz)
            enddo
         enddo

      end if
   end if
   !
   return
end
