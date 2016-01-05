!***********************************************************************
subroutine leggi_pareti_pom(ti)
   !***********************************************************************
   !     read input data from pom simulation, used in nesting procedure
   use myarrays_nesting
   use myarrays_buffer_bodyforce
   use mysending
   !
   use scala3
   use tipologia
   use orl
   !
   use mpi


   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,n,isc
      
   real ti
   real yyy
      
   real var1,var2,var3,var4,var(nscal)
      
   logical read_plane
   !-----------------------------------------------------------------------

   read_plane = .true.

   do while(read_plane)
      if(myid==0)write(*,*)'READ PLANE NESTING'

      ntke = ntke + 1

      ! read side 1
      if(infout1 /= 0)then

         read(  81,*)ti_pom_new

         do k=0,jz+1
            do j=0,jy+1
               read(  81,*)var1,var2,var3,var4,(var(n),n=1,nscal)
	    
               if(k.ge.kparasta.and.k.le.kparaend)then
                  un1(j,k) = var1*index_out1(j,k)
                  vn1(j,k) = var2*index_out1(j,k)
                  wn1(j,k) = var3*index_out1(j,k)
                  tkepom1(j,k,ntke)=var4*index_out1(j,k)
                  do isc = 1,nscal
                     rhovn1(isc,j,k) = var(isc)*index_rho1(j,k)
                  end do
               endif
            end do
         end do

      end if

      ! read side 2
      if(infout2 /= 0)then

         read(  82,*)ti_pom_new

         do k=0,jz+1
            do j=0,jy+1
               read(   82,*)var1,var2,var3,var4,(var(n),n=1,nscal)
               if(k.ge.kparasta.and.k.le.kparaend)then
                  un2(j,k) = var1*index_out2(j,k)
                  vn2(j,k) = var2*index_out2(j,k)
                  wn2(j,k) = var3*index_out2(j,k)
                  tkepom2(j,k,ntke)=var4*index_out2(j,k)
                  do isc = 1,nscal
                     rhovn2(isc,j,k) = var(isc)*index_rho2(j,k)
                  end do
               endif
            enddo
         enddo
 
      end if

      ! read side 5
      if(infout5 /= 0)then
         read(  85,*)ti_pom_new
	 
         do j=1,jy
            do i=1,jx
	 
               read( 85,*)var1,var2,var3,var4,(var(n),n=1,nscal)
	 
               un5(i,j) = var1*index_out5(i,j)
               vn5(i,j) = var2*index_out5(i,j)
               wn5(i,j) = var3*index_out5(i,j)
               tkepom5(i,j,ntke)=var4*index_out5(i,j)
               do isc = 1,nscal
                  rhovn5(isc,i,j) = var(isc)*index_rho5(i,j)
               end do
            end do
         end do
	 
      end if
	 
      ! read side 6
      if(infout6 /= 0)then
         read(  86,*)ti_pom_new
	 
         do j=1,jy
            do i=1,jx
	 
               read( 86,*)var1,var2,var3,var4,(var(n),n=1,nscal)
	 
               un6(i,j) = var1*index_out6(i,j)
               vn6(i,j) = var2*index_out6(i,j)
               wn6(i,j) = var3*index_out6(i,j)
               tkepom6(i,j,ntke)=var4*index_out6(i,j)
               do isc = 1,nscal
                  rhovn6(isc,i,j) = var(isc)*index_rho6(i,j)
               end do
            end do
         end do

      end if

      !cc         ti_pom_new = ti_pom_new * 3600.


      if(ti > ti_pom_new)then
         if(infout1 /= 0)then
            uo1=un1
            vo1=vn1
            wo1=wn1
            rhovo1=rhovn1
         end if
	     
         if(infout2 /= 0)then
            uo2=un2
            vo2=vn2
            wo2=wn2
            rhovo2=rhovn2
         end if
	    
         if(infout5 /= 0)then
            uo5=un5
            vo5=vn5
            wo5=wn5
            rhovo5=rhovn5
         end if
	   
         if(infout6 /= 0)then
            uo6=un6
            vo6=vn6
            wo6=wn6
            rhovo6=rhovn6
         end if

         ti_pom_old=ti_pom_new
      else
	    
         read_plane = .false.

      end if
	 	 
	 
   end do !do while

8753 format(6e14.8)
      
   return
end
