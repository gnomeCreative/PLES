!***********************************************************************
subroutine inflow(myid,nproc,lett,bodyforce,area1,area2,area3,area4,area5,area6,kparasta,kparaend)
   !***********************************************************************
   use mysettings_boundary
   use myarrays_buffer_bodyforce
   use myarrays_velo3
   use myarrays2
   !
   use scala3
   use tipologia
   use orl
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   ! array declaration
   character*20 filepn
   integer i,j,k,isc,npn,myid,lett,ierr,nproc
   integer l,numero_celle_bloccate
   integer bodyforce,kparasta,kparaend
   integer nnx,nny,nnz

   real area1,area2,area3,area4,area5,area6

   real parziale1,parziale2
   real parziale3,parziale4
   real parziale5,parziale6

   integer cont1,cont2,cont3,cont4,cont5,cont6

   real disct1A(jy,jz)
   real disct1B,disct1,disct2,a1,a2,a3
   real xib,yib,zib,xv,yv,zv,fattore
   integer kk,pianoturbo,npiani

   real massa,massatot
      
   real val_tke1,val_tke2
   real val_tke5,val_tke6
   !-----------------------------------------------------------------------
   !     read the tke

   !     side 1
   if(ibodybuffer1==1)then
      open(121,file='tke1.dat',status='unknown')
      do k=1,jz
         do j=1,jy
            read(121,*)val_tke1
            if(k.ge. kparasta-1 .and. k.le. kparaend+1)then
               tke1(j,k) = val_tke1
            end if
         end do
      end do
      close(121)
   end if

   !     side 2
   if(ibodybuffer2==1)then
      open(121,file='tke1.dat',status='unknown')
      do k=1,jz
         do j=1,jy
            read(122,*)val_tke2
            if(k.ge. kparasta-1 .and. k.le. kparaend+1)then
               tke2(j,k) = val_tke2
            end if
         end do
      end do
      close(122)
   end if
      
   !     side 5
   if(myid==0)then
      if(ibodybuffer5==1)then
         open(125,file='tke5.dat',status='unknown')
         do j=1,jy
            do i=1,jx
               read(125,*)val_tke5
               tke5(i,j) = val_tke5
            end do
         end do
         close(125)
      end if
   end if

   !     side 6
   if(myid==nproc-1)then
      if(ibodybuffer6==1)then
         open(126,file='tke6.dat',status='unknown')
         do j=1,jy
            do i=1,jx
               read(126,*)val_tke6
               tke6(i,j) = val_tke6
            end do
         end do
         close(126)
      end if
   end if
      
   !-----------------------------------------------------------------------
   !    OUTFLOW
   !-----------------------------------------------------------------------
   ! in case of partial outflow
   ! if no outflow index = 0
   ! if    outflow index = 1
   !
   !-----------------------------------------------------------------------
   ! inflow

   !     planes at constant csi
   if(infout1.eq.0)then
	
      filepn='piano1.dat'
      open(10011,file=filepn,status='unknown')
      read(10011,*)nny,nnz !jy,jz
      read(10011,*)npn1
      read(10011,*)npn

      do k=1,jz
         do j=1,jy
            read(10011,100)up1(j,k),&
               vp1(j,k),&
               wp1(j,k),&
               (rhovp1(isc,j,k),isc=1,nscal)
     
            up1(j,k)=up1(j,k)*real(index_out1(j,k))
            vp1(j,k)=vp1(j,k)*real(index_out1(j,k))
            wp1(j,k)=wp1(j,k)*real(index_out1(j,k))
            do isc=1,nscal
               rhovp1(isc,j,k)=rhovp1(isc,j,k)*real(index_out1(j,k))
            end do
         end do
      end do
	

	
	
      if(myid.eq.0) then
         print*,myid,'inflow plane face 1'
      endif
   endif
100 format(10e13.5)

   if(infout2.eq.0)then
	
      filepn='piano2.dat'
      open(10012,file=filepn,status='unknown')
      read(10012,*)nny,nnz !jy,jz
      read(10012,*)npn2
      read(10012,*)npn

      do k=1,jz
         do j=1,jy
            read(10012,101)up2(j,k),&
               vp2(j,k),&
               wp2(j,k),&
               (rhovp2(isc,j,k),isc=1,nscal)
     
            up2(j,k)=up2(j,k)*real(index_out2(j,k))
            vp2(j,k)=vp2(j,k)*real(index_out2(j,k))
            wp2(j,k)=wp2(j,k)*real(index_out2(j,k))
            do isc=1,nscal
               rhovp2(isc,j,k)=rhovp2(isc,j,k)*real(index_out2(j,k))
            end do
         enddo
      enddo
      if(myid.eq.0) then
         print*,myid,'inflow plane face 2'
      endif
   endif
101 format(10e13.5)

   !     planes at constant eta
   if(infout3.eq.0)then

      filepn='piano3.dat'
      open(10013,file=filepn,status='unknown')
      read(10013,*)nnx,nnz !jx,jz
      read(10013,*)npn3
      read(10013,*)npn
      do k=1,jz
         do i=1,jx
            read(10013,102)up3(i,k),&
               vp3(i,k),&
               wp3(i,k),&
               (rhovp3(isc,i,k),isc=1,nscal)

            up3(i,k)=up3(i,k)*real(index_out3(i,k))
            vp3(i,k)=vp3(i,k)*real(index_out3(i,k))
            wp3(i,k)=wp3(i,k)*real(index_out3(i,k))
            do isc=1,nscal
               rhovp3(isc,i,k)=rhovp3(isc,i,k)*real(index_out3(i,k))
            end do
         enddo
      enddo
      if(myid.eq.0) then
         print*,myid,'inflow plane face 3'
      endif
   endif
102 format(10e13.5)
 
   if(infout4.eq.0)then

      filepn='piano4.dat'
      open(10014,file=filepn,status='unknown')
      read(10014,*)nnx,nnz !jx,jz
      read(10014,*)npn4
      read(10014,*)npn
	
      do k=1,jz
         do i=1,jx
            read(10014,103)up4(i,k),&
               vp4(i,k),&
               wp4(i,k),&
               (rhovp4(isc,i,k),isc=1,nscal)

            up4(i,k)=up4(i,k)*real(index_out4(i,k))
            vp4(i,k)=vp4(i,k)*real(index_out4(i,k))
            wp4(i,k)=wp4(i,k)*real(index_out4(i,k))
            do isc=1,nscal
               rhovp4(isc,i,k)=rhovp4(isc,i,k)*real(index_out4(i,k))
            end do
         enddo
      enddo
      if(myid.eq.0) then
         print*,myid,'inflow plane face 4'
      endif
   endif
103 format(10e13.5)

   !     planes at constant eta
   if(infout5.eq.0)then

      filepn='piano5.dat'
      open(10015,file=filepn,status='unknown')
      read(10015,*)nnx,nny !jx,jy
      read(10015,*)npn5
      read(10015,*)npn
	
      do j=1,jy
         do i=1,jx
            read(10015,104)up5(i,j),&
               vp5(i,j),&
               wp5(i,j),&
               (rhovp5(isc,i,j),isc=1,nscal)

            up5(i,j)=up5(i,j)*real(index_out5(i,j))
            vp5(i,j)=vp5(i,j)*real(index_out5(i,j))
            wp5(i,j)=wp5(i,j)*real(index_out5(i,j))
            do isc=1,nscal
               rhovp5(isc,i,j)=rhovp5(isc,i,j)*real(index_out5(i,j))
            end do
                           if (myid==nproc-1) then
               write(*,*)i*j,up5(i,j),&
               vp5(i,j),&
               wp5(i,j),&
               (rhovp5(isc,i,j),isc=1,nscal)
                end if
         enddo
      enddo
      if(myid.eq.0) then
         print*,myid,' inflow plane face 5'
      endif
   endif
104 format(10e13.5)

   if(infout6.eq.0)then

      filepn='piano6.dat'
      open(10016,file=filepn,status='unknown')
      read(10016,*)nnx,nny !jx,jy
      read(10016,*)npn6
      read(10016,*)npn
	
      do j=1,jy
         do i=1,jx
            read(10016,105)up6(i,j),&
               vp6(i,j),&
               wp6(i,j),&
               (rhovp6(isc,i,j),isc=1,nscal)

            up6(i,j)=up6(i,j)*real(index_out6(i,j))
            vp6(i,j)=vp6(i,j)*real(index_out6(i,j))
            wp6(i,j)=wp6(i,j)*real(index_out6(i,j))
            do isc=1,nscal
               rhovp6(isc,i,j)=rhovp6(isc,i,j)*real(index_out6(i,j))
            end do
         enddo
      enddo
      if(myid.eq.0) then
         print*,myid,'inflow plane face 6'
      endif
   endif
105 format(10e13.5)
   !-----------------------------------------------------------------------
   return
end
