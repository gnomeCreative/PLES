!***********************************************************************
subroutine facce(myid,nproc,kparasta,kparaend, &
   area1,area2,area3,area4,area5,area6)
   !***********************************************************************
   !     find area of cells at the sides

   use myarrays2
   use myarrays_metri3
   use scala3
   use tipologia
   use output_module, only: info_run_file
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     variables declaration
   integer ierr,myid,nproc,kparasta,kparaend
   integer i,j,k    !,jx,jy,jz,n1,n2,n3
   integer status(mpi_status_size)
   real area1,area2,area3,area4,area5,area6
   real area1_loc, area2_loc, area3_loc
   real area4_loc, area5_loc, area6_loc
   real dxdeta,dydeta,dzdeta,dxdzet,dydzet,dzdzet
   real dxdcsi,dydcsi,dzdcsi
   real ra, sump
   !
   !-----------------------------------------------------------------------
   ! area face 1
   i = 0
   area1     = 0.
   area1_loc = 0.

   do k=kparasta,kparaend
      do j=1,n2
         areola1(j,k) = 0.
      end do
   end do

   do k=kparasta,kparaend
      do j=1,n2

         dxdeta = .5*( x(i,j,k)  + x(i,j,k-1)   )- &
            .5*( x(i,j-1,k)+ x(i,j-1,k-1) )

         dydeta = .5*( y(i,j,k)  + y(i,j,k-1)   )- &
            .5*( y(i,j-1,k)+ y(i,j-1,k-1) )

         dzdeta = .5*( z(i,j,k)  + z(i,j,k-1)   )- &
            .5*( z(i,j-1,k)+ z(i,j-1,k-1) )

         dxdzet = .5*( x(i,j,k) + x(i,j-1,k)    )- &
            .5*( x(i,j,k-1)+x(i,j-1,k-1)  )

         dydzet = .5*( y(i,j,k) + y(i,j-1,k)    )- &
            .5*( y(i,j,k-1)+y(i,j-1,k-1)  )

         dzdzet = .5*( z(i,j,k) + z(i,j-1,k)    )- &
            .5*( z(i,j,k-1)+z(i,j-1,k-1)  )

         gg221(j,k) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
         gg331(j,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
         gg231(j,k) = dxdeta*dxdzet + dydeta*dydzet + dzdeta*dzdzet

         areola1(j,k)=sqrt(gg221(j,k)*gg331(j,k)-gg231(j,k)*gg231(j,k) )

      end do
   end do

   do k=kparasta,kparaend
      do j=1,n2
         area1_loc = area1_loc + areola1(j,k)
      end do
   end do
   !
   ! areas sum, the sum is made by myid=0
   !
   if (myid.eq.0) then
      area1 = area1_loc
      sump  = 0.
      ra    = 0.
   end if

   do i=1,nproc-1

      if (myid.eq.0) then

         call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
         sump = sump + ra

      else if(myid.eq.i) then

         call mpi_send(area1_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)

      end if

   end do

   if (myid.eq.0) then
      area1 = area1 + sump
      write(*,*)'area face 1 =',area1
   end if
   !
   ! myid 0 send the computed total area to all procs
   ! mpi broadcast

   call mpi_bcast(area1,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

   write(info_run_file,*)'PE=',myid,'area1=',area1

   !
   !-----------------------------------------------------------------------
   ! area face 2

   i = n1
   area2     = 0.
   area2_loc = 0.


   do k=kparasta,kparaend
      do j=1,n2
         areola2(j,k) = 0.
      end do
   end do

   do k=kparasta,kparaend
      do j=1,n2

         dxdeta = .5*( x(i,j,k)  + x(i,j,k-1)   )- &
            .5*( x(i,j-1,k)+ x(i,j-1,k-1) )

         dydeta = .5*( y(i,j,k)  + y(i,j,k-1)   )- &
            .5*( y(i,j-1,k)+ y(i,j-1,k-1) )

         dzdeta = .5*( z(i,j,k)  + z(i,j,k-1)   )- &
            .5*( z(i,j-1,k)+ z(i,j-1,k-1) )

         dxdzet = .5*( x(i,j,k) + x(i,j-1,k)    )- &
            .5*( x(i,j,k-1)+x(i,j-1,k-1)  )

         dydzet = .5*( y(i,j,k) + y(i,j-1,k)    )- &
            .5*( y(i,j,k-1)+y(i,j-1,k-1)  )

         dzdzet = .5*( z(i,j,k) + z(i,j-1,k)    )- &
            .5*( z(i,j,k-1)+z(i,j-1,k-1)  )

         gg221(j,k) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
         gg331(j,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
         gg231(j,k) = dxdeta*dxdzet + dydeta*dydzet + dzdeta*dzdzet

         areola2(j,k)=sqrt(gg221(j,k)*gg331(j,k)-gg231(j,k)*gg231(j,k))

      end do
   end do

   do k=kparasta,kparaend
      do j=1,n2
         area2_loc = area2_loc + areola2(j,k)
      end do
   end do


   ! sum to have total area

   if (myid.eq.0) then
      area2 = area2_loc
      sump  = 0.
      ra    = 0.
   end if

   do i=1,nproc-1

      if (myid.eq.0) then
         call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
         sump = sump + ra

      else if(myid.eq.i) then
         call mpi_send(area2_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)
 
      end if

   end do

   if (myid.eq.0) then
      area2 = area2 +sump
      write(*,*)'area face 2 =',area2
   end if

   call mpi_bcast(area2,1,MPI_REAL_SD,0,mpi_comm_world,ierr)
      
   write(info_run_file,*)'PE=',myid,'area2=',area2

   !
   !-----------------------------------------------------------------------
   ! area face 5

   if (myid.eq.0) then

      k = 0 
      area5     = 0.
      area5_loc = 0.

      do j=1,n2
         do i=1,n1
            areola5(i,j) = 0.
         end do
      end do

      do j=1,n2
         do i=1,n1

            dxdcsi = .5*( x(i,j,k)   + x(i,j-1,k)   )- &
               .5*( x(i-1,j,k) + x(i-1,j-1,k) )
            dydcsi = .5*( y(i,j,k)   + y(i,j-1,k)   )- &
               .5*( y(i-1,j,k) + y(i-1,j-1,k) )
            dzdcsi = .5*( z(i,j,k)   + z(i,j-1,k)   )- &
               .5*( z(i-1,j,k) + z(i-1,j-1,k) )
            dxdeta = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
               .5*( x(i,j-1,k) + x(i-1,j-1,k) )
            dydeta = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
               .5*( y(i,j-1,k) + y(i-1,j-1,k) )
            dzdeta = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
               .5*( z(i,j-1,k) + z(i-1,j-1,k) )

            gg112(i,j) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
            gg222(i,j) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
            gg122(i,j) = dxdcsi*dxdeta + dydcsi*dydeta + dzdcsi*dzdeta

            areola5(i,j)=sqrt(gg112(i,j)*gg222(i,j)-gg122(i,j)*gg122(i,j))

         end do
      end do
      
      do j=1,n2
         do i=1,n1
            area5_loc = area5_loc + areola5(i,j)
         end do
      end do

      area5 = area5_loc
     
      write(info_run_file,*)'area face 5 =',area5

   end if ! close if PE=0

   call mpi_bcast(area5,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

   write(info_run_file,*)'PE=',myid,'area5=',area5

   !-----------------------------------------------------------------------
   ! area face 6

   if (myid.eq.(nproc-1)) then

      k = n3
      area6     = 0.
      area6_loc = 0.

      do j=1,n2
         do i=1,n1
            areola6(i,j) = 0.
         end do
      end do

      do j=1,n2
         do i=1,n1

            dxdcsi = .5*( x(i,j,k)   + x(i,j-1,k)   )- &
               .5*( x(i-1,j,k) + x(i-1,j-1,k) )
            dydcsi = .5*( y(i,j,k)   + y(i,j-1,k)   )- &
               .5*( y(i-1,j,k) + y(i-1,j-1,k) )
            dzdcsi = .5*( z(i,j,k)   + z(i,j-1,k)   )- &
               .5*( z(i-1,j,k) + z(i-1,j-1,k) )
            dxdeta = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
               .5*( x(i,j-1,k) + x(i-1,j-1,k) )
            dydeta = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
               .5*( y(i,j-1,k) + y(i-1,j-1,k) )
            dzdeta = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
               .5*( z(i,j-1,k) + z(i-1,j-1,k) )

            gg112(i,j) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
            gg222(i,j) = dxdeta*dxdeta + dydeta*dydeta + dzdeta*dzdeta
            gg122(i,j) = dxdcsi*dxdeta + dydcsi*dydeta + dzdcsi*dzdeta

            areola6(i,j)=sqrt(gg112(i,j)*gg222(i,j)-gg122(i,j)*gg122(i,j))

         end do
      end do
      
      do j=1,n2
         do i=1,n1
            area6_loc = area6_loc + areola6(i,j)
         end do
      end do

      area6 = area6_loc

      write(info_run_file,*)'area della faccia 6 =',area6

   end if ! close if PE=nproc-1

   call mpi_bcast(area6,1,MPI_REAL_SD,nproc-1,mpi_comm_world,ierr)

   write(info_run_file,*)'PE=',myid,'area6=',area6

   !
   !-----------------------------------------------------------------------
   ! area face 3

   j = 0
   area3     = 0.
   area3_loc = 0.

   do k=kparasta,kparaend
      do i=1,n1
         areola3(i,k) = 0.
      end do
   end do

   do k=kparasta,kparaend
      do i=1,n1

         dxdcsi = .5*( x(i,j,k)   + x(i,j,k-1)   )- &
            .5*( x(i-1,j,k) + x(i-1,j,k-1) )
         dydcsi = .5*( y(i,j,k)   + y(i,j,k-1)   )- &
            .5*( y(i-1,j,k) + y(i-1,j,k-1) )
         dzdcsi = .5*( z(i,j,k)   + z(i,j,k-1)   )- &
            .5*( z(i-1,j,k) + z(i-1,j,k-1) )
         dxdzet = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
            .5*( x(i,j,k-1) + x(i-1,j,k-1) )
         dydzet = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
            .5*( y(i,j,k-1) + y(i-1,j,k-1) )
         dzdzet = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
            .5*( z(i,j,k-1) + z(i-1,j,k-1) )
      
         gg113(i,k) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
         gg333(i,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
         gg133(i,k) = dxdcsi*dxdzet + dydcsi*dydzet + dzdcsi*dzdzet

         areola3(i,k)=sqrt(gg113(i,k)*gg333(i,k)-gg133(i,k)*gg133(i,k))

      end do
   end do

   do k=kparasta,kparaend
      do i=1,n1
         area3_loc = area3_loc + areola3(i,k)
      end do
   end do
    
   ! call to sum total area

   if (myid.eq.0) then
      area3 = area3_loc
      sump  = 0.
      ra    = 0.
   end if

   do i=1,nproc-1

      if (myid.eq.0) then
         call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
         sump = sump + ra

      else if(myid.eq.i) then
         call mpi_send(area3_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)

      end if

   end do

   if (myid.eq.0) then
      area3 = area3 + sump
      write(*,*)'area face 3 =',area3
   end if

   call mpi_bcast(area3,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

   write(info_run_file,*)'PE=',myid,'area3=',area3

   !
   !-----------------------------------------------------------------------
   ! area face 4

   j = n2
   area4     = 0.
   area4_loc = 0.

   do k=kparasta,kparaend
      do i=1,n1
         areola4(i,k) = 0.
      end do
   end do

   do k=kparasta,kparaend
      do i=1,n1

         dxdcsi = .5*( x(i,j,k)   + x(i,j,k-1)   )- &
            .5*( x(i-1,j,k) + x(i-1,j,k-1) )
         dydcsi = .5*( y(i,j,k)   + y(i,j,k-1)   )- &
            .5*( y(i-1,j,k) + y(i-1,j,k-1) )
         dzdcsi = .5*( z(i,j,k)   + z(i,j,k-1)   )- &
            .5*( z(i-1,j,k) + z(i-1,j,k-1) )
         dxdzet = .5*( x(i,j,k)   + x(i-1,j,k)   )- &
            .5*( x(i,j,k-1) + x(i-1,j,k-1) )
         dydzet = .5*( y(i,j,k)   + y(i-1,j,k)   )- &
            .5*( y(i,j,k-1) + y(i-1,j,k-1) )
         dzdzet = .5*( z(i,j,k)   + z(i-1,j,k)   )- &
            .5*( z(i,j,k-1) + z(i-1,j,k-1) )
      
         gg113(i,k) = dxdcsi*dxdcsi + dydcsi*dydcsi + dzdcsi*dzdcsi
         gg333(i,k) = dxdzet*dxdzet + dydzet*dydzet + dzdzet*dzdzet
         gg133(i,k) = dxdcsi*dxdzet + dydcsi*dydzet + dzdcsi*dzdzet

         areola4(i,k)=sqrt(gg113(i,k)*gg333(i,k)-gg133(i,k)*gg133(i,k))

      end do
   end do

   do k=kparasta,kparaend
      do i=1,n1
         area4_loc = area4_loc + areola4(i,k)
      end do
   end do

   ! call to know the total area

   if (myid.eq.0) then
      area4 = area4_loc
      sump  = 0.
      ra    = 0.
   end if

   do i=1,nproc-1
 
      if (myid.eq.0) then
         call mpi_recv(ra,1,MPI_REAL_SD,i,i,mpi_comm_world,status,ierr)
         sump = sump + ra

      else if (myid.eq.i) then
         call mpi_send(area4_loc,1,MPI_REAL_SD,0,i,mpi_comm_world,ierr)
      end if

   end do

   if (myid.eq.0) then
      area4 = area4 + sump
      write(*,*)'area face 4 =',area4
   end if

   call mpi_bcast(area4,1,MPI_REAL_SD,0,mpi_comm_world,ierr)

   write(info_run_file,*)'PE=',myid,'area4=',area4


end subroutine
