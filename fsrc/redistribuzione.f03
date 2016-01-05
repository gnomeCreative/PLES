!***********************************************************************
subroutine redistribuzione(bodyforce,area1,area2,area5,area6,kparasta,kparaend,myid,nproc)
   !***********************************************************************
   ! mass distribution on boundary to obtain a flow field with
   ! zero divergence. Used for nesting
   !
   use myarrays_nesting
   use myarrays_velo3
   use myarrays2
   !
   use scala3
   use tipologia
   use orl
   !
   use mpi


   implicit none
   !
   !-----------------------------------------------------------------------
   !     variables declaration
   integer i,j,k
   integer bodyforce
   integer kparasta,kparaend,myid,nproc
   integer ierr ,status(mpi_status_size)
   real area1,area2,area5,area6
   real areatot
   real massa1,massa2,massa3,massa4,massa5,massa6
   real massatot
   real massa1tot,massa2tot
   real massa3tot,massa4tot
   real massa5tot,massa6tot
   real bilancio

   real ar1,ar2,ar3,ar4,ar5,ar6
   real artot
   real ar1tot,ar2tot
   real ar3tot,ar4tot
   real ar5tot,ar6tot
   !
   !-----------------------------------------------------------------------
   ! put to zero fluxes at faces 3 and 4
   !
   vc(:, 0,kparasta:kparaend)=0.
   vc(:,jy,kparasta:kparaend)=0.
	 
   !-----------------------------------------------------------------------
   ! mass check
   !     faces 1 and 2
   ar1 = 0.
   massa1=0.
   do k=kparasta,kparaend
      do j=1,jy
         massa1=massa1+uc(0,j,k)*index_out1(j,k)
         if(uc(0,j,k).lt.0.)ar1 = ar1 + areola1(j,k)*index_out1(j,k)
      end do
   end do
      
   ar2 = 0.
   massa2=0.
   do k=kparasta,kparaend
      do j=1,jy
         massa2=massa2-uc(jx,j,k)*index_out2(j,k)
         if(uc(jx,j,k).gt.0.)ar2 = ar2 + areola2(j,k)*index_out2(j,k)
      end do
   end do

   !     solid walls faces 3 and 4
   massa3=0.
   massa4=0.

   !     faces 5 and 6

   ar5 = 0.
   massa5=0.
   if(myid.eq.0)then
      do j=1,jy
         do i=1,jx
            massa5=massa5+wc(i,j,0)*index_out5(i,j)
            if(wc(i,j,0).lt.0.)ar5 = ar5 + areola5(i,j)*index_out5(i,j)
         end do
      end do
   end if
      
   ar6 = 0.
   massa6=0.
   if(myid.eq.nproc-1)then
      do j=1,jy
         do i=1,jx
            massa6=massa6-wc(i,j,jz)*index_out6(i,j)
            if(wc(i,j,jz).gt.0.)ar6 = ar6 + areola6(i,j)*index_out6(i,j)
         end do
      end do
   end if
      
   !     global mass at the faces
   call MPI_ALLREDUCE(massa1,massa1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa2,massa2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa3,massa3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa4,massa4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa5,massa5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa6,massa6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

   !     global mass at the faces
   call MPI_ALLREDUCE(ar1,ar1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ar2,ar2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ar3,ar3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ar4,ar4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ar5,ar5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ar6,ar6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

   if(myid.eq.0)then
      bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot
     
      artot = ar1tot + ar2tot + ar3tot + ar4tot + ar5tot + ar6tot
      write(*,*)'mass balance before redistribution: ',bilancio
      write(*,*)massa1tot,massa2tot,massa3tot,massa4tot,massa5tot,massa6tot
   end if
      
   call MPI_BCAST(bilancio,1,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(artot,1,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)
   !-----------------------------------------------------------------------
   ! excluding ibm and inflow
   areatot = 0.
      
   if(bodyforce == 0)then
      if(infout1 /= 0 .and. massa1tot<0.)areatot = areatot + area1  !goes out
      if(infout2 /= 0 .and. massa2tot<0.)areatot = areatot + area2
      if(infout5 /= 0 .and. massa5tot<0.)areatot = areatot + area5
      if(infout6 /= 0 .and. massa6tot<0.)areatot = areatot + area6
   else
      if(infout1/=0 .and. massa1tot<0.)areatot=areatot+area_bagnata1
      if(infout2/=0 .and. massa2tot<0.)areatot=areatot+area_bagnata2
      if(infout5/=0 .and. massa5tot<0.)areatot=areatot+area_bagnata5
      if(infout6/=0 .and. massa6tot<0.)areatot=areatot+area_bagnata6
   end if
      


   !-----------------------------------------------------------------------
   ! fluxes adjustment with mass defeact

   !     face 1
   if(infout1 /= 0)then
      do k=kparasta,kparaend
         do j=1,jy
            if(uc(0,j,k).lt.0.)then
               uc(0,j,k)=uc(0,j,k)-bilancio*areola1(j,k)/artot
               uc(0,j,k)=uc(0,j,k)*index_out1(j,k)
               ucp1(j,k)=uc(0,j,k)
            end if
         end do
      end do
   end if
   !     face 2
   if(infout2 /= 0)then
      do k=kparasta,kparaend
         do j=1,jy
            if(uc(jx,j,k).gt.0.)then
               uc(jx,j,k)=uc(jx,j,k)+bilancio*areola2(j,k)/artot
               uc(jx,j,k)=uc(jx,j,k)*index_out2(j,k)
               ucp2(j,k)=uc(jx,j,k)
            end if
         end do
      end do
   end if
   !     face 5
   if(infout5 /= 0)then
      if(myid.eq.0)then
         do j=1,jy
            do i=1,jx
               if(wc(i,j,0).lt.0.)then
                  wc(i,j,0)=wc(i,j,0)-bilancio*areola5(i,j)/artot
                  wc(i,j,0)=wc(i,j,0)*index_out5(i,j)
                  wcp5(i,j)=wc(i,j,0)
               end if
            end do
         end do
      end if
   end if
   !     face 6
   if(infout6 /= 0)then
      if(myid.eq.nproc-1)then
         do j=1,jy
            do i=1,jx
               if(wc(i,j,jz).gt.0.)then
                  wc(i,j,jz)=wc(i,j,jz)+bilancio*areola6(i,j)/artot
                  wc(i,j,jz)=wc(i,j,jz)*index_out6(i,j)
                  wcp6(i,j)=wc(i,j,jz)
               end if
            end do
         end do
      end if
   end if
   !
   !-----------------------------------------------------------------------
   ! mass check after redistribution
   !     faces 1 and 2
   massa1=0.
   do k=kparasta,kparaend
      do j=1,jy
         massa1=massa1+uc(0,j,k)
      end do
   end do

   massa2=0.
   do k=kparasta,kparaend
      do j=1,jy
         massa2=massa2-uc(jx,j,k)
      end do
   end do

   !     solid walls faces 3 and 4
   massa3=0.
   massa4=0.

   !     faces 5 and 6
   massa5=0.
   if(myid.eq.0)then
      do j=1,jy
         do i=1,jx
            massa5=massa5+wc(i,j,0)
         end do
      end do
   end if

   massa6=0.
   if(myid.eq.nproc-1)then
      do j=1,jy
         do i=1,jx
            massa6=massa6-wc(i,j,jz)
         end do
      end do
   end if

   !     global mass at the faces
   call MPI_ALLREDUCE(massa1,massa1tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa2,massa2tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa3,massa3tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa4,massa4tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa5,massa5tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(massa6,massa6tot,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)

   if(myid.eq.0)then
      bilancio = massa1tot + massa2tot + massa3tot + massa4tot + massa5tot + massa6tot
      write(*,*)'mass balance after redistribution: ',bilancio
   end if

   ! chicco     NOTA andrebbero forse cambiate le velocita' con i
   ! chicco     coseni direttori in seguito alla ridistribuzione massa

   return
end
