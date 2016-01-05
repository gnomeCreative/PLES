!***********************************************************************
subroutine flu_turbo(kparasta,kparaend)
   !***********************************************************************
   ! compute explicit terms for turbulence model like
   ! d/dcsi(Uf*annit), periodic version
   !
   use myarrays_velo3
   use myarrays_metri3
                         
   use scala3
   use period
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,kper
   !
   integer ierr,myid,nproc
   integer ncolperproc,kparasta,kparaend,m
   integer kparaendp
   !-----------------------------------------------------------------------
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   !-----------------------------------------------------------------------
   ! compute controvariant flux Uf (jordan notation)
   !
   !     side left and right, periodic
   do k=kparasta,kparaend
      do j=1,jy
         !
         cgra1(0 ,j,k)=csx(0 ,j,k)*.5*(gra1(jx,j,k)+gra1(1,j,k))+ &
            csy(0 ,j,k)*.5*(gra2(jx,j,k)+gra2(1,j,k))+ &
            csz(0 ,j,k)*.5*(gra3(jx,j,k)+gra3(1,j,k))- &
            cgra1(0,j,k)
         cgra1(0 ,j,k)=(1-ip)*cgra1(0 ,j,k)

         cgra1(jx,j,k)=csx(jx,j,k)*.5*(gra1(jx,j,k)+gra1(1,j,k))+ &
            csy(jx,j,k)*.5*(gra2(jx,j,k)+gra2(1,j,k))+ &
            csz(jx,j,k)*.5*(gra3(jx,j,k)+gra3(1,j,k))- &
            cgra1(jx,j,k)
         cgra1(jx,j,k)=(1-ip)*cgra1(jx,j,k)
      !
      end do
   end do
   !
   !     into the field
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx-1
            !
            cgra1(i,j,k)=csx(i,j,k)*.5*(gra1(i+1,j,k)+gra1(i,j,k))+ &
               csy(i,j,k)*.5*(gra2(i+1,j,k)+gra2(i,j,k))+ &
               csz(i,j,k)*.5*(gra3(i+1,j,k)+gra3(i,j,k))- &
               cgra1(i,j,k)
         !
         end do
      end do
   end do
   !
   ! compute controvariant flux Vf (jordan notation)
   !
   !     sides bottom and upper, not periodic
   do k=kparasta,kparaend
      do i=1,jx

         cgra2(i,0,k)=etx(i,0,k)*.5*(gra1(i,jy,k)+gra1(i,1,k))+ &
            ety(i,0,k)*.5*(gra2(i,jy,k)+gra2(i,1,k))+ &
            etz(i,0,k)*.5*(gra3(i,jy,k)+gra3(i,1,k))- &
            cgra2(i,0,k)
         cgra2(i,0,k)=(1-jp)*cgra2(i,0,k)

         cgra2(i,jy,k)=etx(i,jy,k)*.5*(gra1(i,jy,k)+gra1(i,1,k))+ &
            ety(i,jy,k)*.5*(gra2(i,jy,k)+gra2(i,1,k))+ &
            etz(i,jy,k)*.5*(gra3(i,jy,k)+gra3(i,1,k))- &
            cgra2(i,jy,k)
         cgra2(i,jy,k)=(1-jp)*cgra2(i,jy,k)

      end do
   end do

   !     into the field
   do k=kparasta,kparaend
      do j=1,jy-1
         do i=1,jx
            !
            cgra2(i,j,k)=etx(i,j,k)*.5*(gra1(i,j+1,k)+gra1(i,j,k))+ &
               ety(i,j,k)*.5*(gra2(i,j+1,k)+gra2(i,j,k))+ &
               etz(i,j,k)*.5*(gra3(i,j+1,k)+gra3(i,j,k))- &
               cgra2(i,j,k)
         !
         end do
      end do
   end do
   !
   ! compute controvariant flux Wf (jordan notation)
   !
   do j=1,jy
      do i=1,jx
         !
         if (myid.eq.0) then

            cgra3(i,j,0 )=ztx(i,j,0 ) &
               *.5*(gra1(i,j,1)+gra1_appoggio(i,j,jz))+ &
               zty(i,j,0 ) &
               *.5*(gra2(i,j,1)+gra2_appoggio(i,j,jz))+ &
               ztz(i,j,0 ) &
               *.5*(gra3(i,j,1)+gra3_appoggio(i,j,jz))- &
               cgra3(i,j,0)
            cgra3(i,j,0 )=(1-kp)*cgra3(i,j,0 )

         else if (myid.eq.nproc-1) then

            cgra3(i,j,jz)=ztx(i,j,jz) &
               *.5*(gra1_appoggio(i,j,1)+gra1(i,j,jz))+ &
               zty(i,j,jz) &
               *.5*(gra2_appoggio(i,j,1)+gra2(i,j,jz))+ &
               ztz(i,j,jz) &
               *.5*(gra3_appoggio(i,j,1)+gra3(i,j,jz))- &
               cgra3(i,j,jz)
            cgra3(i,j,jz)=(1-kp)*cgra3(i,j,jz)

         endif
      !
      end do
   end do
   !
   !     into the field
   if(myid.eq.nproc-1)then
      kparaendp=kparaend-1
   else
      kparaendp=kparaend
   endif
   do k=kparasta,kparaendp
      do j=1,jy
         do i=1,jx
            !
            cgra3(i,j,k)=ztx(i,j,k)*.5*(gra1(i,j,k+1)+gra1(i,j,k))+ &
               zty(i,j,k)*.5*(gra2(i,j,k+1)+gra2(i,j,k))+ &
               ztz(i,j,k)*.5*(gra3(i,j,k+1)+gra3(i,j,k))- &
               cgra3(i,j,k)
         !
         end do
      end do
   end do
                                        
   !hicco ??? come mai questo vecchio commento senza poi il passaggio tra procs
   ! rendo visibili cgra1 cgra2 cgra3 a tutti i PEs



   return
end
