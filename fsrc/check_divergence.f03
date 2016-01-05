!***********************************************************************
subroutine check_divergence(tipo)
   !***********************************************************************
   use mysending, only: deepl, deepr, kparaend, kparasta, myid, nproc
   use mysettings, only: lett
   use myarrays_velo3, only: uc, uc1_orl, uc2_orl, vc, vc3_orl, vc4_orl, wc, wc5_orl, wc6_orl
   use myarrays_metri3, only: giac
   use myarrays_ibm, only: bodyforce
   !
   use scala3, only: jx, jy, jz, n1, n2
   use tipologia
   use orl, only: infout1, infout2, infout3, infout4, infout5, infout6
   !
   use mpi

   implicit none

       
   !-----------------------------------------------------------------------
   !     array declarations
   integer i,j,k
   integer ierr,status(MPI_STATUS_SIZE)
   integer,intent(in) :: tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   !integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

   real div,divg
   real adiv,adivg
   real divmax,divgmax
   real divmax_loc,divgmax_loc

   real, allocatable :: uc_prov(:,:,:)
   real, allocatable :: vc_prov(:,:,:)
   real, allocatable :: wc_prov(:,:,:)

   !-----------------------------------------------------------------------
   !      call MPI_SENDRECV(wc(1,1,kparaend),jx*jy,
   !     >  		MPI_REAL_SD,rightpem,51+myid,
   !     >  		wc(1,1,kparasta-1),jx*jy,
   !     >  		MPI_REAL_SD,leftpem,50+myid,
   !     >  		MPI_COMM_WORLD,status,ierr)

   allocate(uc_prov(0:n1,1:n2,kparasta  :kparaend))
   allocate(vc_prov(1:n1,0:n2,kparasta  :kparaend))
   allocate(wc_prov(1:n1,1:n2,kparasta-1:kparaend))

    divmax=0.
    divgmax=0.

   do k=kparasta,kparaend
      do j=1,jy
         do i=0,jx
            uc_prov(i,j,k) = uc(i,j,k)
         end do
      end do
   end do
      
   do k=kparasta,kparaend
      do j=0,jy
         do i=1,jx
            vc_prov(i,j,k) = vc(i,j,k)
         end do
      end do
   end do
      
   do k=kparasta-1,kparaend
      do j=1,jy
         do i=1,jx
            wc_prov(i,j,k) = wc(i,j,k)
         end do
      end do
   end do
      
      
   if(lett==1)then
      
      if(infout1 == 0)then
         do k=kparasta,kparaend
            do j=1,jy
               uc_prov(0,j,k) = uc1_orl(j,k)
            end do
         end do
      end if

      if(infout2 == 0)then
         do k=kparasta,kparaend
            do j=1,jy
               uc_prov(jx,j,k) = uc2_orl(j,k)
            end do
         end do
      end if


      if(infout3==0)then      
         do k=kparasta,kparaend
            do i=1,jx
               vc_prov(i,0,k) = vc3_orl(i,k)
            end do
         end do
      end if

      if(infout4==0)then      
         do k=kparasta,kparaend
            do i=1,jx
               vc_prov(i,jy,k) = vc4_orl(i,k)
            end do
         end do
      end if
     
      if(infout5==0 .and. myid==0)then
         do j=1,jy
            do i=1,jx
               wc_prov(i,j,0) = wc5_orl(i,j)
            end do
         end do
      end if

      if(infout6==0 .and. myid==nproc-1)then
         do j=1,jy
            do i=1,jx
               wc_prov(i,j,jz) = wc6_orl(i,j)
            end do
         end do
      end if
      
   end if ! on lett
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     without IBM bodyforce=0

   if(bodyforce==0)then
      divmax_loc=0.
      divgmax_loc=0.
      
      do k=kparasta,kparaend
         do j=2,jy-1
            do i=2,jx-1 !1,jx
               !
               if(k==1 .or. k==jz)cycle
               div=uc_prov(i,j,k)-uc_prov(i-1,j,k)+vc_prov(i,j,k)-vc_prov(i,j-1,k)+wc_prov(i,j,k)-wc_prov(i,j,k-1)
     
               !
               divg=div/giac(i,j,k)
               !
               adiv=abs(div)
               adivg=abs(divg)
               !
               divmax_loc=max(divmax_loc,adiv)
               divgmax_loc=max(divgmax_loc,adivg)
            end do
         end do
      end do
      
      !     find the max between all procs
      call MPI_REDUCE(divmax_loc,divmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(divgmax_loc,divgmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      !     all procs know divmax
      !
      if (myid==0) then
         write(*,*)myid,'divmax  inside field',divmax
         write(*,*)myid,'divgmax inside field',divgmax
      endif      
      !.......................................................................
      divmax_loc=0.
      divgmax_loc=0.
      
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               div=uc_prov(i,j,k)-uc_prov(i-1,j,k)+vc_prov(i,j,k)-vc_prov(i,j-1,k)+wc_prov(i,j,k)-wc_prov(i,j,k-1)
     
               !
               divg=div/giac(i,j,k)
               !
               adiv=abs(div)
               adivg=abs(divg)
               !
               divmax_loc=max(divmax_loc,adiv)
               divgmax_loc=max(divgmax_loc,adivg)
            end do
         end do
      end do
      !     find the max between all procs
      !
      call MPI_REDUCE(divmax_loc,divmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(divgmax_loc,divgmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !


      !
      !     all procs know divmax
      !

      if (myid==0) then
         write(*,*)myid,'divmax  all',divmax
         write(*,*)myid,'divgmax all',divgmax
      endif         

   end if
   !-----------------------------------------------------------------------

   !     with IBM bodyforce=1
   if(bodyforce==1)then

      divmax_loc=0.
      divgmax_loc=0.

      do k=kparasta,kparaend
         do j=2,jy-1
            do i=2,jx-1

               if(k==1 .or. k==jz)cycle

               !       only for fluid cells tipo=2
               if(tipo(i,j,k)==2)then
                  div=uc_prov(i,j,k)-uc_prov(i-1,j,k)+vc_prov(i,j,k)-vc_prov(i,j-1,k)+wc_prov(i,j,k)-wc_prov(i,j,k-1)
                  !
                  divg=div/giac(i,j,k)
                  !
                  adiv=abs(div)
                  adivg=abs(divg)
                  !
                  divmax_loc=max(divmax_loc,adiv)
                  divgmax_loc=max(divgmax_loc,adivg)
               end if
            end do
         end do
      end do
      
      !     find the max between all procs
      call MPI_REDUCE(divmax_loc,divmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(divgmax_loc,divgmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      !     all procs know divmax
      !
      if (myid==0) then
         write(*,*)myid,'divmax  inside field fluid',divmax
         write(*,*)myid,'divgmax inside field fluid',divgmax
      endif        
      
      !.......................................................................
      divmax_loc=0.
      divgmax_loc=0.

      do k=kparasta,kparaend
         do j=2,jy-1
            do i=2,jx-1

               if(k==1 .or. k==jz)cycle

               !       only for fluid cells tipo=2
               if(tipo(i,j,k)==1)then
                  div=uc_prov(i,j,k)-uc_prov(i-1,j,k)+vc_prov(i,j,k)-vc_prov(i,j-1,k)+wc_prov(i,j,k)-wc_prov(i,j,k-1)
                  !
                  divg=div/giac(i,j,k)
                  !
                  adiv=abs(div)
                  adivg=abs(divg)
                  !
                  divmax_loc=max(divmax_loc,adiv)
                  divgmax_loc=max(divgmax_loc,adivg)
               end if
            end do
         end do
      end do
      
      !     find the max between all procs
      call MPI_REDUCE(divmax_loc,divmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(divgmax_loc,divgmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      !     all procs know divmax
      !
      if (myid==0) then
         write(*,*)myid,'divmax  on ib',divmax
         write(*,*)myid,'divgmax on ib',divgmax
      endif       
      
      !.......................................................................
      divmax_loc=0.
      divgmax_loc=0.

      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx

               div=uc_prov(i,j,k)-uc_prov(i-1,j,k)+vc_prov(i,j,k)-vc_prov(i,j-1,k)+wc_prov(i,j,k)-wc_prov(i,j,k-1)
               !
               divg=div/giac(i,j,k)
               !
               adiv=abs(div)
               adivg=abs(divg)
               !
               divmax_loc=max(divmax_loc,adiv)
               divgmax_loc=max(divgmax_loc,adivg)

            end do
         end do
      end do
      
      !     find the max between all procs
      call MPI_REDUCE(divmax_loc,divmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(divgmax_loc,divgmax,1,MPI_REAL_SD,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      !
      !     all procs know divmax
      !
      if (myid==0) then
         write(*,*)myid,'divmax  all',divmax
         write(*,*)myid,'divgmax all',divgmax
      endif       
      
   end if
   !-----------------------------------------------------------------------
   !

   !-----------------------------------------------------------------------
   deallocate(uc_prov)
   deallocate(vc_prov)
   deallocate(wc_prov)



end
