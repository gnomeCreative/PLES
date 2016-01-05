!***********************************************************************
subroutine condi3
   !***********************************************************************
   ! set boundary conditions on du dv dw
   ! at faces 5 and 6 "avanti" and "indietro"
   ! as in Kim and Moin
   !
   use myarrays_velo3
   use mysending
   use scala3

   implicit none
   !-----------------------------------------------------------------------
   ! array declaration
   integer ierr
   integer i,j
   !-----------------------------------------------------------------------
   !
   !     face 5 avanti
   if (myid.eq.0) then
      !
      do j=1,jy
         do i=1,jx
            !
            delu(i,j,0)= &
               1.875*gra1(i,j,1)-1.25*gra1(i,j,2)+.375*gra1(i,j,3)
            delu(i,j,0)=-delu(i,j,0)
            !
            delv(i,j,0)= &
               1.875*gra2(i,j,1)-1.25*gra2(i,j,2)+.375*gra2(i,j,3)
            delv(i,j,0)=-delv(i,j,0)
            !
            delw(i,j,0)= &
               1.875*gra3(i,j,1)-1.25*gra3(i,j,2)+.375*gra3(i,j,3)
            delw(i,j,0)=-delw(i,j,0)
         !
         enddo
      enddo

   endif
   !
   !     face 6 indietro
   if (myid.eq.nproc-1) then
      !
      do j=1,jy
         do i=1,jx
            !
            delu(i,j,jz+1)= &
               .375*gra1(i,j,jz-2)-1.25*gra1(i,j,jz-1)+1.875*gra1(i,j,jz)
            delu(i,j,jz+1)=-delu(i,j,jz+1)
            !
            delv(i,j,jz+1)= &
               .375*gra2(i,j,jz-2)-1.25*gra2(i,j,jz-1)+1.875*gra2(i,j,jz)
            delv(i,j,jz+1)=-delv(i,j,jz+1)
            !
            delw(i,j,jz+1)= &
               .375*gra3(i,j,jz-2)-1.25*gra3(i,j,jz-1)+1.875*gra3(i,j,jz)
            delw(i,j,jz+1)=-delw(i,j,jz+1)
         !
         enddo
      enddo
   !
   endif

   return
end
