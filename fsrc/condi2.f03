!***********************************************************************
subroutine condi2
   !***********************************************************************
   ! set boundary conditions on du dv dw
   ! at faces 3 and 4 "sotto" and "sopra"
   ! as in Kim and Moin
   !
   use myarrays_velo3
   use mysending
   use scala3

   implicit none
   !-----------------------------------------------------------------------
   ! array declaration
   integer i,k
   !-----------------------------------------------------------------------
   do k=kparasta,kparaend
      do i=1,jx
         !     face 3 sotto
         delu(i,0,k)= &
            1.875*gra1(i,1,k)-1.25*gra1(i,2,k)+.375*gra1(i,3,k)
         delu(i,0,k)=-delu(i,0,k)
         !
         delv(i,0,k)= &
            1.875*gra2(i,1,k)-1.25*gra2(i,2,k)+.375*gra2(i,3,k)
         delv(i,0,k)=-delv(i,0,k)
         !
         delw(i,0,k)= &
            1.875*gra3(i,1,k)-1.25*gra3(i,2,k)+.375*gra3(i,3,k)
         delw(i,0,k)=-delw(i,0,k)
         !
         !     face 4 sopra
         delu(i,jy+1,k)= &
            .375*gra1(i,jy-2,k)-1.25*gra1(i,jy-1,k)+1.875*gra1(i,jy,k)
         delu(i,jy+1,k)=-delu(i,jy+1,k)
         !
         delv(i,jy+1,k)= &
            .375*gra2(i,jy-2,k)-1.25*gra2(i,jy-1,k)+1.875*gra2(i,jy,k)
         delv(i,jy+1,k)=-delv(i,jy+1,k)
         !
         delw(i,jy+1,k)= &
            .375*gra3(i,jy-2,k)-1.25*gra3(i,jy-1,k)+1.875*gra3(i,jy,k)
         delw(i,jy+1,k)=-delw(i,jy+1,k)
      !
      enddo
   enddo
   !
   return
end
