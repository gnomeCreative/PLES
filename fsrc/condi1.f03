!***********************************************************************
subroutine condi1
   !************************************************************************
      ! set boundary conditions on du dv dw
      ! at faces 1 and 2 "sinistra" and "destra"
      ! as in Kim and Moin
      !
      use myarrays_velo3, only: delu, delv, delw, gra1, gra2, gra3
      use mysending, only: kparaend, kparasta
      use scala3, only: jx, jy

      implicit none
      !-----------------------------------------------------------------------
      ! array declaration
      integer j,k
      !-----------------------------------------------------------------------
      !
      do k=kparasta,kparaend
         do j=1,jy
            !     face 1 sinistra
            delu(0,j,k)= &
               1.875*gra1(1,j,k)-1.25*gra1(2,j,k)+.375*gra1(3,j,k)
            delu(0,j,k)=-delu(0,j,k)
            !
            delv(0,j,k)= &
               1.875*gra2(1,j,k)-1.25*gra2(2,j,k)+.375*gra2(3,j,k)
            delv(0,j,k)=-delv(0,j,k)
            !
            delw(0,j,k)= &
               1.875*gra3(1,j,k)-1.25*gra3(2,j,k)+.375*gra3(3,j,k)
            delw(0,j,k)=-delw(0,j,k)
            !
            !     face 2 destra
            delu(jx+1,j,k)= &
               .375*gra1(jx-2,j,k)-1.25*gra1(jx-1,j,k)+1.875*gra1(jx,j,k)
            delu(jx+1,j,k)=-delu(jx+1,j,k)
            !
            delv(jx+1,j,k)= &
               .375*gra2(jx-2,j,k)-1.25*gra2(jx-1,j,k)+1.875*gra2(jx,j,k)
            delv(jx+1,j,k)=-delv(jx+1,j,k)
            !
            delw(jx+1,j,k)= &
               .375*gra3(jx-2,j,k)-1.25*gra3(jx-1,j,k)+1.875*gra3(jx,j,k)
            delw(jx+1,j,k)=-delw(jx+1,j,k)
         !
         enddo
      enddo
      !
      return

   end subroutine condi1
