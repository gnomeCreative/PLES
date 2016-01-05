!***********************************************************************
subroutine cellep(xx1,yy1,zz1,xx2,yy2,zz2,kgridparasta)
   !***********************************************************************
   ! periodic points out of the domain
   !
   use myarrays_metri3
   use mysending
   use scala3
   use period

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer kgridparasta
   integer kgpsta,kgpend

   real :: xx1(-8:n1+8,-8:n2+8,n3-8:n3)
   real :: yy1(-8:n1+8,-8:n2+8,n3-8:n3)
   real :: zz1(-8:n1+8,-8:n2+8,n3-8:n3)

   real :: xx2(-8:n1+8,-8:n2+8,0:8)
   real :: yy2(-8:n1+8,-8:n2+8,0:8)
   real :: zz2(-8:n1+8,-8:n2+8,0:8)
   !
   !-----------------------------------------------------------------------
   ! csi direction, adiacent blocks to faces sn and dx (1 and 2)
   !
   do i=1,8
      do j=0,jy
         do k=kparasta-1,kparaend !+1 !0,jz
            x(-i,j,k)=x(0,j,k)-(x(jx,j,k)-x(jx-i,j,k))
            y(-i,j,k)=y(jx-i,j,k)
            z(-i,j,k)=z(jx-i,j,k)
            !
            x(jx+i,j,k)=x(jx,j,k)-(x(0,j,k)-x(i,j,k))
            y(jx+i,j,k)=y(i,j,k)
            z(jx+i,j,k)=z(i,j,k)
         end do
      end do
   end do

   if(myid.eq.0)then
      do i=1,8
         do j=0,jy
            do k=n3-8,n3
               xx1(-i,j,k)=xx1(0,j,k)-(xx1(jx,j,k)-xx1(jx-i,j,k))
               yy1(-i,j,k)=yy1(jx-i,j,k)
               zz1(-i,j,k)=zz1(jx-i,j,k)

               xx1(jx+i,j,k)=xx1(jx,j,k)-(xx1(0,j,k)-xx1(i,j,k))
               yy1(jx+i,j,k)=yy1(i,j,k)
               zz1(jx+i,j,k)=zz1(i,j,k)
            end do
         end do
      end do
   end if

   if(myid.eq.nproc-1)then
      do i=1,8
         do j=0,jy
            do k=0,8
               xx2(-i,j,k)=xx2(0,j,k)-(xx2(jx,j,k)-xx2(jx-i,j,k))
               yy2(-i,j,k)=yy2(jx-i,j,k)
               zz2(-i,j,k)=zz2(jx-i,j,k)

               xx2(jx+i,j,k)=xx2(jx,j,k)-(xx2(0,j,k)-xx2(i,j,k))
               yy2(jx+i,j,k)=yy2(i,j,k)
               zz2(jx+i,j,k)=zz2(i,j,k)
            end do
         end do
      end do
   end if
   !-----------------------------------------------------------------------
   ! zita direction, blocks adiacent to faces in and av (5 and 6)
   ! fill the corner
   if(myid.eq.0)then
      do k=1,8
         do i=-8,jx+8
            do j=0,jy
               x(i,j,-k)=xx1(i,j,jz-k)
               y(i,j,-k)=yy1(i,j,jz-k)
               z(i,j,-k)=z(i,j,0)-(zz1(i,j,jz)-zz1(i,j,jz-k))
            end do
         end do
      end do
   end if

   if(myid.eq.nproc-1)then
      do k=1,8
         do i=-8,jx+8
            do j=0,jy
               x(i,j,jz+k)=xx2(i,j,k)
               y(i,j,jz+k)=yy2(i,j,k)
               z(i,j,jz+k)=z(i,j,jz)-(zz2(i,j,0)-zz2(i,j,k))
            end do
         end do
      end do
   end if
   !-----------------------------------------------------------------------
   ! eta direction, blocks adiacent to st and sp (3 and 4)
   !
   !     k grid parallel start -> kgpsta
   !     k grid parallel end   -> kgpend
   if(myid.eq.0)then
      kgpsta=kgridparasta-8
      kgpend=kparaend !+1
   elseif(myid.eq.nproc-1)then
      kgpsta=kparasta-1
      kgpend=kparaend+8
   else
      kgpsta=kparasta-1
      kgpend=kparaend  !-1
   end if


   do j=1,8
      do i=-8,jx+8
         do k=kgpsta,kgpend  !-8,jz+8
            x(i,-j,k)=x(i,jy-j,k)
            y(i,-j,k)=y(i,0,k)-(y(i,jy,k)-y(i,jy-j,k))
            z(i,-j,k)=z(i,jy-j,k)
            !
            x(i,jy+j,k)=x(i,j,k)
            y(i,jy+j,k)=y(i,jy,k)-(y(i,0,k)-y(i,j,k))
            z(i,jy+j,k)=z(i,j,k)
         end do
      end do
   end do
   !
   ! blocco di spigolo su angolo csi=0 , zita=0     e angolo csi=jx zita=jz
   !
   !      do i=1,8
   !      do k=1,8
   !      do j=0,jy
   !      x(-i,j,-k)=x(0,j,0)-(x(jx,j,jz)-x(jx-i,j,jz-k))
   !      y(-i,j,-k)=y(jx-i,j,jz-k)
   !      z(-i,j,-k)=z(0,j,0)-(z(jx,j,jz)-z(jx-i,j,jz-k))
   !c
   !      x(jx+i,j,jz+k)=x(jx,j,jz)-(x(0,j,0)-x(i,j,k))
   !      y(jx+i,j,jz+k)=y(i,j,k)
   !      z(jx+i,j,jz+k)=z(jx,j,jz)-(z(0,j,0)-z(i,j,k))
   !      end do
   !      end do
   !
   ! blocco di spigolo su angolo csi=0 zita=jz e angolo csi=jx zita=0
   !
   !      do i=1,8
   !      do k=1,8
   !      do j=0,jy
   !      x(-i,j,jz+k)=x(0,j,jz)-(x(jx,j,0)-x(jx-i,j,k))
   !      y(-i,j,jz+k)=y(jx-i,j,k)
   !      z(-i,j,jz+k)=z(0,j,jz)-(z(jx,j,0)-z(jx-i,j,k))
   !c
   !      x(jx+i,j,-k)=x(jx,j,0)-(x(0,j,jz)-x(i,j,jz-k))
   !      y(jx+i,j,-k)=y(i,j,jz-k)
   !      z(jx+i,j,-k)=z(jx,j,0)-(z(0,j,jz)-z(i,j,jz-k))
   return
end
