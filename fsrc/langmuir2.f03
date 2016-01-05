!***********************************************************************
      subroutine langmuir2(myid,nproc,kparasta,kparaend)
!***********************************************************************
! for langmuir circulation, compute the controvariant term
      use myarrays_WB
      use myarrays_LC
      use myarrays_velo3
      use myarrays_metri3
      !
      use scala3
      use period
      !
      use mpi

      implicit none

!-----------------------------------------------------------------------
!     array declaration
      integer i,j,k
      integer myid,nproc,kparasta,kparaend
      integer kparastal,kparaendl

      real uinter, winter
!-----------------------------------------------------------------------
! sides right and left for ucs
!
      do k=kparasta,kparaend
      do j=1,jy

      ucs(0,j,k)=u_drift(0,j,k)*csx(0,j,k)+ &
                 w_drift(0,j,k)*csz(0,j,k)

      ucs(jx,j,k)=u_drift(jx+1,j,k)*csx(jx,j,k)+ &
                  w_drift(jx+1,j,k)*csz(jx,j,k)

      end do
      end do
!
! into the field
!
      do k=kparasta,kparaend
      do j=1,jy
      do i=ip,jx-ip

      uinter=.5*(u_drift(i,j,k)+u_drift(i+1,j,k))
      winter=.5*(w_drift(i,j,k)+w_drift(i+1,j,k))

      ucs(i,j,k)=uinter*csx(i,j,k)+ &
                 winter*csz(i,j,k)

      end do
      end do
      end do
!
!-----------------------------------------------------------------------
!
! sides back and front for wcs
!
      if (myid.eq.0) then
      do i=1,jx
      do j=1,jy

      wcs(i,j,0)=u_drift(i,j,0)*ztx(i,j,0)+ &
                 w_drift(i,j,0)*ztz(i,j,0)

      end do
      end do
      end if

      if (myid.eq.nproc-1) then
      do i=1,jx
      do j=1,jy

      wcs(i,j,jz)=u_drift(i,j,jz+1)*ztx(i,j,jz)+ &
                  w_drift(i,j,jz+1)*ztz(i,j,jz)
      end do
      end do
      end if

!
! into the field
!
      if (myid.eq.0) then
      kparastal=kp
      kparaendl=kparaend
      else if (myid.eq.nproc-1) then
      kparastal=kparasta
      kparaendl=kparaend-kp
      else
      kparastal=kparasta
      kparaendl=kparaend
      end if

      do k=kparastal,kparaendl
      do j=1,jy
      do i=1,jx

      uinter=.5*(u_drift(i,j,k)+u_drift(i,j,k+1))
      winter=.5*(w_drift(i,j,k)+w_drift(i,j,k+1))

      wcs(i,j,k)=uinter*ztx(i,j,k)+ &
                 winter*ztz(i,j,k)

      end do
      end do
      end do
!
      end subroutine
