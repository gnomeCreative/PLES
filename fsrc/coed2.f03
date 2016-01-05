!***********************************************************************
      subroutine coed2(i,k,del,aa,rh,kparasta,kparaend,isc)
!***********************************************************************
! compute the coefficents for band tridiagonal matrix for the solution
! along eta of the scalar eq.
!
      use myarrays_metri3
      use myarrays_density
      use scala3
      use period

      implicit none
!-----------------------------------------------------------------------
!     array declaration
      integer i,j,k,jj,isc
      integer kparasta,kparaend
      real aa(3,*),rh(*)
      real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !0:n3+1)
!-----------------------------------------------------------------------
!
!
      do jj=1,jp
!     computation on bottom side
      j=1
!
      aa(1,j)=0.
!
      aa(2,j)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)+ &
              3.*akaptV(isc,i,j-1,k)*g22(i,j-1,k)
      aa(2,j)=1.+aa(2,j)*dt/giac(i,j,k)/2.
!
      aa(3,j)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)+ &
              akaptV(isc,i,j-1,k)*g22(i,j-1,k)/3.
      aa(3,j)=-dt*aa(3,j)/giac(i,j,k)/2.
!
      rh(j)=del(i,j,k)
!
!     computation on upper side
      j=jy
!
      aa(1,j)=akaptV(isc,i,j+1,k)*g22(i,j,k)/3.+ &
           .5*(akaptV(isc,i,j,k)+akaptV(isc,i,j-1,k))*g22(i,j-1,k)
      aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
!
      aa(2,j)=akaptV(isc,i,j+1,k)*g22(i,j,k)*3.+ &
           .5*(akaptV(isc,i,j,k)+akaptV(isc,i,j-1,k))*g22(i,j-1,k)
      aa(2,j)=1.+dt*aa(2,j)/giac(i,j,k)/2.
!
      aa(3,j)=0.     
      rh(j)=del(i,j,k)
!
      enddo
!
!     computation into the field
      do j=1+jp,jy-jp
         !
         aa(1,j)=.5*(akaptV(isc,i,j,k) &
            +akaptV(isc,i,j-1,k))*g22(i,j-1,k)
      aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
!
      aa(3,j)=.5*(akaptV(isc,i,j,k)+akaptV(isc,i,j+1,k))*g22(i,j,k)
      aa(3,j)=-dt*aa(3,j)/giac(i,j,k)/2.
!
      aa(2,j)=1.-aa(1,j)-aa(3,j)
!
      rh(j)=del(i,j,k)
!
      enddo
!
      return
      end
