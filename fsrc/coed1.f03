!***********************************************************************
      subroutine coed1(j,k,del,aa,rh,kparasta,kparaend,isc)
!***********************************************************************
! compute the coefficents for band tridiagonal matrix for the solution
! along csi of the scalar eq.
!
      use myarrays_metri3
      use myarrays_velo3
      use myarrays_density
      use scala3
      use period

      implicit none
!-----------------------------------------------------------------------
!     array declaration
      integer i,j,k,ii,isc
      integer kparasta,kparaend
      real aa(3,*),rh(*)
      real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1) !0:n3+1)
!-----------------------------------------------------------------------
!
      do ii=1,ip
!     computation on left side
      i=1
!
      aa(1,i)=0.
!
      aa(2,i)=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)+ &
              3.*akapt(isc,i-1,j,k)*g11(i-1,j,k)
      aa(2,i)=1.+aa(2,i)*dt/giac(i,j,k)/2.
!
      aa(3,i)=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))*g11(i,j,k)+ &
              akapt(isc,i-1,j,k)*g11(i-1,j,k)/3.
      aa(3,i)=-dt*aa(3,i)/giac(i,j,k)/2.      
!
!
      rh(i)=rhs(i,j,k)
!
!     computation on right side
      i=jx
!
      aa(1,i)=akapt(isc,i+1,j,k)*g11(i,j,k)/3.+ &
              .5*(akapt(isc,i,j,k)+akapt(isc,i-1,j,k))*g11(i-1,j,k)
      aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
      !
      aa(2,i)=akapt(isc,i+1,j,k)*g11(i,j,k)*3.+ &
         .5*(akapt(isc,i,j,k)+akapt(isc,i-1,j,k))*g11(i-1,j,k)
      aa(2,i)=1.+dt*aa(2,i)/giac(i,j,k)/2.
!
      aa(3,i)=0.      
!
      rh(i)=rhs(i,j,k)
!
      enddo
!
!     computation into the field
      do i=1+ip,jx-ip
!
      aa(1,i)=.5*(akapt(isc,i,j,k)+akapt(isc,i-1,j,k))*g11(i-1,j,k)
      aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
!
      aa(3,i)=.5*(akapt(isc,i+1,j,k)+akapt(isc,i,j,k))*g11(i,j,k)
      aa(3,i)=-dt*aa(3,i)/giac(i,j,k)/2.
!
      aa(2,i)=1.-aa(1,i)-aa(3,i)     
!
      rh(i)=rhs(i,j,k)
!
      enddo
!
      return
      end
