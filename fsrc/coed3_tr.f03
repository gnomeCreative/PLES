!***********************************************************************
subroutine coed3_tr(j,k,aka_tr,del,aa,rh,myid,g33_tr,giac_tr, &
   iparasta,iparaend )
!***********************************************************************
! compute the coefficent for the banded tridiagonal matrix
! for solution in zita
      use scala3

      implicit none
!
!-----------------------------------------------------------------------
!     array declaration
      integer i,j,k,kk,myid
      integer iparasta,iparaend
      real aa(3,*),rh(*)
      real aka_tr(0:n3+1,0:n2+1,iparasta-1:iparaend+1)
      real del(0:n3+1,0:n2+1,iparasta:iparaend)
      real g33_tr(0:n3,n2,iparasta:iparaend) !n1)
      real giac_tr(n3,n2,iparasta:iparaend) !n1)
!-----------------------------------------------------------------------
!
!     inside the field
!
      do i=1,jz
!
      aa(1,i)=.5*(aka_tr(i,j,k)+aka_tr(i-1,j,k))*g33_tr(i-1,j,k)
      aa(1,i)=-dt*aa(1,i)/giac_tr(i,j,k)/2.
!
      aa(3,i)=.5*(aka_tr(i,j,k)+aka_tr(i+1,j,k))*g33_tr(i,j,k)
      aa(3,i)=-dt*aa(3,i)/giac_tr(i,j,k)/2.
!
      aa(2,i)=1.-aa(1,i)-aa(3,i)
!
      rh(i)=del(i,j,k)
!
      enddo
!
      return
      end
