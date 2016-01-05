!***********************************************************************
      subroutine coef2_par(i,k,del,aa,rh,kparasta,kparaend)
!***********************************************************************
! compute the coefficent of the band tridiagonal matrix for solution
! along eta for equations in u,v,w
!
      use myarrays_metri3
      use myarrays_velo3
      use scala3
      use period

      implicit none
!-----------------------------------------------------------------------
!     array declaration
      integer kparasta,kparaend
      integer i,j,k,jj
      real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1),aa(3,*),rh(*)
!-----------------------------------------------------------------------
!
      do jj=1,jp
!
!     side bottom
      j=1
!
      aa(1,j)=0.
!
      aa(2,j)=.5*(annitV(i,j,k)+annitV(i,j+1,k))*g22(i,j,k) &
             + 3.*annitV(i,j-1,k)*g22(i,j-1,k)
      aa(2,j)=1.+aa(2,j)*dt/giac(i,j,k)/2.
!
      aa(3,j)=.5*(annitV(i,j,k)+annitV(i,j+1,k))*g22(i,j,k) &
             + annitV(i,j-1,k)*g22(i,j-1,k)/3.
      aa(3,j)=-dt*aa(3,j)/giac(i,j,k)/2.
!
      rh(j)=del(i,j,k)+ &
       annitV(i,j-1,k)*g22(i,j-1,k)* &
       del(i,j-1,k)*8.*dt/giac(i,j,k)/6.
!
!     side upper
      j=jy
!
      aa(1,j)=annitV(i,j+1,k)*g22(i,j,k)/3.+ &
            .5*(annitV(i,j,k)+annitV(i,j-1,k))*g22(i,j-1,k)
      aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
!
      aa(2,j)=annitV(i,j+1,k)*g22(i,j,k)*3.+ &
            .5*(annitV(i,j,k)+annitV(i,j-1,k))*g22(i,j-1,k)
      aa(2,j)=1.+dt*aa(2,j)/giac(i,j,k)/2.
!
      aa(3,j)=0.
!
      rh(j)=del(i,j,k)+ &
            annitV(i,j+1,k)*g22(i,j,k) &
           *del(i,j+1,k)*8.*dt/giac(i,j,k)/6.
!
      enddo
!
!     into the field
      do j=1+jp,jy-jp
!
      aa(1,j)=.5*(annitV(i,j,k) &
             +annitV(i,j-1,k))*g22(i,j-1,k)
      aa(1,j)=-dt*aa(1,j)/giac(i,j,k)/2.
      !
      aa(3,j)=.5*(annitV(i,j,k) &
         +annitV(i,j+1,k))*g22(i,j,k)
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
