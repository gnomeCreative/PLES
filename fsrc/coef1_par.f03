!***********************************************************************
subroutine coef1_par(j,k,del,aa,rh,kparasta,kparaend)
   !***********************************************************************
   ! compute the coefficent of the band tridiagonal matrix for solution
   ! along csi for equation u,v,  w
   !
   use myarrays_metri3
   use myarrays_velo3
   use scala3
   use period

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k,ii
   integer kparasta,kparaend
   real del(0:n1+1,0:n2+1,kparasta-1:kparaend+1),aa(3,*),rh(*)
   !-----------------------------------------------------------------------
   !
   do ii=1,ip
      !     side left
      i=1
      !
      aa(1,i)=0.
      !
      aa(2,i)=.5*(annit(i,j,k)+annit(i+1,j,k))*g11(i,j,k)+ &
         3.*annit(i-1,j,k)*g11(i-1,j,k)
      aa(2,i)=1.+aa(2,i)*dt/giac(i,j,k)/2.
      !
      aa(3,i)=.5*(annit(i,j,k)+annit(i+1,j,k))*g11(i,j,k)+ &
         annit(i-1,j,k)*g11(i-1,j,k)/3.
      aa(3,i)=-dt*aa(3,i)/giac(i,j,k)/2.
      !
      rh(i)=rhs(i,j,k)+ &
         annit(i-1,j,k)*g11(i-1,j,k)* &
         del(i-1,j,k)*8.*dt/giac(i,j,k)/6.
      !
      !     side right
      i=jx
      !
      aa(1,i)=annit(i+1,j,k)*g11(i,j,k)/3.+ &
         .5*(annit(i,j,k)+annit(i-1,j,k))*g11(i-1,j,k)
      aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
      !
      aa(2,i)=annit(i+1,j,k)*g11(i,j,k)*3.+ &
         .5*(annit(i,j,k)+annit(i-1,j,k))*g11(i-1,j,k)
      aa(2,i)=1.+dt*aa(2,i)/giac(i,j,k)/2.
      !
      aa(3,i)=0.
      !
      rh(i)=rhs(i,j,k)+ &
         annit(i+1,j,k)*g11(i,j,k)* &
         del(i+1,j,k)*8.*dt/giac(i,j,k)/6.
   !
   enddo
   !
   !     into the field
   do i=1+ip,jx-ip
      !
      aa(1,i)=.5*(annit(i,j,k)+annit(i-1,j,k))*g11(i-1,j,k)
      aa(1,i)=-dt*aa(1,i)/giac(i,j,k)/2.
      !
      aa(3,i)=.5*(annit(i+1,j,k)+annit(i,j,k))*g11(i,j,k)
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
