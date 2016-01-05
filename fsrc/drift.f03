!***********************************************************************
subroutine drift(myid,nproc,kparasta,kparaend)
   !***********************************************************************
   ! compute Stokes drift velocity ( Noh, 2003)

   use myarrays_metri3
   use myarrays_LC
   !
   use scala3
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer i,j,k
   integer myid,nproc,kparasta,kparaend
   real pig,g,Aus,Aws
   !-----------------------------------------------------------------------
   !
   pig = 3.14159
   g   = 9.81

   !      h_0   = 1.  ! wave amplitude
   !      lamb= 40.   ! wave length
   !      h_0    = .24
   !      lamb = 12.
      
   Aus=(pig*h_0/lamb)*(pig*h_0/lamb)*(sqrt((g*lamb)/(2.*pig)))* &
      cos(betaw)
   Aws=(pig*h_0/lamb)*(pig*h_0/lamb)*(sqrt((g*lamb)/(2.*pig)))*  &
      sin(betaw)
      
   do k=kparasta-1,kparaend+1
      do j=0,n2+1
         do i=0,n1+1
            u_drift(i,j,k) = Aus*exp(-4*pig* &
               ( y(i,n2,k)-y(i,j,k) )/lamb)
            w_drift(i,j,k) = Aws*exp(-4*pig* &
               ( y(i,n2,k)-y(i,j,k) )/lamb)
         end do
      end do
   end do


   return
end subroutine

