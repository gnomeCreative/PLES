!**********************************************************************
subroutine buffvect2d(rbuff,var,n)
   !**********************************************************************
   !
   ! aggiorna il vettore rbuffd per lo scambio di variabili

   use scala3

   implicit none
   !
   double precision var(n2)
   double precision rbuff((n1+2)*(n2+2)*40)
   integer j,n,m

   do j=1,jy

      m=(n-1)*jy+j
      var(j)=rbuff(m)

   enddo
   !
   return
end
