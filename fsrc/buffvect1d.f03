!**********************************************************************
subroutine buffvect1d(sbuff,var,n)
   !**********************************************************************
   !
   ! aggiorna il vettore sbuffd per lo scambio di variabili

   use scala3

   implicit none
   !
   double precision var(n2)
   double precision sbuff((n1+2)*(n2+2)*40)
   integer j,n,m

   do j=1,jy

      m=(n-1)*jy+j
      sbuff(m)=var(j)

   enddo
   !
   return
end
