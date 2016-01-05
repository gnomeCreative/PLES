!-------------------------------------------------------------
! dichiarazione area di common e dimensionamenti per codice 3d
!-------------------------------------------------------------

module scala3

! SONO DEFINITI:
! i parametri n1,n2,n3 massimi di dimensione
! il numero di scalari che si considerano nscal
! il numero di celle in x,y,z
!
      integer n1,n2,n3
      integer nscal
      integer jx,jy,jz
      integer potenziale
      integer single_or_double
      real dt,dt_start,re
!
      parameter (n1=32,n2=32,n3=32)
      parameter (nscal=1) !MUST NOT BE LESS THAN 1
!
      parameter (single_or_double=1)
!
!     definizione di singola o doppia precisione per chiamate MPI
!     viene definito un datatype MPI_REAL_SD contenuto nell'area di
!     common tipologia.h      
!     singola precisione = 0
!     doppia  precisione = 1
!     IMPORTANTE se doppia compilare con -qrealsize=8 o flag equivalente 
!     per ridefinire i real in real*8      
!
      !common/tervi/jx,jy,jz
      !common/tempo/dt,dt_start,re
      !common/pot/potenziale

end module scala3

