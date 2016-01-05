module subgrid

   use scala3
   !
   ! common per variabili subgrid cartesiani eddy-viscosity
   ! viene dichiarata anche la media, max e rms della costante
   ! del modello - tutti valori mediati sui piani di omogeneità xz
   !
   double precision sub(n2),sub11(n2),sub22(n2),sub33(n2)
   double precision         sub12(n2),sub13(n2),sub23(n2)
   double precision sus(n2),sus11(n2),sus22(n2),sus33(n2)
   double precision         sus12(n2),sus13(n2),sus23(n2)
   double precision subrho11(n2),subrho22(n2),subrho33(n2)
   double precision susrho11(n2),susrho22(n2),susrho33(n2)
   real c11(n2),c22(n2),c33(n2)
   !
   common/subgrid1/sub,sub11,sub22,sub33,sub12,sub23,sub13
   common/subgrid2/sus,sus11,sus22,sus33,sus12,sus23,sus13
   common/subgrid3/subrho11,subrho22,subrho33
   common/subgrid4/susrho11,susrho22,susrho33
   common/costante/c11,c22,c33

end module subgrid
