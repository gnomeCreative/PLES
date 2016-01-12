module subgrid

   !
   ! common per variabili subgrid cartesiani eddy-viscosity
   ! viene dichiarata anche la media, max e rms della costante
   ! del modello - tutti valori mediati sui piani di omogeneità xz
   !
   double precision,allocatable :: sub(:),sub11(:),sub22(:),sub33(:)
   double precision,allocatable ::         sub12(:),sub13(:),sub23(:)
   double precision,allocatable :: sus(:),sus11(:),sus22(:),sus33(:)
   double precision,allocatable ::         sus12(:),sus13(:),sus23(:)
   double precision,allocatable :: subrho11(:),subrho22(:),subrho33(:)
   double precision,allocatable :: susrho11(:),susrho22(:),susrho33(:)
   real,allocatable :: c11(:),c22(:),c33(:)
   !
   !common/subgrid1/sub,sub11,sub22,sub33,sub12,sub23,sub13
   !common/subgrid2/sus,sus11,sus22,sus33,sus12,sus23,sus13
   !common/subgrid3/subrho11,subrho22,subrho33
   !common/subgrid4/susrho11,susrho22,susrho33
   !common/costante/c11,c22,c33

end module subgrid
