module period

    use iso_c_binding
   !
   ! common relativo alla periodicità
   ! indici dichiarati in visco.in
   ! che definiscono le BC periodiche
   !
   integer(kind=c_int),bind(C) :: ip,jp,kp
   !logical(kind=c_bool),bind(C) :: csi_periodic,eta_periodic,zita_periodic
   !
   !common/periodicity/ip,jp,kp

end module period
