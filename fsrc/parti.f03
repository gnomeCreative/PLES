module parti
   !--------------------------------------------------------------
   ! dichiarazione common per le particelle
   !
   !       integer ixpa(kpart),izpa(kpart)
   !       real xpart(kpart,3),upart(kpart,3),vp(kpart,3)
   !       real diam(kpart),ropart(kpart)
   !       real deltax,deltaz
   !       real ro1,ro2,ro3,di1,di2,di3
   real alx,alz,ti
   !
   !       common/parti/ixpa,izpa,xpart,upart,vp
   !       common/caratt/diam,ropart
   !       common/del/deltax,deltaz
   !       common/incaratt/ro1,ro2,ro3,di1,di2,di3
   common/dominio/alx,alz
   common/tempopa/ti

end module parti
