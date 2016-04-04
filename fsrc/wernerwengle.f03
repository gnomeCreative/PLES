subroutine wernerwengle(distanza,ustar,vtan)
    ! wall function with Werner-Wengle:
    ! if y<11.8       u+ = y+              (linear)
    !    y>11.8       u+ = 8.3*y+**1/7     (exponential)
    !-----------------------------------------------------------------------

    use scala3, only: re

    implicit none
    !---------------------------------------------------------------
    real,intent(out) :: ustar
    real,intent(in) :: distanza
    real,intent(in) :: vtan

    !---------------------------------------------------------------
    ! parameter Werner-Wengle
    real,parameter :: wwa=8.3
    real,parameter :: wwb=1.0/7.0
    real,parameter :: const1=0.5*(1.0-wwb)*wwa**((1.0+wwb)/(1.0-wwb))
    real,parameter :: const2=(1.0+wwb)/wwa
    real,parameter :: const3=wwa**(2.0/(1.0-wwb))
    real,parameter :: const4=2.0/(1.0+wwb)

    !-----------------------------------------------------------------------
    ! compute wall shear stress with Werner-Wengle power law (1/7)
    ! for reference, see:
    ! http://aerojet.engr.ucdavis.edu/fluenthelp/html/ug/node516.htm
    !
    !---------------------------------------------------------------
    real :: visco,delta_y,inv_delta_y,ratio
    real :: tau_w,u_limit


    visco=1.0/re

    delta_y=2.0*distanza
    ratio=visco/delta_y
    inv_delta_y=1.0/delta_y

    ! limit velocity for viscous sublayer
    u_limit=0.5*visco*inv_delta_y*const3

    if (vtan<u_limit) then
        tau_w=2.0*visco*vtan*inv_delta_y
    else
        tau_w=(const1*ratio**(1.0+wwb)+const2*vtan*ratio**wwb)**const4
    end if

    ustar=sqrt(tau_w)


!    !-----------------------------------------------------------------------
!    !    Alessandro : I believe this old method is wrong
!    !
!    !---------------------------------------------------------------
!    real d_centro_plus,dycell,rycell
!    real visco,sub
!    real dvtan,vtankr
!    real taupow,tausub
!
!
!    visco=1.0/re
!
!    d_centro_plus=distanza/visco
!    dycell=2.0*distanza
!    rycell=1.0/dycell
!
!    ! sub bring the function to linear profile if necessary,
!    !  if the node is in the viscous sublayer
!    !
!    vtankr=0.5*visco*rycell*const3
!    dvtan=vtankr-vtan
!    sub=max(sign(1.0,dvtan),0.0)
!
!    tausub=2.0*visco*vtan*rycell
!
!    taupow=(vtan/(wwa*(d_centro_plus)**wwb))**(7.0/4.0)
!
!    if (taupow<tausub) then
!        sub=1.0
!    end if
!
!    ustar=sqrt(sub*tausub+(1.0-sub)*taupow)


    return
end
