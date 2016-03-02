subroutine wernerwengle(distanza,ustar,vtan)
    ! wall function with Werner-Wengle:
    ! if y<11.8       u+ = y+              (linear)
    !    y>11.8       u+ = 8.3*y+**1/7     (exponential)
    !-----------------------------------------------------------------------

    ! for reference, see:
    ! http://aerojet.engr.ucdavis.edu/fluenthelp/html/ug/node516.htm

    use scala3, only: re

    implicit none
    !---------------------------------------------------------------
    real,intent(inout) :: ustar
    real,intent(in) :: distanza
    real,intent(in) :: vtan

    !---------------------------------------------------------------
    ! parameter Werner-Wengle
    real,parameter :: aaa=8.3
    real,parameter :: bbb=1.0/7.0
    real,parameter :: const1=0.5*(1.0-bbb)*aaa**((1.0+bbb)/(1.0-bbb))
    real,parameter :: const2=(1.0+bbb)/aaa
    real,parameter :: const3=aaa**(2.0/(1.0-bbb))
    real,parameter :: const4=2.0/(1.0+bbb)

    !---------------------------------------------------------------
    real d_centro_plus,dycell,rycell
    real visco,sub
    real dvtan,vtankr
    real taupow,tausub

    !-----------------------------------------------------------------------
    !    compute wall shear stress with WW POWER LAW 1/7
    !
    !distanza=dist_ib_parete(l)+dist_pp_ib(l)

    visco=1.0/re

    d_centro_plus=distanza/visco
    dycell=2.0*distanza
    rycell=1.0/dycell

    ! sub bring the function to linear profile if necessary,
    !  if the node is in the viscous sublayer
    !
    vtankr=0.5*visco*rycell*const3
    dvtan=vtankr-vtan
    sub=max(sign(1.0,dvtan),0.0)

    tausub=2.0*visco*vtan*rycell

    taupow=(vtan/(aaa*(d_centro_plus)**bbb))**(7.0/4.0)

    if (taupow<tausub) then
        sub=1.0
    end if

    ustar=sqrt(sub*tausub+(1.0-sub)*taupow)

    return
end
