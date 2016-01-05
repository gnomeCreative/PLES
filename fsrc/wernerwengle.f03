!***********************************************************************
subroutine wernerwengle(l,MN,vtan, &
    dist_ib_parete,dist_pp_ib,ustar)
    !***********************************************************************
    !     wall function with Werner - Wengle:
    !          se y<11.8       u+ = y+              linear
    !             y>11.8       u+ = 8.3*y+**1/7     exponential
    !-----------------------------------------------------------------------
    use myarrays2
    !
    use scala3
    use period

    implicit none
    !---------------------------------------------------------------
    !     array declaration
    !
    integer MN,MP
    integer l
    integer ptx,pty,ptz

    real aaa,bbb
    real const1,const2,const3,const4

    real dist_pp_ib(MN),dist_ib_parete(MN)

    real d_centro,d_centro_plus,dycell,rycell
    real visco,sub
    real vtan,vnor,dvtan,vtankr
    real taupow,tausub
    real ustar(MN)

    real distanza,distanza_ib_parete
    real distanza_ib_solida,coefparete
    real denominatore

    !-----------------------------------------------------------------------
    !     only for check
    ptx=1
    pty=1
    ptz=1 !25
    !
    !-----------------------------------------------------------------------
    !     parameter Werner-Wengle

    aaa=8.3
    bbb=0.1428571429

    const1 = 0.5*(1. - bbb)*aaa**((1. + bbb) /(1. - bbb))
    const2 = (1. + bbb) / aaa
    const3 = aaa ** (2. / (1. - bbb))
    const4 = 2. / (1. + bbb)
    !-----------------------------------------------------------------------
    !    compute wall shear stress with WW POWER LAW 1/7
    !
    d_centro=dist_ib_parete(l)+dist_pp_ib(l)

    visco=1./re

    d_centro_plus=d_centro/(1./re)
    dycell=2.*d_centro
    rycell = 1. / dycell
    !....................................................................
    ! sub bring the function to linear profile if necessary,
    !  if the node is in the viscous sublayer
    !
    vtankr = 0.5 * visco * rycell * const3
    dvtan  = vtankr - vtan
    sub    = MAX (SIGN(1.,dvtan),0.)
    !..................................................................
    tausub   = 2. * visco * vtan * rycell

    taupow  =  ( vtan/(aaa*(d_centro/visco)**bbb) )**(7./4.)

    if(taupow .lt. tausub)then
        sub=1.
    end if

    ustar(l)=sqrt(sub*tausub+(1.-sub)*taupow)

    return
end
