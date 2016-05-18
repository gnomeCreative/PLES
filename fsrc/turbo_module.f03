module turbo_module

    !  array for turbulence model
    !-----------------------------------------------------------------------

    ! THIS IS THE OLD TURBO2_DATA

    real,allocatable :: s11f(:,:,:),smods11f(:,:,:)
    real,allocatable :: s12f(:,:,:),smods12f(:,:,:)
    real,allocatable :: s13f(:,:,:),smods13f(:,:,:)
    real,allocatable :: s21f(:,:,:),smods21f(:,:,:)
    real,allocatable :: s22f(:,:,:),smods22f(:,:,:)
    real,allocatable :: s23f(:,:,:),smods23f(:,:,:)
    real,allocatable :: s31f(:,:,:),smods31f(:,:,:)
    real,allocatable :: s32f(:,:,:),smods32f(:,:,:)
    real,allocatable :: s33f(:,:,:),smods33f(:,:,:)
    !
    real,allocatable :: rho11f(:,:,:),smodrho11f(:,:,:)
    real,allocatable :: rho22f(:,:,:),smodrho22f(:,:,:)
    real,allocatable :: rho33f(:,:,:),smodrho33f(:,:,:)
    real,allocatable :: smodf(:,:,:)
    real,allocatable :: fil11(:,:,:),fil12(:,:,:),fil13(:,:,:)
    real,allocatable :: fil21(:,:,:),fil22(:,:,:),fil23(:,:,:)
    real,allocatable :: fil31(:,:,:),fil32(:,:,:),fil33(:,:,:)
    real,allocatable :: fipp11(:,:,:),fipp12(:,:,:),fipp13(:,:,:)
    real,allocatable :: fipp21(:,:,:),fipp22(:,:,:),fipp23(:,:,:)
    real,allocatable :: fipp31(:,:,:),fipp32(:,:,:),fipp33(:,:,:)
    real,allocatable :: m11(:,:,:),m12(:,:,:),m13(:,:,:)
    real,allocatable :: m21(:,:,:),m22(:,:,:),m23(:,:,:)
    real,allocatable :: m31(:,:,:),m32(:,:,:),m33(:,:,:)
    real,allocatable :: mrho11(:,:,:),mrho22(:,:,:),mrho33(:,:,:)
    real,allocatable :: am11(:,:,:),am12(:,:,:),am13(:,:,:)
    real,allocatable :: am21(:,:,:),am22(:,:,:),am23(:,:,:)
    real,allocatable :: am31(:,:,:),am32(:,:,:),am33(:,:,:)
    real,allocatable :: amrho11(:,:,:),amrho22(:,:,:)
    real,allocatable :: amrho33(:,:,:)
    !
    real,allocatable :: lmf11(:,:,:),lmf12(:,:,:),lmf13(:,:,:)
    real,allocatable :: lmf21(:,:,:),lmf22(:,:,:),lmf23(:,:,:)
    real,allocatable :: lmf31(:,:,:),lmf32(:,:,:),lmf33(:,:,:)
    !
    real,allocatable :: rhof(:,:,:),uf(:,:,:),vf(:,:,:),wf(:,:,:)
    real,allocatable :: ucof(:,:,:),vcof(:,:,:),wcof(:,:,:)
    real,allocatable :: uucof(:,:,:),uvcof(:,:,:),uwcof(:,:,:)

    real,allocatable :: vucof(:,:,:),vvcof(:,:,:),vwcof(:,:,:)
    real,allocatable :: wucof(:,:,:),wvcof(:,:,:),wwcof(:,:,:)
    real,allocatable :: rhoucof(:,:,:),rhovcof(:,:,:)
    real,allocatable :: rhowcof(:,:,:)
    !
    real,allocatable :: l11(:,:,:),l12(:,:,:),l13(:,:,:)
    real,allocatable :: l21(:,:,:),l22(:,:,:),l23(:,:,:)
    real,allocatable :: l31(:,:,:),l32(:,:,:),l33(:,:,:)
    real,allocatable :: lrho11(:,:,:),lrho22(:,:,:),lrho33(:,:,:)
    !
    real,allocatable :: ass11(:,:,:),ass12(:,:,:),ass13(:,:,:)
    real,allocatable :: ass21(:,:,:),ass22(:,:,:),ass23(:,:,:)
    real,allocatable :: ass31(:,:,:),ass32(:,:,:),ass33(:,:,:)
    real,allocatable :: assrho11(:,:,:),assrho22(:,:,:)
    real,allocatable :: assrho33(:,:,:)
    !

    real,allocatable :: m11m(:,:,:),m12m(:,:,:),m13m(:,:,:)
    real,allocatable :: m21m(:,:,:),m22m(:,:,:),m23m(:,:,:)
    real,allocatable :: m31m(:,:,:),m32m(:,:,:),m33m(:,:,:)
    !
    real,allocatable :: al11(:,:,:),al12(:,:,:),al13(:,:,:)
    real,allocatable :: al21(:,:,:),al22(:,:,:),al23(:,:,:)
    real,allocatable :: al31(:,:,:),al32(:,:,:),al33(:,:,:)
    real,allocatable :: alrho11(:,:,:),alrho22(:,:,:)
    real,allocatable :: alrho33(:,:,:)
    !
    real,allocatable :: rhofb(:,:,:)
    real,allocatable :: ucofb(:,:,:),vcofb(:,:,:),wcofb(:,:,:)
    real,allocatable :: lmfb11(:,:,:),lmfb12(:,:,:),lmfb13(:,:,:)
    real,allocatable :: lmfb21(:,:,:),lmfb22(:,:,:),lmfb23(:,:,:)
    real,allocatable :: lmfb31(:,:,:),lmfb32(:,:,:),lmfb33(:,:,:)
    !
    real,allocatable :: uucofb(:,:,:),uvcofb(:,:,:),uwcofb(:,:,:)
    real,allocatable :: vucofb(:,:,:),vvcofb(:,:,:),vwcofb(:,:,:)
    real,allocatable :: wucofb(:,:,:),wvcofb(:,:,:),wwcofb(:,:,:)
    real,allocatable :: rhoucofb(:,:,:),rhovcofb(:,:,:)
    real,allocatable :: rhowcofb(:,:,:)
    !
    real,allocatable :: ufucof(:,:,:),ufvcof(:,:,:),ufwcof(:,:,:)
    real,allocatable :: vfucof(:,:,:),vfvcof(:,:,:),vfwcof(:,:,:)
    real,allocatable :: wfucof(:,:,:),wfvcof(:,:,:),wfwcof(:,:,:)
    real,allocatable :: rhofucof(:,:,:),rhofvcof(:,:,:)
    real,allocatable :: rhofwcof(:,:,:)
    !
    real,allocatable :: bl11(:,:,:),bl12(:,:,:),bl13(:,:,:)
    real,allocatable :: bl21(:,:,:),bl22(:,:,:),bl23(:,:,:)
    real,allocatable :: bl31(:,:,:),bl32(:,:,:),bl33(:,:,:)
    real,allocatable :: blrho11(:,:,:),blrho22(:,:,:)
    real,allocatable :: blrho33(:,:,:)
    !
    real,allocatable :: sgs11(:,:,:),sgs12(:,:,:),sgs13(:,:,:)
    real,allocatable :: sgs21(:,:,:),sgs22(:,:,:),sgs23(:,:,:)
    real,allocatable :: sgs31(:,:,:),sgs32(:,:,:),sgs33(:,:,:)
    real,allocatable :: sgsrho11(:,:,:),sgsrho22(:,:,:)
    real,allocatable :: sgsrho33(:,:,:)
    !
    real,allocatable :: piano1(:,:),piano2(:,:)
    real,allocatable :: piano3(:,:),piano4(:,:)
    real,allocatable :: piano5(:,:),piano6(:,:)
    real,allocatable :: piano7(:,:),piano8(:,:)
    real,allocatable :: piano9(:,:),piano10(:,:)
    real,allocatable :: piano11(:,:),piano12(:,:),piano13(:,:)
    real,allocatable :: piano14(:,:),piano15(:,:)
    real,allocatable :: piano16(:,:),piano17(:,:)
    real,allocatable :: piano18(:,:),piano19(:,:)
    real,allocatable :: piano20(:,:),piano21(:,:)
    real,allocatable :: piano22(:,:),piano23(:,:)
    real,allocatable :: piano24(:,:),piano25(:,:)
    real,allocatable :: piano26(:,:),piano27(:,:)
    real,allocatable :: piano28(:,:),piano29(:,:)
    real,allocatable :: piano30(:,:)


    real,allocatable :: appo1(:,:,:)
    real,allocatable :: appo2(:,:,:)
    real,allocatable :: appo1rho(:,:,:)
    real,allocatable :: appo2rho(:,:,:)
    real,allocatable :: appo1_piano(:,:,:)
    real,allocatable :: appo2_piano(:,:,:)
    real,allocatable :: appo1rho_piano(:,:,:)
    real,allocatable :: appo2rho_piano(:,:,:)

    real,allocatable :: alalm(:,:,:)
    real,allocatable :: alamm(:,:,:)
    real,allocatable :: alalmrho(:,:,:)
    real,allocatable :: alammrho(:,:,:)


    !      real,allocatable :: piano1dp(:,:)
    !      real,allocatable :: piano2dp(:,:)
    !      real,allocatable :: piano3dp(:,:)
    !      real,allocatable :: piano4dp(:,:)

    ! THIS IS THE OLD TURBO3BIS

    real,allocatable :: smod(:,:,:)
    real,allocatable :: smodV(:,:,:)
    real,allocatable :: smodH(:,:,:)

    real,allocatable :: s11(:,:,:)
    real,allocatable :: s12(:,:,:)
    real,allocatable :: s13(:,:,:)
    real,allocatable :: s21(:,:,:)
    real,allocatable :: s22(:,:,:)
    real,allocatable :: s23(:,:,:)
    real,allocatable :: s31(:,:,:)
    real,allocatable :: s32(:,:,:)
    real,allocatable :: s33(:,:,:)
    !
    real,allocatable :: uco(:,:,:)
    real,allocatable :: vco(:,:,:)
    real,allocatable :: wco(:,:,:)
    real,allocatable :: uuco(:,:,:)
    real,allocatable :: uvco(:,:,:)
    real,allocatable :: uwco(:,:,:)
    real,allocatable :: vuco(:,:,:)
    real,allocatable :: vvco(:,:,:)
    real,allocatable :: vwco(:,:,:)
    real,allocatable :: wuco(:,:,:)
    real,allocatable :: wvco(:,:,:)
    real,allocatable :: wwco(:,:,:)
    !
    real,allocatable :: rhofl(:,:,:)
    real,allocatable :: rho11(:,:,:,:)
    real,allocatable :: rho22(:,:,:,:)
    real,allocatable :: rho33(:,:,:,:)
    real,allocatable :: rhouco(:,:,:)
    real,allocatable :: rhovco(:,:,:)
    real,allocatable :: rhowco(:,:,:)
    !
    real,allocatable :: apcsx(:,:,:)
    real,allocatable :: apcsy(:,:,:)
    real,allocatable :: apcsz(:,:,:)
    real,allocatable :: apetx(:,:,:)
    real,allocatable :: apety(:,:,:)
    real,allocatable :: apetz(:,:,:)
    real,allocatable :: apztx(:,:,:)
    real,allocatable :: apzty(:,:,:)
    real,allocatable :: apztz(:,:,:)
    !
    real,allocatable :: ap11(:,:,:)
    real,allocatable :: ap12(:,:,:)
    real,allocatable :: ap13(:,:,:)
    real,allocatable :: ap21(:,:,:)
    real,allocatable :: ap22(:,:,:)
    real,allocatable :: ap23(:,:,:)
    real,allocatable :: ap31(:,:,:)
    real,allocatable :: ap32(:,:,:)
    real,allocatable :: ap33(:,:,:)
    real,allocatable :: pp0(:,:,:)
    real,allocatable :: pp1(:,:,:)
    real,allocatable :: pp2(:,:,:)
    real,allocatable :: pp3(:,:,:)
    real,allocatable :: pc1(:,:,:)
    real,allocatable :: pc2(:,:,:)
    real,allocatable :: pc3(:,:,:)
    real,allocatable :: p0a(:,:,:)
    real,allocatable :: p0b(:,:,:)

    real,allocatable :: sbuff1(:)
    real,allocatable :: rbuff1(:)

contains

    subroutine initialize_turbo()

        use scala3, only: n1,n2,n3,nscal
        use mysending, only: kparasta,kparaend,myid,nproc
        use mysettings, only: inmod,inmodrho,nsgs

        implicit none

        if (myid==0) write(*,*) myid,'start allocation for turbo_statico'

        if (inmod .or. inmodrho .or. nsgs>=2) then

            allocate (m21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (ucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (l11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (l33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (ass11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (ass33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (piano4(0:n1+1,0:n2+1))
            allocate (piano5(0:n1+1,0:n2+1))
            allocate (piano6(0:n1+1,0:n2+1))
            allocate (piano7(0:n1+1,0:n2+1))
            allocate (piano8(0:n1+1,0:n2+1))
            allocate (piano9(0:n1+1,0:n2+1))
            allocate (piano10(0:n1+1,0:n2+1))
            allocate (piano11(0:n1+1,0:n2+1))
            allocate (piano12(0:n1+1,0:n2+1))
            allocate (piano13(0:n1+1,0:n2+1))
            allocate (piano14(0:n1+1,0:n2+1))
            allocate (piano15(0:n1+1,0:n2+1))
            allocate (piano16(0:n1+1,0:n2+1))
            allocate (piano17(0:n1+1,0:n2+1))
            allocate (piano18(0:n1+1,0:n2+1))
            allocate (piano19(0:n1+1,0:n2+1))
            allocate (piano20(0:n1+1,0:n2+1))
            allocate (piano21(0:n1+1,0:n2+1))
            allocate (piano22(0:n1+1,0:n2+1))

            allocate (uco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uuco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (uvco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (uwco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vuco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vvco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vwco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wuco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wvco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wwco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        end if

        if (inmod .or. nsgs >=2) then

            allocate (lmf11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lmf33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (uvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (uwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (vwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (wvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (m11m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m12m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m13m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m21m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m22m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m23m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m31m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m32m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m33m(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        end if

        if (inmodrho .or. nsgs>=2) then
            allocate (rhof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhofl(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        end if
        !
        allocate (piano1(0:n1+1,0:n2+1))
        allocate (piano2(0:n1+1,0:n2+1))
        allocate (piano3(0:n1+1,0:n2+1))
        !
        !     now turbo_statico matrix
        !
        allocate(smod(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(smodV(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate(smodH(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        !

        allocate (apcsx(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apcsy(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apcsz(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apetx(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apety(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apetz(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apztx(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apzty(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (apztz(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        !
        allocate (pp0(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pp1(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pp2(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pp3(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pc1(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pc2(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (pc3(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (p0a(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (p0b(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
        allocate (ap33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        allocate (rbuff1((n1+2)*(n2+2)*40))
        allocate (sbuff1((n1+2)*(n2+2)*40))

        !-----------------------------------------------------------------------
        !***********************************************************************
        !-----------------------------------------------------------------------
        !     allocation for dynamic procedure

        if (nsgs>=2) then


            allocate (s11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s12f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s13f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s21f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s23f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s31f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s32f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods12f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods13f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods21f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods23f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods31f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods32f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smods33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (rho11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodrho11f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodrho22f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodrho33f(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (smodf(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fil33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (fipp33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (m11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (m32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !      allocate (m33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (mrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (mrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (mrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (am33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (amrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (amrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (amrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))


            allocate (rhoucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhovcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhowcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !

            allocate (lrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (lrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

            allocate (sgs11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgs33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgsrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgsrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (sgsrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !
            allocate (   appo1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))
            allocate (   appo2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))
            allocate (appo1rho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))
            allocate (appo2rho(0:n1+1,0:n2+1,kparasta-1:kparaend+1)) !0:n3+1))

            if (myid==0) then
                allocate(   appo1_piano(0:n1+1,0:n2+1,n3:n3)) !0:n3+1))
                allocate(   appo2_piano(0:n1+1,0:n2+1,n3:n3)) !0:n3+1))
                allocate(appo1rho_piano(0:n1+1,0:n2+1,n3:n3)) !0:n3+1))
                allocate(appo2rho_piano(0:n1+1,0:n2+1,n3:n3)) !0:n3+1))
            elseif (myid==nproc-1) then
                allocate(   appo1_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
                allocate(   appo2_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
                allocate(appo1rho_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
                allocate(appo2rho_piano(0:n1+1,0:n2+1,1:1)) !0:n3+1))
            end if

            allocate(   alalm(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(   alamm(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(alalmrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate(alammrho(0:n1+1,0:n2+1,kparasta-1:kparaend+1))


            if (inmod) then

                allocate (assrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (assrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (assrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (al11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (al33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (alrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (alrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (alrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (lmfb11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (lmfb33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

                allocate (rhofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (ucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (uucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (uvcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (uwcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vvcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vwcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wvcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wwcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhoucofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhovcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhowcofb(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

                allocate (ufucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (ufvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (ufwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vfucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vfvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (vfwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wfucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wfvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (wfwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhofucof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhofvcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (rhofwcof(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                !
                allocate (bl11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (bl33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (blrho11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (blrho22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
                allocate (blrho33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            end if

            !      allocate (piano22(0:n1+1,0:n2+1))
            allocate (piano23(0:n1+1,0:n2+1))
            allocate (piano24(0:n1+1,0:n2+1))
            allocate (piano25(0:n1+1,0:n2+1))
            allocate (piano26(0:n1+1,0:n2+1))
            allocate (piano27(0:n1+1,0:n2+1))
            allocate (piano28(0:n1+1,0:n2+1))
            allocate (piano29(0:n1+1,0:n2+1))
            allocate (piano30(0:n1+1,0:n2+1))
            !
            !
            !     allocation for turbo3bis:
            !
            allocate (s11(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s12(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s13(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s21(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s22(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s23(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s31(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s32(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (s33(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            !
            !      allocate (rhofl(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho11(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho22(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rho33(nscal,0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhouco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhovco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
            allocate (rhowco(0:n1+1,0:n2+1,kparasta-1:kparaend+1))

        end if ! dynamic

    end subroutine initialize_turbo

    subroutine execute_turbo(ktime,i_rest,in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kpstamg,kpendmg)

        use mysettings, only: nsgs,pran,prsc
        use scala3, only: re,jx,jy,n1,n2,nscal
        use mysending, only: kparasta,kparaend,myid,nproc
        use myarrays_metri3, only: annit,annitV,annit_piano,annitV_piano
        use myarrays_velo3, only: akapt,akaptV,akapt_piano,akaptV_piano

        implicit none

        integer,intent(in) :: ktime,i_rest
        integer,intent(in) :: kpstamg(0:4),kpendmg(0:4)
        integer,intent(in) :: in_dx1(n1,n2,kpstamg(1):kpendmg(1))
        integer,intent(in) :: in_sn1(n1,n2,kpstamg(1):kpendmg(1))
        integer,intent(in) :: in_sp1(n1,n2,kpstamg(1):kpendmg(1))
        integer,intent(in) :: in_st1(n1,n2,kpstamg(1):kpendmg(1))
        integer,intent(in) :: in_av1(n1,n2,kpstamg(1):kpendmg(1))
        integer,intent(in) :: in_in1(n1,n2,kpstamg(1):kpendmg(1))

        integer i,j,k,isc

        ! -----------------------------------------------------------------

        ! this external loop is a temporary solution to have the eddy diffusivity computed
        ! for each scalar without reynolds analogy

        if (nsgs==1) then
            ! compute the eddy viscosity and diffusivity with smagorinksy
            ! model with fixed constant
            call turbo_statico(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kpstamg,kpendmg)

        else if (nsgs==2 .or. nsgs==3) then
            !   compute the eddy viscosity and diffusivity with smagorinksy
            !   model with dynamic procedure for the constant
            call turbo_lagrdin(ktime,in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1,kpstamg,kpendmg)

        else if (nsgs==0) then
            ! if DNS:
            do k=kparasta-1,kparaend+1 !0,jz+1
                do j=0,jy+1
                    do i=0,jx+1
                        ! eddy viscosity
                        annit(i,j,k) =1./re
                        annitV(i,j,k)=1./re
                        ! eddy diffusivity
                        do isc=1,nscal
                            akapt(isc,i,j,k) =1./re/pran(isc)
                            akaptV(isc,i,j,k)=1./re/pran(isc)
                        end do
                    end do
                end do
            end do
            if (myid==0 .or. myid==nproc-1) then
                annit_piano  = 1./re
                annitV_piano = 1./re
            end if

            if (myid==0 .or. myid==nproc-1) then
                do isc=1,nscal
                    akapt_piano  = 1./re/pran(isc)
                    akaptV_piano = 1./re/pran(isc)
                end do
            end if

        end if

    end subroutine execute_turbo

subroutine periodic(r1,r2,r3)
    !***********************************************************************
    ! set periodic boundary condition on model matrix
    !
    use scala3
    use period
    use mysending, only: kparasta,kparaend
    !
    use mpi

    implicit none
    !
    !-----------------------------------------------------------------------
    !     array declaration
    integer i,j,k

    real r1(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
    real r2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
    real r3(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
    !-----------------------------------------------------------------------
    !
    !      periodicity on csi
    !
    if (ip==0) then
        do k=kparasta,kparaend
            do j=1,jy
                r1(0,j,k)=r1(jx,j,k)
                r2(0,j,k)=r2(jx,j,k)
                r3(0,j,k)=r3(jx,j,k)
                !
                r1(jx+1,j,k)=r1(1,j,k)
                r2(jx+1,j,k)=r2(1,j,k)
                r3(jx+1,j,k)=r3(1,j,k)
            end do
        end do
    end if
    !
    !      periodicity on eta
    !
    if (jp==0) then
        do k=kparasta,kparaend
            do i=1,jx
                r1(i,0,k)=r1(i,jy,k)
                r2(i,0,k)=r2(i,jy,k)
                r3(i,0,k)=r3(i,jy,k)
                !
                r1(i,jy+1,k)=r1(i,1,k)
                r2(i,jy+1,k)=r2(i,1,k)
                r3(i,jy+1,k)=r3(i,1,k)
            end do
        end do
    end if
    !
    !      periodicity on zeta made in a next step
    !
    return
end

!***********************************************************************
end module turbo_module
!***********************************************************************
