!***********************************************************************
module turbo2_data
    !***********************************************************************
    ! array for turbulence model
    !-----------------------------------------------------------------------

    use iso_c_binding

    real,bind(C) :: cost,costH,costV
       
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

contains

    subroutine initialize_turbo2()

        use mysending
        use scala3

        implicit none

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

    end subroutine initialize_turbo2

!***********************************************************************
end module turbo2_data
!***********************************************************************
