!***********************************************************************
module myarrays_wallmodel
    !***********************************************************************
    !      for wall modeling
    !-----------------------------------------------------------------------
    use iso_c_binding

    real,allocatable :: tauw3x(:,:),tauw4x(:,:)
    real,allocatable :: tauw3z(:,:),tauw4z(:,:)
    real,allocatable :: punto_wfp3(:,:,:,:)
    real,allocatable :: punto_wfp4(:,:,:,:)
    real,allocatable :: u_t(:,:,:)

    integer,allocatable ::att_mod_par(:,:,:)
    real,allocatable :: utangente(:,:,:)

    integer,bind(C) :: wfp1,wfp2,wfp3,wfp4,wfp5,wfp6
    integer :: eseguo34
    integer,bind(C) :: att_wm_sgs,rough
    real,bind(C) :: z0

    contains

    subroutine initialize_wallmodel(coef_wall)

        use scala3, only: jx
        use mysending, only: kparasta,kparaend

        implicit none

        integer,intent(in) :: coef_wall

        if (coef_wall == 0) then
            wfp1 = 0
            wfp2 = 0
            wfp3 = 0
            wfp4 = 0
            wfp5 = 0
            wfp6 = 0
            att_wm_sgs=0
        else if (coef_wall==1) then
            allocate(att_mod_par(1:jx,2,kparasta:kparaend))
            allocate(u_t(1:jx,2,kparasta:kparaend))
            allocate(utangente(1:jx,2,kparasta:kparaend))
            allocate(punto_wfp3(3,3,jx,kparasta:kparaend))
            allocate(punto_wfp4(3,3,jx,kparasta:kparaend))
            att_mod_par=0
            utangente=1.
            u_t=1.

            if (wfp3==1) then
                att_mod_par(:,1,:)=1 ! first index=1 : face 3
            end if

            if (wfp4==1) then
                att_mod_par(:,2,:)=1 ! second index=2 : face 4
            end if

        end if

    end subroutine initialize_wallmodel


!***********************************************************************
end module myarrays_wallmodel
!***********************************************************************
