! Bingwei Long  10/20/2018

module mod_vfunc_smpl

  use mod_obj_smplchn
  use potwrap
  use testwrap
  use delta_epmod
  use testwrap
  use chengdu
  implicit none

  integer :: V1_BASEN_sprbl_mvfs = 4,     &
    &        V2_BASEN_sprbl_mvfs = 4,    &
    &        V3_BASEN_sprbl_mvfs = 4


  real(NER) :: V1_BND_sprbl_mvfs = 1.0E-4_NER,  &
    &          V2_BND_sprbl_mvfs = 1.0E-4_NER,  &
    &          V3_BND_sprbl_mvfs = 1.0E-4_NER
contains

  subroutine SetUpVfunc_sprbl_smpl(chn)

    class(obj_smplchn), intent(inout) :: chn

    ! NLO pot does *not* vanishes
    chn%vanishing(1) = .false.
    chn%npar(1:4) = NPAR_SPR1S0_SMPL

    chn%pVfunc_V0 => VLO_sprbl_epwrap

    chn%pVfunc_V1 => VNLO_sprbl_epwrap
    chn%V1_cheb_base_n = V1_BASEN_sprbl_mvfs
    chn%V1_cheb_bnd = V1_BND_sprbl_mvfs

    chn%pVfunc_V2 => VN2LO_sprbl_epwrap
    chn%V2_cheb_base_n = V2_BASEN_sprbl_mvfs
    chn%V2_cheb_bnd = V2_BND_sprbl_mvfs

    chn%pVfunc_V3 => VN3LO_sprbl_epwrap
    chn%V3_cheb_base_n = V3_BASEN_sprbl_mvfs
    chn%V3_cheb_bnd = V3_BND_sprbl_mvfs

  end subroutine SetUpVfunc_sprbl_smpl

  subroutine SetUpVfunc_mmwly_smpl(chn)

    class(obj_smplchn), intent(inout) :: chn

    select case(chn%chnid)
      case (CHNID_1S0_PP)
        chn%vanishing(1) = .false.
        chn%npar(1:3) = (/1, 3, 6/)
        chn%pVfunc_V0 => VLO_withc0_epwrap
        chn%pVfunc_V1 => VNLO_MMWLY_epwrap
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER
!         chn%pVfunc_V2 => VN2LO_MMWLY_epwrap
! debug fit without TPE
        chn%pVfunc_V2 => VN2LO_MMWLY_SHRG_epwrap
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-7_NER

      case (CHNID_3S1_PP)
        chn% V2_cheb_pertBD_n = 5
        chn%V2_cheb_pertBD_bnd = 1.0E-6_NER
        chn%vanishing(1) = .true.
        chn%npar(1:3) = (/1, 1, 4/)
        chn%pVfunc_V0 => VLO_withc0_epwrap
! debug fit with Vbar
!         chn%pVfunc_V0 => VbarLO_withc0_epwrap

        chn%pVfunc_V1 => zero_V

!         chn%vlvs_separated(2) = .true.
!         chn%pVfunc_V2 => VN2LO_withc0_epwrap
!         chn%pVfunc_VL2 =>  TPE0_epwrap
!         chn%pVfunc_VS2 =>  VN2LO_withc0_SHRG_epwrap
! debug fit without TPE
        chn%pVfunc_V2 => VN2LO_withc0_SHRG_epwrap
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-6_NER

      case (CHNID_3P0_PP)
        chn%vanishing(1) = .true.
        chn%npar(1:3) = (/1, 1, 3/)
        chn%pVfunc_V0 => VLO_withc0_epwrap
        chn%pVfunc_V1 => zero_V
!         chn%pVfunc_V2 => VN2LO_withc0_epwrap
! debug fit without TPE
        chn%pVfunc_V2 => VN2LO_withc0_SHRG_epwrap
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-5_NER

      case (CHNID_3P2_PP, CHNID_3P1_PP, CHNID_1P1_PP)
        chn%vanishing(0) = .true.
        chn%npar(1:3) = (/0, 0, 1/)
        chn%pVfunc_V0 => zero_V
        chn%pVfunc_V1 => OPE_epwrap
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-3_NER


      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%vanishing(0) = .true.
        chn%npar(1:3) = (/0, 0, 0/)
        chn%pVfunc_V0 => zero_V
        chn%pVfunc_V1 => OPE_epwrap
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER
        chn%pVfunc_V2 => zero_V

      case default ! F and G waves
        chn%vanishing(0) = .true.
        chn%npar(1:3) = (/0, 0, 0/)
        chn%pVfunc_V0 => zero_V
        chn%pVfunc_V1 => OPE_epwrap
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER
        chn%pVfunc_V2 => zero_V
        chn%pVfunc_V2 => zero_V
    end select

  end subroutine SetUpVfunc_mmwly_smpl

! setup pVfunc with chnegdu for getting phsae shifts to verify the values of parameters in chengdu

  subroutine SetUpVfunc_mmwly_chengdu_smpl(chn, pot)

    class(obj_smplchn), intent(inout) :: chn
    character(len = 20), intent(in)   :: pot

    call chengdu_dispatch(pot)
!     debug
!     print *, 'pot = ', pot
    chn%pVfunc_V0 => LO_chengdu_epwrap
    chn%pVfunc_V1 => NLO_chengdu_epwrap
    chn%pVfunc_V2 => N2LO_chengdu_epwrap
    chn%npar(1:3) = (/0, 0, 0/)

    select case(chn%chnid)
      case (CHNID_1S0_PP)
        chn%vanishing(1) = .false.
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-7_NER

      case (CHNID_3S1_PP)
        chn% V2_cheb_pertBD_n = 5
        chn%V2_cheb_pertBD_bnd = 1.0E-6_NER
        chn%vanishing(1) = .true.
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-6_NER

      case (CHNID_3P0_PP)
        chn%vanishing(1) = .true.
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-5_NER

      case (CHNID_3P2_PP, CHNID_3P1_PP, CHNID_1P1_PP)
        chn%vanishing(0) = .true.
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER
        chn%V2_cheb_base_n = 4
        chn%V2_cheb_bnd = 1.0E-3_NER

      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%vanishing(0) = .true.
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER

      case default ! F and G waves
        chn%vanishing(0) = .true.
        chn%V1_cheb_base_n = 4
        chn%V1_cheb_bnd = 1.0E-3_NER
    end select

  end subroutine SetUpVfunc_mmwly_chengdu_smpl

  subroutine SetUpVfunc_withc0_smpl(chn)

    class(obj_smplchn), intent(inout) :: chn

    ! NLO pot vanishes
    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL
    else
      chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap
    ! debug
    chn%pVfunc_V1 => zero_V

    chn%pVfunc_V2 => VN2LO_withc0_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-5_NER

    chn%pVfunc_V3 => VN3LO_withc0_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine SetUpVfunc_withc0_smpl

!   subroutine SetUpVfunc_withc0_DR_smpl(chn)

!     class(obj_smplchn), intent(inout) :: chn

!     ! NLO pot vanishes
!     chn%vanishing(1) = .true.
!     if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
!       chn%npar(1:4) = (/1,1,3,6/)
!     else
!       chn%npar(1:4) = (/1,1,4,10/)
!     end if

!     chn%pVfunc_V0 => VLO_withc0_epwrap
!     ! debug
!     chn%pVfunc_V1 => zero_V

!     chn%pVfunc_V2 => VN2LO_withc0_epwrap
!     chn%V2_cheb_base_n = 4
!     chn%V2_cheb_bnd = 1.0E-3_NER

!     chn%pVfunc_V3 => VN3LO_withc0_DR_epwrap
!     chn%V3_cheb_base_n = 4
!     chn%V3_cheb_bnd = 1.0E-3_NER

!   end subroutine SetUpVfunc_withc0_DR_smpl

  ! for 3p0 and 3s1

  subroutine SetUp_withc0_del_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    ! PC_delta  =  305.0_NER
    ! PC_delta_Ha = 1.450_NER

    ! NLO pot vanishes
    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL
    else
      chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap

    chn%vlvs_separated(2) = .true.
    chn%pVfunc_V2  =>  VN2LO_withc0_del_epwrap
    chn%pVfunc_VL2 =>  TPE0_del_epwrap
    chn%pVfunc_VS2 =>  VN2LO_withc0_SHRG_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%vlvs_separated(3) = .true.
    chn%pVfunc_V3  =>  VN3LO_withc0_del_epwrap
    chn%pVfunc_VL3 =>  TPE1_del_epwrap
    chn%pVfunc_VS3 =>  VN3LO_withc0_SHRG_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine SetUp_withc0_del_smpl

  subroutine setUP_withc0_sfr_smpl(chn)

    class(obj_smplchn), intent(inout) :: chn

    ! NLO pot vanishes
    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL
    else
      chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap

    chn%vlvs_separated(2) = .true.
    chn%pVfunc_V2  =>  VN2LOSFR_withc0_epwrap
    chn%pVfunc_VL2 =>  TPESFR_epwrap
    chn%pVfunc_VS2 =>  VN2LO_withc0_SHRG_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%vlvs_separated(3) = .true.
    chn%pVfunc_V3  =>  VN3LOSFR_withc0_epwrap
    chn%pVfunc_VL3 =>  TPESFR1_epwrap
    chn%pVfunc_VS3 =>  VN3LO_withc0_SHRG_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER
  end subroutine setUP_withc0_sfr_smpl

  subroutine setUP_withc0_sfr2_smpl(chn)

    class(obj_smplchn), intent(inout) :: chn

    ! NLO pot vanishes
    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL
    else
      chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap
    ! debug
!     chn%pVfunc_V1 => zero_V

    chn%pVfunc_V2 => VN2LOSFR_withc0_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V3 => VN3LOSFR_withc0_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER
  end subroutine setUP_withc0_sfr2_smpl

  subroutine SetUp_withc0_SFRdel_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL
    else
      chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap

    chn%vlvs_separated(2) = .true.
    chn%pVfunc_V2  =>  VN2LOSFR_withc0_del_epwrap
    chn%pVfunc_VL2 =>  TPE0SFR_del_epwrap
    chn%pVfunc_VS2 =>  VN2LO_withc0_SHRG_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%vlvs_separated(3) = .true.
    chn%pVfunc_V3  =>  VN3LOSFR_withc0_del_epwrap
    chn%pVfunc_VL3 =>  TPE1SFR_del_epwrap
    chn%pVfunc_VS3 =>  VN3LO_withc0_SHRG_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine SetUp_withc0_SFRdel_smpl

  subroutine SetUp_withc0_SFRdelc_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = (/1,1,3,6/)
    else
      chn%npar(1:4) = (/1,1,4,10/)
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap

    chn%vlvs_separated(2) = .true.
    chn%pVfunc_V2  =>  VN2LOSFR_withc0_del_epwrap
    chn%pVfunc_VL2 =>  TPE0SFR_del_epwrap
    chn%pVfunc_VS2 =>  VN2LO_withc0_SHRG_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%vlvs_separated(3) = .true.
    chn%pVfunc_V3  =>  VN3LOSFR_withc0_delc_epwrap
    chn%pVfunc_VL3 =>  TPE1SFR_del_epwrap
    chn%pVfunc_VS3 =>  VN3LO_withc0_deltac_SHRG_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine SetUp_withc0_SFRdelc_smpl

  subroutine SetUp_withc0_SFR2del_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL
    else
      chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap

    chn%pVfunc_V2 =>  VN2LOSFR_withc0_del_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V3 => VN3LOSFR_withc0_del_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine SetUp_withc0_SFR2del_smpl

  ! Perturbative pion -- suppressed scheme and deltaless

  subroutine SetUp_pope_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => OPE_epwrap
    chn%V1_cheb_base_n = 4
    chn%V1_cheb_bnd = 1.0E-3_NER

    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 5/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
      chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 4/)
      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 1/)
      case default
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 0/)
    end select

    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V3 => VN3LO_pope_smpl_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V4 => VN4LO_pope_smpl_epwrap
    chn%V4_cheb_base_n = 4
    chn%V4_cheb_bnd = 1.0E-3_NER

  end subroutine SetUp_pope_smpl

    subroutine SetUp_pope_dr_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => OPE_epwrap
    chn%V1_cheb_base_n = 2
    chn%V1_cheb_bnd = 0.05_NER

    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 5/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
      chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 4/)
      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 1/)
      case default
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 0/)
    end select

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V3 => VN3LO_pope_smpl_epwrap
    chn%V3_cheb_base_n = 5
    chn%V3_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V4 => VN4LO_pope_smpl_epwrap
    chn%V4_cheb_base_n = 5
    chn%V4_cheb_bnd = 1.0E-6_NER

  end subroutine SetUp_pope_dr_smpl

  ! Perturbative pion -- suppressed scheme and deltaless
  subroutine SetUp_pope_del_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => OPE_epwrap
    chn%V1_cheb_base_n = 4
    chn%V1_cheb_bnd = 1.0E-3_NER

    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 5/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
      chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 4/)
      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 1/)
      case default
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 0/)
    end select

    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V3 => VN3LO_pope_del_smpl_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V4 => VN4LO_pope_del_smpl_epwrap
    chn%V4_cheb_base_n = 4
    chn%V4_cheb_bnd = 1.0E-3_NER

  end subroutine SetUp_pope_del_smpl

   ! Perturbative pion -- suppressed scheme and deltaless with SFR

  subroutine SetUp_pope_del_SFR_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => OPE_epwrap
!     chn%V1_cheb_base_n = 2
    chn%V1_cheb_base_n = 2
!     chn%V1_cheb_bnd = 1.0E-3_NER
    chn%V1_cheb_bnd = 0.05_NER

!     chn%V2_cheb_base_n = 5
!     chn%V2_cheb_bnd = 1.0E-6_NER

    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 5/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
      chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 2, 4/)
      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 1/)
      case default
        chn%pVfunc_V2 => zero_V
        chn%npar(1:5) = (/0, 0, 0, 0, 0/)
    end select

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V3 => VN3LO_pope_del_SFR_smpl_epwrap
    chn%V3_cheb_base_n = 5
    chn%V3_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V4 => VN4LO_pope_del_SFR_smpl_epwrap
    chn%V4_cheb_base_n = 5
    chn%V4_cheb_bnd = 1.0E-6_NER

  end subroutine SetUp_pope_del_SFR_smpl

! the next sub for pope in smpl TPE lower than OPE only one order

  subroutine SetUp_del_SFR_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 => VNLO_del_SFR_smpl_epwrap
    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%npar(1:4) = (/ 0, 1, 2, 5/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
        chn%npar(1:4) = (/ 0, 1, 2, 4/)
      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%npar(1:4) = (/ 0, 0, 0, 1/)
      case default
        chn%npar(1:4) = (/ 0, 0, 0, 0/)
    end select

    chn%V1_cheb_base_n = 4
    chn%V1_cheb_bnd = 0.05_NER

    chn%pVfunc_V2 => VN2LO_del_SFR_smpl_epwrap
    chn%V2_cheb_base_n = 2
    chn%V2_cheb_bnd = 0.05_NER

    chn%pVfunc_V3 => VN3LO_del_SFR_smpl_epwrap
    chn%V3_cheb_base_n = 2
    chn%V3_cheb_bnd = 0.05_NER

  end subroutine SetUp_del_SFR_smpl
!
!the following code for 3P0

  subroutine setUP_withc0_sup3p0_smpl(chn)

    class(obj_smplchn), intent(inout) :: chn


    chn%vanishing(1) = .false.
    chn%npar(1:5) = (/ 0, 1, 2, 4, 6 /)


    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 =>  VLO_withc0_epwrap

    chn%pVfunc_V2 => VN2LO_supper3p0_smpl_epwrap

    chn%pVfunc_V3 => VN3LO_supper3p0_smpl_epwrap

    chn%pVfunc_V4 => VN4LO_supper3p0_smpl_epwrap

  end subroutine setUP_withc0_sup3p0_smpl

  subroutine setUP_withc0_3p0_smpl(chn)
    class(obj_smplchn),intent(inout)  ::  chn

    chn%vanishing(1) = .false.
    chn%npar(1:5) = (/ 0, 1, 3, 5, 5 /)

    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 =>  VLO_withc0_epwrap

    chn%pVfunc_V2 => VN2LO_per3p0_smpl_epwrap

    chn%pVfunc_V3 => VN3LO_per3p0_smpl_epwrap

  end subroutine setUP_withc0_3p0_smpl

  subroutine setUP_withc0_high3p0_smpl(chn)
    class(obj_smplchn),intent(inout)  ::  chn

    chn%vanishing(1) = .false.
    chn%npar(1:5) = (/ 0, 1, 3, 6, 6 /)

    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 =>  VLO_withc0_epwrap

    chn%pVfunc_V2 => VN2LO_per3p0_smpl_epwrap

    chn%pVfunc_V3 => VN3LO_per3p0_high_smpl_epwrap

  end subroutine setUP_withc0_high3p0_smpl

!   subroutine setUP_withc0_high3s1_smpl(chn)
!     class(obj_smplchn),intent(inout)  ::  chn

!     chn%vanishing(1) = .false.
!     chn%npar(1:5) = (/ 0, 1, 4, 10, 10 /)

!     chn%pVfunc_V0 => zero_V

!     chn%pVfunc_V1 =>  VLO_withc0_epwrap

!     chn%pVfunc_V2 => VN2LO_per3p0_smpl_epwrap

!     chn%pVfunc_V3 => VN3LO_per3p0_high_smpl_epwrap

!   end subroutine setUP_withc0_high3s1_smpl

  subroutine setUP_withc0_highwot3p0_smpl(chn)
    class(obj_smplchn),intent(inout)  ::  chn

    chn%vanishing(1) = .false.
    chn%npar(1:5) = (/ 0, 1, 3, 6, 6 /)

    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 =>  VLO_withc0_epwrap

    chn%pVfunc_V2 => VN2LO_perwot3p0_smpl_epwrap

    chn%pVfunc_V3 => VN3LO_perwot3p0_high_smpl_epwrap

  end subroutine setUP_withc0_highwot3p0_smpl

  subroutine setUP_withc0_3p0NLO_smpl(chn)
    class(obj_smplchn),intent(inout)  ::  chn

    chn%vanishing(1) = .false.
    chn%npar(1:5) = (/ 0, 1, 1, 1, 1 /)

    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 =>  VLO_withc0_epwrap

    chn%pVfunc_V2 => zero_V

    chn%pVfunc_V3 => zero_V

    chn%pVfunc_V4 => zero_V

  end subroutine setUP_withc0_3p0NLO_smpl

  subroutine setUP_withc0_m3p0_smpl(chn)
    class(obj_smplchn),intent(inout)  ::  chn

    chn%vanishing(1) = .false.
    chn%npar(1:5) = (/ 0, 6, 6, 6, 6 /)

    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 =>  VNLO_withc0_m3p0_epwrap

    chn%pVfunc_V2 => zero_V

    chn%pVfunc_V3 => zero_V

    chn%pVfunc_V4 => zero_V

  end subroutine setUP_withc0_m3p0_smpl

!#####################################
! add contact term V4 to N3LO
!#####################################

  subroutine SetUp_pope_del_SFR_smpl_test(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => OPE_epwrap
!     chn%V1_cheb_base_n = 2
    chn%V1_cheb_base_n = 2
!     chn%V1_cheb_bnd = 1.0E-3_NER
    chn%V1_cheb_bnd = 0.05_NER

!     chn%V2_cheb_base_n = 5
!     chn%V2_cheb_bnd = 1.0E-6_NER

    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 4, 10/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
      chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 3, 6/)
    end select

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V3 => VN3LO_pope_del_SFR_smpl_V3V4_epwrap
    chn%V3_cheb_base_n = 5
    chn%V3_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V4 => VN4LO_pope_del_SFR_smpl_V3V4_epwrap
    chn%V4_cheb_base_n = 5
    chn%V4_cheb_bnd = 1.0E-6_NER

  end subroutine SetUp_pope_del_SFR_smpl_test

  subroutine SetUp_pope_SFR_smpl_test(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => OPE_epwrap
!     chn%V1_cheb_base_n = 2
    chn%V1_cheb_base_n = 2
!     chn%V1_cheb_bnd = 1.0E-3_NER
    chn%V1_cheb_bnd = 0.05_NER

!     chn%V2_cheb_base_n = 5
!     chn%V2_cheb_bnd = 1.0E-6_NER

    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 4, 10/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
      chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 3, 6/)
    end select

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V3 => VN3LO_pope_SFR_smpl_V3V4_epwrap
    chn%V3_cheb_base_n = 5
    chn%V3_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V4 => VN4LO_pope_SFR_smpl_V3V4_epwrap
    chn%V4_cheb_base_n = 5
    chn%V4_cheb_bnd = 1.0E-6_NER

  end subroutine SetUp_pope_SFR_smpl_test

  subroutine SetUp_pope_DR_smpl_test(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => OPE_epwrap
!     chn%V1_cheb_base_n = 2
    chn%V1_cheb_base_n = 2
!     chn%V1_cheb_bnd = 1.0E-3_NER
    chn%V1_cheb_bnd = 0.05_NER

!     chn%V2_cheb_base_n = 5
!     chn%V2_cheb_bnd = 1.0E-6_NER

    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 4, 10/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
      chn%pVfunc_V2 => VN2LO_smplper_epwrap
        chn%npar(1:5) = (/0, 0, 1, 3, 6/)
    end select

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V3 => VN3LO_pope_DR_smpl_V3V4_epwrap
    chn%V3_cheb_base_n = 5
    chn%V3_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V4 => VN4LO_pope_DR_smpl_V3V4_epwrap
    chn%V4_cheb_base_n = 5
    chn%V4_cheb_bnd = 1.0E-6_NER

  end subroutine SetUp_pope_DR_smpl_test
end module mod_vfunc_smpl



