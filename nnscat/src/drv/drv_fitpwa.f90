! Bingwei Long Oct 06, 2018
! Bingwei Long Jun 15, 2018
! Bingwei Long Jan 27, 2018

module drv_fitpwa
! Driver subroutines used to obtain partial-wave phase shifts

  ! use mod_obj_dibcpld
  use mod_obj_1s0SYLV
  use mod_obj_1s0LY
  use mod_obj_lbwdib
  use mod_obj_sprbl1s0
  use mod_obj_cpld
  use mod_obj_withc0_sngl
  use mod_obj_pope_sngl
  use mod_obj_suppope_sngl
  use mod_obj_pope_cpld
  use mod_obj_suppope_cpld

  use mod_vfunc

  use mod_obj_smplchn
  use mod_vfunc_smpl

  use drv_pwphase,    only :                                    &
    & printout_basic_phycnst, write_to_file_basic_phycnst,      &
    & get_phase_as_value_for_klist_smpl, get_phase_Tmtrx_as_value_for_klist_smpl, &
    & get_phase_as_value_for_klist_forBD3s1_smpl, get_phase_as_value_for_klist_forSCL1s0_smpl
  use eft_potspwd
  use eft_phaseconv
  use util_io
  implicit none

  ! uptoQn is the order:
  !   Leading order (LO) -> uptoQn = 0
  !   Next-leading order (NLO) -> uptoQn = 1, and so on.
  ! phaselst(ii, 1) is LO for k = klist(ii), phaselst(ii, 2) is NLO correction, etc.
  abstract interface
    subroutine GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
      & num_paras, klist, numk, phaselst)
      import :: NER
      real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
      integer, intent(in)         :: Nmesh, uptoQn, chnid, num_paras, numk
      real(NER), intent(out)      :: phaselst(:, :)
    end subroutine GetSngl_fitpwa
    subroutine GetCpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
      & num_paras, klist, numk, phaselst)
      import  :: NER, triplet_phs_cpld
      real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
      integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
      type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)
    end subroutine GetCpld_fitpwa
    subroutine Get3s1BD_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
      & num_paras, klist, numk, phaselst, paralstout)
      import  :: NER, triplet_phs_cpld
      real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
      integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
      type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)
      real(NER), intent(out)    :: paralstout(:)
    end subroutine Get3s1BD_fitpwa
    subroutine GetSLBEPS_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
      & num_paras, klist, numk, phaselst, Tmtrx)
      import  :: NER
      real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
      integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
      real(NER), intent(out)      :: phaselst(:, :)
      real(NER), intent(out)      :: Tmtrx(:, :)
    end subroutine GetSLBEPS_fitpwa
    subroutine GetSnglTmtrx_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
      & num_paras, klist, numk, Tmtrx)
      import :: NER
      real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
      integer, intent(in)         :: Nmesh, uptoQn, chnid, num_paras, numk
      real(NER), intent(out)      :: Tmtrx(:, :)
    end subroutine GetSnglTmtrx_fitpwa
  end interface

contains

!\\\\---- 1S0 ----////
  ! Compute 1S0 phase shifts using obj_lbwdib
  ! uptoQn = 0: num_paras = 2
  ! uptoQn = 1: num_paras = 5
  ! uptoQn = 2: num_paras = 9
  ! Output: phaselst(ii, n) = phase at klist(ii) at the nth order
  subroutine Get1s0dib_Long_fitpwa(Nmesh, Mambda, uptoQn, paralst,             &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype, chnid
    class(obj_lbwdib), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    chnid = CHNID_1S0_PP

    allocate(obj_lbwdib::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda, klist,   &
      & numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Get1s0dib_Long_fitpwa

  ! Compute 1S0 phase shifts using obj_lbwdib
  ! Slightly different interface than Get1s0dib_Long_fitpwa
  ! uptoQn = 0: num_paras = 2
  ! uptoQn = 1: num_paras = 5
  ! uptoQn = 2: num_paras = 9
  ! uptoQn = 3: num_paras = 14
  subroutine Get1s0dib_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_lbwdib), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_lbwdib::ptrchn)
    call create_channel(ptrchn, CHNID_1S0_PP, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda, klist,   &
      & numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Get1s0dib_fitpwa

  ! Compute 1S0 phase shifts using obj_sprbl1s0
  ! uptoQn = 0: num_paras = 2
  subroutine GetSprbl1s0_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,         &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_sprbl1s0), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_sprbl1s0::ptrchn)
    call create_channel(ptrchn, CHNID_1S0_PP, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda, klist,   &
      & numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetSprbl1s0_fitpwa

  ! Compute 1S0 phase shifts using obj_smplchn
  ! uptoQn = 0: num_paras = 2
  ! uptoQn = 1: num_paras = 5
  ! uptoQn = 2: num_paras = 9
  ! uptoQn = 3: num_paras = 14
  subroutine Getspr1s0_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUpVfunc_sprbl_smpl(ptrchn)

    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Getspr1s0_smpl_fitpwa

  ! Compute LO 1S0 phase shifts using obj_1s0LY
  ! num_paras = 1
  subroutine get1s0LYLO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist, numk, &
    & phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_1s0LY), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_1S0_PP

    allocate(obj_1s0LY::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get1s0LYLO_fitpwa

  ! Compute LO 1S0 phase shifts using obj_1s0SYLV
  ! num_paras = 3
  subroutine get1s0SYLVLO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist,     &
    & numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_1s0SYLV), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_1S0_PP

    allocate(obj_1s0SYLV::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    ! Warning: remember to release the object pointer by doing the following
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get1s0SYLVLO_fitpwa

  ! Compute LO 1S0 phase shifts using obj_1s0LY, with LO being pure Yukawa
  subroutine get1s0YukawaLO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist,   &
    & numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    real(NER)                   :: PC_CgA, modparalst(1:5)
    class(obj_1s0LY), pointer   :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_1S0_PP

    allocate(obj_1s0LY::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    PC_CgA = PC_gA*PC_gA*PC_mN/(8.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi)
    modparalst(1) = paralst(1) - PC_CgA
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get1s0YukawaLO_fitpwa

  !\\\\---- Generic sub for Non-Perturbative pion in coupled-channels----////
  ! Applicable for 3S1-3D1, 3P2-3F2, etc., and uptoQn <= 3
  ! uptoQn = 0: num_paras = 1
  ! uptoQn = 1: NLO correction vanishes
  ! uptoQn = 2: num_paras = 4
  ! uptoQn = 3: num_paras = 7
  subroutine GetNPCpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetNPCpld_fitpwa

  ! subroutine GetNPCpld_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
  !   & num_paras, klist, numk, phaselst)

  !   real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
  !   integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
  !   type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

  !   integer :: regtype
  !   class(obj_smplchn), pointer :: ptrchn

  !   regtype = REGTYPE_GAUSSIAN
  !   allocate(obj_smplchn::ptrchn)
  !   call create_smplchn(ptrchn, chnid, uptoQn)
  !   ptrchn%vanishing(1) = .true.
  !   ptrchn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL

  !   ptrchn%pVfunc_V0 => Vfunc_VLO_withc0_chn

  !   ptrchn%pVfunc_V2 => Vfunc_VN2LO_withc0_chn
  !   ptrchn%V2_cheb_base_n = 4
  !   ptrchn%V2_cheb_bnd = 1.0E-3_NER

  !   ptrchn%pVfunc_V3 => Vfunc_VN3LO_withc0_chn
  !   ptrchn%V3_cheb_base_n = 4
  !   ptrchn%V3_cheb_bnd = 1.0E-3_NER

  !   call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
  !     & klist, numk, paralst, num_paras, phaselst)

  !   call ptrchn%erase()
  !   deallocate(ptrchn)

  ! end subroutine GetNPCpld_smpl_fitpwa

  !\\\\---- Generic sub for Non-Perturbative pion in uncoupled channels----////
  ! Without C0
  ! Applicable for 1P1, 3P1; 1D2, and uptoQn <= 3
  ! * For 1P1 and 3P1
  !   - uptoQn = 0: num_paras = 0
  !   - uptoQn = 1: NLO correction vanishes
  !   - uptoQn = 2: num_paras = 1
  !   - uptoQn = 3: num_paras = 2
  ! * For 1D2, num_paras = 0 for all orders

  subroutine GetNPWOC0Sngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,       &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer :: regtype
    class(obj_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetNPWOC0Sngl_fitpwa

  !\\\\---- Generic sub for Non-Perturbative pion in uncoupled channels----////
  ! With C0
  ! Applicable for 3P0; 3D2
  ! uptoQn = 0: num_paras = 1
  ! uptoQn = 1: NLO correction vanishes
  ! uptoQn = 2: num_paras = 3
  ! uptoQn = 3: num_paras = 5

  subroutine GetNPWithC0Sngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,     &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer :: regtype
    class(obj_withc0_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_withc0_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetNPWithC0Sngl_fitpwa

  subroutine GetNPWithC0Sngl_del_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,     &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer :: regtype
    class(obj_withc0_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_withc0_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_VDelTPE0_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VDelTPE1_sngl
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetNPWithC0Sngl_del_fitpwa

! SFR deltaful

  subroutine GetNPWithC0Sngl_SFRdel_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,     &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer :: regtype
    class(obj_withc0_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_withc0_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_VDelTPE0sfr_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VDelTPE1sfr_sngl
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetNPWithC0Sngl_SFRdel_fitpwa

  !\\\\---- OPE-ladder for coupled-channels----////
  ! cpld phase shifts with only OPE ladders, (uptoQn)-pion-exchanges
  ! No free parameters are asked as inputs;
  ! chnid is pre-defined ID integer number for a given channel (see eft_potspwd).
  ! Eg. CHNID_3P2_PP, CHNID_3D3_PP, etc.
  subroutine GetOPELadderCpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, klist,      &
    & numk, phaselst)

    real(NER), intent(in)       :: Mambda, klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, chnid, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype, ii
    class(obj_pope_cpld), pointer :: ptrchn
    complex(NEC)  :: ladderTQn(2, 2, MAX_NUMBER_ORDER)

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_pope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%init(Nmesh, regtype, Mambda)
    do ii = 1, numk
      call ptrchn%update_k(klist(ii))
      call ptrchn%get_opeladder_pope(ladderTQn)
      ptrchn%TQn(1:2, 1:2, 1:uptoQn+1) = ladderTQn(1:2, 1:2, 1:uptoQn+1)
      ptrchn%TQn_updated = .true.
      call ptrchn%convert_TQn_to_phase_quiet()
      call ptrchn%output_incremental_phase(phaselst(ii, 1:ptrchn%uptoQn+1))
    end do
    call ptrchn%release()
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetOPELadderCpld_fitpwa

  !!! WARNING:
  ! * The subs for the original manuscript (Unsuppressed scenario) are renamed as
  !   GetPOPECpld_fitpwa_premod and GetPOPESngl_fitpwa_premod.
  ! * The modified subs, however, inherited the previous names.
  ! * While the calling interface, e.g., type of input variables, stay the same
  !   as before, the number of fitting paras at each order has changed. See
  !   below for more detail.
  ! * The modified sub can now do uptoQn = 4
  ! %%%% Suppressed-scenario (or natural) Perturbative OPE %%%%
  ! For coupled channels
  ! uptoQn <= 4
  ! For 3P2-3F2:
  !   * When uptoQn <= 1, num_paras = 0
  !   * When uptoQn = 2, num_paras = 1 with paralst(1)
  !   * When uptoQn = 3, num_paras = 2 with paralst(1) and (2).
  !     - paralst(1) should be the same paralst(1) as for uptoQn = 2.
  !   * When uptoQn = 4, num_paras = 5 with paralst(1) and (2).
  !     - paralst(1) should be the same paralst(1) as for uptoQn = 2.
  !     - paralst(2) should be the same paralst(1) as for uptoQn = 3.
  !     - Therefore, there are 3 free physics parameters, instead of 5.
  !     - paralst(5) affects only the mixing angle E2.
  ! For 3D3-3G3:
  !   * When uptoQn <= 3, num_paras = 0
  !   * When uptoQn = 4, num_paras = 1 with paralst(1)
  ! For 3F4-3H4, 3G5-3I5, and higher waves:
  !   * num_paras = 0
  subroutine GetPOPECpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_SUPpope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPOPECpld_fitpwa

  ! Same as GetPOPECpld_fitpwa, kept as legacy
  subroutine GetAltPOPECpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    call GetPOPECpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
      & num_paras, klist, numk, phaselst)

  end subroutine GetAltPOPECpld_fitpwa

  !\\\\---- Generic peripheral PW for uncoupled-channels----////
  ! uptoQn <= 4
  ! For 3P0, 1P1, and 3P1:
  !   * When uptoQn <= 1, num_paras = 0
  !   * When uptoQn = 2, num_paras = 1 with paralst(1)
  !   * When uptoQn = 3, num_paras = 2 with paralst(1) and (2). However,
  !     paralst(1) should be the same paralst(1) as for uptoQn = 2. Ie., there
  !     are two free physics parameter rather than three.
  !   * When uptoQn = 4, num_paras = 4 with paralst(1), (2), (3) and (4).
  ! For 1D2 and 3D2:
  !   * When uptoQn <= 3, num_paras = 0
  !   * When uptoQn = 4, num_paras = 1 with paralst(1)
  ! For 1F3, 3F3, 1G4, 3G4 and higher waves (^1L_L and ^3L_L with L>=3):
  !   * num_paras = 0
  subroutine GetPOPESngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,         &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPOPESngl_fitpwa

 subroutine GetPOPEDeltaSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,         &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)

      ptrchn%pVfunc_VL2 => Vfunc_VDelTPE0_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VDelTPE1_sngl

    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPOPEDeltaSngl_fitpwa

  !\\\\---- altope for coupled-channels----////
  ! uptoQn <= 3
  ! For 3P2-3F2:
  !   * When uptoQn <= 1, num_paras = 0
  !   * When uptoQn = 2, num_paras = 1 with paralst(1)
  !   * When uptoQn = 3, num_paras = 4 with paralst(1), (2), (3) and (4).
  !      - paralst(4) affects only the mixing angle E2.
  !      - paralst(1) should be the same paralst(1) as for uptoQn = 2.
  !        Ie., there are three free physics parameters, instead of four.
  ! For 3D3-3G3:
  !   * When uptoQn <= 2, num_paras = 0
  !   * When uptoQn = 3, num_paras = 1 with paralst(1)
  ! For 3F4-3H4, 3G5-3I5, and higher waves:
  !   * num_paras = 0
  subroutine GetPOPECpld_fitpwa_premod(Nmesh, Mambda, chnid, uptoQn,           &
    & paralst, num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_pope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_pope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda, klist,   &
      & numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPOPECpld_fitpwa_premod

  !\\\\---- Generic peripheral PW for uncoupled-channels----////
  ! uptoQn <= 3
  ! For 3P0, 1P1, and 3P1:
  !   * When uptoQn <= 1, num_paras = 0
  !   * When uptoQn = 2, num_paras = 1 with paralst(1)
  !   * When uptoQn = 3, num_paras = 3 with paralst(1), (2) and (3). However,
  !     paralst(1) should be the same paralst(1) as for uptoQn = 2. Ie., there
  !     are two free physics parameter rather than three.
  ! For 1D2 and 3D2:
  !   * When uptoQn <= 2, num_paras = 0
  !   * When uptoQn = 3, num_paras = 1 with paralst(1)
  ! For 1F3, 3F3, 1G4, 3G4 and higher waves (^1L_L and ^3L_L with L>=3):
  !   * num_paras = 0
  subroutine GetPOPESngl_fitpwa_premod(Nmesh, Mambda, chnid, uptoQn, paralst,  &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_pope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_pope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPOPESngl_fitpwa_premod

!\\\\---- 3S1-3D1 ----////
  subroutine get3s1LO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist, numk,   &
    & phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    ! See eft_phaseconv for definition of type(triplet_phs_cpld)
    type(triplet_phs_cpld), intent(out)  :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_3S1_PP

    allocate(obj_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3s1LO_fitpwa

  ! Bookmark 03/19/2018; get3s1res_fitpwa to be worked on
  ! subroutine get3s1res_fitpwa(Nmesh, Mambda, paralst, uptoQn, klist, numk, phaselst)

  !     real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
  !     integer, intent(in)         :: Nmesh, uptoQn, numk
  !     real(NER), intent(out)      :: phaselst(:)

  !     integer :: regtype, chnid, num_paras
  !     class(obj_cpld), pointer :: ptrchn
  !     integer :: ii
  !     integer, dimension(1:10)  :: uflag

  !     regtype = REGTYPE_GAUSSIAN
  !     chnid = CHNID_3S1_PP
  !     num_paras = 2

  !     allocate(obj_withc0_sngl::ptrchn)
  !     call create_channel(ptrchn, chnid, uptoQn)
  !     call ptrchn%load_inputs(num_paras, paralst)
  !     call ptrchn%init(Nmesh, regtype, Mambda, Mambda*RLOM)
  !     do ii = 1, numk
  !         call ptrchn%update_k(klist(ii))
  !         call ptrchn%get_residual_TQn()
  !         call ptrchn%convert_TQn_to_phase_quiet(uflag)
  !         call ptrchn%output_total_phase(phaselst(ii))
  !     end do
  !     call ptrchn%release()
  !     call ptrchn%erase()
  !     deallocate(ptrchn)

  ! end subroutine get3s1res_fitpwa

!\\\\---- 1P1 ----////
  ! Compute 1P1 phase shifts with C0 at LO
  subroutine get1p1withC0_fitpwa(Nmesh, Mambda, paralst, num_paras, klist,     &
    & numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_withc0_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_1P1_PP

    allocate(obj_withc0_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get1p1withC0_fitpwa

  ! Compute 1P1 phase shifts w/o C0 at LO
  subroutine get1p1woC0_fitpwa(Nmesh, Mambda, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, klist(:)
    integer, intent(in)         :: Nmesh, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_sngl), pointer :: ptrchn
    real(NER)   :: paralst(1:2)
    integer     :: num_paras

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_1P1_PP

    allocate(obj_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get1p1woC0_fitpwa

!\\\\---- 3P0 ----////
  subroutine get3p0LO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist, numk,   &
    & phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_withc0_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    chnid = CHNID_3P0_PP
    uptoQn = 0

    allocate(obj_withc0_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3p0LO_fitpwa

  ! 3P0 phase shifts by get_residual_TQn
  ! paralst must have at least two elements
  subroutine get3p0res_fitpwa(Nmesh, Mambda, paralst, uptoQn, klist, numk,     &
    & phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, num_paras
    class(obj_withc0_sngl), pointer :: ptrchn
    integer :: ii
    ! integer, dimension(1:10)  :: uflag

    regtype = REGTYPE_GAUSSIAN
    chnid = CHNID_3P0_PP
    num_paras = 2

    allocate(obj_withc0_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%load_inputs(num_paras, paralst)
    call ptrchn%init(Nmesh, regtype, Mambda)
    do ii = 1, numk
      call ptrchn%update_k(klist(ii))
      call ptrchn%get_residual_TQn()
      call ptrchn%convert_TQn_to_phase_quiet()
      call ptrchn%output_total_phase(phaselst(ii))
    end do
    call ptrchn%release()
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3p0res_fitpwa

!\\\\---- 3P1 ----////
  ! Compute 3P1 phase shifts with C0 at LO
  subroutine get3p1withC0_fitpwa(Nmesh, Mambda, paralst, num_paras, klist,     &
    & numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_withc0_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_3P1_PP

    allocate(obj_withc0_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3p1withC0_fitpwa

  ! Compute 3P1 phase shifts w/o C0 at LO
  subroutine get3p1woC0_fitpwa(Nmesh, Mambda, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, klist(:)
    integer, intent(in)         :: Nmesh, numk
    real(NER), intent(out)      :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_sngl), pointer :: ptrchn
    real(NER)   :: paralst(1:2)
    integer     :: num_paras

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_3P1_PP

    allocate(obj_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3p1woC0_fitpwa

!\\\\---- 3P2-3F2 ----////
  subroutine get3p2LO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist, numk,   &
    & phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, num_paras, numk
    type(triplet_phs_cpld), intent(out)  :: phaselst(:)

    integer :: regtype, chnid, uptoQn
    class(obj_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    uptoQn = 0
    chnid = CHNID_3P2_PP

    allocate(obj_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
      & paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3p2LO_fitpwa

  ! Frozen
  ! num_paras = 2
  ! 3P2 dibaryon, energy-dependent C
  ! subroutine getdib3p2LO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist,      &
  !   & numk, phaselst)

  !   real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
  !   integer, intent(in)         :: Nmesh, num_paras, numk
  !   type(triplet_phs_cpld), intent(out)  :: phaselst(:)

  !   integer :: regtype, chnid, uptoQn
  !   class(obj_dibcpld), pointer :: ptrchn

  !   regtype = REGTYPE_GAUSSIAN
  !   uptoQn = 0
  !   chnid = CHNID_3P2_PP

  !   allocate(obj_dibcpld::ptrchn)
  !   call create_channel(ptrchn, chnid, uptoQn)
  !   call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
  !     & paralst, num_paras, phaselst)
  !   call ptrchn%erase()
  !   deallocate(ptrchn)

  ! end subroutine getdib3p2LO_fitpwa

  ! Frozen
  ! Next-to-leading order phase shifts of 3P2
  ! num_paras = 5
  ! where paralst(1) and paralst(2) = the two paras used to call getdib3p2LO_fitpwa
  ! Therefore, only paralst(3, 4, 5) are actually free parameters.
  ! paralst(1) & (2) give sort of "reference point"
  ! subroutine getdib3p2_NLO_fitpwa(Nmesh, Mambda, paralst, num_paras, klist,    &
  !   & numk, phaselst)

  !   real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
  !   integer, intent(in)         :: Nmesh, num_paras, numk
  !   type(triplet_phs_cpld), intent(out)  :: phaselst(:)

  !   integer :: regtype, chnid, uptoQn
  !   class(obj_dibcpld), pointer :: ptrchn

  !   regtype = REGTYPE_GAUSSIAN
  !   uptoQn = 1
  !   chnid = CHNID_3P2_PP

  !   allocate(obj_dibcpld::ptrchn)
  !   call create_channel(ptrchn, chnid, uptoQn)
  !   call ptrchn%get_phase_for_klist(Nmesh, regtype, Mambda, klist, numk,       &
  !     & paralst, num_paras, phaselst)
  !   call ptrchn%erase()
  !   deallocate(ptrchn)

  ! end subroutine getdib3p2_NLO_fitpwa

!----\\\\ Perturbative Pion ////----
! LO is always zero with the following subroutine.

  ! 3P2-3F2 phase shifts by obj_pope_cpld
  ! When uptoQn <= 2, num_paras = 0 and paralst is empty
  ! uptoQn = 3, num_paras = 1, let paralst(1) be C1
  ! uptoQn = 4, num_paras = 2; however, there is only free parameter.
  ! paralst(1) = C1 which is always the parameter used for uptoQn=3, so the
  ! actual free parameter for uptoQn = 4 is paralst(2) = C2. Stated differently,
  ! when uptoQn = 4, we fit to only ONE phase shift.
  subroutine get3p2_pope_fitpwa(Nmesh, Mambda, uptoQn, paralst, num_paras,     &
    & klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype, chnid
    class(obj_pope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    chnid = CHNID_3P2_PP

    allocate(obj_pope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3p2_pope_fitpwa

!\\\\---- 3D3-3G3 ----////
  ! 3D3-3G3 phase shifts by obj_pope_cpld
  ! When uptoQn <= 4, num_paras = 0 and paralst is empty
  subroutine get3d3_pope_fitpwa(Nmesh, Mambda, uptoQn, paralst, num_paras,     &
    & klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype, chnid
    class(obj_pope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    chnid = CHNID_3D3_PP

    allocate(obj_pope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3d3_pope_fitpwa

!\\\\---- 3D2 ----////
  ! 3D2 phase shifts by obj_pope_sngl
  ! When uptoQn <= 4, num_paras = 0 and paralst is empty
  subroutine get3d2_pope_fitpwa(Nmesh, Mambda, uptoQn, paralst, num_paras,     &
    & klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype, chnid
    class(obj_pope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    chnid = CHNID_3D2_PP

    allocate(obj_pope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3d2_pope_fitpwa

!\\\\---- 3F3 ----////
  ! 3F3 phase shifts by obj_pope_sngl
  ! When uptoQn <= 6, num_paras = 0 and paralst is empty

  subroutine get3f3_pope_fitpwa(Nmesh, Mambda, uptoQn, paralst, num_paras,     &
    & klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype, chnid
    class(obj_pope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    chnid = CHNID_3F3_PP

    allocate(obj_pope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,          &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine get3f3_pope_fitpwa

! \\\\for deltaSFR in smpl////

  subroutine Getsfrdeltacpld_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUp_pope_del_SFR_smpl(ptrchn)

    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Getsfrdeltacpld_smpl_fitpwa

  subroutine Getsfrdeltasngl_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUp_pope_del_SFR_smpl(ptrchn)
    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Getsfrdeltasngl_smpl_fitpwa

  subroutine Get1s0_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, pot, amplitude)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER)                   :: phaselst(numk,uptoQn+1)
    real(NER), intent(out)      :: amplitude(:, :)
    integer :: regtype
    character(len=20), intent(in)       :: pot
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUpVfunc_mmwly_chengdu_smpl(ptrchn, pot)
    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    amplitude(1,1:uptoQn+1) = real(ptrchn%TQn(1,1,1:uptoQn+1))
!     print *, phaselst(1, uptoQn+1), sum(ptrchn%TQn(1,1,1:uptoQn+1))
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Get1s0_smpl_fitpwa

  subroutine Get_SLBEPS_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
    & num_paras, klist, numk, phaselst, Tmtrx)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)
    real(NER), intent(out)  :: Tmtrx(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUpVfunc_mmwly_smpl(ptrchn)

    call get_phase_Tmtrx_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst, Tmtrx)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Get_SLBEPS_fitpwa

! \\\\for 3s1 rest para with binding energy in smpl////

  subroutine Get3s1_withbd_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst, paralstout)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)
    real(NER),intent(out)                   :: paralstout(:)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUpVfunc_mmwly_smpl(ptrchn)

    call get_phase_as_value_for_klist_forBD3s1_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst, paralstout)

    call ptrchn%fit_erase()
    deallocate(ptrchn)

  end subroutine Get3s1_withbd_smpl_fitpwa

! \\\\for 1s0 rest para with scattering length in smpl////

  subroutine Get1s0_withscl_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUpVfunc_mmwly_smpl(ptrchn)

    call get_phase_as_value_for_klist_forSCL1s0_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    call ptrchn%fit_erase()
    deallocate(ptrchn)

  end subroutine Get1s0_withscl_smpl_fitpwa

end module drv_fitpwa
