! pray  08/16/2019
! Bingwei Long  01/25/2019
! long-range potentia tree:
!   V_long-rang
!       |
!       + -> V_ope
!       |
!       + -> V_tpe
!             |
!             + -> V_tpedeltaless
!             |        |
!             |        + -> V_NLOdeltaless
!             |        |
!             |        + -> V_NNLOdeltaless

module potwrap

  use eft_potspwd
  use ope_pwd
  use vtpe0_pwd
  use vtpe1_pwd
  use mod_deltafullNLO_pwd,  only:TPEdel_jpm,TPEdel_jmm,TPEdel_jpp,TPEdel_j0j,TPEdel_j1j
  use mod_deltafullN2LO_pwd,  only:TPEdel_sub_jpm,TPEdel_sub_jmm,TPEdel_sub_jpp,TPEdel_sub_j0j,TPEdel_sub_j1j, &
  &Select_cis2PE_vdel
  use VTPE_SFR_pwd,  only:TPESFR_jpm,TPESFR_jmm,TPESFR_j0j,TPESFR_j1j,TPESFR_jpp
  use VTPE1_SFR_pwd, only:TPESFR1_jpm,TPESFR1_jmm,TPESFR1_j0j,TPESFR1_j1j,TPESFR1_jpp,Select_cisSFR_vtpe1,set_cisSFR_vtpe1
  implicit none

  integer, parameter, private  :: NCH = 2

  ! * WARNING: potval is always 2X2, even for uncoupled channels
  ! * For uncoupled channels, only potval(1, 1) is meaningful
  ! * Normalization
  ! If V were so weak as to validate the Born approximation, the relation of
  ! V_1S0 and the scattering length 'a' would have been
  ! < L, p | V | L, p> = 2/pi * a + ... for p -> 0
  ! * Sign
  ! For attractive V, sign of V is minus
  abstract interface
    subroutine gnrc_eft_pot(L, S, J, regtype, Mambda, para, p1, p2, potval)
      import                  :: NER, NCH
      implicit none
      integer, intent(in)     :: L, S, J, regtype
      real(NER), intent(in)   :: Mambda, p1, p2, para(:)
      real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    end subroutine gnrc_eft_pot
  end interface

contains

  subroutine zero_V(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    potval = 0.0_NER

  end subroutine zero_V

  subroutine OPE_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = OPE_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = OPE_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = OPE_jpp(j, p1, p2)
        else
          potval(1, 1) = OPE_jmm(j, p1, p2)
          potval(2, 1) = OPE_jpm(j, p1, p2)
          potval(1, 2) = OPE_jpm(j, p2, p1)
          potval(2, 2) = OPE_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine OPE_epwrap

  subroutine TPE0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = VTPE0_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = VTPE0_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = VTPE0_jpp(j, p1, p2)
        else
          potval(1, 1) = VTPE0_jmm(j, p1, p2)
          potval(2, 1) = VTPE0_jpm(j, p1, p2)
          potval(1, 2) = VTPE0_jpm(j, p2, p1)
          potval(2, 2) = VTPE0_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine TPE0_epwrap

  subroutine TPE1_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = VTPE1_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = VTPE1_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = VTPE1_jpp(j, p1, p2)
        else
          potval(1, 1) = VTPE1_jmm(j, p1, p2)
          potval(2, 1) = VTPE1_jpm(j, p1, p2)
          potval(1, 2) = VTPE1_jpm(j, p2, p1)
          potval(2, 2) = VTPE1_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine TPE1_epwrap



  ! Polynomial potential
  ! vab(1, 1) = reg(p1)*reg(p2)*(p1*p2)^L

  subroutine VPNQ0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    select case (L)
      case (0)
        potval(1, 1) = 1.0_NER
      case (1)
        potval(1, 1) = p1*p2
      case default
        potval(1, 1) = (p1*p2)**L
    end select

    potval(1, 1) = potval(1, 1)*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine VPNQ0_epwrap

  ! Polynomial potential
  ! vab(1, 1) = (p1^2 + p2^2)*reg(p1)*reg(p2)*(p1*p2)^L

  subroutine VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: pref

    potval = 0.0_NER

    select case (L)
      case (0)
        pref = 1.0_NER
      case (1)
        pref = p1*p2
      case default
        pref = (p1*p2)**L
    end select

    potval(1, 1) = (p1*p1 + p2*p2)*regltr_PP(regtype, p1, p2, Mambda)*pref

  end subroutine VPNQ2_diag_epwrap

  subroutine VPNQ2_off_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: pref

    potval = 0.0_NER

    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      select case (L)
        case (0)
          pref = 1.0_NER
        case (1)
          pref = p1*p2
        case default
          pref = (p1*p2)**L
      end select
      pref = pref * regltr_PP(regtype, p1, p2, Mambda)
      potval(1, 2) = p2*p2 * pref
      potval(2, 1) = p1*p1 * pref
    end if

  end subroutine VPNQ2_off_epwrap


! ###################################################
!   contact interaction at N5LO for suppope
!  vab(1, 1) = (p1*p1*p2*p2*)*reg(p1)*reg(p2)*p1*p2^L
! ###################################################

  subroutine VPNQ5_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: pref

    potval = 0.0_NER

    select case (L)
      case (0)
        pref = 1.0_NER
      case (1)
        pref = p1*p2
      case default
        pref = (p1*p2)**L
    end select

    potval(1, 1) = (p1*p1*p2*p2)*regltr_PP(regtype, p1, p2, Mambda)*pref

  end subroutine VPNQ5_epwrap

!#############################################
!VPNQ5_off =  (p1*p2)^L*reg(p1)*reg(p2)*|     0          p1*p1*p2*p2 |
!                                       | p1*p1*p2*p2        0       |
! ###########################################

  subroutine VPNQ5_off_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: pref

    potval = 0.0_NER

    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      select case (L)
        case (0)
          pref = 1.0_NER
        case (1)
          pref = p1*p2
        case default
          pref = (p1*p2)**L
      end select
      pref = pref * regltr_PP(regtype, p1, p2, Mambda)
      potval(1, 2) = p2*p2*p1*p1 * pref
      potval(2, 1) = potval(1, 2)
    end if

  end subroutine VPNQ5_off_epwrap
!#####################################################
! VPNQ5_diag = (p1*p2)^L*reg(p1)*reg(p2)*|         0                  0      |
!                                        |         0            p1*p1*p2*p2  |

  subroutine VPNQ5_diag_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: pref

    potval = 0.0_NER

    select case (L)
      case (0)
        pref = 1.0_NER
      case (1)
        pref = p1*p2
      case default
        pref = (p1*p2)**L
    end select

    potval(2, 2) = (p1*p1*p2*p2)*regltr_PP(regtype, p1, p2, Mambda)*pref

  end subroutine VPNQ5_diag_epwrap

  ! Applicable to perturbative channel with deltaless potential;
  subroutine VN3LO_pope_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH), vs(1:NCH, 1:NCH)

    potval = 0.0_NER

    call TPE0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe)
    select case (L)
      case(1)
        call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
        vs(1,1) = Cs(2) * vs(1,1)
        potval = vtpe + vs
      case default
        potval = vtpe
    end select

  end subroutine VN3LO_pope_smpl_epwrap

  subroutine VN4LO_pope_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH), vs(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    real(NER) :: vq(1:NCH, 1:NCH)
    potval = 0.0_NER

    call TPE1_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vq)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    select case (L)
      case (1)
        if(triqn_is_coupled(L, S, J))  then  ! 3p2-3f2
          call VPNQ2_off_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vl)
          potval = vtpe + Cs(3) * vq + Cs(4) * vs + Cs(5) * vl
        else
          potval(1, 1) = vtpe(1, 1) + Cs(3) * vq(1, 1) + Cs(4) * vs(1, 1)
        end if
      case (2)
        potval = vtpe + Cs(1) * vq
      case default
        potval = vtpe
    end select

  end subroutine VN4LO_pope_smpl_epwrap

  ! Perturbative channels for both deltaless and deltaful
  subroutine VN2LO_smplper_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vs(1:NCH, 1:NCH)

    potval = 0.0_NER

    select case (L)
      case(1)
        call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
        vs(1,1) = Cs(1) * vs(1,1)
        potval = vs
      case default
        potval = 0.0_NER
    end select
  end subroutine VN2LO_smplper_epwrap

  ! Applicable to 3S1 - 3D1;

  subroutine VLO_withc0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vope(1:NCH, 1:NCH), vs(1:NCH, 1:NCH)

    potval = 0.0_NER

    call OPE_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vope)
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    vs(1, 1) = Cs(1) * vs(1, 1)

    ! debug
    ! vope = 0.0_NER

    potval = vope + vs

  end subroutine VLO_withc0_epwrap

  ! Contact of VN2LO_withc0

  subroutine VN2LO_withc0_SHRG_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vs0, vs2, vs2off

    potval = 0.0_NER

    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs0)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2)

    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      call VPNQ2_off_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2off)
      potval = Cs(2)*vs0 + Cs(3)*vs2 + Cs(4)*vs2off
    else
      potval = Cs(2)*vs0 + Cs(3)*vs2
    end if

  end subroutine VN2LO_withc0_SHRG_epwrap

  ! Contact of VN3LO_withc0

  subroutine VN3LO_withc0_SHRG_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vs0, vs2, vs2off

    potval = 0.0_NER

    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs0)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2)

    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      call VPNQ2_off_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2off)
      potval = Cs(5)*vs0 + Cs(6)*vs2 + Cs(7)*vs2off
    else
      potval = Cs(4)*vs0 + Cs(5)*vs2
    end if

  end subroutine VN3LO_withc0_SHRG_epwrap
! contact test terms


  subroutine VN2LO_withc0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPE0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vl)
    call VN2LO_withc0_SHRG_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN2LO_withc0_epwrap

  subroutine VN3LO_withc0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPE1_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vl)
    call VN3LO_withc0_SHRG_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN3LO_withc0_epwrap

!==================================================================
! Minimally modified Weinberg 1S0 (Long & Yang PRC '11 - '12)
! VLO = OPE + C0 for 1S0
! VNLO = C1 + D0*(p^2 + p'^2) for 1S0
! VN2LO = TPE0 + C2 + D1*(p^2 + p'^2) + E0*p^2*p'^2 for 1S0
! =================================================================

  subroutine VLO_MMWLY_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    call VLO_withc0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

  end subroutine VLO_MMWLY_epwrap

  subroutine VNLO_MMWLY_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vs(1:NCH, 1:NCH), voff(1:NCH,1:NCH)

    call VPNQ0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vs)
    potval = para(2)*potval + para(3)*vs

  end subroutine VNLO_MMWLY_epwrap

  subroutine VN2LO_MMWLY_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vs(1:NCH, 1:NCH), v5(1:NCH,1:NCH), vtpe0(1:NCH,1:NCH)

    call VPNQ0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vs)
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, para, p1, p2, v5)
    call TPE0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vtpe0)

    potval = para(4)*potval + para(5)*vs + para(6)*v5 + vtpe0

  end subroutine VN2LO_MMWLY_epwrap

  subroutine VN2LO_MMWLY_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vs(1:NCH, 1:NCH), v5(1:NCH,1:NCH)

    call VPNQ0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vs)
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, para, p1, p2, v5)
    potval = para(4)*potval + para(5)*vs + para(6)*v5

  end subroutine VN2LO_MMWLY_SHRG_epwrap

!============================================================
!
! "Short-range" part of 1S0 separable pots
!
!============================================================

  subroutine VLO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: kps, vope(1:NCH, 1:NCH)
    real(NER) :: CgA

    potval = 0.0_NER
    CgA = PC_gA*PC_gA*PC_mN/(8.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi)

    kps = PC_mN*para(1)
    potval(1, 1) = PC_mN*para(2)/(sqrt(p1*p1 + kps)*sqrt(p2*p2 + kps))
    potval(1, 1) = (potval(1, 1) - CgA) * regltr_PP(regtype, p1, p2, Mambda)

  end subroutine VLO_sprbl_SHRG_epwrap

  subroutine VNLO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: kps, y0, Delta1, y1, denom1, denom2, C0

    potval = 0.0_NER

    ! para(1)=Delta, para(2)=y, para(3)=Delta1, para(4)=y1, para(5)=C0
    kps = PC_mN*para(1)
    y0 = para(2)
    denom1 = p1*p1 + kps
    denom2 = p2*p2 + kps
    Delta1 = para(3)
    y1 = para(4)
    C0 = para(5)

    potval(1, 1) = PC_mN/(sqrt(denom1)*sqrt(denom2))      &
      & *(y1 - PC_mN*y0/2.0_NER*(1.0_NER/denom1 + 1.0_NER/denom2)*Delta1) + C0

    potval(1, 1) = potval(1, 1)*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine VNLO_sprbl_SHRG_epwrap

  subroutine VN2LO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vtpe0
    real(NER) :: kps, y0, Delta1, y1, Delta2, y2, denom1, denom2, C1, D0

    potval = 0.0_NER

    ! para(1)=Delta, para(2)=y,
    ! para(3)=Delta1, para(4)=y1, para(5)=C0
    ! para(6)=Delta2, para(7)=y2, para(8)=C1, para(9)=D0
    kps = PC_mN*para(1)
    y0 = para(2)
    denom1 = p1*p1 + kps
    denom2 = p2*p2 + kps
    Delta1 = para(3)
    y1 = para(4)
    Delta2 = para(6)
    y2 = para(7)
    C1 = para(8)
    D0 = para(9)

    potval(1, 1) = PC_mN/(sqrt(denom1)*sqrt(denom2))*(y2 - PC_mN/2.0_NER       &
      & *(1.0_NER/denom1 + 1.0_NER/denom2)*(y0*Delta2 + y1*Delta1)             &
      & + PC_mN**2/8.0_NER*(3.0_NER/denom1**2 + 2.0_NER/(denom1*denom2)        &
      & + 3.0_NER/denom2**2)*y0*Delta1**2) + C1 + D0*(p1*p1 + p2*p2)

    potval(1, 1) = potval(1, 1)*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine VN2LO_sprbl_SHRG_epwrap

  subroutine VN3LO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vtpe1
    real(NER) :: kps, y0, Delta1, y1, Delta2, y2, Delta3, y3, denom1, denom2,  &
      & C2, D1, E0

    potval = 0.0_NER

    ! para(1)  = Delta,  para(2)  = y,
    ! para(3)  = Delta1, para(4)  = y1, para(5)   = C0
    ! para(6)  = Delta2, para(7)  = y2, para(8)   = C1, para(9)  = D0
    ! para(10) = Delta3, para(11) = y3, paras(12) = C2, para(13) = D1,
    ! para(14) = E0
    kps = PC_mN*para(1)
    y0 = para(2)
    denom1 = p1*p1 + kps
    denom2 = p2*p2 + kps
    Delta1 = para(3)
    y1 = para(4)
    Delta2 = para(6)
    y2 = para(7)
    Delta3 = para(10)
    y3 = para(11)
    C2 = para(12)
    D1 = para(13)
    E0 = para(14)

    potval(1, 1) = PC_mN/(sqrt(denom1)*sqrt(denom2))*(y3 - PC_mN/2.0_NER       &
      & *(1.0_NER/denom1 + 1.0_NER/denom2)*(y0*Delta3 + y1*Delta2 + y2*Delta1) &
      & + PC_mN**2/8.0_NER*(3.0_NER/denom1**2 + 2.0_NER/(denom1*denom2)        &
      &      + 3.0_NER/denom2**2)*(y1*Delta1**2 + 2.0_NER*y0*Delta1*Delta2)    &
      & - PC_mN**3/16.0_NER*(5.0_NER/denom1**3 + 3.0_NER/(denom1**2*denom2)    &
      &      + 3.0_NER/(denom1*denom2**2) + 5.0_NER/denom2**3)*y0*Delta1**3)   &
      & + C2 + D1*(p1*p1 + p2*p2) + E0*p1*p1*p2*p2

    potval(1, 1) = potval(1, 1)*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine VN3LO_sprbl_SHRG_epwrap

!============================================================
!
! Separable 1S0 pots plus deltaless pion-exchanges
!
!============================================================

  subroutine VLO_sprbl_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vope(1:NCH, 1:NCH)

    call VLO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)
    call OPE_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vope)
    potval = potval + vope

  end subroutine VLO_sprbl_epwrap

  subroutine VNLO_sprbl_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    call VNLO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

  end subroutine VNLO_sprbl_epwrap

  subroutine VN2LO_sprbl_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe0(1:NCH, 1:NCH)

    call TPE0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vtpe0)
    call VN2LO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)
    potval(1, 1) = potval(1, 1) + vtpe0(1, 1)

  end subroutine VN2LO_sprbl_epwrap

  subroutine VN3LO_sprbl_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe1(1:NCH, 1:NCH)

    call TPE1_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vtpe1)
    call VN3LO_sprbl_SHRG_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)
    potval(1, 1) = potval(1, 1) + vtpe1(1, 1)

  end subroutine VN3LO_sprbl_epwrap


!============================================================
!
! Relativistic corrections
!
!============================================================

  subroutine RelOPE_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = rel_OPE_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = rel_OPE_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = rel_OPE_jpp(j, p1, p2)
        else
          potval(1, 1) = rel_OPE_jmm(j, p1, p2)
          potval(2, 1) = rel_OPE_jpm(j, p1, p2)
          potval(1, 2) = rel_OPE_jpm(j, p2, p1)
          potval(2, 2) = rel_OPE_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine RelOPE_epwrap

  subroutine VN2LO_withc0_RelOPE_epwrap(L, S, J, regtype, Mambda, para, p1,    &
    & p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: tpe0, relope

    call VN2LO_withc0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, tpe0)
    call RelOPE_epwrap(L, S, J, regtype, Mambda, para, p1, p2, relope)

    potval = tpe0 + relope

  end subroutine VN2LO_withc0_RelOPE_epwrap

!===================================================================
!
! Short-range part of energy-dependent dibaryon pots
!
!===================================================================

  subroutine VLO_LBWDIB_SHRG_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, Ecm, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER)   :: CgA, C0, tmppara(1:1)

    call VPNQ0_epwrap(L, S, J, regtype, Mambda, tmppara, p1, p2, potval)
    CgA = PC_gA*PC_gA*PC_mN/(8.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi)
    ! C0 = y/(E + Delta) - PC_CgA
    ! Delta = paras(1), y = paras(2)
    C0 = para(2)/(Ecm + para(1)) - CgA
    potval = C0 * potval

  end subroutine VLO_LBWDIB_SHRG_epwrap

!===================================================================
!
! Energy-dependent dibaryon pots (with deltaless pion-exchg)
!
!===================================================================

  subroutine VLO_LBWDIB_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, Ecm, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vope(1:NCH, 1:NCH)

    call VLO_LBWDIB_SHRG_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)
    call OPE_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vope)
    potval = potval + vope

  end subroutine VLO_LBWDIB_epwrap

!===================================================================
!
! Short-range part of energy-dependent 1.5-dibaryon pots
!
!===================================================================

  subroutine VLO_SYLV_SHRG_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, Ecm, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER)   :: CgA, C0, tmppara(1:1)

    call VPNQ0_epwrap(L, S, J, regtype, Mambda, tmppara, p1, p2, potval)
    CgA = PC_gA*PC_gA*PC_mN/(8.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi)
    ! C0 = y/(E + Delta) - CgA + z
    ! Delta = para(1), y = para(2), z = para(3)
    C0 = para(2)/(Ecm + para(1)) - CgA + para(3)
    potval = C0 * potval

  end subroutine VLO_SYLV_SHRG_epwrap

!===================================================================
!
! Energy-dependent 1.5-dibaryon pots (with deltaless pion-exchg)
!
!===================================================================

  subroutine VLO_SYLV_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, Ecm, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vope(1:NCH, 1:NCH)

    call VLO_SYLV_SHRG_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)
    call OPE_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vope)
    potval = potval + vope

  end subroutine VLO_SYLV_epwrap

!===================================================================
!
! for perturbative 3P0 with C0
!
!===================================================================

  subroutine VN2LO_per3p0_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vs(1:NCH, 1:NCH), vtpe(1:NCH, 1:NCH)

    potval = 0.0_NER
    vtpe = 0.0_NER
    vs = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call TPE0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    potval = vtpe + Cs(2)*potval + Cs(3)*vs
  end subroutine VN2LO_per3p0_smpl_epwrap

  subroutine VN2LO_perwot3p0_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vs(1:NCH, 1:NCH)

    potval = 0.0_NER
    vs = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    potval =+ Cs(2)*potval + Cs(3)*vs
  end subroutine VN2LO_perwot3p0_smpl_epwrap

  subroutine VN3LO_per3p0_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vs(1:NCH, 1:NCH), vtpe(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call TPE1_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    potval = vtpe + Cs(4)*potval + Cs(5)*vs
  end subroutine VN3LO_per3p0_smpl_epwrap

  subroutine VN3LO_per3p0_high_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vs(1:NCH, 1:NCH), vtpe(1:NCH, 1:NCH), vst(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call TPE1_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vst)
    potval = vtpe + Cs(4)*potval + Cs(5)*vs + Cs(6)*vst
  end subroutine VN3LO_per3p0_high_smpl_epwrap

  subroutine VN3LO_perwot3p0_high_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vs(1:NCH, 1:NCH), vst(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vst)
    potval = + Cs(4)*potval + Cs(5)*vs + Cs(6)*vst
  end subroutine VN3LO_perwot3p0_high_smpl_epwrap

  subroutine VN2LO_supper3p0_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = Cs(2)*potval
  end subroutine VN2LO_supper3p0_smpl_epwrap

  subroutine VN3LO_supper3p0_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vs(1:NCH, 1:NCH), vtpe(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call TPE0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    potval = vtpe + Cs(3)*potval + Cs(4)*vs
  end subroutine VN3LO_supper3p0_smpl_epwrap

  subroutine VN4LO_supper3p0_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vs(1:NCH, 1:NCH), vtpe(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call TPE1_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs)
    potval = vtpe + Cs(5)*potval + Cs(6)*vs
  end subroutine VN4LO_supper3p0_smpl_epwrap

end module potwrap

