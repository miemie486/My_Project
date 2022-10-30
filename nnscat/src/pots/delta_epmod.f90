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
!             |
!             + -> V_tpedeltafull
!                       |
!                       + -> V_NLOdeltafull
!                       |
!                       + -> V_NNLOdeltafull
!  pray  08/16/2019
module delta_epmod

  use potwrap

  implicit none

  integer, parameter, private  :: NCH = 2
contains

! contact test terms

  subroutine VN3LO_withc0_deltac_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vs0, vs2, vs2off,vct5,vs5off,vs5

    potval = 0.0_NER

    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs0)
    call VPNQ5_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs5)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs2)
    call VPNQ5_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct5)

    if (L /= j .and. j /= 0) then  ! is coupled channel?
      call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs2off)
        call VPNQ5_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs5off)
      potval = Cs(5)*vs0 + Cs(6)*vs2 + Cs(7)*vs2off + Cs(8)*vs5 + Cs(9)*vct5 +&
      & Cs(10)*vs5off
    else
      potval = Cs(4)*vs0 + Cs(5)*vs2 + Cs(6)*vs5
    end if

  end subroutine VN3LO_withc0_deltac_SHRG_epwrap

  subroutine VN2LOsfr_withc0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPESFR_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    call VN2LO_withc0_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN2LOsfr_withc0_epwrap

  subroutine VN3LOsfr_withc0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPESFR1_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    call VN3LO_withc0_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN3LOsfr_withc0_epwrap

  subroutine VN2LOSFR_withc0_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPE0SFR_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    call VN2LO_withc0_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs
  end subroutine VN2LOSFR_withc0_del_epwrap

  subroutine VN3LOSFR_withc0_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPE1SFR_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    call VN3LO_withc0_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN3LOSFR_withc0_del_epwrap

  subroutine VN3LOSFR_withc0_delc_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPE1SFR_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    call VN3LO_withc0_deltac_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN3LOSFR_withc0_delc_epwrap

  subroutine TPE0SFR_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vt

    call TPESFR_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vt)
    call TPE0_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    potval = vl +  vt
  end subroutine TPE0SFR_del_epwrap

  subroutine TPE1SFR_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vt

    call TPESFR1_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vt)
    call TPE1_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    potval = vl + vt

  end subroutine TPE1SFR_del_epwrap


! ==========================================================================
!  add SFR deltafull and deltaless TPE potential
!  add some sub for Born approximation . AUG/12/2019 pray
!  TPE0_delonly: NLOdeltafull potential
!  TPE1_delonly: NNLOdeltafull potential; both of them from Krebs/Epelbuam07
!  Eur. Phys. J. A 32, 127-137(2007)
! ==========================================================================

  subroutine TPE0_delonly_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = TPEdel_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = TPEdel_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = TPEdel_jpp(j, p1, p2)
        else
          potval(1, 1) = TPEdel_jmm(j, p1, p2)
          potval(2, 1) = TPEdel_jpm(j, p1, p2)
          potval(1, 2) = TPEdel_jpm(j, p2, p1)
          potval(2, 2) = TPEdel_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine TPE0_delonly_epwrap

  subroutine TPE1_delonly_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = TPEdel_sub_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = TPEdel_sub_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = TPEdel_sub_jpp(j, p1, p2)
        else
          potval(1, 1) = TPEdel_sub_jmm(j, p1, p2)
          potval(2, 1) = TPEdel_sub_jpm(j, p1, p2)
          potval(1, 2) = TPEdel_sub_jpm(j, p2, p1)
          potval(2, 2) = TPEdel_sub_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine TPE1_delonly_epwrap

  subroutine OPEplusNLOdeltaonly_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl(1:NCH, 1:NCH)

    potval = 0.0_NER
      call TPE0_delonly_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
      call OPE_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl)
    potval = potval + vl

  end subroutine OPEplusNLOdeltaonly_epwrap

    subroutine OPEplusTPE_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl(1:NCH, 1:NCH)

    potval = 0.0_NER
      call TPE0_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
      call OPE_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl)
    potval = potval + vl

  end subroutine OPEplusTPE_epwrap
!  ======================================================================
!  TPESFR and TPESFR1 correspond respecyively to V_tpe0 and V_tpe1 caculated
!  with spectral functoin regularization (SFR).
!  using Epelbuam04:Eur. Phys. J. A 19, 125-137(2004)
!  ======================================================================

  subroutine TPESFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = TPESFR_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = TPESFR_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = TPESFR_jpp(j, p1, p2)
        else
          potval(1, 1) = TPESFR_jmm(j, p1, p2)
          potval(2, 1) = TPESFR_jpm(j, p1, p2)
          potval(1, 2) = TPESFR_jpm(j, p2, p1)
          potval(2, 2) = TPESFR_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine TPESFR_epwrap

  subroutine TPESFR1_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = TPESFR1_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = TPESFR1_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = TPESFR1_jpp(j, p1, p2)
        else
          potval(1, 1) = TPESFR1_jmm(j, p1, p2)
          potval(2, 1) = TPESFR1_jpm(j, p1, p2)
          potval(1, 2) = TPESFR1_jpm(j, p2, p1)
          potval(2, 2) = TPESFR1_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_PP(regtype, p1, p2, Mambda)

  end subroutine TPESFR1_epwrap

!=================================================================
! set potential structure for Born approximation

!=================================================================

! || OPE + TPE0SFR ||

  subroutine OPEplusTPESFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl(1:NCH, 1:NCH)

    potval = 0.0_NER
      call TPESFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
      call OPE_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl)
    potval = potval + vl

  end subroutine OPEplusTPESFR_epwrap

! || OPE + TPE0SFR +TPE1SFR  ||

  subroutine OPEplusAllTPESFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl0(1:NCH, 1:NCH),vl1(1:NCH, 1:NCH)

    potval = 0.0_NER
      call TPESFR1_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl1)
      call TPESFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
      call OPE_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl0)
    potval = potval + vl0 + vl1

  end subroutine OPEplusAllTPESFR_epwrap

! ||OPE + TPE0SFR + NLOdeltaSFR||

  subroutine OPEplusTPE1deltaSFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl(1:NCH, 1:NCH),vl1(1:NCH, 1:NCH),vl2(1:NCH, 1:NCH)

    potval = 0.0_NER
      call TPESFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
      call TPE0_delonly_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl1)
      call OPE_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl)

    potval = potval + vl + vl1

  end subroutine OPEplusTPE1deltaSFR_epwrap

! ||OPE + TPE0SFR + NLOdeltaSFR + NNLOdeltaSFR||

  subroutine OPEplusallTPEdeltaSFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl(1:NCH, 1:NCH),vl1(1:NCH, 1:NCH),vl2(1:NCH, 1:NCH),vl3(1:NCH, 1:NCH)

    potval = 0.0_NER
      call OPE_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
      call TPESFR_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl)
      call TPE0_delonly_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl1)
      call TPE1_delonly_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl2)
      call TPESFR1_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl3)

    potval = potval + vl + vl1 + vl2 + vl3

  end subroutine OPEplusallTPEdeltaSFR_epwrap

! ==================================================================
! set SFRpotential structure for suppope power counting
!
! ==================================================================

  subroutine VN3LO_pope_del_SFR_smpl_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: v1(1:NCH, 1:NCH), v2(1:NCH, 1:NCH), v3(1:NCH, 1:NCH), v4(1:NCH, 1:NCH)

    potval = 0.0_NER

    call TPESFR_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, v1)
    call TPE0_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, v2)
    v3 = v1 + v2
    select case (L)
      case(1)
        call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, v4)
        v4(1,1) = Cs(2) * v4(1,1)
        potval = v3 + v4
      case default
        potval = v3
    end select

  end subroutine VN3LO_pope_del_SFR_smpl_epwrap

  subroutine VN4LO_pope_del_SFR_smpl_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH), vs(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    real(NER) :: vq(1:NCH, 1:NCH), vdel(1:NCH, 1:NCH)
    potval = 0.0_NER

    call TPESFR1_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe)
    call TPE1_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vdel)
    vtpe = vtpe + vdel
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vq)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    select case (L)
      case (1)
        if(L /= j .and. j /= 0) then  ! 3p2-3f2
          call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
          potval = vtpe + Cs(3) * vq + Cs(4) * vs + Cs(5) * vl
        else
          potval(1, 1) = vtpe(1, 1) + Cs(3) * vq(1, 1) + Cs(4) * vs(1, 1)
        end if
      case (2)
        potval = vtpe + Cs(1) * vq
      case default
        potval = vtpe
    end select

  end subroutine VN4LO_pope_del_SFR_smpl_epwrap
! #############################################################
! the following structures for pope
! #############################################################

  subroutine VNLO_del_SFR_smpl_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: v1(1:NCH, 1:NCH)

    potval = 0.0_NER
    call OPE_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)
    call VN2LO_smplper_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, v1)

    potval = potval + v1

  end subroutine VNLO_del_SFR_smpl_epwrap

  subroutine VN2LO_del_SFR_smpl_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: v1(1:NCH, 1:NCH), v2(1:NCH, 1:NCH)

    potval = 0.0_NER
    call TPESFR_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)
    call TPE0_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, v1)
    select case (L)
      case(1)
        call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, v2)
        v2(1,1) = Cs(2) * v2(1,1)
        potval = potval + v1 + v2
      case default
        potval = potval + v1
      end select

  end subroutine VN2LO_del_SFR_smpl_epwrap

  subroutine VN3LO_del_SFR_smpl_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vdel1(1:NCH, 1:NCH), vsfr1(1:NCH, 1:NCH), vct1(1:NCH, 1:NCH)
    real(NER) :: vq(1:NCH, 1:NCH), vs(1:NCH, 1:NCH)

    potval = 0.0_NER
    call TPESFR1_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vsfr1)
    call TPE1_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vdel1)
    vsfr1 = vsfr1 + vdel1
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct1)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vq)
    select case(L)
      case(1)
        if(L /= j .and. j /= 0) then  ! 3p2-3f2
          call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
          potval = vsfr1 + Cs(3) * vct1 + Cs(4) * vq + Cs(5) * vs
        else
          potval(1,1) = vsfr1(1,1) + Cs(3) * vct1(1,1) + Cs(4) * vq(1,1)
        end if
      case(2)
        potval = vsfr1 + Cs(1) * vct1
      case default
        potval = vsfr1
    end select

  end subroutine VN3LO_del_SFR_smpl_epwrap

!#################################################
!# VN3LO = V_TPE^0 + V_vt^3 + V_vt^4
!# for P-waves
!##############################################

  subroutine VN3LO_pope_del_SFR_smpl_V3V4_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH), vdel(1:NCH, 1:NCH)
    real(NER) :: vct3(1:NCH, 1:NCH), vct4(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)

    potval = 0.0_NER

    call TPESFR_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe)
    call TPE0_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vdel)
    vtpe = vtpe + vdel
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct3)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct4)

      if(L /= j .and. j /= 0) then  ! 3p2-3f2
        call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
          potval = vtpe + Cs(2) * vct3 + Cs(3) * vct4 + Cs(4) * vl
      else
          potval(1, 1) = vtpe(1, 1) + Cs(2) * vct3(1, 1) + Cs(3) * vct4(1, 1)
      end if

  end subroutine VN3LO_pope_del_SFR_smpl_V3V4_epwrap

  subroutine VN4LO_pope_del_SFR_smpl_V3V4_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe1(1:NCH, 1:NCH), vdel1(1:NCH, 1:NCH)
    real(NER) :: vct3(1:NCH, 1:NCH), vct4(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    real(NER) :: vct5(1:NCH, 1:NCH), vl1(1:NCH, 1:NCH), vctf(1:NCH, 1:NCH)

    potval = 0.0_NER

    call TPESFR1_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe1)
    call TPE1_delonly_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vdel1)
    vtpe1 = vtpe1 + vdel1
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct3)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct4)
    call VPNQ5_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vctf)
    call VPNQ5_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct5)
      if(L /= j .and. j /= 0) then  ! 3p2-3f2
        call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
        call VPNQ5_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl1)
          potval = vtpe1 + Cs(5) * vct3 + Cs(6) * vct4 + Cs(7) * vct5          &
          & + Cs(8) * vctf + Cs(9) * vl +Cs(10) * vl1
      else
          potval(1, 1) = vtpe1(1, 1) + Cs(4) * vct3(1, 1) + Cs(5) * vct4(1, 1) &
          & + Cs(6) * vctf(1,1)
      end if

  end subroutine VN4LO_pope_del_SFR_smpl_V3V4_epwrap

  subroutine VN3LO_pope_SFR_smpl_V3V4_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH)
    real(NER) :: vct3(1:NCH, 1:NCH), vct4(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    potval = 0.0_NER

    call TPESFR_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct3)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct4)

      if(L /= j .and. j /= 0) then  ! 3p2-3f2
        call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
          potval = vtpe + Cs(2) * vct3 + Cs(3) * vct4 + Cs(4) * vl
      else
          potval(1, 1) = vtpe(1, 1) + Cs(2) * vct3(1, 1) + Cs(3) * vct4(1, 1)
      end if

  end subroutine VN3LO_pope_SFR_smpl_V3V4_epwrap

  subroutine VN4LO_pope_SFR_smpl_V3V4_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe1(1:NCH, 1:NCH)
    real(NER) :: vct3(1:NCH, 1:NCH), vct4(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    real(NER) :: vct5(1:NCH, 1:NCH), vl1(1:NCH, 1:NCH), vctf(1:NCH, 1:NCH)

    potval = 0.0_NER

    call TPESFR1_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe1)
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct3)
    call VPNQ5_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vctf)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct4)
    call VPNQ5_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct5)
      if(L /= j .and. j /= 0) then  ! 3p2-3f2
        call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
        call VPNQ5_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl1)
          potval = vtpe1 + Cs(5) * vct3 + Cs(6) * vct4 + Cs(7) * vct5          &
          & + Cs(8) * vctf + Cs(9) * vl +Cs(10) * vl1
      else
          potval(1, 1) = vtpe1(1, 1) + Cs(4) * vct3(1, 1) + Cs(5) * vct4(1, 1) &
          & + Cs(6) * vctf(1,1)
      end if

  end subroutine VN4LO_pope_SFR_smpl_V3V4_epwrap

  subroutine VN3LO_pope_DR_smpl_V3V4_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH)
    real(NER) :: vct3(1:NCH, 1:NCH), vct4(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    potval = 0.0_NER

    call TPE0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct3)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct4)

      if(L /= j .and. j /= 0) then  ! 3p2-3f2
        call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
          potval = vtpe + Cs(2) * vct3 + Cs(3) * vct4 + Cs(4) * vl
      else
          potval(1, 1) = vtpe(1, 1) + Cs(2) * vct3(1, 1) + Cs(3) * vct4(1, 1)
      end if

  end subroutine VN3LO_pope_DR_smpl_V3V4_epwrap

  subroutine VN4LO_pope_DR_smpl_V3V4_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe1(1:NCH, 1:NCH)
    real(NER) :: vct3(1:NCH, 1:NCH), vct4(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    real(NER) :: vct5(1:NCH, 1:NCH), vl1(1:NCH, 1:NCH), vctf(1:NCH, 1:NCH)

    potval = 0.0_NER

    call TPE1_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe1)
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct3)
    call VPNQ5_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vctf)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct4)
    call VPNQ5_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vct5)
      if(L /= j .and. j /= 0) then  ! 3p2-3f2
        call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
        call VPNQ5_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl1)
          potval = vtpe1 + Cs(5) * vct3 + Cs(6) * vct4 + Cs(7) * vct5          &
          & + Cs(8) * vctf + Cs(9) * vl +Cs(10) * vl1
      else
          potval(1, 1) = vtpe1(1, 1) + Cs(4) * vct3(1, 1) + Cs(5) * vct4(1, 1) &
          & + Cs(6) * vctf(1,1)
      end if

  end subroutine VN4LO_pope_DR_smpl_V3V4_epwrap

!============================================================
!
! Deltaful pion-exchanges
!
!============================================================

  subroutine TPE0_del_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = TPEdel_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = TPEdel_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = TPEdel_jpp(j, p1, p2)
        else
          potval(1, 1) = TPEdel_jmm(j, p1, p2)
          potval(2, 1) = TPEdel_jpm(j, p1, p2)
          potval(1, 2) = TPEdel_jpm(j, p2, p1)
          potval(2, 2) = TPEdel_jpp(j, p1, p2)
        end if
      end if
    end if
    call TPE0_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl)
    potval = potval*regltr_PP(regtype, p1, p2, Mambda)
    potval = vl + potval

  end subroutine TPE0_del_epwrap

  subroutine TPE1_del_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vl(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = TPEdel_sub_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = TPEdel_sub_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = TPEdel_sub_jpp(j, p1, p2)
        else
          potval(1, 1) = TPEdel_sub_jmm(j, p1, p2)
          potval(2, 1) = TPEdel_sub_jpm(j, p1, p2)
          potval(1, 2) = TPEdel_sub_jpm(j, p2, p1)
          potval(2, 2) = TPEdel_sub_jpp(j, p1, p2)
        end if
      end if
    end if
    potval = potval*regltr_PP(regtype, p1, p2, Mambda)
    call TPE1_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vl)
    potval = vl + potval

  end subroutine TPE1_del_epwrap

!============================================================
!
! Deltaful pion-exchanges plus withc0 short-range pots
!
!============================================================

  subroutine VN2LO_withc0_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPE0_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    call VN2LO_withc0_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN2LO_withc0_del_epwrap

  subroutine VN3LO_withc0_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER), dimension(1:NCH, 1:NCH) :: vl, vs

    call TPE1_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
    call VN3LO_withc0_SHRG_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    potval = vl + vs

  end subroutine VN3LO_withc0_del_epwrap



!============================================================
!
! Deltaful pion-exchanges plus 1s0 separable short-range pots
!
!============================================================

  subroutine VN2LO_sprbl_del_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe0(1:NCH, 1:NCH)

    call TPE0_del_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vtpe0)
    call VN2LO_sprbl_SHRG_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
    potval(1, 1) = potval(1, 1) + vtpe0(1, 1)

  end subroutine VN2LO_sprbl_del_epwrap

  subroutine VN3LO_sprbl_del_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe1(1:NCH, 1:NCH)

    call TPE1_del_epwrap(L, S, j, regtype, Mambda, para, p1, p2, vtpe1)
    call VN3LO_sprbl_SHRG_epwrap(L, S, j, regtype, Mambda, para, p1, p2, potval)
    potval(1, 1) = potval(1, 1) + vtpe1(1, 1)

  end subroutine VN3LO_sprbl_del_epwrap

  ! Perturbative channels for only deltaful

  subroutine VN3LO_pope_del_smpl_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH), vs(1:NCH, 1:NCH)

    potval = 0.0_NER

    call TPE0_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe)
    select case (L)
      case(1)
        call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
        vs(1,1) = Cs(2) * vs(1,1)
        potval = vtpe + vs
      case default
        potval = vtpe
    end select

  end subroutine VN3LO_pope_del_smpl_epwrap

  ! Perturbative channels for only deltaful

  subroutine VN4LO_pope_del_smpl_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vtpe(1:NCH, 1:NCH), vs(1:NCH, 1:NCH), vl(1:NCH, 1:NCH)
    real(NER) :: vq(1:NCH, 1:NCH)
    potval = 0.0_NER

    call TPE1_del_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vtpe)
    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vq)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    select case (L)
      case (1)
        if(L /= j .and. j /= 0) then  ! 3p2-3f2
          call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
          potval = vtpe + Cs(3) * vq + Cs(4) * vs + Cs(5) * vl
        else
          potval(1, 1) = vtpe(1, 1) + Cs(3) * vq(1, 1) + Cs(4) * vs(1, 1)
        end if
      case (2)
        potval = vtpe + Cs(1) * vq
      case default
        potval = vtpe
    end select

  end subroutine VN4LO_pope_del_smpl_epwrap
end module delta_epmod
