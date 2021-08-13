! Bingwei Long  01/05/2020
! - Complex valued version of potwrap

module cmplx_epmod

  use eft_potspwd
  use ope_pwd

  integer, parameter, private  :: NCH = 2

  abstract interface
    subroutine gnrc_eft_cmplx_pot(L, S, J, regtype, Mambda, para, p1, p2, potval)
      import                  :: NER, NEC, NCH
      implicit none
      integer, intent(in)       :: L, S, J, regtype
      real(NER), intent(in)     :: Mambda, para(:)
      complex(NEC), intent(in)  :: p1, p2
      complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)
    end subroutine gnrc_eft_cmplx_pot
  end interface

contains

  subroutine zero_V_cmplx(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER

  end subroutine zero_V_cmplx

  subroutine OPE_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    if (S == 0) then
      potval(1, 1) = OPE_cmplx_j0j(j, p1, p2)
    else
      if (j == L) then
        potval(1, 1) = OPE_cmplx_j1j(j, p1, p2)
      else
        if (j == 0) then
          potval(1, 1) = OPE_cmplx_jpp(j, p1, p2)
        else
          potval(1, 1) = OPE_cmplx_jmm(j, p1, p2)
          potval(2, 1) = OPE_cmplx_jpm(j, p1, p2)
          potval(1, 2) = OPE_cmplx_jpm(j, p2, p1)
          potval(2, 2) = OPE_cmplx_jpp(j, p1, p2)
        end if
      end if
    end if

    potval = potval*regltr_gaussian_cmplx_PP(p1, p2, Mambda)

  end subroutine OPE_cmplx_epwrap

  subroutine VPNQ0_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    select case (L)
      case (0)
        potval(1, 1) = 1.0_NER
      case (1)
        potval(1, 1) = p1*p2
      case default
        potval(1, 1) = (p1*p2)**L
    end select

    potval(1, 1) = potval(1, 1)*regltr_gaussian_cmplx_PP(p1, p2, Mambda)

  end subroutine VPNQ0_cmplx_epwrap

  subroutine VPNQ2_diag_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

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

    potval(1, 1) = (p1*p1 + p2*p2)*regltr_gaussian_cmplx_PP(p1, p2, Mambda)*pref

  end subroutine VPNQ2_diag_cmplx_epwrap

  subroutine VPNQ2_off_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

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
      pref = pref * regltr_gaussian_cmplx_PP(p1, p2, Mambda)
      potval(1, 2) = p2*p2 * pref
      potval(2, 1) = p1*p1 * pref
    end if

  end subroutine VPNQ2_off_cmplx_epwrap

  subroutine VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    complex(NEC) :: vope(1:NCH, 1:NCH), vs(1:NCH, 1:NCH)

    potval = 0.0_NER

    call OPE_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vope)
    call VPNQ0_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vs)
    vs(1, 1) = para(1) * vs(1, 1)

    potval = vope + vs

  end subroutine VLO_withc0_cmplx_epwrap

  subroutine VLO_MMWLY_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

  end subroutine VLO_MMWLY_cmplx_epwrap

  subroutine VNLO_MMWLY_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    complex(NEC) :: vs(1:NCH, 1:NCH), voff(1:NCH,1:NCH)

    call VPNQ0_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)
    call VPNQ2_diag_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vs)
    potval = para(2)*potval + para(3)*vs

  end subroutine VNLO_MMWLY_cmplx_epwrap


!===================================================================
!
! Short-range part of energy-dependent dibaryon pots
!
!===================================================================

  subroutine VLO_LBWDIB_SHRG_cmplx_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: Ecm, p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    complex(NEC)  :: C0
    real(NER)     :: CgA, tmppara(1:1)

    call VPNQ0_cmplx_epwrap(L, S, J, regtype, Mambda, tmppara, p1, p2, potval)
    CgA = PC_gA*PC_gA*PC_mN/(8.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi)
    ! C0 = y/(E + Delta) - PC_CgA
    ! Delta = paras(1), y = paras(2)
    C0 = para(2)/(Ecm + para(1)) - CgA
    potval = C0 * potval

  end subroutine VLO_LBWDIB_SHRG_cmplx_epwrap

!===================================================================
!
! Energy-dependent dibaryon pots (with deltaless pion-exchg)
!
!===================================================================

  subroutine VLO_LBWDIB_cmplx_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)

    integer, intent(in)       :: L, S, J, regtype
    real(NER), intent(in)     :: Mambda, para(:)
    complex(NEC), intent(in)  :: Ecm, p1, p2
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    complex(NEC) :: vope(1:NCH, 1:NCH)

    call VLO_LBWDIB_SHRG_cmplx_epwrap(L, S, J, regtype, Mambda, para, Ecm, p1, p2, potval)
    call OPE_cmplx_epwrap(L, S, J, regtype, Mambda, para, p1, p2, vope)
    potval = potval + vope

  end subroutine VLO_LBWDIB_cmplx_epwrap

end module cmplx_epmod
