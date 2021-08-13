! new power counting test
! set potential with special counterterms

module testwrap

  use potwrap
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

!   abstract interface
!     subroutine gnrc_eft_pot(L, S, J, regtype, Mambda, para, p1, p2, potval)
!       import                  :: NER, NCH
!       implicit none
!       integer, intent(in)     :: L, S, J, regtype
!       real(NER), intent(in)   :: Mambda, p1, p2, para(:)
!       real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
!     end subroutine gnrc_eft_pot
!   end interface

contains

! ########################
! KSW power counting
! #######################

  subroutine VLO_ksw_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval(1, 1) = Cs(1)*potval(1, 1)

  end subroutine VLO_ksw_smpl_epwrap

  subroutine VNLO_ksw_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: vs0(1:NCH, 1:NCH), vs2(1:NCH, 1:NCH)
    real(NER)               :: vs2off(1:NCH, 1:NCH), vs(1:NCH, 1:NCH)

    potval = 0.0_NER
    call OPE_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs0)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2)

    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      call VPNQ2_off_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2off)
!       vs = Cs(2)*vs0 + Cs(3)*vs2 + Cs(4)*vs2off
! turn off offterms
        vs = Cs(2)*vs0 + Cs(3)*vs2
    else
      vs = Cs(2)*vs0 + Cs(3)*vs2
    end if
    potval = potval + vs

  end subroutine VNLO_ksw_smpl_epwrap

  subroutine VN2LO_ksw_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      potval(1, 1) = Cs(4)*potval(1, 1)
    else
      potval(1, 1) = Cs(4)*potval(1, 1)
    end if
  end subroutine VN2LO_ksw_smpl_epwrap

  subroutine VN3LO_ksw_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      potval(1, 1) = Cs(5)*potval(1, 1)
    else
      potval(1, 1) = Cs(5)*potval(1, 1)
    end if
  end subroutine VN3LO_ksw_smpl_epwrap

  subroutine VN4LO_ksw_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    if (triqn_is_coupled(L, S, J))  then  ! is coupled channel?
      potval(1, 1) = Cs(6)*potval(1, 1)
    else
      potval(1, 1) = Cs(6)*potval(1, 1)
    end if
  end subroutine VN4LO_ksw_smpl_epwrap

! ########################
! sepble power counting
! LO = OPE + C0 +C2(pp+p'p')
! ########################

  subroutine VLO_sepble_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: Vope(1:NCH, 1:NCH), vs2(1:NCH, 1:NCH),vs0(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs0)
    call OPE_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, Vope)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2)
    potval(1, 1) = Vope(1, 1)+ Cs(1)*vs0(1, 1) + Cs(2)*vs2(1, 1)
  end subroutine VLO_sepble_smpl_epwrap

  subroutine VNLO_sepble_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: Vs4(1:NCH, 1:NCH), vs2(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2)
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs4)
    potval(1, 1) = Cs(3)*potval(1, 1) + Cs(4)*vs2(1, 1)+ Cs(5)*Vs4(1, 1)
  end subroutine VNLO_sepble_smpl_epwrap

  subroutine VN2LO_sepble_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER)               :: Vs4(1:NCH, 1:NCH), vs2(1:NCH, 1:NCH), vs6(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2)
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs4)
    call VPNQ6_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs6)
    potval(1, 1) = Cs(6)*potval(1, 1) + Cs(7)*vs2(1, 1)+ Cs(8)*Vs4(1, 1)+ Cs(9)*Vs6(1, 1)
  end subroutine VN2LO_sepble_smpl_epwrap

! ###################################################
!  vab(1, 1) = (p1*p1+p2*p2)(p1*p1*p2*p2)*reg(p1)*reg(p2)*p1*p2^L
! ###################################################

  subroutine VPNQ6_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

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

    potval(1, 1) = (p1*p1+p2*p2)*(p1*p1*p2*p2)*regltr_PP(regtype, p1, p2, Mambda)*pref

  end subroutine VPNQ6_epwrap

  subroutine VNLO_withc0_m3p0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    real(NER) :: vope(1:NCH, 1:NCH), vtpe0(1:NCH, 1:NCH), vtpe1(1:NCH, 1:NCH)
    real(NER) :: vs1(1:NCH, 1:NCH), vs2(1:NCH, 1:NCH), vs3(1:NCH, 1:NCH)

    potval = 0.0_NER
    call OPE_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vope)
    call TPE0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe0)
    call TPE1_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vtpe1)
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs1)
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs2)
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, vs3)

    potval(1, 1) = vope(1, 1) + vtpe0(1, 1) + vtpe1(1, 1) + (Cs(1)+ Cs(2)+ Cs(4))*vs1(1, 1) + &
    & (Cs(3) + Cs(5))*vs2(1, 1) + Cs(6)*vs3(1, 1)

  end subroutine VNLO_withc0_m3p0_epwrap

  subroutine VLO_withc011_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VLO_withc0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval(1, 2) = 0.0_NER
    potval(2, 1) = 0.0_NER
    potval(2, 2) = 0.0_NER
  end subroutine VLO_withc011_smpl_epwrap

  subroutine VNLO_withc0offdig_smpl_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VLO_withc0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval(1, 1) = 0.0_NER
  end subroutine VNLO_withc0offdig_smpl_epwrap


  subroutine reVPNQ2_off_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ2_off_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = potval/(Mambda**2)

  end subroutine reVPNQ2_off_epwrap

  subroutine reVPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = potval/(Mambda**2)

  end subroutine reVPNQ2_diag_epwrap

  subroutine reTPE0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call TPE0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = potval/mambda

  end subroutine reTPE0_epwrap

  subroutine reVPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = potval/mambda

  end subroutine reVPNQ0_epwrap

  subroutine reVPNQ0_SCL_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ0_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = potval/(mambda**2)

  end subroutine reVPNQ0_SCL_epwrap

  subroutine reVPNQ2_SCL_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ2_diag_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = potval/(Mambda**4)

  end subroutine reVPNQ2_SCL_diag_epwrap

  subroutine reVPNQ5_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)      :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call VPNQ5_epwrap(L, S, J, regtype, Mambda, Cs, p1, p2, potval)
    potval = potval/(Mambda**4)

  end subroutine reVPNQ5_epwrap
end module testwrap
