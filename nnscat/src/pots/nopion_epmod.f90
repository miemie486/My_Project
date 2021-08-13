!pray  potwrap for pionless

module nopion_epmod

  use potwrap
  implicit none

  integer, parameter, private  :: NCH = 2
contains

  subroutine VLO_withc0_pionless_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    potval(1, 1) = Cs(1) * exp(-(p1**2 + p2**2)/Mambda**2)

  end subroutine VLO_withc0_pionless_epwrap

  subroutine VLO_withc0_pionless_NP_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vs(1:NCH, 1:NCH)

    potval = 0.0_NER
! print *, exp(-(p1**2 + p2**2)/Mambda**2)
!     call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
!     potval  = Cs(1) * vs

    potval(1,1) = Cs(1) * exp(-(p1**2 + p2**2)/Mambda**2)

  end subroutine VLO_withc0_pionless_NP_epwrap


  subroutine VLO_withc0_pionless_NPLO_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    real(NER) :: vs(1:NCH, 1:NCH), vl(1:NCH,1:NCH)

    potval = 0.0_NER

    call VPNQ0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vs)
    call VPNQ2_diag_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    if(L /= j .and. j /= 0) then
      call VPNQ2_off_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, vl)
      potval = Cs(1) * Vs + potval * Cs(2) + Cs(3) * vl
    else
      potval = Cs(1) * vs + Cs(2) * potval
    end if
  end subroutine VLO_withc0_pionless_NPLO_epwrap

end module nopion_epmod
