! Bingwei Long, 07/15/2018: Cleaning up
! Bingwei 12/22/2010
! major revision on 04/11/2012
! partial-wave decomposition (PWD) of delta-less, leading two-pion exchange

! Reference
! KBW97: Kaiser, et al., NPA 625, 758 (1997)
! i^L convention in partial-wave decomposition followed that of KBW97

module vtpe0_pwd

  use nneft_type,     only: NER, NEC, PI_NE
  use nneft_phyconst, only: PC_mN, PC_gA, PC_fpi, PC_mpi

  ! lgndr_table_pwo, mshlg_pwo, and wghtslg_pwo are renamed in this module
  ! in order to shorten expression in partial-wave projection subs
  use ope_pwd,        only: ONEOVEREIGHTPI_pwo, initmshlg_pwo, Nlg_pwo,        &
  & mshlg_inited_pwo, PLtbl => lgndr_table_pwo, mshz => mshlg_pwo,             &
  & wghtz => wghtslg_pwo
  implicit none

  contains

  ! WC_VTPE0 and VT_NLO: Leading TPE found in Eqs.(14) & (15) of KBW97
  ! - polynomials in q^2, including ln(mpi) terms, are dropped
  ! - U = V + (4*I - 3)*W   (cf. ope_pwd.f90)
  ! - L + s + I = odd
  ! - VS_VTPE0 = -q^2*VT_VTPE0 is implemented in projection subs like
  !   VTPE0_j0j, VTPE0_j1j, VTPE0_jpp, etc.

  ! overnorm_pwo = 2.0*mN/pi * 1/8pi
  ! WC_VTPE0 = ovalnorm_pwo*(384*pi^2*fpi^4)^{-1}*[4mpi^2*(5*gA^4 - 4*gA^2 - 1)
  !             + q^2*(23*gA^4 - 10*gA^2 - 1) + 48*gA^4*mpi^4/(4mpi^2 + q^2)]
  !             *L(q)
  pure function WC_VTPE0(qsqr, Lq)

    real(NER) :: WC_VTPE0
    real(NER), intent(in)   :: qsqr, Lq

    real(NER)   :: w, q, a, b, c, pref
    real(NER)   :: ovalnorm_pwo, mpi_sqr_quad, PC_mpi_sqr, PC_mpi_4th,         &
    & PC_gA_sqr, PC_gA_4th, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_mpi_sqr      = PC_mpi * PC_mpi
    PC_mpi_4th      = PC_mpi_sqr * PC_mpi_sqr
    PC_gA_sqr       = PC_gA * PC_gA
    PC_gA_4th       = PC_gA_sqr * PC_gA_sqr
    PC_fpi_4th      = PC_fpi**4
    mpi_sqr_quad    = 4.0_NER * PC_mpi_sqr

    a = mpi_sqr_quad * (5.0_NER*PC_gA_4th - 4.0_NER*PC_gA_sqr - 1.0_NER)
    b = 23.0_NER*PC_gA_4th - 10.0_NER*PC_gA_sqr - 1.0_NER
    c = 48.0_NER*PC_gA_4th*PC_mpi_4th
    pref = ovalnorm_pwo / (384.0_NER*PI_NE*PI_NE*PC_fpi_4th)

    WC_VTPE0 = pref * (a + b*qsqr + c/(mpi_sqr_quad + qsqr)) * Lq

  end function WC_VTPE0

  ! VT_VTPE0 = ovalnorm_pwo * 3*gA^4/(64*pi^2*fpi^4)*L(q)
  pure function VT_VTPE0(Lq)

    real(NER) :: VT_VTPE0
    real(NER), intent(in)   :: Lq

    real(NER)   :: w, q
    real(NER)   :: pref
    real(NER)   :: ovalnorm_pwo, PC_gA_4th, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_gA_4th       = PC_gA**4
    PC_fpi_4th      = PC_fpi**4

    pref = ovalnorm_pwo*3.0_NER*PC_gA_4th / (64.0_NER*PI_NE*PI_NE*PC_fpi_4th)
    VT_VTPE0 = pref * Lq

  end function VT_VTPE0

  pure subroutine get_Lq_from_qsqr(Lq, qsqr)

    real(NER), intent(out)  :: Lq
    real(NER), intent(in)   :: qsqr

    real(NER)   :: q, w

    q   = sqrt(qsqr)
    w   = sqrt(4.0_NER * PC_mpi * PC_mpi + qsqr)
    Lq  = w/q * log((w+q)/PC_mpi*0.5_NER)

  end subroutine get_Lq_from_qsqr

  ! <j0j|VTPE0|j0j>
  function VTPE0_j0j(j, p, k)

    real(NER)               :: VTPE0_j0j
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, Lq, z, ut, uc, us, psqr, ksqr, pk

    if (j < 0) then
      VTPE0_j0j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    select case (j)
      case (0)
        do ii = 1, Nlg_pwo
          z       = mshz(ii)
          qsqr    = psqr + ksqr - 2.0_NER*pk*z
          call get_Lq_from_qsqr(Lq, qsqr)
          uc      = WC_VTPE0(qsqr, Lq)
          ut      = VT_VTPE0(Lq)
          us      = -qsqr*ut
          regsum  = regsum + wghtz(ii)*(uc - 2.0_NER*us)
        end do
      case (1)
        do ii = 1, Nlg_pwo
          z       = mshz(ii)
          qsqr    = psqr + ksqr - 2.0_NER*pk*z
          call get_Lq_from_qsqr(Lq, qsqr)
          uc      = -3.0_NER*WC_VTPE0(qsqr, Lq)
          ut      = VT_VTPE0(Lq)
          us      = -qsqr*ut
          regsum  = regsum + wghtz(ii)*(uc - 2.0_NER*us)*z
        end do
      case default
        if (mod(j, 2) == 1) then
          do ii = 1, Nlg_pwo     ! j == odd
            z       = mshz(ii)
            qsqr    = psqr + ksqr - 2.0_NER*pk*z
            call get_Lq_from_qsqr(Lq, qsqr)
            uc      = -3.0_NER*WC_VTPE0(qsqr, Lq)
            ut      = VT_VTPE0(Lq)
            us      = -qsqr*ut
            regsum  = regsum + wghtz(ii)*(uc - 2.0_NER*us)*PLtbl(j, ii)
          end do
        else
          do ii = 1, Nlg_pwo     ! j == even
            z       = mshz(ii)
            qsqr    = psqr + ksqr - 2.0_NER*pk*z
            call get_Lq_from_qsqr(Lq, qsqr)
            uc      = WC_VTPE0(qsqr, Lq)
            ut      = VT_VTPE0(Lq)
            us      = -qsqr*ut
            regsum  = regsum + wghtz(ii)*(uc - 2.0_NER*us)*PLtbl(j, ii)
          end do
        end if
    end select

    VTPE0_j0j = -regsum

  end function VTPE0_j0j

  ! <j1j|VTPE0|j1j>
  function VTPE0_j1j(j, p, k)

    real(NER)               :: VTPE0_j1j
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, Lq, z, pk, psqr, ksqr, uc, ut

    if (j < 1) then
      VTPE0_j1j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    if (mod(j, 2) == 1) then
      do ii = 1, Nlg_pwo        ! j == odd => uc = wc
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        call get_Lq_from_qsqr(Lq, qsqr)
        uc   = WC_VTPE0(qsqr, Lq)
        ut   = VT_VTPE0(Lq)
        regsum = regsum + wghtz(ii)*((uc + 4.0_NER*pk*z*ut)*PLtbl(j,ii)        &
        &                 - 2.0_NER*pk*ut*(PLtbl(j-1, ii) + PLtbl(j+1, ii)))
      end do
    else
      do ii = 1, Nlg_pwo        ! j == even => uc = -3.0_NER*wc
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        call get_Lq_from_qsqr(Lq, qsqr)
        uc   = -3.0_NER*WC_VTPE0(qsqr, Lq)
        ut   = VT_VTPE0(Lq)
        regsum = regsum + wghtz(ii)*((uc + 4.0_NER*pk*z*ut)*PLtbl(j, ii)       &
        &                 - 2.0_NER*pk*ut*(PLtbl(j-1, ii) + PLtbl(j+1, ii)))
      end do
    end if

    VTPE0_j1j = -regsum

  end function VTPE0_j1j

  ! <j+1, 1j|VTPE0|j+1, 1j>
  function VTPE0_jpp(j, p, k)

    real(NER)               :: VTPE0_jpp
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, Lq, z, pk, psqr, ksqr, pref1, pref2, uc, ut, us

    if (j < 0) then
      VTPE0_jpp = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum  = 0.0_NER
    pref1   = 2.0_NER/(2.0_NER*real(j, NER)+1.0_NER)
    pref2   = 1.0_NER/(2.0_NER*real(j, NER)+1.0_NER)
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    if (mod(j, 2) == 1) then
      do ii = 1, Nlg_pwo        ! j == odd => uc = -3.0*wc
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        call get_Lq_from_qsqr(Lq, qsqr)
        uc   = -3.0_NER*WC_VTPE0(qsqr, Lq)
        ut   = VT_VTPE0(Lq)
        us   = -qsqr*ut
        regsum = regsum + wghtz(ii)*(pref1*pk*ut*PLtbl(j, ii)         &
        & + (uc + us - pref2*(psqr + ksqr)*ut)*PLtbl(j+1, ii))
      end do
    else
      do ii = 1, Nlg_pwo        ! j == even => uc = wc
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        call get_Lq_from_qsqr(Lq, qsqr)
        uc   = WC_VTPE0(qsqr, Lq)
        ut   = VT_VTPE0(Lq)
        us   = -qsqr*ut
        regsum = regsum + wghtz(ii)*(pref1*pk*ut*PLtbl(j, ii)         &
        & + (uc + us - pref2*(psqr + ksqr)*ut)*PLtbl(j+1, ii))
      end do
    end if

    VTPE0_jpp = -regsum

  end function VTPE0_jpp

  ! <j-1, 1j|VTPE0|j-1, 1j>
  function VTPE0_jmm(j, p, k)

    real(NER)               :: VTPE0_jmm
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, Lq, z, pk, psqr, ksqr, pref1, pref2, uc, ut, us

    if (j < 0) then
      VTPE0_jmm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    pref1 = 2.0_NER/(2.0_NER*real(j, NER)+1.0_NER)
    pref2 = 1.0_NER/(2.0_NER*real(j, NER)+1.0_NER)
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    if (mod(j, 2) == 1) then
      do ii = 1, Nlg_pwo        ! j == odd => uc = -3.0*wc
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        call get_Lq_from_qsqr(Lq, qsqr)
        uc   = -3.0_NER*WC_VTPE0(qsqr, Lq)
        ut   = VT_VTPE0(Lq)
        us   = -qsqr*ut
        regsum = regsum + wghtz(ii)*(-pref1*pk*ut*PLtbl(j, ii)        &
        &      + (uc + us + pref2*(psqr + ksqr)*ut)*PLtbl(j-1, ii))
      end do
    else
      do ii = 1, Nlg_pwo        ! j == even => uc = wc
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        call get_Lq_from_qsqr(Lq, qsqr)
        uc   = WC_VTPE0(qsqr, Lq)
        ut   = VT_VTPE0(Lq)
        us   = -qsqr*ut
        regsum = regsum + wghtz(ii)*(-pref1*pk*ut*PLtbl(j, ii)        &
        &      + (uc + us + pref2*(psqr + ksqr)*ut)*PLtbl(j-1, ii))
      end do
    end if

    VTPE0_jmm = -regsum

  end function VTPE0_jmm

  ! <j+1, 1j|VTPE0|j-1, 1j>
  function VTPE0_jpm(j, p, k)

    real(NER)               :: VTPE0_jpm
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, Lq, z, pk, psqr, ksqr, pref1, pref2, ut

    if (j < 0) then
      VTPE0_jpm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    do ii = 1, Nlg_pwo
      z = mshz(ii)
      qsqr = psqr + ksqr - 2.0_NER*pk*z
      call get_Lq_from_qsqr(Lq, qsqr)
      ut   = VT_VTPE0(Lq)
      regsum = regsum + wghtz(ii)*ut*(-4.0_NER*pk*PLtbl(j, ii)        &
      &                 + 2.0_NER*psqr*PLtbl(j-1, ii)                 &
      &                 + 2.0_NER*ksqr*PLtbl(j+1, ii))
    end do

    VTPE0_jpm = sqrt(j*(j+1.0_NER))/(2.0_NER*real(j, NER)+1.0_NER)*regsum

  end function VTPE0_jpm

end module vtpe0_pwd
