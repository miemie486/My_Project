! Bingwei Long, 07/15/2018: Cleaning up
! Bingwei 12/22/2010 Jerry 2/1/2011
! major revision on 04/11/2012
! partial-wave decomposition (PWD) of next-to-leading two-pion exchange


! Reference
! KBW97: Kaiser, et al., NPA 625, 758 (1997)
! i^L convention in partial-wave decomposition followed that of KBW97

module vtpe1_pwd

  use feb_type,       only: standard_error_unit
  use nneft_type,     only: NER, NEC, PI_NE
  use nneft_phyconst, only: PC_mN, PC_gA, PC_fpi, PC_mpi, strct_cis,           &
  &                         PC_Pavon_c1, PC_Pavon_c3, PC_Pavon_c4,             &
  &                         PC_RSE2_c1, PC_RSE2_c3, PC_RSE2_c4,                &
  &                         PC_RSE3_c1, PC_RSE3_c3, PC_RSE3_c4,                &
  &                         PC_RSE4_c1, PC_RSE4_c3, PC_RSE4_c4,                &
  &                         PC_default_c1, PC_default_c3, PC_default_c4,       &
  &                         PC_DELTAF1_c1, PC_DELTAF1_c3, PC_DELTAF1_c4,       &
  &                         PC_DELTAF2_c1, PC_DELTAF2_c3, PC_DELTAF2_c4,       &
  &                         PC_DELTAF3_c1, PC_DELTAF3_c3, PC_DELTAF3_c4,       &
  &                         PC_DELTAF4_c1, PC_DELTAF4_c3, PC_DELTAF4_c4,       &
  &                         PC_DELTAF5_c1, PC_DELTAF5_c3, PC_DELTAF5_c4,       &
  &                         PC_DELTAF6_c1, PC_DELTAF6_c3, PC_DELTAF6_c4,       &
  &                         PC_DELTAF7_c1, PC_DELTAF7_c3, PC_DELTAF7_c4,       &
  &                         PC_DELTAF8_c1, PC_DELTAF8_c3, PC_DELTAF8_c4,       &
  &                         PC_DELTAF9_c1, PC_DELTAF9_c3, PC_DELTAF9_c4,       &
  &                         PC_DELTAF10_c1, PC_DELTAF10_c3, PC_DELTAF10_c4,    &
  &                         PC_DELTAF11_c1, PC_DELTAF11_c3, PC_DELTAF11_c4,    &
  &                         PC_DELTAF12_c1, PC_DELTAF12_c3, PC_DELTAF12_c4
  ! lgndr_table_pwo, mshlg_pwo, and wghtslg_pwo are renamed in this module
  ! in order to shorten expression in partial-wave projection subs
  use ope_pwd,        only: ONEOVEREIGHTPI_pwo, initmshlg_pwo, Nlg_pwo,        &
  & mshlg_inited_pwo, PLtbl => lgndr_table_pwo, mshz => mshlg_pwo,             &
  & wghtz => wghtslg_pwo

 implicit none

  type (strct_cis) :: cis2PE

contains

  subroutine set_cis2PE_vtpe1(c1, c3, c4)

    real(NER), intent(in) :: c1, c3, c4

    cis2PE%c1 = c1
    cis2PE%c3 = c3
    cis2PE%c4 = c4

  end subroutine set_cis2PE_vtpe1

  ! Select ci's values according to set_n.
  ! set_n = 1 ~ 4
  !     1 : the default set used before May 10, 2018
  ! See nneft_phyconst.f90 for source of the predefined ci's:
  !   PC_Pavon_c?, PC_RSE2_c?, PC_RSE3_c?, PC_RSE4_c?
  subroutine Select_cis2PE_vtpe1(set_n)

    integer, intent(in) :: set_n

    select case (set_n)
      case (1)
        call set_cis2PE_vtpe1(PC_Pavon_c1, PC_Pavon_c3, PC_Pavon_c4)
      case (2)
        call set_cis2PE_vtpe1(PC_RSE2_c1, PC_RSE2_c3, PC_RSE2_c4)
      case (3)
        call set_cis2PE_vtpe1(PC_RSE3_c1, PC_RSE3_c3, PC_RSE3_c4)
      case (4)
        call set_cis2PE_vtpe1(PC_RSE4_c1, PC_RSE4_c3, PC_RSE4_c4)
      case (5)
        call set_cis2PE_vtpe1(PC_DELTAF1_c1, PC_DELTAF1_c3, PC_DELTAF1_c4)
      case (6)
        call set_cis2PE_vtpe1(PC_DELTAF2_c1, PC_DELTAF2_c3, PC_DELTAF2_c4)
      case (7)
        call set_cis2PE_vtpe1(PC_DELTAF3_c1, PC_DELTAF3_c3, PC_DELTAF3_c4)
      case (8)
        call set_cis2PE_vtpe1(PC_DELTAF4_c1, PC_DELTAF4_c3, PC_DELTAF4_c4)
      case (9)
        call set_cis2PE_vtpe1(PC_DELTAF5_c1, PC_DELTAF5_c3, PC_DELTAF5_c4)
      case (10)
        call set_cis2PE_vtpe1(PC_DELTAF6_c1, PC_DELTAF6_c3, PC_DELTAF6_c4)
      case (11)
        call set_cis2PE_vtpe1(PC_DELTAF7_c1, PC_DELTAF7_c3, PC_DELTAF7_c4)
      case (12)
        call set_cis2PE_vtpe1(PC_DELTAF8_c1, PC_DELTAF8_c3, PC_DELTAF8_c4)
      case (13)
        call set_cis2PE_vtpe1(PC_DELTAF9_c1, PC_DELTAF9_c3, PC_DELTAF9_c4)
      case (14)
        call set_cis2PE_vtpe1(PC_DELTAF10_c1, PC_DELTAF10_c3, PC_DELTAF10_c4)
      case (15)
        call set_cis2PE_vtpe1(PC_DELTAF11_c1, PC_DELTAF11_c3, PC_DELTAF11_c4)
      case (16)
        call set_cis2PE_vtpe1(PC_DELTAF12_c1, PC_DELTAF12_c3, PC_DELTAF12_c4)
      case default
        write (standard_error_unit, '(a)')  &
        & 'Select_cis2PE_vtpe1: unknown id of ci set. Default set pre-defined will be used'
!       case default
!         write (standard_error_unit, '(a)')  &
!         & 'Select_cis2PE_vtpe1: unknown id of ci set. Default set pre-defined will be used'
        ! call set_cis2PE_vtpe1(PC_default_c1, PC_default_c3, PC_default_c4)
    end select

  end subroutine Select_cis2PE_vtpe1

  ! Return the current values of cis2PE
  subroutine get_cis2PE_vtpe1(cisout)

    type (strct_cis), intent(out) :: cisout

    cisout%c1 = cis2PE%c1
    cisout%c3 = cis2PE%c3
    cisout%c4 = cis2PE%c4

  end subroutine get_cis2PE_vtpe1

  ! VC_TPE1, WC_TPE1, VT_TPE1, WT_TPE1, VSO_TPE1, and WSO_TPE1:
  ! next-to-leading TPE found in Eqs.(17) - (23) of KBW97
  ! - polynomials in q^2, including ln(mpi), mpi, and mpi^3 terms, are dropped
  ! - U = V + (4*I - 3)*W   (cf. ope_pwd.f90)
  ! - L + s + I = odd
  pure function VC_TPE1(qsqr, aq)

    real(NER) :: VC_TPE1
    real(NER), intent(in) :: qsqr, aq

    real(NER)   :: q, a, d, f, g
    real(NER)   :: ovalnorm_pwo, PC_mpi_sqr, PC_mpi_4th, PC_gA_sqr, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_mpi_sqr      = PC_mpi * PC_mpi
    PC_mpi_4th      = PC_mpi_sqr * PC_mpi_sqr
    PC_gA_sqr       = PC_gA * PC_gA
    PC_fpi_4th      = PC_fpi**4

    a = ovalnorm_pwo * 3.0_NER * PC_gA_sqr / (16.0_NER * PI_NE * PC_fpi_4th)
    d = PC_gA_sqr * PC_mpi_4th * PC_mpi / (16.0_NER * PC_mN)
    f = 2.0_NER * PC_mpi_sqr * (2.0_NER * cis2PE%c1 - cis2PE%c3)
    g = cis2PE%c3 + 3.0_NER * PC_gA_sqr / (16.0_NER * PC_mN)
    VC_TPE1 = a * (-d/(4.0_NER*PC_mpi_sqr + qsqr)         &
    &         + (f - g*qsqr)*(2.0_NER*PC_mpi_sqr + qsqr) * aq)

  end function VC_TPE1

  pure function WC_TPE1(qsqr, aq)

    real(NER) :: WC_TPE1
    real(NER), intent(in) :: qsqr, aq

    real(NER)   :: q, a, b, c, d
    real(NER)   :: ovalnorm_pwo, PC_mpi_sqr, PC_mpi_4th, PC_gA_sqr, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_mpi_sqr      = PC_mpi * PC_mpi
    PC_mpi_4th      = PC_mpi_sqr * PC_mpi_sqr
    PC_gA_sqr       = PC_gA * PC_gA
    PC_fpi_4th      = PC_fpi**4

    a = ovalnorm_pwo * PC_gA_sqr / (128.0_NER * PI_NE * PC_mN * PC_fpi_4th)
    b = 3.0_NER * PC_gA_sqr * PC_mpi_4th * PC_mpi
    c = 4.0_NER * PC_mpi_sqr * (1.0_NER - PC_gA_sqr)
    d = 2.0_NER - 3.0_NER * PC_gA_sqr
    WC_TPE1 = a * (-b/(4.0_NER*PC_mpi_sqr + qsqr)                       &
    &         + (c + d*qsqr)*(2.0_NER*PC_mpi_sqr + qsqr) * aq)

  end function WC_TPE1

  pure function VT_TPE1(qsqr, aq)

    real(NER) :: VT_TPE1
    real(NER), intent(in) :: qsqr, aq

    real(NER)   :: q, a
    real(NER)   :: ovalnorm_pwo, PC_mpi_sqr, PC_gA_4th, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_mpi_sqr      = PC_mpi * PC_mpi
    PC_gA_4th       = PC_gA**4
    PC_fpi_4th      = PC_fpi**4

    a = -9.0_NER*ovalnorm_pwo*PC_gA_4th / (512.0_NER*PI_NE*PC_mN*PC_fpi_4th)
    VT_TPE1 = a * (2.0_NER*PC_mpi_sqr + qsqr) * aq

  end function VT_TPE1

  pure function WT_TPE1(qsqr, aq)

    real(NER) :: WT_TPE1
    real(NER), intent(in) :: qsqr, aq
    real(NER)   :: q, a, b, c
    real(NER)   :: ovalnorm_pwo, PC_mpi_sqr, PC_gA_sqr, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_mpi_sqr      = PC_mpi * PC_mpi
    PC_gA_sqr       = PC_gA * PC_gA
    PC_fpi_4th      = PC_fpi**4

    a = ovalnorm_pwo * PC_gA_sqr / (32.0_NER*PI_NE*PC_fpi_4th)
    b = (cis2PE%c4 + 0.25_NER/PC_mN)*4.0_NER*PC_mpi_sqr             &
    &   - 1.25_NER*PC_gA_sqr*PC_mpi_sqr/PC_mN
    c = cis2PE%c4 + 0.25_NER/PC_mN - 0.375_NER*PC_gA_sqr/PC_mN
    WT_TPE1 = a * (b + c * qsqr) * aq

  end function WT_TPE1

  pure function VSO_TPE1(qsqr, aq)

    real(NER) :: VSO_TPE1
    real(NER), intent(in) :: qsqr, aq

    real(NER)   :: q, a
    real(NER)   :: ovalnorm_pwo, PC_mpi_sqr, PC_gA_4th, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_mpi_sqr      = PC_mpi * PC_mpi
    PC_gA_4th       = PC_gA**4
    PC_fpi_4th      = PC_fpi**4

    a = 3.0_NER * ovalnorm_pwo*PC_gA_4th / (64.0_NER*PI_NE*PC_mN*PC_fpi_4th)
    VSO_TPE1 = a * (2.0_NER*PC_mpi_sqr + qsqr) * aq

  end function VSO_TPE1

  pure function WSO_TPE1(qsqr, aq)

    real(NER) :: WSO_TPE1
    real(NER), intent(in) :: qsqr, aq

    real(NER)   :: q, a
    real(NER)   :: ovalnorm_pwo, PC_mpi_sqr, PC_gA_sqr, PC_fpi_4th

    ovalnorm_pwo    = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_mpi_sqr      = PC_mpi * PC_mpi
    PC_gA_sqr       = PC_gA * PC_gA
    PC_fpi_4th      = PC_fpi**4

    a = ovalnorm_pwo*PC_gA_sqr*(1.0_NER - PC_gA_sqr)          &
    &    /(64.0_NER*PI_NE*PC_mN*PC_fpi_4th)
    WSO_TPE1 = a * (4.0_NER*PC_mpi_sqr + qsqr) * aq

  end function WSO_TPE1

  ! <j0j|VTPE1|j0j>  <Lsj|
  function VTPE1_j0j(j, p, k)

    real(NER) :: VTPE1_j0j
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer :: ii
    real(NER) :: regsum, qsqr, q, aq, z, ut, uc, us, psqr, ksqr, pk

    if (j < 0) then
      VTPE1_j0j = 0.0_NER
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
          z = mshz(ii)
          qsqr = psqr + ksqr - 2.0_NER*pk*z
          q  = sqrt(qsqr)
          aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
          uc = VC_TPE1(qsqr, aq) + WC_TPE1(qsqr, aq)
          ut = VT_TPE1(qsqr, aq) + WT_TPE1(qsqr, aq)
          us = -qsqr*ut
          regsum = regsum + wghtz(ii)*(uc - 2.0_NER*us)
        end do
      case (1)
        do ii = 1, Nlg_pwo
          z = mshz(ii)
          qsqr = psqr + ksqr - 2.0_NER*pk*z
          q  = sqrt(qsqr)
          aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
          uc = VC_TPE1(qsqr, aq) - 3.0_NER*WC_TPE1(qsqr, aq)
          ut = VT_TPE1(qsqr, aq) - 3.0_NER*WT_TPE1(qsqr, aq)
          us = -qsqr*ut
          regsum = regsum + wghtz(ii)*(uc - 2.0_NER*us)*z
        end do
      case default
        if (mod(j, 2) == 1) then
          do ii = 1, Nlg_pwo     ! j == odd => u = v - 3w
            z = mshz(ii)
            qsqr = psqr + ksqr - 2.0_NER*pk*z
            q  = sqrt(qsqr)
            aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
            uc = VC_TPE1(qsqr, aq) - 3.0_NER*WC_TPE1(qsqr, aq)
            ut = VT_TPE1(qsqr, aq) - 3.0_NER*WT_TPE1(qsqr, aq)
            us = -qsqr*ut
            regsum = regsum + wghtz(ii)*(uc - 2.0_NER*us)*PLtbl(j, ii)
          end do
        else
          do ii = 1, Nlg_pwo     ! j == even => u = v + w
            z = mshz(ii)
            qsqr = psqr + ksqr - 2.0_NER*pk*z
            q  = sqrt(qsqr)
            aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
            uc = VC_TPE1(qsqr, aq) + WC_TPE1(qsqr, aq)
            ut = VT_TPE1(qsqr, aq) + WT_TPE1(qsqr, aq)
            us = -qsqr*ut
            regsum = regsum + wghtz(ii)*(uc - 2.0_NER*us)*PLtbl(j, ii)
          end do
        end if
    end select
    VTPE1_j0j = -regsum

  end function VTPE1_j0j

  ! <j1j|VTPE1|j1j>
  function VTPE1_j1j(j, p, k)

    real(NER) :: VTPE1_j1j
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer :: ii
    real(NER) :: regsum, qsqr, q, aq, z, pk, psqr, ksqr, uc, ut, uso

    if (j < 1) then
      VTPE1_j1j = 0.0_NER
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
      do ii = 1, Nlg_pwo        ! j == odd => u = v + w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        uc = VC_TPE1(qsqr, aq) + WC_TPE1(qsqr, aq)
        ut = VT_TPE1(qsqr, aq) + WT_TPE1(qsqr, aq)
        uso = VSO_TPE1(qsqr, aq) + WSO_TPE1(qsqr, aq)
        regsum = regsum + wghtz(ii)*((uc + 4.0_NER*pk*z*(ut - uso))       &
        & * PLtbl(j, ii)                                                  &
        & + 2.0_NER*pk*(uso - ut)*(PLtbl(j-1, ii) + PLtbl(j+1, ii)))
      end do
    else
      do ii = 1, Nlg_pwo        ! j == even => u = v - 3w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        uc = VC_TPE1(qsqr, aq) - 3.0_NER*WC_TPE1(qsqr, aq)
        ut = VT_TPE1(qsqr, aq) - 3.0_NER*WT_TPE1(qsqr, aq)
        uso = VSO_TPE1(qsqr, aq) - 3.0_NER*WSO_TPE1(qsqr, aq)
        regsum = regsum                                                     &
        & + wghtz(ii)*((uc + 4.0_NER*pk*z*(ut - uso))*PLtbl(j, ii)          &
        & + 2.0_NER*pk*(uso - ut)*(PLtbl(j-1, ii) + PLtbl(j+1, ii)))
      end do
    end if

    VTPE1_j1j = -regsum

  end function VTPE1_j1j

  ! <j+1, 1j|VTPE1|j+1, 1j>
  function VTPE1_jpp(j, p, k)

    real(NER) :: VTPE1_jpp
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer   :: ii
    real(NER) :: regsum, qsqr, q, aq, z, pk, psqr, ksqr, jr, pref1, pref2,  &
    & uc, ut, us, uso

    if (j < 0) then
      VTPE1_jpp = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    jr = real(j, NER)
    regsum = 0.0_NER
    pref1 = 2.0_NER/(2.0_NER*jr+1.0_NER)
    pref2 = 1.0_NER/(2.0_NER*jr+1.0_NER)
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    if (mod(j, 2) == 1) then
      do ii = 1, Nlg_pwo        ! j == odd => u = v - 3w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        uc = VC_TPE1(qsqr, aq) - 3.0_NER*WC_TPE1(qsqr, aq)
        ut = VT_TPE1(qsqr, aq) - 3.0_NER*WT_TPE1(qsqr, aq)
        uso = VSO_TPE1(qsqr, aq) - 3.0_NER*WSO_TPE1(qsqr, aq)
        us = -qsqr*ut
        regsum = regsum                                                     &
        & + wghtz(ii)*((2.0_NER*uso + pref1*ut)*pk*PLtbl(j, ii)             &
        & + (uc + us - 2.0_NER*pk*z*uso - pref2*(psqr + ksqr)*ut)           &
        &   *PLtbl(j+1, ii))
      end do
    else
      do ii = 1, Nlg_pwo        ! j == even => u = v + w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        uc = VC_TPE1(qsqr, aq) + WC_TPE1(qsqr, aq)
        ut = VT_TPE1(qsqr, aq) + WT_TPE1(qsqr, aq)
        uso = VSO_TPE1(qsqr, aq) + WSO_TPE1(qsqr, aq)
        us   = -qsqr*ut
        regsum = regsum                                                     &
        & + wghtz(ii)*((2.0_NER*uso + pref1*ut)*pk*PLtbl(j, ii)             &
        & + (uc + us - 2.0_NER*pk*z*uso - pref2*(psqr + ksqr)*ut)           &
        &   *PLtbl(j+1, ii))
      end do
    end if

    VTPE1_jpp = -regsum

  end function VTPE1_jpp


  ! <j-1, 1j|VTPE1|j-1, 1j>
  function VTPE1_jmm(j, p, k)

    real(NER) :: VTPE1_jmm
    real(NER), intent(in) :: p, k
    integer, intent(in) :: j

    integer   :: ii
    real(NER) :: regsum, qsqr, q, aq, z, pk, psqr, ksqr, pref1, pref2, uc,  &
    & ut, us, uso, jr

    if (j < 0) then
      VTPE1_jmm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    jr = real(j, NER)
    regsum = 0.0_NER
    pref1 = 2.0_NER/(2.0_NER*jr+1.0_NER)
    pref2 = 1.0_NER/(2.0_NER*jr+1.0_NER)
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    if (mod(j, 2) == 1) then
      do ii = 1, Nlg_pwo        ! j == odd => u = v - 3w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        uc = VC_TPE1(qsqr, aq) - 3.0_NER*WC_TPE1(qsqr, aq)
        ut = VT_TPE1(qsqr, aq) - 3.0_NER*WT_TPE1(qsqr, aq)
        uso = VSO_TPE1(qsqr, aq) - 3.0_NER*WSO_TPE1(qsqr, aq)
        us = -qsqr*ut
        regsum = regsum                                                     &
        & + wghtz(ii)*((2.0_NER*uso - pref1*ut)*pk*PLtbl(j, ii)             &
        & + (uc + us - 2.0_NER*pk*z*uso + pref2*(psqr + ksqr)*ut)           &
        &    *PLtbl(j-1, ii))
      end do
    else
      do ii = 1, Nlg_pwo        ! j == even => u = v + w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        uc = VC_TPE1(qsqr, aq) + WC_TPE1(qsqr, aq)
        ut = VT_TPE1(qsqr, aq) + WT_TPE1(qsqr, aq)
        uso = VSO_TPE1(qsqr, aq) + WSO_TPE1(qsqr, aq)
        us   = -qsqr*ut
        regsum = regsum                                                     &
        & + wghtz(ii)*((2.0_NER*uso - pref1*ut)*pk*PLtbl(j, ii)             &
        & + (uc + us - 2.0_NER*pk*z*uso + pref2*(psqr + ksqr)*ut)           &
        &    *PLtbl(j-1, ii))
      end do
    end if

    VTPE1_jmm = -regsum

  end function VTPE1_jmm

  ! <j+1, 1j|VTPE1|j-1, 1j>
  function VTPE1_jpm(j, p, k)

    real(NER) :: VTPE1_jpm
    real(NER), intent(in) :: p, k
    integer, intent(in) :: j

    integer :: ii
    real(NER) :: regsum, qsqr, q, aq, z, pk, psqr, ksqr, ut, jr

    if (j < 0) then
      VTPE1_jpm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    jr = real(j, NER)
    if (mod(j, 2) == 1) then
      do ii = 1, Nlg_pwo      ! j == odd => u = v - 3w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        ut = VT_TPE1(qsqr, aq) - 3.0_NER*WT_TPE1(qsqr, aq)
        regsum = regsum + wghtz(ii)*ut*(-4.0_NER*pk*PLtbl(j, ii)            &
        & + 2.0_NER*psqr*PLtbl(j-1, ii) + 2.0_NER*ksqr*PLtbl(j+1, ii))
      end do
    else
      do ii = 1, Nlg_pwo      ! j == even => u = v + w
        z = mshz(ii)
        qsqr = psqr + ksqr - 2.0_NER*pk*z
        q  = sqrt(qsqr)
        aq = atan(q/PC_mpi*0.5_NER)/q*0.5_NER
        ut = VT_TPE1(qsqr, aq) + WT_TPE1(qsqr, aq)
        regsum = regsum + wghtz(ii)*ut*(-4.0_NER*pk*PLtbl(j, ii)            &
        & + 2.0_NER*psqr*PLtbl(j-1, ii) + 2.0_NER*ksqr*PLtbl(j+1, ii))
      end do
    end if
    VTPE1_jpm = sqrt(jr*(jr+1.0_NER))/(2.0_NER*jr+1.0_NER)*regsum

  end function VTPE1_jpm

end module vtpe1_pwd
