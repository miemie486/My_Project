! Bingwei Long 07/15/2018: cleaning up and improving comments
! Major revision on 04/11/2012
! Bingwei 12/16/2010

! Remarks:
! - Partial-wave decomposition of OPE
! - Normalization of potentials must conform with eft/nneft_lsesv. See also
!   doc/physnotes.pdf
! - Units in MeV or MeV^-1

! Reference
! KBW97: Kaiser, et al., NPA 625, 758 (1997)
! EGM00: Epelbaum, Glockle, and Meissner, NPA 671, 295 (2000)
! i^L convention in partial-wave decomposition followed that of KBW97
! which, compared with EGM00, leads to opposite signs of S-D off-diagonal terms

module ope_pwd

  use nneft_type,       only: NER, NEC, PI_NE, standard_output_unit
  use nneft_phyconst,   only: PC_mN, PC_gA, PC_fpi, PC_mpi
  use util_gauleg,      only: feb_glabsc_std
  use util_spcfnc,      only: lgndr_array
  implicit none

  real(NER), parameter    :: ONEOVEREIGHTPI_pwo  = 1.0_NER/(8.0_NER*PI_NE)

  ! mshlg_pwo & wghtslg_pwo:
  ! Gaussian-Legendre mesh to be used in z integrals, where z = cos\theta
  ! Nlg_pwo : Upper bound of mshlg_pwo() & wghtslg_pwo()
  ! mshlg_inited_pwo: If .true., mshlg_pwo and wghtslg_pwo are initialized
  ! - The mesh will *not* be automatically initialized; therefore, whatever
  ! subroutine tries to use it should check the initialization flag:
  ! mshlg_inited_pwo. If mshlg_inited_pwo == .false., call initmshlg_pwo
  logical             :: mshlg_inited_pwo = .false.
  integer, parameter  :: Nlg_pwo = 100
  real(NER) :: mshlg_pwo(1:Nlg_pwo), wghtslg_pwo(1:Nlg_pwo)

  ! lgndr_table_pwo: Pre-calculated table of Legendre functions up to lgrank_pwo
  ! - lgndr_table_pwo will be initialized by initmshlg_pwo
  integer, parameter  :: lgrank_pwo = 20
  real(NER) :: lgndr_table_pwo(0:lgrank_pwo, 1:Nlg_pwo)

  ! Frozen!!!
  ! ! As a fail-safe, when p or k < EPS_pk, except for S-waves, potential value
  ! ! set to zero.
  ! real(NER), parameter, private :: EPS_pk = 1.0E-4_NER

  contains

  ! Initializing mesh points and Legendre-function table that will be used in z
  ! integrals
  subroutine initmshlg_pwo()

    integer                             :: ii
    real(NER), dimension(0:lgrank_pwo)  :: pn

    if (mshlg_inited_pwo) then
      return
    end if
    call feb_glabsc_std(Nlg_pwo, mshlg_pwo, wghtslg_pwo)
    do ii = 1, Nlg_pwo
      call lgndr_array(lgrank_pwo, mshlg_pwo(ii), pn)
      lgndr_table_pwo(0:lgrank_pwo, ii) = pn(0:lgrank_pwo)
    end do
    mshlg_inited_pwo = .true.

    ! debug
    write(standard_output_unit, '(a)')    &
      & "-->initmshlg_pwo: initialized successfully"

  end subroutine initmshlg_pwo

  ! UT = VT + (4*I - 3)WT
  ! L + s + I = odd
  ! WT_OPE = 2.0*mN/pi * 1/8pi * gA^2/(4 fpi^2) * (mpi^2 + q^2)
  pure function WT_OPE(qsqr)

    real(NER)             :: WT_OPE
    real(NER), intent(in) :: qsqr

    real(NER) :: pref
    ! overall normalization factor of the partial-wave potential, adjustable
    ! according to the form of the LS eqn it has the 1/8pi factor
    real(NER)   :: ovalnorm_pwo, PC_fpi_sqr, PC_gA_sqr, PC_mpi_sqr

    ovalnorm_pwo = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_fpi_sqr = PC_fpi * PC_fpi
    PC_gA_sqr = PC_gA * PC_gA
    PC_mpi_sqr = PC_mpi * PC_mpi

    pref = ovalnorm_pwo * 0.25_NER * PC_gA_sqr / PC_fpi_sqr
    WT_OPE = pref / (PC_mpi_sqr + qsqr)

  end function WT_OPE

  pure function WT_cmplx_OPE(qsqr)

    complex(NEC)              :: WT_cmplx_OPE
    complex(NEC), intent(in)  :: qsqr

    real(NER) :: pref
    ! overall normalization factor of the partial-wave potential, adjustable
    ! according to the form of the LS eqn it has the 1/8pi factor
    real(NER)   :: ovalnorm_pwo, PC_fpi_sqr, PC_gA_sqr, PC_mpi_sqr

    ovalnorm_pwo = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_fpi_sqr = PC_fpi * PC_fpi
    PC_gA_sqr = PC_gA * PC_gA
    PC_mpi_sqr = PC_mpi * PC_mpi

    pref = ovalnorm_pwo * 0.25_NER * PC_gA_sqr / PC_fpi_sqr
    WT_cmplx_OPE = pref / (PC_mpi_sqr + qsqr)

  end function WT_cmplx_OPE

  ! <j0j, p|OPE|j0j, k>
  ! <Lsj|
  ! j = even --> I = 1 and (4*I - 3) = 1
  ! j = odd  --> I = 0 and (4*I - 3) = -3
  ! Symmetric
  ! Vj0j(p, k) = \int dz q^2*UT(q^2) P_j(z) ,
  !   where q^2 = p^2 + k^2 - 2*p*k*z
  function OPE_j0j(j, p, k)

    real(NER)               :: OPE_j0j
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, z, ut, psqr, ksqr, pk

    if (j < 0) then
      OPE_j0j = 0.0_NER
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
      z      = mshlg_pwo(ii)
      qsqr   = psqr + ksqr - 2.0_NER*pk*z
      regsum = regsum + qsqr*WT_OPE(qsqr)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)
    end do
    if (mod(j, 2) == 1) regsum = -3.0_NER*regsum
    OPE_j0j = regsum

  end function OPE_j0j

  ! <j1j, p|OPE|j1j, k>
  ! Symmetric
  ! j = even --> I = 0 and (4*I - 3) = -3
  ! j = odd  --> I = 1 and (4*I - 3) = 1
  function OPE_j1j(j, p, k)

    real(NER)               :: OPE_j1j
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, z, pk, psqr, ksqr, pref

    if (j < 1) then
      OPE_j1j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr   = p*p
    ksqr   = k*k
    pk     = p*k
    do ii = 1, Nlg_pwo
      z = mshlg_pwo(ii)
      qsqr    = psqr + ksqr - 2.0_NER*pk*z
      regsum  = regsum + ((psqr + 2.0_NER*pk*z + ksqr)*lgndr_table_pwo(j, ii) &
      & - 2.0_NER*pk*(lgndr_table_pwo(j-1, ii)+lgndr_table_pwo(j+1, ii)))     &
      & *WT_OPE(qsqr)*wghtslg_pwo(ii)
    end do

    if (mod(j, 2) == 1) then
      OPE_j1j = -regsum
    else
      OPE_j1j = 3.0_NER*regsum
    endif

  end function OPE_j1j

  ! For jpp, jpm, and jmm:
  ! j = even --> I = 1 and (4*I - 3) = 1
  ! j = odd  --> I = 0 and (4*I - 3) = -3

  ! <(j+1)1j, p|OPE|(j+1)1j, k>
  ! Symmetric
  function OPE_jpp(j, p, k)

    real(NER)               :: OPE_jpp
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, z, pk, psqr, ksqr, pref, pref1, pref2

    if (j < 0) then
      OPE_jpp = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    if (mod(j, 2) == 1) then
      pref = -3.0_NER
    else
      pref = 1.0_NER
    end if

    pref1 = 2.0_NER/(2.0_NER*j+1.0_NER)
    pref2 = 1.0_NER/(2.0_NER*j+1.0_NER)
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    do ii = 1, Nlg_pwo
      z = mshlg_pwo(ii)
      qsqr = psqr + ksqr - 2.0_NER*pk*z
      regsum = regsum + (pref1*lgndr_table_pwo(j, ii)*pk               &
      & - pref2*lgndr_table_pwo(j+1, ii)*(psqr + ksqr))*WT_OPE(qsqr)   &
      & *wghtslg_pwo(ii)
    end do

    OPE_jpp = -pref*regsum

  end function OPE_jpp

  ! <(j-1)1j, p|OPE|(j-1)1j, k>
  function OPE_jmm(j, p, k)

    real(NER)               :: OPE_jmm
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, z, pk, psqr, ksqr, pref, pref1, pref2
    real(NER), parameter :: minus_two_thirds = -2.0_NER/3.0_NER,      &
    & one_third = 1.0_NER/3.0_NER

    if (j < 1) then
      OPE_jmm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    if (mod(j, 2) == 1) then
      pref = -3.0_NER
    else
      pref = 1.0_NER
    end if

    psqr    = p*p
    ksqr    = k*k
    pk      = p*k
    if (j == 1) then
      do ii = 1, Nlg_pwo
        z       = mshlg_pwo(ii)
        qsqr    = psqr + ksqr - 2.0_NER*pk*z
        regsum  = regsum + (minus_two_thirds*z*pk + one_third*(psqr + ksqr))  &
        & *WT_OPE(qsqr)*wghtslg_pwo(ii)
      end do
    else
      pref1 = 2.0_NER/(2.0_NER*j+1.0_NER)
      pref2 = 1.0_NER/(2.0_NER*j+1.0_NER)
      do ii = 1, Nlg_pwo
        z       = mshlg_pwo(ii)
        qsqr    = psqr + ksqr - 2.0_NER*pk*z
        regsum  = regsum + (-pref1*lgndr_table_pwo(j, ii)*pk                  &
        & + pref2*lgndr_table_pwo(j-1, ii)*(psqr + ksqr))*WT_OPE(qsqr)        &
        & *wghtslg_pwo(ii)
      end do
    end if

    OPE_jmm = -pref*regsum

  end function OPE_jmm

  ! <(j+1)1j, p|OPE|(j-1)1j, k>  I = 1 when j = even, I = 0 when j = odd
  function OPE_jpm(j, p, k)

    real(NER)               :: OPE_jpm
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum, qsqr, z, pk, psqr, ksqr, pref

    if (j < 1) then
      OPE_jpm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    if (mod(j, 2) == 1) then
      pref = -3.0_NER
    else
      pref = 1.0_NER
    end if
    pref = pref*sqrt(j*(j+1.0_NER))/(2.0_NER*j+1.0_NER)
    psqr    = p*p
    ksqr    = k*k
    pk      = p*k
    do ii = 1, Nlg_pwo
      z       = mshlg_pwo(ii)
      qsqr    = psqr + ksqr - 2.0_NER*pk*z
      regsum  = regsum  &
        & + (-4.0_NER*pk*lgndr_table_pwo(j, ii)         &
        & + 2.0_NER*psqr*lgndr_table_pwo(j-1, ii)       &
        & + 2.0_NER*ksqr*lgndr_table_pwo(j+1, ii))*WT_OPE(qsqr)*wghtslg_pwo(ii)
    end do
    OPE_jpm = pref*regsum

  end function OPE_jpm

  function yukawa_1s0_mpisqr(mpisqr, p, k)

    real(NER)               :: yukawa_1s0_mpisqr
    real(NER), intent(in)   :: p, k, mpisqr

    integer     :: ii
    real(NER)   :: regsum, qsqr, z, pref
    real(NER)   :: ovalnorm_pwo, PC_fpi_sqr, PC_gA_sqr

    ovalnorm_pwo = 2.0_NER*PC_mN/PI_NE * ONEOVEREIGHTPI_pwo
    PC_fpi_sqr = PC_fpi * PC_fpi
    PC_gA_sqr = PC_gA * PC_gA
    pref = ovalnorm_pwo * mpisqr * 0.25_NER * PC_gA_sqr / PC_fpi_sqr
    regsum = 0.0_NER
    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
      qsqr   = p*p + k*k - 2.0_NER*p*k*z
      regsum = regsum - wghtslg_pwo(ii)/(mpisqr + qsqr)
    end do

    yukawa_1s0_mpisqr = regsum*pref

  end function yukawa_1s0_mpisqr

  function yukawa_1S0(p, k)

    real(NER)               :: yukawa_1S0
    real(NER), intent(in)   :: p, k

    integer     :: ii
    real(NER)   :: regsum, qsqr, z
    real(NER)   :: PC_mpi_sqr

    PC_mpi_sqr = PC_mpi * PC_mpi
    regsum = 0.0_NER
    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
      qsqr   = p*p + k*k - 2.0_NER*p*k*z
      regsum = regsum - PC_mpi_sqr*WT_OPE(qsqr)*wghtslg_pwo(ii)
    end do

    yukawa_1S0 = regsum

  end function yukawa_1S0


!===================================================================
!
! Complex valued OPE PW decomposition
!
!===================================================================

  ! Symmetric even for complex p, k
  ! Vj0j(p, k) = \int dz q^2*UT(q^2) P_j(z) ,
  !   where q^2 = p^2 + k^2 - 2*p*k*z
  function OPE_cmplx_j0j(j, p, k)

    complex(NEC)            :: OPE_cmplx_j0j
    complex(NEC), intent(in):: p, k
    integer, intent(in)     :: j

    integer     :: ii
    complex(NEC):: regsum, qsqr, z, ut, psqr, ksqr, pk

    if (j < 0) then
      OPE_cmplx_j0j = 0.0_NER
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
      z      = mshlg_pwo(ii)
      qsqr   = psqr + ksqr - 2.0_NER*pk*z
      regsum = regsum + qsqr*WT_cmplx_OPE(qsqr)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)
    end do
    if (mod(j, 2) == 1) regsum = -3.0_NER*regsum
    OPE_cmplx_j0j = regsum

  end function OPE_cmplx_j0j

  function OPE_cmplx_j1j(j, p, k)

    complex(NEC)            :: OPE_cmplx_j1j
    complex(NEC), intent(in):: p, k
    integer, intent(in)     :: j

    integer     :: ii
    complex(NEC):: regsum, qsqr, z, pk, psqr, ksqr, pref

    if (j < 1) then
      OPE_cmplx_j1j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr   = p*p
    ksqr   = k*k
    pk     = p*k
    do ii = 1, Nlg_pwo
      z = mshlg_pwo(ii)
      qsqr    = psqr + ksqr - 2.0_NER*pk*z
      regsum  = regsum + ((psqr + 2.0_NER*pk*z + ksqr)*lgndr_table_pwo(j, ii) &
      & - 2.0_NER*pk*(lgndr_table_pwo(j-1, ii)+lgndr_table_pwo(j+1, ii)))     &
      & *WT_cmplx_OPE(qsqr)*wghtslg_pwo(ii)
    end do

    if (mod(j, 2) == 1) then
      OPE_cmplx_j1j = -regsum
    else
      OPE_cmplx_j1j = 3.0_NER*regsum
    endif

  end function OPE_cmplx_j1j

  function OPE_cmplx_jpp(j, p, k)

    complex(NEC)            :: OPE_cmplx_jpp
    complex(NEC), intent(in):: p, k
    integer, intent(in)     :: j

    integer     :: ii
    complex(NEC):: regsum, qsqr, z, pk, psqr, ksqr, pref, pref1, pref2

    if (j < 0) then
      OPE_cmplx_jpp = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    if (mod(j, 2) == 1) then
      pref = -3.0_NER
    else
      pref = 1.0_NER
    end if

    pref1 = 2.0_NER/(2.0_NER*j+1.0_NER)
    pref2 = 1.0_NER/(2.0_NER*j+1.0_NER)
    psqr = p*p
    ksqr = k*k
    pk   = p*k
    do ii = 1, Nlg_pwo
      z = mshlg_pwo(ii)
      qsqr = psqr + ksqr - 2.0_NER*pk*z
      regsum = regsum + (pref1*lgndr_table_pwo(j, ii)*pk               &
      & - pref2*lgndr_table_pwo(j+1, ii)*(psqr + ksqr))*WT_cmplx_OPE(qsqr)   &
      & *wghtslg_pwo(ii)
    end do

    OPE_cmplx_jpp = -pref*regsum

  end function OPE_cmplx_jpp

  function OPE_cmplx_jmm(j, p, k)

    complex(NEC)            :: OPE_cmplx_jmm
    complex(NEC), intent(in):: p, k
    integer, intent(in)     :: j

    integer     :: ii
    complex(NEC):: regsum, qsqr, z, pk, psqr, ksqr, pref, pref1, pref2
    real(NER), parameter :: minus_two_thirds = -2.0_NER/3.0_NER,      &
    & one_third = 1.0_NER/3.0_NER

    if (j < 1) then
      OPE_cmplx_jmm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    if (mod(j, 2) == 1) then
      pref = -3.0_NER
    else
      pref = 1.0_NER
    end if

    psqr    = p*p
    ksqr    = k*k
    pk      = p*k
    if (j == 1) then
      do ii = 1, Nlg_pwo
        z       = mshlg_pwo(ii)
        qsqr    = psqr + ksqr - 2.0_NER*pk*z
        regsum  = regsum + (minus_two_thirds*z*pk + one_third*(psqr + ksqr))  &
        & *WT_cmplx_OPE(qsqr)*wghtslg_pwo(ii)
      end do
    else
      pref1 = 2.0_NER/(2.0_NER*j+1.0_NER)
      pref2 = 1.0_NER/(2.0_NER*j+1.0_NER)
      do ii = 1, Nlg_pwo
        z       = mshlg_pwo(ii)
        qsqr    = psqr + ksqr - 2.0_NER*pk*z
        regsum  = regsum + (-pref1*lgndr_table_pwo(j, ii)*pk                  &
        & + pref2*lgndr_table_pwo(j-1, ii)*(psqr + ksqr))*WT_cmplx_OPE(qsqr)        &
        & *wghtslg_pwo(ii)
      end do
    end if

    OPE_cmplx_jmm = -pref*regsum

  end function OPE_cmplx_jmm

  function OPE_cmplx_jpm(j, p, k)

    complex(NEC)            :: OPE_cmplx_jpm
    complex(NEC), intent(in):: p, k
    integer, intent(in)     :: j

    integer     :: ii
    complex(NEC):: regsum, qsqr, z, pk, psqr, ksqr, pref

    if (j < 1) then
      OPE_cmplx_jpm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    if (mod(j, 2) == 1) then
      pref = -3.0_NER
    else
      pref = 1.0_NER
    end if
    pref = pref*sqrt(j*(j+1.0_NER))/(2.0_NER*j+1.0_NER)
    psqr    = p*p
    ksqr    = k*k
    pk      = p*k
    do ii = 1, Nlg_pwo
      z       = mshlg_pwo(ii)
      qsqr    = psqr + ksqr - 2.0_NER*pk*z
      regsum  = regsum  &
        & + (-4.0_NER*pk*lgndr_table_pwo(j, ii)         &
        & + 2.0_NER*psqr*lgndr_table_pwo(j-1, ii)       &
        & + 2.0_NER*ksqr*lgndr_table_pwo(j+1, ii))*WT_cmplx_OPE(qsqr)*wghtslg_pwo(ii)
    end do
    OPE_cmplx_jpm = pref*regsum

  end function OPE_cmplx_jpm

end module ope_pwd
