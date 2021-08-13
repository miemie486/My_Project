! Bingwei Long
! 06/26/2013
! 05/19/2012
! 04/09/2012
! 06/18/2011

! Long-range potentials are "rescaled" in the sense that no explicit dependence
! on mN appears in the LS eqn (see nneft_lsesv.f90)
! For coupled channels, the quantum number with lower angular momentum is used.
! e.g., 3S1, 3F2, etc.
module eft_potspwd

  use nneft_type
  use nneft_phyconst
  use ope_pwd
  use vtpe0_pwd
  use vtpe1_pwd
  use rel_ope_pwd
  implicit none


  integer, parameter ::       &
    & CHNID_UNDEFINED = 0,    &
    & CHNID_1S0_PP = 100,     &
    & CHNID_3S1_PP = 301,     &
    & CHNID_1P1_PP = 111,     &
    & CHNID_3P0_PP = 310,     &
    & CHNID_3P1_PP = 311,     &
    & CHNID_3P2_PP = 312,     &
    & CHNID_1D2_PP = 122,     &
    ! & CHNID_3D1_PP = 321,      &
    & CHNID_3D2_PP = 322,     &
    & CHNID_3D3_PP = 323,     &
    & CHNID_1F3_PP = 133,     &
    ! & CHNID_3F2_PP = 332,     &
    & CHNID_3F3_PP = 333,     &
    & CHNID_3F4_PP = 334,     &
    & CHNID_1G4_PP = 144,     &
    ! & CHNID_3G3_PP = 343,     &
    & CHNID_3G4_PP = 344,     &
    & CHNID_3G5_PP = 345

  integer, parameter ::       &
    & REGTYPE_GAUSSIAN = 1,   &
    & REGTYPE_SHARP    = 2,   &
    & REGTYPE_NONE     = 100

  integer, parameter :: CHTYPE_SNGL = 10, CHTYPE_CPLD = 20
  integer, parameter :: ALGID_1S0_DIBARYON = 110, ALGID_1S0_NORMAL = 120
  integer, parameter :: POTS_DWAVE_PERT = 100, POTS_WPC = 200, POTS_J2_PERT = 300
  character(len = 20), parameter :: orb_L_id = 'spdfghijklmnopqrstuv'

  integer, parameter               ::       &
    & PHASE_INCREMENTAL = 10,               &
    & PHASE_TOTAL       = 20

  integer, parameter      ::                &
    & OUTOPT_PHASE      = 10,               &
    & OUTOPT_B0B2       = 20

  integer, parameter  ::    &
    & MAX_NCH           = 2

  integer, parameter               ::       &
    & MAX_NUMBER_POWFAC = 10,               &
    & MAX_NUMBER_ARRAY  = 10,               &
    & MAX_NUMBER_ORDER  = 10,               &
    & MAX_NUMBER_INPUTS = 40,               &
    & DEFAULT_N         = 100,              &
    & DEFAULT_REGTYPE   = REGTYPE_GAUSSIAN, &
    & DEFAULT_UPTOQN    = 0,                &
    & MAX_NUMBER_CS     = 30

  real(NER), parameter             ::       &
    & DEFAULT_LAMBDA = 1400.0_NER,          &
    & DEFAULT_MAMBDA = 800.0_NER,           &
    & DEFAULT_k      = 10.0_NER

  ! the contact part of OPE_1S0
  ! real(NER) :: PC_CgA = PC_gA*PC_gA*PC_mN/(8.0*PI_NE*PI_NE*PC_fpi*PC_fpi)

  type :: lsj_symbol

    integer :: L, S, J

  end type lsj_symbol

  ! option structure passed to functions that generate potential elements

  type :: potential_options

    integer          :: regtype = REGTYPE_GAUSSIAN, pots_type = POTS_WPC, chnid = CHNID_UNDEFINED
    real(NER)        :: Mambda
    type(lsj_symbol) :: lsj

  end type potential_options


  contains

  ! Convert threshold T(p, 0; 0)/0 to \alpha in fm^3
  ! both T and p are in MeV
  function t2alpha_fm3(T, k)

    real(NER)                   :: t2alpha_fm3
    complex(NEC), intent(in)    :: T
    real(NER), intent(in)       :: k

    t2alpha_fm3 = real(T) / (k*k) * PC_mN * PC_hbarc ** 3

  end function t2alpha_fm3

  ! Compute Tlab from kcm (applicable for both relativisitic and non-rel.)
  function k2Tlab(k)

    real(NER)               :: k2Tlab
    real(NER), intent(in)   :: k

    k2Tlab = 2.0_NER * k * k / PC_mN

  end function k2Tlab

  ! Compute kcm from Tlab (applicable for both relativisitic and non-rel.)
  function Tlab2k(Tl)

    real(NER)               :: Tlab2k
    real(NER), intent(in)   :: Tl

    Tlab2k = sqrt(PC_mN * Tl * 0.5_NER)

  end function Tlab2k

  ! Compute Ecm from kcm (relativistic)
  function k2Ecm(k)

    real(NER)               :: k2Ecm
    real(NER), intent(in)   :: k
    k2Ecm = 2.0_NER * (sqrt(k*k + PC_mN * PC_mN) - PC_mN)

  end function k2Ecm

  ! Compute Ecm from kcm (nonrelativistic)
  function k2Ecm_nonrel(k)

    real(NER)               :: k2Ecm_nonrel
    real(NER), intent(in)   :: k

    k2Ecm_nonrel = k*k/PC_mN

  end function k2Ecm_nonrel

  function L_from_chnid(chnid)

    integer             :: L_from_chnid
    integer, intent(in) :: chnid

    type(lsj_symbol):: lsj
    logical         :: succ_flag

    call convert_chnid_to_lsj(chnid, lsj, succ_flag)
    if (succ_flag) then
      L_from_chnid = lsj%L
    else
      L_from_chnid = 0
    end if

  end function L_from_chnid


  ! chnid is a three-digit integer -> 1(3)LJ. e.g., 3S1 corresponds to 301

  subroutine convert_lsj_to_chnid(lsj, chnid)

    type(lsj_symbol), intent(in) :: lsj
    integer, intent(out)         :: chnid

    if (.not. is_lsj_valid(lsj)) then
      chnid = CHNID_UNDEFINED
      return
    end if
    if (lsj%J .gt. 10 .or. lsj%L .gt. 10) then
      chnid = CHNID_UNDEFINED
      return
    end if
    chnid = lsj%J + lsj%L*10
    if (lsj%S .eq. 1) then
      chnid = chnid + 300
    else
      chnid = chnid + 100
    end if

  end subroutine convert_lsj_to_chnid


  subroutine convert_chnid_to_lsj(chnid, lsj, succ_flag)

    integer, intent(in)             :: chnid
    type(lsj_symbol), intent(out)   :: lsj
    logical, intent(out)            :: succ_flag

    integer :: pc

    pc = mod(chnid, 100)
    select case (int(chnid/100))
      case (1)
        lsj%S = 0
      case (3)
        lsj%S = 1
      case default
        succ_flag = .false.
        return
    end select
    lsj%L = int(pc/10)
    lsj%J = mod(pc, 10)
    succ_flag = .true.

  end subroutine convert_chnid_to_lsj

  subroutine convert_chnid_to_text(chnid, chnstr, succ_flag)

    integer, intent(in)             :: chnid
    character(len = *), intent(out) :: chnstr
    logical, intent(out)            :: succ_flag

    type(lsj_symbol)  :: lsj

    call convert_chnid_to_lsj(chnid, lsj, succ_flag)
    if (.not. succ_flag) return
    call convert_lsj_to_text(lsj, chnstr)

  end subroutine convert_chnid_to_text

  function is_lsj_valid(lsj)

    logical :: is_lsj_valid
    type(lsj_symbol), intent(in) :: lsj

    if (lsj%S .lt. 0 .or. lsj%S .gt. 1) then
      is_lsj_valid = .false.
      return
    end if
    if (lsj%L .eq. 0) then
      if (lsj%J .eq. lsj%S) then
        is_lsj_valid = .true.
      else
        is_lsj_valid = .true.
      end if
      return
    end if
    if (lsj%J .gt. lsj%L+lsj%S .or. lsj%J .lt. lsj%L-lsj%S) then
      is_lsj_valid = .false.
    else
      is_lsj_valid = .true.
    end if

  end function is_lsj_valid

  ! string ID is in lower case, e.g., 1s0, 3d3, etc.
  subroutine convert_lsj_to_text(lsj, text)

    type(lsj_symbol), intent(in)    :: lsj
    character(len = *), intent(out) :: text

    text = ''
    if (.not. is_lsj_valid(lsj)) then
      text = ''
      return
    end if
    if (lsj%J .gt. 10 .or. lsj%L .gt. 10) then
      text = ''
      return
    end if
    if (lsj%S .eq. 1) then
      text(1:1) = '3'
    else
      text(1:1) = '1'
    end if
    text(2:2) = orb_L_id(lsj%L+1:lsj%L+1)
    write (text(3:3), '(i1)') lsj%J

  end subroutine convert_lsj_to_text


  subroutine convert_text_to_lsj(text, lsj, succ_flag)

    character(len = *), intent(in)  :: text
    type(lsj_symbol), intent(out)   :: lsj
    logical, intent(out)            :: succ_flag

    integer :: tmpi, ios

    read (text(1:1), '(i1)', iostat = ios) tmpi
    if (ios /= 0) then
      succ_flag = .false.
      return
    end if
    if (tmpi .eq. 3) then
      lsj%S = 1
    else
      lsj%S = 0
    end if
    lsj%L = scan(orb_L_id, text(2:2)) - 1
    read (text(3:3), '(i1)', iostat = ios) lsj%J
    if (ios /= 0) then
      succ_flag = .false.
      return
    end if
    succ_flag = is_lsj_valid(lsj)

  end subroutine convert_text_to_lsj


  subroutine convert_text_to_chnid(text, chnid, succ_flag)

    character(len = *), intent(in)  :: text
    integer, intent(out)            :: chnid
    logical, intent(out)            :: succ_flag

    type(lsj_symbol)                :: lsj

    call convert_text_to_lsj(text, lsj, succ_flag)
    if (succ_flag) then
      call convert_lsj_to_chnid(lsj, chnid)
    else
      chnid = CHNID_UNDEFINED
    end if

  end subroutine convert_text_to_chnid


  ! Whether chnid is a coupled or single channel

  function type_of_channel_PP(chnid)

    integer :: type_of_channel_PP
    integer, intent(in) :: chnid

    logical :: succeeded
    type(lsj_symbol) :: lsj

    call convert_chnid_to_lsj(chnid, lsj, succeeded)

    if (lsj_is_coupled(lsj)) then
      type_of_channel_PP = CHTYPE_CPLD
    else
      type_of_channel_PP = CHTYPE_SNGL
    end if

  end function type_of_channel_PP

  logical function lsj_is_coupled(lsj)

    type(lsj_symbol), intent(in) :: lsj

    lsj_is_coupled = .false.
    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a, i10, 1x, i10, 1x, i10, a)')  &
        & 'lsj_is_coupled: ', lsj%L, lsj%S, lsj%J, " are not legit."
      lsj_is_coupled = .false.
    else
      if (lsj%J == 0) then
        lsj_is_coupled = .false.
      else
        if (lsj%L == lsj%J .or. lsj%J == 0) then
          lsj_is_coupled = .false.
        else
          lsj_is_coupled = .true.
        end if
      end if
    end if

  end function lsj_is_coupled

  logical function triqn_is_coupled(L, S, J)

    integer, intent(in) :: L, S, J

    type(lsj_symbol) :: lsj
    lsj%L = L
    lsj%S = S
    lsj%J = J

    triqn_is_coupled = lsj_is_coupled(lsj)

  end function triqn_is_coupled

  function regltr_half_PP(regtype, p, Mambda)

    integer,   intent(in)   :: regtype
    real(NER), intent(in)   :: p, Mambda
    real(NER)               :: regltr_half_PP

    select case (regtype)
      case (REGTYPE_GAUSSIAN)
        regltr_half_PP = regltr_gaussian_half_PP(p, Mambda)
      case (REGTYPE_SHARP)
        regltr_half_PP = 1.0
      case default
        write (standard_error_unit, '(a)')  &
          'regltr_half_PP: Warning: unknown regtype. REGTYPE_SHARP is assumed.'
        regltr_half_PP = 1.0
    end select

  end function regltr_half_PP

  function regltr_PP(regtype, p1, p2, Mambda)

    integer,   intent(in)   :: regtype
    real(NER), intent(in)   :: p1, p2, Mambda
    real(NER)               :: regltr_PP

    select case (regtype)
      case (REGTYPE_GAUSSIAN)
        regltr_PP = regltr_gaussian_PP(p1, p2, Mambda)
      case (REGTYPE_SHARP)
        regltr_PP = 1.0_NER
      case default
        write (standard_error_unit, '(a)')  &
          'regltr_PP: Warning: unknown regtype. REGTYPE_SHARP is assumed.'
        regltr_PP = 1.0_NER
    end select

  end function regltr_PP

  ! Gaussian regulator, Mambda is the cutoff parameter
  pure function regltr_gaussian_PP(p1, p2, Mambda)

    real(NER), intent(in)   :: p1, p2, Mambda
    real(NER)               :: regltr_gaussian_PP

    regltr_gaussian_PP =    &
      & exp(-(p1**GAUSSIAN_POWER + p2**GAUSSIAN_POWER)/Mambda**GAUSSIAN_POWER)

  end function regltr_gaussian_PP

  pure function regltr_gaussian_cmplx_PP(p1, p2, Mambda)

    real(NER), intent(in)   :: Mambda
    complex(NEC), intent(in):: p1, p2
    complex(NEC)            :: regltr_gaussian_cmplx_PP

    regltr_gaussian_cmplx_PP =    &
      & exp(-(p1**PC_GAUSS_CMPLX + p2**PC_GAUSS_CMPLX)/Mambda**PC_GAUSS_CMPLX)

  end function regltr_gaussian_cmplx_PP

  pure function regltr_gaussian_half_PP(p, Mambda)

    real(NER), intent(in)   :: p, Mambda
    real(NER)               :: regltr_gaussian_half_PP

    regltr_gaussian_half_PP = exp(-(p/Mambda)**GAUSSIAN_POWER)

  end function regltr_gaussian_half_PP


end module eft_potspwd
