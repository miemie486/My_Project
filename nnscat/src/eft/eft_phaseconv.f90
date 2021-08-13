! Bingwei Long 07/07/2012
! Conversion of T-matrix to phase shifts and vice versa

module eft_phaseconv

  use nneft_type
  implicit none

  ! d1 = theta_(j-1); d2 = theta_(j+1); e is the mixing angle

  type triplet_phs_cpld
    real(NER) :: d1
    real(NER) :: d2
    real(NER) :: e
  end type triplet_phs_cpld

  integer, parameter :: UNITARITY_SAFE = 0, UNITARITY_BROKEN = 1
  real(NER), parameter, private ::    &
    & SMATRIX_UNITARITY_VIOLATION_LIMIT = 1.0E-3_NER
  character(len=32), parameter, private :: phase_format_string = 'f15.10'

  interface write_phsn_to_text
    module procedure write_phsn_to_text_sngl, write_phsn_to_text_cpld
  end interface write_phsn_to_text
  interface t2d1d2e_np
    module procedure t2d1d2e_np_series, t2d1d2e_np_struct
  end interface t2d1d2e_np
  interface get_nth_delta_cpld
    module procedure get_nth_delta_cpld_series, get_nth_delta_cpld_struct
  end interface get_nth_delta_cpld
  interface write_incremental_phsn_to_text
    module procedure write_incremental_phsn_to_text_cpld, write_incremental_phsn_to_text_sngl
  end interface write_incremental_phsn_to_text

contains


  !------------------------------------------------------------------------|
  ! Utility subroutines that converts phase shifts of the uncoupled and    |
  ! coupled channels to strings.                                           |
  !------------------------------------------------------------------------|

  ! Print a series of phases in one line of text

  subroutine write_phsn_to_text_sngl(num, phsn, outstr)

    integer, intent(in)                 :: num
    real(NER), dimension(:), intent(in) :: phsn
    character(len = *), intent(out)     :: outstr

    integer             :: ii
    character(len = 512):: buff
    character(len = 128):: format_str

    outstr = ''
    format_str = '(2x, '//phase_format_string//')'
    do ii = 1, num
      write (buff, format_str) phsn(ii)
      outstr = trim(outstr)//trim(buff)
    end do

  end subroutine write_phsn_to_text_sngl

  subroutine write_phsn_to_text_cpld(num, phsn, outstr)

    integer, intent(in)                                 :: num
    type(triplet_phs_cpld), dimension(:), intent(in)    :: phsn
    character(len = *), intent(out)                     :: outstr

    integer             :: ii
    character(len = 512):: buff
    character(len = 128):: format_str

    outstr = ''
    format_str = '(3(2x, '//phase_format_string//'))'
    do ii = 1, num
      write (buff, format_str) phsn(ii)%d1, phsn(ii)%d2, phsn(ii)%e
      outstr = trim(outstr)//trim(buff)
    end do

  end subroutine write_phsn_to_text_cpld


  ! Input array phsn are expressed as LO, NLO correction, NNLO correction, and
  ! so on. The result text, however, is total phases at each order.

  subroutine write_incremental_phsn_to_text_sngl(num, phsn, outstr)

    integer, intent(in)                 :: num
    real(NER), dimension(:), intent(in) :: phsn
    character(len = *), intent(out)     :: outstr

    integer             :: ii
    real(NER), dimension(1:num) :: tmp_phs

    tmp_phs(1) = phsn(1)
    do ii = 2, num
      tmp_phs(ii) = tmp_phs(ii-1) + phsn(ii)
    end do
    call write_phsn_to_text_sngl(num, tmp_phs, outstr)

  end subroutine write_incremental_phsn_to_text_sngl

  subroutine write_incremental_phsn_to_text_cpld(num, phsn, outstr)

    integer, intent(in)                                 :: num
    type(triplet_phs_cpld), dimension(:), intent(in)    :: phsn
    character(len = *), intent(out)                     :: outstr

    integer             :: ii
    type(triplet_phs_cpld), dimension(1:num)    :: tmp_triplet

    tmp_triplet(1)%d1 = phsn(1)%d1
    tmp_triplet(1)%d2 = phsn(1)%d2
    tmp_triplet(1)%e  = phsn(1)%e
    do ii = 2, num
      tmp_triplet(ii)%d1 = tmp_triplet(ii-1)%d1 + phsn(ii)%d1
      tmp_triplet(ii)%d2 = tmp_triplet(ii-1)%d2 + phsn(ii)%d2
      tmp_triplet(ii)%e  = tmp_triplet(ii-1)%e  + phsn(ii)%e
    end do
    call write_phsn_to_text_cpld(num, tmp_triplet, outstr)

  end subroutine write_incremental_phsn_to_text_cpld

  subroutine trplt_phs_to_array(d12e, phsa)

    type(triplet_phs_cpld), intent(in)  :: d12e
    real(NER), dimension(:), intent(out):: phsa

    phsa(1) = d12e%d1
    phsa(2) = d12e%d2
    phsa(3) = d12e%e

  end subroutine trplt_phs_to_array

  subroutine array_to_trplt_phs(phsa, d12e)

    real(NER), dimension(:), intent(in) :: phsa
    type(triplet_phs_cpld), intent(out) :: d12e

    d12e%d1 = phsa(1)
    d12e%d2 = phsa(2)
    d12e%e  = phsa(3)

  end subroutine array_to_trplt_phs

  pure function factorial(n)

    integer             :: factorial
    integer, intent(in) :: n

    integer :: ii

    factorial = 1
    if (n .le. 1) return
    select case (n)
      case (2)
        factorial = 2
      case (3)
        factorial = 6
      case (4)
        factorial = 24
      case (5)
        factorial = 120
      case default
        factorial = 120
        do ii = 6, n
          factorial = factorial * ii
        end do
    end select

  end function factorial

! (ix)^n/n!, the n-th term of Taylor expansion of e^(ix)

  pure function expi_taylor(n, x)

    complex(NEC)            :: expi_taylor
    integer,    intent(in)  :: n
    real(NER),  intent(in)  :: x

    integer :: ii

    expi_taylor = 1.0_NER
    if (n .le. 0) return
    select case (n)
      case (1)
        expi_taylor = IMGUNT_NE*x
      case (2)
        expi_taylor = -x*x*0.5_NER
      case default
        expi_taylor = -x*x*0.5_NER
        do ii = 3, n
          expi_taylor = expi_taylor*(IMGUNT_NE*x)/ii
        end do
    end select

  end function expi_taylor

! the n-th term of Taylor expansion of cos(x)

  pure function cos_taylor(n, x)

    real(NEC)               :: cos_taylor
    integer,    intent(in)  :: n
    real(NER),  intent(in)  :: x

    if (mod(n, 2) .eq. 1) then
      cos_taylor = 0.0_NER
      return
    end if
    cos_taylor = x**n/factorial(n)
    if (mod(n, 4) .eq. 2) cos_taylor = -cos_taylor

  end function cos_taylor

! the n-th term of Taylor expansion of sin(x)

  pure function sin_taylor(n, x)

    real(NEC)               :: sin_taylor
    integer,    intent(in)  :: n
    real(NER),  intent(in)  :: x

    if (mod(n, 2) .eq. 0) then
      sin_taylor = 0.0_NER
      return
    end if
    sin_taylor = x**n/factorial(n)
    if (mod(n, 4) .eq. 3) sin_taylor = -sin_taylor

  end function sin_taylor


! exp(i(x1 + x2 + ...)) = 1 + f1 + f2 + ...
! f1 = ix1; f2 = -1/2(x1)^2 + ix2, f3 = ...

  subroutine get_nth_expi_prdt_subtracted(n, xk, fn)

    integer,                 intent(in)  :: n
    real(NER), dimension(:), intent(in)  :: xk
    complex(NEC),            intent(out) :: fn

    complex(NEC), dimension(n, n) :: s_n
    integer :: ii, jj

    do ii = 1, n-1
      do jj = 1, n
        if (ii*jj .le. n) s_n(ii, jj) = expi_taylor(jj, xk(ii))
      end do
    end do
    select case (n)
      case (0, 1)
        fn = 0.0_NER
      case (2)
        fn = s_n(1, 2)
      case (3)
        fn = s_n(1, 1)*s_n(2, 1) + s_n(1, 3)
      case (4)
        fn = s_n(2, 2) + s_n(1, 1)*s_n(3, 1) + s_n(1, 2)*s_n(2, 1) + s_n(1, 4)
      case (5)
        fn = s_n(1, 1)*s_n(4, 1) + s_n(2, 1)*s_n(3, 1) + s_n(1, 2)*s_n(3, 1) + s_n(1, 1)*s_n(2, 2) + s_n(1, 3)*s_n(2, 1) + s_n(1, 5)
      case (6)
        fn = s_n(1, 1)*s_n(5, 1) + s_n(2, 1)*s_n(4, 1) + s_n(1, 2)*s_n(4, 1) + s_n(3, 2) + s_n(1, 1)*s_n(2, 1)*s_n(3, 1) + s_n(1, 3)*s_n(3, 1) + s_n(2, 3) + s_n(1, 2)*s_n(2, 2) + s_n(1, 4)*s_n(2, 1) + s_n(1, 6)
      case default
        fn = 0.0_NER
        write (standard_error_unit, '(a)') 'get_nth_expi_prdt_subtracted: n is too high. fn is set zero.'
    end select

  end subroutine  get_nth_expi_prdt_subtracted

! s = exp(2*i delta/180*pi)

  subroutine get_s_from_delta_sngl(delta, s)

    real(NER),      intent(in)  :: delta
    complex(NEC),   intent(out) :: s

    s = exp(IMGUNT_NE*delta*PIONINETY_NE)

  end subroutine get_s_from_delta_sngl

! phs(1) = d1, phs(2) = d2, phs(3) = e
! phases are in degrees
! s = { cos(2e)*e^(2id1), isin(2e)*e^(i(d1+d2)) }
!     { isin(2e)*e^(i(d1+d2)), cos(2e)*e^(2id2) }

  subroutine get_s_from_delta_cpld(delta, sab)

    real(NER),      dimension(:),       intent(in)  :: delta
    complex(NEC),   dimension(2, 2),    intent(out) :: sab

    real(NER) :: cos_2e, sin_2e

    cos_2e = cos(delta(3)*PIONINETY_NE)
    sin_2e = sin(delta(3)*PIONINETY_NE)
    sab(1, 1) = cos_2e*exp(IMGUNT_NE*delta(1)*PIONINETY_NE)
    sab(2, 2) = cos_2e*exp(IMGUNT_NE*delta(2)*PIONINETY_NE)
    sab(1, 2) = IMGUNT_NE*sin_2e*exp(IMGUNT_NE*(delta(1)+delta(2))*PIOOEZ_NE)
    sab(2, 1) = sab(1, 2)

  end subroutine get_s_from_delta_cpld

  ! violation = s^+ s - 1, where s^+ = s^\dagger

  subroutine get_delta_from_s_sngl(s, delta, violation, flag)

    complex(NEC), intent(in)    :: s
    real(NER),    intent(out)   :: delta
    real(NER),    intent(out)   :: violation
    integer,      intent(out)   :: flag

    real(NER), parameter :: local_EPS = SMATRIX_UNITARITY_VIOLATION_LIMIT

    violation = real(conjg(s)*s) - 1.0_NER
    if (abs(violation) .lt. local_EPS) then
      flag = UNITARITY_SAFE
    else
      flag = UNITARITY_BROKEN
    end if
    delta = real(-0.5_NER*IMGUNT_NE*log(s)) * OEZOPI_NE

  end subroutine get_delta_from_s_sngl

  subroutine get_delta_from_s_cpld(s, delta, violation, flag)

    complex(NEC), dimension(1:2, 1:2), intent(in)   :: s
    real(NER),    dimension(:),        intent(out)  :: delta
    complex(NEC), dimension(1:2, 1:2), intent(out)  :: violation
    integer,                           intent(out)  :: flag

    real(NER), parameter :: local_EPS = SMATRIX_UNITARITY_VIOLATION_LIMIT
    real(NER)   :: ee, d1, d2

    violation = matmul(transpose(conjg(s)), s)
    violation(1, 1) = violation(1, 1) - 1.0_NER
    violation(2, 2) = violation(2, 2) - 1.0_NER
    if ((abs(violation(1, 1)) + abs(violation(1, 2)) + abs(violation(2, 1)) + abs(violation(2, 2))) .lt. local_EPS) then
      flag = UNITARITY_SAFE
    else
      flag = UNITARITY_BROKEN
    end if

    d1 = real(-0.5_NER*IMGUNT_NE*log(s(1, 1)/abs(s(1, 1))))
    d2 = real(-0.5_NER*IMGUNT_NE*log(s(2, 2)/abs(s(2, 2))))
    ee = 0.5_NER*asin(aimag(s(1, 2)*exp(-IMGUNT_NE*(d1+d2))))
    delta(1) = d1 * OEZOPI_NE
    delta(2) = d2 * OEZOPI_NE
    delta(3) = ee * OEZOPI_NE

  end subroutine get_delta_from_s_cpld

  ! The nth order of s, however, with the part linear to delta_n subtracted
  ! Expansion of the S-matrix and delta are expressed incrementally.
  ! s = s0 + s^(1) + s^(2) + ... + s^(n) + ...
  ! delta = delta_0 + delta_1 + ... delta_n + ...
  ! s^(n) = sbtrct_n(delta_0, delta_1, ..., delta_(n-1)) + s0*(2i delta_n)
  ! in array delta(ii): delta(1) = delta_0, delta(2) = delta_1, and so on

  subroutine get_nth_s_subtracted_sngl(n, delta, sbtrct_n)

    integer,                 intent(in)  :: n
    real(NER), dimension(:), intent(in)  :: delta
    complex(NEC),            intent(out) :: sbtrct_n

    real(NER), dimension(1:n-1) :: xn

    sbtrct_n = 0.0_NER
    if (n .le. 1) return

    xn(1:n-1) = delta(2:n)*PIONINETY_NE
    call get_nth_expi_prdt_subtracted(n, xn, sbtrct_n)
    sbtrct_n = exp(IMGUNT_NE*delta(1)*PIONINETY_NE)*sbtrct_n

  end subroutine get_nth_s_subtracted_sngl

! s^(n) = sbtrct_n(delta_0, delta_1, ..., delta_(n-1))
! + { {cos(2e_0)e^(2id1_0)*(2id1_n) - sin(2e_0)e^(2id1_0)*(2e_n), cos(2e_0)e^(i(d1_0+d2_0))*(2ie_n) - sin(2e_0)e^(i(d1_0+d2_0))*(d1_n+d2_n)}, {diag., cos(2e_0)e^(2id2_0)*(2id2_n) - sin(2e_0)e^(2id2_0)*(2e_n)}

  subroutine get_nth_s_subtracted_cpld(n, delta, sbtrct_n)

    integer,                 intent(in)  :: n
    real(NER), dimension(:), intent(in)  :: delta
    complex(NEC), dimension(1:2, 1:2), intent(out) :: sbtrct_n

    complex(NER)    :: fn1, fn2
    real(NER), dimension(1:n-1) :: xn
    real(NER), dimension(1:3*n) :: phs

    sbtrct_n = 0.0_NER
    if (n .le. 1) return

    phs(1:3*n) = delta(1:3*n) * PIOOEZ_NE
    xn(1:n-1) = 2.0_NER*(phs(4:3*n:3) + phs(6:3*n:3))
    call get_nth_expi_prdt_subtracted(n, xn, fn1)
    xn(1:n-1) = 2.0_NER*(phs(4:3*n:3) - phs(6:3*n:3))
    call get_nth_expi_prdt_subtracted(n, xn, fn2)
    sbtrct_n(1, 1) = 0.5_NER*exp(2.0_NER*IMGUNT_NE*phs(1)) * (exp(2.0_NER*IMGUNT_NE*phs(3))*fn1 + exp(-2.0_NER*IMGUNT_NE*phs(3))*fn2)

    xn(1:n-1) = phs(4:3*n:3) + phs(5:3*n:3) + 2.0_NER*phs(6:3*n:3)
    call get_nth_expi_prdt_subtracted(n, xn, fn1)
    xn(1:n-1) = phs(4:3*n:3) + phs(5:3*n:3) - 2.0_NER*phs(6:3*n:3)
    call get_nth_expi_prdt_subtracted(n, xn, fn2)
    sbtrct_n(1, 2) = 0.5_NER*exp(IMGUNT_NE*(phs(1)+phs(2))) * (exp(2.0_NER*IMGUNT_NE*phs(3))*fn1 - exp(-2.0_NER*IMGUNT_NE*phs(3))*fn2)
    sbtrct_n(2, 1) = sbtrct_n(1, 2)

    xn(1:n-1) = 2.0_NER*(phs(5:3*n:3) + phs(6:3*n:3))
    call get_nth_expi_prdt_subtracted(n, xn, fn1)
    xn(1:n-1) = 2.0_NER*(phs(5:3*n:3) - phs(6:3*n:3))
    call get_nth_expi_prdt_subtracted(n, xn, fn2)
    sbtrct_n(2, 2) = 0.5_NER*exp(2.0_NER*IMGUNT_NE*phs(2))*(exp(2.0_NER*IMGUNT_NE*phs(3))*fn1 + exp(-2.0_NER*IMGUNT_NE*phs(3))*fn2)

  end subroutine get_nth_s_subtracted_cpld

  ! Non-Perturbative conversion

  subroutine t2d1d2e_np_series(k, T, phs, flag)

    real(NER),                         intent(in)   :: k
    complex(NEC), dimension(1:2, 1:2), intent(in)   :: T
    real(NER), dimension(:),           intent(out)  :: phs
    integer,                           intent(out)  :: flag

    complex(NEC), dimension(1:2, 1:2) :: s, violation

    s = - IMGUNT_NE*PI_NE*k * T
    s(1, 1) = s(1, 1) + 1.0_NER
    s(2, 2) = s(2, 2) + 1.0_NER
    call get_delta_from_s_cpld(s, phs, violation, flag)

  end subroutine t2d1d2e_np_series

  ! Non-Perturbative conversion

  subroutine t2d1d2e_np_struct(k, T, triplet, flag)

    real(NER),                         intent(in)   :: k
    complex(NEC), dimension(1:2, 1:2), intent(in)   :: T
    type(triplet_phs_cpld),            intent(out)  :: triplet
    integer,                           intent(out)  :: flag

    real(NER), dimension(1:3)           :: phs

    call t2d1d2e_np_series(k, T, phs, flag)
    call array_to_trplt_phs(phs, triplet)

  end subroutine t2d1d2e_np_struct

  ! violation \equiv 2* (s0^+ * s1)

  subroutine get_nlo_delta_from_s1_sngl(delta0, s1, delta1, violation, flag)

    real(NER),     intent(in)   :: delta0
    complex(NEC),  intent(in)   :: s1
    real(NER),     intent(out)  :: delta1
    real(NER),     intent(out)  :: violation
    integer,       intent(out)  :: flag

    real(NER), parameter            :: local_EPS = SMATRIX_UNITARITY_VIOLATION_LIMIT
    complex(NEC)   :: s0, s0cs1

    call get_s_from_delta_sngl(delta0, s0)
    s0cs1 = conjg(s0)*s1
    violation = 2.0_NER*real(s0cs1)
    ! If unitarity is preserved perfectly in the perturbation-theory sense
    ! s0^+ * s1 will be purely imaginary
    if (abs(violation) .lt. local_EPS) then
      flag = UNITARITY_SAFE  ! (s0^+ s1) + (s1^+ s0) = 0 => unitarity satisfied
    else
      flag = UNITARITY_BROKEN
    end if
    delta1 = aimag(s0cs1)*NINETYOPI_NE

  end subroutine get_nlo_delta_from_s1_sngl

  subroutine get_s1_from_nlo_delta_sngl(delta0, delta1, s1)

    real(NER),    intent(in)   :: delta0, delta1
    complex(NEC), intent(out)  :: s1

    complex(NEC)    :: s0

    call get_s_from_delta_sngl(delta0, s0)
    s1 = s0*(IMGUNT_NE*delta1*PIONINETY_NE)

  end subroutine get_s1_from_nlo_delta_sngl

  ! violation \equiv (S0^\dagger S1 + S1^\dagger S0)

  subroutine get_nlo_delta_from_s1_cpld(delta0, s1, delta1, violation, flag)

    real(NER),    dimension(:),    intent(in)   :: delta0
    complex(NEC), dimension(2, 2), intent(in)   :: s1
    real(NER),    dimension(:),    intent(out)  :: delta1
    complex(NEC), dimension(2, 2), intent(out)  :: violation
    integer,                       intent(out)  :: flag

    real(NER), parameter            :: local_EPS = SMATRIX_UNITARITY_VIOLATION_LIMIT
    complex(NEC), dimension(2, 2)   :: s0, tmp
    real(NER)    :: cse0, r1, r2, e1
    complex(NEC) :: x, ex2d1c, ex2d2c, exd1d2c

    call get_s_from_delta_cpld(delta0, s0)
    ! violation = matmul(transpose(conjg(s0)), s1)
    tmp = matmul(transpose(conjg(s0)), s1)
    violation = tmp + transpose(conjg(tmp))
    ! if unitarity is preserved perfectly in the perturbation-theory sense
    ! s0^+ s1 will be purely anti-hermitian
    if ((abs(violation(1, 1)) + abs(violation(1, 2)) + abs(violation(2, 1)) + abs(violation(2, 2))) .lt. local_EPS) then
      flag = UNITARITY_SAFE  ! (s0^+ s1) + (s1^+ s0) = 0 => unitarity satisfied
    else
      flag = UNITARITY_BROKEN
    end if

    cse0   = cos(delta0(3)*PIONINETY_NE)
    ex2d1c = exp(-IMGUNT_NE*delta0(1)*PIONINETY_NE)
    ex2d2c = exp(-IMGUNT_NE*delta0(2)*PIONINETY_NE)
    exd1d2c = exp(-IMGUNT_NE*(delta0(1)+delta0(2))*PIOOEZ_NE)

    r1   = aimag(s1(1, 1)*ex2d1c)/cse0       ! 2*d_1^(1) in radian
    r2   = aimag(s1(2, 2)*ex2d2c)/cse0       ! 2*d_2^(1) in radian
    e1   = aimag(s1(1, 2)*exd1d2c)/cse0      ! 2*epsilon^(1)  in radian

    delta1(1) = r1 * NINETYOPI_NE
    delta1(2) = r2 * NINETYOPI_NE
    delta1(3) = e1 * NINETYOPI_NE

  end subroutine get_nlo_delta_from_s1_cpld

  subroutine get_s1_from_nlo_delta_cpld(delta0, delta1, s1)

    real(NER), dimension(:),       intent(in)   :: delta0, delta1
    complex(NEC), dimension(2, 2), intent(out)  :: s1

    real(NER)       :: r1, r2, e1, cse0, sne0
    complex(NEC)    :: x, ex2d1, ex2d2, exd1d2

    cse0  = cos(delta0(3)*PIONINETY_NE)
    sne0  = sin(delta0(3)*PIONINETY_NE)
    ex2d1 = exp(IMGUNT_NE*delta0(1)*PIONINETY_NE)
    ex2d2 = exp(IMGUNT_NE*delta0(2)*PIONINETY_NE)
    exd1d2 = exp(IMGUNT_NE*(delta0(1)+delta0(2))*PIOOEZ_NE)
    r1  = delta1(1)*PIOOEZ_NE
    r2  = delta1(2)*PIOOEZ_NE
    e1  = delta1(3)*PIOOEZ_NE

    s1(1, 1) = 2.0_NER*ex2d1*(-sne0*e1 + IMGUNT_NE*cse0*r1)
    s1(2, 2) = 2.0_NER*ex2d2*(-sne0*e1 + IMGUNT_NE*cse0*r2)
    s1(1, 2) = exd1d2*(-(r1+r2)*sne0 + 2.0_NER*IMGUNT_NE*cse0*e1)
    s1(2, 1) = s1(1, 2)

  end subroutine get_s1_from_nlo_delta_cpld

  subroutine get_nth_s_sngl(n, delta, snn)

    integer,                 intent(in)  :: n
    real(NER), dimension(:), intent(in)  :: delta
    complex(NEC),            intent(out) :: snn

    complex(NEC)    :: sn1

    call get_nth_s_subtracted_sngl(n, delta, snn)
    call get_s1_from_nlo_delta_sngl(delta(1), delta(n+1), sn1)
    snn = snn + sn1

  end subroutine get_nth_s_sngl

  subroutine get_nth_s_cpld(n, delta, snn)

    integer,                            intent(in)  :: n
    real(NER),    dimension(:),         intent(in)  :: delta
    complex(NEC), dimension(1:2, 1:2),  intent(out) :: snn

    complex(NEC), dimension(1:2, 1:2) :: sn1

    call get_nth_s_subtracted_cpld(n, delta, snn)
    call get_s1_from_nlo_delta_cpld(delta(1:3), delta(3*n+1:3*n+3), sn1)
    snn = snn + sn1

  end subroutine get_nth_s_cpld

  ! kcot = ik - 2/(pi*T)

  subroutine convert_t_to_kcot_sngl(k, T, kcot)

    real(NER),               intent(in)  :: k
    complex(NEC),            intent(in)  :: T
    real(NER),               intent(out) :: kcot

    ! phs = t2delta_np(T, k, flag)
    ! kcot = k/tan(phs/180.0*PI_NE)
    kcot = real(IMGUNT_NE*k - 2.0_NER/(PI_NE*T))

  end subroutine convert_t_to_kcot_sngl

  subroutine convert_phase_to_kcot_sngl(k, phs, kcot)

    real(NER), intent(in)  :: k
    real(NER), intent(in)  :: phs
    real(NER), intent(out) :: kcot

    ! phs = t2delta_np(T, k, flag)
    kcot = k/tan(phs/180.0*PI_NE)
    ! kcot = real(IMGUNT_NE*k - 2.0_NER/(PI_NE*T))

  end subroutine convert_phase_to_kcot_sngl

  ! phsn(1) is the LO phase delta_0, phsn(2) has NLO
  ! For nonpertbative LO, T_0 has complete unitarity.
  ! For perturbative LO, T_0 (and delta_0) must be zero.
  ! violation is amount of unitarity violation, defined as

  subroutine get_nth_delta_sngl(n, T_n, phsn, k, flag, violation, delta_n)

    integer,   intent(in)     :: n
    complex(NEC), intent(in)  :: T_n
    real(NER), intent(in)     :: phsn(:)
    real(NER),    intent(in)  :: k
    integer,      intent(out) :: flag
    real(NER),    intent(out) :: violation, delta_n

    real(NER)       :: local_EPS = SMATRIX_UNITARITY_VIOLATION_LIMIT
    complex(NEC)    :: snn, sbtrct_n, sn1

    if (n <= 0) then
      write (standard_error_unit, '(a)')                  &
        & 'get_nth_delta_sngl: n must be greater than 0.'
      flag = UNITARITY_BROKEN
    else
      ! snn = -i*k*pi* Tn
      snn = -IMGUNT_NE*k*PI_NE*T_n
      if (n == 1) then
        sn1 = snn
      else
        call get_nth_s_subtracted_sngl(n, phsn(1:n), sbtrct_n)
        ! snn = sbtrct_n + sn1
        ! sn1 = s0^+ * (2*i*delta_n)
        sn1 = snn - sbtrct_n
      end if
      call get_nlo_delta_from_s1_sngl(phsn(1), sn1, delta_n, violation, flag)
    end if

  end subroutine get_nth_delta_sngl

  subroutine get_nth_delta_cpld_series(n, T_n, phsn, k, flag, violation, delta_n)

    integer,                            intent(in)  :: n
    complex(NEC), dimension(1:2, 1:2),  intent(in)  :: T_n
    real(NER),    dimension(:),         intent(in)  :: phsn
    real(NER),                          intent(in)  :: k
    integer,                            intent(out) :: flag
    complex(NEC), dimension(1:2, 1:2),  intent(out) :: violation
    real(NER),    dimension(:),         intent(out) :: delta_n

    real(NER)                           :: local_EPS = SMATRIX_UNITARITY_VIOLATION_LIMIT
    complex(NEC), dimension(1:2, 1:2)   :: snn, sbtrct_n, sn1

    if (n <= 0) then
      write (standard_error_unit, '(a)')                  &
        & 'get_nth_delta_sngl: n must be greater than 0.'
      flag = UNITARITY_BROKEN
    else
      ! snn = -i*k*pi* Tn
      snn = -IMGUNT_NE*k*PI_NE*T_n
      if (n .eq. 1) then
        sn1 = snn
      else
        call get_nth_s_subtracted_cpld(n, phsn(1:3*n), sbtrct_n)
        ! snn = sbtrct_n + sn1
        sn1 = snn - sbtrct_n
      end if
      call get_nlo_delta_from_s1_cpld(phsn(1:3), sn1, delta_n, violation, flag)
    end if

  end subroutine get_nth_delta_cpld_series

  subroutine get_nth_delta_cpld_struct(n, T_n, trplt_array, k, flag, violation, trplt_n)

    integer,                            intent(in)  :: n
    complex(NEC), dimension(1:2, 1:2),  intent(in)  :: T_n
    type(triplet_phs_cpld), dimension(:), intent(in)  ::  trplt_array
    real(NER),                          intent(in)  :: k
    integer,                            intent(out) :: flag
    complex(NEC), dimension(1:2, 1:2),  intent(out) :: violation
    type(triplet_phs_cpld), intent(out)             :: trplt_n

    integer :: ii
    real(NER), dimension(1:3*n) :: lower_order_phs
    real(NER), dimension(1:3)   :: phs_n

    do ii = 1, n
      lower_order_phs(3*ii-2) = trplt_array(ii)%d1
      lower_order_phs(3*ii-1) = trplt_array(ii)%d2
      lower_order_phs(3*ii) = trplt_array(ii)%e
    end do
    call get_nth_delta_cpld_series(n, T_n, lower_order_phs, k, flag, violation, phs_n)
    call array_to_trplt_phs(phs_n, trplt_n)

  end subroutine get_nth_delta_cpld_struct

  ! e^(2i delta) = 1 - i*\pi*k*T

  function t2delta_np(T, k, violation, flag)

    real(NER)                   :: t2delta_np
    complex(NEC), intent(in)    :: T
    real(NER), intent(in)       :: k
    real(NER), intent(out)      :: violation
    integer, intent(out)        :: flag

    complex(NEC)    :: s

    s = -k*PI_NE*IMGUNT_NE * T + 1.0_NER
    call get_delta_from_s_sngl(s, t2delta_np, violation, flag)

  end function t2delta_np

  subroutine get_delta_np_sngl(T, k, violation, flag, delta)

    complex(NEC), intent(in)    :: T
    real(NER), intent(in)       :: k
    real(NER), intent(out)      :: violation, delta
    integer, intent(out)        :: flag

    complex(NEC)    :: s

    s = -k*PI_NE*IMGUNT_NE * T + 1.0_NER
    call get_delta_from_s_sngl(s, delta, violation, flag)

  end subroutine get_delta_np_sngl


  ! T1 = -delta1*exp(2idelta0)*2/(pi*k)
  ! delta1 and delta0 are in degrees

  function delta2t_nlo(delta1, delta0, k)
    complex(NEC)           :: delta2t_nlo
    real(NER), intent(in)  :: delta1, delta0, k

    delta2t_nlo = -delta1/90.0_NER * exp(IMGUNT_NE * delta0*PIONINETY_NE) / k

  end function delta2t_nlo

end module eft_phaseconv
