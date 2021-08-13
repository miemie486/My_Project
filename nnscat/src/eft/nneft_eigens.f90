! Bingwei Long 08/31/2015
! Bingwei Long 11/09/2012

module nneft_eigens

  use nneft_type
  implicit none
  external sgeev, dgeev

  integer, parameter :: eigens_SUCCEEDED = 0, eigens_FAILED = -1

  abstract interface
    function realfunc_onevar_eigens(q)
      import      :: NER
      implicit none
      real(NER)   :: realfunc_onevar_eigens
      real(NER), intent(in)   :: q
    end function
    function realfunc_twovars_eigens(q, p)
      import      :: NER
      real(NER)   :: realfunc_twovars_eigens
      real(NER), intent(in)   :: q, p
    end function
  end interface

  private :: realfunc_onevar_eigens, realfunc_twovars_eigens

  contains

  ! eta*Gamma_i = arrA_ij Gamma_j
  ! Return the eta that is closest to 1
  ! This is a wrapper for LAPACK subroutines
  subroutine get_eigens_matrix_cmplx_modeigens(arrA, N, eta, Gamma, flag)

    complex(NEC), intent(in)    :: arrA(:, :)
    integer, intent(in)         :: N
    complex(NEC), intent(out)   :: eta, Gamma(:)
    integer, intent(out)        :: flag

    integer, parameter  :: buff_length_para = 8
    integer :: ii, jj, INFO, LDdummy, LDVR, LWORK, id_eigen_unit
    real(NER)   :: reg, factor, current_min, dist, RWORK(1:2*N)
    complex(NEC)   :: W(1:N), dummy(1:N, 1:N), K(1:N, 1:N), VR(1:N, 1:N), WORK(1:buff_length_para*N*N)

    LDdummy = 1
    LDVR = N
    LWORK = buff_length_para*N*N
    K(1:N, 1:N) = arrA(1:N, 1:N)

    if (NEC .eq. DPC) then
      call zgeev('N', 'V', N, K, N, W, dummy, LDdummy, VR, LDVR, WORK, LWORK, RWORK, INFO)
    else
      call cgeev('N', 'V', N, K, N, W, dummy, LDdummy, VR, LDVR, WORK, LWORK, RWORK, INFO)
    end if
    if (INFO .ne. 0) then
      flag = eigens_FAILED
      return
    end if

    id_eigen_unit = 1
    current_min = abs(W(1) - 1.0_NER)
    do ii = 2, N
      dist = abs(W(ii) - 1.0_NER)
      if (dist .lt. current_min) then
        current_min = dist
        id_eigen_unit = ii
      end if
    end do
    eta = W(id_eigen_unit)
    Gamma(1:N) = VR(1:N, id_eigen_unit)
    flag = eigens_SUCCEEDED

  end subroutine get_eigens_matrix_cmplx_modeigens

! with real type

  subroutine get_eigens_matrix_real_modeigens(arrA, N, eta, Gamma, flag)

    real(NER), intent(in)    :: arrA(:, :)
    integer, intent(in)         :: N
    real(NER), intent(out)   :: eta, Gamma(:)
    integer, intent(out)        :: flag

    integer, parameter  :: buff_length_para = 8
    integer :: ii, jj, INFO, LDdummy, LDVR, LWORK, id_eigen_unit
    real(NER)   :: reg, factor, current_min, dist, RWORK(1:2*N)
    complex(NEC)   :: W(1:N), dummy(1:N, 1:N), K(1:N, 1:N), VR(1:N, 1:N), WORK(1:buff_length_para*N*N)

    LDdummy = 1
    LDVR = N
    LWORK = buff_length_para*N*N
    K(1:N, 1:N) = arrA(1:N, 1:N)

    if (NEC .eq. DPC) then
      call zgeev('N', 'V', N, K, N, W, dummy, LDdummy, VR, LDVR, WORK, LWORK, RWORK, INFO)
    else
      call cgeev('N', 'V', N, K, N, W, dummy, LDdummy, VR, LDVR, WORK, LWORK, RWORK, INFO)
    end if
    if (INFO .ne. 0) then
      flag = eigens_FAILED
      return
    end if

    id_eigen_unit = 1
    current_min = abs(W(1) - 1.0_NER)
    do ii = 2, N
      dist = abs(W(ii) - 1.0_NER)
      if (dist .lt. current_min) then
        current_min = dist
        id_eigen_unit = ii
      end if
    end do
    eta = W(id_eigen_unit)
    Gamma(1:N) = VR(1:N, id_eigen_unit)
    flag = eigens_SUCCEEDED

  end subroutine get_eigens_matrix_real_modeigens


  ! eta*Gamma(p) = 1/f(p)*\int^\Lambda dq*funcM(p, q; B)*Gamma(q)
  subroutine get_eigens_modeigens(funcf, funcM, N, msh, wght, eta, Gamma)

    procedure(realfunc_onevar_eigens)      :: funcf
    procedure(realfunc_twovars_eigens)     :: funcM
    integer, intent(in)         :: N
    real(NER), intent(in)       :: msh(:), wght(:)
    real(NER), intent(out)      :: Gamma(:)
    complex(NEC), intent(out)   :: eta

    integer :: ii, jj, flag
    real(NER)   :: pii, pjj, vecf(1:N), matM(1:N, 1:N)

    do jj = 1, N
      pjj = msh(jj)
      vecf(jj) = funcf(pjj)
      do ii = 1, N
        pii = msh(ii)
        matM(ii, jj) = funcM(pii, pjj)
      end do
    end do
    call eigen_invf_mat_modeigens(N, vecf, matM, msh, wght, eta, Gamma, flag)
    if (flag == eigens_FAILED) then
      write(standard_error_unit, '(a)')   &
        'get_eigens_modeigens: eigen_invf_mat_modeigens failed.'
      eta = -1.0_NER
    end if

  end subroutine get_eigens_modeigens

  ! eta*Gamma(p) = 1/f(p, B)*\int^\Lambda dq*potM(p, q; B)*Gamma(q)
  ! return the eigen value, eta, that is closest to 1
  subroutine eigen_invf_mat_modeigens(N, f, potM, msh, wght, eta, Gamma, flag)

    integer, intent(in)             :: N
    real(NER), intent(in)           :: f(:), potM(:, :), msh(:), wght(:)
    complex(NEC), intent(out)       :: eta
    real(NER), intent(out)          :: Gamma(:)
    integer, intent(out)            :: flag

    integer, parameter  :: buff_length_para = 16
    integer :: ii, jj, INFO, LDdummy, LDVR, LWORK, id_eigen_unit
    real(NER)   :: current_min, dist
    real(NER)   :: WR(1:N), WI(1:N), dummy(1:N, 1:N), K(1:N, 1:N), VR(1:N, 1:N), WORK(1:buff_length_para*N*N)

    LDdummy = 1
    LDVR = N
    LWORK = buff_length_para*N*N
    do jj = 1, N
      do ii = 1, N
        K(ii, jj) = potM(ii, jj)*wght(jj)/f(ii)
      end do
    end do
    if (NER .eq. DP) then
      call dgeev('N', 'V', N, K, N, WR, WI, dummy, LDdummy, VR, LDVR, WORK, LWORK, INFO)
    else
      call sgeev('N', 'V', N, K, N, WR, WI, dummy, LDdummy, VR, LDVR, WORK, LWORK, INFO)
    end if
    if (INFO .ne. 0) then
      flag = eigens_FAILED
      return
    end if

    ! look for the eigen value that is closest to 1
    id_eigen_unit = 1
    current_min = abs(cmplx(WR(1) - 1.0_NER, WI(1), NEC))
    do ii = 2, N
      dist = abs(cmplx(WR(ii) - 1.0_NER, WI(ii), NEC))
      if (dist .lt. current_min) then
        current_min = dist
        id_eigen_unit = ii
      end if
    end do
    eta = cmplx(WR(id_eigen_unit), WI(id_eigen_unit))
    Gamma(1:N) = VR(1:N, id_eigen_unit)
    flag = eigens_SUCCEEDED

  end subroutine eigen_invf_mat_modeigens


  ! eta*Gamma(p) = 1/(-mN*B - p^2)*\int^\infty_0 dq*q^2*v(p, q; B)*Gamma(q)
  ! return the eigen value, eta, that is closest to 1
  subroutine eigen_mat_modeigens(N, B, mN, VM, msh, wght, qsqr, eta, Gamma, flag)

    integer, intent(in)             :: N
    real(NER), intent(in)           :: B, mN, VM(:, :), msh(:), wght(:), qsqr(:)
    complex(NEC), intent(out)       :: eta
    real(NER), intent(out)          :: Gamma(:)
    integer, intent(out)            :: flag

    integer     :: ii, jj
    real(NER)   :: mB, f(1:N), potM(1:N, 1:N)

    mB = B*mN
    if (mB < 0.0_NER) then
      write (standard_error_unit, '(a)') '(solve_eigen_modeigens): mB must be positive'
      flag = eigens_FAILED
      return
    end if
    do jj = 1, N
      f(jj) = 1.0_NER/(-mB - qsqr(jj))
      do ii = 1, N
        potM(ii, jj) = VM(ii, jj)*qsqr(jj)
      end do
    end do

    call eigen_invf_mat_modeigens(N, f, potM, msh, wght, eta, Gamma, flag)

  end subroutine eigen_mat_modeigens


  subroutine matrx_elmnt_eigen_modeigens(mN, N, Gamma, potM, msh, wght, elem)

    real(NER), intent(in)   :: mN
    integer, intent(in)     :: N
    real(NER), intent(in) :: Gamma(:), msh(:), wght(:), potM(:, :)

    real(NER), intent(out)  :: elem

    integer     :: ii, jj
    real(NER)   :: tmpsum1, tmpsum2, ppi, ppj

    tmpsum1 = 0.0
    tmpsum2 = 0.0
    do ii = 1, N
      ppi = msh(ii)
      tmpsum1 = tmpsum1 + Gamma(ii)*Gamma(ii)*ppi*ppi*wght(ii)
      do jj = 1, N
        ppj = msh(jj)
        tmpsum2 = tmpsum2 + Gamma(ii)*Gamma(jj)*potM(ii, jj)*ppi*ppi*ppj*ppj*wght(ii)*wght(jj)
      end do
    end do

    elem = tmpsum2/(mN*tmpsum1)

  end subroutine matrx_elmnt_eigen_modeigens


end module nneft_eigens
