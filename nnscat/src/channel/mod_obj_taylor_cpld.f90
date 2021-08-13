! Bingwei Long 11/05/2013
! Class obj_taylor_cpld (Calculate perturbative insertions in the coupled channels
! by numerical Taylor expansion of the full iteration.)

! mod_obj_taylor_sngl is now obsolete, 07/14/2018
module mod_obj_taylor_cpld

  use mod_obj_cpld
  use util_cheb
  implicit none

  type, extends(obj_cpld) :: obj_taylor_cpld

    type(VM_cpld)   :: VM_VTOT

  contains

    procedure   :: init                         => init_taylor_cpld
    procedure   :: release                      => release_taylor_cpld

    procedure   :: GetSinglePertV_default   &
      & => GetSinglePertV_default_taylor_cpld

    procedure   :: GetSinglePertV_input_ab  &
      & => GetSinglePertV_input_ab_taylor_cpld

    procedure   :: GetDoublePertV_input_ab  &
      & => GetDoublePertV_input_ab_taylor_cpld

    procedure   :: GetDoublePertV_default   &
      & => GetDoublePertV_default_taylor_cpld

  end type obj_taylor_cpld


contains


  subroutine init_taylor_cpld(self, N, regtype, Mambda)

    class(obj_taylor_cpld), intent(inout) :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    call init_channel(self, N, regtype, Mambda)
    call self%create_Tk(self%Tk)
    call self%create_VM(self%VM_VTOT)

  end subroutine init_taylor_cpld


  subroutine release_taylor_cpld(self)

    class(obj_taylor_cpld), intent(inout) :: self

    call self%release_Tk(self%Tk)
    call self%release_VM(self%VM_VTOT)
    call release_channel(self)

  end subroutine release_taylor_cpld

  ! n insertions of VM_pert. a and b are the boundary for the dimensionless
  ! parameter x, by which the Taylor expansion is defined (see
  ! cheb_taylor_real). cheb_base+n+1 is the actual number of insertions
  ! calculated by cheb_taylor subroutines. TQn(i) is the value with (i-1)
  ! insertions.
  ! Note: TQ(1) is expected to be zero if the accuracy is infinite, so its
  ! returned value serves as a metric for actual numerical accuracy.
  subroutine GetSinglePertV_input_ab_taylor_cpld(self, VM_lo, n, VM_pert, TQn, a, b, cheb_base_n)

    class(obj_taylor_cpld), intent(inout)           :: self
    type(VM_cpld), intent(inout)                    :: VM_lo, VM_pert
    complex(NEC), dimension(:, :, :), intent(out)   :: TQn
    integer, intent(in)                             :: n, cheb_base_n
    real(NER), intent(in)                           :: a, b

    real(NER), dimension(6, 1:cheb_base_n+n+5) :: g
    integer :: ii

    call cheb_taylor(a, b, 6, cheb_base_n+n+1, g, func_amp)
    do ii = 1, n+1
      TQn(1, 1, ii) = cmplx(g(1, ii), g(2, ii))
      TQn(2, 2, ii) = cmplx(g(3, ii), g(4, ii))
      TQn(1, 2, ii) = cmplx(g(5, ii), g(6, ii))
      TQn(2, 1, ii) = TQn(1, 2, ii)
    end do

  contains

    subroutine func_amp(x, fnval)

      real(NER), intent(in) :: x
      real(NER), dimension(:), intent(out) :: fnval

      call get_allTQn_channel(self)
      call self%b_times_x_plus_a(VM_lo, VM_pert, x, self%VM_VTOT)
      call self%get_Tk_from_VM(self%VM_VTOT, self%Tk)
      fnval(1) = real(self%Tk%ptr(1, 1, self%N+1))
      fnval(2) = aimag(self%Tk%ptr(1, 1, self%N+1))
      fnval(3) = real(self%Tk%ptr(2, 2, self%N+1))
      fnval(4) = aimag(self%Tk%ptr(2, 2, self%N+1))
      fnval(5) = real(self%Tk%ptr(1, 2, self%N+1))
      fnval(6) = aimag(self%Tk%ptr(1, 2, self%N+1))

    end subroutine func_amp

  end subroutine GetSinglePertV_input_ab_taylor_cpld

  subroutine GetSinglePertV_default_taylor_cpld(self, VM_lo, n, VM_pert, TQn)

    class(obj_taylor_cpld), intent(inout)     :: self
    type(VM_cpld), intent(inout)       :: VM_lo, VM_pert
    integer, intent(in)                     :: n
    complex(NEC), dimension(:, :, :), intent(out) :: TQn

    real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER, CHEB_B_PERT = 0.05_NER
    ! real(NEC), parameter    :: CHEB_A_PERT = -0.1_NER, CHEB_B_PERT = 0.1_NER
    integer, parameter      :: CHEB_BASE_N = 2

    call GetSinglePertV_input_ab_taylor_cpld(self, VM_lo, n, VM_pert, TQn, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)

  end subroutine GetSinglePertV_default_taylor_cpld

  ! TQn(1:2, 1:2, i, j) = (i-1) insertions of v1 and (j-1) insertions of v2
  ! Note: TQ(1, 1) is expected to be zero if the accuracy is infinite, so its
  ! returned value serves as a metric for actual numerical accuracy.
  subroutine GetDoublePertV_input_ab_taylor_cpld(self, VM_lo, n1, VM_1, n2, VM_2, TQn,  a1, b1, a2, b2, cheb_base_n1, cheb_base_n2)

    class(obj_taylor_cpld), intent(inout)                 :: self
    type(VM_cpld), intent(inout)                   :: VM_lo, VM_1, VM_2
    complex(NEC), dimension(:, :, :, :), intent(out)    :: TQn
    integer, intent(in)                                 :: n1, n2, cheb_base_n1, cheb_base_n2
    real(NER), intent(in)                               :: a1, b1, a2, b2

    real(NER), dimension(1:6, 1:cheb_base_n1+n1+5, 1:cheb_base_n2+n2+5) :: g
    integer :: ii, jj

    call cheb_taylor(a1, b1, cheb_base_n1+n1+1, a2, b2, cheb_base_n2+n2+1, 6, g, func_amp)
    forall (ii = 1:n1+1, jj = 1:n2+1)
      TQn(1, 1, ii, jj) = cmplx(g(1, ii, jj), g(2, ii, jj))
      TQn(2, 2, ii, jj) = cmplx(g(3, ii, jj), g(4, ii, jj))
      TQn(1, 2, ii, jj) = cmplx(g(5, ii, jj), g(6, ii, jj))
      TQn(2, 1, ii, jj) = TQn(1, 2, ii, jj)
    end forall

  contains

    subroutine func_amp(x, y, rslts)

      real(NER), intent(in)                   :: x, y
      real(NER), dimension(:), intent(out)    :: rslts

      call get_allTQn_channel(self)
      call self%c_times_y_plus_b_times_x_plus_a(VM_lo, VM_1, x, VM_2, y, self%VM_VTOT)
      call self%get_Tk_from_VM(self%VM_VTOT, self%Tk)
      rslts(1) = real(self%Tk%ptr(1, 1, self%N+1))
      rslts(2) = aimag(self%Tk%ptr(1, 1, self%N+1))
      rslts(3) = real(self%Tk%ptr(2, 2, self%N+1))
      rslts(4) = aimag(self%Tk%ptr(2, 2, self%N+1))
      rslts(5) = real(self%Tk%ptr(1, 2, self%N+1))
      rslts(6) = aimag(self%Tk%ptr(1, 2, self%N+1))

    end subroutine func_amp

  end subroutine GetDoublePertV_input_ab_taylor_cpld

  ! Generally, the smaller CHEB_A_PERT and CHEB_B_PERT and the larger
  ! CHEB_BASE_N, the more precisions we get
  subroutine GetDoublePertV_default_taylor_cpld(self, VM_lo, n1, VM_1, n2, VM_2, TQn)

    class(obj_taylor_cpld), intent(inout)           :: self
    type(VM_cpld), intent(inout)                    :: VM_lo, VM_1, VM_2
    complex(NEC), dimension(:, :, :, :), intent(out):: TQn
    integer, intent(in)                             :: n1, n2

    ! real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER, CHEB_B_PERT = 0.05_NER
    real(NEC), parameter    :: CHEB_A_PERT = -0.001_NER, CHEB_B_PERT = -CHEB_A_PERT
    integer, parameter      :: CHEB_BASE_N = 3

    call self%GetDoublePertV_input_ab(VM_lo, n1, VM_1, n2, VM_2, TQn, CHEB_A_PERT, CHEB_B_PERT, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N, CHEB_BASE_N)

  end subroutine GetDoublePertV_default_taylor_cpld


end module mod_obj_taylor_cpld
