! Bingwei Long 10/23/2013
! Class obj_taylor_sngl (Calculate perturbative insertions in the single channels
! by numerical Taylor expansion of the full iteration.)

! mod_obj_taylor_sngl is now obsolete, 07/14/2018
module mod_obj_taylor_sngl

  use mod_obj_sngl
  use util_cheb
  implicit none

  type, extends(obj_sngl) :: obj_taylor_sngl

    type(VM_sngl)   :: VM_VTOT

  contains

    procedure   :: init                     => init_taylor_sngl
    procedure   :: release                  => release_taylor_sngl
    procedure   :: GetSinglePertV_default  => GetSinglePertV_default_taylor_sngl
    procedure   :: GetSinglePertV_input_ab => GetSinglePertV_input_ab_taylor_sngl
    procedure   :: GetDoublePertV_input_ab => GetDoublePertV_input_ab_taylor_sngl
    procedure   :: GetDoublePertV_default   => GetDoublePertV_default_taylor_sngl

  end type obj_taylor_sngl


contains


  subroutine init_taylor_sngl(self, N, regtype, Mambda)

    class(obj_taylor_sngl), intent(inout) :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    call init_channel(self, N, regtype, Mambda)
    call self%create_Tk(self%Tk)
    call self%create_VM(self%VM_VTOT)

  end subroutine init_taylor_sngl

  subroutine release_taylor_sngl(self)

    class(obj_taylor_sngl), intent(inout) :: self

    call self%release_Tk(self%Tk)
    call self%release_VM(self%VM_VTOT)
    call release_channel(self)

  end subroutine release_taylor_sngl

  ! n insertions of VM_pert. a and b are the boundary for the dimensionless
  ! parameter x, by which the Taylor expansion is defined (see
  ! cheb_taylor_real). cheb_base+n+1 is the actual number of insertions
  ! calculated by cheb_taylor subroutines. TQn(i) is the value with (i-1)
  ! insertions.
  subroutine GetSinglePertV_input_ab_taylor_sngl(self, VM_lo, n, VM_pert, TQn, a, b, cheb_base_n)

    class(obj_taylor_sngl), intent(inout)   :: self
    type(VM_sngl), intent(inout)            :: VM_lo, VM_pert
    complex(NEC), dimension(:), intent(out) :: TQn
    integer, intent(in)                     :: n, cheb_base_n
    real(NER), intent(in)                   :: a, b

    complex(NEC), dimension(1:cheb_base_n+n+5) :: g

    call cheb_taylor(a, b, cheb_base_n+n+1, g, func_totamp)
    TQn(1:n+1) = g(1:n+1)

  contains

    function func_totamp(x)

      real(NER), intent(in) :: x
      complex(NEC) :: func_totamp

      call get_allTQn_channel(self)
      call self%b_times_x_plus_a(VM_lo, VM_pert, x, self%VM_VTOT)
      call self%get_Tk_from_VM(self%VM_VTOT, self%Tk)
      func_totamp = self%Tk%ptr(self%N+1)

    end function func_totamp

  end subroutine GetSinglePertV_input_ab_taylor_sngl


  subroutine GetSinglePertV_default_taylor_sngl(self, VM_lo, n, VM_pert, TQn)

    class(obj_taylor_sngl), intent(inout)   :: self
    type(VM_sngl), intent(inout)            :: VM_lo, VM_pert
    integer, intent(in)                     :: n
    complex(NEC), dimension(:), intent(out) :: TQn

    real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER, CHEB_B_PERT = 0.05_NER
    integer, parameter      :: CHEB_BASE_N = 4

    call GetSinglePertV_input_ab_taylor_sngl(self, VM_lo, n, VM_pert, TQn, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)

  end subroutine GetSinglePertV_default_taylor_sngl


  ! TQn(i, j) = (i-1) insertions of v1 and (j-1) insertions of v2
  subroutine GetDoublePertV_input_ab_taylor_sngl(self, VM_lo, n1, VM_1, n2, VM_2, TQn,  a1, b1, cheb_base_n1, a2, b2, cheb_base_n2)

    class(obj_taylor_sngl), intent(inout)       :: self
    type(VM_sngl), intent(inout)                :: VM_lo, VM_1, VM_2
    complex(NEC), dimension(:, :), intent(out)  :: TQn
    integer, intent(in)                         :: n1, n2, cheb_base_n1, cheb_base_n2
    real(NER), intent(in)                       :: a1, b1, a2, b2

    real(NER), dimension(1:2, 1:cheb_base_n1+n1+5, 1:cheb_base_n2+n2+5) :: g
    integer :: ii, jj

    call cheb_taylor(a1, b1, cheb_base_n1+n1+1, a2, b2, cheb_base_n2+n2+1, 2, g, func_amp)
    forall (ii = 1:n1+1, jj = 1:n2+1)
      TQn(ii, jj) = cmplx(g(1, ii, jj), g(2, ii, jj))
    end forall

  contains

    subroutine func_amp(x, y, rslts)

      real(NER), intent(in)                   :: x, y
      real(NER), dimension(:), intent(out)    :: rslts

      call get_allTQn_channel(self)
      call self%c_times_y_plus_b_times_x_plus_a(VM_lo, VM_1, x, VM_2, y,       &
      & self%VM_VTOT)
      call self%get_Tk_from_VM(self%VM_VTOT, self%Tk)
      rslts(1) = real(self%Tk%ptr(self%N+1))
      rslts(2) = aimag(self%Tk%ptr(self%N+1))

    end subroutine func_amp

  end subroutine GetDoublePertV_input_ab_taylor_sngl


  subroutine GetDoublePertV_default_taylor_sngl(self, VM_lo, n1, VM_1, n2, VM_2, TQn)

    class(obj_taylor_sngl), intent(inout)       :: self
    type(VM_sngl), intent(inout)                :: VM_lo, VM_1, VM_2
    complex(NEC), dimension(:, :), intent(out)  :: TQn
    integer, intent(in)                         :: n1, n2

    real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER, CHEB_B_PERT = 0.05_NER
    integer, parameter      :: CHEB_BASE_N = 2

    call self%GetDoublePertV_input_ab(VM_lo, n1, VM_1, n2, VM_2, TQn, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)

  end subroutine GetDoublePertV_default_taylor_sngl


end module mod_obj_taylor_sngl
