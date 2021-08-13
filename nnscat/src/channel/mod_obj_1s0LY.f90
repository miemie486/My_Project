! Bingwei Long 07/10/2013
! Class: 1S0 with BwL & CJY counting, PRC 86, 024001 (2012)


module mod_obj_1s0LY

  use mod_obj_withc0_sngl
  implicit none

  type, extends(obj_withc0_sngl)  :: obj_1s0LY

  contains

    procedure   :: init                 => init_1s0LY
    procedure   :: release              => release_1s0LY
    procedure   :: get_num_paras        => get_num_paras_1s0LY
    procedure   :: update_k             => update_k_1s0LY
    procedure   :: get_C1D0_by_phase    => get_C1D0_by_phase_1s0LY
    procedure   :: get_allTQn           => get_allTQn_1s0LY

  end type obj_1s0LY


contains


  subroutine init_1s0LY(self, N, regtype, Mambda)

    class(obj_1s0LY), intent(inout)  :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    if (self%chnid /= CHNID_1S0_PP) &
      write (standard_error_unit, '(a)')  &
        'init_1s0LY: Warning: chnid should be CHNID_1S0_PP.'
    call init_withc0_sngl(self, N, regtype, Mambda)
    if (self%uptoQn > 1) call self%create_TM(self%Tpq)

  end subroutine init_1s0LY

  subroutine release_1s0LY(self)

    class(obj_1s0LY), intent(inout)  :: self

    if (self%Tpq%created) call self%release_TM(self%Tpq)
    call release_withc0_sngl(self)

  end subroutine release_1s0LY

  subroutine get_num_paras_1s0LY(self, num_cnttrm)

    class(obj_1s0LY), intent(inout)  :: self
    integer, intent(out)                :: num_cnttrm

    select case (self%uptoQn)
      case (0)
        ! C0
        num_cnttrm = 1
      case (1)
        ! C0, C1, D0
        num_cnttrm = 3
      case (2)
        ! C0, C1, D0, C2, D1, E0
        num_cnttrm = 6
      case default
        write (standard_error_unit, '(a, i2)')  &
          'get_num_paras_1s0LY: I don''t know how to handle order ', self%uptoQn
        num_cnttrm = 0
    end select

  end subroutine get_num_paras_1s0LY


  ! inparas(2) = kcm_1
  ! inparas(3) = phase_1
  ! inparas(4) = kcm_2
  ! inparas(5) = phase_2

  subroutine get_C1D0_by_phase_1s0LY(self)

    class(obj_1s0LY), intent(inout)  :: self

    integer         :: tmp_uptoQn, uflag(1:MAX_NUMBER_ORDER)
    real(NER)       :: a(1:2, 1:2), b(1:2)
    complex(NEC)    :: CT, DT

    call self%update_k(self%inparas(2))
    tmp_uptoQn = self%uptoQn
    self%uptoQn = 0
    call self%get_allTQn()
    call self%convert_TQn_to_phase_quiet()
    b(1) = real(delta2t_nlo(self%inparas(3)-self%phsn(1), self%phsn(1), self%k))
    call self%get_CDT(self%Tk, CT, DT)
    a(1, 1) = real(CT)
    a(1, 2) = real(DT)

    call self%update_k(self%inparas(4))
    call self%get_allTQn()
    call self%convert_TQn_to_phase_quiet()
    b(2) = real(delta2t_nlo(self%inparas(5)-self%phsn(1), self%phsn(1), self%k))
    call self%get_CDT(self%Tk, CT, DT)
    a(2, 1) = real(CT)
    a(2, 2) = real(DT)

    ! C1 = Cs(2), D0 = Cs(3)
    self%Cs(3) = (b(2)*a(1, 1) - b(1)*a(2, 1))/(a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))
    self%Cs(2) = (b(1)-a(1, 2)*self%Cs(3))/a(1, 1)

    self%uptoQn = tmp_uptoQn

  end subroutine get_C1D0_by_phase_1s0LY


  subroutine update_k_1s0LY(self, k)

    class(obj_1s0LY), intent(inout)  :: self
    real(NER), intent(in)               :: k

    call update_k_withc0_sngl(self, k)
    self%Tpq%updated = .false.

  end subroutine update_k_1s0LY


  subroutine get_allTQn_1s0LY(self)

    class(obj_1s0LY), intent(inout)  :: self

    integer         :: N
    real(NER)       :: VL2
    complex(NEC)    :: t2, L1, L2, CT, DT, ET, xi0, xi2, xi4

    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        & 'get_allTQn_1s0LY: k is not defined. Nothing will be done.'
      return
    end if
    N = self%N
    if (.not. self%VMVLO%inited) call self%init_VMVLO(self%Cs(1))

    if (self%uptoQn > 1) then
      call self%get_Tpq_from_VM(self%VMVLO, self%Tpq)
      self%Tk%ptr(1:N+1) = self%Tpq%ptr(1:N+1, N+1)
      self%Tk%updated = .true.
    else
      call self%get_Tk_from_VM(self%VMVLO, self%Tk)
    end if
    self%TQn(1) = self%Tk%ptr(N+1)
    if (self%uptoQn == 0) then
      self%TQn_updated = .true.
      return
    end if

    if (self%uptoQn == 1) then
      call self%get_CDT(self%Tk, CT, DT)
    else
      call self%get_CDET(self%Tk, CT, DT, ET)
    end if
    ! C1 = self%Cs(2), D0 = self%Cs(3)
    self%TQn(2) = self%Cs(2)*CT + self%Cs(3)*DT
    if (self%uptoQn == 1) then
      self%TQn_updated = .true.
      return
    end if

    VL2 = self%vs_vl2%VM%ptr(N+1, N+1)
    call self%get_1loop_os_VMTk(self%vs_vl2%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl2%VM, self%Tk, L2)
    ! C2 = self%Cs(4), D1 = self%Cs(5), E0 = self%Cs(6)
    call self%get_CD_xis(self%Tpq, self%Cs(2), self%Cs(3), xi0, xi2, xi4)
    self%TQn(3) = VL2 + 2.0_NER*L1 + L2 + (xi0 + self%Cs(4))*CT + (xi2 + self%Cs(5))*DT + (xi4 + self%Cs(6))*ET
    self%TQn_updated = .true.

  end subroutine get_allTQn_1s0LY

end module mod_obj_1s0LY
