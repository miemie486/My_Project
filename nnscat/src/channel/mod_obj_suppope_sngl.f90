! Bingwei Long 10/03/2018
! Bingwei Long 08/22/2018
! Class obj_SUPpope_sngl (Uncoupled channels where OPE is perturbative.)
! Suppressed scenario

module mod_obj_suppope_sngl

  use mod_obj_sngl
  use eft_phaseconv
  implicit none

  type, extends(obj_sngl) :: obj_SUPpope_sngl

    ! vs_dummy_lo will not be initialized by any Vfunc, but is supposed to be
    ! set zero when created
    type(vstruct_sngl)  :: vs_dummy_lo
    ! , vs_VQ4
    type(VM_sngl)       :: VM_VTOT ! To be used by Taylor methods

  contains

    procedure   :: init                 => init_SUPpope_sngl
    procedure   :: release              => release_SUPpope_sngl
    procedure   :: get_num_paras        => get_num_paras_SUPpope_sngl
    procedure   :: update_k             => update_k_SUPpope_sngl
    procedure   :: get_allTQn           => get_allTQn_SUPpope_sngl
    procedure   :: get_othersTQn        => get_othersTQn_SUPpope_sngl
    procedure   :: get_opeladder_pope   => get_opeladder_SUPpope_sngl

    procedure   :: convert_TQn_to_phase_quiet                             &
      & => convert_TQn_to_phase_quiet_SUPpope_sngl

  end type obj_SUPpope_sngl


contains


  subroutine init_SUPpope_sngl(self, N, regtype, Mambda)

    class(obj_SUPpope_sngl), intent(inout) :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    call init_channel(self, N, regtype, Mambda)
    call self%create_Tk(self%Tk)
    call self%create_VM(self%VM_VTOT)
    call self%create_vstruct(self%vs_vl0)
    call self%create_vstruct(self%vs_dummy_lo)
    call self%init_symmtr_vstruct(self%vs_vl0, Vfunc_OPE_sngl)

  end subroutine init_SUPpope_sngl

  subroutine release_SUPpope_sngl(self)

    class(obj_SUPpope_sngl), intent(inout) :: self

    call self%release_vstruct(self%vs_dummy_lo)
    call self%release_vstruct(self%vs_vl0)
    call self%release_Tk(self%Tk)
    call self%release_VM(self%VM_VTOT)
    call release_channel(self)

  end subroutine release_SUPpope_sngl

  subroutine get_num_paras_SUPpope_sngl(self, num_cnttrm)

    class(obj_SUPpope_sngl), intent(inout) :: self
    integer, intent(out)                :: num_cnttrm

    integer :: L

    if (self%uptoQn > 4) then
      write (standard_error_unit, '(a, i2)')  &
      & 'get_num_paras_SUPpope_sngl: Don''t know how to handle order ', self%uptoQn
      num_cnttrm = 0
      return
    end if
    L = self%lsj%L
    if (L > 2) then
      num_cnttrm = 0
      return
    end if

    select case (self%uptoQn)
      case (0, 1)
        num_cnttrm = 0
      case (2)
        if (L == 1) then
          num_cnttrm = 1
        else
          num_cnttrm = 0
        end if
      case (3)
        if (L == 1) num_cnttrm = 2
      case (4)
        if (L == 1) then
          num_cnttrm = 4
        else
          num_cnttrm = 1
        end if
    end select

  end subroutine get_num_paras_SUPpope_sngl

  subroutine update_k_SUPpope_sngl(self, k)

    class(obj_SUPpope_sngl), intent(inout) :: self
    real(NER), intent(in)               :: k

    call update_k_sngl(self, k)
    call self%update_symmtr_edge_vstruct(self%vs_vl0)
    if (self%vstcVQ2%VM%inited) call self%update_symmtr_edge_vstruct(self%vstcVQ2)
    if (self%vstcVQ3%VM%inited) call self%update_symmtr_edge_vstruct(self%vstcVQ3)
    if (self%vstcVQ4%VM%inited) call self%update_symmtr_edge_vstruct(self%vstcVQ4)
    self%Tk%updated = .false.

  end subroutine update_k_SUPpope_sngl

  subroutine get_allTQn_SUPpope_sngl(self)

    class(obj_SUPpope_sngl), intent(inout) :: self

    complex(NEC)    :: tmpTQn(1:MAX_NUMBER_ORDER)
    complex(NEC)    :: TQn_ladder(1:MAX_NUMBER_ORDER)

    call self%get_opeladder_pope(TQn_ladder)
    call self%get_othersTQn(tmpTQn)
    self%TQn(1:self%uptoQn+1) = TQn_ladder(1:self%uptoQn+1)   &
      & + tmpTQn(1:self%uptoQn+1)
    self%TQn_updated = .true.

  end subroutine get_allTQn_SUPpope_sngl

  ! Calculate pure OPE ladders: Vpi + Vpi G0 Vpi + Vpi G0 Vpi G0 Vpi + ...
  subroutine get_opeladder_SUPpope_sngl(self, restTQn)

    class(obj_SUPpope_sngl), intent(inout) :: self
    complex(NEC), intent(inout)            :: restTQn(:)

    real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER,   &
      &                        CHEB_B_PERT = -CHEB_A_PERT
    integer, parameter      :: CHEB_BASE_N = 4

    call get_allTQn_channel(self)
    self%TQn(1) = 0.0_NER
    if (self%uptoQn == 0) return
    call self%GetSinglePertV(self%vs_dummy_lo%VM, self%uptoQn, self%vs_vl0%VM, &
      & self%VM_VTOT, restTQn, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
    self%TQn(1) = 0.0_NER

  end subroutine get_opeladder_SUPpope_sngl

  ! Calculate all the diagrams except for the pure OPE ladder. Fill
  ! in all elements of TQn(1:uptoQn+1)
  subroutine get_othersTQn_SUPpope_sngl(self, restTQn)

    class(obj_SUPpope_sngl), intent(inout) :: self
    complex(NEC), intent(inout)         :: restTQn(:)


    real(NEC), parameter    :: CHEB_A_PERT = -1.0E-6_NER, CHEB_B_PERT = -CHEB_A_PERT
    integer, parameter      :: CHEB_BASE_N = 5
    real(NEC), parameter    :: CHEB_A_2VQ2 = -1.0E-6_NER,   &
      &                        CHEB_B_2VQ2 = -CHEB_A_2VQ2
    integer, parameter      :: CHEB_BASE_2VQ2 = 5
    real(NER)   :: VQamp
    complex(NEC):: tmpTQn_2d(1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER), &
      & two_ope_one_VQ2, two_VQ2, tmpTQn_1d(1:MAX_NUMBER_ORDER), one_ope_one_VQ3

    restTQn(1:self%uptoQn+1) = 0.0_NER
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)') 'get_othersTQn_SUPpope_sngl: k is not defined.'
      return
    end if
    if (self%uptoQn <= 1) return

    if (self%Lsj%L == 1) then

      call Vfunc_VQ2_SUPpope_sngl(self, self%k, self%k, VQamp)
      restTQn(3) = VQamp
      if (self%uptoQn == 2) return

      if (.not. self%vstcVQ2%VM%inited) then
        if (.not. self%vstcVQ2%VM%created) call self%create_vstruct(self%vstcVQ2)
        call self%init_symmtr_vstruct(self%vstcVQ2, Vfunc_VQ2_SUPpope_sngl)
        call self%update_symmtr_edge_vstruct(self%vstcVQ2)
      end if
      call self%GetDoublePertV(self%vs_dummy_lo%VM, 2, self%vs_vl0%VM, 1,      &
        & self%vstcVQ2%VM, self%VM_VTOT, tmpTQn_2d, CHEB_A_PERT, CHEB_B_PERT,  &
        & CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
      call Vfunc_VQ3_SUPpope_sngl(self, self%k, self%k, VQamp)
      restTQn(4) = VQamp + tmpTQn_2d(2, 2)
      if (self%uptoQn == 3) return

      two_ope_one_VQ2 = tmpTQn_2d(3, 2)
      call self%GetSinglePertV(self%vs_dummy_lo%VM, 2, self%vstcVQ2%VM,        &
        & self%VM_VTOT, tmpTQn_1d, CHEB_A_2VQ2, CHEB_B_2VQ2, CHEB_BASE_2VQ2)
      two_VQ2 = tmpTQn_1d(3)
      if (.not. self%vstcVQ3%VM%inited) then
          if (.not. self%vstcVQ3%VM%created)        &
            & call self%create_vstruct(self%vstcVQ3)
          call self%init_symmtr_vstruct(self%vstcVQ3, Vfunc_VQ3_SUPpope_sngl)
          call self%update_symmtr_edge_vstruct(self%vstcVQ3)
      end if
      call self%GetDoublePertV(self%vs_dummy_lo%VM, 1, self%vs_vl0%VM, 1,      &
        & self%vstcVQ3%VM, self%VM_VTOT, tmpTQn_2d, CHEB_A_PERT, CHEB_B_PERT,  &
        & CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
      one_ope_one_VQ3 = tmpTQn_2d(2, 2)
      call Vfunc_VQ4_SUPpope_sngl(self, self%k, self%k, VQamp)
      restTQn(5) = two_ope_one_VQ2 + two_VQ2 + one_ope_one_VQ3 + VQamp
      if (self%uptoQn == 4) return

    else

      if (self%uptoQn == 2) return

      call Vfunc_VQ3_SUPpope_sngl(self, self%k, self%k, VQamp)
      restTQn(4) = VQamp
      if (self%uptoQn == 3) return

      if (.not. self%vstcVQ3%VM%inited) then
          if (.not. self%vstcVQ3%VM%created)        &
            & call self%create_vstruct(self%vstcVQ3)
          call self%init_symmtr_vstruct(self%vstcVQ3, Vfunc_VQ3_SUPpope_sngl)
          call self%update_symmtr_edge_vstruct(self%vstcVQ3)
      end if
      call self%GetDoublePertV(self%vs_dummy_lo%VM, 1, self%vs_vl0%VM, 1,      &
        & self%vstcVQ3%VM, self%VM_VTOT, tmpTQn_2d, CHEB_A_PERT, CHEB_B_PERT,  &
        & CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
      one_ope_one_VQ3 = tmpTQn_2d(2, 2)
      call Vfunc_VQ4_SUPpope_sngl(self, self%k, self%k, VQamp)
      restTQn(5) = one_ope_one_VQ3 + VQamp
      if (self%uptoQn == 4) return
    end if

    write (standard_error_unit, '(a)')  &
      & 'get_othersTQn_SUPpope_sngl: don''t know to do uptoQn > 4'

  end subroutine get_othersTQn_SUPpope_sngl

  subroutine convert_TQn_to_phase_quiet_SUPpope_sngl(self)

    class(obj_SUPpope_sngl), intent(inout)   :: self

    integer               :: ii

    if (.not. self%TQn_updated) then
      write (standard_error_unit, '(a)')  &
      & 'convert_TQn_to_phase_quiet_SUPpope_sngl: Warning: TQn is not updated.'
      return
    end if

    self%uflag(1) = UNITARITY_SAFE
    self%phsn(1) = 0.0_NER
    do ii = 1, self%uptoQn
      call get_nth_delta_sngl(ii, self%TQn(ii+1), self%phsn(1:ii), self%k,     &
        & self%uflag(ii+1), self%univio(ii+1), self%phsn(ii+1))
    end do
    self%phsn_updated = .true.
    call convert_TQn_to_phase_quiet_channel(self)

  end subroutine convert_TQn_to_phase_quiet_SUPpope_sngl

  subroutine Vfunc_VQ2_SUPpope_sngl(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER)   :: vs

    v = 0.0_NER
    if (osngl%lsj%L == 1) then
      call Vfunc_VS0_sngl(osngl, p1, p2, vs)
      v = osngl%Cs(1)*vs
    end if

  end subroutine Vfunc_VQ2_SUPpope_sngl

  subroutine Vfunc_VQ3_SUPpope_sngl(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER)   :: vl, vs

    call osngl%pvfunc_VL2(osngl, p1, p2, vl)
    select case (osngl%lsj%L)
      case (1)
        call Vfunc_VS0_sngl(osngl, p1, p2, vs)
        v = vl + osngl%Cs(2)*vs
      case default
        v = vl
    end select

  end subroutine Vfunc_VQ3_SUPpope_sngl

  subroutine Vfunc_VQ4_SUPpope_sngl(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER)   :: vl, vs

    call osngl%pVfunc_VL3(osngl, p1, p2, vl)
    select case (osngl%lsj%L)
      case (1)
        call Vfunc_VS0_sngl(osngl, p1, p2, vs)
        v = vl + (osngl%Cs(3) + osngl%Cs(4)*(p1*p1 + p2*p2))*vs
      case (2)
        call Vfunc_VS0_sngl(osngl, p1, p2, vs)
        v = vl + osngl%Cs(1)*vs
      case default
        v = vl
    end select

  end subroutine Vfunc_VQ4_SUPpope_sngl

end module mod_obj_suppope_sngl
