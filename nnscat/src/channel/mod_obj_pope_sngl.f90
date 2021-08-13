! Bingwei Long 10/03/2018
! Bingwei Long 10/23/2013
! Class obj_pope_sngl (Uncoupled channels where OPE is perturbative.)
! Unsuppressed scenario for uncoupled channel

module mod_obj_pope_sngl

  use mod_obj_sngl
  use eft_phaseconv
  implicit none

  type, extends(obj_sngl) :: obj_pope_sngl

    ! vs_dummy_lo is assumed to be set zero when created
    type(vstruct_sngl)  :: vs_dummy_lo
    ! , vs_VQ4
    type(VM_sngl)       :: VM_VTOT ! To be used by Taylor methods

  contains

    procedure   :: init                 => init_pope_sngl
    procedure   :: release              => release_pope_sngl
    ! procedure   :: get_paras_filename   => get_paras_filename_pope_sngl
    procedure   :: get_num_paras        => get_num_paras_pope_sngl
    procedure   :: update_k             => update_k_pope_sngl
    procedure   :: get_allTQn           => get_allTQn_pope_sngl
    procedure   :: get_othersTQn        => get_othersTQn_pope_sngl
    procedure   :: get_opeladder_pope   => get_opeladder_pope_sngl

    procedure   :: convert_TQn_to_phase_quiet       &
    & => convert_TQn_to_phase_quiet_pope_sngl

  end type obj_pope_sngl


contains


  subroutine init_pope_sngl(self, N, regtype, Mambda)

    class(obj_pope_sngl), intent(inout) :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    call init_channel(self, N, regtype, Mambda)
    call self%create_Tk(self%Tk)
    call self%create_VM(self%VM_VTOT)
    call self%create_vstruct(self%vs_vl0)
    call self%create_vstruct(self%vs_dummy_lo)
    call self%init_symmtr_vstruct(self%vs_vl0, self%pVfunc_VL0)

  end subroutine init_pope_sngl

  subroutine release_pope_sngl(self)

    class(obj_pope_sngl), intent(inout) :: self

    call self%release_vstruct(self%vs_dummy_lo)
    call self%release_vstruct(self%vs_vl0)
    call self%release_Tk(self%Tk)
    call self%release_VM(self%VM_VTOT)
    call release_channel(self)

  end subroutine release_pope_sngl

  ! subroutine get_paras_filename_pope_sngl(self, fname)

  !   class(obj_pope_sngl), intent(inout) :: self
  !   character(len = *), intent(out)     :: fname

  !   fname = 'run_'//self%chnstr//'_Lambda.in'

  ! end subroutine get_paras_filename_pope_sngl

  subroutine get_num_paras_pope_sngl(self, num_cnttrm)

    class(obj_pope_sngl), intent(inout) :: self
    integer, intent(out)                :: num_cnttrm

    integer :: L

    if (self%uptoQn > 3) then
      write (standard_error_unit, '(a, i2)')  &
      & 'get_num_paras_pope_sngl: Don''t know how to handle order ', self%uptoQn
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
        if (L == 1) then
          num_cnttrm = 3
        else
          num_cnttrm = 1
        end if
    end select

  end subroutine get_num_paras_pope_sngl

  subroutine update_k_pope_sngl(self, k)

    class(obj_pope_sngl), intent(inout) :: self
    real(NER), intent(in)               :: k

    call update_k_sngl(self, k)
    call self%update_symmtr_edge_vstruct(self%vs_vl0)
    if (self%vstcVQ2%VM%inited) call self%update_symmtr_edge_vstruct(self%vstcVQ2)
    if (self%vstcVQ3%VM%inited) call self%update_symmtr_edge_vstruct(self%vstcVQ3)
    self%Tk%updated = .false.

  end subroutine update_k_pope_sngl

  subroutine get_allTQn_pope_sngl(self)

    class(obj_pope_sngl), intent(inout) :: self

    complex(NEC)    :: tmpTQn(1:MAX_NUMBER_ORDER)
    complex(NEC)    :: TQn_ladder(1:MAX_NUMBER_ORDER)

    call self%get_opeladder_pope(TQn_ladder)
    call self%get_othersTQn(tmpTQn)
    self%TQn(1:self%uptoQn+1) = TQn_ladder(1:self%uptoQn+1)   &
      & + tmpTQn(1:self%uptoQn+1)
    self%TQn_updated = .true.

  end subroutine get_allTQn_pope_sngl

  ! Calculate pure OPE ladders: Vpi + Vpi G0 Vpi + Vpi G0 Vpi G0 Vpi + ...
  subroutine get_opeladder_pope_sngl(self, restTQn)

    class(obj_pope_sngl), intent(inout) :: self
    complex(NEC), intent(inout)         :: restTQn(:)

    real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER,   &
      &                        CHEB_B_PERT = -CHEB_A_PERT
    integer, parameter      :: CHEB_BASE_N = 4

    call get_allTQn_channel(self)
    self%TQn(1) = 0.0_NER
    if (self%uptoQn == 0) return
    call self%GetSinglePertV(self%vs_dummy_lo%VM, self%uptoQn, self%vs_vl0%VM, &
      & self%VM_VTOT, restTQn, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
    self%TQn(1) = 0.0_NER

  end subroutine get_opeladder_pope_sngl

  ! Calculate all the diagrams except for the pure OPE ladder. Fill
  ! in all elements of TQn(1:uptoQn+1)
  subroutine get_othersTQn_pope_sngl(self, restTQn)

    class(obj_pope_sngl), intent(inout) :: self
    complex(NEC), intent(inout)         :: restTQn(:)

    real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER, CHEB_B_PERT = 0.05_NER
    integer, parameter      :: CHEB_BASE_N = 2
    real(NER)   :: VQamp
    complex(NEC):: loopval, tmpTQn(1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER)

    restTQn(1:self%uptoQn+1) = 0.0_NER
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)') 'get_othersTQn_pope_sngl: k is not defined.'
      return
    end if
    if (self%uptoQn <= 1) return
    call Vfunc_VQ2_pope_sngl(self, self%k, self%k, VQamp)
    restTQn(3) = VQamp
    if (self%uptoQn == 2) return

    if (.not. self%vstcVQ2%VM%inited) then
      if (.not. self%vstcVQ2%VM%created) call self%create_vstruct(self%vstcVQ2)
      call self%init_symmtr_vstruct(self%vstcVQ2, Vfunc_VQ2_pope_sngl)
      call self%update_symmtr_edge_vstruct(self%vstcVQ2)
    end if
    call self%Get1LoopVMG0VM_OS(self%vstcVQ2%VM, self%vs_vl0%VM, loopval)
    restTQn(4) = loopval
    call self%Get1LoopVMG0VM_OS(self%vs_vl0%VM, self%vstcVQ2%VM, loopval)
    restTQn(4) = loopval + restTQn(4)
    call Vfunc_VQ3_pope_sngl(self, self%k, self%k, VQamp)
    restTQn(4) = VQamp + restTQn(4)
    if (self%uptoQn == 3) return

    write (standard_error_unit, '(a)') 'get_othersTQn_pope_sngl: don''t know to do uptoQn > 3'

  end subroutine get_othersTQn_pope_sngl

  subroutine convert_TQn_to_phase_quiet_pope_sngl(self)

    class(obj_pope_sngl), intent(inout)   :: self

    integer               :: ii

    if (.not. self%TQn_updated) then
      write (standard_error_unit, '(a)')  &
      & 'convert_TQn_to_phase_quiet_pope_sngl: Warning: TQn is not updated.'
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

  end subroutine convert_TQn_to_phase_quiet_pope_sngl

  subroutine Vfunc_VQ2_pope_sngl(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER)   :: vl, vs

! debug
!   print *, osngl%Cs(1),osngl%Cs(2),osngl%Cs(3)
    call osngl%pVfunc_VL2(osngl, p1, p2, vl)
    if (osngl%lsj%L == 1) then
      call Vfunc_VS0_sngl(osngl, p1, p2, vs)
      v = vl + osngl%Cs(1)*vs
    else
      v = vl
    end if

  end subroutine Vfunc_VQ2_pope_sngl

  subroutine Vfunc_VQ3_pope_sngl(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER)   :: vl, vs

! debug
!   print *, osngl%Cs(1),osngl%Cs(2),osngl%Cs(3)
    call osngl%pVfunc_VL3(osngl, p1, p2, vl)
    select case (osngl%lsj%L)
      case (1)
        call Vfunc_VS0_sngl(osngl, p1, p2, vs)
        v = vl + (osngl%Cs(2) + osngl%Cs(3)*(p1*p1 + p2*p2))*vs
      case (2)
        call Vfunc_VS0_sngl(osngl, p1, p2, vs)
        v = vl + osngl%Cs(1)*vs
      case default
        v = vl
    end select

  end subroutine Vfunc_VQ3_pope_sngl

end module mod_obj_pope_sngl
