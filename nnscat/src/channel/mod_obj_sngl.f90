! Bingwei Long 07/14/2018
! Bingwei Long 03/04/2018
! Bingwei Long 08/31/2015
! Bingwei Long 07/01/2013
! Class obj_sngl (Uncoupled channels)

module mod_obj_sngl

  use mod_obj_channel
  use nneft_lsesv,      only: lsesv_full_sngl, lsesv_half_sngl
  use nneft_loops,      only: get_1loop_TM1G0TM2_slp, get_1loop_Tk1G0Tk2_slp,  &
    & loop_1loop_os_VMTk_sngl, loop_1loop_VkVk_sngl, loop_1loop_VkTk_sngl,     &
    & loop_2loops_TkVMTk_sym_sngl, loop_2loops_Vk1TpqVk2_asym_sngl
  use mod_toypot
  use nneft_eigens
  implicit none

  type TM_sngl
    logical :: created = .false., updated = .false.
    complex(NEC), dimension(:, :), pointer :: ptr => null()
  end type TM_sngl

  type Tk_sngl
    logical :: created = .false., updated = .false.
    complex(NEC), dimension(:), pointer :: ptr => null()
  end type Tk_sngl

  ! Wrapper structure for single-channel VM array and the potential that acts on it.
  type vstruct_sngl

    type(VM_sngl)                           :: VM
    procedure(Vfunc_sngl), pointer, nopass  :: pVfunc => null()

  end type vstruct_sngl


  ! Uncoupled channels where C starts to appear at O(Q^2)
  !
  ! VS =  q^L*p^L*[ (C0+C1+...) + (D0+D1+...)*(q^2 + p^2) + (E0+E1+...)*q^2*p^2]
  !
  type, extends(obj_channel)  :: obj_sngl

    ! TQn and phsn are incremental
    complex(NEC), dimension(1:MAX_NUMBER_ORDER) :: TQn
    real(NER), dimension(1:MAX_NUMBER_ORDER)    :: phsn
    type(vstruct_sngl)  :: vs_vl0, vs_vl2, vs_vl3
    type(vstruct_sngl)  :: vstcVQ1, vstcVQ2, vstcVQ3, vstcVQ4
    type(Tk_sngl)       :: Tk
    type(TM_sngl)       :: Tpq
    complex(NEC)        :: CT_threshold, DT_threshold
    logical             :: CDT_threshold_set = .false.
    procedure(Vfunc_sngl), pointer, nopass  :: pVfunc_VL0 => Vfunc_OPE_sngl
    procedure(Vfunc_sngl), pointer, nopass  :: pVfunc_VL2 => Vfunc_VTPE0_sngl
    procedure(Vfunc_sngl), pointer, nopass  :: pVfunc_VL3 => Vfunc_VTPE1_sngl

  contains

    procedure   :: init                         => init_sngl
    procedure   :: release                      => release_sngl
    procedure   :: get_num_paras                => get_num_paras_sngl

    procedure   :: update_k                     => update_k_sngl
    procedure   :: convert_TQn_to_phase_quiet   => convert_TQn_to_phase_quiet_sngl
    procedure   :: output_phase_to_string       => output_phase_to_string_sngl
    procedure   :: output_total_phase           => output_total_phase_sngl
    procedure   :: output_incremental_phase     => output_incremental_phase_sngl

    ! Methods acting on type(VM_sngl), using attributes of self
    procedure   :: create_VM            => create_VM_objsngl
    procedure   :: release_VM           => release_VM_objsngl
    procedure   :: add_2VMs             => add_2VMs_objsngl
    procedure   :: b_times_x_plus_a     => b_times_x_plus_a_sngl
    procedure   :: c_times_y_plus_b_times_x_plus_a => c_times_y_plus_b_times_x_plus_a_sngl
    procedure   :: init_symmtr_VM       => init_symmtr_VM_sngl
    procedure   :: update_symmtr_edge_VM => update_symmtr_edge_VM_sngl

    ! Methods acting on type(vstruct_sngl), using attributes of self
    procedure   :: create_vstruct       => create_vstruct_sngl
    procedure   :: release_vstruct      => release_vstruct_sngl
    procedure   :: init_symmtr_vstruct  => init_symmtr_vstruct_sngl
    procedure   :: update_symmtr_edge_vstruct => update_symmtr_edge_vstruct_sngl

    ! Methods acting on type(TM_sngl), using attributes of self
    procedure   :: create_TM            => create_TM_sngl
    procedure   :: release_TM           => release_TM_sngl

    ! Methods acting on type(Tk_sngl), using attributes of self
    procedure   :: create_Tk            => create_Tk_sngl
    procedure   :: release_Tk           => release_Tk_sngl

    procedure   :: get_Tk_from_VM       => get_Tk_from_VM_sngl
    procedure   :: get_Tpq_from_VM      => get_Tpq_from_VM_sngl
    ! Frozen
    ! procedure   :: get_kcot             => get_kcot_sngl
    ! The following method is frozen
    ! procedure   :: get_Tk_nonsing     => get_Tk_nonsing_sngl
    procedure   :: get_allTQn           => get_allTQn_sngl
    procedure   :: get_residual_TQn     => get_residual_TQn_sngl

    procedure   :: get_eigens           => get_eigens_sngl
    ! procedure   :: get_eigens_pole      => get_eigens_pole_sngl
    procedure   :: get_TMG0TM           => get_TMG0TM_sngl
    procedure   :: get_TkG0Tk           => get_TkG0Tk_sngl
    procedure   :: get_1loop_os_VMTk    => get_1loop_os_VMTk_sngl
    procedure   :: Get1LoopVMG0VM_OS    => Get1LoopVMG0VM_OS_sngl
    procedure   :: get_2loops_TkVMTk    => get_2loops_TkVMTk_sngl
    procedure   :: get_1loop_powTk      => get_1loop_powTk_sngl
    procedure   :: get_2loops_powTMpow  => get_2loops_powTMpow_sngl
    procedure   :: get_loop_Ik          => get_loop_Ik_sngl
    procedure   :: get_CT               => get_CT_sngl
    procedure   :: get_CDT              => get_CDT_sngl
    procedure   :: get_CDET             => get_CDET_sngl
    ! To be implemented
    ! procedure   :: get_CDEFT             => get_CDEFT_sngl
    procedure   :: get_CD_xis           => get_CD_xis_sngl
    procedure   :: get_1C1Dxi           => get_1C1Dxi_sngl
    procedure   :: get_2Cxi             => get_2Cxi_sngl
    procedure   :: set_CDTthreshold     => set_CDTthreshold_sngl

    ! Methods for Taylor expansion
    procedure   :: GetSinglePertV       => GetSinglePertV_sngl
    procedure   :: GetDoublePertV       => GetDoublePertV_sngl

    ! Frozen
    ! procedure   :: get_B0               => get_B0_sngl
    ! procedure   :: get_B2               => get_B2_sngl

  end type obj_sngl

  ! Abstract interface for single-channel potentials
  abstract interface
    subroutine Vfunc_sngl(os, p1, p2, potval)
      import                          :: obj_sngl, NER
      implicit none
      class(obj_sngl), intent(inout)  :: os
      real(NER), intent(in)           :: p1, p2
      real(NER), intent(out)          :: potval
    end subroutine
    function realfunc_onevar_sngl(q)
      import      :: NER
      implicit none
      real(NER)   :: realfunc_onevar_sngl
      real(NER), intent(in)   :: q
    end function
    function realfunc_twovars_sngl(q, p)
      import      :: NER
      real(NER)   :: realfunc_twovars_sngl
      real(NER), intent(in)   :: q, p
    end function
  end interface

  private :: realfunc_twovars_sngl, realfunc_onevar_sngl

contains

  ! According to self%uptoQn, decide whether to create and initialize vs_vl0,
  ! vs_vl2, and vs_vl3
  subroutine init_sngl(self, N, regtype, Mambda)

    class(obj_sngl), intent(inout)   :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    call init_channel(self, N, regtype, Mambda)
    call self%create_Tk(self%Tk)
    call self%create_vstruct(self%vs_vl0)
    call self%init_symmtr_vstruct(self%vs_vl0, self%pVfunc_VL0)
    if (self%uptoQn < 2) return
    call self%create_vstruct(self%vs_vl2)
    call self%init_symmtr_vstruct(self%vs_vl2, self%pVfunc_VL2)
    if (self%uptoQn < 3) return
    call self%create_vstruct(self%vs_vl3)
    call self%init_symmtr_vstruct(self%vs_vl3, self%pVfunc_VL3)
    if (self%uptoQn > 3) then
      write (standard_error_unit, '(a)')  &
        & 'init_sngl: I do not know how to handle orders that are higher than 3. Order 3 is assumed.'
      self%uptoQn = 3
    end if

  end subroutine init_sngl

  subroutine release_sngl(self)

    class(obj_sngl), intent(inout)    :: self

    call self%release_Tk(self%Tk)
    self%CDT_threshold_set = .false.
    if (associated(self%vstcVQ1%VM%ptr)) call self%release_vstruct(self%vstcVQ1)
    if (associated(self%vstcVQ2%VM%ptr)) call self%release_vstruct(self%vstcVQ2)
    if (associated(self%vstcVQ3%VM%ptr)) call self%release_vstruct(self%vstcVQ3)
    if (associated(self%vs_vl0%VM%ptr)) call self%release_vstruct(self%vs_vl0)
    if (associated(self%vs_vl2%VM%ptr)) call self%release_vstruct(self%vs_vl2)
    if (associated(self%vs_vl3%VM%ptr)) call self%release_vstruct(self%vs_vl3)
    call release_channel(self)

  end subroutine release_sngl

  subroutine get_num_paras_sngl(self, num_cnttrm)

    class(obj_sngl), intent(inout)  :: self
    integer, intent(out)            :: num_cnttrm

    select case (self%uptoQn)
      case (0, 1)
        num_cnttrm = 0
      case (2)
        num_cnttrm = 1
      case (3)
        num_cnttrm = 2
      case default
        write (standard_error_unit, '(a, i2)')  &
          & 'get_num_paras_sngl: I don''t know how to handle order ', self%uptoQn
        num_cnttrm = 0
    end select

  end subroutine get_num_paras_sngl


  subroutine update_k_sngl(self, k)

    class(obj_sngl), intent(inout)    :: self
    real(NER), intent(in)             :: k

    call update_k_channel(self, k)
    self%TQn_updated  = .false.
    self%phsn_updated = .false.
    call self%update_symmtr_edge_vstruct(self%vs_vl0)
    if (self%vs_vl2%VM%inited) then
      call self%update_symmtr_edge_vstruct(self%vs_vl2)
    end if
    if (self%vs_vl3%VM%inited) then
      call self%update_symmtr_edge_vstruct(self%vs_vl3)
    end if
    self%Tk%updated = .false.

  end subroutine update_k_sngl


  ! Convert TQn to phase shifts. uflag are the unitarity flags for each order
  ! Assuming that LO is nonperturbative.
  ! In perturbative OPE calculations, this method needs to be overriden.
  subroutine convert_TQn_to_phase_quiet_sngl(self)

    class(obj_sngl), intent(inout)        :: self
    ! integer, dimension(:), intent(out)  :: uflag

    integer               :: ii
    character(len = 512)  :: msg
    logical               :: is_broken

    if (.not. self%TQn_updated) then
      write (standard_error_unit, '(a)')  &
        'convert_TQn_to_phase_quiet_sngl: Warning: TQn is not updated.'
      return
    end if

    call get_delta_np_sngl(self%TQn(1), self%k, self%univio(1), self%uflag(1), self%phsn(1))
    do ii = 1, self%uptoQn
      call get_nth_delta_sngl(ii, self%TQn(ii+1), self%phsn(1:ii), self%k,     &
        & self%uflag(ii+1), self%univio(ii+1), self%phsn(ii+1))
    end do
    self%phsn_updated = .true.
    call convert_TQn_to_phase_quiet_channel(self)

  end subroutine convert_TQn_to_phase_quiet_sngl

  subroutine output_phase_to_string_sngl(self, outstr)

    class(obj_sngl), intent(inout)  :: self
    character(len = *), intent(out) :: outstr

    call output_phase_to_string_channel(self, outstr)
    call write_incremental_phsn_to_text(self%uptoQn+1, self%phsn, outstr)

  end subroutine output_phase_to_string_sngl

  subroutine output_total_phase_sngl(self, out_phs)

    class(obj_sngl), intent(inout)  :: self
    class(*), intent(inout)         :: out_phs

    call output_total_phase_channel(self, out_phs)
    select type (out_phs)
      type is (real(NER))
        out_phs = sum(self%phsn(1:self%uptoQn+1))
      class default
        write (standard_error_unit, '(a)')  &
          'output_total_phase_sngl: The type of output phase must be real. The returned value is un-defined.'
    end select

  end subroutine output_total_phase_sngl

  subroutine output_incremental_phase_sngl(self, out_phs)

    class(obj_sngl), intent(inout)  :: self
    class(*), intent(inout)         :: out_phs(:)

    call output_incremental_phase_channel(self, out_phs)
    select type (out_phs)
      type is (real(NER))
        out_phs(1:self%uptoQn+1) = self%phsn(1:self%uptoQn+1)
      class default
        write (standard_error_unit, '(a)')  &
          'output_incremental_phase_sngl: The type of output phase must be real. The returned value is un-defined.'
    end select

  end subroutine output_incremental_phase_sngl

  subroutine create_VM_objsngl(self, VM)

    class(obj_sngl), intent(inout)  :: self
    type(VM_sngl), intent(inout)    :: VM

    call self%create_VM_sngl(VM)

  end subroutine create_VM_objsngl


  subroutine release_VM_objsngl(self, VM)

    class(obj_sngl), intent(inout)  :: self
    type(VM_sngl), intent(inout)    :: VM

    call self%release_VM_sngl(VM)

  end subroutine release_VM_objsngl


  subroutine add_2VMs_objsngl(self, VM1, VM2, VMdest)

    class(obj_sngl), intent(inout)  :: self
    type(VM_sngl), intent(inout)    :: VM1, VM2
    type(VM_sngl), intent(inout)    :: VMdest

    call self%add_2VMs_sngl(VM1, VM2, VMdest)

  end subroutine add_2VMs_objsngl


  ! VMc = VMa + x * VMb
  ! If VMa is not inited, then VMc = x * VMb

  subroutine b_times_x_plus_a_sngl(self, VMa, VMb, x, VMrslt)

    class(obj_sngl), intent(inout)  :: self
    type(VM_sngl), intent(inout)    :: VMa, VMb
    real(NEC), intent(in)           :: x
    type(VM_sngl), intent(inout)    :: VMrslt

    integer :: n

    n = self%N
    if (VMa%inited) then
      VMrslt%ptr(1:n+1, 1:n+1) = VMa%ptr(1:n+1, 1:n+1) + x* VMb%ptr(1:n+1, 1:n+1)
    else
      VMrslt%ptr(1:n+1, 1:n+1) = x* VMb%ptr(1:n+1, 1:n+1)
    end if
    VMrslt%inited = .true.

  end subroutine b_times_x_plus_a_sngl


  subroutine c_times_y_plus_b_times_x_plus_a_sngl(self, VMa, VMb, x, VMc, y, VMrslt)

    class(obj_sngl), intent(inout)  :: self
    type(VM_sngl), intent(inout)    :: VMa, VMb, VMc
    real(NEC), intent(in)           :: x, y
    type(VM_sngl), intent(inout)    :: VMrslt

    integer :: n

    n = self%N
    if (VMa%inited) then
      VMrslt%ptr(1:n+1, 1:n+1) = VMa%ptr(1:n+1, 1:n+1) + x* VMb%ptr(1:n+1, 1:n+1) + y* VMc%ptr(1:n+1, 1:n+1)
    else
      VMrslt%ptr(1:n+1, 1:n+1) = x* VMb%ptr(1:n+1, 1:n+1) + y* VMc%ptr(1:n+1, 1:n+1)
    end if
    VMrslt%inited = .true.

  end subroutine c_times_y_plus_b_times_x_plus_a_sngl


  ! Initialize VM(1:N, 1:N)

  subroutine init_symmtr_VM_sngl(self, Vfunc, VM)

    class(obj_sngl), intent(inout)  :: self
    procedure(Vfunc_sngl)           :: Vfunc
    type(VM_sngl), intent(inout)    :: VM

    integer     :: ii, jj
    real(NER)   :: ppi, ppj

    call self%check_warning_inited('init_symmtr_VM_sngl')
    if (.not. VM%created) then
      write (standard_error_unit, '(a)') 'init_symmtr_VM_sngl: VM has not been created.'
      return
    end if
    do ii = 1, self%N
      ppi = self%msh(ii)
      do jj = 1, ii - 1
        ppj = self%msh(jj)
        call Vfunc(self, ppi, ppj, VM%ptr(ii, jj))
        VM%ptr(jj, ii) = VM%ptr(ii, jj)
      end do
      call Vfunc(self, ppi, ppi, VM%ptr(ii, ii))  ! set up the diagonal elements
    end do
    VM%inited = .true.

  end subroutine init_symmtr_VM_sngl


  subroutine update_symmtr_edge_VM_sngl(self, Vfunc, VM)

    class(obj_sngl), intent(inout)  :: self
    procedure(Vfunc_sngl)           :: Vfunc
    type(VM_sngl), intent(inout)    :: VM

    integer     :: jj, N
    real(NER)   :: ppj, k

    if (.not. VM%created) then
      write (standard_error_unit, '(a)')  &
        'update_symmtr_edge_VM_sngl: VM has not been created. Nothing will be done.'
      return
    end if
    if (.not. (self%k_defined)) then
      write(standard_error_unit, '(a)')   &
        'update_symmtr_edge_VM_sngl: k is not defined. Nothing will be done.'
      return
    end if
    k = self%k
    N = self%N
    do jj = 1, N
      ppj = self%msh(jj)
      call Vfunc(self, k, ppj, VM%ptr(N+1, jj))
      VM%ptr(jj, N+1) = VM%ptr(N+1, jj)
    end do
    call Vfunc(self, k, k, VM%ptr(N+1, N+1))

  end subroutine update_symmtr_edge_VM_sngl


  subroutine create_vstruct_sngl(self, vs)

    class(obj_sngl), intent(inout)      :: self
    type(vstruct_sngl), intent(inout)   :: vs

    call self%create_VM(vs%VM)

  end subroutine create_vstruct_sngl


  subroutine release_vstruct_sngl(self, vs)

    class(obj_sngl), intent(inout)      :: self
    type(vstruct_sngl), intent(inout)   :: vs

    call self%release_VM(vs%VM)
    vs%pVfunc => null()

  end subroutine release_vstruct_sngl


  subroutine init_symmtr_vstruct_sngl(self, vs, Vfunc)

    class(obj_sngl), intent(inout)      :: self
    type(vstruct_sngl), intent(inout)   :: vs
    procedure(Vfunc_sngl)               :: Vfunc

    vs%pVfunc => Vfunc
    call self%init_symmtr_VM(vs%pVfunc, vs%VM)

  end subroutine init_symmtr_vstruct_sngl


  subroutine update_symmtr_edge_vstruct_sngl(self, vs)

    class(obj_sngl), intent(inout)      :: self
    type(vstruct_sngl), intent(inout)   :: vs

    call self%update_symmtr_edge_VM(vs%pVfunc, vs%VM)

  end subroutine update_symmtr_edge_vstruct_sngl


  subroutine create_TM_sngl(self, TM)

    class(obj_sngl), intent(inout)  :: self
    type(TM_sngl), intent(inout)    :: TM

    call self%check_warning_inited('create_TM_sngl')
    if (TM%created)   &
      write (standard_error_unit, '(a)')  &
        'create_TM_sngl: Warning! Trying to create a TM that has already been created'
    if (associated(TM%ptr)) deallocate(TM%ptr)
    allocate(TM%ptr(1:self%N+1, 1:self%N+1))
    TM%created = .true.

  end subroutine create_TM_sngl


  subroutine release_TM_sngl(self, TM)

    class(obj_sngl), intent(inout)  :: self
    type(TM_sngl), intent(inout)    :: TM

    call self%check_warning_inited('release_TM_sngl')
    if (.not. associated(TM%ptr)) then
      write (standard_error_unit, '(a)') 'release_TM_sngl: No need to release'
    else
      deallocate(TM%ptr)
    endif
    TM%created = .false.
    TM%updated = .false.

  end subroutine release_TM_sngl


  subroutine get_Tpq_from_VM_sngl(self, VM, TM)

    class(obj_sngl), intent(inout)  :: self
    type(VM_sngl), intent(in)       :: VM
    type(TM_sngl), intent(inout)    :: TM

    integer                         :: flag

    if (.not. self%k_defined) then
      write(standard_error_unit, '(a)')   &
        'get_Tpq_from_VM_sngl: k not defined. I will not do anything.'
      return
    end if
    if (.not. VM%inited) then
      write(standard_error_unit, '(a)')   &
        'get_Tpq_from_VM_sngl: VM not initialized. I will not do anything.'
      return
    end if
    if (.not. TM%created) then
      write(standard_error_unit, '(a)')   &
        'get_Tpq_from_VM_sngl: TM not created. I will not do anything.'
      return
    end if
    call lsesv_full_sngl(self%k, self%N, VM%ptr, self%msh, self%wght, self%lsqr, TM%ptr, flag)
    TM%updated = .true.

  end subroutine get_Tpq_from_VM_sngl


  subroutine create_Tk_sngl(self, Tk)

    class(obj_sngl), intent(inout)  :: self
    type(Tk_sngl), intent(inout)    :: Tk

    call self%check_warning_inited('create_Tk_sngl')
    if (Tk%created)   &
      write (standard_error_unit, '(a)')  &
        'create_Tk_sngl: Warning! Trying to create a Tk that has already been created'
    if (associated(Tk%ptr)) deallocate(Tk%ptr)
    allocate(Tk%ptr(1:self%N+1))
    Tk%created = .true.
    Tk%updated = .false.

  end subroutine create_Tk_sngl


  subroutine release_Tk_sngl(self, Tk)

    class(obj_sngl), intent(inout)  :: self
    type(Tk_sngl), intent(inout)    :: Tk

    call self%check_warning_inited('release_Tk_sngl')
    if (.not. associated(Tk%ptr)) then
      write (standard_error_unit, '(a)') 'release_Tk_sngl: No need to release'
    else
      deallocate(Tk%ptr)
    endif
    Tk%created = .false.
    Tk%updated = .false.

  end subroutine release_Tk_sngl


  subroutine get_Tk_from_VM_sngl(self, VM, Tk)

    class(obj_sngl), intent(inout)  :: self
    type(VM_sngl), intent(in)       :: VM
    type(Tk_sngl), intent(inout)    :: Tk

    integer                         :: flag

    if (.not. self%k_defined) then
      write(standard_error_unit, '(a)')   &
        'get_Tk_from_VM_sngl: k not defined. I will not do anything.'
      return
    end if
    if (.not. VM%inited) then
      write(standard_error_unit, '(a)')   &
        'get_Tk_from_VM_sngl: VM not initialized. I will not do anything.'
      return
    end if
    if (.not. Tk%created) then
      write(standard_error_unit, '(a)')   &
        'get_Tk_from_VM_sngl: Tk not created. I will not do anything.'
      return
    end if
    call lsesv_half_sngl(self%k, self%N, VM%ptr, self%msh, self%wght, self%lsqr, Tk%ptr, flag)
    Tk%updated = .true.

  end subroutine get_Tk_from_VM_sngl


  ! eta*Gamma(p) = 1/f(p)*\int^\Lambda dq*funcM(p, q; B)*q^2*Gamma(q)
  subroutine get_eigens_sngl(self, funcf, funcM, eta, Gamma)

    class(obj_sngl), intent(inout)      :: self
    procedure(realfunc_onevar_sngl)      :: funcf
    procedure(realfunc_twovars_sngl)     :: funcM
    real(NER), intent(out)      :: Gamma(:)
    complex(NEC), intent(out)   :: eta

    integer :: ii, jj, N, flag
    real(NER)   :: pii, pjj, vecf(1:self%N), matM(1:self%N, 1:self%N)

    N = self%N
    do jj = 1, N
      pjj = self%msh(jj)
      vecf(jj) = funcf(pjj)
      do ii = 1, N
        pii = self%msh(ii)
        matM(ii, jj) = funcM(pii, pjj)*pjj*pjj*self%half_regfac(ii)*self%half_regfac(jj)
      end do
    end do

    call eigen_invf_mat_modeigens(N, vecf, matM, self%msh, self%wght, eta, Gamma, flag)
    if (flag == eigens_FAILED) then
      write(standard_error_unit, '(a)')   &
        'get_eigens_sngl: eigen_invf_mat_modeigens failed.'
      eta = -1.0_NER
    end if

  end subroutine get_eigens_sngl


  ! Frozen
  ! ! eta*Gamma(p) = \int^\Lambda dq*M(p, q; B)/[s(q) - i0]*q^2*Gamma(q)
  ! ! s(q) has a zero at q*, where ds/dq = a0^(-1)

  ! subroutine get_eigens_pole_sngl(self, fs, lstar, a0, fM, eta, Gamma)

  !   class(obj_sngl), intent(inout)  :: self
  !   real(NER), intent(in)           :: lstar, a0
  !   procedure(realfunc_onevar_sngl)      :: fs
  !   procedure(realfunc_twovars_sngl)     :: fM
  !   complex(NEC), intent(out)      :: Gamma(:)
  !   complex(NEC), intent(out)   :: eta

  !   integer     :: ii, jj, N, flag
  !   real(NER)   :: pii, pjj, vecs(1:self%N), matM(1:self%N+1, 1:self%N+1)

  !   N = self%N
  !   do jj = 1, N
  !     pjj = self%msh(jj)
  !     vecs(jj) = fs(pjj)
  !     do ii = 1, N
  !       pii = self%msh(ii)
  !       matM(ii, jj) = fM(pii, pjj)*pjj*pjj*self%half_regfac(ii)*self%half_regfac(jj)
  !     end do
  !     matM(N+1, jj) = fM(lstar, pjj)*pjj*pjj*self%half_regulator(lstar)*self%half_regfac(jj)
  !     matM(jj, N+1) = fM(pjj, lstar)*lstar*lstar*self%half_regulator(lstar)*self%half_regfac(jj)
  !   end do
  !   matM(N+1, N+1) = fM(lstar, lstar)*lstar*lstar*self%half_regulator(lstar)*self%half_regulator(lstar)

  !   call eigen_invsM_pole_mat_modeigens(N, vecs, lstar, a0, matM, self%msh, self%wght, eta, Gamma, flag)
  !   if (flag == eigens_FAILED) then
  !     write(standard_error_unit, '(a)')   &
  !       'get_eigens_sngl: eigen_invf_mat_modeigens failed.'
  !     eta = -1.0_NER
  !   end if

  ! end subroutine get_eigens_pole_sngl


  ! Frozen, Jan 30/2018

  ! ! T(q) = V(q) + \int^Lambda_0 dl M(q, l) T(l), where q > 0 and M(q, l) does
  ! ! not have singularity for l > 0
  ! subroutine get_Tk_nonsing_sngl(self, fV, fM, arrV, arrM, arrT)

  !     class(obj_sngl), intent(inout)  :: self
  !     procedure(realfunc_onevar_sngl)      :: fV
  !     procedure(realfunc_twovars_sngl)     :: fM
  !     real(NER), intent(out)          :: arrV(:), arrM(:, :)
  !     real(NER), intent(out)          :: arrT(:)

  !     integer     :: ii, jj, N, flag
  !     real(NER)   :: pjj, pii

  !     N = self%N
  !     do jj = 1, N
  !         pjj = self%msh(jj)
  !         arrV(jj) = fV(pjj)
  !         do ii = 1, N
  !             pii = self%msh(ii)
  !             ! arrM(ii, jj) = fM(pii, pjj)*self%half_regfac(jj)*self%half_regfac(jj)
  !             arrM(ii, jj) = fM(pii, pjj)
  !         end do
  !     end do
  !     call lsesv_nonsing_half_sngl(N, arrV, arrM, self%wght, arrT, flag)
  !     if (flag == lsesv_FAILED) then
  !         write (standard_error_unit, '(a)') 'get_Tk_nonsing_sngl: lsesv_FAILED'
  !     end if

  ! end subroutine get_Tk_nonsing_sngl

  ! Frozen
  ! ! Calculate kcot, assuming that Cs have been setup properly.
  ! subroutine get_kcot_sngl(self, k, kcot)

  !     class(obj_sngl), intent(inout)          :: self
  !     real(NER), intent(in)                   :: k
  !     real(NER), intent(out)                  :: kcot

  !     integer, dimension(1:MAX_NUMBER_ORDER)  :: uflag

  !     call self%update_k(k)
  !     call self%get_allTQn()
  !     call self%convert_TQn_to_phase_quiet(uflag)
  !     call convert_phase_to_kcot_sngl(k, sum(self%phsn(1:self%uptoQn+1)), kcot)

  ! end subroutine get_kcot_sngl

  subroutine get_allTQn_sngl(self)

    class(obj_sngl), intent(inout)   :: self

    real(NER)       :: VL2, VL3
    complex(NEC)    :: t2, t3, L1, L2, CT

    call get_allTQn_channel(self)
    call self%get_Tk_from_VM(self%vs_vl0%VM, self%Tk)
    self%TQn(1) = self%Tk%ptr(self%N+1)
    if (self%uptoQn == 0) then
      self%TQn_updated = .true.
      return
    end if

    self%TQn(2) = 0.0_NER
    if (self%uptoQn == 1)  then
      self%TQn_updated = .true.
      return
    end if

    call self%get_1loop_os_VMTk(self%vs_vl2%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl2%VM, self%Tk, L2)
    call self%get_CT(self%Tk, CT)
    VL2 = self%vs_vl2%VM%ptr(self%N+1, self%N+1)
    ! C2 = self%Cs(1)
    self%TQn(3) = VL2 + 2.0_NER*L1 + L2 + self%Cs(1)*CT
    if (self%uptoQn == 2)  then
      self%TQn_updated = .true.
      return
    end if

    call self%get_1loop_os_VMTk(self%vs_vl3%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl3%VM, self%Tk, L2)
    VL3 = self%vs_vl3%VM%ptr(self%N+1, self%N+1)
    ! C3 = self%Cs(2)
    self%TQn(4) = VL3 + 2.0_NER*L1 + L2 + self%Cs(2)*CT
    self%TQn_updated = .true.

  end subroutine get_allTQn_sngl

  ! TQ(1) = fully iterated LO
  ! uptoQn = 0 : TQ(2) = self%Cs(2) * CT
  ! uptoQn = 1 : TQ(2) = C1 * CT + self%Cs(2) * DT, where C1 is chosen such that
  !   TQ(2) = 0 when k = 0
  subroutine get_residual_TQn_sngl(self)

    class(obj_sngl), intent(inout)   :: self

    integer :: saved_uptoQn
    real(NER)       :: C1
    complex(NEC)    :: CT, DT

    saved_uptoQn = self%uptoQn
    self%uptoQn = 0
    call self%get_allTQn()
    self%TQn(1) = self%Tk%ptr(self%N+1)

    if (saved_uptoQn == 0) then
      call self%get_CT(self%Tk, CT)
      self%TQn(2) = self%Cs(2)*CT
    else
      if (.not. self%CDT_threshold_set) then
        call self%set_CDTthreshold()
        C1 = - self%Cs(2) * real(self%DT_threshold/self%CT_threshold)
      end if
      call self%get_CDT(self%Tk, CT, DT)
      self%TQn(2) = C1*CT + self%Cs(2)*DT
    end if

    self%uptoQn = saved_uptoQn

  end subroutine get_residual_TQn_sngl


!-------------------------------------!
!                                     !
! Loop integrals for single channels  !
!                                     !
!-------------------------------------!

  subroutine get_1loop_os_VMTk_sngl(self, VM, Tk, loopval)

    class(obj_sngl), intent(inout)    :: self
    type(VM_sngl), intent(in)           :: VM
    type(Tk_sngl), intent(in)           :: Tk
    complex(NEC), intent(out)           :: loopval

    call self%check_warning_inited('get_1loop_os_VMTk_sngl')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_1loop_os_VMTk_sngl: k not defined. Nothing will be done.'
      return
    end if
    if (.not. VM%inited) &
      write (standard_error_unit, '(a)')  &
        'get_1loop_os_VMTk_sngl: Warning! VM not initialized.  Nothing will be done.'
    if (.not. Tk%updated) &
      write (standard_error_unit, '(a)')  &
        'get_1loop_os_VMTk_sngl: Warning! Tk not defined. Nothing will be done.'
    loopval = loop_1loop_os_VMTk_sngl(self%k, self%N, self%msh, self%wght, self%lsqr, VM%ptr, Tk%ptr)

  end subroutine get_1loop_os_VMTk_sngl

  subroutine Get1LoopVMG0VM_OS_sngl(self, VMa, VMb, loopval)

    class(obj_sngl), intent(inout) :: self
    type(VM_sngl), intent(in)   :: VMa, VMb
    complex(NEC), intent(out)   :: loopval

    call self%check_warning_inited('Get1LoopVMVM_OS_sngl')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
      & 'Get1LoopVMVM_OS_sngl: k not defined. Nothing will be done.'
      return
    end if
    if ((.not. VMa%inited) .or. (.not. VMb%inited)) &
    & write (standard_error_unit, '(a)')            &
    & 'Get1LoopVMVM_OS_sngl: Warning! VM not initialized.  Nothing will be done.'
    loopval = loop_1loop_VkVk_sngl(self%k, self%N, self%msh,  &
    & self%wght, self%lsqr, VMa%ptr(self%N+1, 1:self%N+1), VMb%ptr(1:self%N+1, self%N+1))

  end subroutine Get1LoopVMG0VM_OS_sngl

  subroutine get_TMG0TM_sngl(self, TMa, TMb, rsltTM)

    class(obj_sngl), intent(inout) :: self
    complex(NEC), intent(in)   :: TMa(:, :), TMb(:, :)
    complex(NEC), intent(out)  :: rsltTM(:, :)

    call self%check_warning_inited('get_TMG0TM_sngl')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
      & 'get_TMG0TM_sngl: k not defined. Nothing will be done.'
      return
    end if
    call get_1loop_TM1G0TM2_slp(self%k, self%N, self%msh,           &
      & self%wght, self%lsqr, TMa, TMb, rsltTM)

  end subroutine get_TMG0TM_sngl

  subroutine get_TkG0Tk_sngl(self, Tka, Tkb, val)

    class(obj_sngl), intent(inout)  :: self
    complex(NEC), intent(in)        :: Tka(:), Tkb(:)
    complex(NEC), intent(out)       :: val

    call self%check_warning_inited('get_TkG0Tk_sngl')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
      & 'get_TkG0Tk_sngl: k not defined. Nothing will be done.'
      return
    end if
    call get_1loop_Tk1G0Tk2_slp(self%k, self%N, self%msh,         &
      & self%wght, self%lsqr, Tka, Tkb, val)

  end subroutine get_TkG0Tk_sngl

  subroutine get_2loops_TkVMTk_sngl(self, VM, Tk, loopval)

    class(obj_sngl), intent(inout) :: self
    type(VM_sngl), intent(in)   :: VM
    type(Tk_sngl), intent(in)   :: Tk
    complex(NEC), intent(out)   :: loopval

    call self%check_warning_inited('get_2loops_TkVMTk_sngl')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_2loops_TkVMTk_sngl: k not defined. Nothing will be done.'
      return
    end if
    if (.not. VM%inited) &
      write (standard_error_unit, '(a)')  &
        'get_2loops_TkVMTk_sngl: Warning! VM not initialized. Nothing will be done.'
    if (.not. Tk%updated) &
      write (standard_error_unit, '(a)')  &
        'get_2loops_TkVMTk_sngl: Warning! Tk not defined. Nothing will be done.'
    loopval = loop_2loops_TkVMTk_sym_sngl(self%k, self%N, self%msh, self%wght, self%lsqr, VM%ptr, Tk%ptr)

  end subroutine get_2loops_TkVMTk_sngl


  subroutine get_1loop_powTk_sngl(self, pow_arrayID, Tk, loopval)

    class(obj_sngl), intent(inout) :: self
    integer, intent(in)            :: pow_arrayID
    type(Tk_sngl), intent(in)      :: Tk
    complex(NEC), intent(out)      :: loopval

    real(NER), dimension(:), pointer    :: ptrVk

    call self%check_warning_inited('get_1loop_powTk_sngl')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_1loop_powTk_sngl: k not defined. Nothing will be done.'
      return
    end if
    if (.not. Tk%updated) &
      write (standard_error_unit, '(a)')  &
        'get_1loop_powTk_sngl: Warning! Tk not defined. Nothing will be done.'
    call self%set_ptr_to_powfac(pow_arrayID, ptrVk)
    loopval = loop_1loop_VkTk_sngl(self%k, self%N, self%msh,      &
      & self%wght, self%lsqr, ptrVk, Tk%ptr)

  end subroutine get_1loop_powTk_sngl

  subroutine get_2loops_powTMpow_sngl(self, n_pow_arrayID, m_pow_arrayID, TM, loopval)

    class(obj_sngl), intent(inout) :: self
    integer, intent(in)            :: n_pow_arrayID, m_pow_arrayID
    type(TM_sngl), intent(in)      :: TM
    complex(NEC), intent(out)      :: loopval

    real(NER), dimension(:), pointer    :: n_ptrVk, m_ptrVk

    call self%check_warning_inited('get_2loops_powTMpow')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_2loops_powTkpow_sngl: k not defined. Nothing will be done.'
      return
    end if
    if (.not. TM%updated) &
      write (standard_error_unit, '(a)')  &
        'get_2loops_powTMpow_sngl: Warning! TM not defined. Nothing will be done.'
    call self%set_ptr_to_powfac(n_pow_arrayID, n_ptrVk)
    call self%set_ptr_to_powfac(m_pow_arrayID, m_ptrVk)
    loopval = loop_2loops_Vk1TpqVk2_asym_sngl(self%k, self%N,     &
      & self%msh, self%wght, self%lsqr, n_ptrVk, m_ptrVk, TM%ptr)

  end subroutine get_2loops_powTMpow_sngl


  subroutine get_loop_Ik_sngl(self, n_pow_arrayID, m_pow_arrayID, TM, loopval)

    class(obj_sngl), intent(inout) :: self
    integer, intent(in)            :: n_pow_arrayID, m_pow_arrayID
    type(TM_sngl), intent(in)      :: TM
    complex(NEC), intent(out)      :: loopval

    complex(NEC)    :: loopTM

    call self%get_2loops_powTMpow(n_pow_arrayID, m_pow_arrayID, TM, loopTM)
    loopval = self%loop_powpow(n_pow_arrayID, m_pow_arrayID) + loopTM

  end subroutine get_loop_Ik_sngl


  ! L : angular momentum
  ! one insertion of C p^L q^L into the LO wave function = C * CT

  subroutine get_CT_sngl(self, Tk, CT)

    class(obj_sngl), intent(inout) :: self
    type(Tk_sngl), intent(in)      :: Tk
    complex(NEC), intent(out)      :: CT

    complex(NEC)                        :: F
    real(NER)                           :: regfac, regsqr
    real(NER), dimension(:), pointer    :: ptr_powfac

    if (.not. Tk%updated) then
      write (standard_error_unit, '(a)')  &
        'get_CT_sngl: Tk has not been updated. Nothing will be done.'
      CT = 0.0_NER
      return
    end if
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_CT_sngl: k is not defined. Nothing will be done.'
      return
    end if
    call self%set_ptr_to_powfac(1, ptr_powfac)
    regfac = ptr_powfac(self%N+1)
    regsqr = regfac*regfac
    call self%get_1loop_powTk(1, Tk, F)
    CT = regsqr + 2.0_NER*F*regfac + F*F

  end subroutine get_CT_sngl


  ! L : angular momentum
  ! one insertion of C p^L q^L into the LO wave function = C * CT
  ! one insertion of D p^L q^L(p^2 + q^2) into the LO wave function = D * DT

  subroutine get_CDT_sngl(self, Tk, CT, DT)

    class(obj_sngl), intent(inout) :: self
    type(Tk_sngl), intent(in)      :: Tk
    complex(NEC), intent(out)      :: CT, DT

    complex(NEC)                        :: F, F2
    real(NER)                           :: ksqr, regfac, regsqr
    real(NER), dimension(:), pointer    :: ptr_powfac

    if (.not. Tk%updated) then
      write (standard_error_unit, '(a)')  &
        'get_CDT_sngl: Tk has not been updated.'
      CT = 0.0_NER
      DT = 0.0_NER
      return
    end if
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_CDT_sngl: k is not defined. Nothing will be done.'
      return
    end if
    call self%set_ptr_to_powfac(1, ptr_powfac)
    regfac = ptr_powfac(self%N+1)
    regsqr = regfac*regfac
    ksqr   = self%k*self%k
    call self%get_1loop_powTk(1, Tk, F)
    call self%get_1loop_powTk(2, Tk, F2)
    CT = regsqr + 2.0_NER*F*regfac + F*F
    DT = 2.0_NER*(ksqr*regsqr + (ksqr*F + F2)*regfac + F2*F)

  end subroutine get_CDT_sngl


  ! L : angular momentum
  ! one insertion of C p^L q^L into the LO wave function = C * CT
  ! one insertion of D p^L q^L(p^2 + q^2) into the LO wave function = D * DT
  ! one insertion of E p^L q^L(p^2q^2) into the LO wave function = E * ET

  subroutine get_CDET_sngl(self, Tk, CT, DT, ET)

    class(obj_sngl), intent(inout)  :: self
    type(Tk_sngl), intent(in)       :: Tk
    complex(NEC), intent(out)       :: CT, DT, ET

    complex(NEC)                        :: F, F2
    real(NER)                           :: ksqr, regfac, regsqr
    real(NER), dimension(:), pointer    :: ptr_powfac

    if (.not. Tk%updated) then
      write (standard_error_unit, '(a)')  &
        'get_CDET_sngl: Tk has not been updated.'
      CT = 0.0_NER
      DT = 0.0_NER
      ET = 0.0_NER
      return
    end if
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_CDET_sngl: k is not defined. Nothing will be done.'
      return
    end if
    call self%set_ptr_to_powfac(1, ptr_powfac)
    regfac = ptr_powfac(self%N+1)
    regsqr = regfac*regfac
    ksqr   = self%k*self%k
    call self%get_1loop_powTk(1, Tk, F)
    call self%get_1loop_powTk(2, Tk, F2)
    CT = regsqr + 2.0_NER*F*regfac + F*F
    DT = 2.0_NER*(ksqr*regsqr + (ksqr*F + F2)*regfac + F2*F)
    ET = ksqr*ksqr*regsqr + 2.0_NER*ksqr*F2*regfac + F2*F2

  end subroutine get_CDET_sngl


  ! Two insertions of (C + D(p^2+q^2))*(pq)^L = xi0*CT + xi2*DT + xi4*ET

  subroutine get_CD_xis_sngl(self, TM, C, D, xi0, xi2, xi4)

    class(obj_sngl), intent(inout)  :: self
    type(TM_sngl), intent(in)       :: TM
    real(NER), intent(in)           :: C, D
    complex(NEC), intent(out)       :: xi0, xi2, xi4

    complex(NEC)    :: I0, I2, I4, J00, J02, J22

    I0 = self%loop_powpow(1, 1)
    I2 = self%loop_powpow(1, 2)
    I4 = self%loop_powpow(2, 2)
    call self%get_2loops_powTMpow(1, 1, TM, J00)
    call self%get_2loops_powTMpow(1, 2, TM, J02)
    call self%get_2loops_powTMpow(2, 2, TM, J22)
    xi0 = D*D*(I4 + J22) + 2.0_NER*C*D*(I2 + J02) + C*C*(I0 + J00)
    xi2 = D*D*(I2 + J02) + C*D*(I0 + J00)
    xi4 = D*D*(I0 + J00)

  end subroutine get_CD_xis_sngl

  ! One insertion of (pq)^L and one insertion of (p^2 + q^2)*(pq)^L
  ! = xi0*CT + xi2*DT
  subroutine get_1C1Dxi_sngl(self, TM, xi0, xi2)

    class(obj_sngl), intent(inout)  :: self
    type(TM_sngl), intent(in)       :: TM
    complex(NEC), intent(out)       :: xi0, xi2

    complex(NEC)    :: I2, J02, J00, I0
    ! complex(NEC)    :: I0, I2, I4, J00, J02, J22

    I0 = self%loop_powpow(1, 1)
    I2 = self%loop_powpow(1, 2)
    ! I4 = self%loop_powpow(2, 2)
    call self%get_2loops_powTMpow(1, 1, TM, J00)
    call self%get_2loops_powTMpow(1, 2, TM, J02)
    ! call self%get_2loops_powTMpow(2, 2, TM, J22)
    ! xi0 = D*D*(I4 + J22) + 2.0_NER*C*D*(I2 + J02) + C*C*(I0 + J00)
    xi0 = 2.0_NER*(I2 + J02)
    ! xi2 = D*D*(I2 + J02) + C*D*(I0 + J00)
    xi2 = I0 + J00
    ! xi4 = D*D*(I0 + J00)

  end subroutine get_1C1Dxi_sngl

  ! Two insertions of (pq)^L = xi0*CT
  subroutine get_2Cxi_sngl(self, TM, xi0)

    class(obj_sngl), intent(inout)  :: self
    type(TM_sngl), intent(in)       :: TM
    complex(NEC), intent(out)       :: xi0

    complex(NEC)    :: I0, J00

    I0 = self%loop_powpow(1, 1)
    call self%get_2loops_powTMpow(1, 1, TM, J00)
    xi0 = (I0 + J00)

  end subroutine get_2Cxi_sngl

  subroutine set_CDTthreshold_sngl(self)

    class(obj_sngl), intent(inout)  :: self

    integer         :: save_uptoQn
    real(NER)       :: zero_k, save_k
    real(NER), parameter    :: local_EPS = 1.0E-6_NER

    save_uptoQn = self%uptoQn
    self%uptoQn = 0
    save_k = self%k
    zero_k = local_EPS*self%Mambda
    call self%update_k(zero_k)
    call self%get_allTQn()
    call self%get_CDT(self%Tk, self%CT_threshold, self%DT_threshold)
    self%CDT_threshold_set = .true.
    self%uptoQn = save_uptoQn
    self%k = save_k

  end subroutine set_CDTthreshold_sngl

  ! Frozen
  ! ! B0 = k/f0(k) * d f0(k)/dk where f0(k) = pi/2*CT(k)/(1 + 2ik t) and t = - pi/2*T
  ! subroutine get_B0_sngl(self, k, B0, err_B0)

  !     class(obj_sngl), intent(inout)  :: self
  !     real(NER), intent(in)           :: k
  !     real(NER), intent(out)          :: B0, err_B0

  !     real, parameter :: local_epsilon = 1.0E-2_NER
  !     real(NER)       :: deri_of_f0, value_of_f0, h, err
  !     complex(NEC)    :: t

  !     value_of_f0 = f0(k)
  !     h = self%Mambda * 1.0E-3_NER
  !     deri_of_f0 = dfridr(f0, k, h, err)
  !     B0 = k*deri_of_f0/value_of_f0
  !     err_B0 = k*err/value_of_f0

  ! contains

  !     function f0(mu)

  !         real(NER), intent(in)   :: mu
  !         real(NER)               :: f0

  !         integer         :: save_uptoQn
  !         complex(NEC)    :: t, CT, tmp_f0

  !         save_uptoQn = self%uptoQn
  !         self%uptoQn = 0
  !         call self%update_k(mu)
  !         call self%get_allTQn()
  !         t = -PI_NE*0.5*self%TQn(1)
  !         call self%get_CT(self%Tk, CT)
  !         ! tmp_f0 = PI_NE*0.5*CT/(t*t)
  !         tmp_f0 = PI_NE*0.5*mu*CT/(1.0_NER + 2.0_NER*IMGUNT_NE*mu*t)
  !         if (aimag(tmp_f0)/abs(tmp_f0) > local_epsilon) then
  !             write (standard_error_unit, '(a, 2x, f9.2, 2x, f9.3)') 'get_B0_woc0_sngl: Warning! f0 is not real.', self%Mambda, mu
  !         end if
  !         f0 = real(tmp_f0)
  !         self%uptoQn = save_uptoQn

  !     end function f0

  ! end subroutine get_B0_sngl


  ! ! B2 = k/f2(k) * d f2(k)/dk where
  ! ! f2(k) = pi/2*[DT(k) - DT(0.0)/CT(0.0)*CT(k)]/(1 + 2ik t) and t = - pi/2*T
  ! subroutine get_B2_sngl(self, k, B2, err_B2)

  !     class(obj_sngl), intent(inout)  :: self
  !     real(NER), intent(in)           :: k
  !     real(NER), intent(out)          :: B2, err_B2

  !     real, parameter :: local_epsilon = 1.0E-2_NER
  !     integer         :: save_uptoQn
  !     real(NER)       :: deri_of_f2, value_of_f2, h, err, zero_k, save_k
  !     complex(NEC)    :: t

  !     if (.not. self%CDT_threshold_set) then
  !         call self%set_CDTthreshold()
  !     end if

  !     value_of_f2 = f2(k)
  !     h = self%Mambda * 1.0E-3_NER
  !     deri_of_f2 = dfridr(f2, k, h, err)
  !     B2 = k*deri_of_f2/value_of_f2
  !     err_B2 = k*err/value_of_f2

  ! contains

  !     function f2(mu)

  !         real(NER), intent(in)   :: mu
  !         real(NER)               :: f2
  !         integer         :: save_uptoQn
  !         complex(NEC)    :: t, CT, DT, tmp_f2

  !         save_uptoQn = self%uptoQn
  !         self%uptoQn = 0
  !         call self%update_k(mu)
  !         call self%get_allTQn()
  !         t = -PI_NE*0.5*self%TQn(1)
  !         call self%get_CDT(self%Tk, CT, DT)
  !         ! tmp_f2 = PI_NE*0.5*(DT - self%DT_threshold*CT/self%CT_threshold)/(t*t)
  !         tmp_f2 = PI_NE*0.5*mu*(DT - self%DT_threshold*CT/self%CT_threshold)/(1.0_NER + 2.0_NER*IMGUNT_NE*mu*t)
  !         if (aimag(tmp_f2)/abs(tmp_f2) > local_epsilon) then
  !             write (standard_error_unit, '(a, 2x, f9.2, 2x, f9.3)') 'get_B2_woc0_sngl: Warning! f2 is not real.', self%Mambda, mu
  !         end if
  !         f2 = real(tmp_f2)
  !         self%uptoQn = save_uptoQn

  !     end function f2

  ! end subroutine  get_B2_sngl


  ! n insertions of VM_pert. a and b are the boundary for the dimensionless
  ! parameter x, by which the Taylor expansion is defined (see
  ! cheb_taylor_real). cheb_base+n+1 is the actual number of insertions
  ! calculated by cheb_taylor subroutines. TQn(i) is the value with (i-1)
  ! insertions.
  ! VM_tmp : working memory to hold VM that will be iterated to all orders.
  subroutine GetSinglePertV_sngl(self, VM_lo, n, VM_pert, VM_tmp, TQn, a, b,   &
    & cheb_base_n)

    class(obj_sngl), intent(inout)          :: self
    type(VM_sngl), intent(inout)            :: VM_lo, VM_pert, VM_tmp
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
      call self%b_times_x_plus_a(VM_lo, VM_pert, x, VM_tmp)
      call self%get_Tk_from_VM(VM_tmp, self%Tk)
      func_totamp = self%Tk%ptr(self%N+1)

    end function func_totamp

  end subroutine GetSinglePertV_sngl

  ! TQn(i, j) = (i-1) insertions of v1 and (j-1) insertions of v2
  subroutine GetDoublePertV_sngl(self, VM_lo, n1, VM_1, n2, VM_2, VM_tmp, TQn, &
    & a1, b1, cheb_base_n1, a2, b2, cheb_base_n2)

    class(obj_sngl), intent(inout)              :: self
    type(VM_sngl), intent(inout)                :: VM_lo, VM_1, VM_2, VM_tmp
    complex(NEC), dimension(:, :), intent(out)  :: TQn
    integer, intent(in)                         :: n1, n2, cheb_base_n1,       &
    & cheb_base_n2
    real(NER), intent(in)                       :: a1, b1, a2, b2

    real(NER), dimension(1:2, 1:cheb_base_n1+n1+5, 1:cheb_base_n2+n2+5) :: g
    integer :: ii, jj

    call cheb_taylor(a1, b1, cheb_base_n1+n1+1, a2, b2, cheb_base_n2+n2+1, 2, g, func_amp)
    forall (ii = 1:n1+1, jj = 1:n2+1)
      TQn(ii, jj) = cmplx(g(1, ii, jj), g(2, ii, jj), NEC)
    end forall

  contains

    subroutine func_amp(x, y, rslts)

      real(NER), intent(in)                   :: x, y
      real(NER), dimension(:), intent(out)    :: rslts

      call get_allTQn_channel(self)
      call self%c_times_y_plus_b_times_x_plus_a(VM_lo, VM_1, x, VM_2, y, VM_tmp)
      call self%get_Tk_from_VM(VM_tmp, self%Tk)
      rslts(1) = real(self%Tk%ptr(self%N+1))
      rslts(2) = aimag(self%Tk%ptr(self%N+1))

    end subroutine func_amp

  end subroutine GetDoublePertV_sngl



!--------------------------------------------!
!                                            !
! Potential functions acting on obj_sngl     !
!                                            !
!--------------------------------------------!

  ! subroutine Vfunc_zero_sngl(osngl, p1, p2, potval)

  !   class(obj_sngl), intent(inout)    :: osngl
  !   real(NER), intent(in)          :: p1, p2
  !   real(NER), intent(out)         :: potval

  !   call osngl%check_warning_inited('Vfunc_zero_sngl')
  !   potval = 0.0_NER

  ! subroutine Vfunc_zero_sngl

  ! Vfunc_OPE_sngl needs regularization, so osngl must have been initialized.
  subroutine Vfunc_OPE_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval

    call osngl%check_warning_inited('Vfunc_OPE_sngl')
    if (osngl%lsj%S == 0) then
      potval = OPE_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval = OPE_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval = OPE_jpp(osngl%lsj%j, p1, p2)
        else
          potval = OPE_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_OPE_sngl

  ! ! Same as OPE except that the contact part in 1S0 is removed
  ! subroutine Vfunc_Yukawa_sngl(osngl, p1, p2, potval)

  !     class(obj_sngl), intent(inout)    :: osngl
  !     real(NER), intent(in)          :: p1, p2
  !     real(NER), intent(out)         :: potval

  !     call osngl%check_warning_inited('Vfunc_Yukawa_sngl')
  !     if (osngl%lsj%S == 0) then
  !         if (osngl%lsj%j == 0) then
  !             potval =
  !         else
  !             potval = OPE_j0j(osngl%lsj%j, p1, p2)
  !         end if
  !     else
  !         if (osngl%lsj%j == osngl%lsj%L) then
  !             potval = OPE_j1j(osngl%lsj%j, p1, p2)
  !         else
  !             if (osngl%lsj%L == osngl%lsj%j+1) then
  !                 potval = OPE_jpp(osngl%lsj%j, p1, p2)
  !             else
  !                 potval = OPE_jmm(osngl%lsj%j, p1, p2)
  !             end if
  !         end if
  !     end if

  !     if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
  !         potval = potval * osngl%regulator(p1, p2)

  ! end subroutine Vfunc_Yukawa_sngl


  ! Singular attractive potential
  subroutine Vfunc_sa_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval

    call pwd_coord(p1, p2, osngl%lsj%L, singattrpot_toy, potval)

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_sa_sngl


  subroutine Vfunc_VTPE0_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval

    call osngl%check_warning_inited('Vfunc_VTPE0_sngl')
    if (osngl%lsj%S == 0) then
      potval = VTPE0_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval = VTPE0_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval = VTPE0_jpp(osngl%lsj%j, p1, p2)
        else
          potval = VTPE0_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_VTPE0_sngl


  subroutine Vfunc_VTPE1_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: potval

    call osngl%check_warning_inited('Vfunc_VTPE1_sngl')
    if (osngl%lsj%S == 0) then
      potval = VTPE1_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval = VTPE1_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval = VTPE1_jpp(osngl%lsj%j, p1, p2)
        else
          potval = VTPE1_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_VTPE1_sngl


  ! potval = reg(p1)*reg(p2)*(p1*p2)^L

  subroutine Vfunc_VS0_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: potval

    integer :: L

    L = osngl%lsj%L
    select case (L)
      case (0)
        potval = 1.0_NER
      case (1)
        potval = p1*p2
      case default
        potval = (p1*p2)**L
    end select

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_VS0_sngl


  ! potval = reg(p1)*reg(p2)*(p1*p2)^L*(p1^2 + p2^2)

  subroutine Vfunc_VS1_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: potval

    call Vfunc_VS0_sngl(osngl, p1, p2, potval)
    potval = potval * (p1*p1 + p2*p2)

  end subroutine Vfunc_VS1_sngl


end module mod_obj_sngl
