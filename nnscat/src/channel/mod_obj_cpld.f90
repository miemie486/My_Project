! Bingwei Long 03/06/2018
! Bingwei Long 08/10/2013
! Class obj_cpld (coupled channels)

module mod_obj_cpld

  use mod_obj_channel
  use nneft_lsesv,    only: lsesv_2ch_half_cpld, lsesv_2ch_full_cpld
  use nneft_loops,    only: loop_1loop_os_VMTk_cpld, loop_1loop_VkTk_ab_cpld,  &
    & loop_2loops_TkVMTk_sym_cpld
  implicit none

  type VM_cpld
    logical :: created = .false., inited = .false.
    real(NER), dimension(:, :, :, :), pointer :: ptr => null()
  end type VM_cpld
  type TM_cpld
    logical :: created = .false., updated = .false.
    complex(NEC), dimension(:, :, :, :), pointer :: ptr => null()
  end type TM_cpld
  type Tk_cpld
    logical :: created = .false., updated = .false.
    complex(NEC), dimension(:, :, :), pointer :: ptr => null()
  end type Tk_cpld

  ! Wrapper structure for coupled-channel VM array and the potential that acts on it.
  type vstruct_cpld

    type(VM_cpld)                           :: VM
    procedure(Vfunc_cpld), pointer, nopass  :: pVfunc => null()

  end type vstruct_cpld

  ! Coupled channels; by default, power counting in Long & Yang ('12) is used.
  ! C appears at O(1), D at O(Q^2), etc.
  !
  !          ( (C0+C1+...) + (D0+D1+...)*(q^2 + p^2), (E0+E1+...)*p^2)
  ! VS = q^L (                                                       ) p^L
  !          ( (E0+E1+...)*q^2                      , 0              )
  !
  type, extends(obj_channel)  :: obj_cpld

    complex(NEC)              :: TQn(2, 2, MAX_NUMBER_ORDER)
    type(triplet_phs_cpld)    :: phsn(1:MAX_NUMBER_ORDER)

    ! Pre-defined long-range potential structures
    type(vstruct_cpld)  :: vs_vl0, vs_vl2, vs_vl3

    complex(NEC) :: CT_thrsh(1:2, 1:2), DT_thrsh(1:2, 1:2), ET_thrsh(1:2, 1:2)
    logical             :: CDET_threshold_set = .false.

    type(Tk_cpld)   :: Tk
    type(VM_cpld)   :: VMVLO
    procedure(Vfunc_cpld), pointer, nopass  :: pVfunc_VL0 => Vfunc_OPE_cpld
    procedure(Vfunc_cpld), pointer, nopass  :: pVfunc_VL2 => Vfunc_VTPE0_cpld
    procedure(Vfunc_cpld), pointer, nopass  :: pVfunc_VL3 => Vfunc_VTPE1_cpld

  contains

    procedure   :: init                        => init_cpld
    procedure   :: release                     => release_cpld
    procedure   :: update_k                    => update_k_cpld
    procedure   :: convert_TQn_to_phase_quiet  => convert_TQn_to_phase_quiet_cpld
    procedure   :: output_phase_to_string      => output_phase_to_string_cpld
    procedure   :: output_total_phase          => output_total_phase_cpld
    procedure   :: output_incremental_phase    => output_incremental_phase_cpld

    procedure   :: create_VM             => create_VM_cpld
    procedure   :: release_VM            => release_VM_cpld
    procedure   :: add_2VMs              => add_2VMs_cpld
    procedure   :: b_times_x_plus_a      => b_times_x_plus_a_cpld
    procedure   :: c_times_y_plus_b_times_x_plus_a => c_times_y_plus_b_times_x_plus_a_cpld
    procedure   :: init_symmtr_VM        => init_symmtr_VM_cpld
    procedure   :: update_symmtr_edge_VM => update_symmtr_edge_VM_cpld

    procedure   :: create_vstruct        => create_vstruct_cpld
    procedure   :: release_vstruct       => release_vstruct_cpld
    procedure   :: init_symmtr_vstruct   => init_symmtr_vstruct_cpld
    procedure   :: update_symmtr_edge_vstruct => update_symmtr_edge_vstruct_cpld

    procedure   :: create_TM             => create_TM_cpld
    procedure   :: release_TM            => release_TM_cpld

    procedure   :: create_Tk             => create_Tk_cpld
    ! procedure   :: FillTkByVM
    procedure   :: release_Tk            => release_Tk_cpld

    procedure   :: get_num_paras         => get_num_paras_cpld
    ! Frozen
    ! procedure   :: change_a_C            => change_a_C_cpld
    procedure   :: init_VMVLO            => init_VMVLO_cpld
    procedure   :: update_edge_VMVLO     => update_edge_VMVLO_cpld

    procedure   :: get_Tk_from_VM        => get_Tk_from_VM_cpld
    procedure   :: get_Tpq_from_VM       => get_Tpq_from_VM_cpld
    procedure   :: get_allTQn            => get_allTQn_cpld

    procedure   :: get_1loop_os_VMTk     => get_1loop_os_VMTk_cpld
    ! Frozen
    ! procedure   :: Get1LoopVMVM_OS       => Get1LoopVMVM_OS_cpld
    procedure   :: get_2loops_TkVMTk     => get_2loops_TkVMTk_cpld
    procedure   :: get_1loop_powTk       => get_1loop_powTk_cpld
    procedure   :: get_CDET              => get_CDET_cpld
    procedure   :: set_CDETthreshold     => set_CDETthreshold_cpld

    procedure   :: GetSinglePertV        => GetSinglePertV_cpld
    procedure   :: GetDoublePertV        => GetDoublePertV_cpld
  end type obj_cpld

  integer, parameter  :: npar_cpld(1:4) = (/1, 1, 4, 7/)

  ! Abstract interface for coupled-channel potentials
  abstract interface
    subroutine Vfunc_cpld(ocpld, p1, p2, potval)
      import  :: obj_cpld, NER
      implicit none
      class(obj_cpld), intent(inout)  :: ocpld
      real(NER), intent(in)           :: p1, p2
      real(NER), intent(out)          :: potval(1:2, 1:2)
    end subroutine
  end interface

contains

  subroutine init_cpld(self, N, regtype, Mambda)

    class(obj_cpld), intent(inout)  :: self
    integer, intent(in)             :: N, regtype
    real(NER), intent(in)           :: Mambda

    call init_channel(self, N, regtype, Mambda)
    call self%create_Tk(self%Tk)
    call self%create_vstruct(self%vs_vl0)
    call self%init_symmtr_vstruct(self%vs_vl0, self%pVfunc_VL0)
    call self%create_VM(self%VMVLO)
    call self%create_VM_sngl(self%VM_poly0)
    call self%init_VMpoly0(self%VM_poly0)

    if (self%uptoQn < 2) return
    call self%create_vstruct(self%vs_vl2)
    call self%init_symmtr_vstruct(self%vs_vl2, self%pVfunc_VL2)

    if (self%uptoQn < 3) return
    call self%create_vstruct(self%vs_vl3)
    call self%init_symmtr_vstruct(self%vs_vl3, self%pVfunc_VL3)
    if (self%uptoQn > 3)   &
      write (standard_error_unit, '(a)')  &
        'init_cpld: I do not know how to handle orders that are higher than 3. Order 3 is assumed.'

  end subroutine init_cpld

  subroutine release_cpld(self)

    class(obj_cpld), intent(inout)    :: self

    call self%release_Tk(self%Tk)
    call self%release_VM(self%VMVLO)
    call self%release_VM_sngl(self%VM_poly0)
    self%CDET_threshold_set = .false.
    if (associated(self%vs_vl0%VM%ptr)) call self%release_vstruct(self%vs_vl0)
    if (associated(self%vs_vl2%VM%ptr)) call self%release_vstruct(self%vs_vl2)
    if (associated(self%vs_vl3%VM%ptr)) call self%release_vstruct(self%vs_vl3)
    call release_channel(self)

  end subroutine release_cpld

  subroutine get_num_paras_cpld(self, num_cnttrm)

    class(obj_cpld), intent(inout)  :: self
    integer, intent(out)            :: num_cnttrm

    if (self%uptoQn > 3) then
      write (standard_error_unit, '(a, i2)')  &
        & 'get_num_paras_cpld: don''t know how to handle order ',   &
        & self%uptoQn
      num_cnttrm = 0
    else
      num_cnttrm = npar_cpld(self%uptoQn + 1)
    end if

    ! select case (self%uptoQn)
    !   case (0, 1)
    !     num_cnttrm = 1
    !   case (2)
    !     num_cnttrm = 4
    !   case (3)
    !     num_cnttrm = 7
    !   case default
    !     write (standard_error_unit, '(a, i2)')  &
    !       'get_num_paras_cpld: I don''t know how to handle order ', self%uptoQn
    !     num_cnttrm = 0
    ! end select

  end subroutine get_num_paras_cpld

  ! Frozen
  ! ! If self%Cs(1) is changed, VMVLO needs to be re-initialized.
  ! subroutine change_a_C_cpld(self, n, Cn)

  !     class(obj_cpld), intent(inout)   :: self
  !     integer, intent(in)                     :: n
  !     real(NER), intent(in)                   :: Cn

  !     call change_a_C_channel(self, n, Cn)
  !     if (n == 1) self%VMVLO%inited = .false.

  ! end subroutine change_a_C_cpld

  subroutine init_VMVLO_cpld(self, C0)

    class(obj_cpld), intent(inout)   :: self
    real(NER), intent(in)            :: C0

    integer :: ii, jj, locN

    if (.not. self%VMVLO%created) then
      write (standard_error_unit, '(a)')  &
        'init_VMVLO_cpld: VMVLO has yet been created. Nothing will be done.'
      return
    end if
    locN = self%N
    do ii = 1, locN
      do jj = 1, ii - 1
        self%VMVLO%ptr(1:2, 1:2, ii, jj) = self%vs_vl0%VM%ptr(1:2, 1:2, ii, jj)
        self%VMVLO%ptr(1, 1, ii, jj) = self%VMVLO%ptr(1, 1, ii, jj) + C0*self%VM_poly0%ptr(ii, jj)
        self%VMVLO%ptr(1:2, 1:2, jj, ii) = transpose(self%VMVLO%ptr(1:2, 1:2, ii, jj))
      end do
      self%VMVLO%ptr(1:2, 1:2, ii, ii) = self%vs_vl0%VM%ptr(1:2, 1:2, ii, ii)
      self%VMVLO%ptr(1, 1, ii, ii) = self%VMVLO%ptr(1, 1, ii, ii) + C0*self%VM_poly0%ptr(ii, ii)
    end do
    self%VMVLO%inited = .true.

  end subroutine init_VMVLO_cpld

  subroutine update_edge_VMVLO_cpld(self, k, C0)

    class(obj_cpld), intent(inout)          :: self
    real(NER), intent(in)                   :: k
    real(NER), intent(in)                   :: C0

    integer :: locN, ii

    if (.not. self%VMVLO%created) then
      write (standard_error_unit, '(a)')  &
        'update_edge_VMVLO_cpld: VMVLO has yet been created. Nothing will be done.'
      return
    end if
    locN = self%N
    do ii = 1, locN
      self%VMVLO%ptr(1:2, 1:2, ii, locN+1) = self%vs_vl0%VM%ptr(1:2, 1:2, ii, locN+1)
      self%VMVLO%ptr(1, 1, ii, locN+1) = self%VMVLO%ptr(1, 1, ii, locN+1) + C0*self%VM_poly0%ptr(ii, locN+1)
      self%VMVLO%ptr(1:2, 1:2, locN+1, ii) = transpose(self%VMVLO%ptr(1:2, 1:2, ii, locN+1))
    end do
    self%VMVLO%ptr(1:2, 1:2, locN+1, locN+1) = self%vs_vl0%VM%ptr(1:2, 1:2, locN+1, locN+1)
    self%VMVLO%ptr(1, 1, locN+1, locN+1) = self%VMVLO%ptr(1, 1, locN+1, locN+1) + C0*self%VM_poly0%ptr(locN+1, locN+1)

  end subroutine update_edge_VMVLO_cpld

  subroutine update_k_cpld(self, k)

    class(obj_cpld), intent(inout)    :: self
    real(NER), intent(in)             :: k

    call update_k_channel(self, k)
    self%TQn_updated  = .false.
    self%phsn_updated = .false.
    call self%update_symmtr_edge_vstruct(self%vs_vl0)
    call self%update_edge_VMpoly0(self%VM_poly0)
    ! C0 = self%Cs(1)
    call self%update_edge_VMVLO(k, self%Cs(1))
    if (self%vs_vl2%VM%inited) then
      call self%update_symmtr_edge_vstruct(self%vs_vl2)
    end if
    if (self%vs_vl3%VM%inited) then
      call self%update_symmtr_edge_vstruct(self%vs_vl3)
    end if
    self%Tk%updated = .false.

  end subroutine update_k_cpld

  subroutine set_CDETthreshold_cpld(self)

    class(obj_cpld), intent(inout)  :: self

    integer         :: save_uptoQn
    real(NER)       :: zero_k, save_k
    real(NER), parameter    :: local_EPS = 1.0E-6_NER

    save_uptoQn = self%uptoQn
    self%uptoQn = 0
    save_k = self%k
    zero_k = local_EPS*self%Mambda
    call self%update_k(zero_k)
    call self%get_allTQn()
    call self%get_CDET(self%Tk, self%CT_thrsh, self%DT_thrsh, self%ET_thrsh)
    self%CDET_threshold_set = .true.
    self%uptoQn = save_uptoQn
    self%k = save_k

  end subroutine set_CDETthreshold_cpld

  ! Convert TQn to phase shifts. uflag are the unitarity flags for each order
  ! Assuming that LO is nonperturbative.
  ! In perturbative OPE calculations, this method needs to be overriden.
  subroutine convert_TQn_to_phase_quiet_cpld(self)

    class(obj_cpld), intent(inout)      :: self
    ! integer, dimension(:), intent(out)  :: uflag

    integer                         :: ii
    complex(NEC), dimension(2, 2)   :: violation

    if (.not. self%TQn_updated) then
      write (standard_error_unit, '(a)')  &
        'convert_TQn_to_phase_quiet_cpld: Warning: TQn is not updated.'
      return
    end if

    call t2d1d2e_np(self%k, self%TQn(1:2, 1:2, 1), self%phsn(1), self%uflag(1))
    ! for 3S1, adjust d1 and mixing angle
    if (self%chnid .eq. CHNID_3S1_PP) then
      if (self%k/PC_mpi .lt. 1.4_NER .and. self%phsn(1)%d1 .lt. 0.0_NER) then
        self%phsn(1)%d1 = self%phsn(1)%d1 + 180.0_NER
        self%phsn(1)%e = -self%phsn(1)%e
      end if
    end if
    ! for 3P2, adjust d1
    if (self%chnid .eq. CHNID_3P2_PP) then
      if (self%phsn(1)%d1 .lt. 0.0_NER) then
        self%phsn(1)%d1 = self%phsn(1)%d1 + 180.0_NER
        self%phsn(1)%e = -self%phsn(1)%e
      end if
    end if
    do ii = 1, self%uptoQn
      call get_nth_delta_cpld(ii, self%TQn(1:2, 1:2, ii+1), self%phsn(1:ii), self%k, self%uflag(ii+1), violation, self%phsn(ii+1))
      self%univio(ii+1) = abs(violation(1, 1)) + abs(violation(1, 2))          &
        & + abs(violation(2, 1)) + abs(violation(2, 2))
    end do
    self%phsn_updated = .true.
    call convert_TQn_to_phase_quiet_channel(self)

  end subroutine convert_TQn_to_phase_quiet_cpld


  subroutine output_phase_to_string_cpld(self, outstr)

    class(obj_cpld), intent(inout)  :: self
    character(len = *), intent(out) :: outstr

    call output_phase_to_string_channel(self, outstr)
    call write_incremental_phsn_to_text(self%uptoQn+1, self%phsn, outstr)

  end subroutine output_phase_to_string_cpld


  subroutine output_total_phase_cpld(self, out_phs)

    class(obj_cpld), intent(inout)  :: self
    class(*), intent(inout)         :: out_phs

    integer :: ii

    call output_total_phase_channel(self, out_phs)
    select type (out_phs)
      type is (triplet_phs_cpld)
        out_phs%d1 = 0.0_NER
        out_phs%d2 = 0.0_NER
        out_phs%e = 0.0_NER
        do ii = 1, self%uptoQn+1
          out_phs%d1 = out_phs%d1 + self%phsn(ii)%d1
          out_phs%d2 = out_phs%d2 + self%phsn(ii)%d2
          out_phs%e = out_phs%e + self%phsn(ii)%e
        end do
      class default
        write (standard_error_unit, '(a)')  &
          'output_total_phase_cpld: The type of output phase must be triplet_phs_cpld. The returned value is un-defined.'
    end select

  end subroutine output_total_phase_cpld

  subroutine output_incremental_phase_cpld(self, out_phs)

    class(obj_cpld), intent(inout)  :: self
    class(*), intent(inout)         :: out_phs(:)

    integer :: ii

    call output_incremental_phase_channel(self, out_phs)
    select type (out_phs)
      type is (triplet_phs_cpld)
        do ii = 1, self%uptoQn+1
          out_phs(ii)%d1 = self%phsn(ii)%d1
          out_phs(ii)%d2 = self%phsn(ii)%d2
          out_phs(ii)%e = self%phsn(ii)%e
        end do
      class default
        write (standard_error_unit, '(a)')  &
          'output_incremental_phase_cpld: The type of output phase must be triplet_phs_cpld. The returned value is un-defined.'
    end select

  end subroutine output_incremental_phase_cpld

  ! create_VM_cpld will initialize all elements to zero
  subroutine create_VM_cpld(self, VM)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(inout)    :: VM

    call self%check_warning_inited('create_VM_cpld')
    if (VM%created)   &
      write (standard_error_unit, '(a)') 'create_VM_cpld: Warning! Trying to create a VM that has already been created'
    if (associated(VM%ptr)) deallocate(VM%ptr)
    allocate(VM%ptr(1:2, 1:2, 1:self%N+1, 1:self%N+1))
    VM%ptr(1:2, 1:2, 1:self%N+1, 1:self%N+1) = 0.0_NER
    VM%created = .true.
    VM%inited  = .false.

  end subroutine create_VM_cpld

  subroutine release_VM_cpld(self, VM)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(inout)    :: VM

    if (.not. associated(VM%ptr)) then
      write (standard_error_unit, '(a)')  &
        'release_VM_cpld: Warning: No need to release.'
    else
      deallocate(VM%ptr)
    endif
    VM%created = .false.
    VM%inited  = .false.

  end subroutine release_VM_cpld


  ! Add VM1 and VM2 up to fill VMdest. If k_defined is true, also act on the
  ! edges of the VMs.

  subroutine add_2VMs_cpld(self, VM1, VM2, VMdest)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(inout)    :: VM1, VM2, VMdest

    integer :: ii, jj

    call self%check_warning_inited('add_2VMs_cpld')
    do ii = 1, self%N
      do jj = 1, self%N
        VMdest%ptr(1:2, 1:2, ii, jj) = VM1%ptr(1:2, 1:2, ii, jj) + VM2%ptr(1:2, 1:2, ii, jj)
      end do
    end do
    if (self%k_defined) then
      do ii = 1, self%N
        VMdest%ptr(1:2, 1:2, ii, self%N+1) = VM1%ptr(1:2, 1:2, ii, self%N+1) + VM2%ptr(1:2, 1:2, ii, self%N+1)
        VMdest%ptr(1:2, 1:2, self%N+1, ii) = VM1%ptr(1:2, 1:2, self%N+1, ii) + VM2%ptr(1:2, 1:2, self%N+1, ii)
      end do
      VMdest%ptr(1:2, 1:2, self%N+1, self%N+1) = VM1%ptr(1:2, 1:2, self%N+1, self%N+1) + VM2%ptr(1:2, 1:2, self%N+1, self%N+1)
    end if
    VMdest%inited = .true.

  end subroutine add_2VMs_cpld


  ! Return VMa + x*VMb
  ! If VMa is not inited, then return x*VMb
  subroutine b_times_x_plus_a_cpld(self, VMa, VMb, x, VMrslt)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(inout)    :: VMa, VMb
    real(NEC), intent(in)           :: x
    type(VM_cpld), intent(inout)    :: VMrslt

    integer :: n

    n = self%N
    if (VMa%inited) then
      VMrslt%ptr(1:2, 1:2, 1:n+1, 1:n+1) = VMa%ptr(1:2, 1:2, 1:n+1, 1:n+1) + x* VMb%ptr(1:2, 1:2, 1:n+1, 1:n+1)
    else
      VMrslt%ptr(1:2, 1:2, 1:n+1, 1:n+1) = x* VMb%ptr(1:2, 1:2, 1:n+1, 1:n+1)
    end if
    VMrslt%inited = .true.

  end subroutine b_times_x_plus_a_cpld

  ! Return VMa + x*VMb + y*VMc
  ! If VMa is not inited, then return x*VMb + y*VMc
  subroutine c_times_y_plus_b_times_x_plus_a_cpld(self, VMa, VMb, x, VMc, y, VMrslt)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(inout)    :: VMa, VMb, VMc
    real(NEC), intent(in)           :: x, y
    type(VM_cpld), intent(inout)    :: VMrslt

    integer :: n

    n = self%N
    if (VMa%inited) then
      VMrslt%ptr(1:2, 1:2, 1:n+1, 1:n+1) = VMa%ptr(1:2, 1:2, 1:n+1, 1:n+1) + x* VMb%ptr(1:2, 1:2, 1:n+1, 1:n+1) + y* VMc%ptr(1:2, 1:2, 1:n+1, 1:n+1)
    else
      VMrslt%ptr(1:2, 1:2, 1:n+1, 1:n+1) = x* VMb%ptr(1:2, 1:2, 1:n+1, 1:n+1) + y* VMc%ptr(1:2, 1:2, 1:n+1, 1:n+1)
    end if
    VMrslt%inited = .true.

  end subroutine c_times_y_plus_b_times_x_plus_a_cpld


  subroutine init_symmtr_VM_cpld(self, Vfunc, VM)

    class(obj_cpld), intent(inout)  :: self
    procedure(Vfunc_cpld)           :: Vfunc
    type(VM_cpld), intent(inout)    :: VM

    integer     :: ii, jj
    real(NER)   :: ppi, ppj
    real(NER), dimension(2, 2)  :: vab

    call self%check_warning_inited('init_symmtr_VM_cpld')
    if (.not. VM%created) then
      write (standard_error_unit, '(a)') 'init_symmtr_VM_cpld: VM has not been created.'
      return
    end if
    do ii = 1, self%N
      ppi = self%msh(ii)
      do jj = 1, ii - 1
        ppj = self%msh(jj)
        call Vfunc(self, ppi, ppj, vab)
        VM%ptr(1:2, 1:2, ii, jj) = vab
        VM%ptr(1:2, 1:2, jj, ii) = transpose(vab)
      end do
      ! set up the diagonal elements
      call Vfunc(self, ppi, ppi, VM%ptr(1:2, 1:2, ii, ii))
    end do
    VM%inited = .true.

  end subroutine init_symmtr_VM_cpld


  subroutine update_symmtr_edge_VM_cpld(self, Vfunc, VM)

    class(obj_cpld), intent(inout)  :: self
    procedure(Vfunc_cpld)           :: Vfunc
    type(VM_cpld), intent(inout)    :: VM

    integer     :: jj, N
    real(NER)   :: ppj, k
    real(NER), dimension(2, 2)  :: vab

    if (.not. VM%created) then
      write (standard_error_unit, '(a)')  &
        'update_symmtr_edge_VM_cpld: VM has not been created. Nothing will be done.'
      return
    end if
    if (.not. (self%k_defined)) then
      write(standard_error_unit, '(a)')   &
        'update_symmtr_edge_VM_cpld: k is not defined. Nothing will be done.'
      return
    end if
    k = self%k
    N = self%N
    do jj = 1, N
      ppj = self%msh(jj)
      call Vfunc(self, k, ppj, vab)
      VM%ptr(1:2, 1:2, N+1, jj) = vab
      VM%ptr(1:2, 1:2, jj, N+1) = transpose(vab)
    end do
    call Vfunc(self, k, k, VM%ptr(1:2, 1:2, N+1, N+1))

  end subroutine update_symmtr_edge_VM_cpld


  subroutine create_vstruct_cpld(self, vs)

    class(obj_cpld), intent(inout)      :: self
    type(vstruct_cpld), intent(inout)   :: vs

    call self%create_VM(vs%VM)

  end subroutine create_vstruct_cpld


  subroutine release_vstruct_cpld(self, vs)

    class(obj_cpld), intent(inout)      :: self
    type(vstruct_cpld), intent(inout)   :: vs

    call self%release_VM(vs%VM)
    vs%pVfunc => null()

  end subroutine release_vstruct_cpld


  subroutine init_symmtr_vstruct_cpld(self, vs, Vfunc)

    class(obj_cpld), intent(inout)      :: self
    type(vstruct_cpld), intent(inout)   :: vs
    procedure(Vfunc_cpld)               :: Vfunc

    vs%pVfunc => Vfunc
    call self%init_symmtr_VM(vs%pVfunc, vs%VM)

  end subroutine init_symmtr_vstruct_cpld


  subroutine update_symmtr_edge_vstruct_cpld(self, vs)

    class(obj_cpld), intent(inout)      :: self
    type(vstruct_cpld), intent(inout)   :: vs

    call self%update_symmtr_edge_VM(vs%pVfunc, vs%VM)

  end subroutine update_symmtr_edge_vstruct_cpld


  subroutine create_TM_cpld(self, TM)

    class(obj_cpld), intent(inout)  :: self
    type(TM_cpld), intent(inout)    :: TM

    call self%check_warning_inited('create_TM_cpld')
    if (TM%created)   &
      write (standard_error_unit, '(a)')  &
        'create_TM_cpld: Warning! Trying to create a TM that has already been created'
    if (associated(TM%ptr)) deallocate(TM%ptr)
    allocate(TM%ptr(1:2, 1:2, 1:self%N+1, 1:self%N+1))
    TM%created = .true.

  end subroutine create_TM_cpld


  subroutine release_TM_cpld(self, TM)

    class(obj_cpld), intent(inout)  :: self
    type(TM_cpld), intent(inout)    :: TM

    call self%check_warning_inited('release_TM_cpld')
    if (.not. associated(TM%ptr)) then
      write (standard_error_unit, '(a)') 'release_TM_cpld: No need to release'
    else
      deallocate(TM%ptr)
    endif
    TM%created = .false.
    TM%updated = .false.

  end subroutine release_TM_cpld


  subroutine create_Tk_cpld(self, Tk)

    class(obj_cpld), intent(inout)  :: self
    type(Tk_cpld), intent(inout)    :: Tk

    call self%check_warning_inited('create_Tk_cpld')
    if (Tk%created)   &
      write (standard_error_unit, '(a)')  &
        'create_Tk_cpld: Warning! Trying to create a Tk that has already been created'
    if (associated(Tk%ptr)) deallocate(Tk%ptr)
    allocate(Tk%ptr(1:2, 1:2, 1:self%N+1))
    Tk%created = .true.
    Tk%updated = .false.

  end subroutine create_Tk_cpld


  subroutine release_Tk_cpld(self, Tk)

    class(obj_cpld), intent(inout)  :: self
    type(Tk_cpld), intent(inout)    :: Tk

    call self%check_warning_inited('release_Tk_cpld')
    if (.not. associated(Tk%ptr)) then
      write (standard_error_unit, '(a)') 'release_Tk_cpld: No need to release'
    else
      deallocate(Tk%ptr)
    endif
    Tk%created = .false.
    Tk%updated = .false.

  end subroutine release_Tk_cpld


  subroutine get_Tk_from_VM_cpld(self, VM, Tk)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(in)       :: VM
    type(Tk_cpld), intent(inout)    :: Tk

    integer                         :: flag

    if (.not. self%k_defined) then
      write(standard_error_unit, '(a)')   &
        'get_Tk_from_VM_cpld: k not defined. I will not do anything.'
      return
    end if
    if (.not. VM%inited) then
      write(standard_error_unit, '(a)')   &
        'get_Tk_from_VM_cpld: VM not initialized. I will not do anything.'
      return
    end if
    if (.not. Tk%created) then
      write(standard_error_unit, '(a)')   &
        'get_Tk_from_VM_cpld: Tk not created. I will not do anything.'
      return
    end if
    call lsesv_2ch_half_cpld(self%k, self%N, VM%ptr, self%msh, self%wght, self%lsqr, Tk%ptr, flag)
    Tk%updated = .true.

  end subroutine get_Tk_from_VM_cpld

  subroutine get_Tpq_from_VM_cpld(self, VM, TM)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(in)       :: VM
    type(TM_cpld), intent(inout)    :: TM

    integer                         :: flag

    if (.not. self%k_defined) then
      write(standard_error_unit, '(a)')   &
        'get_Tpq_from_VM_cpld: k not defined. I will not do anything.'
      return
    end if
    if (.not. VM%inited) then
      write(standard_error_unit, '(a)')   &
        'get_Tpq_from_VM_cpld: VM not initialized. I will not do anything.'
      return
    end if
    if (.not. TM%created) then
      write(standard_error_unit, '(a)')   &
        'get_Tpq_from_VM_cpld: TM not created. I will not do anything.'
      return
    end if
    call lsesv_2ch_full_cpld(self%k, self%N, VM%ptr, self%msh, self%wght, self%lsqr, TM%ptr, flag)
    TM%updated = .true.

  end subroutine get_Tpq_from_VM_cpld

  subroutine get_allTQn_cpld(self)

    class(obj_cpld), intent(inout)   :: self

    integer                         :: N
    real(NER), dimension(2, 2)      :: VL
    complex(NEC), dimension(2, 2)   :: L1, L2, CT, DT, ET

    call get_allTQn_channel(self)
    if (.not. self%VMVLO%inited) call self%init_VMVLO(self%Cs(1))
    call self%get_Tk_from_VM(self%VMVLO, self%Tk)
    N = self%N
    self%TQn(1:2, 1:2, 1) = self%Tk%ptr(1:2, 1:2, N+1)
    if (self%uptoQn == 0) then
      self%TQn_updated = .true.
      return
    end if

    self%TQn(1:2, 1:2, 2) = 0.0_NER
    if (self%uptoQn == 1) then
      self%TQn_updated = .true.
      return
    end if

    VL(1:2, 1:2) = self%vs_vl2%VM%ptr(1:2, 1:2, N+1, N+1)
    call self%get_1loop_os_VMTk(self%vs_vl2%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl2%VM, self%Tk, L2)
    call self%get_CDET(self%Tk, CT, DT, ET)
    self%TQn(1:2, 1:2, 3) = VL + L1 + transpose(L1) + L2 + self%Cs(2)*CT + self%Cs(3)*DT  &
      + self%Cs(4)*ET
    if (self%uptoQn == 2) then
      self%TQn_updated = .true.
      return
    end if

    VL(1:2, 1:2) = self%vs_vl3%VM%ptr(1:2, 1:2, N+1, N+1)
    call self%get_1loop_os_VMTk(self%vs_vl3%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl3%VM, self%Tk, L2)
    self%TQn(1:2, 1:2, 4) = VL + L1 + transpose(L1) + L2 + self%Cs(5)*CT + self%Cs(6)*DT  &
      + self%Cs(7)*ET
    if (self%uptoQn == 3) then
      self%TQn_updated = .true.
      return
    end if

    write (standard_error_unit, '(a)')  &
    & 'get_allTQn_cpld: don''t know how to do uptoQn > 3'

  end subroutine get_allTQn_cpld

  subroutine get_1loop_os_VMTk_cpld(self, VM, Tk, loopval)

    class(obj_cpld), intent(inout) :: self
    type(VM_cpld), intent(in)   :: VM
    type(Tk_cpld), intent(in)   :: Tk
    complex(NEC), dimension(1:2, 1:2), intent(out) :: loopval

    call self%check_warning_inited('get_1loop_os_VMTk_cpld')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_1loop_os_VMTk_cpld: k not defined. Nothing will be done.'
      return
    end if
    if (.not. VM%inited) &
      write (standard_error_unit, '(a)')  &
        'get_1loop_os_VMTk_cpld: Warning! VM not initialized.  Nothing will be done.'
    if (.not. Tk%updated) &
      write (standard_error_unit, '(a)')  &
        'get_1loop_os_VMTk_cpld: Warning! Tk not defined. Nothing will be done.'
    call loop_1loop_os_VMTk_cpld(self%k, self%N, self%msh, self%wght, self%lsqr, VM%ptr, Tk%ptr, loopval)

  end subroutine get_1loop_os_VMTk_cpld

  ! Frozen
  ! subroutine Get1LoopVMVM_OS_cpld(self, VMa, VMb, loopval)

  !   class(obj_cpld), intent(inout) :: self
  !   type(VM_cpld), intent(in)   :: VMa, VMb
  !   complex(NEC), intent(out)   :: loopval(1:2, 1:2)

  !   call self%check_warning_inited('Get1LoopVMVM_OS_cpld')
  !   if (.not. self%k_defined) then
  !     write (standard_error_unit, '(a)')  &
  !       'Get1LoopVMVM_OS_cpld: k not defined. Nothing will be done.'
  !     return
  !   end if
  !   if ((.not. VMa%inited) .or. (.not. VMb%inited)) &
  !     write (standard_error_unit, '(a)')  &
  !       'Get1LoopVMVM_OS_cpld: Warning! VM not initialized.  Nothing will be done.'
    ! call Get1LoopVMVMOSCpld_loop(self%k, self%N, self%msh, self%wght, self%lsqr, VMa%ptr, VMa%ptr, loopval)

  ! end subroutine Get1LoopVMVM_OS_cpld

  subroutine get_2loops_TkVMTk_cpld(self, VM, Tk, loopval)

    class(obj_cpld), intent(inout) :: self
    type(VM_cpld), intent(in)   :: VM
    type(Tk_cpld), intent(in)   :: Tk
    complex(NEC), dimension(1:2, 1:2), intent(out) :: loopval

    call self%check_warning_inited('get_2loops_TkVMTk_cpld')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_2loops_TkVMTk_cpld: k not defined. Nothing will be done.'
      return
    end if
    if (.not. VM%inited) &
      write (standard_error_unit, '(a)')  &
        'get_2loops_TkVMTk_cpld: Warning! VM not initialized. Nothing will be done.'
    if (.not. Tk%updated) &
      write (standard_error_unit, '(a)')  &
        'get_2loops_TkVMTk_cpld: Warning! Tk not defined. Nothing will be done.'
    call loop_2loops_TkVMTk_sym_cpld(self%k, self%N, self%msh, self%wght, self%lsqr, VM%ptr, Tk%ptr, loopval)

  end subroutine get_2loops_TkVMTk_cpld

  subroutine get_1loop_powTk_cpld(self, pow_arrayID, Tk, loopval)

    class(obj_cpld), intent(inout) :: self
    integer, intent(in)            :: pow_arrayID
    type(Tk_cpld), intent(in)      :: Tk
    complex(NEC), dimension(1:2, 1:2), intent(out)  :: loopval

    integer :: aa, bb
    real(NER), dimension(:), pointer    :: ptrVk

    call self%check_warning_inited('get_1loop_powTk_cpld')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_1loop_powTk_cpld: k not defined. Nothing will be done.'
      return
    end if
    if (.not. Tk%updated) &
      write (standard_error_unit, '(a)')  &
        'get_1loop_powTk_cpld: Warning! Tk not defined. Nothing will be done.'

    call self%set_ptr_to_powfac(pow_arrayID, ptrVk)
    do aa = 1, 2
      do bb = 1, 2
        loopval(aa, bb) = loop_1loop_VkTk_ab_cpld(self%k, self%N, self%msh, self%wght, self%lsqr, ptrVk, Tk%ptr, aa, bb)
      end do
    end do

  end subroutine get_1loop_powTk_cpld

  subroutine get_CDET_cpld(self, Tk, CT, DT, ET)

    class(obj_cpld), intent(inout)  :: self
    type(Tk_cpld), intent(in)       :: Tk
    complex(NEC), dimension(2, 2), intent(out) :: CT, DT, ET

    complex(NEC), dimension(2, 2)       :: F, F2
    real(NER)                           :: ksqr, regfac, regsqr
    real(NER), dimension(:), pointer    :: ptr_powfac

    if (.not. Tk%updated) then
      write (standard_error_unit, '(a)')  &
        'get_CDET_cpld: Tk has not been updated.'
      CT = 0.0
      DT = 0.0
      ET = 0.0
      return
    end if
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        'get_CDET_cpld: k is not defined. Nothing will be done.'
      return
    end if

    call self%set_ptr_to_powfac(1, ptr_powfac)
    regfac = ptr_powfac(self%N+1)
    regsqr = regfac*regfac
    ksqr   = self%k*self%k

    call self%get_1loop_powTk(1, Tk, F)
    call self%get_1loop_powTk(2, Tk, F2)

    CT(1, 1) = regsqr + 2.0*regfac*F(1,1) + F(1,1)*F(1,1)
    CT(1, 2) = regfac*F(1,2) + F(1,1)*F(1,2)
    CT(2, 1) = CT(1, 2)
    CT(2, 2) = F(1, 2)*F(1, 2)
    DT(1, 1) = 2.0*(regsqr*ksqr + regfac*(ksqr*F(1,1) + F2(1,1)) + F2(1,1)*F(1,1))
    DT(1, 2) = regfac*(ksqr*F(1,2) + F2(1,2)) + F(1,1)*F2(1,2) + F(1,2)*F2(1,1)
    DT(2, 1) = DT(1, 2)
    DT(2, 2) = 2.0*F(1,2)*F2(1,2)
    ET(1, 1) = 2.0*(regfac*F2(2,1) + F(1,1)*F2(2,1))
    ET(1, 2) = regsqr*ksqr + regfac*(ksqr*F(1,1) + F2(2,2)) + F(1,2)*F2(2,1) + F(1,1)*F2(2,2)
    ET(2, 1) = ET(1, 2)
    ET(2, 2) = 2.0*(regfac*ksqr*F(1,2) + F(1,2)*F2(2,2))

  end subroutine get_CDET_cpld

  ! n insertions of VM_pert. a and b are the boundary for the dimensionless
  ! parameter x, by which the Taylor expansion is defined (see
  ! cheb_taylor_real). cheb_base+n+1 is the actual number of insertions
  ! calculated by cheb_taylor subroutines. TQn(i) is the value with (i-1)
  ! insertions.
  ! Note: TQ(1) is expected to be zero if the accuracy is infinite, so its
  ! returned value serves as a metric for actual numerical accuracy.
  subroutine GetSinglePertV_cpld(self, VM_lo, n, VM_pert, VM_tmp, TQn, a, b, cheb_base_n)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(inout)    :: VM_lo, VM_pert, VM_tmp
    complex(NEC), intent(out)       :: TQn(:, :, :)
    integer, intent(in)             :: n, cheb_base_n
    real(NER), intent(in)           :: a, b

    real(NER) :: g(6, 1:cheb_base_n+n+5)
    integer   :: ii

    call cheb_taylor(a, b, 6, cheb_base_n+n+1, g, func_amp)
    do ii = 1, n+1
      TQn(1, 1, ii) = cmplx(g(1, ii), g(2, ii), NEC)
      TQn(2, 2, ii) = cmplx(g(3, ii), g(4, ii), NEC)
      TQn(1, 2, ii) = cmplx(g(5, ii), g(6, ii), NEC)
      TQn(2, 1, ii) = TQn(1, 2, ii)
    end do

  contains

    subroutine func_amp(x, fnval)

      real(NER), intent(in) :: x
      real(NER), dimension(:), intent(out) :: fnval

      call get_allTQn_channel(self)
      call self%b_times_x_plus_a(VM_lo, VM_pert, x, VM_tmp)
      call self%get_Tk_from_VM(VM_tmp, self%Tk)
      fnval(1) = real(self%Tk%ptr(1, 1, self%N+1))
      fnval(2) = aimag(self%Tk%ptr(1, 1, self%N+1))
      fnval(3) = real(self%Tk%ptr(2, 2, self%N+1))
      fnval(4) = aimag(self%Tk%ptr(2, 2, self%N+1))
      fnval(5) = real(self%Tk%ptr(1, 2, self%N+1))
      fnval(6) = aimag(self%Tk%ptr(1, 2, self%N+1))

    end subroutine func_amp

  end subroutine GetSinglePertV_cpld

  ! TQn(1:2, 1:2, i, j) = (i-1) insertions of v1 and (j-1) insertions of v2
  ! Note: TQ(1, 1) is expected to be zero if the accuracy is infinite, so its
  ! returned value serves as a metric for actual numerical accuracy.
  subroutine GetDoublePertV_cpld(self, VM_lo, n1, VM_1, n2, VM_2, VM_tmp, TQn, &
    & a1, b1, cheb_base_n1, a2, b2, cheb_base_n2)

    class(obj_cpld), intent(inout)  :: self
    type(VM_cpld), intent(inout)    :: VM_lo, VM_1, VM_2, VM_tmp
    complex(NEC), intent(out)       :: TQn(:, :, :, :)
    integer, intent(in)             :: n1, n2, cheb_base_n1, cheb_base_n2
    real(NER), intent(in)           :: a1, b1, a2, b2

    real(NER) :: g(1:6, 1:cheb_base_n1+n1+5, 1:cheb_base_n2+n2+5)
    integer   :: ii, jj

    call cheb_taylor(a1, b1, cheb_base_n1+n1+1, a2, b2, cheb_base_n2+n2+1, 6, g, func_amp)
    forall (ii = 1:n1+1, jj = 1:n2+1)
      TQn(1, 1, ii, jj) = cmplx(g(1, ii, jj), g(2, ii, jj), NEC)
      TQn(2, 2, ii, jj) = cmplx(g(3, ii, jj), g(4, ii, jj), NEC)
      TQn(1, 2, ii, jj) = cmplx(g(5, ii, jj), g(6, ii, jj), NEC)
      TQn(2, 1, ii, jj) = TQn(1, 2, ii, jj)
    end forall

  contains

    subroutine func_amp(x, y, rslts)

      real(NER), intent(in)                   :: x, y
      real(NER), dimension(:), intent(out)    :: rslts

      call get_allTQn_channel(self)
      call self%c_times_y_plus_b_times_x_plus_a(VM_lo, VM_1, x, VM_2, y, VM_tmp)
      call self%get_Tk_from_VM(VM_tmp, self%Tk)
      rslts(1) = real(self%Tk%ptr(1, 1, self%N+1))
      rslts(2) = aimag(self%Tk%ptr(1, 1, self%N+1))
      rslts(3) = real(self%Tk%ptr(2, 2, self%N+1))
      rslts(4) = aimag(self%Tk%ptr(2, 2, self%N+1))
      rslts(5) = real(self%Tk%ptr(1, 2, self%N+1))
      rslts(6) = aimag(self%Tk%ptr(1, 2, self%N+1))

    end subroutine func_amp

  end subroutine GetDoublePertV_cpld


!--------------------------------------------!
!                                            !
! Potential functions for coupled channels   !
!                                            !
!--------------------------------------------!

  ! < p1 | OPE | p2 >
  subroutine Vfunc_OPE_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)              :: ocpld
    real(NER), intent(in)                       :: p1, p2
    real(NER), dimension(1:2, 1:2), intent(out) :: vab

    vab = 0.0
    if (ocpld%lsj%j .eq. 0 .or. ocpld%lsj%L .ne. ocpld%lsj%j-1) return
    vab(1, 1) = OPE_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) = OPE_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) = OPE_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) = OPE_jpp(ocpld%lsj%j, p1, p2)
    if (ocpld%regtype .ne. REGTYPE_NONE .and. ocpld%regtype .ne. REGTYPE_SHARP)   &
      vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_OPE_cpld

  ! < p1 | VTPE0 | p2 >
  ! delta-less TPE_0
  subroutine Vfunc_VTPE0_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)                 :: ocpld
    real(NER), intent(in)                       :: p1, p2
    real(NER), dimension(1:2, 1:2), intent(out) :: vab

    vab = 0.0_NER
    if (ocpld%lsj%j .eq. 0 .or. ocpld%lsj%L .ne. ocpld%lsj%j-1) return
    vab(1, 1) = VTPE0_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) = VTPE0_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) = VTPE0_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) = VTPE0_jpp(ocpld%lsj%j, p1, p2)
    if (ocpld%regtype .ne. REGTYPE_NONE .and. ocpld%regtype .ne. REGTYPE_SHARP)   &
      vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_VTPE0_cpld

  ! < p1 | VTPE1 | p2 >
  ! delta-less TPE_1
  subroutine Vfunc_VTPE1_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)                 :: ocpld
    real(NER), intent(in)                       :: p1, p2
    real(NER), dimension(1:2, 1:2), intent(out) :: vab

    vab = 0.0_NER
    if (ocpld%lsj%j .eq. 0 .or. ocpld%lsj%L .ne. ocpld%lsj%j-1) return
    vab(1, 1) = VTPE1_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) = VTPE1_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) = VTPE1_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) = VTPE1_jpp(ocpld%lsj%j, p1, p2)
    if (ocpld%regtype .ne. REGTYPE_NONE .and. ocpld%regtype .ne. REGTYPE_SHARP)   &
      vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_VTPE1_cpld

  ! vab(1, 1) = reg(p1)*reg(p2)*(p1*p2)^L
  subroutine Vfunc_VS0_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)                 :: ocpld
    real(NER), intent(in)                       :: p1, p2
    real(NER), dimension(1:2, 1:2), intent(out) :: vab

    integer :: L

    vab = 0.0_NER
    L = ocpld%lsj%L
    select case (L)
      case (0)
        vab(1, 1) = 1.0_NER
      case (1)
        vab(1, 1) = p1*p2
      case default
        vab(1, 1) = (p1*p2)**L
    end select
    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
      vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_VS0_cpld

  ! vab(1, 2) = reg(p1)*reg(p2)*(p1*p2)^L * (p2^2)
  ! vab(2, 1) = reg(p1)*reg(p2)*(p1*p2)^L * (p1^2)
  ! < p1 | VSoff | p2>
  subroutine Vfunc_VSoff_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)                 :: ocpld
    real(NER), intent(in)                       :: p1, p2
    real(NER), dimension(1:2, 1:2), intent(out) :: vab

    integer :: L
    real(NER)   :: pref

    vab = 0.0_NER
    L = ocpld%lsj%L
    vab(1, 2) = p2*p2
    vab(2, 1) = p1*p1
    pref = 0.0_NER
    select case (L)
      case (0)
        pref = 1.0_NER
      case (1)
        pref = p1*p2
      case default
        pref = (p1*p2)**L
    end select
    vab = vab*pref
    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
      vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_VSoff_cpld

end module mod_obj_cpld
