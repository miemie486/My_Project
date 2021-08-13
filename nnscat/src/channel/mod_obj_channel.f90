! Bingwei Long 06/30/2013
! Class channel
! * The default units are MeV and MeV^(-1), unless noted otherwise.

module mod_obj_channel

  use nneft_type
  use util_gauleg,    only: feb_tanabsc
  use eft_potspwd
  use nneft_loops
  use eft_phaseconv
  use util_io
  use util_cheb
  use nneft_phyconst, only: PC_delta_cutoff

  implicit none

  ! private ::  &
  !   & DEFAULT_N, DEFAULT_LAMBDA, DEFAULT_MAMBDA, DEFAULT_UPTOQN, DEFAULT_k

  ! integer, parameter      ::                &
  !   & PHASE_INCREMENTAL = 10,               &
  !   & PHASE_TOTAL       = 20

  ! integer, parameter      ::                &
  !   & MAX_NUMBER_POWFAC = 10,               &
  !   & MAX_NUMBER_ARRAY  = 10,               &
  !   & MAX_NUMBER_ORDER  = 10,               &
  !   & MAX_NUMBER_INPUTS = 40,               &
  !   & DEFAULT_N         = 100,              &
  !   & DEFAULT_REGTYPE   = REGTYPE_GAUSSIAN, &
  !   & DEFAULT_UPTOQN    = 0,                &
  !   & MAX_NUMBER_CS     = 30

  ! ! In unit of GeV
  ! real(NER), parameter    ::                &
  !   & DEFAULT_LAMBDA = 1.4_NER,             &
  !   & DEFAULT_MAMBDA = 0.8_NER,             &
  !   & BIG_LAMBDA     = 3.0_NER,             &
  !   & LOGMAP_LAMBDA  = 2.4_NER,             &
  !   & DEFAULT_k      = 0.01_NER

  ! ! In unit of MeV
  ! real(NER), parameter    ::                &
  !   & DEFAULT_LAMBDA = 1400.0_NER,          &
  !   & DEFAULT_MAMBDA = 800.0_NER,           &
  !   & BIG_LAMBDA     = 3000.0_NER,          &
  !   & LOGMAP_LAMBDA  = 2400.0_NER,          &
  !   & DEFAULT_k      = 10.0_NER

  type :: aux_real_array
    logical                             :: created = .false., inited = .false.
    real(NER), dimension(:), pointer    :: ptr => null()
  end type aux_real_array

  type VM_sngl
    logical :: created = .false., inited = .false.
    real(NER), dimension(:, :), pointer :: ptr => null()
  end type VM_sngl

  type :: obj_channel

    ! uptoQn is in the Long_Yang convention: LO is 0, NLO is 1, and so on
    integer             ::                            &
      & chnid   = CHNID_UNDEFINED,                    &
      & NCH     = 1,                                  &
      & N       = DEFAULT_N,                          &
      & Nln     = DEFAULT_N,                          &
      & uptoQn  = DEFAULT_UPTOQN,                     &
      & regtype = DEFAULT_REGTYPE,                    &
      & uflag(1:MAX_NUMBER_ORDER)

    real(NER)           ::                            &
      & k          = DEFAULT_k,                       &
      & Mambda     = DEFAULT_MAMBDA,                  &
      & Mambda_sqr = DEFAULT_MAMBDA*DEFAULT_MAMBDA,   &
      & Cs(1:MAX_NUMBER_CS),                          &
      & inparas(1:MAX_NUMBER_INPUTS),                 &
      & univio(1:MAX_NUMBER_ORDER)

    ! Ecm is independent of k
    complex(NEC)        :: Ecm

    type(lsj_symbol)    :: lsj

    character(len = 3)  :: chnstr

    logical             ::                            &
      & created          = .false.,                   &
      & mesh_inited      = .false.,                   &
      & k_defined        = .false.,                   &
      & TQn_updated      = .false.,                   &
      & phsn_updated     = .false.,                   &
      & inpar_set        = .false.,                   &
      & log_mapped       = .false.,                   &
      & unitarity_warning_enabled = .true.

    real(NER), allocatable ::                         &
      & msh(:),                                       &
      & wght(:),                                      &
      & lsqr(:),                                      &
      & half_regfac(:)

    ! See create_powfac_channel for comments on %powfac
    type(aux_real_array), dimension(1:MAX_NUMBER_POWFAC), private   :: powfac
    ! VM_poly0 = reg(p)*reg(q)*(p*q)^L
    type(VM_sngl)   :: VM_poly0

  contains

    procedure   :: init                         => init_channel
    procedure   :: check_warning_inited         => check_warning_inited_channel
    procedure   :: init_variables               => init_variables_channel
    procedure   :: init_mesh                    => init_mesh_channel
    procedure   :: update_k                     => update_k_channel
    procedure   :: release                      => release_channel
    procedure   :: erase                        => erase_channel

    procedure, private   :: create_powfac       => create_powfac_channel
    procedure, private   :: init_powfac         => init_powfac_channel
    procedure, private   :: update_powfac       => update_powfac_channel
    procedure   :: set_ptr_to_powfac            => set_ptr_to_powfac_channel

    procedure   :: create_VM_sngl               => create_VM_sngl_channel
    procedure   :: release_VM_sngl              => release_VM_sngl_channel
    procedure   :: add_2VMs_sngl                => add_2VMs_sngl_channel

    procedure   :: init_VMpoly0                 => init_VMpoly0_channel
    procedure   :: update_edge_VMpoly0          => update_edge_VMpoly0_channel

    procedure   :: half_regulator               => half_regulator_channel
    procedure   :: regulator                    => regulator_channel

    procedure   :: get_paras_filename           => get_paras_filename_channel
    procedure   :: get_writeout_string          => get_writeout_string_channel
    procedure   :: get_num_paras                => get_num_paras_channel
    procedure   :: read_inputfile               => read_inputfile_channel
    procedure   :: load_inputs                  => load_inputs_channel
    ! change_a_C is frozen
    ! procedure   :: change_a_C                   => change_a_C_channel

    procedure   :: loop_powpow                  => loop_powpow_channel
    procedure   :: loop_powpow_cmplxE           => loop_powpow_cmplxE_channel

    ! get_kcot and get_ere_paras are frozen
    ! procedure   :: get_kcot                     => get_kcot_channel
    ! procedure   :: get_ere_paras                => get_ere_paras_channel

    procedure   :: get_allTQn                   => get_allTQn_channel
    procedure   :: convert_TQn_to_phase_quiet             &
      & => convert_TQn_to_phase_quiet_channel
    procedure   :: output_phase_to_string                 &
      & => output_phase_to_string_channel
    procedure   :: output_total_phase                     &
      & => output_total_phase_channel
    procedure   :: output_incremental_phase               &
      & => output_incremental_phase_channel
    procedure   :: convert_TQn_to_phase_as_formatted_text &
      & => convert_TQn_to_phase_as_formatted_text_channel
    procedure   :: get_output_formatted_text_for_klist    &
      & => get_output_formatted_text_for_klist_channel
    procedure   :: get_phase_formatted_text_for_klist     &
      & => get_phase_formatted_text_for_klist_channel
    procedure   :: get_phase_for_klist                    &
      & => get_phase_for_klist_channel
    procedure   :: get_allorders_phase_for_klist          &
      & => get_allorders_phase_for_klist_channel
    procedure   :: GenerateUnitarityMsg                   &
      & => GenerateUnitarityMsg_channel

  end type obj_channel


contains


  ! Check and set chnid, chnstr, lsj, and uptoQn. If successful, set self%created
  ! to be true. The order of create_channel and init_channel is interchangeable.
  subroutine create_channel(self, chnid, uptoQn)

    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: chnid, uptoQn

    logical                             :: succeeded

    if (uptoQn < 0) then
      write (standard_error_unit, '(a)')  &
        & 'create_channel: uptoQn cannot be negative. Zeroth order is assumed.'
      self%uptoQn = 0
    else
      if (uptoQn > MAX_NUMBER_ORDER) then
        write (standard_error_unit, '(a, i2, a)')  &
          & 'create_channel: uptoQn cannot be larger than ', MAX_NUMBER_ORDER, &
          & '. Zeroth order is assumed.'
        self%uptoQn = 0
      else
        self%uptoQn = uptoQn
      end if
    end if
    call convert_chnid_to_lsj(chnid, self%lsj, succeeded)
    if (.not. succeeded) then
      write(standard_error_unit, '(a)')   &
        & 'create_channel: Conversion of chnid to lsj failed. l is assumed to be 0.'
      self%lsj%l = 0
      self%chnstr = 'NUL'
    else
      self%chnid = chnid
      call convert_lsj_to_text(self%lsj, self%chnstr)
    end if
    self%created = .true.

  end subroutine create_channel


  ! Call init_variables()
  ! Call init_mesh()
  ! The order of create_channel and init_channel is interchangeable.
  subroutine init_channel(self, N, regtype, Mambda)

    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    call self%init_variables(N, regtype, Mambda)
    call self%init_mesh()

  end subroutine init_channel

  ! Set N, regtype, Mambda, and Lambda.
  ! Alert if Mambda is larger than Lambda.
  subroutine init_variables_channel(self, N, regtype, Mambda)

    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    self%N               = N
    self%regtype         = regtype
    ! ! Converting from MeV to GeV
    ! self%Mambda          = Mambda*1.0E-3_NER
    ! self%Mambda_sqr      = Mambda*Mambda
    ! self%Lambda          = Lambda*1.0E-3_NER

    self%Mambda          = Mambda
    self%Mambda_sqr      = Mambda*Mambda
!     PC_delta_cutoff      = Mambda
    PC_delta_cutoff      = 700.0_NER
!     PC_delta_cutoff      = 5000.0_NER
!     PC_delta_cutoff      = 10000.0_NER

  end subroutine init_variables_channel

  ! Allocate memory for msh, wght, lsgr, and half_regfac
  subroutine init_mesh_channel(self)

    class(obj_channel), intent(inout)     :: self

    integer     :: ii, N

    N = self%N
    if (allocated(self%msh)) deallocate(self%msh)
    if (allocated(self%wght)) deallocate(self%wght)
    if (allocated(self%lsqr)) deallocate(self%lsqr)
    if (allocated(self%half_regfac)) deallocate(self%half_regfac)
    allocate(self%msh(1:N), self%wght(1:N), self%lsqr(1:N), self%half_regfac(1:N))
    ! When self%mesh_inited is set, mesh, wght, lsqr, and half_regfac must
    ! have been allocated, but their initial values have yet been set.

    ! TanMesh used, value of Lambda is set but has no impact
    call feb_tanabsc(self%N,                                                   &
      & self%Mambda*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)),         &
      & self%msh, self%wght)

    self%mesh_inited = .true.
    do ii = 1, self%N
      self%lsqr(ii)        = self%msh(ii)*self%msh(ii)
      self%half_regfac(ii) = self%half_regulator(self%msh(ii))
    end do

  end subroutine init_mesh_channel

  ! Frozen
  ! Allocate memory for msh, wght, lsgr, and half_regfac. msh has two
  ! sections. One is linear mapping and the other is logarithmic mapping.
  ! subroutine init_logmesh_channel(self)

  !     class(obj_channel), intent(inout)     :: self

  !     integer     :: ii, N, N1, N2

  !     N = self%N
  !     if (allocated(self%msh)) deallocate(self%msh)
  !     if (allocated(self%wght)) deallocate(self%wght)
  !     if (allocated(self%lsqr)) deallocate(self%lsqr)
  !     if (allocated(self%half_regfac)) deallocate(self%half_regfac)
  !     allocate(self%msh(1:N), self%wght(1:N), self%lsqr(1:N), self%half_regfac(1:N))
  !     ! When self%mesh_inited is set, mesh, wght, lsqr, and half_regfac must
  !     ! have been allocated, but their initial values have yet been set.
  !     if (self%Lambda < BIG_LAMBDA) then
  !         call set_mesh_gnr(self%N, 0.0_NER, self%Lambda, self%msh, self%wght)
  !     else
  !         N1 = int(N*0.7)
  !         N2 = N - N1
  !         call set_mesh_gnr(N1, 0.0_NER, LOGMAP_LAMBDA, self%msh, self%wght)
  !         call feb_logabsc_2pt(N2, self%msh(N1+1:N), self%wght(N1+1:N), LOGMAP_LAMBDA, self%Lambda)
  !     end if
  !     self%mesh_inited = .true.
  !     do ii = 1, self%N
  !         self%lsqr(ii)        = self%msh(ii)*self%msh(ii)
  !         self%half_regfac(ii) = self%half_regulator(self%msh(ii))
  !     end do

  ! end subroutine init_logmesh_channel

  ! Release arrays or other objects allocated during self's existence.
  ! Reset some of the attributes to its default value.
  subroutine release_channel(self)

    class(obj_channel), intent(inout)   :: self

    integer :: ii

    if (allocated(self%msh)) deallocate(self%msh)
    if (allocated(self%wght)) deallocate(self%wght)
    if (allocated(self%lsqr)) deallocate(self%lsqr)
    if (allocated(self%half_regfac)) deallocate(self%half_regfac)
    do ii = 1, MAX_NUMBER_POWFAC
      if (associated(self%powfac(ii)%ptr)) deallocate(self%powfac(ii)%ptr)
      self%powfac(ii)%inited  = .false.
      self%powfac(ii)%created = .false.
    end do
    self%N              = DEFAULT_N
    self%Nln            = DEFAULT_N
    self%regtype        = DEFAULT_REGTYPE
    self%k              = DEFAULT_k
    self%Mambda         = DEFAULT_MAMBDA
    self%Mambda_sqr     = DEFAULT_MAMBDA*DEFAULT_MAMBDA
    self%k_defined      = .false.
    self%TQn_updated    = .false.
    self%phsn_updated   = .false.
    self%mesh_inited    = .false.
    self%inpar_set = .false.
    self%log_mapped     = .false.
    self%unitarity_warning_enabled = .true.

  end subroutine release_channel


  ! Call %release()
  ! Reset %chnid, %uptoQn, and %created to their default values.
  subroutine erase_channel(self)

    class(obj_channel), intent(inout)   :: self

    ! If self still has things allocated, release them
    if (self%mesh_inited) call self%release()
    self%chnid   = CHNID_UNDEFINED
    self%uptoQn  = DEFAULT_UPTOQN
    self%created = .false.

  end subroutine erase_channel

  ! Check whether self is initialized.
  ! If not, show a warning message with invoking_source shown in the front
  subroutine check_warning_inited_channel(self, invoking_source)

    class(obj_channel), intent(in)    :: self
    character(len = *), intent(in)      :: invoking_source

    if (.not. self%mesh_inited) write (standard_error_unit, '(a, a)')          &
      & invoking_source, ': Mesh not initialized.'

  end subroutine check_warning_inited_channel

  ! Update the current k and change the powfac arrays accordingly if they have
  ! been initialized.

  subroutine update_k_channel(self, k)

    class(obj_channel), intent(inout)   :: self
    real(NER), intent(in)               :: k

    integer :: ii

    self%k = k
    self%k_defined = .true.
    do ii = 1, MAX_NUMBER_POWFAC
      if (self%powfac(ii)%inited) call self%update_powfac(ii)
    end do

  end subroutine update_k_channel


  ! powfac(1) = ppi**L*half_regfac
  ! powfac(2) = ppi**(L+2)*half_regfac
  ! powfac(3) = ppi**(L+4)*half_regfac
  ! ...

  ! The array powfac has been designed as such that children of obj_channel
  ! do not need to know how powfac is maintained.
  ! powfac is not initialized until it is cited through set_ptr_to_powfac.
  ! Therefore, set_ptr_to_powfac has to detect whether the required powfac is
  ! properly initialized or up-to-date with the current k.

  ! Note that update_k_channel will also update the powfac arrays when k is
  ! updated.

  ! Methods of obj_channel utilizing powfac will create/initialize powfac when
  ! powfac for the first time is used.
  subroutine create_powfac_channel(self, array_ID)

    class(obj_channel), intent(inout)       :: self
    integer, intent(in)                     :: array_ID

    if (.not. self%mesh_inited) then
      write (standard_error_unit, '(a)')      &
        & 'create_powfac_channel: Mesh not initialized'
      return
    end if
    if (associated(self%powfac(array_ID)%ptr)) then
      write (standard_error_unit, '(a)')      &
        & 'create_powfac_channel: Inconsistency warning 1'
      deallocate(self%powfac(array_ID)%ptr)
    end if
    allocate(self%powfac(array_ID)%ptr(1:self%N+1))
    self%powfac(array_ID)%created = .true.

  end subroutine create_powfac_channel

  subroutine init_powfac_channel(self, array_ID)

    class(obj_channel), intent(inout)       :: self
    integer, intent(in)                     :: array_ID

    integer     :: ii, power
    real(NER)   :: ppi

    if (.not. self%mesh_inited) then
      write (standard_error_unit, '(a)')      &
        & 'init_powfac_channel: Mesh not initialized'
      return
    end if
    if (self%powfac(array_ID)%inited) return
    if (.not. self%powfac(array_ID)%created)  &
      & call self%create_powfac(array_ID)
    power = self%lsj%l + (array_ID-1)*2
    do ii = 1, self%N
      ppi = self%msh(ii)
      self%powfac(array_ID)%ptr(ii) = self%half_regfac(ii) * ppi**power
    end do
    self%powfac(array_ID)%inited = .true.

  end subroutine init_powfac_channel

  ! Update powfac for current k
  subroutine update_powfac_channel(self, array_ID)

    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: array_ID

    integer :: power

    if (.not. self%mesh_inited) then
      write (standard_error_unit, '(a)')      &
        & 'update_powfac_channel: Mesh not initialized'
      return
    end if
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')      &
        & 'update_powfac_channel: k is not defined'
      return
    end if
    if (.not. self%powfac(array_ID)%inited)   &
      & call self%init_powfac(array_ID)
    power = self%lsj%l + (array_ID-1)*2
    self%powfac(array_ID)%ptr(self%N+1) =     &
      & self%half_regulator(self%k) * self%k**power

  end subroutine update_powfac_channel

  subroutine set_ptr_to_powfac_channel(self, array_ID, ptr_out)

    class(obj_channel), intent(inout)               :: self
    integer, intent(in)                             :: array_ID
    real(NER), dimension(:), pointer, intent(out)   :: ptr_out

    if (.not. self%powfac(array_ID)%inited) then
      call self%init_powfac(array_ID)
      if (self%k_defined) call self%update_powfac(array_ID)
    end if
    ptr_out => self%powfac(array_ID)%ptr

  end subroutine set_ptr_to_powfac_channel


  ! Methods to manipulate VM_sngl
  ! created: whether the memory of VM has been allocated.
  ! inited: whether VM has been initialized with certain potenital

  ! create_VM_sngl_channel will initialize all elements to zero
  subroutine create_VM_sngl_channel(self, VM)

    class(obj_channel), intent(inout)   :: self
    type(VM_sngl), intent(inout)        :: VM

    call self%check_warning_inited('create_VM_sngl_channel')
    if (VM%created) write (standard_error_unit, '(a)')      &
      & 'create_VM_sngl_channel: Warning! Trying to create a VM that has already been created'
    if (associated(VM%ptr)) deallocate(VM%ptr)
    allocate(VM%ptr(1:self%N+1, 1:self%N+1))
    VM%ptr(1:self%N+1, 1:self%N+1) = 0.0_NER
    VM%created = .true.
    VM%inited  = .false.

  end subroutine create_VM_sngl_channel

  subroutine release_VM_sngl_channel(self, VM)

    class(obj_channel), intent(inout)   :: self
    type(VM_sngl), intent(inout)        :: VM

    if (.not. associated(VM%ptr)) then
      write (standard_error_unit, '(a)')  &
        & 'release_VM_sngl_channel: Warning: No need to release.'
    else
      deallocate(VM%ptr)
    endif
    VM%created = .false.
    VM%inited  = .false.

  end subroutine release_VM_sngl_channel

  ! Add VM1 and VM2 up to fill VMdest. If k_defined is true, also act on the
  ! edges of the VMs.
  subroutine add_2VMs_sngl_channel(self, VM1, VM2, VMdest)

    class(obj_channel), intent(inout)  :: self
    type(VM_sngl), intent(inout)       :: VM1, VM2, VMdest

    integer :: ii, jj

    call self%check_warning_inited('add_2VMs_sngl_channel')
    do ii = 1, self%N
      do jj = 1, self%N
        VMdest%ptr(ii, jj) = VM1%ptr(ii, jj) + VM2%ptr(ii, jj)
      end do
    end do
    if (self%k_defined) then
      do ii = 1, self%N
        VMdest%ptr(ii, self%N+1) = VM1%ptr(ii, self%N+1) + VM2%ptr(ii, self%N+1)
        VMdest%ptr(self%N+1, ii) = VM1%ptr(self%N+1, ii) + VM2%ptr(self%N+1, ii)
      end do
      VMdest%ptr(self%N+1, self%N+1) = VM1%ptr(self%N+1, self%N+1) +           &
        & VM2%ptr(self%N+1, self%N+1)
    end if
    VMdest%inited = .true.

  end subroutine add_2VMs_sngl_channel

  ! (N+1)X(N+1) VMpoly0(i, j) = powfac(i; 1) * powfac(j; 1)
  subroutine init_VMpoly0_channel(self, VMpoly0)

    class(obj_channel), intent(inout)   :: self
    type(VM_sngl), intent(inout)        :: VMpoly0

    integer :: ii, jj
    real(NER), dimension(:), pointer    :: ptr

    if (.not. VMpoly0%created) then
      write (standard_error_unit, '(a)')  &
        'init_VMpoly0_channel: VM not created. Nothing will be done.'
      return
    end if
    if (.not. self%mesh_inited) then
      write (standard_error_unit, '(a)')  &
        'init_VMpoly0_channel: Mesh not initialized. Nothing will be done.'
      return
    end if
    call self%set_ptr_to_powfac(1, ptr)
    do ii = 1, self%N
      do jj = 1, ii - 1
        VMpoly0%ptr(ii, jj) = ptr(ii) * ptr(jj)
        VMpoly0%ptr(jj, ii) = VMpoly0%ptr(ii, jj)
      end do
      VMpoly0%ptr(ii, ii) = ptr(ii)*ptr(ii)
    end do
    VMpoly0%inited = .true.

  end subroutine init_VMpoly0_channel

  subroutine update_edge_VMpoly0_channel(self, VMpoly0)

    class(obj_channel), intent(inout)   :: self
    type(VM_sngl), intent(inout)        :: VMpoly0

    integer :: ii, locN
    real(NER), dimension(:), pointer :: ptr

    call self%set_ptr_to_powfac(1, ptr)
    locN = self%N
    do ii = 1, locN
      VMpoly0%ptr(ii, locN+1) = ptr(ii)*ptr(locN+1)
      VMpoly0%ptr(locN+1, ii) = VMpoly0%ptr(ii, locN+1)
    end do
    VMpoly0%ptr(locN+1, locN+1) = ptr(locN+1)*ptr(locN+1)

  end subroutine update_edge_VMpoly0_channel

  subroutine get_paras_filename_channel(self, fname)

    class(obj_channel), intent(inout)   :: self
    character(len = *), intent(out)     :: fname

    fname = 'inputs_'//self%chnstr//'.in'

  end subroutine get_paras_filename_channel

  ! A characteristic string that will be used as a part of output file names
  ! Its maximum length is 10.
  ! (To be used by get_output_for_Lambda_list in drv_pwphase)
  subroutine get_writeout_string_channel(self, outstr)

    class(obj_channel), intent(inout)   :: self
    character(len = 10), intent(out)    :: outstr

    outstr = self%chnstr

  end subroutine get_writeout_string_channel

  subroutine get_num_paras_channel(self, num_cnttrm)

    class(obj_channel), intent(inout)   :: self
    integer, intent(out)                :: num_cnttrm

    num_cnttrm = 0

  end subroutine get_num_paras_channel

  ! Format of .in :
  ! Mambda <separator> para(1) <separator> para(2) ...
  subroutine read_inputfile_channel(self, funt, num_cnttrms, nLmbds,           &
    & Lambda_cnttrms, succeeded)

    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: funt
    integer, intent(out)                :: nLmbds, num_cnttrms
    real(NER), dimension(:, :), intent(out) :: Lambda_cnttrms
    logical, intent(out)                    :: succeeded

    integer                 :: flag
    character(len = 512)    :: fname

    call self%get_num_paras(num_cnttrms)
    call self%get_paras_filename(fname)
    call read_ncols_io(funt, trim(fname), num_cnttrms+1, Lambda_cnttrms, nLmbds, flag)
    if (flag == read_ncols_FAIL) then
      write (standard_error_unit, '(a)')  &
        'read_inputfile_channel: couldn''t open or read '//trim(fname)
      succeeded = .false.
    else
      succeeded = .true.
    end if

  end subroutine read_inputfile_channel

  ! Load inputs onto self%inparas. And, by default, take all input
  ! parameters as C's. load_inputs is expected to be overriden.
  subroutine load_inputs_channel(self, num_inputs, inputs)

    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: num_inputs
    real(NER), dimension(:), intent(in) :: inputs

    self%inparas(1:num_inputs) = inputs(1:num_inputs)
    self%Cs(1:num_inputs) = self%inparas(1:num_inputs)
    self%inpar_set = .true.

  end subroutine load_inputs_channel

  ! Frozen
  ! subroutine change_a_C_channel(self, n, Cn)

  !     class(obj_channel), intent(inout)   :: self
  !     integer, intent(in)                 :: n
  !     real(NER), intent(in)               :: Cn

  !     self%Cs(n) = Cn

  ! end subroutine change_a_C_channel

  function half_regulator_channel(self, p)

    real(NER)                       :: half_regulator_channel
    class(obj_channel), intent(in)  :: self
    real(NER), intent(in)           :: p

    call self%check_warning_inited('half_regulator_channel')
    half_regulator_channel = regltr_half_PP(self%regtype, p, self%Mambda)

  end function half_regulator_channel

  function regulator_channel(self, p1, p2)

    real(NER)                       :: regulator_channel
    class(obj_channel), intent(in)  :: self
    real(NER), intent(in)           :: p1, p2

    call self%check_warning_inited('regulator_channel')
    regulator_channel = regltr_PP(self%regtype, p1, p2, self%Mambda)

  end function regulator_channel

  ! \int dl l^2/(k^2 - l^2)*powfac(n)*powfac(m)
  function loop_powpow_channel(self, n_powID, m_powID)

    complex(NEC)                        :: loop_powpow_channel
    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: n_powID, m_powID

    real(NER), dimension(:), pointer    :: ptr_n, ptr_m

    call self%check_warning_inited('loop_powpow_channel')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)') 'loop_powpow_channel: k not defined.'
      return
    end if
    call self%set_ptr_to_powfac(n_powID, ptr_n)
    call self%set_ptr_to_powfac(m_powID, ptr_m)
    loop_powpow_channel = loop_1loop_VkVk_sngl(self%k, self%N, self%msh,       &
      & self%wght, self%lsqr, ptr_n, ptr_m)

  end function loop_powpow_channel

  function loop_powpow_cmplxE_channel(self, n_powID, m_powID, mE)

    complex(NEC)                        :: loop_powpow_cmplxE_channel
    class(obj_channel), intent(inout)   :: self
    integer, intent(in)                 :: n_powID, m_powID
    complex(NEC), intent(in)            :: mE

    real(NER), dimension(:), pointer    :: ptr_n, ptr_m

    call self%check_warning_inited('loop_powpow_cmplxE_channel')
    call self%set_ptr_to_powfac(n_powID, ptr_n)
    call self%set_ptr_to_powfac(m_powID, ptr_m)
    loop_powpow_cmplxE_channel = loop_1loop_VkVk_cmplxE_sngl(mE, self%N,       &
      & self%msh, self%wght, self%lsqr, ptr_n, ptr_m)

  end function loop_powpow_cmplxE_channel

  ! Checks several prerequisites before overriden method calculates real
  ! amplitudes.
  subroutine get_allTQn_channel(self)

    class(obj_channel), intent(inout)   :: self

    call self%check_warning_inited('get_allTQn_channel')
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        & 'get_allTQn_channel: Warning: k is not defined.'
    end if

  end subroutine get_allTQn_channel


  ! Frozen
  ! ! A place holder for methods that calculate the k-matrix element for a given
  ! ! k. The method will be called by get_ere_paras_channel to calculate the
  ! ! effective range expansion parameters.
  ! subroutine get_kcot_channel(self, k, kcot)

  !     class(obj_channel), intent(inout)   :: self
  !     real(NER), intent(in)               :: k
  !     real(NER), intent(out)              :: kcot

  !     call self%update_k(k)
  !     kcot = 0.0

  ! end subroutine get_kcot_channel


  ! Frozen
  ! ! Taylor expands the k-matrix element provided by get_kcot_channel.
  ! ! unit_k is the value of 1 unit of CM momentum. num_paras is the order of
  ! ! sought effective range expansion with 1 meaning that only the scattering
  ! ! is desired.
  ! subroutine get_ere_paras_channel(self, unit_k, num_paras, ereparas)

  !     class(obj_channel), intent(inout)   :: self
  !     real(NER), intent(in)               :: unit_k
  !     integer, intent(in)                 :: num_paras
  !     real(NER), dimension(:), intent(out):: ereparas

  !     real(NER), parameter    :: CHEB_A = 0.01_NER, CHEB_B = 0.1_NER
  !     integer, parameter      :: CHEB_BASE = 4
  !     integer     :: ii
  !     real(NER)   :: k, kcot
  !     real(NER), dimension(1:num_paras+CHEB_BASE+1) :: tmp_paras

  !     call cheb_taylor(CHEB_A, CHEB_B, num_paras + CHEB_BASE, tmp_paras, func_km)
  !     do ii = 1, num_paras
  !         !ereparas(ii) = tmp_paras(ii)*(PC_hbarc)**(2*ii-1)
  !         ereparas(ii) = tmp_paras(ii)/(unit_k)**(2*ii-2)
  !     end do

  ! contains

  !     function func_km(x)

  !         real(NER), intent(in)   :: x
  !         real(NER)               :: func_km

  !         real(NER)   :: y

  !         if (x < 0.0_NER) then
  !             write (standard_error_unit, '(a)') 'get_ere_paras_channel: k must be positive.'
  !             y = 0.01
  !         else
  !             y = x
  !         end if
  !         k = sqrt(y)*unit_k
  !         call self%get_kcot(k, kcot)
  !         func_km = kcot

  !     end function func_km

  ! end subroutine get_ere_paras_channel

  ! Display warning message in case unitarity broken
  subroutine convert_TQn_to_phase_quiet_channel(self)

    class(obj_channel), intent(inout)   :: self

    character(len = 1024) :: msg
    logical               :: is_broken

!debug
!       print *, self%inparas(1)
    if (self%unitarity_warning_enabled) then
      call self%GenerateUnitarityMsg(is_broken, msg)
      if (is_broken) then
        write(standard_error_unit, '(e10.4e1, x, e10.4e1, x, a3, x, a)')       &
          & self%k, self%Mambda, self%chnstr, trim(msg)
      end if
    end if

  end subroutine convert_TQn_to_phase_quiet_channel

  ! Check phsn_updated
  subroutine output_phase_to_string_channel(self, outstr)

    class(obj_channel), intent(inout)   :: self
    character(len = *), intent(out)     :: outstr

    if (.not. self%phsn_updated) then
      write (standard_error_unit, '(a)')  &
        'output_phase_to_string_channel: Warning: phsn is not updated.'
      return
    end if

  end subroutine output_phase_to_string_channel

  ! Check phsn_updated
  ! The derived class should return the sum of phsn(:)
  subroutine output_total_phase_channel(self, out_phs)

    class(obj_channel), intent(inout)   :: self
    class(*), intent(inout)             :: out_phs

    if (.not. self%phsn_updated) then
      write (standard_error_unit, '(a)')  &
        'output_total_phase_channel: Warning: phsn is not updated.'
      return
    end if

  end subroutine output_total_phase_channel

  subroutine output_incremental_phase_channel(self, out_phs)

    class(obj_channel), intent(inout)   :: self
    class(*), intent(inout)             :: out_phs(:)

    if (.not. self%phsn_updated) then
      write (standard_error_unit, '(a)')  &
        'output_incremental_phase_channel: Warning: phsn is not updated.'
      return
    end if

  end subroutine output_incremental_phase_channel

  ! Generate a warning message, according to unitarity check flags given by
  ! self%uflag(:).
  subroutine GenerateUnitarityMsg_channel(self, is_broken, msg)

    class(obj_channel), intent(inout)   :: self
    logical, intent(out)                :: is_broken
    character(len = *), intent(out)     :: msg

    integer                 :: ii
    character(len = 128)    :: minibuff

    minibuff = ''
    msg = '### uni. vio. '
    is_broken = .false.
    do ii = 1, self%uptoQn+1
      if (self%uflag(ii) == UNITARITY_BROKEN) then
        is_broken = .true.
        write (minibuff, '(1x, a, i1, x, e10.4e1)')  'Q^', ii-1, self%univio(ii)
        msg = trim(msg)//trim(minibuff)
      end if
    end do

  end subroutine GenerateUnitarityMsg_channel

  subroutine convert_TQn_to_phase_as_formatted_text_channel(self, outstr)

    class(obj_channel), intent(inout)   :: self
    character(len = *), intent(out)     :: outstr

    character(len = 512)    :: tmpstr, longbuff
    logical                 :: redflag

    call self%convert_TQn_to_phase_quiet()
    call self%output_phase_to_string(tmpstr)
    call self%GenerateUnitarityMsg(redflag, longbuff)
    if (redflag) then
      outstr = trim(tmpstr)//trim(longbuff)
    else
      outstr = trim(tmpstr)
    end if

  end subroutine convert_TQn_to_phase_as_formatted_text_channel

  ! If ouput_option is OUTOPT_PHASE call get_phase_formatted_text_for_klist_channel
  ! If other output is wanted, this method needs to be overridden
  subroutine get_output_formatted_text_for_klist_channel(self, N, regtype,     &
    & Mambda, kk, numk, inputs, num_inputs, output_option, output_text)

    class(obj_channel), intent(inout)       :: self
    integer, intent(in)                     :: regtype, N, numk, num_inputs
    real(NER), intent(in)                   :: Mambda
    real(NER), dimension(1:), intent(in)    :: kk, inputs
    integer, intent(in)                     :: output_option
    character(len=*), dimension(1:), intent(inout) :: output_text

    if (output_option == OUTOPT_PHASE) then
      call self%get_phase_formatted_text_for_klist(N, regtype, Mambda, kk,     &
        & numk, inputs, num_inputs, output_text)
    else
      return
    end if

  end subroutine get_output_formatted_text_for_klist_channel


  ! self must have been created but not initialized.
  subroutine get_phase_formatted_text_for_klist_channel(self, N, regtype,      &
    & Mambda, kk, numk, inputs, num_inputs, phase_out_text)

    class(obj_channel), intent(inout)       :: self
    integer, intent(in)                     :: regtype, N, numk, num_inputs
    real(NER), intent(in)                   :: Mambda
    real(NER), dimension(1:), intent(in)    :: kk, inputs
    character(len=*), dimension(1:), intent(inout) :: phase_out_text

    integer :: ii

    call self%load_inputs(num_inputs, inputs)
    call self%init(N, regtype, Mambda)
    do ii = 1, numk
      call self%update_k(kk(ii))
      call self%get_allTQn()
      call self%convert_TQn_to_phase_as_formatted_text(phase_out_text(ii))
    end do
    call self%release()

  end subroutine get_phase_formatted_text_for_klist_channel

  ! Note the order of methods called before the amplitudes can be harvested.
  subroutine get_phase_for_klist_channel(self, N, regtype, Mambda, kk, numk,   &
    & inputs, num_inputs, out_phase_lst)

    class(obj_channel), intent(inout)       :: self
    integer, intent(in)                     :: regtype, N, numk, num_inputs
    real(NER), intent(in)                   :: Mambda
    real(NER), intent(in)                   :: kk(:), inputs(:)
    !!! Attn: "class(*), intent(out)" may cause problem !!!
    class(*), intent(inout)                 :: out_phase_lst(:)

    integer :: ii
    ! integer, dimension(1:MAX_NUMBER_ORDER)  :: uflag

    call self%load_inputs(num_inputs, inputs)
    call self%init(N, regtype, Mambda)
    do ii = 1, numk
      call self%update_k(kk(ii))
      call self%get_allTQn()
      call self%convert_TQn_to_phase_quiet()
      call self%output_total_phase(out_phase_lst(ii))
    end do
    call self%release()

  end subroutine get_phase_for_klist_channel

  ! output_phase_lst(ii, n) is the phase shift correction for uptoQn=(n-1) at
  ! kk(ii)
  ! Eg., output_phase_lst(ii, 1) is the LO, output_phase_lst(ii, 2) is the
  ! first CORRECTION, and so on
  subroutine get_allorders_phase_for_klist_channel(self, N, regtype, Mambda,   &
      & kk, numk, inputs, num_inputs, out_phase_lst)

    class(obj_channel), intent(inout)       :: self
    integer, intent(in)                     :: regtype, N, numk, num_inputs
    real(NER), intent(in)                   :: Mambda
    real(NER), intent(in)                   :: kk(:), inputs(:)
    !!! Attn: "class(*), intent(out)" may cause problem !!!
    class(*), intent(inout)                 :: out_phase_lst(:, :)

    integer :: ii

    call self%load_inputs(num_inputs, inputs)
    call self%init(N, regtype, Mambda)
    do ii = 1, numk
      call self%update_k(kk(ii))
      call self%get_allTQn()
      call self%convert_TQn_to_phase_quiet()
      call self%output_incremental_phase(out_phase_lst(ii, 1:self%uptoQn+1))
    end do
    call self%release()

  end subroutine get_allorders_phase_for_klist_channel

end module mod_obj_channel
