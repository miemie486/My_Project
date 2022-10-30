
! Bingwei Long & Rui Peng  07/29/2019
! - Optimized the initialization of VMs, especially for strict perturbative
!   theory
! Bingwei Long  03/17/2019
! - switch from Vfunc_smplchn to gnrc_eft_pot

! Bingwei Long  10/20/2018

! Class: SiMPLified CHaNnel
! * This class has simpler structure than obj_channel so that it is easy to
!   easy to learn
! * The default units are MeV and MeV^(-1), unless noted otherwise.

module mod_obj_smplchn

  use nneft_type
  use util_gauleg,    only: feb_tanabsc
  use eft_potspwd
  use potwrap, only: gnrc_eft_pot, TPE0_epwrap, VPNQ0_epwrap, VPNQ2_diag_epwrap, &
      & VPNQ2_off_epwrap, VN2LO_MMWLY_epwrap, VLO_withc0_epwrap, zero_V, VPNQ5_epwrap, &
      & VNLO_MMWLY_epwrap, VN2LO_MMWLY_SHRG_epwrap
  use testwrap, only: reVPNQ2_off_epwrap, reVPNQ2_diag_epwrap, reVPNQ0_epwrap, reTPE0_epwrap, &
      & reVPNQ0_SCL_epwrap, reVPNQ2_SCL_diag_epwrap, reVPNQ5_epwrap
  use nneft_lsesv
  use eft_phaseconv
  use util_io
  use util_cheb
  use nneft_phyconst, only:PC_delta_cutoff, PC_mN
  use eft_h2wfs
  use util_rootfinder, only:zbrent
  implicit none

  integer, parameter  ::                              &
    & NPAR_SPR1S0_SMPL(1:4)      = (/2, 5, 9, 14/),   &
    & NPAR_LBWDIB_SMPL(1:4)      = (/2, 5, 9, 14/),   &
    & NPAR_WITHC0_SNGL_SMPL(1:4) = (/1, 1, 3, 5/),    &
    & NPAR_WOC0_SNGL_SMPL(1:4)   = (/0, 0, 1, 2/),    &
    & NPAR_WITHC0_CPLD_SMPL(1:4) = (/1, 1, 4, 7/),    &
    & NPAR_LY1S0_SMPL(1:3)       = (/1, 3, 6/)

    logical             ::                            &
!     pray allpert_BDE initial flag
      & allpert_BDE_initial = .false.,                &
      & pert_SCL_initial(2) = .false.

    !pray save 3s1 pert_binding energy
    real(NER)  ::                         &
      & allpert_BDE(4) = 0.0_NER,         &
      & NLOpert_SCL(2) = 0.0_NER,         &
      & N2LOpert_SCL(5) = 0.0_NER

    real(NER),parameter     ::            &
      & empirical_scatlength = -23.740_NER

  type :: obj_smplchn

    ! uptoQn is in the Long_Yang convention: LO is 0, NLO is 1, and so on
    integer             ::                            &
      & chnid   = CHNID_UNDEFINED,                    &
      & NCH     = 1,                                  &
      & N       = DEFAULT_N,                          &
      & totN    = DEFAULT_N,                          &
      & uptoQn  = DEFAULT_UPTOQN,                     &
      & regtype = DEFAULT_REGTYPE,                    &
      & npar(1:MAX_NUMBER_ORDER) = 0,                 &
      & uflag(1:MAX_NUMBER_ORDER) = UNITARITY_SAFE

    real(NER)           ::                            &
      & k          = DEFAULT_k,                       &
      & Mambda     = DEFAULT_MAMBDA,                  &
      & input(1:MAX_NUMBER_INPUTS),                   &
      & potpara(1:MAX_NUMBER_CS),                     &
      & miscpara(1:MAX_NUMBER_INPUTS),                &
      & univio(1:MAX_NUMBER_ORDER)

    ! Ecm is independent of k
    complex(NEC)        :: Ecm
    type(lsj_symbol)    :: lsj
    character(len = 3)  :: chnstr

    logical             ::                            &
      & created          = .false.,                   &
      & k_defined        = .false.,                   &
      & TQn_updated      = .false.,                   &
      & phsn_updated     = .false.,                   &
      & unitarity_warning_enabled = .true.,           &
      & VMinited(0:8)    = .false.,                   &
      & VLMinited(0:8)   = .false.,                   &
      ! vanishing(i) = .true. => V_{Q^i} = 0
      & vanishing(0:8)   = .false.,                   &
      & vlvs_separated(0:8) = .false.


    real(NER), allocatable ::                         &
      & msh(:),                                       &
      & wght(:),                                      &
      & lsqr(:)

    integer   ::                          &
      & V1_cheb_base_n = 4,               &
      & V2_cheb_base_n = 4,               &
      & V3_cheb_base_n = 4,               &
      & V4_cheb_base_n = 4,               &
      & V2_cheb_pertBD_n = 4
    real(NER) ::                          &
      & V1_cheb_bnd = 1.0E-3_NER,         &
      & V2_cheb_bnd = 1.0E-3_NER,         &
      & V3_cheb_bnd = 1.0E-3_NER,         &
      & V4_cheb_bnd = 1.0E-6_NER,         &
      & V2_cheb_pertBD_bnd = 1.0E-5_NER

    ! VMLO(1:NCH*(1+N), 1:NCH*(1+N))
    real(NER), allocatable ::             &
      & VMLO(:, :),                       &
      & VM1(:, :),                        &
      & VM2(:, :),                        &
      & VM3(:, :),                        &
      & VM4(:, :),                        &
      & VLM2(:, :),                       &
      & VLM3(:, :),                       &
      & VLM4(:, :),                       &
      & VMtmp(:, :)

    ! Tk(1:NCH*(1+N), 1:NCH)
    complex(NEC), allocatable :: Tk(:, :)
    complex(NEC)  :: TQn(1:MAX_NCH, 1:MAX_NCH, 1:MAX_NUMBER_ORDER)
    real(NER)     :: phsn(1:MAX_NUMBER_ORDER)
    type(triplet_phs_cpld)    :: trpln(1:MAX_NUMBER_ORDER)

    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_V0 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_V1 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_V2 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_V3 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_V4 => null()

    ! procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VL0 => null()
    ! procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VL1 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VL2 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VL3 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VL4 => null()

    ! procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VS0 => null()
    ! procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VS1 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VS2 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VS3 => null()
    procedure(gnrc_eft_pot), pointer, nopass  :: pVfunc_VS4 => null()

  contains

    procedure   :: init                         => init_smplchn
    procedure   :: init_variables               => init_variables_smplchn
    procedure   :: init_mesh                    => init_mesh_smplchn
    procedure   :: init_pots                    => init_pots_smplchn
    procedure   :: update_k                     => update_k_smplchn
    procedure   :: release                      => release_smplchn
    procedure   :: reset_VLMinited              => reset_VLMinited_smplchn
    procedure   :: erase                        => erase_smplchn
    procedure   :: fit_erase                    => fit_erase_smplchn

    procedure   :: get_num_paras                => get_num_paras_smplchn
    procedure   :: get_paras_filename           => get_paras_filename_smplchn
    procedure   :: read_inputfile               => read_inputfile_smplchn
    procedure   :: load_inputs                  => load_inputs_smplchn
    procedure   :: half_regulator               => half_regulator_smplchn
    procedure   :: regulator                    => regulator_smplchn

    procedure   :: init_VM                      => init_VM_smplchn
    procedure   :: init_VM_vlvs                 => init_VM_vlvs_smplchn
    procedure   :: update_edge_VM               => update_edge_VM_smplchn
    procedure   :: get_Tk_from_VM               => get_Tk_from_VM_smplchn
    procedure   :: get_allTQn                   => get_allTQn_smplchn
    procedure   :: GetSinglePertV               => GetSinglePertV_smplchn
    procedure   :: GetDoublePertV               => GetDoublePertV_smplchn

    procedure   :: convert_TQn_to_phase         => convert_TQn_to_phase_smplchn
    procedure   :: convert_phsn_to_text         => convert_phsn_to_text_smplchn
    procedure   :: get_incremental_phsn_value   &
      &              => get_incremental_phsn_value_smplchn
    procedure   :: GenerateUnitarityMsg         => GenerateUnitarityMsg_smplchn

  end type obj_smplchn

  private :: loc_series

contains

  ! Check and set chnid, chnstr, lsj, and uptoQn. If successful, set self%created
  ! to be true. The order of create_smplchn and init_smplchn is interchangeable.
  subroutine create_smplchn(self, chnid, uptoQn)

    class(obj_smplchn), intent(inout)   :: self
    integer, intent(in)                 :: chnid, uptoQn

    logical                             :: succeeded

    if (uptoQn < 0) then
      write (standard_error_unit, '(a)')  &
        & 'create_smplchn: uptoQn cannot be negative. Zeroth order is assumed.'
      self%uptoQn = 0
    else
      if (uptoQn > MAX_NUMBER_ORDER) then
        write (standard_error_unit, '(a, i2, a)')  &
          & 'create_smplchn: uptoQn cannot be larger than ', MAX_NUMBER_ORDER, &
          & '. Zeroth order is assumed.'
        self%uptoQn = 0
      else
        self%uptoQn = uptoQn
      end if
    end if
    call convert_chnid_to_lsj(chnid, self%lsj, succeeded)
    if (.not. succeeded) then
      write(standard_error_unit, '(a)')   &
        & 'create_smplchn: Conversion of chnid to lsj failed. l is assumed to be 0.'
      self%lsj%l = 0
      self%chnstr = 'NUL'
    else
      self%chnid = chnid
      call convert_lsj_to_text(self%lsj, self%chnstr)
    end if
    select case (type_of_channel_PP(chnid))
      case (CHTYPE_SNGL)
        self%NCH = 1
      case (CHTYPE_CPLD)
        self%NCH = 2
      case default
        write(standard_error_unit, '(a, i3, a)')   &
          & 'create_smplchn: ', chnid, ' is neither uncoupled or coupled channel'
        self%NCH = 1
    end select

    self%created = .true.

  end subroutine create_smplchn

  ! Call init_variables()
  ! Call init_mesh()
  ! The order of create_smplchn and init_smplchn is *not* interchangeable.
  subroutine init_smplchn(self, N, regtype, Mambda)

    class(obj_smplchn), intent(inout)   :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    if (.not. self%created) then
      write(standard_error_unit, '(a, i3, a)')   &
        & 'init_smplchn: WARNING! channel has not yet been created'
    end if
    call self%init_variables(N, regtype, Mambda)
    call self%init_mesh()

    if (allocated(self%VMtmp)) deallocate(self%VMtmp)
    allocate(self%VMtmp(1:self%totN, 1:self%totN))

    if (allocated(self%Tk)) deallocate(self%Tk)
    allocate(self%Tk(1:self%totN, 1:self%NCH))

  end subroutine init_smplchn

  ! Set N, regtype, Mambda, and Lambda.
  ! Alert if Mambda is larger than Lambda.
  subroutine init_variables_smplchn(self, N, regtype, Mambda)

    class(obj_smplchn), intent(inout)   :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    self%N               = N
    self%totN            = self%NCH*(self%N+1)
    self%regtype         = regtype
    self%Mambda          = Mambda
!     PC_delta_cutoff      = Mambda
!     PC_delta_cutoff      = 700.0_NER
!     PC_delta_cutoff      = 5000.0_NER
    PC_delta_cutoff      = 10000.0_NER
    ! Quick and dirty


  end subroutine init_variables_smplchn

  ! Allocate memory for msh, wght, and lsgr

  subroutine init_mesh_smplchn(self)

    class(obj_smplchn), intent(inout)     :: self

    integer     :: ii, N

    N = self%N
    if (allocated(self%msh)) deallocate(self%msh)
    if (allocated(self%wght)) deallocate(self%wght)
    if (allocated(self%lsqr)) deallocate(self%lsqr)
    allocate(self%msh(1:N), self%wght(1:N), self%lsqr(1:N))
    ! TanMesh used
    call feb_tanabsc(self%N,                                                   &
      & self%Mambda*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)),         &
      & self%msh, self%wght)
    do ii = 1, self%N
      self%lsqr(ii)        = self%msh(ii)*self%msh(ii)
    end do

  end subroutine init_mesh_smplchn

  ! Allocate memory for self%VMLO(1, 2, 3) and fill them out with
  ! If necessary, call self%pVfunc_*

  subroutine init_pots_smplchn(self, lbl_V)

    class(obj_smplchn), intent(inout)     :: self
    integer, intent(in) :: lbl_V

    if (self%VMinited(lbl_V)) return
    select case (lbl_V)
      case (0)
        if (allocated(self%VMLO)) deallocate(self%VMLO)
        allocate(self%VMLO(1:self%totN, 1:self%totN))
        call self%init_VM(self%pVfunc_V0, self%VMLO)
      case (1)
        if (allocated(self%VM1)) deallocate(self%VM1)
        allocate(self%VM1(1:self%totN, 1:self%totN))
        call self%init_VM(self%pVfunc_V1, self%VM1)
      case (2)
        if (allocated(self%VM2)) deallocate(self%VM2)
        allocate(self%VM2(1:self%totN, 1:self%totN))
        if (self%vlvs_separated(2)) then
          if (.not. self%VLMinited(2)) then
            if (allocated(self%VLM2)) deallocate(self%VLM2)
            allocate(self%VLM2(1:self%totN, 1:self%totN))
            call self%init_VM(self%pVfunc_VL2, self%VLM2)
            self%VLMinited(2) = .true.
          end if
          call self%init_VM_vlvs(self%pVfunc_VS2, self%VLM2, self%VM2)
        else
          call self%init_VM(self%pVfunc_V2, self%VM2)
        end if
      case (3)
        if (allocated(self%VM3)) deallocate(self%VM3)
        allocate(self%VM3(1:self%totN, 1:self%totN))
        if (self%vlvs_separated(3)) then
          if (.not. self%VLMinited(3)) then
            if (allocated(self%VLM3)) deallocate(self%VLM3)
            allocate(self%VLM3(1:self%totN, 1:self%totN))
            call self%init_VM(self%pVfunc_VL3, self%VLM3)
            self%VLMinited(3) = .true.
          end if
          call self%init_VM_vlvs(self%pVfunc_VS3, self%VLM3, self%VM3)
        else
          call self%init_VM(self%pVfunc_V3, self%VM3)
        end if
      case (4)
        if (allocated(self%VM4)) deallocate(self%VM4)
        allocate(self%VM4(1:self%totN, 1:self%totN))
        if (self%vlvs_separated(4)) then
          if (.not. self%VLMinited(4)) then
            if (allocated(self%VLM4)) deallocate(self%VLM4)
            allocate(self%VLM4(1:self%totN, 1:self%totN))
            call self%init_VM(self%pVfunc_VL4, self%VLM4)
            self%VLMinited(4) = .true.
          end if
          call self%init_VM_vlvs(self%pVfunc_VS4, self%VLM4, self%VM4)
        else
          call self%init_VM(self%pVfunc_V4, self%VM4)
        end if
      case default
        self%VMinited(lbl_V) = .false.
        return
    end select
    self%VMinited(lbl_V) = .true.

  end subroutine init_pots_smplchn

  ! Reset some attributes to its default value.
  ! Do *not* reset the following attributes
  !   - uptoQn, chnid, NCH
  !   - VLMinited(:), vanishing(:), and vlvs_separated(:)
  !   - pVfuncV*

  subroutine release_smplchn(self)

    class(obj_smplchn), intent(inout)   :: self

    self%N              = DEFAULT_N
    self%totN           = DEFAULT_N
    self%regtype        = DEFAULT_REGTYPE
    self%k              = DEFAULT_k
    self%Mambda         = DEFAULT_MAMBDA
    self%k_defined      = .false.
    self%TQn_updated    = .false.
    self%phsn_updated   = .false.
    self%unitarity_warning_enabled = .true.
    self%VMinited(0:8)  = .false.

  end subroutine release_smplchn

  subroutine reset_VLMinited_smplchn(self)

    class(obj_smplchn), intent(inout)   :: self

    self%VLMinited(0:8)  = .false.

  end subroutine reset_VLMinited_smplchn

  ! Call %release()
  ! Reset %chnid, %uptoQn, and %created to their default values.

  subroutine erase_smplchn(self)

    class(obj_smplchn), intent(inout)   :: self

    call self%release()

    if (allocated(self%msh)) deallocate(self%msh)
    if (allocated(self%wght)) deallocate(self%wght)
    if (allocated(self%lsqr)) deallocate(self%lsqr)
    if (allocated(self%VMLO)) deallocate(self%VMLO)
    if (allocated(self%VM1)) deallocate(self%VM1)
    if (allocated(self%VM2)) deallocate(self%VM2)
    if (allocated(self%VM3)) deallocate(self%VM3)
    if (allocated(self%VM4)) deallocate(self%VM4)
    if (allocated(self%VLM2)) deallocate(self%VLM2)
    if (allocated(self%VLM3)) deallocate(self%VLM3)
    if (allocated(self%VLM4)) deallocate(self%VLM4)
    if (allocated(self%VMtmp)) deallocate(self%VMtmp)
    if (allocated(self%Tk)) deallocate(self%Tk)

    self%chnid   = CHNID_UNDEFINED
    self%uptoQn  = DEFAULT_UPTOQN
    self%created = .false.

  end subroutine erase_smplchn

! erase for fitting

  subroutine fit_erase_smplchn(self)

    class(obj_smplchn), intent(inout)   :: self

    call self%release()

    if (allocated(self%msh)) deallocate(self%msh)
    if (allocated(self%wght)) deallocate(self%wght)
    if (allocated(self%lsqr)) deallocate(self%lsqr)
    if (allocated(self%VMLO)) deallocate(self%VMLO)
    if (allocated(self%VM1)) deallocate(self%VM1)
    if (allocated(self%VM2)) deallocate(self%VM2)
    if (allocated(self%VM3)) deallocate(self%VM3)
    if (allocated(self%VM4)) deallocate(self%VM4)
!     if (allocated(self%VLM2)) deallocate(self%VLM2)
!     if (allocated(self%VLM3)) deallocate(self%VLM3)
!     if (allocated(self%VLM4)) deallocate(self%VLM4)
    if (allocated(self%VMtmp)) deallocate(self%VMtmp)
    if (allocated(self%Tk)) deallocate(self%Tk)

    self%chnid   = CHNID_UNDEFINED
    self%uptoQn  = DEFAULT_UPTOQN
    self%created = .false.

  end subroutine fit_erase_smplchn

  subroutine get_paras_filename_smplchn(self, fname)

    class(obj_smplchn), intent(inout)   :: self
    character(len = *), intent(out)     :: fname

    fname = 'inputs_'//self%chnstr//'.in'

  end subroutine get_paras_filename_smplchn

  subroutine get_num_paras_smplchn(self, num_cnttrm)

    class(obj_smplchn), intent(inout)   :: self
    integer, intent(out)                :: num_cnttrm

    num_cnttrm = self%npar(self%uptoQn + 1)

  end subroutine get_num_paras_smplchn

  ! Format of .in :
  ! Mambda <separator> para(1) <separator> para(2) ...

  subroutine read_inputfile_smplchn(self, funt, num_cnttrms, nLmbds,           &
    & Lambda_cnttrms, succeeded)

    class(obj_smplchn), intent(inout)   :: self
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
        'read_inputfile_smplchn: couldn''t open or read '//trim(fname)
      succeeded = .false.
    else
      succeeded = .true.
    end if

  end subroutine read_inputfile_smplchn

  ! Load inputs onto self%para. By default, take all input
  ! parameters as C's.

  subroutine load_inputs_smplchn(self, num_inputs, inputs)

    class(obj_smplchn), intent(inout)   :: self
    integer, intent(in)                 :: num_inputs
    real(NER), dimension(:), intent(in) :: inputs

    self%input(1:num_inputs) = inputs(1:num_inputs)
    self%potpara(1:num_inputs) = self%input(1:num_inputs)

  end subroutine load_inputs_smplchn

  function half_regulator_smplchn(self, p)

    real(NER)                       :: half_regulator_smplchn
    class(obj_smplchn), intent(in)  :: self
    real(NER), intent(in)           :: p

    half_regulator_smplchn = regltr_half_PP(self%regtype, p, self%Mambda)

  end function half_regulator_smplchn

  function regulator_smplchn(self, p1, p2)

    real(NER)                       :: regulator_smplchn
    class(obj_smplchn), intent(in)  :: self
    real(NER), intent(in)           :: p1, p2

    regulator_smplchn = regltr_PP(self%regtype, p1, p2, self%Mambda)

  end function regulator_smplchn

  ! Display warning message in case unitarity broken

  subroutine convert_TQn_to_phase_smplchn(self)

    class(obj_smplchn), intent(inout)   :: self

    integer               :: ii
    character(len = 1024) :: msg
    logical               :: is_broken
    complex(NEC)          :: violation(2, 2)

    if (.not. self%TQn_updated) then
      write (standard_error_unit, '(a)')  &
        'convert_TQn_to_phase_smplchn: Warning: TQn is not updated.'
      return
    end if

    if (self%NCH == 1) then
      call get_delta_np_sngl(self%TQn(1, 1, 1), self%k, self%univio(1),        &
        & self%uflag(1), self%phsn(1))
      do ii = 1, self%uptoQn
        call get_nth_delta_sngl(ii, self%TQn(1, 1, ii+1), self%phsn(1:ii),     &
          & self%k, self%uflag(ii+1), self%univio(ii+1), self%phsn(ii+1))
      end do
    else
      call t2d1d2e_np(self%k, self%TQn(1:2, 1:2, 1), self%trpln(1),            &
        & self%uflag(1))
      ! for 3S1, adjust d1 and mixing angle
      if (self%chnid .eq. CHNID_3S1_PP) then
        if (self%k/PC_mpi .lt. 1.4_NER .and. self%trpln(1)%d1 .lt. 0.0_NER) then
          self%trpln(1)%d1 = self%trpln(1)%d1 + 180.0_NER
          self%trpln(1)%e = -self%trpln(1)%e
        end if
      end if
      ! for 3P2, adjust d1
      if (self%chnid .eq. CHNID_3P2_PP) then
        if (self%trpln(1)%d1 .lt. 0.0_NER) then
          self%trpln(1)%d1 = self%trpln(1)%d1 + 180.0_NER
          self%trpln(1)%e = -self%trpln(1)%e
        end if
      end if
      do ii = 1, self%uptoQn
        call get_nth_delta_cpld(ii, self%TQn(1:2, 1:2, ii+1),                  &
          & self%trpln(1:ii), self%k, self%uflag(ii+1), violation,             &
          & self%trpln(ii+1))
        self%univio(ii+1) = abs(violation(1, 1)) + abs(violation(1, 2))        &
          & + abs(violation(2, 1)) + abs(violation(2, 2))
      end do
    end if
    self%phsn_updated = .true.

    if (self%unitarity_warning_enabled) then
      call self%GenerateUnitarityMsg(is_broken, msg)
      if (is_broken) then
        write(standard_error_unit, '(e10.4e1, x, e10.4e1, x, a3, x, a)')       &
          & self%k, self%Mambda, self%chnstr, trim(msg)
      end if
    end if

  end subroutine convert_TQn_to_phase_smplchn

  subroutine convert_phsn_to_text_smplchn(self, outstr)

    class(obj_smplchn), intent(inout)  :: self
    character(len = *), intent(out) :: outstr

    select case (self%NCH)
      case (1)
        call write_incremental_phsn_to_text_sngl(self%uptoQn+1, self%phsn, outstr)
      case (2)
        call write_incremental_phsn_to_text_cpld(self%uptoQn+1, self%trpln, outstr)
      case default
        write(standard_error_unit, '(a)')   &
          & 'convert_phsn_to_text_smplchn: type of channel not recognized'
    end select

  end subroutine convert_phsn_to_text_smplchn

  ! out_phs(n) is the phase shift correction for uptoQn=(n-1)
  ! Eg., out_phs(1) is the LO, out_phs(2) is the first CORRECTION, and so on

  subroutine get_incremental_phsn_value_smplchn(self, out_phs)

    class(obj_smplchn), intent(inout)  :: self
    class(*), intent(inout)         :: out_phs(:)

    integer :: ii

    select type (out_phs)
      type is (triplet_phs_cpld)
        if (self%NCH /= 2) then
          write(standard_error_unit, '(a)')  &
            & 'get_incremental_phsn_value_smplchn: out_phs is not for coupled chn'
          return
        end if
        do ii = 1, self%uptoQn+1
          out_phs(ii)%d1 = self%trpln(ii)%d1
          out_phs(ii)%d2 = self%trpln(ii)%d2
          out_phs(ii)%e  = self%trpln(ii)%e
        end do
      type is (real(NER))
        if (self%NCH /= 1) then
          write(standard_error_unit, '(a)')  &
            & 'get_incremental_phsn_value_smplchn: out_phs is not for uncoupled chn'
          return
        end if
        out_phs(1:self%uptoQn+1) = self%phsn(1:self%uptoQn+1)
      class default
        write (standard_error_unit, '(a)')  &
          'get_incremental_phsn_value_smplchn: type of out_phs not recognized.'
    end select

  end subroutine get_incremental_phsn_value_smplchn

  ! Generate a warning message, according to unitarity check flags given by
  ! self%uflag(:).

  subroutine GenerateUnitarityMsg_smplchn(self, is_broken, msg)

    class(obj_smplchn), intent(inout)   :: self
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

  end subroutine GenerateUnitarityMsg_smplchn

  ! VM presumed real and symmetric

  subroutine init_VM_smplchn(self, Vfunc, VM)

    class(obj_smplchn), intent(inout)  :: self
    procedure(gnrc_eft_pot)            :: Vfunc
    real(NER), intent(inout)           :: VM(:, :)

    integer     :: ii, jj
    ! Vfunc only outputs 2X2 matrix
    real(NER)   :: ppi, ppj, potval(1:MAX_NCH, 1:MAX_NCH),                     &
      & tt(1:MAX_NCH, 1:MAX_NCH)

    do ii = 1, self%N
      ppi = self%msh(ii)
      do jj = 1, ii - 1
        ppj = self%msh(jj)
        call Vfunc(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype,           &
          & self%Mambda, self%potpara, ppi, ppj, potval)
        call SetLargeVM(VM, self%NCH, self%N+1, ii, jj, potval)
        tt = transpose(potval)
        call SetLargeVM(VM, self%NCH, self%N+1, jj, ii, tt)
      end do
      call Vfunc(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype,             &
        & self%Mambda, self%potpara, ppi, ppi, potval)
      call SetLargeVM(VM, self%NCH, self%N+1, ii, ii, potval)
    end do

  end subroutine init_VM_smplchn

  ! VM presumed real and symmetric

  subroutine init_VM_vlvs_smplchn(self, VSfunc, VLM, VMtot)

    class(obj_smplchn), intent(inout) :: self
    procedure(gnrc_eft_pot)           :: VSfunc
    real(NER), intent(inout)          :: VLM(:, :), VMtot(:, :)

    call self%init_VM(VSfunc, VMtot)
    VMtot(1:self%totN, 1:self%totN) = VMtot(1:self%totN, 1:self%totN) +        &
      & VLM(1:self%totN, 1:self%totN)

  end subroutine init_VM_vlvs_smplchn

  ! VM presumed real and symmetric

  subroutine update_edge_VM_smplchn(self, Vfunc, VM)

    class(obj_smplchn), intent(inout)  :: self
    procedure(gnrc_eft_pot)            :: Vfunc
    real(NER), intent(inout)           :: VM(:, :)

    integer     :: jj, N
    real(NER)   :: ppj, k, potval(1:MAX_NCH, 1:MAX_NCH),            &
      & tt(1:MAX_NCH, 1:MAX_NCH)

    if (.not. (self%k_defined)) then
      write(standard_error_unit, '(a)')                             &
        'update_edge_VM_smplchn: k is not defined. Nothing will be done.'
      return
    end if
    k = self%k
    N = self%N

    do jj = 1, N
      ppj = self%msh(jj)
      call Vfunc(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype,  &
          & self%Mambda, self%potpara, k, ppj, potval)
      call SetLargeVM(VM, self%NCH, N+1, N+1, jj, potval)
      tt = transpose(potval)
      call SetLargeVM(VM, self%NCH, N+1, jj, N+1, tt)
    end do

    call Vfunc(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype, self%Mambda,  &
      & self%potpara, k, k, potval)
    call SetLargeVM(VM, self%NCH, N+1, N+1, N+1, potval)

  end subroutine update_edge_VM_smplchn

  subroutine update_k_smplchn(self, k)

    class(obj_smplchn), intent(inout)   :: self
    real(NER), intent(in)               :: k

    self%k = k
    self%k_defined = .true.

    if (self%VMinited(0)) call self%update_edge_VM(self%pVfunc_V0, self%VMLO)
    if (self%VMinited(1)) call self%update_edge_VM(self%pVfunc_V1, self%VM1)
    if (self%VMinited(2)) call self%update_edge_VM(self%pVfunc_V2, self%VM2)
    if (self%VMinited(3)) call self%update_edge_VM(self%pVfunc_V3, self%VM3)
    if (self%VMinited(4)) call self%update_edge_VM(self%pVfunc_V4, self%VM4)

  end subroutine update_k_smplchn

  subroutine get_Tk_from_VM_smplchn(self, VM, Tk)

    class(obj_smplchn), intent(inout) :: self
    real(NER), intent(in)             :: VM(:, :)
    complex(NEC), intent(inout)       :: Tk(:, :)

    integer                           :: flag

    if (.not. self%k_defined) &
      & write(standard_error_unit, '(a)')   &
        & 'get_Tk_from_VM_smplchn: Warning! k not defined.'
    call lsesv_nch_half_cpld(self%NCH, self%k, self%N, VM, self%msh,           &
      & self%wght, self%lsqr, Tk, flag)

  end subroutine get_Tk_from_VM_smplchn

  ! * uptoQn <= 4.
  ! * Must rewrite for uptoQn > 4
  ! * Checks several prerequisites before overridden method calculates real
  ! amplitudes.

  subroutine get_allTQn_smplchn(self)

    class(obj_smplchn), intent(inout)   :: self

    integer :: NCH, nins_VQ2
    real(NER) :: tmpVQ(1:MAX_NCH, 1:MAX_NCH)
    complex(NEC), dimension(1:MAX_NCH, 1:MAX_NCH, 1:MAX_NUMBER_ORDER)  ::      &
      & TQV1, TQV2, TQV3, TQV4
    complex(NEC), dimension(1:MAX_NCH, 1:MAX_NCH, 1:MAX_NUMBER_ORDER,          &
      & 1:MAX_NUMBER_ORDER)  ::  TQV1V2, TQV1V3

    ! Important to set zero the following arrays
    TQV1 = 0.0_NER
    TQV2 = 0.0_NER
    TQV3 = 0.0_NER
    TQV4 = 0.0_NER
    TQV1V2 = 0.0_NER
    TQV1V3 = 0.0_NER
    NCH = self%NCH

    if (.not. self%k_defined) then
      write(standard_error_unit, '(a)')                       &
        & 'get_allTQn_smplchn: Warning! k not defined, using k = 0.1 MeV'
      call self%update_k(0.1_NER)
    end if

    if (.not. self%vanishing(0)) then
      if (.not. self%VMinited(0)) then
        call self%init_pots(0)
        call self%update_edge_VM(self%pVfunc_V0, self%VMLO)
      end if
      call self%get_Tk_from_VM(self%VMLO, self%Tk)
      call GetTQFromTk(self%Tk, NCH, self%N+1, self%N+1,      &
        & self%TQn(1:NCH, 1:NCH, 1))
      ! debug
      ! print *, "k = ", self%k
      ! print *, "smpl TQ0 = ", self%TQn(1:NCH, 1:NCH, 1)
    else
      self%TQn(1:NCH, 1:NCH, 1) = 0.0_NER
    end if
    if (self%uptoQn == 0) then
      self%TQn_updated = .true.
      return
    end if

    if (.not. self%vanishing(1)) then
      if (self%uptoQn >= 2 .or. (.not. self%vanishing(0)))  then
        if (.not. self%VMinited(1)) then
          call self%init_pots(1)
          call self%update_edge_VM(self%pVfunc_V1, self%VM1)
        end if
        call self%GetSinglePertV(self%VMLO, self%uptoQn, self%VM1,             &
          & self%VMtmp, TQV1, -self%V1_cheb_bnd, self%V1_cheb_bnd,             &
          & self%V1_cheb_base_n)
      end if
      if (self%vanishing(0)) then
        call self%pVfunc_V1(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype,  &
          self%Mambda, self%potpara, self%k, self%k, tmpVQ)
        TQV1(1:NCH, 1:NCH, 2) = tmpVQ(1:NCH, 1:NCH)
      end if
    end if
    self%TQn(1:NCH, 1:NCH, 2) = TQV1(1:NCH, 1:NCH, 2)
    if (self%uptoQn == 1) then
      self%TQn_updated = .true.
      return
    end if

    if (mod(self%uptoQn, 2) == 1)  then
      nins_VQ2 = (self%uptoQn - 1) / 2
    else
      nins_VQ2 = self%uptoQn / 2
    end if
    if (.not. self%vanishing(2)) then
      if (nins_VQ2 >= 2 .or. (.not. self%vanishing(0))) then
        if (.not. self%VMinited(2)) then
          call self%init_pots(2)
          call self%update_edge_VM(self%pVfunc_V2, self%VM2)
        end if
        call self%GetSinglePertV(self%VMLO, nins_VQ2, self%VM2, self%VMtmp,    &
          & TQV2, -self%V2_cheb_bnd, self%V2_cheb_bnd, self%V2_cheb_base_n)
      end if
      if (self%vanishing(0)) then
        call self%pVfunc_V2(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype,  &
          self%Mambda, self%potpara, self%k, self%k, tmpVQ)
        TQV2(1:NCH, 1:NCH, 2) = tmpVQ(1:NCH, 1:NCH)
      end if
    end if
    self%TQn(1:NCH, 1:NCH, 3) = TQV1(1:NCH, 1:NCH, 3) + TQV2(1:NCH, 1:NCH, 2)
!     debug
! print *, 'TQn = ', self%TQn(1, 1, 3), 'TQV1 = ', TQV1(1, 1, 3)

    if (self%uptoQn == 2) then
      self%TQn_updated = .true.
      return
    end if

    if (.not. self%vanishing(3)) then
      if (.not. self%vanishing(0)) then
        if (.not. self%VMinited(3)) then
          call self%init_pots(3)
          call self%update_edge_VM(self%pVfunc_V3, self%VM3)
        end if
        call self%GetSinglePertV(self%VMLO, 1, self%VM3, self%VMtmp, TQV3,     &
          & -self%V3_cheb_bnd, self%V3_cheb_bnd, self%V3_cheb_base_n)
      else
        call self%pVfunc_V3(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype,  &
          self%Mambda, self%potpara, self%k, self%k, tmpVQ)
        TQV3(1:NCH, 1:NCH, 2) = tmpVQ(1:NCH, 1:NCH)
      end if
    end if
    if (.not. (self%vanishing(1) .or. self%vanishing(2))) then
      if (.not. self%VMinited(1)) then
        call self%init_pots(1)
        call self%update_edge_VM(self%pVfunc_V1, self%VM1)
      end if
      if (.not. self%VMinited(2)) then
        call self%init_pots(2)
        call self%update_edge_VM(self%pVfunc_V2, self%VM2)
      end if
      call self%GetDoublePertV(self%VMLO, 2, self%VM1, 1, self%VM2,            &
        & self%VMtmp, TQV1V2,                                                  &
        & -self%V2_cheb_bnd, self%V2_cheb_bnd, self%V2_cheb_base_n,            &
        & -self%V2_cheb_bnd, self%V2_cheb_bnd, self%V2_cheb_base_n)
    end if
    self%TQn(1:NCH, 1:NCH, 4) = TQV1(1:NCH, 1:NCH, 4)                          &
      & + TQV1V2(1:NCH, 1:NCH, 2, 2) + TQV3(1:NCH, 1:NCH, 2)

    if (self%uptoQn == 3) then
      self%TQn_updated = .true.
      return
    end if

    if (.not. self%vanishing(4)) then
      if (.not. self%vanishing(0)) then
        if (.not. self%VMinited(4)) then
          call self%init_pots(4)
          call self%update_edge_VM(self%pVfunc_V4, self%VM4)
        end if
        call self%GetSinglePertV(self%VMLO, 1, self%VM4, self%VMtmp, TQV4,     &
          & -self%V4_cheb_bnd, self%V4_cheb_bnd, self%V4_cheb_base_n)
      else
        call self%pVfunc_V4(self%lsj%L, self%lsj%S, self%lsj%j, self%regtype,  &
          self%Mambda, self%potpara, self%k, self%k, tmpVQ)
        TQV4(1:NCH, 1:NCH, 2) = tmpVQ(1:NCH, 1:NCH)
      end if
    end if
    if (.not. (self%vanishing(1) .or. self%vanishing(3))) then
      if (.not. self%VMinited(1)) call self%init_pots(1)
      if (.not. self%VMinited(3)) then
        call self%init_pots(3)
        call self%update_edge_VM(self%pVfunc_V3, self%VM3)
      end if
      call self%GetDoublePertV(self%VMLO, 1, self%VM1, 1, self%VM3,            &
        & self%VMtmp, TQV1V3,                                                  &
        & -self%V2_cheb_bnd, self%V2_cheb_bnd, self%V2_cheb_base_n,            &
        & -self%V3_cheb_bnd, self%V3_cheb_bnd, self%V3_cheb_base_n)
    end if
    self%TQn(1:NCH, 1:NCH, 5) = TQV1(1:NCH, 1:NCH, 5)                          &
      & + TQV1V2(1:NCH, 1:NCH, 3, 2) + TQV1V3(1:NCH, 1:NCH, 2, 2)              &
      & + TQV2(1:NCH, 1:NCH, 3) + TQV4(1:NCH, 1:NCH, 2)

    if (self%uptoQn == 4) then
      self%TQn_updated = .true.
      return
    end if

    self%TQn_updated = .true.

  end subroutine get_allTQn_smplchn

  ! Mapping {a, b} -> 2*((aa-1)*NCH + bb) - 1

  pure function loc_series(aa, bb, NCH)

    integer, intent(in) :: aa, bb, NCH
    integer :: loc_series

    loc_series = 2*((aa-1)*NCH + bb) - 1
  end function loc_series

  ! n insertions of VM_pert. a and b are the boundary for the dimensionless
  ! parameter x, by which the Taylor expansion is defined (see
  ! cheb_taylor_real). cheb_base+n+1 is the actual number of insertions
  ! calculated by cheb_taylor subroutines. TQTayS(i) is the value with (i-1)
  ! insertions.

  subroutine GetSinglePertV_smplchn(self, VMlo, n, VM_pert, VMtmp, TQTayS, a,  &
    & b, cheb_base_n)

    class(obj_smplchn), intent(inout)   :: self
    ! real(NER), intent(in), dimension(:, :), pointer :: VMlo, VM_pert
    real(NER), intent(in), dimension(:, :)  :: VMlo, VM_pert
    real(NER), intent(out)              :: VMtmp(:, :)
    complex(NEC), intent(out)           :: TQTayS(:, :, :)
    integer, intent(in)                 :: n, cheb_base_n
    real(NER), intent(in)               :: a, b

    real(NER) :: g(MAX_NCH*MAX_NCH*2, cheb_base_n+n+5)
    integer   :: ii, aa, bb, NCH, totN, tmploc
    complex(NEC)  :: tmpTQ(1:MAX_NCH, 1:MAX_NCH)

    NCH = self%NCH
    totN = self%totN

    if (NCH == 2) then
      call cheb_taylor_multi(a, b, 3*2, cheb_base_n+n+1, g, func_amp)
    else
      call cheb_taylor_multi(a, b, 2, cheb_base_n+n+1, g, func_amp)
    end if

    do ii = 1, n+1
      if (NCH == 2) then
        TQTayS(1, 1, ii) = cmplx(g(1, ii), g(2, ii), NEC)
        TQTayS(2, 2, ii) = cmplx(g(3, ii), g(4, ii), NEC)
        TQTayS(1, 2, ii) = cmplx(g(5, ii), g(6, ii), NEC)
        TQTayS(2, 1, ii) = TQTayS(1, 2, ii)
      else
        TQTayS(1, 1, ii) = cmplx(g(1, ii), g(2, ii), NEC)
      end if
    end do

  contains

    subroutine func_amp(x, fnval)

      real(NER), intent(in)   :: x
      real(NER), intent(out)  :: fnval(:)

      integer       :: aa_in, bb_in, loc

      if (self%vanishing(0)) then
        VMtmp(1:totN, 1:totN) = x*VM_pert(1:totN, 1:totN)
      else
        VMtmp(1:totN, 1:totN) = VMlo(1:totN, 1:totN) + x*VM_pert(1:totN, 1:totN)
      end if
      call self%get_Tk_from_VM(VMtmp, self%Tk)
      call GetTQFromTk(self%Tk, NCH, self%N+1, self%N+1, tmpTQ)

      if (NCH == 2) then
        fnval(1) = real(tmpTQ(1, 1))
        fnval(2) = aimag(tmpTQ(1, 1))
        fnval(3) = real(tmpTQ(2, 2))
        fnval(4) = aimag(tmpTQ(2, 2))
        fnval(5) = real(tmpTQ(1, 2))
        fnval(6) = aimag(tmpTQ(1, 2))
      else
        fnval(1) = real(tmpTQ(1, 1))
        fnval(2) = aimag(tmpTQ(1, 1))
      end if

    end subroutine func_amp

  end subroutine GetSinglePertV_smplchn

  ! TQn(1:2, 1:2, i, j) = (i-1) insertions of v1 and (j-1) insertions of v2
  ! Note: TQ(1, 1) is expected to be zero if the accuracy is infinite, so its
  ! returned value serves as a metric for actual numerical accuracy.

  subroutine GetDoublePertV_smplchn(self, VMlo, n1, VM1, n2, VM2, VMtmp,       &
    & TQTayS, a1, b1, cheb_base_n1, a2, b2, cheb_base_n2)

    class(obj_smplchn), intent(inout) :: self
    real(NER), intent(in)             :: VMlo(:, :), VM1(:, :), VM2(:, :)
    real(NER), intent(out)            :: VMtmp(:, :)
    complex(NEC), intent(out)         :: TQTayS(:, :, :, :)
    integer, intent(in)               :: n1, n2, cheb_base_n1, cheb_base_n2
    real(NER), intent(in)             :: a1, b1, a2, b2

    real(NER) :: g(MAX_NCH*MAX_NCH*2, cheb_base_n1+n1+5, cheb_base_n2+n2+5)
    integer   :: ii, jj, aa, bb, NCH, totN, tmploc
    complex(NEC)  :: tmpTQ(1:MAX_NCH, 1:MAX_NCH)

    NCH = self%NCH
    totN = self%totN

    if (NCH == 2) then
      call cheb_taylor_multi_xy(a1, b1, cheb_base_n1+n1+1, a2, b2,             &
      & cheb_base_n2+n2+1, 3*2, g, func_amp)
    else
      call cheb_taylor_multi_xy(a1, b1, cheb_base_n1+n1+1, a2, b2,             &
      & cheb_base_n2+n2+1, 2, g, func_amp)
    end if

    do ii = 1, n1+1
      do jj = 1, n2+1
        if (NCH == 2) then
          TQTayS(1, 1, ii, jj) = cmplx(g(1, ii, jj), g(2, ii, jj), NEC)
          TQTayS(2, 2, ii, jj) = cmplx(g(3, ii, jj), g(4, ii, jj), NEC)
          TQTayS(1, 2, ii, jj) = cmplx(g(5, ii, jj), g(6, ii, jj), NEC)
          TQTayS(2, 1, ii, jj) = TQTayS(1, 2, ii, jj)
        else
          TQTayS(1, 1, ii, jj) = cmplx(g(1, ii, jj), g(2, ii, jj), NEC)
        end if
      end do
    end do

  contains

    subroutine func_amp(x, y, fnval)

      real(NER), intent(in)  :: x, y
      real(NER), intent(out) :: fnval(:)

      integer       :: aa_in, bb_in, loc

      if (self%vanishing(0)) then
        VMtmp(1:totN, 1:totN) = x*VM1(1:totN, 1:totN) + y*VM2(1:totN, 1:totN)
      else
        VMtmp(1:totN, 1:totN) = VMlo(1:totN, 1:totN) + x*VM1(1:totN, 1:totN)   &
        & + y*VM2(1:totN, 1:totN)
      end if
      call self%get_Tk_from_VM(VMtmp, self%Tk)
      call GetTQFromTk(self%Tk, NCH, self%N+1, self%N+1, tmpTQ)

      if (NCH == 2) then
        fnval(1) = real(tmpTQ(1, 1))
        fnval(2) = aimag(tmpTQ(1, 1))
        fnval(3) = real(tmpTQ(2, 2))
        fnval(4) = aimag(tmpTQ(2, 2))
        fnval(5) = real(tmpTQ(1, 2))
        fnval(6) = aimag(tmpTQ(1, 2))
      else
        fnval(1) = real(tmpTQ(1, 1))
        fnval(2) = aimag(tmpTQ(1, 1))
      end if

    end subroutine func_amp

  end subroutine GetDoublePertV_smplchn

! get cheby parameters from file 'cheby_paras.in'

  subroutine get_cheby_paras(chn)
    class(obj_smplchn), intent(inout) :: chn
    integer   ::  chebyn(1:4), flag
    real(NER) ::  bnd(1:4)
    open(unit=57,file="cheby_paras.in")
    if (chn%uptoQn == 0) return
    if (chn%uptoQn > 0) then
      do
        read(57,*,iostat=flag)chebyn(1:chn%uptoQn)
        read(57,*,iostat=flag)bnd(1:chn%uptoQn)
      if (flag/=0) exit
      end do
    end if

    chn%V1_cheb_base_n = chebyn(1)
    chn%V1_cheb_bnd = bnd(1)
    if (chn%uptoQn == 1) then
      close(57)
      return
    end if
    chn%V2_cheb_base_n = chebyn(2)
    chn%V2_cheb_bnd = bnd(2)
    if (chn%uptoQn == 2) then
      close(57)
      return
    end if
    chn%V3_cheb_base_n = chebyn(3)
    chn%V3_cheb_bnd = bnd(3)
    if (chn%uptoQn == 3) then
      close(57)
      return
    end if
    chn%V4_cheb_base_n = chebyn(4)
    chn%V4_cheb_bnd = bnd(4)
    close(57)
  end subroutine get_cheby_paras

  subroutine get_BDE_smplchn(chn, uptoQn, x, BDE)

    class(obj_smplchn), intent(inout) :: chn
    real(NER), intent(out)            :: BDE
    real(NER), intent(in)             :: x
    integer, intent(in)               :: uptoQn

    real(NER)          :: Rwf(2*chn%N)
    real(NER)          :: Rmsh(chn%N), Rwght(chn%N), RMRatio
!  get Cmsh and Cwght
    RMRatio = 0.5_NER
    call feb_tanabsc(chn%N, chn%Mambda*RMRatio, Rmsh, Rwght)
!     Cmsh = chn%msh
!     Cwght = chn%wght
    call get_h2wfs_real(chn%Mambda, chn%potpara, chn%regtype,  &
      &  chn%N, Rmsh, Rwght, pc_mN, pot, BDE, Rwf)

  contains
    subroutine pot(L, S, j, regtype, Mambda, para, p1, p2, potval)

      integer, intent(in)     :: L, S, J, regtype
      real(NER), intent(in)   :: Mambda, p1, p2, para(:)
      real(NER), intent(out)  :: potval(1:MAX_NCH, 1:MAX_NCH)

      real(NER)               :: V0(1:MAX_NCH, 1:MAX_NCH), V2(1:MAX_NCH, 1:MAX_NCH)

      potval = 0.0_NER
      call chn%pVfunc_V0(L, S, j, regtype, Mambda, para, p1, p2, V0)
      select case(uptoQn)
        case(2)
          call chn%pVfunc_V2(L, S, j, regtype, Mambda, para, p1, p2, V2)
          potval = V0 + x * V2
!           potval = V0 + x**2 * V2
        case default
          potval = V0
      end select
    end subroutine pot

  end subroutine get_BDE_smplchn

  subroutine get_pertBDE_smplchn(chn, n, a, b, cheb_base_n, pert_BDE)

    class(obj_smplchn), intent(inout) :: chn
    real(NER), intent(in)                   :: a, b
    integer, intent(in)                     :: n, cheb_base_n
    real(NER), intent(out)            :: pert_BDE(2)

    real(NER)     :: g(cheb_base_n+n+5)

    g = 0.0_NER
    call cheb_taylor_real(a, b, cheb_base_n+n+1, g, BDE)

    pert_BDE = g(1:2)

  contains
    function BDE(x)
      real(NER), intent(in)  :: x
      real(NER)              :: BDE

      call get_BDE_smplchn(chn, chn%uptoQn, x, BDE)
    end function BDE
  end subroutine get_pertBDE_smplchn

  subroutine get_pertpara_smplchn(chn, xd, xu, rslt)

    class(obj_smplchn), intent(inout) :: chn
    real(NER),intent(in)  :: xd, xu
    real(NER),intent(out) :: rslt
    real(NER)             :: tol

!     tol =1.0E-12_NER
    tol =1.0E-9_NER
    rslt = zbrent(pert_E, xd, xu, tol)

  contains
    function pert_E(cs)
      real(NER),intent(in)  :: cs
      real(NER)             :: pert_E

      real(NER)          :: tmp_pertBDE(2)

      chn%potpara(chn%npar(chn%uptoQn)+1) = cs

      call get_pertBDE_smplchn(chn, chn%uptoQn, chn%V2_cheb_pertBD_bnd, -chn%V2_cheb_pertBD_bnd, &
      &    chn%V2_cheb_pertBD_n, tmp_pertBDE)
      pert_E = tmp_pertBDE(2)
    end function
  end subroutine get_pertpara_smplchn

  subroutine reset_para_smplchn(chn)

    class(obj_smplchn), intent(inout) :: chn

    select case (chn%chnid)
      case(CHNID_3S1_PP)
        if (chn%uptoqn >= 2) then
          call get_pertpara2_smplchn(chn)
          return
        end if
      case(CHNID_1S0_PP)
        if (chn%uptoqn == 1) then
          call get_NLOpertSCL_smplchn(chn)
        end if
        if (chn%uptoqn == 2) then
          call get_N2LOpertSCL_smplchn(chn)
        end if
!         debug
!         print *, 'NLOpert_scl = ', NLOpert_SCL
!         print *, 'N2LOpert_scl = ', N2LOpert_SCL
        call get_N2LO1s0para_smplchn(chn)
      case default
        write(standard_output_unit, '(a)')    &
          & "-->reset_para_smplchn: reset paras 3s1 and 1s0 only"
        return

      chn%pVfunc_V1 => VNLO_MMWLY_epwrap
!       chn%pVfunc_V2 => VN2LO_MMWLY_epwrap
! debug fit without TPE
        chn%pVfunc_V2 => VN2LO_MMWLY_SHRG_epwrap
    end select

  end subroutine reset_para_smplchn


  subroutine get_allpertBDE_smplchn(chn)

    class(obj_smplchn), intent(inout) :: chn
    real(NER)                   :: a, b
    integer                     :: n, cheb_base_n
    real(NER)                   :: pert_BDE(4)

    real(NER),allocatable       :: g(:)

    if (allpert_BDE_initial) then
      return
    end if
    n = chn%uptoqn
    a = chn%V2_cheb_pertBD_bnd
    b = -chn%V2_cheb_pertBD_bnd
    cheb_base_n = chn%V2_cheb_pertBD_n

    allocate(g(1:cheb_base_n+n+5))

    chn%pVfunc_V2 => TPE0_epwrap
!     chn%pVfunc_V2 => reTPE0_epwrap
    g = 0.0_NER
    call cheb_taylor_real(a, b, cheb_base_n+n+1, g, BDE)
    pert_BDE(1) = g(2)

    chn%pVfunc_V2 => VPNQ0_epwrap
!     chn%pVfunc_V2 => reVPNQ0_epwrap
    g = 0.0_NER
    call cheb_taylor_real(a, b, cheb_base_n+n+1, g, BDE)
    pert_BDE(2) = g(2)

!     chn%pVfunc_V2 => VPNQ2_diag_epwrap
    chn%pVfunc_V2 => reVPNQ2_diag_epwrap
    g = 0.0_NER
    call cheb_taylor_real(a, b, cheb_base_n+n+1, g, BDE)
    pert_BDE(3) = g(2)

!     chn%pVfunc_V2 => VPNQ2_off_epwrap
    chn%pVfunc_V2 => reVPNQ2_off_epwrap
    g = 0.0_NER
    call cheb_taylor_real(a, b, cheb_base_n+n+1, g, BDE)
    pert_BDE(4) = g(2)

    allpert_BDE = pert_BDE
    allpert_BDE_initial = .true.

    ! debug
    write(standard_output_unit, '(a)')    &
      & "-->allpert_BDE_initial: initialized successfully"
    print *, '==============================='
!     print *, 'allpert_BDE =', pert_BDE

  contains
    function BDE(x)
      real(NER), intent(in)  :: x
      real(NER)              :: BDE

      call get_BDE_smplchn(chn, chn%uptoQn, x, BDE)
    end function BDE
  end subroutine get_allpertBDE_smplchn

  subroutine get_pertpara2_smplchn(chn)

    class(obj_smplchn), intent(inout) :: chn
    real(NER)          ::  C, D, E, mambda


    call get_allpertBDE_smplchn(chn)
    mambda = chn%mambda
    D = chn%potpara(chn%npar(chn%uptoQn)+2)
    E = chn%potpara(chn%npar(chn%uptoQn)+3)

!     c = -(allpert_BDE(1) + mambda**2*D*allpert_BDE(3) + mambda**2*E*allpert_BDE(4))/(allpert_BDE(2))
!     debug  fit without TPE
    c = -(mambda**2*D*allpert_BDE(3) + mambda**2*E*allpert_BDE(4))/(allpert_BDE(2))

    chn%potpara(chn%npar(chn%uptoQn)+1) = C

  end subroutine get_pertpara2_smplchn

! get 1s0 scattering length

  subroutine get_1s0sclength_smplchn(chn, uptoQn, sclength)

    class(obj_smplchn), intent(inout) :: chn
    integer,intent(in)                :: uptoQn
    real(NER),intent(out)             :: sclength(1:uptoqn+1)
    real(NER)      :: amplitude(1:uptoQn+1)
    integer        ::  ii

    chn%VMinited = .false.
    call chn%get_allTQn()
    amplitude(1:uptoQn+1) = real(chn%TQn(1,1,1:uptoQn+1))
!     amplitude(1:uptoQn+1) = abs(chn%TQn(1,1,1:uptoQn+1))
!     debug
! print*, 'TQn =', chn%TQn(1, 1, 1: uptoqn+1)
    do ii = 1, uptoQn+1
      sclength(ii) = amplitude(ii)*PI_NE/2.0_NER*PC_default_hbarc
    end do
  end subroutine get_1s0sclength_smplchn

  subroutine get_NLOpertSCL_smplchn(chn)

    class(obj_smplchn), intent(inout) :: chn
    integer                     :: n
    real(NER)                   :: NLOpert_sclength(2)

    if (pert_SCL_initial(1) .or. chn%uptoqn /= 1) then
      return
    end if

    chn%k = 1.0E-12_NER
!     chn%k_defined = .true.
    call chn%update_k(chn%k)
    n = chn%uptoqn
!     chn%pVfunc_V1 => VPNQ0_epwrap
    chn%pVfunc_V1 => reVPNQ0_SCL_epwrap
    call get_1s0sclength_smplchn(chn, n, NLOpert_sclength)
    NLOpert_SCL(1) = NLOpert_sclength(2)

!     chn%pVfunc_V2 => VPNQ2_diag_epwrap
    chn%pVfunc_V1 => reVPNQ2_SCL_diag_epwrap
    call get_1s0sclength_smplchn(chn, n, NLOpert_sclength)
    NLOpert_SCL(2) = NLOpert_sclength(2)

!     chn%k_defined = .false.
    pert_SCL_initial(1) = .true.

    ! debug
    write(standard_output_unit, '(a)')    &
      & "-->pert_SCL_initial(1): initialized successfully"
    print *, '==============================='
    print *, 'NLOpert_SCL =', NLOpert_SCL

  end subroutine get_NLOpertSCL_smplchn

  subroutine get_N2LOpertSCL_smplchn(chn)

    class(obj_smplchn), intent(inout) :: chn
    integer                     :: n
    real(NER)                   :: N2LOpert_sclength(3)

    if (pert_SCL_initial(1) .or. chn%uptoqn /= 2) then
      return
    end if

    n = chn%uptoqn
    chn%k = 1.0E-12_NER
!     chn%k_defined = .true.
    call chn%update_k(chn%k)

    chn%pVfunc_V2 => zero_V
    call get_1s0sclength_smplchn(chn, n, N2LOpert_sclength)
    N2LOpert_SCL(5) = N2LOpert_sclength(3)

    chn%pVfunc_V2 => TPE0_epwrap
!     chn%pVfunc_V2 => reTPE0_epwrap
    call get_1s0sclength_smplchn(chn, n, N2LOpert_sclength)
    N2LOpert_SCL(1) = N2LOpert_sclength(3) - N2LOpert_SCL(5)

!     chn%pVfunc_V2 => VPNQ0_epwrap
    chn%pVfunc_V2 => reVPNQ0_epwrap
    call get_1s0sclength_smplchn(chn, n, N2LOpert_sclength)
    N2LOpert_SCL(2) = N2LOpert_sclength(3) - N2LOpert_SCL(5)

!     chn%pVfunc_V2 => reVPNQ2_SCL_diag_epwrap
    chn%pVfunc_V2 => reVPNQ2_diag_epwrap
    call get_1s0sclength_smplchn(chn, n, N2LOpert_sclength)
    N2LOpert_SCL(3) = N2LOpert_sclength(3) - N2LOpert_SCL(5)

    chn%pVfunc_V2 => VPNQ5_epwrap
    chn%pVfunc_V2 => reVPNQ5_epwrap
    call get_1s0sclength_smplchn(chn, n, N2LOpert_sclength)
    N2LOpert_SCL(4) = N2LOpert_sclength(3) - N2LOpert_SCL(5)

    pert_SCL_initial(2) = .true.
    ! debug
    write(standard_output_unit, '(a)')    &
      & "-->pert_SCL_initial(2): initialized successfully"
!     print *, '==============================='
!     print *, 'N2LOpert_SCL =', N2LOpert_SCL

  end subroutine get_N2LOpertSCL_smplchn

  subroutine get_N2LO1s0para_smplchn(chn)

    class(obj_smplchn), intent(inout) :: chn
    real(NER)        :: c1, d0, c2, d1, e0, mambda
    real(NER)        :: xu, xd, tol

    mambda = chn%mambda
    if(chn%uptoQn == 1) then
      d0 = chn%potpara(chn%npar(chn%uptoQn)+2)
      c1 = -NLOpert_SCL(2)*d0*(mambda**4)/NLOpert_SCL(1)/(mambda**2)
      chn%potpara(chn%npar(chn%uptoQn)+1) = c1
!       debug
print *, 'k =', chn%k
print *, 'NLO_scl = ',  NLOpert_SCL
    print *, 'c1 / d0 = ', c1, d0
    end if
    if(chn%uptoqn == 2) then
      tol = 1.0E-9_NER
      c2 = chn%potpara(chn%npar(chn%uptoQn)+1)
      xd = abs(c2) * 10.0_NER
!       xu = abs(c2) * 1.0E-2_NER
      xu = -xd
      d1 = chn%potpara(chn%npar(chn%uptoQn)+2)
      e0 = chn%potpara(chn%npar(chn%uptoQn)+3)
!       debug
      print *, 'c2 by initin =', c2
      c2 = zbrent(N2LOpert_SCL2, xd, xu, tol)
!       debug
      print *, 'c2 by zbrent =', c2
!       c2 = -(N2LOpert_SCL(1) + N2LOpert_SCL(3)*d1*(Mambda**2) + N2LOpert_SCL(4)*e0*(Mambda**4) &
!         & + N2LOpert_SCL(5))/N2LOpert_SCL(2)/mambda
! !       debug
      print *, 'c2 by linefu =', c2
      chn%potpara(chn%npar(chn%uptoQn)+1) = c2
!       debug
! print *, 'k =', chn%k
print *, 'N2LO_scl = ',  N2LOpert_SCL
    end if
    contains
      function N2LOpert_SCL2(x)
      real(NER), intent(in)      :: x
      real(NER)      :: N2LOpert_SCL2

      real(NER)      :: N2LOpert_sclength(3), k_init, k_SCL

!       chn%pVfunc_V2 => VN2LO_MMWLY_epwrap
! debug fit without TPE
        chn%pVfunc_V2 => VN2LO_MMWLY_SHRG_epwrap
      chn%potpara(chn%npar(chn%uptoQn)+1) = x

      k_init = chn%k
      k_SCL  = 1.0E-12_NER
!     chn%k_defined = .true.
      call chn%update_k(k_SCL)
      call get_1s0sclength_smplchn(chn, chn%uptoqn, N2LOpert_sclength)
      N2LOpert_SCL2 = N2LOpert_sclength(3)
      chn%k = k_init
!       chn%k_defined = .false.

      end function N2LOpert_SCL2
  end subroutine get_N2LO1s0para_smplchn

  subroutine get_LO1s0para_smplchn(chn, xu, xd, C0)

    class(obj_smplchn), intent(inout) :: chn
    real(NER), intent(out)             :: C0
    real(NER), intent(in)        :: xu, xd
    real(NER)                    :: tol

    tol = 1.0E-9_NER
    C0 = zbrent(LO_SCL2, xu, xd, tol)

    contains
      function LO_SCL2(x)
      real(NER), intent(in)      :: x
      real(NER)      :: LO_SCL2

      real(NER)      :: LO_sclength(1)

      chn%pVfunc_V0 => VLO_withc0_epwrap
      chn%potpara(1) = x

      call get_1s0sclength_smplchn(chn, chn%uptoqn, LO_sclength)
      LO_SCL2 = LO_sclength(1) - empirical_scatlength
!       debug
! print*, x, LO_sclength

      end function LO_SCL2
  end subroutine get_LO1s0para_smplchn

!   subroutine get_1s0para_c2_smplchn(d1, e0)

!     real(NER), intent(in)       :: d1, e0
!     class(obj_smplchn), pointer :: chn

!     allocate(obj_smplchn::ptrchn)


!   end subroutine get_1s0para_c2_smplchn
end module mod_obj_smplchn
