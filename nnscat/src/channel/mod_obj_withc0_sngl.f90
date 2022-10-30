! Bingwei Long 07/10/2013
! Class: single channels where C is present at LO and tensor OPE applies. E.g., 3P0

module mod_obj_withc0_sngl

  use mod_obj_sngl
  implicit none

  ! Uncoupled channels where C starts to appear at O(1), D at O(Q^2), etc.
  type, extends(obj_sngl)    :: obj_withc0_sngl

    ! VM_powfac will be used to assemble VLO: VLO = OPE + C * VM_VS0
    type(VM_sngl)   :: VMVLO

  contains

    procedure   :: init             => init_withc0_sngl
    procedure   :: release          => release_withc0_sngl
    ! Frozen
    ! procedure   :: change_a_C       => change_a_C_withc0_sngl
    procedure   :: init_VMVLO       => init_VMVLO_withc0_sngl
    procedure   :: update_k         => update_k_withc0_sngl
    procedure   :: update_edge_VMVLO  => update_edge_VMVLO_withc0_sngl
    procedure   :: get_allTQn       => get_allTQn_withc0_sngl
    procedure   :: get_num_paras    => get_num_paras_withc0_sngl

  end type obj_withc0_sngl

  integer, parameter    :: npar_withc0_sngl(1:4) = (/1, 1, 3, 5/)

contains

  ! Call init_sngl. Create and initialize VM_powfac.
  ! Create VMVLO without initializing it. Instead, initialization of VMVLO will be
  ! deferred until it's used for the first time, e.g., by get_allTQn.

  subroutine init_withc0_sngl(self, N, regtype, Mambda)

    class(obj_withc0_sngl), intent(inout)   :: self
    integer, intent(in)                     :: N, regtype
    real(NER), intent(in)                   :: Mambda

    call init_sngl(self, N, regtype, Mambda)
    call self%create_VM(self%VMVLO)
    call self%create_VM_sngl(self%VM_poly0)
    call self%init_VMpoly0(self%VM_poly0)

  end subroutine init_withc0_sngl


  subroutine release_withc0_sngl(self)

    class(obj_withc0_sngl), intent(inout)  :: self

    call self%release_VM(self%VMVLO)
    call self%release_VM(self%VM_poly0)
    call release_sngl(self)

  end subroutine release_withc0_sngl


  ! Frozen
  ! ! If self%Cs(1) is changed, VMVLO needs to be re-initialized.
  ! subroutine change_a_C_withc0_sngl(self, n, Cn)

  !     class(obj_withc0_sngl), intent(inout)   :: self
  !     integer, intent(in)                     :: n
  !     real(NER), intent(in)                   :: Cn

  !     ! call self%obj_woc0_sngl%change_a_C(n, Cn)
  !     call change_a_C_channel(self, n, Cn)
  !     if (n == 1) then
  !         self%VMVLO%inited = .false.
  !         self%CDT_threshold_set = .false.
  !     end if

  ! end subroutine change_a_C_withc0_sngl


  ! self%VMVLO = self%vs_vl0 + C0 * self%VM_poly0
  subroutine init_VMVLO_withc0_sngl(self, C0)

    class(obj_withc0_sngl), intent(inout)   :: self
    real(NER), intent(in)                   :: C0

    integer :: ii, jj, locN

    locN = self%N
    do ii = 1, locN
      do jj = 1, ii - 1
        self%VMVLO%ptr(ii, jj) = self%vs_vl0%VM%ptr(ii, jj) + C0*self%VM_poly0%ptr(ii, jj)
        self%VMVLO%ptr(jj, ii) = self%VMVLO%ptr(ii, jj)
      end do
      self%VMVLO%ptr(ii, ii) = self%vs_vl0%VM%ptr(ii, ii) + C0*self%VM_poly0%ptr(ii, ii)
    end do
    self%VMVLO%inited = .true.

  end subroutine init_VMVLO_withc0_sngl


  ! In addition to calling update_k_woc0_sngl, update the edge of VMVLO.
  ! Note that when VMVLO's edge is being updated, VMVLO is not necessarily
  ! initialized.

  subroutine update_k_withc0_sngl(self, k)

    class(obj_withc0_sngl), intent(inout)   :: self
    real(NER), intent(in)                   :: k

    call update_k_sngl(self, k)
    call self%update_edge_VMpoly0(self%VM_poly0)
    ! C0 = self%Cs(1)
    call self%update_edge_VMVLO(k, self%Cs(1))

  end subroutine update_k_withc0_sngl

  subroutine update_edge_VMVLO_withc0_sngl(self, k, C0)

    class(obj_withc0_sngl), intent(inout)   :: self
    real(NER), intent(in)                   :: k
    real(NER), intent(in)                   :: C0

    integer :: locN, ii
    real(NER), dimension(:), pointer :: ptr

    locN = self%N
    call self%set_ptr_to_powfac(1, ptr)
    do ii = 1, locN
      self%VMVLO%ptr(ii, locN+1) = self%vs_vl0%VM%ptr(ii, locN+1) + C0*self%VM_poly0%ptr(ii, locN+1)
      self%VMVLO%ptr(locN+1, ii) = self%VMVLO%ptr(ii, locN+1)
    end do
    self%VMVLO%ptr(locN+1, locN+1) = self%vs_vl0%VM%ptr(locN+1, locN+1) + C0*self%VM_poly0%ptr(locN+1, locN+1)

  end subroutine update_edge_VMVLO_withc0_sngl


  ! get_allTQn does not have any prerequisite.

  subroutine get_allTQn_withc0_sngl(self)

    class(obj_withc0_sngl), intent(inout)   :: self

    integer         :: N
    real(NER)       :: VL
    complex(NEC)    :: L1, L2, CT, DT

    if (.not. self%VMVLO%inited) call self%init_VMVLO(self%Cs(1))

    N = self%N
    call self%get_Tk_from_VM(self%VMVLO, self%Tk)
    self%TQn(1) = self%Tk%ptr(N+1)
    if (self%uptoQn .eq. 0) then
      self%TQn_updated = .true.
      return
    end if

    self%TQn(2) = 0.0_NER
    if (self%uptoQn .eq. 1) then
      self%TQn_updated = .true.
      return
    end if

    VL = self%vs_vl2%VM%ptr(N+1, N+1)
    call self%get_1loop_os_VMTk(self%vs_vl2%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl2%VM, self%Tk, L2)
    call self%get_CDT(self%Tk, CT, DT)
    ! C2 = self%Cs(2), D2 = self%Cs(3)
    self%TQn(3) = VL + 2.0_NER*L1 + L2 + self%Cs(2)*CT + self%Cs(3)*DT
    if (self%uptoQn .eq. 2) then
      self%TQn_updated = .true.
      return
    end if

    VL = self%vs_vl3%VM%ptr(N+1, N+1)
    call self%get_1loop_os_VMTk(self%vs_vl3%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl3%VM, self%Tk, L2)
    ! C3 = self%Cs(4), D3 = self%Cs(5)
    self%TQn(4) = VL + 2.0_NER*L1 + L2 + self%Cs(4)*CT + self%Cs(5)*DT
    self%TQn_updated = .true.

  end subroutine get_allTQn_withc0_sngl


  subroutine get_num_paras_withc0_sngl(self, num_cnttrm)

    class(obj_withc0_sngl), intent(inout)   :: self
    integer, intent(out)                    :: num_cnttrm

    select case (self%uptoQn)
      case (0, 1)
        num_cnttrm = 1
      case (2)
        num_cnttrm = 3
      case (3)
        num_cnttrm = 5
      case default
        write (standard_error_unit, '(a, i2)')  &
          'get_num_paras_withc0_sngl: I don''t know how to handle order ', self%uptoQn
        num_cnttrm = 0
    end select

  end subroutine get_num_paras_withc0_sngl


end module mod_obj_withc0_sngl
