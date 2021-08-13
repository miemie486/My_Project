! Bingwei Long 11/20/2017
! Implementation of the 1S0 paper by Sanchez, Yang, Long, and van Kolck, arXiv:1704.08524

module mod_obj_1s0SYLV

  use mod_obj_1s0LY
  implicit none

  type, extends(obj_1s0LY) :: obj_1s0SYLV

  contains

    ! procedure   :: get_paras_filename   => get_paras_filename_1s0SYLV
    ! procedure   :: get_writeout_string  => get_writeout_string_1s0SYLV
    procedure   :: get_num_paras        => get_num_paras_1s0SYLV
    ! Frozen
    ! procedure   :: change_a_C           => change_a_C_1s0SYLV
    procedure   :: SetEdepCs            => SetEdepCs_1s0SYLV
    procedure   :: update_k             => update_k_1s0SYLV
    ! procedure   :: get_allTQn           => get_allTQn_1s0SYLV

  end type obj_1s0SYLV

contains

  ! subroutine get_paras_filename_1s0SYLV(self, fname)

  !   class(obj_1s0SYLV), intent(inout)   :: self
  !   character(len = *), intent(out)     :: fname

  !   select case (self%uptoQn)
  !     case (0)
  !       fname = 'run_3halves_Cs_Q0.in'
  !     ! case (1)
  !     !     fname = 'run_3halves_Cs_Q1.in'
  !     ! case (2)
  !     !     fname = 'run_3halves_Cs_Q2.in'
  !     ! case (3)
  !     !     fname = 'run_3halves_Cs_Q3.in'
  !     case default
  !       write (standard_error_unit, '(a, i2)')  &
  !         'get_paras_filename_1s0SYLV: I don''t know how to handle order ', self%uptoQn
  !       fname = ''
  !   end select

  ! end subroutine get_paras_filename_1s0SYLV

  ! subroutine get_writeout_string_1s0SYLV(self, outstr)

  !   class(obj_1s0SYLV), intent(inout)   :: self
  !   character(len = 10), intent(out)        :: outstr

  !   outstr = '3halves'

  ! end subroutine get_writeout_string_1s0SYLV

  subroutine get_num_paras_1s0SYLV(self, num_cnttrm)

    class(obj_1s0SYLV), intent(inout)   :: self
    integer, intent(out)                :: num_cnttrm

    select case (self%uptoQn)
      case (0)
        num_cnttrm = 3
      ! case (1)
      !     num_cnttrm = 5
      ! case (2)
      !     num_cnttrm = 9
      case default
        write (standard_error_unit, '(a, i2)')  &
          'get_num_paras_1s0SYLV: I don''t know how to handle order ', self%uptoQn
        num_cnttrm = 0
    end select

  end subroutine get_num_paras_1s0SYLV


  ! Frozen
  ! ! If self%Cs(1) or self%Cs(2) or self%Cs(3) is changed, VMVLO needs to be re-initialized.
  ! subroutine change_a_C_1s0SYLV(self, n, Cn)

  !     class(obj_1s0SYLV), intent(inout)   :: self
  !     integer, intent(in)                     :: n
  !     real(NER), intent(in)                   :: Cn

  !     call change_a_C_channel(self, n, Cn)
  !     if (n==1 .or. n==2 .or. n==3) self%VMVLO%inited = .false.

  ! end subroutine change_a_C_1s0SYLV

  subroutine SetEdepCs_1s0SYLV(self, E)

    class(obj_1s0SYLV), intent(inout)   :: self
    real(NER), intent(in)               :: E

    real(NER)   :: PC_CgA

    ! Delta = inparas(1), y = inparas(2), z = inparas(3)
    ! Cs(1) = y/(Ecm + Delta) + z - CgA
    PC_CgA = PC_gA*PC_gA*PC_mN/(8.0*PI_NE*PI_NE*PC_fpi*PC_fpi)
    self%Cs(1) = self%inparas(2)/(E + self%inparas(1)) - PC_CgA + self%inparas(3)
    if (self%uptoQn == 0) return

  end subroutine SetEdepCs_1s0SYLV

  ! Since the LO contact interaction is energy-dependent, update_k is called to
  ! initialize VMVLO.
  subroutine update_k_1s0SYLV(self, k)

    class(obj_1s0SYLV), intent(inout)  :: self
    real(NER), intent(in)               :: k

    call self%SetEdepCs(k2Ecm(k))
    self%VMVLO%inited = .false.
    call update_k_1s0LY(self, k)

    ! call update_k_sngl(self, k)
    ! call self%update_edge_VMpoly0(self%VM_poly0)
    ! call self%init_VMVLO(C0)
    ! call self%update_edge_VMVLO(k, C0)
    ! self%Tpq%updated = .false.

  end subroutine update_k_1s0SYLV


  ! subroutine get_allTQn_1s0SYLV(self)

  !     class(obj_1s0SYLV), intent(inout)   :: self

  !     integer         :: N
  !     ! real(NER)       :: VL2, C1, C2, D2, Ecm
  !     ! complex(NEC)    :: t2, L1, L2, CT, DT, xi0

  !     N = self%N
  !     call get_allTQn_sngl(self)

  !     if (self%uptoQn > 1) then
  !         call self%get_Tpq_from_VM(self%VMVLO, self%Tpq)
  !         self%Tk%ptr(1:N+1) = self%Tpq%ptr(1:N+1, N+1)
  !         self%Tk%updated = .true.
  !     else
  !         call self%get_Tk_from_VM(self%VMVLO, self%Tk)
  !     end if
  !     self%TQn(1) = self%Tk%ptr(N+1)
  !     if (self%uptoQn == 0) then
  !         self%TQn_updated = .true.
  !         return
  !     end if

  !     if (self%uptoQn > 0) then
  !         write (standard_error_unit, '(a, i2)')  &
  !                 'get_allTQn_1s0SYLV: don''t know how to handle order ', self%uptoQn
  !     end if

    ! if (self%uptoQn == 1) then
    !     call self%get_CT(self%Tk, CT)
    ! else
    !     call self%get_CDT(self%Tk, CT, DT)
    ! end if
    ! Ecm = k2Ecm(self%k)
    ! C1 = self%Cs(3)/(Ecm + self%Cs(1))**2 + self%Cs(4)/(Ecm + self%Cs(1)) + self%Cs(5)
    ! self%TQn(2) = C1*CT
    ! if (self%uptoQn == 1) then
    !     self%TQn_updated = .true.
    !     return
    ! end if

    ! C2 = self%Cs(3)*self%Cs(3)/(self%Cs(2)*(Ecm + self%Cs(1))**3) + self%Cs(6)/(Ecm + self%Cs(1))**2    &
    !     + self%Cs(7)/(Ecm + self%Cs(1)) + self%Cs(8)
    ! D2 = self%Cs(9)
    ! call self%get_C_xis(self%Tpq, C1, xi0)
    ! VL2 = self%vs_vl2%VM%ptr(N+1, N+1)
    ! call self%get_1loop_os_VMTk(self%vs_vl2%VM, self%Tk, L1)
    ! call self%get_2loops_TkVMTk(self%vs_vl2%VM, self%Tk, L2)
    ! self%TQn(3) = VL2 + 2.0*L1 + L2 + (xi0 + C2)*CT + D2*DT
    ! self%TQn_updated = .true.

  ! end subroutine get_allTQn_1s0SYLV


end module mod_obj_1s0SYLV
