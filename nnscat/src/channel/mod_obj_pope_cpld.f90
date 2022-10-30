! Bingwei Long 10/03/2018
! Bingwei Long 03/30/2018
! Unsuppressed scenario for coupled channel

module mod_obj_pope_cpld

  use mod_obj_cpld
  implicit none

  type, extends(obj_cpld)    :: obj_pope_cpld

    ! vs_dummy_lo is assumed to be set zero when created
    type(vstruct_cpld)  :: vs_dummy_lo, vs_VQ2, vs_VQ3, vs_VQ4

    ! Working memory to be used by Taylor expansion
    type(VM_cpld)             :: VM_VTOT


  contains

    procedure   :: init                 => init_pope_cpld
    procedure   :: release              => release_pope_cpld
    procedure   :: update_k             => update_k_pope_cpld
    procedure   :: get_num_paras        => get_num_paras_pope_cpld
    procedure   :: get_allTQn           => get_allTQn_pope_cpld
    procedure   :: get_opeladder_pope   => get_opeladder_pope_cpld
    procedure   :: get_othersTQn        => get_othersTQn_pope_cpld
    procedure   :: convert_TQn_to_phase_quiet   &
      & => convert_TQn_to_phase_quiet_pope_cpld

  end type obj_pope_cpld

contains

  subroutine init_pope_cpld(self, N, regtype, Mambda)

    class(obj_pope_cpld), intent(inout) :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    call init_channel(self, N, regtype, Mambda)
    call self%create_Tk(self%Tk)
    call self%create_VM(self%VM_VTOT)
    call self%create_vstruct(self%vs_vl0)
    call self%init_symmtr_vstruct(self%vs_vl0, self%pVfunc_VL0)
    call self%create_vstruct(self%vs_dummy_lo)

  end subroutine init_pope_cpld

  subroutine release_pope_cpld(self)

    class(obj_pope_cpld), intent(inout) :: self

    if (associated(self%vs_VQ2%VM%ptr)) call self%release_vstruct(self%vs_VQ2)
    ! call release_pope_cpld(self)
    call self%release_vstruct(self%vs_vl0)
    call self%release_vstruct(self%vs_dummy_lo)
    if (associated(self%vs_VQ3%VM%ptr)) call self%release_vstruct(self%vs_VQ3)
    if (associated(self%vs_VQ4%VM%ptr)) call self%release_vstruct(self%vs_VQ4)
    call self%release_Tk(self%Tk)
    call self%release_VM(self%VM_VTOT)
    call release_channel(self)

  end subroutine release_pope_cpld

  subroutine update_k_pope_cpld(self, k)

    class(obj_pope_cpld), intent(inout) :: self
    real(NER), intent(in)               :: k

    call update_k_channel(self, k)
    call self%update_symmtr_edge_vstruct(self%vs_vl0)
    if (self%vs_VQ3%VM%inited) call self%update_symmtr_edge_vstruct(self%vs_VQ3)
    if (self%vs_VQ4%VM%inited) call self%update_symmtr_edge_vstruct(self%vs_VQ4)
    self%Tk%updated = .false.
    ! call update_k_pope_cpld(self, k)
    if (self%vs_VQ2%VM%inited) call self%update_symmtr_edge_vstruct(self%vs_VQ2)

  end subroutine update_k_pope_cpld

  subroutine get_num_paras_pope_cpld(self, num_cnttrm)

    class(obj_pope_cpld), intent(inout)  :: self
    integer, intent(out)                    :: num_cnttrm

    integer :: L

    if (self%uptoQn > 3) then
      write (standard_error_unit, '(a, i2)')  &
        & 'get_num_paras_pope_cpld: Don''t know how to handle order ', self%uptoQn
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
          num_cnttrm = 4
        else
          num_cnttrm = 1
        end if
    end select

  end subroutine get_num_paras_pope_cpld

  subroutine get_allTQn_pope_cpld(self)

    class(obj_pope_cpld), intent(inout) :: self

    complex(NEC)  :: tmpTQn(2, 2, MAX_NUMBER_ORDER),    &
      & TQn_ladder(2, 2, MAX_NUMBER_ORDER)

    call self%get_opeladder_pope(TQn_ladder)
    call self%get_othersTQn(tmpTQn)
    self%TQn(1:2, 1:2, 1:self%uptoQn+1) = TQn_ladder(1:2, 1:2, 1:self%uptoQn+1) &
      & + tmpTQn(1:2, 1:2, 1:self%uptoQn+1)
    self%TQn_updated = .true.

  end subroutine get_allTQn_pope_cpld

  subroutine get_opeladder_pope_cpld(self, rsltTQn)

    class(obj_pope_cpld), intent(inout) :: self
    complex(NEC), intent(inout)         :: rsltTQn(:, :, :)

    real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER, CHEB_B_PERT = 0.05_NER
    ! real(NEC), parameter    :: CHEB_A_PERT = -0.1_NER, CHEB_B_PERT = 0.1_NER
    integer, parameter      :: CHEB_BASE_N = 2

    call get_allTQn_channel(self)
    self%TQn(1:2, 1:2, 1) = 0.0_NER
    if (self%uptoQn == 0) return
    ! if (self%uptoQn == 1) then
    !   self%TQn(1:2, 1:2, 2) = self%vs_vl0%VM%ptr(1:2, 1:2, self%N+1, self%N+1)
    !   return
    ! end if
    call self%GetSinglePertV(self%vs_dummy_lo%VM, self%uptoQn, self%vs_vl0%VM, &
      & self%VM_VTOT, rsltTQn, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
    self%TQn(1:2, 1:2, 1) = 0.0_NER

  end subroutine get_opeladder_pope_cpld

  subroutine get_othersTQn_pope_cpld(self, restTQn)

    class(obj_pope_cpld), intent(inout) :: self
    complex(NEC), intent(inout)         :: restTQn(:, :, :)

    real(NER)   :: vsab(2, 2), vlab(2, 2), VQamp(2, 2)
    complex(NEC):: tmpTQn(2, 2, 1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER),       &
    &   tmpVQ2VQ2(2, 2, 1:MAX_NUMBER_ORDER), loopval(2, 2)

    real(NEC), parameter    :: CHEB_A_PERT = -1.0E-6_NER, CHEB_B_PERT = -CHEB_A_PERT
    ! real(NEC), parameter    :: CHEB_A_PERT = -0.1_NER, CHEB_B_PERT = 0.1_NER
    integer, parameter      :: CHEB_BASE_N = 5

    restTQn(1:2, 1:2, 1:self%uptoQn+1) = 0.0_NER
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')    &
        & 'get_othersTQn_pope_cpld: k is not defined.'
      return
    end if
    if (self%uptoQn <= 1) return
    call Vfunc_VQ2_pope_cpld(self, self%k, self%k, VQamp)
    restTQn(1:2, 1:2, 3) = VQamp
    if (self%uptoQn == 2) return

    if (.not. self%vs_VQ2%VM%inited) then
      if (.not. self%vs_VQ2%VM%created) call self%create_vstruct(self%vs_VQ2)
      call self%init_symmtr_vstruct(self%vs_VQ2, Vfunc_VQ2_pope_cpld)
      call self%update_symmtr_edge_vstruct(self%vs_VQ2)
    end if
    call Vfunc_VQ3_pope_cpld(self, self%k, self%k, VQamp)
    call self%GetDoublePertV(self%vs_dummy_lo%VM, 1, self%vs_vl0%VM, 1,        &
      & self%vs_VQ2%VM, self%VM_VTOT, tmpTQn, CHEB_A_PERT, CHEB_B_PERT,        &
      & CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)

    restTQn(1:2, 1:2, 4) = VQamp + tmpTQn(1:2, 1:2, 2, 2)
    if (self%uptoQn == 3) return

    ! if (.not. self%vs_VQ3%VM%inited) then
    !     if (.not. self%vs_VQ3%VM%created) call self%create_vstruct(self%vs_VQ3)
    !     call self%init_symmtr_vstruct(self%vs_VQ3, Vfunc_VQ3_pope_cpld)
    !     call self%update_symmtr_edge_vstruct(self%vs_VQ3)
    ! end if
    ! restTQn(1:2, 1:2, 5) = tmpTQn(1:2, 1:2, 3, 2) + tmpVQ2VQ2(1:2, 1:2, 3) + VQamp
    ! if (self%uptoQn == 4) return

    write (standard_error_unit, '(a)')  &
    & 'get_othersTQn_pope_cpld: don''t know how to do uptoQn > 3'

  end subroutine get_othersTQn_pope_cpld

  subroutine convert_TQn_to_phase_quiet_pope_cpld(self)

    class(obj_pope_cpld), intent(inout)      :: self

    integer               :: ii
    complex(NEC)          :: violation(2, 2)
    character(len = 512)  :: msg
    logical               :: is_broken

    if (.not. self%TQn_updated) then
      write (standard_error_unit, '(a)')  &
      & 'convert_TQn_to_phase_quiet_pope_cpld: Warning: TQn is not updated.'
      return
    end if

    self%phsn(1)%d1 = 0.0_NER
    self%phsn(1)%d2 = 0.0_NER
    self%phsn(1)%e = 0.0_NER
    self%uflag(1) = UNITARITY_SAFE
    if (self%uptoqn == 0) return
    call get_nth_delta_cpld(1, self%TQn(1:2, 1:2, 2), self%phsn(1:1), self%k,  &
      & self%uflag(2), violation, self%phsn(2))
    if (self%chnid .eq. CHNID_3P2_PP) then
      if (self%phsn(2)%d1 < 0.0_NER) then
        self%phsn(2)%d1 = self%phsn(2)%d1 + 180.0_NER
        self%phsn(2)%e = -self%phsn(2)%e
      end if
    end if
    do ii = 2, self%uptoQn
      call get_nth_delta_cpld(ii, self%TQn(1:2, 1:2, ii+1), self%phsn(1:ii),   &
        & self%k, self%uflag(ii+1), violation, self%phsn(ii+1))
    end do
    self%phsn_updated = .true.

    call convert_TQn_to_phase_quiet_channel(self)

  end subroutine convert_TQn_to_phase_quiet_pope_cpld

  subroutine Vfunc_VQ2_pope_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)              :: ocpld
    real(NER), intent(in)                       :: p1, p2
    real(NER), dimension(1:2, 1:2), intent(out) :: vab

    real(NER)   :: vlab(2, 2), vsab(2, 2)

    call ocpld%pVfunc_VL2(ocpld, p1, p2, vlab)
    if (ocpld%lsj%L == 1) then
      call Vfunc_VS0_cpld(ocpld, p1, p2, vsab)
      vab = vlab + ocpld%Cs(1)*vsab
    else
      vab = vlab
    end if

  end subroutine Vfunc_VQ2_pope_cpld

  subroutine Vfunc_VQ3_pope_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout) :: ocpld
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out) :: vab(1:2, 1:2)

    real(NER)   :: vlab(2, 2), vs0ab(2, 2), vsoffab(2, 2)

    call ocpld%pVfunc_VL3(ocpld, p1, p2, vlab)
    select case (ocpld%lsj%L)
      case (1)
        call Vfunc_VS0_cpld(ocpld, p1, p2, vs0ab)
        call Vfunc_VSoff_cpld(ocpld, p1, p2, vsoffab)
        vab = vlab + (ocpld%Cs(2) + ocpld%Cs(3)*(p1*p1 + p2*p2))*vs0ab      &
          & + ocpld%Cs(4)*vsoffab
      case (2)
        call Vfunc_VS0_cpld(ocpld, p1, p2, vs0ab)
        vab = vlab + ocpld%Cs(1)*vs0ab
      case default
        vab = vlab
    end select

  end subroutine Vfunc_VQ3_pope_cpld

  ! subroutine Vfunc_VQ4_pope_cpld(self, p1, p2, vab)

  !   class(obj_cpld), intent(inout)                 :: self
  !   real(NER), intent(in)                       :: p1, p2
  !   real(NER), dimension(1:2, 1:2), intent(out) :: vab

  !   real(NER)   :: vlab(2, 2), vsab(2, 2)

  !   vlab(1:2, 1:2) = 0.0_NER
  !   if (self%lsj%L == 1) then
  !     call Vfunc_VS0_cpld(self, p1, p2, vsab)
  !     vab = vlab + (self%Cs(3) + self%Cs(4)*(p1*p1 + p2*p2))*vsab
  !     call Vfunc_VSoff_cpld(self, p1, p2, vsab)
  !     vab = vab + self%Cs(5)*vsab
  !   else
  !     vab = vlab
  !   end if

  ! end subroutine Vfunc_VQ4_pope_cpld

end module mod_obj_pope_cpld
