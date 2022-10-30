! Bingwei Long  07/29/2018
! Implementation of Long & Yang 1S0 with Taylor expansion method

module mod_obj_tlrLY

  use mod_obj_withc0_sngl
  implicit none

  type, extends(obj_withc0_sngl)  :: obj_tlrLY

    type(VM_sngl)       :: VM_VTOT ! To be used by Taylor methods

  contains

    procedure   :: init                 => init_tlrLY
    procedure   :: release              => release_tlrLY
    procedure   :: get_num_paras        => get_num_paras_tlrLY
    procedure   :: update_k             => update_k_tlrLY
    procedure   :: get_allTQn           => get_allTQn_tlrLY

  end type obj_tlrLY

contains

  subroutine init_tlrLY(self, N, regtype, Mambda)

    class(obj_tlrLY), intent(inout)  :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    if (self%chnid /= CHNID_1S0_PP)       &
      write (standard_error_unit, '(a)')  &
        & 'init_tlrLY: Warning: chnid should be CHNID_1S0_PP.'
    call init_withc0_sngl(self, N, regtype, Mambda)
    if (self%uptoQn > 1) call self%create_TM(self%Tpq)
    if (self%uptoQn > 2) call self%create_VM(self%VM_VTOT)

  end subroutine init_tlrLY

  subroutine release_tlrLY(self)

    class(obj_tlrLY), intent(inout)  :: self

    if (associated(self%VM_VTOT%ptr)) call self%release_VM(self%VM_VTOT)
    if (self%Tpq%created) call self%release_TM(self%Tpq)
    call release_withc0_sngl(self)

  end subroutine release_tlrLY

  subroutine get_num_paras_tlrLY(self, num_cnttrm)

    class(obj_tlrLY), intent(inout)  :: self
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
      case (3)
        ! C0, C1, D0, C2, D1, E0, C3, D2, E1, F0
        num_cnttrm = 10
      case default
        write (standard_error_unit, '(a, i2)')  &
          & 'get_num_paras_tlrLY: I don''t know how to handle order ', self%uptoQn
        num_cnttrm = 0
    end select

  end subroutine get_num_paras_tlrLY

  subroutine update_k_tlrLY(self, k)

    class(obj_tlrLY), intent(inout) :: self
    real(NER), intent(in)           :: k

    call update_k_withc0_sngl(self, k)
    if (self%vstcVQ1%VM%inited)                           &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ1)
    if (self%vstcVQ2%VM%inited)                           &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ2)
    if (self%vstcVQ3%VM%inited)                           &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ3)
    self%Tk%updated = .false.

  end subroutine update_k_tlrLY

  ! Frozen  08/02/2018
  ! Complete Taylor expansion for Q1 to Q3
  ! subroutine get_allTQn_tlrLY(self)

  !   class(obj_tlrLY), intent(inout)   :: self

  !   integer :: N
  !   real(NER), parameter    :: A1 = -0.1E-1_NER, B1 = -A1,  &
  !     &                        A2 = -0.1E-2_NER, B2 = -A2,  &
  !     &                        A3 = -0.1E-4_NER, B3 = -A3
  !   integer, parameter      :: N1 = 40, N2 = 40, N3 = 50
  !   complex(NEC)            :: tmp_2d(1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER), &
  !     &                        tmp_VQ2(1:MAX_NUMBER_ORDER),                    &
  !     &                        tmp_VQ3(1:MAX_NUMBER_ORDER)

  !   if (.not. self%k_defined) then
  !     write (standard_error_unit, '(a)')  &
  !       & 'get_allTQn_tlrLY: k is not defined. Nothing will be done.'
  !     return
  !   end if
  !   N = self%N
  !   if (.not. self%VMVLO%inited) call self%init_VMVLO(self%Cs(1))

  !   if (self%uptoQn == 0) then
  !     call self%get_Tk_from_VM(self%VMVLO, self%Tk)
  !     self%TQn(1) = self%Tk%ptr(N+1)
  !     self%TQn_updated = .true.
  !     return
  !   end if

  !   if (.not. self%vstcVQ1%VM%inited) then
  !     if (.not. self%vstcVQ1%VM%created) call self%create_vstruct(self%vstcVQ1)
  !     call self%init_symmtr_vstruct(self%vstcVQ1, Vfunc_VQ1_tlrLY)
  !     call self%update_symmtr_edge_vstruct(self%vstcVQ1)
  !   end if
  !   call self%GetSinglePertV(self%VMVLO, self%uptoQn, self%vstcVQ1%VM,         &
  !     & self%VM_VTOT, self%TQn, A1, B1, N1)
  !   if (self%uptoQn == 1) then
  !     self%TQn_updated = .true.
  !     return
  !   end if

  !   if (.not. self%vstcVQ2%VM%inited) then
  !     if (.not. self%vstcVQ2%VM%created) call self%create_vstruct(self%vstcVQ2)
  !       call self%init_symmtr_vstruct(self%vstcVQ2, Vfunc_VQ2_tlrLY)
  !       call self%update_symmtr_edge_vstruct(self%vstcVQ2)
  !   end if
  !   call self%GetSinglePertV(self%VMVLO, int(self%uptoQn/2), self%vstcVQ2%VM,  &
  !     & self%VM_VTOT, tmp_VQ2, A2, B2, N2)
  !   self%TQn(3) = self%TQn(3) + tmp_VQ2(2)
  !   if (self%uptoQn == 2) then
  !     self%TQn_updated = .true.
  !     return
  !   end if

  !   if (.not. self%vstcVQ3%VM%inited) then
  !     if (.not. self%vstcVQ3%VM%created) call self%create_vstruct(self%vstcVQ3)
  !       call self%init_symmtr_vstruct(self%vstcVQ3, Vfunc_VQ3_tlrLY)
  !       call self%update_symmtr_edge_vstruct(self%vstcVQ3)
  !   end if
  !   call self%GetDoublePertV(self%VMVLO, 1, self%vstcVQ1%VM,                   &
  !     & 1, self%vstcVQ2%VM, self%VM_VTOT, tmp_2d, A1, B1, N1, A2,              &
  !     & B2, N2)
  !   call self%GetSinglePertV(self%VMVLO, self%uptoQn-2, self%vstcVQ3%VM,       &
  !     & self%VM_VTOT, tmp_VQ3, A3, B3, N3)
  !   ! debug
  !   ! self%TQn(4) = self%TQn(4) + tmp_2d(2, 2)
  !   self%TQn(4) = self%TQn(4) + tmp_2d(2, 2) + tmp_VQ3(2)
  !   if (self%uptoQn == 3) then
  !     self%TQn_updated = .true.
  !     return
  !   end if

  !   write (standard_error_unit, '(a, I1.1)')  &
  !     & 'get_allTQn_tlrLY: cannot do order ', self%uptoQn
  !   return

  ! end subroutine get_allTQn_tlrLY

  ! Hybrid method, using Taylor expansion to calculate multiple insertions of
  ! VQ1 and combination of VQ1 and VQ2, and using get_?loop_* subs to calculate
  ! single insertion of VQ3
  ! As of 08/01/2018 11:55pm, uptoQn = 3 can be done for Mambda <= 1600 without
  ! breaking unitarity
  subroutine get_allTQn_tlrLY(self)

    class(obj_tlrLY), intent(inout)   :: self

    integer :: N
    real(NER), parameter    :: A1 = -0.1E-1_NER, B1 = -A1,  &
      &                        A2 = -0.1E-2_NER, B2 = -A2
    integer, parameter      :: N1 = 10, N2 = 20
    complex(NEC)            :: tmp_2d(1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER)

    real(NER)       :: V3, VL2
    complex(NEC)    :: L1, L2, CT, DT, ET, xi0, xi2, xi4

    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        & 'get_allTQn_tlrLY: k is not defined. Nothing will be done.'
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

    if (self%uptoQn <= 2) then
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
      self%TQn(3) = VL2 + 2.0_NER*L1 + L2 + (xi0 + self%Cs(4))*CT              &
        & + (xi2 + self%Cs(5))*DT + (xi4 + self%Cs(6))*ET
      self%TQn_updated = .true.
      return
    end if

    if (.not. self%vstcVQ1%VM%inited) then
      if (.not. self%vstcVQ1%VM%created) call self%create_vstruct(self%vstcVQ1)
      call self%init_symmtr_vstruct(self%vstcVQ1, Vfunc_VQ1_tlrLY)
      call self%update_symmtr_edge_vstruct(self%vstcVQ1)
    end if
    if (.not. self%vstcVQ2%VM%inited) then
      if (.not. self%vstcVQ2%VM%created) call self%create_vstruct(self%vstcVQ2)
      call self%init_symmtr_vstruct(self%vstcVQ2, Vfunc_VQ2_tlrLY)
      call self%update_symmtr_edge_vstruct(self%vstcVQ2)
    end if
    if (.not. self%vstcVQ3%VM%inited) then
      if (.not. self%vstcVQ3%VM%created) call self%create_vstruct(self%vstcVQ3)
      call self%init_symmtr_vstruct(self%vstcVQ3, Vfunc_VQ3_tlrLY)
      call self%update_symmtr_edge_vstruct(self%vstcVQ3)
    end if

    call self%GetSinglePertV(self%VMVLO, self%uptoQn, self%vstcVQ1%VM,         &
      & self%VM_VTOT, self%TQn, A1, B1, N1)
    call self%GetDoublePertV(self%VMVLO, 1, self%vstcVQ1%VM,                   &
      & 1, self%vstcVQ2%VM, self%VM_VTOT, tmp_2d, A1, B1, N1, A2,              &
      & B2, N2)

    self%TQn(3) = self%TQn(3) + tmp_2d(1, 2)
    V3 = self%vstcVQ3%VM%ptr(N+1, N+1)
    call self%get_Tk_from_VM(self%VMVLO, self%Tk)
    call self%get_1loop_os_VMTk(self%vstcVQ3%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vstcVQ3%VM, self%Tk, L2)
    self%TQn(4) = self%TQn(4) + tmp_2d(2, 2) + V3 + 2.0_NER*L1 + L2
    if (self%uptoQn == 3) then
      self%TQn_updated = .true.
      return
    end if

    write (standard_error_unit, '(a, I1.1)')  &
      & 'get_allTQn_tlrLY: cannot do order ', self%uptoQn
    return

  end subroutine get_allTQn_tlrLY

  ! Potentials at subleading orders
  subroutine Vfunc_VQ1_tlrLY(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER) :: vs

    call Vfunc_VS0_sngl(osngl, p1, p2, vs)
    v = (osngl%Cs(2) + osngl%Cs(3)*(p1*p1 + p2*p2))*vs

  end subroutine Vfunc_VQ1_tlrLY

  subroutine Vfunc_VQ2_tlrLY(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER) :: vs, vl

    call Vfunc_VTPE0_sngl(osngl, p1, p2, vl)
    call Vfunc_VS0_sngl(osngl, p1, p2, vs)
    v = (osngl%Cs(4) + osngl%Cs(5)*(p1*p1 + p2*p2) + osngl%Cs(6)*p1*p1*p2*p2)*vs
    v = v + vl

  end subroutine Vfunc_VQ2_tlrLY

  subroutine Vfunc_VQ3_tlrLY(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER) :: vs, vl

    call Vfunc_VTPE1_sngl(osngl, p1, p2, vl)
    call Vfunc_VS0_sngl(osngl, p1, p2, vs)
    v = vs*(osngl%Cs(7) + osngl%Cs(9)*p1*p1*p2*p2                     &
      & + (osngl%Cs(8) + osngl%Cs(10)*p1*p1*p2*p2)*(p1*p1 + p2*p2))
    v = v + vl

  end subroutine Vfunc_VQ3_tlrLY

end module mod_obj_tlrLY
