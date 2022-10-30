! Bingwei Long 07/27/2013
! Implementation of the 1s0 dibaryon paper by Long, PRC 88, 014002 (2013)

module mod_obj_lbwdib

  use mod_obj_withc0_sngl
  use eft_potspwd
  implicit none

  type, extends(obj_withc0_sngl) :: obj_lbwdib

    type(VM_sngl)       :: VM_VTOT ! To be used by Taylor methods
    type(TM_sngl)       :: TMa, TMb

  contains

    procedure   :: init                 => init_lbwdib
    procedure   :: release              => release_lbwdib
    procedure   :: get_num_paras        => get_num_paras_lbwdib
    procedure   :: load_inputs          => load_inputs_lbwdib
    ! Frozen
    ! procedure   :: change_a_C           => change_a_C_lbwdib
    procedure   :: SetEdepCs            => SetEdepCs_lbwdib
    procedure   :: update_k             => update_k_lbwdib
    procedure   :: get_allTQn           => get_allTQn_lbwdib

  end type obj_lbwdib

  integer, parameter  :: npar_lbwdib(1:4) = (/2, 5, 9, 14/)

contains

  subroutine init_lbwdib(self, N, regtype, Mambda)

    class(obj_lbwdib), intent(inout)  :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    if (self%chnid /= CHNID_1S0_PP) &
      write (standard_error_unit, '(a)')  &
        'init_lbwdib: Warning: chnid should be CHNID_1S0_PP.'
    call init_withc0_sngl(self, N, regtype, Mambda)
    if (self%uptoQn > 1) call self%create_TM(self%Tpq)
    if (self%uptoQn > 2) call self%create_VM(self%VM_VTOT)
    call self%create_TM(self%TMa)
    call self%create_TM(self%TMb)

  end subroutine init_lbwdib

  subroutine release_lbwdib(self)

    class(obj_lbwdib), intent(inout)  :: self

    if (associated(self%VM_VTOT%ptr)) call self%release_VM(self%VM_VTOT)
    if (self%Tpq%created) call self%release_TM(self%Tpq)
    if (self%TMa%created) call self%release_TM(self%TMa)
    if (self%TMb%created) call self%release_TM(self%TMb)
    call release_withc0_sngl(self)

  end subroutine release_lbwdib

  subroutine get_num_paras_lbwdib(self, num_cnttrm)

    class(obj_lbwdib), intent(inout)  :: self
    integer, intent(out)              :: num_cnttrm

    if (self%uptoQn > 3) then
      write (standard_error_unit, '(a, i2)')  &
        & 'get_num_paras_lbwdib: I don''t know how to handle order ',   &
        & self%uptoQn
      num_cnttrm = 0
    else
      num_cnttrm = npar_lbwdib(self%uptoQn + 1)
    end if

  end subroutine get_num_paras_lbwdib

  subroutine load_inputs_lbwdib(self, num_inputs, inputs)

    class(obj_lbwdib), intent(inout)    :: self
    integer, intent(in)                 :: num_inputs
    real(NER), dimension(:), intent(in) :: inputs

    self%inparas(1:num_inputs) = inputs(1:num_inputs)
    ! ! Converting from MeV to GeV
    ! self%inparas(1) = inputs(1)*1.0E-3_NER
    ! if (self%uptoQn > 0) then
    !   self%inparas(3) = inputs(3)*1.0E-3_NER
    !   self%inparas(5) = inputs(5)*1.0E3_NER
    ! end if
    ! if (self%uptoQn > 1) then
    !   self%inparas(6) = inputs(6)*1.0E-3_NER
    !   self%inparas(8) = inputs(8)*1.0E3_NER
    !   self%inparas(9) = inputs(9)*1.0E9_NER
    ! end if
    ! if (self%uptoQn > 2) then
    !   self%inparas(10) = inputs(10)*1.0E-3_NER
    !   self%inparas(12) = inputs(12)*1.0E3_NER
    !   self%inparas(13) = inputs(13)*1.0E9_NER
    !   self%inparas(14) = inputs(14)*1.0E15_NER
    ! end if

    if (self%uptoQn >= 2) then
    ! Cs(4) = D0 = self%inparas(9)
      self%Cs(4) = self%inparas(9)
      if (self%uptoQn >= 3) then
        ! Cs(6) = D1, Cs(7) = E0
        self%Cs(6) = self%inparas(13)
        self%Cs(7) = self%inparas(14)
      end if
    end if
    self%inpar_set = .true.

  end subroutine load_inputs_lbwdib

  ! The naming scheme of Delta1, y1, etc. is not self-explanatory.
  subroutine SetEdepCs_lbwdib(self, Ecm)

    class(obj_lbwdib), intent(inout)    :: self
    real(NER), intent(in)               :: Ecm

    real(NER)   :: PC_CgA

    PC_CgA = PC_gA*PC_gA*PC_mN/(8.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi)
    ! C0 = y/(E + Delta) - PC_CgA
    ! Delta = inparas(1), y = inparas(2)
    self%Cs(1) = self%inparas(2)/(Ecm + self%inparas(1)) - PC_CgA
    if (self%uptoQn == 0) return

    ! Delta1 = inparas(3), y1 = inparas(4), C0 = inparas(5)
    ! Delta1/(Ecm + Delta)^2 + y1/(Ecm + Delta) + C0
    ! Note that: Delta /= Delta + Delta1 + Delta2 + ...
    self%Cs(2) = self%inparas(3)/(Ecm + self%inparas(1))**2   &
      & + self%inparas(4)/(Ecm + self%inparas(1))             &
      & + self%inparas(5)
    if (self%uptoQn == 1) return

    ! Delta2 = inparas(6), y2 = inparas(7), C1 = inparas(8), D0 = inparas(9)
    ! (y^{-1} Delta1^2)/(Ecm + Delta)^3 + Delta2/(Ecm + Delta)^2
    !   + y2/(Ecm + Delta) + C1
    self%Cs(3) = self%inparas(3)*self%inparas(3)          &
      & /(self%inparas(2)*(Ecm + self%inparas(1))**3)     &
      &  + self%inparas(6)/(Ecm + self%inparas(1))**2     &
      &  + self%inparas(7)/(Ecm + self%inparas(1))        &
      &  + self%inparas(8)
    if (self%uptoQn == 2) return

    ! Bookmark 0807/2018: Check the following formula
    self%Cs(5) = self%inparas(3)**3/self%inparas(2)**2    &
      &      /(Ecm + self%inparas(1))**4                  &
      & + self%inparas(3)/self%inparas(2)                 &
      &    *(2.0_NER*self%inparas(6) - self%inparas(4)    &
      &       /self%inparas(2)*self%inparas(3))           &
      &   /(Ecm + self%inparas(1))**3                     &
      & + self%inparas(10)/(Ecm + self%inparas(1))**2     &
      & + self%inparas(11)/(Ecm + self%inparas(1))        &
      & + self%inparas(12)

  end subroutine SetEdepCs_lbwdib

  ! Since the LO contact interaction is energy-dependent, update_k is called to
  ! initialize VMVLO.
  subroutine update_k_lbwdib(self, k)

    class(obj_lbwdib), intent(inout)  :: self
    real(NER), intent(in)             :: k

    real(NER) :: Ecm

    ! ! Converting from MeV to GeV
    ! self%k = k*1.0E-3_NER

    self%k = k

    Ecm = k2Ecm(self%k)
    call self%SetEdepCs(Ecm)
    self%VMVLO%inited = .false.
    self%Tpq%updated = .false.
    if (self%vstcVQ1%VM%inited)                             &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ1)
    if (self%vstcVQ2%VM%inited)                             &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ2)
    if (self%vstcVQ3%VM%inited)                             &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ3)
    call update_k_withc0_sngl(self, self%k)

  end subroutine update_k_lbwdib

  subroutine get_allTQn_lbwdib(self)

    class(obj_lbwdib), intent(inout)   :: self

    integer         :: N, ii, jj
    real(NER)       :: VL, V3
    complex(NEC)    :: L1, L2, CT, DT, ET, xi0, xi0_1C1D, xi2_1C1D,   &
      &                TC1triple, TC1D0, TC1C2, TC3D1E0, TVL3, TC1VL2

    real(NER), parameter    :: A1 = -0.1E-5_NER, B1 = -A1,            &
      &                        A2 = -0.1E-5_NER, B2 = -A2
    integer, parameter      :: N1 = 10, N2 = 10
    complex(NEC)  :: tmp_VSQ1VL2(1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER)

    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')  &
        & 'get_allTQn_lbwdib: k is not defined. Nothing will be done.'
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

    select case (self%uptoQn)
      case (1)
        call self%get_CT(self%Tk, CT)
      case (2)
        call self%get_CDT(self%Tk, CT, DT)
      case (3)
        call self%get_CDET(self%Tk, CT, DT, ET)
    end select
    self%TQn(2) = self%Cs(2)*CT
    if (self%uptoQn == 1) then
      self%TQn_updated = .true.
      return
    end if

    call self%get_2Cxi(self%Tpq, xi0)
    VL = self%vs_vl2%VM%ptr(N+1, N+1)
    call self%get_1loop_os_VMTk(self%vs_vl2%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl2%VM, self%Tk, L2)

    self%TQn(3) = VL + 2.0_NER*L1 + L2                      &
      & + (self%Cs(2)*self%Cs(2)*xi0 + self%Cs(3))*CT       &
      & + self%Cs(4)*DT
    if (self%uptoQn == 2) then
      self%TQn_updated = .true.
      return
    end if

    TC1triple = self%Cs(2)**3*xi0*xi0*CT
    TC1C2 = 2.0_NER*self%Cs(2)*self%Cs(3)*xi0*CT
    call self%get_1C1Dxi(self%Tpq, xi0_1C1D, xi2_1C1D)

    TC1D0 = self%Cs(2)*self%Cs(4)*(xi0_1C1D*CT + xi2_1C1D*DT)
    TC3D1E0 = self%Cs(5)*CT + self%Cs(6)*DT + self%Cs(7)*ET
    call self%get_1loop_os_VMTk(self%vs_vl3%VM, self%Tk, L1)
    call self%get_2loops_TkVMTk(self%vs_vl3%VM, self%Tk, L2)

    TVL3 = self%vs_vl3%VM%ptr(N+1, N+1) + 2.0_NER*L1 + L2

    if (.not. self%vstcVQ1%VM%inited) then
      if (.not. self%vstcVQ1%VM%created) call self%create_vstruct(self%vstcVQ1)
      call self%init_symmtr_vstruct(self%vstcVQ1, Vfunc_VSQ1_lbwdib)
      call self%update_symmtr_edge_vstruct(self%vstcVQ1)
    end if
    call self%GetDoublePertV(self%VMVLO, 1, self%vstcVQ1%VM,                   &
      & 1, self%vs_vl2%VM, self%VM_VTOT, tmp_VSQ1VL2, A1, B1, N1, A2,          &
      & B2, N2)
    self%Tk%updated = .false.
    TC1VL2 = self%Cs(2)*tmp_VSQ1VL2(2, 2)

    self%TQn(4) = TC1triple + TC1C2 + TC1D0 + TC3D1E0 + TVL3 + TC1VL2
    if (self%uptoQn == 3) then
      self%TQn_updated = .true.
      return
    end if

    write (standard_error_unit, '(a, I1.1)')  &
      & 'get_allTQn_lbwdib: cannot do order ', self%uptoQn
    return

  end subroutine get_allTQn_lbwdib

  ! Potentials at subleading orders
  subroutine Vfunc_VSQ1_lbwdib(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER) :: vs

    call Vfunc_VS0_sngl(osngl, p1, p2, vs)
    v = vs

  end subroutine Vfunc_VSQ1_lbwdib

  subroutine Vfunc_VSQ2_lbwdib(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER) :: vs, vl

    call Vfunc_VS0_sngl(osngl, p1, p2, vs)
    v = (p1*p1 + p2*p2)*vs

  end subroutine Vfunc_VSQ2_lbwdib

  subroutine Vfunc_VSQ3_lbwdib(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER) :: vs, vl

    call Vfunc_VS0_sngl(osngl, p1, p2, vs)
    v = vs*p1*p1*p2*p2

  end subroutine Vfunc_VSQ3_lbwdib

end module mod_obj_lbwdib
