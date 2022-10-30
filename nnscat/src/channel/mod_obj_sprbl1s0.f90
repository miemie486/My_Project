! Bingwei Long 08/14/2018 11:57pm

module mod_obj_sprbl1s0

  use mod_obj_lbwdib
  implicit none

  ! Uncoupled channels where C starts to appear at O(1), D at O(Q^2), etc.
  type, extends(obj_lbwdib)    :: obj_sprbl1s0

    ! type(vstruct_sngl)  :: vstc_VS0

  contains

    procedure   :: init                 => init_sprbl1s0
    procedure   :: get_num_paras        => get_num_paras_sprbl1s0
    procedure   :: load_inputs          => load_inputs_sprbl1s0
    procedure   :: update_k             => update_k_sprbl1s0

  end type obj_sprbl1s0

  integer, parameter  :: npar_sprbl1s0(1:4) = (/2, 5, 9, 14/)

contains

  subroutine init_sprbl1s0(self, N, regtype, Mambda)

    class(obj_sprbl1s0), intent(inout)  :: self
    integer, intent(in)                 :: N, regtype
    real(NER), intent(in)               :: Mambda

    real(NER) :: PC_CgA

    self%pVfunc_VL0 => Vfunc_VL0_sprbl1s0
    PC_CgA = PC_gA*PC_gA*PC_mN/(8.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi)
    self%Cs(1) = -PC_CgA
    call init_lbwdib(self, N, regtype, Mambda)

  end subroutine init_sprbl1s0

  subroutine get_num_paras_sprbl1s0(self, num_cnttrm)

    class(obj_sprbl1s0), intent(inout)  :: self
    integer, intent(out)              :: num_cnttrm

    if (self%uptoQn > 3) then
      write (standard_error_unit, '(a, i2)')  &
        & 'get_num_paras_sprbl1s0: I don''t know how to handle order ',   &
        & self%uptoQn
      num_cnttrm = 0
    else
      num_cnttrm = npar_sprbl1s0(self%uptoQn + 1)
    end if

  end subroutine get_num_paras_sprbl1s0

  subroutine load_inputs_sprbl1s0(self, num_inputs, inputs)

    class(obj_sprbl1s0), intent(inout)    :: self
    integer, intent(in)                 :: num_inputs
    real(NER), dimension(:), intent(in) :: inputs

    self%inparas(1:num_inputs) = inputs(1:num_inputs)

    if (self%uptoQn >= 1) then
      ! C1
      self%Cs(2) = self%inparas(5)
    end if

    if (self%uptoQn >= 2) then
      ! C2
      self%Cs(3) = self%inparas(8)
      ! Cs(4) = D0 = self%inparas(9)
      self%Cs(4) = self%inparas(9)
      if (self%uptoQn >= 3) then
        ! C3
        self%Cs(5) = self%inparas(12)
        ! Cs(6) = D1, Cs(7) = E0
        self%Cs(6) = self%inparas(13)
        self%Cs(7) = self%inparas(14)
      end if
    end if
    self%inpar_set = .true.

  end subroutine load_inputs_sprbl1s0

  subroutine update_k_sprbl1s0(self, k)

    class(obj_sprbl1s0), intent(inout)  :: self
    real(NER), intent(in)             :: k

    self%k = k
    self%Tpq%updated = .false.
    if (self%vstcVQ1%VM%inited)                             &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ1)
    if (self%vstcVQ2%VM%inited)                             &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ2)
    if (self%vstcVQ3%VM%inited)                             &
      & call self%update_symmtr_edge_vstruct(self%vstcVQ3)
    call update_k_withc0_sngl(self, self%k)

  end subroutine update_k_sprbl1s0

  subroutine Vfunc_VL0_sprbl1s0(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: potval

    real(NER) :: kps, vope

    kps = PC_mN*osngl%inparas(1)
    potval = PC_mN*osngl%inparas(2)/(sqrt(p1*p1 + kps)*sqrt(p2*p2 + kps))
    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)    &
      & potval = potval * osngl%regulator(p1, p2)

    call Vfunc_OPE_sngl(osngl, p1, p2, vope)
    potval = potval + vope

  end subroutine Vfunc_VL0_sprbl1s0


end module mod_obj_sprbl1s0
