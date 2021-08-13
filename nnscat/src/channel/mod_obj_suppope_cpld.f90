! Bingwei Long 10/03/2018
! Bingwei Long 08/22/2018
! Class obj_SUPpope_cpld (coupled channels where OPE is perturbative.)
! Power counting suggested by the referee

module mod_obj_suppope_cpld

  use mod_obj_pope_cpld
  implicit none

  type, extends(obj_pope_cpld)    :: obj_SUPpope_cpld

  contains

    procedure   :: get_num_paras        => get_num_paras_SUPpope_cpld
    procedure   :: get_othersTQn        => get_othersTQn_SUPpope_cpld

  end type obj_SUPpope_cpld

contains

  subroutine get_num_paras_SUPpope_cpld(self, num_cnttrm)

    class(obj_SUPpope_cpld), intent(inout)  :: self
    integer, intent(out)                    :: num_cnttrm

    integer :: L

    if (self%uptoQn > 4) then
      write (standard_error_unit, '(a, i2)')  &
        & 'get_num_paras_SUPpope_cpld: Don''t know how to handle order ',      &
        & self%uptoQn
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
          num_cnttrm = 2
        else
          num_cnttrm = 0
        end if
      case (4)
        if (L == 1) then
          num_cnttrm = 5
        else
          num_cnttrm = 1
        end if
    end select

  end subroutine get_num_paras_SUPpope_cpld

  subroutine get_othersTQn_SUPpope_cpld(self, restTQn)

    class(obj_SUPpope_cpld), intent(inout) :: self
    complex(NEC), intent(inout)         :: restTQn(:, :, :)

    real(NER)   :: vsab(2, 2), vlab(2, 2), VQamp(2, 2)
    complex(NEC):: tmpTQn_2d(2, 2, 1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER),    &
      & tmpTQn_1d(2, 2, 1:MAX_NUMBER_ORDER), two_ope_one_VQ2(2, 2),            &
      & two_VQ2(2, 2), one_ope_one_VQ3(2, 2)

    real(NEC), parameter    :: CHEB_A_PERT = -1.0E-5_NER, CHEB_B_PERT = -CHEB_A_PERT
    integer, parameter      :: CHEB_BASE_N = 5

    restTQn(1:2, 1:2, 1:self%uptoQn+1) = 0.0_NER
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')    &
      &   'get_othersTQn_SUPpope_cpld: k is not defined.'
      return
    end if
    if (self%uptoQn <= 1) return

    if (self%Lsj%L == 1) then

      call Vfunc_VQ2_SUPpope_cpld(self, self%k, self%k, VQamp)
      restTQn(1:2, 1:2, 3) = VQamp(1:2, 1:2)
      if (self%uptoQn == 2) return

      if (.not. self%vs_VQ2%VM%inited) then
        if (.not. self%vs_VQ2%VM%created) call self%create_vstruct(self%vs_VQ2)
        call self%init_symmtr_vstruct(self%vs_VQ2, Vfunc_VQ2_SUPpope_cpld)
        call self%update_symmtr_edge_vstruct(self%vs_VQ2)
      end if
      call self%GetDoublePertV(self%vs_dummy_lo%VM, 2, self%vs_vl0%VM, 1,      &
        & self%vs_VQ2%VM, self%VM_VTOT, tmpTQn_2d, CHEB_A_PERT, CHEB_B_PERT,   &
        & CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
      call Vfunc_VQ3_SUPpope_cpld(self, self%k, self%k, VQamp)
      restTQn(1:2, 1:2, 4) = VQamp(1:2, 1:2) + tmpTQn_2d(1:2, 1:2, 2, 2)
      if (self%uptoQn == 3) return

      two_ope_one_VQ2 = tmpTQn_2d(1:2, 1:2, 3, 2)
      call self%GetSinglePertV(self%vs_dummy_lo%VM, 2, self%vs_VQ2%VM,         &
        & self%VM_VTOT, tmpTQn_1d, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
      two_VQ2 = tmpTQn_1d(1:2, 1:2, 3)
      if (.not. self%vs_VQ3%VM%inited) then
          if (.not. self%vs_VQ3%VM%created)        &
            & call self%create_vstruct(self%vs_VQ3)
          call self%init_symmtr_vstruct(self%vs_VQ3, Vfunc_VQ3_SUPpope_cpld)
          call self%update_symmtr_edge_vstruct(self%vs_VQ3)
      end if
      call self%GetDoublePertV(self%vs_dummy_lo%VM, 1, self%vs_vl0%VM, 1,      &
        & self%vs_VQ3%VM, self%VM_VTOT, tmpTQn_2d, CHEB_A_PERT, CHEB_B_PERT,   &
        & CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
      one_ope_one_VQ3 = tmpTQn_2d(1:2, 1:2, 2, 2)
      call Vfunc_VQ4_SUPpope_cpld(self, self%k, self%k, VQamp)
      restTQn(1:2, 1:2, 5) = two_ope_one_VQ2 + two_VQ2 + one_ope_one_VQ3 + VQamp
      if (self%uptoQn == 4) return

    else

      if (self%uptoQn == 2) return

      call Vfunc_VQ3_SUPpope_cpld(self, self%k, self%k, VQamp)
      restTQn(1:2, 1:2, 4) = VQamp(1:2, 1:2)
      if (self%uptoQn == 3) return

      if (.not. self%vs_VQ3%VM%inited) then
          if (.not. self%vs_VQ3%VM%created)        &
            & call self%create_vstruct(self%vs_VQ3)
          call self%init_symmtr_vstruct(self%vs_VQ3, Vfunc_VQ3_SUPpope_cpld)
          call self%update_symmtr_edge_vstruct(self%vs_VQ3)
      end if
      call self%GetDoublePertV(self%vs_dummy_lo%VM, 1, self%vs_vl0%VM, 1,      &
        & self%vs_VQ3%VM, self%VM_VTOT, tmpTQn_2d, CHEB_A_PERT, CHEB_B_PERT,   &
        & CHEB_BASE_N, CHEB_A_PERT, CHEB_B_PERT, CHEB_BASE_N)
      one_ope_one_VQ3 = tmpTQn_2d(1:2, 1:2, 2, 2)
      call Vfunc_VQ4_SUPpope_cpld(self, self%k, self%k, VQamp)
      restTQn(1:2, 1:2, 5) = one_ope_one_VQ3 + VQamp
      if (self%uptoQn == 4) return
    end if

    write (standard_error_unit, '(a)')  &
      & 'get_othersTQn_SUPpope_cpld: don''t know to do uptoQn > 4'

  end subroutine get_othersTQn_SUPpope_cpld

  subroutine Vfunc_VQ2_SUPpope_cpld(self, p1, p2, vab)

    class(obj_cpld), intent(inout)                 :: self
    real(NER), intent(in)                       :: p1, p2
    real(NER), dimension(1:2, 1:2), intent(out) :: vab

    real(NER)   :: vsab(2, 2)

    if (self%lsj%L == 1) then
      call Vfunc_VS0_cpld(self, p1, p2, vsab)
      vab = self%Cs(1)*vsab
    else
      vab(1:2, 1:2) = 0.0_NER
    end if

  end subroutine Vfunc_VQ2_SUPpope_cpld

  subroutine Vfunc_VQ3_SUPpope_cpld(self, p1, p2, vab)

    class(obj_cpld), intent(inout) :: self
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out) :: vab(1:2, 1:2)

    real(NER)   :: vlab(2, 2), vs0ab(2, 2)

    call self%pVfunc_VL2(self, p1, p2, vlab)
    select case (self%lsj%L)
      case (1)
        call Vfunc_VS0_cpld(self, p1, p2, vs0ab)
        vab = vlab + self%Cs(2)*vs0ab
      case default
        vab = vlab
    end select

  end subroutine Vfunc_VQ3_SUPpope_cpld

  subroutine Vfunc_VQ4_SUPpope_cpld(self, p1, p2, vab)

    class(obj_cpld), intent(inout) :: self
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out) :: vab(1:2, 1:2)

    real(NER)   :: vlab(2, 2), vs0ab(2, 2), vsoffab(2, 2)

    call self%pVfunc_VL3(self, p1, p2, vlab)
    select case (self%lsj%L)
      case (1)
        call Vfunc_VS0_cpld(self, p1, p2, vs0ab)
        call Vfunc_VSoff_cpld(self, p1, p2, vsoffab)
        vab = vlab + (self%Cs(3) + self%Cs(4)*(p1*p1 + p2*p2))*vs0ab           &
          & + self%Cs(5)*vsoffab
      case (2)
        call Vfunc_VS0_cpld(self, p1, p2, vs0ab)
        vab = vlab + self%Cs(1)*vs0ab
      case default
        vab = vlab
    end select

  end subroutine Vfunc_VQ4_SUPpope_cpld

end module mod_obj_suppope_cpld
