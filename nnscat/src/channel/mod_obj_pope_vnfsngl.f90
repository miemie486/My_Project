! Bingwei Long  07/17/2018

! Template to calculate uncoupled-channel phase shifts using perturbative
! scheme. Power counting is similar but slightly different from arXiv:1807.04407
! (Shaowei Wu and Bingwei Long )
module mod_obj_pope_vnfsngl

  use mod_obj_pope_sngl
  use mod_vfunc
  implicit none

  type, extends(obj_pope_sngl)    :: obj_pope_vnfsngl

  contains

    procedure   :: get_num_paras        => get_num_paras_pope_vnfsngl
    procedure   :: get_othersTQn        => get_othersTQn_pope_vnfsngl

  end type obj_pope_vnfsngl

contains

  subroutine Vfunc_VQ2_pope_vnfsngl(osngl, p1, p2, v)

    class(obj_sngl), intent(inout) :: osngl
    real(NER), intent(in)       :: p1, p2
    real(NER), intent(out)      :: v

    real(NER)   :: vl, vs

    ! Here is where Vfunc_VNF_sngl enters the calculation
    call Vfunc_VNF_sngl(osngl, p1, p2, vl)
    v = vl

  end subroutine Vfunc_VQ2_pope_vnfsngl

  ! Return the number of contact terms, i.e., the number of V_CT coupling
  ! constants (cf. arXiv:1807.04407), for the order given by self%uptoQn
  ! Do *not* change the following code unless you know what you are doing.
  subroutine get_num_paras_pope_vnfsngl(self, num_cnttrm)

    class(obj_pope_vnfsngl), intent(inout)  :: self
    integer, intent(out)                    :: num_cnttrm

    integer :: L

    num_cnttrm = 0

  end subroutine get_num_paras_pope_vnfsngl

  subroutine get_othersTQn_pope_vnfsngl(self, restTQn)

    class(obj_pope_vnfsngl), intent(inout) :: self
    complex(NEC), intent(inout)         :: restTQn(:)

    ! real(NEC), parameter    :: CHEB_A_PERT = -0.05_NER, CHEB_B_PERT = 0.05_NER
    ! integer, parameter      :: CHEB_BASE_N = 2
    real(NER)   :: VQamp
    complex(NEC):: loopval, tmpTQn(1:MAX_NUMBER_ORDER, 1:MAX_NUMBER_ORDER)

    restTQn(1:self%uptoQn+1) = 0.0_NER
    if (.not. self%k_defined) then
      write (standard_error_unit, '(a)')                &
      & 'get_othersTQn_pope_vnfsngl: k is not defined.'
      return
    end if
    if (self%uptoQn <= 1) return
    call Vfunc_VQ2_pope_vnfsngl(self, self%k, self%k, VQamp)
    restTQn(3) = VQamp
    if (self%uptoQn == 2) return

    write (standard_error_unit, '(a)')                  &
    & 'get_othersTQn_pope_vnfsngl: don''t know to do uptoQn > 2'

  end subroutine get_othersTQn_pope_vnfsngl

end module mod_obj_pope_vnfsngl
