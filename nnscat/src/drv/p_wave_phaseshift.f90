!this module provide phase shift  used in fitting programs

module p_wave_phaseshift

  use mod_obj_suppope_sngl
  use mod_obj_suppope_cpld
  use mod_vfunc

  implicit none


contains

  !sub for 1p1 3p0 3p1 1d2 3d1 3d2 using deltafull potentials

  subroutine GetPS_sngl_delta(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_VDelTPE0sfr_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VDelTPE1sfr_sngl

    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
  !debug
!       print *, ptrchn%inparas(1), ptrchn%phsn(2), ptrchn%phsn(3), ptrchn%phsn(4)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_sngl_delta

    !sub for 1p1 3p0 3p1 1d2 3d1 using deltaless potentials

  subroutine GetPS_sngl_deltaless(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_sfrTPE0_sngl
      ptrchn%pVfunc_VL3 => Vfunc_sfrTPE1_sngl

    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_sngl_deltaless

! sub for 3p2 3d3 using deltafull potentials

  subroutine GetPS_cpld_delta(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld),  intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_VDelTPE0sfr_cpld
      ptrchn%pVfunc_VL3 => Vfunc_VDelTPE1sfr_cpld

    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_cpld_delta

  ! sub for 3p2 3d3 using deltaless potentials

  subroutine GetPS_cpld_deltaless(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld),  intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_sfrTPE0_cpld
      ptrchn%pVfunc_VL3 => Vfunc_sfrTPE1_cpld

    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_cpld_deltaless

  subroutine Get_3p0_wotpe_delta(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => V0sngl
      ptrchn%pVfunc_VL3 => V0sngl

    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Get_3p0_wotpe_delta

  subroutine GetPS_1p1_delta(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_VDelTPE0sfr_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VDelTPE1sfr_sngl

    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_1p1_delta



end module p_wave_phaseshift
