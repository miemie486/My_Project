module pw_phase

  use mod_obj_suppope_sngl
  use mod_obj_suppope_cpld
  use mod_obj_sngl
  use mod_obj_cpld
  use mod_vfunc_rel_ope
  implicit none
  

contains
  subroutine GetPS_1p1_rel(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    
    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_vtpe0_rel_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VTPE1_sngl
  	
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_1p1_rel

  subroutine GetPS_3p1_rel(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    
    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_vtpe0_rel_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VTPE1_sngl
  	
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_3p1_rel

  subroutine GetPS_3p2_rel(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld),  intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    
    allocate(obj_SUPpope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_vtpe0_rel_cpld
      ptrchn%pVfunc_VL3 => Vfunc_VTPE1_cpld
  	
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_3p2_rel


  subroutine GetPS_1d2_rel(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    
    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_vtpe0_rel_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VTPE1_sngl
    
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_1d2_rel


  subroutine GetPS_3d2_rel(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_sngl), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    
    allocate(obj_SUPpope_sngl::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_vtpe0_rel_sngl
      ptrchn%pVfunc_VL3 => Vfunc_VTPE1_sngl
    
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_3d2_rel


    subroutine GetPS_3d3_rel(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)
    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld),  intent(out)  :: phaselst(:, :)

    integer                 :: regtype
    class(obj_SUPpope_cpld), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    
    allocate(obj_SUPpope_cpld::ptrchn)
    call create_channel(ptrchn, chnid, uptoQn)
      ptrchn%pVfunc_VL2 => Vfunc_vtpe0_rel_cpld
      ptrchn%pVfunc_VL3 => Vfunc_VTPE1_cpld
    
    call ptrchn%get_allorders_phase_for_klist(Nmesh, regtype, Mambda,   &
      & klist, numk, paralst, num_paras, phaselst)
    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine GetPS_3d3_rel

end module pw_phase
