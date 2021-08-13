

module mod_vfunc_rel_ope

  use rel_ope_pwd
  use mod_obj_sngl
  use mod_obj_cpld
  implicit none


contains



subroutine Vfunc_rel_only_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval



    if (osngl%lsj%S == 0) then
      potval =rel_OPE_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval =rel_OPE_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval =rel_OPE_jpp(osngl%lsj%j, p1, p2)
        else
          potval =rel_OPE_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)    &
    & potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_rel_only_sngl


  subroutine Vfunc_vtpe0_rel_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval
    
    real(NER)        :: vll

    if (osngl%lsj%S == 0) then
      potval =rel_OPE_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval =rel_OPE_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval =rel_OPE_jpp(osngl%lsj%j, p1, p2)
        else
          potval =rel_OPE_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)    &
    & potval = potval * osngl%regulator(p1, p2)

    call Vfunc_VTPE0_sngl(osngl, p1, p2, vll)
     potval=potval+vll

  end subroutine Vfunc_vtpe0_rel_sngl





 subroutine Vfunc_rel_only_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)

    vab = 0.0
    if (ocpld%lsj%j == 0 .or. ocpld%lsj%L /= ocpld%lsj%j-1) return
    vab(1, 1) =rel_OPE_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) =rel_OPE_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) =rel_OPE_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) =rel_OPE_jpp(ocpld%lsj%j, p1, p2)
   

    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
    & vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_rel_only_cpld


  subroutine Vfunc_vtpe0_rel_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)
    real(NER)        :: vll(1:2,1:2)
    vab = 0.0
    if (ocpld%lsj%j == 0 .or. ocpld%lsj%L /= ocpld%lsj%j-1) return
    vab(1, 1) =rel_OPE_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) =rel_OPE_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) =rel_OPE_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) =rel_OPE_jpp(ocpld%lsj%j, p1, p2)
   

    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
    & vab = vab * ocpld%regulator(p1, p2)

    call Vfunc_VTPE0_cpld(ocpld, p1, p2, vll)
     vab=vab+vll

  end subroutine Vfunc_vtpe0_rel_cpld



end module mod_vfunc_rel_ope
