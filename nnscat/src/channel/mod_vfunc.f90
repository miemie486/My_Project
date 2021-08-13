! Bingwei Long 07/14/2018

! More Vfunc_sngl's and Vfunc_cpld's that are *not* already included in
! mod_obj_sngl or mod_obj_cpld

module mod_vfunc

  use mod_deltafullNLO_pwd,  only:TPEdel_jpm,TPEdel_jmm,TPEdel_jpp,TPEdel_j0j,TPEdel_j1j
  use mod_deltafullN2LO_pwd,  only:TPEdel_sub_jpm,TPEdel_sub_jmm,TPEdel_sub_jpp,TPEdel_sub_j0j,TPEdel_sub_j1j
  use mod_obj_sngl
  use mod_obj_cpld
  use VTPE_SFR_pwd, only:TPESFR_j0j,TPESFR_j1j,TPESFR_jpm,TPESFR_jpp,TPESFR_jmm
  use VTPE1_SFR_pwd, only:TPESFR1_j0j,TPESFR1_j1j,TPESFR1_jpm,TPESFR1_jpp,TPESFR1_jmm

  implicit none

contains

  ! Pot. subs. for obj_cpld and obj_sngl
  ! Unlike VNF_j0j, VNF_j1j, etc. the following subs are coupled with obj_cpld and
  ! obj_sngl, so they can access information like regulator, Lsj, and so on.

  ! < p1 | VNF | p2 > for coupled channels
  subroutine Vfunc_VDelTPE0_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)
    real(NER)                       :: val_deltaless(1:2,1:2)

    vab = 0.0
    call Vfunc_VDelTPE0_only_cpld(ocpld, p1, p2, vab)
    call Vfunc_VTPE0_cpld(ocpld, p1, p2,val_deltaless)

    vab = vab + val_deltaless

  end subroutine Vfunc_VDelTPE0_cpld

  subroutine Vfunc_VDelTPE1_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)
    real(NER)                       :: val_deltaless(1:2,1:2)

    vab = 0.0
    call Vfunc_VDelTPE1_only_cpld(ocpld, p1, p2, vab)
    call Vfunc_VTPE1_cpld(ocpld, p1, p2,val_deltaless)

    vab = vab + val_deltaless


  end subroutine Vfunc_VDelTPE1_cpld

! SFR cpld potential

  subroutine Vfunc_sfrTPE0_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)

    vab = 0.0
    if (ocpld%lsj%j == 0 .or. ocpld%lsj%L /= ocpld%lsj%j-1) return
    vab(1, 1) = TPESFR_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) = TPESFR_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) = TPESFR_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) = TPESFR_jpp(ocpld%lsj%j, p1, p2)
    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
    & vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_sfrTPE0_cpld

  subroutine Vfunc_sfrTPE1_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)

    vab = 0.0
    if (ocpld%lsj%j == 0 .or. ocpld%lsj%L /= ocpld%lsj%j-1) return
    vab(1, 1) = TPESFR1_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) = TPESFR1_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) = TPESFR1_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) = TPESFR1_jpp(ocpld%lsj%j, p1, p2)
    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
    & vab = vab * ocpld%regulator(p1, p2)

  end subroutine Vfunc_sfrTPE1_cpld

! Vsfr deltaless + deltafull cpld

  subroutine Vfunc_VDelTPE0sfr_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)
    real(NER)                       :: val_deltaless(1:2,1:2)

    vab = 0.0
    call Vfunc_VDelTPE0_only_cpld(ocpld, p1, p2, vab)
    call Vfunc_sfrTPE0_cpld(ocpld, p1, p2,val_deltaless)

    vab = vab + val_deltaless

  end subroutine Vfunc_VDelTPE0sfr_cpld

  subroutine Vfunc_VDelTPE1sfr_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)
    real(NER)                       :: val_deltaless(1:2,1:2)

    vab = 0.0
    call Vfunc_VDelTPE1_only_cpld(ocpld, p1, p2, vab)
    call Vfunc_sfrTPE1_cpld(ocpld, p1, p2,val_deltaless)

    vab = vab + val_deltaless


  end subroutine Vfunc_VDelTPE1sfr_cpld

  subroutine Vfunc_VDelTPE0_only_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)


    vab = 0.0
    if (ocpld%lsj%j == 0 .or. ocpld%lsj%L /= ocpld%lsj%j-1) return

    vab(1, 1) = TPEdel_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) = TPEdel_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) = TPEdel_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) = TPEdel_jpp(ocpld%lsj%j, p1, p2)

    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
    & vab = vab * ocpld%regulator(p1, p2)


  end subroutine Vfunc_VDelTPE0_only_cpld

  subroutine Vfunc_VDelTPE1_only_cpld(ocpld, p1, p2, vab)

    class(obj_cpld), intent(inout)  :: ocpld
    real(NER), intent(in)           :: p1, p2
    real(NER), intent(out)          :: vab(1:2, 1:2)


    vab = 0.0
    if (ocpld%lsj%j == 0 .or. ocpld%lsj%L /= ocpld%lsj%j-1) return

    vab(1, 1) = TPEdel_sub_jmm(ocpld%lsj%j, p1, p2)
    vab(2, 1) = TPEdel_sub_jpm(ocpld%lsj%j, p1, p2)
    vab(1, 2) = TPEdel_sub_jpm(ocpld%lsj%j, p2, p1)
    vab(2, 2) = TPEdel_sub_jpp(ocpld%lsj%j, p1, p2)

    if (ocpld%regtype /= REGTYPE_NONE .and. ocpld%regtype /= REGTYPE_SHARP)    &
    & vab = vab * ocpld%regulator(p1, p2)




  end subroutine Vfunc_VDelTPE1_only_cpld

  ! < p1 | VNF | p2 > for uncoupled channels
  ! subroutine Vfunc_VNF_sngl(osngl, p1, p2, potval)

  !   class(obj_sngl), intent(inout)    :: osngl
  !   real(NER), intent(in)          :: p1, p2
  !   real(NER), intent(out)         :: potval

  !   real(NER) :: val_deltaless

  !   if (osngl%lsj%S == 0) then
  !     potval = TPEdel_sub_j0j(osngl%lsj%j, p1, p2)
  !   else
  !     if (osngl%lsj%j == osngl%lsj%L) then
  !       potval = TPEdel_sub_j1j(osngl%lsj%j, p1, p2)
  !     else
  !       if (osngl%lsj%L == osngl%lsj%j+1) then
  !         potval = TPEdel_sub_jpp(osngl%lsj%j, p1, p2)
  !       else
  !         potval = TPEdel_sub_jmm(osngl%lsj%j, p1, p2)
  !       end if
  !     end if
  !   end if

  !   if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)    &
  !   & potval = potval * osngl%regulator(p1, p2)



  !end subroutine Vfunc_VNF_sngl
!######################################################################
!VTPE0_deltaless_DR + VTPE0_deltaful_SFR  in sngl chn

  subroutine Vfunc_VDelTPE0_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval
    real(NER)                       :: val_deltaless

    call Vfunc_VDelTPE0_only_sngl(osngl, p1, p2, potval)
    call Vfunc_VTPE0_sngl(osngl, p1, p2,val_deltaless)
    potval=potval+val_deltaless
  end subroutine Vfunc_VDelTPE0_sngl


  subroutine Vfunc_VDelTPE1_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval
    real(NER)                       :: val_deltaless

    call Vfunc_VDelTPE1_only_sngl(osngl, p1, p2, potval)
    call Vfunc_VTPE1_sngl(osngl,p1,p2,val_deltaless)
    potval=potval+val_deltaless

  end subroutine Vfunc_VDelTPE1_sngl
! ################################################################
! vsfr sngl

  subroutine Vfunc_sfrTPE0_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval

    if (osngl%lsj%S == 0) then
      potval = TPESFR_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval = TPESFR_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval = TPESFR_jpp(osngl%lsj%j, p1, p2)
        else
          potval = TPESFR_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_sfrTPE0_sngl


  subroutine Vfunc_sfrTPE1_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval

    if (osngl%lsj%S == 0) then
      potval = TPESFR1_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval = TPESFR1_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval = TPESFR1_jpp(osngl%lsj%j, p1, p2)
        else
          potval = TPESFR1_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_sfrTPE1_sngl

! Vsfr deltaless + deltafull sngl

  subroutine Vfunc_VDelTPE0sfr_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval
    real(NER)                       :: val_deltaless

    call Vfunc_VDelTPE0_only_sngl(osngl, p1, p2, potval)
    call Vfunc_sfrTPE0_sngl(osngl, p1, p2,val_deltaless)
    potval=potval+val_deltaless
  end subroutine Vfunc_VDelTPE0sfr_sngl


  subroutine Vfunc_VDelTPE1sfr_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval
    real(NER)                       :: val_deltaless

    call Vfunc_VDelTPE1_only_sngl(osngl, p1, p2, potval)
    call Vfunc_sfrTPE1_sngl(osngl,p1,p2,val_deltaless)
    potval=potval+val_deltaless

  end subroutine Vfunc_VDelTPE1sfr_sngl

  subroutine Vfunc_VDelTPE0_only_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval


    if (osngl%lsj%S == 0) then
      potval = TPEdel_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval = TPEdel_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval = TPEdel_jpp(osngl%lsj%j, p1, p2)
        else
          potval = TPEdel_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)


  end subroutine Vfunc_VDelTPE0_only_sngl


  subroutine Vfunc_VDelTPE1_only_sngl(osngl, p1, p2, potval)

    class(obj_sngl), intent(inout)    :: osngl
    real(NER), intent(in)          :: p1, p2
    real(NER), intent(out)         :: potval
    real(NER)                       :: val_deltaless

    if (osngl%lsj%S == 0) then
      potval = TPEdel_sub_j0j(osngl%lsj%j, p1, p2)
    else
      if (osngl%lsj%j == osngl%lsj%L) then
        potval = TPEdel_sub_j1j(osngl%lsj%j, p1, p2)
      else
        if (osngl%lsj%L == osngl%lsj%j+1) then
          potval = TPEdel_sub_jpp(osngl%lsj%j, p1, p2)
        else
          potval = TPEdel_sub_jmm(osngl%lsj%j, p1, p2)
        end if
      end if
    end if

    if (osngl%regtype /= REGTYPE_NONE .and. osngl%regtype /= REGTYPE_SHARP)   &
      potval = potval * osngl%regulator(p1, p2)

  end subroutine Vfunc_VDelTPE1_only_sngl

subroutine V0cpld(self,p1,p2,potval)

class(obj_cpld), intent(inout)    :: self
 real(NER), intent(in)          :: p1, p2
 real(NER), intent(out)         :: potval(1:2, 1:2)

 potval=0.0_NER

end subroutine V0cpld

subroutine V0sngl(self,p1,p2,potval)

class(obj_sngl), intent(inout)    :: self
 real(NER), intent(in)          :: p1, p2
 real(NER), intent(out)         :: potval

 potval=0.0_NER

end subroutine V0sngl

!test

! subroutine V0cpld(self,p1,p2,potval)

! class(obj_cpld), intent(inout)    :: self
!  real(NER), intent(in)          :: p1, p2
!  real(NER), intent(out)         :: potval(1:2, 1:2)

!  potval=0.0_NER

! end subroutine V0cpld

! subroutine V0sngl(self,p1,p2,potval)

! class(obj_sngl), intent(inout)    :: self
!  real(NER), intent(in)          :: p1, p2
!  real(NER), intent(out)         :: potval

!  potval=0.0_NER

! end subroutine V0sngl

end module mod_vfunc
