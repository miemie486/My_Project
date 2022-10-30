! pray 24/12/2018
!this program fitting nonperterbative OPE in 3p0

program fit_nonperdelta_3p0

  use drv_fitting
  implicit none


  integer               :: npar_Qn(1:10)
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname, inpar_fname, out_fname, errmsg
  logical               :: succeeded
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0  C3  D1                                                  '  &
    & /)
  integer :: regtype
  class(obj_smplchn), allocatable   :: chn
  real(NER)                   :: t0, t1

  call cpu_time(t0)

  fitphs_fname = 'fit3p0_phs.in'
  inpar_fname = 'fit3p0_inpar.in'
  out_fname = 'fit3p0_result.out'
  npar_Qn(1:4) = npar_withc0_sngl(1:4)
  initial_fraction = 1.0E-3_NER
  tolerance = 1.0E-7_NER
  errmsg = 'fit_delta_3p0: Something went wrong'

  call init_smpl
!   call PrintOutFitSngl_fitting(GetNPWithC0Sngl_del_fitpwa,                         &
!     & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
!     & initial_fraction, tolerance, succeeded)
  call PrintOutFitSngl_fitting(GetNPSngl_smpl_fitdelta3p0,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not. succeeded) then
    write(standard_error_unit, '(a)') errmsg
  else
    call chn%erase()
    deallocate(chn)
    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if

contains

  subroutine init_smpl()

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_smplchn::chn)

    chn%vanishing(1) = .true.
    chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL

    chn%pVfunc_V0 => VLO_withc0_epwrap

    chn%vlvs_separated(2) = .true.
    chn%pVfunc_V2  =>  VN2LOSFR_withc0_del_epwrap
    chn%pVfunc_VL2 =>  TPE0SFR_del_epwrap
    chn%pVfunc_VS2 =>  VN2LO_withc0_SHRG_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%vlvs_separated(3) = .true.
    chn%pVfunc_V3  =>  VN3LOSFR_withc0_del_epwrap
    chn%pVfunc_VL3 =>  TPE1SFR_del_epwrap
    chn%pVfunc_VS3 =>  VN3LO_withc0_SHRG_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine init_smpl

  subroutine GetNPSngl_smpl_fitdelta3p0(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    use mod_obj_smplchn
    use drv_pwphase
    use mod_vfunc_smpl
    implicit none

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    call create_smplchn(chn, chnid, uptoQn)
    call SetUp_withc0_SFRdel_smpl(chn)

    call get_phase_as_value_for_klist_smpl(chn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

  end subroutine GetNPSngl_smpl_fitdelta3p0

end program fit_nonperdelta_3p0
