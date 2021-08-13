! pray  9/18/2019
! using smpl and SFR potential
! number of parameter(1,1,4,7)

program fit_deltasmpl_3s1


  use drv_fitting
  use potwrap
  ! use mod_vfunc_smpl
  ! use drv_pwphase
  implicit none

  integer               :: npar_Qn(1:10), chnid, uptoQn
  real(NER)             :: initial_fraction, tolerance
  character(len = 512)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname, errmsg
  character(len = 512)  :: inpar_fname, out_fname
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0  E0                                                      ', &
    & '# Mambda  C0  C2  D0  E0  C3  D1  E1                                          '  &
    & /)
  logical               :: succeeded
  integer :: regtype
  class(obj_smplchn), allocatable   :: chn
  real(NER)             :: t0, t1

  call cpu_time(t0)
  fitphs1_fname = 'fit3s1_phs1.in'
  fitphs2_fname = 'fit3s1_phs2.in'
  fitmxa_fname = 'fit3s1_mxa.in'
  inpar_fname = 'fit3s1_inpar.in'
  out_fname = 'fit3s1_result.out'
  npar_Qn(1:4) = NPAR_WITHC0_CPLD_SMPL
  initial_fraction = 1.0E-4_NER
  tolerance = 1.0E-5_NER
  errmsg = 'fit_nonperdelta_3s1: Something went wrong'

  call init_smpl()

  call PrintOutFitCpld_fitting(GetNPCpldSfr_smpl_fitdelta3s1, fitphs1_fname,      &
    & fitphs2_fname, fitmxa_fname, inpar_fname, npar_Qn, out_fname,            &
    & parname_Qn, initial_fraction, tolerance, succeeded)
  if (.not.succeeded) then
    write(standard_error_unit, '(a)') trim(errmsg)
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
    chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL

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

  subroutine GetNPCpldSfr_smpl_fitdelta3s1(Nmesh, Mambda, chnid, uptoQn, paralst, &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    call create_smplchn(chn, chnid, uptoQn)
    call get_phase_as_value_for_klist_smpl(chn, Nmesh, regtype, Mambda, klist, &
      & numk, paralst, num_paras, phaselst)

  end subroutine GetNPCpldSfr_smpl_fitdelta3s1

end program fit_deltasmpl_3s1