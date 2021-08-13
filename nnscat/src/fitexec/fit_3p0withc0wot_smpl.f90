! pray 11/8/2019  for 3p0 with C0  without TPE unsuppressed-per

program fit_3p0withc0wot_smpl

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
  real(NER)             :: t0, t1

  call cpu_time(t0)

  call init_smpl()

  fitphs_fname = 'fit3p0_phs.in'
  inpar_fname = 'fit3p0_inpar.in'
  out_fname = 'fit3p0_result.out'
  npar_Qn(1:5) = chn%npar(1:5)
  initial_fraction = 1.0E-3_NER
  tolerance = 1.0E-7_NER
  errmsg = 'fit_3p0withc0wot_smpl: Something went wrong'


  call PrintOutFitSngl_fitting(GetNPSngl_smpl_fit3p0withc0,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not. succeeded) write(standard_error_unit, '(a)') trim(errmsg)

    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'

contains

  subroutine init_smpl()

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_smplchn::chn)

    chn%vanishing(1) = .false.

!     chn%npar(1:5) = (/ 0, 1, 3, 5, 5 /)

    chn%npar(1:5) = (/ 0, 1, 3, 6, 6 /)

    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 =>  VLO_withc0_epwrap

    chn%pVfunc_V2 => VN2LO_perwot3p0_smpl_epwrap
    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

!     chn%pVfunc_V3 => VN3LO_per3p0_smpl_epwrap
    chn%pVfunc_V3 => VN3LO_perwot3p0_high_smpl_epwrap
    chn%V3_cheb_base_n = 8
    chn%V3_cheb_bnd = 1.0E-9_NER



  end subroutine init_smpl

  subroutine GetNPSngl_smpl_fit3p0withc0(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    call create_smplchn(chn, chnid, uptoQn)

    call get_phase_as_value_for_klist_smpl(chn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

  end subroutine GetNPSngl_smpl_fit3p0withc0

end program fit_3p0withc0wot_smpl