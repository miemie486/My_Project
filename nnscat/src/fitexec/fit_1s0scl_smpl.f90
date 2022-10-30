! Pray  10/21/2020  for 1s0 with scattering length

program fit_1s0scl_smpl

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
  real(NER)                   :: t0, t1
  call cpu_time(t0)

  fitphs_fname = 'fit1s0_phs.in'
  inpar_fname = 'fit1s0_inpar.in'
  out_fname = 'fit1s0_result.out'
  npar_Qn(1:4) = (/ 1, 3, 6, 9/)
  initial_fraction = 1.0E-3_NER
  tolerance = 1.0E-7_NER
  errmsg = 'fit_1s0scl_smpl: Something went wrong'

  call PrintOutFit1s0_SCL_fitting(Get1s0_withscl_smpl_fitpwa,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not. succeeded) then
    write(standard_error_unit, '(a)') trim(errmsg)
  else
    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if

end program fit_1s0scl_smpl