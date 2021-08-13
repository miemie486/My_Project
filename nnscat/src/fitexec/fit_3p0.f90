! Bingwei Long  10/15/2018  for DR potential

program fit_3p0

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

  fitphs_fname = 'fit3p0_phs.in'
  inpar_fname = 'fit3p0_inpar.in'
  out_fname = 'fit3p0_result.out'
  npar_Qn(1:4) = npar_withc0_sngl(1:4)
  initial_fraction = 1.0E-3_NER
  tolerance = 1.0E-7_NER
  errmsg = 'fit_3p0: Something went wrong'

  call PrintOutFitSngl_fitting(GetNPWithC0Sngl_fitpwa,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not. succeeded) write(standard_error_unit, '(a)') trim(errmsg)

end program fit_3p0
