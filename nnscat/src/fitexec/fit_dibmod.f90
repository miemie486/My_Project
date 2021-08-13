! Bingwei Long  10/28/2018

program fit_dibmod

  use drv_fitting
  implicit none

  integer               :: npar_Qn(1:10)
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname, inpar_fname, out_fname, errmsg
  logical :: succeeded
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  Delta  y                                                            ', &
    & '# Mambda  Delta  y  Delta1  y1  C0                                            ', &
    & '# Mambda  Delta  y  Delta1  y1  C0  Delta2  y2  C1  D0                        ', &
    & '# Mambda  Delta  y  Delta1  y1  C0  Delta2  y2  C1  D0  Delta2  y2  C2  D1  E0'  &
    & /)

  fitphs_fname = 'fitlbwdib_phs.in'
  inpar_fname = 'fitlbwdib_inpar.in'
  out_fname = 'fitlbwdib_result.out'
  npar_Qn(1:4) = npar_lbwdib(1:4)
  initial_fraction = 1.0E-3_NER
  tolerance = 1.0E-7_NER
  errmsg = 'fit_dibmod: Something went wrong'

  call PrintOutFitSngl_fitting(GetNPWithC0Sngl_fitpwa,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not. succeeded) write(standard_error_unit, '(a)') errmsg

end program fit_dibmod