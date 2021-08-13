program fit_supdeltasfr_psngl

  use drv_fitting
  implicit none

  integer               :: npar_Qn(1:10)
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname, inpar_fname, out_fname, errmsg
  logical               :: succeeded
  character(len = 78), dimension(1:5), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0  C3  D1                                                  '  &
    & /)
  real(NER)                   :: t0, t1

  call cpu_time(t0)

  fitphs_fname = 'fitsngl_phs.in'
  inpar_fname = 'fitsngl_inpar.in'
  out_fname = 'fitsngl_result.out'
  npar_Qn(1:5) = (/0,0,1,2,4/)
  initial_fraction = 1.0E-3_NER
  tolerance = 1.0E-7_NER
  errmsg = 'fit_supdeltasfr_psngl: Something went wrong'

  call PrintOutFitSngl_fitting(Getsfrdeltasngl_smpl_fitpwa,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not.succeeded) then
      write(standard_error_unit, '(a)') trim(errmsg)
    else
      call cpu_time(t1)
      print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if

end program fit_supdeltasfr_psngl