! pray 09/16/2020
!this program fitting nonperterbative 3s1 with binding energy

program fit_3s1bd_smpl

  use drv_fitting
  implicit none


  integer               :: npar_Qn(1:10)
  real(NER)             :: initial_fraction, tolerance
  character(len = 512)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname, errmsg
  character(len = 512)  :: inpar_fname, out_fname
  logical               :: succeeded
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0  C3  D1                                                  '  &
    & /)

  real(NER)                   :: t0, t1
  call cpu_time(t0)

      npar_Qn(1:4) = (/ 1, 1, 4, 7/)
      fitphs1_fname = 'fit3s1_phs1.in'
      fitphs2_fname = 'fit3s1_phs2.in'
      fitmxa_fname = 'fit3s1_mxa.in'
      inpar_fname = 'fit3s1_inpar.in'
      out_fname = 'fit3s1_result.out'
      initial_fraction = 1.0E-4_NER
!       tolerance = 1.0E-8_NER
      tolerance = 1.0E-12_NER
!       initial_fraction = 1.0E-3_NER
!       tolerance = 1.0E-8_NER
      errmsg = 'fit_3s1bd_smpl: Something went wrong'


  call PrintOutFit3s1_BD_fitting(Get3s1_withbd_smpl_fitpwa, fitphs1_fname,      &
    & fitphs2_fname, fitmxa_fname, inpar_fname, npar_Qn, out_fname,            &
    & parname_Qn, initial_fraction, tolerance, succeeded)
  if (.not.succeeded) then
    write(standard_error_unit, '(a)') trim(errmsg)
  else
    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if

end program fit_3s1bd_smpl