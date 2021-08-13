! Bingwei Long  10/19/2018

! Using Getspr1s0_smpl_fitpwa to fit parameters for spr at each order
! If necessary, change Getspr1s0_smpl_fitpwa to reflect the algorithm wished for
program fitspr_smpl

  use drv_fitting
  implicit none

  integer               :: npar_Qn(1:10)
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname, inpar_fname, out_fname, errmsg
  logical               :: succeeded
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  Delta  y                                                            ', &
    & '# Mambda  Delta  y  Delta1  y1  C0                                            ', &
    & '# Mambda  Delta  y  Delta1  y1  C0  Delta2  y2  C1  D0                        ', &
    & '# Mambda  Delta  y  Delta1  y1  C0  Delta2  y2  C1  D0  Delta3  y3  C2  D1  E0'  &
    & /)

  fitphs_fname = 'fitspr_phs.in'
  inpar_fname = 'fitspr_inpar.in'
  out_fname = 'fitspr_result.out'
  npar_Qn(1:4) = NPAR_SPR1S0_SMPL
  initial_fraction = 1.0E-4_NER
  tolerance = 1.0E-10_NER
  errmsg = 'fitspr_smpl: Something went wrong'

  ! V1_BASEN_sprbl_mvfs = 2
  ! V1_BND_sprbl_mvfs = 1.0E-3_NER
  ! V2_BASEN_sprbl_mvfs = 2
  ! V2_BND_sprbl_mvfs = 1.0E-4_NER
  ! V3_BASEN_sprbl_mvfs = 2
  ! V3_BND_sprbl_mvfs = 1.0E-4_NER

  call PrintOutFitSngl_fitting(Getspr1s0_smpl_fitpwa,                          &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not. succeeded) write(standard_error_unit, '(a)') trim(errmsg)

end program fitspr_smpl
