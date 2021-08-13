! 27/7/2020
! fit for MMWLY power counting
! add scattering length in 1s0 andbinding energy in 3s1 for fitting

program fit_mmwly_smpl

  use drv_fitting
  implicit none


  integer               :: npar_Qn(1:10), ios, cmdcnt, ii, N
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname, inpar_fname, out_fname, errmsg
  character(len = 512)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname
  logical               :: succeeded
  character(len = 3)          :: cmdstr
  character(len = 1024)       :: cmdbuff, tmp
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0  C3  D1                                                  '  &
    & /)
  integer :: regtype
  class(obj_smplchn), allocatable   :: chn
  type(lsj_symbol)            :: lsj
  real(NER)                   :: t0, t1

  call cpu_time(t0)
  cmdbuff = ''
  cmdcnt = command_argument_count()
  do ii = 1, cmdcnt
    call get_command_argument(ii, tmp)
    cmdbuff = trim(cmdbuff)//' '//trim(tmp)
  end do
  read (cmdbuff, *, IOSTAT=ios) N, cmdstr

  select case (trim(cmdstr))
    case('1s0')
      npar_Qn(1:3) = (/1, 3, 6/)
    case('3s1')
      npar_Qn(1:3) = (/1, 1, 4/)
    case ('3p0')
      npar_Qn(1:3) =  (/1, 1, 3/)
    case ('1d2', '3d2', '3d3')
      npar_Qn(1:3) = (/ 0, 0, 0/)
    case ('1p1', '3p1')
      npar_Qn(1:3) = (/ 0, 0, 1/)
    case ('3p2')
      npar_Qn(1:3) = (/ 0, 0, 1/)
    case default
      npar_Qn(1:3) = (/ 0, 0, 0/)
  end select

  call convert_text_to_lsj(trim(cmdstr), lsj, succeeded)
  if (.not.succeeded) then
    print *, 'convert chn_symbol to lsj not succeeded'
  end if

    inpar_fname = 'fit'//trim(cmdstr)//'_inpar.in'
    out_fname = 'fit'//trim(cmdstr)//'_result.out'
    errmsg = 'fit_mmwly_smpl: Something went wrong'

  if (lsj%L /= lsj%j .and. lsj%j /= 0) then
    fitphs1_fname = 'fit'//trim(cmdstr)//'_phs1.in'
    fitphs2_fname = 'fit'//trim(cmdstr)//'_phs2.in'
    fitmxa_fname = 'fit'//trim(cmdstr)//'_mxa.in'
    initial_fraction = 1.0E-4_NER
    tolerance = 1.0E-5_NER
  else
    fitphs1_fname = 'fit'//trim(cmdstr)//'_phs.in'
    fitphs2_fname = ''
    fitmxa_fname = ''
    initial_fraction = 1.0E-3_NER
    tolerance = 1.0E-7_NER
  end if

  call PrintOutFitSLBEPS_fitting(Get_SLBEPS_fitpwa, fitphs1_fname,      &
    & fitphs2_fname, fitmxa_fname, inpar_fname, npar_Qn, out_fname,            &
    & parname_Qn, initial_fraction, tolerance, succeeded)
  if (.not.succeeded) then
    write(standard_error_unit, '(a)') trim(errmsg)
  else
    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if

end program fit_mmwly_smpl