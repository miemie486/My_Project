
program fit_supsfr_pwave_contacttest

  use drv_fitting
  use mod_vfunc_smpl
  implicit none

  integer               :: npar_Qn(1:10), ios, cmdcnt, ii, N
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname
  character(len = 512)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname, errmsg
  character(len = 512)  :: inpar_fname, out_fname
  character(len = 3)          :: cmdstr
  character(len = 1024)       :: cmdbuff, tmp
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
  cmdbuff = ''
  cmdcnt = command_argument_count()
  do ii = 1, cmdcnt
    call get_command_argument(ii, tmp)
    cmdbuff = trim(cmdbuff)//' '//trim(tmp)
  end do
  read (cmdbuff, *, IOSTAT=ios) N, cmdstr
  select case (trim(cmdstr))
!     case ('1d2', '3d2', '3d3')
!       npar_Qn(1:5) = (/ 0, 0, 0, 0, 1/)
    case ('1p1', '3p1', '3p0')
      npar_Qn(1:5) = (/ 0, 0, 1, 3, 6/)
    case ('3p2')
      npar_Qn(1:5) = (/ 0, 0, 1, 4, 10/)
!     case default
!       npar_Qn(1:5) = (/ 0, 0, 0, 0, 0/)
  end select
select case(trim(cmdstr))
  case('3p2')
    fitphs1_fname = 'fitcpld_phs1.in'
    fitphs2_fname = 'fitcpld_phs2.in'
    fitmxa_fname = 'fitcpld_mxa.in'
    inpar_fname = 'fitcpld_inpar.in'
    out_fname = 'fitcpld_result.out'
    initial_fraction = 1.0E-3_NER
    tolerance = 1.0E-8_NER
    errmsg = 'fit_supsfr_pwave_contacttest: Something went wrong'

    call PrintOutFitcpld_fitting(Getsfrcpld_smpl_test_fitpwa,                         &
    & fitphs1_fname, fitphs2_fname, fitmxa_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
    if (.not.succeeded) then
        write(standard_error_unit, '(a)') trim(errmsg)
      else
        call cpu_time(t1)
        print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
    end if

  case default
    fitphs_fname = 'fitsngl_phs.in'
    inpar_fname = 'fitsngl_inpar.in'
    out_fname = 'fitsngl_result.out'
    initial_fraction = 1.0E-3_NER
    tolerance = 1.0E-7_NER
    errmsg = 'fit_supsfr_pwave_contacttest: Something went wrong'

    call PrintOutFitSngl_fitting(Getsfrsngl_smpl_test_fitpwa,                         &
      & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
      & initial_fraction, tolerance, succeeded)
    if (.not.succeeded) then
        write(standard_error_unit, '(a)') trim(errmsg)
      else
        call cpu_time(t1)
        print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
    end if
end select


contains

  subroutine Getsfrsngl_smpl_test_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUp_pope_SFR_smpl_test(ptrchn)
    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Getsfrsngl_smpl_test_fitpwa

  subroutine Getsfrcpld_smpl_test_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUp_pope_SFR_smpl_test(ptrchn)

    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Getsfrcpld_smpl_test_fitpwa

end program fit_supsfr_pwave_contacttest