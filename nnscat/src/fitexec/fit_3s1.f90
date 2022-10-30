! Bingwei Long  10/13/2018

program fit_3s1

  use drv_fitting
  use util_datetime
  implicit none

  external GetNPCpld_smpl_fit3s1

  integer, parameter    :: MAX_NUMK = 100, MAX_NMBD = 50, MAX_NPAR = 30
  integer               :: Nmesh, Nphs1, Nphs2, Nmxa, nMmbd, ii, funt, npar_known
  integer               :: chnid, npar_unknown, uptoQn, setn_cis
  integer               :: cmdcnt, ios
  real(NER)             :: kcm_fitphs1(MAX_NUMK, 2), kcm_fitphs2(MAX_NUMK, 2)
  real(NER)             :: kcm_fitmxa(MAX_NUMK, 2)
  real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
  real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
  real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
  real(NER)             :: minchisqr
  character(len = 512)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname
  character(len = 512)  :: inpar_fname, out_fname
  character(len = 80)   :: timestamp, chnstr
  character(len = 1024) :: cmdbuff, cmdline, tmp
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0  E0                                                      ', &
    & '# Mambda  C0  C2  D0  E0  C3  D1  E1                                          '  &
    & /)
  logical               :: succeeded
  real(NER)                   :: t0, t1

  call cpu_time(t0)

  fitphs1_fname = 'fit3s1_phs1.in'
  fitphs2_fname = 'fit3s1_phs2.in'
  fitmxa_fname = 'fit3s1_mxa.in'
  inpar_fname = 'fit3s1_inpar.in'
  out_fname = 'fit3s1_result.out'
  chnid = CHNID_3S1_PP

  call get_command(cmdline)
  cmdbuff = ''
  cmdcnt = command_argument_count()
  do ii = 1, cmdcnt
    call get_command_argument(ii, tmp)
    cmdbuff = trim(cmdbuff)//' '//trim(tmp)
  end do
  read (cmdbuff, *, IOSTAT=ios) Nmesh, uptoQn
  if  (ios /= 0) then
    write (standard_error_unit, '(a)')  &
      & 'fit_3s1: Something wrong with the arguments'
    write (standard_error_unit, '(a)')  &
      & 'Usage: fit_3s1 Nmesh uptoQn [setn_cis]'
    stop
  end if
  if (uptoQn >= 3) then
    read (cmdbuff, *, IOSTAT=ios) Nmesh, uptoQn, setn_cis
    if  (ios /= 0) then
      write (standard_error_unit, '(a)')  &
        & 'fit_3s1: for uptoQn = 3, setn_cis is needed'
      write (standard_error_unit, '(a)')  &
        & 'Usage: fit_3s1 Nmesh uptoQn [setn_cis]'
      stop
    end if
!     if (setn_cis < 0 .or. setn_cis > 4) then
!       write (standard_error_unit, '(a)') 'fit_3s1: setn_cis = 1 - 4'
!       stop
!     else
!       call Select_cis2PE_vtpe1(setn_cis)
!     end if
    call Select_cis2PE_vtpe1(setn_cis)
  end if

  call read_fitphs_fitting(fitphs1_fname, 1, kcm_fitphs1, Nphs1)
  call read_fitphs_fitting(fitphs2_fname, 1, kcm_fitphs2, Nphs2)
  call read_fitphs_fitting(fitmxa_fname, 1, kcm_fitmxa, Nmxa)

  if (uptoQn == 0) then
    npar_known = 0
  else
    npar_known = NPAR_WITHC0_CPLD_SMPL(uptoQn)
  end if
  npar_unknown = NPAR_WITHC0_CPLD_SMPL(uptoQn+1) - npar_known
  call read_inpar_fitting(inpar_fname, npar_known+npar_unknown, Mambda_inpar,  &
    & nMmbd)

  do ii = 1, nMmbd
    call FitCpld_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,            &
      & Mambda_inpar(ii, 2:npar_known+1), npar_known,                          &
      & kcm_fitphs1, Nphs1, kcm_fitphs2, Nphs2, kcm_fitmxa, Nmxa,              &
      & Mambda_inpar(ii, npar_known+2:npar_known+1+npar_unknown),              &
      &   npar_unknown, 1.0E-3_NER, 1.0E-8_NER,                                &
      & minchisqr, GetNPCpld_smpl_fit3s1, 1000)
      ! & minchisqr, GetNPCpld_fitpwa, 1000)
  end do

  funt = get_free_unit_io()
  open (funt, file = out_fname, status = 'replace')

  call get_timestring_hms(timestamp)
  write (funt, '(a)') '# '//trim(timestamp)
  write (funt, '(a)') '# '//trim(cmdline)
  call convert_chnid_to_text(chnid, chnstr, succeeded)
  if (succeeded) write (funt, '(a)') '# CHN = '//trim(chnstr)
  write (funt, '(a, I1)') '# uptoQn = ', uptoQn
  if (uptoQn >= 3) write (funt, '(a, I1, x, a, I1)') '# setn_cis = ', setn_cis
  write (funt, '(a, I5)') '# Nmesh = ', Nmesh
  call write_to_file_basic_phycnst(funt)

  write (funt, '(a)') '# fit to: kcm  phs1'
  do ii = 1, Nphs1
    write (funt, '(a, x, f9.3, x, f8.4)') '# ', kcm_fitphs1(ii, 1),            &
      & kcm_fitphs1(ii, 2)
  end do

  write (funt, '(a)') '# fit to: kcm  phs2'
  do ii = 1, Nphs2
    write (funt, '(a, x, f9.3, x, f8.4)') '# ', kcm_fitphs2(ii, 1),            &
      & kcm_fitphs2(ii, 2)
  end do

  write (funt, '(a)') '# fit to: kcm  mxa'
  do ii = 1, Nmxa
    write (funt, '(a, x, f9.3, x, f8.4)') '# ', kcm_fitmxa(ii, 1),             &
      & kcm_fitmxa(ii, 2)
  end do

  write (funt, '(a)') trim(parname_Qn(uptoQn+1))
  do ii = 1, nMmbd
    write (funt, '(f7.2, 14(2x, e26.17e3))') Mambda_inpar(ii, 1),              &
      & Mambda_inpar(ii, 2:npar_known+npar_unknown+1)
  end do
  close(funt)

  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'

end program fit_3s1

subroutine GetNPCpld_smpl_fit3s1(Nmesh, Mambda, chnid, uptoQn, paralst,        &
    & num_paras, klist, numk, phaselst)

    use mod_obj_smplchn
    use drv_pwphase
    use mod_vfunc_smpl
    implicit none

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    integer :: regtype
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    ptrchn%vanishing(1) = .true.
    ptrchn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL

    ptrchn%pVfunc_V0 => VLO_withc0_epwrap

    ptrchn%pVfunc_V2 => VN2LO_withc0_epwrap
    ptrchn%V2_cheb_base_n = 4
    ptrchn%V2_cheb_bnd = 1.0E-3_NER

    ptrchn%pVfunc_V3 => VN3LO_withc0_epwrap
    ptrchn%V3_cheb_base_n = 4
    ptrchn%V3_cheb_bnd = 1.0E-3_NER

    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

    call ptrchn%erase()
    deallocate(ptrchn)

end subroutine GetNPCpld_smpl_fit3s1
