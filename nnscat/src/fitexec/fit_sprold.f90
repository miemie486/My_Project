! Bingwei Long  10/19/2018

program fit_spr

  use drv_fitting
  use util_datetime
  implicit none

  integer, parameter    :: MAX_NUMK = 20, MAX_NMBD = 50, MAX_NPAR = 30
  integer               :: Nmesh, numk, nMmbd, ii, funt, npar_known
  integer               :: npar_unknown, uptoQn, setn_cis
  real(NER)             :: kcm_fitphs(MAX_NUMK, 2)
  real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
  real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
  real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
  real(NER)             :: minchisqr
  character(len = 512)  :: fitphs_fname, inpar_fname, out_fname, timestamp
  character(len = 1024) :: cmdbuff, cmdline, tmp
  integer               :: cmdcnt, ios, chnid
  character(len = 10)   :: chnstr
  logical               :: succeeded
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  Delta  y                                                            ', &
    & '# Mambda  Delta  y  Delta1  y1  C0                                            ', &
    & '# Mambda  Delta  y  Delta1  y1  C0  Delta2  y2  C1  D0                        ', &
    & '# Mambda  Delta  y  Delta1  y1  C0  Delta2  y2  C1  D0  Delta2  y2  C2  D1  E0'  &
    & /)

  fitphs_fname = 'fitspr_phs.in'
  inpar_fname = 'fitspr_inpar.in'
  out_fname = 'fitspr_result.out'
  chnid = CHNID_1S0_PP

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
      & 'fit_spr: Something wrong with the arguments'
    write (standard_error_unit, '(a)')  &
      & 'Usage: fit_spr Nmesh uptoQn [setn_cis]'
    stop
  end if
  if (uptoQn >= 3) then
    read (cmdbuff, *, IOSTAT=ios) Nmesh, uptoQn, setn_cis
    if  (ios /= 0) then
      write (standard_error_unit, '(a)')  &
        & 'fit_spr: for uptoQn = 3, setn_cis is needed'
      write (standard_error_unit, '(a)')  &
        & 'Usage: fit_spr Nmesh uptoQn [setn_cis]'
      stop
    end if
    if (setn_cis < 0 .or. setn_cis > 4) then
      write (standard_error_unit, '(a)') 'fit_spr: setn_cis = 1 - 4'
      stop
    else
      call Select_cis2PE_vtpe1(setn_cis)
    end if
  end if

  call read_fitphs_fitting(fitphs_fname, 1, kcm_fitphs, numk)

  if (uptoQn == 0) then
    npar_known = 0
  else
    npar_known = npar_sprbl1s0(uptoQn)
  end if
  npar_unknown = npar_sprbl1s0(uptoQn+1) - npar_known
  call read_inpar_fitting(inpar_fname, npar_known+npar_unknown, Mambda_inpar,  &
    & nMmbd)

  do ii = 1, nMmbd
    call FitSngl_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,     &
      & Mambda_inpar(ii, 2:npar_known+1), npar_known,                          &
      & kcm_fitphs, numk,                                                      &
      & Mambda_inpar(ii, npar_known+2:npar_known+1+npar_unknown),              &
      &   npar_unknown, 1.0E-4_NER, 1.0E-8_NER,                                &
      & minchisqr, GetSprbl1s0_fitpwa, 1000)
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

  write (funt, '(a)') '# fit to: kcm  phs'
  do ii = 1, numk
    write (funt, '(a, x, f9.3, x, f8.4)') '# ', kcm_fitphs(ii, 1),             &
      & kcm_fitphs(ii, 2)
  end do

  write (funt, '(a)') trim(parname_Qn(uptoQn+1))
  do ii = 1, nMmbd
    write (funt, '(f7.2, 14(2x, e26.17e3))') Mambda_inpar(ii, 1),              &
      & Mambda_inpar(ii, 2:npar_known+npar_unknown+1)
  end do
  close(funt)

end program fit_spr
