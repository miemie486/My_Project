! Bingwei Long 07/26/2018

! Driver routines for phase shift fitting

module drv_fitting

  use drv_fitpwa
  use util_datetime
  use util_gadgets
  use mod_deltafullN2LO_pwd
  use mod_obj_smplchn, only:allpert_BDE_initial, allpert_BDE
  ! use omp_lib

  implicit none

  ! external subroutines from MINUIT
  external MNINIT, MNEXCM, MNPARM, MNSETI, MNPOUT

  integer, parameter  :: MAXNUM_CALLMINUIT = 500

contains

  subroutine ParseCommandLine_fitting(cmdline, N, chnstr, uptoQn, setn_cis, succeeded)

    integer, intent(out)  :: N, uptoQn, setn_cis
    character(len = *), intent(out) :: cmdline, chnstr
    logical, intent(out)  :: succeeded

    integer :: cmdcnt, ii, ios
    character(len = 1024) :: cmdbuff, tmp

    succeeded = .true.
    call get_command(cmdline)
    cmdbuff = ''
    cmdcnt = command_argument_count()
    do ii = 1, cmdcnt
      call get_command_argument(ii, tmp)
      cmdbuff = trim(cmdbuff)//' '//trim(tmp)
    end do
    read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn
    if  (ios /= 0) then
      write (standard_error_unit, '(a)')  &
        & 'ParseCommandLine: Something wrong with the arguments'
      write (standard_error_unit, '(a)')  &
        & 'Usage: fit_something Nmesh chn uptoQn [setn_cis]'
      succeeded = .false.
    end if
    if (uptoQn >= 3) then
      read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, setn_cis
      if  (ios /= 0) then
        write (standard_error_unit, '(a)')  &
          & 'ParseCommandLine: for uptoQn = 3, setn_cis is needed'
        write (standard_error_unit, '(a)')  &
          & 'Usage: fit_something Nmesh chn uptoQn [setn_cis]'
        succeeded = .false.
      end if
    end if

  end subroutine ParseCommandLine_fitting

  ! Input:
  !   numphs : number of phases at each kcm (1 for uncoupled channel)
  ! Output:
  !   fitphs(:, 1) : list of kcm to be fitted at
  ! Format of input file named fname
  ! kcm1        delta_A1, delta_B1 (if any), ...
  ! kcm2        delta_A2, delta_B2 (if any), ...
  ! ...
  subroutine read_fitphs_fitting(fname, numphs, fitphs, numrow)

    character(len = *), intent(in)  :: fname
    integer, intent(in)             :: numphs
    real(NER), intent(out)          :: fitphs(:, :)
    integer, intent(out)            :: numrow
    integer :: flag, funt

    funt = get_free_unit_io()
    call read_ncols_io(funt, fname, numphs+1, fitphs, numrow, flag)

    if (flag /= 0) then
      write (standard_error_unit, '(a)')    &
        & 'read_fitphs_fitting: couldn''t open "'//trim(fname)//'"'
      return
    end if

  end subroutine read_fitphs_fitting

  ! Format of input file
  ! Mambda1    inpar1(1)    inpar1(2)   ...
  ! Input:
  !   num_inpar: how many input parameters to be read for each row
  ! Output:
  !   Mambda_inpar(:, 1): Mambda
  !   Mambda_inpar(:, 2): the first input parameter
  !   Mambda_inpar(:, 3): the second input parameter
  ! ...
  subroutine read_inpar_fitting(fname, num_inpar, Mambda_inpar, numrow)

    character(len = *), intent(in)  :: fname
    integer, intent(in)             :: num_inpar
    real(NER), intent(out)          :: Mambda_inpar(:, :)
    integer, intent(out)            :: numrow

    integer :: funt, flag

    funt = get_free_unit_io()
    call read_ncols_io(funt, trim(fname), num_inpar+1, Mambda_inpar, numrow,   &
      & flag)
    if (flag == read_ncols_FAIL) then
      write (standard_error_unit, '(a)')  &
        & 'read_inpar_fitting: couldn''t open or read '//trim(fname)
    end if

  end subroutine read_inpar_fitting

  ! subroutine write_results_fitting(funt, parname_Qn, nMmbd, Mambda_paras)

  ! end subroutine write_results_fitting

  ! Input:
  ! uptoQn : the order to be fitted at
  ! unknpara: suggested initial values for paras to be fitted
  ! Output:
  ! unknpara: resulting values
  ! minchsqr: chi square for resulting values
  subroutine  Fit1s0lbwdib_fitting(Nmesh, Mambda, uptoQn, knpara, klist, numk, &
    & fitphase, unknpara, minchisqr)

    integer, intent(in)       :: Nmesh, uptoQn, numk
    real(NER), intent(in)     :: Mambda, knpara(:), klist(:), fitphase(:)
    real(NER), intent(inout)  :: unknpara(:)
    real(NER), intent(out)    :: minchisqr

    ! npar_uptoQn: total number of inparas up to this order
    ! npar_uptoPrev: total number of inparas up to previous order
    integer   :: errflag, jj, ivarbl, npar_uptoQn, npar_uptoPrev, npar_unknown
    real(NER) :: argmn(1:5), parup, parlow, parerr
    character(len = 4)  :: parname(1:MAX_NUMBER_ORDER)

    npar_uptoQn   = npar_lbwdib(uptoQn+1)
    if (uptoQn > 0) then
      npar_uptoPrev = npar_lbwdib(uptoQn)
    else
      npar_uptoPrev = 0
    end if
    npar_unknown  = npar_uptoQn - npar_uptoPrev

    call MNINIT (5, 6, 7)
    call MNSETI ('Fit1s0lbwdib_fitting')

    ! -1 no output except from SHOW commands
    ! 0 minimum output (no starting values or intermediate results)
    ! 1 default value, normal output
    ! 2 additional output giving intermediate results.
    ! 3 maximum output,showing progress of minimizations.
    argmn(1) = -1.0_NER
    call MNEXCM (chisqr, 'SET PRINTOUT', argmn, 1, errflag, function_zero)

    do jj = 1, npar_unknown
      write (parname(jj), '(a, I1.1)') 'g', jj
      ! initial steps: abs(unknpara(jj)*1.0E-3
      call MNPARM (jj, trim(parname(jj)), unknpara(jj),                        &
        & abs(unknpara(jj)*1.0E-3_NER), 0.0_NER, 0.0_NER, errflag)
    end do
    ! Strategy level = 0, 1, 2 (higher number means more derivative computes)
    argmn(1) = 1.0_NER
    call MNEXCM (chisqr, 'SET STRATEGY', argmn, 1, errflag, function_zero)
    ! max number of calls
    argmn(1) = 1000
    ! tolerance
    argmn(2) = abs(unknpara(jj)*1.0E-8_NER)
    call MNEXCM (chisqr, 'MIGRAD', argmn, 2, errflag, function_zero)
    do jj = 1, npar_unknown
      call MNPOUT (jj, trim(parname(jj)), unknpara(jj), parerr, parup, parlow, &
        & ivarbl)
    end do

  contains

    subroutine chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

      integer, intent(in)     :: NPAR, IFLAG
      real(NER), intent(in)   :: XVAL(1:NPAR)
      real(NER), intent(out)  :: FVAL
      real(NER), intent(out)  :: GRAD(1:NPAR)

      external FUTIL

      integer   :: ii
      real(NER) :: predct_phase(1:numk, 1:uptoQn+1), paralst(1:npar_uptoQn)

      paralst(1:npar_uptoPrev) = knpara(1:npar_uptoPrev)
      paralst(npar_uptoPrev+1:npar_uptoQn) = XVAL(1:NPAR)
      call Get1s0dib_Long_fitpwa(Nmesh, Mambda, uptoQn, paralst, npar_uptoQn,  &
        & klist, numk, predct_phase)

      FVAL = 0.0_NER
      do ii = 1, numk
        FVAL = FVAL + (sum(predct_phase(ii, 1:uptoQn+1)) - fitphase(ii))**2/   &
          & abs(klist(ii))
      end do

    end subroutine chisqr

  end subroutine Fit1s0lbwdib_fitting

  ! kphs_lst(i, 1) = the ith kcm, kphs_lst(i, 2) = the ith phase shift
  ! initfrac: initial step = initfrac*initvalue, advised value: 1.0E-3
  ! tolerance: advised value: 1.0E-7 or larger
  ! Nmaxcall: max number of calls ~ 1000

  subroutine  FitSngl_fitting(Nmesh, Mambda, chnid, uptoQn, knwpara, Nknw,     &
    & kphs_lst, numk, unknwpara, Nunknw, initfrac, tolerance, minchisqr,       &
    & GetPhaseSngl, Nmaxcall)

    integer, intent(in)       :: Nmesh, chnid, uptoQn, numk, Nknw, Nunknw
    integer, intent(in)       :: Nmaxcall
    real(NER), intent(in)     :: Mambda, initfrac, tolerance
    real(NER), intent(in)     :: knwpara(:), kphs_lst(:, :)
    real(NER), intent(inout)  :: unknwpara(:)
    real(NER), intent(out)    :: minchisqr
    procedure(GetSngl_fitpwa) :: GetPhaseSngl

    ! npar_uptoQn: total number of inparas up to this order
    integer   :: errflag, jj, ivarbl, npar_uptoQn
    real(NER) :: argmn(1:5), parup, parlow, parerr
    character(len = 4)  :: parname(1:MAX_NUMBER_ORDER)
    real(NER) :: predct_phase(1:numk, 1:uptoQn+1), paralst(1:Nknw + Nunknw)

    npar_uptoQn = Nknw + Nunknw

    call MNINIT (5, 6, 7)
    call MNSETI ('FitSngl_fitting')

    ! -1 no output except from SHOW commands
    ! 0 minimum output (no starting values or intermediate results)
    ! 1 default value, normal output
    ! 2 additional output giving intermediate results.
    ! 3 maximum output,showing progress of minimizations.
    argmn(1) = -1.0_NER
    call MNEXCM (chisqr, 'SET PRINTOUT', argmn, 1, errflag, function_zero)

    do jj = 1, Nunknw
      write (parname(jj), '(a, I1.1)') 'g', jj
      ! initial steps: abs(unknwpara(jj)*initfrac
      call MNPARM (jj, trim(parname(jj)), unknwpara(jj),                       &
        & abs(unknwpara(jj)*initfrac), 0.0_NER, 0.0_NER, errflag)
    end do
    ! Strategy level = 0, 1, 2 (higher number means more derivative computes)
    argmn(1) = 1.0_NER
    call MNEXCM (chisqr, 'SET STRATEGY', argmn, 1, errflag, function_zero)
    ! max number of calls
    argmn(1) = real(Nmaxcall)
    ! tolerance
    argmn(2) = abs(unknwpara(jj)*tolerance)
    call MNEXCM (chisqr, 'MIGRAD', argmn, 2, errflag, function_zero)
    do jj = 1, Nunknw
      call MNPOUT (jj, trim(parname(jj)), unknwpara(jj), parerr, parup, parlow, &
        & ivarbl)
    end do

  contains

    subroutine chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

      integer, intent(in)     :: NPAR, IFLAG
      real(NER), intent(in)   :: XVAL(1:NPAR)
      real(NER), intent(out)  :: FVAL
      real(NER), intent(out)  :: GRAD(1:NPAR)

      external FUTIL

      integer   :: ii

      if (Nknw > 0) paralst(1:Nknw) = knwpara(1:Nknw)
      paralst(Nknw+1:npar_uptoQn) = XVAL(1:NPAR)
      call GetPhaseSngl(Nmesh, Mambda, chnid, uptoQn, paralst, npar_uptoQn,    &
        & kphs_lst(1:numk, 1), numk, predct_phase)

      FVAL = 0.0_NER
      do ii = 1, numk
!         FVAL = FVAL + (sum(predct_phase(ii, 1:uptoQn+1)) - kphs_lst(ii, 2))**2 &
!           &               /abs(kphs_lst(ii, 1))
        FVAL = FVAL + (sum(predct_phase(ii, 1:uptoQn+1)) - kphs_lst(ii, 2))**2
      end do

    end subroutine chisqr

  end subroutine FitSngl_fitting

  ! kphs1_lst(i, 1) = the ith kcm, kphs1_lst(i, 2) = the ith phase shift (lower L)
  ! kphs2_lst(i, 1) = the ith kcm, kphs2_lst(i, 2) = the ith phase shift (higher L)
  ! kmxa_lst(i, 1) = the ith kcm, kmxa_lst(i, 2) = the ith mixing angle
  ! initfrac: initial step = initfrac*initvalue, recommended value: 1.0E-3
  ! tolerance: recommended value: 1.0E-7 or larger

  subroutine  FitCpld_fitting(Nmesh, Mambda, chnid, uptoQn, knwpara, Nknw,     &
    & kphs1_lst, Nphs1, kphs2_lst, Nphs2, kmxa_lst, Nmxa, unknwpara, Nunknw,   &
    & initfrac, tolerance, minchisqr, GetPhaseCpld, Nmaxcall)

    integer, intent(in)       :: Nmesh, chnid, uptoQn
    integer, intent(in)       :: Nknw, Nunknw, Nphs1, Nphs2, Nmxa, Nmaxcall
    real(NER), intent(in)     :: Mambda, initfrac, tolerance
    real(NER), intent(in)     :: knwpara(:), kphs1_lst(:, :), kphs2_lst(:, :), &
      & kmxa_lst(:, :)
    real(NER), intent(inout)  :: unknwpara(:)
    real(NER), intent(out)    :: minchisqr
    procedure(GetCpld_fitpwa) :: GetPhaseCpld

    integer   :: errflag, jj, ivarbl, npar_uptoQn
    real(NER) :: argmn(1:5), parup, parlow, parerr, klist(1:Nphs1+Nphs2+Nmxa)
    character(len = 4)  :: parname(1:MAX_NUMBER_ORDER)
    real(NER) :: paralst(1:Nknw+Nunknw)
    type(triplet_phs_cpld)  :: predct_phase(1:Nphs1+Nphs2+Nmxa, 1:uptoQn+1)

    npar_uptoQn = Nknw + Nunknw
    klist(1:Nphs1) = kphs1_lst(1:Nphs1, 1)
    if (Nphs2 > 0) klist(Nphs1+1:Nphs1+Nphs2) = kphs2_lst(1:Nphs2, 1)
    if (Nmxa > 0) klist(Nphs1+Nphs2+1:Nphs1+Nphs2+Nmxa) = kmxa_lst(1:Nmxa, 1)

    call MNINIT (5, 6, 7)
    call MNSETI ('FitCpld_fitting')
!debug
print *, Nphs1,Nphs2,Nmxa
    argmn(1) = -1.0_NER
    call MNEXCM (chisqr, 'SET PRINTOUT', argmn, 1, errflag, function_zero)

    do jj = 1, Nunknw
      write (parname(jj), '(a, I1.1)') 'g', jj
      call MNPARM (jj, trim(parname(jj)), unknwpara(jj),                       &
        & abs(unknwpara(jj)*initfrac), 0.0_NER, 0.0_NER, errflag)
    end do
    ! Strategy level = 0, 1, 2 (higher number means more derivative computes)
    argmn(1) = 1.0_NER
    call MNEXCM (chisqr, 'SET STRATEGY', argmn, 1, errflag, function_zero)
    ! max number of calls
    argmn(1) = real(Nmaxcall)
    ! tolerance
    argmn(2) = abs(unknwpara(jj)*tolerance)
    call MNEXCM (chisqr, 'MIGRAD', argmn, 2, errflag, function_zero)
    do jj = 1, Nunknw
      call MNPOUT (jj, trim(parname(jj)), unknwpara(jj), parerr, parup,        &
        & parlow, ivarbl)
    end do

  contains

    subroutine chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

      integer, intent(in)     :: NPAR, IFLAG
      real(NER), intent(in)   :: XVAL(1:NPAR)
      real(NER), intent(out)  :: FVAL
      real(NER), intent(out)  :: GRAD(1:NPAR)

      external FUTIL

      integer   :: ii, jj
      real(NER) :: tmpsum

      if (Nknw > 0) paralst(1:Nknw) = knwpara(1:Nknw)
      paralst(Nknw+1:npar_uptoQn) = XVAL(1:NPAR)
      call GetPhaseCpld(Nmesh, Mambda, chnid, uptoQn, paralst, npar_uptoQn,    &
        & klist, Nphs1+Nphs2+Nmxa, predct_phase)

      FVAL = 0.0_NER
      do ii = 1, Nphs1
        tmpsum = 0.0_NER
        do jj = 1, uptoQn+1
          tmpsum = tmpsum + predct_phase(ii, jj)%d1
        end do
!         FVAL = FVAL + (tmpsum - kphs1_lst(ii, 2))**2 /abs(kphs1_lst(ii, 1))
        FVAL = FVAL + (tmpsum - kphs1_lst(ii, 2))**2
      end do
      do ii = 1, Nphs2
        tmpsum = 0.0_NER
        do jj = 1, uptoQn+1
          tmpsum = tmpsum + predct_phase(Nphs1+ii, jj)%d2
        end do
!         FVAL = FVAL + (tmpsum - kphs2_lst(ii, 2))**2 /abs(kphs2_lst(ii, 1))
        FVAL = FVAL + (tmpsum - kphs2_lst(ii, 2))**2
      end do
      do ii = 1, Nmxa
        tmpsum = 0.0_NER
        do jj = 1, uptoQn+1
          tmpsum = tmpsum + predct_phase(Nphs1+Nphs2+ii, jj)%e
        end do
!         FVAL = FVAL + (tmpsum - kmxa_lst(ii, 2))**2 /abs(kmxa_lst(ii, 1))
        FVAL = FVAL + (tmpsum - kmxa_lst(ii, 2))**2
      end do

     end subroutine chisqr

  end subroutine FitCpld_fitting

  subroutine PrintOutFitSngl_fitting(GetPhaseSngl, fitphs_fname, inpar_fname,  &
    & nparlst, out_fname, parname_Qn, initfrac, tolerance, successful)

    procedure(GetSngl_fitpwa)       :: GetPhaseSngl
    character(len = *), intent(in)  :: fitphs_fname, inpar_fname, out_fname
    integer, intent(in)             :: nparlst(:)
    character(len = *), intent(in)  :: parname_Qn(:)
    real(NER), intent(in)           :: initfrac, tolerance
    logical, intent(out)            :: successful

    integer, parameter    :: MAX_NUMK = 100, MAX_NMBD = 50, MAX_NPAR = 30
    integer               :: Nmesh, chnid, numk, nMmbd, ii, funt
    integer               :: uptoQn, npar_known, npar_unknown, setn_cis
    real(NER)             :: kcm_fitphs(MAX_NUMK, 2)
    real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
    real(NER)             :: minchisqr
    character(len = 512)  :: cmdline, timestamp
    character(len = 10)   :: chnstr
    logical               :: succeeded

    successful = .true.
    call ParseCommandLine_fitting(cmdline, Nmesh, chnstr, uptoQn, setn_cis,    &
      & succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if
    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if

    if (uptoQn >= 3) then
      call Select_cis2PE_vtpe1(setn_cis)
      call Select_cisSFR_vtpe1(setn_cis)
      call Select_cis2PE_vdel(setn_cis)
    end if

    call read_fitphs_fitting(fitphs_fname, 1, kcm_fitphs, numk)
    if (uptoQn == 0) then
      npar_known = 0
    else
      npar_known = nparlst(uptoQn)
    end if
    npar_unknown = nparlst(uptoQn+1) - npar_known
    call read_inpar_fitting(inpar_fname, npar_known+npar_unknown,              &
      & Mambda_inpar, nMmbd)

    do ii = 1, nMmbd

      call FitSngl_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,          &
        & Mambda_inpar(ii, 2:npar_known+1), npar_known,                        &
        & kcm_fitphs, numk,                                                    &
        & Mambda_inpar(ii, npar_known+2:npar_known+1+npar_unknown),            &
        &   npar_unknown, initfrac, tolerance,                                 &
        & minchisqr, GetPhaseSngl, 1000)
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
      write (funt, '(a, x, f9.3, x, f8.4)') '# ', kcm_fitphs(ii, 1),           &
        & kcm_fitphs(ii, 2)
    end do

    write (funt, '(a)') trim(parname_Qn(uptoQn+1))
    do ii = 1, nMmbd
      write (funt, '(f7.2, 14(2x, e26.17e3))') Mambda_inpar(ii, 1),            &
        & Mambda_inpar(ii, 2:npar_known+npar_unknown+1)
    end do

    close(funt)

  end subroutine PrintOutFitSngl_fitting


  ! Specializing in invoking deltaful TPEs
  ! Nonperturbative LO
  subroutine  PrintOutFitCpld_fitting(GetPhaseCpld, fitphs1_fname,             &
    & fitphs2_fname,fitmxa_fname, inpar_fname, nparlst, out_fname, parname_Qn, &
    & initfrac, tolerance, successful)

    procedure(GetCpld_fitpwa)       :: GetPhaseCpld
    integer, parameter    :: MAX_NUMK = 100, MAX_NMBD = 50, MAX_NPAR = 30
    integer               :: Nmesh, Nphs1, Nphs2, Nmxa, nMmbd, ii, funt, npar_known
    integer               :: chnid, npar_unknown, uptoQn, setn_cis
    integer               :: cmdcnt, ios
    integer, intent(in)   :: nparlst(:)

    real(NER), intent(in) :: initfrac, tolerance
    real(NER)             :: kcm_fitphs1(MAX_NUMK, 2), kcm_fitphs2(MAX_NUMK, 2)
    real(NER)             :: kcm_fitmxa(MAX_NUMK, 2)
    real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
    real(NER)             :: minchisqr

    character(len = *)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname
    character(len = *)  :: inpar_fname, out_fname
    character(len = 80)   :: timestamp, chnstr
    character(len = 1024) :: cmdbuff, cmdline, tmp
    character(len = *), intent(in)  :: parname_Qn(:)

    logical               :: succeeded
    logical, intent(out)  :: successful

    successful = .true.
    call ParseCommandLine_fitting(cmdline, Nmesh, chnstr, uptoQn, setn_cis,    &
      & succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if
    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if

    if (uptoQn >= 3) then
      if (setn_cis >= 5) then
        call Select_cis2PE_vdel(setn_cis)
        call Select_cisSFR_vtpe1(setn_cis)
        call Select_cis2PE_vtpe1(setn_cis)
      else
        call Select_cis2PE_vdel(5)
        ! select fifth group of cis as initial value for vdel
        call Select_cis2PE_vtpe1(setn_cis)
      end if
    end if
    call read_fitphs_fitting(fitphs1_fname, 1, kcm_fitphs1, Nphs1)
    call read_fitphs_fitting(fitphs2_fname, 1, kcm_fitphs2, Nphs2)
    call read_fitphs_fitting(fitmxa_fname, 1, kcm_fitmxa, Nmxa)

    if (uptoQn == 0) then
      npar_known = 0
    else
      npar_known = nparlst(uptoQn)
    end if
    npar_unknown = nparlst(uptoQn+1) - npar_known
    call read_inpar_fitting(inpar_fname, npar_known+npar_unknown,              &
      & Mambda_inpar, nMmbd)

    do ii = 1, nMmbd
      call FitCpld_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,          &
        & Mambda_inpar(ii, 2:npar_known+1), npar_known,                        &
        & kcm_fitphs1, Nphs1, kcm_fitphs2, Nphs2, kcm_fitmxa, Nmxa,            &
        & Mambda_inpar(ii, npar_known+2:npar_known+1+npar_unknown),            &
        &   npar_unknown, initfrac, tolerance,                                 &
        & minchisqr, GetPhaseCpld, MAXNUM_CALLMINUIT)
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

  end subroutine PrintOutFitCpld_fitting

  ! fitting amplitude
  ! k_phase(i, 1) = the ith kcm, k_phase(i, 2) = the ith phase
  ! initfrac: initial step = initfrac*initvalue, advised value: 1.0E-3
  ! tolerance: advised value: 1.0E-7 or larger
  ! Nmaxcall: max number of calls ~ 1000

  subroutine  FitSnglTmtrx_fitting(Nmesh, Mambda, chnid, uptoQn, knwpara, Nknw,     &
    & k_phase, numk, unknwpara, Nunknw, initfrac, tolerance, minchisqr,       &
    & GetTmtrxSngl, Nmaxcall)

    integer, intent(in)       :: Nmesh, chnid, uptoQn, numk, Nknw, Nunknw
    integer, intent(in)       :: Nmaxcall
    real(NER), intent(in)     :: Mambda, initfrac, tolerance
    real(NER), intent(in)     :: knwpara(:), k_phase(:, :)
    real(NER), intent(inout)  :: unknwpara(:)
    real(NER), intent(out)    :: minchisqr
    procedure(GetSnglTmtrx_fitpwa) :: GetTmtrxSngl

    ! npar_uptoQn: total number of inparas up to this order
    integer   :: errflag, jj, ivarbl, npar_uptoQn
    real(NER) :: argmn(1:5), parup, parlow, parerr
    character(len = 4)  :: parname(1:MAX_NUMBER_ORDER)
    real(NER) :: paralst(1:Nknw + Nunknw)
    real(NER) :: predct_phase(1:numk, 1:uptoQn+1)

    npar_uptoQn = Nknw + Nunknw

    call MNINIT (5, 6, 7)
    call MNSETI ('FitSnglTmtrx_fitting')

    ! -1 no output except from SHOW commands
    ! 0 minimum output (no starting values or intermediate results)
    ! 1 default value, normal output
    ! 2 additional output giving intermediate results.
    ! 3 maximum output,showing progress of minimizations.
    argmn(1) = -1.0_NER
    call MNEXCM (chisqr, 'SET PRINTOUT', argmn, 1, errflag, function_zero)

    do jj = 1, Nunknw
      write (parname(jj), '(a, I1.1)') 'g', jj
      ! initial steps: abs(unknwpara(jj)*initfrac
      call MNPARM (jj, trim(parname(jj)), unknwpara(jj),                       &
        & abs(unknwpara(jj)*initfrac), 0.0_NER, 0.0_NER, errflag)
    end do
    ! Strategy level = 0, 1, 2 (higher number means more derivative computes)
    argmn(1) = 1.0_NER
    call MNEXCM (chisqr, 'SET STRATEGY', argmn, 1, errflag, function_zero)
    ! max number of calls
    argmn(1) = real(Nmaxcall)
    ! tolerance
    argmn(2) = abs(unknwpara(jj)*tolerance)
    call MNEXCM (chisqr, 'MIGRAD', argmn, 2, errflag, function_zero)
    do jj = 1, Nunknw
      call MNPOUT (jj, trim(parname(jj)), unknwpara(jj), parerr, parup, parlow, &
        & ivarbl)
    end do

  contains

    subroutine chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

      integer, intent(in)     :: NPAR, IFLAG
      real(NER), intent(in)   :: XVAL(1:NPAR)
      real(NER), intent(out)  :: FVAL
      real(NER), intent(out)  :: GRAD(1:NPAR)

      complex(NER),parameter  :: im = (0, 1)
      real(NER),parameter     :: empirical_scatlength = -23.740_NER
      complex(NER)            :: predct(1: numk), tmtrx(1: numk)
      real(NER)               :: phase1(1: numk)

      external FUTIL

      integer   :: ii

      if (Nknw > 0) paralst(1:Nknw) = knwpara(1:Nknw)
      paralst(Nknw+1:npar_uptoQn) = XVAL(1:NPAR)
      call GetTmtrxSngl(Nmesh, Mambda, chnid, uptoQn, paralst, npar_uptoQn,    &
        & k_phase(1:numk, 1), numk, predct_phase)

      phase1(1) = sum(predct_phase(1, 1:uptoQn+1))
      phase1(1) = phase1(1) / 180.0_NER * PI_NE
      predct(1) = -sin(phase1(1))*exp(im*phase1(1))/k_phase(1, 1)
      FVAL = (real(predct(1))*PC_default_hbarc - empirical_scatlength)**2
      do ii = 2, numk
        phase1(ii) = sum(predct_phase(ii, 1:uptoQn+1))
        FVAL = FVAL + abs((phase1(ii) - k_phase(ii,2)))**2
      end do
    end subroutine chisqr

  end subroutine FitSnglTmtrx_fitting



  subroutine PrintOutFitSnglTmtrx_fitting(GetTmtrxSngl, fitphs_fname, inpar_fname,  &
    & nparlst, out_fname, parname_Qn, initfrac, tolerance, successful)

    procedure(GetSnglTmtrx_fitpwa)       :: GetTmtrxSngl
    character(len = *), intent(in)  :: fitphs_fname, inpar_fname, out_fname
    integer, intent(in)             :: nparlst(:)
    character(len = *), intent(in)  :: parname_Qn(:)
    real(NER), intent(in)           :: initfrac, tolerance
    logical, intent(out)            :: successful

    integer, parameter    :: MAX_NUMK = 100, MAX_NMBD = 50, MAX_NPAR = 30
    integer               :: Nmesh, chnid, numk, nMmbd, ii, funt
    integer               :: uptoQn, npar_known, npar_unknown, setn_cis
    real(NER)             :: kcm_fittmtrx(MAX_NUMK, 2)
    real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
    real(NER)             :: minchisqr
    character(len = 512)  :: cmdline, timestamp
    character(len = 10)   :: chnstr
    logical               :: succeeded

    successful = .true.
    call ParseCommandLine_fitting(cmdline, Nmesh, chnstr, uptoQn, setn_cis,    &
      & succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if
    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if

    if (uptoQn >= 3) then
      call Select_cis2PE_vtpe1(setn_cis)
      call Select_cisSFR_vtpe1(setn_cis)
      call Select_cis2PE_vdel(setn_cis)
    end if

    call read_fitphs_fitting(fitphs_fname, 1, kcm_fittmtrx, numk)
    if (uptoQn == 0) then
      npar_known = 0
    else
      npar_known = nparlst(uptoQn)
    end if
    npar_unknown = nparlst(uptoQn+1) - npar_known
    call read_inpar_fitting(inpar_fname, npar_known+npar_unknown,              &
      & Mambda_inpar, nMmbd)

    do ii = 1, nMmbd

      call FitSnglTmtrx_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,          &
        & Mambda_inpar(ii, 2:npar_known+1), npar_known,                        &
        & kcm_fittmtrx, numk,                                                    &
        & Mambda_inpar(ii, npar_known+2:npar_known+1+npar_unknown),            &
        &   npar_unknown, initfrac, tolerance,                                 &
        & minchisqr, GetTmtrxSngl, MAXNUM_CALLMINUIT)
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
      write (funt, '(a, x, f9.3, x, f8.4)') '# ', kcm_fittmtrx(ii, 1),           &
        & kcm_fittmtrx(ii, 2)
    end do

    write (funt, '(a)') trim(parname_Qn(uptoQn+1))
    do ii = 1, nMmbd
      write (funt, '(f7.2, 14(2x, e26.17e3))') Mambda_inpar(ii, 1),            &
        & Mambda_inpar(ii, 2:npar_known+npar_unknown+1)
    end do

    close(funt)

  end subroutine PrintOutFitSnglTmtrx_fitting

! merge cpld and sngl together
! fit 1s0 and 3s1 with scattering length and binding energy repectively
! choose the fitting data by chnid
! SLBEPS = Scattering Length Binding Energy Phase Shifts
! pray 25/7/2020 add for mmwly PC

  subroutine  Fit_SLBEPS_fitting(Nmesh, Mambda, chnid, uptoQn, knwpara, Nknw,     &
    & kphs1_lst, Nphs1, kphs2_lst, Nphs2, kmxa_lst, Nmxa, unknwpara, Nunknw,   &
    & initfrac, tolerance, minchisqr, GetSLBEPSfitpwa, Nmaxcall)

    integer, intent(in)       :: Nmesh, chnid, uptoQn
    integer, intent(in)       :: Nknw, Nunknw, Nphs1, Nphs2, Nmxa, Nmaxcall
    real(NER), intent(in)     :: Mambda, initfrac, tolerance
    real(NER), intent(in)     :: knwpara(:), kphs1_lst(:, :), kphs2_lst(:, :), &
      & kmxa_lst(:, :)
    real(NER), intent(inout)  :: unknwpara(:)
    real(NER), intent(out)    :: minchisqr
    procedure(GetSLBEPS_fitpwa) :: GetSLBEPSfitpwa

    integer   :: errflag, jj, ivarbl, npar_uptoQn
    real(NER) :: argmn(1:5), parup, parlow, parerr, klist(1:Nphs1+Nphs2+Nmxa)
    character(len = 4)  :: parname(1:MAX_NUMBER_ORDER)
    real(NER) :: paralst(1:Nknw+Nunknw)
    real(NER) :: predct_phase(1:Nphs1*3+Nphs2*3+Nmxa*3, 1:uptoQn+1)
    real(NER)             :: predct_Tmtrx(1:Nphs1+Nphs2+Nmxa, 1:uptoQn+1)
    character(len = 10)   :: chnstr
    logical               :: succeeded
    real(NER)             :: kfor_sc(1:1), prphsc(1:3,1:uptoQn+1), prTmxsc(1:1,1:uptoQn+1)
    type(lsj_symbol)      :: lsj
    real(NER)             :: predct_SC, predct_BD

    kfor_sc(1:1) = 1E-9_NER
    npar_uptoQn = Nknw + Nunknw
    klist(1:Nphs1) = kphs1_lst(1:Nphs1, 1)
    if (Nphs2 > 0) klist(Nphs1+1:Nphs1+Nphs2) = kphs2_lst(1:Nphs2, 1)
    if (Nmxa > 0) klist(Nphs1+Nphs2+1:Nphs1+Nphs2+Nmxa) = kmxa_lst(1:Nmxa, 1)

    call MNINIT (5, 6, 7)
    call MNSETI ('Fit_SLBEPS_fitting')

    argmn(1) = -1.0_NER
    call MNEXCM (chisqr, 'SET PRINTOUT', argmn, 1, errflag, function_zero)

    do jj = 1, Nunknw
      write (parname(jj), '(a, I1.1)') 'g', jj
      call MNPARM (jj, trim(parname(jj)), unknwpara(jj),                       &
        & abs(unknwpara(jj)*initfrac), 0.0_NER, 0.0_NER, errflag)
    end do
    ! Strategy level = 0, 1, 2 (higher number means more derivative computes)
    argmn(1) = 1.0_NER
    call MNEXCM (chisqr, 'SET STRATEGY', argmn, 1, errflag, function_zero)
    ! max number of calls
    argmn(1) = real(Nmaxcall)
    ! tolerance
    argmn(2) = abs(unknwpara(jj)*tolerance)
    call MNEXCM (chisqr, 'MIGRAD', argmn, 2, errflag, function_zero)
    do jj = 1, Nunknw
      call MNPOUT (jj, trim(parname(jj)), unknwpara(jj), parerr, parup,        &
        & parlow, ivarbl)
    end do

  contains

    subroutine chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

      integer, intent(in)     :: NPAR, IFLAG
      real(NER), intent(in)   :: XVAL(1:NPAR)
      real(NER), intent(out)  :: FVAL
      real(NER), intent(out)  :: GRAD(1:NPAR)

      complex(NER),parameter  :: im = (0, 1)
      real(NER),parameter     :: empirical_scatlength = -23.740_NER

      external FUTIL

      integer   :: ii, jj

      if (Nknw > 0) paralst(1:Nknw) = knwpara(1:Nknw)
      paralst(Nknw+1:npar_uptoQn) = XVAL(1:NPAR)
      call GetSLBEPSfitpwa(Nmesh, Mambda, chnid, uptoQn, paralst, npar_uptoQn,    &
        & klist, Nphs1+Nphs2+Nmxa, predct_phase, predct_Tmtrx)

      FVAL = 0.0_NER
      call convert_chnid_to_text(chnid, chnstr, succeeded)
      call convert_chnid_to_lsj(chnid, lsj, succeeded)

      select case(trim(chnstr))
        case('1s0')
          call GetSLBEPSfitpwa(Nmesh, Mambda, chnid, uptoQn, paralst, npar_uptoQn,    &
            & kfor_sc, 1, prphsc, prTmxsc)
          predct_SC = real(sum(prTmxsc(1,1:uptoQn+1)))*PI_NE/2.0_NER*PC_default_hbarc
!           debug
!             print*, "predct_SC = ", predct_SC
          FVAL = (real(predct_SC) - empirical_scatlength)**2
!         case('3s1')
!           get 3s1 binding energy and FVAL
        case default
          FVAL = 0.0_NER
      end select
      if (lsj%L /= lsj%j .and. lsj%j /= 0) then  !cpld
        do ii = 1, Nphs1
          FVAL = FVAL + (sum(predct_phase(ii, 1:uptoQn+1))-kphs1_lst(ii, 2))**2
        end do
        do ii = 1, Nphs2
          FVAL = FVAL + (sum(predct_phase(Nphs1+Nphs2+Nmxa+NPhs1+ii, 1:uptoQn+1))-kphs2_lst(ii, 2))**2
        end do
        do ii = 1, Nmxa
          FVAL = FVAL + (sum(predct_phase(Nphs1*2+Nphs2*2+Nmxa*2+Nphs2+NPhs1+ii, 1:uptoQn+1))-kmxa_lst(ii, 2))**2
        end do
      else
        do ii = 1, Nphs1
          FVAL = FVAL + (sum(predct_phase(ii, 1:uptoQn+1))-kphs1_lst(ii, 2))**2
        end do
      end if
     end subroutine chisqr
  end subroutine Fit_SLBEPS_fitting

  subroutine  PrintOutFitSLBEPS_fitting(GetSLBEPSfitpwa, fitphs1_fname,             &
    & fitphs2_fname,fitmxa_fname, inpar_fname, nparlst, out_fname, parname_Qn, &
    & initfrac, tolerance, successful)

    procedure(GetSLBEPS_fitpwa) :: GetSLBEPSfitpwa
    integer, parameter    :: MAX_NUMK = 100, MAX_NMBD = 50, MAX_NPAR = 30
    integer               :: Nmesh, Nphs1, Nphs2, Nmxa, nMmbd, ii, funt, npar_known
    integer               :: chnid, npar_unknown, uptoQn, setn_cis
    integer               :: cmdcnt, ios
    integer, intent(in)   :: nparlst(:)

    real(NER), intent(in) :: initfrac, tolerance
    real(NER)             :: kcm_fitphs1(MAX_NUMK, 2), kcm_fitphs2(MAX_NUMK, 2)
    real(NER)             :: kcm_fitmxa(MAX_NUMK, 2)
    real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
    real(NER)             :: minchisqr

    character(len = *)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname
    character(len = *)  :: inpar_fname, out_fname
    character(len = 80)   :: timestamp, chnstr
    character(len = 1024) :: cmdbuff, cmdline, tmp
    character(len = *), intent(in)  :: parname_Qn(:)

    logical               :: succeeded
    logical, intent(out)  :: successful

    successful = .true.
    call ParseCommandLine_fitting(cmdline, Nmesh, chnstr, uptoQn, setn_cis,    &
      & succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if
    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if
    if (uptoQn >= 3) then
      if (setn_cis >= 5) then
        call Select_cis2PE_vdel(setn_cis)
        call Select_cisSFR_vtpe1(setn_cis)
        call Select_cis2PE_vtpe1(setn_cis)
      else
        call Select_cis2PE_vdel(5)
        ! select fifth group of cis as initial value for vdel
        call Select_cis2PE_vtpe1(setn_cis)
      end if
    end if
    if (trim(fitphs2_fname) == '') then
      call read_fitphs_fitting(fitphs1_fname, 1, kcm_fitphs1, Nphs1)
      kcm_fitmxa = 0.0_NER
      kcm_fitphs2 = kcm_fitmxa
      Nmxa = 0
      Nphs2 = Nmxa
    else
      call read_fitphs_fitting(fitphs1_fname, 1, kcm_fitphs1, Nphs1)
      call read_fitphs_fitting(fitphs2_fname, 1, kcm_fitphs2, Nphs2)
      call read_fitphs_fitting(fitmxa_fname, 1, kcm_fitmxa, Nmxa)
    end if
    if (uptoQn == 0) then
      npar_known = 0
    else
      npar_known = nparlst(uptoQn)
    end if
    npar_unknown = nparlst(uptoQn+1) - npar_known
    call read_inpar_fitting(inpar_fname, npar_known+npar_unknown,              &
      & Mambda_inpar, nMmbd)

    do ii = 1, nMmbd
      call Fit_SLBEPS_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,          &
        & Mambda_inpar(ii, 2:npar_known+1), npar_known,                        &
        & kcm_fitphs1, Nphs1, kcm_fitphs2, Nphs2, kcm_fitmxa, Nmxa,            &
        & Mambda_inpar(ii, npar_known+2:npar_known+1+npar_unknown),            &
        &   npar_unknown, initfrac, tolerance,                                 &
        & minchisqr, GetSLBEPSfitpwa, MAXNUM_CALLMINUIT)
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
    if (trim(fitphs2_fname) /= '') then
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
    end if
    write (funt, '(a)') trim(parname_Qn(uptoQn+1))
    do ii = 1, nMmbd
      write (funt, '(f7.2, 14(2x, e26.17e3))') Mambda_inpar(ii, 1),              &
        & Mambda_inpar(ii, 2:npar_known+npar_unknown+1)
    end do

    close(funt)

  end subroutine PrintOutFitSLBEPS_fitting

  subroutine  Fit3s1_BD_fitting(Nmesh, Mambda, chnid, uptoQn, knwpara, Nknw,     &
    & kphs1_lst, Nphs1, kphs2_lst, Nphs2, kmxa_lst, Nmxa, unknwpara, Nunknw,   &
    & initfrac, tolerance, minchisqr, Get3s1BD_Cpld, Nmaxcall)

    integer, intent(in)       :: Nmesh, chnid, uptoQn
    integer, intent(in)       :: Nknw, Nunknw, Nphs1, Nphs2, Nmxa, Nmaxcall
    real(NER), intent(in)     :: Mambda, initfrac, tolerance
    real(NER), intent(in)     :: kphs1_lst(:, :), kphs2_lst(:, :), &
      & kmxa_lst(:, :)
    real(NER), intent(inout)  :: unknwpara(:), knwpara(:)
    real(NER), intent(out)    :: minchisqr
    procedure(Get3s1BD_fitpwa) :: Get3s1BD_Cpld

    integer   :: errflag, jj, ivarbl, npar_uptoQn
    real(NER) :: argmn(1:5), parup, parlow, parerr, klist(1:Nphs1+Nphs2+Nmxa)
    character(len = 4)  :: parname(1:MAX_NUMBER_ORDER)
    real(NER) :: paralst(1:Nknw+Nunknw), paralstout(1:Nknw+Nunknw)
    type(triplet_phs_cpld)  :: predct_phase(1:Nphs1+Nphs2+Nmxa, 1:uptoQn+1)

    npar_uptoQn = Nknw + Nunknw
    klist(1:Nphs1) = kphs1_lst(1:Nphs1, 1)
    if (Nphs2 > 0) klist(Nphs1+1:Nphs1+Nphs2) = kphs2_lst(1:Nphs2, 1)
    if (Nmxa > 0) klist(Nphs1+Nphs2+1:Nphs1+Nphs2+Nmxa) = kmxa_lst(1:Nmxa, 1)

    call MNINIT (5, 6, 7)
    call MNSETI ('Fit3s1_BD_fitting')
!debug
print *, 'Nphs1,Nphs2,Nmxa =', Nphs1,Nphs2,Nmxa
    argmn(1) = -1.0_NER
    call MNEXCM (chisqr, 'SET PRINTOUT', argmn, 1, errflag, function_zero)

    do jj = 1, Nunknw
      write (parname(jj), '(a, I1.1)') 'g', jj
      call MNPARM (jj, trim(parname(jj)), unknwpara(jj),                       &
        & abs(unknwpara(jj)*initfrac), 0.0_NER, 0.0_NER, errflag)
    end do
    ! Strategy level = 0, 1, 2 (higher number means more derivative computes)
    argmn(1) = 1.0_NER
    call MNEXCM (chisqr, 'SET STRATEGY', argmn, 1, errflag, function_zero)
    ! max number of calls
    argmn(1) = real(Nmaxcall)
    ! tolerance
    argmn(2) = abs(unknwpara(jj)*tolerance)
    call MNEXCM (chisqr, 'MIGRAD', argmn, 2, errflag, function_zero)
    do jj = 1, Nunknw
      call MNPOUT (jj, trim(parname(jj)), unknwpara(jj), parerr, parup,        &
        & parlow, ivarbl)
    end do

  contains

    subroutine chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

      integer, intent(in)     :: NPAR, IFLAG
      real(NER), intent(in)   :: XVAL(1:NPAR)
      real(NER), intent(out)  :: FVAL
      real(NER), intent(out)  :: GRAD(1:NPAR)

      external FUTIL

      integer   :: ii, jj
      real(NER) :: tmpsum

      if (Nknw > 0) paralst(1:Nknw) = knwpara(1:Nknw)
      paralst(Nknw+1:npar_uptoQn) = XVAL(1:NPAR)
      call Get3s1BD_Cpld(Nmesh, Mambda, chnid, uptoQn, paralst, npar_uptoQn,    &
        & klist, Nphs1+Nphs2+Nmxa, predct_phase, paralstout)
!       debug
!       print*,'paralstout = ', paralstout(1:4)
!       knwpara(Nknw) = paralstout(Nknw)

      FVAL = 0.0_NER
      do ii = 1, Nphs1
        tmpsum = 0.0_NER
        do jj = 1, uptoQn+1
          tmpsum = tmpsum + predct_phase(ii, jj)%d1
        end do
!         FVAL = FVAL + (tmpsum - kphs1_lst(ii, 2))**2 /abs(kphs1_lst(ii, 1))
        FVAL = FVAL + (tmpsum - kphs1_lst(ii, 2))**2
      end do
      do ii = 1, Nphs2
        tmpsum = 0.0_NER
        do jj = 1, uptoQn+1
          tmpsum = tmpsum + predct_phase(Nphs1+ii, jj)%d2
        end do
!         FVAL = FVAL + (tmpsum - kphs2_lst(ii, 2))**2 /abs(kphs2_lst(ii, 1))
        FVAL = FVAL + (tmpsum - kphs2_lst(ii, 2))**2
      end do
      do ii = 1, Nmxa
        tmpsum = 0.0_NER
        do jj = 1, uptoQn+1
          tmpsum = tmpsum + predct_phase(Nphs1+Nphs2+ii, jj)%e
        end do
!         FVAL = FVAL + (tmpsum - kmxa_lst(ii, 2))**2 /abs(kmxa_lst(ii, 1))
        FVAL = FVAL + (tmpsum - kmxa_lst(ii, 2))**2
      end do

     end subroutine chisqr

  end subroutine Fit3s1_BD_fitting

  subroutine  PrintOutFit3s1_BD_fitting(Get3s1BD_Cpld, fitphs1_fname,             &
    & fitphs2_fname,fitmxa_fname, inpar_fname, nparlst, out_fname, parname_Qn, &
    & initfrac, tolerance, successful)

    procedure(Get3s1BD_fitpwa)       :: Get3s1BD_Cpld
    integer, parameter    :: MAX_NUMK = 100, MAX_NMBD = 50, MAX_NPAR = 30
    integer               :: Nmesh, Nphs1, Nphs2, Nmxa, nMmbd, ii, funt, npar_known
    integer               :: chnid, npar_unknown, uptoQn, setn_cis
    integer               :: cmdcnt, ios
    integer, intent(in)   :: nparlst(:)

    real(NER), intent(in) :: initfrac, tolerance
    real(NER)             :: kcm_fitphs1(MAX_NUMK, 2), kcm_fitphs2(MAX_NUMK, 2)
    real(NER)             :: kcm_fitmxa(MAX_NUMK, 2)
    real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
    real(NER)             :: minchisqr, Mambda_intext

    character(len = *)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname
    character(len = *)  :: inpar_fname, out_fname
    character(len = 80)   :: timestamp, chnstr
    character(len = 1024) :: cmdbuff, cmdline, tmp
    character(len = *), intent(in)  :: parname_Qn(:)

    logical               :: succeeded
    logical, intent(out)  :: successful

    successful = .true.
    call ParseCommandLine_fitting(cmdline, Nmesh, chnstr, uptoQn, setn_cis,    &
      & succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if
    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if

    if (uptoQn >= 3) then
      if (setn_cis >= 5) then
        call Select_cis2PE_vdel(setn_cis)
        call Select_cisSFR_vtpe1(setn_cis)
        call Select_cis2PE_vtpe1(setn_cis)
      else
        call Select_cis2PE_vdel(5)
        ! select fifth group of cis as initial value for vdel
        call Select_cis2PE_vtpe1(setn_cis)
      end if
    end if
    call read_fitphs_fitting(fitphs1_fname, 1, kcm_fitphs1, Nphs1)
    call read_fitphs_fitting(fitphs2_fname, 1, kcm_fitphs2, Nphs2)
    call read_fitphs_fitting(fitmxa_fname, 1, kcm_fitmxa, Nmxa)

    if (uptoQn == 0) then
      npar_known = 0
    else
      npar_known = nparlst(uptoQn)
    end if
    npar_unknown = nparlst(uptoQn+1) - npar_known
    call read_inpar_fitting(inpar_fname, npar_known+npar_unknown,              &
      & Mambda_inpar, nMmbd)

    do ii = 1, nMmbd
      call Fit3s1_BD_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,          &
        & Mambda_inpar(ii, 2:npar_known+2), npar_known+1,                        &
        & kcm_fitphs1, Nphs1, kcm_fitphs2, Nphs2, kcm_fitmxa, Nmxa,            &
        & Mambda_inpar(ii, npar_known+3:npar_known+1+npar_unknown),            &
        &   npar_unknown-1, initfrac, tolerance,                                 &
        & minchisqr, Get3s1BD_Cpld, MAXNUM_CALLMINUIT)
      Mambda_intext = Mambda_inpar(ii, 1)
!       Mambda_inpar(ii, 3) = -(allpert_BDE(1) + Mambda_intext**2*Mambda_inpar(ii, 4)*allpert_BDE(3) +  &
!         & Mambda_intext**2*Mambda_inpar(ii, 5)*allpert_BDE(4))/(allpert_BDE(2))
! debug fit without TPE
      Mambda_inpar(ii, 3) = -(Mambda_intext**2*Mambda_inpar(ii, 4)*allpert_BDE(3) +  &
        & Mambda_intext**2*Mambda_inpar(ii, 5)*allpert_BDE(4))/(allpert_BDE(2))
      allpert_BDE = 0.0_NER
      allpert_BDE_initial = .false.
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

  end subroutine PrintOutFit3s1_BD_fitting


  subroutine PrintOutFit1s0_SCL_fitting(GetPhaseSngl, fitphs_fname, inpar_fname,  &
    & nparlst, out_fname, parname_Qn, initfrac, tolerance, successful)

    procedure(GetSngl_fitpwa)       :: GetPhaseSngl
    character(len = *), intent(in)  :: fitphs_fname, inpar_fname, out_fname
    integer, intent(in)             :: nparlst(:)
    character(len = *), intent(in)  :: parname_Qn(:)
    real(NER), intent(in)           :: initfrac, tolerance
    logical, intent(out)            :: successful

    integer, parameter    :: MAX_NUMK = 100, MAX_NMBD = 50, MAX_NPAR = 30
    integer               :: Nmesh, chnid, numk, nMmbd, ii, funt
    integer               :: uptoQn, npar_known, npar_unknown, setn_cis
    real(NER)             :: kcm_fitphs(MAX_NUMK, 2)
    real(NER)             :: Mambda(MAX_NMBD), Mambda_inpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: knownpar(MAX_NMBD, MAX_NPAR)
    real(NER)             :: unknpara(MAX_NMBD, MAX_NPAR)
    real(NER)             :: minchisqr, Mambda_intext
    character(len = 512)  :: cmdline, timestamp
    character(len = 10)   :: chnstr
    logical               :: succeeded

    successful = .true.
    call ParseCommandLine_fitting(cmdline, Nmesh, chnstr, uptoQn, setn_cis,    &
      & succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if
    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (.not. succeeded) then
      successful = .false.
      return
    end if

    if (uptoQn >= 3) then
      call Select_cis2PE_vtpe1(setn_cis)
      call Select_cisSFR_vtpe1(setn_cis)
      call Select_cis2PE_vdel(setn_cis)
    end if

    call read_fitphs_fitting(fitphs_fname, 1, kcm_fitphs, numk)
    if (uptoQn == 0) then
      npar_known = 0
    else
      npar_known = nparlst(uptoQn)
    end if
    npar_unknown = nparlst(uptoQn+1) - npar_known
    call read_inpar_fitting(inpar_fname, npar_known+npar_unknown,              &
      & Mambda_inpar, nMmbd)
! debug
    print*, "nparlst =", nparlst(1:4), '////', npar_known, npar_unknown
    do ii = 1, nMmbd

      call FitSngl_fitting(Nmesh, Mambda_inpar(ii, 1), chnid, uptoQn,          &
        & Mambda_inpar(ii, 2:npar_known+2), npar_known+1,                        &
        & kcm_fitphs, numk,                                                    &
        & Mambda_inpar(ii, npar_known+3:npar_known+1+npar_unknown),            &
        &   npar_unknown-1, initfrac, tolerance,                                 &
        & minchisqr, GetPhaseSngl, 1000)
        Mambda_intext = Mambda_inpar(ii, 1)
        if (uptoQn == 1) then
          Mambda_inpar(ii, 3) = -NLOpert_SCL(2)*Mambda_inpar(ii, 4)*(Mambda_intext**2)/NLOpert_SCL(1)
          NLOpert_SCL = 0.0_NER
          pert_SCL_initial(1) = .false.
        end if
        if (uptoQn == 2) then
          Mambda_inpar(ii, 5) = -(N2LOpert_SCL(3)*Mambda_inpar(ii, 6)*(Mambda_intext**2) + N2LOpert_SCL(4)*Mambda_inpar(ii, 7)*(Mambda_intext**4) &
          & + N2LOpert_SCL(5))/N2LOpert_SCL(2)/Mambda_intext
! debug fit without TPE
!           Mambda_inpar(ii, 5) = -(N2LOpert_SCL(1) + N2LOpert_SCL(3)*Mambda_inpar(ii, 6)*(Mambda_intext**2) + N2LOpert_SCL(4)*Mambda_inpar(ii, 7)*(Mambda_intext**4) &
!           & + N2LOpert_SCL(5))/N2LOpert_SCL(2)/Mambda_intext
          N2LOpert_SCL = 0.0_NER
          pert_SCL_initial(2) = .false.
        end if
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
      write (funt, '(a, x, f9.3, x, f8.4)') '# ', kcm_fitphs(ii, 1),           &
        & kcm_fitphs(ii, 2)
    end do

    write (funt, '(a)') trim(parname_Qn(uptoQn+1))
    do ii = 1, nMmbd
      write (funt, '(f7.2, 14(2x, e26.17e3))') Mambda_inpar(ii, 1),            &
        & Mambda_inpar(ii, 2:npar_known+npar_unknown+1)
    end do

    close(funt)

  end subroutine PrintOutFit1s0_SCL_fitting



end module drv_fitting

