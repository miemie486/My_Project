! pray 05/29/2020
!this program fitting KSW in 3s1 and 1s0

program fit_ksw_smpl

  use drv_fitting
  use testwrap
  implicit none


  integer               :: npar_Qn(1:10), ios, cmdcnt, ii, N
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname
  character(len = 512)  :: fitphs1_fname, fitphs2_fname, fitmxa_fname, errmsg
  character(len = 512)  :: inpar_fname, out_fname
  character(len = 3)          :: cmdstr
  character(len = 1024)       :: cmdbuff, tmp
  logical               :: succeeded
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0  C3  D1                                                  '  &
    & /)
  integer :: regtype
  class(obj_smplchn), allocatable   :: chn
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
    case ('1s0')
    npar_Qn(1:5) = (/ 1, 3, 4, 5, 6/)
    fitphs_fname = 'fit1s0_phs.in'
    inpar_fname = 'fit1s0_inpar.in'
    out_fname = 'fit1s0_result.out'
!   npar_Qn(1:4) = npar_withc0_sngl(1:4)
    initial_fraction = 1.0E-3_NER
    tolerance = 1.0E-7_NER
    errmsg = 'fit_ksw_smpl: Something went wrong'

  call init_smpl
  call PrintOutFitSngl_fitting(GetNPSngl_smpl_fitDR,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)
  if (.not. succeeded) then
    write(standard_error_unit, '(a)') errmsg
  else
    call chn%erase()
    deallocate(chn)
    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if
case default
      npar_Qn(1:5) = (/ 1, 3, 4, 5, 6/)
      fitphs1_fname = 'fit3s1_phs1.in'
      fitphs2_fname = 'fit3s1_phs2.in'
      fitmxa_fname = 'fit3s1_mxa.in'
      inpar_fname = 'fit3s1_inpar.in'
      out_fname = 'fit3s1_result.out'
      initial_fraction = 1.0E-4_NER
      tolerance = 1.0E-5_NER
      errmsg = 'fit_ksw_smpl: Something went wrong'

  call init_smpl()

  call PrintOutFitCpld_fitting(GetNPcpld_smpl_fitDR, fitphs1_fname,      &
    & fitphs2_fname, fitmxa_fname, inpar_fname, npar_Qn, out_fname,            &
    & parname_Qn, initial_fraction, tolerance, succeeded)
  if (.not.succeeded) then
    write(standard_error_unit, '(a)') trim(errmsg)
  else
    call chn%erase()
    deallocate(chn)
    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if
end select
contains

  subroutine init_smpl()

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_smplchn::chn)
    ! NLO pot vanishes
    chn%vanishing(0) = .false.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:5) = (/1, 3, 4, 5, 6/)
    else
      chn%npar(1:5) = (/1, 3, 4, 5, 6/)
    end if

    chn%pVfunc_V0 => VLO_ksw_smpl_epwrap


    chn%pVfunc_V1 => VNLO_ksw_smpl_epwrap
    chn%V1_cheb_base_n = 2
    chn%V1_cheb_bnd = 0.005_NER

    chn%pVfunc_V2 => VN2LO_ksw_smpl_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V3 => VN3LO_ksw_smpl_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V4 => VN4LO_ksw_smpl_epwrap
    chn%V4_cheb_base_n = 4
    chn%V4_cheb_bnd = 1.0E-3_NER

  end subroutine init_smpl

  subroutine GetNPSngl_smpl_fitDR(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    use mod_obj_smplchn
    use drv_pwphase
    use mod_vfunc_smpl
    implicit none

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

!     integer :: regtype
!     class(obj_smplchn), pointer :: chn

!     regtype = REGTYPE_GAUSSIAN
    call create_smplchn(chn, chnid, uptoQn)
!     call SetUpVfunc_withc0_smpl(chn)

    call get_phase_as_value_for_klist_smpl(chn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

  end subroutine GetNPSngl_smpl_fitDR

  subroutine GetNPcpld_smpl_fitDR(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

    use mod_obj_smplchn
    use drv_pwphase
    use mod_vfunc_smpl
    implicit none

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

!     integer :: regtype
!     class(obj_smplchn), pointer :: chn

!     regtype = REGTYPE_GAUSSIAN
    call create_smplchn(chn, chnid, uptoQn)
!     call SetUpVfunc_withc0_smpl(chn)

    call get_phase_as_value_for_klist_smpl(chn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

  end subroutine GetNPcpld_smpl_fitDR

end program fit_ksw_smpl