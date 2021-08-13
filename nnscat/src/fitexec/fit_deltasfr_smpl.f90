! 9/26/2019
!for pope  using smpl  both sngl and cpld chn

program fit_deltasfr_smpl

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
    case ('1d2', '3d2', '3d3')
      npar_Qn(1:4) = (/ 0, 0, 0, 1/)
    case ('1p1', '3p1', '3p0')
      npar_Qn(1:4) = (/ 0, 0, 1, 3/)
    case ('3p2')
      npar_Qn(1:4) = (/ 0, 0, 1, 4/)
    case default
      npar_Qn(1:4) = (/ 0, 0, 0, 0/)
  end select
  call init_smpl
select case(trim(cmdstr))
  case('3p2', '3d3')
    fitphs1_fname = 'fitcpld_phs1.in'
    fitphs2_fname = 'fitcpld_phs2.in'
    fitmxa_fname = 'fitcpld_mxa.in'
    inpar_fname = 'fitcpld_inpar.in'
    out_fname = 'fitcpld_result.out'
!   npar_Qn(1:4) =  chn%npar(1:4)
    initial_fraction = 1.0E-4_NER
    tolerance = 1.0E-5_NER
    errmsg = 'fit_deltasfr_cpld: Something went wrong'


    call PrintOutFitCpld_fitting(GetNPCpldSfr_smpl_fitdeltacpld, fitphs1_fname,      &
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
  case default
    fitphs_fname = 'fitsngl_phs.in'
    inpar_fname = 'fitsngl_inpar.in'
    out_fname = 'fitsngl_result.out'
!   npar_Qn(1:4) = npar(1:4)
    initial_fraction = 1.0E-3_NER
    tolerance = 1.0E-7_NER
    errmsg = 'fit_deltasfr_sngl: Something went wrong'
!   debug
    print *, npar_Qn
!    print *, chn%chnid
    call PrintOutFitSngl_fitting(GetNPSngl_smpl_fitdeltasngl,                         &
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
end select
contains

  subroutine init_smpl()

    regtype = REGTYPE_GAUSSIAN

    allocate(obj_smplchn::chn)

    chn%vanishing(0) = .true.
    chn%pVfunc_V0 => zero_V

    chn%pVfunc_V1 => VNLO_del_SFR_smpl_epwrap
    select case(chn%chnid)
      case (CHNID_3P2_PP)
        chn%npar(1:4) = (/ 0, 1, 2, 5/)
      case (CHNID_3P0_PP, CHNID_3P1_PP, CHNID_1P1_PP)
        chn%npar(1:4) = (/ 0, 1, 2, 4/)
      case (CHNID_1D2_PP, CHNID_3D2_PP, CHNID_3D3_PP)
        chn%npar(1:4) = (/ 0, 0, 0, 1/)
      case default
        chn%npar(1:4) = (/ 0, 0, 0, 0/)
    end select

    chn%V1_cheb_base_n = 4
    chn%V1_cheb_bnd = 0.05_NER

    chn%pVfunc_V2 => VN2LO_del_SFR_smpl_epwrap
    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%pVfunc_V3 => VN3LO_del_SFR_smpl_epwrap
    chn%V3_cheb_base_n = 5
    chn%V3_cheb_bnd = 1.0E-6_NER

  end subroutine init_smpl

  subroutine GetNPSngl_smpl_fitdeltasngl(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, phaselst)

!     use mod_obj_smplchn
!     use drv_pwphase
!     use mod_vfunc_smpl
!     implicit none

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    real(NER), intent(out)  :: phaselst(:, :)

    call create_smplchn(chn, chnid, uptoQn)
!     call SetUp_del_SFR_smpl(chn)

    call get_phase_as_value_for_klist_smpl(chn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, phaselst)

  end subroutine GetNPSngl_smpl_fitdeltasngl

  subroutine GetNPCpldSfr_smpl_fitdeltacpld(Nmesh, Mambda, chnid, uptoQn, paralst, &
    & num_paras, klist, numk, phaselst)

    real(NER), intent(in)   :: Mambda, paralst(:), klist(:)
    integer, intent(in)     :: Nmesh, uptoQn, chnid, num_paras, numk
    type(triplet_phs_cpld), intent(out)      :: phaselst(:, :)

    call create_smplchn(chn, chnid, uptoQn)
    call get_phase_as_value_for_klist_smpl(chn, Nmesh, regtype, Mambda, klist, &
      & numk, paralst, num_paras, phaselst)

  end subroutine GetNPCpldSfr_smpl_fitdeltacpld
end program fit_deltasfr_smpl