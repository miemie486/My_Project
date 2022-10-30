! Bingwei Long 10/21/2018
! Calculate nonperturbative channels
! Using obj_smplchn
! Using deltaless TPEs
program get_NP_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  integer                     :: N, chnid, uptoQn = 0, setn_cis
  class(obj_smplchn), allocatable :: chn
  character(len = 3)          :: chnstr
  character(len = 10)         :: cmdstr
  logical                     :: fxk_is_needed = .false., successful
  real(NER)                   :: t0, t1

  call cpu_time(t0)

  call parse_cmdln_getprintout(N, cmdstr, uptoQn, fxk_is_needed, setn_cis, successful)

  if (.not. successful) stop
  if (setn_cis > 0) call Select_cis2PE_vtpe1(setn_cis)
  if (setn_cis > 4) call Select_cis2PE_vdel(setn_cis)

  select case (trim(cmdstr))
    case ('spr')
      chnstr = '1s0'
    case ('ros')
      chnstr = '3s1'
    case ('rop')
      chnstr = '3p0'
    case ('del3s1')
      chnstr = '3s1'
    case ('del3p0')
      chnstr = '3p0'
    case ('per3p0')
      chnstr = '3p0'
    case default
      chnstr = cmdstr
  end select

  call convert_text_to_chnid(chnstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_NP_smpl: converting '//chnstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::chn)
  call create_smplchn(chn, chnid, uptoQn)

  select case (trim(cmdstr))
    case ('1s0')
      call SetUpVfunc_mmwly_smpl(chn)
    case ('3s1', '3p0')
      call SetUpVfunc_withc0_smpl(chn)
    case ('del3s1', 'del3p0')
      call SetUp_withc0_del_smpl(chn)
    case ('per3p0')
      call SetUpVfunc_pert3P0(chn)
    case ('spr')

      ! V1_BASEN_sprbl_mvfs = 2
      ! V1_BND_sprbl_mvfs = 1.0E-3_NER
      ! V2_BASEN_sprbl_mvfs = 2
      ! V2_BND_sprbl_mvfs = 1.0E-4_NER
      ! V3_BASEN_sprbl_mvfs = 2
      ! V3_BND_sprbl_mvfs = 1.0E-4_NER

      call SetUpVfunc_sprbl_smpl(chn)

    case ('ros', 'rop')
      call SetUpVfunc_withc0_relope(chn)

    case default
      write (standard_error_unit, '(a)')  &
        & 'get_NP_smpl: do not understand '//trim(cmdstr)
      stop
  end select

  call get_phase_as_text_for_MmbdkcmList_smpl(chn, N, regtype,    &
    & fxk_is_needed, preamble)
  call chn%erase()
  deallocate(chn)

  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'

contains

  subroutine preamble(funt)

    integer, intent(in) :: funt

    integer             :: ii
    character(len=1024) :: cmdline

    call get_command(cmdline)
    write(funt, '(a)') '# '//trim(cmdline)
    do ii = 1, chn%npar(chn%uptoQn+1)
      write(funt, '(a, i2, a, e19.10e3)') '# Channel para (', ii, ') = ', chn%potpara(ii)
    end do

  end subroutine preamble

  subroutine SetUpVfunc_withc0_relope(chn)

    class(obj_smplchn), intent(inout) :: chn

    ! NLO pot vanishes
    chn%vanishing(1) = .true.
    if (type_of_channel_PP(chn%chnid) == CHTYPE_SNGL) then
      chn%npar(1:4) = NPAR_WITHC0_SNGL_SMPL
    else
      chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL
    end if

    chn%pVfunc_V0 => VLO_withc0_epwrap

    chn%pVfunc_V2 => VN2LO_withc0_RelOPE_epwrap
    chn%V2_cheb_base_n = 4
    chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V3 => VN3LO_withc0_epwrap
    chn%V3_cheb_base_n = 4
    chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine SetUpVfunc_withc0_relope

  subroutine SetUpVfunc_pert3P0(chn)

    class(obj_smplchn), intent(inout) :: chn

    ! LO vanishes
    chn%vanishing(0) = .true.
    chn%vanishing(2) = .true.
    chn%npar(1:4) = (/0, 1, 1, 1/)

    chn%pVfunc_V0 => zero_V
    chn%pVfunc_V1 => VLO_withc0_epwrap

    ! chn%pVfunc_V2 => VN2LO_withc0_RelOPE_epwrap
    ! chn%V2_cheb_base_n = 4
    ! chn%V2_cheb_bnd = 1.0E-3_NER

    ! chn%pVfunc_V3 => VN3LO_withc0_epwrap
    ! chn%V3_cheb_base_n = 4
    ! chn%V3_cheb_bnd = 1.0E-3_NER

  end subroutine SetUpVfunc_pert3P0

end program get_NP_smpl
