! pray 2019/12/16

program get_NP_3s1_smpl

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

  call  parse_cmdln_getprintout(N, cmdstr, uptoQn, fxk_is_needed, setn_cis, successful)
  if (.not. successful) stop
  if (setn_cis > 0) call Select_cis2PE_vtpe1(setn_cis)
  if (setn_cis > 4) call Select_cis2PE_vdel(setn_cis)

  chnstr = cmdstr

  call convert_text_to_chnid(chnstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_NP_3s1_smpl: converting '//chnstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::chn)
  call create_smplchn(chn, chnid, uptoQn)

!   chn%potpara(1) = -4.5126679738286763E-003_NER   ! cuoff = 400 MeV
!   chn%potpara(1) = -5.4945264707960580E-003_NER   ! 350
  call SetUp_3s1_offdigper_smpl(chn)

  call get_phase_as_text_for_MmbdkcmList_smpl(chn, N, regtype,    &
    & fxk_is_needed, preamble)
    print *, chn%potpara(1), chn%potpara(2)
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


  subroutine SetUp_3s1_offdigper_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%npar(1:5) = (/1, 1, 1, 1, 1/)
    chn%pVfunc_V0 => VLO_withc011_smpl_epwrap
    chn%pVfunc_V1 => VNLO_withc0offdig_smpl_epwrap
    chn%V1_cheb_base_n = 5
    chn%V1_cheb_bnd = 1.0E-2_NER
    chn%pVfunc_V2 => zero_V
    chn%pVfunc_V3 => zero_V
    chn%pVfunc_V4 => zero_V

  end subroutine SetUp_3s1_offdigper_smpl

end program get_NP_3s1_smpl
