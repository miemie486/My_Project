! pray 2020/08/02
! using the pot form chengdu to investigate the values of parameters
! getting phase shifts with chengdu
! notice :: there is no Mambda initial, mambda in sub: get_phase_as_text_for_chengdu_smpl

program get_phase_withCD_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN, N_mambda = 40
  integer                     :: N, chnid, uptoQn = 0, setn_cis, kk
  class(obj_smplchn), allocatable :: chn
  character(len = 3)          :: chnstr
  character(len = 10)         :: cmdstr
  character(len = 20)         :: pot(N_mambda)
  logical                     :: fxk_is_needed = .false., successful
  real(NER)                   :: t0, t1, Mambda(N_mambda)

  call cpu_time(t0)

  call parse_cmdln_getprintout(N, cmdstr, uptoQn, fxk_is_needed, setn_cis, successful)
  if (.not. successful) stop
  if (setn_cis > 0) call Select_cis2PE_vtpe1(setn_cis)
  if (setn_cis > 4) call Select_cis2PE_vdel(setn_cis)

    chnstr = trim(cmdstr)

  call convert_text_to_chnid(chnstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_phase_withCD_smpl: converting '//chnstr//' wasn''t successful'
    stop
  end if

  pot(1: N_mambda) = (/'chengdu_MMWLY_350 ','chengdu_MMWLY_400 ','chengdu_MMWLY_600 ', &
      &'chengdu_MMWLY_800 ','chengdu_MMWLY_1600','chengdu_MMWLY_3200','chengdu_MMWLY_1200','chengdu_MMWLY_2000','chengdu_MMWLY_2400','chengdu_MMWLY_2800', &
      &'chengdu_MMWLY_3600','chengdu_MMWLY_4000','chengdu_MMWLY_4400','chengdu_MMWLY_4800','chengdu_MMWLY_500 ','chengdu_MMWLY_1800','chengdu_MMWLY_2200', &
      &'chengdu_MMWLY_900 ','chengdu_MMWLY_1300','chengdu_MMWLY_1400','chengdu_MMWLY_6400','chengdu_VLbar_1200','chengdu_VLbar_1300','chengdu_VLbar_1400', &
      &'chengdu_VLbar_1600','chengdu_WOTPE_1000','chengdu_WOTPE_1200','chengdu_WOTPE_1300','chengdu_WOTPE_1400','chengdu_WOTPE_1600','chengdu_WOTPE_2000', &
      &'chengdu_WOTPE_2200','chengdu_WOTPE_2400','chengdu_WOTPE_800 ','chengdu_WOTPE_2800','chengdu_WOTPE_3200','chengdu_WOTPE_4000','chengdu_WOTPE_4800', &
      &'chengdu_WOTPE_400 ','chengdu_WOTPE_600 ' /)
  Mambda(1: N_mambda) = (/350.0_NER, 400.0_NER, 600.0_NER, 800.0_NER, 1600.0_NER, 3200.0_NER, 1200.0_NER, 2000.0_NER, 2400.0_NER, 2800.0_NER, 3600.0_NER, 4000.0_NER, &
      & 4400.0_NER, 4800.0_NER, 500.0_NER, 1800.0_NER, 2200.0_NER, 900.0_NER, 1300.0_NER, 1400.0_NER, 6400.0_NER, 1200.0_NER, 1300.0_NER, 1400.0_NER, 1600.0_NER,     &
      & 1000.0_NER, 1200.0_NER, 1300.0_NER, 1400.0_NER, 1600.0_NER, 2000.0_NER, 2200.0_NER, 2400.0_NER, 800.0_NER, 2800.0_NER, 3200.0_NER, 4000.0_NER, 4800.0_NER,    &
      & 4000.0_NER, 600.0_NER/)
  do kk = 1, N_mambda
    allocate(obj_smplchn::chn)
    call create_smplchn(chn, chnid, uptoQn)
    call SetUpVfunc_mmwly_chengdu_smpl(chn, pot(kk))
    call get_phase_as_text_for_chengdu_smpl(chn, N, regtype,    &
      & fxk_is_needed, preamble, Mambda(kk), pot(kk))
!     call get_phase_as_text_for_MmbdkcmList_smpl(chn, N, regtype,    &
!       & fxk_is_needed, preamble)
    call chn%erase()
    deallocate(chn)
  end do
  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'

contains

  subroutine preamble(funt)

    integer, intent(in) :: funt

    integer             :: ii
    character(len=1024) :: cmdline

    call get_command(cmdline)
    write(funt, '(a)') '# '//trim(cmdline)
    write(funt, *) '# '//trim(pot(kk))

  end subroutine preamble

end program get_phase_withCD_smpl
