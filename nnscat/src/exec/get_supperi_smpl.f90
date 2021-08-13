program get_supperi_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  integer                     :: N, chnid, uptoQn = 0, setn_cis
  class(obj_smplchn), allocatable :: chn
  character(len = 3)          :: cmdstr
  logical                     :: fxk_is_needed = .false., successful

  call parse_cmdln_getprintout(N, cmdstr, uptoQn, fxk_is_needed, setn_cis, successful)
  if (.not. successful) stop

  if (uptoQn >= 4) then
    call Select_cis2PE_vtpe1(setn_cis)
    call Select_cis2PE_vdel(setn_cis)
  end if

  call convert_text_to_chnid(cmdstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_vdelta_smpl: converting '//cmdstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::chn)
  call create_smplchn(chn, chnid, uptoQn)
  call SetUp_pope_smpl(chn)

  call get_phase_as_text_for_MmbdkcmList_smpl(chn, N, regtype, fxk_is_needed,  &
    & preamble)
  call chn%erase()
  deallocate(chn)

contains

  subroutine preamble(funt)

    integer, intent(in) :: funt

    character(len=1024) :: cmdline

    call get_command(cmdline)
    write(funt, '(a)') '# '//trim(cmdline)
    call write_to_file_basic_phycnst(funt)

  end subroutine preamble

end program get_supperi_smpl
