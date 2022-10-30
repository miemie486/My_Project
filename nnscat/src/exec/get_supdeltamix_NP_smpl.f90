!  Pary 12/20/2018
!  Calculate channels where OPE is nonperturbative using obj_smplchn

program get_supdeltamix_NP_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  use nneft_phyconst,  only: PC_hA
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  integer                     :: N, chnid, uptoQn = 0, setn_cis
  class(obj_smplchn), pointer :: ptr_chn => null()
  character(len = 3)          :: cmdstr
  logical                     :: fxk_is_needed = .false., successful

  call  parse_cmdln_getprintout(N, cmdstr, uptoQn, fxk_is_needed, setn_cis, successful)
  if (.not. successful) stop

  if (uptoQn >= 3) then
      call Select_cis2PE_vtpe1(setn_cis)
      call Select_cis2PE_vdel(setn_cis)
  end if

  call convert_text_to_chnid(cmdstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_vdelta_smpl: converting '//cmdstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::ptr_chn)
  call create_smplchn(ptr_chn, chnid, uptoQn)

  select case (cmdstr)
!      case ('1s0')
!        call setUpvfunc_1s0_withdelta_smpl(ptr_chn)
     case ('3s1', '3p0')
       call SetUp_withc0_del_smpl(ptr_chn)
     ! case ('1d2', '3d2', '3d3','1f3','3f3', '3f4')
     !   call setUpvfunc_withsupdelta_highpwd_smpl
     case default
      write (standard_error_unit, '(a)')  &
        & 'get_vdelta_smpl: do not understand '//trim(cmdstr)
      stop
   end select

  call get_phase_as_text_for_MmbdkcmList_smpl(ptr_chn, N, regtype,    &
    & fxk_is_needed, preamble)
  call ptr_chn%erase()
  deallocate(ptr_chn)

contains

  subroutine preamble(funt)

    integer, intent(in) :: funt

    character(len=1024) :: cmdline

    call get_command(cmdline)
    write(funt, '(a)') '# '//trim(cmdline)
    call write_to_file_basic_phycnst(funt)

  end subroutine preamble

end program get_supdeltamix_NP_smpl
