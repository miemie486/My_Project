! Bingwei Long 07/17/2018
! Calculate phase shifts for uncoupled channels with L >= 2
!
! Usage: get_vnfperi chn uptoQn fxkopt
! chn = '1d2', '3d2', '1f3', '3f3', ...
! uptoQn = 2
! fxkopt = 'fxk' or 'nofxk'
! Input files:
!   * inputs_chn.in, e.g., inputs_1d2.in
!   * run_klist.in
! E.g.
! $ get_vnfperi 1d2 2 nofxk
program get_vnfperi

  use mod_obj_pope_vnfsngl
  use drv_pwphase
  implicit none

  integer, parameter          :: N = 100, regtype = REGTYPE_GAUSSIAN
  integer                     :: ii, ios, chnid, uptoQn = 0, cmdcnt
  class(obj_channel), pointer :: ptr_chn => null()
  character(len = 3)          :: chnstr, my_chnstr
  character(len = 10)         :: opt
  character(len = 1024)       :: cmdbuff, tmp
  logical                     :: fxk_is_needed = .false., succeeded
  type(lsj_symbol)    :: lsj

  cmdbuff = ''
  cmdcnt = command_argument_count()
  do ii = 1, cmdcnt
    call get_command_argument(ii, tmp)
    cmdbuff = trim(cmdbuff)//' '//trim(tmp)
  end do
  read (cmdbuff, *, IOSTAT=ios) chnstr, uptoQn, opt
  if  (ios /= 0) then
    write (standard_error_unit, '(a)')                                      &
      & 'get_vnfperi: Something wrong with the arguments'
    write (standard_error_unit, '(a)')                                      &
      & 'Usage: get_vnfperi chn uptoQn opt'
    stop
  end if
  select case (trim(opt))
    case ('fxk')
      fxk_is_needed = .true.
    case ('nofxk')
      fxk_is_needed = .false.
    case default
      write (standard_error_unit, '(a)')                                    &
      & 'get_vnfperi: unknown option '//trim(opt)
  end select

  my_chnstr = chnstr(1:3)
  call convert_text_to_lsj(my_chnstr, lsj, succeeded)
  if (.not. succeeded) then
    write (standard_error_unit, '(a, a3, a)')                               &
    & 'allocate_create_peripheral_class: Conversion of ', my_chnstr,        &
    & ' to an channel IDfailed. Nothing will be done.'
    return
  end if
  call convert_lsj_to_chnid(lsj, chnid)
  if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
    allocate(obj_pope_vnfsngl::ptr_chn)
  else
    write (standard_error_unit, '(a)')                                      &
      & 'get_vnfperi: '//my_chnstr//' is not an uncoupled channel.'
  end if
  call create_channel(ptr_chn, chnid, uptoQn)

  call get_output_for_Lambda_list(ptr_chn, N, regtype, OUTOPT_PHASE, fxk_is_needed)
  call ptr_chn%erase()
  deallocate(ptr_chn)

end program get_vnfperi
