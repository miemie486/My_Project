! Bingwei Long 07/24/2018
! Bingwei Long 08/02/2013
! Calculate channels where OPE is nonperturbative

program get_low

  use drv_pwphase
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  integer                     :: N, ii, ios, chnid, uptoQn = 0, cmdcnt, setn_cis
  class(obj_channel), pointer :: ptr_chn => null()
  character(len = 1024)       :: cmdbuff, tmp
  character(len = 3)          :: chnstr
  character(len = 10)         :: opt
  logical                     :: fxk_is_needed = .false., successful

  cmdbuff = ''
  cmdcnt = command_argument_count()
  do ii = 1, cmdcnt
    call get_command_argument(ii, tmp)
    cmdbuff = trim(cmdbuff)//' '//trim(tmp)
  end do
  read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt
  if  (ios /= 0) then
    write (standard_error_unit, '(a)')  &
      & 'get_low: Something wrong with the arguments'
    write (standard_error_unit, '(a)')  &
      & 'Usage: get_low N chn uptoQn fxkopt [setn_cis]'
    stop
  end if
  select case (trim(opt))
    case ('fxk')
      fxk_is_needed = .true.
    case ('nofxk')
      fxk_is_needed = .false.
    case default
      write (standard_error_unit, '(a)')  &
        & 'get_low: unknown option '//trim(opt)
      stop
  end select

  if (uptoQn == 3) then
    read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt, setn_cis
    if  (ios /= 0) then
      write (standard_error_unit, '(a)')  &
        & 'get_low: for uptoQn = 3, setn_cis is needed'
      write (standard_error_unit, '(a)')  &
        & 'Usage: get_low N chn uptoQn opt [setn_cis]'
      stop
    end if
    if (setn_cis < 0 .or. setn_cis > 4) then
      write (standard_error_unit, '(a)') 'get_low: setn_cis = 1 - 4'
      stop
    else
      call Select_cis2PE_vtpe1(setn_cis)
    end if
  end if

  call AllocateNonPertCHN(ptr_chn, chnstr, uptoQn, successful)
  if (.not. successful) stop
  call get_output_for_Lambda_list(ptr_chn, N, regtype, OUTOPT_PHASE,           &
    & fxk_is_needed)
  call ptr_chn%erase()
  deallocate(ptr_chn)

end program get_low
