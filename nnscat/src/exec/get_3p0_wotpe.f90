! get 3p0 without TPEP
program get_3p0_wotpe

  use drv_pwphase
  use mod_vfunc

  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  integer                     :: N, ii, ios, chnid, uptoQn = 0, cmdcnt, setn_cis
  class(obj_channel), pointer :: ptr_chn => null()
  character(len = 3)          :: chnstr
  character(len = 10)         :: opt
  character(len = 1024)       :: cmdbuff, tmp
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
      & 'get_supperi: Something wrong with the arguments'
    write (standard_error_unit, '(a)')  &
      & 'Usage: get_supperi N chn uptoQn opt [setn_cis]'
    stop
  end if

  select case (trim(opt))
    case ('fxk')
      fxk_is_needed = .true.
    case ('nofxk')
      fxk_is_needed = .false.
    case default
      write (standard_error_unit, '(a)')  &
        & 'get_supperi: unknown option '//trim(opt)
  end select

  if (uptoQn == 4)  then
    read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt
    if  (ios /= 0) then
      write (standard_error_unit, '(a)')  &
        & 'get_supperi: for uptoQn = 4, setn_cis is needed'
      write (standard_error_unit, '(a)')  &
        & 'Usage: get_supperi N chn uptoQn opt [setn_cis]'
      stop
    end if
    ! if (setn_cis < 0 .or. setn_cis > 4) then
    !   write (standard_error_unit, '(a)') 'get_supperi: setn_cis = 1 - 4'
    !   stop
    ! else
!       call Select_cis2PE_vtpe1(setn_cis)
!       call Select_cis2PE_vdel(setn_cis)
    ! end if
  end if

  call allocate_create_supperi_class(ptr_chn, chnstr, uptoQn, successful)
  if (.not. successful) stop
  ! Set N2LO and N3LO long-range potentials
  ! You can replace Vfunc_VTPE0_* or Vfunc_VTPE1_* to other long-range Vfunc
     select type (ptr_chn)
   class is (obj_SUPpope_sngl)
      ptr_chn%pVfunc_VL2 => V0sngl
      ptr_chn%pVfunc_VL3 => V0sngl


  end select
  call get_output_for_Lambda_list(ptr_chn, N, regtype, OUTOPT_PHASE, fxk_is_needed)
  call ptr_chn%erase()
  deallocate(ptr_chn)

end program get_3p0_wotpe



