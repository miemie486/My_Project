! for deltaless SFR potential 9/11/2019 pray

program get_supsfr_chn


  use drv_pwphase
  use mod_vfunc
  use VTPE1_SFR_pwd, only: Select_cisSFR_vtpe1
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  integer                     :: N, ii, ios, chnid, uptoQn = 0, cmdcnt, setn_cis
  class(obj_channel), pointer :: ptr_chn => null()
  character(len = 3)          :: chnstr
  character(len = 10)         :: opt
  character(len = 1024)       :: cmdbuff, tmp
  logical                     :: fxk_is_needed = .false., successful
  real(NER)                   :: t0, t1

  call cpu_time(t0)

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
    read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt, setn_cis
    if  (ios /= 0) then
      write (standard_error_unit, '(a)')  &
        & 'get_supperi: for uptoQn = 4, setn_cis is needed'
      write (standard_error_unit, '(a)')  &
        & 'Usage: get_supperi N chn uptoQn opt [setn_cis]'
      stop
    end if
      call Select_cisSFR_vtpe1(setn_cis)
  end if

  call allocate_create_supperi_class(ptr_chn, chnstr, uptoQn, successful)
  if (.not. successful) stop
  ! Set N2LO and N3LO long-range potentials
  ! You can replace Vfunc_VTPE0_* or Vfunc_VTPE1_* to other long-range Vfunc
  select type (ptr_chn)
   class is (obj_SUPpope_sngl)
      ptr_chn%pVfunc_VL2 => Vfunc_sfrTPE0_sngl
      ptr_chn%pVfunc_VL3 => Vfunc_sfrTPE1_sngl
   class is (obj_SUPpope_cpld)
      ptr_chn%pVfunc_VL2 => Vfunc_sfrTPE0_cpld
      ptr_chn%pVfunc_VL3 => Vfunc_sfrTPE1_cpld
  end select
  call get_output_for_Lambda_list(ptr_chn, N, regtype, OUTOPT_PHASE, fxk_is_needed)
  call ptr_chn%erase()
  deallocate(ptr_chn)

  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'

end program get_supsfr_chn