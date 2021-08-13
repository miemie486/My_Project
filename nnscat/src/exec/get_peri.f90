! Bingwei Long 05/01/2018
! Calculate phase shifts for peripheral channels where OPE is perturbative


! Usage: get_peri Nmesh chn uptoQn fxkopt
! Nmesh: number of mesh points, normally 50 ~ 120
! chn = '1p1', '3p1', '3p2', '1d2', '3d2', '3d3' ...
! uptoQn = 1, 2, or 3
! fxkopt = fxk or nofxk
! Input files:
!   * inputs_chn.in, e.g., inputs_1d2.in
!   * run_klist.in
program get_peri

  use drv_pwphase
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
      & 'get_peri: Something wrong with the arguments'
    write (standard_error_unit, '(a)')  &
      & 'Usage: get_peri N chn uptoQn opt [setn_cis]'
    stop
  end if
  select case (trim(opt))
    case ('fxk')
      fxk_is_needed = .true.
    case ('nofxk')
      fxk_is_needed = .false.
    case default
          write (standard_error_unit, '(a)')  &
            & 'get_peri: unknown option '//trim(opt)
  end select
  if (uptoQn == 3)  then
    read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt, setn_cis
    if  (ios /= 0) then
      write (standard_error_unit, '(a)')  &
        & 'get_peri: for uptoQn = 3, setn_cis is needed'
      write (standard_error_unit, '(a)')  &
        & 'Usage: get_peri N chn uptoQn opt [setn_cis]'
      stop
    end if
    if (setn_cis < 0 .or. setn_cis > 4) then
      write (standard_error_unit, '(a)') 'get_peri: setn_cis = 1 - 4'
      stop
    else
      call Select_cis2PE_vtpe1(setn_cis)
    end if
  end if

  call allocate_create_peripheral_class(ptr_chn, chnstr, uptoQn, successful)
  if (.not. successful) stop
  ! Set N2LO and N3LO long-range potentials
  ! You can replace Vfunc_VTPE0_* or Vfunc_VTPE1_* to other long-range Vfunc
  select type (ptr_chn)
    class is (obj_sngl)
      ptr_chn%pVfunc_VL2 => Vfunc_VTPE0_sngl
      ptr_chn%pVfunc_VL3 => Vfunc_VTPE1_sngl
    class is (obj_cpld)
      ptr_chn%pVfunc_VL2 => Vfunc_VTPE0_cpld
      ptr_chn%pVfunc_VL3 => Vfunc_VTPE1_cpld
  end select
  call get_output_for_Lambda_list(ptr_chn, N, regtype, OUTOPT_PHASE, fxk_is_needed)
  call ptr_chn%erase()
  deallocate(ptr_chn)

end program get_peri
