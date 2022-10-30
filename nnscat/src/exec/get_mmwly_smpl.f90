! pray 25/7/2020
! using mmwly power countting

program get_mmwly_smpl

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

!   select case (trim(cmdstr))
!     case ('1s0')
!       chnstr = '1s0'
!     case ('ros')
!       chnstr = '3s1'
!     case ('rop')
!       chnstr = '3p0'
!     case ('3s1')
!       chnstr = '3s1'
!     case ('3p0')
!       chnstr = '3p0'
!     case ('3p0')
!       chnstr = '3p0'
!     case default
      chnstr = cmdstr
!   end select

  call convert_text_to_chnid(chnstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_mmwly_smpl: converting '//chnstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::chn)
  call create_smplchn(chn, chnid, uptoQn)
  call SetUpVfunc_mmwly_smpl(chn)

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

end program get_mmwly_smpl
