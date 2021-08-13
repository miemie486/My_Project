! Bingwei Long 10/21/2018
! Calculate nonperturbative channels
! Using obj_smplchn
! Using deltaless TPEs
program get_3p0_smpl

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
!   if (setn_cis > 4) call Select_cis2PE_vdel(setn_cis)

  select case (trim(cmdstr))
    case ('per3p0','hper3p0')
      chnstr = '3p0'
    case ('sup3p0','3p0','m3p0','wot3p0')
      chnstr = '3p0'
  end select

  call convert_text_to_chnid(chnstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_3p0_smpl: converting '//chnstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::chn)
  call create_smplchn(chn, chnid, uptoQn)

  select case (trim(cmdstr))
    case ('sup3p0')
      call setUP_withc0_sup3p0_smpl(chn)

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%V3_cheb_base_n = 8
    chn%V3_cheb_bnd = 1.0E-6_NER

    chn%V4_cheb_base_n = 8
    chn%V4_cheb_bnd = 1.0E-7_NER
    case ('per3p0')
      call setUP_withc0_3p0_smpl(chn)

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%V3_cheb_base_n = 8
    chn%V3_cheb_bnd = 1.0E-5_NER

    case ('hper3p0')
      call setUP_withc0_high3p0_smpl(chn)

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%V3_cheb_base_n = 8
    chn%V3_cheb_bnd = 1.0E-9_NER

    case ('wot3p0')
      call setUP_withc0_highwot3p0_smpl(chn)

    chn%V2_cheb_base_n = 5
    chn%V2_cheb_bnd = 1.0E-6_NER

    chn%V3_cheb_base_n = 8
    chn%V3_cheb_bnd = 1.0E-9_NER

    case ('3p0')
      call setUP_withc0_3p0NLO_smpl(chn)
      chn%V1_cheb_base_n = 4
      chn%V1_cheb_bnd = 0.05_NER

    case('m3p0')
      call setUP_withc0_m3p0_smpl(chn)
    case default
      write (standard_error_unit, '(a)')  &
        & 'get_3p0_smpl: do not understand '//trim(cmdstr)
      stop
  end select

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

end program get_3p0_smpl
