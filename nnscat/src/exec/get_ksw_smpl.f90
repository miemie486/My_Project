! pray 05/29/2020
! Calculate ksw in Lambda = 300-400MeV in 3s1 and 1s0
! Using obj_smplchn
! Using deltaless TPEs

program get_ksw_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  use testwrap
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
!   debug
 print *,  setn_cis
!  print *, uptoQn, setn_cis
  if (.not. successful) stop
   if (setn_cis > 0) call Select_cis2PE_vtpe1(setn_cis)
!   if (setn_cis > 4) call Select_cis2PE_vdel(setn_cis)

  chnstr = cmdstr

  call convert_text_to_chnid(chnstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_ksw_smpl: converting '//chnstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::chn)
    call create_smplchn(chn, chnid, uptoQn)

    call setUP_ksw_smpl(chn)

! #######################
 ! set cheby paras
! #######################

!   select case(chnstr)

!     case ('3s1')
!     chn%V2_cheb_base_n = 5
!     chn%V2_cheb_bnd = 1.0E-6_NER

!     chn%V3_cheb_base_n = 8
!     chn%V3_cheb_bnd = 1.0E-9_NER

!     case ('1s0')
!     chn%V2_cheb_base_n = 5
!     chn%V2_cheb_bnd = 1.0E-6_NER

!     chn%V3_cheb_base_n = 8
!     chn%V3_cheb_bnd = 1.0E-9_NER

!     case default
!       write (standard_error_unit, '(a)')  &
!         & 'get_ksw_smpl: do not understand '//trim(cmdstr)
!       stop
!   end select

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

  subroutine SetUp_ksw_smpl(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(0) = .false.

    select case(chn%chnid)
      case (CHNID_3s1_PP)
        chn%npar(1:5) = (/1, 3, 4, 5, 6/)
      case (CHNID_1s0_PP)
        chn%npar(1:5) = (/1, 3, 4, 5, 6/)
    end select

    call get_cheby_paras(chn)
!   debug
!     print *, chn%V1_cheb_base_n, chn%V1_cheb_bnd
    chn%pVfunc_V0 => VLO_ksw_smpl_epwrap

    chn%pVfunc_V1 => VNLO_ksw_smpl_epwrap
!     chn%V1_cheb_base_n = 2
!     chn%V1_cheb_bnd = 0.005_NER

    chn%pVfunc_V2 => VN2LO_ksw_smpl_epwrap
!     chn%V2_cheb_base_n = 4
!     chn%V2_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V3 => VN3LO_ksw_smpl_epwrap
!     chn%V3_cheb_base_n = 4
!     chn%V3_cheb_bnd = 1.0E-3_NER

    chn%pVfunc_V4 => VN4LO_ksw_smpl_epwrap
!     chn%V4_cheb_base_n = 4
!     chn%V4_cheb_bnd = 1.0E-3_NER

  end subroutine SetUp_ksw_smpl



end program get_ksw_smpl