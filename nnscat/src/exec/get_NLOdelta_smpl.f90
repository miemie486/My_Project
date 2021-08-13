program get_NLOdelta_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  use nneft_phyconst,  only: Delta_paras, PC_delta_cutoff
  use mod_deltafullNLO_pwd
  use mod_deltafullN2LO_pwd, only:set_cis2PE_vdel
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  integer                     :: N, chnid, uptoQn = 0, setn_cis
  class(obj_smplchn), pointer :: ptr_chn => null()
  character(len = 3)          :: cmdstr
  logical                     :: fxk_is_needed = .false., successful
  real(NER)                   :: t0, t1
  type(Delta_paras)           ::delparas
  call cpu_time(t0)

  call  parse_cmdln_getprintout(N, cmdstr, uptoQn, fxk_is_needed, setn_cis, successful)
  if (.not. successful) stop
    !       PC_gA = 1.27_NER
!    Krebs07 ha = 1.34
!     Delparas%del = 5000.0_NER
    delparas%ha=1.34_NER
!     delparas%ha=5.34_NER
  call set_2PEparas_vdel( delparas%ha, Delparas%del)
!   call set_cisSFR_vtpe1(-0.57E-3_NER, -0.79E-3_NER, 1.33E-3_NER)
  call set_cis2PE_vdel(-0.57E-3_NER, -0.25E-3_NER, -0.79E-3_NER, 1.33E-3_NER, delparas%ha, Delparas%del)
!   call set_cisSFR_vtpe1(-1.42E-3_NER, -5.58E-3_NER, 3.50E-3_NER)

! call set_cisSFR_vtpe1(-0.57E-3_NER, -3.87E-3_NER, 2.89E-3_NER)
! call set_cis2PE_vdel(-0.57E-3_NER, 2.84E-3_NER, -3.87E-3_NER, 2.89E-3_NER, delparas%ha, Delparas%del)



!       Delparas%ha = 1.05_NER
!       call set_2PEparas_vdel( Delparas%ha, Delparas%del)
!       call set_cisSFR_vtpe1(-0.57E-3_NER, -1.87E-3_NER, 1.87E-3_NER)
!       call set_cis2PE_vdel(-0.57E-3_NER, 0.83E-3_NER, -1.87E-3_NER, 1.87E-3_NER, Delparas%ha, Delparas%del)

      call set_cisSFR_vtpe1(-0.81E-3_NER, -4.70E-3_NER, 3.40E-3_NER)
!       Delparas%ha = 0.53_NER
!       call set_2PEparas_vdel( 0.53_NER, Delparas%del)

!   if (uptoQn >= 4) then
!     call Select_cis2PE_vtpe1(setn_cis)
!     call Select_cis2PE_vdel(setn_cis)
!   end if

    ! PC_delta_cutoff      = Mambda
!     PC_delta_cutoff      = 700.0_NER

  call convert_text_to_chnid(cmdstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_vdelta_smpl: converting '//cmdstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::ptr_chn)
  call create_smplchn(ptr_chn, chnid, uptoQn)

  ptr_chn%vanishing(0) = .true.
  select case (cmdstr)
     case ('1d2', '3d2', '3d3', '1f3', '3f3', '3f4')
       call setdw(ptr_chn)
     ! case ('1d2', '3d2', '3d3','1f3','3f3', '3f4')
     !   call setUpvfunc_withsupdelta_highpwd_smpl
     case default
      write (standard_error_unit, '(a)')  &
        & 'get_NLOdelta_smpl: do not understand '//trim(cmdstr)
      stop
   end select

  call get_phase_as_text_for_MmbdkcmList_smpl(ptr_chn, N, regtype,    &
    & fxk_is_needed, preamble)
  call ptr_chn%erase()
  deallocate(ptr_chn)

  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'

contains


  subroutine preamble(funt)

    integer, intent(in) :: funt

    character(len=1024) :: cmdline

    call get_command(cmdline)
    write(funt, '(a)') '# '//trim(cmdline)
    call write_to_file_basic_phycnst(funt)

  end subroutine preamble

  subroutine setdw(chn)

    class(obj_smplchn),intent(inout)  :: chn

    chn%vanishing(1) = .false.

    chn%npar(1:5) = (/0, 0, 0, 0, 1/)
!     else
!       chn%npar(1:4) = NPAR_WITHC0_CPLD_SMPL

    chn%pVfunc_V0 => zero_V

!     chn%pVfunc_V1 =>TPE0_delonly_epwrap
!     chn%pVfunc_V1 => OPEplusNLOdelta_epwrap
!     chn%pVfunc_V1 => OPEplusNLOdeltaonly_epwrap
!     chn%pVfunc_V1 => OPEplusTPE_epwrap
!     chn%pVfunc_V1 => OPE_epwrap
!     chn%pVfunc_V1 => OPEplusTPESFR_epwrap
    chn%pVfunc_V1 => OPEplusAllTPESFR_epwrap
!       chn%pVfunc_V1 => OPEplusTPE1deltaSFR_epwrap
!       chn%pVfunc_V1 => OPEplusallTPEdeltaSFR_epwrap

!         chn%pVfunc_V1 => TPESFR1_epwrap
!       chn%pVfunc_V1 => TPE0_delonly_epwrap
!       chn%pVfunc_V1 => TPE1_delonly_epwrap
    chn%V1_cheb_base_n = 4
    chn%V1_cheb_bnd = 1.0E-3_NER

  end subroutine setdw
end program get_NLOdelta_smpl