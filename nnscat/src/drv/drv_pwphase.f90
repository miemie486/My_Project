! Created by Bingwei Long in 07/28/2013
! Driver and utility subroutines that compute the phase shifts for a given channel

! July 28, 2015:
! - Removed the 3-character limit on and to rename to command_str the input parameter chnstr of allocate_create_channel_class

! Lambda, Mambda, kcm, Tlab, and input parameters in inputs_chn.in are in unit
! of MeV or MeV^-1

module drv_pwphase

  ! use mod_obj_withc0_residual_sngl
  ! use mod_obj_1s0_mpistar_sngl
  use mod_obj_1s0SYLV
  use mod_obj_1s0LY
  ! use mod_obj_tlrLY
  use mod_obj_lbwdib
  use mod_obj_sprbl1s0
  use mod_obj_cpld
  use mod_obj_withc0_sngl
  use mod_obj_pope_sngl
  use mod_obj_suppope_sngl
  use mod_obj_pope_cpld
  use mod_obj_suppope_cpld

  use mod_obj_smplchn

  use vtpe1_pwd,      only : get_cis2PE_vtpe1
  use mod_deltafullN2LO_pwd,   only : get_deltaparas_N2LO
  use eft_potspwd
  use util_io
  use util_datetime
  use mod_vfunc_smpl
  implicit none

  real(NER)   :: Mambda_1s0
  integer, parameter  :: NMAX_PHYCNST = 100
!   set values of numks_max, nLmbds_max and num_cnttrms_max
!   numks_max : the max number of k in klist could be calculated in the same time
!   nLmbds_max : the max number of Mambda could be calculated in the same time
!   num_cnttrms_max : the max number of counterterm's parameters
  integer, parameter  :: numks_max = 500, nLmbds_max = 50,   &
    & num_cnttrms_max = 50

  private     :: NMAX_PHYCNST

  ! prefix_callback is a callback function to be used to print out some preamble
  abstract interface
    subroutine prefix_callback(unt)
      integer, intent(in) :: unt
    end subroutine prefix_callback
  end interface

contains

  ! Return a list of strings that show the current values of physical constants,
  ! excluding counterterms, used by nnscat
  subroutine printout_basic_phycnst(strlst, nitem)

    character(len = *), intent(inout) :: strlst(:)
    integer, intent(out)              :: nitem

    integer             :: ii
    real(NER)           :: valuelst(1:NMAX_PHYCNST)
    type(strct_cis)     :: tmpcis
    type(Delta_paras)   :: tmpdelpara
    character(len = 64) :: namelst(1:NMAX_PHYCNST)
    character(len = 20) :: formatstr(1:NMAX_PHYCNST)

    call get_cis2PE_vtpe1(tmpcis)
    call get_deltaparas_N2LO(tmpdelpara)

    nitem = 15
    namelst(1:nitem) = (/ &
      & 'GAUSSIAN POWER                ',   &
      & 'gA                            ',   &
      & 'Average mN                    ',   &
      & 'fpi                           ',   &
      & 'Average mpi                   ',   &
      & 'hbarc                         ',   &
      & 'c1_2PE                        ',   &
      & 'c3_2PE                        ',   &
      & 'c4_2PE                        ',   &
      & 'c1_del                        ',   &
      & 'c2_del                        ',   &
      & 'c3_del                        ',   &
      & 'c4_del                        ',   &
      & 'hA_del                        ',   &
      & 'delta_del                     '    /)

    formatstr(1:nitem) = (/ &
      &  'f7.1               ',              &
      & ('f15.10             ', ii = 1, nitem-1)  /)

    valuelst(1:nitem)  = (/ &
      & real(GAUSSIAN_POWER, NER),          &
      & PC_gA,                              &
      & PC_mN,                              &
      & PC_fpi,                             &
      & PC_mpi,                             &
      & PC_hbarc,                           &
      & tmpcis%c1,                          &
      & tmpcis%c3,                          &
      & tmpcis%c4,                          &
      & tmpdelpara%c1,                      &
      & tmpdelpara%c2,                      &
      & tmpdelpara%c3,                      &
      & tmpdelpara%c4,                      &
      & tmpdelpara%ha,                      &
      & tmpdelpara%del      /)

    do ii = 1, nitem
      write (strlst(ii), '(a24, a, '//trim(formatstr(ii))//')')      &
        & '# '//namelst(ii), ' = ', valuelst(ii)
    end do

  end subroutine printout_basic_phycnst

  subroutine write_to_file_basic_phycnst(funt)

    integer, intent(in)   :: funt

    integer :: nPC, ii
    character(len = 80)   :: PClst(1:NMAX_PHYCNST)

    call printout_basic_phycnst(PClst, nPC)
    do ii = 1, nPC
      write (funt, '(a)') trim(PClst(ii))
    end do

  end subroutine write_to_file_basic_phycnst

  ! Allocate ptr_chn with the proper channel class according to chnstr
  ! When chnstr == 'dib', obj_dibaryon will be invoked.
  ! 'SAS', singular attractive S-wave
  ! 'RAS', regular attractive S-wave
  subroutine allocate_create_channel_class(ptr_chn, command_str, uptoQn, is_successful)

    class(obj_channel), pointer, intent(inout)  :: ptr_chn
    character(len = *), intent(in)              :: command_str
    integer, intent(in)                         :: uptoQn
    logical, intent(out)                        :: is_successful

    integer         :: chnid
    logical         :: succeeded, is_dibaryon
    type(lsj_symbol):: lsj
    character(len = 3) :: my_chnstr

    is_successful = .false.
    if (associated(ptr_chn)) then
      write (standard_error_unit, '(a)')  &
      & 'allocate_create_channel_class: Input pointer must be null. Nothing will be done.'
        return
    end if
    is_dibaryon = .false.
    select case (command_str)
      case ('dib')
        is_dibaryon = .true.
        my_chnstr = '1s0'
      case ('SAS')
        my_chnstr = '3s1'
      case ('RAS')
        my_chnstr = '1s0'
      case ('uni')
        my_chnstr = '1s0'
      case default
        my_chnstr = command_str(1:3)
    end select
    call convert_text_to_lsj(my_chnstr, lsj, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a, a3, a)')  &
      & 'allocate_create_channel_class: Conversion of ', my_chnstr, ' to an channel ID failed. Nothing will be done.'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case (command_str)
      ! case ('uni')
      !   allocate(obj_1s0_mpistar_sngl::ptr_chn)
      ! case ('SAS', 'RAS')
      !   allocate(obj_withc0_residual_sngl::ptr_chn)
      case default
        if (.not. is_dibaryon) then
          select case (chnid)
            case (CHNID_3P0_PP)
              allocate(obj_withc0_sngl::ptr_chn)
            case (CHNID_1P1_PP, CHNID_3P1_PP)
              allocate(obj_sngl::ptr_chn)
            case (CHNID_1S0_PP)
              allocate(obj_1s0LY::ptr_chn)
            case (CHNID_3S1_PP)
              allocate(obj_cpld::ptr_chn)
            case default
              write (standard_error_unit, '(a, a3)')  &
              & 'allocate_create_channel_class: Don''t know how to handle ', command_str
              return
          end select
        else
          allocate(obj_lbwdib::ptr_chn)
        end if
    end select

    is_successful = .true.
    call create_channel(ptr_chn, chnid, uptoQn)

  end subroutine allocate_create_channel_class

  ! Allocate and create channel for nonperturbative partial waves
  ! '1s0' : Long and Yang's 1s0
  ! 'lbw' : Long's dibaryon 1s0
  ! 'spr' : Long's sprbl1s0
  ! 'syl' : Sanchez, Yang, Long and van Kolck's 1s0
  ! 'tly' : obj_tlrLY
  subroutine AllocateNonPertCHN(ptrchn, command_str, uptoQn, is_successful)

    class(obj_channel), pointer, intent(inout)  :: ptrchn
    character(len = *), intent(in)              :: command_str
    integer, intent(in)                         :: uptoQn
    logical, intent(out)                        :: is_successful

    integer         :: chnid
    logical         :: succeeded
    type(lsj_symbol):: lsj
    character(len = 3) :: my_chnstr

    is_successful = .false.
    if (associated(ptrchn)) then
      write (standard_error_unit, '(a)')  &
        & 'AllocateNonPertCHN: Input pointer must be null. Nothing will be done.'
      return
    end if

    select case (command_str)
      case ('spr')
        my_chnstr = '1s0'
      case ('lbw')
        my_chnstr = '1s0'
      case ('syl')
        my_chnstr = '1s0'
      case ('tly')
        my_chnstr = '1s0'
      case default
        my_chnstr = command_str(1:3)
    end select
    call convert_text_to_lsj(my_chnstr, lsj, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a, a3, a)')                 &
        & 'AllocateNonPertCHN: Conversion of ', my_chnstr,      &
        & ' to an channel ID failed. Nothing will be done.'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case (command_str)
      case ('spr')
        allocate(obj_sprbl1s0::ptrchn)
      case ('lbw')
        allocate(obj_lbwdib::ptrchn)
      case ('syl')
        allocate(obj_1s0SYLV::ptrchn)
      ! case ('tly')
      !   allocate(obj_tlrLY::ptrchn)
      case default
        select case (chnid)
          case (CHNID_3P0_PP)
            allocate(obj_withc0_sngl::ptrchn)
          case (CHNID_1P1_PP, CHNID_3P1_PP)
            allocate(obj_sngl::ptrchn)
          case (CHNID_1S0_PP)
            allocate(obj_1s0LY::ptrchn)
          case (CHNID_3S1_PP, CHNID_3P2_PP)
            allocate(obj_cpld::ptrchn)
          case default
            write (standard_error_unit, '(a, a3)')  &
              & 'AllocateNonPertCHN: Don''t know how to handle ', command_str
            return
        end select
    end select

    is_successful = .true.
    call create_channel(ptrchn, chnid, uptoQn)

  end subroutine AllocateNonPertCHN

  subroutine allocate_create_peripheral_class(ptr_chn, command_str, uptoQn,    &
    & is_successful)

    class(obj_channel), pointer, intent(inout)  :: ptr_chn
    character(len = *), intent(in)              :: command_str
    integer, intent(in)                         :: uptoQn
    logical, intent(out)                        :: is_successful

    integer             :: chnid
    logical             :: succeeded
    type(lsj_symbol)    :: lsj
    character(len = 3)  :: my_chnstr

    is_successful = .false.
    if (associated(ptr_chn)) then
      write (standard_error_unit, '(a)')  &
        & 'allocate_create_peripheral_class: Input pointer must be null. Nothing will be done.'
        return
    end if

    my_chnstr = command_str(1:3)
    call convert_text_to_lsj(my_chnstr, lsj, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a, a3, a)')  &
        & 'allocate_create_peripheral_class: Conversion of ', my_chnstr,       &
        & ' to an channel ID failed. Nothing will be done.'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
      allocate(obj_pope_sngl::ptr_chn)
    else
      allocate(obj_pope_cpld::ptr_chn)
    end if

    is_successful = .true.
    call create_channel(ptr_chn, chnid, uptoQn)

  end subroutine allocate_create_peripheral_class

  subroutine allocate_create_supperi_class(ptr_chn, command_str, uptoQn,       &
    & is_successful)

    class(obj_channel), pointer, intent(inout)  :: ptr_chn
    character(len = *), intent(in)              :: command_str
    integer, intent(in)                         :: uptoQn
    logical, intent(out)                        :: is_successful

    integer             :: chnid
    logical             :: succeeded
    type(lsj_symbol)    :: lsj
    character(len = 3)  :: my_chnstr

    is_successful = .false.
    if (associated(ptr_chn)) then
      write (standard_error_unit, '(a)')  &
        & 'allocate_create_supperi_class: Input pointer must be null. Nothing will be done.'
        return
    end if

    my_chnstr = command_str(1:3)
    call convert_text_to_lsj(my_chnstr, lsj, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a, a3, a)')  &
        & 'allocate_create_supperi_class: Conversion of ', my_chnstr,          &
        & ' to an channel ID failed. Nothing will be done.'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
      allocate(obj_SUPpope_sngl::ptr_chn)
    else
      allocate(obj_SUPpope_cpld::ptr_chn)
    end if

    is_successful = .true.
    call create_channel(ptr_chn, chnid, uptoQn)

  end subroutine allocate_create_supperi_class

  ! channel must have been created.
  subroutine get_output_for_Lambda_list(channel, N, regtype, output_option,    &
    & fxk_is_needed)

    class(obj_channel), intent(inout)   :: channel
    integer, intent(in)                 :: N, regtype
    integer, intent(in)                 :: output_option
    logical, intent(in)                 :: fxk_is_needed

!     integer, parameter  :: numks_max = 2000, nLmbds_max = 2000,   &
!       & num_cnttrms_max = 100
    integer             :: funt, numks, nLmbds, num_cnttrms, ii
    character(len = 10) :: outstr
    logical             :: succeeded

    real(NER)             ::                                &
      & kk(1:numks_max),                                    &
      & Mambda(1:nLmbds_max),                               &
      & Mambda_cnttrms(1:nLmbds_max, 1:num_cnttrms_max+1)
    character(len = 256)  ::                                &
      & text_output_k(1:numks_max),                         &
      & text_output_lmbd_k(1:nLmbds_max, 1:numks_max)

    Mambda_cnttrms = 0.0_NER
    funt = get_free_unit_io()
    call open_klist(channel%chnstr, funt, numks, kk, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a)')  &
      & 'get_output_for_Lambda_list: Couldn''t open the klist file. Nothing will be done.'
      return
    end if
    if (numks > numks_max) then
      write (standard_error_unit, '(a)')  &
      & 'get_output_for_Lambda_list: the number of k''s exceeds the limit.'
      return
    end if

    call channel%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
    if  (.not. succeeded) then
      write (standard_error_unit, '(a)')  &
      & 'get_output_for_Lambda_list: Couldn''t open the counterterm/Lambda file. Nothing will be done.'
      return
    end if
    if (nLmbds > nLmbds_max .or. nLmbds == 0) then
      write (standard_error_unit, '(a)')  &
      & 'get_output_for_Lambda_list: the number of Lambda''s is zero or exceeds the limit. Nothing will be done.'
      return
    end if
    if (num_cnttrms > num_cnttrms_max) then
      write (standard_error_unit, '(a)')  &
      & 'get_output_for_Lambda_list: the number of counterterms exceeds the limit. Nothing will be done.'
      return
    end if
    !debug
    outstr = channel%chnstr
    Mambda(1:nLmbds) = Mambda_cnttrms(1:nLmbds, 1)
    do ii = 1, nLmbds
      call channel%get_output_formatted_text_for_klist(N, regtype, Mambda(ii), &
        & kk, numks, Mambda_cnttrms(ii, 2:num_cnttrms+1), num_cnttrms,         &
        & output_option, text_output_k)
      call write_fxLmbd_without_prefix(outstr, funt, kk, numks, Mambda(ii),    &
        & text_output_k)
      text_output_lmbd_k(ii, 1:numks) = text_output_k(1:numks)
    end do

    call channel%get_writeout_string(outstr)
    if (fxk_is_needed) then
      call write_fxk_without_prefix(outstr, funt, kk, numks, Mambda, nLmbds,   &
        & text_output_lmbd_k)
    end if

  end subroutine get_output_for_Lambda_list

  !================================
  ! Subs working with obj_smplchn
  !================================

  ! phaselst(ii, n) is the phase shift correction for uptoQn=(n-1) at kk(ii)
  ! Eg., phaselst(ii, 1) is the LO, phaselst(ii, 2) is the
  ! first CORRECTION, and so on

  subroutine get_phase_as_value_for_klist_smpl(chn, N, regtype, Mambda, kk,    &
    & numk, inputs, num_inputs, phaselst)

    class(obj_smplchn), intent(inout) :: chn
    integer, intent(in)               :: regtype, N, numk, num_inputs
    real(NER), intent(in)             :: Mambda
    real(NER), intent(in)             :: kk(:), inputs(:)
    class(*), intent(inout)             :: phaselst(:, :)

    integer :: ii

    call chn%load_inputs(num_inputs, inputs)
    call chn%init(N, regtype, Mambda)
    do ii = 1, numk
      call chn%update_k(kk(ii))
      call chn%get_allTQn()
      call chn%convert_TQn_to_phase()
      call chn%get_incremental_phsn_value(phaselst(ii, 1:chn%uptoQn+1))
    end do

    call chn%release()

  end subroutine get_phase_as_value_for_klist_smpl

  subroutine get_phase_as_value_for_klist_forBD3s1_smpl(chn, N, regtype, Mambda, kk,    &
    & numk, inputs, num_inputs, phaselst, paralstout)

    class(obj_smplchn), intent(inout) :: chn
    integer, intent(in)               :: regtype, N, numk, num_inputs
    real(NER), intent(in)             :: Mambda
    real(NER), intent(in)             :: kk(:), inputs(:)
    class(*), intent(inout)           :: phaselst(:, :)
    real(NER),intent(out)             :: paralstout(:)

    integer :: ii

    call chn%load_inputs(num_inputs, inputs)
    call chn%init(N, regtype, Mambda)
    call reset_para_smplchn(chn)
    paralstout(1:chn%npar(chn%uptoQn+1)) = chn%potpara(1:chn%npar(chn%uptoQn+1))
    call SetUpVfunc_mmwly_smpl(chn)
    do ii = 1, numk
      call chn%update_k(kk(ii))
      call chn%get_allTQn()
      call chn%convert_TQn_to_phase()
      call chn%get_incremental_phsn_value(phaselst(ii, 1:chn%uptoQn+1))
    end do

    call chn%release()

  end subroutine get_phase_as_value_for_klist_forBD3s1_smpl

  subroutine get_phase_as_value_for_klist_forSCL1s0_smpl(chn, N, regtype, Mambda, kk,    &
    & numk, inputs, num_inputs, phaselst)

    class(obj_smplchn), intent(inout) :: chn
    integer, intent(in)               :: regtype, N, numk, num_inputs
    real(NER), intent(in)             :: Mambda
    real(NER), intent(in)             :: kk(:), inputs(:)
    class(*), intent(inout)             :: phaselst(:, :)

    integer :: ii

    call chn%load_inputs(num_inputs, inputs)
    call chn%init(N, regtype, Mambda)
!     debug
    print *, '==============================='
    print *, 'paralst_in =', inputs
    call reset_para_smplchn(chn)
    print *, 'paralst_out =', chn%potpara(1:6)
    call SetUpVfunc_mmwly_smpl(chn)
    do ii = 1, numk
      call chn%update_k(kk(ii))
      call chn%get_allTQn()
      call chn%convert_TQn_to_phase()
      call chn%get_incremental_phsn_value(phaselst(ii, 1:chn%uptoQn+1))
    end do

    call chn%release()

  end subroutine get_phase_as_value_for_klist_forSCL1s0_smpl

  subroutine get_phase_Tmtrx_as_value_for_klist_smpl(chn, N, regtype, Mambda, kk,    &
    & numk, inputs, num_inputs, phaselst, Tmtrx)

    class(obj_smplchn), intent(inout) :: chn
    integer, intent(in)               :: regtype, N, numk, num_inputs
    real(NER), intent(in)             :: Mambda
    real(NER), intent(in)             :: kk(:), inputs(:)
    real(NER), intent(inout)        :: phaselst(:, :)
    real(NER), intent(out)            :: Tmtrx(:, :)
    real(NER), allocatable                   :: phasesngl(:, :)
    type(triplet_phs_cpld), allocatable      :: phasecpld(:, :)
    type(lsj_symbol)      :: lsj
    logical               :: succeeded

    integer :: ii

    call chn%load_inputs(num_inputs, inputs)
    call chn%init(N, regtype, Mambda)
    call convert_chnid_to_lsj(chn%chnid, lsj, succeeded)
!     debug
!     print*, "---------------------"
!     print*, 'potpara =', inputs(1:4)
    allocate(phasesngl(1:numk, 1:chn%uptoQn+1))
    allocate(phasecpld(1:numk, 1:chn%uptoQn+1))

    phaselst = 0.0_NER
    do ii = 1, numk
      call chn%update_k(kk(ii))
      call chn%get_allTQn()
      call chn%convert_TQn_to_phase()
      Tmtrx(ii,1:chn%uptoQn+1) = real(chn%TQn(1,1,1:chn%uptoQn+1))
      if (lsj%L /= lsj%j .and. lsj%j /= 0) then
        call chn%get_incremental_phsn_value(phasecpld(ii, 1:chn%uptoQn+1))
        phaselst(ii, 1:chn%uptoQn+1) = phasecpld(ii, 1:chn%uptoQn+1)%d1
        phaselst(ii+numk, 1:chn%uptoQn+1) = phasecpld(ii, 1:chn%uptoQn+1)%d2
        phaselst(ii+numk*2, 1:chn%uptoQn+1) = phasecpld(ii, 1:chn%uptoQn+1)%e
      else
        call chn%get_incremental_phsn_value(phasesngl(ii, 1:chn%uptoQn+1))
        phaselst(ii, 1:chn%uptoQn+1) = phasesngl(ii, 1:chn%uptoQn+1)
        phaselst(ii+numk, 1:chn%uptoQn+1) = 0.0_NER
        phaselst(ii+numk*2, 1:chn%uptoQn+1) = 0.0_NER
      end if
    end do
    deallocate(phasesngl, phasecpld)

    call chn%release()

  end subroutine get_phase_Tmtrx_as_value_for_klist_smpl

  ! smplchn must have been created

  subroutine get_phase_as_text_for_klist_smpl(smplchn, N, regtype, Mambda, kk, &
    & numk, inputs, num_inputs, phase_out_text)

    class(obj_smplchn), intent(inout) :: smplchn
    integer, intent(in)               :: regtype, N, numk, num_inputs
    real(NER), intent(in)             :: Mambda
    real(NER), intent(in)             :: kk(:), inputs(:)
    character(len=*), intent(inout)   :: phase_out_text(:)

    integer :: ii

    call smplchn%load_inputs(num_inputs, inputs)
    call smplchn%init(N, regtype, Mambda)
!     debug use vbar to calculate phase shift in 3s1
!     call init_Rmesh(smplchn%N, smplchn%msh, smplchn%wght)
    do ii = 1, numk
      call smplchn%update_k(kk(ii))
      call smplchn%get_allTQn()
      call smplchn%convert_TQn_to_phase()
      call smplchn%convert_phsn_to_text(phase_out_text(ii))
    end do
    call smplchn%release()
    call smplchn%reset_VLMinited()
!     debug
!     call relase_spu_flag

  end subroutine get_phase_as_text_for_klist_smpl

  ! smplchn must have been created.

  subroutine get_phase_as_text_for_MmbdkcmList_smpl(smplchn, N, regtype,       &
    & fxk_is_needed, preamble_func)

    class(obj_smplchn), intent(inout)   :: smplchn
    integer, intent(in)                 :: N, regtype
    logical, intent(in)                 :: fxk_is_needed
    procedure(prefix_callback)          :: preamble_func

!     integer, parameter  :: numks_max = 1000, nLmbds_max = 1000,   &
!       & num_cnttrms_max = 100
    integer             :: funt, numks, nLmbds, num_cnttrms, ii
    character(len = 10) :: outstr
    logical             :: succeeded

    real(NER)             ::                                &
      & kk(1:numks_max),                                    &
      & Mambda(1:nLmbds_max),                               &
      & Mambda_cnttrms(1:nLmbds_max, 1:num_cnttrms_max+1)
    character(len = 256)  ::                                &
      & text_output_k(1:numks_max),                         &
      & text_output_lmbd_k(1:nLmbds_max, 1:numks_max)

    Mambda_cnttrms = 0.0_NER
    funt = get_free_unit_io()
    call open_klist(smplchn%chnstr, funt, numks, kk, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a)')  &
      & 'get_phase_as_text_for_MmbdkcmList_smpl: Couldn''t open the klist file. Nothing will be done.'
      return
    end if
    if (numks > numks_max) then
      write (standard_error_unit, '(a)')  &
      & 'get_phase_as_text_for_MmbdkcmList_smpl: the number of k''s exceeds the limit.'
      return
    end if
    call smplchn%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
    if  (.not. succeeded) then
      write (standard_error_unit, '(a)')  &
      & 'get_phase_as_text_for_MmbdkcmList_smpl: Couldn''t open the counterterm/Lambda file. Nothing will be done.'
      return
    end if
    if (nLmbds > nLmbds_max .or. nLmbds == 0) then
      write (standard_error_unit, '(a)')  &
      & 'get_phase_as_text_for_MmbdkcmList_smpl: the number of Lambda''s is zero or exceeds the limit. Nothing will be done.'
      return
    end if
    if (num_cnttrms > num_cnttrms_max) then
      write (standard_error_unit, '(a)')  &
      & 'get_phase_as_text_for_MmbdkcmList_smpl: the number of counterterms exceeds the limit. Nothing will be done.'
      return
    end if
    Mambda(1:nLmbds) = Mambda_cnttrms(1:nLmbds, 1)
    do ii = 1, nLmbds
      call get_phase_as_text_for_klist_smpl(smplchn, N, regtype, Mambda(ii),   &
        & kk, numks, Mambda_cnttrms(ii, 2:num_cnttrms+1), num_cnttrms,         &
        & text_output_k)
      call write_fxLmbd_with_prefix(smplchn%chnstr, funt, kk, numks,           &
        & Mambda(ii), text_output_k, preamble_func)
      text_output_lmbd_k(ii, 1:numks) = text_output_k(1:numks)
    end do

    ! call write_fxLmbd_with_prefix(smplchn%chnstr, funt, kk, numks, Mambda,  &
      ! & nLmbds, text_output_lmbd_k, preamble_func)
    if (fxk_is_needed) then
      call write_fxk_with_prefix(smplchn%chnstr, funt, kk, numks, Mambda,   &
        & nLmbds, text_output_lmbd_k, preamble_func)
    end if

  end subroutine get_phase_as_text_for_MmbdkcmList_smpl

  ! using potential in chengdu as input instead

  subroutine get_phase_as_text_for_chengdu_smpl(smplchn, N, regtype,       &
    & fxk_is_needed, preamble_func, Mambda_in ,pot)

    class(obj_smplchn), intent(inout)   :: smplchn
    integer, intent(in)                 :: N, regtype
    logical, intent(in)                 :: fxk_is_needed
    real(NER), intent(in)               :: Mambda_in
    character(len = 20), intent(in)     :: pot
    procedure(prefix_callback)          :: preamble_func

!     integer, parameter  :: numks_max = 1000, nLmbds_max = 1000,   &
!       & num_cnttrms_max = 100
    integer             :: funt, numks, nLmbds, num_cnttrms, ii
    character(len = 10) :: outstr
    logical             :: succeeded

    real(NER)             ::                                &
      & kk(1:numks_max),                                    &
      & Mambda(1:nLmbds_max),                               &
      & Mambda_cnttrms(1:nLmbds_max, 1:num_cnttrms_max+1)
    character(len = 256)  ::                                &
      & text_output_k(1:numks_max)

    funt = get_free_unit_io()
    call open_klist(smplchn%chnstr, funt, numks, kk, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a)')  &
      & 'get_phase_as_text_for_chengdu_smpl: Couldn''t open the klist file. Nothing will be done.'
      return
    end if
    if (numks > numks_max) then
      write (standard_error_unit, '(a)')  &
      & 'get_phase_as_text_for_chengdu_smpl: the number of k''s exceeds the limit.'
      return
    end if

    nLmbds = 1
    Mambda_cnttrms(1:nLmbds, 1) = Mambda_in
    num_cnttrms = 0
    Mambda(1:nLmbds) = Mambda_cnttrms(1:nLmbds, 1)
    ii = 1
      call get_phase_as_text_for_klist_smpl(smplchn, N, regtype, Mambda(ii),   &
        & kk, numks, Mambda_cnttrms(ii, 2:num_cnttrms+1), num_cnttrms,         &
        & text_output_k)
    call write_fxLmbd_with_prefix_chengdu(smplchn%chnstr, funt, kk, numks,           &
      & pot, text_output_k, preamble_func)

  end subroutine get_phase_as_text_for_chengdu_smpl

  !================================!
  ! Utility subroutines            !
  !================================!

  ! If setn_cis is not set by command line, -999 will be returned
  subroutine parse_cmdln_getprintout(N, chnstr, uptoQn, fxk_is_needed,         &
    & setn_cis, succeeded)

    integer, intent(out)  :: N, uptoQn, setn_cis
    character(len = *), intent(out) :: chnstr
    logical, intent(out)  :: fxk_is_needed, succeeded

    integer     :: cmdcnt, ii, ios
    character(len = 1024)   :: cmdbuff, tmp
    character(len = 10)     :: opt

    succeeded = .true.
    cmdbuff = ''
    cmdcnt = command_argument_count()
    do ii = 1, cmdcnt
      call get_command_argument(ii, tmp)
      cmdbuff = trim(cmdbuff)//' '//trim(tmp)
    end do

    read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt
!     if  (ios /= 0) then
!       read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt, setn_cis
!       if (ios /= 0) then
!         write (standard_error_unit, '(a)')  &
!           & 'parse_cmdln_getprintout: Something wrong with the arguments'
!         write (standard_error_unit, '(a)')  &
!           & 'Usage: get_something N chn uptoQn fxkopt [setn_cis]'
!         write (standard_error_unit, '(a)')  &
!           & 'fxkopt = nofxk or fxk'
!         succeeded = .false.
!         return
!       end if
!     else
!       setn_cis = -999
!     end if
  if (uptoQn >= 3)  then
      read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt, setn_cis
      if (ios /= 0) then
        write (standard_error_unit, '(a)')  &
          & 'parse_cmdln_getprintout: Something wrong with the arguments'
        write (standard_error_unit, '(a)')  &
          & 'Usage: get_something N chn uptoQn fxkopt [setn_cis]'
        write (standard_error_unit, '(a)')  &
          & 'fxkopt = nofxk or fxk'
        succeeded = .false.
        return
      end if
  end if

    select case (trim(opt))
      case ('fxk')
        fxk_is_needed = .true.
      case ('nofxk')
        fxk_is_needed = .false.
      case default
        write (standard_error_unit, '(a)')  &
          & 'parse_cmdln_getprintout: unknown option '//trim(opt)
        succeeded = .false.
        return
    end select

    ! if (uptoQn >= 3) then
    !   read (cmdbuff, *, IOSTAT=ios) N, chnstr, uptoQn, opt, setn_cis
    !   if  (ios /= 0) then
    !     write (standard_error_unit, '(a)')  &
    !       & 'parse_cmdln_getprintout: for uptoQn = 3, setn_cis is needed'
    !     write (standard_error_unit, '(a)')  &
    !       & 'Usage: get_something N chn uptoQn opt [setn_cis]'
    !     succeeded = .false.
    !   end if
    ! end if

  end subroutine parse_cmdln_getprintout

  ! output fxk files
  subroutine write_fxk_with_prefix(chnstr, funt, kk, numks, Mambda, nLmbds,    &
    & text_output_lmbd_k, preamble)

    character(len = *), intent(in)          :: chnstr
    integer, intent(in)                     :: funt, numks, nLmbds
    real(NER), dimension(1:), intent(in)    :: kk, Mambda
    character(len = *), dimension(1:, 1:), intent(in)   :: text_output_lmbd_k
    procedure(prefix_callback)              :: preamble

    character(len = 64)    :: fname, timestamp
    integer :: ii, jj

    call get_timestring_hms(timestamp)
    do jj = 1, numks
      write (fname, '(i3.3)') int(kk(jj))
      fname = adjustl(fname)
      fname = 'run_'//trim(chnstr)//'_fxk_'//trim(fname)//'.out'
      open (funt, file=fname, status='replace')
      write (funt, '(a)') '# '//trim(timestamp)
      call preamble(funt)
      write (funt, '(a, f6.2)') '# kcm = ', kk(jj)
      write (funt, '(a)') '# Mambda(i),   Mambda(1)/Mambda(i) ...'
      do ii = 1, nLmbds
        write (funt, '(f8.1, 2x, f8.5, 4x, a)') Mambda(ii),                    &
          & Mambda(1)/Mambda(ii), trim(text_output_lmbd_k(ii, jj))
      end do
      close(funt)
    end do

  end subroutine write_fxk_with_prefix

  subroutine write_fxk_without_prefix(chnstr, funt, kk, numks, Mambda, nLmbds, text_output_lmbd_k)

    character(len = *), intent(in)          :: chnstr
    integer, intent(in)                     :: funt, numks, nLmbds
    real(NER), dimension(1:), intent(in)    :: kk, Mambda
    character(len = *), dimension(1:, 1:), intent(in)   :: text_output_lmbd_k

    call write_fxk_with_prefix(chnstr, funt, kk, numks, Mambda, nLmbds,        &
      & text_output_lmbd_k, no_preamble)

  end subroutine write_fxk_without_prefix

  subroutine write_fxLmbd_with_prefix(chnstr, funt, kk, numks, Mambda,         &
    & text_output_k, preamble)

    character(len = *), intent(in)          :: chnstr
    integer, intent(in)                     :: funt, numks
    real(NER), intent(in)                   :: kk(:), Mambda
    character(len = *), intent(in)          :: text_output_k(:)
    procedure(prefix_callback)              :: preamble

    character(len = 64)    :: fname, timestamp
    integer :: jj

    call get_timestring_hms(timestamp)
    write (fname, '(i6.4)') int(Mambda)
    fname = adjustl(fname)
    fname = 'run_'//trim(chnstr)//'_lmbd_'//trim(fname)//'.out'
    open (funt, file=fname, status='replace')
    write (funt, '(a)') '# '//trim(timestamp)
    call preamble(funt)
    write(funt, '(a, f8.1)') '# Mambda = ', Mambda
    call write_to_file_basic_phycnst(funt)
    write(funt, '(a)') '# kcm, Tlab, LO, NLO, N2LO, N3LO, ...'
    do jj = 1, numks
      write (funt, '(2(f10.3, 3x), a)') kk(jj), k2Tlab(kk(jj)),  &
        & trim(text_output_k(jj))
    end do
    close(funt)

  end subroutine write_fxLmbd_with_prefix

  subroutine write_fxLmbd_with_prefix_chengdu(chnstr, funt, kk, numks, pot,         &
    & text_output_k, preamble)

    character(len = *), intent(in)          :: chnstr
    integer, intent(in)                     :: funt, numks
    real(NER), intent(in)                   :: kk(:)
    character(len = 20), intent(in)         :: pot
    character(len = *), intent(in)          :: text_output_k(:)
    procedure(prefix_callback)              :: preamble

    character(len = 64)    :: fname, timestamp
    integer :: jj

    call get_timestring_hms(timestamp)
    fname = trim(chnstr)//trim(pot)//'_phase.out'
    open (funt, file=fname, status='replace')
    write (funt, '(a)') '# '//trim(timestamp)
    call preamble(funt)
    write(funt, '(a, f8.1)') '# pot = '//trim(pot)
    call write_to_file_basic_phycnst(funt)
    write(funt, '(a)') '# kcm, Tlab, LO, NLO, N2LO, N3LO, ...'
    do jj = 1, numks
      write (funt, '(2(f10.3, 3x), a)') kk(jj), k2Tlab(kk(jj)),  &
        & trim(text_output_k(jj))
    end do
    close(funt)

  end subroutine write_fxLmbd_with_prefix_chengdu

  ! subroutine write_fxLmbd_with_prefix(chnstr, funt, kk, numks, Mambda, nLmbds, &
  !   & text_output_lmbd_k, preamble)

  !   character(len = *), intent(in)          :: chnstr
  !   integer, intent(in)                     :: funt, numks, nLmbds
  !   real(NER), dimension(1:), intent(in)    :: kk, Mambda
  !   character(len = *), dimension(1:, 1:), intent(in)   :: text_output_lmbd_k
  !   procedure(prefix_callback)              :: preamble

  !   character(len = 64)    :: fname, timestamp
  !   integer :: ii, jj

  !   call get_timestring_hms(timestamp)
  !   do ii = 1, nLmbds
  !     write (fname, '(i6.4)') int(Mambda(ii))
  !     fname = adjustl(fname)
  !     fname = 'run_'//trim(chnstr)//'_lmbd_'//trim(fname)//'.out'
  !     open (funt, file=fname, status='replace')
  !     write (funt, '(a)') '# '//trim(timestamp)
  !     call preamble(funt)
  !     write(funt, '(a, f8.1)') '# Mambda = ', Mambda(ii)
  !     call write_to_file_basic_phycnst(funt)
  !     write(funt, '(a)') '# kcm, Tlab, LO, NLO, N2LO, N3LO, ...'
  !     do jj = 1, numks
  !       write (funt, '(2(f10.3, 3x), a)') kk(jj), k2Tlab(kk(jj)),  &
  !         & trim(text_output_lmbd_k(ii, jj))
  !     end do
  !     close(funt)
  !   end do

  ! end subroutine write_fxLmbd_with_prefix

  subroutine write_fxLmbd_without_prefix(chnstr, funt, kk, numks, Mambda,      &
    & text_output_k)

    character(len = *), intent(in)          :: chnstr
    integer, intent(in)                     :: funt, numks
    real(NER), intent(in)                   :: kk(:), Mambda
    character(len = *), intent(in)          :: text_output_k(:)

    call write_fxLmbd_with_prefix(chnstr, funt, kk, numks, Mambda,             &
      & text_output_k, no_preamble)

  end subroutine write_fxLmbd_without_prefix

  subroutine open_klist(chnstr, funt, numks, kk, lgc_succeeded)

    character(len = *), intent(in)  :: chnstr
    integer, intent(in)             :: funt
    real(NER), dimension(1:), intent(out) :: kk
    integer, intent(out)            :: numks
    logical, intent(out)            :: lgc_succeeded

    integer                 :: flag
    character(len = 512)    :: fname

    fname = 'run_klist.in'
    call read_1col_io(funt, fname, kk, numks, flag)
    if (flag /= 0) then
      write (standard_error_unit, '(a)') 'couldn''t open '//trim(fname)
      lgc_succeeded = .false.
      return
    end if
    lgc_succeeded = .true.

  end subroutine open_klist

  subroutine open_tlablist(chnstr, funt, numks, kk, lgc_succeeded)

    character(len = *), intent(in)  :: chnstr
    integer, intent(in)             :: funt
    real(NER), dimension(1:), intent(out) :: kk
    integer, intent(out)            :: numks
    logical, intent(out)            :: lgc_succeeded

    integer                 :: flag, ii
    character(len = 512)    :: fname

    fname = 'run_'//trim(chnstr)//'_tlablist.in'
    call read_1col_io(funt, fname, kk, numks, flag)
    if (flag .ne. 0) then
      write (standard_error_unit, '(a)') 'couldn''t open '//trim(fname)
      lgc_succeeded = .false.
      return
    end if
    lgc_succeeded = .true.
    do ii = 1, numks
      kk(ii) = tlab2k(kk(ii))
    end do

  end subroutine open_tlablist

  subroutine no_preamble(funt)

    integer, intent(in) :: funt

  end subroutine no_preamble


end module drv_pwphase
