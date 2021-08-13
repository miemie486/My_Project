! Bingwei Long 08/12/2013
! Driver routines for solving for the values of counterterms

! Frozen on 03/10/2018
module drv_solve

  use drv_pwphase
  use util_rootfinder
  use util_io
  implicit none

  real(NER), parameter   :: goal_default_drv_solve = 1.0E-5,  &
    & rel_tol_default_drv_solve = 1.0E-7
  integer, parameter      :: read_file_FAIL = 10, read_file_SUCCESS = 20

  ! Used by diff_phase_C0
  class(obj_channel), pointer :: chn_diff_phase_C0 => null()
  real(NER)                   :: k_diff_phase_C0
  real(NER)                   :: phase_diff_phase_C0

contains


  ! Assuming chn_diff_phase_C0 has already been properly initialized.

  function diff_phase_C0(C0)

    real(NER)               :: diff_phase_C0
    real(NER), intent(in)   :: C0

    integer, dimension(1:2) :: uflag

    call chn_diff_phase_C0%change_a_C(1, C0)
    call chn_diff_phase_C0%update_k(k_diff_phase_C0)
    call chn_diff_phase_C0%get_allTQn()
    call chn_diff_phase_C0%convert_TQn_to_phase_quiet(uflag)
    select type (chn_diff_phase_C0)
      class is (obj_sngl)
        diff_phase_C0 = chn_diff_phase_C0%phsn(1) - phase_diff_phase_C0
      class is (obj_cpld)
        diff_phase_C0 = chn_diff_phase_C0%phsn(1)%d1 - phase_diff_phase_C0
    end select

  end function diff_phase_C0


  ! Solve for C0, the value for the leading order counterterm, to fit the
  ! input phase shift at the given. If C0 is not found or there is an error,
  ! num_C0s will return 0.

  subroutine solve_for_C0(chnstr, N, regtype, Mambda, rLtM, k, phase,          &
    & C0_start, C0_end, num_segs, rel_tol, goal, num_C0s, out_C0)

    character(len = *), intent(in)  :: chnstr
    integer, intent(in)             :: N, regtype, num_segs
    real(NER), intent(in)           :: Mambda, rLtM, k, phase, C0_start,       &
      & C0_end, rel_tol, goal
    integer, intent(out)            :: num_C0s
    real(NER), intent(out)          :: out_C0(:)

    integer                 :: jj, bn
    real(NER)               :: C0
    character(len = 30)     :: known_chns
    real(NER), pointer      :: xb1(:), xb2(:)
    class(obj_channel), pointer :: chn

    known_chns = '1s0_3s1_3p0_3p2'
    if (index(known_chns, chnstr) == 0) then
      write (standard_error_unit, '(a, a)')  &
        'solve_for_C0: Do not know how to solve for C0 of ', chnstr
      num_C0s = 0
      return
    end if
    if (associated(chn_diff_phase_C0)) then
      if (chn_diff_phase_C0%created) call chn_diff_phase_C0%erase()
      deallocate(chn_diff_phase_C0)
    end if
    call AllocateNonPertCHN(chn, chnstr, 0)
    chn_diff_phase_C0 => chn
    call chn%init(N, regtype, Mambda, Mambda*rLtM)
    k_diff_phase_C0 = k
    phase_diff_phase_C0 = phase
    call zbrak_pointer(diff_phase_C0, C0_start, C0_end, num_segs, xb1, xb2, bn)
    num_C0s = 0
    do jj = 1, bn
      C0 = zbrent(diff_phase_C0, xb1(jj), xb2(jj), rel_tol*abs(C0_start-C0_end))
      if (abs(diff_phase_C0(C0)/phase) < goal) then
        num_C0s = num_C0s + 1
        out_C0(num_C0s) = C0
      end if
    end do
    deallocate(xb1, xb2)
    nullify(chn_diff_phase_C0)
    call chn%erase()
    deallocate(chn)

  end subroutine solve_for_C0


  ! Format of the input file
  ! 1st line:
  ! Tlab    Phase
  ! Lambda1
  ! Lambda2
  ! Lambda3
  ! ...

  subroutine read_input_file_slvc0(input_fname, k, phase, num_Mambda,          &
    & Mambda_list, C0_start_list, C0_end_list, flag)

    character(len = *), intent(in)      :: input_fname
    real(NER), intent(out)              :: k, phase
    real(NER), dimension(:), intent(out)    :: Mambda_list, C0_start_list,     &
      & C0_end_list
    integer, intent(out)                    :: num_Mambda, flag

    integer :: funt, ii, ios
    logical :: is_end_of_file, succeeded
    real(NER), dimension(1:3)  :: para_buff

    funt = get_free_unit_io()
    flag = read_file_FAIL
    open (funt, file=input_fname, status='old', iostat = ios, action = 'read')
    if (ios .ne. 0) then
      close (funt)
      return
    end if

    call read_next_useful_line_as_real_array(funt, 2, para_buff, succeeded)
    if (.not. succeeded) then
      write (standard_error_unit, '(a)')  &
        & 'read_input_file_slvc0: The first line must be: k   phase   . Nothing will be done.'
      close (funt)
      return
    end if
    k     = Tlab2k(para_buff(1))
    phase = para_buff(2)
    ii = 0
    do
      call read_next_useful_line_as_real_array(funt, 3, para_buff, succeeded)
      if (is_end_of_file) exit
      if (.not. succeeded) cycle
      ii = ii + 1
      Mambda_list(ii)   = para_buff(1)
      C0_start_list(ii) = para_buff(2)
      C0_end_list(ii)   = para_buff(3)
    end do

    close (funt)
    num_Mambda = ii
    if (ii == 0) then
      write (standard_error_unit, '(a)')  &
        & 'read_input_file_slvc0: No Mambdas have been read.'
    else
      flag = read_file_SUCCESS
    end if

  end subroutine read_input_file_slvc0


  subroutine write_output_file_slvc0(output_fname, k, phase, nrows, C0, Mambda)

    character(len = *), intent(in)  :: output_fname
    real(NER), intent(in)           :: k, phase
    integer, intent(in)             :: nrows
    real(NER), dimension(:), intent(in) :: C0, Mambda

    integer :: funt, ii
    real(NER), dimension(1:nrows, 1:2)  :: matrix

    do ii = 1, nrows
      matrix(ii, 1) = Mambda(ii)
      matrix(ii, 2) = C0(ii)
    end do

    funt=get_free_unit_io()
    open (funt, file = output_fname, status = 'replace')
    write (funt, '(a, f8.3, 2x, a, f9.4)') '# k = ', k, 'phase = ', phase
    write (funt, '(a)') '# Mambda    C0'
    call write_ncols_io(funt, 2, nrows, matrix, '*')
    close (funt)

  end subroutine write_output_file_slvc0


  ! Maximum of 512 Mambdas

  subroutine get_C0_for_Mambda_list(chnstr, N, regtype, rLtM, input_fname,     &
    & output_fname, num_segs, rel_tol, goal)

    character(len = 3), intent(in)  :: chnstr
    character(len = *), intent(in)  :: input_fname, output_fname
    integer, intent(in)             :: N, regtype, num_segs
    real(NER), intent(in)           :: rLtM, rel_tol, goal

    integer :: num_Mmbds, ii, file_flag, num_C0s, cnt
    real(NER)   :: k, phase
    real(NER), dimension(1:512) :: C0_1, C0_2, C0, Mmbd, Mmbd_out
    real(NER), dimension(1:10)  :: tmp_C0

    call read_input_file_slvc0(input_fname, k, phase, num_Mmbds, Mmbd, C0_1,   &
      & C0_2, file_flag)
    if (file_flag == read_file_FAIL) then
      write (standard_error_unit, '(a)')  &
        'get_C0_for_Mambda_list: Couldn''t read '//trim(input_fname)
      return
    end if

    cnt = 1
    do ii = 1, num_Mmbds
      call solve_for_C0(chnstr, N, regtype, Mmbd(ii), rLtM, k, phase,          &
        & C0_1(ii), C0_2(ii), num_segs, rel_tol, goal, num_C0s, tmp_C0)
      if (num_C0s > 0) then
        C0(cnt:cnt+num_C0s-1) = tmp_C0(1:num_C0s)
        Mmbd_out(cnt:cnt+num_C0s-1) = Mmbd(ii)
        cnt = cnt + num_C0s
      end if
    end do

    call write_output_file_slvc0(output_fname, k, phase, cnt, C0, Mmbd_out)

  end subroutine get_C0_for_Mambda_list


end module drv_solve
