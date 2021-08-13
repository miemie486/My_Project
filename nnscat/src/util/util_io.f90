! Bingwei 02/27/2011
! - add support for comment lines, ignoring un-recognizable line

! Bingwei 12/28/2010
! Input/output subroutines

module util_io

  use nneft_type
  implicit none

  ! the maximum size of a line in data files
  private :: line_buff_size

  integer, parameter  ::        &
  &   line_buff_size = 1024,    &
  &   read_ncols_SUCCESS = 0,   &
  &   read_ncols_FAIL    = -1


contains

  ! Search from MIN_UNIT_NUMBER to MAX_UNIT_NUMBER for an available file unit
  ! number and return it.
  function get_free_unit_io()

    integer :: get_free_unit_io

    integer :: unit
    logical :: is_connected
    integer, parameter :: &
    & MIN_UNIT_NUMBER = 10,     &
    & MAX_UNIT_NUMBER = 512

    do unit = MIN_UNIT_NUMBER, MAX_UNIT_NUMBER
      inquire(unit = unit, opened = is_connected)
      if (.not. is_connected) then
        get_free_unit_io = unit
        return
      end if
    end do
    unit = MAX_UNIT_NUMBER
    write(standard_error_unit, '(a)')  &
    & 'get_free_unit_io: No available unit found. 99 is returned.'

  end function get_free_unit_io


  ! Read a line from the given file unit and convert it to a given number of
  ! real parameters. Ignore lines that a) have '#' as the first non-white
  ! character, or b) are not correctly formatted. When succeeded is set false,
  ! the end of file must have been reached.

  subroutine read_next_useful_line_as_real_array(funt, num_para, para_lst, succeeded)

    integer, intent(in)     :: funt, num_para
    real(NER), dimension(:), intent(out)    :: para_lst
    logical, intent(out)    :: succeeded

    integer :: ios
    character (len = line_buff_size) :: buff

    do
      read (funt, '(a)', iostat = ios) buff
      if (ios /= 0) then
        if (ios < 0) then ! end of file detected
          succeeded = .false.
          exit
        else
          ! something went wrong, ignore the current line
          cycle
        end if
      end if
      buff = trim(buff)
      if (buff(1:1) == '#') then
        ! the current line is a comment, ignore
        cycle
      end if
      read (buff, *, iostat = ios) para_lst(1:num_para)
      if (ios /= 0) then
        ! the current line can't be read properly, ignore
        cycle
      else
        succeeded = .true.
        exit
      end if
    end do

  end subroutine read_next_useful_line_as_real_array


  subroutine read_ncols_io(unt, fname, numcols, matrix, numrow, flag)

    integer, intent(in) :: unt, numcols
    character(len=*), intent(in) :: fname
    real(NER), dimension(:, :), intent(out) :: matrix
    integer, intent(out) :: numrow, flag

    integer :: ios, ii
    character (len = line_buff_size) :: buff
    logical :: succeeded_reading_a_line

    open (unt, file=fname, status='old', iostat = ios, action = 'read')
    if (ios .ne. 0) then
      flag = read_ncols_FAIL
      close (unt)
      return
    end if
    ii = 1
    do
      call read_next_useful_line_as_real_array(unt, numcols, matrix(ii, 1:numcols), succeeded_reading_a_line)
      if (.not. succeeded_reading_a_line) exit
      ii = ii + 1
    end do
    close (unt)
    numrow = ii - 1
    flag = read_ncols_SUCCESS

    ! if (numrow > 0) then
    !   flag = read_ncols_SUCCESS
    ! else
    !   flag = read_ncols_FAIL
    ! end if

  end subroutine read_ncols_io


  ! read a single column of data, return number of data, flag of success
  ! flag = read_ncols_SUCCESS -> success; flag != 0 -> fail to open the specified file

  subroutine read_1col_io(unt, fname, col, numrow, flag)

    integer, intent(in) :: unt
    character (len=*), intent(in) :: fname
    real(NER), dimension(:), intent(out) :: col
    integer, intent(out) :: numrow, flag

    real(NER), dimension(:, :), allocatable :: para_buff

    allocate(para_buff(size(col), 1:1))
    call read_ncols_io(unt, fname, 1, para_buff, numrow, flag)
    if (numrow > 0) col(1:numrow) = para_buff(1:numrow, 1)
    deallocate(para_buff)

  end subroutine read_1col_io


  subroutine write_ncols_io(unt, ncols, nrows, matrix, format_string)

    integer, intent(in)         :: unt, ncols, nrows
    real(NER), dimension(:, :)  :: matrix
    character(len=*), intent(in):: format_string

    integer :: ii, jj

    do ii = 1, nrows
      write (unt, format_string) (matrix(ii, jj), jj = 1, ncols)
    end do

  end subroutine write_ncols_io

end module util_io
