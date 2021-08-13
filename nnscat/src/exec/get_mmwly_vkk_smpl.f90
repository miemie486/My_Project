! get 3p0 vkk with MMWLY_800

program get_mmwly_vkk_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  real, parameter             :: pr_hbarc = 197.3286_NER
  integer                     :: N, chnid, uptoQn = 0, setn_cis, ii, jj, L, S, J, funt1, num_cnttrms, nLmbds, unit
  class(obj_smplchn), allocatable :: chn
  character(len = 3)          :: chnstr
  character(len = 10)         :: cmdstr
  logical                     :: fxk_is_needed = .false., successful
  real(NER)                   :: t0, t1, tmpVQ0(1:2,1:2), tmpVQ2(1:2,1:2),Mambda_cnttrms(1:1000, 1:100)
  real(NER), allocatable      :: pp(:)
  logical                     :: succeeded
  character(len = 20)         :: cmambda
  character(len = 64)      :: fname

  call cpu_time(t0)

  call parse_cmdln_getprintout(N, cmdstr, uptoQn, fxk_is_needed, setn_cis, successful)
  if (.not. successful) stop
  if (setn_cis > 0) call Select_cis2PE_vtpe1(setn_cis)
  if (setn_cis > 4) call Select_cis2PE_vdel(setn_cis)

  chnstr = cmdstr


  allocate(pp(1:N))
  call convert_text_to_chnid(chnstr, chnid, successful)
  if (.not. successful) then
    write  (standard_error_unit, '(a)')  &
      & 'get_mmwly_smpl: converting '//chnstr//' wasn''t successful'
    stop
  end if

  allocate(obj_smplchn::chn)
  call create_smplchn(chn, chnid, uptoQn)
  call SetUpVfunc_mmwly_smpl(chn)

  L = 1
  S = 1
  J = 0
  Mambda_cnttrms = 0.0_NER
  call chn%read_inputfile(funt1, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
  if  (.not. succeeded) then
    write (standard_error_unit, '(a)')  &
    & 'Couldn''t open the counterterm/Lambda file. Nothing will be done.'
    return
  end if
  do ii = 1, nLmbds
    write (cmambda, '(i4.4)') int(Mambda_cnttrms(ii, 1))
    fname = trim(cmambda)//'_tmtrx.out'
    open(unit = 57+ii, file=fname)
    write(57+ii,'(a)') '#  K  V0  V2'
    do jj = 1, N
      pp(jj) = Mambda_cnttrms(ii, 1)*jj/N
      call chn%pVfunc_V0(L, S, j, regtype, Mambda_cnttrms(ii, 1), Mambda_cnttrms(ii, 2:num_cnttrms+1), pp(jj), pp(jj), tmpVQ0)
      call chn%pVfunc_V2(L, S, j, regtype, Mambda_cnttrms(ii, 1), Mambda_cnttrms(ii, 2:num_cnttrms+1), pp(jj), pp(jj), tmpVQ2)
      write(57+ii,*) pp(jj), tmpVQ0(1,1)*pr_hbarc, tmpVQ2(1,1)*pr_hbarc
    end do
    close(57+ii)
  end do
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

end program get_mmwly_vkk_smpl
