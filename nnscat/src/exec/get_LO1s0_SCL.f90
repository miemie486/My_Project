! Rui 10/27/2020
! calculate 1s0 LO scattering length

program get_LO1s0_SCL

  use mod_obj_smplchn
  use potwrap
  use util_io, only: get_free_unit_io
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  class(obj_smplchn), allocatable :: ptr_chn

  logical                     :: succeeded
  real(NER)                   :: C0, Mambda, sclength(1), Mambda_cnttrms(1:50, 1:20), fc
  integer                     :: jj, N, chnid, uptoqn, nLmbds, funt, num_cnttrms, kk
  real(NER)                   :: t0, t1

  call cpu_time(t0)

  allocate(obj_smplchn::ptr_chn)

  ptr_chn%pVfunc_V0 => VLO_withc0_epwrap
  N = 200
  uptoqn = 0
  ptr_chn%NCH = 1
  ptr_chn%npar(1:3) = (/1, 3, 6/)
  ptr_chn%k = 1.0E-12_NER
  ptr_chn%k_defined = .true.

  call convert_text_to_chnid('1s0', chnid, succeeded)
  call create_smplchn(ptr_chn, chnid, uptoQn)
  funt = get_free_unit_io()

  call ptr_chn%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
  do jj = 1, nLmbds
    Mambda = Mambda_cnttrms(jj, 1)
    ptr_chn%potpara(1) = Mambda_cnttrms(jj, 2)
    C0 =  Mambda_cnttrms(jj, 2)
    call ptr_chn%init(N, regtype, Mambda)
    call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
    print *, 'Mambda = ', Mambda
    print*, '//////////////////////////////////////'
    print *, 'init_C0 = ', ptr_chn%potpara(1)
    print *, 'LO_sclength1 = ', sclength(1)
!     do kk = 1, 5
!       fc = 10.0_NER**(-5-kk)
!       print *, 'fc = ', fc
!       ptr_chn%potpara(1) = C0*(1.0_NER-fc)
!       call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!     print *, 'C0 = ', ptr_chn%potpara(1)
!       print *, 'LO_sclength2 = ', sclength(1)
!       print *, '====================================='
!     end do

  end do
  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
end program get_LO1s0_SCL