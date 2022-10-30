! pray 21/12/2020.

program get_3p0_wfs

  use vbar_3p0_wfs
  use util_gauleg
  use util_io, only: get_free_unit_io
  use mod_obj_smplchn

  implicit none

  integer, parameter :: N = 100
  integer   :: ii, jj, NPi, funt, num_cnttrms, nLmbds, kk, ll
  real(NER) :: Mambda, CMRatio, Mambda_cnttrms(1:50, 1:20)
  real(NER) :: para(1:10)
  real(NEC)  :: f, h2BE, norm, eta(100), tmp_e(100)
  real(NEC)  :: Rmsh(N), Rwght(N), Cwf(2*N), wfvals(2*N), plst(N)
  class(obj_smplchn), allocatable :: chn
  logical    :: succeeded
  real(NER)                   :: t0, t1

  call cpu_time(t0)


  print*, 'Nmesh =', N
  funt = get_free_unit_io()
  allocate(obj_smplchn::chn)
  chn%chnstr = '3p0'
  chn%npar(1) = 1
  chn%uptoQn = 0
  call chn%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
!   print *, 'num_cnttrms =', num_cnttrms
!   print *, 'nLmbds= ',nLmbds
do kk = 1, nLmbds
  Mambda = Mambda_cnttrms(kk, 1)
  para(1:num_cnttrms) = Mambda_cnttrms(kk, 2:num_cnttrms+1)
!   print *, 'para = ', para(1:num_cnttrms)
  call feb_tanabsc(N, Mambda*0.5_NER, Rmsh, Rwght)

  call get_h2wfs_real_wrap(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, PC_mN, VLO_withc0_epwrap,      &
    &  h2BE, Cwf)

  print *, '===================================================================='
  print *, 'Mambda = ', Mambda
  print *, 'para = ', para(1:num_cnttrms)
  print *, "h2BE1 = ", h2BE

  call get_h2wfs_real_wrap(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, PC_mN, VbarLO_3p0_withc0_chengdu_epwrap,      &
    &  h2BE, Cwf)
  print *, "h2BE2 = ", h2BE

!   call get_eta_3p0_real(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, PC_mN, VLO_withc0_epwrap, h2BE, Cwf, tmp_e, eta)

  call get_eta_3p0_real(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, PC_mN, VbarLO_3p0_withc0_chengdu_epwrap, h2BE, Cwf, tmp_e, eta)

  do ll = 1, 100
    print *, "E =",  tmp_e(ll), "eta = ",  eta(ll)
  end do

  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
end do

end program get_3p0_wfs