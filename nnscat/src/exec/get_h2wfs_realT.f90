!  pray 2020/08/20

program get_h2wfs_realT

  use eft_tmtrx
  use util_gauleg
  use util_io, only: get_free_unit_io
  use mod_obj_smplchn
  implicit none

  integer, parameter :: N = 100
  integer   :: ii, jj, NPi, funt, num_cnttrms, nLmbds, kk, ll
  real(NER) :: Mambda, CMRatio, phi, Mambda_cnttrms(1:50, 1:20)
  real(NER) :: msh(N), wght(N), extrapara(NUM_PARA_TWOBODY), para(1:20)
  real(NEC)  :: f, h2BE, norm, eta(100), tmp_e(100)
  real(NEC)  :: Rmsh(N), Rwght(N), Cwf(2*N), wfvals(2*N), plst(N)
  class(obj_smplchn), allocatable :: chn
  logical    :: succeeded
  real(NER)                   :: t0, t1

  call cpu_time(t0)


!   Mambda = 400.0_NER
  CMRatio = 0.5_NER
  phi = 5.0_NER
  ! phi = 0.0_NER
!   f = exp(-IMGUNT_NE*phi/180.0*PI_NE)
  f = 1.0_NER
  print*, 'Nmesh =', N
  funt = get_free_unit_io()
  allocate(obj_smplchn::chn)
  chn%chnstr = '3s1'
  chn%npar(1) = 1
  chn%uptoQn = 0
  call chn%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
!   print *, 'num_cnttrms =', num_cnttrms
!   print *, 'nLmbds= ',nLmbds
do kk = 1, nLmbds
  Mambda = Mambda_cnttrms(kk, 1)
  para(1:num_cnttrms) = Mambda_cnttrms(kk, 2:num_cnttrms+1)
!   print *, 'para = ', para(1:num_cnttrms)
  call feb_tanabsc(N, Mambda*CMRatio, msh, wght)
  do ii = 1, N
    Rmsh(ii) = f*msh(ii)
    Rwght(ii) = f*wght(ii)
  end do
  call get_h2wfs_real_potwrap(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, VLO_withc0_epwrap,      &
    & extrapara, h2BE, Cwf)
!   call get_h2wfs_vbd_real_potwrap(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, VLO_withc0_epwrap,      &
!     & extrapara, h2BE, Cwf)
!   call get_eta_no_real(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, PC_mN, VLO_withc0_epwrap, h2BE, Cwf, tmp_e, eta)
  NPi = N
  plst(1:N) = Rmsh(1:N)

!   call get_generic_h2wfvals_real_potwrap(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, NPi, plst,    &
!     & VLO_withc0_epwrap, extrapara, h2BE, wfvals)

  norm = 0.0_NER
  do ii = 1, N
    norm = norm + (Cwf(ii)**2 + Cwf(ii+N)**2)*Rmsh(ii)*Rmsh(ii)*Rwght(ii)
  end do

  print *, '===================================================================='
  print *, 'Mambda = ', Mambda
  print *, 'para = ', para(1:num_cnttrms)
  print *, "h2BE1 = ", h2BE
!   do ll = 1, 100
!     print *, "E =",  tmp_e(ll), "eta = ",  eta(ll)
!   end do
!   print *, "norm = ", norm

  call init_Rmesh(N, Rmsh, Rwght)
!   call get_h2wfs_vbar_real(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, PC_mN, VbarLO_withc0_epwrap,      &
!     &  h2BE, Cwf)
!   print *, "h2BE2_vbar =", h2BE

  call get_h2wfs_vbd_real_potwrap(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, VbarLO_withc0_epwrap,      &
    & extrapara, h2BE, Cwf)
!   call get_eta_no_real(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, PC_mN, VbarLO_withc0_epwrap, h2BE, Cwf, tmp_e, eta)

!   do ll = 1, 100
!     print *, "E =",  tmp_e(ll), "eta = ",  eta(ll)
!   end do
!   do ll = 1, N
!     print *, "Rmsh =" , Rmsh(ll),"Rwght =" , Rwght(ll)
!   end do
  print *, "h2BE2 =", h2BE

!   call get_h2wfs_real_potwrap(Mambda, para, REGTYPE_GAUSSIAN, N, Rmsh, Rwght, VbarLO_withc0_epwrap,      &
!     & extrapara, h2BE, Cwf)
!   print *, "h2BE2 =", h2BE
  call relase_spu_flag()


  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
end do

end program get_h2wfs_realT