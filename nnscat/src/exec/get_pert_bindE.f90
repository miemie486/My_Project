! Rui 09/12/2020
! calculate perterbative binding energy in 3s1 at N2LO
! run it with inputfile 'inputs_3s1.in'

program get_pert_bindE

  use mod_obj_smplchn
  use util_io, only: get_free_unit_io
  use potwrap
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  class(obj_smplchn), pointer :: ptr_chn => null()

  real(NER)                   :: a, pert_BDE(2), Mambda_cnttrms(1:50, 1:20)
  integer                     :: jj, cheb_n, funt, num_cnttrms, nLmbds, kk
  logical    :: succeeded
  character(len = 20)         :: cmambda

  real(NER)   :: x(20), BDE(20), y(20)
  real(NER)                   :: t0, t1

  call cpu_time(t0)

  do kk = 1, 20
  	x(kk) = 0.001 + kk*(0.004_NER/20)
  end do

!   a = 1.0E-5_NER
!   cheb_n = 4
!   a = 1.0E-8_NER
!   cheb_n = 6

!   a = 1.0E-9_NER
!   cheb_n = 7

!   a = 1.0E-7_NER
!   cheb_n = 5
  a = 1.0E-6_NER
  cheb_n = 5
!   a = 1.0E-10_NER
!   cheb_n = 6

  funt = get_free_unit_io()
  allocate(obj_smplchn::ptr_chn)

  ptr_chn%pVfunc_V0 => VLO_withc0_epwrap
!   ptr_chn%pVfunc_V0 => VbarLO_withc0_epwrap
!   ptr_chn%pVfunc_V2 => VN2LO_withc0_epwrap
  ptr_chn%pVfunc_V2 => VN2LO_withc0_SHRG_epwrap

  ptr_chn%N = 200
!   ptr_chn%N = 100
!   ptr_chn%N = 150
  ptr_chn%regtype = regtype
  ptr_chn%chnstr = '3s1'
  ptr_chn%uptoqn = 2
  ptr_chn%NCH = 2
  ptr_chn%npar(1:3) = (/1, 1, 4/)
  ptr_chn%V2_cheb_pertBD_bnd = a
  ptr_chn%V2_cheb_pertBD_n = cheb_n

  call ptr_chn%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
  print*, 'nLmbds =', nLmbds
  print*, 'Nmesh = ', ptr_chn%N
  print*, 'cheb_n and bnd = ', cheb_n, 'and', a
  do kk = 1, nLmbds
  	ptr_chn%mambda = Mambda_cnttrms(kk, 1)
    print*, 'mambda =', ptr_chn%mambda
    ptr_chn%potpara(1: num_cnttrms) = Mambda_cnttrms(kk, 2:num_cnttrms+1)
    print*, 'paras =', ptr_chn%potpara(1:4)
!     debug
!     call reset_para_smplchn(ptr_chn)
!     print*, 'reset paras =', ptr_chn%potpara(1:4)
!     ptr_chn%pVfunc_V2 => VN2LO_withc0_epwrap
    ptr_chn%pVfunc_V2 => VN2LO_withc0_SHRG_epwrap
    call get_pertBDE_smplchn(ptr_chn, ptr_chn%uptoqn, a, -a, cheb_n, pert_BDE)
    allpert_BDE_initial = .false.
    print*, 'pert_BDE=', pert_BDE
    print*, '---------------------------------------'
!     debug
!     write (cmambda, '(i4.4)') int(ptr_chn%mambda)
!     open(unit=int(57)+kk, file=trim(cmambda)//'.out')
!     do jj =1, 20
!       call get_BDE_smplchn(ptr_chn, ptr_chn%uptoqn, x(jj), BDE(jj))
!       y(jj) = pert_BDE(1) + pert_BDE(2)*x(jj)
!     	write (int(57)+kk,*) x(jj), BDE(jj), y(jj)
!     end do
!     close(unit=int(57)+kk)
  end do

  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
end program get_pert_bindE