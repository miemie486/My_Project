! Rui 10/16/2020
! calculate perterbative scatteringlength in 1s0 at NLO and N2LO
! run it with inputfile 'inputs_1s0.in'

program get_pert_SCL

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
  real(NER),allocatable       :: sclength(:)

  real(NER)                   :: t0, t1, k

  call cpu_time(t0)



  funt = get_free_unit_io()
  allocate(obj_smplchn::ptr_chn)

        ptr_chn%vanishing(1) = .false.
        ptr_chn%npar(1:3) = (/1, 3, 6/)
        ptr_chn%pVfunc_V0 => VLO_withc0_epwrap
        ptr_chn%pVfunc_V1 => VNLO_MMWLY_epwrap
        ptr_chn%V1_cheb_base_n = 4
        ptr_chn%V1_cheb_bnd = 1.0E-3_NER
!         ptr_chn%pVfunc_V2 => VN2LO_MMWLY_epwrap
! debug fit without TPE
        ptr_chn%pVfunc_V2 => VN2LO_MMWLY_SHRG_epwrap
        ptr_chn%V2_cheb_base_n = 4
        ptr_chn%V2_cheb_bnd = 1.0E-7_NER

!   ptr_chn%N = 200
  ptr_chn%N = 100
!   ptr_chn%N = 150
  ptr_chn%regtype = regtype
  ptr_chn%chnstr = '1s0'
  ptr_chn%NCH = 1
  ptr_chn%npar(1:3) = (/1, 3, 6/)
      print *, 'V2_cheb_base_n/bnd =', ptr_chn%V2_cheb_base_n, ' / ', ptr_chn%V2_cheb_bnd

!   call ptr_chn%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
  k = 1.0E-12_NER
!   k = 1.0E-11_NER
!   do jj = 1, 1
  do jj = 2, 2
!   do jj = 0, 0
!     print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++'
!     print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++'
!     jj = 0
    ptr_chn%uptoqn = jj
    call ptr_chn%read_inputfile(funt, num_cnttrms, nLmbds, Mambda_cnttrms, succeeded)
    allocate(sclength(1:ptr_chn%uptoqn+1))
!     do kk = 1, 1
    do kk = 1, nLmbds
  	  ptr_chn%mambda = Mambda_cnttrms(kk, 1)
!       print *, 'num_cnttrms =', num_cnttrms
      print*, 'mambda =', ptr_chn%mambda
      ptr_chn%potpara(1: num_cnttrms) = Mambda_cnttrms(kk, 2:num_cnttrms+1)
      print*, 'paras =', ptr_chn%potpara(1:ptr_chn%npar(jj+1))
      call create_smplchn(ptr_chn, CHNID_1S0_PP, jj)
      call ptr_chn%init(ptr_chn%N, ptr_chn%regtype, ptr_chn%mambda)
      call ptr_chn%update_k(k)

!     debug
        ptr_chn%pVfunc_V1 => VNLO_MMWLY_epwrap
!         ptr_chn%pVfunc_V2 => VN2LO_MMWLY_epwrap
! debug fit without TPE
        ptr_chn%pVfunc_V2 => VN2LO_MMWLY_SHRG_epwrap
      call get_1s0sclength_smplchn(ptr_chn, ptr_chn%uptoQn, sclength)
      call get_NLOpertSCL_smplchn(ptr_chn)
!       print *, 'NLO_scl = ',  NLOpert_SCL
print *, 'TQn(3) =', ptr_chn%TQn(1,1,3)
      print*, 'sclength without resetparas=', sclength
      print*, '---------------------------------------'

      call reset_para_smplchn(ptr_chn)
      print*, 'reset paras =', ptr_chn%potpara(1:ptr_chn%npar(jj+1))
        ptr_chn%pVfunc_V1 => VNLO_MMWLY_epwrap
!         ptr_chn%pVfunc_V2 => VN2LO_MMWLY_epwrap
! debug fit without TPE
        ptr_chn%pVfunc_V2 => VN2LO_MMWLY_SHRG_epwrap
      call ptr_chn%update_k(k)
      call get_1s0sclength_smplchn(ptr_chn, ptr_chn%uptoQn, sclength)
      pert_SCL_initial(jj) = .false.
print *, 'TQn(3) =', ptr_chn%TQn(1,1,3)
      print*, 'sclength=', sclength
      print*, '---------------------------------------'
    end do
    deallocate(sclength)
  end do

  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
end program get_pert_SCL