!  pray 2020/10/29
! using the pot form chengdu to investigate the values of parameters
! getting binding energy and scattering length with chengdu

program get_bdscl_withCD_smpl

  use drv_pwphase
  use mod_vfunc_smpl
  implicit none

  integer, parameter          :: regtype = REGTYPE_GAUSSIAN, N_mambda = 15
  integer                     :: N, chnid, uptoQn = 0, setn_cis, kk
  class(obj_smplchn), allocatable :: chn
  character(len = 3)          :: chnstr
  character(len = 20)         :: pot(N_mambda)
  logical                     :: fxk_is_needed = .false., successful
  real(NER)                   :: t0, t1, Mambda(N_mambda), pert_BDE(2), BDE1, scl(3)

  call cpu_time(t0)


  pot(1: N_mambda) = (/'chengdu_MMWLY_350 ','chengdu_MMWLY_400 ','chengdu_MMWLY_600 ', &
      & 'chengdu_MMWLY_800 ','chengdu_MMWLY_1600','chengdu_MMWLY_3200','chengdu_MMWLY_1200', &
      & 'chengdu_MMWLY_2000','chengdu_MMWLY_2400','chengdu_MMWLY_2800','chengdu_MMWLY_500 ', &
      & 'chengdu_MMWLY_1800','chengdu_MMWLY_2200','chengdu_MMWLY_4000','chengdu_MMWLY_4800'/)
  Mambda(1: N_mambda) = (/350.0_NER, 400.0_NER, 600.0_NER, 800.0_NER, 1600.0_NER, &
      & 3200.0_NER, 1200.0_NER, 2000.0_NER, 2400.0_NER, 2800.0_NER, 500.0_NER, 1800.0_NER, &
      & 2200.0_NER, 4000.0_NER, 4800.0_NER/)


  do kk = 11, N_mambda
    allocate(obj_smplchn::chn)
    chn%N = 100
    chn%Mambda = Mambda(kk)
    chn%regtype = regtype
    chn%NCH = 1
    chn%npar(1:3) = (/1, 3, 6/)
    chn%potpara(1: 6) = 0.0_NER
      chnstr = '1s0'
      chn%k = 1.0E-11_NER
      chn%k_defined = .true.

      call convert_text_to_chnid(chnstr, chnid, successful)

      call create_smplchn(chn, chnid, 2)
      call chn%init(chn%N, chn%regtype, chn%mambda)
      call SetUpVfunc_mmwly_chengdu_smpl(chn, pot(kk))
        print *, "Mambda = ", Mambda(kk)
      call get_1s0sclength_smplchn(chn, chn%uptoQn, scl)
        print*, 'scattering length = ', scl

    call chn%release()
    call chn%erase()
    deallocate(chn)
  end do

!   do kk = 1, 10
!     allocate(obj_smplchn::chn)
!     chn%N = 100
!     chn%Mambda = Mambda(kk)
!       chnstr = '3s1'
!       call convert_text_to_chnid(chnstr, chnid, successful)
!       call create_smplchn(chn, chnid, 0)
!       call SetUpVfunc_mmwly_chengdu_smpl(chn, pot(kk))

!       call get_BDE_smplchn(chn, 0, 0.0_NER, BDE1)
!         print *, "Mambda = ", Mambda(kk)
!         print *, 'LO_BDE = ', BDE1

!       call create_smplchn(chn, chnid, 2)
!       call get_pertBDE_smplchn(chn, chn%uptoqn, chn%V2_cheb_pertBD_bnd, -chn%V2_cheb_pertBD_bnd, &
!         & chn%V2_cheb_pertBD_n, pert_BDE)
!         print *, 'N2LO_BDE = ', pert_BDE
!     call chn%release()
!     call chn%erase()
!     deallocate(chn)
!   end do


  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'

end program get_bdscl_withCD_smpl
