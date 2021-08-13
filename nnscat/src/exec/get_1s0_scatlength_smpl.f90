!pray  get 1s0 scatlength with MMWLY potential
!empirical scatlength = -23.740(20)
!numpara = (/ 1, 3, 6/), uptoQn = 0,1,2

program get_1s0_scatlength_smpl

  use drv_fitpwa
  use mod_obj_smplchn
  implicit none

  real(NER)   :: Mambda(14), paralst(1:6), k(1:1), scatlength
  integer     :: uptoQn, ii, jj
  character(len = 20)  :: pot(14)

    PC_gA = 1.29_NER
    PC_mN = 939.0_NER
    PC_mpi = 138.0_NER
    PC_fpi = 92.4_NER

  pot(1: 14) = (/'chengdu_MMWLY_350 ','chengdu_MMWLY_400 ','chengdu_MMWLY_600 ', &
      & 'chengdu_MMWLY_800 ','chengdu_MMWLY_1600','chengdu_MMWLY_3200','chengdu_MMWLY_1200',&
      & 'chengdu_MMWLY_2000','chengdu_MMWLY_2400','chengdu_MMWLY_2800', &
      & 'chengdu_MMWLY_3600','chengdu_MMWLY_4000','chengdu_MMWLY_4400','chengdu_MMWLY_4800'/)
  Mambda(1: 14) = (/350.0_NER, 400.0_NER, 600.0_NER, 800.0_NER, 1600.0_NER, 3200.0_NER, 1200.0_NER, 2000.0_NER, &
    & 2400.0_NER, 2800.0_NER, 3600.0_NER, 4000.0_NER, 4400.0_NER, 4800.0_NER/)


do jj = 1, 3
  uptoQn = jj -1
  print *, 'uptoQn = ', uptoQn
  do ii = 1, 14
    print *, 'pot = ', pot(ii)
    call get_scatlength(Mambda(ii), pot(ii), uptoQn, k, scatlength)
    k(1) = 1E-12_NER
    print *, 'k1=', k(1), 'scatlength1 = ' , scatlength
    k(1) = 1E-9_NER
    call get_scatlength(Mambda(ii), pot(ii), uptoQn, k, scatlength)
    print *, 'k2=', k(1), 'scatlength2 = ' , scatlength
    k(1) = 1E-6_NER
    call get_scatlength(Mambda(ii), pot(ii), uptoQn, k, scatlength)
    print *, 'k3=', k(1), 'scatlength3 = ' , scatlength
  end do
end do
contains

	subroutine get_scatlength(Mambda, pot, uptoQn, k, scatlength)

		real(NER), intent(out) :: scatlength

    real(NER)                   :: Mambda, paralst(1:6)
    real(NER), intent(in)       :: k(1:1)
    integer                     :: num_para, Nmesh, chnid, numk
    integer, intent(in)         :: uptoQn
    real(NER)                   :: amplitude(1,uptoQn+1)
    type(lsj_symbol)            :: lsj
    character(len=20), intent(in)       :: pot

    lsj%L = 0
    lsj%J = 0
    lsj%S = 0
    call convert_lsj_to_chnid(lsj, chnid)

!     select case(uptoQn)
!        case(0)
!          num_para = 1
!        case(1)
!          num_para = 3
!        case(2)
!          num_para = 6
!        case default
!          print *, 'uptoQn need to be 0 or 1'
!     end select

    Nmesh = 100
    paralst= 0.0_NER
    num_para = 0
    numk = 1

    call Get1s0_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst, num_para, k, numk, pot, amplitude)

    scatlength = real(sum(amplitude(1,1:uptoQn+1)))*PI_NE/2.0_NER*PC_default_hbarc
    print *, "amplitude" ,sum(amplitude(1,1:uptoQn+1))
	end subroutine get_scatlength

end program get_1s0_scatlength_smpl