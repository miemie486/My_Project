
program show_mtrx

  use nneft_type
  use chengdu
  implicit none

  integer       ::  uptoQn, L, S, J
  real(NER)     ::  pout, pin, potval(1:2, 1:2)

  ! Magnitude of in/out momenta
  pout = 2.5_NER
  pin = 2.5_NER
  ! QM # of partial waves
  L = 0
  S = 1
  J = 1
  ! Label of EFT order: 0 -> LO, 1 -> NLO, etc.
  uptoQn = 0

  ! One can call directly chengdu routines to generate mtrx elements
  call chengdu_DLSPR_350(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  ! Or, firstly set the name of the potential,
  call chengdu_dispatch("chengdu_DLSPR_400")
  ! secondly, call the generic routine
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_MMWLY_400(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_MMWLY_600(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_MMWLY_800(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_3200(L, S, J, uptoQn, pout, pin, potval)
  print *, potval


  uptoQn = 0
  L = 0
  S = 0
  J = 0

  call chengdu_dispatch("chengdu_MMWLY_400")
  print *, "chengdu_MMWLY_400, 1S0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_600")
  print *, "chengdu_MMWLY_600, 1S0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_800")
  print *, "chengdu_MMWLY_800, 1S0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_1600")
  print *, "chengdu_MMWLY_1600, 1S0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_3200")
  print *, "chengdu_MMWLY_3200, 1S0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_MMWLY_800(L, S, J, uptoQn, pout, pin, potval)
  print *, potval




  uptoQn = 0
  L = 1
  S = 1
  J = 0

  call chengdu_dispatch("chengdu_MMWLY_400")
  print *, "chengdu_MMWLY_400, 3p0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_600")
  print *, "chengdu_MMWLY_600, 3p0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_800")
  print *, "chengdu_MMWLY_800, 3p0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_1600")
  print *, "chengdu_MMWLY_1600, 3p0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_3200")
  print *, "chengdu_MMWLY_3200, 3p0"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval





  uptoQn = 0
  L = 0
  S = 1
  J = 1

  call chengdu_dispatch("chengdu_MMWLY_400")
  print *, "chengdu_MMWLY_400, 3s1"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_600")
  print *, "chengdu_MMWLY_600, 3s1"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_800")
  print *, "chengdu_MMWLY_800, 3s1"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_1600")
  print *, "chengdu_MMWLY_1600, 3s1"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_dispatch("chengdu_MMWLY_3200")
  print *, "chengdu_MMWLY_3200, 3s1"
  call chengdu_hotpot(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

end program show_mtrx
