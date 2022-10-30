
program show_mtrx

  use nneft_type
  use chengdu
  implicit none

  integer       ::  uptoQn, L, S, J
  real(NER)     ::  pout, pin, potval(1:2, 1:2)

  ! Magnitude of in/out momenta
  pout = 11.0_NER
  pin = 200.0_NER
  ! QM # of partial waves
  L = 0
  S = 1
  J = 1
  ! Label of EFT order: 0 -> LO, 1 -> NLO, etc.
  uptoQn = 0

  call chengdu_DLSPR(L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_generic(1, L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_generic(2, L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_generic(3, L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_generic(4, L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  uptoQn = 1
  L = 0
  S = 0
  J = 0

  call chengdu_DLSPR_generic(1, L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_generic(2, L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_generic(3, L, S, J, uptoQn, pout, pin, potval)
  print *, potval

  call chengdu_DLSPR_generic(4, L, S, J, uptoQn, pout, pin, potval)
  print *, potval


end program show_mtrx
