
program get_LO3s1para_BD
  use eft_h2wfs
  use nneft_type

  real(NER)   :: x1 , x2 , tol
  real(NER)             :: mambda, rslt
  integer               :: Nmesh

  Nmesh = 100
  tol =1.0E-9_NER

!   mambda = 350.0_NER
!   x1 = -5.2E-003_NER
!   x2 = -5.6E-003_NER

!   print*, '--------------------------'
!   print*, 'Nmesh =', Nmesh
!   print*, '--------------------------'
!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 400.0_NER
!   x1 = -4.1E-003_NER
!   x2 = -4.6E-003_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 600.0_NER
!   x1 = -1.1E-003_NER
!   x2 = -1.6E-003_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 800.0_NER
!   x1 = 2.8E-003_NER
!   x2 = 3.1E-003_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 1200.0_NER
!   x1 = -1.2E-002_NER
!   x2 = -1.3E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 1600.0_NER
!   x1 = -4.6E-003_NER
!   x2 = -4.7E-003_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 2000.0_NER
!   x1 = -2.1E-003_NER
!   x2 = -2.3E-003_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 2400.0_NER
!   x1 = 6.1E-004_NER
!   x2 = 6.4E-004_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 2800.0_NER
!   x1 = 1.2E-002_NER
!   x2 = 1.4E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 3200.0_NER
!   x1 = -1.2E-002_NER
!   x2 = -1.3E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 4000.0_NER
!   x1 = -0.3E-002_NER
!   x2 = -0.45E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 4800.0_NER
!   x1 = -0.1E-002_NER
!   x2 = -0.2E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 6400.0_NER
!   x1 = -0.1E-001_NER
!   x2 = -0.3E-001_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 1000.0_NER
!   x1 = 0.4E-001_NER
!   x2 = 0.5E-001_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 1400.0_NER
!   x1 = -0.6E-002_NER
!   x2 = -0.7E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 1800.0_NER
!   x1 = -0.2E-002_NER
!   x2 = -0.4E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

!   mambda = 2200.0_NER
!   x1 = -0.7E-003_NER
!   x2 = -0.1E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

  mambda = 900.0_NER
  x1 = 0.80E-002_NER
  x2 = 1.00E-001_NER

  call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
  print*,'mambda = ', mambda
  print *, '3s1_LOpara(1) = ', rslt

!   mambda = 500.0_NER
!   x1 = -0.2E-002_NER
!   x2 = -0.3E-002_NER

!   call find_LOpara(mambda, Nmesh, x1, x2, tol, rslt)
!   print*,'mambda = ', mambda
!   print *, '3s1_LOpara(1) = ', rslt

end program get_LO3s1para_BD