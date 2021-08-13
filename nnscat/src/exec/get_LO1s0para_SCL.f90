! Rui 10/21/2020
! calculate 1s0 LO parameter C0

program get_LO1s0para_SCL

  use mod_obj_smplchn
  use potwrap
  implicit none

!   integer, parameter          :: regtype = REGTYPE_SHARP
  integer, parameter          :: regtype = REGTYPE_GAUSSIAN
  class(obj_smplchn), pointer :: ptr_chn => null()

  logical                     :: succeeded
  real(NER)                   :: xu, xd, C0, Mambda, sclength(1)
  integer                     :: jj, N, chnid, uptoqn

  real(NER)                   :: t0, t1

  call cpu_time(t0)

  allocate(obj_smplchn::ptr_chn)

  ptr_chn%pVfunc_V0 => VLO_withc0_epwrap
!   ptr_chn%N = 200
  N = 100
!   ptr_chn%N = 150
!   ptr_chn%chnstr = '1s0'
  uptoqn = 0
  ptr_chn%NCH = 1
  ptr_chn%npar(1:3) = (/1, 3, 6/)
  ptr_chn%k = 1.0E-12_NER
  ptr_chn%k_defined = .true.

  call convert_text_to_chnid('1s0', chnid, succeeded)
  call create_smplchn(ptr_chn, chnid, uptoQn)
!   call ptr_chn%load_inputs(num_inputs, inputs)

!   Mambda = 1000.0_NER
!   ptr_chn%potpara(1) = -0.33181124678628828E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda = ', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.99_NER
!   xd = ptr_chn%potpara(1)*1.01_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='

!   Mambda = 1300.0_NER
!   ptr_chn%potpara(1) = -0.31103986808504844E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.999_NER
!   xd = ptr_chn%potpara(1)*1.001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='

!   Mambda = 400.0_NER
!   ptr_chn%potpara(1) = -0.46345112093629202E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.999_NER
!   xd = ptr_chn%potpara(1)*1.001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='

!   Mambda = 800.0_NER
!   ptr_chn%potpara(1) = -0.35387941376915466E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.999_NER
!   xd = ptr_chn%potpara(1)*1.001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='

!   Mambda = 4000.0_NER
!   ptr_chn%potpara(1) = -0.26070755680477078E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.999_NER
!   xd = ptr_chn%potpara(1)*1.001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='

!   Mambda = 4800.0_NER
!   ptr_chn%potpara(1) = -0.25624472971879611E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.9999_NER
!   xd = ptr_chn%potpara(1)*1.0001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='

!   Mambda = 6400.0_NER
!   ptr_chn%potpara(1) = -0.25051121960230978E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.9995_NER
!   xd = ptr_chn%potpara(1)*1.0005_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='


!   Mambda = 500.0_NER
!   ptr_chn%potpara(1) = -0.41940507430956253E-002_NER
! !   ptr_chn%potpara(1) = -0.13340507430956253E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.9999_NER
!   xd = ptr_chn%potpara(1)*1.0001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
! !   print *, '====================================='


!   Mambda = 1000.0_NER
!   ptr_chn%potpara(1) = -0.33181124698890589E-002_NER
! !   ptr_chn%potpara(1) = -0.13340507430956253E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.9999_NER
!   xd = ptr_chn%potpara(1)*1.0001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='


!   Mambda = 1400.0_NER
!   ptr_chn%potpara(1) = -0.30600398721280560E-002_NER
! !   ptr_chn%potpara(1) = -0.13340507430956253E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.9999_NER
!   xd = ptr_chn%potpara(1)*1.0001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='


!   Mambda = 1800.0_NER
!   ptr_chn%potpara(1) = -0.29117205777111277E-002_NER
! !   ptr_chn%potpara(1) = -0.13340507430956253E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.9999_NER
!   xd = ptr_chn%potpara(1)*1.0001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
!   print *, '====================================='


!   Mambda = 2200.0_NER
!   ptr_chn%potpara(1) = -0.28144571255542625E-002_NER
! !   ptr_chn%potpara(1) = -0.13340507430956253E-002_NER
!   call ptr_chn%init(N, regtype, Mambda)
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'Mambda =', Mambda
!   print *, 'init_C0 = ', ptr_chn%potpara(1)
!   print *, 'LO_sclength1 = ', sclength(1)
!   xu =  ptr_chn%potpara(1)*0.9999_NER
!   xd = ptr_chn%potpara(1)*1.0001_NER
!   call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
!   print *, "C0 = ", C0
!   ptr_chn%potpara(1) = C0
!   call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
!   print *, 'LO_sclength2 = ', sclength(1)
  print *, '====================================='


  Mambda = 900.0_NER
  ptr_chn%potpara(1) = -0.34088524933945379E-002_NER
  call ptr_chn%init(N, regtype, Mambda)
  call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
  print *, 'Mambda =', Mambda
  print *, 'init_C0 = ', ptr_chn%potpara(1)
  print *, 'LO_sclength1 = ', sclength(1)
  xu =  ptr_chn%potpara(1)*0.99_NER
  xd = ptr_chn%potpara(1)*1.01_NER
  call get_LO1s0para_smplchn(ptr_chn, xu, xd, C0)
  print *, "C0 = ", C0
  ptr_chn%potpara(1) = C0
  call get_1s0sclength_smplchn(ptr_chn, uptoqn, sclength(1))
  print *, 'LO_sclength2 = ', sclength(1)
  print *, '====================================='


  call cpu_time(t1)
  print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
end program get_LO1s0para_SCL