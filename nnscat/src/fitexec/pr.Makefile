include ../../common.Makefile


KEYEXETGTS = \
  fit_3s1bd_smpl\
  fit_mmwly_smpl\
  fit_ksw_smpl \
  fit_1s0scl_smpl


.PHONY : keyexes

keyexes : $(KEYEXETGTS)

fit_3s1bd_smpl : drv_fitting.o

fit_mmwly_smpl : drv_fitting.o

fit_ksw_smpl : drv_fitting.o

fit_1s0scl_smpl : drv_fitting.o