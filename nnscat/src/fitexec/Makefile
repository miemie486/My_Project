# Bingwei Long 08/30/2020
# Bingwei Long  07/27/2018

include ../../common.Makefile

KEYTGTS = \
  fit_lbwdib \
  fit_3s1 \
  fit_3p0 \
  fitspr_smpl

.PHONY : keyexes

keyexes : $(KEYTGTS)

fit_lbwdib : fit_lbwdib.f90 drv_fitting.o util_datetime.o

fitspr_smpl : fitspr_smpl.f90 drv_fitting.o

fit_3s1 : fit_3s1.f90 drv_fitting.o util_datetime.o drv_pwphase.o \
	mod_obj_smplchn.o mod_vfunc_smpl.o

fit_3p0 : fit_3p0.f90 drv_fitting.o util_datetime.o
