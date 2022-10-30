# Bingwei Long  08/29/2020
# Bingwei Long  Dec 22/2017

include ../../common.Makefile

KEYEXETGTS = \
  get_low.exe    \
  get_supperi.exe \
  get_peri.exe \
  get_NP_smpl.exe \
  show_mtrx.exe

.PHONY : keyexes

keyexes : $(KEYEXETGTS)

get_low.exe : drv_pwphase.o

get_NP_smpl.exe	: drv_pwphase.o mod_vfunc_smpl.o

get_peri.exe : drv_pwphase.o

get_supperi.exe : drv_pwphase.o

show_mtrx.exe :


