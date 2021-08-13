# Bingwei Long  08/29/2020
# Bingwei Long  Dec 21/2017

include ../../common.Makefile

EFTOBJS = nneft_loops.o \
          nneft_lsesv.o \
          nneft_eigens.o \
          eft_tmtrx.o \
          eft_phaseconv.o \
          eft_h2wfs.o \
          eft_tmtrx_C.o # This one is C interface

# formfactor.o \
# wave_func_deute.o \
# get_3s1bindenergy.o \
# util_optstruct.o \

.PHONY : eftmods

eftmods : $(EFTOBJS)

eft_h2wfs.o : chengdu.o nneft_type.o util_rootfinder.o nneft_lsesv.o \
  nneft_eigens.o util_gauleg.o

nneft_lsesv.o : nneft_type.o

nneft_eigens.o : nneft_type.o

eft_tmtrx.o : chengdu.o nneft_lsesv.o nneft_loops.o util_gauleg.o \
  eft_potspwd.o eft_h2wfs.o

eft_phaseconv.o : nneft_type.o

nneft_loops.o : nneft_type.o eft_potspwd.o nneft_lsesv.o

wave_func_deute.o : util_gauleg.o nneft_phyconst.o rel_ope_pwd.o vtpe0_pwd.o \
  vtpe1_pwd.o potwrap.o nneft_eigens.o chengdu.o

formfactor.o : wave_func_deute.o

eft_tmtrx_C.o : eft_tmtrx.o chengdu.o nneft_lsesv.o util_gauleg.o eft_potspwd.o

get_3s1bindenergy.o : wave_func_deute.o
