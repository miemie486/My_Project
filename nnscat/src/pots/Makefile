# Bingwei Long  08/30/2020
# Bingwei Long  Feb 17/2019

include ../../common.Makefile

POTSOBJS = eft_potspwd.o \
          potwrap.o \
          chengdu.o \
          ope_pwd.o \
          vtpe0_pwd.o \
          vtpe1_pwd.o \
          mod_toypot.o \
          mod_deltafullNLO_pwd.o \
          mod_deltafullN2LO_pwd.o \
          rel_ope_pwd.o \
          VTPE_SFR_pwd.o \
          VTPE1_SFR_pwd.o \
          delta_epmod.o \
          nopion_epmod.o \
          cmplx_epmod.o \
          testwrap.o \
          vbar_wfs.o \
          vbar_3p0_wfs.o

potsmods : $(POTSOBJS)

eft_potspwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o rel_ope_pwd.o \
  vtpe0_pwd.o vtpe1_pwd.o util_gadgets.o
potwrap.o : eft_potspwd.o ope_pwd.o vtpe0_pwd.o vtpe1_pwd.o \
  mod_deltafullNLO_pwd.o mod_deltafullN2LO_pwd.o VTPE_SFR_pwd.o VTPE1_SFR_pwd.o
cmplx_epmod.o : eft_potspwd.o ope_pwd.o
delta_epmod.o : potwrap.o
chengdu.o : potwrap.o nopion_epmod.o cmplx_epmod.o testwrap.o vbar_wfs.o vbar_3p0_wfs.o
nopion_epmod.o : potwrap.o
ope_pwd.o : nneft_type.o nneft_phyconst.o util_spcfnc.o util_gauleg.o
vtpe0_pwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o rel_ope_pwd.o
vtpe1_pwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o feb_type.o
mod_toypot.o : ope_pwd.o util_gauleg.o feb_type.o
mod_deltafullNLO_pwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o util_gauleg.o
mod_deltafullN2LO_pwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o util_gauleg.o mod_deltafullNLO_pwd.o
rel_ope_pwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o
VTPE_SFR_pwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o util_gauleg.o
VTPE1_SFR_pwd.o : nneft_type.o nneft_phyconst.o ope_pwd.o feb_type.o
testwrap.o : potwrap.o
vbar_wfs.o : nneft_type.o util_rootfinder.o nneft_lsesv.o nneft_eigens.o util_gauleg.o
vbar_3p0_wfs.o : nneft_type.o util_rootfinder.o nneft_lsesv.o nneft_eigens.o util_gauleg.o

