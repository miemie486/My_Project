# Bingwei Long  Aug 03/2019
# Bingwei Long  Feb 17/2019

FFOPTS = -c $< -o $@
LNKOPTS = -o $@ *.o

# Use the following defs for gfortran
FC = gfortran
GNUOPTS = -ffree-line-length-none -O2
FFOPTS += $(GNUOPTS)
ECOPTS += $(GNUOPTS)

# Use the following defs for ifort
# FC = ifort
# IFORTOPTS = -O2
# FFOPTS += $(IFORTOPTS)
# ECOPTS += $(IFORTOPTS)

.PHONY : all
all : defmods utlmods potsmods show_mtrx

show_mtrx : show_mtrx.o
	$(FC) -o show_mtrx *.o

show_mtrx.o : show_mtrx.f90 chengdu.o nneft_type.o
	$(FC) -c show_mtrx.f90 -o show_mtrx.o -ffree-line-length-none -O2

DEFOBJS = feb_type.o \
          nneft_phyconst.o \
          nneft_type.o

defmods : $(DEFOBJS)

feb_type.o : feb_type.f90
	$(FC) $(FFOPTS)

nneft_type.o : nneft_type.f90 feb_type.o
	$(FC) $(FFOPTS)

nneft_phyconst.o : nneft_phyconst.f90 nneft_type.o
	$(FC) $(FFOPTS)

UTLOBJS = util_gauleg.o \
          util_spcfnc.o \
          util_nr.o \
          util_gadgets.o


utlmods : $(UTLOBJS)

# ****WARNING!****: Must put the .f90 source as the first prerequisite!

util_gauleg.o : util_gauleg.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_gadgets.o : util_gadgets.f90 nneft_type.o nneft_phyconst.o
	$(FC) $(FFOPTS)

util_spcfnc.o : util_spcfnc.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_nr.o : util_nr.f90 nneft_type.o
	$(FC) $(FFOPTS)

POTSOBJS = eft_potspwd.o \
          potwrap.o \
          testwrap.o \
          chengdu.o \
          ope_pwd.o \
          vtpe0_pwd.o \
          vtpe1_pwd.o \
          rel_ope_pwd.o \
          mod_deltafullNLO_pwd.o \
          mod_deltafullN2LO_pwd.o \
          VTPE_SFR_pwd.o \
          VTPE1_SFR_pwd.o \
          delta_epmod.o \
          nopion_epmod.o \
          cmplx_epmod.o

potsmods : $(POTSOBJS)

eft_potspwd.o : eft_potspwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o rel_ope_pwd.o vtpe0_pwd.o vtpe1_pwd.o util_gadgets.o
	$(FC) $(FFOPTS)

testwrap.o : testwrap.f90 potwrap.o
	$(FC) $(FFOPTS)

nopion_epmod.o : nopion_epmod.f90 potwrap.o
	$(FC) $(FFOPTS)

delta_epmod.o : delta_epmod.f90 potwrap.o
	$(FC) $(FFOPTS)

cmplx_epmod.o : cmplx_epmod.f90 eft_potspwd.o ope_pwd.o
	$(FC) $(FFOPTS)

chengdu.o : chengdu.f90 potwrap.o nopion_epmod.o cmplx_epmod.o
	$(FC) $(FFOPTS)

ope_pwd.o : ope_pwd.f90 nneft_type.o nneft_phyconst.o util_spcfnc.o util_gauleg.o
	$(FC) $(FFOPTS)

vtpe0_pwd.o : vtpe0_pwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o rel_ope_pwd.o
	$(FC) $(FFOPTS)

vtpe1_pwd.o : vtpe1_pwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o feb_type.o
	$(FC) $(FFOPTS)

mod_deltafullNLO_pwd.o : mod_deltafullNLO_pwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o util_gauleg.o
	$(FC) $(FFOPTS)

mod_deltafullN2LO_pwd.o : mod_deltafullN2LO_pwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o util_gauleg.o mod_deltafullNLO_pwd.o vtpe1_pwd.o
	$(FC) $(FFOPTS)

VTPE_SFR_pwd.o : VTPE_SFR_pwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o util_gauleg.o
	$(FC) $(FFOPTS)

VTPE1_SFR_pwd.o : VTPE1_SFR_pwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o feb_type.o
	$(FC) $(FFOPTS)

rel_ope_pwd.o : rel_ope_pwd.f90 nneft_type.o nneft_phyconst.o ope_pwd.o
	$(FC) $(FFOPTS)

potwrap.o : potwrap.f90 eft_potspwd.o ope_pwd.o vtpe0_pwd.o vtpe1_pwd.o mod_deltafullNLO_pwd.o mod_deltafullN2LO_pwd.o VTPE_SFR_pwd.o VTPE1_SFR_pwd.o
	$(FC) $(FFOPTS)

.PHONY : clean
clean :
	rm -f ./*.mod
	rm -f ./*.o

