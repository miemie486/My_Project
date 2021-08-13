# Bingwei Long,  Dec 22/2017

include ../Makefile.common

OBJS =   mod_obsv.o \
			util_mathrle.o \
			mod_fit.o
# Frozen
# drv_solve.o \

mods : $(OBJS)

# ****WARNING!****: Must put the .f90 source as the first prerequisite!

mod_obsv.o : mod_obsv.f90 util_mathrle.o nneft_type.o 
	$(FC) $(FFOPTS)

# Frozen
# drv_solve.o : drv_solve.f90 util_rootfinder.o drv_pwphase.o util_io.o
# 	$(FC) $(FFOPTS)

util_mathrle.o : util_mathrle.f90 nneft_type.o 
	$(FC) $(FFOPTS)

mod_fit.o : mod_fit.f90 mod_obsv.o drv_fitpwa.o
	$(FC) $(FFOPTS)
