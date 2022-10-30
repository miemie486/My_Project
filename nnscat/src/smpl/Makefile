# Bingwei Long,  10/20/2018

include ../../common.Makefile

SMPLOBJS =  \
	mod_obj_smplchn.o \
	mod_vfunc_smpl.o \

#	mod_obj_zeronlo.o

smplmods : $(SMPLOBJS)

# ****WARNING!****: Must put the .f90 source as the first prerequisite!

mod_obj_smplchn.o : mod_obj_smplchn.f90 nneft_type.o eft_potspwd.o util_gauleg.o nneft_lsesv.o eft_phaseconv.o util_io.o util_cheb.o potwrap.o nneft_phyconst.o eft_h2wfs.o util_rootfinder.o testwrap.o
	$(FC) $(FFOPTS)

mod_vfunc_smpl.o : mod_vfunc_smpl.f90 mod_obj_smplchn.o potwrap.o delta_epmod.o testwrap.o chengdu.o
	$(FC) $(FFOPTS)
# mod_obj_zeronlo.o : mod_obj_zeronlo.f90 mod_obj_smplchn.o
# 	$(FC) $(FFOPTS)
