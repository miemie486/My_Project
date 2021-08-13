# Bingwei Long,  Dec 22/2017

include ../../common.Makefile

DRVOBJS = drv_pwphase.o \
          drv_fitpwa.o \
          drv_fitting.o \
          p_wave_phaseshift.o


# Frozen
# drv_solve.o \

drvmods : $(DRVOBJS)

# ****WARNING!****: Must put the .f90 source as the first prerequisite!


drv_fitpwa.o : drv_fitpwa.f90 mod_obj_1s0SYLV.o eft_potspwd.o eft_phaseconv.o util_io.o mod_obj_withc0_sngl.o mod_obj_1s0LY.o mod_obj_cpld.o mod_obj_pope_cpld.o mod_obj_lbwdib.o mod_obj_sprbl1s0.o mod_obj_suppope_sngl.o mod_obj_suppope_cpld.o drv_pwphase.o mod_obj_smplchn.o mod_vfunc_smpl.o mod_vfunc.o mod_obj_pope_sngl.o
	$(FC) $(FFOPTS)

# Frozen
# drv_solve.o : drv_solve.f90 util_rootfinder.o drv_pwphase.o util_io.o
# 	$(FC) $(FFOPTS)

drv_pwphase.o : drv_pwphase.f90 mod_obj_smplchn.o mod_obj_sngl.o mod_obj_cpld.o mod_obj_withc0_sngl.o mod_obj_1s0LY.o mod_obj_lbwdib.o mod_obj_sprbl1s0.o mod_obj_pope_sngl.o mod_obj_pope_cpld.o mod_obj_1s0SYLV.o eft_potspwd.o util_io.o util_datetime.o mod_obj_suppope_sngl.o mod_obj_suppope_cpld.o vtpe1_pwd.o mod_deltafullN2LO_pwd.o mod_vfunc_smpl.o
	$(FC) $(FFOPTS)

drv_fitting.o : drv_fitting.f90 drv_fitpwa.o util_gadgets.o util_datetime.o mod_deltafullN2LO_pwd.o mod_obj_smplchn.o
	$(FC) $(FFOPTS)

p_wave_phaseshift.o : p_wave_phaseshift.f90 mod_obj_suppope_sngl.o mod_obj_suppope_cpld.o mod_vfunc.o
	$(FC) $(FFOPTS)

