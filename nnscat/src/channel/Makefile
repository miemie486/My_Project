# Bingwei Long,  Dec 20/2017

include ../../common.Makefile

CHNOBJS =  \
	mod_obj_channel.o \
	mod_obj_smplchn.o \
	mod_obj_sngl.o \
	mod_obj_1s0LY.o \
	mod_obj_withc0_sngl.o \
	mod_obj_lbwdib.o \
	mod_obj_pope_sngl.o \
	mod_obj_pope_cpld.o \
	mod_obj_suppope_cpld.o \
	mod_obj_suppope_sngl.o \
	mod_obj_cpld.o \
	mod_obj_1s0SYLV.o \
	mod_vfunc.o \
	mod_obj_sprbl1s0.o

# mod_vfunc_rel_ope.o
#	mod_obj_withc0_residual_sngl.o \
#	mod_obj_1s0_mpistar_sngl.o \
# mod_obj_dibcpld.o \
# mod_obj_tlrLY.o \
# mod_obj_woc0_taylor_sngl.o \
# mod_obj_withc0_cpld.o \
# mod_obj_taylor_sngl.o \
# mod_obj_taylor_cpld.o \
# mod_obj_altpope_cpld.o \
# mod_obj_pope_vnfsngl.o \

chnmods : $(CHNOBJS)

# ****WARNING!****: Must put the .f90 source as the first prerequisite!

mod_obj_smplchn.o : mod_obj_smplchn.f90 nneft_type.o eft_potspwd.o util_gauleg.o nneft_loops.o eft_phaseconv.o util_io.o util_cheb.o
	$(FC) $(FFOPTS)

mod_obj_sprbl1s0.o : mod_obj_sprbl1s0.f90 mod_obj_lbwdib.o
	$(FC) $(FFOPTS)

mod_obj_suppope_cpld.o : mod_obj_suppope_cpld.f90 mod_obj_pope_cpld.o
	$(FC) $(FFOPTS)

mod_obj_suppope_sngl.o : mod_obj_suppope_sngl.f90 eft_phaseconv.o mod_obj_sngl.o
	$(FC) $(FFOPTS)

# mod_obj_dibcpld.o : mod_obj_dibcpld.f90 mod_obj_cpld.o
# 	$(FC) $(FFOPTS)

mod_obj_1s0SYLV.o : mod_obj_1s0SYLV.f90 mod_obj_1s0LY.o
	$(FC) $(FFOPTS)

mod_obj_pope_cpld.o : mod_obj_pope_cpld.f90 mod_obj_cpld.o
	$(FC) $(FFOPTS)

mod_obj_cpld.o : mod_obj_cpld.f90 mod_obj_channel.o nneft_lsesv.o nneft_loops.o
	$(FC) $(FFOPTS)

mod_obj_pope_sngl.o : mod_obj_pope_sngl.f90 eft_phaseconv.o mod_obj_sngl.o
	$(FC) $(FFOPTS)

# mod_obj_withc0_residual_sngl.o : mod_obj_withc0_residual_sngl.f90 mod_obj_withc0_sngl.o
# 	$(FC) $(FFOPTS)

mod_obj_lbwdib.o : mod_obj_lbwdib.f90 mod_obj_withc0_sngl.o util_gadgets.o eft_potspwd.o
	$(FC) $(FFOPTS)

mod_obj_1s0LY.o : mod_obj_1s0LY.f90 mod_obj_withc0_sngl.o
	$(FC) $(FFOPTS)

# mod_obj_tlrLY.o : mod_obj_tlrLY.f90 mod_obj_withc0_sngl.o
# 	$(FC) $(FFOPTS)

# mod_obj_1s0_mpistar_sngl.o : mod_obj_1s0_mpistar_sngl.f90 mod_obj_1s0LY.o
# 	$(FC) $(FFOPTS)

mod_obj_withc0_sngl.o : mod_obj_withc0_sngl.f90 mod_obj_sngl.o
	$(FC) $(FFOPTS)

mod_obj_sngl.o : mod_obj_sngl.f90 mod_obj_channel.o nneft_loops.o nneft_lsesv.o mod_toypot.o nneft_eigens.o
	$(FC) $(FFOPTS)

mod_obj_channel.o : mod_obj_channel.f90 nneft_type.o eft_potspwd.o util_gauleg.o nneft_loops.o eft_phaseconv.o util_io.o util_cheb.o nneft_phyconst.o
	$(FC) $(FFOPTS)

mod_vfunc.o : mod_vfunc.f90 mod_deltafullNLO_pwd.o mod_obj_sngl.o mod_obj_cpld.o mod_deltafullN2LO_pwd.o VTPE_SFR_pwd.o VTPE1_SFR_pwd.o
	$(FC) $(FFOPTS)

# mod_vfunc_rel_ope.o : mod_vfunc_rel_ope.f90 rel_ope_pwd.o mod_obj_sngl.o mod_obj_cpld.o
# 	$(FC) $(FFOPTS)

# mod_obj_pope_vnfsngl.o : mod_obj_pope_vnfsngl.f90 mod_obj_pope_sngl.o mod_vfunc.o
	# $(FC) $(FFOPTS)

# mod_obj_altpope_cpld.o : mod_obj_altpope_cpld.f90 mod_obj_pope_cpld.o
# 	$(FC) $(FFOPTS)

# mod_obj_taylor_cpld.o : mod_obj_taylor_cpld.f90 mod_obj_cpld.o util_cheb.o
#   $(FC) $(FFOPTS)

# mod_obj_woc0_taylor_sngl.o : mod_obj_woc0_taylor_sngl.f90 mod_obj_taylor_sngl.o
#   $(FC) $(FFOPTS)

# mod_obj_taylor_sngl.o : mod_obj_taylor_sngl.f90 mod_obj_sngl.o util_cheb.o
#   $(FC) $(FFOPTS)


