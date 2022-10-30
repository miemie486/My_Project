include ../../common.Makefile


KEYEXETGTS = \
  get_phase_withCD_smpl  \
  get_mtrx_data \
  get_h2wfs_realT \
  get_1s0_scatlength_smpl \
  get_LO3s1para_BD \
  get_pert_bindE \
  get_mmwly_smpl \
  get_pert_SCL \
  get_LO1s0para_SCL \
  get_LO1s0_SCL \
  get_bdscl_withCD_smpl \
  get_3p0_wfs

.PHONY : keyexes

keyexes : $(KEYEXETGTS)

get_phase_withCD_smpl : drv_pwphase.o mod_vfunc_smpl.o nneft_phyconst.o

get_mtrx_data : nneft_type.o chengdu.o util_gauleg.o

get_h2wfs_realT : eft_tmtrx.o util_gauleg.o util_io.o mod_obj_smplchn.o

get_1s0_scatlength_smpl : drv_fitpwa.o mod_obj_smplchn.o

get_LO3s1para_BD : eft_h2wfs.o nneft_type.o

get_pert_bindE : mod_obj_smplchn.o util_io.o potwrap.o

get_mmwly_smpl : drv_pwphase.o mod_vfunc_smpl.o

get_pert_SCL : mod_obj_smplchn.o util_io.o potwrap.o

get_LO1s0para_SCL : mod_obj_smplchn.o potwrap.o

get_LO1s0_SCL : mod_obj_smplchn.o potwrap.o util_io.o

get_bdscl_withCD_smpl : drv_pwphase.o mod_vfunc_smpl.o

get_3p0_wfs : vbar_3p0_wfs.o util_gauleg.o util_io.o mod_obj_smplchn.o
