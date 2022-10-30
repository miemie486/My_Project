# Bingwei Long  08/30/2020
# Bingwei Long  Dec 21/2017

include ../../common.Makefile

UTLOBJS = util_cheb.o \
          util_gauleg.o \
          util_gadgets.o \
          util_inter_extra.o \
          util_io.o \
          util_nr.o \
          util_rootfinder.o \
          util_spcfnc.o \
          util_datetime.o\
          util_cbcspn.o

# util_mesh.o \
# util_bscint.o \

utlmods : $(UTLOBJS)

util_gauleg.o : nneft_type.o
util_inter_extra.o : nneft_type.o
# util_mesh.o : util_mesh.f90 nneft_type.o util_gauleg.o
# util_bscint.o : util_bscint.f90 nneft_type.o
util_gadgets.o : nneft_type.o nneft_phyconst.o
util_io.o : nneft_type.o
util_spcfnc.o : nneft_type.o
util_rootfinder.o : nneft_type.o
util_nr.o : nneft_type.o
util_cbcspn.o : nneft_type.o
util_cheb.o	: nneft_type.o util_nr.o
util_datetime.o :
