# Bingwei Long 08/30/2020

# To compile your code:
# run make mods in nnscat/
# $ make mods
# $ cd src/exec
# $ make -f sample.Makefile get_h2wfs

# Include common.Makefile in root dir
include ../../common.Makefile

# List modules needed by the executable. No longer necessary to write the
# compiling and liking commands. Those commands are included in the pattern rule
# section of common.Makefile
get_h2wfs : eft_tmtrx.o util_gauleg.o

get_cmplx : nneft_lsesv.o nneft_loops.o util_gauleg.o cmplx_epmod.o eft_tmtrx.o
