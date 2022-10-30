# Bingwei Long  08/29/2020
# Bingwei Long  Dec 20/2016

# Choose compiler according to local env $EFTFCC; if undefined, use gfortran
ifdef TEGFC
FC = $(TEGFC)
else
FC = gfortran
endif

# Defining directories
DMOD=../../mod
DOBJ=../../obj
DBIN=../../bin
DDEF=../def
DCHN=../channel
DSMP=../smpl
DEFT=../eft
DDRV=../drv
DEXE=../exec
DUTL=../util
DFEX=../fitexec

#VPATH is a special variable that tells "make" which path to look for prerequisites
VPATH = $(DBIN):$(DOBJ):$(DDEF):$(DCHN):$(DSMP):$(DEFT):$(DDRV):$(DEXE):$(DUTL):$(DFEX)

#======= Compile flags/options =======#

# Compiler options used to build *modules* and shared by both gfortran and ifort
# -o $(DOBJ)/$@ says to put .o file in $(DOBJ). $@ is an automatic make variable
# that refers to the file name in a rule
# $< is the name of the first prerequisite
#
# ****WARNING!****: Must put the .f90 source as the first prerequisite!
#
# FFOPTS = -c -o $(DOBJ)/$@ $<
# FFOPTS = -c -O2 -o $(DOBJ)/$@ $<
FFOPTS = -c -O3 -o $(DOBJ)/$@ $<

# Compiler options used to build *executables* and shared by both gfortran and ifort
# ECOPTS = -c -o $(DBIN)/$@ $<
# ECOPTS = -c -O2 -o $(DBIN)/$@ $<
ECOPTS = -O3 -o $(DBIN)/$@ $<

# Compiler options for gfortran
# -I$(DMOD) suggests a path to look for the .mod if it is absent in the PWD
# -J$(DMOD) says where to put the .mod when compiling a module
# GNUOPTS = -I$(DMOD) -J$(DMOD) -ffree-line-length-none -fopenmp
GNUOPTS = -I$(DMOD) -J$(DMOD) -ffree-line-length-none
ifdef TEG_NOOMP_FLAG
else
GNUOPTS += -fopenmp
endif

# Compiler options for ifort
# -I $(DMOD) suggests a path to look for the .mod if it is absent in the PWD
# -module $(DMOD) says where to put the .mod when compiling a module
# IFORTOPTS= -I $(DMOD) -module $(DMOD) -qopenmp
IFORTOPTS= -I $(DMOD) -module $(DMOD)
ifdef TEG_NOOMP_FLAG
else
IFORTOPTS += -qopenmp
endif

ifeq ($(FC), ifort)
FFOPTS += $(IFORTOPTS)
ECOPTS += $(IFORTOPTS)
else
FFOPTS += $(GNUOPTS)
ECOPTS += $(GNUOPTS)
endif

#======= Link options =======#

LNKOPTS = $(DOBJ)/*.o

ifeq ($(FC), ifort)
LIBS += -mkl -lminuitif
# LIBS += -mkl
else
LIBS += -llapack -lblas -lminuitgf
# LIBS += -L/usr/lib -llapack -lblas -lminuitgf
# LIBS += -llapack -lblas
endif
ifeq ($(FC), ifort)
LNKOPTS += -qopenmp
else
LNKOPTS += -fopenmp
endif

#======= Pattern rules =======#

# Regular modules
%.o: %.f90
	$(FC) $(FFOPTS)

# Executables
% : %.f90
	$(FC) $(ECOPTS) $(LNKOPTS) $(LIBS)

%.exe : %.f90
	$(FC) $(ECOPTS) $(LNKOPTS) $(LIBS)
