Package that generates nucleon-nucleon scattering phase shifts based on a renormalization-group invariant power counting. The code is written in Fortran 2003. The object-oriented facility "class" is utilized so that the code can be easily modified to accommodate new modifications made to the power counting.

# External libraries
- LAPACK (normally to be installed by the package manager of given Linux distro)
- MINUIT (to be installed by the source tar ball located in extlibs/)


# make existing modules
- To compile modules in the nnscat package, navigate to the root dir of nnscat
$ cd nnscat/
$ make mods

- To make essential executables
$ make exec

- To make executables listed in src/exec/misc.Makefile and src/fitexec/misc.Makefile
$ make all

- To remove obj and binaries
$ make clean


# Compile your own executables
Suppose your own executable source is "AProg.f90", which uses nneft_type and
nneft_lsesv.
1. It has to be in a directory one level beneath nnscat/src/, e.g., nnscat/src/exec
2. Write a Makefile of your own(!!!), e.g., my.Makefile

Inside nnscat/src/exec/my.Makefile
-------------------

include ../../common.Makefile

AProg.exe : nneft_type.o nneft_lsesv.o

-------------------
Note that you do not need to add compiling commands; they are already included
as pattern rules in common.Makefile. There is an example Makefile:
src/exec/sample.Makefile

To compile,

$ cd nnscat/
$ make mods; cd src/exec; make -f my.Makefile AProg.exe; cd -


# Compile your own modules
Suppose your module is in "AMod.f90", which uses nneft_type and nneft_lsesv.
1. It has to be one level under nnscat/src/, e.g., in nnscat/src/eft/

2. The best practice is to write a simple Makefile of your own, e.g.,
my.Makefile. This way, your module can be kept away from main packages before
it is fully accepted.

Inside my.Makefile:
-------------------

include ../../common.Makefile

AMod.o : nneft_type.o nneft_lsesv.o

-------------------
Note that you do not need to add compiling commands; they are already included
as pattern rules in common.Makefile.

To compile

$ cd nnscat/
$ make mods; cd src/eft; make -f my.Makefile AMod.o;cd -


# Env vars
- To invoke openMP, set TEG_OMP_FLAG to any value
$ export TEG_OMP_FLAG=1


# Misc

- Compilation with ifort has not been well tested
