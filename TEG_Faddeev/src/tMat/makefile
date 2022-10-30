include ../makefile.common
NNOBJS = $(NNPATH)/obj/*.o

# ****************************************************
# Targets

# allobj: tMat.o tMatGen.o
allobj: tMat.o

# tMatGen.o: gaussIntegral.h tMatGen.h def.h
tMat.o: gaussIntegral.h tMat.h def.h msg.h msg.o $(NNOBJS)
