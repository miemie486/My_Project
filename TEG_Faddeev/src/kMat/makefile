include ../makefile.common

# allobj: spline.o matrixK.o 3nj.o gAndG.o gaussIntegral.o channels.o kMatGen.o
allobj: fstr.o gAndG.o channels.o kMatGen.o uMat.o

fstr.o: fstr.h csvReader.h csvReader.o
gAndG.o: gAndG.h 3nj.h def.h
channels.o: channels.h

kMatGen.o : kMatGen.cpp kMatGen.h kmatinh.tpp kmathmgns.tpp def.h gaussIntegral.h  gaussIntegral.o \
	spline.h spline.o gAndG.h gAndG.o tMat.h tMat.o channels.h channels.o \
	matOperator.h matOperator.o msg.o msg.h \

uMat.o : uMat.cpp uMat.h kMatGen.h kMatGen.o fstr.h fstr.o msg.h msg.o
