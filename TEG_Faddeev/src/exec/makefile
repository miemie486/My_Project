include ../makefile.common

NNOBJS = $(NNPATH)/obj/*.o

# ****************************************************
# Targets
.PHONY : allbin

# allbin: scatlen.exe eigenH3.exe 
allbin: scatlen.exe eigenH3.exe showdetk.exe elasui.exe elasuidelta.exe

eigenH3.exe : kMatGen.h kmathmgns.tpp kMatGen.o channels.h channels.o fstr.h fstr.o msg.h msg.o $(NNOBJS)

showdetk.exe : kMatGen.h kMatGen.o channels.h channels.o fstr.h fstr.o msg.h msg.o $(NNOBJS)

scatlen.exe : kMatGen.h kmathmgns.tpp kMatGen.o uMat.h uMat.o msg.h msg.o $(NNOBJS)

elasui.exe : kMatGen.h kmathmgns.tpp kMatGen.o msg.h msg.o $(NNOBJS)

elasuidelta.exe : kMatGen.h kmathmgns.tpp kMatGen.o msg.h msg.o delta.h $(NNOBJS)

# ****************************************************
# Obsolete
# ****************************************************

listalpha: $(OBJS_eigenH3) $(NNOBJS)
	$(CXX) -o $(BINDIR)/$@.exe  listalpha.cpp $^ $(LIB) $(INC)

sublH3: $(OBJS_eigenH3) $(NNOBJS)
	$(CXX) -o $(BINDIR)/$@.exe  sublH3.cpp $^ $(LIB) $(INC)

testFile: $(OBJS_scat) $(NNOBJS)
	$(CXX) -o $(BINDIR)/$@  testFile.cpp $^ $(LIB) $(INC)

