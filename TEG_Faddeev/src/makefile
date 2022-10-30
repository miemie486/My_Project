# "-" at the beginning of each line asks make to ignore error message and not to stop
# "@" asks make not to print the commands that are being executed.

MK_NNPATH=../../nnscat
MK_UTIL=util
MK_TMAT=tMat
MK_KMAT=kMat
MK_EXEC=exec
MK_BIN=../bin
MK_OBJ=../obj

.PHONY : nothing all objs allTest clean cobjbin cbin cobj cleanTest

nothing :
	@-echo "No default make option."
	@-echo "Try make all, make objs, make allTest, etc."

execs = scatlen eigenH3 elasui elasuidelta

$(execs) :
	@-make objs
	@-cd $(MK_EXEC);make $@.exe

all :
	@-make objs
	@-cd $(MK_EXEC);make allbin --keep-going
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""

objs :
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo "****////   NOTE: Make objs.    ////****"
	@-export TEG_OMP_FLAG=1
	@-cd $(MK_NNPATH); make eft --keep-going
	@-cd $(MK_UTIL); make allobj --keep-going
	@-cd $(MK_TMAT); make allobj --keep-going
	@-cd $(MK_KMAT); make allobj --keep-going
	@-echo ""

clean:
	@-echo "****////   Cleaning  ////****"
	@-rm -f $(MK_BIN)/* $(MK_OBJ)/*
	@-cd $(MK_NNPATH); make clean --keep-going
	@-echo "****////   Cleaning Done  ////****"

cobjbin:
	@-echo "****////   Cleaning  ////****"
	@-rm -f $(MK_BIN)/* $(MK_OBJ)/*

cbin:
	@-echo "****////   Cleaning  ////****"
	@-rm -f $(MK_BIN)/*

cobj:
	@-echo "****////   Cleaning  ////****"
	@-rm -f $(MK_OBJ)/*

allTest:
	@-echo ""
	@-echo "****////   NOTE: Run all tests   ////****"
	@-echo ""
	@-cd ./test/testTemp; make allTest
	@-echo ""
	@-echo "****////   NOTE: Run all tests   ////****"

	@-echo ""


cleanTest:
	@-cd ./test/testTemp; make clean
