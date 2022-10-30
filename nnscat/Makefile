# Bingwei Long  08/29/2020
# Bingwei Long  Dec 21/2017
# Shaowei Wu Jan 17/2018

.PHONY : noting eft mods clean keys all

# "-" at the beginning of each line asks make to ignore error message and not to stop
# "@" asks make not to print the commands that are being executed.
nothing :
	@-echo "No default make option."
	@-echo "Try make eft, make mods, etc."

eft :
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""
	@-cd src/def/; make --keep-going
	@-cd src/util/; make --keep-going
	@-cd src/pots/; make --keep-going
	@-cd src/eft/; make --keep-going
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""

mods :
	@-make eft
	@-cd src/channel/; make --keep-going
	@-cd src/smpl/; make --keep-going
	@-cd src/drv/; make --keep-going

exec :
	@-make mods
	@-cd src/exec/; make --keep-going
	@-cd src/fitexec/; make --keep-going
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""

all :
	@-make exec
	@-cd src/exec/; make -f misc.Makefile --keep-going
	@-cd src/fitexec/; make -f misc.Makefile --keep-going
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""

clean :
	rm -f mod/*.mod
	rm -f obj/*.o
	rm -f bin/*
