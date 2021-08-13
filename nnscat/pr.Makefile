.PHONY : noting eft mods clean keys all

# "-" at the beginning of each line asks make to ignore error message and not to stop
# "@" asks make not to print the commands that are being executed.
nothing :
	@-cd src/exec/; make -f pr.Makefile --keep-going
	@-cd src/fitexec/; make -f pr.Makefile --keep-going
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""

all :
	@-make mods
	@-cd src/exec/; make -f pr.Makefile --keep-going
	@-cd src/fitexec/; make -f pr.Makefile --keep-going
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""

