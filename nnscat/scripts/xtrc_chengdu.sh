# DEFf90="feb_type.f90 nneft_phyconst.f90 nneft_type.f90"

AUXDIR=../../chengdu/

rm chengdu/*.Makefile
rm chengdu/*.f90
cp chengdu.Makefile chengdu/
cd src/def/
cp feb_type.f90 nneft_phyconst.f90 nneft_type.f90 $AUXDIR
cd ../util/
cp util_gauleg.f90 util_spcfnc.f90 util_nr.f90 util_gadgets.f90 $AUXDIR
cd ../pots/
cp eft_potspwd.f90 potwrap.f90 testwrap.f90 chengdu.f90 ope_pwd.f90 vtpe0_pwd.f90 \
vtpe1_pwd.f90 rel_ope_pwd.f90 mod_deltafullNLO_pwd.f90 \
mod_deltafullN2LO_pwd.f90 VTPE_SFR_pwd.f90 VTPE1_SFR_pwd.f90 \
*_epmod.f90 $AUXDIR
cd ../exec/
cp show_mtrx.f90 $AUXDIR
cd ../../doc
cp chengdu.README ../chengdu/
cd ../chengdu/
tar -cvzf chengdupots.tar.gz chengdu.Makefile *.f90

# tar -cvf chengdupots.tar chengdu.Makefile
# cd src/def/
# tar -rvf ../../chengdupots.tar feb_type.f90 nneft_phyconst.f90 nneft_type.f90
# cd ../util/
# tar -rvf ../../chengdupots.tar util_gauleg.f90 util_spcfnc.f90 util_nr.f90 util_gadgets.f90
# cd ../pots/
# tar -rvf ../../chengdupots.tar eft_potspwd.f90 potwrap.f90 chengdu.f90 ope_pwd.f90 vtpe0_pwd.f90 vtpe1_pwd.f90 rel_ope_pwd.f90 mod_deltafullNLO_pwd.f90 mod_deltafullN2LO_pwd.f90
# cd ../exec/
# tar -rvf ../../chengdupots.tar show_mtrx.f90
# cd ../../doc
# tar -rvf ../chengdupots.tar chengdu.README
