#!/bin/sh

CHEM_MECH=$1

echo $CHEM_MECH
cp -f CHEM_MECH_READY/inSPACK_$CHEM_MECH inSPACK


mkdir $CHEM_MECH
# run SPACK to create fortran routines to compile with ccatt-brams core

rm -f non_zero.dat  *.f90* *.log *.mod *.o
SPACK > spack.log

# run convert program to rewrite spack routines for optimization
cd ../Post-OPT
Make_post
cp post.x ../Application
cd ../Application

post.x > post.log


# rename the files for compilation
mv chem_spack_fexprod.f90   chem_spack_fexprod.f90_orig
mv chem_spack_fexprod_s.f90 chem_spack_fexprod.f90

mv chem_spack_fexloss.f90   chem_spack_fexloss.f90_orig
mv chem_spack_fexloss_s.f90 chem_spack_fexloss.f90

mv chem_spack_jacdchemdc.f90   chem_spack_jacdchemdc.f90_orig
mv chem_spack_jacdchemdc_s.f90 chem_spack_jacdchemdc.f90

mv chem_spack_fexchem.f90   chem_spack_fexchem.f90_orig
mv chem_spack_fexchem_s.f90 chem_spack_fexchem.f90


rm -f *_orig 

mv chem_spack*.f90 chem1*.f90 $CHEM_MECH
cp ../Others_specific_routines/$CHEM_MECH/*.f90  $CHEM_MECH


