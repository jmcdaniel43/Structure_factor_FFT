#!/bin/bash


module load intel/19.0.5

COMPILER=mpiifort
MKLLIB="-L${MKLROOT}/lib/intel64"
include="../shared_source/"

#OPT="-qopenmp -static"
OPT="-qopenmp -static -check bounds -check uninit -check format -warn declarations -traceback"  #-warn unused 
FFLAGS="-O3 -mkl=sequential $MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl "

# turn this on if you want to profile the code ...
#FFLAGS="-O3 -profile-functions -mkl=sequential $MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl "

echo MKLROOT: $MKLROOT
echo MKLLIB: $MKLLIB
echo FFLAGS: $FFLAGS

export KMP_AFFINITY=verbose,none

$COMPILER $FFLAGS -c ${include}glob_v.f90 ${include}routines.f90 ${include}pme.f90 structure_factor.f90 ${include}main_structure.f90
$COMPILER $FFLAGS *.o $MKLLIB  -o main_structure
