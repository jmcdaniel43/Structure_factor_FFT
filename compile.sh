#!/bin/bash


module load intel/12.1.4
module load openmpi/1.6
module load mkl/10.0


#export PATH=$PATH:/opt/intel/mkl/10.0.3.020/include/
#export LD_LIBRARY_PATH=/opt/intel/composerxe-2011.2.137/mkl/lib/intel64/:$LD_LIBRARY_PATH

OMPINCLUDE=/usr/local/pacerepov1/openmpi/1.6/intel-12.1.4/lib/openmpi/
MKLLIB=/usr/local/pacerepov1/intel/mkl/10.0.5.25/lib/em64t/



#OPT="-openmp -static -check bounds -check uninit -check format -warn declarations -traceback"  #-warn unused 

OPT="-openmp -static"

ifort $OPT -c glob_v.f90 routines.f90 pme.f90 structure_factor.f90 main_structure.f90  -I$MKLLIB -I$OMPINCLUDE -I./mkl_modules
ifort $OPT *.o -L$MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$OMPINCLUDE -lompi_dbg_msgq -o main_structure
