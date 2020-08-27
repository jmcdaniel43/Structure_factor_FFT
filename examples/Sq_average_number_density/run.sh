#!/bin/bash

# calculate scattering structure factor with PME approach...
../../Sq_Average_FFT_source/main_structure nvt_300K.gro traj.gro ../../atomic_form_factors/ calculation_parameters.pmt > Sq_number.out

# now use a little perl script to average the scattering and plot isotropic scattering as a function of wavevector magnitude
../../perl_scripts/average_Sq.pl Sq_number.out > Sq_mag
