#!/bin/sh
#
# Runs a benchmark suite consisting of 12 of the longer testscripts
#
./dalton-benchmark -f -c options.example -t cc_geoopt2,cc_rsp_twophotdirect,dft_optimize_nosym,dft_properties_sym,energy_solv,geoopt_prop3,prop_roa,prop_vibana,rsp_excipolar,rsp_mnfphosph,walk_polar2,walk_solvmag
