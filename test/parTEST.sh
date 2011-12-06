#!/bin/bash
#
# -------------------------------------------------------------------------------
# Run all DALTON2011 parallel tests using MPI
#
# This is a sample script for running the the parallel tests with PBS
# assuming at least four CPU cores per node.
# Modify to fit your own setup.
#
# -- Hans Jørgen Aa. Jensen, Dec. 2010
# -------------------------------------------------------------------------------
#
# Ask PBS for 2 nodes.
#
#PBS -l nodes=2
#PBS -l walltime=00:20:00
#PBS -q express
#
# cd to the directory where qsub was issued.
#
cd $PBS_O_WORKDIR
#
#
./TEST -keep -benchmark -param "-N 8" longpar
#
