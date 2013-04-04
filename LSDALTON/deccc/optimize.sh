#!/bin/bash 

# Script for DEC-MP2 geometry optimization.
#
# Requires the starting files:
# lsdalton.x    (to run the gradient calculation)
# rungeoopt.x   (to run the optimizer and get new geometry)
# geoopt.inp    (Input file for optimizer)
# grad.inp      (Input file for Dalton/DEC MP2 gradient calculation)
# MOLECULE.INP  (Starting geometry)

# Copy starting geometry to standalone reference file
cp MOLECULE.INP runopt.initMOLECULE

# Set control parameters
control=1
criteria=`grep "Geometry converged" LSDALTON.OUT | wc -l`
i="0"

# Start cycle
while [ ${criteria} -ne ${control} ] 
do
i=$[$i+1]
echo 'Starting step:' $i

# Delete old info
rm C* cmo* dens.restart fock.restart *info lcm* lsitem vdens* *~

# Run MP2 gradient calculation using geometry of MOLECULE.INP --> produces/updates runopt.history
cp grad.inp DALTON.INP
./lsdalton.x

# Save gradient output, just in case
mv LSDALTON.OUT grad.out$i

# Run optimizer using runopt.history. Next geometry is stored in MOLECULE.OUT
cp geoopt.inp DALTON.INP
./rungeoopt.x
cp MOLECULE.OUT MOLECULE.INP

# Grep to check for convergence
criteria=`grep "Geometry converged" LSDALTON.OUT | wc -l`
done

echo 'Geometry optimization converged! Number of steps: ' $i
