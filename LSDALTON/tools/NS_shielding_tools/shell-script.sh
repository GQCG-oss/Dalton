#!/bin/bash 

## Input file and Python Path 
echo " Write the name of xyz input_file (Without .xyz extension) "
read Input_file
path=`pwd`
PYTHON="/usr/bin/python"
Script_Path_order_xyz="sort-xyz.py"
Script_Path_xyz_shielding_input="xyz2mol-shield.py"

## Sorting of atoms in xyz file 
echo ""
echo " ------  sorting atoms in xyzfile ------"
#$PYTHON $Script_Path_order_xyz amylose4.xyz
$PYTHON $Script_Path_order_xyz $Input_file.xyz

## Creating .mol file 
echo""
echo " ------ Construction of shielding input  (.mol) file ------"  
#$PYTHON $Script_Path_xyz_shielding_input amylose4.xyz
$PYTHON $Script_Path_xyz_shielding_input $Input_file.xyz

## Visualize sorted xyz file
echo ""
echo " ------ intiallising jmol to visulaise xyz file ------ "
PATH_TO_JMOLE=/Users/ckumar/Documents/jmol-14.2.7_2014.10.13/jmol.sh

Jmol_Input=$Input_file'ordered.xyz'
$PATH_TO_JMOLE $Jmol_Input

## Total No. of selected atoms for NMR shielding calculations
echo ""
echo " ------ Enter total number of selected Atoms ------ "
read Num_selected_atoms

## Serial no. of Selected atoms in sorted xyz file
FILENAME=$Input_file.mol
x=0
for i in `seq 1 $Num_selected_atoms`
do
echo "Enter Serial Number of selected-atoms"
read Serial_no_of_selected_atom
x=$((2*($Serial_no_of_selected_atom)+4))
echo "$x"
x1=$((($x)-1))
z=$x1
y=$x
sed -i.bk  "$z,$y s/NONNMR/NMR/"   $FILENAME > $FILENAME.bk
done 
rm $FILENAME.bk
