#!/usr/bin/python

#----------------------------------------------------------
# Title:   xyz2mol.py
# 
# Purpose: convert xyz file to another xyz file with ordered atoms 
#          (For ex. - All the hydrogen will be continously placed in the xyz file) 
# Usage:   $ ./sortxyz.py filename.xyz
#          $ python sortxyz.py filename.xyz
# Author:  This script is modification of Pablo Baudin's code. (Chandan Kumar)
#
# Date:    June 2015
#----------------------------------------------------------

import sys
import numpy as np

atoms=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',
       'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
       'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
       'I','Xe',]

verbose=False

# Check that user provided an xyz file in input
if (len(sys.argv)!=2):
    print 'You need to provide only one argument, the .xyz file:\n'
    print '$ ./sortxyz.py file.xyz'
    print '$ python sortxyz.py file.xyz'
    exit()

name = (sys.argv[1]).split('.')[0]

if (verbose):
    print '\nGenerating MOLECULE.INP file \n... \n'


# extract content of xyz file
xyzfile = open(sys.argv[1], "r")

lines = []
for line in xyzfile:
    lines.append(line)
xyzfile.close()

# Get number of atoms:
try:
    natoms = int(lines[0])
except:
    print 'file: '+sys.argv[1]+' does not fit the required format'
    exit()

coord=np.zeros((natoms,3))

o=2
labels=[]
Atomtypes=[]
for i in range(natoms):
    line = lines[i+o].split()
    atom = line[0]
    if (atom not in labels):
        Atomtypes.append(atom)

    labels.append(atom)

    for x in range(3):
        coord[i,x] = float(line[x+1])
NumAtoms=str(lines[0])
title=str(lines[1])
#print 'Atomtypes'
#print Atomtypes 
#print 'labels'
#print labels

if (verbose):
    print 'Print reordered in terminal and file '+name+'.mol\n'
    print titl1
output = open(name+'ordered.xyz','w')

#output.write(str(len(lines)))
output.write(NumAtoms)
output.write(title)

for atype in Atomtypes:
    try:
        charge=str(atoms.index(atype)+1)
        natomtype=str(labels.count(atype)) 
    except ValueError:
        print "Atom "+atype+" is not regognized by the script,"
        print "you might have to modify the script to include the full periodic table."
        # remove wrong output file:
        output.close()
        os.remove('./'+name+'.mol')
        exit()
#    output.write('Charge='+charge+' Atoms='+natomtype+' SubSystem='+Shield+'\n')
    if (verbose):
        print 'Charge='+charge+' Atoms='+natomtype
    for i in range(natoms):
	if(labels[i] ==atype): 
             output.write((labels[i]+'  {:10.6f}  {:10.6f}  {:10.6f} \n').format(coord[i,0],coord[i,1],coord[i,2]))
	     if (verbose):
                    print (labels[i]+'  {:10.6f}  {:10.6f}  {:10.6f}').format(coord[i,0],coord[i,1],coord[i,2])
output.close()
print ""
print "Another xyz file with sorted atoms "+name+"ordered.xyz has been created"


