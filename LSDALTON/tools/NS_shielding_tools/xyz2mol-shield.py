#!/usr/bin/env python

#----------------------------------------------------------
# Title:   xyz2mol.py
# 
# Purpose: convert xyz file to DALTON MOLECULE.INP format for Nuclei-selected shielding calculations
          
# Usage:   $ ./xyz2mol-shield.py filename.xyz
#          $ python xyz2mol-shield.py filename.xyz
#
# Author:  This script is modification of Pablo Baudin's code. (Chandan Kumar)
# Date:    June 2015
#----------------------------------------------------------

import sys
import numpy as np
import pickle 

atoms=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',
       'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
       'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
       'I','Xe',]

verbose=False

# Check that user provided an xyz file in input
if (len(sys.argv)!=2):
    print 'You need to provide only one argument, the .xyz file:\n'
    print '$ ./xyz2mol-sheild.py file.xyz'
    print '$ python xyz2mol-sheild.py file.xyz'
    exit()

basis = raw_input('Enter basis set [cc-pVTZ]:\n')
if (basis==''):
    basis = 'cc-pVTZ\n'
else:
    basis+='\n'
totch = raw_input('Enter total charge of the molecule [0]:\n')
if (totch==''):
    totch = '0'
titl1 = raw_input('Enter title line 1 []:\n')
titl1+='\n'
titl2 = raw_input('Enter title line 2 []:\n')
titl2+='\n'

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
#print 'natoms'
#print natoms
# Check for title:
title=str(lines[1])
if (title!=''):
    if (titl1=='\n'):
        titl1=title
    elif (titl2=='\n'):
        titl2=title

coord=np.zeros((natoms,3))


o=2
labels=[]
listedatom=[]
for i in range(natoms):
    line = lines[i+o].split()
    atom = line[0]
    if (atom not in labels):
        listedatom.append(atom)

    labels.append(atom)

    for x in range(3):
        coord[i,x] = float(line[x+1])
#print 'listedatom'
#print listedatom 
#print 'labels'
#print labels
NumAtoms=str(natoms)
# Print MOLECULE,INP file:
if (verbose):
    print 'Print MOLECULE.INP in terminal and file '+name+'.mol\n'
    print 'BASIS'
    print basis
    print titl1
    print titl2
    print 'listedatom='+natoms+' Nosymmetry SubSystems Angstrom Charge='+totch  
output = open(name+'.mol','w')
output.write('BASIS\n')
output.write(basis)
output.write(titl1)
output.write(titl2)
output.write('Atomtypes=' +NumAtoms+ ' Angstrom SubSystems Charge='+totch+ '\n') 
for atype in listedatom:
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
    if (verbose):
        print 'Charge='+charge+' Atoms='+natomtype
    for i in range(natoms):
	  if(labels[i]==atype): 
                        output.write('Charge='+charge+' Atoms=1 SubSystem=NONNMR \n')
                        output.write((labels[i]+'  {:10.6f}  {:10.6f}  {:10.6f} \n').format(coord[i,0],coord[i,1],coord[i,2]))
    if (verbose):
         print (labels[i]+'  {:10.6f}  {:10.6f}  {:10.6f}').format(coord[i,0],coord[i,1],coord[i,2])
output.close()
print "(LS)DALTON molecule file: "+name+".mol has been created"


