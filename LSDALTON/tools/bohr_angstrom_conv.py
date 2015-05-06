#!/usr/bin/python
#
#----------------------------------------------------------
# Title:   bohr_angstrom_conv.py
# 
# Purpose: Script for conversion of xyz file in Bohr to Aangstrom
#
# Author:  Pablo Baudin
# Date:    Mar. 2015
#----------------------------------------------------------

import sys
import numpy as np


# factor of conversion:
bohr2angstrom = 0.529177249
angstrom2bohr = 1.889725989

# Check that user provided an xyz file in input
if (len(sys.argv)!=2):
    print 'You need to provide only one argument the .xyz file \n'
    print '$ ./bohr_angstrom_conv.py file.xyz'
    print '$ python bohr_angstrom_conv.py file.xyz'
    exit()


# Get basis total charge and title
conv = raw_input('For bohr to angstrom conversion type 1 \n'
        'for angstrom to bohr type 0 [1]:\n')

if (conv==''):
    conv = bohr2angstrom
    ext="_angs.xyz"
    title="(Angstrom)"
    print "Bohr to angstrom conversion is performed"
else:
    try:
        if (int(conv)==1):
            conv = bohr2angstrom
            ext="_angs.xyz"
            title="(Angstrom)"
            print "Bohr to angstrom conversion is performed"
        elif (int(conv)==0):
            conv = angstrom2bohr
            ext="_bohr.xyz"
            title="(Bohr)"
            print "Angstrom to bohr conversion is performed"
        else:
            print "Wrong input, allowed are: <empty>, 0 or 1"
            exit()
    except:
        print "Wrong input, allowed are: <empty>, 0 or 1"
        exit()


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

# Check for title:
title+=str(lines[1])

coord=np.zeros((natoms,3))

# Open output file:
name = (sys.argv[1]).split('.')[0]+ext
output = open(name,"w")

output.write(str(natoms)+"\n")
output.write(title)

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
        # read and convert to abgstroms
        coord[i,x] = float(line[x+1])*conv

    # print to output file
    output.write((labels[i]+'  {:10.6f}  {:10.6f}  {:10.6f} \n').format(coord[i,0],coord[i,1],coord[i,2]))


