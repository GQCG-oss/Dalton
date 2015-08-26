#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#***********************************************************#
# coord2dalmol.py, v0.1, Mikael Johansson 31.1.2008         #
#                                                           #
# Convert Turbomole coord to Dalton MOLECULE.INP format     #
#                                                           #
# Version history:                                          #
#   v0.1 : First release                                    #
#                                                           #
# Check for new versions and contact information at:        #
# http://www.iki.fi/~mpjohans/                              #
#***********************************************************#

import sys

VERSION='v0.1'
DATE='31.1.2008'
AAU=0.5291772108 # Ångström per a.u. (CODATA 2002)
ATOMS=[
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18  <-IUPAC column nr.
'Q' , #                                                                                     IUPAC row nr.
'H' ,                                                                                'He', # 1
'Li','Be',                                                  'B' ,'C' ,'N' ,'O' ,'F' ,'Ne', # 2
'Na','Mg',                                                  'Al','Si','P' ,'S' ,'Cl','Ar', # 3
'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', # 4
'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe', # 5
'Cs','Ba',                                                                                 # 6
          'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Ly',      # lanthanides
               'Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', # 6 cont.
'Fr','Ra',                                                                                 # 7
          'Ac','Th','Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      # actinides
               'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg']                                    # 7 cont.

class AtomClass:
    def __init__(self):
        self.type=None    # atom type, ex. 'H'
        self.charge=None  # atom type number (charge of nucleus), ex. 1
        self.x=None       # x-coordinate
        self.y=None         
        self.z=None

class ATypeClass:
    def __init__(self):
        self.type=None    # atom type, ex. 'H'
        self.charge=None  # atom type number (charge of nucleus), ex. 1
        self.count=None   # number of atoms of this type

        
def usage():
    print "coord2dalmol "+VERSION+" "+DATE+"\n"
    print "Converts a Turbomole coord file to Dalton MOLECULE.INP style."
    print "The output is for the moment dumped on stdout:"
    print "Usage: "+sys.argv[0]+" <cooordfile>"
    print "            <coord>: Turbomole coord file"

def init():
    args=sys.argv
        
    if len(args)<>2:
        sys.stderr.write("Error: exactly one argument expected\n")
        usage()
        sys.exit(1)
    
    coordfile=open(args[1],'r')
    dalfile=sys.stdout # For possible future feature enhancement

    return coordfile, dalfile
    

#*******************
# def GetAtoms()
#
# will return list 'atoms' of following structure:
# atoms[i]=[x-coord, y-coord, z-coord, atomtype(char), atomncharge]
# The atoms are sorted according to charge, from H, q=1.
#*********************************
def GetAtoms(coordfile):

    found=0
    line=coordfile.readline()
    while line:
        if line.find("$coord")<>-1:
            found=1
            break
        line=coordfile.readline()

    if not(found):
        sys.stderr.write("Error: unable to find keyword $coord\n")
        sys.exit(2)
    

    column=coordfile.readline().split()
    atoms=[]

    AN=0
    while column[0][0]<>'$':
        
        column[3]=column[3].capitalize()
        if ATOMS.count(column[3])==0:
            sys.stderr.write("Warning: Unknown atomtype found: "+column[3]+", replaced with a hydrogen\n")
            column[3]='H'
        atoms.append(AtomClass())
        atoms[AN].type=column[3]
        atoms[AN].charge=ATOMS.index(column[3])
        atoms[AN].x=float(column[0])
        atoms[AN].y=float(column[1])
        atoms[AN].z=float(column[2])
        column=coordfile.readline().split()
        AN=AN+1

    atoms.sort(lambda x, y: cmp(x.charge,y.charge))
    # Put a dummy atom at end for easier looping without bound checks:
    atoms.append(AtomClass())
    atoms[AN].type="Xx"

    return atoms

#************
# def GetAtomTypes()
#*******************

def GetAtomTypes(atoms):

    atypelist=[]
    atypecount=0
    aindex=0
    for atype in ATOMS:
        acount=0
        while atoms[aindex+acount].type==atype:
            acount=acount+1
        if acount>0:
            atypelist.append(ATypeClass())
            atypelist[atypecount].type=atype
            atypelist[atypecount].charge=ATOMS.index(atype)
            atypelist[atypecount].count=acount
            atypecount=atypecount+1
        aindex=aindex+acount

    return atypelist


#********
# def MakeFormula(atoms)
#**********************

def MakeFormula(atoms):

    return "DummyMolecule"


#***********
# def PutMolHeader()
# 
# defines the first lines of the output file.
#**************************

def PutMolHeader(formula, atomtypes, generators, charge, molfile):
    molfile.write("ATOMDF\n")
    molfile.write("%s\n" % (formula))
    molfile.write("created by coord2dalmol.py\n")
    molfile.write("Atomtypes=%d Generators=%s Charge=%d\n" %(atomtypes, generators, charge))


#***********
# def PutAtoms()
#**************

def PutAtoms(atoms, atypelist,  molfile):

    aindex=0
    for tindex in range(0,len(atypelist)):
        molfile.write("Charge=%0.1f  Atoms=%d  Bas=Turbomole-SVP  Aux=Ahlrichs-Coulomb-Fit\n" %
            (atypelist[tindex].charge, atypelist[tindex].count))
            
        while (atoms[aindex].type==atypelist[tindex].type):
            molfile.write("%s %20.14f  %20.14f  %20.14f\n" %
                ((atoms[aindex].type+`aindex+1`).ljust(8), atoms[aindex].x, atoms[aindex].y, atoms[aindex].z))
            aindex=aindex+1
  

def outit(coordfile, molfile):
    coordfile.close()
    molfile.close()
   
   
#*****
#************
# Main
#*************
#***********

coordfile, molfile=init()    

atoms=GetAtoms(coordfile)
formula=MakeFormula(atoms)
atypelist=GetAtomTypes(atoms)

generators="0"
charge=0
PutMolHeader(formula, len(atypelist), generators, charge, molfile)

PutAtoms(atoms, atypelist, molfile)

outit(coordfile, molfile)

