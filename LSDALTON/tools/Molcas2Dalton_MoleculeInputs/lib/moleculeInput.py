#!/usr/bin/env python

import sys
import string

Atomic_NUMBERS={
"H": 1,
"He": 2,
"Li": 3,
"Be": 4,
"B": 5,
"C": 6,
"N": 7,
"O": 8,
"F": 9,
"Ne": 10,
"Na": 11,
"Mg": 12,
"Al": 13,
"Si": 14,
"P": 15,
"S": 16,
"Cl": 17,
"Ar": 18,
"K": 19,
"Ca": 20,
"Sc": 21,
"Ti": 22,
"V": 23,
"Cr": 24,
"Mn": 25,
"Fe": 26,
"Co": 27,
"Ni": 28,
"Cu": 29,
"Zn": 30,
"Ga": 31,
"Ge": 32,
"As": 33,
"Se": 34,
"Br": 35,
"Kr": 36,
"Rb": 37,
"Sr": 38,
"Y": 39,
"Zr": 40,
"Nb": 41,
"Mo": 42,
"Tc": 43,
"Ru": 44,
"Rh": 45,
"Pd": 46,
"Ag": 47,
"Cd": 48,
"In": 49,
"Sn": 50,
"Sb": 51,
"Te": 52,
"I": 53,
"Xe": 54,
"Cs": 55,
"Ba": 56,
"La": 57,
"Hf": 72,
"Ta": 73,
"W": 74,
"Re": 75,
"Os": 76,
"Ir": 77,
"Pt": 78,
"Au": 79,
"Hg": 80,
"Tl": 81,
"Pb": 82,
"Bi": 83,
"Po": 84,
"At": 85,
"Rn": 86,
"Fr": 87,
"Ra": 88,
"Ac": 89,
"Ku": 104,
"Ha": 105
}

class atomInfos:
    def __init__(self,atomSymbol,xCoord,yCoord,zCoord,atomCharge=None,unitDistances=None):
        self.atomSymbol = atomSymbol
        if atomCharge is None: 
            self.atomCharge = float(Atomic_NUMBERS[atomSymbol])
        else:
            self.atomCharge = float(atomCharge)
        self.unitDistances = unitDistances
        self.xCoord = float(xCoord)
        self.yCoord = float(yCoord)
        self.zCoord = float(zCoord)
    def setSymbol(self, symbol):
        self.atomSymbol = symbol
    def setCharge(self, charge):
        self.atomCharge = float(charge)
    def setAtomCoord(self, x,y,z):
        self.xCoord = float(x)
        self.yCoord = float(y)
        self.zCoord = float(z)
    def getContent_atomCoord(self):
        s = ''
        s += '{0:s}     {1: 7.14f}     {2: 7.14f}     {3: 7.14f}\n'.format(self.atomSymbol, self.xCoord, self.yCoord, self.zCoord)
        return s
    def print_atomCoord(self):
        print '{0}'.format(self.getContent_atomCoord())


class groupSameAtoms:
    def __init__(self):
        self.atomsSymbol = None
        self.atomTypeCharge = None
        self.nbAtomsInGroup = 0
        self.listAtomsCoord = []
    def addAtomInfo(self,atom):
        assert isinstance(atom,atomInfos), 'Trying to add something which is not an atomInfos object to a groupSameAtoms object'
        self.atomsSymbol = atom.atomSymbol
        self.atomTypeCharge = atom.atomCharge
        self.nbAtomsInGroup += 1
        self.listAtomsCoord.append(atom)
    def setAtomTypeCharge(self, charge):
        self.atomTypeCharge = float(charge)
    def getContent_DALTON_groupSameAtoms(self):
        s = ''
#        s += 'Charge={0:g} Atoms={1:d}\n'.format(self.atomTypeCharge, self.nbAtomsInGroup)
        s += 'Charge={0} Atoms={1:d}\n'.format(self.atomTypeCharge, self.nbAtomsInGroup)
        for a in self.listAtomsCoord:
            s += a.getContent_atomCoord()
        return s
    def print_DALTON_groupSameAtoms(self):
        print '{0}'.format(self.getContent_DALTON_groupSameAtoms())

class moleculeInfo:
    def __init__(self,name,comments="",unitDistances=None):
        self.moleculeName = name
        self.nbAtomsInMolecule = 0
        self.unitDistances = unitDistances
        self.comments = comments
        self.listAtoms          = []
        self.listGroupSameAtoms = []
    def setUnitDistances(self,unitDistances):
        self.unitDistances = unitDistances
    def addGroupSameAtomInfo(self, groupAtoms):
        assert isinstance(groupAtoms,groupSameAtoms), 'Trying to add something which is not an groupSameAtoms object to a molecule object'
        self.nbAtomsInMolecule += groupAtoms.nbAtomsInGroup
        self.listGroupSameAtoms.append(groupAtoms)
        for a in groupAtoms.listAtomsCoord:
            self.addAtomInfo(a)
    def addAtomInfo(self,atom):
        self.nbAtomsInMolecule += 1
        assert isinstance(atom,atomInfos), 'Trying to add something which is not an atomInfos object to a molecule object'
        self.listAtoms.append(atom)
        if not atom.unitDistances is None: # hoping that all atoms with unitDistance defined have actually the same :)
            self.setUnitDistances(atom.unitDistances)
    def getContent_DALTON_Molecule(self):
        s=''
        for group in self.listGroupSameAtoms:
            s += '{0}'.format(group.getContent_DALTON_groupSameAtoms())
        return s
    def print_DALTON_Molecule(self):
        print 'Name: {0}'.format(self.moleculeName)
        print 'nb. atoms: {0}'.format(self.nbAtomsInMolecule)
        print 'Unit: {0}'.format(self.unitDistances)
        print 'comments: {0}'.format(self.comments)
        print '{0}'.format(self.getContent_DALTON_Molecule())
    def getContent_MOLCAS_Molecule(self):
        s=''
        for atom in self.listAtoms:
            s += '{0}'.format(atom.getContent_atomCoord())
        return s
    def print_MOLCAS_Molecule(self):
        print '{0}'.format(self.getContent_MOLCAS_Molecule())
    def create_groupsSameAtomsDALTON(self): # for creating DALTON inputs
        listUniqueSymbols = []
        groups = []
        # find the list of unique symbols used in the molecule
        for a in self.listAtoms:
            if not any(a.atomSymbol==s for s in listUniqueSymbols):
                listUniqueSymbols.append(a.atomSymbol)
        # create a group of atoms for each unique atom type in the molecule
        for s in listUniqueSymbols:
            newGroup = groupSameAtoms()
            for a in self.listAtoms:
                if a.atomSymbol==s:
                    newGroup.addAtomInfo(a)
            groups.append(newGroup)    
        self.listGroupSameAtoms = groups
        return groups



class moleculeInput:
    def __init__(self,molecule,comments="",usingSymmetry="Nosymmetry",regBasis=None,auxBasis=None):
        assert isinstance(molecule,moleculeInfo), 'Trying to add something which is not a moleculeInfo object to a moleculeInput object'
        self.molecule = molecule
        self.regBasis = regBasis
        self.auxBasis = auxBasis
        self.comments = comments
        self.usingSymmetry = usingSymmetry
    def getContent_DALTON_MoleculeInput(self):
        s = 'BASIS\n'
        if self.regBasis is None:
            s += 'regbasis Aux=auxbasis\n'
        else:
            if self.auxBasis is None:
                s += '{0} Aux=auxbasis\n'.format(self.regBasis)                
            else:
                s += '{0} Aux={1}\n'.format(self.regBasis,self.auxBasis)                
        s += '{0}\n'.format(self.molecule.moleculeName)
        s += '{0} {1}\n'.format(self.comments, self.molecule.comments)
        assert not self.molecule.unitDistances is None, 'Trying to print molecule infos without knowing the unitDistance variable'
        groups = self.molecule.create_groupsSameAtomsDALTON()
        s += 'Atomtypes={0:d} {1} {2}\n'.format(len(self.molecule.listGroupSameAtoms), self.molecule.unitDistances, self.usingSymmetry)
        s += self.molecule.getContent_DALTON_Molecule()
        return s
    def print_DALTON_MoleculeInput(self):
        print '{0}'.format(self.getContent_DALTON_MoleculeInput())
    def print_MOLCAS_MoleculeInput(self):
        print '{0}'.format(self.getContent_MOLCAS_MoleculeInput())
    def getContent_MOLCAS_MoleculeInput(self):
        s=''
        s += '{0:d}\n'.format(self.molecule.nbAtomsInMolecule)
        if self.molecule.unitDistances is None:
            s += '{0}\n'.format(self.molecule.moleculeName)
        else:
            s += '{0} (in {1})\n'.format(self.molecule.moleculeName,self.molecule.unitDistances)
        s += self.molecule.getContent_MOLCAS_Molecule()
        s += '\n'
        return s
            
