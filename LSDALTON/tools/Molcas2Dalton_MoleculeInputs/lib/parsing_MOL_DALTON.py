#!/usr/bin/env python

import sys
import string
from pyparsing     import *
from libFiles      import *  # library to get full path of files, checking their existence, ...
from moleculeInput import * # (universal) class definition to store data from molecule input files


#------------------------------------------------------------------------
# Open a DALTON input file (XXX.mol) and extract information about inputs of 
# the molecule to investigate.
#------------------------------------------------------------------------
def parsing_MOL_DALTON(pathToFile):
    # --- Open/Read the file (XXX...XXX.mol)   
    if is_filePresent(pathToFile):
        infile = open(pathToFile,'r')
        file_str = infile.read()
        infile.close()   
        
        # --- Define Grammars for the parser (pyparsing)
        integer     = Word(nums)        ### '0123456789'
        StrangeName = Word(printables)  ### '0123456789abc...wxyzABC...WXYZ!"#$%&\\' ()*+,-./:;<=>?@[\\]^_`{|}~'
        number  = Word(alphanums)       ### 'abc...wyzABC...WXYZ0123456789'
        name    = Word(alphas)          ### 'abc...wxyzABC...WXYZ'

        endLine = Literal("\n")
        end     = Literal("\n").suppress()  # go to the end of the line, and suppress it, same as EOL?
        EOL     = LineEnd().suppress()      # go to the end of the line, and suppress it, same as end?
        all     = SkipTo(end)               # go to the end of the line, match the ENTIRE next line

        element = oneOf( """H He Li Be B C N O F Ne Na Mg Al Si P S Cl
            Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge
            As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag
            Cd In Sn Sb Te I Xe Cs Ba Lu Hf Ta W Re Os
            Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Lr Rf
            Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus
            Uuo La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm
            Yb Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No""" )
        
        get_infos =  Literal("BASIS").suppress() + EOL +  StrangeName.setResultsName("regBase") + Optional(Literal("Aux=")+StrangeName.setResultsName("auxBase")) + EOL + StrangeName.setResultsName("moleculeName") + all.setResultsName("comments") + EOL + Literal("Atomtypes=").suppress() + integer.setResultsName("nbAtomsTypes") + Optional(name).setResultsName("unitDistances") + Optional(name).setResultsName("hasSymmetry")
        
        coordinate = Combine((Optional(Literal("-"))+Optional(integer)+Literal(".")+integer))
        AtomCoordinates = element.setResultsName("atomAbrev") + coordinate.setResultsName("xCoord") +coordinate.setResultsName("yCoord") +coordinate.setResultsName("zCoord") + EOL

        #AtomCoordinates.ignore(get_infos)
        get_AtomTypeInfos  =  Literal("Charge=").suppress() + StrangeName.setResultsName("chargeAtom") + Literal("Atoms=").suppress() + integer.setResultsName("nbAtoms") + EOL
        get_sameAtomsCoord = AtomCoordinates.setResultsName("atomCoordinates")
        get_AtomsInfos = get_AtomTypeInfos + OneOrMore(get_sameAtomsCoord)
        logEntry = get_infos | get_AtomTypeInfos | get_sameAtomsCoord
        atomsEntry = get_AtomTypeInfos | get_sameAtomsCoord

        groupEntry = get_AtomTypeInfos.setResultsName("groupInfo") + OneOrMore(get_sameAtomsCoord).setResultsName("listSameAtoms")

        # --- EXTRACT THE DATA
        regBase      = None
        auxBase      = None
        moleculeName = None
        comments     = ""
        nbAtomsTypes = None
        unitDistances= None
        hasSymmetry  = "Nosymmetry"
        for tokens in get_infos.searchString(file_str):
            #print tokens.dump()
            if tokens.regBase:        regBase = tokens.regBase
            if tokens.auxBase:        auxBase = tokens.auxBase
            if tokens.moleculeName:   moleculeName = tokens.moleculeName
            if tokens.comments:       comments = tokens.comments
            if tokens.nbAtomsTypes:   nbAtomsTypes = tokens.nbAtomsTypes
            if tokens.unitDistances:  unitDistances = tokens.unitDistances
            if tokens.hasSymmetry:    hasSymmetry = tokens.hasSymmetry

        MyMolecule      = moleculeInfo(moleculeName,"",unitDistances) 

        charge     = None
        nbSameAtoms= None
        atomSymbol = None
        xCoord     = None
        yCoord     = None
        zCoord     = None
        #data = atomsEntry.searchString(file_str)
        #listOfGroups = []

#         test = groupEntry.searchString(file_str)
#         print test[0]
        
        

            


        for tokens in groupEntry.searchString(file_str):
            ##print tokens.dump()
            groupOfAtoms = groupSameAtoms()
            charge      = tokens[0]
            nbSameAtoms = int(tokens[1])
            for a in range(0,nbSameAtoms):
                atomSymbol  = tokens[0+a*4+2]
                xCoord      = tokens[1+a*4+2]
                yCoord      = tokens[2+a*4+2]
                zCoord      = tokens[3+a*4+2]
                atom = atomInfos(atomSymbol,xCoord,yCoord,zCoord,charge)
                groupOfAtoms.addAtomInfo(atom)
            MyMolecule.addGroupSameAtomInfo(groupOfAtoms)
        MyMoleculeInput = moleculeInput(MyMolecule,comments,hasSymmetry,regBase,auxBase)
        return MyMoleculeInput
