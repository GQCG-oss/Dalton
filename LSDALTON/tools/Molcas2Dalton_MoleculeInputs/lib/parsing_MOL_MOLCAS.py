#!/usr/bin/env python

import sys
import string
sys.path.append("../../pyparsing-1.5.5/")
from pyparsing     import *
sys.path.append("../../")
from libFiles      import *  # library to get full path of files, checking their existence, ...
from moleculeInput import * # (universal) class definition to store data from molecule input files




#------------------------------------------------------------------------
# Open a MOLCAS input file and extract information about inputs of 
# the molecule to investigate.
#------------------------------------------------------------------------
def parsing_MOL_MOLCAS(pathToFile):
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
        

        coordinate = Combine((Optional(Literal("-"))+Optional(integer)+Literal(".")+integer))
        AtomCoordinates = element.setResultsName("atomAbrev") + coordinate.setResultsName("xCoord") +coordinate.setResultsName("yCoord") +coordinate.setResultsName("zCoord") + EOL
        atomsEntry = AtomCoordinates.setResultsName("atomCoordinates")

        get_infos = StringStart() + integer.setResultsName("nbTotAtoms")  + Optional(StrangeName.setResultsName("moleculeName") + Optional(Literal("(in") + name.setResultsName("unitDistances") + Literal(")"))) + AtomCoordinates.setResultsName("firstAtomInfo")
        # --- EXTRACT THE DATA
        moleculeName = None
        nbTotAtoms = None
        unitDistances= None
        moleculeNameDist = None
        for tokens in get_infos.searchString(file_str):
            #print tokens.dump()
            if tokens.moleculeName:   moleculeName = tokens.moleculeName
            if tokens.unitDistances:  unitDistances = tokens.unitDistances
            if tokens.nbTotAtoms:     nbTotAtoms = tokens.nbTotAtoms
           
        if moleculeName == None:
            moleculeName = getFileNameExtension(pathToFile)[1]
        if unitDistances == None:
            unitDistances = "Angstrom" # default unit in Molcas
            
        MyMolecule      = moleculeInfo(moleculeName,"",unitDistances) 
        #MyMolecule.print_DALTON_Molecule()
        atomSymbol = None
        xCoord     = None
        yCoord     = None
        zCoord     = None
        for tokens in atomsEntry.searchString(file_str):
            atomSymbol = tokens[0]
            xCoord     = tokens[1]
            yCoord     = tokens[2]
            zCoord     = tokens[3]
            #print xCoord," ",yCoord," ",zCoord," \n"
            atom = atomInfos(atomSymbol,xCoord,yCoord,zCoord)
            MyMolecule.addAtomInfo(atom)

        MyMoleculeInput = moleculeInput(MyMolecule,"")
        return MyMoleculeInput


#myMol = parsing_MOL_MOLCAS("./BeF202H4.xyz")
#myMol = parsing_MOL_MOLCAS("./amylose_16.xyz")
#print myMol.getContent_DALTON_MoleculeInput()
