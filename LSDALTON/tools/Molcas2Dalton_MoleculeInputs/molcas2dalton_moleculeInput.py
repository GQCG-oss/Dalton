#!/usr/bin/env python

import sys
sys.path.append("./lib")
sys.path.append("./lib/pyparsing-1.5.5/")
from libFiles import *
from parsing_MOL_DALTON import *
from parsing_MOL_MOLCAS import *


# directory containing the files to convert
pathToInputFile = sys.argv[1]

# directory containing the files converted
pathToConvertedInputFile = ensure_dir(sys.argv[2]+"/")


listFiles =  get_fullpath_fileList(pathToInputFile,'*.xyz')

nbMolecule = 0
nbconverted = 0
listConvertedMolecule = []
liste = ''
for f in listFiles: # f > fullpath/nom.erxtension
    nbMolecule += 1
    (path, nom, extension) = getFileNameExtension(f)
    MolcasMoleculeInput = parsing_MOL_MOLCAS(f)
    content = MolcasMoleculeInput.getContent_DALTON_MoleculeInput()
    
    doConvert = True
    print "molecule: ",MolcasMoleculeInput.molecule.moleculeName
    i = 0
    for group in MolcasMoleculeInput.molecule.listGroupSameAtoms:
        charge = group.atomTypeCharge
        i +=1
        # if charge > 30.0:
        #     doConvert = False
        #     print "\t charge atom({0:d}): {1:.1f}  (NOT CONVERTED)".format(i,charge)
        # else:
        #     print "\t charge atom({0:d}): {1:.1f}".format(i,charge)
    if doConvert:
        writeToFile(content,pathToConvertedInputFile+'/'+nom+'.mol')
        nbconverted += 1
        listConvertedMolecule.append(nom+" ")
    #print "{0:d} molecules converted out of {1:d}".format(nbconverted,nbMolecule)
    print "{0:d} molecules converted".format(nbconverted)
print liste.join(listConvertedMolecule)
