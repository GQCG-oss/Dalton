#!/usr/bin/env python

import os, sys, glob, fnmatch, re
import math
sys.path.append("./pyparsing-1.5.5/")
from pyparsing import *


### Sort a dictionnary made of dictionaries
### example:
# # main = {"molecule1": {"p1": -23.3,
# #                       "p2":  34.2},
# #         "molecule2": {"p1":  17.4,
# #                       "p2": -22.2},
# #         "molecule3": {"p1":  -2.3,
# #                       "p2":  38.1},
# #         "molecule4": {"p1":   0.0,
# #                       "p2":  87.2}}
def sortDictOfDicts(dictOfDicts,str):
    ### str = name of the key to sort
    ### dictOfDicts = input dictionary to sort into a list
    items = sorted(dictOfDicts.items(), key = lambda tup: tup[1][str])
#     print type(items)
#     for item in items:
#         print item[0], item[1][str]
    return items

### Sort a dictionnary by its key
def sortDictByKey(mydict):
    items = sorted(mydict.iterkeys())
    return items

### Sort a dictionnary by its values
def sortDictKeysByValues(mydict):
    list_Keys = []
    for key, value in sorted(mydict.iteritems(), key=lambda (k,v): (v,k)):
        list_Keys.append(key)
    return list_Keys


### Function that replace some character in a string to be used in .tex files
def correctSpelling(str):
    newstr = str.replace("_", " ")
    return newstr


### Gaussian (or normal) distribution function
def normalDistrib(mu,std):
    def f(x):
        variance = (1.0*std)**2
        z = -0.5*(x-mu)**2/variance
        C = 1.0*math.sqrt(2.0*math.pi)*std
        expo = math.exp(z)
        return 1.0*expo/C
    return f

### print a dictionnary in a pretty way one line per item, with tab indentation for the depth
def pretty(d, indent=0):
    for key in sorted(d.iterkeys()):
        value = d[key]
        print '\t' * indent + str(key)
        if isinstance(value, dict):
            pretty(value, indent+1)
        else:
            print '\t' * (indent+1) + str(value)
# def pretty(d, indent=0):
#     for key, value in d.iteritems():
#         print '\t' * indent + str(key)
#         if isinstance(value, dict):
#             pretty(value, indent+1)
#         else:
#             print '\t' * (indent+1) + str(value)


### Example: get_fullpath_fileList('./','*.txt')
### ['/Users/patrime/Documents/CODE_DEVELOPMENT/ExtractionRESULTS_Python/parsingDALTON/titi.txt', '/Users/patrime/Documents/CODE_DEVELOPMENT/ExtractionRESULTS_Python/parsingDALTON/toto.txt']
def get_fullpath_oneFile(filename):
    return os.path.abspath(filename)

def get_fullpath_fileList(directory, extension):
    # print directory+extension
    fileList =  glob.glob(directory+extension)
    pathFilelist = []
    for file in fileList:
        pathFilelist.append(get_fullpath_oneFile(file))
    return pathFilelist



### Split name+extension of a filename
def splitFileName(filename):
    (nom, extension) = os.path.splitext(filename)
    return (nom, extension)

### Example: get_filenameList('./','*.txt')
### ['titi', 'toto']
def get_filenameList(directory, extension):
    filelist = []
    filenamelist = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, extension):
            filelist.append(file)
    for file in filelist:
        (nom, extension) = splitFileName(file)
        filenamelist.append(nom)
    return filenamelist


def is_filePresent(*argv):
    # args = sys.argv[1:]
    args = argv
    if len(args) != 1:
        print 'USAGE: python is_filePresent.py <path/file.extension>'
        sys.exit(-1)
    # filePathIn = sys.argv[1]
    filePathIn = argv[0]
    pathExist = os.path.isdir(os.path.dirname(filePathIn))
    if pathExist:
        fileExist = os.path.isfile(filePathIn)
        if fileExist:
            check = True
        else:
            print 'This file is NOT found: ', os.path.basename(filePathIn)
            check = False
    else:
        print 'This path is NOT found: ', os.path.dirname(filePathIn)
        check = False
    return check

def are_filesPresent(*argv):
    check = True
    for file in argv:
        check = check and is_filePresent(file)
    return check

def writeToFile(content,filename):
    f = open(filename, 'w')
    f.write(content)
    f.close()

def ensure_dir(dir):
    d = os.path.dirname(dir)
    if not os.path.exists(d):
        os.makedirs(d)
    return os.path.abspath(dir)

def getFileNameExtension(fullpath):
    path = os.path.dirname(fullpath)
    nameExtension = os.path.basename(fullpath) # nom.extension
    (nom, extension) = splitFileName(nameExtension)
    return (path, nom, extension)
