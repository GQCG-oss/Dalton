#!/usr/local/bin/python
#
#  This program is a script that reads in a basisfile in a turbomole format
#  and converts it into an EMSL dalton format
# 
#  Run command: 
#  python turbo2dal.py input output
# 
#  Author: Y. M. Wang 
#  Date: Okt 2015
#
#  Feel free to make modifications
#
# Requirements:
#   https://pypi.python.org/pypi/periodictable
#   
# Linux or Mac OS:  
#   sudo pip install periodic
#   # 
#


import sys
import re
from periodic.table import element
import datetime
import numpy as np

ifile = open(sys.argv[1])
ofile = open(sys.argv[2],'w+')

nmol = 20
data = [] 
atoms = [[] for i in range(nmol)]
orbs = [[[] for i in range(7)] for j in range(nmol) ]
counter = [[0]*7 for i in range(nmol)]

for str in ifile:
   array = str.split()
   data.append(array)

switch = 0
iter = 0
natoms = 0
c = True

for str in data:
   m = re.match("\*",str[0])
   
   if m:
      switch = switch + 1
      b = re.match('\$end',data[iter+1][0])
      c = True

   if switch == 2:
      switch = 0

   if (m and c) and switch == 1 and not b:
      c = False
      atoms[0].append(data[iter+1][0])
      atoms[1].append(data[iter+2][2])
      basis = data[iter+1][1]
      natoms = natoms + 1

   if (switch == 0):
      b = re.match('\*',data[iter+1][0])
      if not b:
         s = re.match('s',data[iter+1][1])
         p = re.match('p',data[iter+1][1])
         d = re.match('d',data[iter+1][1])
         f = re.match('f',data[iter+1][1])
         g = re.match('g',data[iter+1][1])
         h = re.match('h',data[iter+1][1])
         i = re.match('i',data[iter+1][1])

         if s:
            orbs[natoms-1][0].append(data[iter+2][0])
	    counter[natoms-1][0] += 1
	 elif p:
	    orbs[natoms-1][1].append(data[iter+2][0])
	    counter[natoms-1][1] += 1
	 elif d:
	    orbs[natoms-1][2].append(data[iter+2][0])
	    counter[natoms-1][2] += 1
	 elif f:
	    orbs[natoms-1][3].append(data[iter+2][0])
	    counter[natoms-1][3] += 1
	 elif g:
	    orbs[natoms-1][4].append(data[iter+2][0])
	    counter[natoms-1][4] += 1
         elif h:
	    orbs[natoms-1][5].append(data[iter+2][0])
	    counter[natoms-1][5] += 1
	 elif i:
	    orbs[natoms-1][6].append(data[iter+2][0])
	    counter[natoms-1][6] += 1

   iter = iter + 1

format = "%a %b %d %H:%M:%S %Y"
today = datetime.datetime.today()

ofile.write('!  %s  EMSL  Basis Set Exchange Library   %s\n' % (basis,today) )
ofile.write('! Elements                             References\n')
ofile.write('! --------                             ----------\n')
ofile.write('! H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr: F. Weigend, A. Kohn, C. Hattig, Efficient use of the correlation consistent basis sets in resolution of the identity MP2 calculations, The Journal of Chemical Physics 116, 3175 (2002).\n')
ofile.write('! Li Be Na Mg: Christof Haettig, Optimization of auxiliary basis sets for RI-MP2 and RI-CC2 calculations: Core-valence and quintuple-? basis sets for H to Ar and QZVPP basis sets for Li to Kr, Physical Chemistry Chemical Physics 7, 59 (2005).\n')
ofile.write('!\n')
ofile.write('\n')
ofile.write('\n')
ofile.write('\n')
ofile.write('\n')
ofile.write('! Basis = %s\n' % (basis))

orbnames = ['s','p','d','f','g','h','i']
for n in range(natoms):
    string = element('%s' % (atoms[0][n]))
    atomname = (string.name).upper()
    data = atoms[1][n]
    b = re.findall('\d+\w',data[1:-1])
    configuration = ','.join(b)
    ofile.write('! %s       (%s)\n' % (atomname,configuration))
    ofile.write('! %s       (%s)\n' % (atomname,configuration))
    for i in range(7):
       if counter[n][i] > 0:
	  ofile.write('! %s functions\n' % (orbnames[i]))
          ofile.write('H    %s    %s\n' % (counter[n][i],counter[n][i]))
	  a = np.identity(counter[n][i])
	  for j in range(counter[n][i]):
	     ofile.write('      %s' % (orbs[n][i][j]))
	     for k in range(counter[n][i]):
	        ofile.write('      %.6f' % (a[j][k]))
		if k == 5:
		   ofile.write('\n')
	     ofile.write('\n')
	  

ifile.close
ofile.close
