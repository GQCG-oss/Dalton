#!/usr/bin/env python

import sys

charge2elem = { 0:  'X',  1:  'H',  2: 'He',  3: 'Li',  4: 'Be',  5:  'B',
                6:  'C',  7:  'N',  8:  'O',  9:  'F', 10: 'Ne', 11: 'Na',
               12: 'Mg', 13: 'Al', 14: 'Si', 15:  'P', 16:  'S', 17: 'Cl',
               18: 'Ar', 19:  'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23:  'V',
               24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
               30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
               36: 'Kr', 37: 'Rb', 38: 'Sr', 39:  'Y', 40: 'Zr', 41: 'Nb',
               42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag',
               48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53:  'I',
               54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr',
               60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
               66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
               72: 'Hf', 73: 'Ta', 74:  'W', 75: 'Re', 76: 'Os', 77: 'Ir',
               78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi',
               84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra'}

au2aa = 0.5291772108

fin = open('{}'.format(sys.argv[1]), 'r')

line = fin.readline().split()
unit = line[0]

line = fin.readline().split()

nsites = int(line[0])
mulorder = int(line[1])
polorder = int(line[2])
lexlst = int(line[3])

if len(line) == 5:
    elems = int(line[4])
elif len(line) == 4:
    elems = 0
else:
    exit('Error reading line')

exlists = []
elements = []
coords = []
Q0s = []
Q1s = []
Q2s = []
Q3s = []
isoalphas = []
alphas = []

if mulorder == 3:
    pad = 20
elif mulorder == 2:
    pad = 10
elif mulorder == 1:
    pad = 4
elif mulorder == 0:
    pad = 1

for line in fin:
    line = line.split()
    if not line:
        continue
    exlist = [int(n) for n in line[0:lexlst]]
    exlists.append(exlist)
    coord = [float(x) for x in line[lexlst+elems:lexlst+elems+3]]
    coords.append(coord)
    if elems == 1:
        element = int(line[lexlst])
        elements.append(element)
    elif elems == 0:
        elements.append(0)
    if mulorder >= 0:
        Q0 = float(line[lexlst+elems+3])
        Q0s.append(Q0)
    if mulorder >= 1:
        Q1 = [float(x) for x in line[lexlst+elems+4:lexlst+elems+7]]
        Q1s.append(Q1)
    if mulorder >= 2:
        Q2 = [float(x) for x in line[lexlst+elems+7:lexlst+elems+13]]
        Q2s.append(Q2)
    if mulorder >= 3:
        Q3 = [float(x) for x in line[lexlst+elems+13:lexlst+elems+23]]
        Q3s.append(Q3)
    if polorder == 1:
        isoalpha = float(line[lexlst+elems+3+pad])
        isoalphas.append(isoalpha)
    elif polorder == 2:
        alpha = [float(x) for x in line[lexlst+elems+pad+3:lexlst+elems+pad+9]]
        alphas.append(alpha)

fin.close()

if unit == 'AU':
    coords = [[comp * au2aa for comp in coord] for coord in coords]

ndec = len(str(nsites)) + 1

fout = open('new_{}'.format(sys.argv[1]), 'w')
body = 'coordinates\n'
body += '{}\n'.format(nsites)
body += 'AA\n'
for i, coord in enumerate(coords):
    body += '{0} {1[0]:12.6f} {1[1]:12.6f} {1[2]:12.6f}\n'.format(charge2elem[elements[i]], coord)

if mulorder >= 0:
    body += 'monopoles\n'
    body += '{}\n'.format(nsites)
    for i, Q0 in enumerate(Q0s):
        body += '{0:{1}d} {2:12.6f}\n'.format(i+1, ndec, Q0)

if mulorder >= 1:
    body += 'dipoles\n'
    body += '{}\n'.format(nsites)
    for i, Q1 in enumerate(Q1s):
        body += '{0:{1}d} {2[0]:12.6f} {2[1]:12.6f} {2[2]:12.6f}\n'.format(i+1, ndec, Q1)

if mulorder >= 2:
    body += 'quadrupoles\n'
    body += '{}\n'.format(nsites)
    for i, Q2 in enumerate(Q2s):
        body += '{0:{1}d} {2[0]:12.6f} {2[1]:12.6f} {2[2]:12.6f} {2[3]:12.6f} {2[4]:12.6f} {2[5]:12.6f}\n'.format(i+1, ndec, Q2)

if mulorder >= 3:
    body += 'octopoles\n'
    body += '{}\n'.format(nsites)
    for i, Q3 in enumerate(Q3s):
        body += '{0:{1}d} {2[0]:12.6f} {2[1]:12.6f} {2[2]:12.6f} {2[3]:12.6f} {2[4]:12.6f} {2[5]:12.6f} {2[6]:12.6f} {2[7]:12.6f} {2[8]:12.6f} {2[9]:12.6f}\n'.format(i+1, ndec, Q3)

if polorder == 1:
    body += 'isoalphas\n'
    body += '{}\n'.format(nsites)
    for i, isoalpha in enumerate(isoalphas):
        body += '{0:{1}d} {2:12.6f}\n'.format(i+1, ndec, isoalpha)

if polorder == 2:
    body += 'alphas\n'
    body += '{}\n'.format(nsites)
    for i, alpha in enumerate(alphas):
        body += '{0:{1}d} {2[0]:12.6f} {2[1]:12.6f} {2[2]:12.6f} {2[3]:12.6f} {2[4]:12.6f} {2[5]:12.6f}\n'.format(i+1, ndec, alpha)

if exlists:
    body += 'exlists\n'
    body += '{}\n'.format(lexlst)
    for i, exlist in enumerate(exlists):
        for ex in exlist:
            body += ' {0:{1}}'.format(ex, ndec)
        body += '\n'

body += '\n'

fout.write(body)
