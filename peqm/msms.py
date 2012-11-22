#!/usr/bin/env python

import os
import sys
import numpy as np

au2aa = 0.5291772108

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

elem2vdw = { 'H': 1.20, 'He': 1.40, 'Li': 2.20, 'Be': 1.90,  'B': 1.80,
             'C': 1.70,  'N': 1.60,  'O': 1.55,  'F': 1.50, 'Ne': 1.54,
            'Na': 2.40, 'Mg': 2.20, 'Al': 2.10, 'Si': 2.10,  'P': 1.95,
             'S': 1.80, 'Cl': 1.80, 'Ar': 1.88, 'Br': 1.90,  'X': 1.50}

def centroid(a, b, c):
    """Calculate the centroid of a triangle"""
    xc = (a[0] + b[0] + c[0]) / 3.0
    yc = (a[1] + b[1] + c[1]) / 3.0
    zc = (a[2] + b[2] + c[2]) / 3.0
    return [xc, yc, zc]

def area(a, b, c):
    """Calculate the area of a triangle"""
    d = np.array([[a[0], b[0], c[0]], [a[1], b[1], c[1]], [1.0, 1.0, 1.0]])
    e = np.array([[a[1], b[1], c[1]], [a[2], b[2], c[2]], [1.0, 1.0, 1.0]])
    f = np.array([[a[2], b[2], c[2]], [a[0], b[0], c[0]], [1.0, 1.0, 1.0]])
    d = np.linalg.det(d)**2
    e = np.linalg.det(e)**2
    f = np.linalg.det(f)**2
    return 0.5 * np.sqrt(d + e + f)

def write_surface(centroids, areas):
    """Write surface file for PE module."""
    fout = open('surface.dat', 'w')
    fout.write('NPOINTS\n')
    fout.write('{}\n'.format(len(areas)))
    output = ''
    fout.write('coordinates\n')
    for area, centroid in zip(areas, centroids):
        output += ('{0[0]:12.6f}, {0[1]:12.6f}, {0[2]:12.6f}, '.format(centroid) +
                   '{0:12.6f}\n'.format(area))
    fout.write(output)
    fout.close()


if __name__ == "__main__":

    fvert = open(sys.argv[1], 'r')
    fface = open(sys.argv[2], 'r')

    verts = []
    for line in fvert:
        if line[0] == '#':
            continue
        data = line.split()
        verts.append([float(vert) for vert in data[0:3]])
    fvert.close()

    faces = []
    for line in fface:
        data = line.split()
        if line[0] == '#':
            continue
        elif len(line.split()) == 4:
            nfaces = data[0]
            nspheres = data[1]
            continue
        faces.append([(int(face) - 1) for face in data[0:3]])
    fface.close()

    centroids = []
    for face in faces:
        a = verts[face[0]]
        b = verts[face[1]]
        c = verts[face[2]]
        centroids.append(centroid(a, b, c))

    areas = []
    for face in faces:
        a = verts[face[0]]
        b = verts[face[1]]
        c = verts[face[2]]
        areas.append(area(a, b, c))
    write_surface(centroids, areas)
