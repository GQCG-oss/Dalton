#!/usr/bin/env python

import os
import sys
import argparse as ap
import math
import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt

parser = ap.ArgumentParser(description='Cube analyzer 0.1',
                           epilog='Have a nice day :-)',
                           usage='%(prog)s [options]',
                           fromfile_prefix_chars='@')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('-i', dest='inputfiles', nargs='+', metavar='INPUT_FILE',
                    default=[],
                    help='''Specify the names of the cube files''')
parser.add_argument('-o', dest='outputfile', metavar='OUTPUT_FILE',
                    default='',
                    help='''Specify the name of the output file
                            [default: %(default)s]''')
parser.add_argument('-ref', dest='refidx', default=0, type=int,
                    metavar='REFCUBE',
                    help='''Specify the index of the reference cube to use in
                            different analysis. [default: %(default)s]''')
parser.add_argument('-log', dest='logfile', metavar='LOG_FILE',
                    default='cube_analyzer.log',
                    help='''Specify the name of the log file
                            [default: %(default)s]''')
parser.add_argument('-add', dest='addlist', nargs='+', default=[], type=int,
                    metavar=('CUBE1', 'CUBE2'),
                    help='''Specify which cubes to add. The cubes are numbered
                            according to the input order starting from 0''')
parser.add_argument('-sub', dest='sublist', nargs='+', default=[], type=int,
                    metavar=('CUBE1', 'CUBE2'),
                    help='''Specify which cubes to subtract. The cubes are
                            numbered according to the input order starting
                            from 0''')
parser.add_argument('-pln', dest='plnpts', nargs=9, default=[], type=float,
                    metavar=('X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2', 'X3', 'Y3', 'Z3'),
                    help='''Define plane cartesian coordinates''')
parser.add_argument('-mae', dest='mae', action='store_true', default=False,
                    help='''Calculate MAE with REFCUBE as the reference.''')
parser.add_argument('-vdw', dest='vdw', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a vdw analysis. RMSD is calculated relative to
                            a reference in volumes between MIN times vdw radius
                            and MAX times vdw radius in STEP steps''')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()
args = parser.parse_args()
if not args.inputfiles:
    parser.print_help()
    sys.exit()

charge2radius = { 1.0: 1.20,  2.0: 1.40,  3.0: 2.20,  4.0: 1.90,  5.0: 1.80,
                  6.0: 1.70,  7.0: 1.60,  8.0: 1.55,  9.0: 1.50, 10.0: 1.54,
                 11.0: 2.40, 12.0: 2.20, 13.0: 2.10, 14.0: 2.10, 15.0: 1.95,
                 16.0: 1.80, 17.0: 1.80, 18.0: 1.88, 19.0: 1.90,  0.0: 1.00}

aa2au = 1.8897261249935897

class Cube(object):

    def __init__(self):
        pass

    def read(self, filename):
        self.filename = filename
        self.readheader()
        self.readgrid()

    def readheader(self):
        fcube = open(self.filename, 'r')
        line = fcube.readline()
        line = fcube.readline()
        line = fcube.readline().split()
        self.natoms = int(line[0])
        self.origo = np.array([float(coord) for coord in line[1:4]])
        line = fcube.readline().split()
        self.xpoints = int(line[0])
        self.xstep = float(line[1])
        line = fcube.readline().split()
        self.ypoints = int(line[0])
        self.ystep = float(line[2])
        line = fcube.readline().split()
        self.zpoints = int(line[0])
        self.zstep = float(line[3])
        self.charges = []
        self.coords = []
        for i in xrange(self.natoms):
            line = fcube.readline().split()
            self.charges.append(float(line[1]))
            self.coords.append([float(coord) for coord in line[2:5]])
        self.grid = np.empty((self.xpoints, self.ypoints, self.zpoints))
        fcube.close()

    def copyheader(self, cube):
        self.natoms = cube.natoms
        self.origo = cube.origo
        self.xpoints = cube.xpoints
        self.xstep = cube.xstep
        self.ypoints = cube.ypoints
        self.ystep = cube.ystep
        self.zpoints = cube.zpoints
        self.zstep = cube.zstep
        self.charges = cube.charges
        self.coords = cube.coords
        self.grid = np.empty((self.xpoints, self.ypoints, self.zpoints))

    def readgrid(self):
        fcube = open(self.filename, 'r')
        line = fcube.readline()
        line = fcube.readline()
        line = fcube.readline().split()
        line = fcube.readline().split()
        line = fcube.readline().split()
        line = fcube.readline().split()
        for i in xrange(self.natoms):
            line = fcube.readline().split()
        z = 0
        for x in xrange(self.xpoints):
            for y in xrange(self.ypoints):
                for i in xrange(int(math.ceil(float(self.zpoints) / 6.0))):
                    line = fcube.readline().split()
                    for value in line:
                        self.grid[x][y][z] = float(value)
                        z += 1
                z = 0
        fcube.close()

    def writecube(self, filename):
        self.filename = filename
        fcube = open(self.filename, 'w')
        header = ''
        header += 'CUBE file\n'
        header += 'Generated by the CUBE Analyzer 0.1\n'
        header += '{0:5d}'.format(self.natoms)
        header += '{0[0]:12.6f}{0[1]:12.6f}{0[2]:12.6f}\n'.format(self.origo)
        header += '{0:5d}'.format(self.xpoints)
        header += '{0:12.6f}{1:12.6f}{1:12.6f}\n'.format(self.xstep, 0.0)
        header += '{0:5d}'.format(self.ypoints)
        header += '{1:12.6f}{0:12.6f}{1:12.6f}\n'.format(self.zstep, 0.0)
        header += '{0:5d}'.format(self.zpoints)
        header += '{1:12.6f}{1:12.6f}{0:12.6f}\n'.format(self.zstep, 0.0)
        for charge, coord in zip(self.charges, self.coords):
            header += '{0:5d}'.format(int(charge))
            header += '{0:12.6f}'.format(charge)
            header += '{0[0]:12.6f}{0[1]:12.6f}{0[2]:12.6f}\n'.format(coord)
        fcube.write(header)
        grid = ''
        for x in xrange(self.xpoints):
            for y in xrange(self.ypoints):
                i = 0
                for z in xrange(self.zpoints):
                    i += 1
                    grid += '{0:13.5e}'.format(self.grid[x][y][z])
                    if i == 6:
                        grid += '\n'
                        i = 0
                grid += '\n'
            fcube.write(grid)
            grid = ''
        fcube.close()

def plane_analysis(cube, pt1, pt2, pt3):
    fpln = open('test.dat', 'w')
    plnlog = ''
    plnfmt = '{0[0]:12.6f}{0[1]:12.6f}{0[2]:12.6f}{1:12.6f}\n'
    plane = []
    a = (pt1[1] * (pt2[2] - pt3[2]) + 
         pt2[1] * (pt3[2] - pt1[2]) + 
         pt3[1] * (pt1[2] - pt2[2]))
    b = (pt1[2] * (pt2[0] - pt3[0]) + 
         pt2[2] * (pt3[0] - pt1[0]) + 
         pt3[2] * (pt1[0] - pt2[0]))
    c = (pt1[0] * (pt2[1] - pt3[1]) + 
         pt2[0] * (pt3[1] - pt1[1]) + 
         pt3[0] * (pt1[1] - pt2[1]))
    d = - (pt1[0] * (pt2[1] * pt3[2] - pt3[1] * pt2[2]) +
           pt2[0] * (pt3[1] * pt1[2] - pt1[1] * pt3[2]) +
           pt3[0] * (pt1[1] * pt2[2] - pt2[1] * pt1[2]))
    point = [0.0, 0.0, 0.0]
    for x in xrange(cube.xpoints):
        for y in xrange(cube.ypoints):
            for z in xrange(cube.zpoints):
                point[0] = cube.origo[0] + x * cube.xstep
                point[1] = cube.origo[1] + y * cube.ystep
                point[2] = cube.origo[2] + z * cube.zstep
                if inplane(point, a, b, c, d):
                    plane.append([point, cube.grid[x][y][z]])
                    plnlog += plnfmt.format(point, cube.grid[x][y][z])
        plnlog += '\n'
    fpln.write(plnlog)
    fpln.close()

def inplane(point, a, b, c, d):
    s = a * point[0] + b * point[1] + c * point[2] + d
    if s < 0.1 and s > -0.1:
        return True
    else:
        return False

def vdw_analysis(cube, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum is less than minimum')
    rmsds = []
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        rmsds.append(0.0)
        grdpts.append(0)
        shells.append([round(inner, 4), round(outer, 4)])
        inner += step
        outer += step
    point = [0.0, 0.0, 0.0]
    for x in xrange(cube.xpoints):
        for y in xrange(cube.ypoints):
            for z in xrange(cube.zpoints):
                point[0] = cube.origo[0] + x * cube.xstep
                point[1] = cube.origo[1] + y * cube.ystep
                point[2] = cube.origo[2] + z * cube.zstep
                for idx, shell in enumerate(shells):
                    inner = shell[0]
                    outer = shell[1]
                    include = False
                    for center, charge in zip(cube.coords, cube.charges):
                        radius = charge2radius[charge] * aa2au
                        if (involume(point, center, inner * radius,
                                     outer * radius) and not
                            overlap(point, cube.coords, cube.charges, inner)):
                            include = True
                            break
                        else:
                            continue
                    if include:
                        rmsds[idx] += (cube.grid[x][y][z] - refcub.grid[x][y][z])**2
                        grdpts[idx] += 1
    fvdw = open('{}.log'.format(cube.filename[:-5]), 'w')
    vdw = 'Reference: {}\n'.format(refcub.filename)
    vd += '{}\n'.format(cube.filename)
    vdw += ' Points  Volume   Midpoint    RMSD\n'
    for idx, shell in enumerate(shells):
        inner = shell[0]
        outer = shell[1]
        vdw += '{0:6d} '.format(grdpts[idx])
        vdw += '{0:5.2f}-{1:<5.2f} '.format(inner, outer)
        vdw += '{0:5.2f} '.format(round(inner + 0.5 * step, 4))
        vdw += '{0:12.4e}\n'.format(math.sqrt(rmsds[idx] / grdpts[idx]))
    fvdw.write(vdw)
    fvdw.close()

def overlap(point, coords, charges, vdwfac):
    """"Return True if point is inside other vdw sphere"""
    for coord, charge in zip(coords, charges):
        if inside(point, coord, vdwfac * charge2radius[charge] * aa2au):
            return True
        else:
            return False

def involume(point, center, inner, outer):
    """return True if point is inside volume"""
    r2 = ((point[0] - center[0])**2 +
          (point[1] - center[1])**2 +
          (point[2] - center[2])**2)
    return (r2 < outer**2 and r2 >= inner**2)

def inside(point, center, radius):
    """Return True if point is inside sphere"""
    r2 = ((point[0] - center[0])**2 +
          (point[1] - center[1])**2 +
          (point[2] - center[2])**2)
    return r2 < radius**2

def outside(point, center, radius):
    """Return True if point is outside or on sphere"""
    r2 = ((point[0] - center[0])**2 +
          (point[1] - center[1])**2 +
          (point[2] - center[2])**2)
    return r2 >= radius**2

def mae_analysis(cube, refcub):
    mae = 0.0
    for x in xrange(refcub.xpoints):
        for y in xrange(refcub.ypoints):
            for z in xrange(refcub.zpoints):
                mae += abs(cube.grid[x][y][z] - refcub.grid[x][y][z])
    return mae / (refcub.xpoints * refcub.ypoints * refcub.zpoints)


if __name__ == "__main__":

    cubelist = []
    for cubefile in args.inputfiles:
        cube = Cube()
        cube.read(cubefile)
        cubelist.append(cube)
    
    newcube = Cube()
    newcube.copyheader(cubelist[0])

    if args.addlist:
        for index in args.addlist:
            newcube.grid += cubelist[index].grid
    
    if args.sublist:
        for index in args.sublist:
            newcube.grid -= cubelist[index].grid

    if args.mae:
        print('MAEs:')
        refcub = cubelist[args.refidx]
        print('Reference: {}'.format(refcub.filename))
        for cube in cubelist:
            if cube == refcub:
                continue
            mae = mae_analysis(cube, refcub)
            print('{0}: {1:13.5e}'.format(cube.filename, mae))

    if args.vdw:
        refcub = cubelist[args.refidx]
        for cube in cubelist:
            if cube == refcub:
                continue
            vdw_analysis(cube, refcub, *args.vdw)

    if args.plnpts:
        plane_analysis(cube, args.plnpts[0:3], args.plnpts[3:6], args.plnpts[6:9])

    if args.outputfile:
        newcube.writecube(args.outputfile)

