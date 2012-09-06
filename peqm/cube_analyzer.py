#!/usr/bin/env python

import os
import sys
import argparse as ap
import math
import numpy as np
import time
import multiprocessing as mp

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
        self.grid = np.zeros((self.xpoints, self.ypoints, self.zpoints))
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
        self.grid = np.zeros((self.xpoints, self.ypoints, self.zpoints))

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
        header += '{1:12.6f}{0:12.6f}{1:12.6f}\n'.format(self.ystep, 0.0)
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
                if grid[-1] != '\n':
                    grid += '\n'
            fcube.write(grid)
            grid = ''
        fcube.close()

def vdw_analysis(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum is less than minimum')
    rmsds = []
    for i in range(len(cubelist)):
        rmsds.append([])
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            rmsds[i].append(0.0)
        grdpts.append(0)
        shells.append([round(inner, 4), round(outer, 4)])
        inner += step
        outer += step
    points = []
    point = [0.0, 0.0, 0.0]
    for x in xrange(refcub.xpoints):
        point[0] = refcub.origo[0] + x * refcub.xstep
        for y in xrange(refcub.ypoints):
            point[1] = refcub.origo[1] + y * refcub.ystep
            for z in xrange(refcub.zpoints):
                point[2] = refcub.origo[2] + z * refcub.zstep
                points.append([list(point), x, y, z])
    out_queue = mp.Queue()
    nprocs = mp.cpu_count()
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, rmsds, grdpts, refcub, cubelist, out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for rmsd, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                rmsds[j][i] += rmsd[j][i]
    for ic, cube in enumerate(cubelist):
        fvdw = open('{}.log'.format(cube.filename[:-5]), 'w')
        vdw = 'Reference: {}\n'.format(refcub.filename)
        vdw += '{}\n'.format(cube.filename)
        vdw += ' Points   vdW Volume   Midpoint    RMSD\n'
        for ish, shell in enumerate(shells):
            inner = shell[0]
            outer = shell[1]
            vdw += '{0:8d} '.format(grdpts[ish])
            vdw += '{0:5.2f}-{1:<5.2f} '.format(inner, outer)
            vdw += '{0:5.2f} '.format(round(inner + 0.5 * step, 4))
            vdw += '{0:12.4e}\n'.format(math.sqrt(rmsds[ic][ish] / grdpts[ish]))
        fvdw.write(vdw)
        fvdw.close()

def worker(points, shells, rmsds, grdpts, refcub, cubelist, out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    rmsds[ic][ish] += (cube.grid[x][y][z] -
                                       refcub.grid[x][y][z])**2
                grdpts[ish] += 1
                break
    out_queue.put((rmsds, grdpts))
 
def inorout(point, shell, coords, charges):
    radii = [charge2radius[charge] * aa2au for charge in charges]
    for center, radius in zip(coords, radii):
        r2 = ((point[0] - center[0])**2 +
              (point[1] - center[1])**2 +
              (point[2] - center[2])**2)
        if r2 < (shell[1] * radius)**2 and r2 >= (shell[0] * radius)**2:
            for coord, rad in zip(coords, radii):
                r2 = ((point[0] - coord[0])**2 +
                      (point[1] - coord[1])**2 +
                      (point[2] - coord[2])**2)
                if r2 > (shell[0] * rad)**2:
                    continue
                else:
                    return False
            return True
        else:
            continue
    return False

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
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        vdw_analysis(cubana, refcub, *args.vdw)

    if args.outputfile:
        newcube.writecube(args.outputfile)

