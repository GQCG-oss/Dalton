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
parser.add_argument('-mae', dest='mae', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a MAE analysis. MAE is calculated relative to
                            the reference cube in volumes between MIN times
                            vdw radius and MAX times vdw radius in STEP steps.''')
parser.add_argument('-mre', dest='mre', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a MRE analysis. MRE is calculated relative to
                            the reference cube in volumes between MIN times
                            vdw radius and MAX times vdw radius in STEP steps.''')
parser.add_argument('-print-shells', dest='printshells', nargs=3, default=[],
                    type=float, metavar=('MIN','MAX','STEP'),
                    help='''Print difference values relative to the reference
                            cube divided in volumes between MIN times vdw radius
                            and MAX times vdw radius in STEP steps.''')
parser.add_argument('-rmsd', dest='rmsd', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a RMSD analysis. RMSD is calculated relative to
                            the reference cube in volumes between MIN times
                            vdw radius and MAX times vdw radius in STEP steps''')
parser.add_argument('-nrmsd', dest='nrmsd', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a normalized RMSD analysis. RMSD is calculated
                            relative to the reference cube in volumes between
                            MIN times vdw radius and MAX times vdw radius in
                            STEP steps''')
parser.add_argument('-ncores', dest='ncores', default=1, type=int,
                    metavar='NCORES',
                    help='''Number of cores [default: %(default)s]''')

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

def rmsd_analysis(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
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
    nprocs = args.ncores
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=rmsd_worker,
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
        frmsd = open('{}_rmsd.log'.format(cube.filename[:-5]), 'w')
        rmsd = 'Reference: {}\n'.format(refcub.filename)
        rmsd += '{}\n'.format(cube.filename)
        rmsd += ' Points   vdW Volume   Midpoint    RMSD\n'
        for ish, shell in enumerate(shells):
            inner = shell[0]
            outer = shell[1]
            rmsd += '{0:8d} '.format(grdpts[ish])
            rmsd += '{0:5.2f}-{1:<5.2f} '.format(inner, outer)
            rmsd += '{0:5.2f} '.format(round(inner + 0.5 * step, 4))
            rmsd += '{0:12.4e}\n'.format(math.sqrt(rmsds[ic][ish] / grdpts[ish]))
        frmsd.write(rmsd)
        frmsd.close()

def rmsd_worker(points, shells, rmsds, grdpts, refcub, cubelist, out_queue):
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

def nrmsd_analysis(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    nrmsds = []
    for i in range(len(cubelist)):
        nrmsds.append([])
    grdpts = []
    maxvals = []
    minvals = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            nrmsds[i].append(0.0)
        grdpts.append(0)
        maxvals.append(0.0)
        minvals.append(0.0)
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
    nprocs = args.ncores
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=nrmsd_worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, nrmsds, grdpts, maxvals, minvals,
                                refcub, cubelist, out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for nrmsd, pts, maxval, minval in results:
        for i in xrange(len(shells)):
            maxvals[i] = max(maxvals[i], maxval[i])
            minvals[i] = min(minvals[i], minval[i])
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                nrmsds[j][i] += nrmsd[j][i]
    for ic, cube in enumerate(cubelist):
        fnrmsd = open('{}_nrmsd.log'.format(cube.filename[:-5]), 'w')
        nrmsd = 'Reference: {}\n'.format(refcub.filename)
        nrmsd += '{}\n'.format(cube.filename)
        nrmsd += ' Points   vdW Volume   Midpoint    NRMSD\n'
        for ish, shell in enumerate(shells):
            inner = shell[0]
            outer = shell[1]
            nrmsd += '{0:8d} '.format(grdpts[ish])
            nrmsd += '{0:5.2f}-{1:<5.2f} '.format(inner, outer)
            nrmsd += '{0:5.2f} '.format(round(inner + 0.5 * step, 4))
            nrmsd += '{0:12.4e}\n'.format(math.sqrt(nrmsds[ic][ish] / grdpts[ish])
                                            / (maxvals[ish] - minvals[ish]))
        fnrmsd.write(nrmsd)
        fnrmsd.close()

def nrmsd_worker(points, shells, nrmsds, grdpts, maxvals, minvals,
                 refcub, cubelist, out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    nrmsds[ic][ish] += (cube.grid[x][y][z] -
                                        refcub.grid[x][y][z])**2
                    if refcub.grid[x][y][z] > maxvals[ish]:
                        maxvals[ish] = refcub.grid[x][y][z]
                    if refcub.grid[x][y][z] < minvals[ish]:
                        minvals[ish] = refcub.grid[x][y][z]
                grdpts[ish] += 1
                break
    out_queue.put((nrmsds, grdpts, maxvals, minvals))

def mae_analysis(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    maes = []
    for i in range(len(cubelist)):
        maes.append([])
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            maes[i].append(0.0)
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
    nprocs = args.ncores 
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=mae_worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, maes, grdpts, refcub, cubelist, out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for mae, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                maes[j][i] += mae[j][i]
    for ic, cube in enumerate(cubelist):
        fmae = open('{}_mae.log'.format(cube.filename[:-5]), 'w')
        mae = 'Reference: {}\n'.format(refcub.filename)
        mae += '{}\n'.format(cube.filename)
        mae += ' Points   vdW Volume   Midpoint    MAE\n'
        for ish, shell in enumerate(shells):
            inner = shell[0]
            outer = shell[1]
            mae += '{0:8d} '.format(grdpts[ish])
            mae += '{0:5.2f}-{1:<5.2f} '.format(inner, outer)
            mae += '{0:5.2f} '.format(round(inner + 0.5 * step, 4))
            mae += '{0:12.4e}\n'.format(maes[ic][ish] / grdpts[ish])
        fmae.write(mae)
        fmae.close()

def mae_worker(points, shells, maes, grdpts, refcub, cubelist, out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    maes[ic][ish] += np.abs(cube.grid[x][y][z] -
                                       refcub.grid[x][y][z])
                grdpts[ish] += 1
                break
    out_queue.put((maes, grdpts))

def mre_analysis(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    mres = []
    for i in range(len(cubelist)):
        mres.append([])
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            mres[i].append(0.0)
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
    nprocs = args.ncores 
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=mre_worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, mres, grdpts, refcub, cubelist, out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for mre, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                mres[j][i] += mre[j][i]
    for ic, cube in enumerate(cubelist):
        fmre = open('{}_mre.log'.format(cube.filename[:-5]), 'w')
        mre = 'Reference: {}\n'.format(refcub.filename)
        mre += '{}\n'.format(cube.filename)
        mre += ' Points   vdW Volume   Midpoint    MAE\n'
        for ish, shell in enumerate(shells):
            inner = shell[0]
            outer = shell[1]
            mre += '{0:8d} '.format(grdpts[ish])
            mre += '{0:5.2f}-{1:<5.2f} '.format(inner, outer)
            mre += '{0:5.2f} '.format(round(inner + 0.5 * step, 4))
            mre += '{0:12.4e}\n'.format(mres[ic][ish] / grdpts[ish])
        fmre.write(mre)
        fmre.close()

def mre_worker(points, shells, mres, grdpts, refcub, cubelist, out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    mres[ic][ish] += np.abs(1.0 - (cube.grid[x][y][z] /
                                                   refcub.grid[x][y][z]))
                grdpts[ish] += 1
                break
    out_queue.put((mres, grdpts))

def print_shells(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    vals = []
    rels = []
    for i in range(len(cubelist)):
        vals.append([])
        rels.append([])
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            vals[i].append([])
            rels[i].append([])
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
    nprocs = args.ncores 
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=print_shell_worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, vals, rels, grdpts, refcub, cubelist,
                                out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for val, rel, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                for k in xrange(pts[i]):
                    vals[j][i].append(val[j][i][k])
                    rels[j][i].append(rel[j][i][k])
    for ic, cube in enumerate(cubelist):
        for ish, shell in enumerate(shells):
            fprt = open('{0}_shell_{1}.log'.format(cube.filename[:-5], ish), 'w')
            prt = 'Reference: {0}\n'.format(refcub.filename)
            prt += '{0}\n'.format(cube.filename)
            prt += 'Inner-outer radii: {0[0]:5.2f}-{0[1]:<5.2f}\n'.format(shell)
            prt += 'Shell midpoint: {0:5.2f}\n'.format(round(shell[0] + 0.5 * step, 4))
            prt += 'Number of grid points: {0:9d}\n'.format(grdpts[ish])
            fprt.write(prt)
            prt = ''
            for i in xrange(grdpts[ish]):
                prt += '{0:12.4e} {1:12.4e}\n'.format(vals[ic][ish][i],
                                                      rels[ic][ish][i])
            fprt.write(prt)
            fprt.close()

def print_shell_worker(points, shells, vals, rels, grdpts, refcub, cubelist,
                       out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    vals[ic][ish].append(cube.grid[x][y][z]
                                         - refcub.grid[x][y][z])
                    rels[ic][ish].append(abs((refcub.grid[x][y][z]
                                              - cube.grid[x][y][z])
                                             / refcub.grid[x][y][z]))
                grdpts[ish] += 1
                break
    out_queue.put((vals, rels, grdpts))

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
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        mae_analysis(cubana, refcub, *args.mae)
    
    if args.mre:
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        mre_analysis(cubana, refcub, *args.mre)

    if args.rmsd:
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        rmsd_analysis(cubana, refcub, *args.rmsd)

    if args.nrmsd:
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        nrmsd_analysis(cubana, refcub, *args.nrmsd)

    if args.printshells:
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        print_shells(cubana, refcub, *args.printshells)

    if args.outputfile:
        newcube.writecube(args.outputfile)

