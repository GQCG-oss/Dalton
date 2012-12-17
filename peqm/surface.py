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

class Sphere(object):

# Adapted from code by The Little Grasshopper (http://prideout.net/blog/?p=44)

    def __init__(self, center=[0.0, 0.0, 0.0], radius=1.0, detail=1):
        self.center = center
        self.radius = radius
        self.detail = detail
        self.icosahedron()
        self.add_points()
        self.calculate()

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = [centroid[0] for centroid in self.centroids]
        y = [centroid[1] for centroid in self.centroids]
        z = [centroid[2] for centroid in self.centroids]
        ax.scatter(x, y, z, s=3, color='red')
        plt.show()

    def add_points(self):
        i = 1
        while i <= self.detail:
            i += 1
            self.subdivide()

    def icosahedron(self):
        """Construct an icosahedron with radius r"""
        phi = 26.56505
        phia = np.pi * phi / 180.0
        self.vertices = []
        self.vertices.append((self.center[0] + 0.0,
                              self.center[1] + 0.0,
                              self.center[2] + self.radius))
        theta = 0.0
        for i in xrange(1, 6):
            self.vertices.append(
                (self.center[0] + self.radius * np.cos(theta) * np.cos(phia),
                 self.center[1] + self.radius * np.sin(theta) * np.cos(phia),
                 self.center[2] + self.radius * np.sin(phia)))
            theta = theta + np.pi * 72.0 / 180.0
        theta = np.pi * 36.0 / 180.0
        for i in xrange(6, 11):
            self.vertices.append(
                (self.center[0] + self.radius * np.cos(theta) * np.cos(-phia),
                 self.center[1] + self.radius * np.sin(theta) * np.cos(-phia),
                 self.center[2] + self.radius * np.sin(-phia)))
            theta = theta + np.pi * 72.0 / 180.0
        self.vertices.append((self.center[0] + 0.0,
                              self.center[1] + 0.0,
                              self.center[2] - self.radius))
        self.faces = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 5, 1),
                      (11, 6, 7), (11, 7, 8), (11, 8, 9), (11, 9, 10),
                      (11, 10, 6), (1, 2, 6), (2, 3, 7), (3, 4, 8), (4, 5, 9),
                      (5, 1, 10), (6, 7, 2), (7, 8, 3), (8, 9, 4), (9, 10, 5),
                      (10, 6, 1)]

    def subdivide(self):
        """Subdivide each triangle into four triangles, pushing vertices to the
           sphere."""
        for faceIndex in xrange(len(self.faces)):
            face = self.faces[faceIndex]
            # Create three new vertices at the midpoints of each edge:
            a, b, c = (np.array(self.vertices[vertIndex])-np.array(self.center)
                       for vertIndex in face)
            d = np.array(((a + b) / np.linalg.norm(a + b)) * self.radius)
            e = np.array(((b + c) / np.linalg.norm(b + c)) * self.radius)
            f = np.array(((c + a) / np.linalg.norm(c + a)) * self.radius)
            self.vertices.append(tuple(self.center + d))
            self.vertices.append(tuple(self.center + e))
            self.vertices.append(tuple(self.center + f))
            # Split the current triangle into four smaller triangles:
            i = len(self.vertices) - 3
            j, k = i + 1, i + 2
            self.faces.append((i, j, k))
            self.faces.append((face[0], i, k))
            self.faces.append((i, face[1], j))
            self.faces[faceIndex] = (k, j, face[2])

    def calculate(self):
        """Calculate the areas and centroids of all triangles"""
        areas = []
        centroids = []
        for faceIndex in xrange(len(self.faces)):
            a, b, c = (np.array(self.vertices[vertIndex])
                       for vertIndex in self.faces[faceIndex])
            # calculate triangle area
            d = np.array([[a[0], b[0], c[0]], [a[1], b[1], c[1]], [1.0, 1.0, 1.0]])
            e = np.array([[a[1], b[1], c[1]], [a[2], b[2], c[2]], [1.0, 1.0, 1.0]])
            f = np.array([[a[2], b[2], c[2]], [a[0], b[0], c[0]], [1.0, 1.0, 1.0]])
            d = np.linalg.det(d)**2
            e = np.linalg.det(e)**2
            f = np.linalg.det(f)**2
            areas.append(0.5 * np.sqrt(d+e+f))
            # calculate triangle centroid
            xc = (1.0/3.0) * (a[0] + b[0] + c[0])
            yc = (1.0/3.0) * (a[1] + b[1] + c[1])
            zc = (1.0/3.0) * (a[2] + b[2] + c[2])
            centroids.append([xc, yc, zc])
        self.areas = areas
        self.centroids = centroids


class MolecularSurface(object):

    def __init__(self, elems=['X'], coords=[[0.0, 0.0, 0.0]], allcoords=[[0.0, 0.0, 0.0]],
                 detail=2, vdwfactor=1.0, vdwadd=1.0):
        self.elems = elems
        self.coords = coords
        self.allcoords = allcoords
        self.detail = detail
        self.vdwfactor = vdwfactor
        self.vdwadd = vdwadd
        vdws = []
        for elem in elems:
            vdws.append(self.vdwfactor * elem2vdw[elem])
        self.vdws = vdws
        self.add_midbonds()
        self.create_surface()

    def add_midbonds(self):
        midbonds = []
        elems = []
        vdws = []
        for i, (coord1, elem1) in enumerate(zip(self.coords, self.elems)):
            for j, (coord2, elem2) in enumerate(zip(self.coords, self.elems)):
                if j <= i:
                    continue
                elif self.bonded(coord1, elem1, coord2, elem2):
                    midbonds.append(self.midpoint(coord1, coord2))
                    elems.extend('X')
                    vdws.append((elem2vdw[elem1] + elem2vdw[elem2]) / 2.0)
        self.coords.extend(midbonds)
        self.elems.extend(elems)
        self.vdws.extend(vdws)

    def bonded(self, coord1, elem1, coord2, elem2):
        dist = np.linalg.norm(np.array(coord2) - np.array(coord1))
        if dist <= elem2vdw[elem1] + elem2vdw[elem2]:
            return True
        return False

    def midpoint(self, coord1, coord2):
        return [(comp1 + comp2) / 2.0 for comp1, comp2 in zip(coord1, coord2)]

    def create_surface(self):
        spheres = []
        for coord, radius in zip(self.coords, self.vdws):
            sphere = Sphere(coord, (radius + self.vdwadd), self.detail)
            spheres.append(sphere)
        self.spheres = spheres
        spheres = []
        for coord, radius in zip(self.coords, self.vdws):
            sphere = Sphere(coord, radius, self.detail)
            spheres.append(sphere)
        self.innerspheres = spheres
        self.remove_overlap()
        areas = []
        centroids = []
        for sphere in self.innerspheres:
            areas.extend(sphere.areas)
            centroids.extend(sphere.centroids)
        self.areas = areas
        self.centroids = centroids

    def inside(self, point, center, radius):
        """Return True if point is (almost) inside sphere"""
        r2 = ((point[0] - center[0])**2 +
              (point[1] - center[1])**2 +
              (point[2] - center[2])**2)
        return r2 <= (radius)**2

    def remove_overlap(self):
        """Remove vertices from spheres that are overlapping with other
        spheres"""
        for sphere, insphere in zip(self.spheres, self.innerspheres):
            for other in self.spheres:
                if other.center == sphere.center:
                    continue
                elif np.linalg.norm(np.array(sphere.center) - np.array(other.center)) > 5.0:
                    continue
                areas = []
                centroids = []
                inareas = []
                incentroids = []
                for inarea, incentroid, area, centroid in zip(insphere.areas, insphere.centroids,
                                                              sphere.areas, sphere.centroids):
                    if self.inside(centroid, other.center, other.radius):
                        continue
                    areas.append(area)
                    centroids.append(centroid)
                    inareas.append(inarea)
                    incentroids.append(incentroid)
                insphere.centroids = incentroids
                insphere.areas = inareas
                sphere.centroids = centroids
                sphere.areas = areas
            for center in self.allcoords:
                if center in self.coords:
                    continue
                elif np.linalg.norm(np.array(sphere.center) - np.array(center)) > 5.0:
                    continue
                areas = []
                centroids = []
                inareas = []
                incentroids = []
                for inarea, incentroid, area, centroid in zip(insphere.areas, insphere.centroids,
                                                              sphere.areas, sphere.centroids):
                    if self.inside(centroid, center, 2.0):
                        continue
                    areas.append(area)
                    centroids.append(centroid)
                    inareas.append(inarea)
                    incentroids.append(incentroid)
                insphere.centroids = incentroids
                insphere.areas = inareas
                sphere.centroids = centroids
                sphere.areas = areas

    def write_surface(self):
        """Write input files for PE module."""
        fout = open('surface.dat', 'w')
        fout.write('NPOINTS\n')
        fout.write('{}\n'.format(len(self.areas)))
        output = ''
        fout.write('coordinates\n')
        for area, centroid in zip(self.areas, self.centroids):
            output += ('{0[0]:12.6f}, {0[1]:12.6f}, {0[2]:12.6f}, '.format([float(comp)/au2aa for comp in centroid]) +
                       '{0:12.6f}\n'.format(float(area)/au2aa**2))
        fout.write(output)
        fout.close()


if __name__ == "__main__":

    fsurf = open(sys.argv[1], 'r')
    felems = open(sys.argv[2], 'r')
    fall = open(sys.argv[3], 'r')
    level = int(sys.argv[4])

    surfcoords = [[float(coord) * au2aa for coord in line.split()] for line in fsurf]
    elems = [charge2elem[int(line)] for line in felems]
    allcoords = [[float(coord) * au2aa for coord in line.split()] for line in fall]

    fsurf.close()
    felems.close()
    fall.close()

    mol = MolecularSurface(elems=elems, coords=surfcoords, allcoords=allcoords, detail=level, vdwfactor=1.2, vdwadd=1.4)
    mol.write_surface()

