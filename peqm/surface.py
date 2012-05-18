#!/usr/bin/env python
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

# TODO: tuples vs. lists?

class Sphere(object):

# Adapted from code by The Little Grasshopper (http://prideout.net/blog/?p=44)

    def __init__(self, center=[0.0, 0.0, 0.0], radius=1.0):
        self.center = center
        self.radius = radius
        self.area = 4.0 * np.pi * radius**2
        self.icosahedron()

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = [vertice[0] for vertice in self.vertices]
        y = [vertice[1] for vertice in self.vertices]
        z = [vertice[2] for vertice in self.vertices]
        ax.scatter(x, y, z, s=3, color='red')
        plt.show()

    def add_points(self, level=1):
        i = 1
        while i <= level:
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
        self.faces = [(2, 1, 0), (3, 2, 0), (4, 3, 0), (5, 4, 0), (1, 5, 0),
                      (11, 6, 7), (11, 7, 8), (11, 8, 9), (11, 9, 10),
                      (11, 10, 6), (1, 2, 6), (2, 3, 7), (3, 4, 8), (4, 5, 9),
                      (5, 1, 10), (2, 7, 6), (3, 8, 7), (4, 9, 8), (5, 10, 9),
                      (1, 6, 10)]
        self.area = (4.0 * np.pi * self.radius**2) / len(self.vertices)

    def subdivide(self):
        """Subdivide each triangle into four triangles, pushing vertices to the
           sphere of radius r"""
        for faceIndex in xrange(len(self.faces)):
            # Create three new vertices at the midpoints of each edge:
            face = self.faces[faceIndex]
            a, b, c = (np.array(self.vertices[vertIndex])-np.array(self.center)
                       for vertIndex in face)
            d = (((a + b) / np.linalg.norm(a + b)) * self.radius)
            e = (((b + c) / np.linalg.norm(b + c)) * self.radius)
            f = (((a + c) / np.linalg.norm(a + c)) * self.radius)
            self.vertices.append(tuple(np.array(self.center) + d))
            self.vertices.append(tuple(np.array(self.center) + e))
            self.vertices.append(tuple(np.array(self.center) + f))
            # Split the current triangle into four smaller triangles:
            i = len(self.vertices) - 3
            j, k = i + 1, i + 2
            self.faces.append((i, j, k))
            self.faces.append((face[0], i, k))
            self.faces.append((i, face[1], j))
            self.faces[faceIndex] = (k, j, face[2])
        self.area = (4.0 * np.pi * self.radius**2) / len(self.vertices)


class MolecularSurface(object):

    elem2vdw = { 'H': 1.20, 'He': 1.40, 'Li': 2.20, 'Be': 1.90,  'B': 1.80,
                 'C': 1.70,  'N': 1.60,  'O': 1.55,  'F': 1.50, 'Ne': 1.54,
                'Na': 2.40, 'Mg': 2.20, 'Al': 2.10, 'Si': 2.10,  'P': 1.95,
                 'S': 1.80, 'Cl': 1.80, 'Ar': 1.88, 'Br': 1.90,  'X': 0.00}

    def __init__(self, elems, coords, level, vdwfactor):
        self.elems = elems
        self.coords = coords
        self.level = level
        self.vdwfactor = vdwfactor
        vdws = []
        for elem in elems:
            vdws.append(self.elem2vdw[elem])
        self.vdws = vdws
        self.create_surface()

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = [vertice[0] for vertice in self.vertices]
        y = [vertice[1] for vertice in self.vertices]
        z = [vertice[2] for vertice in self.vertices]
        ax.scatter(x, y, z, s=1, color='red')
        plt.show()

    def create_surface(self):
        spheres = []
        for coord, radius in zip(self.coords, self.vdws):
            sphere = Sphere(coord, self.vdwfactor*radius)
            sphere.add_points(self.level)
            spheres.append(sphere)
        self.spheres = spheres
        for sphere in self.spheres:
            self.remove_overlap(sphere)
        areas = []
        vertices = []
        for sphere in self.spheres:
            areas.append(sphere.area)
            vertices.extend(sphere.vertices)
        self.areas = areas
        self.vertices = vertices

    def move(self, other, center):
        """Move sphere to center"""
        for index, vertice in enumerate(other.vertices):
            other.vertices.pop(index)
            other.vertices.insert(index, (center[0] - vertice[0],
                                          center[1] - vertice[1],
                                          center[2] - vertice[2]))

    def inside(self, point, center, radius):
        """Return True if point is inside sphere at center with radius"""
        r2 = ((point[0] - center[0])**2 +
              (point[1] - center[1])**2 +
              (point[2] - center[2])**2)
        return r2 <= radius**2 + 0.1 * radius

    def remove_overlap(self, sphere):
        """Remove vertices from spheres that are overlapping with other
        spheres"""
        for sph in self.spheres:
            if sph.center == sphere.center:
                continue
            vertices = []
            for index, vertice in enumerate(sphere.vertices):
                if self.inside(vertice, sph.center, sph.radius):
                    continue
                vertices.append(vertice)
            sphere.vertices = list(set(vertices))


elems = ['C', 'C', 'C', 'C', 'C', 'O', 'O', 'N', 'N', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']

coords=[[-3.323000, -1.195000,  0.000000],
        [-2.351000, -0.029000,  0.000000],
        [ 0.000000,  0.661000,  0.000000],
        [ 1.369000, -0.024000,  0.000000],
        [ 3.807000,  0.313000,  0.000000],
        [-2.732000,  1.143000,  0.000000],
        [ 1.486000, -1.250000,  0.000000],
        [-1.028000, -0.352000,  0.000000],
        [ 2.437000,  0.813000,  0.000000],
        [-2.833000, -2.174000,  0.000000],
        [-3.967000, -1.118000,  0.883000],
        [-3.967000, -1.118000, -0.883000],
        [-0.700000, -1.312000,  0.000000],
        [-0.099000,  1.311000, -0.880000],
        [-0.099000,  1.311000,  0.880000],
        [ 2.278000,  1.811000,  0.000000],
        [ 4.488000,  1.166000,  0.000000],
        [ 3.995000, -0.300000,  0.887000],
        [ 3.995000, -0.300000, -0.887000]]

mol = MolecularSurface(elems, coords, 2, 1.0)
mep = open('mep.inp', 'w')
mep.write('{}\n'.format(len(mol.vertices)))
mep.write('AA\n')
coords = ''
for vertice, area in zip(mol.vertices, mol.areas):
    coords += '{0[0]:12.6f} {0[1]:12.6f} {0[2]:12.6f} {1:12.6f}\n'.format(vertice, area)
mep.write(coords)

