#!/usr/bin/env python
try:
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    print('Import from matplotlib failed! Plotting will fail.')
try:
    import matplotlib.pyplot as plt
except ImportError:
    print('Import from matplotlib failed! Plotting will fail.')
import numpy as np

# TODO: tuples vs. lists?

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

    elem2vdw = { 'H': 1.20, 'He': 1.40, 'Li': 2.20, 'Be': 1.90,  'B': 1.80,
                 'C': 1.70,  'N': 1.60,  'O': 1.55,  'F': 1.50, 'Ne': 1.54,
                'Na': 2.40, 'Mg': 2.20, 'Al': 2.10, 'Si': 2.10,  'P': 1.95,
                 'S': 1.80, 'Cl': 1.80, 'Ar': 1.88, 'Br': 1.90,  'X': 1.00}

    def __init__(self, elems=['X'], coords=[[0.0, 0.0, 0.0]], detail=2,
                 vdwfactor=1.0):
        self.elems = elems
        self.coords = coords
        self.detail = detail
        self.vdwfactor = vdwfactor
        vdws = []
        for elem in elems:
            vdws.append(1.17 * self.elem2vdw[elem])
#            vdws.append(2.0 * self.elem2vdw[elem])
        self.vdws = vdws
        self.create_surface()

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = [centroid[0] for centroid in self.centroids]
        y = [centroid[1] for centroid in self.centroids]
        z = [centroid[2] for centroid in self.centroids]
        ax.scatter(x, y, z, s=1, color='red')
        plt.show()

    def create_surface(self):
        spheres = []
        for coord, radius in zip(self.coords, self.vdws):
            sphere = Sphere(coord, self.vdwfactor*radius, self.detail)
            spheres.append(sphere)
        self.spheres = spheres
        self.remove_overlap()
        areas = []
        centroids = []
        for sphere in self.spheres:
            areas.extend(sphere.areas)
            centroids.extend(sphere.centroids)
        self.areas = areas
        self.centroids = centroids

    def inside(self, point, center, radius):
        """Return True if point is (almost) inside sphere"""
        r2 = ((point[0] - center[0])**2 +
              (point[1] - center[1])**2 +
              (point[2] - center[2])**2)
        return r2 <= (radius + 0.01 * radius)**2

    def remove_overlap(self):
        """Remove vertices from spheres that are overlapping with other
        spheres"""
        for sphere in self.spheres:
            for other in self.spheres:
                if other.center == sphere.center:
                    continue
                areas = []
                centroids = []
                for area, centroid in zip(sphere.areas, sphere.centroids):
                    if self.inside(centroid, other.center, other.radius):
                        continue
                    areas.append(area)
                    centroids.append(centroid)
                sphere.areas = areas
                sphere.centroids = centroids

    def write_surface(self):
        """Write input files for PE module."""
        fout = open('surface.dat', 'w')
        fout.write('NPOINTS\n')
        fout.write('{}\n'.format(len(self.areas)))
        output = ''
        fout.write('coordinates\n')
        for area, centroid in zip(self.areas, self.centroids):
            output += ('{0[0]:12.6f} {0[1]:12.6f} {0[2]:12.6f}'.format(centroid) +
                       '{0:12.6f}\n'.format(area))
        fout.write(output)
        fout.close()

if __name__ == "__main__":
    elems = ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'O', 'H', 'H', 'O', 
            'H', 'H']  
#    elems = ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H']
 
    coords=[[-0.145335, -0.546770,  0.000607],
            [ 1.274009, -0.912471, -0.000167],
            [ 1.630116, -2.207690, -0.000132],
            [-0.560104,  0.608977,  0.000534],
            [-0.871904, -1.386459,  0.001253],
            [ 2.004448, -0.101417, -0.000710],
            [ 0.879028, -3.000685,  0.000484],
            [ 2.675323, -2.516779, -0.000673],
            [-3.328551, -0.103229, -0.000415],
            [-2.503795,  0.413221,  0.000339],
            [-4.039214,  0.546729, -0.000849],
            [ 1.742297,  2.341361, -0.000745],
            [ 0.841678,  1.971807, -0.000820],
            [ 1.632558,  3.298301,  0.004154]]



#    elems = ['C', 'C', 'C', 'C', 'C', 'O', 'O', 'N', 'N', 'H', 'H', 'H',
#            'H', 'H', 'H', 'H', 'H', 'H', 'H']

#    coords=[[-3.323000, -1.195000,  0.000000],
#            [-2.351000, -0.029000,  0.000000],
#            [ 0.000000,  0.661000,  0.000000],
#            [ 1.369000, -0.024000,  0.000000],
#            [ 3.807000,  0.313000,  0.000000],
#            [-2.732000,  1.143000,  0.000000],
#            [ 1.486000, -1.250000,  0.000000],
#            [-1.028000, -0.352000,  0.000000],
#            [ 2.437000,  0.813000,  0.000000],
#            [-2.833000, -2.174000,  0.000000],
#            [-3.967000, -1.118000,  0.883000],
#            [-3.967000, -1.118000, -0.883000],
#            [-0.700000, -1.312000,  0.000000],
#            [-0.099000,  1.311000, -0.880000],
#            [-0.099000,  1.311000,  0.880000],
#            [ 2.278000,  1.811000,  0.000000],
#            [ 4.488000,  1.166000,  0.000000],
#            [ 3.995000, -0.300000,  0.887000],
#            [ 3.995000, -0.300000, -0.887000]]

    mol = MolecularSurface(elems, coords, detail=1)
#    mol = MolecularSurface(elems, coords, detail=2)
    mol.write_surface()
    mol.plot()
