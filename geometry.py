#!/usr/bin/env python

from collections import namedtuple
from math import cos, sin, sqrt
from time import time


import matplotlib.path as mpltPath
import matplotlib.pyplot as plt
from numba import jit, njit
import numba
import numpy as np


# Proximity = NewType('Proximity', Tuple[np.ndarray, float])
Proximity = namedtuple('Proximity', ['vec', 'dist'])  # np.ndarray, float


def dist_vert_vert(a, b):
    return sqrt(dist2_vert_vert(a, b))


def dist2_vert_vert(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def do_edges_cross(e1, e2):
    a1, b1 = e1
    a2, b2 = e2
    orientation_a2 = orientation(a1, b1, a2)
    orientation_b2 = orientation(a1, b1, b2)
    orientation_a1 = orientation(a2, b2, a1)
    orientation_b1 = orientation(a2, b2, b1)

    if (orientation_a2 * orientation_b2 == -1
            and orientation_a1 * orientation_b1 == -1):
        return True

    return False


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def edge_vec(e):
    """@returns: Vector pointing from e[0] to e[1]"""
    return [e[1][0] - e[0][0], e[1][1] - e[0][1]]


def is_in_polygon(polygon, point):
    return is_inside_sm(polygon, point)


def norm(v):
    return sqrt(norm2(v))


def norm2(v):
    return dot(v, v)


def orientation(a, b, c):
    """Returns 0 if points are colinear, 1 if clockwise, -1 if counter-clockwise
       Does so by computing determinant of matrix row vectors a, b, c, extended w/ 1.
    """
    def sgn(x):
        return 1 if x > 0 else (-1 if x < 0 else 0)

    det = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])
    return sgn(det)


def proj_vec_onto(v, dir):
    return vec_smult(dir, dot(v, dir) / norm2(dir))


# TODO: Take edges into account
def proximity_to_verts(verts, node):
    """Return vector pointing from node to nearest vert."""
    min_dist = dist_vert_vert(verts[0], node)
    min_vec  = verts[0] - node
    for v in verts:
        cur_dist = dist_vert_vert(v, node)
        if  cur_dist < min_dist:
            min_dist = cur_dist
            min_vec = v - node
    return Proximity(min_vec, min_dist)


def proximity_to_edge(a, b, node):
    dist2 = dist2_vert_vert(a, b)
    if dist2 == 0.0:
        return dist_vert_vert(a, node)

    t = dot(node - a, b - a) / dist2
    if t < 0.0:
        return a - node, dist_vert_vert(a, node)
    elif t > 1.0:
        return b - node, dist_vert_vert(b, node)
    else:
        min_vec = a + t * (b - a) - node
        min_dist = dist_vert_vert(min_vec, node)
        return Proximity(min_vec, min_dist)


def proximity_to_polygon(poly, node):
    edges = verts_to_closed_polygon_edges(poly)
    e0 = edges[0]
    min_vec, min_dist = proximity_to_edge(e0[0], e0[1], node)
    for e in edges[1:]:
        cur_vec, cur_dist = proximity_to_edge(e0[0], e0[1], node)
        if cur_dist < min_dist:
            min_vec = cur_vec
            min_dist = cur_dist
    return Proximity(min_vec, min_dist)


def reflect_across_edge(e, vert):
    parallel = edge_vec(e)
    normal = rotate90(parallel)
    delta = vert - e[0]
    projected = proj_vec_onto(delta, normal)
    result = vert - vec_smult(projected, 2)
    return result


def rotate_around_point(p, radians, v):
    c = cos(radians)
    s = sin(radians)
    return [c * v[0] - s * v[1], s * v[0] + c * v[1]]


def rotate90(vert):
    if type(vert).__name__ == 'tuple':
        raise Exception('Look at backtrace; change to list')
    result = [-vert[1], vert[0]]
    return result


def scale_around_point(self, center, sfactor, v):
    return center + sfactor * (v - center)


def vec_smult(vec, scalar):
    return [vec[0] * scalar, vec[1] * scalar]


def verts_to_closed_polygon_edges(verts):
    return (
            [(verts[k], verts[k+1]) for k in range(0, len(verts) - 1)]
            + [(verts[-1], verts[0])]
            )


# ----------------------------------------
# From https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
#   - No license file. Will remove from this repo upon requrest.
# Found on stackoverflow at https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
#  - Posted by @Mehdi
#
@jit(nopython=True)
def is_inside_sm(polygon, point):
    length = len(polygon)-1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if point[0] > F: # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif point[0] == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (point[0]==polygon[jj][0] or (dy==0 and (point[0]-polygon[ii][0])*(point[0]-polygon[jj][0])<=0)):
                return 2

        ii = jj
        jj += 1

    #print 'intersections =', intersections
    return intersections & 1


# From https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
# Found on stackoverflow at https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
#  - Posted by @Mehdi
#
@njit(parallel=True)
def is_inside_sm_parallel(points, polygon):
    ln = len(points)
    D = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        D[i] = is_inside_sm(polygon,points[i])
    return D


# Found on stackoverflow at https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
#  - Posted by @Mehdi
def test_is_inside_sm_parallel():
    np.random.seed(2)

    time_parallelpointinpolygon=[]
    time_mpltPath=[]
    time_ray_tracing_numpy_numba=[]
    time_is_inside_sm_parallel=[]
    time_is_inside_postgis_parallel=[]
    n_points=[]

    for i in range(1, 10000002, 1000000):
        n_points.append(i)

        lenpoly = 100
        polygon = [[np.sin(x)+0.5,np.cos(x)+0.5] for x in np.linspace(0,2*np.pi,lenpoly)]
        polygon = np.array(polygon)
        N = i
        points = np.random.uniform(-1.5, 1.5, size=(N, 2))

        start_time = time()
        inside4=is_inside_sm_parallel(points,polygon)
        time_is_inside_sm_parallel.append(time()-start_time)

    plt.plot(n_points,time_is_inside_sm_parallel,label='is_inside_sm_parallel')
    plt.xlabel("N points")
    plt.ylabel("time (sec)")
    plt.legend(loc = 'best')
    plt.show()


if __name__ == '__main__':
    test_is_inside_sm_parallel()
