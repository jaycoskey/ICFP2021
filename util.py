#!/usr/bin/env python

from math import sqrt
from time import time


import matplotlib.path as mpltPath
import matplotlib.pyplot as plt
from numba import jit, njit
import numba
import numpy as np


def do_edges_cross(e1, e2):
    a1, b1 = e1
    a2, b2 = e2
    orientation_a2 = geom_orientation(a1, b1, a2)
    orientation_b2 = geom_orientation(a1, b1, b2)
    orientation_a1 = geom_orientation(a2, b2, a1)
    orientation_b1 = geom_orientation(a2, b2, b1)

    if (orientation_a2 * orientation_b2 == -1
            and orientation_a1 * orientation_b1 == -1):
        return True

    return False


def geom_vec_dist_edge_vertex(a, b, node):
    dist2 = geom_dist2_vertex_vertex(a, b)
    if dist2 == 0.0:
        return geom_dist_vertex_vertex(a, node)

    t = geom_dot(node - a, b - a) / dist2
    if t < 0.0:
        return a - node, geom_dist_vertex_vertex(a, node)
    elif t > 1.0:
        return b - node, geom_dist_vertex_vertex(b, node)
    else:
        min_vec = a + t * (b - a) - node
        min_dist = geom_dist_vertex_vertex(min_vec, node)
        return min_vec, min_dist


def geom_vec_dist_polygon_vertex(poly, node):
    edges = verts_to_closed_polygon_edges(poly)
    e0 = edges[0]
    min_vec, min_dist = geom_vec_dist_edge_vertex(e0[0], e0[1], node)
    for e in edges[1:]:
        cur_vec, cur_dist = geom_vec_dist_edge_vertex(e0[0], e0[1], node)
        if cur_dist < min_dist:
            min_vec = cur_vec
            min_dist = cur_dist
    return min_vec, min_dist

def geom_dist_vertex_vertex(a, b):
    return sqrt(geom_dist2_vertex_vertex(a, b))


def geom_dist2_vertex_vertex(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def geom_dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def geom_norm(v):
    return sqrt(geom_dot(v, v))


def geom_orientation(a, b, c):
    """Returns 0 if points are colinear, 1 if clockwise, -1 if counter-clockwise
       Does so by computing determinant of matrix row vectors a, b, c, extended w/ 1.
    """
    det = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])
    return sgn(det)


# TODO: Take edges into account
def geom_proximity(verts, edges, node):
    """Return vector pointing from node to nearest vert."""
    min_dist = geom_distance_vertex_vertex(verts[0], node)
    min_vec  = verts[0] - node
    for v in verts:
        cur_dist = geom_distance_vertex_vertex(v, node)
        if  cur_dist < min_dist:
            min_dist = cur_dist
            min_vec = v - node
    return min_vec


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


def verts_to_closed_polygon_edges(verts):
    return (
            [(verts[k], verts[k+1]) for k in range(0, len(verts) - 1)]
            + [(verts[-1], verts[0])]
            )


def sgn(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


if __name__ == '__main__':
    test_is_inside_sm_parallel()
