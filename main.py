#!/usr/bin/env python

from collections import defaultdict
import copy
from dataclasses import dataclass
from pprint import pprint
import math
import os
import sys
from typing import Any, Dict, List, Optional, Set, Tuple
from typing import NewType

import json
from matplotlib import pyplot as plt
import networkx as nx
from networkx.algorithms.components import connected_components
from networkx.classes.graph import Graph
import numpy as np
# import numpy.typing as npt

import http.client as httplib
import requests
# import urllib3
from urllib.parse import urlencode

import geometry as geom
from geometry import Proximity


# Abbreviations used
#   * cc   = connected component
#   * cp   = cut point
#   * fig  = figure
#   * rot  = rotation
#   * soln = solution
#   * vec  = vector
#   * vert = vertex


@dataclass
class Topology:
    graph = None

    cc_graph = None
    cut_points = None
    conn_components = None
    flip_triples = None

    def __init__(self, prob):
        self.graph:Graph = nx.convert.from_edgelist(prob.fig_edges)  # Nodes are vertex ordinals
        assert(self.graph is not None)

        self.cut_points:List[np.int64] = list(nx.articulation_points(self.graph))
        assert(self.cut_points is not None)
        cp2nbrs:Map[np.int64, Set[np.int64]] = {
                cp:set(self.graph.neighbors(cp))
                for cp in self.cut_points
                }

        # ============================================================ 
        # cc_graph: Closed Component Graph & rot_components
        # ============================================================ 
        sliced_graph = self.graph.copy()
        sliced_graph.remove_nodes_from(self.cut_points)
        self.conn_components = list(map(frozenset, connected_components(sliced_graph)))
        ccid_to_cpid:Dict[Int, Set[Int]] = defaultdict(set)
        cpid_to_ccid:Dict[Int, Set[Int]] = defaultdict(set)

        for cpid in range(len(self.cut_points)):
            for ccid in range(len(self.conn_components)):
                if any((v in cp2nbrs[self.cut_points[cpid]] for v in self.conn_components[ccid])):
                    # print(f'INFO: Adding cut point {cpid} to a rot component')
                    ccid_to_cpid[ccid].add(-1 * cpid)  # Use negative numbers for cut points, to avoid ID collision
                    cpid_to_ccid[-1 * cpid].add(ccid)

        # Make the connected component graph (from the original figure, less cut points)
        self.cc_graph = nx.Graph()  # Nodes have type="component" (connected component) or type="cutpoint" (cutpoints)
        for ccid in range(len(self.conn_components)):
            self.cc_graph.add_node(ccid, type='component')
        for cpid in range(len(self.cut_points)):
            self.cc_graph.add_node(-1 * cpid, type='cutpoint')

        for cpid, ccids in cpid_to_ccid.items():
            for ccid in ccids:
                self.cc_graph.add_edge(cpid, ccid)

        # ============================================================ 
        # flip_graph: Closed Component Graph & rot_components
        # TODO: Find flip components based on cut components with two cut points
        # ============================================================ 
        self.flip_triples = []
        for e in prob.fig_edges:
            g_less_e = self.graph.copy()
            g_less_e.remove_edge(e[0], e[1])
            g_less_e.remove_node(e[0])
            g_less_e.remove_node(e[1])
            ccs = list(connected_components(g_less_e))
            if len(ccs) > 1:  # Found a cut capped edge (see README)
                self.flip_triples.append((ccs[0], e, ccs[1]))

    # Note: This doesn't represent all the proper subsets of the graph that can be rotated without changing edge lengths.
    #       For example, if the figure is a graph that forms a polyline, then the pose can be rotated around every vertex,
    #       but this method would only return the end-points.
    def rot_component_ids(self):
        for cc, attrs in self.cc_graph.nodes(data=True):
            if attrs['type'] == 'component' and self.cc_graph.degree(cc) == 1:
                yield cc


# TODO: Keep attempt iterations distinct from solution?
class PoseProb:
    PROBLEM_DIR = './problems'
    SOLUTION_DIR = './solutions'

    def __init__(self, prob_id):
        """Loads problem from local file.
        Either solve() or read_pose() is called later.
        """
        self.prob_id = prob_id

        prob = self._get_prob(prob_id)
        self.hole = np.array(prob[0])
        self.fig_verts = np.array(prob[1])
        self.fig_edges = np.array(prob[2])
        self.epsilon = int(prob[3])

        # Populated by init_topology(), which is called by solve()
        self.topology = None

        # Populated by solve or read_pose()
        self.pose_verts = None  # Best vertex coords found thus far.

    def __repr__(self):
        raise NotImplementedError

    def __str__(self):
        return (f'{self.prob_id}: '
                + f'hole.shape={self.hole.shape}; '
                + f'fig_verts.shape={self.fig_verts.shape}; '
                + f'fig_edges.shape={self.fig_edges.shape}; '
                + f'epsilon%={100*self.epsilon/1_000_000}'
                )

    def _get_prob(self, prob_id):
        """Reads problem from local JSON file
        Returns hole, fig_verts, fig_edges, epsilon
        Called by constructor
        Only one problem file per id
        """
        in_path = os.path.join(PoseProb.PROBLEM_DIR, str(prob_id) + '.problem')
        with open(in_path, 'r') as f:
            problem = json.load(f)
            hole:List[List[int]]      = problem['hole']
            fig_verts:List[List[int]] = problem['figure']['vertices']
            fig_edges:List[List[int]] = problem['figure']['edges']
            epsilon:int               = problem['epsilon']

            return hole, fig_verts, fig_edges, epsilon

    def centroid(self):
        total = np.array([0, 0])
        for v in self.pose_verts:
            total += v
        return (1 / len(self.pose_verts)) * total  # Does this round down to ints?

    def dislikes(self) -> int:
        result = 0
        for h in self.hole:
            distances = [geom.dist_vert_vert(h, v) for v in self.pose_verts]
            result += min(distances)
        return result

    def display(self):
        def vflip(x):
            return -x

        plt.title(f'Problem {self.prob_id}')
        hole_x = self.hole[:,0].flatten()
        hole_y = self.hole[:,1].flatten()
        hole_y = vflip(hole_y)
        plt.fill(hole_x, hole_y, fill=False,
                linestyle='-', color='blue')

        fig_x = self.fig_verts[:,0].flatten()
        fig_y = self.fig_verts[:,1].flatten()
        fig_y = vflip(fig_y)
        plt.plot(fig_x[self.fig_edges.T], fig_y[self.fig_edges.T],
                linestyle='-', color='red', markerfacecolor='red', marker='o')

        if self.pose_verts is not None:
            pose_x = self.pose_verts[:,0].flatten()
            pose_y = self.pose_verts[:,1].flatten()
            pose_y = vflip(pose_y)
            plt.plot(pose_x[self.fig_edges.T], pose_y[self.fig_edges.T],
                    linestyle='-', color='green', markerfacecolor='green', marker='o')

        print('INFO: Displaying pose. Close matplotlib window to proceed.')
        plt.show()

 
    def get_pose_nudge(self, vids):
        """
        Sometimes a figure needs a nudge (i.e., translation) to get it to fit in the hole.
        @returns Recommended translation for pose
            In case of conflict (e.g., vectors pointing in opposing directions), return None.
        """
        proximities = []
        max_vec = None
        max_dist = None
        for vid in vids:
            # Ignore vertices that fit in the hole
            if self.is_in_hole(self.pose_verts[vid]):
                continue

            proximity = self.hole_proximity_of_vert(self.pose_verts[vid])
            proximities.append(proximity)
            if max_dist is None or proximity.dist > max_dist:
                max_vec = proximity.vec
                max_dist = proximity.dist

        if len(proximities) == 0:
            return None

        # Q: Are all the proximities within 90 degrees of max_vec?
        if all(map(lambda prox: geom.dot(prox.vec, max_vec) >= 0, proximities)):
            return Proximity(max_vec, max_dist)
        else:
            return None

    def get_pose_twist(self, pose_verts, center, alpha=1.0):  # -> Optional[float]:
        """
        Sometimes a figure needs a slight twist (i.e., rotation) to get it to fit in the hole.
        If the top of a region needs to be moved left (as indicated by the Promitiy field) ,
            and the bottom needs to be moved right, then rotating the verts might get it to fit in the hole.
            Compare with how torque is computed with moments.

        This is intended to be applied to the entire pose.
        A separate method will be used for rotation component.

        @result Angle in radians
            In case of conflict (e.g., opposing rotations), return None
        """
        angles = []
        for v in pose_verts:
            # Ignore vertices that fit in the hole
            if self.is_in_hole(v):
                continue

            # TODO: Upgrade to proximity_to_polygon
            proximity = geom.proximity_to_verts(pose_verts, v)

            # Find angle from Law of Cosines
            a_vec = v - center
            c_vec = proximity.vec
            b_vec = a_vec + c_vec
            a = geom.norm(a_vec)
            b = geom.norm(b_vec)
            c = geom.norm(c_vec)
            angle = math.acos((a**2 + b**2 - c**2) / (2 * a * b))
            angles.append(angle)

        if len(angles) == 0:
            return None

        if all(map(lambda x: x >= 0, angles)):
            return alpha * max(angles)
        elif all(map(lambda x: x <= 0, angles)):
            return alpha * min(angles)
        else:
            return None

    def get_rot_component_rotation(self, ccid):
        """
        Suggest how much to rotate a rot component.
        Find all verts of the rot component outside the hole.
        @returns Recommended angle
          * If they are all within 90 degrees of the average (from the cut point),
            then return 180 degrees. (TODO: Improve algorithm)
          * Otherwise, return None
        """
        raise NotImplementedError

    def golf_score(self, pose_verts) -> Tuple[int, float]:
        """Returns (0, 0) for poses that fit in the hole.
        Otherwise, returns (count, dist),
            where count is the number of pose vertices that fall outside the hole,
            and dist is the largest proximity distance over pose vertices outisde the hole.
        """
        outside_count = 0
        max_outside_dist = 0.0
        for v in pose_verts:
            if self.is_in_hole(v):
                continue

            outside_count += 1
            prox = self.hole_proximity_of_vert(v)
            if prox.dist > max_outside_dist:
                max_outside_dist = prox.dist
        return (outside_count, max_outside_dist)

    def is_rot_component_in_hole(self, ccid):
        verts = map(lambda vid: self.pose_verts[vid], self.conn_components[ccid])
        return all(lambda v: self.is_in_hole(v), verts)

    def hole_proximity_of_vert(self, v):
        return geom.proximity_to_polygon(self.hole, v)

    def init_topology(self):
        """Initializes info helpful in solving problem.
        Example: Cut points in graph, around which components can be freely rotated.
        """
        self.topology = Topology(prob=self)

    def is_edgelist_in_hole(self, verts) -> bool:
        # Step 1: Test whether verts[0] lies in the hole
        if not self.is_in_hole(verts[0]):
            return False

        # Step 2: Test whether any edges cross hole boundary
        hole_edges = geom.verts_to_closed_polygon_edges(self.hole)
        fig_edges = [(self.fig_verts[a], self.fig_verts[b]) for a,b in self.fig_edges]
        for hole_edge in hole_edges:
            for fig_edge in fig_edges:
                if geom.do_edges_cross(hole_edge, fig_edge):
                    return False

        return True

    def is_in_hole(self, v):
        return geom.is_in_polygon(self.hole, v)

    # TODO: Test cases: figure {vertex,edge} {intersects,overlaps} hole {vertex,edge}
    def is_soln(self, verts) -> bool:
        """Solution criteria (from spec PDF):
        (a) Conn: All graph manipulations preserve vertex connectedness
        (b) Eps:  Edges can only be compressed or stretched within epsilon threshold: abs((d' / d) - 1) <= eps/1_000_000
        (c) Fit:  Every point of pose must lie on or in hole boundary
        """
        # (a) Connectedness: For now, assume we haven't goofed by breaking connectedness
        #     TODO (optional): Add sanity check

        # (b) Epsilon
        for ei in range(len(self.fig_edges)):
            old_dist = geom.dist_vert_vert(
                            self.fig_verts[self.fig_edges[ei][0]],
                            self.fig_verts[self.fig_edges[ei][1]]
                            )
            new_dist = geom.dist_vert_vert(
                            verts[self.fig_edges[ei][0]],
                            verts[self.fig_edges[ei][1]]
                            )
            distortion = abs((new_dist / old_dist) - 1)
            if distortion > (self.epsilon / 1_000_000):
                # print(f'INFO: is_soln: Failed Epsilon criterion')
                return False

        # (c) Fit
        if not self.is_edgelist_in_hole(verts):
            # print(f'INFO: is_soln: Failed Fit criterion')
            return False

        return True

    def read_pose(self, prob_id) -> None:
        """Reads solution from local JSON file
        Allows for multiple solutions per problem
        """
        pass

    def shrink_pose_verts(self):
        shrink_factor = (1_000_000 - (0.95 * self.epsilon)) / 1_000_000
        shrink_test_verts = self.pose_verts.copy()
        self.xform_shrink(shrink_test_verts, shrink_factor)
        self.pose_verts = shrink_test_verts

    def solve(self, max_iters=10) -> bool:
        """Attempt to solve the problem, using the problem definition and data computed in self.topology
        Things to try (iteratively), not necessarily in this order:
          * Translate figure as a whole.
          * Rotate "free" cut components around their cut points.
          * Flip portions of the figure over "creases" (cut edges that do not contain a cut point).
          * Squash & stretch.
        @param max_iters: Number of times to iterate through transformation attempts
        """
        if self.pose_verts == None:
            self.pose_verts = self.fig_verts.copy()

        # Is the original problem already solved?
        if self.is_soln(self.pose_verts):
            return True

        if self.topology is None:
            self.init_topology()
        assert(self.topology)

        cur_golf_score = self.golf_score(self.pose_verts)
        # ============================================================ 
        # Iterative transforms
        # TODO: Refactor to reduce repetition
        # ============================================================ 
        did_find_solution = True

        for iter_count in range(1, max_iters + 1):
            ########## Rotate rot components ########## TODO: Refactor to reduce repetition
            for ccid in self.topology.rot_component_ids():
                angle_test_verts = self.pose_verts.copy()
                for angle in range(30, 360, 30):
                    rads = math.radians(angle)
                    angle_test_verts = self.pose_verts.copy()
                    self.xform_rotate_component(angle_test_verts, ccid, rads)
                    angle_golf_score = self.golf_score(angle_test_verts)
                    if angle_golf_score < cur_golf_score:
                        og = f'cur_golf_score'
                        ng = f'angle_golf_score'
                        # print(f'INFO: Rot component rotation improved golf score from {og} to {ng}')
                        print('.', end='')
                        cur_golf_score = angle_golf_score
                        self.pose_verts = angle_test_verts
                        if cur_golf_score == (0, 0.0):
                            return True

            # ========== Flip flip components ========== TODO: Refactor to reduce repetition
            for flip_trip in self.topology.flip_triples:
                left_comp, edge, right_comp = flip_trip

                left_test_verts = self.pose_verts.copy()
                self.xform_flip_component(left_test_verts, edge, left_comp)
                left_golf_score = self.golf_score(left_test_verts)
                if left_golf_score < cur_golf_score:
                    og = f'cur_golf_score'
                    ng = f'left_golf_score'
                    # print(f'INFO: Flip improved golf score from {og} to {ng}')
                    print('.', end='')
                    cur_golf_score = left_golf_score
                    self.pose_verts = left_test_verts
                    if cur_golf_score == (0, 0.0):
                        return True

                right_test_verts = self.pose_verts.copy()
                self.xform_flip_component(right_test_verts, edge, right_comp)
                right_golf_score = self.golf_score(right_test_verts)
                if right_golf_score < cur_golf_score:
                    og = f'cur_golf_score'
                    ng = f'right_golf_score'
                    # print(f'INFO: Flip improved golf score from {og} to {ng}')
                    print('.', end='')
                    cur_golf_score = right_golf_score
                    self.pose_verts = right_test_verts
                    if cur_golf_score == (0, 0.0):
                        return True

            # ========== Nudge pose ========== TODO: Refactor to reduce repetition
            nudge_test_verts = self.pose_verts.copy()
            nudge = self.get_pose_nudge(range(len(self.pose_verts)))
            if nudge is not None:
                self.xform_translate(nudge_test_verts, nudge.vec)
                nudge_golf_score = self.golf_score(nudge_test_verts)
                if nudge_golf_score < cur_golf_score:
                    og = f'cur_golf_score'
                    ng = f'nudge_golf_score'
                    # print(f'INFO: Nudge improved golf score from {og} to {ng}')
                    print('.', end='')
                    cur_golf_score = nudge_golf_score
                    self.pose_verts = nudge_test_verts
                    if cur_golf_score == (0, 0.0):
                        return True

            # ========== Twist pose ========== TODO: Refactor to reduce repetition
            twist_test_verts = self.pose_verts.copy()

            angle = self.get_pose_twist(twist_test_verts, self.centroid())
            if angle is not None:
                self.xform_twist(twist_test_verts, self.centroid(), angle)
                twist_golf_score = self.golf_score(twist_test_verts)
                if twist_golf_score < cur_golf_score:
                    og = f'cur_golf_score'
                    ng = f'nudge_golf_score'
                    # print(f'INFO: Nudge improved golf score from {og} to {ng}')
                    print('.', end='')
                    cur_golf_score = twist_golf_score
                    self.pose_verts = twist_test_verts
                    if cur_golf_score == (0, 0.0):
                        return True

            # ========== TODO: Stretch ========== 

        return False

    def submit(self, do_submit=True) -> None:
        print(f'INFO: Entering submit for problem #{self.prob_id}')
        host = 'poses.live'
        path = f'/api/problems/{self.prob_id}/solutions'
        url = f'https://{host}{path}'
        api_token = os.environ['ICFP_2021_API_TOKEN']
        headers = {
                'User-Agent': 'python',
                'Content-Type': 'application/x-www-form-urlencoded',
                'Authorization': 'Bearer ' + api_token
                }
        vertices = [[int(coord) for coord in point] for point in self.pose_verts]
        content = { "vertices": vertices }
        print(f'Unencoded content={content}')
        payload = urlencode(content)

        net_pkg = 'httplib'
        # net_pkg = 'requests'

        # Add 1 & 2 to identifiers to satisfy mypy.
        if do_submit:
            print(f'Submitting solution for #{self.prob_id} ... ', end='')
            if net_pkg == 'httplib':
                conn1 = httplib.HTTPSConnection(host)
                conn1.request('POST', path, payload, headers)
                response1 = conn1.getresponse()
                data1 = response1.read()
            elif net_pkg == 'requests':
                response2 = requests.post(f'https://{host}{path}', data=payload, json=headers)
            print('Done!')
    
            if net_pkg == 'httplib':
                print(f'\tResponse status: {response1.status} ({response1.reason})')
                print(f'\tResponse data  : {str(data1)}')
            elif net_pkg == 'requests':
                print(f'\tResponse status: {response2.status_code} ({response2.reason})')
        else:
            print(f'INFO: Skipping submission')

    def write_pose(self, fname) -> None:
        if not os.path.exists(PoseProb.SOLUTION_DIR):
            os.makedirs(PoseProb.SOLUTION_DIR)

        with open(f'solutions/{self.prob_id}.solution', 'w') as f:
            json.dump({ 'vertices': self.pose_verts.tolist() }, f)

    def xform_flip_component(self, pose_verts, edge, fc) -> None:
        """Flip the specified FlipComponent (@fc) across its flip line (@param edge)"""
        for vi in fc:
            edge2 = (pose_verts[edge[0]], pose_verts[edge[1]])
            pose_verts[vi] = geom.reflect_across_edge(edge2, pose_verts[vi])

    def xform_rotate_component(self, pose_verts, ccid, radians) -> None:
        """Rotate the specified rot_component (@rc) clockwise about its cut point through the given angle (@angle)
        Since we're using integral coordinates, we round off to the nearest integer.
        """
        ccg = self.topology.cc_graph
        rot_component = ccg[ccid]
        rot_component_nbrs = ccg.neighbors(ccid)
        cpid = list(rot_component_nbrs)[0]  # In cc_graph, rot components have exactly one neighbor: a cut point
        cp = pose_verts[cpid]
        for vid in self.topology.conn_components[ccid]:
            v = pose_verts[vid]
            v = geom.rotate_around_point(cp, radians, v)

    def xform_shrink(self, pose_verts, sfactor):
        """Axial shrinkabe, centered at the pose's centroid()"""
        for v in pose_verts:
            v = scale_around_point(self.centroid(), sfactor, v)

    def xform_stretch(self, verts, base, dir, sfactor) -> None:
        """
        Axial stretching: Move each vert v in the direction dir by an amount determined by (v-base) and sfactor.
        @param sfactor:
          - Values less than 1 contract v along dir toward base
          - Values greater than 1 extend v along dir away from base
        """
        for v in verts:
            v += (1 - sfactor) * geom.proj_vec_onto(v - base, dir)

    def xform_translate(self, pose_verts, translation_vec) -> None:
        for v in pose_verts:
            # TODO: Audit choices of types to make lines like the following concise
            v = [v[0] + translation_vec[0], v[1] + translation_vec[1]]

    def xform_twist(self, pose_verts, center, radians) -> None:
        for v in pose_verts:
            v = geom.rotate_around_point(center, radians, v)


def test_init_topology() -> None:
    prob3 = PoseProb(3)
    assert(len(prob3.fig_verts) == 36)
    prob3.init_topology()
    assert(len(prob3.topology.graph) == 36)
    assert(len(prob3.topology.cut_points) == 16)


def test_display() -> None:
    def hshift5(p):
        return np.array([p[0] + 5, p[1]])

    prob2 = PoseProb(2)

    prob2.pose_verts = prob2.fig_verts.copy()
    # np_hshift5 = np.vectorize(hshift5)  # Vectorize isn't very efficient, but that's OK.
    # prob2.pose_verts = np_hshift5(prob2.pose_verts)
    for k in range(len(prob2.pose_verts)):
        prob2.pose_verts[k][0] += 5
    prob2.display()


def test_is_soln() -> None:
    prob6mod = PoseProb(6)

    assert(not prob6mod.is_soln(prob6mod.fig_verts))  # Original figure is not a solution
    # Remove the dent from the hole, chaning if from pac-man-like to roundish
    prob6mod.hole[-2][0] = prob6mod.hole[-1][0]
    # The original is a solution for the modified hole
    assert(prob6mod.is_soln(prob6mod.fig_verts))


if __name__ == '__main__':
    test_init_topology()
    # test_display()
    test_is_soln()
    first = 1
    last = 132  # As of 48 hours into contest, 132 available
    prob_ids:List[int] = list(range(first, last + 1))
    excluded_prob_ids = [58]  # TODO: Fix math.acos() range error arising with figure #58
    shrink_prob_ids:List[int] = []  # Attempts that might work after contraction toward centroid

    prob_ids = [95, 118, 119, 125]
    shrink_prob_ids = [66, 68, 74, 75, 89, 90, 102, 105, 106, 107, 110, 121, 123, 128]
    # Note: 129 would be solved by one more xform---a very slight vertical squash

    for prob_id in prob_ids:
        if prob_id in excluded_prob_ids:
            continue
        prob = PoseProb(prob_id)
        print(prob)
        if prob.solve():
            print('INFO: Found a solution!')
            if prob_id in shrink_prob_ids:
                prob.shrink_pose_verts()
            prob.display()
            prob.write_pose(os.path.join(PoseProb.SOLUTION_DIR, f'{prob.prob_id}.solution'))
            prob.submit()
