#!/usr/bin/env python

from collections import defaultdict
import copy
from dataclasses import dataclass
import os
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
import urllib3

import geom


Proximity = NewType('Proximity', Tuple[np.ndarray, float])


# Abbreviations used
#   * cc   = connected component
#   * cp   = cut point
#   * fig  = figure
#   * rot  = rotation
#   * soln = solution
#   * vec  = vector
#   * vert = vertex


@dataclass
class PoseTools:
    graph = None

    cut_points = None
    # rot_components = None

    flip_components = None

    vertnum2is_in_hole = None
    vertnum2proximity = None

    def __init__(self, prob):
        self.graph:Graph = nx.convert.from_edgelist(prob.fig_edges)  # Nodes are vertex ordinals
        assert(self.graph is not None)

        self.cut_points:List[np.int64] = list(nx.articulation_points(self.graph))
        assert(self.cut_points is not None)
        cp2nbrs:Map[np.int64, Set[np.int64]] = {
                cp:set(self.graph.neighbors(cp))
                for cp in self.cut_points
                }

        ############################################################
        # cc_graph: Closed Component Graph & rot_components
        ############################################################
        sliced_graph = self.graph.copy()
        sliced_graph.remove_nodes_from(self.cut_points)
        self.conn_components = list(map(frozenset, connected_components(sliced_graph)))
        ccid_to_cpid:Dict[Int, Set[Int]] = defaultdict(set)
        cpid_to_ccid:Dict[Int, Set[Int]] = defaultdict(set)

        for cpid in range(len(self.cut_points)):
            for ccid in range(len(self.conn_components)):
                if any((v in cp2nbrs[self.cut_points[cpid]] for v in self.conn_components[ccid])):
                    # print(f'INFO: Adding cut point {cpid} to a rot component')
                    ccid_to_cpid[ccid].add(-cpid)
                    cpid_to_ccid[-cpid].add(ccid)

        # Make the connected component graph (from the original figure, less cut points)
        self.cc_graph = nx.Graph()  # Nodes have type="component" (connected component) or type="cutpoint" (cutpoints)
        for ccid in range(len(self.conn_components)):
            self.cc_graph.add_node(ccid, type='component')
        for cpid in range(len(self.cut_points)):
            self.cc_graph.add_node(-cpid, type='cutpoint')

        for cpid, ccids in cpid_to_ccid.items():
            for ccid in ccids:
                self.cc_graph.add_edge(cpid, ccid)

        ############################################################
        # flip_graph: Closed Component Graph & rot_components
        # TODO: Find flip components based on cut components with two cut points
        ############################################################
        self.flip_triples = []
        for e in prob.fig_edges:
            g_less_e = self.graph.copy()
            g_less_e.remove_edge(e[0], e[1])
            g_less_e.remove_node(e[0])
            g_less_e.remove_node(e[1])
            ccs = list(connected_components(g_less_e))
            if len(ccs) > 1:  # Found a cut capped edge (see README)
                self.flip_triples.append((ccs[0], e, ccs[1]))

        self.vertnum2is_in_hole:Map[int, bool] = {
                k:geom.is_inside_polygon(prob.hole, prob.fig_verts[k])
                for k in range(len(prob.fig_verts))
                }

        self.vertnum2proximity:Dict[int, Proximity] = {
                k:prob.hole_proximity_vert(prob.fig_verts[k])
                for k in range(len(prob.fig_verts))
                }

        # self.edge2proximity = {e:self.hole_proximity_edge(e) for e in self.graph.edges()}

    def rot_component_ids(self):
        for cc, attrs in self.cc_graph.nodes(data=True):
            if attrs['type'] == 'component' and self.cc_graph.degree(cc) == 1:
                yield cc


# TODO: Keep attempt iterations distinct from solution?
class PoseProb:
    PROBLEM_DIR = './problems'
    SOLUTION_DIR = './solutions'

    def __init__(self, id):
        """Loads problem from local file.
        Either solve() or read_pose() is called later.
        """
        self.id = id

        prob = self._get_prob(id)
        self.hole = np.array(prob[0])
        self.fig_verts = np.array(prob[1])
        self.fig_edges = np.array(prob[2])
        self.epsilon = int(prob[3])

        # Populated by init_tools()
        self.tools = None

        # Populated by solve or read_pose()
        self.pose_verts = None
        self.pose_is_valid = None
        self.pose_dislikes = None

    def __repr__(self):
        raise NotImplementedError

    def __str__(self):
        return (f'{self.id}: '
                + f'hole.shape={self.hole.shape}; '
                + f'fig_verts.shape={self.fig_verts.shape}; '
                + f'fig_edges.shape={self.fig_edges.shape}; '
                + f'epsilon%={100*self.epsilon/1_000_000}'
                )

    def _get_prob(self, id):
        """Reads problem from local JSON file
        Returns hole, fig_verts, fig_edges, epsilon
        Called by constructor
        Only one problem file per id
        """
        in_path = os.path.join(PoseProb.PROBLEM_DIR, str(id) + '.problem')
        with open(in_path, 'r') as f:
            problem = json.load(f)
            hole:List[List[int]]      = problem['hole']
            fig_verts:List[List[int]] = problem['figure']['vertices']
            fig_edges:List[List[int]] = problem['figure']['edges']
            epsilon:int               = problem['epsilon']

            return hole, fig_verts, fig_edges, epsilon

    def dislikes(self) -> int:
        result = 0
        for h in self.hole:
            distances = [geom.dist_vert_vert(h, v) for v in self.pose_verts]
            result += min(distances)
        return result

    def display(self):
        def vflip(x):
            return -x

        plt.title(f'Problem {self.id}')
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

        plt.show()

    def get_nudge_dist(self, verts):
        """
        Each vertex crossing the hole boundary requires a nudge (a translation) to get it to fit the hole.
        @returns the largest vector needed for any of the vertices in the collection (TODO: Exact type)
        """
        raise NotImplementedError

    def get_twist_angle(self, center, verts) -> float:
        """
        If the top of a region needs to be moved left, and the bottom needs to be moved right,
        then rotating the verts might get it to fit in the hole.
        """
        raise NotImplementedError

    def hole_proximity_vert(self, node):
        return geom.vec_dist_polygon_vert(self.hole, node)

    def init_tools(self):
        """Initializes info helpful in solving problem.
        Example: Cut points in graph, around which components can be freely rotated.
        """
        self.tools = PoseTools(prob=self)

    def is_edgelist_in_hole(self, verts) -> bool:
        # Step 1: Test whether verts[0] lies inside hole
        if not geom.is_inside_polygon(self.hole, verts[0]):
            return False

        # Step 2: Test whether any edges cross hole boundary
        hole_edges = geom.verts_to_closed_polygon_edges(self.hole)
        fig_edges = [(self.fig_verts[a], self.fig_verts[b]) for a,b in self.fig_edges]
        for hole_edge in hole_edges:
            for fig_edge in fig_edges:
                if geom.do_edges_cross(hole_edge, fig_edge):
                    return False

        return True

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

    def pose_verts_to_json(self):
        result = {'vertices': json.dumps(self.pose_verts)}
        # print(f'INFO: pose_verts_to_json: result={result}')
        return result

    def read_pose(self, id) -> None:
        """Reads solution from local JSON file
        Allows for multiple solutions per problem
        """
        pass

    def solve(self) -> bool:
        """Attempt to solve the problem, using the problem definition and data computed in self.tools
        Things to try (iteratively), not necessarily in this order:
          * Translate figure as a whole.
          * Rotate "free" cut components around their cut points.
          * Flip portions of the figure over "creases" (cut edges that do not contain a cut point).
          * Squash & stretch.
        """
        # Is the original problem already solved?
        if self.is_soln(self.fig_verts):
            self.pose_verts = self.fig_verts.copy()
            return True

        if self.tools is None:
            self.init_tools()
        assert(self.tools)

        # Temporary question: Is there a cut component that, with the cut points, lives entirely within the hole?
        # TODO: Also check edges.
        # for ccid in self.tools.rot_component_ids():
        #     if all(geom.is_inside_polygon(self.hole, self.fig_verts[node_i]):
        #         print(f'INFO: {self.id}: Found cut component with verts entirely within hole')

        # print('INFO: Failed to find solution')

        self.pose_verts = self.fig_verts.copy()

        # TODO: Iterative transforms

        return self.is_soln(self.pose_verts)

    def submit(self) -> None:
        host = 'poses.live'
        url = f'/problems/{self.id}/solutions'
        api_token = os.environ['ICFP_2021_API_TOKEN']
        headers = {
                'User-Agent': 'python',
                'Content-Type': 'application/x-www-form-urlencoded',
                'Authorization': api_token,
                }
        content = urllib3.urlencode(self.pose_verts_to_json())

        conn = httplib.HTTPSConnection(host)
        conn.request('POST', url, content, headers)
        response = conn.getresponse()
        data = response.read()
        print('Solution for #{self.id} submitted')
        print(f'Response status: {response.status}; reason={response.reason}')
        print('Data: ')
        print(data)

    def write_pose(self, fname) -> None:
        if not os.path.exists(PoseProb.SOLUTION_DIR):
            os.makedirs(PoseProb.SOLUTION_DIR)

        with open(f'solutions/{self.id}.solution', 'w') as f:
            json.dump({ 'vertices': self.pose_verts.tolist() }, f)

    def xform_flip_component(self, edge, fc) -> None:
        """Flip the specified FlipComponent (@fc) across its flip line (@param edge)"""
        for vi in fc:
            self.pose_verts[vi] = geom.reflect_across_edge(edge, self.pose_verts[vi])
            
    def xform_rotate_component(self, rc, angle) -> None:
        """Rotate the specified rot_component (@rc) clockwise about its cut point through the given angle (@angle)
        Since we're using integral coordinates, we round off to the nearest integer.
        """
        raise NotImplementedError

    def xform_stretch(self, verts, base, dir, sfactor) -> None:
        """
        Move each specified vert v in the direction dir by an amount determined by (v-base) and sfactor.
        @param sfactor:
          - Values less than 1 contract v along dir toward base
          - Values greater than 1 extend v along dir away from base
        """
        for v in verts:
            v += (1 - sfactor) * proj_vec_onto(v - base, dir)


def test_init_tools() -> None:
    prob3 = PoseProb(3)
    assert(len(prob3.fig_verts) == 36)
    prob3.init_tools()
    assert(len(prob3.tools.graph) == 36)
    assert(len(prob3.tools.cut_points) == 16)


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
    test_init_tools()
    # test_display()
    test_is_soln()
    for id in range(1, 10+1):
       prob = PoseProb(id)
       print(prob)
       prob.init_tools()
       prob.solve()
       # prob.display()
       if prob.pose_verts is not None:
           prob.write_pose(os.path.join(PoseProb.SOLUTION_DIR, f'{prob.id}.solution'))
