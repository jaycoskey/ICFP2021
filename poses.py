#!/usr/bin/env python

from dataclasses import dataclass
import os

import json
from matplotlib import pyplot as plt
import networkx as nx
from networkx.algorithms.components import connected_components
import numpy as np

from util import is_inside_sm, do_edges_cross, geom_vec_dist_polygon_vertex, verts_to_closed_polygon_edges


@dataclass
class PoseTools:
    graph = None
    cut_points = None
    cut_components = None
    vertnum2is_in_hole = None
    vertnum2proximity = None

    def __init__(self, prob):
        self.graph = nx.convert.from_edgelist(prob.fig_edges)
        assert(self.graph is not None)
        self.cut_points = list(nx.articulation_points(self.graph))
        assert(self.cut_points is not None)
        cp2nbrs = {cp:set(self.graph.neighbors(cp))
                      for cp in self.cut_points
                      }

        cut_graph = self.graph.copy()
        cut_graph.remove_nodes_from(self.cut_points)
        self.cut_components = connected_components(cut_graph)

        # Convert cut_components from generate to list, allowing the cc's to be individually referenceable.
        self.cut_components = [cc for cc in self.cut_components]

        # Add each cutpoint to each of its neighboring cut components
        for cp in self.cut_points:
            for cc in self.cut_components:
                if any((node in cp2nbrs[cp] for node in cc)):
                    # print(f'INFO: Adding cut point {cp} to a cut component')
                    cc.add(cp)

        self.vertnum2is_in_hole = {k:is_inside_sm(prob.hole, prob.fig_verts[k])
                                      for k in range(len(prob.fig_verts))
                                      }
        self.vertnum2proximity = {k:prob.hole_proximity_vert(prob.fig_verts[k])
                                  for k in range(len(prob.fig_verts))
                                  }
        # self.edge2proximity = {e:self.hole_proximity_edge(e) for e in self.graph.edges()}


# TODO: Keep attempt iterations distinct from solution?
class PoseProb:
    problem_dir = './problems'

    @classmethod
    def download(id):
        raise NotImplementedError

    def __init__(self, id):
        """Loads problem from local file.
        Either solve() or read_soln() is called later.
        """
        self.id = id

        prob = self._get_prob(id)
        self.hole = np.array(prob[0])
        self.fig_verts = np.array(prob[1])
        self.fig_edges = np.array(prob[2])
        self.epsilon = int(prob[3])

        # Populated by init_tools()
        self.tools = None

        # Populated by solve or read_soln()
        # soln_score can be used to determine whether soln_verts is actually valid---not just an attempt
        self.soln_verts = None
        self.soln_score = None

    def __repr__(self):
        pass

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
        in_path = os.path.join(PoseProb.problem_dir, str(id) + '.problem')
        with open(in_path, 'r') as f:
            problem = json.load(f)
            hole = problem['hole']
            fig_verts=problem['figure']['vertices']
            fig_edges=problem['figure']['edges']
            epsilon=problem['epsilon']
            return hole, fig_verts, fig_edges, epsilon

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

        if self.soln_verts is not None:
            soln_x = self.soln_verts[:,0].flatten()
            soln_y = self.soln_verts[:,1].flatten()
            soln_y = vflip(soln_y)
            plt.plot(soln_x[self.fig_edges.T], soln_y[self.fig_edges.T],
                    linestyle='-', color='green', markerfacecolor='green', marker='o')

        plt.show()

    def hole_proximity_vert(self, node):
        return geom_vec_dist_polygon_vertex(self.hole, node)

    def init_tools(self):
        """Initializes info helpful in solving problem.
        Example: Cut points in graph, around which components can be freely rotated.
        """
        self.tools = PoseTools(prob=self)

    def is_edgelist_in_hole(self, verts):
        # Step 1: Test whether verts[0] lies inside hole
        if not is_inside_sm(self.hole, verts[0]):
            return False

        # Step 2: Test whether any edges cross hole boundary
        hole_edges = verts_to_closed_polygon_edges(self.hole)
        fig_edges = [(self.fig_verts[a], self.fig_verts[b]) for a,b in self.fig_edges]
        for hole_edge in hole_edges:
            for fig_edge in fig_edges:
                if do_edges_cross(hole_edge, fig_edge):
                    return False

        return True

    # TODO: Test cases: figure {vertex,edge} {intersects,overlaps} hole {vertex,edge}
    def is_soln(self, verts):
        done = self.is_edgelist_in_hole(verts)
        if done:
            return True
        # TODO: Add test for epsilon distortion

    def read_soln(self, id):
        """Reads solution from local JSON file
        Allows for multiple solutions per problem
        """
        pass

    def solve(self):
        # Is the original problem already solved?
        if self.is_soln(self.fig_verts):
            self.soln_verts = self.fig_verts.copy()
            return

        # Temporary question: Is there a cut component that, with the cut points, lives entirely within the hole?
        # TODO: Also check edges.
        for cc in self.tools.cut_components:
            if all(is_inside_sm(self.hole, self.fig_verts[node_i]) for node_i in cc):
                print(f'INFO: {self.id}: Found cut component with verts entirely within hole')

        # print('INFO: Failed to find solution')

    def score(self, verts):
        pass


    def write_soln(self):
        pass


def test_init_tools():
    prob3 = PoseProb(3)
    assert(len(prob3.fig_verts) == 36)
    prob3.init_tools()
    assert(len(prob3.tools.graph) == 36)
    assert(len(prob3.tools.cut_points) == 16)


def test_display():
    def hshift5(p):
        return np.array([p[0] + 5, p[1]])

    prob2 = PoseProb(2)

    prob2.soln_verts = prob2.fig_verts.copy()
    # np_hshift5 = np.vectorize(hshift5)  # Vectorize isn't very efficient, but that's OK.
    # prob2.soln_verts = np_hshift5(prob2.soln_verts)
    for k in range(len(prob2.soln_verts)):
        prob2.soln_verts[k][0] += 5
    prob2.display()


def test_is_soln():
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
       prob.write_soln()
