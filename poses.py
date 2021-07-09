#!/usr/bin/env python

import os

import json
from matplotlib import pyplot as plt
import networkx as nx
import numpy as np


# Note: See README for TODO list

# TODO
# def display():

POINT_AT_INFINITY = (1_000_000, 1_000_000)


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

        # Populated by analyze()
        self.graph = None
        self.cut_points = None

        # Populated by solve or read_soln()
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

    def analyze(self):
        """Adds info helpful to solve method
        Example: Cut points in graph, where parts can be freely rotated.
        """
        self.graph = nx.convert.from_edgelist(self.fig_edges)
        self.cut_points = list(nx.articulation_points(self.graph))
        assert(self.graph is not None)
        assert(self.cut_points is not None)

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
        # print(f'y:{type(y)}={y}')
        plt.plot(fig_x[self.fig_edges.T], fig_y[self.fig_edges.T],
                linestyle='-', color='red', markerfacecolor='red', marker='o')
        plt.show()

    # TODO: Test cases: figure {vertex,edge} {intersects,overlaps} hole {vertex,edge}
    def is_soln(self, verts):
        # Step 1: Test whether verts[0] lies inside hole
        vert0 = verts[0]

        # Step 2: Test whether any edges cross hole boundary
        edges = [(self.fig_verts[a], self.fig_verts[b]) for a,b in self.fig_edges]
        # TODO: Use is_inside_sm_parallel
        raise NotImplementedError

    def read_soln(self, id):
        """Reads solution from local JSON file
        Allows for multiple solutions per problem
        """
        pass

    def solve(self):
        pass

    def score(self, verts):
        pass


    def write_soln(self):
        pass


def test_analyze():
    prob3 = PoseProb(3)
    assert(len(prob3.fig_verts) == 36)
    assert(prob3.cut_points is None)
    prob3.analyze()
    assert(len(prob3.graph) == 36)
    assert(len(prob3.cut_points) == 16)


if __name__ == '__main__':
    test_analyze()
    for id in range(1, 10+1):
       prob = PoseProb(id)
       print(prob)
       prob.display()
       # print(f'{prob.id}: Is fig_verts a solution?: {prob.is_soln(prob.fig_verts)}')
       # prob.analyze()
       # prob.solve()
       # prob.display()
       # prob.write_soln()
