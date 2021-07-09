#!/usr/bin/env python

import os

import json
import networkx as nx
import numpy as np

from util import dict_merge

# Note: See README for TODO list

# TODO
# def display():

class PoseProb:
    problem_dir = './problems'

    @classmethod
    def download(prob_num): 
        raise NotImplementedError

    def __init__(self, prob_num):
        """Loads problem from local file.
        Either solve() or read_soln() is called later.
        """
        self.prob_num = prob_num

        prob = self._get_prob(prob_num)
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
        return (f'{self.prob_num}: '
                + f'hole.shape={self.hole.shape}; '
                + f'fig_verts.shape={self.fig_verts.shape}; '
                + f'fig_edges.shape={self.fig_edges.shape}; '
                + f'epsilon={self.epsilon}'
                )

    def _get_prob(self, prob_num): 
        """Reads problem from local JSON file
        Returns hole, fig_verts, fig_edges, epsilon
        Called by constructor
        Only one problem file per prob_num
        """
        in_path = os.path.join(PoseProb.problem_dir, str(prob_num) + '.problem')
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
        # vert_count = len(self.fig_verts)
        # adj_dict = dict_merge(
        #                {a:b for (a,b) in self.fig_edges},
        #                {b:a for (a,b) in self.fig_edges}
        #                )
        # adj_matrix = [[1 if (a,b) for (a,b) in adj_dict.items())
        self.graph = nx.convert.from_edgelist(self.fig_edges)
        self.cut_points = list(nx.articulation_points(self.graph))
        assert(self.graph is not None)
        assert(self.cut_points is not None)

    def display(self):
        pass

    def read_soln(self, prob_num): 
        """Reads solution from local JSON file
        Allows for multiple solutions per prob_num
        """
        pass

    def solve(self):
        pass

    def score(self, verts):
        pass


    def write_soln(self):
        pass


def test():
    prob3 = PoseProb(3)
    assert(len(prob3.fig_verts) == 36)
    assert(prob3.cut_points is None)
    prob3.analyze()
    assert(len(prob3.graph) == 36)
    assert(len(prob3.cut_points) == 16)


if __name__ == '__main__':
    test() 
    prob_num = 1
    prob = PoseProb(prob_num)
    print(prob)
    prob.analyze()
    prob.solve()
    # prob.display()
    prob.write_soln()
