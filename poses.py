#!/usr/bin/env python

import json
import numpy as np
import os


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
        prob = self._get_prob(prob_num)
        self.hole = prob[0]
        self.fig_verts = prob[1]
        self.fig_edges = prob[2]
        self.epsilon = prob[3]

        self.soln_verts = None
        self.soln_score = None

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
        pass

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

if __name__ == '__main__':
    prob_num = 1
    prob = PoseProb(prob_num)
    prob.solve()
    # prob.display()
    prob.write_soln()
