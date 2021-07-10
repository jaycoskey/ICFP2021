#!/usr/bin/env python

import unittest

from util import test_is_inside_sm_parallel
from poses import test_init_tools


class TestPose(unittest.TestCase):
    def test_util(self):
        test_is_inside_sm_parallel()

    def test_poseprob(self):
        test_init_tools()


if __name__ == '__main__':
    unittest.main()
