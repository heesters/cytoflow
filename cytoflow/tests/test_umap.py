#!/usr/bin/env python3.8
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2025
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Tests for cytoflow.operations.umap
'''

import unittest, sys
from unittest import mock

import pandas as pd

import cytoflow as flow
import cytoflow.utility as util
from .test_base import ImportedDataTest


class TestUMAP(ImportedDataTest):

    def setUp(self):
        # ImportedDataTest thins the data set to a few hundred events, which
        # keeps these tests fast.
        super().setUp()
        self.op = flow.UMAPOp(name = "UMAP",
                              channels = ["V2-A", "V2-H", "Y2-A", "Y2-H"],
                              scale = {"V2-A" : "log",
                                       "V2-H" : "log",
                                       "Y2-A" : "log",
                                       "Y2-H" : "log"},
                              sample = 0.5,
                              random_state = 0)

    def testEstimate(self):
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)
        self.assertIn("UMAP_1", ex2.channels)
        self.assertIn("UMAP_2", ex2.channels)
        self.assertNotIn("UMAP_3", ex2.channels)

        self.assertIsInstance(ex2.data.index, pd.RangeIndex)
        self.assertFalse(ex2.data[["UMAP_1", "UMAP_2"]].isna().any().any())

    def testEstimateBy(self):
        self.op.by = ["Dox"]
        self.op.estimate(self.ex)
        self.assertEqual(len(self.op._umap), 3)

        ex2 = self.op.apply(self.ex)

        self.assertIn("UMAP_1", ex2.channels)
        self.assertIn("UMAP_2", ex2.channels)

    def testEstimateSubset(self):
        self.op.estimate(self.ex, subset = "Dox == 10.0")

        ex2 = self.op.apply(self.ex)

        self.assertIn("UMAP_1", ex2.channels)
        self.assertIn("UMAP_2", ex2.channels)

    def testComponents(self):
        self.op.n_components = 3
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)
        self.assertIn("UMAP_1", ex2.channels)
        self.assertIn("UMAP_2", ex2.channels)
        self.assertIn("UMAP_3", ex2.channels)

    def testParameters(self):
        self.op.n_neighbors = 5
        self.op.min_dist = 0.0
        self.op.metric = "cosine"
        self.op.estimate(self.ex)
        model = self.op._umap[True]
        self.assertEqual(model.n_neighbors, 5)
        self.assertEqual(model.min_dist, 0.0)
        self.assertEqual(model.metric, "cosine")
        self.op.apply(self.ex)

    def testReproducible(self):
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        op2 = self.op.clone_traits()
        op2.estimate(self.ex)
        ex3 = op2.apply(self.ex)

        pd.testing.assert_frame_equal(ex2.data, ex3.data)

    def testApplyWithoutEstimate(self):
        with self.assertRaises(util.CytoflowOpError):
            self.op.apply(self.ex)

    def testBadNeighbors(self):
        self.op.n_neighbors = 1
        with self.assertRaises(util.CytoflowOpError):
            self.op.estimate(self.ex)

    def testTooFewEvents(self):
        self.op.sample = 0.01
        with self.assertRaises(util.CytoflowOpError):
            self.op.estimate(self.ex)

    def testNameCollision(self):
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self.op.apply(ex2)

    def testMissingUMAP(self):
        # simulate umap-learn not being installed
        with mock.patch.dict(sys.modules, {'umap' : None}):
            with self.assertRaisesRegex(util.CytoflowOpError, "umap-learn"):
                self.op.estimate(self.ex)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
