#!/usr/bin/env python3.8
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2022
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
Tests for PhenographOp (graph-based clustering)
'''
import unittest, random

import numpy as np
import pandas as pd

import cytoflow as flow
import cytoflow.utility as util
import cytoflow.operations.phenograph as pg
from cytoflow.operations.phenograph import _jaccard_graph


def make_blobs(sizes, seed = 0, tubes = None):
    """
    An experiment with one well-separated 2D blob per entry in ``sizes``.
    The channel ``Truth`` says which blob each event came from.

    If ``tubes`` is a dict of name --> function(data) returning a boolean mask,
    the events are split among the values of a ``Tube`` condition.
    """
    rng = np.random.RandomState(seed)
    centers = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0], [100.0, 100.0]])

    xy = np.vstack([rng.normal(centers[i], 1.0, size = (n, 2))
                    for i, n in enumerate(sizes)])
    truth = np.repeat(np.arange(len(sizes)), sizes)

    # shuffle, so the blobs aren't contiguous
    order = rng.permutation(len(truth))
    data = pd.DataFrame({'X' : xy[order, 0],
                         'Y' : xy[order, 1],
                         'Truth' : truth[order].astype(float)})

    ex = flow.Experiment()
    for c in ['X', 'Y', 'Truth']:
        ex.add_channel(c)

    if tubes is None:
        ex.add_events(data, {})
    else:
        ex.add_condition('Tube', 'category')
        for name, mask in tubes.items():
            ex.add_events(data[mask(data)], {'Tube' : name})

    return ex


class TestPhenograph(unittest.TestCase):

    def setUp(self):
        # blobs are made in order of decreasing size, so blob i should
        # end up as cluster i
        self.sizes = [400, 300, 200]
        self.ex = make_blobs(self.sizes)
        self.op = flow.PhenographOp(name = "PG",
                                    channels = ["X", "Y"],
                                    scale = {"X" : "linear",
                                             "Y" : "linear"},
                                    num_neighbors = 15,
                                    resolution = 0.01)

    def assertRecovered(self, ex):
        """Each blob is exactly one cluster, numbered by size"""
        np.testing.assert_array_equal(ex['PG'].astype(int).values,
                                      ex['Truth'].astype(int).values)

    def testJaccardGraph(self):
        # two triangles, far apart.  with k = 2, each point's neighbors are
        # the other two points in its triangle; two neighbors share one
        # neighbor out of three --> 1/3
        x = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.1],
                      [50.0, 50.0], [51.0, 50.0], [50.0, 51.1]])
        edges, weights, _ = _jaccard_graph(x, 2)

        self.assertEqual(sorted(map(tuple, edges.tolist())),
                         [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)])
        np.testing.assert_allclose(weights, 1.0 / 3.0)

    def testJaccardGraphChunks(self):
        # the chunked, sparse computation should agree with brute force
        x = np.random.RandomState(1).normal(size = (500, 3))
        k = 10

        chunk = pg._CHUNK
        try:
            pg._CHUNK = 64
            edges, weights, nn = _jaccard_graph(x, k)
        finally:
            pg._CHUNK = chunk

        neighbors = [set(n) for n in nn.kneighbors(return_distance = False)]
        expected = {}
        for i in range(len(x)):
            for j in neighbors[i]:
                jac = len(neighbors[i] & neighbors[j]) / len(neighbors[i] | neighbors[j])
                key = (min(i, j), max(i, j))
                expected[key] = expected.get(key, 0.0) + jac / 2.0
        expected = {e : w for e, w in expected.items() if w > 0}

        actual = {tuple(e) : w for e, w in zip(edges.tolist(), weights)}
        self.assertEqual(set(actual), set(expected))
        for e in expected:
            self.assertAlmostEqual(actual[e], expected[e])

    def testEstimate(self):
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        self.assertRecovered(ex2)
        self.assertIsInstance(ex2.data.index, pd.RangeIndex)
        self.assertEqual(ex2.data['PG'].dtype.name, 'category')

    def testLouvain(self):
        self.op.method = "louvain"
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        self.assertRecovered(ex2)

    def testEstimate1D(self):
        self.op.channels = ["X"]
        self.op.scale = {"X" : "linear"}
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        # blobs 0 and 2 have the same X, so they can't be told apart -- but
        # no cluster should have events from both them and blob 1.  (in 1D the
        # graph is a chain, which modularity is happy to cut into pieces.)
        self.assertTrue((ex2['PG'].astype(int) >= 0).all())
        is_blob_1 = (ex2['Truth'] == 1).groupby(ex2['PG'].astype(int)).nunique()
        self.assertTrue((is_blob_1 == 1).all())
        
    def testStatistic(self):
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        stat = ex2.statistics['PG']
        self.assertEqual(list(stat.index.names), ['PG'])
        self.assertEqual(list(stat.columns), ['X', 'Y', 'Proportion'])
        self.assertEqual(len(stat), 3)

        np.testing.assert_allclose(stat['Proportion'].values,
                                   [4 / 9, 3 / 9, 2 / 9])
        np.testing.assert_allclose(stat[['X', 'Y']].values,
                                   [[0, 0], [100, 0], [0, 100]],
                                   atol = 0.5)

    def testResolution(self):
        # higher resolution --> more clusters
        self.op.resolution = 5.0
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)
        self.assertGreater(len(ex2['PG'].unique()), 3)

        # and they're numbered by size
        counts = ex2.data['PG'].astype(int).value_counts()
        sizes = [counts[i] for i in range(len(counts))]
        self.assertEqual(sizes, sorted(sizes, reverse = True))

    def testDeterministic(self):
        # a resolution where the answer isn't obvious, and a sample
        self.op.resolution = 5.0
        self.op.sample = 0.5
        self.op.random_state = 7

        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        # disturb the global generators -- shouldn't matter
        random.seed(12345)
        np.random.seed(12345)

        op2 = flow.PhenographOp(name = "PG",
                                channels = ["X", "Y"],
                                scale = {"X" : "linear", "Y" : "linear"},
                                num_neighbors = 15,
                                resolution = 5.0,
                                sample = 0.5,
                                random_state = 7)
        op2.estimate(self.ex)
        ex3 = op2.apply(self.ex)

        self.assertGreater(len(ex2['PG'].unique()), 3)
        np.testing.assert_array_equal(ex2['PG'].astype(int).values,
                                      ex3['PG'].astype(int).values)

    def testGlobalRandomStateUntouched(self):
        # estimating must not consume or reset python's global generator,
        # nor leave igraph using our seeded one
        random.seed(99)
        expected = random.Random(99).random()
        self.op.estimate(self.ex)
        self.assertEqual(random.random(), expected)

    def testSample(self):
        # cluster a tenth of the events; all of them still get assigned
        self.op.sample = 0.1
        self.op.num_neighbors = 10
        self.op.estimate(self.ex)
        self.assertEqual(len(self.op._labels[True]), 90)

        ex2 = self.op.apply(self.ex)
        self.assertRecovered(ex2)

    def testApplyNewEvents(self):
        # apply to events that estimate() never saw
        self.op.estimate(self.ex)

        ex_new = make_blobs([50, 60, 70], seed = 42)
        ex2 = self.op.apply(ex_new)
        self.assertRecovered(ex2)

    def testApplyMajorityVote(self):
        # new events between the clusters go with their nearest neighbors
        self.op.estimate(self.ex)

        ex_new = flow.Experiment()
        for c in ['X', 'Y', 'Truth']:
            ex_new.add_channel(c)
        ex_new.add_events(pd.DataFrame({'X' : [40.0, 60.0, 2.0],
                                        'Y' : [0.0, 0.0, 60.0],
                                        'Truth' : [0.0, 1.0, 2.0]}), {})
        ex2 = self.op.apply(ex_new)
        self.assertEqual(list(ex2['PG'].astype(int)), [0, 1, 2])

    def testApplyGraphLabels(self):
        # with sample = 1, events that were in the graph get the graph's
        # clusters, not a vote.
        self.op.resolution = 5.0
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)
        np.testing.assert_array_equal(ex2['PG'].astype(int).values,
                                      self.op._labels[True])

    def testInvalidScale(self):
        # events that the scale can't handle are in cluster -1
        self.op.scale = {"X" : "log", "Y" : "linear"}
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        negative = (self.ex['X'] <= 0).values
        self.assertTrue(negative.any())
        self.assertTrue((ex2['PG'].astype(int).values[negative] == -1).all())
        self.assertTrue((ex2['PG'].astype(int).values[~negative] >= 0).all())

    def testFewEvents(self):
        # fewer events than neighbors
        ex = make_blobs([5, 4])
        self.op.num_neighbors = 30
        self.op.estimate(ex)
        ex2 = self.op.apply(ex)
        self.assertEqual(len(ex2), 9)
        self.assertTrue((ex2['PG'].astype(int) >= 0).all())

    def testBadParams(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.PhenographOp(name = "PG", channels = []).estimate(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            flow.PhenographOp(name = "PG", channels = ["Nope"]).estimate(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            self.op.apply(self.ex)   # no estimate

        self.op.sample = 0.001       # --> no events
        with self.assertRaises(util.CytoflowOpError):
            self.op.estimate(self.ex)

        with self.assertRaises(Exception):
            self.op.method = "kmeans"

    def testPlot(self):
        self.op.estimate(self.ex)
        self.op.default_view().plot(self.ex)

        ex2 = self.op.apply(self.ex)
        self.op.default_view().plot(ex2)

    def testPlot1D(self):
        self.op.channels = ["X"]
        self.op.scale = {"X" : "linear"}
        self.op.estimate(self.ex)
        self.op.default_view().plot(self.ex)


class TestPhenographBy(unittest.TestCase):

    def setUp(self):
        # tube A has blobs 1 and 2 and half of blob 0; tube B has the other
        # half of blob 0 and blob 3, which is its biggest.
        self.ex = make_blobs([300, 200, 100, 400],
                             tubes = {'A' : lambda d: d['Truth'].isin([1, 2]) | ((d['Truth'] == 0) & (d['X'] > 0)),
                                      'B' : lambda d: (d['Truth'] == 3) | ((d['Truth'] == 0) & (d['X'] <= 0))})

        self.op = flow.PhenographOp(name = "PG",
                                    channels = ["X", "Y"],
                                    scale = {"X" : "linear",
                                             "Y" : "linear"},
                                    num_neighbors = 15,
                                    resolution = 0.01,
                                    by = ["Tube"])

    def testEstimateBy(self):
        self.op.estimate(self.ex)
        self.assertEqual(len(self.op._centers), 2)

        ex2 = self.op.apply(self.ex)
        self.assertEqual(len(ex2), len(self.ex))

        # one clustering per tube: three clusters in A, two in B.
        a = ex2.data[ex2.data['Tube'] == 'A']
        b = ex2.data[ex2.data['Tube'] == 'B']
        self.assertEqual(sorted(a['PG'].astype(int).unique()), [0, 1, 2])
        self.assertEqual(sorted(b['PG'].astype(int).unique()), [0, 1])

        # each blob is one cluster, and in each tube they're numbered by size.
        for tube, blobs_by_size in [(a, [1.0, 0.0, 2.0]), (b, [3.0, 0.0])]:
            for cluster, blob in enumerate(blobs_by_size):
                self.assertEqual(tube[tube['PG'].astype(int) == cluster]['Truth'].unique().tolist(),
                                 [blob])

    def testStatisticBy(self):
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)

        stat = ex2.statistics['PG']
        self.assertEqual(list(stat.index.names), ['Tube', 'PG'])
        self.assertEqual(list(stat.index),
                         [('A', 0), ('A', 1), ('A', 2), ('B', 0), ('B', 1)])
        self.assertFalse(stat.isna().any().any())
        self.assertAlmostEqual(stat.loc['A', 'Proportion'].sum(), 1.0)
        self.assertAlmostEqual(stat.loc['B', 'Proportion'].sum(), 1.0)
        np.testing.assert_allclose(stat.loc[('B', 0), ['X', 'Y']].values.astype(float),
                                   [100, 100], atol = 0.5)

    def testApplyMissingGroup(self):
        self.op.estimate(self.ex, subset = "Tube == 'A'")
        with self.assertRaises(util.CytoflowOpError):
            self.op.apply(self.ex)

    def testPlotBy(self):
        self.op.estimate(self.ex)
        ex2 = self.op.apply(self.ex)
        dv = self.op.default_view()
        for v in dv.enum_plots(ex2):
            self.op.default_view().plot(ex2, plot_name = v)


if __name__ == "__main__":
#     import sys;sys.argv = ['', 'TestPhenograph.testEstimate']
    unittest.main()
