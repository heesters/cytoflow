'''
Tests for the cross entropy test (CrossEntropyTestOp)
'''

import unittest, warnings

import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util
from cytoflow.operations.cross_entropy import (pointwise_cross_entropy, 
                                               _conditional_probabilities, _adjust)


def brute_force(high, low, perplexity):
    """The paper's equations, with dense matrices and all pairs."""
    n = len(high)
    d2 = ((high[:, None, :] - high[None, :, :]) ** 2).sum(axis = 2)
    cond = np.zeros((n, n))
    for i in range(n):
        others = np.delete(np.arange(n), i)
        cond[i, others] = _conditional_probabilities(d2[i, others][None, :], perplexity)[0]
    p = (cond + cond.T) / (2 * n)
    
    e2 = ((low[:, None, :] - low[None, :, :]) ** 2).sum(axis = 2)
    w = 1 / (1 + e2)
    np.fill_diagonal(w, 0)
    q = w / w.sum()
    
    h = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j and p[i, j] > 0:
                h[i] -= n * p[i, j] * np.log(n * q[i, j])
    return h, p, q


class TestPointwiseCrossEntropy(unittest.TestCase):
    
    def setUp(self):
        self.rng = np.random.default_rng(0)
        
    def testPerplexity(self):
        d2 = self.rng.uniform(0, 10, (50, 40))
        p = _conditional_probabilities(d2, 12.0)
        np.testing.assert_allclose(p.sum(axis = 1), 1.0)
        perp = np.exp(-(p * np.log(p)).sum(axis = 1))
        np.testing.assert_allclose(perp, 12.0, rtol = 1e-3)
        
    def testAgainstBruteForce(self):
        high = self.rng.normal(size = (120, 5))
        low = self.rng.normal(size = (120, 2))
        ref, p, q = brute_force(high, low, 10.0)
        self.assertAlmostEqual(p.sum(), 1.0)
        self.assertAlmostEqual(q.sum(), 1.0)
        # all neighbors -> exactly the dense computation
        ce = pointwise_cross_entropy(high, low, perplexity = 10.0, num_neighbors = 119)
        np.testing.assert_allclose(ce, ref, rtol = 1e-4)
        # the default (3 * perplexity neighbors) is a close approximation
        ce = pointwise_cross_entropy(high, low, perplexity = 10.0)
        np.testing.assert_allclose(ce, ref, rtol = 0.05)
        
    def testGoodEmbeddingHasLowerCrossEntropy(self):
        high = np.vstack([self.rng.normal(0, 1, (150, 4)), self.rng.normal(8, 1, (150, 4))])
        good = high[:, :2] + self.rng.normal(0, 0.1, (300, 2))
        bad = self.rng.normal(size = (300, 2))
        self.assertLess(pointwise_cross_entropy(high, good, 15).mean(),
                        pointwise_cross_entropy(high, bad, 15).mean())
        
    def testDuplicatePoints(self):
        high = np.vstack([np.zeros((30, 3)), self.rng.normal(size = (100, 3))])
        low = self.rng.normal(size = (130, 2))
        self.assertTrue(np.all(np.isfinite(pointwise_cross_entropy(high, low, 10))))
        
    def testAdjust(self):
        p = [0.01, 0.04, 0.03, 0.5]
        np.testing.assert_allclose(_adjust(p, 'bonferroni'), [0.04, 0.16, 0.12, 1.0])
        np.testing.assert_allclose(_adjust(p, 'holm'), [0.04, 0.09, 0.09, 0.5])
        np.testing.assert_allclose(_adjust(p, 'none'), p)


class TestCrossEntropyTestOp(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(3)
        n = 600
        
        def cells(frac_b):
            # two populations with different *structure*: a round blob, and
            # a long thin filament.  (the test compares the distributions of
            # cross entropy, which depend on the local structure.)
            nb = int(n * frac_b)
            a = rng.normal([200, 200, 200], 40, (n - nb, 3))
            t = rng.uniform(0, 1, nb)
            b = np.column_stack([900 + 600 * t, 800 + 500 * t, 300 + 50 * t]) \
                + rng.normal(0, 4, (nb, 3))
            return np.vstack([a, b])
        
        # samples 1 and 2 have the same composition; sample 3 is different
        ex = flow.Experiment()
        for c in ["A", "B", "C"]:
            ex.add_channel(c)
        ex.add_condition("Sample", "int")
        for s, frac in ((1, 0.2), (2, 0.2), (3, 0.8)):
            ex.add_events(pd.DataFrame(cells(frac), columns = ["A", "B", "C"]),
                          {"Sample" : s})
            
        from sklearn.manifold import TSNE
        emb = TSNE(perplexity = 30, random_state = 0, init = 'pca').fit_transform(ex.data[["A", "B", "C"]].values)
        ex.add_channel("tSNE_1", emb[:, 0])
        ex.add_channel("tSNE_2", emb[:, 1])
        cls.ex = ex
        
    def tearDown(self):
        plt.close('all')
        
    def _op(self, **kwargs):
        return flow.CrossEntropyTestOp(name = "CE", by = "Sample",
                                       channels = ["A", "B", "C"],
                                       embedding = ["tSNE_1", "tSNE_2"],
                                       perplexity = 30, **kwargs)
        
    def testTest(self):
        op = self._op()
        op.estimate(self.ex)
        ex2 = op.apply(self.ex)
        
        pairs = ex2.statistics["CE"]
        self.assertEqual(list(pairs.index.names), ["group1", "group2"])
        self.assertEqual(len(pairs), 3)
        
        # same composition: not different.  different composition: different
        self.assertFalse(pairs.loc[(1, 2), "reject"])
        self.assertTrue(pairs.loc[(1, 3), "reject"])
        self.assertTrue(pairs.loc[(2, 3), "reject"])
        self.assertLess(pairs.loc[(1, 2), "D"], pairs.loc[(1, 3), "D"])
        
        # the statistics are what scipy gives for the stored distributions
        res = scipy.stats.ks_2samp(op._cross_entropy[1], op._cross_entropy[3])
        self.assertEqual(pairs.loc[(1, 3), "D"], res.statistic)
        self.assertEqual(pairs.loc[(1, 3), "p"], res.pvalue)
        self.assertGreaterEqual(pairs.loc[(1, 3), "p_adj"], pairs.loc[(1, 3), "p"])
        
        groups = ex2.statistics["CE_Groups"]
        self.assertEqual(list(groups["n"]), [600, 600, 600])
        
        m = op.distance_matrix()
        self.assertEqual(m.loc[1, 3], m.loc[3, 1])
        self.assertEqual(m.loc[1, 1], 0)
        
    def testSubsetAndMaxEvents(self):
        op = self._op(max_events = 900)
        op.estimate(self.ex, subset = "Sample != 2")
        self.assertEqual(sorted(op._cross_entropy.keys()), [1, 3])
        self.assertEqual(sum(len(v) for v in op._cross_entropy.values()), 900)
        self.assertEqual(len(op._pairs), 1)
        
    def testDeterministic(self):
        op1 = self._op(max_events = 900); op1.estimate(self.ex)
        op2 = self._op(max_events = 900); op2.estimate(self.ex)
        pd.testing.assert_frame_equal(op1._pairs, op2._pairs)
        
    def testViews(self):
        op = self._op()
        with self.assertRaises(util.CytoflowViewError):
            op.default_view().plot(self.ex)
        op.estimate(self.ex)
        op.default_view().plot(self.ex)
        op.default_view(dendrogram = True).plot(self.ex)
        
    def testErrors(self):
        with self.assertRaises(util.CytoflowOpError):
            self._op().apply(self.ex)                       # no estimate
        for kwargs in [dict(by = ""), dict(by = "nope"), dict(channels = ["A"]),
                       dict(embedding = ["tSNE_1"]), dict(embedding = ["A", "tSNE_1"]),
                       dict(channels = ["A", "nope"]), dict(perplexity = 0.5)]:
            op = self._op()
            op.trait_set(**kwargs)
            with self.assertRaises(util.CytoflowOpError):
                op.estimate(self.ex)
        op = self._op()
        with self.assertRaises(util.CytoflowOpError):
            op.estimate(self.ex, subset = "Sample == 1")    # one group


if __name__ == "__main__":
    unittest.main()
