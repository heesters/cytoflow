'''
Tests for ProliferationOp and CellCycleOp
'''

import unittest, warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util
from cytoflow.operations.proliferation import (fit_generations, proliferation_metrics, LOG2)
from cytoflow.operations.cell_cycle import fit_cell_cycle


def dye(rng, fracs, n = 8000, mu0 = 5.0, sigma = 0.07):
    g = rng.choice(len(fracs), size = n, p = fracs)
    return 10 ** (mu0 - g * LOG2 + rng.normal(0, sigma, n)), g


class TestProliferationMetrics(unittest.TestCase):
    
    def testHandComputed(self):
        # N = [100, 200, 400] -> precursors P = [100, 100, 100]
        m = proliferation_metrics([100, 200, 400])
        self.assertAlmostEqual(m["Percent_divided"], 200 / 3)
        self.assertAlmostEqual(m["Division_index"], (0 * 100 + 1 * 100 + 2 * 100) / 300)
        self.assertAlmostEqual(m["Proliferation_index"], (1 * 100 + 2 * 100) / 200)
        self.assertAlmostEqual(m["Expansion_index"], 700 / 300)
        self.assertAlmostEqual(m["Replication_index"], 600 / 200)
        
    def testUndivided(self):
        m = proliferation_metrics([500, 0, 0])
        self.assertEqual(m["Percent_divided"], 0)
        self.assertEqual(m["Division_index"], 0)
        self.assertEqual(m["Expansion_index"], 1)
        self.assertTrue(np.isnan(m["Proliferation_index"]))
        self.assertTrue(np.isnan(m["Replication_index"]))
        
    def testOneDivisionEach(self):
        # every cell divided exactly once: 2 daughters per precursor
        m = proliferation_metrics([0, 1000])
        self.assertAlmostEqual(m["Percent_divided"], 100)
        self.assertAlmostEqual(m["Division_index"], 1)
        self.assertAlmostEqual(m["Proliferation_index"], 1)
        self.assertAlmostEqual(m["Expansion_index"], 2)
        self.assertAlmostEqual(m["Replication_index"], 2)


class TestFitGenerations(unittest.TestCase):
    
    def setUp(self):
        self.rng = np.random.default_rng(0)
        
    def testRecoverFractions(self):
        for fracs in ([0.1, 0.15, 0.25, 0.3, 0.15, 0.05], [0.6, 0.4], [1.0]):
            x, _ = dye(self.rng, fracs, n = 20000)
            fit = fit_generations(np.log10(x))
            w = fit['weights']
            for g, f in enumerate(fracs):
                self.assertAlmostEqual(w[g], f, delta = 0.02)
            self.assertAlmostEqual(w[len(fracs):].sum(), 0, delta = 0.02)
            self.assertAlmostEqual(fit['mu0'], 5.0, delta = 0.02)
            self.assertAlmostEqual(fit['sigma'], 0.07, delta = 0.01)
            
    def testFixedPeak(self):
        # nearly everything has divided: generation 0 can't be found from 
        # the data, but it can be given
        fracs = [0.0, 0.02, 0.18, 0.5, 0.3]
        x, _ = dye(self.rng, fracs, n = 20000)
        fit = fit_generations(np.log10(x), undivided_peak = 5.0)
        np.testing.assert_allclose(fit['weights'][:5], fracs, atol = 0.02)
        
    def testFitSpacing(self):
        g = self.rng.choice(4, size = 20000, p = [0.25, 0.25, 0.25, 0.25])
        x = 5.0 - g * LOG2 * 1.08 + self.rng.normal(0, 0.05, 20000)
        fit = fit_generations(x, fit_spacing = True)
        self.assertAlmostEqual(fit['spacing'], 1.08 * LOG2, delta = 0.005)
        np.testing.assert_allclose(fit['weights'][:4], 0.25, atol = 0.02)
        
    def testTooFew(self):
        with self.assertRaises(util.CytoflowOpError):
            fit_generations(np.ones(10))


class TestProliferationOp(unittest.TestCase):
    
    def setUp(self):
        rng = np.random.default_rng(1)
        self.fracs = {"none" : [1.0], 
                      "low" : [0.5, 0.3, 0.2], 
                      "high" : [0.05, 0.1, 0.2, 0.35, 0.3]}
        ex = flow.Experiment()
        ex.add_channel("CFSE")
        ex.add_condition("Stim", "category")
        for stim, fracs in self.fracs.items():
            x, _ = dye(rng, fracs)
            ex.add_events(pd.DataFrame({"CFSE" : x}), {"Stim" : stim})
        self.ex = ex
        
    def tearDown(self):
        plt.close('all')
        
    def testByWithControl(self):
        op = flow.ProliferationOp(name = "Prolif", channel = "CFSE", by = ["Stim"],
                                  control_subset = "Stim == 'none'")
        op.estimate(self.ex)
        ex2 = op.apply(self.ex)
        stat = ex2.statistics["Prolif"]
        
        for stim, fracs in self.fracs.items():
            for g, f in enumerate(fracs):
                self.assertAlmostEqual(stat.loc[stim, "Generation_{}".format(g)], 100 * f, delta = 3)
            truth = proliferation_metrics(fracs)
            for k, v in truth.items():
                if np.isnan(v):
                    continue
                self.assertAlmostEqual(stat.loc[stim, k], v, delta = 0.05 * max(1.0, abs(v)), 
                                       msg = "{} {}".format(stim, k))
            self.assertAlmostEqual(stat.loc[stim, "Undivided_peak"], 1e5, delta = 5e3)
            
        self.assertAlmostEqual(stat.loc["none", "Percent_divided"], 0, delta = 2)
        
        # the per-event generation
        gen = ex2.data.groupby("Stim", observed = True)["Prolif"]
        self.assertGreater((gen.get_group("none") == 0).mean(), 0.97)
        self.assertGreater((gen.get_group("high") >= 3).mean(), 0.55)
        
    def testPooledAndViews(self):
        op = flow.ProliferationOp(name = "Prolif", channel = "CFSE")
        with self.assertRaises(util.CytoflowViewError):
            op.default_view().plot(self.ex)
        op.estimate(self.ex, subset = "Stim == 'low'")
        ex2 = op.apply(self.ex)
        self.assertEqual(list(ex2.statistics["Prolif"].index), ["All"])
        op.default_view().plot(ex2)
        
        op2 = flow.ProliferationOp(name = "P2", channel = "CFSE", by = ["Stim"], 
                                   undivided_peak = 1e5)
        op2.estimate(self.ex)
        op2.default_view().plot(self.ex, title = "Proliferation")
        op2.default_view(subset = "Stim == 'high'").plot(self.ex)
        
    def testNonPositive(self):
        self.ex.data.loc[self.ex.data.index[:2000], "CFSE"] = -5.0
        op = flow.ProliferationOp(name = "Prolif", channel = "CFSE", by = ["Stim"])
        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter("always")
            op.estimate(self.ex)
        self.assertTrue(any("signal <= 0" in str(x.message) for x in w))
        op.apply(self.ex)
        
    def testErrors(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.ProliferationOp(name = "P", channel = "CFSE").apply(self.ex)
        for kwargs in [dict(channel = ""), dict(channel = "nope"), dict(by = ["nope"]),
                       dict(num_generations = 0), dict(undivided_peak = -1.0),
                       dict(control_subset = "Stim == 'nothing'")]:
            op = flow.ProliferationOp(name = "P", channel = "CFSE")
            op.trait_set(**kwargs)
            with self.assertRaises(util.CytoflowOpError):
                op.estimate(self.ex)


def dna(rng, g1, s, g2, n = 20000, m = 50000.0, cv = 0.04, r = 2.0):
    n1, ns = int(n * g1), int(n * s)
    pos = m * (1 + rng.uniform(0, 1, ns) * (r - 1))
    return np.concatenate([rng.normal(m, cv * m, n1),
                           rng.normal(pos, cv * pos),
                           rng.normal(r * m, cv * r * m, n - n1 - ns)])


class TestCellCycle(unittest.TestCase):
    
    def setUp(self):
        self.rng = np.random.default_rng(2)
        
    def tearDown(self):
        plt.close('all')
        
    def testRecoverPhases(self):
        for g1, s, g2, cv in [(0.6, 0.25, 0.15, 0.04), (0.45, 0.4, 0.15, 0.03), 
                              (0.7, 0.1, 0.2, 0.06), (0.85, 0.1, 0.05, 0.04)]:
            fit = fit_cell_cycle(dna(self.rng, g1, s, g2, cv = cv))
            self.assertAlmostEqual(fit["G1"], g1, delta = 0.03)
            self.assertAlmostEqual(fit["S"], s, delta = 0.03)
            self.assertAlmostEqual(fit["G2M"], g2, delta = 0.03)
            self.assertAlmostEqual(fit["G1_CV"], cv, delta = 0.005)
            self.assertAlmostEqual(fit["G2_G1_ratio"], 2.0, delta = 0.03)
            self.assertAlmostEqual(fit["G1_mean"], 50000, delta = 500)
            self.assertLess(fit["RMSE"], 0.05)
            
    def testScaleInvariant(self):
        # the same cells on an instrument with a 100x smaller range
        x = dna(self.rng, 0.6, 0.25, 0.15)
        a, b = fit_cell_cycle(x), fit_cell_cycle(x / 100.0)
        for k in ("G1", "S", "G2M", "G1_CV"):
            self.assertAlmostEqual(a[k], b[k], delta = 0.005)
            
    def testSubG1AndAggregates(self):
        x = np.concatenate([dna(self.rng, 0.6, 0.25, 0.15, n = 18000),
                            self.rng.uniform(5000, 30000, 1000),      # debris / apoptotic
                            self.rng.normal(150000, 5000, 1000)])     # aggregates
        fit = fit_cell_cycle(x)
        self.assertAlmostEqual(fit["SubG1"], 0.05, delta = 0.01)
        self.assertAlmostEqual(fit["AboveG2M"], 0.05, delta = 0.01)
        self.assertAlmostEqual(fit["G1"], 0.6, delta = 0.04)
        
    def testOp(self):
        ex = flow.Experiment()
        ex.add_channel("DNA")
        ex.add_condition("Drug", "category")
        truth = {"ctrl" : (0.6, 0.25, 0.15), "arrest" : (0.2, 0.15, 0.65)}
        for drug, (g1, s, g2) in truth.items():
            ex.add_events(pd.DataFrame({"DNA" : dna(self.rng, g1, s, g2, n = 10000)}), {"Drug" : drug})
            
        op = flow.CellCycleOp(name = "CC", channel = "DNA", by = ["Drug"])
        with self.assertRaises(util.CytoflowOpError):
            op.apply(ex)
        op.estimate(ex)
        ex2 = op.apply(ex)
        stat = ex2.statistics["CC"]
        for drug, (g1, s, g2) in truth.items():
            self.assertAlmostEqual(stat.loc[drug, "G1"], 100 * g1, delta = 4)
            self.assertAlmostEqual(stat.loc[drug, "S"], 100 * s, delta = 4)
            self.assertAlmostEqual(stat.loc[drug, "G2M"], 100 * g2, delta = 4)
            self.assertAlmostEqual(stat.loc[drug, "G1"] + stat.loc[drug, "S"] + stat.loc[drug, "G2M"], 100)
            
        # the per-event phases roughly agree with the statistic
        phases = ex2.data[ex2.data["Drug"] == "arrest"]["CC"].value_counts(normalize = True)
        self.assertAlmostEqual(phases["G2M"], 0.65, delta = 0.08)
        self.assertEqual(set(ex2.data["CC"].unique()) - {"G1", "S", "G2M", "SubG1", "AboveG2M"}, set())
        
        op.default_view().plot(ex2)
        op.default_view(subset = "Drug == 'ctrl'").plot(ex2, title = "Control")
        
        # pooled
        op2 = flow.CellCycleOp(name = "CC2", channel = "DNA")
        op2.estimate(ex, subset = "Drug == 'ctrl'")
        self.assertEqual(list(op2.apply(ex).statistics["CC2"].index), ["All"])
        
    def testErrors(self):
        with self.assertRaises(util.CytoflowOpError):
            fit_cell_cycle(np.ones(50))
        with self.assertRaises(util.CytoflowOpError):
            fit_cell_cycle(np.ones(5000))
        # one peak only: a warning, not nonsense
        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter("always")
            fit = fit_cell_cycle(self.rng.normal(50000, 2000, 5000))
        self.assertTrue(any("Only one peak" in str(x.message) for x in w))
        self.assertGreater(fit["G1"], 0.95)


if __name__ == "__main__":
    unittest.main()
