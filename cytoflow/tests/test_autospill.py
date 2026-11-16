'''
Tests for AutoSpillOp (a port of the AutoSpill R package)
'''

import os, shutil, tempfile, unittest, warnings

import numpy as np
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util
from cytoflow.operations.autospill import autospill, robust_linear_fit

CHANNELS = ["B1-A", "B2-A", "R1-A", "V1-A"]


def make_spillover(rng, n):
    s = np.eye(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                s[i, j] = rng.uniform(0, 0.25) * (0.5 ** abs(i - j))
    return s


def make_control(rng, s, i, n_events = 10000, outliers = True):
    n = len(s)
    true = np.zeros((n_events, n))
    pos = rng.random(n_events) < 0.6
    true[:, i] = np.where(pos, rng.lognormal(8, 1.0, n_events), rng.normal(0, 30, n_events))
    raw = true @ s
    raw += rng.normal(0, 40, raw.shape) + rng.normal(0, 1, raw.shape) * np.sqrt(np.abs(raw)) * 1.5
    if outliers:
        out = rng.random(n_events) < 0.01
        raw[out] += rng.lognormal(6, 1.5, (out.sum(), n))
    return raw


class TestAutoSpillAlgorithm(unittest.TestCase):
    
    def setUp(self):
        self.rng = np.random.default_rng(1)
        
    def testRobustFit(self):
        import statsmodels.api as sm
        x = self.rng.normal(0, 1, 2000)
        y = 2 + 0.3 * x + self.rng.standard_t(2, 2000)
        ref = sm.RLM(y, sm.add_constant(x), M = sm.robust.norms.HuberT()).fit().params
        np.testing.assert_allclose(robust_linear_fit(x, y), ref, atol = 2e-3)
        
    def testRobustFitDegenerate(self):
        self.assertEqual(robust_linear_fit([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]), (2.0, 0.0))
        self.assertEqual(robust_linear_fit([0.0, 1.0, 2.0], [1.0, 3.0, 5.0]), (1.0, 2.0))
        
    def testRecoverMatrix(self):
        for n in (2, 5, 8):
            s = make_spillover(self.rng, n)
            data = [make_control(self.rng, s, i) for i in range(n)]
            est, conv, err = autospill(data)
            np.testing.assert_allclose(est, s, atol = 2e-3)
            np.testing.assert_allclose(np.diag(est), 1.0)
            self.assertLess(np.abs(err).max(), 1e-2)
            self.assertEqual(conv["scale"].iloc[0], "linear")
            self.assertEqual(conv["scale"].iloc[-1], "arcsinh")
            
    def testCompensatedControlsAreFlat(self):
        s = make_spillover(self.rng, 4)
        data = [make_control(self.rng, s, i, outliers = False) for i in range(4)]
        est, _, _ = autospill(data)
        comp = data[0] @ np.linalg.inv(est)
        for j in range(1, 4):
            self.assertLess(abs(np.polyfit(comp[:, 0], comp[:, j], 1)[0]), 2e-3)
            
    def testNoConvergence(self):
        # two identical "controls" can't be separated
        s = make_spillover(self.rng, 2)
        d = make_control(self.rng, s, 0)
        d[:, 1] = d[:, 0]
        with self.assertRaises(util.CytoflowOpError):
            autospill([d, d])


class TestAutoSpillOp(unittest.TestCase):
    
    def setUp(self):
        self.rng = np.random.default_rng(2)
        self.dir = tempfile.mkdtemp()
        n = len(CHANNELS)
        self.s = make_spillover(self.rng, n)
        ranges = {c : 262144.0 for c in CHANNELS + ["FSC-A"]}
        
        def write(name, data):
            path = os.path.join(self.dir, name)
            fsc = self.rng.normal(50000, 5000, (len(data), 1))
            util.write_fcs(path, CHANNELS + ["FSC-A"], ranges, np.hstack([data, fsc]),
                           compat_chn_names = False, compat_negative = False)
            return path
            
        self.controls = {c : write("control_{}.fcs".format(i), make_control(self.rng, self.s, i, 5000))
                         for i, c in enumerate(CHANNELS)}
        
        # a "sample": every color, independent
        true = self.rng.lognormal(7, 1, (5000, n))
        self.true = true
        self.sample = write("sample.fcs", true @ self.s)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.ex = flow.ImportOp(tubes = [flow.Tube(file = self.sample)]).apply()
        
    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors = True)
        plt.close('all')
        
    def _op(self):
        return flow.AutoSpillOp(controls = dict(self.controls))
        
    def testEstimateAndApply(self):
        op = self._op()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            op.estimate(self.ex)
        
        est = np.array([[op.spillover[f][t] for t in CHANNELS] for f in CHANNELS])
        np.testing.assert_allclose(est, self.s, atol = 5e-3)
        
        ex2 = op.apply(self.ex)
        np.testing.assert_allclose(ex2.data[CHANNELS].values, self.true, 
                                   rtol = 0.05, atol = 60)
        self.assertTrue(ex2.metadata[CHANNELS[0]]['compensated'])
        # FSC is untouched
        np.testing.assert_array_equal(ex2.data["FSC-A"].values, self.ex.data["FSC-A"].values)
        
    def testSubsetAndHistory(self):
        ex = flow.ThresholdOp(name = "Big", channel = "FSC-A", threshold = 50000).apply(self.ex)
        op = self._op()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            op.estimate(ex, subset = "Big == True")
        est = np.array([[op.spillover[f][t] for t in CHANNELS] for f in CHANNELS])
        np.testing.assert_allclose(est, self.s, atol = 1e-2)
        
    def testViews(self):
        op = self._op()
        with self.assertRaises(util.CytoflowViewError):
            op.default_view(convergence = True).plot(self.ex)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            op.estimate(self.ex)
        op.default_view().plot(self.ex)
        op.default_view(convergence = True).plot(self.ex)
        
    def testSwappedControlsWarn(self):
        op = self._op()
        op.controls[CHANNELS[0]], op.controls[CHANNELS[1]] = \
            op.controls[CHANNELS[1]], op.controls[CHANNELS[0]]
        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter("always")
            try:
                op.estimate(self.ex)
            except util.CytoflowOpError:
                pass
        self.assertTrue(any("Is it the right file" in str(x.message) for x in w))
        
    def testErrors(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.AutoSpillOp(controls = {CHANNELS[0] : self.controls[CHANNELS[0]]}).estimate(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.AutoSpillOp(controls = {CHANNELS[0] : self.controls[CHANNELS[0]],
                                         CHANNELS[1] : "/no/such/file.fcs"}).estimate(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.AutoSpillOp(controls = {CHANNELS[0] : self.controls[CHANNELS[0]],
                                         CHANNELS[1] : self.controls[CHANNELS[0]]}).estimate(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.AutoSpillOp(controls = {"nope" : self.controls[CHANNELS[0]],
                                         CHANNELS[1] : self.controls[CHANNELS[1]]}).estimate(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op().apply(self.ex)


if __name__ == "__main__":
    unittest.main()
