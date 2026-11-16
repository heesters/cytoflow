'''
Tests for ScaleSettingsOp, the arcsinh scale and per-channel scale parameters.
'''

import unittest

import numpy as np
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util

from .test_base import ImportedDataTest


class TestArcsinhScale(ImportedDataTest):
    
    def tearDown(self):
        plt.close('all')
    
    def testRoundTrip(self):
        s = util.scale_factory("arcsinh", self.ex, channel = "Y2-A", cofactor = 5)
        x = np.array([-500.0, 0.0, 3.0, 1000.0, 1e5])
        np.testing.assert_allclose(s(x), np.arcsinh(x / 5))
        np.testing.assert_allclose(s.inverse(s(x)), x)
        self.assertAlmostEqual(s(5.0), np.arcsinh(1.0))
        
    def testPlots(self):
        flow.HistogramView(channel = "Y2-A", scale = "arcsinh").plot(self.ex)
        flow.ScatterplotView(xchannel = "V2-A", xscale = "arcsinh",
                             ychannel = "Y2-A", yscale = "arcsinh").plot(self.ex)
        flow.DensityView(xchannel = "V2-A", xscale = "arcsinh",
                         ychannel = "Y2-A", yscale = "arcsinh").plot(self.ex)
        
    def testGate(self):
        ex = flow.PolygonOp(name = "P", xchannel = "V2-A", ychannel = "Y2-A",
                            xscale = "arcsinh", yscale = "arcsinh",
                            vertices = [(0, 0), (1000, 0), (1000, 20000), (0, 20000)]).apply(self.ex)
        self.assertGreater(ex["P"].sum(), 0)
        
    def testBadCofactor(self):
        with self.assertRaises(util.CytoflowError):
            util.scale_factory("arcsinh", self.ex, channel = "Y2-A", cofactor = 0)(1.0)


class TestScaleSettings(ImportedDataTest):
    
    def tearDown(self):
        plt.close('all')

    def testApply(self):
        op = flow.ScaleSettingsOp(settings = {"Y2-A" : {"logicle" : {"W" : 1.5, "M" : 5.0},
                                                        "arcsinh" : {"cofactor" : 5.0}}})
        ex2 = op.apply(self.ex)
        
        s = util.scale_factory("logicle", ex2, channel = "Y2-A")
        self.assertEqual(s.W, 1.5)
        self.assertEqual(s.M, 5.0)
        self.assertEqual(util.scale_factory("arcsinh", ex2, channel = "Y2-A").cofactor, 5.0)
        
        # other channels and the original experiment are unaffected
        self.assertEqual(util.scale_factory("arcsinh", ex2, channel = "V2-A").cofactor, 150.0)
        self.assertNotIn('scale_params', self.ex.metadata["Y2-A"])
        
        # explicit parameters still win
        self.assertEqual(util.scale_factory("logicle", ex2, channel = "Y2-A", W = 0.7).W, 0.7)
        
        # and views use them
        flow.HistogramView(channel = "Y2-A", scale = "logicle").plot(ex2)
        
    def testMerge(self):
        ex2 = flow.ScaleSettingsOp(settings = {"Y2-A" : {"logicle" : {"W" : 1.5}}}).apply(self.ex)
        ex3 = flow.ScaleSettingsOp(settings = {"Y2-A" : {"logicle" : {"M" : 5.0}}}).apply(ex2)
        self.assertEqual(ex3.metadata["Y2-A"]["scale_params"]["logicle"], {"W" : 1.5, "M" : 5.0})
        self.assertEqual(ex2.metadata["Y2-A"]["scale_params"]["logicle"], {"W" : 1.5})
        
    def testExplicitLogicleW(self):
        # this used to be ignored (the cached W didn't depend on the set value)
        self.assertEqual(util.scale_factory("logicle", self.ex, channel = "Y2-A", W = 1.5).W, 1.5)

    def testErrors(self):
        for settings in [{"XXX" : {"logicle" : {"W" : 1.0}}},
                         {"Y2-A" : {"nope" : {"W" : 1.0}}},
                         {"Y2-A" : {"logicle" : {"nope" : 1.0}}},
                         {"Y2-A" : {"logicle" : {"r" : 2.0}}},
                         {"Y2-A" : {"arcsinh" : {"cofactor" : -1.0}}},
                         {}]:
            with self.assertRaises(util.CytoflowOpError):
                flow.ScaleSettingsOp(settings = settings).apply(self.ex)


if __name__ == "__main__":
    unittest.main()
