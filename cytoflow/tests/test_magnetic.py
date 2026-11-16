'''
Tests for "magnetic" (per-group adjusted) polygon and ellipse gates.
'''

import unittest

import numpy as np
import pandas as pd

import cytoflow as flow
import cytoflow.utility as util

SHIFTS = [(0, 0), (300, -200), (-250, 150)]


class TestMagneticGates(unittest.TestCase):
    
    def setUp(self):
        rng = np.random.default_rng(0)
        ex = flow.Experiment()
        ex.add_channel("X")
        ex.add_channel("Y")
        ex.add_channel("Truth")
        ex.add_condition("Tube", "int")
        for t, (dx, dy) in enumerate(SHIFTS):
            pop = rng.normal([1000 + dx, 1000 + dy], 80, size = (3000, 2))
            bg = rng.uniform(0, 3000, size = (1000, 2))
            d = np.vstack([pop, bg])
            ex.add_events(pd.DataFrame({"X" : d[:, 0], "Y" : d[:, 1],
                                        "Truth" : np.r_[np.ones(3000), np.zeros(1000)]}),
                          {"Tube" : t})
        self.ex = ex
        self.verts = [(800, 800), (1200, 800), (1200, 1200), (800, 1200)]
        
    def _recall(self, ex, name):
        return {t : float((g[name] & (g["Truth"] > 0.5)).sum() / (g["Truth"] > 0.5).sum())
                for t, g in ex.data.groupby("Tube")}
        
    def testPolygon(self):
        fixed = flow.PolygonOp(name = "P", xchannel = "X", ychannel = "Y",
                               vertices = self.verts).apply(self.ex)
        mag = flow.PolygonOp(name = "P", xchannel = "X", ychannel = "Y",
                             vertices = self.verts, adjust_by = ["Tube"],
                             max_adjust = 1.0).apply(self.ex)
        
        self.assertLess(self._recall(fixed, "P")[1], 0.2)
        for r in self._recall(mag, "P").values():
            self.assertGreater(r, 0.95)
            
        # the shifts follow the populations (relative to tube 0)
        shifts = mag.statistics["P"]
        for t, (dx, dy) in enumerate(SHIFTS):
            self.assertAlmostEqual(shifts.loc[t, "x_shift"] - shifts.loc[0, "x_shift"], dx, delta = 30)
            self.assertAlmostEqual(shifts.loc[t, "y_shift"] - shifts.loc[0, "y_shift"], dy, delta = 30)
            
    def testFreehand(self):
        mag = flow.FreehandOp(name = "P", xchannel = "X", ychannel = "Y",
                              vertices = self.verts, adjust_by = ["Tube"],
                              max_adjust = 1.0).apply(self.ex)
        for r in self._recall(mag, "P").values():
            self.assertGreater(r, 0.95)
            
    def testEllipse(self):
        mag = flow.EllipseOp(name = "E", xchannel = "X", ychannel = "Y",
                             center = (1000, 1000), width = 500, height = 500,
                             adjust_by = ["Tube"], max_adjust = 1.0).apply(self.ex)
        for r in self._recall(mag, "E").values():
            self.assertGreater(r, 0.95)
            
    def testMaxAdjust(self):
        mag = flow.PolygonOp(name = "P", xchannel = "X", ychannel = "Y",
                             vertices = self.verts, adjust_by = ["Tube"],
                             max_adjust = 0.25).apply(self.ex)
        # the gate is 400 wide, so it can't move more than 100
        self.assertLessEqual(mag.statistics["P"].abs().max().max(), 100.0 + 1e-9)
        
    def testNoAdjustment(self):
        # without adjust_by, nothing changes and no statistic is added
        ex = flow.PolygonOp(name = "P", xchannel = "X", ychannel = "Y",
                            vertices = self.verts).apply(self.ex)
        self.assertNotIn("P", ex.statistics)
        
    def testBadCondition(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.PolygonOp(name = "P", xchannel = "X", ychannel = "Y",
                           vertices = self.verts, adjust_by = ["Nope"]).apply(self.ex)


if __name__ == "__main__":
    unittest.main()
