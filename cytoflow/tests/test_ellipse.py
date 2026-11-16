#!/usr/bin/env python3.8
# coding: latin-1

# (c) Brian Teague 2018-2026
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

import unittest
import numpy as np
import pandas as pd

import cytoflow as flow
import cytoflow.utility as util
from cytoflow.operations.ellipse import (ellipse_contains,
                                         _selector_to_scaled,
                                         _scaled_to_selector)
from .test_base import ImportedDataSmallTest


class Test(ImportedDataSmallTest):
    def setUp(self):
        super().setUp()
        self.gate = flow.EllipseOp(name = "Ellipse",
                                   xchannel = "V2-A",
                                   ychannel = "Y2-A",
                                   xscale = "logicle",
                                   yscale = "logicle",
                                   center = (150.0, 5000.0),
                                   width = 0.3,
                                   height = 0.15,
                                   angle = 20.0)

    def testGate(self):
        ex2 = self.gate.apply(self.ex)

        self.assertIsInstance(ex2.data.index, pd.RangeIndex)
        self.assertEqual(ex2.data["Ellipse"].dtype, bool)

        # compute the expected membership by hand, in scaled space
        xscale = util.scale_factory("logicle", self.ex, channel = "V2-A")
        yscale = util.scale_factory("logicle", self.ex, channel = "Y2-A")
        x = xscale(self.ex.data["V2-A"]) - xscale(150.0)
        y = yscale(self.ex.data["Y2-A"]) - yscale(5000.0)
        t = np.deg2rad(20.0)
        u = x * np.cos(t) + y * np.sin(t)
        v = -x * np.sin(t) + y * np.cos(t)
        expected = (u / 0.15) ** 2 + (v / 0.075) ** 2 <= 1.0

        n = ex2.data.groupby("Ellipse").size()[True]
        self.assertGreater(n, 0)
        self.assertLess(n, len(ex2.data))
        self.assertEqual(n, expected.sum())
        np.testing.assert_array_equal(ex2.data["Ellipse"].values,
                                      expected.values)

    def testLinear(self):
        self.gate.xscale = self.gate.yscale = "linear"
        self.gate.center = (0.0, 0.0)
        self.gate.width = 2.0
        self.gate.height = 1.0
        self.gate.angle = 90.0

        self.ex.data.loc[0, ['V2-A', 'Y2-A']] = (0.0, 0.99)     # in
        self.ex.data.loc[1, ['V2-A', 'Y2-A']] = (0.49, 0.0)     # in
        self.ex.data.loc[2, ['V2-A', 'Y2-A']] = (0.0, 1.0)      # on the edge
        self.ex.data.loc[3, ['V2-A', 'Y2-A']] = (0.51, 0.0)     # out
        self.ex.data.loc[4, ['V2-A', 'Y2-A']] = (0.0, 1.01)     # out

        ex2 = self.gate.apply(self.ex)
        self.assertEqual(list(ex2.data["Ellipse"].values[0:5]),
                         [True, True, True, False, False])

    def testLogMaskedEvents(self):
        # events outside the log scale's domain are never in the gate
        self.gate.xscale = self.gate.yscale = "log"
        self.gate.width = 2.0
        self.gate.height = 2.0
        self.ex.data.loc[0, ['V2-A', 'Y2-A']] = (-10.0, 5000.0)
        ex2 = self.gate.apply(self.ex)
        self.assertFalse(ex2.data["Ellipse"].iloc[0])
        self.assertGreater(ex2.data["Ellipse"].sum(), 0)

    def testBadParams(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.EllipseOp(name = "Ellipse", xchannel = "V2-A",
                           ychannel = "Y2-A", center = (1.0, 1.0),
                           width = 0.0, height = 1.0).apply(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            flow.EllipseOp(xchannel = "V2-A", ychannel = "Y2-A",
                           center = (1.0, 1.0),
                           width = 1.0, height = 1.0).apply(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            flow.EllipseOp(name = "Ellipse", xchannel = "XXX",
                           ychannel = "Y2-A", center = (1.0, 1.0),
                           width = 1.0, height = 1.0).apply(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            flow.EllipseOp(name = "Ellipse", xchannel = "V2-A",
                           ychannel = "Y2-A", xscale = "log",
                           center = (-1.0, 1.0),
                           width = 1.0, height = 1.0).apply(self.ex)

    def testContains(self):
        xy = np.array([[0.0, 0.0], [1.9, 0.0], [0.0, 0.9], [2.1, 0.0],
                       [np.nan, 0.0]])
        np.testing.assert_array_equal(ellipse_contains(xy, (0, 0), 4, 2, 0),
                                      [True, True, True, False, False])
        np.testing.assert_array_equal(ellipse_contains(xy, (0, 0), 4, 2, 90),
                                      [True, False, True, False, False])

    def testSelectorConversion(self):
        for (w, h, a, aspect) in [(2.0, 1.0, 0.0, 1.0),
                                  (2.0, 1.0, 30.0, 2.5),
                                  (1.0, 3.0, -70.0, 0.3),
                                  (0.5, 0.5, 10.0, 4.0)]:
            hw, hh, rot = _scaled_to_selector(w, h, a, aspect)
            self.assertTrue(-45 < rot <= 45)
            w2, h2, a2 = _selector_to_scaled(hw, hh, rot, aspect)

            # the same ellipse?
            xy = np.random.default_rng(0).uniform(-3, 3, size = (2000, 2))
            np.testing.assert_array_equal(ellipse_contains(xy, (0, 0), w, h, a),
                                          ellipse_contains(xy, (0, 0), w2, h2, a2))

    def testPlot(self):
        self.gate.default_view().plot(self.ex)

    def testPlotParams(self):
        self.gate.default_view().plot(self.ex, patch_props = {'color' : 'blue'})

    def testDensityPlot(self):
        self.gate.default_view(density = True).plot(self.ex)

    def testInteractive(self):
        v = self.gate.default_view(interactive = True)
        v.plot(self.ex)

        # the existing gate is loaded into the selector
        self.assertTrue(v._widget._selection_completed)
        self.assertTrue(np.allclose(v._widget.center, v._scaled_center()))

        # and reading it back gives the same ellipse
        center, width, height, angle = v._get_selection()
        self.assertTrue(np.allclose(center, self.gate.center))
        self.assertAlmostEqual(width, self.gate.width)
        self.assertAlmostEqual(height, self.gate.height)
        self.assertAlmostEqual(angle, self.gate.angle)

        v.interactive = False
        self.assertIsNone(v._widget)

    def testInteractiveEmpty(self):
        v = flow.EllipseOp(name = "Ellipse",
                           xchannel = "V2-A",
                           ychannel = "Y2-A").default_view(interactive = True)
        v.plot(self.ex)
        self.assertFalse(v._widget._selection_completed)
        self.assertIsNone(v._patch)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
