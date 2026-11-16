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
from cytoflow.operations._selectors import simplify_path, EasyLassoSelector
from .test_base import ImportedDataSmallTest


class Test(ImportedDataSmallTest):
    def setUp(self):
        super().setUp()
        self.gate = flow.FreehandOp(name = "Freehand",
                                    xchannel = "V2-A",
                                    ychannel = "Y2-A",
                                    vertices = [(-95.86, 12436.45),
                                                (116.29, 22530.75),
                                                (767.63, 4873.08),
                                                (101.64, 939.38),
                                                (-266.93, 2914.59)])

    def testGate(self):
        ex2 = self.gate.apply(self.ex)

        # the same as the equivalent polygon gate (see test_polygon.py)
        self.assertEqual(ex2.data.groupby("Freehand").size()[True], 4126)
        self.assertIsInstance(ex2.data.index, pd.RangeIndex)

        self.assertEqual(ex2.history[-1].id, 'cytoflow.operations.freehand')
        self.assertIsInstance(self.gate, flow.PolygonOp)

    def testPlot(self):
        self.gate.default_view().plot(self.ex)

    def testPlotParams(self):
        self.gate.default_view().plot(self.ex, patch_props = {'color' : 'blue'})

    def testDensityPlot(self):
        self.gate.default_view(density = True).plot(self.ex)

    def testInteractive(self):
        v = self.gate.default_view(interactive = True)
        v.plot(self.ex)
        self.assertIsInstance(v._widget, EasyLassoSelector)
        self.assertIsNotNone(v._patch)
        v.interactive = False
        self.assertIsNone(v._widget)

    def testSimplify(self):
        # a densely-sampled circle gets simplified to <= 100 vertices
        t = np.linspace(0, 2 * np.pi, 2000)
        xy = np.column_stack((200 * np.cos(t), 200 * np.sin(t)))
        idx = simplify_path(xy, tolerance = 1.0, max_vertices = 100)
        self.assertLessEqual(len(idx), 100)
        self.assertGreater(len(idx), 10)
        self.assertEqual(idx[0], 0)
        self.assertEqual(idx[-1], len(xy) - 1)
        self.assertTrue(np.all(np.diff(idx) > 0))

        idx = simplify_path(xy, tolerance = 1.0, max_vertices = 20)
        self.assertLessEqual(len(idx), 20)

        # collinear points collapse to the endpoints
        xy = np.column_stack((np.arange(10.0), np.zeros(10)))
        self.assertEqual(list(simplify_path(xy)), [0, 9])

        self.assertEqual(len(simplify_path(np.zeros((2, 2)))), 2)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
