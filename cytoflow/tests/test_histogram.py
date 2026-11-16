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
Created on Mar 5, 2018

@author: brian
'''

import unittest
import cytoflow as flow

from .test_base import View1DTestBase, ImportedDataTest
import cytoflow.utility as util

class TestHistogram(View1DTestBase, unittest.TestCase):

    def setUp(self):
        super().setUp()
        self.view = flow.HistogramView(channel = "B1-A")

    # Histogram params
        
    def testNumBins(self):
        self.view.plot(self.ex, num_bins = 20)
        
    def testStep(self):
        self.view.plot(self.ex, histtype = 'step') 
        
    def testBar(self):
        self.view.plot(self.ex, histtype = 'bar')
        
    def testLineStyle(self):
        self.view.plot(self.ex, histtype = 'step', linestyle = 'solid')
        self.view.plot(self.ex, histtype = 'step', linestyle = 'dashed')
        self.view.plot(self.ex, histtype = 'step', linestyle = 'dashdot')
        self.view.plot(self.ex, histtype = 'step', linestyle = 'dotted')
        
    def testLineWidth(self):
        self.view.plot(self.ex, histtype = 'step', linestyle = 'solid', linewidth = 5)
        
    def testNormed(self):
        self.view.huefacet = "Dox"
        self.view.plot(self.ex, histtype = 'step', density = True)

        
if __name__ == "__main__":
    import sys;sys.argv = ['', 'TestHistogram.testLogicleScale']
    unittest.main()


class TestHistogramOverlays(ImportedDataTest):
    """normalize_to_mode and offset"""
    
    def setUp(self):
        super().setUp()
        self.view = flow.HistogramView(channel = "Y2-A", scale = "logicle",
                                       huefacet = "Dox")
        
    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        
    def _patches(self):
        import matplotlib.pyplot as plt
        return [p for p in plt.gca().patches]

    def testNormalizeToMode(self):
        import matplotlib.pyplot as plt
        self.view.plot(self.ex, normalize_to_mode = True, histtype = 'step')
        ax = plt.gca()
        self.assertEqual(ax.get_ylabel(), '% of max')
        for p in self._patches():
            self.assertAlmostEqual(p.get_path().vertices[:, 1].max(), 100.0)
            
    def testOffset(self):
        self.view.plot(self.ex, normalize_to_mode = True, offset = 0.5, histtype = 'step')
        bottoms = sorted(p.get_path().vertices[:, 1].min() for p in self._patches())
        self.assertEqual(len(bottoms), len(self.ex.conditions["Dox"]))
        for i, b in enumerate(bottoms):
            self.assertAlmostEqual(b, i * 50.0)
            
    def testOffsetCounts(self):
        self.view.plot(self.ex, offset = 0.5)
        
    def testOffsetHorizontal(self):
        self.view.plot(self.ex, offset = 0.5, orientation = 'horizontal', 
                       normalize_to_mode = True)
            
    def testDensityAndMode(self):
        with self.assertRaises(util.CytoflowViewError):
            self.view.plot(self.ex, density = True, normalize_to_mode = True)
