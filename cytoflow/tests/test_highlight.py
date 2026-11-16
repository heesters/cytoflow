'''
Tests for back-gating (the `highlight` attribute of ScatterplotView and
DensityView).
'''

import unittest

import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

import cytoflow as flow
import cytoflow.utility as util

from .test_base import ImportedDataTest


def _n_points(fig):
    return [len(c.get_offsets()) for ax in fig.axes for c in ax.collections
            if isinstance(c, PathCollection)]


class TestHighlight(ImportedDataTest):

    def setUp(self):
        super().setUp()
        self.ex = flow.RangeOp(name = "R", channel = "Y2-A", 
                               low = 2000, high = 20000).apply(self.ex)
        self.ex = flow.QuadOp(name = "Q", xchannel = "V2-A", ychannel = "Y2-A",
                              xthreshold = 100, ythreshold = 1000).apply(self.ex)
        
    def tearDown(self):
        plt.close('all')

    def testScatterCondition(self):
        flow.ScatterplotView(xchannel = "V2-A", ychannel = "Y2-A",
                             highlight = "R").plot(self.ex)
        self.assertEqual(_n_points(plt.gcf()), [len(self.ex), int(self.ex["R"].sum())])

    def testScatterExpression(self):
        flow.ScatterplotView(xchannel = "V2-A", ychannel = "Y2-A",
                             highlight = "Q == 'Q_2' and R").plot(self.ex, highlight_color = 'k')
        expected = int(((self.ex["Q"] == "Q_2") & self.ex["R"]).sum())
        self.assertEqual(_n_points(plt.gcf())[-1], expected)
        
    def testScatterSubsetFacet(self):
        flow.ScatterplotView(xchannel = "V2-A", ychannel = "Y2-A",
                             xfacet = "Dox", subset = "Dox == 10.0",
                             highlight = "R").plot(self.ex)
        expected = int(self.ex.data.query("Dox == 10.0")["R"].sum())
        self.assertEqual(sum(_n_points(plt.gcf())[1::2]), expected)

    def testDensity(self):
        flow.DensityView(xchannel = "V2-A", ychannel = "Y2-A", 
                         highlight = "R").plot(self.ex)
        self.assertEqual(_n_points(plt.gcf()), [int(self.ex["R"].sum())])
        
    def testNotBoolean(self):
        with self.assertRaises(util.CytoflowViewError):
            flow.ScatterplotView(xchannel = "V2-A", ychannel = "Y2-A",
                                 highlight = "Q").plot(self.ex)
            
    def testBadExpression(self):
        with self.assertRaises(util.CytoflowViewError):
            flow.ScatterplotView(xchannel = "V2-A", ychannel = "Y2-A",
                                 highlight = "NotAGate").plot(self.ex)


if __name__ == "__main__":
    unittest.main()
