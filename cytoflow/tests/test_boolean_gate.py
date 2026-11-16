#!/usr/bin/env python3.11
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2025
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util
from cytoflow.operations.boolean_gate import boolean_condition_names
from .test_base import ImportedDataSmallTest


class TestBooleanGate(ImportedDataSmallTest):
    def setUp(self):
        super().setUp()
        self.ex = flow.ThresholdOp(name = "Y2_high",
                                   channel = "Y2-A",
                                   threshold = 300).apply(self.ex)
        self.ex = flow.RangeOp(name = "B1_mid",
                               channel = "B1-A",
                               low = 100,
                               high = 1000).apply(self.ex)
        self.ex = flow.QuadOp(name = "Q",
                              xchannel = "V2-A",
                              xthreshold = 100,
                              ychannel = "Y2-A",
                              ythreshold = 100).apply(self.ex)
        self.y2 = self.ex.data["Y2_high"]
        self.b1 = self.ex.data["B1_mid"]
        self.q = self.ex.data["Q"]

    def apply(self, expression):
        return flow.BooleanGateOp(name = "G",
                                  expression = expression).apply(self.ex)

    def assertGate(self, expression, expected):
        ex = self.apply(expression)
        self.assertEqual(ex.data["G"].dtype.kind, 'b')
        self.assertEqual(len(ex), len(self.ex))
        np.testing.assert_array_equal(ex.data["G"].values,
                                      np.asarray(expected, dtype = bool))
        return ex

    def testAnd(self):
        ex = self.assertGate("Y2_high and B1_mid", self.y2 & self.b1)
        # make sure the test is meaningful
        self.assertGreater(ex.data["G"].sum(), 0)
        self.assertLess(ex.data["G"].sum(), self.y2.sum())

    def testOr(self):
        self.assertGate("Y2_high or B1_mid", self.y2 | self.b1)

    def testNot(self):
        self.assertGate("not Y2_high", ~self.y2)

    def testParentheses(self):
        self.assertGate("not (Y2_high or B1_mid) or (Y2_high and B1_mid)",
                        ~(self.y2 | self.b1) | (self.y2 & self.b1))

    def testMultipleAnd(self):
        self.assertGate("Y2_high and B1_mid and Q == 'Q_2'",
                        self.y2 & self.b1 & (self.q == 'Q_2'))

    def testBitwiseOperators(self):
        self.assertGate("Y2_high & ~B1_mid | (Q != 'Q_1')",
                        (self.y2 & ~self.b1) | (self.q != 'Q_1'))

    def testCategorical(self):
        self.assertGate("Q == 'Q_2'", self.q == 'Q_2')
        self.assertGate("'Q_2' == Q", self.q == 'Q_2')
        self.assertGate("Q != 'Q_2' and Y2_high", (self.q != 'Q_2') & self.y2)

    def testBoolComparison(self):
        self.assertGate("Y2_high == False", ~self.y2)

    def testNumericComparison(self):
        dox = self.ex.data["Dox"]
        self.assertGate("Dox > 1.0", dox > 1.0)
        self.assertGate("1.0 < Dox", dox > 1.0)
        self.assertGate("Dox == 10", dox == 10)

    def testMetadata(self):
        ex = self.apply("Y2_high and B1_mid")
        self.assertEqual(ex.metadata["G"]["type"], "condition")
        self.assertEqual(ex.metadata["G"]["expression"], "Y2_high and B1_mid")
        self.assertEqual(ex.history[-1].expression, "Y2_high and B1_mid")
        self.assertIn("G", boolean_condition_names(ex))
        self.assertNotIn("Q", boolean_condition_names(ex))

        # chaining: use the new gate in another one
        ex2 = flow.BooleanGateOp(name = "H", expression = "not G").apply(ex)
        np.testing.assert_array_equal(ex2.data["H"].values, ~ex.data["G"].values)

    def testBadName(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.BooleanGateOp(expression = "Y2_high").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.BooleanGateOp(name = "Bad-Name", expression = "Y2_high").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.BooleanGateOp(name = "Y2_high", expression = "B1_mid").apply(self.ex)

    def assertRejected(self, expression):
        with self.assertRaises(util.CytoflowOpError) as cm:
            self.apply(expression)
        self.assertEqual(cm.exception.args[0], 'expression')
        return cm.exception.args[1]

    def testEmpty(self):
        self.assertRejected("")
        self.assertRejected("   ")

    def testSyntaxError(self):
        self.assertRejected("Y2_high and")
        self.assertRejected("(Y2_high")

    def testUnknownName(self):
        msg = self.assertRejected("Y2_high and Nope")
        self.assertIn("Nope", msg)
        # the error lists the gates that are available
        self.assertIn("B1_mid", msg)

    def testChannelIsNotACondition(self):
        self.assertRejected("Y2_high and V2_A")

    def testNonBooleanAlone(self):
        self.assertRejected("Q")
        self.assertRejected("Y2_high and Dox")

    def testUnknownCategory(self):
        msg = self.assertRejected("Q == 'Q_9'")
        self.assertIn("Q_2", msg)

    def testBadComparisons(self):
        self.assertRejected("Q < 'Q_2'")
        self.assertRejected("Y2_high > 1")
        self.assertRejected("Y2_high == 'yes'")
        self.assertRejected("Q == Y2_high")
        self.assertRejected("1 < Dox < 5")
        self.assertRejected("Q in ['Q_1', 'Q_2']")
        self.assertRejected("Q is 'Q_1'")

    def testUnsafe(self):
        for expression in ["__import__('os')",
                           "__import__('os').system('true')",
                           "Y2_high.__class__",
                           "Y2_high.sum()",
                           "Y2_high[0]",
                           "open('/etc/passwd')",
                           "lambda: Y2_high",
                           "[x for x in Y2_high]",
                           "Y2_high if B1_mid else B1_mid",
                           "(Y2_high := B1_mid)",
                           "Y2_high + B1_mid",
                           "f'{Y2_high}'",
                           "'Q_1'"]:
            with self.subTest(expression = expression):
                self.assertRejected(expression)

    def testDeeplyNested(self):
        self.assertRejected("not " * 100000 + "Y2_high")

    def testDefaultView(self):
        op = flow.BooleanGateOp(name = "G", expression = "Y2_high and B1_mid")
        ex = op.apply(self.ex)

        v = op.default_view(xchannel = "B1-A", ychannel = "Y2-A",
                            xscale = "log", yscale = "log")
        self.assertIsInstance(v, flow.ScatterplotView)
        self.assertEqual(v.huefacet, "G")
        v.plot(ex)
        plt.close('all')

        v = op.default_view(channel = "Y2-A", scale = "log")
        self.assertIsInstance(v, flow.HistogramView)
        self.assertEqual(v.huefacet, "G")
        v.plot(ex)
        plt.close('all')


if __name__ == "__main__":
    unittest.main()
