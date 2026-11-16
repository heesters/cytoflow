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
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util
from .test_base import ImportedDataSmallTest


class TestDerivedParameter(ImportedDataSmallTest):

    def setUp(self):
        super().setUp()
        self.y2 = self.ex.data["Y2-A"].values.copy()
        self.b1 = self.ex.data["B1-A"].values.copy()
        self.v2 = self.ex.data["V2-A"].values.copy()

    def apply(self, expression, name = "D"):
        # no events should be dropped, so no warning
        with warnings.catch_warnings():
            warnings.simplefilter("error", util.CytoflowOpWarning)
            return flow.DerivedParameterOp(name = name,
                                           expression = expression).apply(self.ex)

    def assertDerived(self, expression, expected):
        ex = self.apply(expression)
        self.assertIn("D", ex.channels)
        self.assertEqual(len(ex), len(self.ex))
        np.testing.assert_allclose(ex.data["D"].values, expected)
        return ex

    def testRatio(self):
        ex = self.assertDerived("Y2_A / B1_A", self.y2 / self.b1)

        # same answer as RatioOp
        ratio = flow.RatioOp(name = "R",
                             numerator = "Y2-A",
                             denominator = "B1-A").apply(self.ex)
        np.testing.assert_allclose(ex.data["D"].values, ratio.data["R"].values)

    def testArithmetic(self):
        self.assertDerived("(V2_A - Y2_A) * 2", (self.v2 - self.y2) * 2)
        self.assertDerived("-Y2_A + +B1_A", -self.y2 + self.b1)
        self.assertDerived("Y2_A ** 2 / 1e3", self.y2 ** 2 / 1e3)
        self.assertDerived("2 ** 3 * Y2_A", 8 * self.y2)

    def testFunctions(self):
        self.assertDerived("arcsinh(Y2_A / 150)", np.arcsinh(self.y2 / 150))
        self.assertDerived("asinh(Y2_A)", np.arcsinh(self.y2))
        self.assertDerived("abs(Y2_A)", np.abs(self.y2))
        self.assertDerived("exp(Y2_A / 1e5)", np.exp(self.y2 / 1e5))
        self.assertDerived("max(Y2_A, B1_A)", np.maximum(self.y2, self.b1))
        self.assertDerived("min(Y2_A, 0)", np.minimum(self.y2, 0))
        self.assertDerived("log10(abs(Y2_A) + 1)", np.log10(np.abs(self.y2) + 1))
        self.assertDerived("log2(abs(Y2_A) + 1)", np.log2(np.abs(self.y2) + 1))
        self.assertDerived("log(abs(Y2_A) + 1)", np.log(np.abs(self.y2) + 1))
        self.assertDerived("ln(abs(Y2_A) + 1)", np.log(np.abs(self.y2) + 1))
        self.assertDerived("sqrt(abs(Y2_A))", np.sqrt(np.abs(self.y2)))

    def testBacktickNames(self):
        self.assertDerived("`Y2-A` / `B1-A`", self.y2 / self.b1)

    def testExactName(self):
        ex = flow.DerivedParameterOp(name = "Y2_sq",
                                     expression = "`Y2-A` ** 2").apply(self.ex)
        ex = flow.DerivedParameterOp(name = "D",
                                     expression = "Y2_sq / 2").apply(ex)
        np.testing.assert_allclose(ex.data["D"].values, self.y2 ** 2 / 2)

    def testMetadata(self):
        ex = self.apply("Y2_A / B1_A")
        self.assertEqual(ex.metadata["D"]["type"], "channel")
        self.assertEqual(ex.metadata["D"]["expression"], "Y2_A / B1_A")
        self.assertAlmostEqual(ex.metadata["D"]["range"],
                               np.abs(self.y2 / self.b1).max())
        self.assertEqual(ex.history[-1].expression, "Y2_A / B1_A")

    def testNonFiniteDropped(self):
        n_bad = int((self.y2 <= 0).sum())
        self.assertGreater(n_bad, 0)

        with self.assertWarns(util.CytoflowOpWarning):
            ex = flow.DerivedParameterOp(name = "D",
                                         expression = "log10(Y2_A)").apply(self.ex)

        good = self.y2 > 0
        self.assertEqual(len(ex), len(self.ex) - n_bad)
        self.assertTrue(np.isfinite(ex.data["D"]).all())
        np.testing.assert_allclose(ex.data["D"].values, np.log10(self.y2[good]))

        # the other columns are dropped consistently
        np.testing.assert_allclose(ex.data["Y2-A"].values, self.y2[good])
        np.testing.assert_allclose(ex.data["B1-A"].values, self.b1[good])
        self.assertEqual(list(ex.data.index), list(range(len(ex))))

        # and the original experiment is untouched
        self.assertEqual(len(self.ex.data), len(self.y2))
        self.assertNotIn("D", self.ex.data)

    def testDivideByZero(self):
        self.ex.data.loc[0:9, "B1-A"] = 0.0
        with self.assertWarns(util.CytoflowOpWarning):
            ex = flow.DerivedParameterOp(name = "D",
                                         expression = "Y2_A / B1_A").apply(self.ex)
        self.assertEqual(len(ex), len(self.ex) - 10)
        self.assertTrue(np.isfinite(ex.data["D"]).all())

    def testAllNonFinite(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.DerivedParameterOp(name = "D",
                                    expression = "Y2_A / 0").apply(self.ex)

    def testBadName(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.DerivedParameterOp(expression = "Y2_A").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.DerivedParameterOp(name = "Bad-Name", expression = "Y2_A").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.DerivedParameterOp(name = "Dox", expression = "Y2_A").apply(self.ex)

    def assertRejected(self, expression):
        with self.assertRaises(util.CytoflowOpError) as cm:
            flow.DerivedParameterOp(name = "D", expression = expression).apply(self.ex)
        self.assertEqual(cm.exception.args[0], 'expression')
        return cm.exception.args[1]

    def testEmpty(self):
        self.assertRejected("")

    def testSyntaxError(self):
        self.assertRejected("Y2_A /")
        self.assertRejected("`Y2-A")

    def testUnknownChannel(self):
        msg = self.assertRejected("Y2_A / Nope")
        self.assertIn("Nope", msg)
        # the error lists the channels that are available
        self.assertIn("B1-A", msg)
        self.assertRejected("`Nope-A` + 1")

    def testConditionIsNotAChannel(self):
        self.assertRejected("Y2_A * Dox")

    def testNoChannels(self):
        self.assertRejected("1 + 2")

    def testBadFunctions(self):
        self.assertRejected("foo(Y2_A)")
        self.assertRejected("log(Y2_A, 2)")
        self.assertRejected("log(x = Y2_A)")
        self.assertRejected("max(Y2_A)")
        self.assertRejected("log + Y2_A")
        self.assertRejected("np.log(Y2_A)")

    def testUnsafe(self):
        for expression in ["__import__('os')",
                           "__import__('os').system('true')",
                           "Y2_A.__class__",
                           "Y2_A.sum()",
                           "Y2_A[0]",
                           "open('/etc/passwd')",
                           "eval('1')",
                           "lambda: Y2_A",
                           "[x for x in Y2_A]",
                           "Y2_A if Y2_A else B1_A",
                           "(Y2_A := B1_A)",
                           "Y2_A // 2",
                           "Y2_A % 2",
                           "Y2_A < 2",
                           "Y2_A and B1_A",
                           "'a' + Y2_A",
                           "True + Y2_A",
                           "1j * Y2_A",
                           "~Y2_A"]:
            with self.subTest(expression = expression):
                self.assertRejected(expression)

    def testHugeNumbers(self):
        self.assertRejected("Y2_A + 1" + "0" * 400)

        # overflow becomes inf, which is dropped -- here, every event, so
        # it's an error (not a hang or a crash)
        self.assertRejected("(abs(Y2_A) + 2) ** 1e10")

    def testDeeplyNested(self):
        self.assertRejected("-" * 100000 + "Y2_A")

    def testDefaultView(self):
        op = flow.DerivedParameterOp(name = "D", expression = "Y2_A / B1_A")
        ex = op.apply(self.ex)
        v = op.default_view(scale = "log")
        self.assertIsInstance(v, flow.HistogramView)
        self.assertEqual(v.channel, "D")
        v.plot(ex)
        plt.close('all')


if __name__ == "__main__":
    unittest.main()
