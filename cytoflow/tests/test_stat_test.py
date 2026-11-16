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
Tests for cytoflow.operations.stat_test -- check the results against
scipy / statsmodels called directly on the same data.
'''

import unittest
import warnings

import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import cytoflow as flow
import cytoflow.utility as util
from .test_base import ImportedDataTest


class TestStatisticalTest(ImportedDataTest):

    def setUp(self):
        super().setUp()

        # a statistic with known values: three "Dox" groups, with 4, 4 and 3
        # replicates (unequal group sizes on purpose), indexed by Dox and Rep
        self.groups = {0.0 : np.array([10.1, 11.3, 9.8, 10.9]),
                       10.0 : np.array([12.2, 13.1, 12.8, 11.9]),
                       100.0 : np.array([15.3, 14.1, 16.2])}

        idx = []
        values = []
        for dox, v in self.groups.items():
            for i, x in enumerate(v):
                idx.append((dox, i + 1))
                values.append(x)

        self.ex.statistics['Known'] = \
            pd.DataFrame({"Y2-A" : values,
                          "Other" : np.arange(len(values), dtype = float)},
                         index = pd.MultiIndex.from_tuples(idx, names = ["Dox", "Rep"]))

        # and a "real" per-tube statistic, from the data
        self.ex = flow.ChannelStatisticOp(name = "GeoMean",
                                          channel = "Y2-A",
                                          function = flow.geom_mean,
                                          by = ["Dox", "Well"]).apply(self.ex)

    def tearDown(self):
        plt.close('all')

    def _op(self, **kwargs):
        args = dict(name = "T", statistic = "Known", feature = "Y2-A", by = "Dox")
        args.update(kwargs)
        return flow.StatisticalTestOp(**args)

    def testANOVA(self):
        ex = self._op(test = "anova").apply(self.ex)
        stat = ex.statistics["T"]

        expected = scipy.stats.f_oneway(*self.groups.values())
        self.assertEqual(stat.index.names, ["Test"])
        self.assertEqual(stat.index[0], "One-way ANOVA")
        self.assertEqual(stat["statistic"].iloc[0], expected.statistic)
        self.assertEqual(stat["p_value"].iloc[0], expected.pvalue)
        self.assertEqual(stat["df_between"].iloc[0], 2)
        self.assertEqual(stat["df_within"].iloc[0], 11 - 3)
        self.assertEqual(stat["n_groups"].iloc[0], 3)
        self.assertEqual(stat["n_obs"].iloc[0], 11)

        # eta-squared
        allv = np.concatenate(list(self.groups.values()))
        ss_b = sum(len(v) * (v.mean() - allv.mean()) ** 2 for v in self.groups.values())
        ss_t = ((allv - allv.mean()) ** 2).sum()
        self.assertAlmostEqual(stat["effect_size"].iloc[0], ss_b / ss_t, places = 12)

        # also, F = (SSB / dfb) / (SSW / dfw)
        self.assertAlmostEqual(stat["statistic"].iloc[0],
                               (ss_b / 2) / ((ss_t - ss_b) / 8), places = 10)

    def testGroups(self):
        ex = self._op().apply(self.ex)
        desc = ex.statistics["T_Groups"]
        self.assertEqual(desc.index.name, "Dox")
        self.assertEqual(list(desc.index), [0.0, 10.0, 100.0])
        for g, v in self.groups.items():
            self.assertEqual(desc.loc[g, "n"], len(v))
            self.assertAlmostEqual(desc.loc[g, "mean"], v.mean(), places = 12)
            self.assertAlmostEqual(desc.loc[g, "sd"], v.std(ddof = 1), places = 12)
            self.assertAlmostEqual(desc.loc[g, "sem"], scipy.stats.sem(v), places = 12)
            self.assertAlmostEqual(desc.loc[g, "median"], np.median(v), places = 12)

    def testKruskal(self):
        ex = self._op(test = "kruskal").apply(self.ex)
        stat = ex.statistics["T"]

        expected = scipy.stats.kruskal(*self.groups.values())
        self.assertEqual(stat.index[0], "Kruskal-Wallis")
        self.assertEqual(stat["statistic"].iloc[0], expected.statistic)
        self.assertEqual(stat["p_value"].iloc[0], expected.pvalue)
        self.assertEqual(stat["df_between"].iloc[0], 2)
        self.assertTrue(np.isnan(stat["df_within"].iloc[0]))
        self.assertAlmostEqual(stat["effect_size"].iloc[0], expected.statistic / 10, places = 12)

    def testWelch(self):
        ex = self._op(test = "welch", subset = "Dox != 100").apply(self.ex)
        stat = ex.statistics["T"]

        a, b = self.groups[0.0], self.groups[10.0]
        expected = scipy.stats.ttest_ind(a, b, equal_var = False)
        self.assertEqual(stat.index[0], "Welch's t-test")
        self.assertEqual(stat["statistic"].iloc[0], expected.statistic)
        self.assertEqual(stat["p_value"].iloc[0], expected.pvalue)

        # Welch-Satterthwaite degrees of freedom
        va, vb = a.var(ddof = 1) / len(a), b.var(ddof = 1) / len(b)
        df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
        self.assertAlmostEqual(stat["df_between"].iloc[0], df, places = 10)

        d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof = 1) + b.var(ddof = 1)) / 2)
        self.assertAlmostEqual(stat["effect_size"].iloc[0], d, places = 12)
        self.assertEqual(stat["n_obs"].iloc[0], 8)

    def testMannWhitney(self):
        ex = self._op(test = "mannwhitney", subset = "Dox != 0").apply(self.ex)
        stat = ex.statistics["T"]

        a, b = self.groups[10.0], self.groups[100.0]
        expected = scipy.stats.mannwhitneyu(a, b, alternative = 'two-sided')
        self.assertEqual(stat.index[0], "Mann-Whitney U")
        self.assertEqual(stat["statistic"].iloc[0], expected.statistic)
        self.assertEqual(stat["p_value"].iloc[0], expected.pvalue)
        self.assertAlmostEqual(stat["effect_size"].iloc[0],
                               2 * expected.statistic / (len(a) * len(b)) - 1,
                               places = 12)

    def testTukey(self):
        ex = self._op(test = "anova", posthoc = "tukey").apply(self.ex)
        tukey = ex.statistics["T_Tukey"]

        endog = np.concatenate(list(self.groups.values()))
        labels = np.concatenate([[g] * len(v) for g, v in self.groups.items()])
        expected = pairwise_tukeyhsd(endog, labels, alpha = 0.05)

        self.assertEqual(tukey.index.names, ["group1", "group2"])
        self.assertEqual(list(tukey.index),
                         [(0.0, 10.0), (0.0, 100.0), (10.0, 100.0)])
        np.testing.assert_array_equal(tukey["meandiff"].to_numpy(), expected.meandiffs)
        np.testing.assert_array_equal(tukey["p_adj"].to_numpy(), expected.pvalues)
        np.testing.assert_array_equal(tukey["lower"].to_numpy(), expected.confint[:, 0])
        np.testing.assert_array_equal(tukey["upper"].to_numpy(), expected.confint[:, 1])
        np.testing.assert_array_equal(tukey["reject"].to_numpy(), expected.reject)

        # meandiff is group2 - group1
        self.assertAlmostEqual(tukey.loc[(0.0, 100.0), "meandiff"],
                               self.groups[100.0].mean() - self.groups[0.0].mean(),
                               places = 12)

        # cross-check against scipy's implementation of the same test
        sp = scipy.stats.tukey_hsd(*self.groups.values())
        np.testing.assert_allclose(tukey["p_adj"].to_numpy(),
                                   [sp.pvalue[0, 1], sp.pvalue[0, 2], sp.pvalue[1, 2]],
                                   rtol = 1e-6)

        for p, s in zip(tukey["p_adj"], tukey["significance"]):
            self.assertEqual(s, util_stars(p))

    def testTukeyAlpha(self):
        ex = self._op(posthoc = "tukey", alpha = 0.01).apply(self.ex)
        tukey = ex.statistics["T_Tukey"]

        endog = np.concatenate(list(self.groups.values()))
        labels = np.concatenate([[g] * len(v) for g, v in self.groups.items()])
        expected = pairwise_tukeyhsd(endog, labels, alpha = 0.01)
        np.testing.assert_array_equal(tukey["lower"].to_numpy(), expected.confint[:, 0])
        np.testing.assert_array_equal(tukey["reject"].to_numpy(), expected.reject)

    def testReplicatesFromData(self):
        # the per-tube geometric means; wells are the replicates
        ex = flow.StatisticalTestOp(name = "T",
                                    statistic = "GeoMean",
                                    feature = "Y2-A",
                                    by = "Dox",
                                    posthoc = "tukey").apply(self.ex)

        gm = self.ex.statistics["GeoMean"]["Y2-A"]
        groups = [gm.xs(d, level = "Dox").to_numpy() for d in [0.0, 10.0, 100.0]]
        expected = scipy.stats.f_oneway(*groups)
        self.assertEqual(ex.statistics["T"]["statistic"].iloc[0], expected.statistic)
        self.assertEqual(ex.statistics["T"]["p_value"].iloc[0], expected.pvalue)
        self.assertEqual(ex.statistics["T"]["n_obs"].iloc[0], 6)

    def testEvents(self):
        ex = flow.StatisticalTestOp(name = "T",
                                    observations = "events",
                                    channel = "Y2-A",
                                    by = "Dox",
                                    posthoc = "tukey").apply(self.ex)

        groups = [self.ex.data.loc[self.ex.data["Dox"] == d, "Y2-A"].to_numpy()
                  for d in [0.0, 10.0, 100.0]]
        expected = scipy.stats.f_oneway(*groups)
        stat = ex.statistics["T"]
        self.assertEqual(stat["statistic"].iloc[0], expected.statistic)
        self.assertEqual(stat["p_value"].iloc[0], expected.pvalue)
        self.assertEqual(stat["n_obs"].iloc[0], len(self.ex))
        self.assertIn("T_Tukey", ex.statistics)

    def testEventsSubset(self):
        ex = flow.StatisticalTestOp(name = "T",
                                    observations = "events",
                                    channel = "Y2-A",
                                    by = "Dox",
                                    test = "mannwhitney",
                                    subset = "Dox != 100").apply(self.ex)

        a = self.ex.data.loc[self.ex.data["Dox"] == 0.0, "Y2-A"].to_numpy()
        b = self.ex.data.loc[self.ex.data["Dox"] == 10.0, "Y2-A"].to_numpy()
        expected = scipy.stats.mannwhitneyu(a, b, alternative = 'two-sided')
        self.assertEqual(ex.statistics["T"]["statistic"].iloc[0], expected.statistic)
        self.assertEqual(ex.statistics["T"]["p_value"].iloc[0], expected.pvalue)

    def testOriginalUnchanged(self):
        ex = self._op(posthoc = "tukey").apply(self.ex)
        self.assertNotIn("T", self.ex.statistics)
        self.assertIn("T", ex.statistics)
        self.assertIn("T_Groups", ex.statistics)
        self.assertIn("T_Tukey", ex.statistics)
        self.assertNotIn("T_Tukey", self._op(name = "U").apply(self.ex).statistics)

    def testNaN(self):
        self.ex.statistics["Known"].loc[(0.0, 1), "Y2-A"] = np.nan
        with self.assertWarns(util.CytoflowOpWarning):
            ex = self._op().apply(self.ex)
        self.assertEqual(ex.statistics["T"]["n_obs"].iloc[0], 10)

    ### errors

    def testBadName(self):
        with self.assertRaises(util.CytoflowOpError):
            self._op(name = "").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(name = "Bad Name!").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(name = "GeoMean").apply(self.ex)

        ex = self._op(name = "X").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(name = "X").apply(ex)

    def testBadStatistic(self):
        with self.assertRaises(util.CytoflowOpError):
            self._op(statistic = "").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(statistic = "NotAStat").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(feature = "").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(feature = "NotAFeature").apply(self.ex)

    def testBadBy(self):
        with self.assertRaises(util.CytoflowOpError):
            self._op(by = "").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(by = "Well").apply(self.ex)

    def testBadSubset(self):
        with self.assertRaises(util.CytoflowOpError):
            self._op(subset = "Dox ==").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(subset = "Dox == 12345").apply(self.ex)

    def testBadAlpha(self):
        with self.assertRaises(util.CytoflowOpError):
            self._op(alpha = 0.0).apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            self._op(alpha = 1.5).apply(self.ex)

    def testTooFewReplicates(self):
        # one value per group -- no replicates
        ex = flow.ChannelStatisticOp(name = "ByDox",
                                     channel = "Y2-A",
                                     function = flow.geom_mean,
                                     by = ["Dox"]).apply(self.ex)

        with self.assertRaisesRegex(util.CytoflowOpError, "at least two replicates"):
            flow.StatisticalTestOp(name = "T",
                                   statistic = "ByDox",
                                   feature = "Y2-A",
                                   by = "Dox").apply(ex)

        # one group with only one replicate
        with self.assertRaisesRegex(util.CytoflowOpError, "at least two replicates"):
            self._op(subset = "not (Dox == 100 and Rep > 1)").apply(self.ex)

    def testTooFewGroups(self):
        with self.assertRaisesRegex(util.CytoflowOpError, "at least two groups"):
            self._op(subset = "Dox == 0").apply(self.ex)

    def testTwoGroupTests(self):
        for test in ["welch", "mannwhitney"]:
            with self.assertRaisesRegex(util.CytoflowOpError, "exactly two groups"):
                self._op(test = test).apply(self.ex)

    def testBadEvents(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.StatisticalTestOp(name = "T", observations = "events",
                                   by = "Dox").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.StatisticalTestOp(name = "T", observations = "events",
                                   channel = "NotAChannel", by = "Dox").apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.StatisticalTestOp(name = "T", observations = "events",
                                   channel = "Y2-A", by = "NotACondition").apply(self.ex)

    ### views

    def testView(self):
        op = self._op(posthoc = "tukey")
        ex = op.apply(self.ex)
        op.default_view().plot(ex)
        op.default_view(show_ns = True).plot(ex)

    def testViewTwoGroups(self):
        op = self._op(test = "welch", subset = "Dox != 100")
        ex = op.apply(self.ex)
        op.default_view().plot(ex)

    def testViewEvents(self):
        op = flow.StatisticalTestOp(name = "T", observations = "events",
                                    channel = "Y2-A", by = "Dox",
                                    test = "kruskal")
        ex = op.apply(self.ex)
        op.default_view().plot(ex)

    def testViewNotApplied(self):
        op = self._op()
        with self.assertRaises(util.CytoflowViewError):
            op.default_view().plot(self.ex)

    def testTableView(self):
        op = self._op(posthoc = "tukey")
        ex = op.apply(self.ex)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for feature in ["meandiff", "p_adj", "reject", "significance"]:
                flow.TableView(statistic = "T_Tukey",
                               feature = feature,
                               row_facet = "group1",
                               column_facet = "group2").plot(ex)

            flow.TableView(statistic = "T",
                           feature = "p_value",
                           row_facet = "Test").plot(ex)

            flow.TableView(statistic = "T_Groups",
                           feature = "mean",
                           row_facet = "Dox").plot(ex)


def util_stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


if __name__ == "__main__":
    unittest.main()
