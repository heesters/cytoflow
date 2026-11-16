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
GUI tests for the statistical test operation.
'''

import os, unittest, tempfile
import pandas as pd

from cytoflow import geom_mean, geom_sd, geom_sem, ci

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (ChannelStatisticWorkflowOp,
                                             StatisticalTestWorkflowOp,
                                             StatisticalTestWorkflowView)
from cytoflowgui.workflow.subset import CategorySubset, RangeSubset
from cytoflowgui.workflow.views.table import TableWorkflowView
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


class TestStatisticalTest(ImportedDataTest):

    def setUp(self):
        super().setUp()

        self.addTypeEqualityFunc(StatisticalTestWorkflowOp, 'assertHasTraitsEqual')
        self.addTypeEqualityFunc(StatisticalTestWorkflowView, 'assertHasTraitsEqual')

        # the per-tube geometric mean: each Dox concentration has two
        # wells (replicates)
        op = ChannelStatisticWorkflowOp()
        op.name = "GeoMean"
        op.channel = "Y2-A"
        op.function_name = "Geo.Mean"
        op.by = ['Dox', 'Well']

        wi = WorkflowItem(operation = op)
        self.workflow.workflow.append(wi)

        self.op = op = StatisticalTestWorkflowOp()
        op.name = "DoxTest"
        op.statistic = "GeoMean"
        op.feature = "Y2-A"
        op.by = "Dox"
        op.test = "anova"
        op.posthoc = "tukey"

        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")
        wi.views.append(self.wi.default_view)
        self.view = wi.default_view

        self.workflow.workflow.append(wi)
        self.workflow.selected = wi

        self.workflow.wi_waitfor(wi, 'status', "valid")

    def testApply(self):
        self.assertIsNotNone(self.workflow.remote_eval("self.workflow[-1].result"))
        stats = self.workflow.remote_eval("list(self.workflow[-1].result.statistics.keys())")
        self.assertIn("DoxTest", stats)
        self.assertIn("DoxTest_Groups", stats)
        self.assertIn("DoxTest_Tukey", stats)

        # the results are available to views downstream
        self.assertIn("DoxTest_Tukey", self.wi.statistics)
        tukey = self.workflow.remote_eval("self.workflow[-1].result.statistics['DoxTest_Tukey']")
        self.assertEqual(len(tukey), 3)

    def testChangeTest(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.test = "kruskal"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        stat = self.workflow.remote_eval("self.workflow[-1].result.statistics['DoxTest']")
        self.assertEqual(stat.index[0], "Kruskal-Wallis")

    def testNoPosthoc(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.posthoc = "none"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        stats = self.workflow.remote_eval("list(self.workflow[-1].result.statistics.keys())")
        self.assertNotIn("DoxTest_Tukey", stats)

    def testTwoGroupTest(self):
        # three groups -- an error
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.test = "welch"
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.assertIn("exactly two groups", self.workflow.remote_eval("self.workflow[-1].op_error"))

        # pick two
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.subset_list.append(RangeSubset(name = "Dox", values = [1.0, 10.0, 100.0]))
        self.op.subset_list[0].low = 1.0
        self.op.subset_list[0].high = 10.0
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        stat = self.workflow.remote_eval("self.workflow[-1].result.statistics['DoxTest']")
        self.assertEqual(stat["n_groups"].iloc[0], 2)

    def testTooFewReplicates(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.subset_list.append(CategorySubset(name = "Well", values = ["A", "B"]))
        self.op.subset_list[0].selected = ["A"]
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.assertIn("at least two replicates", self.workflow.remote_eval("self.workflow[-1].op_error"))

    def testEvents(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.observations = "events"
        self.op.channel = "Y2-A"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        stat = self.workflow.remote_eval("self.workflow[-1].result.statistics['DoxTest']")
        n = self.workflow.remote_eval("len(self.workflow[-1].result)")
        self.assertEqual(stat["n_obs"].iloc[0], n)

    def testPlot(self):
        self.workflow.wi_sync(self.wi, 'view_error', 'waiting')
        self.wi.current_view = self.view
        self.workflow.wi_waitfor(self.wi, 'view_error', '')

        self.workflow.wi_sync(self.wi, 'view_error', 'waiting')
        self.view.show_ns = True
        self.workflow.wi_waitfor(self.wi, 'view_error', '')

        self.workflow.wi_sync(self.wi, 'view_error', 'waiting')
        self.op.test = "kruskal"
        self.workflow.wi_waitfor(self.wi, 'view_error', '')

    def testTableView(self):
        # the pairwise comparisons, as a group1 x group2 table
        view = TableWorkflowView(statistic = "DoxTest_Tukey",
                                 feature = "p_adj",
                                 row_facet = "group1",
                                 column_facet = "group2")
        self.workflow.wi_sync(self.wi, 'view_error', 'waiting')
        self.wi.views.append(view)
        self.wi.current_view = view
        self.workflow.wi_waitfor(self.wi, 'view_error', '')

        # the overall test (a single row)
        self.workflow.wi_sync(self.wi, 'view_error', 'waiting')
        view.trait_set(statistic = "DoxTest", feature = "p_value",
                       row_facet = "Test", column_facet = "")
        self.workflow.wi_waitfor(self.wi, 'view_error', '')

    def testSerialize(self):
        self.op.subset_list.append(CategorySubset(name = "Well", values = ["A", "B"]))
        self.op.alpha = 0.01
        self.op.observations = "events"
        self.op.channel = "Y2-A"

        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(self.op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)

        self.maxDiff = None
        self.assertEqual(self.op, new_op)

    def testSerializeEmptyOp(self):
        fh, filename = tempfile.mkstemp()
        op = StatisticalTestWorkflowOp()
        try:
            os.close(fh)
            save_yaml(op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)

        self.maxDiff = None
        self.assertEqual(op, new_op)

    def testSerializeView(self):
        self.view.show_ns = True

        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(self.view, filename)
            new_view = load_yaml(filename)
        finally:
            os.unlink(filename)

        self.maxDiff = None
        self.assertEqual(self.view, new_view)

    def testSerializeWorkflowItem(self):
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(self.wi, filename)
            new_wi = load_yaml(filename)
        finally:
            os.unlink(filename)

        self.maxDiff = None
        self.assertEqual(self.wi, new_wi)

    def testNotebook(self):
        code = "import cytoflow as flow\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
            for view in wi.views:
                code = code + view.get_notebook_code(i)

        self.assertIn("default_view", code)

        # the earlier operations' notebook code refers to these directly
        namespace = {"geom_mean" : geom_mean, "geom_sd" : geom_sd,
                     "geom_sem" : geom_sem, "ci" : ci}
        exec(code, namespace)

        n = len(self.workflow.workflow) - 1
        for name in ["DoxTest", "DoxTest_Groups", "DoxTest_Tukey"]:
            nb_stat = namespace['ex_{}'.format(n)].statistics[name]
            remote_stat = self.workflow.remote_eval("self.workflow[-1].result.statistics['{}']".format(name))
            pd.testing.assert_frame_equal(nb_stat, remote_stat)


if __name__ == "__main__":
#     import sys;sys.argv = ['', 'TestStatisticalTest.testNotebook']
    unittest.main()
