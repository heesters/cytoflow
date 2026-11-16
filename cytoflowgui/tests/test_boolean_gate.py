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

import os, unittest, tempfile
import pandas as pd
from traits.util.async_trait_wait import wait_for_condition

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (ThresholdWorkflowOp, RangeWorkflowOp,
                                             QuadWorkflowOp, BooleanGateWorkflowOp)
from cytoflowgui.workflow.operations import ChannelStatisticWorkflowOp
from cytoflowgui.workflow.serialization import load_yaml, save_yaml
from cytoflowgui.op_plugins.boolean_gate import available_conditions_text

class TestBooleanGate(ImportedDataTest):

    def setUp(self):
        super().setUp()

        self.addTypeEqualityFunc(BooleanGateWorkflowOp, 'assertHasTraitsEqual')

        op = ThresholdWorkflowOp()
        op.name = "Y2_high"
        op.channel = "Y2-A"
        op.threshold = 300
        self.append(op)

        op = RangeWorkflowOp()
        op.name = "B1_mid"
        op.channel = "B1-A"
        op.low = 100
        op.high = 1000
        self.append(op)

        op = QuadWorkflowOp()
        op.name = "Q"
        op.xchannel = "V2-A"
        op.xthreshold = 100
        op.ychannel = "Y2-A"
        op.ythreshold = 100
        self.append(op)

        self.op = op = BooleanGateWorkflowOp()
        op.name = "Bool"
        op.expression = "Y2_high and not B1_mid"
        self.wi = self.append(op)

    def append(self, op):
        wi = WorkflowItem(operation = op,
                          status = 'waiting',
                          view_error = "Not yet plotted")

        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        self.workflow.wi_waitfor(wi, 'status', "valid")
        return wi

    def count(self, query):
        return self.workflow.remote_eval("len(self.workflow[-1].result.data.query(\"{}\"))"
                                         .format(query))

    def testApply(self):
        self.assertIsNotNone(self.workflow.remote_eval("self.workflow[-1].result"))
        self.assertIn("Bool", self.wi.conditions)

    def testEvents(self):
        n = self.count("Bool")
        self.assertGreater(n, 0)
        self.assertEqual(n, self.count("Y2_high and not B1_mid"))

    def testChangeExpression(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "Y2_high or Q == 'Q_2'"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

        self.assertEqual(self.count("Bool"),
                         self.count("Y2_high or Q == 'Q_2'"))

    def testBadExpression(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "Y2_high and Nope"
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.workflow.wi_waitfor(self.wi, 'op_error_trait', 'expression')
        wait_for_condition(lambda wi: "Nope" in wi.op_error, 
                           self.wi, 'op_error', 30)

        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "__import__('os')"
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.workflow.wi_waitfor(self.wi, 'op_error_trait', 'expression')

        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "not Y2_high"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertEqual(self.count("Bool"), self.count("not Y2_high"))

    def testAvailableNames(self):
        text = available_conditions_text(self.wi.previous_wi.conditions)
        self.assertIn("Y2_high", text.splitlines())
        self.assertIn("B1_mid", text.splitlines())
        self.assertIn("Q == 'Q_1'", text)
        self.assertIn("Dox (number)", text)
        self.assertNotIn("Bool", text)

        self.assertIn("make some gates", available_conditions_text({}))

    def testSerialize(self):
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
        op = BooleanGateWorkflowOp()
        try:
            os.close(fh)

            save_yaml(op, filename)
            new_op = load_yaml(filename)

        finally:
            os.unlink(filename)

        self.maxDiff = None

        self.assertEqual(op, new_op)

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
        # skip the statistics operations that ImportedDataTest adds: they
        # don't change the data (and their notebook code is tested elsewhere)
        ops = [wi.operation for wi in self.workflow.workflow
               if not isinstance(wi.operation, ChannelStatisticWorkflowOp)]

        code = "import cytoflow as flow\n"
        for i, op in enumerate(ops):
            code = code + op.get_notebook_code(i)

        namespace = {}
        exec(code, namespace)
        nb_data = namespace['ex_{}'.format(len(ops) - 1)].data
        remote_data = self.workflow.remote_eval("self.workflow[-1].result.data")

        pd.testing.assert_frame_equal(nb_data, remote_data)


if __name__ == "__main__":
    unittest.main()
