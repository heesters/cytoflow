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
import numpy as np
import pandas as pd
from traits.util.async_trait_wait import wait_for_condition

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import DerivedParameterWorkflowOp
from cytoflowgui.workflow.operations import ChannelStatisticWorkflowOp
from cytoflowgui.workflow.serialization import load_yaml, save_yaml
from cytoflowgui.op_plugins.derived import available_channels_text

class TestDerivedParameter(ImportedDataTest):

    def setUp(self):
        super().setUp()

        self.addTypeEqualityFunc(DerivedParameterWorkflowOp, 'assertHasTraitsEqual')

        self.op = op = DerivedParameterWorkflowOp()
        op.name = "Derived"
        op.expression = "(Y2_A - `V2-A`) / 2"

        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")

        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        self.workflow.wi_waitfor(wi, 'status', "valid")

    def remote_data(self):
        return self.workflow.remote_eval("self.workflow[-1].result.data")

    def testApply(self):
        self.assertIsNotNone(self.workflow.remote_eval("self.workflow[-1].result"))
        self.assertIn("Derived", self.wi.channels)
        self.assertEqual(self.workflow.remote_eval("self.workflow[-1].result.metadata['Derived']['expression']"),
                         "(Y2_A - `V2-A`) / 2")

    def testValues(self):
        data = self.remote_data()
        np.testing.assert_allclose(data["Derived"],
                                   (data["Y2-A"] - data["V2-A"]) / 2)

    def testChangeExpression(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "arcsinh(B1_A / 150)"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

        data = self.remote_data()
        np.testing.assert_allclose(data["Derived"], np.arcsinh(data["B1-A"] / 150))

    def testNonFiniteDropped(self):
        n = self.workflow.remote_eval("len(self.workflow[-1].result)")

        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "log10(Y2_A)"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

        # (the warning about the dropped events is checked in the library
        # tests: it doesn't reach op_warning when warnings are filtered, eg
        # with pytest -W ignore)
        data = self.remote_data()
        self.assertLess(len(data), n)
        self.assertTrue(np.isfinite(data["Derived"]).all())

    def testBadExpression(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "Y2_A / Nope"
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.workflow.wi_waitfor(self.wi, 'op_error_trait', 'expression')
        wait_for_condition(lambda wi: "Nope" in wi.op_error, 
                           self.wi, 'op_error', 30)

        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "Y2_A.__class__"
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.workflow.wi_waitfor(self.wi, 'op_error_trait', 'expression')

        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.expression = "Y2_A * 2"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

    def testAvailableNames(self):
        text = available_channels_text(self.wi.previous_wi.channels)
        self.assertIn("Y2_A  (or `Y2-A`)", text)
        self.assertIn("log10", text)
        self.assertNotIn("Derived", text)

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
        op = DerivedParameterWorkflowOp()
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

        pd.testing.assert_frame_equal(nb_data, self.remote_data())


if __name__ == "__main__":
    unittest.main()
