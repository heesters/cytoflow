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

import os, unittest, tempfile
import numpy as np
import pandas as pd

from cytoflowgui.tests.test_base import WorkflowTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (ImportWorkflowOp,
                                             CompensationWorkflowOp,
                                             CompensationWorkflowView,
                                             CompensationSpillover)
from cytoflowgui.workflow.serialization import load_yaml, save_yaml

class TestCompensation(WorkflowTest):

    def setUp(self):
        super().setUp()

        self.addTypeEqualityFunc(ImportWorkflowOp, 'assertHasTraitsEqual')
        self.addTypeEqualityFunc(CompensationWorkflowOp, 'assertHasTraitsEqual')
        self.addTypeEqualityFunc(CompensationWorkflowView, 'assertHasTraitsEqual')
        self.addTypeEqualityFunc(CompensationSpillover, 'assertHasTraitsEqual')

        self.cwd = os.path.dirname(os.path.abspath(__file__))

        # a BD FACSAria II file, with a non-trivial SPILL matrix
        from cytoflow import Tube
        import_op = ImportWorkflowOp()
        import_op.tubes = [Tube(file = self.cwd + "/../../cytoflow/tests/data/instruments/BD - FACS Aria II.fcs")]
        wi = WorkflowItem(operation = import_op,
                          status = "waiting",
                          view_error = "Not yet plotted")
        self.workflow.workflow.append(wi)
        import_op.do_estimate = True
        self.workflow.wi_waitfor(wi, 'status', 'valid')

        self.op = op = CompensationWorkflowOp()

        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")
        wi.views.append(self.wi.default_view)
        self.workflow.workflow.append(wi)
        self.workflow.selected = self.wi

        # read the matrix
        op.do_estimate = True
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

    def testEstimate(self):
        self.assertTrue(self.workflow.remote_eval("self.workflow[-1].result is not None"))

        # the matrix came back to the local process, and the table was filled in
        self.assertEqual(len(self.op.spillover), 9)
        self.assertAlmostEqual(self.op.spillover['B515-A']['V545-A'], 0.030264132476262202)
        self.assertEqual(len(self.op.spillover_list), 9 * 8)

        self.assertTrue(self.workflow.remote_eval("self.workflow[-1].result.metadata['B515-A']['compensated']"))

    def testEditValue(self):
        entry = next(s for s in self.op.spillover_list
                     if s.from_channel == 'B515-A' and s.to_channel == 'V545-A')

        before = self.workflow.remote_eval("self.workflow[-1].result.data['V545-A']")

        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        entry.spillover = 0.5
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

        self.assertAlmostEqual(self.workflow.remote_eval("self.workflow[-1].operation.spillover['B515-A']['V545-A']"), 0.5)
        self.assertAlmostEqual(self.op.spillover['B515-A']['V545-A'], 0.5)

        # the diagonal is kept
        self.assertEqual(self.workflow.remote_eval("self.workflow[-1].operation.spillover['B515-A']['B515-A']"), 1.0)

        after = self.workflow.remote_eval("self.workflow[-1].result.data['V545-A']")
        self.assertFalse(np.allclose(before, after))

    def testCSV(self):
        fh, filename = tempfile.mkstemp(suffix = '.csv')
        os.close(fh)
        try:
            with open(filename, 'w') as f:
                f.write("B515-A,V545-A\n1,0.25\n0.1,1\n")

            self.workflow.wi_sync(self.wi, 'status', 'waiting')
            self.op.spillover_file = filename
            self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
            self.assertEqual(len(self.op.spillover_list), 0)

            self.workflow.wi_sync(self.wi, 'status', 'waiting')
            self.op.do_estimate = True
            self.workflow.wi_waitfor(self.wi, 'status', 'valid')

            self.assertEqual(set(self.op.spillover.keys()), {'B515-A', 'V545-A'})
            self.assertAlmostEqual(self.op.spillover['B515-A']['V545-A'], 0.25)
            self.assertEqual(len(self.op.spillover_list), 2)
        finally:
            os.unlink(filename)

    def testPlot(self):
        self.workflow.wi_sync(self.wi, 'view_error', 'waiting')
        self.wi.current_view = self.wi.default_view
        self.workflow.wi_waitfor(self.wi, 'view_error', '')

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

        # the matrix itself is restored, too (it's not compared above
        # because it's a "status" trait)
        self.assertEqual(self.op.spillover, new_op.spillover)

    def testSerializeEmptyOp(self):
        fh, filename = tempfile.mkstemp()
        op = CompensationWorkflowOp()
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
        # edit a value first, to make sure the notebook gets the edited matrix
        entry = next(s for s in self.op.spillover_list
                     if s.from_channel == 'B515-A' and s.to_channel == 'V545-A')
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        entry.spillover = 0.5
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

        code = "import cytoflow as flow\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)

            for view in wi.views:
                code = code + view.get_notebook_code(i)

        ns = {}
        exec(code, ns)
        nb_data = ns['ex_1'].data
        remote_data = self.workflow.remote_eval("self.workflow[-1].result.data")

        pd.testing.assert_frame_equal(nb_data, remote_data)

if __name__ == "__main__":
    unittest.main()
