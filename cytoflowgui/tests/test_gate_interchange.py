#!/usr/bin/env python3.11
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
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

'''
Importing and exporting gates in the GUI's workflow
'''

import os, unittest, tempfile, shutil

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (PolygonWorkflowOp, ThresholdWorkflowOp,
                                             EllipseWorkflowOp, BooleanGateWorkflowOp,
                                             Range2DWorkflowOp)
from cytoflowgui.workflow.serialization import load_yaml, save_yaml
from cytoflowgui.gate_interchange import import_gates, export_gates


class TestGateInterchange(ImportedDataTest):

    def setUp(self):
        super().setUp()
        self.tmp = tempfile.mkdtemp()

        poly = PolygonWorkflowOp()
        poly.name = "Poly"
        poly.xchannel = "Y2-A"
        poly.xscale = "logicle"
        poly.ychannel = "V2-A"
        poly.yscale = "logicle"
        poly.vertices = [(23.411982294776319, 5158.7027015021222),
                         (102.22182270573683, 23124.0588433874530),
                         (510.94519955277201, 23124.0588433874530),
                         (1089.5215641232173, 3800.3424832180476),
                         (340.56382570202402, 801.98947404942271),
                         (65.42597937575897, 1119.3133482602157)]

        thresh = ThresholdWorkflowOp()
        thresh.name = "Thresh"
        thresh.channel = "Y2-A"
        thresh.threshold = 500.0

        both = BooleanGateWorkflowOp()
        both.name = "Both"
        both.expression = "Poly and not Thresh"

        for op in (poly, thresh, both):
            wi = WorkflowItem(operation = op,
                              status = 'waiting',
                              view_error = "Not yet plotted")
            self.workflow.workflow.append(wi)
            self.workflow.wi_waitfor(wi, 'status', 'valid')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors = True)
        super().tearDown()

    def testExportImport(self):
        path = os.path.join(self.tmp, "gates.xml")
        exported, warns = export_gates(self.workflow, path)
        self.assertEqual(exported, ["Poly", "Thresh", "Both"])
        self.assertEqual(warns, [])

        # the names are taken, so the new gates get new ones
        n = len(self.workflow.workflow)
        new_items, warns = import_gates(self.workflow, path)
        self.assertEqual(warns, [])
        self.assertEqual(len(self.workflow.workflow), n + 3)
        self.assertEqual([type(wi.operation) for wi in new_items],
                         [PolygonWorkflowOp, ThresholdWorkflowOp, BooleanGateWorkflowOp])
        self.assertEqual([wi.operation.name for wi in new_items],
                         ["Poly_2", "Thresh_2", "Both_2"])
        self.assertEqual(new_items[0].operation.xscale, "logicle")
        self.assertEqual(len(new_items[0].operation.vertices), 6)
        self.assertEqual(new_items[2].operation.expression, "Poly_2 and not Thresh_2")
        self.assertIsNotNone(new_items[0].default_view)

        self.workflow.wi_waitfor(new_items[-1], 'status', 'valid')
        for name in ("Poly", "Thresh", "Both"):
            self.assertTrue(self.workflow.remote_eval(
                "(self.workflow[-1].result['{0}'] == self.workflow[-1].result['{0}_2']).all()"
                .format(name)))
        self.assertTrue(self.workflow.remote_eval("self.workflow[-1].result['Both_2'].sum() > 100"))

        # and they can be saved and loaded
        flow_path = os.path.join(self.tmp, "workflow.flow")
        save_yaml(self.workflow.workflow, flow_path)
        loaded = load_yaml(flow_path)
        self.assertEqual(loaded[-3].operation.vertices, new_items[0].operation.vertices)

    def testImportFlowJo(self):
        wsp = self.cwd + "/../../cytoflow/tests/data/flowjo/synthetic.wsp"

        with self.assertRaisesRegex(Exception, "specify which 'sample' or 'group'"):
            import_gates(self.workflow, wsp)

        new_items, warns = import_gates(self.workflow, wsp, sample = "RFP_Well_A3.fcs")
        text = "\n".join(warns)
        self.assertIn("CurlyQuad", text)
        self.assertIn("Compensation is NOT imported", text)

        kinds = {wi.operation.name : type(wi.operation) for wi in new_items}
        self.assertIs(kinds["Cells"], PolygonWorkflowOp)
        self.assertIs(kinds["Singlets"], Range2DWorkflowOp)
        self.assertIs(kinds["Singlets_gated"], BooleanGateWorkflowOp)
        self.assertNotIn(EllipseWorkflowOp, kinds.values())

        self.workflow.wi_waitfor(new_items[-1], 'status', 'valid', timeout = 120)
        self.assertTrue(self.workflow.remote_eval(
            "self.workflow[-1].result['Either_gated'].sum() > 100"))


if __name__ == "__main__":
    unittest.main()
