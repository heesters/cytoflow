'''
GUI tests for the proliferation and cell cycle operations.  (The models
themselves are tested in cytoflow/tests/test_proliferation.py; these check
the plumbing, on whatever channels the test data has.)
'''
import os, unittest, tempfile

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (ProliferationWorkflowOp, ProliferationWorkflowView,
                                             CellCycleWorkflowOp, CellCycleWorkflowView)
from cytoflowgui.workflow.subset import CategorySubset
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


class _Base(ImportedDataTest):
    
    def add(self, op):
        self.op = op
        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")
        self.view = view = wi.default_view
        wi.views.append(view)
        wi.current_view = view
        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        op.do_estimate = True
        self.workflow.wi_waitfor(wi, 'status', 'valid', timeout = 120)
        
    def notebook(self):
        code = "import cytoflow as flow\nfrom cytoflow import geom_mean, geom_sd\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
            for view in wi.views:
                code = code + view.get_notebook_code(i)
        ns = {}
        exec(code, ns)
        return ns['ex_{}'.format(len(self.workflow.workflow) - 1)]


class TestProliferation(_Base):

    def setUp(self):
        super().setUp()
        op = ProliferationWorkflowOp(name = "Prolif", channel = "Y2-A", by = ["Well"],
                                     num_generations = 4)
        op.control_subset_list.append(CategorySubset(name = "Well", values = ["A", "B"]))
        op.control_subset_list[0].selected = ["A"]
        self.add(op)

    def testEstimateApplyPlot(self):
        cols = self.workflow.remote_eval("list(self.workflow[-1].result.statistics['Prolif'].columns)")
        self.assertIn("Division_index", cols)
        self.assertIn("Generation_0", cols)
        self.assertTrue(self.workflow.remote_eval("'Prolif' in self.workflow[-1].result.conditions"))
        self.workflow.wi_waitfor(self.wi, 'view_error', '')
        self.assertEqual(self.workflow.remote_eval("self.workflow[-1].operation.control_subset"),
                         self.op.control_subset)
        self.assertNotEqual(self.op.control_subset, "")
        
        # changing a parameter clears the estimate
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.num_generations = 3
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.op.do_estimate = True
        self.workflow.wi_waitfor(self.wi, 'status', 'valid', timeout = 120)

    def testNotebook(self):
        self.assertIn("Prolif", self.notebook().statistics)


class TestCellCycle(_Base):

    def setUp(self):
        super().setUp()
        self.add(CellCycleWorkflowOp(name = "CC", channel = "FSC-A", by = ["Well"]))

    def testEstimateApplyPlot(self):
        cols = self.workflow.remote_eval("list(self.workflow[-1].result.statistics['CC'].columns)")
        self.assertEqual(cols[:3], ["G1", "S", "G2M"])
        self.assertTrue(self.workflow.remote_eval("'CC' in self.workflow[-1].result.conditions"))
        self.workflow.wi_waitfor(self.wi, 'view_error', '')

    def testNotebook(self):
        self.assertIn("CC", self.notebook().statistics)


class TestSerialize(unittest.TestCase):

    def _roundtrip(self, obj):
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(obj, filename)
            return load_yaml(filename)
        finally:
            os.unlink(filename)

    def testProliferation(self):
        op = ProliferationWorkflowOp(name = "P", channel = "Y2-A", by = ["Well"],
                                     num_generations = 5, undivided_peak = 1e5, fit_spacing = True)
        op.control_subset_list.append(CategorySubset(name = "Well", values = ["A", "B"]))
        new_op = self._roundtrip(op)
        self.assertIsInstance(new_op, ProliferationWorkflowOp)
        self.assertEqual((new_op.name, new_op.channel, new_op.by, new_op.num_generations,
                          new_op.undivided_peak, new_op.fit_spacing),
                         ("P", "Y2-A", ["Well"], 5, 1e5, True))
        self.assertEqual(len(new_op.control_subset_list), 1)
        self.assertIsInstance(self._roundtrip(op.default_view()), ProliferationWorkflowView)
        
    def testCellCycle(self):
        op = CellCycleWorkflowOp(name = "CC", channel = "FSC-A", by = ["Well"])
        new_op = self._roundtrip(op)
        self.assertIsInstance(new_op, CellCycleWorkflowOp)
        self.assertEqual((new_op.name, new_op.channel, new_op.by), ("CC", "FSC-A", ["Well"]))
        self.assertIsInstance(self._roundtrip(op.default_view()), CellCycleWorkflowView)


if __name__ == "__main__":
    unittest.main()
