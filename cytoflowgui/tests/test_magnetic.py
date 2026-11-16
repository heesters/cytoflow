'''
GUI tests for "magnetic" (per-sample adjusted) gates
'''
import os, unittest, tempfile

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import PolygonWorkflowOp, EllipseWorkflowOp
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


class TestMagnetic(ImportedDataTest):

    def setUp(self):
        super().setUp()

        self.op = op = PolygonWorkflowOp()
        op.name = "Poly"
        op.xchannel = "Y2-A"
        op.ychannel = "V2-A"
        op.vertices = [(0, 0), (20000, 0), (20000, 2000), (0, 2000)]

        self.wi = wi = WorkflowItem(operation = op, status = 'waiting')
        self.workflow.workflow.append(wi)
        self.workflow.wi_waitfor(wi, 'status', "valid")

    def testAdjustBy(self):
        self.assertFalse(self.workflow.remote_eval("'Poly' in self.workflow[-1].result.statistics"))
        
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.adjust_by = ["Well"]
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        
        cols = self.workflow.remote_eval("list(self.workflow[-1].result.statistics['Poly'].columns)")
        self.assertEqual(cols, ['x_shift', 'y_shift'])
        
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.max_adjust = 0.1
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

    def testNotebook(self):
        self.op.adjust_by = ["Well"]
        code = "import cytoflow as flow\nfrom cytoflow import geom_mean, geom_sd\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
        ns = {}
        exec(code, ns)
        ex = ns['ex_{}'.format(len(self.workflow.workflow) - 1)]
        self.assertIn("Poly", ex.statistics)


class TestMagneticSerialize(unittest.TestCase):
    
    def testSerialize(self):
        self.op = PolygonWorkflowOp(name = "Poly", xchannel = "Y2-A", ychannel = "V2-A",
                                    vertices = [(0, 0), (20000, 0), (20000, 2000), (0, 2000)],
                                    adjust_by = ["Well"], max_adjust = 0.5)
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(self.op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)
        self.assertEqual(new_op.adjust_by, ["Well"])
        self.assertEqual(new_op.max_adjust, 0.5)
        self.assertEqual(new_op.vertices, self.op.vertices)
        
    def testSerializeEllipse(self):
        op = EllipseWorkflowOp(name = "E", xchannel = "Y2-A", ychannel = "V2-A",
                               adjust_by = ["Well"], max_adjust = 0.5)
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)
        self.assertEqual(new_op.adjust_by, ["Well"])
        self.assertEqual(new_op.max_adjust, 0.5)


if __name__ == "__main__":
    unittest.main()
