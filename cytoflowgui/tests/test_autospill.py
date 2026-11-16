'''
GUI tests for the AutoSpill compensation operation
'''
import os, unittest, tempfile

from cytoflow import Tube
from cytoflowgui.tests.test_base import WorkflowTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (ImportWorkflowOp, AutoSpillWorkflowOp, 
                                             BleedthroughChannel)
from cytoflowgui.workflow.serialization import load_yaml, save_yaml

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                    "..", "..", "cytoflow", "tests", "data", "Plate01")

# single-color controls
CONTROLS = {"Y2-A" : "RFP_Well_A3.fcs",
            "V2-A" : "CFP_Well_A4.fcs",
            "B1-A" : "YFP_Well_A7.fcs"}


def make_op():
    op = AutoSpillWorkflowOp()
    op.channels_list = [BleedthroughChannel(channel = c, file = os.path.join(DATA, f))
                        for c, f in CONTROLS.items()]
    return op


class TestAutoSpill(WorkflowTest):

    def setUp(self):
        super().setUp()

        import_op = ImportWorkflowOp()
        import_op.tubes = [Tube(file = os.path.join(DATA, "RFP_Well_A6.fcs"))]
        wi = WorkflowItem(operation = import_op, status = "waiting")
        self.workflow.workflow.append(wi)
        import_op.do_estimate = True
        self.workflow.wi_waitfor(wi, 'status', 'valid')

        self.op = op = make_op()
        self.wi = wi = WorkflowItem(operation = op, status = 'waiting')
        wi.views.append(self.wi.default_view)
        self.workflow.workflow.append(wi)
        self.workflow.selected = self.wi

        op.do_estimate = True
        self.workflow.wi_waitfor(self.wi, 'status', 'valid', timeout = 120)

    def testEstimateApplyAndEdit(self):
        self.assertTrue(self.workflow.remote_eval("self.workflow[-1].result is not None"))

        # the matrix came back to the local process and filled in the table
        self.assertEqual(set(self.op.spillover.keys()), set(CONTROLS.keys()))
        self.assertEqual(len(self.op.spillover_list), 6)
        for c in CONTROLS:
            self.assertEqual(self.op.spillover[c][c], 1.0)
        self.assertTrue(self.workflow.remote_eval(
            "self.workflow[-1].result.metadata['Y2-A']['compensated']"))
        
        # editing a value re-applies with the edited matrix
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        s = next(x for x in self.op.spillover_list 
                 if x.from_channel == "Y2-A" and x.to_channel == "B1-A")
        s.spillover = 0.123
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertEqual(self.workflow.remote_eval(
            "self.workflow[-1].operation.spillover['Y2-A']['B1-A']"), 0.123)

    def testNotebook(self):
        code = "import cytoflow as flow\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
        ns = {}
        exec(code, ns)
        self.assertTrue(ns['ex_1'].metadata['Y2-A']['compensated'])


class TestAutoSpillSerialize(unittest.TestCase):

    def testSerialize(self):
        op = make_op()
        op.cofactor = 500.0
        op.spillover = {"Y2-A" : {"Y2-A" : 1.0, "V2-A" : 0.1},
                        "V2-A" : {"Y2-A" : 0.2, "V2-A" : 1.0}}
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)
        self.assertIsInstance(new_op, AutoSpillWorkflowOp)
        self.assertEqual([(c.channel, c.file) for c in new_op.channels_list],
                         [(c.channel, c.file) for c in op.channels_list])
        self.assertEqual(new_op.cofactor, 500.0)
        self.assertEqual(new_op.spillover, op.spillover)


if __name__ == "__main__":
    unittest.main()
