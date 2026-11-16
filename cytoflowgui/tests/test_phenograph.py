'''
GUI tests for graph clustering (PhenographOp)
'''
import os, unittest, tempfile

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (PhenographWorkflowOp, 
                                             PhenographWorkflowView, tSNEChannel)
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


def make_op():
    op = PhenographWorkflowOp()
    op.name = "Clust"
    op.channels_list = [tSNEChannel(channel = "V2-A", scale = "logicle"),
                        tSNEChannel(channel = "Y2-A", scale = "logicle"),
                        tSNEChannel(channel = "B1-A", scale = "logicle")]
    op.num_neighbors = 15
    op.sample = 0.05
    return op


class TestPhenograph(ImportedDataTest):

    def setUp(self):
        super().setUp()

        self.op = op = make_op()
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

    def testEstimateApplyPlot(self):
        self.assertTrue(self.workflow.remote_eval("'Clust' in self.workflow[-1].result.conditions"))
        self.assertTrue(self.workflow.remote_eval("'Clust' in self.workflow[-1].result.statistics"))
        self.workflow.wi_waitfor(self.wi, 'view_error', '')
        
        # changing a parameter clears the estimate
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.resolution = 0.5
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.op.do_estimate = True
        self.workflow.wi_waitfor(self.wi, 'status', 'valid', timeout = 120)

    def testNotebook(self):
        code = "import cytoflow as flow\nfrom cytoflow import geom_mean, geom_sd\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
            for view in wi.views:
                code = code + view.get_notebook_code(i)
        ns = {}
        exec(code, ns)
        self.assertIn("Clust", ns['ex_{}'.format(len(self.workflow.workflow) - 1)].conditions)


class TestPhenographSerialize(unittest.TestCase):

    def _roundtrip(self, obj):
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(obj, filename)
            return load_yaml(filename)
        finally:
            os.unlink(filename)

    def testSerializeOp(self):
        op = make_op()
        new_op = self._roundtrip(op)
        self.assertIsInstance(new_op, PhenographWorkflowOp)
        self.assertEqual(new_op.channels, op.channels)
        self.assertEqual(new_op.scale, op.scale)
        self.assertEqual((new_op.num_neighbors, new_op.sample, new_op.method), (15, 0.05, "leiden"))
        
    def testSerializeView(self):
        new_view = self._roundtrip(make_op().default_view())
        self.assertIsInstance(new_view, PhenographWorkflowView)


if __name__ == "__main__":
    unittest.main()
