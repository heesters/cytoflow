'''
GUI tests for the cross entropy test
'''
import os, unittest, tempfile

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (CrossEntropyTestWorkflowOp, 
                                             CrossEntropyTestWorkflowView, tSNEChannel)
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


def make_op():
    op = CrossEntropyTestWorkflowOp()
    op.name = "CE"
    op.channels_list = [tSNEChannel(channel = "V2-A", scale = "logicle"),
                        tSNEChannel(channel = "Y2-A", scale = "logicle"),
                        tSNEChannel(channel = "B1-A", scale = "logicle")]
    # any two channels will do as the "embedding" for these tests
    op.xembedding = "FSC-A"
    op.yembedding = "SSC-A"
    op.by = "Well"
    op.perplexity = 20
    op.max_events = 1500
    return op


class TestCrossEntropy(ImportedDataTest):

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
        cols = self.workflow.remote_eval("list(self.workflow[-1].result.statistics['CE'].columns)")
        self.assertEqual(cols, ['D', 'p', 'p_adj', 'reject'])
        self.assertTrue(self.workflow.remote_eval("'CE_Groups' in self.workflow[-1].result.statistics"))
        self.workflow.wi_waitfor(self.wi, 'view_error', '')
        
        # changing a parameter clears the estimate
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.perplexity = 25
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
        self.op.do_estimate = True
        self.workflow.wi_waitfor(self.wi, 'status', 'valid', timeout = 120)

    def testNotebook(self):
        code = "import cytoflow as flow\nfrom cytoflow import geom_mean, geom_sd\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
        ns = {}
        exec(code, ns)
        self.assertIn("CE", ns['ex_{}'.format(len(self.workflow.workflow) - 1)].statistics)


class TestCrossEntropySerialize(unittest.TestCase):

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
        self.assertIsInstance(new_op, CrossEntropyTestWorkflowOp)
        self.assertEqual(new_op.channels, op.channels)
        self.assertEqual(new_op.scale, op.scale)
        self.assertEqual(new_op.embedding, ["FSC-A", "SSC-A"])
        self.assertEqual((new_op.by, new_op.perplexity, new_op.max_events), ("Well", 20, 1500))
        
    def testSerializeView(self):
        op = make_op()
        view = op.default_view()
        view.dendrogram = True
        new_view = self._roundtrip(view)
        self.assertIsInstance(new_view, CrossEntropyTestWorkflowView)
        self.assertTrue(new_view.dendrogram)


if __name__ == "__main__":
    unittest.main()
