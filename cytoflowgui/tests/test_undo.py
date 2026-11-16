'''
Tests for cytoflowgui.workflow.undo
'''
import time, unittest

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import ThresholdWorkflowOp, RangeWorkflowOp
from cytoflowgui.workflow.undo import UndoManager, MERGE_SECONDS


class TestUndo(ImportedDataTest):

    def setUp(self):
        super().setUp()
        
        self.undo = UndoManager(model = self.workflow)

        self.op = op = ThresholdWorkflowOp()
        op.name = "Thresh"
        op.channel = "Y2-A"
        op.threshold = 1000

        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")
        wi.views.append(wi.default_view)
        wi.current_view = wi.default_view
        
        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        self.workflow.wi_waitfor(wi, 'status', "valid")
        
        # start each test with an empty history
        self.undo.clear()
        
    def _pause(self):
        time.sleep(MERGE_SECONDS * 1.5)

    def testUndoRedoParameter(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.threshold = 2000
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertTrue(self.undo.can_undo)
        
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.undo.undo()
        self.assertEqual(self.op.threshold, 1000)
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        
        # the remote process got the old value back, too
        self.assertEqual(self.workflow.remote_eval("self.workflow[-1].operation.threshold"), 1000)
        self.assertTrue(self.undo.can_redo)
        self.assertFalse(self.undo.can_undo)
        
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.undo.redo()
        self.assertEqual(self.op.threshold, 2000)
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertEqual(self.workflow.remote_eval("self.workflow[-1].operation.threshold"), 2000)

    def testSeparateSteps(self):
        self.op.threshold = 2000
        self._pause()
        self.op.threshold = 3000
        self.undo.undo()
        self.assertEqual(self.op.threshold, 2000)
        self.undo.undo()
        self.assertEqual(self.op.threshold, 1000)
        
    def testMergedStep(self):
        # quick successive changes (eg, a gate drag) are one step
        self.op.threshold = 2000
        self.op.channel = "V2-A"
        self.assertEqual(len(self.undo.undo_stack), 1)
        self.undo.undo()
        self.assertEqual(self.op.threshold, 1000)
        self.assertEqual(self.op.channel, "Y2-A")
        
    def testNewChangeClearsRedo(self):
        self.op.threshold = 2000
        self.undo.undo()
        self.assertTrue(self.undo.can_redo)
        self._pause()
        self.op.threshold = 500
        self.assertFalse(self.undo.can_redo)
        
    def testUndoDelete(self):
        self.workflow.workflow.remove(self.wi)
        self.assertNotIn(self.wi, self.workflow.workflow)
        self.undo.undo()
        self.assertIn(self.wi, self.workflow.workflow)
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertEqual(self.workflow.remote_eval("self.workflow[-1].operation.threshold"), 1000)
        
    def testUndoAdd(self):
        op = RangeWorkflowOp(name = "R", channel = "Y2-A", low = 100, high = 1000)
        wi = WorkflowItem(operation = op, status = 'waiting')
        n = len(self.workflow.workflow)
        self.workflow.workflow.append(wi)
        self.workflow.wi_waitfor(wi, 'status', 'valid')
        self._pause()
        self.undo.undo()
        self.assertNotIn(wi, self.workflow.workflow)
        self.assertEqual(self.workflow.remote_eval("len(self.workflow)"), n)
        self.undo.redo()
        self.assertIn(wi, self.workflow.workflow)

    def testViewChange(self):
        view = self.wi.current_view
        old = view.scale
        view.scale = "log"
        self.undo.undo()
        self.assertEqual(view.scale, old)


if __name__ == "__main__":
    unittest.main()
