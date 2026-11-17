#!/usr/bin/env python3.8
# coding: latin-1

import os, unittest, tempfile

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import SampleLabelsWorkflowOp
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


class TestSampleLabels(ImportedDataTest):

    def setUp(self):
        super().setUp()

        self.addTypeEqualityFunc(SampleLabelsWorkflowOp, 'assertHasTraitsEqual')

        self.op = op = SampleLabelsWorkflowOp()
        op.condition = 'Dox'
        op.labels = {'10.0': 'High', '1.0': 'Low'}

        self.wi = wi = WorkflowItem(operation=op,
                                    status='waiting',
                                    view_error="Not yet plotted")
        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        self.workflow.wi_waitfor(wi, 'status', 'valid')

    def testApply(self):
        self.assertTrue(
            self.workflow.remote_eval("self.workflow[-1].result is not None"))

    def testChangeLabels(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.labels = {'10.0': 'HighDose', '1.0': 'LowDose'}
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertTrue(
            self.workflow.remote_eval("self.workflow[-1].result is not None"))

    def testNewCondition(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.new_condition = 'Group'
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertTrue(
            self.workflow.remote_eval(
                "'Group' in self.workflow[-1].result.conditions"))

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
        op = SampleLabelsWorkflowOp()
        try:
            os.close(fh)
            save_yaml(op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)
        self.maxDiff = None
        self.assertEqual(op, new_op)


if __name__ == '__main__':
    unittest.main()
