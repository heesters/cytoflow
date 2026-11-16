'''
GUI tests for the "Axis Transforms" operation (ScaleSettingsOp)
'''
import os, unittest, tempfile

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import ScaleSettingsWorkflowOp
from cytoflowgui.workflow.operations.scale_settings import ChannelScaleSettings
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


class TestScaleSettings(ImportedDataTest):

    def setUp(self):
        super().setUp()

        self.addTypeEqualityFunc(ScaleSettingsWorkflowOp, 'assertHasTraitsEqual')
        self.addTypeEqualityFunc(ChannelScaleSettings, 'assertHasTraitsEqual')

        self.op = op = ScaleSettingsWorkflowOp()
        op.settings_list = [ChannelScaleSettings(channel = "Y2-A", 
                                                 logicle_W = 1.25,
                                                 arcsinh_cofactor = 5)]

        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")
        self.workflow.workflow.append(wi)
        self.workflow.wi_waitfor(wi, 'status', "valid")

    def testApply(self):
        md = self.workflow.remote_eval("self.workflow[-1].result.metadata['Y2-A']['scale_params']")
        self.assertEqual(md, {"logicle" : {"W" : 1.25}, "arcsinh" : {"cofactor" : 5.0}})

    def testChange(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.settings_list[0].logicle_W = 0.75
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        md = self.workflow.remote_eval("self.workflow[-1].result.metadata['Y2-A']['scale_params']")
        self.assertEqual(md["logicle"]["W"], 0.75)
        
    def testAddChannel(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.settings_list.append(ChannelScaleSettings(channel = "V2-A", logicle_M = 5))
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        md = self.workflow.remote_eval("self.workflow[-1].result.metadata['V2-A']['scale_params']")
        self.assertEqual(md, {"logicle" : {"M" : 5.0}})
        
    def testBadChannel(self):
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.settings_list[0].channel = "nope"
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')

    def testSerialize(self):
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(self.op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)
        self.assertEqual(len(new_op.settings_list), 1)
        self.assertEqual(new_op.settings_list[0].logicle_W, 1.25)
        self.assertEqual(new_op.settings_list[0].arcsinh_cofactor, 5.0)
        self.assertIsNone(new_op.settings_list[0].logicle_M)

    def testNotebook(self):
        code = "import cytoflow as flow\nfrom cytoflow import geom_mean, geom_sd\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
        ns = {}
        exec(code, ns)
        self.assertEqual(ns['ex_{}'.format(len(self.workflow.workflow) - 1)]
                         .metadata['Y2-A']['scale_params']['logicle']['W'], 1.25)


if __name__ == "__main__":
    unittest.main()
