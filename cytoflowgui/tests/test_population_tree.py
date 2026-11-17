'''
GUI tests for the population tree
'''
import os, unittest, tempfile

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import (PopulationTreeWorkflowOp, PopulationTreeWorkflowView,
                                             PopulationTreePopulation, populations_from_workflow,
                                             ThresholdWorkflowOp, RangeWorkflowOp, QuadWorkflowOp,
                                             BooleanGateWorkflowOp)
from cytoflowgui.workflow.subset import BoolSubset, CategorySubset
from cytoflowgui.workflow.serialization import load_yaml, save_yaml


class TestPopulationTree(ImportedDataTest):

    def setUp(self):
        super().setUp()
        
        gate = ThresholdWorkflowOp(name = "Bright", channel = "Y2-A", threshold = 1000)
        wi = WorkflowItem(operation = gate, status = 'waiting')
        self.workflow.workflow.append(wi)
        self.workflow.wi_waitfor(wi, 'status', 'valid')

        self.op = op = PopulationTreeWorkflowOp()
        op.name = "Tree"
        op.populations_list = [PopulationTreePopulation(name = "Bright", gate = "Bright", parent = ""),
                               PopulationTreePopulation(name = "Dim", gate = "not Bright", parent = "")]
        op.by = ["Well"]

        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")
        self.view = view = wi.default_view
        wi.views.append(view)
        wi.current_view = view
        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        self.workflow.wi_waitfor(wi, 'status', 'valid')

    def testApplyPlotAndEdit(self):
        self.assertEqual(self.workflow.remote_eval("list(self.workflow[-1].result.statistics['Tree'].columns)"),
                         ['Count', 'Percent_of_parent', 'Percent_of_total', 'Depth'])
        self.assertEqual(self.workflow.remote_eval("list(self.workflow[-1].result.statistics['Tree'].index.names)"),
                         ['Well', 'Population'])
        self.workflow.wi_waitfor(self.wi, 'view_error', '')
        
        # edit a population: re-applies
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.populations_list[1].parent = "Bright"
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')
        self.assertEqual(self.workflow.remote_eval(
            "int(self.workflow[-1].result.statistics['Tree'].xs('Dim', level = 'Population')['Count'].sum())"), 0)
        
        # a bad parent is an error
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.populations_list[1].parent = "Nope"
        self.workflow.wi_waitfor(self.wi, 'status', 'invalid')

    def testNotebook(self):
        code = "import cytoflow as flow\nfrom cytoflow import geom_mean, geom_sd\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
            for view in wi.views:
                code = code + view.get_notebook_code(i)
        ns = {}
        exec(code, ns)
        self.assertIn("Tree", ns['ex_{}'.format(len(self.workflow.workflow) - 1)].statistics)


class TestPopulationsFromWorkflow(unittest.TestCase):
    
    def testInferFromSubsets(self):
        cells = ThresholdWorkflowOp(name = "Cells")
        singlets = RangeWorkflowOp(name = "Singlets")
        pos = ThresholdWorkflowOp(name = "Pos")
        neg = ThresholdWorkflowOp(name = "Neg")
        quad = QuadWorkflowOp(name = "Quad")
        other = ThresholdWorkflowOp(name = "Other")
        
        items = [(cells, []),
                 # drawn on a plot of Cells == True
                 (singlets, [BoolSubset(name = "Cells", selected_t = True)]),
                 (pos, [BoolSubset(name = "Singlets", selected_t = True),
                        BoolSubset(name = "Cells")]),               # unselected: ignored
                 # drawn on "Singlets == False": not a parent
                 (neg, [BoolSubset(name = "Singlets", selected_f = True)]),
                 (quad, []),
                 # two selections: ambiguous, so no parent
                 (other, [BoolSubset(name = "Cells", selected_t = True),
                          BoolSubset(name = "Singlets", selected_t = True)])]
        
        self.assertEqual(populations_from_workflow(items),
                         [("Cells", "Cells", ""),
                          ("Singlets", "Singlets", "Cells"),
                          ("Pos", "Pos", "Singlets"),
                          ("Neg", "Neg", ""),
                          ("Other", "Other", "")])
        
    def testImportedHierarchyWins(self):
        cells = ThresholdWorkflowOp(name = "Cells")
        singlets = RangeWorkflowOp(name = "Singlets")
        gated = BooleanGateWorkflowOp(name = "Singlets_gated", expression = "Singlets and Cells")
        self.assertEqual(populations_from_workflow([(cells, []), (singlets, []), (gated, [])]),
                         [("Cells", "Cells", ""), ("Singlets", "Singlets", "Cells")])


class TestPopulationTreeSerialize(unittest.TestCase):

    def _roundtrip(self, obj):
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
            save_yaml(obj, filename)
            return load_yaml(filename)
        finally:
            os.unlink(filename)

    def testSerialize(self):
        op = PopulationTreeWorkflowOp(name = "Tree", by = ["Well"])
        op.populations_list = [PopulationTreePopulation(name = "A", gate = "G1", parent = ""),
                               PopulationTreePopulation(name = "B", gate = "not G1", parent = "A")]
        new_op = self._roundtrip(op)
        self.assertIsInstance(new_op, PopulationTreeWorkflowOp)
        self.assertEqual(new_op.populations, [("A", "G1", ""), ("B", "not G1", "A")])
        self.assertEqual(new_op.by, ["Well"])
        self.assertIsInstance(self._roundtrip(op.default_view()), PopulationTreeWorkflowView)


if __name__ == "__main__":
    unittest.main()
