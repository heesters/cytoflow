#!/usr/bin/env python3.8
# coding: latin-1
import unittest
import cytoflow as flow


class TestSampleLabels(unittest.TestCase):

    def setUp(self):
        import os
        cwd = os.path.dirname(os.path.abspath(__file__))
        self.ex = flow.ImportOp(
            conditions={'Dox': 'float'},
            tubes=[
                flow.Tube(file=os.path.join(cwd, 'data/Plate01/RFP_Well_A3.fcs'),
                          conditions={'Dox': 10.0}),
                flow.Tube(file=os.path.join(cwd, 'data/Plate01/RFP_Well_A6.fcs'),
                          conditions={'Dox': 1.0}),
            ]).apply()

    def testApply(self):
        op = flow.SampleLabelsOp(
            condition='Dox',
            labels={'10.0': 'High', '1.0': 'Low'})
        ex2 = op.apply(self.ex)
        self.assertIn('Dox', ex2.conditions)
        self.assertIn('High', ex2.data['Dox'].unique())
        self.assertIn('Low', ex2.data['Dox'].unique())

    def testNewCondition(self):
        op = flow.SampleLabelsOp(
            condition='Dox',
            labels={'10.0': 'High', '1.0': 'Low'},
            new_condition='Group')
        ex2 = op.apply(self.ex)
        self.assertIn('Dox', ex2.conditions)
        self.assertIn('Group', ex2.conditions)
        self.assertIn('High', ex2.data['Group'].unique())

    def testUnmappedKept(self):
        """Values not in labels dict keep their original value."""
        op = flow.SampleLabelsOp(
            condition='Dox',
            labels={'10.0': 'High'})
        ex2 = op.apply(self.ex)
        self.assertIn('1.0', ex2.data['Dox'].astype(str).unique())

    def testMissingCondition(self):
        import cytoflow.utility as util
        op = flow.SampleLabelsOp(condition='NonExistent', labels={'a': 'b'})
        with self.assertRaises(util.CytoflowOpError):
            op.apply(self.ex)

    def testEmptyLabels(self):
        import cytoflow.utility as util
        op = flow.SampleLabelsOp(condition='Dox', labels={})
        with self.assertRaises(util.CytoflowOpError):
            op.apply(self.ex)


if __name__ == '__main__':
    unittest.main()
