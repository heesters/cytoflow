'''
Tests for PopulationTreeOp
'''

import os, unittest, warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util


class TestPopulationTree(unittest.TestCase):
    
    def setUp(self):
        # 2 samples x 1000 events.  A = first 600 (sample 1) / 400 (sample 2);
        # B = every 2nd event; C = every 4th; E = nothing
        ex = flow.Experiment()
        ex.add_channel("X")
        ex.add_condition("Sample", "int")
        for s in (1, 2):
            ex.add_events(pd.DataFrame({"X" : np.arange(1000, dtype = float)}), {"Sample" : s})
        i = np.tile(np.arange(1000), 2)
        sample = np.repeat([1, 2], 1000)
        ex.add_condition("A", "bool", pd.Series(i < np.where(sample == 1, 600, 400)))
        ex.add_condition("B", "bool", pd.Series(i % 2 == 0))
        ex.add_condition("C", "bool", pd.Series(i % 4 == 0))
        ex.add_condition("E", "bool", pd.Series(np.zeros(2000, dtype = bool)))
        ex.add_condition("Q", "category", pd.Series(np.where(i % 2 == 0, "Q_1", "Q_2")))
        self.ex = ex
        self.pops = [("Alive", "A", ""),
                     ("Even", "B", "Alive"),
                     ("Fourth", "C", "Even"),
                     ("Odd", "not B", "Alive"),
                     ("Nothing", "E", "Alive"),
                     ("Sub", "B", "Nothing")]
        
    def tearDown(self):
        plt.close('all')
        
    def testPooled(self):
        op = flow.PopulationTreeOp(name = "Tree", populations = self.pops)
        stat = op.apply(self.ex).statistics["Tree"]
        self.assertEqual(list(stat.index.get_level_values("Population")),
                         ["Alive", "Even", "Fourth", "Odd", "Nothing", "Sub"])
        self.assertEqual(list(stat["Count"]), [1000, 500, 250, 500, 0, 0])
        self.assertEqual(list(stat["Depth"]), [0, 1, 2, 1, 1, 2])
        np.testing.assert_allclose(stat["Percent_of_parent"].values[:5], [50, 50, 50, 50, 0])
        self.assertTrue(np.isnan(stat["Percent_of_parent"].values[5]))    # empty parent
        np.testing.assert_allclose(stat["Percent_of_total"].values, [50, 25, 12.5, 25, 0, 0])
        
    def testBy(self):
        op = flow.PopulationTreeOp(name = "Tree", populations = self.pops, by = ["Sample"])
        stat = op.apply(self.ex).statistics["Tree"]
        self.assertEqual(list(stat.index.names), ["Sample", "Population"])
        self.assertEqual(stat.loc[(1, "Alive"), "Count"], 600)
        self.assertEqual(stat.loc[(2, "Alive"), "Count"], 400)
        self.assertEqual(stat.loc[(1, "Fourth"), "Count"], 150)
        self.assertAlmostEqual(stat.loc[(2, "Even"), "Percent_of_parent"], 50.0)
        self.assertAlmostEqual(stat.loc[(2, "Even"), "Percent_of_total"], 20.0)
        
    def testOrderAndExpressions(self):
        # children listed before their parent; a categorical gate
        op = flow.PopulationTreeOp(name = "Tree", 
                                   populations = [("Q1", "Q == 'Q_1'", "Alive"),
                                                  ("Alive", "A", "")])
        stat = op.apply(self.ex).statistics["Tree"]
        self.assertEqual(list(stat.index.get_level_values("Population")), ["Alive", "Q1"])
        self.assertEqual(list(stat["Count"]), [1000, 500])
        
    def testDoesntChangeData(self):
        op = flow.PopulationTreeOp(name = "Tree", populations = self.pops)
        ex2 = op.apply(self.ex)
        self.assertEqual(list(ex2.data.columns), list(self.ex.data.columns))
        self.assertNotIn("Tree", self.ex.statistics)
        
    def testErrors(self):
        for pops in [[],
                     [("A1", "A", ""), ("A1", "B", "")],            # duplicate
                     [("A1", "A", "Nope")],                         # unknown parent
                     [("A1", "A", "A1")],                           # own parent
                     [("A1", "A", "B1"), ("B1", "B", "A1")],        # loop
                     [("A1", "NotACondition", "")],                 # bad gate
                     [("A1", "", "")],
                     [("", "A", "")]]:
            with self.assertRaises(util.CytoflowOpError):
                flow.PopulationTreeOp(name = "Tree", populations = pops).apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.PopulationTreeOp(populations = self.pops).apply(self.ex)
        with self.assertRaises(util.CytoflowOpError):
            flow.PopulationTreeOp(name = "Tree", populations = self.pops, by = ["Nope"]).apply(self.ex)
            
    def testView(self):
        op = flow.PopulationTreeOp(name = "Tree", populations = self.pops)
        op.default_view().plot(self.ex)
        texts = [t.get_text() for t in plt.gca().texts]
        self.assertIn("All events", texts)
        self.assertIn("Fourth", texts)
        self.assertTrue(any("250" in t and "50.0% of parent" in t for t in texts))
        op.default_view(subset = "Sample == 1").plot(self.ex, title = "Sample 1")
        self.assertTrue(any(t.get_text().startswith("600") for t in plt.gca().texts))
        with self.assertRaises(util.CytoflowViewError):
            op.default_view(subset = "Sample == 99").plot(self.ex)
            
    def testTableView(self):
        op = flow.PopulationTreeOp(name = "Tree", populations = self.pops, by = ["Sample"])
        ex2 = op.apply(self.ex)
        flow.TableView(statistic = "Tree", feature = "Count", 
                       row_facet = "Population", column_facet = "Sample").plot(ex2)


class TestInferPopulations(unittest.TestCase):
    
    def testFromBooleanGates(self):
        ops = [flow.PolygonOp(name = "Cells"),
               flow.PolygonOp(name = "Singlets"),
               flow.BooleanGateOp(name = "Singlets_gated", expression = "Singlets and Cells"),
               flow.RangeOp(name = "Pos"),
               flow.BooleanGateOp(name = "Pos_gated", expression = "Pos and Singlets_gated"),
               flow.BooleanGateOp(name = "Either", expression = "Pos or Cells"),
               flow.QuadOp(name = "Quad")]
        pops = flow.infer_populations(ops)
        self.assertEqual(pops, [("Cells", "Cells", ""),
                                ("Singlets", "Singlets", "Cells"),
                                ("Pos", "Pos", "Singlets"),
                                ("Either", "Either", "")])
        
    def testExplicitParents(self):
        ops = [flow.ThresholdOp(name = "A"), flow.ThresholdOp(name = "B"), flow.ThresholdOp(name = "C")]
        pops = flow.infer_populations(ops, parents = {"B" : "A", "C" : "Nope"})
        self.assertEqual(pops, [("A", "A", ""), ("B", "B", "A"), ("C", "C", "")])
        
    def testFlowJoWorkspace(self):
        data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ex = flow.ImportOp(tubes = [flow.Tube(file = os.path.join(data, "Plate01", "RFP_Well_A3.fcs"))]).apply()
            ops = flow.import_flowjo_workspace(os.path.join(data, "flowjo", "synthetic.wsp"), ex,
                                               sample = "RFP_Well_A3.fcs")
            for op in ops:
                ex = op.apply(ex)
        pops = flow.infer_populations(ops)
        parents = {n : p for n, _, p in pops}
        self.assertEqual(parents["Cells"], "")
        self.assertEqual(parents["Singlets"], "Cells")
        
        # the tree's counts are the "_gated" conditions' counts
        stat = flow.PopulationTreeOp(name = "Tree", populations = pops).apply(ex).statistics["Tree"]
        self.assertEqual(stat.loc["Singlets", "Count"], int(ex["Singlets_gated"].sum()))
        for name, _, parent in pops:
            if parent and name + "_gated" in ex.conditions and ex.data[name + "_gated"].dtype == bool:
                self.assertEqual(stat.loc[name, "Count"], int(ex[name + "_gated"].sum()), name)


if __name__ == "__main__":
    unittest.main()
