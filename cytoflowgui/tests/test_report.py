'''
GUI tests for exporting a report (PDF / PowerPoint) of several plots
'''
import os, re, shutil, tempfile, unittest

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import ThresholdWorkflowOp
from cytoflowgui.workflow.views import HistogramWorkflowView, ScatterplotWorkflowView
from cytoflowgui.report_export import list_plots

try:
    import pptx
except ImportError:
    pptx = None


class TestReport(ImportedDataTest):

    def setUp(self):
        super().setUp()
        self.tmp = tempfile.mkdtemp()

        op = ThresholdWorkflowOp(name = "Thresh", channel = "Y2-A", threshold = 1000)
        self.wi = wi = WorkflowItem(operation = op, status = 'waiting',
                                    view_error = "Not yet plotted")
        self.gate_view = wi.default_view
        self.gate_view.scale = "logicle"
        wi.views.append(self.gate_view)
        
        self.hist = HistogramWorkflowView(channel = "V2-A", scale = "logicle", huefacet = "Thresh")
        self.hist.plot_params.title = "My histogram"
        wi.views.append(self.hist)
        self.scatter = ScatterplotWorkflowView(xchannel = "V2-A", ychannel = "Y2-A",
                                               xscale = "logicle", yscale = "logicle")
        wi.views.append(self.scatter)
        
        wi.current_view = self.gate_view
        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        self.workflow.wi_waitfor(wi, 'status', 'valid')
        self.workflow.wi_waitfor(wi, 'view_error', '')
        
    def tearDown(self):
        super().tearDown()
        shutil.rmtree(self.tmp, ignore_errors = True)
        
    def export(self, name, plots = None, **options):
        path = os.path.join(self.tmp, name)
        if plots is None:
            plots = [(self.wi, self.gate_view, "The gate"),
                     (self.wi, self.hist, "Histogram"),
                     (self.wi, self.scatter, "Scatter")]
        self.workflow.export_report(plots, path, **options)
        self.assertTrue(self.workflow.report_event.wait(120))
        return path, self.workflow.report_result

    def testPdf(self):
        path, (ok, message) = self.export("report.pdf", title = "Test", layout = (2, 2), dpi = 80)
        self.assertTrue(ok, message)
        self.assertIn("Saved 3 plots", message)
        with open(path, 'rb') as f:
            self.assertEqual(len(re.findall(rb"/Type\s*/Page[^s]", f.read())), 1)
            
        path, (ok, message) = self.export("vector.pdf", layout = (1, 1))
        self.assertTrue(ok, message)
        with open(path, 'rb') as f:
            self.assertEqual(len(re.findall(rb"/Type\s*/Page[^s]", f.read())), 3)
            
        # the plot on the screen still works afterwards
        self.workflow.wi_waitfor(self.wi, 'view_error', '')
        self.workflow.wi_sync(self.wi, 'view_error', 'waiting')
        self.gate_view.scale = "log"
        self.workflow.wi_waitfor(self.wi, 'view_error', '')
        self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.wi.operation.threshold = 2000
        self.workflow.wi_waitfor(self.wi, 'status', 'valid')

    @unittest.skipIf(pptx is None, "python-pptx isn't installed")
    def testPowerPoint(self):
        path, (ok, message) = self.export("report.pptx", title = "Test", layout = (1, 2), dpi = 80)
        self.assertTrue(ok, message)
        prs = pptx.Presentation(path)
        self.assertEqual(len(prs.slides), 2)
        texts = [s.text_frame.text for s in prs.slides[0].shapes if s.has_text_frame]
        self.assertEqual(texts, ["Test", "The gate", "Histogram"])
        
    def testProblems(self):
        # a view that can't be plotted is reported, the rest are exported
        self.scatter.xchannel = "nope"
        path, (ok, message) = self.export("partial.pdf", dpi = 60)
        self.assertTrue(ok, message)
        self.assertIn("Saved 2 plots", message)
        self.assertIn("Scatter", message)
        
        # a bad file type fails cleanly
        path, (ok, message) = self.export("report.docx")
        self.assertFalse(ok)
        self.assertIn(".pdf or .pptx", message)
        
    def testListPlots(self):
        plots = list_plots(self.workflow.workflow)
        labels = [label for _, _, label in plots]
        self.assertEqual(len(plots), 3)
        self.assertTrue(labels[0].endswith("Thresh: Threshold Selection"), labels[0])
        self.assertTrue(any("My histogram" in l for l in labels), labels)
        self.assertIs(plots[0][0], self.wi)


if __name__ == "__main__":
    unittest.main()
