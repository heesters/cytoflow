'''
Tests for cytoflow.utility.report (PDF / PowerPoint reports)
'''

import os, re, shutil, tempfile, unittest, warnings

import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util

from .test_base import ImportedDataTest

try:
    import pptx
except ImportError:
    pptx = None


def pdf_pages(path):
    with open(path, 'rb') as f:
        return len(re.findall(rb"/Type\s*/Page[^s]", f.read()))


class TestReport(ImportedDataTest):
    
    def setUp(self):
        super().setUp()
        self.tmp = tempfile.mkdtemp()
        ex = flow.ThresholdOp(name = "Bright", channel = "Y2-A", threshold = 1000).apply(self.ex)
        self.ex = ex
        gate = flow.ThresholdOp(name = "G", channel = "Y2-A", threshold = 500)
        self.items = [flow.ReportItem(flow.HistogramView(channel = "Y2-A", scale = "logicle"), ex, "Histogram"),
                      flow.ReportItem(flow.ScatterplotView(xchannel = "V2-A", ychannel = "Y2-A",
                                                           xscale = "logicle", yscale = "logicle",
                                                           huefacet = "Bright"), ex, "Scatter",
                                      plot_kwargs = dict(alpha = 0.5)),
                      flow.ReportItem(flow.DensityView(xchannel = "V2-A", ychannel = "Y2-A",
                                                       xscale = "logicle", yscale = "logicle"), ex, "Density"),
                      flow.ReportItem(gate.default_view(scale = "logicle"), ex, "A gate"),
                      flow.ReportItem(flow.HistogramView(channel = "B1-A", scale = "logicle", 
                                                         xfacet = "Dox"), ex, "Facets")]
        
    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors = True)
        plt.close('all')
        
    def path(self, name):
        return os.path.join(self.tmp, name)
    
    def testVectorPdf(self):
        before = plt.get_fignums()
        problems = flow.export_report(self.items, self.path("vector.pdf"), title = "My report",
                                      layout = (1, 1))
        self.assertEqual(problems, [])
        self.assertEqual(pdf_pages(self.path("vector.pdf")), 5)
        with open(self.path("vector.pdf"), 'rb') as f:
            data = f.read()
        self.assertIn(b"/Font", data)             # real text, not outlines
        self.assertEqual(plt.get_fignums(), before)   # no figures left open
        
    def testGridPdf(self):
        problems = flow.export_report(self.items, self.path("grid.pdf"), title = "My report",
                                      layout = (2, 2), dpi = 100)
        self.assertEqual(problems, [])
        self.assertEqual(pdf_pages(self.path("grid.pdf")), 2)     # 4 + 1 plots
        self.assertGreater(os.path.getsize(self.path("grid.pdf")), 20000)
        
        flow.export_report(self.items, self.path("row.pdf"), layout = (1, 3), 
                           page_size = "Letter landscape", dpi = 60)
        self.assertEqual(pdf_pages(self.path("row.pdf")), 2)
        
    @unittest.skipIf(pptx is None, "python-pptx isn't installed")
    def testPowerPoint(self):
        problems = flow.export_report(self.items, self.path("report.pptx"), title = "My report",
                                      layout = (2, 2), dpi = 100)
        self.assertEqual(problems, [])
        prs = pptx.Presentation(self.path("report.pptx"))
        self.assertEqual(len(prs.slides), 2)
        self.assertAlmostEqual(prs.slide_width / prs.slide_height, 16 / 9, places = 2)
        
        first = prs.slides[0]
        pictures = [s for s in first.shapes if s.shape_type == pptx.enum.shapes.MSO_SHAPE_TYPE.PICTURE]
        texts = [s.text_frame.text for s in first.shapes if s.has_text_frame]
        self.assertEqual(len(pictures), 4)
        self.assertEqual(texts, ["My report", "Histogram", "Scatter", "Density", "A gate"])
        # the pictures are inside the slide
        for p in pictures:
            self.assertGreaterEqual(p.left, 0)
            self.assertLessEqual(p.left + p.width, prs.slide_width)
            self.assertLessEqual(p.top + p.height, prs.slide_height)
        self.assertEqual(len([s for s in prs.slides[1].shapes 
                              if s.shape_type == pptx.enum.shapes.MSO_SHAPE_TYPE.PICTURE]), 1)
        
    def testBadPlotIsReported(self):
        items = self.items[:1] + [flow.ReportItem(flow.HistogramView(channel = "nope"), self.ex, "Broken")]
        for name, layout in (("bad_vector.pdf", (1, 1)), ("bad_grid.pdf", (1, 2))):
            with warnings.catch_warnings(record = True) as w:
                warnings.simplefilter("always")
                problems = flow.export_report(items, self.path(name), layout = layout, dpi = 60)
            self.assertEqual(len(problems), 1)
            self.assertTrue(problems[0].startswith("Broken:"))
            self.assertTrue(any("Broken" in str(x.message) for x in w))
            self.assertEqual(pdf_pages(self.path(name)), 2 if layout == (1, 1) else 1)
            
    def testErrors(self):
        with self.assertRaises(util.CytoflowError):
            flow.export_report([], self.path("x.pdf"))
        with self.assertRaises(util.CytoflowError):
            flow.export_report(self.items, self.path("x.docx"))
        with self.assertRaises(util.CytoflowError):
            flow.export_report(self.items, self.path("x.pdf"), layout = (0, 2))
        with self.assertRaises(util.CytoflowError):
            flow.export_report(self.items, self.path("x.pdf"), page_size = "A0")


if __name__ == "__main__":
    unittest.main()
