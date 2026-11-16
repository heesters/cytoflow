#!/usr/bin/env python3.11
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2026
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Tests for Gating-ML export and import, and FlowJo workspace import.

If FlowKit (https://github.com/whitews/FlowKit) is installed, it is used as
an independent check; otherwise those tests are skipped.  FlowKit is NOT a
dependency of cytoflow.
'''

import unittest, os, re, tempfile, shutil, warnings
import xml.etree.ElementTree as ET

import numpy as np

import cytoflow as flow
import cytoflow.utility as util
import cytoflow.operations.gatingml as gml
from cytoflow.utility.logicle_ext.Logicle import Logicle

try:
    import flowkit as fk
except Exception:  # pragma: no cover
    fk = None

DATA = os.path.dirname(os.path.abspath(__file__)) + "/data/"
FCS = DATA + "Plate01/RFP_Well_A3.fcs"
WSP = DATA + "flowjo/synthetic.wsp"

# events right on the edge of a gate can go either way: cytoflow's logicle
# scale is a (very good) table-based approximation, and curved outlines are
# followed with polygons.  out of 10,000 events:
TOLERANCE = 3

VERTICES = [(23.4, 5158.7), (102.2, 23124.0), (510.9, 23124.0),
            (1089.5, 3800.3), (340.5, 801.9), (65.4, 1119.3)]

SCALES = [('linear', 'linear'), ('log', 'log'), ('logicle', 'logicle'),
          ('arcsinh', 'arcsinh'), ('logicle', 'linear'), ('log', 'arcsinh')]

# reasonable ellipse sizes in each scale's units
SIZES = {'linear' : (600, 12000), 'log' : (1.2, 0.9),
         'logicle' : (0.25, 0.2), 'arcsinh' : (2.5, 2.0)}


def make_ops():
    ops = [flow.ThresholdOp(name = "Thresh", channel = "Y2-A", threshold = 300),
           flow.RangeOp(name = "Range", channel = "B1-A", low = 100, high = 1000),
           flow.Range2DOp(name = "Range2D",
                          xchannel = "Violet", xlow = 10, xhigh = 500,
                          ychannel = "Y2-A", ylow = 100, yhigh = 10000),
           flow.QuadOp(name = "Quad",
                       xchannel = "Violet", xthreshold = 100,
                       ychannel = "Y2-A", ythreshold = 1000)]

    for xscale, yscale in SCALES:
        ops.append(flow.PolygonOp(name = "Poly_{}_{}".format(xscale, yscale),
                                  xchannel = "Violet", ychannel = "Y2-A",
                                  xscale = xscale, yscale = yscale,
                                  vertices = VERTICES))
        ops.append(flow.EllipseOp(name = "Ellipse_{}_{}".format(xscale, yscale),
                                  xchannel = "Violet", ychannel = "Y2-A",
                                  xscale = xscale, yscale = yscale,
                                  center = (200, 4000),
                                  width = SIZES[xscale][0],
                                  height = SIZES[yscale][1],
                                  # a rotation only makes sense if the
                                  # axes have the same units
                                  angle = 25 if xscale == yscale else 0))

    ops.append(flow.FreehandOp(name = "Freehand",
                               xchannel = "Violet", ychannel = "Y2-A",
                               xscale = "logicle", yscale = "log",
                               vertices = VERTICES))

    ops += [flow.BooleanGateOp(name = "Bool_simple",
                               expression = "Thresh and not Range"),
            flow.BooleanGateOp(name = "Bool_nested",
                               expression = "(Thresh or Range) and not (Range2D and Poly_log_log)"),
            flow.BooleanGateOp(name = "Bool_quad",
                               expression = "Quad == 'Quad_2' and Ellipse_logicle_logicle"),
            flow.BooleanGateOp(name = "Bool_not",
                               expression = "not (Thresh or Range)"),
            flow.BooleanGateOp(name = "Bool_alias",
                               expression = "Thresh"),
            flow.BooleanGateOp(name = "Bool_compare",
                               expression = "(Thresh == False) | (Quad != 'Quad_3')"),
            flow.BooleanGateOp(name = "Bool_bool",
                               expression = "Bool_simple or Bool_nested")]
    return ops


def apply_all(ops, ex):
    for op in ops:
        ex = op.apply(ex)
    return ex


def members(ex, name):
    """A gate's membership.  `name` can be a QuadOp category."""
    if name in ex.data:
        return ex[name].values.astype(bool)
    quad = name.rsplit("_", 1)[0]
    return (ex[quad] == name).values


class GateTest(unittest.TestCase):

    def setUp(self):
        # one tube, not sub-sampled, so the events are in the same order
        # as in the FCS file.  and rename a channel.
        self.ex = flow.ImportOp(tubes = [flow.Tube(file = FCS)],
                                channels = {"V2-A" : "Violet",
                                            "Y2-A" : "Y2-A",
                                            "B1-A" : "B1-A",
                                            "FSC-A" : "FSC-A",
                                            "SSC-A" : "SSC-A"}).apply()
        self.tmp = tempfile.mkdtemp()
        self.filename = os.path.join(self.tmp, "gates.xml")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors = True)

    def write(self, body, name = "in.xml"):
        path = os.path.join(self.tmp, name)
        with open(path, 'w') as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n'
                    '<gating:Gating-ML xmlns:gating="{}" xmlns:transforms="{}" '
                    'xmlns:data-type="{}">\n{}\n</gating:Gating-ML>\n'
                    .format(gml.GATING_NS, gml.TRANSFORMS_NS, gml.DATATYPE_NS, body))
        return path

    def do_import(self, func, *args, **kwargs):
        """Returns the ops and the text of the warnings"""
        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter('always')
            ops = func(*args, **kwargs)
        return ops, "\n".join([str(x.message) for x in w
                               if issubclass(x.category, util.CytoflowWarning)])

    def assertSameEvents(self, a, b, msg = None, tolerance = TOLERANCE):
        a = np.asarray(a, dtype = bool)
        b = np.asarray(b, dtype = bool)
        self.assertEqual(len(a), len(b))
        # and make sure the test means something
        self.assertGreater(a.sum(), 50, msg)
        self.assertLess(a.sum(), len(a) - 10, msg)
        self.assertLessEqual((a != b).sum(), tolerance, msg)


class TestExport(GateTest):

    def setUp(self):
        super().setUp()
        self.ops = make_ops()
        self.gated = apply_all(self.ops, self.ex)

    def testRoundTrip(self):
        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter('always')
            exported = flow.export_gatingml(self.ops, self.ex, self.filename)
        self.assertEqual([x for x in w if issubclass(x.category, util.CytoflowWarning)], [])
        self.assertEqual(exported, [op.name for op in self.ops])

        ops, warns = self.do_import(flow.import_gatingml, self.filename, self.ex)
        self.assertEqual(warns, "")
        by_name = {op.name : op for op in ops}

        # same kinds of gates ...
        for op in self.ops:
            expected = flow.PolygonOp if isinstance(op, flow.FreehandOp) else type(op)
            self.assertIs(type(by_name[op.name]), expected, op.name)
            for trait in ('channel', 'xchannel', 'ychannel', 'xscale', 'yscale'):
                if trait in op.trait_names():
                    self.assertEqual(getattr(op, trait), getattr(by_name[op.name], trait))

        # ... with the same parameters ...
        self.assertAlmostEqual(by_name['Thresh'].threshold, 300)
        self.assertAlmostEqual(by_name['Range'].low, 100)
        self.assertAlmostEqual(by_name['Range'].high, 1000)
        self.assertAlmostEqual(by_name['Quad'].xthreshold, 100)
        self.assertAlmostEqual(by_name['Quad'].ythreshold, 1000)
        for xscale, yscale in SCALES:
            p = by_name["Poly_{}_{}".format(xscale, yscale)]
            np.testing.assert_allclose(p.vertices, VERTICES, rtol = 1e-3)
            e = by_name["Ellipse_{}_{}".format(xscale, yscale)]
            np.testing.assert_allclose(e.center, (200, 4000), rtol = 1e-3)
            self.assertAlmostEqual(e.width, SIZES[xscale][0], delta = 1e-4 * e.width)
            self.assertAlmostEqual(e.height, SIZES[yscale][1], delta = 1e-4 * e.height)
            self.assertAlmostEqual(e.angle, 25 if xscale == yscale else 0, places = 4)

        # ... and the same events
        regated = apply_all(ops, self.ex)
        for op in self.ops:
            if isinstance(op, flow.QuadOp):
                self.assertTrue((self.gated[op.name] == regated[op.name]).all())
            else:
                self.assertSameEvents(self.gated[op.name], regated[op.name], op.name)

    @unittest.skipIf(fk is None, "FlowKit isn't installed")
    def testFlowKit(self):
        """Does FlowKit put the same events in the exported gates?"""
        flow.export_gatingml(self.ops, self.ex, self.filename)

        strategy = fk.parse_gating_xml(self.filename)
        sample = fk.Sample(FCS)
        results = strategy.gate_sample(sample)

        # same events, in the same order?
        np.testing.assert_array_equal(
            sample.get_events('raw')[:, sample.pnn_labels.index('V2-A')],
            self.ex['Violet'].values)

        n = 0
        for op in self.ops:
            names = ([op.name + "_{}".format(i) for i in range(1, 5)]
                     if isinstance(op, flow.QuadOp) else [op.name])
            for name in names:
                self.assertSameEvents(members(self.gated, name),
                                      results.get_gate_membership(name),
                                      name)
                n += 1
        self.assertEqual(n, len(self.ops) + 3)

    def testFcsNames(self):
        flow.export_gatingml(self.ops, self.ex, self.filename)
        root = ET.parse(self.filename).getroot()
        self.assertEqual(root.tag, "{{{}}}Gating-ML".format(gml.GATING_NS))
        names = set([el.get("{{{}}}name".format(gml.DATATYPE_NS))
                     for el in root.iter("{{{}}}fcs-dimension".format(gml.DATATYPE_NS))])
        self.assertEqual(names, {"V2-A", "Y2-A", "B1-A"})

        # all four kinds of transform; none on the gates without a scale
        kinds = set([el[0].tag.split('}')[1]
                     for el in root.iter("{{{}}}transformation".format(gml.TRANSFORMS_NS))])
        self.assertEqual(kinds, {'flin', 'flog', 'logicle', 'fasinh'})
        for gate in root.iter("{{{}}}RectangleGate".format(gml.GATING_NS)):
            for dim in gate:
                self.assertIsNone(dim.get("{{{}}}transformation-ref".format(gml.GATING_NS)))

    def testTransformParameters(self):
        """The exported parameters are the ones cytoflow uses"""
        ex = flow.ScaleSettingsOp(settings = {"Violet" : {"logicle" : {"W" : 1.1, "A" : 0.3}},
                                              "Y2-A" : {"arcsinh" : {"cofactor" : 500}}}).apply(self.ex)
        poly = flow.PolygonOp(name = "Poly",
                              xchannel = "Violet", ychannel = "Y2-A",
                              xscale = "logicle", yscale = "arcsinh",
                              vertices = VERTICES)
        flow.export_gatingml([poly], ex, self.filename)
        root = ET.parse(self.filename).getroot()
        t = lambda n: "{{{}}}{}".format(gml.TRANSFORMS_NS, n)

        logicle = next(root.iter(t('logicle')))
        self.assertAlmostEqual(float(logicle.get(t('T'))), 262144.0)
        self.assertAlmostEqual(float(logicle.get(t('W'))), 1.1)
        self.assertAlmostEqual(float(logicle.get(t('M'))), 4.5)
        self.assertAlmostEqual(float(logicle.get(t('A'))), 0.3)

        # asinh(x sinh(M ln 10) / T) == asinh(x / cofactor)
        fasinh = next(root.iter(t('fasinh')))
        T, M, A = [float(fasinh.get(t(n))) for n in "TMA"]
        self.assertAlmostEqual(T / np.sinh(M * np.log(10)), 500.0, places = 6)
        self.assertEqual(A, 0.0)

        coords = [float(c.get("{{{}}}value".format(gml.DATATYPE_NS)))
                  for c in root.iter("{{{}}}coordinate".format(gml.GATING_NS))]
        # the vertices are in the transformed space.  (cytoflow's own
        # logicle scale is binned, so it's only close.)
        scale = util.scale_factory("logicle", ex, channel = "Violet")
        self.assertAlmostEqual(coords[0], scale(VERTICES[0][0]), delta = 1.0 / 4096)
        self.assertAlmostEqual(coords[0], Logicle(262144.0, 1.1, 4.5, 0.3).scale(VERTICES[0][0]))
        self.assertAlmostEqual(coords[1], np.arcsinh(VERTICES[0][1] / 500.0) / (M * np.log(10)))

        # and it comes back the same, with the same settings ...
        ops, warns = self.do_import(flow.import_gatingml, self.filename, ex)
        self.assertEqual(warns, "")
        self.assertEqual(len(ops[0].vertices), len(VERTICES))
        self.assertSameEvents(poly.apply(ex)["Poly"], ops[0].apply(ex)["Poly"])

        # ... and with different settings, it gets extra vertices
        ops, warns = self.do_import(flow.import_gatingml, self.filename, self.ex)
        self.assertIn("follows the original outline", warns)
        self.assertEqual(len(ops[0].vertices), len(VERTICES) * gml._EDGE_SEGMENTS)
        self.assertSameEvents(poly.apply(ex)["Poly"], ops[0].apply(self.ex)["Poly"])

    def testExperimentPerOp(self):
        ex2 = flow.ScaleSettingsOp(settings = {"Violet" : {"logicle" : {"W" : 1.1}}}).apply(self.ex)
        polys = [flow.PolygonOp(name = n,
                                xchannel = "Violet", ychannel = "Y2-A",
                                xscale = "logicle", yscale = "logicle",
                                vertices = VERTICES) for n in ("P1", "P2")]
        flow.export_gatingml(polys, [self.ex, ex2], self.filename)
        root = ET.parse(self.filename).getroot()
        t = lambda n: "{{{}}}{}".format(gml.TRANSFORMS_NS, n)
        ws = sorted([float(el.get(t('W'))) for el in root.iter(t('logicle'))])
        self.assertEqual(len(ws), 3)   # Violet twice, Y2-A once
        self.assertAlmostEqual(ws[-1], 1.1)

        with self.assertRaises(util.CytoflowError):
            flow.export_gatingml(polys, [self.ex], self.filename)

    def testSkipped(self):
        ops = [flow.ThresholdOp(name = "Thresh", channel = "Y2-A", threshold = 300),
               flow.RatioOp(name = "Ratio", numerator = "Y2-A", denominator = "B1-A"),
               flow.ThresholdOp(name = "Unset", channel = "Y2-A"),
               flow.ThresholdOp(name = "BadChannel", channel = "Nope", threshold = 1.0),
               flow.BooleanGateOp(name = "Compare", expression = "Thresh and Dox > 1"),
               flow.BooleanGateOp(name = "Unknown", expression = "Thresh and Other"),
               flow.PolygonOp(name = "NegativeLog",
                              xchannel = "Violet", ychannel = "Y2-A",
                              xscale = "log", yscale = "log",
                              vertices = [(-10, 10), (100, 10), (100, 100)]),
               flow.ThresholdOp(name = "Thresh", channel = "Y2-A", threshold = 500),
               flow.PolygonOp(name = "Magnetic",
                              xchannel = "Violet", ychannel = "Y2-A",
                              vertices = VERTICES, adjust_by = ["Dox"])]
        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter('always')
            exported = flow.export_gatingml(ops, self.ex, self.filename)
        self.assertEqual(exported, ["Thresh", "Magnetic"])
        text = "\n".join([str(x.message) for x in w])
        for name in ("Ratio", "Unset", "BadChannel", "Compare", "Unknown",
                     "NegativeLog", "already a gate"):
            self.assertIn(name, text)
        self.assertIn("un-adjusted", text)

        ops, _ = self.do_import(flow.import_gatingml, self.filename, self.ex)
        self.assertEqual([op.name for op in ops], ["Thresh", "Magnetic"])


class TestImportGatingML(GateTest):

    TRANSFORMS = '''
      <transforms:transformation transforms:id="lgcl">
        <transforms:logicle transforms:T="100000" transforms:W="0.8" transforms:M="4.2" transforms:A="0.2"/>
      </transforms:transformation>
      <transforms:transformation transforms:id="asinh">
        <transforms:fasinh transforms:T="100000" transforms:M="4.0" transforms:A="1.0"/>
      </transforms:transformation>
      <transforms:transformation transforms:id="hyper">
        <transforms:hyperlog transforms:T="100000" transforms:W="0.5" transforms:M="4.5" transforms:A="0"/>
      </transforms:transformation>
      <transforms:transformation transforms:id="lin">
        <transforms:flin transforms:T="10000" transforms:A="1000"/>
      </transforms:transformation>
      <transforms:transformation transforms:id="log">
        <transforms:flog transforms:T="10000" transforms:M="5"/>
      </transforms:transformation>
      <transforms:transformation transforms:id="ratio">
        <transforms:fratio transforms:A="1" transforms:B="0" transforms:C="0">
          <data-type:fcs-dimension data-type:name="V2-A"/>
          <data-type:fcs-dimension data-type:name="Y2-A"/>
        </transforms:fratio>
      </transforms:transformation>
    '''

    @staticmethod
    def dim(name, t = None, lo = None, hi = None, comp = "uncompensated"):
        return ('<gating:dimension gating:compensation-ref="{}"{}{}{}>'
                '<data-type:fcs-dimension data-type:name="{}"/></gating:dimension>'
                .format(comp,
                        ' gating:transformation-ref="{}"'.format(t) if t else '',
                        ' gating:min="{}"'.format(lo) if lo is not None else '',
                        ' gating:max="{}"'.format(hi) if hi is not None else '',
                        name))

    @classmethod
    def polygon(cls, gate_id, tx, ty, vertices, parent = None):
        return ('<gating:PolygonGate gating:id="{}"{}>{}{}{}</gating:PolygonGate>'
                .format(gate_id,
                        ' gating:parent_id="{}"'.format(parent) if parent else '',
                        cls.dim("V2-A", tx), cls.dim("Y2-A", ty),
                        "".join(['<gating:vertex><gating:coordinate data-type:value="{}"/>'
                                 '<gating:coordinate data-type:value="{}"/></gating:vertex>'
                                 .format(x, y) for x, y in vertices])))

    @classmethod
    def ellipsoid(cls, gate_id, tx, ty, mean, cov, d2):
        return ('<gating:EllipsoidGate gating:id="{}">{}{}'
                '<gating:mean><gating:coordinate data-type:value="{}"/>'
                '<gating:coordinate data-type:value="{}"/></gating:mean>'
                '<gating:covarianceMatrix>'
                '<gating:row><gating:entry data-type:value="{}"/><gating:entry data-type:value="{}"/></gating:row>'
                '<gating:row><gating:entry data-type:value="{}"/><gating:entry data-type:value="{}"/></gating:row>'
                '</gating:covarianceMatrix>'
                '<gating:distanceSquare data-type:value="{}"/></gating:EllipsoidGate>'
                .format(gate_id, cls.dim("V2-A", tx), cls.dim("Y2-A", ty),
                        mean[0], mean[1], cov[0][0], cov[0][1], cov[1][0], cov[1][1], d2))

    @staticmethod
    def rectangle(gate_id, dims, parent = None):
        return ('<gating:RectangleGate gating:id="{}"{}>{}</gating:RectangleGate>'
                .format(gate_id,
                        ' gating:parent_id="{}"'.format(parent) if parent else '',
                        "".join(dims)))

    def testHierarchy(self):
        path = self.write(
            # a child before its parent
            self.rectangle("Child", [self.dim("B1-A", lo = 30, hi = 2000)], parent = "Middle")
            + self.rectangle("Middle", [self.dim("Y2-A", lo = 100)], parent = "Top")
            + self.rectangle("Top", [self.dim("V2-A", lo = 0, hi = 800),
                                     self.dim("Y2-A", hi = 30000)])
            + '''<gating:BooleanGate gating:id="NotChild" gating:parent_id="Top">
                   <gating:not><gating:gateReference gating:ref="Child"/></gating:not>
                 </gating:BooleanGate>
                 <gating:BooleanGate gating:id="Or">
                   <gating:or>
                     <gating:gateReference gating:ref="Child"/>
                     <gating:gateReference gating:ref="Middle" gating:use-as-complement="true"/>
                   </gating:or>
                 </gating:BooleanGate>''')

        ops, warns = self.do_import(flow.import_gatingml, path, self.ex)
        self.assertEqual(warns, "")
        self.assertEqual([(type(op).__name__, op.name) for op in ops],
                         [('Range2DOp', 'Top'),
                          ('ThresholdOp', 'Middle'),
                          ('BooleanGateOp', 'Middle_gated'),
                          ('RangeOp', 'Child'),
                          ('BooleanGateOp', 'Child_gated'),
                          ('BooleanGateOp', 'NotChild'),
                          ('BooleanGateOp', 'NotChild_gated'),
                          ('BooleanGateOp', 'Or')])
        by_name = {op.name : op for op in ops}
        self.assertEqual(by_name['Middle_gated'].expression, "Middle and Top")
        self.assertEqual(by_name['Child_gated'].expression, "Child and Middle_gated")
        self.assertEqual(by_name['NotChild'].expression, "not Child_gated")
        self.assertEqual(by_name['NotChild_gated'].expression, "NotChild and Top")
        self.assertEqual(by_name['Or'].expression, "Child_gated or not Middle_gated")
        self.assertEqual(by_name['Top'].xlow, 0)
        self.assertEqual(by_name['Top'].ylow, -np.inf)
        self.assertEqual(by_name['Top'].yhigh, 30000)

        ex = apply_all(ops, self.ex)
        v, y, b = [self.ex[c].values for c in ("Violet", "Y2-A", "B1-A")]
        top = (v >= 0) & (v <= 800) & (y <= 30000)
        middle = top & (y > 100)
        child = middle & (b >= 30) & (b <= 2000)
        self.assertSameEvents(ex['Top'], top, tolerance = 0)
        self.assertSameEvents(ex['Middle_gated'], middle, tolerance = 0)
        self.assertSameEvents(ex['Child_gated'], child, tolerance = 0)
        self.assertSameEvents(ex['NotChild_gated'], top & ~child, tolerance = 0)
        self.assertSameEvents(ex['Or'], child | ~middle, tolerance = 0)

        if fk is not None:
            results = fk.parse_gating_xml(path).gate_sample(fk.Sample(FCS))
            for theirs, ours in (("Top", "Top"), ("Middle", "Middle_gated"),
                                 ("Child", "Child_gated"), ("NotChild", "NotChild_gated"),
                                 ("Or", "Or")):
                self.assertSameEvents(results.get_gate_membership(theirs), ex[ours],
                                      theirs, tolerance = 1)

    def foreign(self):
        """Gates on transforms whose parameters aren't cytoflow's"""
        lgcl = Logicle(100000.0, 0.8, 4.2, 0.2)
        asinh = lambda x: (np.arcsinh(x * np.sinh(4.0 * np.log(10)) / 100000) + np.log(10)) / (5.0 * np.log(10))
        lin = lambda x: (x + 1000) / 11000.0
        log = lambda x: np.log10(x / 10000.0) / 5.0 + 1

        body = self.TRANSFORMS
        body += self.polygon("Poly_logicle", "lgcl", "lgcl",
                             [(lgcl.scale(x), lgcl.scale(y)) for x, y in VERTICES])
        body += self.polygon("Poly_asinh_lin", "asinh", "lin",
                             [(asinh(x), lin(y)) for x, y in VERTICES])
        body += self.polygon("Poly_hyper_log", "hyper", "log",
                             [(0.3, log(900)), (0.3, log(20000)), (0.6, log(20000)), (0.55, log(900))])
        body += self.ellipsoid("Ell_logicle", "lgcl", "lgcl",
                               (lgcl.scale(200.0), lgcl.scale(4000.0)),
                               [[0.004, 0.001], [0.001, 0.003]], 2.5)
        body += self.ellipsoid("Ell_asinh_log", "asinh", "log",
                               (asinh(200.0), log(4000.0)),
                               [[0.004, 0.001], [0.001, 0.003]], 2.5)
        body += self.ellipsoid("Ell_lin_log", "lin", "log",
                               (lin(200.0), log(4000.0)),
                               [[0.0004, 0.0003], [0.0003, 0.003]], 4.0)
        body += self.ellipsoid("Ell_none", None, None,
                               (200.0, 4000.0), [[40000, 300000], [300000, 9000000]], 1.5)
        return self.write(body), lgcl, asinh, lin, log

    def testForeignTransforms(self):
        path, lgcl, asinh, lin, log = self.foreign()
        ops, warns = self.do_import(flow.import_gatingml, path, self.ex)
        by_name = {op.name : op for op in ops}
        self.assertNotIn("not imported", warns)

        # only a shift and a stretch: stay what they are
        self.assertIsInstance(by_name['Ell_lin_log'], flow.EllipseOp)
        self.assertEqual((by_name['Ell_lin_log'].xscale, by_name['Ell_lin_log'].yscale),
                         ('linear', 'log'))
        self.assertIsInstance(by_name['Ell_none'], flow.EllipseOp)
        for name in ('Ell_lin_log', 'Ell_none'):
            self.assertNotIn(name, warns)

        # otherwise, they follow the curved outlines
        for name in ('Ell_logicle', 'Ell_asinh_log'):
            self.assertIsInstance(by_name[name], flow.PolygonOp)
            self.assertEqual(len(by_name[name].vertices), gml._ELLIPSE_VERTICES)
            self.assertIn("'{}'".format(name), warns)
        self.assertEqual(len(by_name['Poly_logicle'].vertices), len(VERTICES) * gml._EDGE_SEGMENTS)
        self.assertEqual((by_name['Poly_asinh_lin'].xscale, by_name['Poly_asinh_lin'].yscale),
                         ('arcsinh', 'linear'))
        self.assertEqual((by_name['Poly_hyper_log'].xscale, by_name['Poly_hyper_log'].yscale),
                         ('logicle', 'log'))

        # check the ellipses against the definition
        ex = apply_all(ops, self.ex)
        v, y = self.ex["Violet"].values, self.ex["Y2-A"].values
        with np.errstate(invalid = 'ignore', divide = 'ignore'):
            spaces = {'lgcl' : (np.vectorize(lgcl.scale)(v), np.vectorize(lgcl.scale)(y)),
                      'asinh_log' : (asinh(v), log(y)),
                      'lin_log' : (lin(v), log(y)),
                      'none' : (v, y)}

        def inside(gx, gy, mean, cov, d2):
            d = np.column_stack((gx - mean[0], gy - mean[1]))
            with np.errstate(invalid = 'ignore'):
                return np.nan_to_num(np.einsum('ij,jk,ik->i', d, np.linalg.inv(cov), d),
                                     nan = np.inf) <= d2

        c1 = [[0.004, 0.001], [0.001, 0.003]]
        self.assertSameEvents(ex['Ell_logicle'],
                              inside(*spaces['lgcl'], (lgcl.scale(200.0), lgcl.scale(4000.0)), c1, 2.5))
        self.assertSameEvents(ex['Ell_asinh_log'],
                              inside(*spaces['asinh_log'], (asinh(200.0), log(4000.0)), c1, 2.5))
        self.assertSameEvents(ex['Ell_lin_log'],
                              inside(*spaces['lin_log'], (lin(200.0), log(4000.0)),
                                     [[0.0004, 0.0003], [0.0003, 0.003]], 4.0),
                              tolerance = 0)
        self.assertSameEvents(ex['Ell_none'],
                              inside(*spaces['none'], (200.0, 4000.0),
                                     [[40000, 300000], [300000, 9000000]], 1.5),
                              tolerance = 0)

    @unittest.skipIf(fk is None, "FlowKit isn't installed")
    def testForeignTransformsFlowKit(self):
        path = self.foreign()[0]
        ops, _ = self.do_import(flow.import_gatingml, path, self.ex)
        ex = apply_all(ops, self.ex)
        results = fk.parse_gating_xml(path).gate_sample(fk.Sample(FCS))
        self.assertEqual(len(ops), 7)
        for op in ops:
            self.assertSameEvents(ex[op.name], results.get_gate_membership(op.name), op.name)

    def testNotImported(self):
        body = self.TRANSFORMS
        body += self.rectangle("Good", [self.dim("Y2-A", "lgcl", lo = 0.5)])
        body += self.rectangle("ThreeD", [self.dim("Y2-A", lo = 1), self.dim("V2-A", lo = 1),
                                          self.dim("B1-A", lo = 1)])
        body += self.rectangle("NoChannel", [self.dim("PE-A", lo = 1)])
        body += self.rectangle("Child", [self.dim("Y2-A", lo = 1)], parent = "NoChannel")
        body += self.rectangle("Orphan", [self.dim("Y2-A", lo = 1)], parent = "Nobody")
        body += self.rectangle("NoTransform", [self.dim("Y2-A", "nope", lo = 1)])
        body += ('<gating:RectangleGate gating:id="Ratio">'
                 '<gating:dimension gating:min="1">'
                 '<data-type:new-dimension data-type:transformation-ref="ratio"/>'
                 '</gating:dimension></gating:RectangleGate>')
        body += self.polygon("OneD", None, None, [(1, 2), (3, 4), (5, 6)]).replace(self.dim("Y2-A"), "")
        body += ('<gating:BooleanGate gating:id="Bool"><gating:and>'
                 '<gating:gateReference gating:ref="Good"/>'
                 '<gating:gateReference gating:ref="ThreeD"/>'
                 '</gating:and></gating:BooleanGate>')
        body += ('<gating:BooleanGate gating:id="Loop"><gating:not>'
                 '<gating:gateReference gating:ref="Loop"/>'
                 '</gating:not></gating:BooleanGate>')
        body += self.rectangle("Good", [self.dim("Y2-A", lo = 5)])
        body += self.rectangle("Comp", [self.dim("Y2-A", lo = 5, comp = "FCS")])

        ops, warns = self.do_import(flow.import_gatingml, self.write(body), self.ex)
        self.assertEqual([op.name for op in ops], ["Good", "Comp"])
        for name, why in (("ThreeD", "3 dimensions"), ("NoChannel", "PE-A"),
                          ("Child", "its parent 'NoChannel'"), ("Orphan", "Nobody"),
                          ("NoTransform", "nope"), ("Ratio", "ratio"), ("OneD", "1 dimension"),
                          ("Bool", "'ThreeD' was not imported"), ("Loop", "depends on itself"),
                          ("Good", "more than one")):
            self.assertRegex(warns, r" - {}: .*{}".format(name, re.escape(why)))
        self.assertIn("Compensation is NOT imported", warns)

        # the bound is converted to data units
        self.assertAlmostEqual(ops[0].threshold, Logicle(100000.0, 0.8, 4.2, 0.2).inverse(0.5))

    def testQuadrants(self):
        divider = ('<gating:divider gating:id="{}" gating:compensation-ref="uncompensated">'
                   '<data-type:fcs-dimension data-type:name="{}"/>{}</gating:divider>')
        quadrant = ('<gating:Quadrant gating:id="{}">'
                    '<gating:position gating:divider_ref="dx" gating:location="{}"/>'
                    '<gating:position gating:divider_ref="dy" gating:location="{}"/>'
                    '</gating:Quadrant>')
        body = ('<gating:QuadrantGate gating:id="Q" gating:parent_id="Top">'
                + divider.format("dx", "V2-A", "<gating:value>100</gating:value>")
                + divider.format("dy", "Y2-A", "<gating:value>1000</gating:value><gating:value>5000</gating:value>")
                + quadrant.format("low_low", 50, 500)
                + quadrant.format("high_mid", 150, 2000)
                + quadrant.format("high_high", 150, 10000)
                + '</gating:QuadrantGate>'
                + self.rectangle("Top", [self.dim("B1-A", lo = 30)]))
        ops, warns = self.do_import(flow.import_gatingml, self.write(body), self.ex)
        self.assertEqual(warns, "")
        self.assertEqual([(type(op).__name__, op.name) for op in ops],
                         [('ThresholdOp', 'Top'),
                          ('Range2DOp', 'low_low'), ('BooleanGateOp', 'low_low_gated'),
                          ('Range2DOp', 'high_mid'), ('BooleanGateOp', 'high_mid_gated'),
                          ('Range2DOp', 'high_high'), ('BooleanGateOp', 'high_high_gated')])
        mid = ops[3]
        self.assertEqual((mid.xlow, mid.xhigh, mid.ylow, mid.yhigh), (100, np.inf, 1000, 5000))

        ex = apply_all(ops, self.ex)
        v, y, b = [self.ex[c].values for c in ("Violet", "Y2-A", "B1-A")]
        self.assertSameEvents(ex['high_mid_gated'],
                              (v >= 100) & (y >= 1000) & (y <= 5000) & (b > 30),
                              tolerance = 0)

    def testQuadOpWithParent(self):
        """An exported QuadOp, to which someone added a parent"""
        quad = flow.QuadOp(name = "Quad",
                           xchannel = "Violet", xthreshold = 100,
                           ychannel = "Y2-A", ythreshold = 1000)
        thresh = flow.ThresholdOp(name = "Top", channel = "B1-A", threshold = 30)
        flow.export_gatingml([quad, thresh], self.ex, self.filename)
        with open(self.filename) as f:
            text = f.read()
        text = text.replace('<gating:QuadrantGate gating:id="Quad">',
                            '<gating:QuadrantGate gating:id="Quad" gating:parent_id="Top">')
        text = text.replace('</gating:Gating-ML>',
                            '<gating:BooleanGate gating:id="Bool"><gating:or>'
                            '<gating:gateReference gating:ref="Quad_1"/>'
                            '<gating:gateReference gating:ref="Quad_4"/>'
                            '</gating:or></gating:BooleanGate></gating:Gating-ML>')
        with open(self.filename, 'w') as f:
            f.write(text)

        ops, warns = self.do_import(flow.import_gatingml, self.filename, self.ex)
        self.assertEqual(warns, "")
        by_name = {op.name : op for op in ops}
        self.assertIsInstance(by_name["Quad"], flow.QuadOp)
        self.assertEqual(by_name["Quad_1_gated"].expression, "(Quad == 'Quad_1') and Top")
        self.assertEqual(by_name["Bool"].expression, "Quad_1_gated or Quad_4_gated")
        ex = apply_all(ops, self.ex)
        v, y, b = [self.ex[c].values for c in ("Violet", "Y2-A", "B1-A")]
        self.assertSameEvents(ex['Bool'],
                              (b > 30) & (((v < 100) & (y > 1000)) | ((v > 100) & (y < 1000))),
                              tolerance = 0)

        # if the name is taken, the quadrants become rectangle gates instead
        ex2 = flow.ThresholdOp(name = "Quad", channel = "Y2-A", threshold = 1).apply(self.ex)
        ops, warns = self.do_import(flow.import_gatingml, self.filename, ex2)
        self.assertEqual(warns, "")
        self.assertEqual([type(op) for op in ops if op.name.startswith("Quad")
                          and not op.name.endswith("gated")], [flow.Range2DOp] * 4)
        ex3 = apply_all(ops, ex2)
        self.assertSameEvents(ex['Bool'], ex3['Bool'], tolerance = 0)

    def testNames(self):
        body = (self.rectangle("CD4+", [self.dim("Y2-A", lo = 1)])
                + self.rectangle("CD4-", [self.dim("Y2-A", hi = 1)])
                + self.rectangle("Y2-A", [self.dim("Y2-A", lo = 2)])
                + self.rectangle("2 much", [self.dim("Y2-A", lo = 3)])
                + self.rectangle("not", [self.dim("Y2-A", lo = 4)])
                + self.rectangle("Dox", [self.dim("Y2-A", lo = 5)], parent = "CD4+"))
        ex = flow.ThresholdOp(name = "Dox", channel = "Y2-A", threshold = 1).apply(self.ex)
        ops, _ = self.do_import(flow.import_gatingml, self.write(body), ex)
        self.assertEqual([op.name for op in ops],
                         ["CD4pos", "CD4neg", "Y2_A_2", "Gate_2_much", "not_gate",
                          "Dox_2", "Dox_2_gated"])
        self.assertEqual(ops[-1].expression, "Dox_2 and CD4pos")
        apply_all(ops, ex)

    def testNotGatingML(self):
        with self.assertRaisesRegex(util.CytoflowError, "FlowJo workspace"):
            flow.import_gatingml(WSP, self.ex)
        with self.assertRaisesRegex(util.CytoflowError, "Couldn't parse"):
            flow.import_gatingml(FCS, self.ex)
        with self.assertRaisesRegex(util.CytoflowError, "Couldn't read"):
            flow.import_gatingml(self.tmp + "/nope.xml", self.ex)

    def testBadTransformParameters(self):
        """Bad logicle parameters used to kill the process"""
        body = '''
          <transforms:transformation transforms:id="bad">
            <transforms:logicle transforms:T="1000" transforms:W="3" transforms:M="4.5" transforms:A="0"/>
          </transforms:transformation>''' + self.polygon("Poly", "bad", "bad", [(0, 0), (1, 0), (1, 1)])
        ops, warns = self.do_import(flow.import_gatingml, self.write(body), self.ex)
        self.assertEqual(ops, [])
        self.assertIn("bad logicle transform parameters", warns)

    @unittest.skipIf(fk is None, "FlowKit isn't installed")
    def testHyperlogInverse(self):
        import flowutils
        x = np.array([-500.0, -1.0, 0.0, 3.0, 100.0, 5000.0, 90000.0])
        g = flowutils.transforms.hyperlog(x.copy(), [0], t = 100000.0, w = 0.5, m = 4.5, a = 0.3)
        np.testing.assert_allclose(gml._HyperlogTransform(100000.0, 0.5, 4.5, 0.3).inverse(g),
                                   x, rtol = 1e-6, atol = 1e-4)


class TestImportFlowJo(GateTest):
    """
    NB: synthetic.wsp was written by hand, following the structure of a FlowJo
    10 workspace.  It was NOT made by FlowJo.
    """

    # FlowJo path -> the cytoflow condition that should have the same events
    POPULATIONS = {
        "Cells" : "Cells",
        "Cells/Singlets" : "Singlets_gated",
        "Cells/Singlets/Y2+" : "Y2pos_gated",
        "Cells/Singlets/Y2-" : "Y2neg_gated",
        "Cells/Singlets/Bright" : "Bright_gated",
        "Cells/Singlets/Bright/Sub" : "Sub_gated",
        "Cells/Singlets/Ell" : "Ell_gated",
        "Cells/Singlets/EllLin" : "EllLin_gated",
        "Cells/Singlets/Outside" : "Outside_gated",
        "Cells/Singlets/Q1: Comp-V2-A- , Comp-Y2-A+" : "Q1_Comp_V2_Aneg_Comp_Y2_Apos_gated",
        "Cells/Singlets/Q2: Comp-V2-A+ , Comp-Y2-A+" : "Q2_Comp_V2_Apos_Comp_Y2_Apos_gated",
        "Cells/Singlets/Q3: Comp-V2-A+ , Comp-Y2-A-" : "Q3_Comp_V2_Apos_Comp_Y2_Aneg_gated",
        "Cells/Singlets/Q4: Comp-V2-A- , Comp-Y2-A-" : "Q4_Comp_V2_Aneg_Comp_Y2_Aneg_gated",
        "Cells/Singlets/NotEll" : "NotEll_gated",
        "Cells/Singlets/Either" : "Either_gated",
        "Cells/Singlets/Both" : "Both_gated"}

    def setUp(self):
        super().setUp()
        self.ops, self.warns = self.do_import(flow.import_flowjo_workspace, WSP, self.ex,
                                              sample = "RFP_Well_A3.fcs")
        self.by_name = {op.name : op for op in self.ops}

    def testSampleAndGroup(self):
        with self.assertRaisesRegex(util.CytoflowError,
                                    "2 samples.*'RFP_Well_A3.fcs', 'CFP_Well_A4.fcs'.*"
                                    "'All Samples', 'Test group', 'Empty group'"):
            flow.import_flowjo_workspace(WSP, self.ex)
        with self.assertRaisesRegex(util.CytoflowError, "no sample 'nope'.*RFP_Well_A3.fcs"):
            flow.import_flowjo_workspace(WSP, self.ex, sample = "nope")
        with self.assertRaisesRegex(util.CytoflowError, "no group 'nope'.*Test group"):
            flow.import_flowjo_workspace(WSP, self.ex, group = "nope")
        with self.assertRaisesRegex(util.CytoflowError, "has no samples"):
            flow.import_flowjo_workspace(WSP, self.ex, group = "Empty group")
        with self.assertRaisesRegex(util.CytoflowError, "not both"):
            flow.import_flowjo_workspace(WSP, self.ex, sample = "a", group = "b")
        with self.assertRaisesRegex(util.CytoflowError, "isn't a FlowJo workspace"):
            flow.export_gatingml([], self.ex, self.filename)
            flow.import_flowjo_workspace(self.filename, self.ex)

        ops, warns = self.do_import(flow.import_flowjo_workspace, WSP, self.ex, group = "Test group")
        self.assertEqual([(type(op), op.name) for op in ops], [(flow.RangeOp, "GroupGate")])

        # by sample ID
        ops, _ = self.do_import(flow.import_flowjo_workspace, WSP, self.ex, sample = "2")
        self.assertEqual([op.name for op in ops], ["GroupGate"])

        # no gates
        ops, warns = self.do_import(flow.import_flowjo_workspace, WSP, self.ex, group = "All Samples")
        self.assertEqual(ops, [])
        self.assertIn("Didn't find any gates", warns)

        self.assertEqual(gml.flowjo_workspace_samples(WSP),
                         (['RFP_Well_A3.fcs', 'CFP_Well_A4.fcs'],
                          ['All Samples', 'Test group', 'Empty group']))

    def testOps(self):
        kinds = {name : type(op) for name, op in self.by_name.items()}
        self.assertIs(kinds["Cells"], flow.PolygonOp)
        self.assertIs(kinds["Singlets"], flow.Range2DOp)
        self.assertIs(kinds["Y2pos"], flow.ThresholdOp)
        self.assertIs(kinds["Y2neg"], flow.RangeOp)
        self.assertIs(kinds["Q1_Comp_V2_Aneg_Comp_Y2_Apos"], flow.Range2DOp)
        for name in self.POPULATIONS.values():
            self.assertIn(name, kinds)
            if name.endswith("_gated"):
                self.assertIs(kinds[name], flow.BooleanGateOp)

        # linear and log axes: the polygon is imported as-is, in data units,
        # and renamed channels are found.
        cells = self.by_name["Cells"]
        self.assertEqual((cells.xchannel, cells.ychannel, cells.xscale, cells.yscale),
                         ("FSC-A", "SSC-A", "linear", "log"))
        np.testing.assert_allclose(cells.vertices,
                                   [(20.0, 300.0), (20.0, 20000.0), (20000.0, 90000.0),
                                    (50000.0, 90000.0), (50000.0, 2000.0), (2000.0, 300.0)])

        # compensated parameters -> the channel; biex -> logicle
        bright = self.by_name["Bright"]
        self.assertEqual((bright.xchannel, bright.ychannel, bright.xscale, bright.yscale),
                         ("Violet", "Y2-A", "logicle", "logicle"))
        self.assertEqual(len(bright.vertices), 6 * gml._EDGE_SEGMENTS)
        np.testing.assert_allclose(bright.vertices[0], VERTICES[0], rtol = 1e-6)
        np.testing.assert_allclose(bright.vertices[gml._EDGE_SEGMENTS], VERTICES[1], rtol = 1e-6)

        self.assertEqual(self.by_name["Singlets_gated"].expression, "Singlets and Cells")
        self.assertEqual(self.by_name["Sub_gated"].expression, "Sub and Bright_gated")
        self.assertEqual(self.by_name["Outside"].expression, "not Outside_region")
        self.assertEqual(self.by_name["NotEll"].expression, "not Ell_gated")
        self.assertEqual(self.by_name["Either"].expression, "Bright_gated or Ell_gated")
        self.assertEqual(self.by_name["Both"].expression, "Bright_gated and Y2pos_gated")

        # the things we said we didn't do
        self.assertRegex(self.warns, r" - Curly: FlowJo 'CurlyQuad' gates aren't supported")
        self.assertRegex(self.warns, r" - CurlyChild: its parent 'Curly' was not imported")
        self.assertIn("Compensation is NOT imported", self.warns)
        self.assertIn("Comp-V2-A", self.warns)
        self.assertIn("biexponential axis was approximated", self.warns)
        self.assertEqual(self.warns.count("not imported:"), 1)

    def testEvents(self):
        """Check against the definitions of the gates, worked out by hand"""
        ex = apply_all(self.ops, self.ex)
        fsc, ssc, v, y, b = [self.ex[c].values for c in ("FSC-A", "SSC-A", "Violet", "Y2-A", "B1-A")]

        # rectangles are in data units
        from matplotlib.path import Path
        cells = Path(np.column_stack(([20.0, 20.0, 20000.0, 50000.0, 50000.0, 2000.0],
                                      np.log10([300.0, 20000.0, 90000.0, 90000.0, 2000.0, 300.0])))) \
                .contains_points(np.column_stack((fsc, np.log10(ssc))))
        self.assertSameEvents(ex["Cells"], cells, tolerance = 1)

        singlets = cells & (fsc >= 30) & (fsc <= 40000) & (ssc >= 400) & (ssc <= 60000)
        self.assertSameEvents(ex["Singlets_gated"], singlets, tolerance = 1)
        self.assertSameEvents(ex["Y2pos_gated"], singlets & (y > 300), tolerance = 1)
        self.assertSameEvents(ex["Y2neg_gated"], singlets & (y <= 300), tolerance = 1)
        self.assertSameEvents(ex["Q3_Comp_V2_Apos_Comp_Y2_Aneg_gated"],
                              singlets & (v >= 100) & (y <= 1000), tolerance = 1)

        # the ellipses are on a 256 x 256 grid over the transformed axes:
        # center (110, 90), semi-axes 60 and 30, rotated by 30 degrees,
        # on Y2-A (logicle) vs B1-A (asinh).
        def inside(gx, gy, cx, cy, a, b, angle):
            t = np.deg2rad(angle)
            dx, dy = gx - cx, gy - cy
            return (((dx * np.cos(t) + dy * np.sin(t)) / a) ** 2 +
                    ((-dx * np.sin(t) + dy * np.cos(t)) / b) ** 2) <= 1

        gy2 = np.vectorize(Logicle(262144.0, 0.5, 4.5, 0.0).scale)(y) * 256
        gb1 = np.arcsinh(b * np.sinh(5 * np.log(10)) / 262144.0) / (5 * np.log(10)) * 256
        ell = inside(gy2, gb1, 110.0, 90.0, 60.0, 30.0, 30.0)
        self.assertSameEvents(ex["Ell"], ell)

        # and FSC-A (linear, 0 to 65536) vs Y2-A
        self.assertSameEvents(ex["EllLin"], inside(fsc / 256.0, gy2, 6.0, 100.0, 50.0, 10.0, 85.0))

        self.assertSameEvents(ex["NotEll_gated"], singlets & ~ell)
        self.assertSameEvents(ex["Outside_gated"], singlets & ~ex["Outside_region"].values, tolerance = 0)
        self.assertSameEvents(ex["Either_gated"],
                              singlets & (ex["Bright"].values | ex["Ell"].values), tolerance = 0)

    def flowkit_workspace(self, text = None):
        # FlowKit doesn't read CurlyQuad gates either
        if text is None:
            with open(WSP) as f:
                text = f.read()
        text, n = re.subn(r'<Population name="Curly".*?</Population>\s*</Subpopulations>\s*</Population>',
                          '', text, flags = re.S)
        self.assertEqual(n, 1)
        path = os.path.join(self.tmp, "flowkit.wsp")
        with open(path, 'w') as f:
            f.write(text)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ws = fk.Workspace(path, fcs_samples = [FCS], filename_as_id = True)
            ws.analyze_samples(sample_id = "RFP_Well_A3.fcs", use_mp = False)
        return ws

    @unittest.skipIf(fk is None, "FlowKit isn't installed")
    def testFlowKit(self):
        """Does FlowKit's FlowJo importer agree about the (synthetic) workspace?"""
        ws = self.flowkit_workspace()
        ex = apply_all(self.ops, self.ex)

        gate_ids = ws.get_gate_ids("RFP_Well_A3.fcs")
        self.assertEqual(len(gate_ids), len(self.POPULATIONS))
        for name, path in gate_ids:
            theirs = ws.get_gate_membership("RFP_Well_A3.fcs", name, gate_path = path)
            ours = ex[self.POPULATIONS["/".join(path[1:] + (name,))]]
            self.assertSameEvents(theirs, ours, name)

    @unittest.skipIf(fk is None, "FlowKit isn't installed")
    def testBiexIsNearlyLogicle(self):
        """The approximation that the importer makes for FlowJo's biex axes"""
        from flowkit._models.transforms._wsp_transforms import WSPBiexTransform
        for neg, width, pos in ((0, -10, 4.42), (0, -100, 4.42), (1, -10, 4.42),
                                (0.5, -50, 5.0), (0, -1000, 4.42)):
            biex = WSPBiexTransform(negative = neg, width = width, positive = pos,
                                    max_value = 262144.0)
            logicle = gml._LogicleTransform(262144.0, np.log10(-width) / 2.0, pos, neg)
            x = np.concatenate((-np.logspace(4, 0, 50), [0.0], np.logspace(0, np.log10(262144), 200)))
            ours = logicle.forward(x)
            theirs = biex.apply(x) / 4096.0
            on_scale = (ours > 0.005) & (theirs > 0.005)
            self.assertGreater(on_scale.sum(), 200)
            self.assertLess(np.abs(ours - theirs)[on_scale].max(), 0.006)

    def testUnsupported(self):
        with open(WSP) as f:
            text = f.read()

        # some other transform on Y2-A: rectangles are fine, nothing else is
        text = text.replace('transforms:logicle', 'transforms:miltenyi')
        # ellipses need a grid whose position we're sure of
        text = text.replace('transforms:minRange="0" transforms:maxRange="65536"',
                            'transforms:minRange="-100" transforms:maxRange="65536"')
        path = os.path.join(self.tmp, "unsupported.wsp")
        with open(path, 'w') as f:
            f.write(text)

        ops, warns = self.do_import(flow.import_flowjo_workspace, path, self.ex, sample = "1")
        names = [op.name for op in ops]
        for name in ("Cells", "Singlets", "Y2pos", "Y2neg", "Q1_Comp_V2_Aneg_Comp_Y2_Apos"):
            self.assertIn(name, names)
        for name, why in (("Bright", "FlowJo's 'miltenyi' axis transform isn't supported"),
                          ("Sub", "its parent 'Bright' was not imported"),
                          ("Ell", "miltenyi"),
                          ("EllLin", "miltenyi"),
                          ("Outside", "miltenyi"),
                          ("NotEll", "the gate 'Ell' was not imported"),
                          ("Either", "the gate 'Bright' was not imported")):
            self.assertRegex(warns, r" - {}: .*{}".format(name, re.escape(why)))
            self.assertNotIn(name, names)

        with open(WSP) as f:
            text = f.read().replace('transforms:minRange="0" transforms:maxRange="65536"',
                                    'transforms:minRange="-100" transforms:maxRange="65536"')
        with open(path, 'w') as f:
            f.write(text)
        ops, warns = self.do_import(flow.import_flowjo_workspace, path, self.ex, sample = "1")
        self.assertRegex(warns, r" - EllLin: .*its linear axis starts at -100")
        self.assertIn("Cells", [op.name for op in ops])
        self.assertIn("Ell", [op.name for op in ops])

    def testMatchingTransforms(self):
        """If cytoflow's scales match FlowJo's axes, gates keep their shapes"""
        # asinh(x sinh(M ln 10) / T) = asinh(x / cofactor)
        k = 5 * np.log(10)
        cofactor = float(262144.0 / np.sinh(k))
        ex = flow.ScaleSettingsOp(settings = {"Y2-A" : {"logicle" : {"W" : 0.5}},
                                              "Violet" : {"logicle" : {"W" : 1.0, "M" : 4.42}},
                                              "B1-A" : {"arcsinh" : {"cofactor" : cofactor}}}).apply(self.ex)
        ops, warns = self.do_import(flow.import_flowjo_workspace, WSP, ex, sample = "RFP_Well_A3.fcs")
        self.assertNotIn("follows the original outline", warns)
        by_name = {op.name : op for op in ops}
        self.assertEqual(len(by_name["Bright"].vertices), 6)
        np.testing.assert_allclose(by_name["Bright"].vertices, VERTICES, rtol = 1e-6)

        ell = by_name["Ell"]
        self.assertIsInstance(ell, flow.EllipseOp)
        self.assertEqual((ell.xchannel, ell.ychannel, ell.xscale, ell.yscale),
                         ("Y2-A", "B1-A", "logicle", "arcsinh"))
        # cytoflow's logicle is the same; the arcsinh scale is M ln 10 times bigger.
        self.assertAlmostEqual(ell.center[0], Logicle(262144.0, 0.5, 4.5, 0.0).inverse(110 / 256.0), places = 3)
        self.assertAlmostEqual(ell.center[1], cofactor * np.sinh(90 / 256.0 * k), places = 3)
        self.assertAlmostEqual(ell.width, 2 * 60 / 256.0, delta = 0.2)
        self.assertGreater(ell.height, 1.0)

        ex2 = apply_all(ops, ex)
        y, b = ex["Y2-A"].values, ex["B1-A"].values
        gx = np.vectorize(Logicle(262144.0, 0.5, 4.5, 0.0).scale)(y) * 256 - 110
        gy = np.arcsinh(b / cofactor) / k * 256 - 90
        t = np.deg2rad(30.0)
        inside = (((gx * np.cos(t) + gy * np.sin(t)) / 60) ** 2 +
                  ((-gx * np.sin(t) + gy * np.cos(t)) / 30) ** 2) <= 1
        self.assertSameEvents(ex2["Ell"], inside)

        if fk is not None:
            ws = self.flowkit_workspace()
            for name, path in ws.get_gate_ids("RFP_Well_A3.fcs"):
                theirs = ws.get_gate_membership("RFP_Well_A3.fcs", name, gate_path = path)
                ours = ex2[self.POPULATIONS["/".join(path[1:] + (name,))]]
                self.assertSameEvents(theirs, ours, name)


if __name__ == "__main__":
    unittest.main()
