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

"""
cytoflow.operations.gatingml
----------------------------

Exchange gates with other flow cytometry software.

`export_gatingml` -- write gate operations as a Gating-ML 2.0 document.

`import_gatingml` -- read the gates in a Gating-ML 2.0 document and return
the equivalent `cytoflow` operations.

`import_flowjo_workspace` -- read the gates of one sample (or group) in a
FlowJo 10 workspace (``.wsp``) and return the equivalent `cytoflow`
operations.

`is_gate` -- is this operation a gate that `export_gatingml` can write?

Coordinate spaces
~~~~~~~~~~~~~~~~~

There are three coordinate spaces to keep straight:

- *data* units -- the numbers in the FCS file (after whatever `cytoflow`
  operations have been applied to the `Experiment`).  `ThresholdOp`,
  `RangeOp`, `Range2DOp` and `QuadOp` are defined entirely in data units, as
  are `PolygonOp.vertices` and `EllipseOp.center`.

- `cytoflow`'s *scaled* space -- the data after an `IScale` (``linear``,
  ``log``, ``logicle``, ``arcsinh``) has been applied.  A polygon's edges are
  straight lines in this space, and an ellipse's width, height and angle are
  measured in it.

- Gating-ML's *transformed* space -- the data after one of the Gating-ML
  transforms.  Gating-ML gate coordinates on a transformed dimension are
  always in this space.

The scales and transforms correspond like this:

=============  =====================================  =========================
cytoflow       Gating-ML                              relationship
=============  =====================================  =========================
``linear``     ``flin(T, A = 0)``                     ``g = s / T``
``log``        ``flog(T, M = 4.5)``                   ``g = (s - log10(T)) / M + 1``
``logicle``    ``logicle(T, W, M, A)``                ``g = s``
``arcsinh``    ``fasinh(T = c sinh(M ln 10), M, 0)``  ``g = s / (M ln 10)``
=============  =====================================  =========================

(``s`` is the `cytoflow` scaled value, ``g`` the Gating-ML one, ``T`` the
channel's range, ``c`` the arcsinh cofactor.)  In every case the two spaces
differ only by a shift and a stretch of each axis, so straight lines stay
straight and ellipses stay ellipses, and the exported gate selects the same
events.

When *importing*, a transform whose parameters don't match what `cytoflow`
would use for that channel (eg, a logicle with a different ``W``) is not just
a shift and a stretch of a `cytoflow` scale.  In that case polygon edges are
sub-divided into many short segments, and ellipses are converted to 128-sided
polygons, so that the `PolygonOp` follows the original (curved) outline.
"""

import ast, keyword, re, warnings
import xml.etree.ElementTree as ET

import numpy as np

import cytoflow.utility as util
from cytoflow.utility.logicle_ext.Logicle import Logicle

from .threshold import ThresholdOp
from .range import RangeOp
from .range2d import Range2DOp
from .quad import QuadOp
from .polygon import PolygonOp
from .ellipse import EllipseOp, _normalize_axes
from .boolean_gate import BooleanGateOp

GATING_NS = "http://www.isac-net.org/std/Gating-ML/v2.0/gating"
TRANSFORMS_NS = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
DATATYPE_NS = "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes"

# the width, in decades, of the exported log and arcsinh transforms.  It's
# arbitrary (it doesn't change which events are in a gate), but 4.5 puts
# a typical channel in [0, 1].
_M = 4.5

# when a polygon's edges aren't straight in cytoflow's scaled space, how many
# segments to split each edge into
_EDGE_SEGMENTS = 32

# how many vertices to use when an ellipse must become a polygon
_ELLIPSE_VERTICES = 128

_LN10 = np.log(10.0)


class _Skip(Exception):
    """Raised when a gate can't be converted; the message says why."""


def is_gate(op):
    """
    Is ``op`` a gate that `export_gatingml` knows how to write?
    """
    return isinstance(op, (ThresholdOp, RangeOp, Range2DOp, QuadOp,
                           PolygonOp, EllipseOp, BooleanGateOp))


###############################################################################
# transforms
###############################################################################

class _Transform(object):
    """
    A Gating-ML (or FlowJo) axis transform, ``g = f(x)``.

    ``scale`` is the name of the most similar `cytoflow` scale.
    """

    scale = 'linear'

    # if set, the transform is only an approximation of the real thing
    approximate = None

    def forward(self, x):
        raise NotImplementedError

    def inverse(self, g):
        raise NotImplementedError

    def affine(self, scale):
        """
        If this transform is a shift and stretch of the `cytoflow` scale
        instance ``scale`` -- ie, ``g = alpha * scale(x) + beta`` -- return
        ``(alpha, beta)``.  Otherwise, return `None`.
        """
        return None


class _IdentityTransform(_Transform):
    def forward(self, x):
        return np.asarray(x, dtype = float)

    def inverse(self, g):
        return np.asarray(g, dtype = float)

    def affine(self, scale):
        return (1.0, 0.0)


class _LinTransform(_Transform):
    """Gating-ML ``flin``: ``g = (x + A) / (T + A)``"""

    def __init__(self, T, A):
        if not (np.isfinite(T) and np.isfinite(A) and T + A > 0):
            raise _Skip("bad linear transform parameters T = {}, A = {}".format(T, A))
        self.T = T
        self.A = A

    def forward(self, x):
        return (np.asarray(x, dtype = float) + self.A) / (self.T + self.A)

    def inverse(self, g):
        return np.asarray(g, dtype = float) * (self.T + self.A) - self.A

    def affine(self, scale):
        return (1.0 / (self.T + self.A), self.A / (self.T + self.A))


class _LogTransform(_Transform):
    """
    Gating-ML ``flog``: ``g = log10(x / T) / M + 1``.

    FlowJo's ``log`` transform, ``g = (log10(x) - log10(offset)) / decades``,
    is the same thing with ``M = decades`` and ``T = offset * 10 ** decades``;
    FlowJo also piles up everything below ``offset`` at ``offset``, which is
    what ``floor`` is for.
    """

    scale = 'log'

    def __init__(self, T, M, floor = None):
        if not (np.isfinite(T) and np.isfinite(M) and T > 0 and M > 0):
            raise _Skip("bad log transform parameters T = {}, M = {}".format(T, M))
        self.T = T
        self.M = M
        self.floor = floor

    def forward(self, x):
        x = np.asarray(x, dtype = float)
        if self.floor is not None:
            x = np.maximum(x, self.floor)
        with np.errstate(invalid = 'ignore', divide = 'ignore'):
            return np.log10(x / self.T) / self.M + 1.0

    def inverse(self, g):
        return self.T * np.power(10.0, (np.asarray(g, dtype = float) - 1.0) * self.M)

    def affine(self, scale):
        return (1.0 / self.M, 1.0 - np.log10(self.T) / self.M)


class _LogicleTransform(_Transform):
    """Gating-ML ``logicle``, which is the same function as `LogicleScale`"""

    scale = 'logicle'

    def __init__(self, T, W, M, A):
        # check here, because the C++ extension terminates the process
        # if it's given bad parameters!
        if not (np.all(np.isfinite([T, W, M, A]))
                and T > 0 and W >= 0 and M > 0 and 2 * W <= M
                and -A <= W and A + W <= M - W):
            raise _Skip("bad logicle transform parameters T = {}, W = {}, "
                        "M = {}, A = {}".format(T, W, M, A))
        self.T, self.W, self.M, self.A = T, W, M, A
        self._logicle = Logicle(float(T), float(W), float(M), float(A))
        self._forward = np.vectorize(self._logicle.scale, otypes = [float])
        self._inverse = np.vectorize(self._logicle.inverse, otypes = [float])

    def forward(self, x):
        return self._forward(np.asarray(x, dtype = float))

    def inverse(self, g):
        return self._inverse(np.asarray(g, dtype = float))

    def affine(self, scale):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            theirs = [scale._T, scale.W, scale.M, scale.A]
        if np.allclose([self.T, self.W, self.M, self.A], theirs,
                       rtol = 1e-6, atol = 1e-9):
            return (1.0, 0.0)
        return None


class _AsinhTransform(_Transform):
    """
    Gating-ML ``fasinh``:
    ``g = (asinh(x sinh(M ln10) / T) + A ln10) / ((M + A) ln10)``
    """

    scale = 'arcsinh'

    def __init__(self, T, M, A):
        if not (np.all(np.isfinite([T, M, A])) and T > 0 and M > 0 and M + A > 0):
            raise _Skip("bad asinh transform parameters T = {}, M = {}, A = {}"
                        .format(T, M, A))
        self.T, self.M, self.A = T, M, A
        self.cofactor = T / np.sinh(M * _LN10)

    def forward(self, x):
        return ((np.arcsinh(np.asarray(x, dtype = float) / self.cofactor) + self.A * _LN10)
                / ((self.M + self.A) * _LN10))

    def inverse(self, g):
        return self.cofactor * np.sinh(np.asarray(g, dtype = float) * (self.M + self.A) * _LN10
                                       - self.A * _LN10)

    def affine(self, scale):
        if np.isclose(self.cofactor, scale.cofactor, rtol = 1e-6, atol = 0):
            return (1.0 / ((self.M + self.A) * _LN10), self.A / (self.M + self.A))
        return None


class _HyperlogTransform(_Transform):
    """
    Gating-ML ``hyperlog``.  Only the inverse (which is all that importing
    a Gating-ML gate needs) is implemented; it has a closed form.
    """

    scale = 'logicle'

    def __init__(self, T, W, M, A):
        if not (np.all(np.isfinite([T, W, M, A]))
                and T > 0 and W > 0 and M > 0 and 2 * W <= M
                and -A <= W and A + W <= M - W):
            raise _Skip("bad hyperlog transform parameters T = {}, W = {}, "
                        "M = {}, A = {}".format(T, W, M, A))
        w = W / (M + A)
        x2 = A / (M + A)
        self.x1 = x2 + w
        x0 = x2 + 2 * w
        self.b = (M + A) * _LN10
        c_a = np.exp(self.b * x0) / w
        f_a = np.exp(self.b * self.x1) + c_a * self.x1
        self.a = T / (np.exp(self.b) + c_a - f_a)
        self.c = c_a * self.a
        self.f = f_a * self.a

    def inverse(self, g):
        g = np.asarray(g, dtype = float)
        negative = g < self.x1
        g = np.where(negative, 2 * self.x1 - g, g)
        x = self.a * np.exp(self.b * g) + self.c * g - self.f
        return np.where(negative, -x, x)


def _silent_scale(name, experiment, channel):
    """Make a cytoflow scale, without the scale's warnings"""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return util.scale_factory(name, experiment, channel = channel)


###############################################################################
# export
###############################################################################

def _q(ns, name):
    return "{{{}}}{}".format(ns, name)


def _f(value):
    return repr(float(value))


def _channel_range(experiment, channel):
    """The `T` parameter cytoflow uses for a channel: see `LogicleScale`"""
    r = experiment.metadata.get(channel, {}).get('range', None)
    try:
        r = float(r)
    except (TypeError, ValueError):
        r = float('nan')
    if not (np.isfinite(r) and r > 0):
        r = float(experiment.data[channel].max())
    if not (np.isfinite(r) and r > 0):
        r = 1.0
    return r


class _Exporter(object):

    def __init__(self, names):
        self.transforms = []      # the <transformation> elements
        self.transform_ids = {}   # (kind, params) -> id
        self.gates = []           # the gate elements
        self.kinds = {}           # exported op name -> 'bool' or 'quad'
        self.ids = set(names)     # every ID that is, or will be, used

    ### transforms

    def transform(self, kind, **params):
        key = (kind, tuple(sorted(params.items())))
        if key not in self.transform_ids:
            tid = "{}_{}".format(kind, len(self.transform_ids) + 1)
            el = ET.Element(_q(TRANSFORMS_NS, 'transformation'),
                            {_q(TRANSFORMS_NS, 'id') : tid})
            ET.SubElement(el, _q(TRANSFORMS_NS, kind),
                          {_q(TRANSFORMS_NS, k) : _f(v) for k, v in params.items()})
            self.transforms.append(el)
            self.transform_ids[key] = tid
        return self.transform_ids[key]

    def axis(self, experiment, channel, scale_name):
        """
        Returns ``(fcs_name, transform_id, to_gml, alpha)``.  ``to_gml`` maps
        a *data* value to the Gating-ML transformed space; ``alpha`` is the
        stretch between cytoflow's scaled space and Gating-ML's.
        """
        if channel not in experiment.channels:
            raise _Skip("channel '{}' isn't in the experiment".format(channel))

        scale = _silent_scale(scale_name, experiment, channel)

        if scale_name == 'linear':
            T = _channel_range(experiment, channel)
            tid = self.transform('flin', T = T, A = 0.0)
            alpha, beta = 1.0 / T, 0.0
        elif scale_name == 'log':
            T = _channel_range(experiment, channel)
            tid = self.transform('flog', T = T, M = _M)
            alpha, beta = 1.0 / _M, 1.0 - np.log10(T) / _M
        elif scale_name == 'logicle':
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                T, W, M, A = scale._T, scale.W, scale.M, scale.A
            tid = self.transform('logicle', T = T, W = W, M = M, A = A)
            alpha, beta = 1.0, 0.0

            # cytoflow's logicle scale uses a lookup table with 4096 bins,
            # which is plenty for plotting and gating but would move the
            # gate's vertices by up to a bin.  compute them exactly instead.
            exact = _LogicleTransform(T, W, M, A)
            def to_gml(x):
                return float(exact.forward(scale.clip(float(x))))
            return self.fcs_name(experiment, channel), tid, to_gml, alpha
        elif scale_name == 'arcsinh':
            tid = self.transform('fasinh',
                                 T = scale.cofactor * np.sinh(_M * _LN10),
                                 M = _M,
                                 A = 0.0)
            alpha, beta = 1.0 / (_M * _LN10), 0.0
        else:
            raise _Skip("Gating-ML has no equivalent of the '{}' scale"
                        .format(scale_name))

        def to_gml(x):
            try:
                s = float(scale(float(x)))
            except util.CytoflowError as e:
                raise _Skip("{} is outside the domain of the {} scale"
                            .format(x, scale_name)) from e
            if not np.isfinite(s):
                raise _Skip("{} is outside the domain of the {} scale"
                            .format(x, scale_name))
            return alpha * s + beta

        return self.fcs_name(experiment, channel), tid, to_gml, alpha

    ### elements

    @staticmethod
    def fcs_name(experiment, channel):
        if channel not in experiment.channels:
            raise _Skip("channel '{}' isn't in the experiment".format(channel))
        return experiment.metadata.get(channel, {}).get('fcs_name', channel)

    @staticmethod
    def dimension(parent, fcs_name, tid = None, lo = None, hi = None, tag = 'dimension',
                  **attrs):
        attrib = {_q(GATING_NS, k) : v for k, v in attrs.items()}
        attrib[_q(GATING_NS, 'compensation-ref')] = 'uncompensated'
        if tid is not None:
            attrib[_q(GATING_NS, 'transformation-ref')] = tid
        if lo is not None and np.isfinite(lo):
            attrib[_q(GATING_NS, 'min')] = _f(lo)
        if hi is not None and np.isfinite(hi):
            attrib[_q(GATING_NS, 'max')] = _f(hi)
        el = ET.SubElement(parent, _q(GATING_NS, tag), attrib)
        ET.SubElement(el, _q(DATATYPE_NS, 'fcs-dimension'), {_q(DATATYPE_NS, 'name') : fcs_name})
        return el

    @staticmethod
    def gate(tag, gate_id):
        return ET.Element(_q(GATING_NS, tag), {_q(GATING_NS, 'id') : gate_id})

    @staticmethod
    def value(parent, tag, value):
        return ET.SubElement(parent, _q(GATING_NS, tag), {_q(DATATYPE_NS, 'value') : _f(value)})

    ### gates

    def add(self, op, experiment):
        if not op.name:
            raise _Skip("it has no name")

        if isinstance(op, ThresholdOp):
            if op.threshold is None:
                raise _Skip("the threshold isn't set")
            el = self.gate('RectangleGate', op.name)
            self.dimension(el, self.fcs_name(experiment, op.channel), lo = op.threshold)
            self.gates.append(el)
            self.kinds[op.name] = 'bool'

        elif isinstance(op, RangeOp):
            if op.low is None or op.high is None:
                raise _Skip("the range isn't set")
            el = self.gate('RectangleGate', op.name)
            self.dimension(el, self.fcs_name(experiment, op.channel), lo = op.low, hi = op.high)
            self.gates.append(el)
            self.kinds[op.name] = 'bool'

        elif isinstance(op, Range2DOp):
            if None in (op.xlow, op.xhigh, op.ylow, op.yhigh):
                raise _Skip("the range isn't set")
            el = self.gate('RectangleGate', op.name)
            self.dimension(el, self.fcs_name(experiment, op.xchannel), lo = op.xlow, hi = op.xhigh)
            self.dimension(el, self.fcs_name(experiment, op.ychannel), lo = op.ylow, hi = op.yhigh)
            self.gates.append(el)
            self.kinds[op.name] = 'bool'

        elif isinstance(op, QuadOp):
            if op.xthreshold is None or op.ythreshold is None:
                raise _Skip("the thresholds aren't set")
            el = self.gate('QuadrantGate', op.name)
            x_id = self.unique_id(op.name + "_x")
            y_id = self.unique_id(op.name + "_y")
            for div_id, channel, threshold in ((x_id, op.xchannel, op.xthreshold),
                                               (y_id, op.ychannel, op.ythreshold)):
                div = self.dimension(el, self.fcs_name(experiment, channel),
                                     tag = 'divider', id = div_id)
                ET.SubElement(div, _q(GATING_NS, 'value')).text = _f(threshold)

            # QuadOp's numbering: 1 = upper left, 2 = upper right,
            # 3 = lower left, 4 = lower right
            for i, (xsign, ysign) in enumerate(((-1, 1), (1, 1), (-1, -1), (1, -1))):
                quadrant = ET.SubElement(el, _q(GATING_NS, 'Quadrant'),
                                         {_q(GATING_NS, 'id') : "{}_{}".format(op.name, i + 1)})
                for div_id, threshold, sign in ((x_id, op.xthreshold, xsign),
                                                (y_id, op.ythreshold, ysign)):
                    # the location is any value on the right side of the divider
                    location = threshold + sign * max(1.0, abs(threshold) * 0.5)
                    ET.SubElement(quadrant, _q(GATING_NS, 'position'),
                                  {_q(GATING_NS, 'divider_ref') : div_id,
                                   _q(GATING_NS, 'location') : _f(location)})
            self.gates.append(el)
            self.kinds[op.name] = 'quad'

        elif isinstance(op, PolygonOp):
            if len(op.vertices) < 3:
                raise _Skip("it has fewer than three vertices")
            xname, xtid, xgml, _ = self.axis(experiment, op.xchannel, op.xscale)
            yname, ytid, ygml, _ = self.axis(experiment, op.ychannel, op.yscale)
            vertices = [(xgml(x), ygml(y)) for x, y in op.vertices]

            el = self.gate('PolygonGate', op.name)
            self.dimension(el, xname, xtid)
            self.dimension(el, yname, ytid)
            for vertex in vertices:
                v = ET.SubElement(el, _q(GATING_NS, 'vertex'))
                self.value(v, 'coordinate', vertex[0])
                self.value(v, 'coordinate', vertex[1])
            self.gates.append(el)
            self.kinds[op.name] = 'bool'

        elif isinstance(op, EllipseOp):
            if not (op.width > 0 and op.height > 0):
                raise _Skip("its width and height aren't set")
            xname, xtid, xgml, xalpha = self.axis(experiment, op.xchannel, op.xscale)
            yname, ytid, ygml, yalpha = self.axis(experiment, op.ychannel, op.yscale)

            # the ellipse, in cytoflow's scaled space, is
            #    (s - c)' inv(cov_s) (s - c) <= 1
            # with cov_s = R diag((w / 2) ^ 2, (h / 2) ^ 2) R'.
            # Gating-ML's space is g = alpha * s + beta on each axis, so
            # cov_g = diag(alpha) cov_s diag(alpha), and distanceSquare = 1
            t = np.deg2rad(op.angle)
            r = np.array([[np.cos(t), -np.sin(t)],
                          [np.sin(t), np.cos(t)]])
            cov_s = r @ np.diag([(op.width / 2.0) ** 2, (op.height / 2.0) ** 2]) @ r.T
            stretch = np.diag([xalpha, yalpha])
            cov_g = stretch @ cov_s @ stretch

            el = self.gate('EllipsoidGate', op.name)
            self.dimension(el, xname, xtid)
            self.dimension(el, yname, ytid)
            mean = ET.SubElement(el, _q(GATING_NS, 'mean'))
            self.value(mean, 'coordinate', xgml(op.center[0]))
            self.value(mean, 'coordinate', ygml(op.center[1]))
            cov = ET.SubElement(el, _q(GATING_NS, 'covarianceMatrix'))
            for i in range(2):
                row = ET.SubElement(cov, _q(GATING_NS, 'row'))
                for j in range(2):
                    self.value(row, 'entry', cov_g[i, j])
            self.value(el, 'distanceSquare', 1.0)
            self.gates.append(el)
            self.kinds[op.name] = 'bool'

        elif isinstance(op, BooleanGateOp):
            try:
                tree = ast.parse(op.expression.strip(), mode = 'eval').body
            except (SyntaxError, ValueError, RecursionError, MemoryError) as e:
                raise _Skip("couldn't parse the expression '{}'".format(op.expression)) from e

            # convert first, so we don't emit half a gate
            tree = self.bool_tree(tree)
            self.bool_emit(tree, op.name, op.name)
            self.kinds[op.name] = 'bool'

        else:
            raise _Skip("{} isn't a gate that Gating-ML can represent"
                        .format(type(op).__name__))

    def unique_id(self, base):
        new_id = base
        i = 1
        while new_id in self.ids:
            i += 1
            new_id = "{}_{}".format(base, i)
        self.ids.add(new_id)
        return new_id

    ### boolean gates

    def bool_tree(self, node):
        """
        Convert a `BooleanGateOp` expression to a tree of
        ``('ref', gate_id, complement)``, ``('and', [...])``, ``('or', [...])``
        and ``('not', child)``
        """
        if isinstance(node, ast.BoolOp):
            return ('and' if isinstance(node.op, ast.And) else 'or',
                    [self.bool_tree(v) for v in node.values])

        elif isinstance(node, ast.BinOp) and isinstance(node.op, (ast.BitAnd, ast.BitOr)):
            return ('and' if isinstance(node.op, ast.BitAnd) else 'or',
                    [self.bool_tree(node.left), self.bool_tree(node.right)])

        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.Not, ast.Invert)):
            child = self.bool_tree(node.operand)
            if child[0] == 'ref':
                return ('ref', child[1], not child[2])
            return ('not', child)

        elif isinstance(node, ast.Name):
            if self.kinds.get(node.id) != 'bool':
                raise _Skip("its expression uses '{}', which isn't an exported "
                            "boolean gate".format(node.id))
            return ('ref', node.id, False)

        elif isinstance(node, ast.Compare) and len(node.ops) == 1:
            left, right, op = node.left, node.comparators[0], node.ops[0]
            if isinstance(left, ast.Constant) and isinstance(right, ast.Name):
                left, right = right, left
            if (isinstance(left, ast.Name) and isinstance(right, ast.Constant)
                    and isinstance(op, (ast.Eq, ast.NotEq))):
                negate = isinstance(op, ast.NotEq)
                kind = self.kinds.get(left.id)
                if (kind == 'bool' and not isinstance(right.value, str)
                        and right.value in (True, False)):
                    return ('ref', left.id, negate != (not right.value))
                elif (kind == 'quad' and isinstance(right.value, str)
                        and right.value in ["{}_{}".format(left.id, i) for i in range(1, 5)]):
                    return ('ref', right.value, negate)

        raise _Skip("its expression contains '{}', which Gating-ML can't "
                    "represent".format(ast.unparse(node)))

    def bool_operand(self, tree, name):
        """Returns a ``(gate_id, complement)`` for ``tree``, adding helper gates as needed"""
        if tree[0] == 'ref':
            return tree[1], tree[2]

        helper = self.unique_id(name + "_part")
        self.bool_emit(tree, helper, name)
        return helper, False

    def bool_emit(self, tree, gate_id, name):
        # a Gating-ML boolean gate is a single and / or / not, whose operands
        # are all references to other gates.  nested sub-expressions must
        # become gates of their own.
        if tree[0] == 'ref':
            # just another name for a gate: a double negative
            operator, operands = 'not', [(tree[1], not tree[2])]
        elif tree[0] == 'not':
            operator, operands = 'not', [self.bool_operand(tree[1], name)]
        else:
            operator, operands = tree[0], [self.bool_operand(c, name) for c in tree[1]]

        el = self.gate('BooleanGate', gate_id)
        op_el = ET.SubElement(el, _q(GATING_NS, operator))
        for ref, complement in operands:
            attrib = {_q(GATING_NS, 'ref') : ref}
            if complement:
                attrib[_q(GATING_NS, 'use-as-complement')] = 'true'
            ET.SubElement(op_el, _q(GATING_NS, 'gateReference'), attrib)
        self.gates.append(el)


def export_gatingml(ops, experiment, filename):
    """
    Write gates to a Gating-ML 2.0 file.  Gating-ML is the ISAC standard for
    exchanging gates; it's supported (to varying degrees) by FlowKit, the
    R/Bioconductor flow cytometry packages (``CytoML``, ``flowCore``),
    FlowJo, Cytobank, and others.  The files this function writes have been
    checked with FlowKit.

    Parameters
    ----------
    ops : list of `IOperation`
        The gates to export.  `ThresholdOp`, `RangeOp` and `Range2DOp` are
        written as a ``RectangleGate``; `PolygonOp` and `FreehandOp` as a
        ``PolygonGate``; `EllipseOp` as an ``EllipsoidGate``; `QuadOp` as a
        ``QuadrantGate`` (whose quadrants have the same names as the `QuadOp`
        categories, eg ``Quad_1``); and `BooleanGateOp` as one or more
        ``BooleanGate`` s.  Anything else is skipped, with a warning.

    experiment : `Experiment`, or a list of `Experiment`
        The experiment the gates are applied to.  It's used to look up the
        channels' names in the FCS files, and the parameters of the gates'
        scales (some of which, like the logicle ``W``, are estimated from
        the data.)  If different gates are applied to different experiments
        (ie, a channel's data changes between them), pass a list with one
        experiment per operation.

    filename : str or path
        The file to write.

    Returns
    -------
    list of str
        The names of the gates that were written.

    Notes
    -----
    - Gate names are the gates' IDs.  `cytoflow` gates aren't nested, so
      no gate has a ``parent_id``; a gate that should only include events
      in some other gate is a `BooleanGateOp`, which is written as a
      ``BooleanGate``.

    - Dimensions are named by the channels' names in the FCS file (``$PnN``),
      and are marked ``uncompensated``: the gates apply to the data *as it
      is in the experiment*.  If you applied compensation (or any other
      operation that changes the data) in `cytoflow`, the other program needs
      to do the same before it applies the gates.

    - `PolygonOp` and `EllipseOp` are written in the Gating-ML transformed
      space that corresponds to their scales -- see the module documentation,
      `cytoflow.operations.gatingml`.  The other gates don't have a scale,
      and are written in untransformed data units.

    - A `BooleanGateOp` expression with nested parentheses becomes several
      ``BooleanGate`` s (the extra ones are named ``<name>_part``.)  Only
      expressions over exported gates can be written; comparisons with
      other conditions (``Dox > 1``) can't.

    - `cytoflow` gates with ``adjust_by`` set ("magnetic" gates) are written
      at their un-adjusted position.

    - There are small differences at the very edges of gates and scales:
      Gating-ML rectangles include their lower bound and exclude their upper
      one, while `RangeOp` includes both and `ThresholdOp` and `QuadOp`
      exclude an event exactly on the threshold; and `cytoflow`'s logicle
      scale piles up events beyond the ends of the scale at its edges, where
      Gating-ML doesn't.
    """

    ops = list(ops)
    if isinstance(experiment, (list, tuple)):
        experiments = list(experiment)
        if len(experiments) != len(ops):
            raise util.CytoflowError("If 'experiment' is a list, it must have "
                                     "one experiment for each operation")
    else:
        experiments = [experiment] * len(ops)

    if any(ex is None for ex in experiments):
        raise util.CytoflowError("No experiment specified")

    exporter = _Exporter([getattr(op, 'name', '') for op in ops])
    exported = []
    skipped = []
    magnetic = []

    for op, ex in zip(ops, experiments):
        name = getattr(op, 'name', '') or type(op).__name__
        try:
            if name in exported:
                raise _Skip("there is already a gate with that name")
            exporter.add(op, ex)
            exported.append(name)
            if getattr(op, 'adjust_by', None):
                magnetic.append(name)
        except _Skip as e:
            skipped.append("{}: {}".format(name, e))

    for prefix, ns in (('gating', GATING_NS),
                       ('transforms', TRANSFORMS_NS),
                       ('data-type', DATATYPE_NS)):
        ET.register_namespace(prefix, ns)

    root = ET.Element(_q(GATING_NS, 'Gating-ML'))
    root.extend(exporter.transforms)
    root.extend(exporter.gates)
    tree = ET.ElementTree(root)
    ET.indent(tree)
    tree.write(filename, encoding = 'UTF-8', xml_declaration = True)

    if skipped:
        warnings.warn("These operations were not exported:\n" + _bullets(skipped),
                      util.CytoflowWarning)
    if magnetic:
        warnings.warn("These gates are adjusted per-sample ('adjust_by'); "
                      "they were exported at their un-adjusted position: {}"
                      .format(", ".join(magnetic)),
                      util.CytoflowWarning)

    return exported


def _bullets(lines):
    return "\n".join([" - " + x for x in lines])


###############################################################################
# import: things shared by Gating-ML and FlowJo
###############################################################################

def _local(tag):
    """An element or attribute name, without its namespace"""
    return tag.rsplit('}', 1)[-1] if isinstance(tag, str) else ''


def _attr(el, name, default = None):
    """Get an attribute, whatever its namespace"""
    for k, v in el.attrib.items():
        if _local(k) == name:
            return v
    return default


def _children(el, name):
    return [c for c in el if _local(c.tag) == name]


def _child(el, name):
    c = _children(el, name)
    return c[0] if c else None


def _float(text, what):
    try:
        value = float(text)
    except (TypeError, ValueError):
        raise _Skip("{} is '{}', which isn't a number".format(what, text))
    if not np.isfinite(value):
        raise _Skip("{} isn't finite".format(what))
    return value


def _parse(filename):
    try:
        return ET.parse(filename).getroot()
    except ET.ParseError as e:
        raise util.CytoflowError("Couldn't parse {}: {}".format(filename, e)) from e
    except OSError as e:
        raise util.CytoflowError("Couldn't read {}: {}".format(filename, e)) from e


class _Dim(object):
    """One dimension of a gate"""
    def __init__(self, name, transform = None, lo = None, hi = None):
        self.name = name             # the name in the FCS file
        self.transform = transform   # a _Transform
        self.lo = lo                 # for rectangles: bounds in DATA units
        self.hi = hi


class _Gate(object):
    """
    A gate from a file, part-way to becoming a `cytoflow` operation. `kind`
    is one of

    - ``rect``: uses `dims` (with their bounds)
    - ``poly``: uses `dims` and `vertices` (in the transformed space)
    - ``ellipse``: uses `dims`, `mean` and `cov` (in the transformed space)
    - ``bool``: uses `operator` and `refs`, a list of (gate id, complement)
    - ``quad``: uses `dims`, `thresholds` (data units) -- becomes a `QuadOp`
    - ``quadrant``: uses `owner` and `value`; a quadrant of a ``quad``
    """
    def __init__(self, gate_id, name, kind = None, parent = None, error = None, **kwargs):
        self.id = gate_id
        self.name = name
        self.kind = kind
        self.parent = parent
        self.error = error     # if set, why the gate can't be imported
        self.complement = False
        self.owner = None
        self.alt_name = None   # another name to try if `name` is taken
        self.__dict__.update(kwargs)


class _Namer(object):
    """Makes unique names that are legal names of cytoflow conditions"""

    def __init__(self, experiment):
        self.used = (set(experiment.data.columns) | set(experiment.conditions)
                     | set([util.sanitize_identifier(c) for c in experiment.channels]))

    @staticmethod
    def clean(name):
        # "CD4+" and "CD4-" shouldn't both become "CD4_"
        name = re.sub(r'\+(?=$|[\s,/])', 'pos', name.strip())
        name = re.sub(r'(?<=\w)-(?=$|[\s,/])', 'neg', name)
        name = util.sanitize_identifier(name) if name and name[0].isalpha() \
               else util.sanitize_identifier("Gate_" + name)
        name = re.sub(r'_+', '_', name).rstrip('_') or "Gate"
        if keyword.iskeyword(name) or name in ('True', 'False', 'None'):
            name += "_gate"
        return name

    def __call__(self, name, alt = None):
        for candidate in (name, alt):
            if candidate and self.clean(candidate) not in self.used:
                self.used.add(self.clean(candidate))
                return self.clean(candidate)
        base = self.clean(name)
        i = 2
        while "{}_{}".format(base, i) in self.used:
            i += 1
        self.used.add("{}_{}".format(base, i))
        return "{}_{}".format(base, i)


def _find_channel(experiment, fcs_name):
    for channel in experiment.channels:
        if experiment.metadata.get(channel, {}).get('fcs_name', None) == fcs_name:
            return channel
    if fcs_name in experiment.channels:
        return fcs_name
    raise _Skip("the experiment doesn't have a channel '{}'".format(fcs_name))


class _Builder(object):
    """Turns a dict of `_Gate` into a list of operations"""

    def __init__(self, gates, experiment):
        self.gates = gates
        self.experiment = experiment
        self.namer = _Namer(experiment)
        self.ops = []
        self.done = {}       # gate id -> its expression (hierarchy included), or None if skipped
        self.raw = {}        # gate id -> name of a quad op
        self.skipped = []
        self.notes = []

    def note(self, msg):
        if msg not in self.notes:
            self.notes.append(msg)

    def build(self):
        for gate_id in self.gates:
            self.resolve(gate_id, ())

        if self.skipped:
            warnings.warn("These gates were not imported:\n" + _bullets(self.skipped),
                          util.CytoflowWarning)
        for msg in self.notes:
            warnings.warn(msg, util.CytoflowWarning)

        return self.ops

    def resolve(self, gate_id, stack):
        """
        Make the operations for a gate (and the gates it depends on); return
        an expression for the gate's events, or None if it was skipped.
        """
        if gate_id in self.done:
            return self.done[gate_id]

        gate = self.gates[gate_id]
        try:
            if gate_id in stack:
                raise _Skip("it depends on itself")
            if gate.error:
                raise _Skip(gate.error)

            # the things it depends on come first
            parent = self.depend(gate.parent, stack + (gate_id,), "its parent") if gate.parent else None
            if gate.kind == 'bool':
                refs = [(self.depend(ref, stack + (gate_id,), "the gate"), complement)
                        for ref, complement in gate.refs]
            elif gate.kind == 'quadrant':
                self.depend(gate.owner, stack + (gate_id,), "its quadrant gate")

            n_ops = len(self.ops)
            try:
                expr = self.make(gate, refs if gate.kind == 'bool' else None)
            except Exception:
                del self.ops[n_ops:]
                raise

            if expr is not None and gate.complement:
                # the events OUTSIDE the shape
                op = self.ops[-1]
                op.name = self.namer(expr + "_region")
                self.ops.append(BooleanGateOp(name = expr,
                                              expression = "not " + op.name))

            if expr is not None and parent:
                # cytoflow gates aren't nested
                name = self.namer((expr if expr.isidentifier() else gate.name) + "_gated")
                self.ops.append(BooleanGateOp(name = name,
                                              expression = "{} and {}".format(expr, parent)))
                expr = name

            self.done[gate_id] = expr if expr is not None else ''

        except _Skip as e:
            self.skipped.append("{}: {}".format(gate.name, e))
            self.done[gate_id] = None

        return self.done[gate_id]

    def depend(self, gate_id, stack, what):
        if gate_id not in self.gates:
            raise _Skip("{} '{}' isn't in the file".format(what, gate_id))
        expr = self.resolve(gate_id, stack)
        if expr is None:
            raise _Skip("{} '{}' was not imported".format(what, self.gates[gate_id].name))
        if expr == '' and what != "its quadrant gate":
            raise _Skip("it refers to the quadrant gate '{}' instead of one of "
                        "its quadrants".format(self.gates[gate_id].name))
        return expr

    def make(self, gate, refs):
        if gate.kind == 'rect':
            return self.make_rect(gate)
        elif gate.kind == 'poly':
            return self.make_polygon(gate, gate.vertices, gate.dims)
        elif gate.kind == 'ellipse':
            return self.make_ellipse(gate)
        elif gate.kind == 'bool':
            return self.make_bool(gate, refs)
        elif gate.kind == 'quad':
            if gate.id in self.namer.used or self.namer.clean(gate.id) != gate.id:
                # the quadrants' names wouldn't match the QuadOp's categories:
                # make a rectangle gate for each quadrant instead.
                self.raw[gate.id] = None
                return None
            name = self.namer(gate.name)
            self.ops.append(QuadOp(name = name,
                                   xchannel = _find_channel(self.experiment, gate.dims[0].name),
                                   ychannel = _find_channel(self.experiment, gate.dims[1].name),
                                   xthreshold = gate.thresholds[0],
                                   ythreshold = gate.thresholds[1]))
            self.raw[gate.id] = name
            return None
        elif gate.kind == 'quadrant':
            if self.raw[gate.owner] is None:
                return self.make_rect(gate)
            return "({} == '{}')".format(self.raw[gate.owner], gate.value)
        else:
            raise _Skip("unknown kind of gate")

    def make_rect(self, gate):
        dims = gate.dims
        if len(dims) == 0:
            raise _Skip("it has no dimensions")
        if len(dims) > 2:
            raise _Skip("it has {} dimensions; only 1- and 2-dimensional "
                        "gates are supported".format(len(dims)))
        for d in dims:
            if d.lo is None and d.hi is None:
                raise _Skip("dimension '{}' has neither a minimum nor a maximum"
                            .format(d.name))
            if d.lo is not None and d.hi is not None and d.lo > d.hi:
                d.lo, d.hi = d.hi, d.lo

        channels = [_find_channel(self.experiment, d.name) for d in dims]
        lo = [d.lo if d.lo is not None else -np.inf for d in dims]
        hi = [d.hi if d.hi is not None else np.inf for d in dims]
        name = self.namer(gate.name, gate.alt_name)

        if len(dims) == 1 and dims[0].hi is None:
            op = ThresholdOp(name = name, channel = channels[0], threshold = lo[0])
        elif len(dims) == 1:
            op = RangeOp(name = name, channel = channels[0], low = lo[0], high = hi[0])
        else:
            op = Range2DOp(name = name,
                           xchannel = channels[0], xlow = lo[0], xhigh = hi[0],
                           ychannel = channels[1], ylow = lo[1], yhigh = hi[1])
        self.ops.append(op)
        return name

    def axes(self, gate, dims):
        if len(dims) != 2:
            raise _Skip("it has {} dimension{}; polygon and ellipse gates must "
                        "have 2".format(len(dims), "" if len(dims) == 1 else "s"))
        channels = [_find_channel(self.experiment, d.name) for d in dims]
        if channels[0] == channels[1]:
            raise _Skip("both of its dimensions are '{}'".format(dims[0].name))
        transforms = [d.transform or _IdentityTransform() for d in dims]
        scales = [_silent_scale(t.scale, self.experiment, c)
                  for t, c in zip(transforms, channels)]
        affine = [t.affine(s) for t, s in zip(transforms, scales)]
        for t, d in zip(transforms, dims):
            if t.approximate:
                self.note("Channel '{}': {}".format(d.name, t.approximate))
        return channels, transforms, affine

    def make_polygon(self, gate, vertices, dims, was = "polygon"):
        channels, transforms, affine = self.axes(gate, dims)
        vertices = np.asarray(vertices, dtype = float)
        if vertices.ndim != 2 or vertices.shape[0] < 3 or vertices.shape[1] != 2:
            raise _Skip("it doesn't have at least three 2-dimensional vertices")

        if None in affine and was == "polygon":
            # the edges aren't straight lines in cytoflow's scaled space.
            # follow them with lots of little segments
            nxt = np.roll(vertices, -1, axis = 0)
            f = np.linspace(0.0, 1.0, _EDGE_SEGMENTS, endpoint = False)
            vertices = (vertices[:, None, :]
                        + f[None, :, None] * (nxt - vertices)[:, None, :]).reshape(-1, 2)

        if None in affine:
            which = [d.name for d, a in zip(dims, affine) if a is None]
            self.note("The {} gate '{}' was drawn with an axis transform on {} "
                      "that is different from cytoflow's {} scale, so its outline "
                      "isn't {} on cytoflow's plots. It was imported as a "
                      "polygon gate with {} vertices that follows the original "
                      "outline."
                      .format(was, gate.name, " and ".join(which),
                              " / ".join(sorted(set([t.scale for t, a in zip(transforms, affine) if a is None]))),
                              "an ellipse" if was == "ellipse" else "made of straight lines",
                              len(vertices)))

        data = np.column_stack([transforms[i].inverse(vertices[:, i]) for i in range(2)])
        if not np.all(np.isfinite(data)):
            raise _Skip("some of its vertices are outside the domain of the "
                        "axis transform")

        name = self.namer(gate.name, gate.alt_name)
        self.ops.append(PolygonOp(name = name,
                                  xchannel = channels[0],
                                  ychannel = channels[1],
                                  xscale = transforms[0].scale,
                                  yscale = transforms[1].scale,
                                  vertices = [(float(x), float(y)) for x, y in data]))
        return name

    def make_ellipse(self, gate):
        channels, transforms, affine = self.axes(gate, gate.dims)
        mean = np.asarray(gate.mean, dtype = float)
        cov = np.asarray(gate.cov, dtype = float)
        if mean.shape != (2,) or cov.shape != (2, 2):
            raise _Skip("it isn't 2-dimensional")
        if not (np.all(np.isfinite(mean)) and np.all(np.isfinite(cov))
                and np.allclose(cov, cov.T, rtol = 1e-6, atol = 1e-12)
                and cov[0, 0] > 0 and np.linalg.det(cov) > 0):
            raise _Skip("its covariance matrix isn't symmetric and positive definite")
        cov = (cov + cov.T) / 2.0

        if None in affine:
            t = np.linspace(0.0, 2 * np.pi, _ELLIPSE_VERTICES, endpoint = False)
            circle = np.column_stack((np.cos(t), np.sin(t)))
            # make the polygon a tiny bit bigger than the one whose vertices
            # are on the ellipse, so it has the same area as the ellipse
            circle *= np.sqrt((2 * np.pi / _ELLIPSE_VERTICES)
                              / np.sin(2 * np.pi / _ELLIPSE_VERTICES))
            vertices = mean + circle @ np.linalg.cholesky(cov).T
            return self.make_polygon(gate, vertices, gate.dims, was = "ellipse")

        # on each axis, g = alpha * s + beta, where s is cytoflow's scaled
        # space.  the center just needs the inverse transform; the shape
        # needs un-stretching.
        alpha = np.array([a[0] for a in affine])
        cov_s = cov / np.outer(alpha, alpha)
        values, vectors = np.linalg.eigh(cov_s)
        width, height, angle = _normalize_axes(vectors @ np.diag(np.sqrt(values)))

        center = [float(transforms[i].inverse(mean[i])) for i in range(2)]
        if not np.all(np.isfinite(center)):
            raise _Skip("its center is outside the domain of the axis transform")

        name = self.namer(gate.name, gate.alt_name)
        self.ops.append(EllipseOp(name = name,
                                  xchannel = channels[0],
                                  ychannel = channels[1],
                                  xscale = transforms[0].scale,
                                  yscale = transforms[1].scale,
                                  center = (center[0], center[1]),
                                  width = width,
                                  height = height,
                                  angle = angle))
        return name

    def make_bool(self, gate, refs):
        terms = [("not " + expr) if complement else expr for expr, complement in refs]
        if gate.operator == 'not':
            if len(terms) != 1:
                raise _Skip("a 'not' gate must have exactly one operand")
            expression = "not ({})".format(terms[0]) if refs[0][1] else "not " + terms[0]
        elif gate.operator in ('and', 'or'):
            if len(terms) < 2:
                raise _Skip("an '{}' gate must have at least two operands"
                            .format(gate.operator))
            expression = " {} ".format(gate.operator).join(terms)
        else:
            raise _Skip("unknown boolean operator '{}'".format(gate.operator))

        name = self.namer(gate.name, gate.alt_name)
        self.ops.append(BooleanGateOp(name = name, expression = expression))
        return name


###############################################################################
# import Gating-ML
###############################################################################

def _gml_transforms(root):
    """Returns a dict of transformation id -> `_Transform`, or a `_Skip` saying why not"""
    transforms = {}
    for el in _children(root, 'transformation'):
        tid = _attr(el, 'id')
        if tid is None or len(el) == 0:
            continue
        t = el[0]
        kind = _local(t.tag)
        p = lambda name: _float(_attr(t, name), "transform parameter " + name)
        try:
            if kind == 'flin':
                transforms[tid] = _LinTransform(p('T'), p('A'))
            elif kind == 'flog':
                transforms[tid] = _LogTransform(p('T'), p('M'))
            elif kind == 'fasinh':
                transforms[tid] = _AsinhTransform(p('T'), p('M'), p('A'))
            elif kind == 'logicle':
                transforms[tid] = _LogicleTransform(p('T'), p('W'), p('M'), p('A'))
            elif kind == 'hyperlog':
                transforms[tid] = _HyperlogTransform(p('T'), p('W'), p('M'), p('A'))
            else:
                raise _Skip("'{}' transforms ('{}') aren't supported".format(kind, tid))
        except _Skip as e:
            transforms[tid] = e
    return transforms


def _gml_dimension(el, transforms, notes):
    """Parse a <dimension> or <divider>.  Bounds stay in the transformed space."""
    fcs = _child(el, 'fcs-dimension')
    if fcs is None:
        if _child(el, 'new-dimension') is not None:
            raise _Skip("it uses a ratio of two parameters, which isn't supported")
        raise _Skip("one of its dimensions has no parameter name")
    name = _attr(fcs, 'name')
    if not name:
        raise _Skip("one of its dimensions has no parameter name")

    comp = _attr(el, 'compensation-ref')
    if comp not in (None, '', 'uncompensated'):
        notes.add(comp)

    tid = _attr(el, 'transformation-ref')
    if tid in (None, ''):
        transform = _IdentityTransform()
    elif tid not in transforms:
        raise _Skip("transformation '{}' isn't in the file".format(tid))
    elif isinstance(transforms[tid], _Skip):
        raise _Skip(str(transforms[tid]))
    else:
        transform = transforms[tid]

    return name, transform


def _gml_bound(el, attr, transform):
    text = _attr(el, attr)
    if text is None:
        return None
    value = float(transform.inverse(_float(text, "the gate's " + attr)))
    if not np.isfinite(value):
        raise _Skip("its {} is outside the domain of the axis transform".format(attr))
    return value


def _gml_coordinates(el, tag = 'coordinate'):
    return [_float(_attr(c, 'value'), "a coordinate") for c in _children(el, tag)]


def _gml_gate(el, transforms, comp_notes):
    """Parse one gate element, returning a list of `_Gate`"""
    kind = _local(el.tag)
    gate_id = _attr(el, 'id')
    gate = _Gate(gate_id, gate_id, parent = _attr(el, 'parent_id') or None)

    try:
        if kind == 'RectangleGate':
            gate.kind = 'rect'
            gate.dims = []
            for d in _children(el, 'dimension'):
                name, transform = _gml_dimension(d, transforms, comp_notes)
                gate.dims.append(_Dim(name,
                                      lo = _gml_bound(d, 'min', transform),
                                      hi = _gml_bound(d, 'max', transform)))

        elif kind == 'PolygonGate':
            gate.kind = 'poly'
            gate.dims = [_Dim(*_gml_dimension(d, transforms, comp_notes))
                         for d in _children(el, 'dimension')]
            gate.vertices = [_gml_coordinates(v) for v in _children(el, 'vertex')]
            if any(len(v) != 2 for v in gate.vertices):
                raise _Skip("its vertices aren't 2-dimensional")

        elif kind == 'EllipsoidGate':
            gate.kind = 'ellipse'
            dims = _children(el, 'dimension')
            if len(dims) != 2:
                raise _Skip("it has {} dimensions; only 2-dimensional ellipsoid "
                            "gates are supported".format(len(dims)))
            gate.dims = [_Dim(*_gml_dimension(d, transforms, comp_notes)) for d in dims]
            mean, cov, dist = (_child(el, 'mean'), _child(el, 'covarianceMatrix'),
                               _child(el, 'distanceSquare'))
            if mean is None or cov is None or dist is None:
                raise _Skip("it's missing its mean, covariance matrix or distance")
            gate.mean = _gml_coordinates(mean)
            d2 = _float(_attr(dist, 'value'), "distanceSquare")
            if not d2 > 0:
                raise _Skip("its distanceSquare isn't positive")
            # (x - mu)' inv(C) (x - mu) <= D^2 is the same as using a
            # covariance of C * D^2 and a distance of 1
            gate.cov = np.array([_gml_coordinates(r, 'entry')
                                 for r in _children(cov, 'row')]) * d2

        elif kind == 'BooleanGate':
            gate.kind = 'bool'
            operators = [c for c in el if _local(c.tag) in ('and', 'or', 'not')]
            if len(operators) != 1:
                raise _Skip("it must have exactly one 'and', 'or' or 'not'")
            gate.operator = _local(operators[0].tag)
            gate.refs = []
            for ref in _children(operators[0], 'gateReference'):
                complement = (_attr(ref, 'use-as-complement', 'false').strip().lower()
                              in ('true', '1'))
                gate.refs.append((_attr(ref, 'ref'), complement))

        elif kind == 'QuadrantGate':
            return _gml_quadrant_gate(el, gate, transforms, comp_notes)

    except _Skip as e:
        gate.error = str(e)

    return [gate]


def _gml_quadrant_gate(el, gate, transforms, comp_notes):
    """
    A QuadrantGate with two dividers of one value each, and four quadrants
    named the way `QuadOp` names them, becomes a `QuadOp`.  Any other
    becomes a rectangle gate for each quadrant.
    """
    quadrants = _children(el, 'Quadrant')
    gates = [_Gate(_attr(q, 'id'), _attr(q, 'id'), parent = gate.parent)
             for q in quadrants]
    try:
        dividers = {}
        order = []
        for d in _children(el, 'divider'):
            name, transform = _gml_dimension(d, transforms, comp_notes)
            values = []
            for v in _children(d, 'value'):
                value = float(transform.inverse(_float(v.text, "a divider value")))
                if not np.isfinite(value):
                    raise _Skip("a divider is outside the domain of the axis transform")
                values.append(value)
            if not values:
                raise _Skip("a divider has no values")
            dividers[_attr(d, 'id')] = (name, transform, sorted(values))
            order.append(_attr(d, 'id'))

        # where is each quadrant?  a (lo, hi) in data units for each divider
        bounds = []
        for q in quadrants:
            b = {}
            for pos in _children(q, 'position'):
                ref = _attr(pos, 'divider_ref')
                if ref not in dividers:
                    raise _Skip("a quadrant refers to a divider that doesn't exist")
                _, transform, values = dividers[ref]
                location = float(transform.inverse(_float(_attr(pos, 'location'),
                                                          "a quadrant's location")))
                idx = int(np.searchsorted(values, location, side = 'right'))
                b[ref] = (values[idx - 1] if idx > 0 else None,
                          values[idx] if idx < len(values) else None)
            bounds.append(b)

        # is it just like a QuadOp?
        if (len(order) == 2 and len(quadrants) == 4
                and all(len(dividers[d][2]) == 1 for d in order)
                and dividers[order[0]][0] != dividers[order[1]][0]
                and all(len(b) == 2 for b in bounds)):
            expected = [(False, True), (True, True), (False, False), (True, False)]
            actual = {q.id : tuple(b[d][0] is not None for d in order)
                      for q, b in zip(gates, bounds)}
            if all(actual.get("{}_{}".format(gate.id, i + 1)) == expected[i]
                   for i in range(4)):
                gate.kind = 'quad'
                gate.dims = [_Dim(dividers[d][0]) for d in order]
                gate.thresholds = [dividers[d][2][0] for d in order]
                # the QuadOp applies to all events: the quadrants get the parent.
                gate.parent = None
                for q in gates:
                    q.kind = 'quadrant'
                    q.owner = gate.id
                    q.value = q.id

        for q, b in zip(gates, bounds):
            q.kind = q.kind or 'rect'
            q.dims = [_Dim(dividers[d][0], lo = b[d][0], hi = b[d][1])
                      for d in order if d in b]
        return ([gate] if gate.kind == 'quad' else []) + gates

    except _Skip as e:
        for q in gates:
            q.error = str(e)
        if not gates:
            gate.error = str(e)
            return [gate]
        return gates


def import_gatingml(filename, experiment):
    """
    Read the gates from a Gating-ML 2.0 file and return the equivalent
    `cytoflow` operations.

    Parameters
    ----------
    filename : str or path
        The Gating-ML file to read.

    experiment : `Experiment`
        The experiment that the gates will be applied to.  It's used to match
        the parameter names in the file to channels (using the channels'
        names in the FCS files, if they've been renamed); to work out
        how the file's axis transforms relate to `cytoflow`'s scales for
        those channels; and to make sure the new gates' names don't collide
        with existing conditions.

    Returns
    -------
    list of `IOperation`
        The operations, in an order in which they can be applied -- ie, a
        `BooleanGateOp` always comes after the gates it refers to.  Apply
        them all to get all of the gates::

            for op in flow.import_gatingml("gates.xml", ex):
                ex = op.apply(ex)

    Notes
    -----
    What becomes what:

    - A ``RectangleGate`` becomes a `ThresholdOp` (one dimension, only a
      minimum), a `RangeOp` (one dimension) or a `Range2DOp` (two).  A missing
      bound is infinite.

    - A ``PolygonGate`` becomes a `PolygonOp`.

    - A two-dimensional ``EllipsoidGate`` becomes an `EllipseOp`.

    - A ``QuadrantGate`` becomes a `QuadOp` if it is exactly like one (as
      written by `export_gatingml`); otherwise, each of its quadrants becomes
      a `RangeOp` or `Range2DOp`.

    - A ``BooleanGate`` becomes a `BooleanGateOp`.

    **Gate hierarchy.**  `cytoflow` gates aren't nested: every gate applies to
    all of the events.  So a gate ``X`` with a ``parent_id`` of ``P`` becomes
    *two* operations: the gate itself, named ``X``, which ignores the parent;
    and a `BooleanGateOp` named ``X_gated`` with the expression ``X and P``
    (or ``X and P_gated``, if ``P`` has a parent too.)  ``X_gated`` is the
    population that other programs call ``X``.  A ``BooleanGate`` that refers
    to ``X`` uses ``X_gated``.

    **Names.**  Gate IDs that aren't legal `cytoflow` names are changed
    (``CD4+`` becomes ``CD4pos``), as are names that are already used in
    ``experiment``.

    **Axis transforms.**  Polygon and ellipse gates get the `cytoflow` scale
    that corresponds to the dimension's transform: ``flin`` (or none) is
    ``linear``, ``flog`` is ``log``, ``logicle`` and ``hyperlog`` are
    ``logicle`` and ``fasinh`` is ``arcsinh``.  If the transform's parameters
    aren't the ones `cytoflow` uses for that channel, the gate's outline
    is curved in `cytoflow`'s scale: polygons get extra vertices to follow it,
    and ellipses become polygons.  (To avoid that, set the channel's scale
    parameters with `ScaleSettingsOp` *before* importing.)  A warning says
    which gates were affected.

    **Not imported** (there is a `CytoflowWarning` that lists the gates
    that were skipped, and why): gates with more than two dimensions; gates
    on ratios of parameters (``fratio``); gates on parameters that aren't in
    ``experiment``; and anything that depends on a gate that was skipped.

    **Compensation** isn't imported: the gates are applied to the
    experiment's channels as they are.  There's a warning if the file's gates
    expect compensated data.
    """

    if experiment is None:
        raise util.CytoflowError("No experiment specified")

    root = _parse(filename)
    if _local(root.tag) != 'Gating-ML':
        if _local(root.tag) == 'Workspace':
            raise util.CytoflowError("{} is a FlowJo workspace, not a Gating-ML "
                                     "file. Use import_flowjo_workspace()."
                                     .format(filename))
        raise util.CytoflowError("{} isn't a Gating-ML file".format(filename))

    transforms = _gml_transforms(root)
    comp_notes = set()
    gates = {}
    duplicates = []
    for el in root:
        if _local(el.tag) in ('RectangleGate', 'PolygonGate', 'EllipsoidGate',
                              'BooleanGate', 'QuadrantGate'):
            for gate in _gml_gate(el, transforms, comp_notes):
                if not gate.id:
                    duplicates.append("(a {}): it has no ID".format(_local(el.tag)))
                elif gate.id in gates:
                    duplicates.append("{}: there is more than one gate with this ID"
                                      .format(gate.id))
                else:
                    gates[gate.id] = gate

    builder = _Builder(gates, experiment)
    builder.skipped.extend(duplicates)
    if comp_notes:
        builder.note("The gates in this file expect compensated data "
                     "(compensation: {}). Compensation is NOT imported: make "
                     "sure the experiment is compensated the same way."
                     .format(", ".join(sorted(comp_notes))))
    return builder.build()


###############################################################################
# import FlowJo workspaces
###############################################################################

def _wsp_transforms(sample_el):
    """Returns a dict of parameter name -> `_Transform`, or a `_Skip`"""
    transforms = {}
    container = _child(sample_el, 'Transformations')
    if container is None:
        return transforms

    for t in container:
        kind = _local(t.tag)
        param = _child(t, 'parameter')
        name = _attr(param, 'name') if param is not None else None
        if not name:
            continue
        p = lambda attr: _float(_attr(t, attr), "transform parameter " + attr)
        try:
            if kind == 'linear':
                # all that matters for polygons is that it's linear.
                # ellipses are stored on a 0 - 256 grid, and for those we
                # need to know where the grid is.
                lo, hi = p('minRange'), p('maxRange')
                transform = _LinTransform(hi, 0.0)
                if lo != 0.0:
                    transform.no_grid = ("its linear axis starts at {}, not 0"
                                         .format(lo))
            elif kind == 'log':
                offset, decades = p('offset'), p('decades')
                if not (offset > 0 and decades > 0):
                    raise _Skip("bad log transform parameters")
                transform = _LogTransform(offset * 10 ** decades, decades, floor = offset)
            elif kind == 'logicle':
                transform = _LogicleTransform(p('T'), p('W'), p('M'), p('A'))
            elif kind == 'fasinh':
                transform = _AsinhTransform(p('T'), p('M'), p('A'))
            elif kind == 'biex':
                # FlowJo's biexponential is computed from a lookup table
                # whose construction isn't documented.  Its parameters map
                # onto the logicle's: the two agree to within about 0.5% of
                # the width of the axis.
                width = p('width')
                if not width < 0:
                    raise _Skip("biex transform with a non-negative width basis")
                transform = _LogicleTransform(T = p('maxRange'),
                                              W = np.log10(-width) / 2.0,
                                              M = p('pos'),
                                              A = p('neg'))
                transform.approximate = (
                    "FlowJo's biexponential axis was approximated with a "
                    "logicle (T = {:g}, W = {:g}, M = {:g}, A = {:g}). The "
                    "two differ by up to about 0.5% of the width of the axis, "
                    "so the edges of polygon and ellipse gates on this channel "
                    "may be very slightly different than in FlowJo."
                    .format(transform.T, transform.W, transform.M, transform.A))
                if _attr(t, 'length') not in (None, '256'):
                    transform.no_grid = ("its biexponential axis has a "
                                         "'length' of {}, not 256"
                                         .format(_attr(t, 'length')))
            else:
                raise _Skip("FlowJo's '{}' axis transform isn't supported"
                            .format(kind))
        except _Skip as e:
            transform = e
        transforms[name] = transform
    return transforms


class _WspSample(object):
    """The things we need to interpret one sample's gates"""

    def __init__(self, sample_el, experiment):
        self.experiment = experiment
        self.transforms = _wsp_transforms(sample_el)
        self.comp = None
        self.comp_used = set()

        matrix = _child(sample_el, 'spilloverMatrix')
        if matrix is not None:
            params_el = _child(matrix, 'parameters')
            params = [_attr(p, 'name') for p in _children(params_el, 'parameter')] \
                     if params_el is not None else []
            self.comp = (matrix.get('prefix', ''), matrix.get('suffix', ''), set(params),
                         matrix.get('name', ''))

    def dimension(self, dim_el):
        """Returns the name of the parameter in the FCS file, and its transform"""
        fcs = _child(dim_el, 'fcs-dimension')
        name = _attr(fcs, 'name') if fcs is not None else None
        if not name:
            raise _Skip("one of its dimensions has no parameter name")

        # compensated parameters are called (prefix)(name)(suffix)
        fcs_name = name
        if self.comp is not None:
            prefix, suffix, params, _ = self.comp
            if (prefix or suffix) and name.startswith(prefix) and name.endswith(suffix):
                stripped = name[len(prefix) : len(name) - len(suffix)]
                if stripped in params:
                    fcs_name = stripped
                    self.comp_used.add(name)

        transform = self.transforms.get(name, self.transforms.get(fcs_name, None))
        return fcs_name, transform

    def gate(self, pop_el, gate):
        """Fill in `gate` from a <Population> element"""
        gate_el = _child(pop_el, 'Gate')
        if gate_el is None or len(gate_el) == 0:
            raise _Skip("it has no gate")
        shape = gate_el[0]
        kind = _local(shape.tag)
        gate.complement = (shape.get('eventsInside', '1').strip() == '0')

        if kind not in ('RectangleGate', 'PolygonGate', 'EllipsoidGate'):
            raise _Skip("FlowJo '{}' gates aren't supported".format(kind))

        dims = [self.dimension(d) for d in _children(shape, 'dimension')]

        if kind == 'RectangleGate':
            # the bounds are in data units, whatever the axis transform
            gate.kind = 'rect'
            gate.dims = []
            for d, (name, _) in zip(_children(shape, 'dimension'), dims):
                lo, hi = _attr(d, 'min'), _attr(d, 'max')
                gate.dims.append(_Dim(name,
                                      lo = _float(lo, "the gate's min") if lo is not None else None,
                                      hi = _float(hi, "the gate's max") if hi is not None else None))
            return

        if len(dims) != 2:
            raise _Skip("it has {} dimensions, not 2".format(len(dims)))
        for name, transform in dims:
            if isinstance(transform, _Skip):
                raise _Skip("{} (channel '{}')".format(transform, name))
        gate.dims = [_Dim(name, transform or _IdentityTransform()) for name, transform in dims]

        if kind == 'PolygonGate':
            # the vertices are in data units; the edges are straight lines
            # in the transformed space
            gate.kind = 'poly'
            vertices = np.array([_gml_coordinates(v) for v in _children(shape, 'vertex')],
                                dtype = float)
            if vertices.ndim != 2 or vertices.shape[1] != 2:
                raise _Skip("its vertices aren't 2-dimensional")
            gate.vertices = np.column_stack([gate.dims[i].transform.forward(vertices[:, i])
                                             for i in range(2)])
            if not np.all(np.isfinite(gate.vertices)):
                raise _Skip("some of its vertices are outside the domain of "
                            "the axis transform")

        else:
            # ellipses are different: two foci and four points on the edge,
            # on a 256 x 256 grid that covers the transformed axes.
            gate.kind = 'ellipse'
            for (name, transform), d in zip(dims, gate.dims):
                if transform is None:
                    raise _Skip("FlowJo stores ellipses in display coordinates, "
                                "and the workspace has no axis transform for "
                                "channel '{}'".format(name))
                if getattr(transform, 'no_grid', None):
                    raise _Skip("FlowJo stores ellipses in display coordinates, "
                                "which can't be interpreted because {} (channel '{}')"
                                .format(transform.no_grid, name))

            foci_el, edge_el = _child(shape, 'foci'), _child(shape, 'edge')
            if foci_el is None or edge_el is None:
                raise _Skip("it doesn't have foci and edge points")
            foci = np.array([_gml_coordinates(v) for v in _children(foci_el, 'vertex')],
                            dtype = float) / 256.0
            edge = np.array([_gml_coordinates(v) for v in _children(edge_el, 'vertex')],
                            dtype = float) / 256.0
            if foci.shape != (2, 2) or edge.ndim != 2 or edge.shape[1] != 2 or len(edge) < 1:
                raise _Skip("it doesn't have two foci and (four) edge points")

            # every point on an ellipse is the same total distance, 2a, from
            # the two foci.  the foci are a distance c from the center, and
            # b^2 = a^2 - c^2
            a = np.mean(np.linalg.norm(edge - foci[0], axis = 1)
                        + np.linalg.norm(edge - foci[1], axis = 1)) / 2.0
            c = np.linalg.norm(foci[1] - foci[0]) / 2.0
            if not a > c:
                raise _Skip("its foci and edge points don't describe an ellipse")
            b = np.sqrt(a ** 2 - c ** 2)
            major = (foci[1] - foci[0]) / (2 * c) if c > 1e-12 else np.array([1.0, 0.0])
            minor = np.array([-major[1], major[0]])

            gate.mean = (foci[0] + foci[1]) / 2.0
            gate.cov = a ** 2 * np.outer(major, major) + b ** 2 * np.outer(minor, minor)


def _wsp_walk(subpops_el, parent, path, sample, gates):
    for el in subpops_el:
        tag = _local(el.tag)
        if tag not in ('Population', 'AndNode', 'OrNode', 'NotNode'):
            continue

        name = el.get('name', '') or tag
        gate_id = "/".join(path + [name])
        gate = _Gate(gate_id, name, parent = parent)
        if path:
            gate.alt_name = "{}_{}".format(path[-1], name)

        try:
            if gate_id in gates:
                raise _Skip("there is more than one population with this name "
                            "under the same parent")
            if tag == 'Population':
                sample.gate(el, gate)
            else:
                # boolean gates refer to other populations by their paths
                gate.kind = 'bool'
                gate.operator = tag[:-4].lower()
                deps = _child(el, 'Dependents')
                gate.refs = [(d.get('name', ''), False)
                             for d in (_children(deps, 'Dependent') if deps is not None else [])]
        except _Skip as e:
            gate.error = str(e)

        if gate_id not in gates:
            gates[gate_id] = gate

        for sub in _children(el, 'Subpopulations'):
            _wsp_walk(sub, gate_id, path + [name], sample, gates)


def import_flowjo_workspace(filename, experiment, sample = None, group = None):
    """
    Read the gates of one sample (or one group) from a FlowJo 10 workspace
    (``.wsp``) and return the equivalent `cytoflow` operations.

    FlowJo keeps a separate gating tree for every sample (and a template for
    every group), but `cytoflow` applies one set of gates to the whole
    experiment.  So, you choose which tree to import.

    Parameters
    ----------
    filename : str or path
        The workspace to read.

    experiment : `Experiment`
        The experiment that the gates will be applied to -- see
        `import_gatingml`.

    sample : str
        The name of the sample whose gates to import (usually the name of
        its FCS file.)

    group : str
        The name of the group whose gates to import.  The group's gates
        are interpreted with the axis transforms of the first sample in
        the group.  Specify either ``sample`` or ``group``.  If the workspace
        only has one sample, you don't need to specify either one; otherwise,
        the error message lists the samples and groups in the workspace.

    Returns
    -------
    list of `IOperation`
        The operations, in an order in which they can be applied.

    Notes
    -----
    What becomes what:

    - Rectangle gates become a `ThresholdOp`, `RangeOp` or `Range2DOp`.
      FlowJo quadrant gates are four rectangle gates, and become four
      `Range2DOp` s.

    - Polygon gates become a `PolygonOp`.

    - Ellipse gates become an `EllipseOp`.

    - Boolean (and / or / not) gates become a `BooleanGateOp`.

    - A gate that selects the events *outside* its outline becomes two
      operations: the outline (``<name>_region``) and a `BooleanGateOp`.

    The **hierarchy** is handled the same way as in `import_gatingml`:
    the population ``Singlets`` under ``Cells`` becomes a gate ``Singlets``
    which applies to *all* of the events, and a `BooleanGateOp`
    ``Singlets_gated`` with the expression ``Singlets and Cells``.  The
    ``_gated`` conditions are the ones that match FlowJo's populations.

    **Coordinates and axis transforms.**  FlowJo saves rectangles and the
    vertices of polygons in data units, so they are imported as-is; the edges
    of a polygon are straight lines on FlowJo's (transformed) axes, and the
    `PolygonOp` gets the corresponding scale: ``linear``, ``log``,
    ``logicle`` (for FlowJo's logicle *and* biexponential axes) or ``arcsinh``.
    FlowJo saves ellipses on a 256 x 256 grid that covers the transformed
    axes.  If FlowJo's axis isn't a shifted and stretched version of
    `cytoflow`'s scale for that channel (as is usually the case for logicle
    and biexponential axes, whose parameters won't match), polygons get extra
    vertices to follow the original edges, and ellipses become polygons.

    Limitations
    -----------
    - FlowJo's biexponential transform is **approximated** by a logicle
      with the corresponding parameters; they differ by up to about 0.5% of
      the axis.  Rectangle gates are not affected.  There's a warning.

    - Polygon and ellipse gates on other FlowJo axis transforms (hyperlog,
      Miltenyi, ...), curly quadrants, and ellipses on axes whose display
      range can't be determined, are skipped, with a warning.

    - **Compensation is not imported.**  Gates on compensated parameters
      (``Comp-FITC-A``) are applied to the experiment's channel of the same
      name (``FITC-A``): compensate the experiment in `cytoflow` first.
      There's a warning.

    - Derived parameters, and per-sample differences in the gates, aren't
      imported.

    - This was written from the published descriptions of the workspace
      format; check the imported gates against FlowJo's event counts.
    """

    if experiment is None:
        raise util.CytoflowError("No experiment specified")

    root = _parse(filename)
    if _local(root.tag) != 'Workspace':
        raise util.CytoflowError("{} isn't a FlowJo workspace (.wsp)".format(filename))

    # the samples
    samples = []     # (name, id, Sample element, SampleNode element)
    sample_list = _child(root, 'SampleList')
    for sample_el in (_children(sample_list, 'Sample') if sample_list is not None else []):
        node = _child(sample_el, 'SampleNode')
        if node is not None:
            samples.append((node.get('name', ''), node.get('sampleID', ''), sample_el, node))

    # and the groups
    groups = {}      # name -> (GroupNode element, [sample ids])
    groups_el = _child(root, 'Groups')
    for group_node in (_children(groups_el, 'GroupNode') if groups_el is not None else []):
        ids = []
        group_el = _child(group_node, 'Group')
        refs = _child(group_el, 'SampleRefs') if group_el is not None else None
        if refs is not None:
            ids = [r.get('sampleID', '') for r in _children(refs, 'SampleRef')]
        groups[group_node.get('name', '')] = (group_node, ids)

    available = ("Samples: {}. Groups: {}."
                 .format(", ".join([repr(s[0]) for s in samples]) or "(none)",
                         ", ".join([repr(g) for g in groups]) or "(none)"))

    if sample is not None and group is not None:
        raise util.CytoflowError("Specify either 'sample' or 'group', not both")

    if group is not None:
        if group not in groups:
            raise util.CytoflowError("There's no group {!r} in {}. {}"
                                     .format(group, filename, available))
        group_node, ids = groups[group]
        members = [s for s in samples if s[1] in ids]
        if not members:
            raise util.CytoflowError("Group {!r} has no samples, so there are "
                                     "no axis transforms to interpret its "
                                     "gates with. {}".format(group, available))
        sample_el = members[0][2]
        tree = group_node
    else:
        if sample is None:
            if len(samples) != 1:
                raise util.CytoflowError("{} has {} samples; specify which "
                                         "'sample' or 'group' to import the "
                                         "gates from. {}"
                                         .format(filename, len(samples), available))
            match = samples
        else:
            match = [s for s in samples if s[0] == sample] or \
                    [s for s in samples if s[1] == str(sample)]
        if not match:
            raise util.CytoflowError("There's no sample {!r} in {}. {}"
                                     .format(sample, filename, available))
        if len(match) > 1:
            raise util.CytoflowError("There are {} samples named {!r} in {}; "
                                     "use the sample ID instead ({})"
                                     .format(len(match), sample, filename,
                                             ", ".join([s[1] for s in match])))
        sample_el = match[0][2]
        tree = match[0][3]

    wsp_sample = _WspSample(sample_el, experiment)
    gates = {}
    for sub in _children(tree, 'Subpopulations'):
        _wsp_walk(sub, None, [], wsp_sample, gates)

    builder = _Builder(gates, experiment)
    if wsp_sample.comp_used:
        builder.note("These gates were drawn on compensated parameters ({}; "
                     "FlowJo compensation matrix '{}'). Compensation is NOT "
                     "imported: the gates are applied to the experiment's "
                     "channels as they are, so compensate the experiment first."
                     .format(", ".join(sorted(wsp_sample.comp_used)),
                             wsp_sample.comp[3]))
    ops = builder.build()

    if not gates:
        warnings.warn("Didn't find any gates for {} in {}"
                      .format("group {!r}".format(group) if group is not None
                              else "the sample", filename),
                      util.CytoflowWarning)
    return ops


def flowjo_workspace_samples(filename):
    """
    List the samples and groups in a FlowJo workspace.

    Returns
    -------
    (list of str, list of str)
        The names of the samples, and the names of the groups.
    """
    root = _parse(filename)
    if _local(root.tag) != 'Workspace':
        raise util.CytoflowError("{} isn't a FlowJo workspace (.wsp)".format(filename))

    sample_list = _child(root, 'SampleList')
    samples = []
    for sample_el in (_children(sample_list, 'Sample') if sample_list is not None else []):
        node = _child(sample_el, 'SampleNode')
        if node is not None:
            samples.append(node.get('name', ''))

    groups_el = _child(root, 'Groups')
    groups = [g.get('name', '')
              for g in (_children(groups_el, 'GroupNode') if groups_el is not None else [])]
    return samples, groups
