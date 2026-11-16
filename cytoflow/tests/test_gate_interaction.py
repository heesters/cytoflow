#!/usr/bin/env python3.8
# coding: latin-1

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
Scripted mouse / keyboard interaction with the ellipse and freehand gate
views, using synthesized matplotlib events.
"""

import unittest
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent, KeyEvent

import cytoflow as flow
from .test_base import ImportedDataSmallTest


def _mouse(canvas, kind, ax, x, y, button = 1):
    X, Y = ax.transData.transform((x, y))
    MouseEvent(kind, canvas, X, Y, button)._process()

def _click(canvas, ax, x, y):
    _mouse(canvas, "button_press_event", ax, x, y)
    _mouse(canvas, "button_release_event", ax, x, y)

def _key(canvas, key):
    KeyEvent("key_press_event", canvas, key, 0, 0)._process()
    KeyEvent("key_release_event", canvas, key, 0, 0)._process()

def _drag(canvas, ax, a, b, n = 8):
    _mouse(canvas, "button_press_event", ax, *a)
    for t in np.linspace(0, 1, n):
        _mouse(canvas, "motion_notify_event", ax,
               a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)
    _mouse(canvas, "button_release_event", ax, *b)


class TestEllipseInteraction(ImportedDataSmallTest):

    def _check_consistent(self, op, view, expect_events = True):
        """The events inside the ellipse on the screen are the gated events."""
        ax = view._ax
        pts = ax.transData.transform(self.ex.data[[op.xchannel, op.ychannel]].values)
        path = view._patch.get_transform().transform_path(view._patch.get_path())
        on_screen = path.contains_points(pts)
        gated = op.apply(self.ex).data[op.name].values

        # a few events right on the boundary may differ, because the
        # patch's outline is a bezier approximation
        if expect_events:
            self.assertGreater(gated.sum(), 0)
        self.assertLessEqual((on_screen != gated).sum(), max(5, 0.005 * gated.sum()))

    def _run(self, scale):
        # with linear scales, the events are all in the bottom-left corner
        events = scale != "linear"
        
        plt.close('all')
        op = flow.EllipseOp(name = "E", xchannel = "V2-A", ychannel = "Y2-A",
                            xscale = scale, yscale = scale)
        v = op.default_view()
        v.plot(self.ex)
        v.interactive = True

        ov = v._overlay  # the selector's axes, in scaled coordinates
        canvas = ov.figure.canvas
        canvas.draw()
        x0, x1 = ov.get_xlim()
        y0, y1 = ov.get_ylim()
        xs, ys = x1 - x0, y1 - y0

        # draw it
        a = (x0 + 0.3 * xs, y0 + 0.55 * ys)
        b = (x0 + 0.7 * xs, y0 + 0.85 * ys)
        _drag(canvas, ov, a, b)
        self.assertAlmostEqual(op.width, 0.4 * xs, delta = 0.02 * xs)
        self.assertAlmostEqual(op.height, 0.3 * ys, delta = 0.02 * ys)
        self.assertEqual(op.angle, 0.0)
        self.assertTrue(np.allclose(v._scaled_center(),
                                    ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2),
                                    atol = 0.01 * max(xs, ys)))
        self._check_consistent(op, v, events)

        # a stray click outside the ellipse leaves it alone
        saved = (op.center, op.width, op.height, op.angle)
        _click(canvas, ov, x0 + 0.05 * xs, y0 + 0.05 * ys)
        self.assertEqual((op.center, op.width, op.height, op.angle), saved)
        self.assertTrue(v._widget._selection_completed)

        # drag the right-hand edge to make it wider
        ex0, ex1, ey0, ey1 = v._widget.extents
        _drag(canvas, ov, (ex1, (ey0 + ey1) / 2), (ex1 + 0.1 * xs, (ey0 + ey1) / 2))
        self.assertAlmostEqual(op.width, 0.5 * xs, delta = 0.02 * xs)

        # drag the inside to move it
        cx, cy = v._widget.center
        _drag(canvas, ov, (cx + 0.01 * xs, cy), (cx - 0.04 * xs, cy - 0.1 * ys))
        self.assertTrue(np.allclose(v._scaled_center(),
                                    (cx - 0.05 * xs, cy - 0.1 * ys),
                                    atol = 0.01 * max(xs, ys)))

        # press 'r', then drag a corner 20 degrees (on the screen) to rotate
        _key(canvas, 'r')
        corner = np.array(v._widget.corners).T[2]
        pc = ov.transData.transform(v._widget.center)
        pk = ov.transData.transform(corner)
        t = np.deg2rad(20)
        d = pk - pc
        pk2 = pc + (d[0] * np.cos(t) - d[1] * np.sin(t),
                    d[0] * np.sin(t) + d[1] * np.cos(t))
        _drag(canvas, ov, tuple(corner),
              tuple(ov.transData.inverted().transform(pk2)), n = 12)
        _key(canvas, 'r')
        self.assertAlmostEqual(v._widget.rotation, 20, delta = 1.5)
        self.assertNotEqual(op.angle, 0.0)
        self._check_consistent(op, v, events)

        # the selector's ellipse is the gate's ellipse
        sel = v._widget._selection_artist
        gate = v._patch.get_transform().transform_path(v._patch.get_path())
        t = np.linspace(0, 2 * np.pi, 60)
        unit = np.column_stack((np.cos(t), np.sin(t)))
        edge = sel.get_transform().transform(unit)
        d = edge - ov.transData.transform(v._widget.center)
        d /= np.hypot(d[:, 0], d[:, 1])[:, None]
        self.assertTrue(gate.contains_points(edge - 1.5 * d).all())
        self.assertFalse(gate.contains_points(edge + 1.5 * d).any())

        # re-plotting loads the gate back into the selector, unchanged
        rotation, extents = v._widget.rotation, v._widget.extents
        saved = (op.center, op.width, op.height, op.angle)
        plt.close('all')
        v2 = op.default_view()
        v2.plot(self.ex)
        v2.interactive = True
        v2._overlay.figure.canvas.draw()
        self.assertAlmostEqual(v2._widget.rotation, rotation, delta = 0.05)
        self.assertTrue(np.allclose(v2._widget.extents, extents, rtol = 1e-3))

        # ... and a stray click doesn't erase it (or reset the rotation)
        _click(v2._overlay.figure.canvas, v2._overlay,
               x0 + 0.02 * xs, y0 + 0.97 * ys)
        self.assertEqual((op.center, op.width, op.height, op.angle), saved)
        self.assertAlmostEqual(v2._widget.rotation, rotation, delta = 0.05)
        self._check_consistent(op, v2, events)

    def testLogicle(self):
        self._run("logicle")

    def testLog(self):
        self._run("log")

    def testLinear(self):
        self._run("linear")


class TestFreehandInteraction(ImportedDataSmallTest):

    def _run(self, scale, density):
        plt.close('all')
        op = flow.FreehandOp(name = "F", xchannel = "V2-A", ychannel = "Y2-A",
                             xscale = scale, yscale = scale)
        v = op.default_view(density = density)
        v.plot(self.ex)
        v.interactive = True
        ax = v._ax
        canvas = ax.figure.canvas
        canvas.draw()

        # draw a densely-sampled loop around the upper population
        bb = ax.bbox
        t = np.linspace(0, 2 * np.pi, 400)
        pix = np.column_stack((bb.x0 + bb.width * (0.5 + 0.3 * np.cos(t)),
                               bb.y0 + bb.height * (0.65 + 0.25 * np.sin(t))))
        data = ax.transData.inverted().transform(pix)
        _mouse(canvas, "button_press_event", ax, *data[0])
        for p in data[1:]:
            _mouse(canvas, "motion_notify_event", ax, *p)
        _mouse(canvas, "button_release_event", ax, *data[-1])

        self.assertGreaterEqual(len(op.vertices), 10)
        self.assertLessEqual(len(op.vertices), 100)

        # the gate matches what's on the screen
        pts = ax.transData.transform(self.ex.data[["V2-A", "Y2-A"]].values)
        path = v._patch.get_transform().transform_path(v._patch.get_path())
        gated = op.apply(self.ex).data["F"].values
        self.assertGreater(gated.sum(), 0)
        self.assertLessEqual((path.contains_points(pts) != gated).sum(), 2)

        # a stray click, or a tiny scribble, doesn't erase the gate
        saved = list(op.vertices)
        _click(canvas, ax, *data[0])
        self.assertEqual(op.vertices, saved)
        _mouse(canvas, "button_press_event", ax, *data[0])
        p = ax.transData.inverted().transform(pix[0] + (2, 1))
        _mouse(canvas, "motion_notify_event", ax, *p)
        _mouse(canvas, "button_release_event", ax, *p)
        self.assertEqual(op.vertices, saved)

        # drawing again replaces the gate
        _mouse(canvas, "button_press_event", ax, *data[0])
        for p in data[1:200]:
            _mouse(canvas, "motion_notify_event", ax, *p)
        _mouse(canvas, "button_release_event", ax, *data[199])
        self.assertNotEqual(op.vertices, saved)

    def testLogicle(self):
        self._run("logicle", False)

    def testLogDensity(self):
        self._run("log", True)


if __name__ == "__main__":
    unittest.main()
