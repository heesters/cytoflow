#!/usr/bin/env python3
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
cytoflow.operations._selectors
------------------------------

Thin subclasses of the `matplotlib.widgets` selectors used by the interactive
gate views.  They add a few ease-of-use behaviours that matter most on a
laptop trackpad (eg, a MacBook), where precise clicks are hard:

- Polygons can be closed by double-clicking, and the last vertex can be
  removed with Backspace/Delete.
- Range and rectangle selectors keep the existing gate (instead of
  collapsing it to nothing) when the user accidentally clicks without
  dragging.
"""

import numpy as np
from matplotlib.widgets import (PolygonSelector, SpanSelector, 
                                RectangleSelector, Cursor)


class EasyPolygonSelector(PolygonSelector):
    """
    A `matplotlib.widgets.PolygonSelector` with a few ease-of-use additions:

    - **Double-click** to close the polygon, instead of having to hit the
      first vertex exactly (hard on a laptop trackpad).
    - **Backspace** / **Delete** removes the most recently placed vertex
      while the polygon is still being drawn.
    - **Escape** (inherited) starts over with a new polygon.
    - A completed polygon's vertices can be dragged to adjust them, and
      a right-click (two-finger click on a trackpad) on a vertex removes it.
    """

    _undo_keys = ('backspace', 'delete')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._swallow_release = False

    def _press(self, event):
        self._swallow_release = False
        if (getattr(event, 'dblclick', False)
                and not self._selection_completed
                and self._active_handle_idx < 0):
            # the first click of the double-click already placed a vertex.
            # close the polygon now, and ignore the release that follows (if
            # there is one) so it doesn't add a duplicate vertex.
            self._swallow_release = True
            verts = list(self.verts)

            # drop any (near-)duplicate trailing vertices
            while len(verts) >= 2 and np.allclose(verts[-1], verts[-2]):
                verts.pop()

            if len(verts) >= 3:
                self.verts = verts
                self.onselect(self.verts)
            return
        super()._press(event)

    def _release(self, event):
        if self._swallow_release:
            self._swallow_release = False
            return
        super()._release(event)

    def _on_key_release(self, event):
        if (event.key in self._undo_keys
                and not self._selection_completed
                and len(self._xys) > 1):
            # self._xys[-1] is the "pending" vertex that follows the cursor
            self._xys.pop(-2)
            self._draw_polygon()
            return

        super()._on_key_release(event)


class _KeepOnClickMixin:
    """
    In interactive mode, matplotlib's span and rectangle selectors throw away
    the current selection when the user clicks (without dragging) outside of
    it.  For a gate, that's almost never what the user wants, so restore the
    previous selection instead.  (Dragging out a new selection still works.)
    """

    def _press(self, event):
        self._extents_before_press = \
            self.extents if self._selection_completed else None
        return super()._press(event)

    def _release(self, event):
        ret = super()._release(event)
        saved = getattr(self, '_extents_before_press', None)
        self._extents_before_press = None
        if saved is not None and not self._selection_completed:
            self.extents = saved
            self._selection_completed = True
            self.set_visible(True)
            self.update()
        return ret


class EasySpanSelector(_KeepOnClickMixin, SpanSelector):
    """A `SpanSelector` that keeps the existing span on a stray click."""


class EasyRectangleSelector(_KeepOnClickMixin, RectangleSelector):
    """A `RectangleSelector` that keeps the existing box on a stray click."""


class GateLineDragger:
    """
    Interactive placement of the line(s) of a threshold or quadrant gate.

    - If there's no gate yet, a crosshair follows the mouse and a click
      places the gate.
    - Once there is a gate, **drag its line** (or, for a quadrant gate, the
      intersection) to move it.  A stray single click elsewhere is ignored,
      so the gate isn't moved by accident; double-click to jump the gate
      to a new location.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes`
    get_x, get_y : callable
        Return the gate's current x / y position (or `None`).
    on_set : callable
        Called as ``on_set(x, y)`` when the user places the gate.  The
        coordinate that wasn't moved is passed as `None`.
    use_x, use_y : bool
        Whether the gate has a vertical (x) and/or horizontal (y) line.
    """

    def __init__(self, ax, get_x, get_y, on_set, use_x = True, use_y = False,
                 grab_range = 15, **line_props):
        self.ax = ax
        self.get_x = get_x
        self.get_y = get_y
        self.on_set = on_set
        self.use_x = use_x
        self.use_y = use_y
        self.grab_range = grab_range
        self._drag = None

        self.cursor = Cursor(ax,
                             horizOn = use_y,
                             vertOn = use_x,
                             useblit = True,
                             **line_props)
        self.cursor.connect_event('button_press_event', self._on_press)
        self.cursor.connect_event('button_release_event', self._on_release)
        self._update_cursor()

    def _has_gate(self):
        return ((not self.use_x or self.get_x() is not None) and
                (not self.use_y or self.get_y() is not None))

    def _update_cursor(self):
        if self._drag:
            self.cursor.vertOn, self.cursor.horizOn = self._drag
            self.cursor.visible = True
        else:
            self.cursor.vertOn, self.cursor.horizOn = self.use_x, self.use_y
            # only show the crosshair if there's no gate yet
            self.cursor.visible = not self._has_gate()

    def _usable(self, event):
        return (self.cursor.active
                and event.inaxes is self.ax
                and event.button == 1
                and event.xdata is not None
                and event.ydata is not None
                and self.ax.figure.canvas.widgetlock.available(self.cursor))

    def _on_press(self, event):
        if not self._usable(event):
            return

        if not self._has_gate() or getattr(event, 'dblclick', False):
            self._drag = (self.use_x, self.use_y)
        else:
            # is the press near one of the existing lines?
            near_x = near_y = False
            gx = self.get_x() if self.use_x else None
            gy = self.get_y() if self.use_y else None
            px, py = self.ax.transData.transform(
                (gx if gx is not None else event.xdata,
                 gy if gy is not None else event.ydata))
            if self.use_x and np.isfinite(px):
                near_x = abs(event.x - px) <= self.grab_range
            if self.use_y and np.isfinite(py):
                near_y = abs(event.y - py) <= self.grab_range
            if not (near_x or near_y):
                return  # a stray click -- leave the gate alone
            self._drag = (near_x, near_y)

        self._update_cursor()
        self.cursor.onmove(event)

    def _on_release(self, event):
        if not self._drag:
            return

        move_x, move_y = self._drag
        self._drag = None

        if (event.inaxes is self.ax
                and event.xdata is not None
                and event.ydata is not None):
            self.on_set(float(event.xdata) if move_x else None,
                        float(event.ydata) if move_y else None)

        self._update_cursor()
        self.cursor.onmove(event)

    def disconnect(self):
        self.cursor.disconnect_events()
        self.cursor.visible = False
