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
from matplotlib.widgets import PolygonSelector, SpanSelector, RectangleSelector


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
