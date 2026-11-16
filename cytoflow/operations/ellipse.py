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
cytoflow.operations.ellipse
---------------------------

Apply an ellipse gate to two channels in an `Experiment`.
`ellipse` has three classes:

`EllipseOp` -- Applies the gate, given a center, width, height and angle.

`ScatterplotEllipseSelectionView` -- an `IView` that allows you to view the
ellipse and/or interactively draw it on a scatterplot.

`DensityEllipseSelectionView` -- an `IView` that allows you to view the
ellipse and/or interactively draw it on a density plot.
"""

from traits.api import (HasStrictTraits, Str, Float, Tuple, provides,
                        Instance, Bool, observe, Any, Dict, Constant)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import cytoflow.utility as util
from cytoflow.views import ISelectionView, ScatterplotView, DensityView

from .i_operation import IOperation
from .base_op_views import Op2DView

from ._selectors import EasyEllipseSelector

@provides(IOperation)
class EllipseOp(HasStrictTraits):
    """
    Apply an ellipse gate to a cytometry experiment.

    The ellipse is defined in *scaled* space -- that is, after the data has
    been transformed by `xscale` and `yscale`.  That way, an ellipse that
    you draw on a plot with (say) logicle-scaled axes looks like an ellipse
    on the plot, and the gate includes exactly the events that you see
    inside it.

    Attributes
    ----------
    name : Str
        The operation name.  Used to name the new metadata field in the
        experiment that's created by `apply`

    xchannel, ychannel : Str
        The names of the x and y channels to apply the gate.

    xscale, yscale : {'linear', 'log', 'logicle'} (default = 'linear')
        The scales applied to the data before drawing the ellipse.

    center : (Float, Float)
        The center of the ellipse, in *data* units (ie, the same units as
        the channels' data and as `PolygonOp.vertices`.)  The center is
        transformed by `xscale` and `yscale` before the gate is applied.

    width, height : Float
        The full width and height of the ellipse (ie, twice the semi-axes),
        in *scaled* units.  For example, if `xscale` is ``log``, a
        `width` of ``1.0`` is one decade.  ``width`` is measured along the
        ellipse's (rotated) x axis, ``height`` along its (rotated) y axis.

    angle : Float (default = 0.0)
        The counter-clockwise rotation of the ellipse, in degrees, in
        scaled space.

    Notes
    -----
    It's much easier to use the interactive view you get from `default_view`
    to draw the ellipse.  Set `ScatterplotEllipseSelectionView.interactive`
    to `True`, then drag out the ellipse's bounding box.  Drag the handles
    to change its size, drag the inside to move it, and press **r** and then
    drag a corner handle to rotate it.  You'll need to do this in a Jupyter
    notebook with ``%matplotlib notebook`` -- see the ``Interactive Plots``
    demo for an example.

    Examples
    --------

    .. plot::
        :context: close-figs

        Make a little data set.

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Create and parameterize the operation.

    .. plot::
        :context: close-figs

        >>> e = flow.EllipseOp(name = "Ellipse",
        ...                    xchannel = "V2-A",
        ...                    ychannel = "Y2-A",
        ...                    xscale = "log",
        ...                    yscale = "log",
        ...                    center = (50.0, 5000.0),
        ...                    width = 1.8,
        ...                    height = 0.7,
        ...                    angle = 10.0)

    Show the default view.

    .. plot::
        :context: close-figs

        >>> df = e.default_view(huefacet = "Dox")
        >>> df.plot(ex)

    .. note::
       If you want to use the interactive default view in a Jupyter notebook,
       make sure you say ``%matplotlib notebook`` in the first cell
       (instead of ``%matplotlib inline`` or similar).  Then call
       `default_view` with ``interactive = True``::

           df = e.default_view(huefacet = "Dox",
                               interactive = True)
           df.plot(ex)

    Apply the gate, and show the result

    .. plot::
        :context: close-figs

        >>> ex2 = e.apply(ex)
        >>> ex2.data.groupby('Ellipse', observed = True).size()
        Ellipse
        False    17194
        True      2806
        dtype: int64

    You can also get (or draw) the ellipse on a density plot instead of a
    scatterplot:

    .. plot::
        :context: close-figs

        >>> df = e.default_view(density = True)
        >>> df.plot(ex)
    """

    # traits
    id = Constant('cytoflow.operations.ellipse')
    friendly_id = Constant("Ellipse")

    name = Str
    xchannel = Str
    ychannel = Str

    center = Tuple(Float, Float)
    width = Float(0.0)
    height = Float(0.0)
    angle = Float(0.0)

    xscale = util.ScaleEnum()
    yscale = util.ScaleEnum()

    _selection_view = Instance('_EllipseSelection', transient = True)

    def apply(self, experiment):
        """Applies the gate to an experiment.

        Parameters
        ----------
        experiment : Experiment
            the old `Experiment` to which this op is applied

        Returns
        -------
        Experiment
            a new 'Experiment`, the same as ``experiment`` but with
            a new column of type `bool` with the same as the operation name.
            The bool is ``True`` if the event's measurement is within the
            ellipse (or on its edge), and ``False`` otherwise.

        Raises
        ------
        CytoflowOpError
            if for some reason the operation can't be applied to this
            experiment. The reason is in the ``args`` attribute.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        # make sure name got set!
        if not self.name:
            raise util.CytoflowOpError('name',
                                       "You have to set the Ellipse gate's name "
                                       "before applying it!")

        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "{} is in the experiment already!"
                                       .format(self.name))

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores."
                                       .format(self.name))

        if not self.xchannel:
            raise util.CytoflowOpError('xchannel',
                                       "Must specify an x channel")

        if not self.ychannel:
            raise util.CytoflowOpError('ychannel',
                                       "Must specify a y channel")

        if not self.xchannel in experiment.channels:
            raise util.CytoflowOpError('xchannel',
                                       "xchannel {0} is not in the experiment"
                                       .format(self.xchannel))

        if not self.ychannel in experiment.channels:
            raise util.CytoflowOpError('ychannel',
                                       "ychannel {0} is not in the experiment"
                                       .format(self.ychannel))

        if not (self.width > 0 and self.height > 0):
            raise util.CytoflowOpError('width',
                                       "The ellipse's width and height must "
                                       "both be > 0")

        if not np.all(np.isfinite(self.center)) or not np.isfinite(self.angle):
            raise util.CytoflowOpError('center',
                                       "The ellipse's center and angle must "
                                       "be finite numbers")

        # just like PolygonOp, the gate is defined in scaled space: transform
        # the data (and the center), then test ellipse membership.
        xscale = util.scale_factory(self.xscale, experiment, channel = self.xchannel)
        yscale = util.scale_factory(self.yscale, experiment, channel = self.ychannel)

        cx = float(xscale(np.array([self.center[0]], dtype = float))[0])
        cy = float(yscale(np.array([self.center[1]], dtype = float))[0])

        if not (np.isfinite(cx) and np.isfinite(cy)):
            raise util.CytoflowOpError('center',
                                       "The ellipse's center {} is outside "
                                       "the domain of the scale"
                                       .format(self.center))

        x = np.asarray(xscale(experiment.data[self.xchannel].values), dtype = float)
        y = np.asarray(yscale(experiment.data[self.ychannel].values), dtype = float)

        in_ellipse = ellipse_contains(np.column_stack((x, y)),
                                      (cx, cy),
                                      self.width,
                                      self.height,
                                      self.angle)

        new_experiment = experiment.clone(deep = False)
        new_experiment.add_condition(self.name, "bool", in_ellipse)
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))

        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns an `IView` that allows a user to view the ellipse or interactively draw it.

        Parameters
        ----------

        density : bool, default = False
            If `True`, return a density plot instead of a scatterplot.
        """

        density = kwargs.pop('density', False)
        if density:
            self._selection_view = DensityEllipseSelectionView(op = self)
        else:
            self._selection_view = ScatterplotEllipseSelectionView(op = self)

        self._selection_view.trait_set(**kwargs)
        return self._selection_view


def ellipse_contains(xy, center, width, height, angle):
    """
    Which points in ``xy`` are inside (or on the edge of) an ellipse?

    Parameters
    ----------
    xy : array, shape (N, 2)
    center : (float, float)
    width, height : float
        The full lengths of the ellipse's axes.
    angle : float
        The counter-clockwise rotation, in degrees.

    Returns
    -------
    numpy.ndarray of bool, shape (N, )
    """
    xy = np.asarray(xy, dtype = float)
    t = np.deg2rad(angle)
    dx = xy[:, 0] - center[0]
    dy = xy[:, 1] - center[1]

    # rotate the points by -angle, so the ellipse is axis-aligned
    u = dx * np.cos(t) + dy * np.sin(t)
    v = -dx * np.sin(t) + dy * np.cos(t)

    with np.errstate(invalid = 'ignore'):
        r = (u / (width / 2.0)) ** 2 + (v / (height / 2.0)) ** 2
        # be a little generous at the edge to avoid floating-point issues
        return np.nan_to_num(r, nan = np.inf) <= 1.0 + 1e-12


def _normalize_axes(m):
    """
    Given a 2x2 matrix ``m`` that maps the unit circle onto an ellipse
    (ie, ``ellipse = center + m @ unit_circle``), return the ellipse's
    ``(width, height, angle)``, choosing the axis closest to the x axis as
    the "width" axis so that ``-45 < angle <= 45``.
    """
    u, s, _ = np.linalg.svd(m)
    angle = np.rad2deg(np.arctan2(u[1, 0], u[0, 0]))
    width, height = 2 * s[0], 2 * s[1]

    # fold into (-90, 90]
    if angle <= -90:
        angle += 180
    elif angle > 90:
        angle -= 180

    # then into (-45, 45], swapping axes
    if angle > 45:
        angle -= 90
        width, height = height, width
    elif angle <= -45:
        angle += 90
        width, height = height, width

    return float(width), float(height), float(angle)

def _rotation(deg):
    t = np.deg2rad(deg)
    return np.array([[np.cos(t), -np.sin(t)],
                     [np.sin(t), np.cos(t)]])

def _selector_to_scaled(half_width, half_height, rotation, aspect):
    """
    Convert an `EllipseSelector`'s shape to a ``(width, height, angle)``
    in scaled (data) space.

    `matplotlib.widgets.EllipseSelector` rotates its ellipse in *display*
    space -- it scales the y axis by ``aspect`` (the display aspect ratio
    of the axes), rotates, then scales it back.
    """
    m = np.diag([1.0, 1.0 / aspect]) \
        @ _rotation(rotation) \
        @ np.diag([half_width, half_height * aspect])
    return _normalize_axes(m)

def _scaled_to_selector(width, height, angle, aspect):
    """
    The inverse of `_selector_to_scaled`: returns
    ``(half_width, half_height, rotation)`` for an `EllipseSelector`.
    """
    # in display-proportional space, the ellipse is a "true" rotated ellipse
    m = np.diag([1.0, aspect]) \
        @ _rotation(angle) \
        @ np.diag([width / 2.0, height / 2.0])
    w, h, rot = _normalize_axes(m)
    return w / 2.0, h / 2.0 / aspect, rot


class _EllipseSelection(Op2DView):
    xfacet = Constant(None)
    yfacet = Constant(None)

    interactive = Bool(False, transient = True)

    # internal state.
    _ax = Any(transient = True)

    # the matplotlib widgets work in data coordinates -- but the axes may
    # have nonlinear (ie, log or logicle) scales.  so we put an invisible,
    # linear set of axes on top of the plot, whose data coordinates are
    # the *scaled* coordinates.  the selector and the ellipse patch live
    # there, so the ellipse looks like an ellipse.
    _overlay = Any(transient = True)
    _xscale = Any(transient = True)
    _yscale = Any(transient = True)
    _lim_cids = Any(transient = True)

    _widget = Instance(EasyEllipseSelector, transient = True)
    _patch = Instance(mpl.patches.Ellipse, transient = True)
    _patch_props = Dict()
    _selector_props = Dict()

    def plot(self, experiment, **kwargs):

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        self._patch_props = kwargs.pop('patch_props',
                                        {'edgecolor' : 'black',
                                         'linewidth' : 2,
                                         'fill' : False})

        self._selector_props = kwargs.pop('selector_props',
                                          {'facecolor' : 'blue',
                                           'edgecolor' : 'blue',
                                           'alpha' : 0.15,
                                           'linewidth' : 2})

        super().plot(experiment, **kwargs)
        self._ax = plt.gca()

        # these are the same scales that EllipseOp.apply uses
        self._xscale = util.scale_factory(self.op.xscale,
                                          experiment,
                                          channel = self.op.xchannel)
        self._yscale = util.scale_factory(self.op.yscale,
                                          experiment,
                                          channel = self.op.ychannel)

        self._make_overlay()
        self._draw_ellipse(None)
        self._interactive(None)

    def _make_overlay(self):
        self._overlay = self._ax.inset_axes([0, 0, 1, 1])
        self._overlay.set_axis_off()
        self._overlay.patch.set_alpha(0.0)
        self._overlay.set_navigate(False)
        self._sync_overlay(None)
        self._lim_cids = (self._ax.callbacks.connect('xlim_changed', self._sync_overlay),
                          self._ax.callbacks.connect('ylim_changed', self._sync_overlay))

    def _sync_overlay(self, _):
        if not self._overlay:
            return

        xlim = self._to_scaled(self._xscale, self._ax.get_xlim())
        ylim = self._to_scaled(self._yscale, self._ax.get_ylim())
        self._overlay.set_xlim(xlim)
        self._overlay.set_ylim(ylim)

    @staticmethod
    def _to_scaled(scale, values):
        return [float(v) for v in scale(np.asarray(values, dtype = float))]

    @staticmethod
    def _from_scaled(scale, values):
        return [float(v) for v in scale.inverse(np.asarray(values, dtype = float))]

    def _scaled_center(self):
        cx = self._to_scaled(self._xscale, [self.op.center[0]])[0]
        cy = self._to_scaled(self._yscale, [self.op.center[1]])[0]
        return cx, cy

    def _has_gate(self):
        return (self.op.width > 0 and self.op.height > 0
                and np.all(np.isfinite(self._scaled_center())))

    @observe('[op.center,op.width,op.height,op.angle]', post_init = True)
    def _draw_ellipse(self, _):
        if not self._overlay:
            return

        if self._patch and self._patch in self._overlay.patches:
            self._patch.remove()
            self._patch = None

        if not self._has_gate():
            return

        self._patch = mpl.patches.Ellipse(self._scaled_center(),
                                          self.op.width,
                                          self.op.height,
                                          angle = self.op.angle,
                                          **self._patch_props)

        self._overlay.add_patch(self._patch)
        plt.draw()

    @observe('interactive', post_init = True)
    def _interactive(self, _):
        if self._overlay and self.interactive:
            if self._widget:
                self._widget.set_active(False)

            # interactive = True lets the user drag the edges of an existing
            # gate (or drag the whole ellipse) to adjust it.
            self._widget = EasyEllipseSelector(
                                self._overlay,
                                onselect = self._onselect,
                                interactive = True,
                                drag_from_anywhere = True,
                                grab_range = 15,
                                spancoords = 'pixels',
                                minspanx = 5,
                                minspany = 5,
                                props = self._selector_props,
                                handle_props = dict(markeredgecolor = 'blue',
                                                    markerfacecolor = 'white'),
                                useblit = True)

            # if there's already a gate, load it into the selector so the
            # user can adjust it instead of having to re-draw it.
            if self._has_gate():
                aspect = self._overlay._get_aspect_ratio()
                hw, hh, rot = _scaled_to_selector(self.op.width,
                                                  self.op.height,
                                                  self.op.angle,
                                                  aspect)
                cx, cy = self._scaled_center()
                self._widget._rotation = np.deg2rad(rot)
                self._widget._set_aspect_ratio_correction()
                self._widget.extents = (cx - hw, cx + hw, cy - hh, cy + hh)
                self._widget._selection_completed = True

        elif self._widget:
            self._widget.set_active(False)
            self._widget.set_visible(False)
            self._widget = None

    def _get_selection(self):
        """
        Returns the selector's ellipse as ``(center, width, height, angle)``,
        with ``center`` in data units and the rest in scaled units; or `None`
        if there's no (or only a degenerate) ellipse.
        """
        x0, x1, y0, y1 = self._widget.extents
        hw, hh = (x1 - x0) / 2.0, (y1 - y0) / 2.0
        if not (hw > 0 and hh > 0):
            return None

        width, height, angle = \
            _selector_to_scaled(hw, hh, self._widget.rotation,
                                self._overlay._get_aspect_ratio())

        cx, cy = self._widget.center
        center = (self._from_scaled(self._xscale, [cx])[0],
                  self._from_scaled(self._yscale, [cy])[0])

        if not np.all(np.isfinite(center)):
            return None

        return center, width, height, angle

    def _onselect(self, _press, _release):
        sel = self._get_selection()

        # ignore a click without a drag
        if sel is None:
            return

        center, width, height, angle = sel
        self.op.trait_set(center = center,
                          width = width,
                          height = height,
                          angle = angle)


@provides(ISelectionView)
class ScatterplotEllipseSelectionView(_EllipseSelection, ScatterplotView):
    """
    Plots, and lets the user interact with, an ellipse gate on a scatterplot.

    Attributes
    ----------
    interactive : bool
        is this view interactive?  Ie, can the user draw (and adjust) the
        ellipse with the mouse?

    Examples
    --------

    In a Jupyter notebook with ``%matplotlib notebook``

    >>> e = flow.EllipseOp(name = "Ellipse",
    ...                    xchannel = "V2-A",
    ...                    ychannel = "Y2-A")
    >>> ev = e.default_view()
    >>> ev.plot(ex2)
    >>> ev.interactive = True
    """

    id = Constant('cytoflow.views.ellipse')
    friendly_id = Constant("Ellipse Selection")

    def plot(self, experiment, **kwargs):
        """
        Plot the default view, and then draw the selection on top of it.

        Parameters
        ----------

        patch_props : Dict
           The properties of the `matplotlib.patches.Ellipse` that is drawn
           on top of the scatterplot.  They're passed
           directly to the `matplotlib.patches.Ellipse` constructor.
           Default: ``{edgecolor : 'black', linewidth : 2, fill : False}``

        selector_props : Dict
           The properties of the interactive selector's ellipse.  They're
           passed directly to the `matplotlib.patches.Ellipse` constructor.
           Default: ``{facecolor : 'blue', edgecolor : 'blue', alpha : 0.15, linewidth : 2}``

        """
        super().plot(experiment, **kwargs)

util.expand_class_attributes(ScatterplotEllipseSelectionView)
util.expand_method_parameters(ScatterplotEllipseSelectionView, ScatterplotEllipseSelectionView.plot)

@provides(ISelectionView)
class DensityEllipseSelectionView(_EllipseSelection, DensityView):
    """
    Plots, and lets the user interact with, an ellipse gate on a density plot.

    Attributes
    ----------
    interactive : bool
        is this view interactive?  Ie, can the user draw (and adjust) the
        ellipse with the mouse?

    Examples
    --------

    In a Jupyter notebook with ``%matplotlib notebook``

    >>> e = flow.EllipseOp(name = "Ellipse",
    ...                    xchannel = "V2-A",
    ...                    ychannel = "Y2-A")
    >>> ev = e.default_view(density = True)
    >>> ev.plot(ex2)
    >>> ev.interactive = True
    """

    id = Constant('cytoflow.views.ellipse_density')
    friendly_id = Constant("Ellipse Selection (Density plot)")

    def plot(self, experiment, **kwargs):
        """
        Plot the default view, and then draw the selection on top of it.

        Parameters
        ----------

        patch_props : Dict
           The properties of the `matplotlib.patches.Ellipse` that is drawn
           on top of the density plot.  They're passed
           directly to the `matplotlib.patches.Ellipse` constructor.
           Default: ``{edgecolor : 'white', linewidth : 2, fill : False}``

        selector_props : Dict
           The properties of the interactive selector's ellipse.  They're
           passed directly to the `matplotlib.patches.Ellipse` constructor.
           Default: ``{facecolor : 'white', edgecolor : 'white', alpha : 0.25, linewidth : 2}``

        """

        patch_props = kwargs.pop('patch_props',
                                 {'edgecolor' : 'white',
                                  'linewidth' : 2,
                                  'fill' : False})

        selector_props = kwargs.pop('selector_props',
                                    {'facecolor' : 'white',
                                     'edgecolor' : 'white',
                                     'alpha' : 0.25,
                                     'linewidth' : 2})

        super().plot(experiment, patch_props = patch_props, selector_props = selector_props, **kwargs)

util.expand_class_attributes(DensityEllipseSelectionView)
util.expand_method_parameters(DensityEllipseSelectionView, DensityEllipseSelectionView.plot)

if __name__ == '__main__':
    import cytoflow as flow
    tube1 = flow.Tube(file = '../../cytoflow/tests/data/Plate01/RFP_Well_A3.fcs',
                      conditions = {"Dox" : 10.0})

    tube2 = flow.Tube(file = '../../cytoflow/tests/data/Plate01/CFP_Well_A4.fcs',
                      conditions = {"Dox" : 1.0})

    ex = flow.ImportOp(conditions = {"Dox" : "float"}, tubes = [tube1, tube2]).apply()

    e = EllipseOp(name = "Ellipse",
                  xchannel = "V2-A",
                  ychannel = "Y2-A")
    v = e.default_view(xscale = "logicle", yscale = "logicle")

    plt.ioff()
    v.plot(ex)
    v.interactive = True
    plt.show()
    print(e.center, e.width, e.height, e.angle)
