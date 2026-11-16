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
cytoflow.operations.freehand
----------------------------

Apply a freehand (lasso) gate to two channels in an `Experiment`.
`freehand` has three classes:

`FreehandOp` -- Applies the gate, given a set of vertices.  It's just a
`PolygonOp` whose default view lets you draw the polygon freehand.

`ScatterplotFreehandSelectionView` -- an `IView` that allows you to view the
gate and/or interactively draw it on a scatterplot.

`DensityFreehandSelectionView` -- an `IView` that allows you to view the
gate and/or interactively draw it on a density plot.
"""

from traits.api import provides, Instance, Bool, observe, Constant

import matplotlib.pyplot as plt

import cytoflow.utility as util
from cytoflow.views import ISelectionView, ScatterplotView, DensityView

from .i_operation import IOperation
from .polygon import PolygonOp, _PolygonSelection

from ._selectors import EasyLassoSelector

@provides(IOperation)
class FreehandOp(PolygonOp):
    """
    Apply a freehand (lasso) gate to a cytometry experiment.

    This is a `PolygonOp` -- the gate is a polygon, and it's applied
    exactly the same way -- but its default view lets you draw the gate
    by dragging the mouse around the population you want, instead of
    clicking each vertex.  The path you draw is simplified to (at most)
    100 vertices.

    Attributes
    ----------
    name : Str
        The operation name.  Used to name the new metadata field in the
        experiment that's created by `apply`

    xchannel, ychannel : Str
        The names of the x and y channels to apply the gate.

    xscale, yscale : {'linear', 'log', 'logicle'} (default = 'linear')
        The scales applied to the data before drawing the gate.

    vertices : List((Float, Float))
        The gate's vertices.  An ordered list of 2-tuples, representing
        the x and y coordinates of the vertices, in data units.

    Notes
    -----
    Set `ScatterplotFreehandSelectionView.interactive` to `True`, then
    press the mouse button and drag around the population.  The gate is
    closed when you let go.  To change the gate, just draw it again.  You'll
    need to do this in a Jupyter notebook with ``%matplotlib notebook`` --
    see the ``Interactive Plots`` demo for an example.

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

    Create and parameterize the operation.  (Usually, you'd draw the
    vertices with the interactive view.)

    .. plot::
        :context: close-figs

        >>> f = flow.FreehandOp(name = "Freehand",
        ...                     xchannel = "V2-A",
        ...                     ychannel = "Y2-A")
        >>> f.vertices = [(23.411982294776319, 5158.7027015021222),
        ...               (102.22182270573683, 23124.058843387455),
        ...               (510.94519955277201, 23124.058843387455),
        ...               (1089.5215641232173, 3800.3424832180476),
        ...               (340.56382570202402, 801.98947404942271),
        ...               (65.42597937575897, 1119.3133482602157)]

    Show the default view.

    .. plot::
        :context: close-figs

        >>> df = f.default_view(huefacet = "Dox",
        ...                     xscale = 'log',
        ...                     yscale = 'log')

        >>> df.plot(ex)

    .. note::
       If you want to use the interactive default view in a Jupyter notebook,
       make sure you say ``%matplotlib notebook`` in the first cell
       (instead of ``%matplotlib inline`` or similar).  Then call
       `default_view` with ``interactive = True``::

           df = f.default_view(huefacet = "Dox",
                               xscale = 'log',
                               yscale = 'log',
                               interactive = True)
           df.plot(ex)

    Apply the gate, and show the result

    .. plot::
        :context: close-figs

        >>> ex2 = f.apply(ex)
        >>> ex2.data.groupby('Freehand', observed = True).size()
        Freehand
        False    15875
        True      4125
        dtype: int64
    """

    # traits
    id = Constant('cytoflow.operations.freehand')
    friendly_id = Constant("Freehand")

    _selection_view = Instance('_FreehandSelection', transient = True)

    def default_view(self, **kwargs):
        """
        Returns an `IView` that allows a user to view the gate or draw it freehand.

        Parameters
        ----------

        density : bool, default = False
            If `True`, return a density plot instead of a scatterplot.
        """

        density = kwargs.pop('density', False)
        if density:
            self._selection_view = DensityFreehandSelectionView(op = self)
        else:
            self._selection_view = ScatterplotFreehandSelectionView(op = self)

        self._selection_view.trait_set(**kwargs)
        return self._selection_view


class _FreehandSelection(_PolygonSelection):

    interactive = Bool(False, transient = True)

    # internal state.
    _widget = Instance(EasyLassoSelector, transient = True)

    @observe('interactive', post_init = True)
    def _interactive(self, _):
        if self._ax and self.interactive:
            if self._widget:
                self._widget.set_active(False)

            props = dict(self._selector_props)
            props.pop('linestyle', None)
            self._widget = EasyLassoSelector(self._ax,
                                             self._onselect,
                                             props = props,
                                             useblit = True)

        elif self._widget:
            self._widget.set_active(False)
            self._widget = None

    def _onselect(self, vertices):
        self.op.vertices = [(float(x), float(y)) for (x, y) in vertices]


@provides(ISelectionView)
class ScatterplotFreehandSelectionView(_FreehandSelection, ScatterplotView):
    """
    Plots, and lets the user draw, a freehand gate on a scatterplot.

    Attributes
    ----------
    interactive : bool
        is this view interactive?  Ie, can the user draw the gate with the
        mouse?

    Examples
    --------

    In a Jupyter notebook with ``%matplotlib notebook``

    >>> f = flow.FreehandOp(xchannel = "V2-A",
    ...                     ychannel = "Y2-A")
    >>> fv = f.default_view()
    >>> fv.plot(ex2)
    >>> fv.interactive = True
    """

    id = Constant('cytoflow.views.freehand')
    friendly_id = Constant("Freehand Selection")

    def plot(self, experiment, **kwargs):
        """
        Plot the default view, and then draw the selection on top of it.

        Parameters
        ----------

        patch_props : Dict
           The properties of the `matplotlib.patches.Patch` that are drawn
           on top of the scatterplot.  They're passed
           directly to the `matplotlib.patches.Patch` constructor.
           Default: ``{edgecolor : 'black', linewidth : 2, fill : False}``

        selector_props : Dict
           The properties of the `matplotlib.lines.Line2D` that is drawn
           while the user drags the mouse.  They're passed
           directly to the `matplotlib.lines.Line2D` constructor.
           Default: ``{color : 'black', linestyle : '-', linewidth : 2, alpha = 0.5}``

        """
        super().plot(experiment, **kwargs)

util.expand_class_attributes(ScatterplotFreehandSelectionView)
util.expand_method_parameters(ScatterplotFreehandSelectionView, ScatterplotFreehandSelectionView.plot)

@provides(ISelectionView)
class DensityFreehandSelectionView(_FreehandSelection, DensityView):
    """
    Plots, and lets the user draw, a freehand gate on a density plot.

    Attributes
    ----------
    interactive : bool
        is this view interactive?  Ie, can the user draw the gate with the
        mouse?

    Examples
    --------

    In a Jupyter notebook with ``%matplotlib notebook``

    >>> f = flow.FreehandOp(xchannel = "V2-A",
    ...                     ychannel = "Y2-A")
    >>> fv = f.default_view(density = True)
    >>> fv.plot(ex2)
    >>> fv.interactive = True
    """

    id = Constant('cytoflow.views.freehand_density')
    friendly_id = Constant("Freehand Selection (Density plot)")

    def plot(self, experiment, **kwargs):
        """
        Plot the default view, and then draw the selection on top of it.

        Parameters
        ----------

        patch_props : Dict
           The properties of the `matplotlib.patches.Patch` that are drawn
           on top of the density plot.  They're passed
           directly to the `matplotlib.patches.Patch` constructor.
           Default: {edgecolor : 'white', linewidth : 2, fill : False}

        selector_props : Dict
           The properties of the `matplotlib.lines.Line2D` that is drawn
           while the user drags the mouse.  They're passed
           directly to the `matplotlib.lines.Line2D` constructor.
           Default: ``{color : 'white', linestyle : '-', linewidth : 2, alpha = 0.5}``

        """

        patch_props = kwargs.pop('patch_props',
                                 {'edgecolor' : 'white',
                                  'linewidth' : 2,
                                  'fill' : False})

        selector_props = kwargs.pop('selector_props',
                                    {'color' : 'white',
                                     'linestyle' : '-',
                                     'linewidth' : 2,
                                     'alpha' : 0.5})

        super().plot(experiment, patch_props = patch_props, selector_props = selector_props, **kwargs)

util.expand_class_attributes(DensityFreehandSelectionView)
util.expand_method_parameters(DensityFreehandSelectionView, DensityFreehandSelectionView.plot)

if __name__ == '__main__':
    import cytoflow as flow
    tube1 = flow.Tube(file = '../../cytoflow/tests/data/Plate01/RFP_Well_A3.fcs',
                      conditions = {"Dox" : 10.0})

    tube2 = flow.Tube(file = '../../cytoflow/tests/data/Plate01/CFP_Well_A4.fcs',
                      conditions = {"Dox" : 1.0})

    ex = flow.ImportOp(conditions = {"Dox" : "float"}, tubes = [tube1, tube2]).apply()

    f = FreehandOp(name = "Freehand",
                   xchannel = "V2-A",
                   ychannel = "Y2-A")
    v = f.default_view(xscale = "logicle", yscale = "logicle")

    plt.ioff()
    v.plot(ex)
    v.interactive = True
    plt.show()
    print(f.vertices)
