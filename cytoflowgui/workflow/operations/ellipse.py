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
cytoflowgui.workflow.operations.ellipse
---------------------------------------

"""

from traits.api import (provides, Instance, Str, observe, Float, Bool, Enum,
                        Constant, Tuple, Property, DelegatesTo)

from cytoflow.operations.ellipse import (EllipseOp, ScatterplotEllipseSelectionView,
                                         DensityEllipseSelectionView, Op2DView,
                                         _EllipseSelection)
import cytoflow.utility as util

from ..views import IWorkflowView, WorkflowView
from ..views.view_base import Data2DPlotParams
from ..views.scatterplot import SCATTERPLOT_MARKERS
from ..serialization import camel_registry, traits_str, cytoflow_class_repr, dedent

from .operation_base import IWorkflowOperation, WorkflowOperation

EllipseOp.__repr__ = cytoflow_class_repr

class EllipsePlotParams(Data2DPlotParams):
    density = Bool(False)

    # density-specific
    gridsize = util.PositiveCInt(50, allow_zero = False)
    smoothed = Bool(False)
    smoothed_sigma = util.PositiveCFloat(1.0, allow_zero = False)

    # scatterplot-specific
    alpha = util.PositiveCFloat(0.25)
    s = util.PositiveCFloat(2)
    marker = Enum(SCATTERPLOT_MARKERS)

class _EllipseSelectionWorkflowView(_EllipseSelection):
    parent = Instance(IWorkflowView)
    subset = DelegatesTo('parent')

    # data flow: user draws (or adjusts) the ellipse.  the remote canvas
    # calls _onselect, which sets the parent's _ellipse.  _ellipse is copied
    # back to the local view (because it's "status = True"). _update_ellipse
    # is called, and the operation is updated.  the operation's _ellipse is
    # sent back to the remote operation (because "apply = True"), where the
    # remote operation is updated.

    # center, width, height and angle are properties (both in the view and
    # in the operation) so that the update happens atomically.

    def _onselect(self, _press, _release):
        sel = self._get_selection()

        # ignore a click without a drag
        if sel is None:
            return

        (cx, cy), width, height, angle = sel
        self.parent._ellipse = (cx, cy, width, height, angle)

@provides(IWorkflowView)
class ScatterplotEllipseSelectionWorkflowView(_EllipseSelectionWorkflowView, ScatterplotEllipseSelectionView):
    pass

@provides(IWorkflowView)
class DensityEllipseSelectionWorkflowView(_EllipseSelectionWorkflowView, DensityEllipseSelectionView):
    pass

@provides(IWorkflowView)
class EllipseSelectionView(WorkflowView, Op2DView):
    id = Constant('cytoflowgui.workflow.operations.ellipseview')

    op = Instance(IWorkflowOperation, fixed = True)
    plot_params = Instance(EllipsePlotParams, ())

    _view = Instance(IWorkflowView)
    _ellipse = Tuple(Float, Float, Float, Float, Float, status = True)

    interactive = Bool(False, transient = True)

    @observe('xchannel,ychannel,xscale,yscale', post_init = True)
    def _reset_ellipse(self, _):
        self._ellipse = (0.0, 0.0, 0.0, 0.0, 0.0)
        self.op._ellipse = (0.0, 0.0, 0.0, 0.0, 0.0)

    def clear_estimate(self):
        # no-op
        return

    def plot(self, experiment, **kwargs):
        density = kwargs.pop('density')
        if density:
            kwargs.pop('alpha')
            kwargs.pop('s')
            kwargs.pop('marker')
            if not self._view or not isinstance(self._view, DensityEllipseSelectionWorkflowView):
                self._view = DensityEllipseSelectionWorkflowView(op = self.op,
                                                                 parent = self,
                                                                 interactive = self.interactive,
                                                                 huescale = self.huescale)
            kwargs['patch_props'] = {'edgecolor' : 'white', 'linewidth' : 3, 'fill' : False}
        else:
            kwargs.pop('gridsize')
            kwargs.pop('smoothed')
            kwargs.pop('smoothed_sigma')
            if not self._view or not isinstance(self._view, ScatterplotEllipseSelectionWorkflowView):
                self._view = ScatterplotEllipseSelectionWorkflowView(op = self.op,
                                                                     parent = self,
                                                                     interactive = self.interactive,
                                                                     huefacet = self.huefacet)
            kwargs['patch_props'] = {'edgecolor' : 'black', 'linewidth' : 3, 'fill' : False}

        self._view.plot(experiment, **kwargs)

    def get_notebook_code(self, idx):
        if self.plot_params.density:
            view = DensityEllipseSelectionWorkflowView()
            plot_params = self.plot_params.clone_traits(copy = "deep")
            plot_params.reset_traits(traits = ['density', 'alpha', 's', 'marker'])
        else:
            view = ScatterplotEllipseSelectionWorkflowView()
            plot_params = self.plot_params.clone_traits(copy = "deep")
            plot_params.reset_traits(traits = ['density', 'gridsize', 'smoothed', 'smoothed_sigma'])

        view.copy_traits(self, view.copyable_trait_names())
        plot_params_str = traits_str(plot_params)

        return dedent("""
        op_{idx}.default_view({traits}{comma}{density}).plot(ex_{prev_idx}{plot_params})
        """
        .format(idx = idx,
                traits = traits_str(view),
                comma = ", " if traits_str(view) and self.plot_params.density else "",
                density = "density = True" if self.plot_params.density else "",
                prev_idx = idx - 1,
                plot_params = ", " + plot_params_str if plot_params_str else ""))

    @observe('interactive', post_init = True)
    def _interactive(self, _):
        if self._view:
            self._view.interactive = self.interactive

    @observe('_ellipse', post_init = True)
    def _update_ellipse(self, _):
        self.op._ellipse = self._ellipse

    @observe('huefacet', post_init = True)
    def _update_huefacet(self, _):
        if self._view:
            self._view.huefacet = self.huefacet

    @observe('huescale', post_init = True)
    def _update_huescale(self, _):
        if self._view:
            self._view.huescale = self.huescale


@provides(IWorkflowOperation)
class EllipseWorkflowOp(WorkflowOperation, EllipseOp):
    name = Str(apply = True)
    xchannel = Str(apply = True)
    ychannel = Str(apply = True)

    # the ellipse is (center x, center y, width, height, angle).  it's a
    # single trait so that updating it from the interactive view is atomic.
    _ellipse = Tuple(Float, Float, Float, Float, Float, apply = True)
    center = Property(Tuple(Float, Float), observe = '_ellipse')
    width = Property(Float, observe = '_ellipse')
    height = Property(Float, observe = '_ellipse')
    angle = Property(Float, observe = '_ellipse')

    # the ellipse is defined in scaled space, so the scales are part of
    # the operation, not just the view.
    xscale = util.ScaleEnum(apply = True)
    yscale = util.ScaleEnum(apply = True)

    def _get_center(self):
        return (self._ellipse[0], self._ellipse[1])

    def _set_center(self, val):
        self._ellipse = (float(val[0]), float(val[1])) + self._ellipse[2:]

    def _get_width(self):
        return self._ellipse[2]

    def _set_width(self, val):
        self._ellipse = self._ellipse[0:2] + (float(val),) + self._ellipse[3:]

    def _get_height(self):
        return self._ellipse[3]

    def _set_height(self, val):
        self._ellipse = self._ellipse[0:3] + (float(val),) + self._ellipse[4:]

    def _get_angle(self):
        return self._ellipse[4]

    def _set_angle(self, val):
        self._ellipse = self._ellipse[0:4] + (float(val),)

    def default_view(self, **kwargs):
        return EllipseSelectionView(op = self, **kwargs)

    def get_notebook_code(self, idx):
        op = EllipseOp()
        op.copy_traits(self, op.copyable_trait_names())

        return dedent("""
        op_{idx} = {repr}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1))


### Serialization
@camel_registry.dumper(EllipseWorkflowOp, 'ellipse', version = 1)
def _dump(op):
    return dict(name = op.name,
                xchannel = op.xchannel,
                ychannel = op.ychannel,
                center = op.center,
                width = op.width,
                height = op.height,
                angle = op.angle,
                xscale = op.xscale,
                yscale = op.yscale)

@camel_registry.loader('ellipse', version = 1)
def _load(data, version):
    data['center'] = tuple(data['center'])
    return EllipseWorkflowOp(**data)

@camel_registry.dumper(EllipseSelectionView, 'ellipse-view', version = 1)
def _dump_view(view):
    return dict(op = view.op,
                huefacet = view.huefacet,
                huescale = view.huescale,
                subset_list = view.subset_list,
                plot_params = view.plot_params,
                current_plot = view.current_plot)

@camel_registry.loader('ellipse-view', version = any)
def _load_view(data, version):
    return EllipseSelectionView(**data)

@camel_registry.dumper(EllipsePlotParams, 'ellipse-params', version = 1)
def _dump_view_params(params):
    return dict(density = params.density,

                # BasePlotParams
                title = params.title,
                xlabel = params.xlabel,
                ylabel = params.ylabel,
                huelabel = params.huelabel,
                col_wrap = params.col_wrap,
                sns_style = params.sns_style,
                sns_context = params.sns_context,
                legend = params.legend,
                sharex = params.sharex,
                sharey = params.sharey,
                despine = params.despine,

                # DataplotParams
                min_quantile = params.min_quantile,
                max_quantile = params.max_quantile,

                # Data2DPlotParams
                xlim = params.xlim,
                ylim = params.ylim,

                # Scatterplot params
                alpha = params.alpha,
                s = params.s,
                marker = params.marker,

                 # Density plot params
                gridsize = params.gridsize,
                smoothed = params.smoothed,
                smoothed_sigma = params.smoothed_sigma )

@camel_registry.loader('ellipse-params', version = any)
def _load_params(data, version):
    return EllipsePlotParams(**data)
