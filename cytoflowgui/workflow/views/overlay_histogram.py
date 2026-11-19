#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.workflow.views.overlay_histogram
--------------------------------------------

"""

from textwrap import dedent

from traits.api import provides, List, Str, Bool, Int, Instance

from cytoflow import OverlayHistogramView
import cytoflow.utility as util

from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr, traits_str
from .view_base import IWorkflowView, WorkflowView, DataPlotParams

OverlayHistogramView.__repr__ = cytoflow_class_repr


class OverlayHistogramPlotParams(DataPlotParams):
    bins = util.PositiveCInt(200, allow_zero=False)
    alpha = util.PositiveCFloat(0.5)
    density = Bool(True)


@provides(IWorkflowView)
class OverlayHistogramWorkflowView(WorkflowView, OverlayHistogramView):
    plot_params = Instance(OverlayHistogramPlotParams, ())

    populations = List(Str, apply=True)

    def get_notebook_code(self, idx):
        view = OverlayHistogramView()
        view.copy_traits(self, view.copyable_trait_names())
        plot_params_str = traits_str(self.plot_params)

        return dedent("""
        {repr}.plot(ex_{idx}{plot_params})
        """.format(repr=repr(view),
                   idx=idx,
                   plot_params=", " + plot_params_str if plot_params_str else ""))


@camel_registry.dumper(OverlayHistogramWorkflowView, 'overlay-histogram', version=1)
def _dump(view):
    return dict(channel=view.channel,
                scale=view.scale,
                populations=view.populations,
                xfacet=view.xfacet,
                yfacet=view.yfacet,
                huefacet=view.huefacet,
                subset_list=view.subset_list,
                plot_params=view.plot_params,
                current_plot=view.current_plot)


@camel_registry.loader('overlay-histogram', version=any)
def _load(data, version):
    return OverlayHistogramWorkflowView(**data)


@camel_registry.dumper(OverlayHistogramPlotParams, 'overlay-histogram-params', version=1)
def _dump_params(params):
    return dict(
        title=params.title,
        xlabel=params.xlabel,
        ylabel=params.ylabel,
        huelabel=params.huelabel,
        col_wrap=params.col_wrap,
        sns_style=params.sns_style,
        sns_context=params.sns_context,
        palette=params.palette,
        legend=params.legend,
        legend_loc=params.legend_loc,
        sharex=params.sharex,
        sharey=params.sharey,
        margin_titles=params.margin_titles,
        despine=params.despine,
        min_quantile=params.min_quantile,
        max_quantile=params.max_quantile,
        bins=params.bins,
        alpha=params.alpha,
        density=params.density)


@camel_registry.loader('overlay-histogram-params', version=any)
def _load_params(data, version):
    return OverlayHistogramPlotParams(**data)
