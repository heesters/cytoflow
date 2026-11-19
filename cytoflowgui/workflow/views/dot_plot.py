#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.workflow.views.dot_plot
------------------------------------

"""

from textwrap import dedent

from traits.api import provides, List, Str, Instance

from cytoflow import DotPlotView
import cytoflow.utility as util

from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr, traits_str
from .view_base import IWorkflowView, WorkflowFacetView, Data2DPlotParams

DotPlotView.__repr__ = cytoflow_class_repr


class DotPlotParams(Data2DPlotParams):
    alpha = util.PositiveCFloat(0.3)
    size = util.PositiveCInt(1, allow_zero=False)


@provides(IWorkflowView)
class DotPlotWorkflowView(WorkflowFacetView, DotPlotView):
    plot_params = Instance(DotPlotParams, ())

    populations = List(Str, apply=True)

    def get_notebook_code(self, idx):
        view = DotPlotView()
        view.copy_traits(self, view.copyable_trait_names())
        plot_params_str = traits_str(self.plot_params)

        return dedent("""
        {repr}.plot(ex_{idx}{plot}{plot_params})
        """.format(repr=repr(view),
                   idx=idx,
                   plot=".subset(" + repr(self.plotfacet) + ", " + repr(self.current_plot) + ")" if self.current_plot else "",
                   plot_params=", " + plot_params_str if plot_params_str else ""))


@camel_registry.dumper(DotPlotWorkflowView, 'dot-plot', version=1)
def _dump(view):
    return dict(xchannel=view.xchannel,
                xscale=view.xscale,
                ychannel=view.ychannel,
                yscale=view.yscale,
                populations=view.populations,
                xfacet=view.xfacet,
                yfacet=view.yfacet,
                huefacet=view.huefacet,
                huescale=view.huescale,
                plotfacet=view.plotfacet,
                subset_list=view.subset_list,
                plot_params=view.plot_params,
                current_plot=view.current_plot)


@camel_registry.loader('dot-plot', version=any)
def _load(data, version):
    return DotPlotWorkflowView(**data)


@camel_registry.dumper(DotPlotParams, 'dot-plot-params', version=1)
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
        xlim=params.xlim,
        ylim=params.ylim,
        alpha=params.alpha,
        size=params.size)


@camel_registry.loader('dot-plot-params', version=any)
def _load_params(data, version):
    return DotPlotParams(**data)
