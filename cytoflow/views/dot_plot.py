#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.views.dot_plot
-----------------------

A scatter plot where every event is colored by which gate it belongs to.

`DotPlotView` -- the `IView` class that makes the plot.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from traits.api import provides, Constant, Str, List, Float, Int, Any

import cytoflow.utility as util
from .i_view import IView
from .base_views import Base2DView


# Okabe-Ito colorblind-safe palette
_OKABE_ITO = [
    '#E69F00', '#56B4E9', '#009E73', '#F0E442',
    '#0072B2', '#D55E00', '#CC79A7',
]

_UNGATED_COLOR = '#cccccc'


@provides(IView)
class DotPlotView(Base2DView):
    """
    A scatter plot whose events are colored by which boolean gate they satisfy.

    For each event the **last** gate in `populations` that is ``True`` "wins";
    events that belong to no listed gate are drawn in light grey.

    Attributes
    ----------
    populations : list of str
        Boolean conditions (gates) to color by.  Order matters: later entries
        draw on top.

    alpha : float
        Point transparency (default 0.3).

    size : int
        Point size in points² (default 1).

    Examples
    --------

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> flow.DotPlotView(xchannel = 'FSC-A',
        ...                  ychannel = 'SSC-A',
        ...                  xscale = 'logicle',
        ...                  yscale = 'logicle',
        ...                  populations = ['Singlets', 'CD4_pos', 'CD8_pos']).plot(ex)
    """

    id = Constant('cytoflow.view.dot_plot')
    friend_id = Constant("Dot Plot")

    populations = List(Str)
    alpha = Float(0.3)
    size = Int(1)

    _color_series = Any(transient=True)

    def plot(self, experiment, **kwargs):
        """
        Plot a dot plot colored by gate membership.
        """
        # Pre-compute per-event colors (Series indexed like experiment.data)
        self._color_series = pd.Series(_UNGATED_COLOR,
                                       index=experiment.data.index,
                                       dtype=object)
        for i, pop in enumerate(self.populations):
            if pop in experiment.conditions and experiment.conditions[pop].dtype.kind == 'b':
                mask = experiment.data[pop].astype(bool)
                self._color_series[mask] = _OKABE_ITO[i % len(_OKABE_ITO)]

        super().plot(experiment, **kwargs)

    def _grid_plot(self, experiment, grid, cmap, **kwargs):

        lim = kwargs.pop('lim')
        xlim = lim[self.xchannel]
        ylim = lim[self.ychannel]

        scale = kwargs.pop('scale')
        xscale = scale[self.xchannel]
        yscale = scale[self.ychannel]

        alpha = self.alpha
        s = self.size

        legend_handles = {}
        for i, pop in enumerate(self.populations):
            if pop in experiment.conditions and experiment.conditions[pop].dtype.kind == 'b':
                from matplotlib.lines import Line2D
                legend_handles[pop] = Line2D([0], [0], marker='o', color='w',
                                             markerfacecolor=_OKABE_ITO[i % len(_OKABE_ITO)],
                                             markersize=6, label=pop)

        grid.map(_dot_plot,
                 self.xchannel,
                 self.ychannel,
                 color_series=self._color_series,
                 alpha=alpha,
                 s=s)

        # add legend manually
        if legend_handles:
            ax = plt.gca()
            from matplotlib.lines import Line2D
            ungated = Line2D([0], [0], marker='o', color='w',
                             markerfacecolor=_UNGATED_COLOR,
                             markersize=6, label='Ungated')
            handles = [ungated] + list(legend_handles.values())
            ax.legend(handles=handles, frameon=True)

        return dict(xlim=xlim,
                    xscale=xscale,
                    ylim=ylim,
                    yscale=yscale)


def _dot_plot(x, y, color_series=None, alpha=0.3, s=1, **kwargs):
    ax = plt.gca()
    kwargs.pop('color', None)
    kwargs.pop('label', None)

    if color_series is None:
        ax.scatter(x.to_numpy(), y.to_numpy(),
                   s=s, alpha=alpha, color=_UNGATED_COLOR,
                   linewidths=0, rasterized=True)
        return ax

    # Reindex to current facet subset
    c = color_series.reindex(x.index)

    unique_colors = c.unique()
    for col in unique_colors:
        mask = (c == col).to_numpy()
        ax.scatter(x.to_numpy()[mask], y.to_numpy()[mask],
                   s=s, alpha=alpha, color=col,
                   linewidths=0, rasterized=True)

    return ax


util.expand_class_attributes(DotPlotView)
util.expand_method_parameters(DotPlotView, DotPlotView.plot)
