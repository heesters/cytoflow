#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.views.overlay_histogram
---------------------------------

Overlay multiple gated populations on one histogram panel.

`OverlayHistogramView` -- the `IView` class that makes the plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted

from traits.api import provides, Constant, Str, List, Float, Int, Bool

import cytoflow.utility as util
from .i_view import IView
from .base_views import BaseDataView


# Okabe-Ito colorblind-safe palette
_OKABE_ITO = [
    '#E69F00', '#56B4E9', '#009E73', '#F0E442',
    '#0072B2', '#D55E00', '#CC79A7', '#000000',
]


@provides(IView)
class OverlayHistogramView(BaseDataView):
    """
    Overlay histograms of multiple gated populations on one panel.

    Each ``population`` is a boolean condition name (eg a gate).  When set,
    only the events where that condition is ``True`` are plotted; they receive
    a unique color from the Okabe-Ito palette.  An "All events" trace is
    always added in light grey.

    Attributes
    ----------
    channel : str
        The channel to histogram.

    scale : {'linear', 'log', 'logicle'}
        Axis scale.

    populations : list of str
        Boolean conditions to overlay.  Each ``True`` subset is one trace.

    alpha : float
        Transparency of each histogram trace (default 0.5).

    bins : int
        Number of histogram bins (default 200).

    density : bool
        Normalise each trace to a probability density (default ``True``).
        This makes populations of different sizes visually comparable.

    Examples
    --------

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> # (after gating with Singlets and CD4_pos)
        >>> flow.OverlayHistogramView(channel = 'BL1-A',
        ...                           scale = 'logicle',
        ...                           populations = ['Singlets', 'CD4_pos']).plot(ex)
    """

    id = Constant('cytoflow.view.overlay_histogram')
    friend_id = Constant("Overlay Histogram")

    channel = Str
    scale = util.ScaleEnum
    populations = List(Str)
    alpha = Float(0.5)
    bins = Int(200)
    density = Bool(True)

    def plot(self, experiment, **kwargs):
        """
        Plot overlapping histogram traces.

        Parameters
        ----------
        title : str
        xlabel : str
        ylabel : str
        """
        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")

        if not self.channel:
            raise util.CytoflowViewError('channel', "Must specify a channel")

        if self.channel not in experiment.data:
            raise util.CytoflowViewError('channel',
                                         "Channel '{}' not in the experiment"
                                         .format(self.channel))

        for pop in self.populations:
            if pop not in experiment.conditions:
                raise util.CytoflowViewError('populations',
                                             "Population '{}' not in experiment"
                                             .format(pop))
            if experiment.conditions[pop].dtype.kind != 'b':
                raise util.CytoflowViewError('populations',
                                             "Population '{}' must be a boolean condition"
                                             .format(pop))

        title = kwargs.pop('title', '')
        xlabel = kwargs.pop('xlabel', util.channel_label(experiment, self.channel))
        ylabel = kwargs.pop('ylabel', 'Density' if self.density else 'Count')

        try:
            import cytoflow
            if cytoflow.RUNNING_IN_GUI:
                import seaborn as sns
                sns.set_style('whitegrid', rc={"xtick.bottom": True, "ytick.left": True})
        except Exception:
            pass

        scale_obj = util.scale_factory(self.scale, experiment, channel=self.channel)
        min_q = kwargs.pop('min_quantile', 0.001)
        max_q = kwargs.pop('max_quantile', 1.0)
        x_min = experiment[self.channel].quantile(min_q)
        x_max = experiment[self.channel].quantile(max_q)
        x_min = scale_obj.clip(x_min)
        x_max = scale_obj.clip(x_max)

        # apply subset
        if self.subset:
            try:
                experiment = experiment.query(self.subset)
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{}' isn't valid".format(self.subset)) from e

        # Build bin edges in scaled space, map back to data space
        bin_edges_scaled = np.linspace(scale_obj(x_min), scale_obj(x_max), self.bins + 1)
        bin_edges = scale_obj.inverse(bin_edges_scaled)
        bin_mids = scale_obj.inverse((bin_edges_scaled[:-1] + bin_edges_scaled[1:]) / 2.0)

        fig, ax = plt.subplots()

        # always plot all events in grey
        x_all = experiment.data[self.channel].values
        x_all_s = scale_obj(x_all)
        valid = ~np.isnan(x_all_s)
        h, _ = np.histogram(x_all_s[valid], bins=bin_edges_scaled)
        if self.density and h.sum() > 0:
            h = h / (h.sum() * np.diff(bin_edges_scaled))
        ax.fill_between(bin_mids, h, alpha=self.alpha * 0.5,
                        color='#aaaaaa', step='mid', label='All events')

        # each population
        for i, pop in enumerate(self.populations):
            color = _OKABE_ITO[i % len(_OKABE_ITO)]
            mask = experiment.data[pop].astype(bool).values
            x_pop_s = x_all_s[mask]
            valid_pop = ~np.isnan(x_pop_s)
            h_pop, _ = np.histogram(x_pop_s[valid_pop], bins=bin_edges_scaled)
            if self.density and h_pop.sum() > 0:
                h_pop = h_pop / (h_pop.sum() * np.diff(bin_edges_scaled))
            ax.fill_between(bin_mids, h_pop, alpha=self.alpha,
                            color=color, step='mid', label=pop)

        ax.set_xscale(scale_obj.name, **scale_obj.get_mpl_params(ax.get_xaxis()))
        ax.set_xlim(x_min, x_max)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title)
        ax.legend(frameon=True)


util.expand_class_attributes(OverlayHistogramView)
util.expand_method_parameters(OverlayHistogramView, OverlayHistogramView.plot)
