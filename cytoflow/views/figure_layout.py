#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.views.figure_layout
-----------------------------

Multi-panel figure layout: arrange several plots in a grid.

`FigureLayoutView` -- the view
`FigurePanel`      -- specification for one panel
"""

import matplotlib.pyplot as plt

from traits.api import HasStrictTraits, Str, Int, List, Any, Constant, provides

import cytoflow.utility as util
from .i_view import IView

# Maps short type names to view class constructors
# populated lazily to avoid circular imports
_VIEW_MAP = None

def _get_view_map():
    global _VIEW_MAP
    if _VIEW_MAP is None:
        from cytoflow.views import (HistogramView, ScatterplotView,
                                     DensityView, Kde1DView, Kde2DView,
                                     Histogram2DView, BarChartView)
        _VIEW_MAP = {
            'Histogram':    HistogramView,
            'Scatterplot':  ScatterplotView,
            'Density':      DensityView,
            'KDE 1D':       Kde1DView,
            'KDE 2D':       Kde2DView,
            'Histogram 2D': Histogram2DView,
            'Bar Chart':    BarChartView,
        }
    return _VIEW_MAP


class FigurePanel(HasStrictTraits):
    """
    Specification for a single panel in a `FigureLayoutView`.

    Attributes
    ----------
    view_type : Str
        The kind of plot: 'Histogram', 'Scatterplot', 'Density', 'KDE 1D',
        'KDE 2D', 'Histogram 2D', or 'Bar Chart'.
    channel_x : Str
        X-axis channel (required for 2D views, also used as the channel for 1D
        views when ``channel`` is not set).
    channel_y : Str
        Y-axis channel (2D views only).
    scale_x : Str
        X-axis scale ('linear', 'log', 'logicle').
    scale_y : Str
        Y-axis scale ('linear', 'log', 'logicle').
    huefacet : Str
        Color facet condition.
    subset : Str
        Subset expression (passed directly to the view).
    title : Str
        Optional panel title override.
    """
    view_type = Str('Histogram')
    channel_x = Str
    channel_y = Str
    scale_x   = Str('logicle')
    scale_y   = Str('logicle')
    huefacet  = Str
    subset    = Str
    title     = Str


@provides(IView)
class FigureLayoutView(HasStrictTraits):
    """
    Arrange multiple plots in a rows × cols grid.

    Each panel is a `FigurePanel` specifying which view to draw and its
    parameters.  All panels operate on the same experiment.

    Attributes
    ----------
    rows : Int
        Number of rows in the grid (default 2).
    cols : Int
        Number of columns in the grid (default 2).
    panels : List(FigurePanel)
        One entry per cell, row-major order.  Cells beyond ``rows × cols``
        are ignored; empty cells are left blank.
    title : Str
        Overall figure title.
    figsize_w : Int
        Width of the figure in inches per column (default 5).
    figsize_h : Int
        Height per row in inches (default 4).

    Examples
    --------
    >>> layout = flow.FigureLayoutView(rows=1, cols=2,
    ...     panels=[
    ...         flow.FigurePanel(view_type='Histogram', channel_x='FSC-A'),
    ...         flow.FigurePanel(view_type='Scatterplot',
    ...                          channel_x='FSC-A', channel_y='SSC-A'),
    ...     ])
    >>> layout.plot(ex)
    """

    id = Constant('cytoflow.view.figure_layout')
    friendly_id = Constant('Figure Layout')

    rows      = Int(2)
    cols      = Int(2)
    panels    = List(FigurePanel)
    title     = Str
    figsize_w = Int(5)
    figsize_h = Int(4)

    def plot(self, experiment, **kwargs):
        """
        Render the figure layout.

        Parameters
        ----------
        experiment : Experiment
            The experiment to plot.
        """
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                "Must specify an experiment")

        if self.rows < 1 or self.cols < 1:
            raise util.CytoflowViewError('rows',
                "rows and cols must be at least 1")

        vm = _get_view_map()

        fig, axes = plt.subplots(
            self.rows, self.cols,
            figsize=(self.cols * self.figsize_w, self.rows * self.figsize_h),
            squeeze=False)

        if self.title:
            fig.suptitle(self.title, fontsize=14, fontweight='bold')

        for idx in range(self.rows * self.cols):
            row = idx // self.cols
            col = idx % self.cols
            ax  = axes[row][col]

            if idx >= len(self.panels):
                ax.set_visible(False)
                continue

            panel = self.panels[idx]
            if not panel.view_type or panel.view_type not in vm:
                ax.set_visible(False)
                continue

            view_cls = vm[panel.view_type]
            view_kwargs = {}

            if panel.channel_x:
                # 1D views use 'channel'; 2D views use 'xchannel'/'ychannel'
                # Try both naming conventions
                try:
                    view = view_cls(channel=panel.channel_x,
                                    scale=panel.scale_x)
                except Exception:
                    view_kwargs['xchannel'] = panel.channel_x
                    view_kwargs['ychannel'] = panel.channel_y or panel.channel_x
                    view_kwargs['xscale'] = panel.scale_x
                    view_kwargs['yscale'] = panel.scale_y
                    view = view_cls(**view_kwargs)
                else:
                    view_kwargs = {}

            else:
                view = view_cls()

            if panel.huefacet and hasattr(view, 'huefacet'):
                view.huefacet = panel.huefacet
            if panel.channel_y and hasattr(view, 'ychannel'):
                view.ychannel = panel.channel_y
            if panel.channel_y and hasattr(view, 'channel_y'):
                view.channel_y = panel.channel_y
            if panel.scale_y and hasattr(view, 'yscale'):
                view.yscale = panel.scale_y
            if panel.subset and hasattr(view, 'subset'):
                view.subset = panel.subset

            # Plot into the axes for this panel
            # We temporarily override plt.gca() by drawing on `ax`
            try:
                # If the view supports an 'ax' kwarg, use it
                try:
                    view.plot(experiment, ax=ax)
                except TypeError:
                    # Fall back: set axes as current and plot
                    plt.sca(ax)
                    view.plot(experiment)
            except Exception as e:
                ax.text(0.5, 0.5, f"Error:\n{e}",
                        ha='center', va='center', transform=ax.transAxes,
                        fontsize=9, color='red', wrap=True)

            if panel.title:
                ax.set_title(panel.title)

        plt.tight_layout()
