#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.views.contour
----------------------

A 2-D density contour plot using fast 2D histogramming.

`ContourView` -- the `IView` class that makes the plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from traits.api import provides, Constant, Bool, Int

import cytoflow.utility as util

from .i_view import IView
from .base_views import Base2DView


@provides(IView)
class ContourView(Base2DView):
    """
    Plots a 2-D density contour map.  Faster than `Kde2DView` because it uses
    ``numpy.histogram2d`` instead of kernel density estimation, making it
    suitable for large datasets (>50 k events).

    Attributes
    ----------
    filled : bool
        If ``True`` (default), draw filled contours.  If ``False``, draw
        contour lines only.

    show_scatter : bool
        If ``True``, draw individual events as a faint scatter underneath the
        contours.  Default ``False``.

    n_levels : int
        Number of contour levels to draw.  Default 8.

    gridsize : int
        Resolution of the underlying 2-D histogram (pixels per axis).
        Default 100.

    Examples
    --------

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

        >>> flow.ContourView(xchannel = 'V2-A',
        ...                  xscale = 'log',
        ...                  ychannel = 'Y2-A',
        ...                  yscale = 'log').plot(ex)
    """

    id = Constant('cytoflow.view.contour')
    friend_id = Constant("Contour Plot")

    filled = Bool(True)
    show_scatter = Bool(False)
    n_levels = Int(8)
    gridsize = Int(100)

    def plot(self, experiment, **kwargs):
        """
        Plot a faceted 2-D density contour map.

        Parameters
        ----------
        alpha : float (default = 0.8)
            Alpha of the contour fill or lines.

        scatter_alpha : float (default = 0.1)
            Alpha of the background scatter (only used when
            ``show_scatter = True``).

        scatter_size : float (default = 1)
            Dot size of the background scatter.

        Notes
        -----
        Other ``kwargs`` are passed to ``matplotlib.axes.Axes.contourf`` /
        ``matplotlib.axes.Axes.contour``.
        """
        super().plot(experiment, **kwargs)

    def _grid_plot(self, experiment, grid, cmap, **kwargs):

        kwargs.setdefault('alpha', 0.8)

        lim = kwargs.pop('lim')
        xlim = lim[self.xchannel]
        ylim = lim[self.ychannel]

        scale = kwargs.pop('scale')
        xscale = scale[self.xchannel]
        yscale = scale[self.ychannel]

        scatter_alpha = kwargs.pop('scatter_alpha', 0.1)
        scatter_size = kwargs.pop('scatter_size', 1)

        grid.map(_contourplot,
                 self.xchannel,
                 self.ychannel,
                 xscale=xscale,
                 yscale=yscale,
                 xlim=xlim,
                 ylim=ylim,
                 filled=self.filled,
                 show_scatter=self.show_scatter,
                 n_levels=self.n_levels,
                 gridsize=self.gridsize,
                 scatter_alpha=scatter_alpha,
                 scatter_size=scatter_size,
                 **kwargs)

        return dict(xlim=xlim,
                    xscale=xscale,
                    ylim=ylim,
                    yscale=yscale)


def _contourplot(x, y, xscale=None, yscale=None,
                 xlim=None, ylim=None,
                 filled=True, show_scatter=False,
                 n_levels=8, gridsize=100,
                 scatter_alpha=0.1, scatter_size=1,
                 **kwargs):
    ax = plt.gca()
    color = kwargs.pop('color', None)
    alpha = kwargs.pop('alpha', 0.8)
    kwargs.pop('label', None)

    xs = xscale(x)
    ys = yscale(y)

    nan_mask = np.isnan(xs) | np.isnan(ys)
    xs = xs[~nan_mask]
    ys = ys[~nan_mask]

    if len(xs) < 10:
        return ax

    xl = (xscale(xlim[0]), xscale(xlim[1])) if xlim else (xs.min(), xs.max())
    yl = (yscale(ylim[0]), yscale(ylim[1])) if ylim else (ys.min(), ys.max())

    H, xedges, yedges = np.histogram2d(xs, ys, bins=gridsize, range=[xl, yl])
    H = H.T  # transpose so rows = y, cols = x

    # Gaussian smooth for nicer contours
    try:
        from scipy.ndimage import gaussian_filter
        H = gaussian_filter(H, sigma=1.5)
    except ImportError:
        pass

    # Build grid in original (data) coordinates
    xc = xscale.inverse((xedges[:-1] + xedges[1:]) / 2.0)
    yc = yscale.inverse((yedges[:-1] + yedges[1:]) / 2.0)
    XX, YY = np.meshgrid(xc, yc)

    if show_scatter:
        ax.scatter(xscale.inverse(xs), yscale.inverse(ys),
                   s=scatter_size, alpha=scatter_alpha,
                   color='#808080', linewidths=0, rasterized=True, zorder=1)

    contour_kwargs = dict(levels=n_levels, alpha=alpha, zorder=2)
    if color is not None:
        contour_kwargs['colors'] = [color]

    if filled:
        ax.contourf(XX, YY, H, **contour_kwargs)
    else:
        ax.contour(XX, YY, H, **contour_kwargs)

    return ax


util.expand_class_attributes(ContourView)
util.expand_method_parameters(ContourView, ContourView.plot)
