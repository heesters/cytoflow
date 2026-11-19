#!/usr/bin/env python3.8
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2022
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
cytoflow.views.scatterplot
--------------------------

A 2-d scatterplot.

`ScatterplotView` -- the `IView` class that makes the plot.
"""

import matplotlib.pyplot as plt
import seaborn as sns
from warnings import warn

from traits.api import provides, Constant, Str, Any
from natsort import natsorted
import numpy as np

import cytoflow.utility as util

from .i_view import IView
from .base_views import Base2DView
from .highlight import highlight_mask, draw_highlight, pop_highlight_kwargs, BACKGROUND_COLOR

try:
    import datashader as _ds
    import datashader.transfer_functions as _tf
    _DATASHADER_AVAILABLE = True
except ImportError:
    _DATASHADER_AVAILABLE = False

@provides(IView)
class ScatterplotView(Base2DView):
    """
    Plots a 2-d scatterplot.  
    
    Attributes
    ----------
    
    huechannel : Str
        If set, color the points using a normed color scale. The norm function
        is set by `huescale`, and the color palette can be changed by passing
        the ``palette`` parameter to `plot`.  
        
    highlight : Str
        "Back-gating": draw these events on top of the others, in a 
        contrasting color, and draw the rest of the events in grey.  Either
        the name of a True/False condition (eg, a gate) or an expression,
        eg ``"Quad == 'Quad_2'"``.
    
    Examples
    --------
    
    Make a little data set.
    
    .. plot::
        :context: close-figs
            
        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()
        
    Plot a density plot
    
    .. plot::
        :context: close-figs
    
        >>> flow.ScatterplotView(xchannel = 'V2-A',
        ...                      xscale = 'log',
        ...                      ychannel = 'Y2-A',
        ...                      yscale = 'log',
        ...                      huefacet = 'Dox').plot(ex)
        
    """
    
    id = Constant('cytoflow.view.scatterplot')
    friend_id = Constant("Scatter Plot")
    huechannel = Str
    highlight = Str
    overlay_condition = Str("")
    
    _highlight_mask = Any(transient = True)
    
    def plot(self, experiment, **kwargs):
        """
        Plot a faceted scatter plot view of a channel
        
        Parameters
        ----------
        
        alpha : float (default = 0.25)
            The alpha blending value, between 0 (transparent) and 1 (opaque).
            
        s : int (default = 2)
            The size in points^2.
            
        marker : a matplotlib marker style, usually a string
            Specfies the glyph to draw for each point on the scatterplot.
            See `matplotlib.markers <http://matplotlib.org/api/markers_api.html#module-matplotlib.markers>`_ for examples.  Default: 'o'
            
        highlight_color, highlight_size, highlight_alpha
            The color, size (in points^2) and alpha of the highlighted events,
            if `highlight` is set.
            
        
        Notes
        -----
        Other ``kwargs`` are passed to `matplotlib.pyplot.scatter <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.scatter.html>`_

        """
        
        if self.huechannel:
            if self.huefacet:
                raise util.CytoflowViewError('huefacet',
                                             "Can't set both 'huefacet' and 'huechannel'!")
        
            if self.huechannel not in experiment.channels:
                raise util.CytoflowViewError('huechannel',
                                             f"Can't find channel {self.huechannel} in the experiment.")
                            
            kwargs.setdefault('palette', 'viridis_r')
            kwargs.setdefault('huelabel', self.huechannel)
            
        self._highlight_mask = highlight_mask(experiment, self.highlight) \
                               if self.highlight else None

        if self.overlay_condition:
            self._overlay_plot(experiment, **kwargs)
        else:
            super().plot(experiment, **kwargs)

    def _overlay_plot(self, experiment, **kwargs):
        """Plot all values of overlay_condition on one axes with per-condition colors."""
        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")
        if self.overlay_condition not in experiment.conditions:
            raise util.CytoflowViewError('overlay_condition',
                                         "Overlay condition '{}' not in experiment"
                                         .format(self.overlay_condition))

        xscale = util.scale_factory(self.xscale, experiment, channel=self.xchannel)
        yscale = util.scale_factory(self.yscale, experiment, channel=self.ychannel)

        subset = kwargs.pop('subset', getattr(self, 'subset', ''))
        if subset:
            try:
                experiment = experiment.query(subset)
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{}' isn't valid".format(subset)) from e

        min_q = kwargs.pop('min_quantile', 0.001)
        max_q = kwargs.pop('max_quantile', 1.0)
        xlim = (experiment[self.xchannel].quantile(min_q),
                experiment[self.xchannel].quantile(max_q))
        ylim = (experiment[self.ychannel].quantile(min_q),
                experiment[self.ychannel].quantile(max_q))

        for k in ('scale', 'lim', 'sns_style', 'sns_context', 'sharex', 'sharey',
                  'height', 'aspect', 'col_wrap', 'margin_titles', 'despine',
                  'legend_loc', 'huelabel', 'palette'):
            kwargs.pop(k, None)

        legend = kwargs.pop('legend', True)
        xlabel = kwargs.pop('xlabel', self.xchannel)
        ylabel = kwargs.pop('ylabel', self.ychannel)
        title = kwargs.pop('title', None)

        alpha = kwargs.pop('alpha', 0.25)
        s = kwargs.pop('s', 2)
        marker = kwargs.pop('marker', 'o')

        values = natsorted(experiment[self.overlay_condition].unique())
        palette = sns.color_palette('tab10', len(values))

        fig, ax = plt.subplots()
        for i, val in enumerate(values):
            subset_data = experiment[experiment[self.overlay_condition] == val]
            ax.scatter(subset_data[self.xchannel], subset_data[self.ychannel],
                       c=[palette[i]], alpha=alpha, s=s, marker=marker,
                       label=str(val))

        ax.set_xscale(xscale.name, **xscale.get_mpl_params(ax.get_xaxis()))
        ax.set_yscale(yscale.name, **yscale.get_mpl_params(ax.get_yaxis()))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(util.channel_label(experiment, xlabel))
        ax.set_ylabel(util.channel_label(experiment, ylabel))
        if title:
            ax.set_title(title)
        if legend:
            ax.legend(title=util.channel_label(experiment, self.overlay_condition),
                      frameon=True, markerscale=5)
        
    def _grid_plot(self, experiment, grid, cmap, **kwargs):

        kwargs.setdefault('alpha', 0.25)
        kwargs.setdefault('s', 2)
        kwargs.setdefault('marker', 'o')
        kwargs.setdefault('antialiased', True)

        lim = kwargs.pop('lim')
        xlim = lim[self.xchannel]
        ylim = lim[self.ychannel]

        scale = kwargs.pop('scale')
        xscale = scale[self.xchannel]
        yscale = scale[self.ychannel]

        highlight_kwargs = pop_highlight_kwargs(kwargs)
        if self._highlight_mask is not None and not self.huefacet and not self.huechannel:
            # draw everything else in grey
            kwargs['color'] = BACKGROUND_COLOR

        # datashader fast path: use for simple density scatter (no per-point
        # hue channel, no back-gating highlight overlay).  Rasterizes all events
        # to a fixed-size RGBA image in O(pixels) instead of O(events).
        _use_datashader = (
            _DATASHADER_AVAILABLE and
            not self.huechannel and
            self._highlight_mask is None
        )

        if self.huechannel:
            if isinstance(cmap, list):
                raise util.CytoflowViewError('palette',
                                             "Must use a continuous palette with huechannel set!")
            kwargs['cmap'] = cmap
            hue_scale = util.scale_factory(self.huescale,
                                           experiment,
                                           channel = self.huechannel)
            norm = hue_scale.norm()
            kwargs['norm'] = norm
            grid.map(_scatterplot, self.xchannel, self.ychannel, self.huechannel, **kwargs)

        elif _use_datashader:
            norm = None
            # Compute limits in transformed (linear display) space
            xlim_t = (_safe_transform(xscale, xlim[0]),
                      _safe_transform(xscale, xlim[1]))
            ylim_t = (_safe_transform(yscale, ylim[0]),
                      _safe_transform(yscale, ylim[1]))
            # Strip scatter-specific kwargs that datashader doesn't need
            ds_kwargs = {k: v for k, v in kwargs.items()
                         if k not in ('alpha', 's', 'marker', 'antialiased', 'color')}
            grid.map(_datashader_scatter, self.xchannel, self.ychannel,
                     xscale=xscale, yscale=yscale,
                     xlim_t=xlim_t, ylim_t=ylim_t, **ds_kwargs)
            # Axes are kept LINEAR (pre-transformed data baked into the image)
            # so return None scales and transformed limits
            return dict(xlim=xlim_t, xscale=None, ylim=ylim_t, yscale=None, norm=None)

        else:
            norm = None
            grid.map(plt.scatter, self.xchannel, self.ychannel, **kwargs)

        if self._highlight_mask is not None:
            grid.map(draw_highlight, self.xchannel, self.ychannel,
                     mask = self._highlight_mask,
                     highlight_kwargs = highlight_kwargs)

        return dict(xlim = xlim,
                    xscale = xscale,
                    ylim = ylim,
                    yscale = yscale,
                    norm = norm)
    
    def _update_legend(self, legend):
        for lh in legend.legend_handles:
            lh.set_alpha(0.8)
            lh.set_sizes([10.0])
            
    
def _scatterplot(x, y, c, **kwargs):
    del kwargs['color']
    plt.scatter(x, y, c = c, **kwargs),


def _safe_transform(scale, value):
    """Apply scale transform to a scalar, returning the value unchanged on error."""
    try:
        result = scale(value)
        if np.isfinite(result):
            return float(result)
    except Exception:
        pass
    return float(value)


def _datashader_scatter(x, y, xscale=None, yscale=None,
                        xlim_t=None, ylim_t=None, **kwargs):
    """Rasterize x/y scatter via datashader in transformed coordinates."""
    import datashader as ds
    import datashader.transfer_functions as tf
    import pandas as pd
    from matplotlib.ticker import FuncFormatter

    ax = plt.gca()

    # Apply scale transforms
    x_t = np.asarray(xscale(x), dtype=float) if xscale else np.asarray(x, dtype=float)
    y_t = np.asarray(yscale(y), dtype=float) if yscale else np.asarray(y, dtype=float)

    # Remove NaN/Inf
    mask = np.isfinite(x_t) & np.isfinite(y_t)
    if not mask.any():
        return

    df = pd.DataFrame({'x': x_t[mask], 'y': y_t[mask]})
    x_min, x_max = xlim_t
    y_min, y_max = ylim_t

    if not (np.isfinite(x_min) and np.isfinite(x_max) and
            np.isfinite(y_min) and np.isfinite(y_max)):
        return
    if x_max <= x_min or y_max <= y_min:
        return

    canvas = ds.Canvas(plot_width=400, plot_height=400,
                       x_range=(x_min, x_max), y_range=(y_min, y_max))
    agg = canvas.points(df, 'x', 'y')

    # Use the facet color if supplied by seaborn, else a default blue
    color = kwargs.pop('color', '#1f77b4')
    if hasattr(color, '__len__') and not isinstance(color, str):
        try:
            import matplotlib.colors as mcolors
            color = mcolors.to_hex(color)
        except Exception:
            color = '#1f77b4'

    img = tf.shade(agg, cmap=['#f0f4ff', color], how='log')
    img_array = np.array(img.to_pil())   # shape (H, W, 4), row 0 at top (y_max)

    ax.imshow(img_array,
              extent=(x_min, x_max, y_min, y_max),
              origin='upper',
              aspect='auto',
              interpolation='nearest',
              zorder=0)

    # Tick formatters: show original data values on the linear axes
    if xscale is not None and xscale.name != 'linear':
        inv_x = xscale.inverse
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda v, _: _fmt_scale_tick(inv_x(v))))
    if yscale is not None and yscale.name != 'linear':
        inv_y = yscale.inverse
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda v, _: _fmt_scale_tick(inv_y(v))))


def _fmt_scale_tick(value):
    """Format a tick value for a scale-transformed axis."""
    try:
        v = float(value)
        if not np.isfinite(v):
            return ''
        if abs(v) >= 10000 or (abs(v) < 0.01 and v != 0):
            return f'{v:.2e}'
        if v == int(v):
            return str(int(v))
        return f'{v:.3g}'
    except Exception:
        return ''


util.expand_class_attributes(ScatterplotView)
util.expand_method_parameters(ScatterplotView, ScatterplotView.plot)
        
if __name__ == '__main__':
    import cytoflow as flow
    tube1 = flow.Tube(file = '../../cytoflow/tests/data/Plate01/RFP_Well_A3.fcs',
                      conditions = {"Dox" : 10.0})
    
    tube2 = flow.Tube(file = '../../cytoflow/tests/data/Plate01/CFP_Well_A4.fcs',
                      conditions = {"Dox" : 1.0})                      

    ex = flow.ImportOp(conditions = {"Dox" : "float"}, tubes = [tube1, tube2])
    
    thresh = flow.ThresholdOp()
    thresh.name = "Y2-A+"
    thresh.channel = 'Y2-A'
    thresh.threshold = 200.0

    ex2 = thresh.apply(ex)
    
    scatter = flow.ScatterplotView()
    scatter.name = "Scatter"
    scatter.xchannel = "FSC-A"
    scatter.ychannel = "SSC-A"
    scatter.xscale = "logicle"
    scatter.yscale = "logicle"
    scatter.huefacet = 'Dox'
    
    plt.ioff()
    scatter.plot(ex2)
    plt.show()