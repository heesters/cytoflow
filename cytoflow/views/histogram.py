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
cytoflow.views.histogram
------------------------

Plots a histogram.

`HistogramView` -- the `IView` class that makes the plot.
"""

from traits.api import Constant, provides, Str
import matplotlib.pyplot as plt

import numpy as np
import math
import bottleneck
from natsort import natsorted

import cytoflow.utility as util
from .i_view import IView
from .base_views import Base1DView

@provides(IView)
class HistogramView(Base1DView):
    """
    Plots a one-channel histogram

    Attributes
    ----------
    overlay_condition : str
        If set, all values of this condition are overlaid on **one** axes
        instead of being split into separate facets.  Each value gets its
        own color; a legend is added automatically.  ``xfacet``, ``yfacet``
        and ``huefacet`` are ignored while this is set.

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
        
    Plot a histogram
    
    .. plot::
        :context: close-figs
    
        >>> flow.HistogramView(channel = 'Y2-A',
        ...                    scale = 'log',
        ...                    huefacet = 'Dox').plot(ex)
    """
    
    # traits
    id = Constant("cytoflow.view.histogram")
    friendly_id = Constant("Histogram")

    overlay_condition = Str("")

    def plot(self, experiment, **kwargs):
        """
        Plot a faceted histogram view of a channel
        
        Parameters
        ----------
        num_bins : int
            The number of bins to plot in the histogram.  Clipped to [100, 1000]
            
        histtype : {'stepfilled', 'step', 'bar'}
            The type of histogram to draw.  ``stepfilled`` is the default, which
            is a line plot with a color filled under the curve.
            
        density: bool
            If `True`, re-scale the histogram to form a probability density
            function, so the area under the histogram is 1.
            
        orientation : {'horizontal', 'vertical'}
            The orientation of the histogram.  ``horizontal`` gives a histogram
            with the intensity on the Y axis and the count on the X axis;
            default is ``vertical``.
        
        linewidth : float
            The width of the histogram line (in points)
                    
        linestyle : ['-' | '--' | '-.' | ':' | "None"]
            The style of the line to plot
            
        alpha : float (default = 0.5)
            The alpha blending value, between 0 (transparent) and 1 (opaque).
            
        normalize_to_mode : bool (default = False)
            Scale each histogram so its tallest bin is 100% ("% of max", 
            like FlowJo's "normalize to mode"). Useful for comparing the 
            *shape* of overlaid histograms with different numbers of events.
            
        offset : float (default = 0)
            Stagger overlaid histograms (one per `huefacet` value) vertically,
            by this fraction of the tallest histogram, so they don't hide 
            each other.  Eg, 0.5.  Usually combined with 
            ``normalize_to_mode``.
    
            
        Notes
        -----
        Other ``kwargs`` are passed to `matplotlib.pyplot.hist <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.hist.html>`_

        
        """
        if kwargs.get('normalize_to_mode', False):
            ylabel = '% of max'
        else:
            ylabel = 'Density' if kwargs.get('density', False) else 'Count'
        if kwargs.get('orientation', 'vertical') == 'vertical':
            kwargs.setdefault('xlabel', self.channel)
            kwargs.setdefault('ylabel', ylabel)
        else:  # flip axis labels
            kwargs.setdefault('xlabel', ylabel)
            kwargs.setdefault('ylabel', self.channel)

        if self.overlay_condition:
            self._overlay_plot(experiment, **kwargs)
        else:
            super().plot(experiment, **kwargs)

    def _overlay_plot(self, experiment, **kwargs):
        """Plot all values of overlay_condition on one axes."""
        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")
        if self.overlay_condition not in experiment.conditions:
            raise util.CytoflowViewError('overlay_condition',
                                         "Overlay condition '{}' not in experiment"
                                         .format(self.overlay_condition))

        import seaborn as sns
        import cytoflow
        if cytoflow.RUNNING_IN_GUI:
            sns.set_style('whitegrid', rc={"xtick.bottom": True, "ytick.left": True})
            sns.set_context('notebook')

        scale = util.scale_factory(self.scale, experiment, channel=self.channel)

        subset = kwargs.pop('subset', getattr(self, 'subset', ''))
        if subset:
            try:
                experiment = experiment.query(subset)
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{}' isn't valid".format(subset)) from e

        lim = kwargs.pop('lim', None)
        if lim is None:
            min_q = kwargs.pop('min_quantile', 0.001)
            max_q = kwargs.pop('max_quantile', 1.0)
            lim = (experiment[self.channel].quantile(min_q),
                   experiment[self.channel].quantile(max_q))
        else:
            kwargs.pop('min_quantile', None)
            kwargs.pop('max_quantile', None)

        # remove remaining unused base-class kwargs
        for k in ('scale', 'sns_style', 'sns_context', 'sharex', 'sharey',
                  'height', 'aspect', 'col_wrap', 'margin_titles', 'despine',
                  'legend_loc', 'huelabel'):
            kwargs.pop(k, None)

        legend = kwargs.pop('legend', True)
        xlabel = kwargs.pop('xlabel', self.channel)
        ylabel = kwargs.pop('ylabel', 'Density' if kwargs.get('density') else 'Count')
        title = kwargs.pop('title', None)

        normalize_to_mode = kwargs.pop('normalize_to_mode', False)
        offset = kwargs.pop('offset', 0) or 0

        scaled_data = scale(experiment[self.channel])
        num_bins = kwargs.pop('num_bins', util.num_hist_bins(scaled_data))
        num_bins = util.num_hist_bins(scaled_data) if num_bins is None else num_bins
        num_bins = max(min(num_bins, 1000), 100)

        xmin = bottleneck.nanmin(scaled_data)
        xmax = bottleneck.nanmax(scaled_data)
        bins = scale.inverse(np.linspace(xmin, xmax, num=int(num_bins), endpoint=True))

        kwargs.setdefault('histtype', 'step')
        kwargs.setdefault('linewidth', 2)
        kwargs.setdefault('density', True)
        kwargs.pop('alpha', None)  # step histograms don't use alpha well

        values = natsorted(experiment[self.overlay_condition].unique())
        palette = sns.color_palette('tab10', len(values))

        fig, ax = plt.subplots()
        tops = []
        for i, val in enumerate(values):
            subset_data = experiment[experiment[self.overlay_condition] == val][self.channel]
            subset_data = subset_data[(subset_data > bins[0]) & (subset_data < bins[-1])]
            if len(subset_data) == 0:
                continue
            if normalize_to_mode:
                counts, _ = np.histogram(subset_data, bins=bins)
                counts = counts.astype(float)
                if counts.max() > 0:
                    counts = counts / counts.max() * 100.0
                n, _, _ = ax.hist(bins[:-1], bins=bins, weights=counts,
                                  color=palette[i], label=str(val),
                                  bottom=i * offset * 100.0 if offset else None,
                                  **kwargs)
            else:
                n, _, _ = ax.hist(subset_data, bins=bins, color=palette[i],
                                  label=str(val), **kwargs)
            tops.append(max(n))

        ax.set_xscale(scale.name, **scale.get_mpl_params(ax.get_xaxis()))
        ax.set_xlim(lim)
        if tops:
            ax.set_ylim(0, 1.05 * max(tops))
        ax.set_xlabel(util.channel_label(experiment, xlabel))
        ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title)
        if legend:
            ax.legend(title=util.channel_label(experiment, self.overlay_condition),
                      frameon=True)

    def _grid_plot(self, experiment, grid, cmap, **kwargs):
                        
        kwargs.setdefault('histtype', 'stepfilled')
        kwargs.setdefault('alpha', 0.5)
        kwargs.setdefault('antialiased', True)

        # estimate a "good" number of bins; see cytoflow.utility.num_hist_bins
        # for a reference.
        scale = kwargs.pop('scale')[self.channel]
        lim = kwargs.pop('lim')[self.channel]
        
        scaled_data = scale(experiment[self.channel])
        num_bins = kwargs.pop('num_bins', util.num_hist_bins(scaled_data))
        num_bins = util.num_hist_bins(scaled_data) if num_bins is None else num_bins
        
        # clip num_bins to (100, 1000)
        num_bins = max(min(num_bins, 1000), 100)
        
        if (self.huefacet 
            and "bins" in experiment.metadata[self.huefacet]
            and experiment.metadata[self.huefacet]["bin_scale"] == self.scale):
            
            # if we color facet by the result of a BinningOp and we don't
            # match the BinningOp bins with the histogram bins, we get
            # gnarly aliasing.
            
            # each color gets at least one bin.  however, if the estimated
            # number of bins for the histogram is much larger than the
            # number of colors, sub-divide each color into multiple bins.
            bins = experiment.metadata[self.huefacet]["bins"]
            scaled_bins = scale(bins)
     
            num_hues = len(experiment[self.huefacet].unique())
            bins_per_hue = math.floor(num_bins / num_hues)

            if bins_per_hue == 1:
                new_bins = scaled_bins
            else:
                new_bins = []
                for idx in range(1, len(scaled_bins)):
                    new_bins = np.append(new_bins, 
                                         np.linspace(scaled_bins[idx - 1],
                                                     scaled_bins[idx],
                                                     bins_per_hue + 1,
                                                     endpoint = False))

            bins = scale.inverse(new_bins)
        else:
            xmin = bottleneck.nanmin(scaled_data)
            xmax = bottleneck.nanmax(scaled_data)
            bins = scale.inverse(np.linspace(xmin, xmax, num=int(num_bins), endpoint = True))
                    
        kwargs.setdefault('bins', bins) 
        kwargs.setdefault('orientation', 'vertical')
        
        if ('linewidth' not in kwargs) or ('linewidth' in kwargs and kwargs['linewidth'] is None):
            kwargs['linewidth'] = 0 if kwargs['histtype'] == "stepfilled" else 2
        
        # if we have a hue facet, the y scaling is frequently wrong.  this
        # will capture the maximum bin count of each call to plt.hist, so 
        # we don't have to compute the histogram multiple times
        count_max = []
        
        normalize_to_mode = kwargs.pop('normalize_to_mode', False)
        offset = kwargs.pop('offset', 0) or 0
        if normalize_to_mode and kwargs.get('density'):
            raise util.CytoflowViewError(None,
                                         "Can't set both 'density' and 'normalize_to_mode'")
        if offset < 0:
            raise util.CytoflowViewError(None, "'offset' must be >= 0")
            
        # how many histograms have been drawn on each axes (for the offset)
        n_drawn = {}
        tops = []
        
        def hist_lims(*args, **kwargs):
            # there's some bug in the above code where we get data that isn't
            # in the range of `bins`, which makes hist() puke.  so get rid
            # of it.
            
            bins = kwargs.get('bins')
            new_args = []
            for x in args:
                x = x[x > bins[0]]
                x = x[x < bins[-1]]
                new_args.append(x)
                
            if normalize_to_mode or offset:
                counts, _ = np.histogram(new_args[0], bins = kwargs["bins"])
                counts = counts.astype(float)
                if normalize_to_mode and counts.max() > 0:
                    counts = counts / counts.max() * 100.0
                kwargs["weights"] = counts
                if offset:
                    ax = plt.gca()
                    i = n_drawn.get(ax, 0)
                    n_drawn[ax] = i + 1
                    kwargs["bottom"] = i * offset * (100.0 if normalize_to_mode else counts_max_all)
                    # draw the later (higher) ones behind
                    kwargs["zorder"] = 10 - i
                    # the count axis doesn't mean much any more
                    if kwargs.get("orientation", "vertical") == "vertical":
                        ax.tick_params(axis = 'y', labelleft = False)
                    else:
                        ax.tick_params(axis = 'x', labelbottom = False)
                n, _, _ = plt.hist(kwargs["bins"][:-1], **kwargs)
                tops.append(max(n) + (kwargs.get("bottom") or 0))
            elif scale.name != "linear" and kwargs.get("density"):
                kwargs["density"] = False
                counts, _ = np.histogram(new_args, bins=kwargs["bins"])
                kwargs["weights"] = counts / np.sum(counts)
                n, _, _ = plt.hist(kwargs["bins"][:-1], **kwargs)
            else:
                n, _, _ = plt.hist(*new_args, **kwargs)

            count_max.append(max(n))
            
        # with an un-normalized offset, stagger by a fraction of the 
        # largest bin count across all the histograms
        counts_max_all = 0
        if offset and not normalize_to_mode:
            for _, data_ijk in grid.facet_data():
                if self.huefacet:
                    groups = [g for _, g in data_ijk.groupby(self.huefacet, observed = True)]
                else:
                    groups = [data_ijk]
                for g in groups:
                    if len(g):
                        c, _ = np.histogram(g[self.channel], bins = kwargs['bins'])
                        counts_max_all = max(counts_max_all, c.max())
                    
        grid.map(hist_lims, self.channel, **kwargs)
        
        top = max(tops) if tops else max(count_max)
        
        ret = {}
        if kwargs['orientation'] == 'vertical':
            ret['xscale'] = scale
            ret['xlim'] = lim
            ret['ylim'] = (0, 1.05 * top)
        else:
            ret['yscale'] = scale
            ret['ylim'] = lim
            ret['xlim'] = (0, 1.05 * top)
            
        return ret

util.expand_class_attributes(HistogramView)
util.expand_method_parameters(HistogramView, HistogramView.plot)
