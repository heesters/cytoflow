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
cytoflow.operations.proliferation
---------------------------------

Dye-dilution proliferation analysis (CFSE, CellTrace Violet, ...).  
`proliferation` has:

`ProliferationOp` -- fits the generations and computes the proliferation
statistics.

`fit_generations` -- the model fit, as a function of a numpy array.

`proliferation_metrics` -- the standard indices, from the number of cells
in each generation.

`ProliferationView` -- plots the histogram with the fitted generations.
"""

import warnings

from traits.api import (HasStrictTraits, Str, List, Dict, Float, Int, Bool, 
                        Any, Constant, Instance, provides)

import numpy as np
import pandas as pd
import scipy.stats
import scipy.ndimage
import matplotlib.pyplot as plt

import cytoflow.views
import cytoflow.utility as util

from .i_operation import IOperation

LOG2 = np.log10(2.0)


def proliferation_metrics(counts):
    """
    The standard proliferation statistics, from the number of cells in each
    generation (``counts[g]`` cells have divided ``g`` times.)
    
    A cell that has divided ``g`` times is one of the :math:`2^g` 
    descendants of a cell of the starting population, so generation ``g`` 
    came from :math:`P_g = N_g / 2^g` *precursors*.  With 
    :math:`P = \\sum_g P_g` (the size of the starting population) and 
    :math:`P_d = \\sum_{g \\ge 1} P_g` (the precursors that divided at 
    least once), the statistics are (as defined by FlowJo's proliferation 
    platform):
    
    - **Percent_divided** = :math:`100\\,P_d / P` -- the percentage of the 
      *starting* cells that divided at least once (also called the 
      precursor frequency.)
    - **Division_index** = :math:`\\sum_g g\\,P_g / P` -- the average 
      number of divisions of *all* the starting cells (including the ones 
      that never divided.)
    - **Proliferation_index** = :math:`\\sum_g g\\,P_g / P_d` -- the 
      average number of divisions of the cells that *did* divide.
    - **Expansion_index** = :math:`\\sum_g N_g / P` -- the fold expansion
      of the whole culture.
    - **Replication_index** = :math:`\\sum_{g \\ge 1} N_g / P_d` -- the 
      fold expansion of the responding cells only.
      
    Returns
    -------
    dict
        The five statistics.  Statistics that divide by the number of 
        responding precursors are not-a-number if no cell divided.
    """
    n = np.asarray(counts, dtype = float)
    g = np.arange(len(n))
    p = n / 2.0 ** g
    p_total = p.sum()
    p_div = p[1:].sum()
    nan = float('nan')
    
    if p_total <= 0:
        return dict(Percent_divided = nan, Division_index = nan, 
                    Proliferation_index = nan, Expansion_index = nan,
                    Replication_index = nan)
        
    return dict(Percent_divided = float(100.0 * p_div / p_total),
                Division_index = float((g * p).sum() / p_total),
                Proliferation_index = float((g * p).sum() / p_div) if p_div > 0 else nan,
                Expansion_index = float(n.sum() / p_total),
                Replication_index = float(n[1:].sum() / p_div) if p_div > 0 else nan)


def _brightest_peak(x):
    """The position of the brightest substantial peak of ``x`` (log units)"""
    lo, hi = np.percentile(x, [0.5, 99.5])
    if hi <= lo:
        return float(np.median(x))
    hist, edges = np.histogram(x, bins = 256, range = (lo, hi))
    smooth = scipy.ndimage.gaussian_filter1d(hist.astype(float), 4)
    centers = (edges[:-1] + edges[1:]) / 2
    peaks = [i for i in range(1, len(smooth) - 1)
             if smooth[i] >= smooth[i - 1] and smooth[i] > smooth[i + 1]
             and smooth[i] >= 0.1 * smooth.max()]
    if not peaks:
        return float(centers[np.argmax(smooth)])
    return float(centers[peaks[-1]])


def fit_generations(x, num_generations = 8, undivided_peak = None, 
                    fit_spacing = False, max_iter = 500, tol = 1e-7):
    """
    Fit a dye-dilution histogram with one Gaussian per generation.
    
    Each division halves the dye, so on a log axis the generations are 
    evenly spaced peaks.  The model is a mixture of Gaussians in 
    :math:`\\log_{10}` of the signal, with means :math:`\\mu_0 - g\\,d` for
    generation :math:`g = 0, 1, ...`, one shared standard deviation, and a 
    free weight for each generation.  It is fitted by 
    expectation-maximization, which is deterministic and uses every event
    (no binning.)
    
    Parameters
    ----------
    x : numpy.ndarray
        The log10 dye signal of each event.
        
    num_generations : int
        The largest number of divisions to model.  Generations dimmer than
        the data (below its 0.5th percentile) are left out.
        
    undivided_peak : float
        :math:`\\mu_0`, the (log10) position of the undivided cells.  If 
        ``None``, it starts at the brightest substantial peak of the 
        histogram and is fitted.
        
    fit_spacing : bool
        If ``False``, :math:`d = \\log_{10} 2` exactly.  If ``True``, 
        :math:`d` is fitted (within 15% of that.)
        
    Returns
    -------
    dict
        ``mu0``, ``spacing``, ``sigma``, ``weights`` (one for each 
        generation that was modelled) and ``means``.
    """
    x = np.asarray(x, dtype = float)
    x = x[np.isfinite(x)]
    if len(x) < 50:
        raise util.CytoflowOpError(None,
                                   "Need at least 50 events to fit generations; got {}"
                                   .format(len(x)))
        
    fixed_mu0 = undivided_peak is not None
    mu0 = float(undivided_peak) if fixed_mu0 else _brightest_peak(x)
    d = LOG2
    
    # leave out generations that would be dimmer than the data
    floor = np.percentile(x, 0.5)
    n_gen = int(min(num_generations, max(0, np.floor((mu0 - floor) / d + 0.5)))) + 1
    g = np.arange(n_gen)
    
    sigma = min(0.1, d / 2.5)
    w = np.full(n_gen, 1.0 / n_gen)
    prev = -np.inf
    
    for _ in range(max_iter):
        means = mu0 - g * d
        # E step (in log space, for numerical stability)
        logp = (np.log(np.maximum(w, 1e-300))[None, :] 
                - 0.5 * ((x[:, None] - means[None, :]) / sigma) ** 2 
                - np.log(sigma))
        m = logp.max(axis = 1, keepdims = True)
        log_norm = m + np.log(np.exp(logp - m).sum(axis = 1, keepdims = True))
        r = np.exp(logp - log_norm)
        ll = float(log_norm.sum())
        
        # M step
        nk = r.sum(axis = 0)
        w = nk / nk.sum()
        if fit_spacing and n_gen > 1:
            # weighted least squares of x on generation number
            sw = nk.sum()
            g_mean = (nk * g).sum() / sw
            x_mean = (r * x[:, None]).sum() / sw
            sgg = (nk * (g - g_mean) ** 2).sum()
            if sgg > 0:
                sgx = (r * (g[None, :] - g_mean) * (x[:, None] - x_mean)).sum()
                d = float(np.clip(-sgx / sgg, 0.85 * LOG2, 1.15 * LOG2))
            if not fixed_mu0:
                mu0 = float(x_mean + d * g_mean)
        elif not fixed_mu0:
            mu0 = float((r * (x[:, None] + g[None, :] * d)).sum() / nk.sum())
        means = mu0 - g * d
        sigma = float(np.sqrt((r * (x[:, None] - means[None, :]) ** 2).sum() / nk.sum()))
        # keep the peaks from merging into one blob, or collapsing
        sigma = float(np.clip(sigma, 0.01, d / 1.5))
        
        if abs(ll - prev) < tol * max(1.0, abs(ll)):
            break
        prev = ll
        
    return dict(mu0 = mu0, spacing = d, sigma = sigma, weights = w, means = mu0 - g * d)


def _posteriors(x, fit):
    means, sigma, w = fit['means'], fit['sigma'], fit['weights']
    logp = (np.log(np.maximum(w, 1e-300))[None, :] 
            - 0.5 * ((x[:, None] - means[None, :]) / sigma) ** 2)
    return logp


@provides(IOperation)
class ProliferationOp(HasStrictTraits):
    """
    Dye-dilution proliferation analysis: find the generations in a CFSE or 
    CellTrace histogram, and compute the proliferation statistics.
    
    Every time a cell divides, its daughters get half of the dye.  On a log
    axis, the generations are therefore evenly spaced peaks.  `estimate`
    fits a Gaussian to each generation (see `fit_generations`); `apply` 
    adds:
    
    - a condition named `name` with each event's generation (0 = undivided,
      1 = divided once, ...), and
    - a statistic named `name` with, for each group of `by`: the 
      percentage of cells in each generation (``Generation_0``, ...), the 
      position of the undivided peak, and ``Percent_divided``, 
      ``Division_index``, ``Proliferation_index``, ``Expansion_index`` and
      ``Replication_index`` (see `proliferation_metrics` for their 
      definitions.)
      
    Gate on single, live cells first.
    
    Attributes
    ----------
    name : Str
        The operation name; names the new condition and statistic.
        
    channel : Str
        The channel with the dye.
        
    num_generations : Int (default = 8)
        The largest number of divisions to look for.
        
    by : List(Str)
        Fit separately for each combination of these conditions' values 
        (eg, each sample.)
        
    undivided_peak : Float
        The position of the undivided cells, in the channel's units.  If 
        it's not set, it is found from an unstimulated control (see
        `control_subset`), or else from the data: the brightest peak of 
        each group.
        
    control_subset : Str
        A subset with *undivided* cells -- eg, ``"Stim == 'none'"`` for an 
        unstimulated control.  Its peak sets the position of generation 0 
        for all the groups.  This is the most reliable way to find it: if 
        nearly all the cells in a sample have divided, the brightest peak 
        in that sample is *not* generation 0.
        
    fit_spacing : Bool (default = False)
        The generations are exactly a factor of 2 apart.  If ``True``, fit 
        the spacing instead (within 15%), which helps if the cytometer's
        amplifiers aren't quite linear.
    """
    
    id = Constant('cytoflow.operations.proliferation')
    friendly_id = Constant("Proliferation")
    
    name = Str
    channel = Str
    num_generations = Int(8)
    by = List(Str)
    undivided_peak = util.FloatOrNone(None)
    control_subset = Str
    fit_spacing = Bool(False)
    
    _fits = Dict(Any, Any, transient = True)       # group -> fit
    
    def _log_data(self, data):
        x = data[self.channel].values.astype(float)
        ok = np.isfinite(x) & (x > 0)
        return np.log10(x[ok]), ok
    
    def estimate(self, experiment, subset = None):
        """
        Fit the generations.
        
        Parameters
        ----------
        experiment : Experiment
        
        subset : Str
            Only use these events (eg, the live single cells.)
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")
        
        if not self.channel:
            raise util.CytoflowOpError('channel', "Must set a channel")
        
        if self.channel not in experiment.channels:
            raise util.CytoflowOpError('channel',
                                       "Channel {} isn't in the experiment".format(self.channel))
            
        if self.num_generations < 1:
            raise util.CytoflowOpError('num_generations', "num_generations must be at least 1")
        
        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "Condition {} isn't in the experiment".format(b))
                
        if self.undivided_peak is not None and not self.undivided_peak > 0:
            raise util.CytoflowOpError('undivided_peak', "undivided_peak must be > 0")
                
        self._fits = {}
        
        if subset:
            try:
                experiment = experiment.query(subset)
            except Exception as e:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' isn't valid".format(subset)) from e
            if len(experiment) == 0:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' returned no events".format(subset))
                
        # where are the undivided cells?
        mu0 = None
        if self.undivided_peak is not None:
            mu0 = float(np.log10(self.undivided_peak))
        elif self.control_subset:
            try:
                control = experiment.query(self.control_subset)
            except Exception as e:
                raise util.CytoflowOpError('control_subset',
                                           "Subset string '{0}' isn't valid"
                                           .format(self.control_subset)) from e
            x, _ = self._log_data(control.data)
            if len(x) < 50:
                raise util.CytoflowOpError('control_subset',
                                           "The control subset has fewer than 50 "
                                           "events with a positive signal")
            # fit the control itself, and use the position of its brightest
            # generation (more accurate than the raw histogram's peak)
            try:
                mu0 = fit_generations(x, num_generations = self.num_generations,
                                      fit_spacing = self.fit_spacing)['mu0']
            except util.CytoflowOpError as e:
                raise util.CytoflowOpError('control_subset', e.args[-1]) from e
            
        if self.by:
            groups = experiment.data.groupby(self.by, observed = True)
        else:
            groups = [((), experiment.data)]
            
        fits = {}
        for group, data in groups:
            key = group if isinstance(group, tuple) else (group,)
            x, ok = self._log_data(data)
            if (~ok).sum() > 0.05 * max(1, len(ok)):
                warnings.warn("{:.0f}% of the events in {} have a signal <= 0 and "
                              "were left out. Is the channel compensated correctly?"
                              .format(100.0 * (~ok).mean(), key if key else "the data"),
                              util.CytoflowOpWarning)
            try:
                fits[key] = fit_generations(x, 
                                            num_generations = self.num_generations,
                                            undivided_peak = mu0,
                                            fit_spacing = self.fit_spacing)
            except util.CytoflowOpError as e:
                raise util.CytoflowOpError('by', "Group {}: {}".format(key, e.args[-1])) from e
            
        # set atomically, for the GUI
        self._fits = fits
        
    def apply(self, experiment):
        """
        Assign each event to a generation, and compute the statistics.
        
        Returns
        -------
        Experiment
            A new experiment with the condition `name` (the event's 
            generation; events with a signal <= 0 are put in the last 
            generation) and the statistic `name`.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")
        
        if not self.name:
            raise util.CytoflowOpError('name', "Must specify a name")
        
        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores.")
            
        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "{} is in the experiment already!".format(self.name))
            
        if self.name in experiment.statistics:
            raise util.CytoflowOpError('name',
                                       "Experiment already has a statistic named {}".format(self.name))
            
        if not self._fits:
            raise util.CytoflowOpError(None, "No model. Did you forget to call estimate()?")
        
        if self.channel not in experiment.channels:
            raise util.CytoflowOpError('channel',
                                       "Channel {} isn't in the experiment".format(self.channel))
            
        max_gen = max(len(f['weights']) for f in self._fits.values())
        generation = np.zeros(len(experiment), dtype = int)
        
        if self.by:
            groups = experiment.data.groupby(self.by, observed = True).indices.items()
        else:
            groups = [((), np.arange(len(experiment)))]
            
        rows = []
        index = []
        for group, idx in groups:
            key = group if isinstance(group, tuple) else (group,)
            if key not in self._fits:
                # a group that wasn't in the subset used for estimate()
                continue
            fit = self._fits[key]
            n_gen = len(fit['weights'])
            x = experiment.data[self.channel].values[idx].astype(float)
            ok = np.isfinite(x) & (x > 0)
            gen = np.full(len(idx), n_gen - 1, dtype = int)
            if ok.any():
                gen[ok] = np.argmax(_posteriors(np.log10(x[ok]), fit), axis = 1)
            generation[idx] = gen
            
            # the model's expected number of cells in each generation: less
            # biased than counting the (hard) assignments when peaks overlap
            counts = np.zeros(max_gen)
            counts[:n_gen] = fit['weights'] * len(idx)
            row = {"Generation_{}".format(g) : 100.0 * counts[g] / counts.sum() 
                   for g in range(max_gen)}
            row.update(proliferation_metrics(counts))
            row["Undivided_peak"] = 10 ** fit['mu0']
            rows.append(row)
            index.append(key)
            
        if self.by:
            idx = pd.MultiIndex.from_tuples(index, names = self.by) if len(self.by) > 1 \
                  else pd.Index([k[0] for k in index], name = self.by[0])
        else:
            idx = pd.Index(["All"], name = "Group")
        stat = pd.DataFrame(rows, index = idx)
        
        new_experiment = experiment.clone(deep = False)
        new_experiment.add_condition(self.name, "int", generation)
        new_experiment.statistics[self.name] = stat
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment
    
    def default_view(self, **kwargs):
        """
        Returns a `ProliferationView`: the histogram of each group with the
        fitted generations.
        """
        v = ProliferationView(op = self)
        v.trait_set(**kwargs)
        return v
    

def plot_fit_grid(groups, by, draw, xlabel, title = None):
    """
    Small multiples: call ``draw(ax, key)`` for each group (at most 16.)
    """
    keys = list(groups)[:16]
    n = len(keys)
    ncols = 1 if n == 1 else (2 if n <= 4 else (3 if n <= 9 else 4))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, squeeze = False,
                             figsize = (4.2 * ncols + 0.8, 2.6 * nrows + 0.8))
    for i, key in enumerate(keys):
        ax = axes[i // ncols][i % ncols]
        draw(ax, key)
        if by and key:
            ax.set_title(", ".join("{} = {}".format(b, k) for b, k in zip(by, key)), 
                         fontsize = 'small')
        if i // ncols == nrows - 1 or i + ncols >= n:
            ax.set_xlabel(xlabel)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].set_axis_off()
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


@provides(cytoflow.views.IView)
class ProliferationView(HasStrictTraits):
    """
    Plots the dye histogram of each group, with the fitted generations.
    
    Attributes
    ----------
    op : Instance(ProliferationOp)
        The operation whose fit to plot.
        
    subset : Str
        Only plot these events.
    """
    
    id = Constant('cytoflow.view.proliferationview')
    friendly_id = Constant("Proliferation")
    
    op = Instance(ProliferationOp)
    subset = Str
    
    def plot(self, experiment, **kwargs):
        """
        Plot the histograms (one panel for each group; at most 16.)
        
        Parameters
        ----------
        title : str
            The figure title
        """
        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")
        
        op = self.op
        if op is None or not op._fits:
            raise util.CytoflowViewError('op', "No model to plot. Did you forget to call estimate()?")
        
        if op.channel not in experiment.channels:
            raise util.CytoflowViewError('op', "Channel {} isn't in the experiment".format(op.channel))
        
        if self.subset:
            try:
                experiment = experiment.query(self.subset)
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' isn't valid"
                                             .format(self.subset)) from e
            if len(experiment) == 0:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' returned no events"
                                             .format(self.subset))
                
        if op.by:
            data = {(k if isinstance(k, tuple) else (k,)) : d
                    for k, d in experiment.data.groupby(op.by, observed = True)}
        else:
            data = {() : experiment.data}
        keys = [k for k in op._fits if k in data]
        if not keys:
            raise util.CytoflowViewError('subset', "No fitted groups to plot")
            
        label = util.channel_label(experiment, op.channel)
        cmap = plt.get_cmap("viridis")
        
        def draw(ax, key):
            fit = op._fits[key]
            x, _ = op._log_data(data[key])
            lo = min(np.percentile(x, 0.2), fit['means'].min() - 3 * fit['sigma'])
            hi = max(np.percentile(x, 99.8), fit['means'].max() + 3 * fit['sigma'])
            bins = np.linspace(lo, hi, 150)
            ax.hist(10 ** x, bins = 10 ** bins, color = '0.8', density = False)
            grid = np.linspace(lo, hi, 600)
            scale = len(x) * (bins[1] - bins[0])
            total = np.zeros_like(grid)
            n_gen = len(fit['weights'])
            for g in range(n_gen):
                comp = fit['weights'][g] * scipy.stats.norm.pdf(grid, fit['means'][g], fit['sigma']) * scale
                total += comp
                ax.plot(10 ** grid, comp, color = cmap(g / max(1, n_gen - 1) * 0.9), lw = 1.2)
                if fit['weights'][g] > 0.01:
                    ax.annotate(str(g), (10 ** fit['means'][g], comp.max()), 
                                textcoords = 'offset points', xytext = (0, 2),
                                ha = 'center', fontsize = 'x-small')
            ax.plot(10 ** grid, total, color = 'k', lw = 1, ls = '--')
            ax.set_xscale('log')
            ax.set_ylabel("Events")
            
        plot_fit_grid(keys, op.by, draw, label, kwargs.get('title'))
        

util.expand_class_attributes(ProliferationOp)
util.expand_method_parameters(ProliferationOp, ProliferationOp.estimate)
