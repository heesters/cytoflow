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
cytoflow.operations.cell_cycle
------------------------------

Cell-cycle analysis from a DNA-content histogram (propidium iodide, DAPI,
Hoechst, ...).  `cell_cycle` has:

`CellCycleOp` -- fits the model and computes the percentage of cells in 
G0/G1, S and G2/M.

`fit_cell_cycle` -- the model fit, as a function of a numpy array.

`CellCycleView` -- plots the histogram with the fitted phases.
"""

import warnings

from traits.api import (HasStrictTraits, Str, List, Dict, Int, Any, Constant, 
                        Instance, provides)

import numpy as np
import pandas as pd
import scipy.optimize
import scipy.ndimage
import scipy.stats
import matplotlib.pyplot as plt

import cytoflow.views
import cytoflow.utility as util

from .i_operation import IOperation
from .proliferation import plot_fit_grid

PHASES = ["G1", "S", "G2M"]

# the S phase is modelled as this many evenly spaced, broadened "slices"
N_S = 12


def _components(x, m, cv, r, a1, a2, s):
    """
    The model's three densities at ``x``: G0/G1, S and G2/M.  ``s`` holds 
    the (non-negative) heights of the S phase distribution at three DNA 
    contents: just after G1, mid-S and just before G2 (a quadratic through
    them gives the weight of each S phase slice.)
    """
    g1 = a1 * scipy.stats.norm.pdf(x, m, cv * m)
    g2 = a2 * scipy.stats.norm.pdf(x, r * m, cv * r * m)
    
    t = (np.arange(N_S) + 0.5) / N_S              # 0 .. 1 through S phase
    pos = m * (1 + t * (r - 1))
    # the quadratic through (0, s0), (0.5, s1), (1, s2)
    w = (s[0] * 2 * (t - 0.5) * (t - 1) - s[1] * 4 * t * (t - 1) + s[2] * 2 * t * (t - 0.5))
    w = np.maximum(w, 0) / N_S
    sp = (w[None, :] * scipy.stats.norm.pdf(x[:, None], pos[None, :], cv * pos[None, :])).sum(axis = 1)
    return g1, sp, g2, float(w.sum())


def fit_cell_cycle(x, bins = 256):
    """
    Fit a DNA-content histogram with a Dean-Jett-Fox style model:
    
    - G0/G1: a Gaussian at the 2N DNA content, ``m``;
    - G2/M: a Gaussian at ``r * m`` (``r`` is fitted, between 1.8 and 2.2);
    - S: cells spread between them according to a quadratic polynomial, 
      *broadened* by the same measurement error as the peaks.
      
    The measurement error has a constant coefficient of variation (CV), so
    the G2/M peak is ``r`` times as wide as the G0/G1 peak.  The model is
    fitted to the histogram by non-linear least squares, starting from the
    histogram's main peak; the fit is deterministic.
    
    Events far outside the model (more than 4 standard deviations below 
    G0/G1, or above G2/M) are debris / apoptotic cells ("sub-G1") or 
    aggregates and polyploid cells; they are counted but not fitted.
    
    Parameters
    ----------
    x : numpy.ndarray
        The DNA signal of each event (linear; use the area, and gate out
        doublets first.)
        
    Returns
    -------
    dict
        ``G1``, ``S``, ``G2M``: the fraction of the *modelled* cells in each
        phase; ``SubG1``, ``AboveG2M``: the fraction of *all* events below
        and above the model; ``G1_mean``, ``G1_CV`` (a fraction), 
        ``G2_G1_ratio``, ``RMSE`` (the root-mean-square error of the fit,
        relative to the height of the G1 peak) and ``params``.
    """
    x = np.asarray(x, dtype = float)
    x = x[np.isfinite(x)]
    if len(x) < 200:
        raise util.CytoflowOpError(None,
                                   "Need at least 200 events to fit a cell cycle; got {}"
                                   .format(len(x)))
        
    # find the G0/G1 peak: the main mode of the (smoothed) histogram
    lo, hi = np.percentile(x, [0.5, 99.5])
    if not hi > lo:
        raise util.CytoflowOpError(None, "All the events have the same DNA content")
    hist, edges = np.histogram(x, bins = bins, range = (lo, hi))
    centers = (edges[:-1] + edges[1:]) / 2
    smooth = scipy.ndimage.gaussian_filter1d(hist.astype(float), 2)
    peak = int(np.argmax(smooth))
    m0 = centers[peak]
    
    # is the main peak G2/M rather than G1?  if there's a substantial peak at 
    # half its DNA content, that one is G1.
    half = np.searchsorted(centers, m0 / 2)
    if 2 < half < len(centers) - 2:
        window = slice(max(0, half - bins // 16), min(len(centers), half + bins // 16))
        if smooth[window].max() > 0.25 * smooth[peak]:
            m0 = centers[window][int(np.argmax(smooth[window]))]
            
    # the peak's width, from its half-maximum on the left (the right side
    # has S phase in it)
    i = int(np.searchsorted(centers, m0))
    i = min(max(i, 1), len(centers) - 1)
    height = smooth[i]
    j = i
    while j > 0 and smooth[j] > height / 2:
        j -= 1
    cv0 = float(np.clip((m0 - centers[j]) / (1.1774 * m0), 0.01, 0.15)) if m0 > 0 else 0.05
    
    # fit the histogram, as a density, inside the model's range
    x_lo = m0 * (1 - 5 * cv0)
    x_hi = 2.2 * m0 * (1 + 5 * cv0)
    inside = x[(x >= x_lo) & (x <= x_hi)]
    if len(inside) < 100:
        raise util.CytoflowOpError(None, "Couldn't find a G0/G1 peak")
    # work in units of the G1 peak's position, so that the density (and the
    # optimizer's tolerances) are of order 1 whatever the instrument's range
    h, e = np.histogram(inside / m0, bins = bins, range = (x_lo / m0, x_hi / m0), 
                        density = True)
    c = (e[:-1] + e[1:]) / 2
    
    def model(p):
        g1, sp, g2, _ = _components(c, p[0], p[1], p[2], p[3], p[4], p[5:8])
        return g1 + sp + g2
    
    frac_g2 = float(np.mean((inside > 1.8 * m0 * (1 - 2 * cv0)) & (inside < 2.2 * m0 * (1 + 2 * cv0))))
    lower = np.array([0.9, 0.005, 1.8, 0, 0, 0, 0, 0])
    upper = np.array([1.1, 0.2, 2.2, 1.5, 1.5, 5, 5, 5])
    
    # a few starting points for the amount of S phase; keep the best fit
    best = None
    for s0 in (0.05, 0.3, 0.6):
        p0 = np.array([1.0, cv0, 2.0, max(0.05, 0.9 - s0 - frac_g2), max(0.02, frac_g2), s0, s0, s0])
        res = scipy.optimize.least_squares(lambda p: model(p) - h, p0, 
                                           bounds = (lower, upper))
        if best is None or res.cost < best.cost:
            best = res
    res = best
    
    # back to the data's units
    p = res.x.copy()
    p[0] *= m0
    m, cv, r, a1, a2 = p[:5]
    _, _, _, s_total = _components(c[:2], m, cv, r, a1, a2, p[5:8])
    total = a1 + a2 + s_total
    if total <= 0:
        raise util.CytoflowOpError(None, "The cell cycle fit failed")
    
    if a2 / total < 0.005 and s_total / total < 0.005:
        warnings.warn("Only one peak was found: there is no sign of S or G2/M "
                      "cells. Is this a DNA stain on a linear scale?",
                      util.CytoflowOpWarning)
    
    below = m * (1 - 4 * cv)
    above = r * m * (1 + 4 * cv)
    
    return dict(G1 = float(a1 / total), 
                S = float(s_total / total), 
                G2M = float(a2 / total),
                SubG1 = float(np.mean(x < below)),
                AboveG2M = float(np.mean(x > above)),
                G1_mean = float(m),
                G1_CV = float(cv),
                G2_G1_ratio = float(r),
                RMSE = float(np.sqrt(np.mean(res.fun ** 2)) / max(h.max(), 1e-300)),
                params = p,
                limits = (float(below), float(above)))


def _assign_phase(x, fit):
    """The most likely phase of each event: an array of strings"""
    p = fit['params']
    g1, sp, g2, _ = _components(x, p[0], p[1], p[2], p[3], p[4], p[5:8])
    phase = np.array(PHASES, dtype = object)[np.argmax(np.column_stack((g1, sp, g2)), axis = 1)]
    below, above = fit['limits']
    phase[x < below] = "SubG1"
    phase[x > above] = "AboveG2M"
    return phase


@provides(IOperation)
class CellCycleOp(HasStrictTraits):
    """
    Cell-cycle analysis: fit the G0/G1, S and G2/M phases of a DNA-content
    histogram.
    
    `estimate` fits the model (see `fit_cell_cycle`) for each group of `by`;
    `apply` adds:
    
    - a statistic named `name` with, for each group: ``G1``, ``S`` and 
      ``G2M`` (the percentage of the cycling cells in each phase; they add 
      up to 100), ``SubG1`` and ``AboveG2M`` (the percentage of *all* the 
      events below G0/G1 -- debris and apoptotic cells -- and above G2/M --
      aggregates and polyploid cells), ``G1_mean``, ``G1_CV`` (in percent),
      ``G2_G1_ratio`` and ``RMSE`` (how well the model fits: the 
      root-mean-square error as a fraction of the height of the G1 peak; 
      more than about 0.05 deserves a look at the plot); and
    - a condition named `name` with each event's most likely phase: 
      ``G1``, ``S``, ``G2M``, ``SubG1`` or ``AboveG2M``.  This is 
      approximate where the phases overlap -- use the statistic for the
      percentages.
      
    The channel must be the DNA stain's **area** on a **linear** scale, and 
    doublets must be gated out first (eg, with a polygon gate on area vs 
    height or width): two G1 cells stuck together look like one G2 cell.
    
    Attributes
    ----------
    name : Str
        The operation name; names the new statistic and condition.
        
    channel : Str
        The DNA stain's channel.
        
    by : List(Str)
        Fit separately for each combination of these conditions' values 
        (eg, each sample.)
    """
    
    id = Constant('cytoflow.operations.cell_cycle')
    friendly_id = Constant("Cell Cycle")
    
    name = Str
    channel = Str
    by = List(Str)
    
    _fits = Dict(Any, Any, transient = True)       # group -> fit
    
    def estimate(self, experiment, subset = None):
        """
        Fit the cell cycle model.
        
        Parameters
        ----------
        experiment : Experiment
        
        subset : Str
            Only use these events (eg, the single cells.)
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")
        
        if not self.channel:
            raise util.CytoflowOpError('channel', "Must set a channel")
        
        if self.channel not in experiment.channels:
            raise util.CytoflowOpError('channel',
                                       "Channel {} isn't in the experiment".format(self.channel))
            
        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "Condition {} isn't in the experiment".format(b))
                
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
                
        if self.by:
            groups = experiment.data.groupby(self.by, observed = True)
        else:
            groups = [((), experiment.data)]
            
        fits = {}
        for group, data in groups:
            key = group if isinstance(group, tuple) else (group,)
            try:
                fits[key] = fit_cell_cycle(data[self.channel].values)
            except util.CytoflowOpError as e:
                raise util.CytoflowOpError('by' if self.by else 'channel', 
                                           "{}{}".format("Group {}: ".format(key) if key else "", 
                                                         e.args[-1])) from e
                
        # set atomically, for the GUI
        self._fits = fits
        
    def apply(self, experiment):
        """
        Compute the cell cycle statistics and assign each event a phase.
        
        Returns
        -------
        Experiment
            A new experiment with the statistic `name` and the condition
            `name`.
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
            
        phase = np.full(len(experiment), "None", dtype = object)
        
        if self.by:
            groups = experiment.data.groupby(self.by, observed = True).indices.items()
        else:
            groups = [((), np.arange(len(experiment)))]
            
        rows = []
        index = []
        for group, idx in groups:
            key = group if isinstance(group, tuple) else (group,)
            if key not in self._fits:
                continue
            fit = self._fits[key]
            x = experiment.data[self.channel].values[idx].astype(float)
            phase[idx] = _assign_phase(np.nan_to_num(x, nan = -np.inf), fit)
            rows.append(dict(G1 = 100 * fit['G1'], S = 100 * fit['S'], G2M = 100 * fit['G2M'],
                             SubG1 = 100 * fit['SubG1'], AboveG2M = 100 * fit['AboveG2M'],
                             G1_mean = fit['G1_mean'], G1_CV = 100 * fit['G1_CV'],
                             G2_G1_ratio = fit['G2_G1_ratio'], RMSE = fit['RMSE']))
            index.append(key)
            
        if self.by:
            idx = pd.MultiIndex.from_tuples(index, names = self.by) if len(self.by) > 1 \
                  else pd.Index([k[0] for k in index], name = self.by[0])
        else:
            idx = pd.Index(["All"], name = "Group")
        stat = pd.DataFrame(rows, index = idx)
        
        new_experiment = experiment.clone(deep = False)
        new_experiment.add_condition(self.name, "category", 
                                     pd.Series(phase).astype("category"))
        new_experiment.statistics[self.name] = stat
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment
    
    def default_view(self, **kwargs):
        """
        Returns a `CellCycleView`: the DNA histogram of each group with the
        fitted phases.
        """
        v = CellCycleView(op = self)
        v.trait_set(**kwargs)
        return v
    
    
@provides(cytoflow.views.IView)
class CellCycleView(HasStrictTraits):
    """
    Plots the DNA histogram of each group, with the fitted G0/G1, S and
    G2/M phases.
    
    Attributes
    ----------
    op : Instance(CellCycleOp)
        The operation whose fit to plot.
        
    subset : Str
        Only plot these events.
    """
    
    id = Constant('cytoflow.view.cellcycleview')
    friendly_id = Constant("Cell Cycle")
    
    op = Instance(CellCycleOp)
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
        colors = {"G1" : "C0", "S" : "C2", "G2M" : "C3"}
        
        def draw(ax, key):
            fit = op._fits[key]
            p = fit['params']
            x = data[key][op.channel].values.astype(float)
            x = x[np.isfinite(x)]
            lo = max(0.0, p[0] * (1 - 8 * p[1]))
            hi = p[2] * p[0] * (1 + 8 * p[1])
            bins = np.linspace(lo, hi, 200)
            ax.hist(x, bins = bins, color = '0.8')
            grid = np.linspace(lo, hi, 800)
            below, above = fit['limits']
            n_model = float(np.sum((x >= p[0] * (1 - 5 * p[1])) & (x <= 2.2 * p[0] * (1 + 5 * p[1]))))
            # the fit is a density over the fitted window
            g1, sp, g2, s_total = _components(grid, p[0], p[1], p[2], p[3], p[4], p[5:8])
            scale = n_model * (bins[1] - bins[0])
            for name, comp in (("G1", g1), ("S", sp), ("G2M", g2)):
                ax.fill_between(grid, comp * scale, color = colors[name], alpha = 0.35, lw = 0)
                ax.plot(grid, comp * scale, color = colors[name], lw = 1.2,
                        label = "{} {:.1f}%".format({"G1" : "G0/G1", "S" : "S", "G2M" : "G2/M"}[name], 
                                                    100 * fit[name]))
            ax.plot(grid, (g1 + sp + g2) * scale, color = 'k', lw = 1, ls = '--')
            ax.set_xlim(lo, hi)
            ax.set_ylabel("Events")
            ax.legend(fontsize = 'x-small', frameon = False)
            
        plot_fit_grid(keys, op.by, draw, label, kwargs.get('title'))
        

util.expand_class_attributes(CellCycleOp)
util.expand_method_parameters(CellCycleOp, CellCycleOp.estimate)
