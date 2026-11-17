#!/usr/bin/env python3
# coding: latin-1

# (c) Brian Teague 2018-2026
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
#
# The algorithm in this module is a Python port of the AutoSpill R package
# (https://github.com/AdrianListon/autospill), which is
#
#   Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#   Software written by Carlos P. Roca, as research funded by the European
#   Union.
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
#   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
#   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
cytoflow.operations.autospill
-----------------------------

Compute a spillover (compensation) matrix from single-color controls with
the *AutoSpill* algorithm, and compensate with it.  `autospill` has:

`AutoSpillOp` -- the `IOperation` (a kind of `CompensationOp`.)

`autospill` -- the algorithm itself, as a function of numpy arrays.

`AutoSpillConvergenceView` -- a diagnostic plot of the algorithm's
convergence.

AutoSpill [1]_ doesn't need the controls to be gated into "positive" and
"negative" populations.  Instead:

1. For each single-color control, it fits a *robust linear regression* of
   every other channel against the control's own channel.  The slopes are
   the first estimate of the spillover matrix.

2. It compensates the controls with that matrix and fits the regressions
   again.  In perfectly compensated controls all the slopes are 0; the
   slopes that are left are the compensation *errors*, which are used to
   correct the matrix.  This is repeated until the errors vanish -- first
   with the data on a linear scale, then on a log-like scale (where errors
   in the dim events are visible.)

Differences from the R package: the scatter gate (step 0 of the R package)
isn't done here -- gate the cells with earlier operations (they are applied
to the controls too) and `subset`; the log-like scale is ``arcsinh`` instead
of FlowJo's bi-exponential; and there is no separate autofluorescence
"dye".

References
----------
.. [1] Roca CP, Burton OT, Gergelits V, et al.  AutoSpill is a principled
   framework that simplifies the analysis of multichromatic flow cytometry
   data.  Nat Commun 12, 2890 (2021).
   https://doi.org/10.1038/s41467-021-23126-8
"""

import os, warnings

from traits.api import (HasStrictTraits, Str, File, Dict, Float, Int, Any,
                        Constant, Instance, provides)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cytoflow.views
import cytoflow.utility as util

from .i_operation import IOperation
from .compensation import CompensationOp
from .import_op import Tube, ImportOp, check_tube


### the algorithm

def robust_linear_fit(x, y, max_iter = 100, k = 1.345, acc = 1e-4):
    """
    Robust simple linear regression of ``y`` on ``x`` (Huber's M-estimator,
    fitted by iteratively re-weighted least squares, with the scale
    re-estimated from the median absolute residual at each step -- the same
    as the default ``MASS::rlm`` in R.)  Falls back to ordinary least
    squares if it doesn't converge.

    Returns
    -------
    (intercept, slope)
    """
    x = np.asarray(x, dtype = float)
    y = np.asarray(y, dtype = float)

    if len(x) < 3 or np.ptp(x) == 0:
        return (float(np.mean(y)) if len(y) else 0.0, 0.0)

    def wls(w):
        sw = w.sum()
        mx = (w * x).sum() / sw
        my = (w * y).sum() / sw
        sxx = (w * (x - mx) ** 2).sum()
        if sxx == 0:
            return my, 0.0
        slope = (w * (x - mx) * (y - my)).sum() / sxx
        return my - slope * mx, slope

    ols = wls(np.ones_like(x))
    inte, slope = ols
    resid = y - (inte + slope * x)

    for _ in range(max_iter):
        scale = np.median(np.abs(resid)) / 0.6745
        if scale == 0:
            return inte, slope

        u = np.abs(resid) / scale
        w = np.where(u <= k, 1.0, k / np.maximum(u, 1e-300))
        inte, slope = wls(w)

        new_resid = y - (inte + slope * x)
        change = np.sqrt(((resid - new_resid) ** 2).sum()
                         / max(1e-20, (resid ** 2).sum()))
        resid = new_resid
        if change < acc:
            return inte, slope

    return ols


def _trim_limits(v, trim):
    """The ``trim`` and ``1 - trim`` extreme values of ``v`` (exclusive.)"""
    n = len(v)
    t = int(round(n * trim))
    if t < 1 or n < 3:
        return -np.inf, np.inf
    s = np.sort(v)
    return s[t - 1], s[n - t]


def _trimmed_fit(proper, other, trim, max_iter):
    """
    Regress ``other`` on ``proper``, after dropping the events with extreme
    values in either.  Returns (intercept, slope, min(proper), max(proper)).
    """
    p_lo, p_hi = _trim_limits(proper, trim)
    o_lo, o_hi = _trim_limits(other, trim)
    keep = (proper > p_lo) & (proper < p_hi) & (other > o_lo) & (other < o_hi)
    x = proper[keep]
    y = other[keep]
    if len(x) < 3:
        return 0.0, 0.0, 0.0, 0.0
    inte, slope = robust_linear_fit(x, y, max_iter = max_iter)
    return inte, slope, float(x.min()), float(x.max())


def _slopes(data, transform, inverse, trim, max_iter):
    """
    The matrix of regression slopes, *in linear units*: ``slopes[i, j]`` is
    the slope of channel ``j`` vs channel ``i`` in control ``i``.  If
    ``transform`` is given, the regression is done on transformed data and
    the slope of the fitted line is converted back to linear units.
    """
    n = len(data)
    slopes = np.eye(n)
    for i in range(n):
        d = data[i] if transform is None else transform(data[i])
        proper = d[:, i]
        for j in range(n):
            if i == j:
                continue
            inte, coef, y1p, y2p = _trimmed_fit(proper, d[:, j], trim, max_iter)

            if transform is None:
                slopes[i, j] = coef
                continue

            # the fitted line's end points, in the transformed and the
            # linear scales
            x1p = inte + coef * y1p
            x2p = inte + coef * y2p
            if y1p == y2p or x1p == x2p:
                slopes[i, j] = 0.0
            else:
                y1, y2 = inverse(np.array([y1p, y2p]))
                x1, x2 = inverse(np.array([x1p, x2p]))
                slopes[i, j] = (x2 - x1) / (y2 - y1) if y2 != y1 else 0.0
    return slopes


def autospill(data,
              cofactor = 150.0,
              trim = 0.01,
              max_iter = 100,
              rlm_max_iter = 100,
              lambda_coarse = 1.0,
              lambda_fine = 0.1,
              history_n = 10,
              threshold_linear = 1e-2,
              threshold_transformed = 1e-4,
              threshold_change = 1e-6):
    """
    Compute a spillover matrix from single-color controls with the
    AutoSpill algorithm.

    Parameters
    ----------
    data : list of numpy.ndarray
        ``data[i]`` is the (events x channels) array of uncompensated data
        from the single-color control for channel ``i``.  All the arrays
        have the same channels, in the same order as the list.

    cofactor : float
        The log-like scale used in the second stage of the refinement is
        ``arcsinh(x / cofactor)``.

    trim : float
        The fraction of events with the most extreme values (at each end,
        in either channel) left out of each regression.

    max_iter : int
        The maximum number of refinement iterations.

    Returns
    -------
    (spillover, convergence, error)
        ``spillover`` is the (channels x channels) spillover matrix:
        ``spillover[i, j]`` is the fraction of channel ``i``'s signal
        detected in channel ``j`` (diagonal = 1).  ``convergence`` is a
        `pandas.DataFrame` with one row per iteration (``scale``,
        ``lambda``, ``delta`` -- the sd of the slope errors --,
        ``delta_max`` and ``delta_change``).  ``error`` is the matrix of
        residual slopes after the last iteration.

    Raises
    ------
    `CytoflowOpError`
        If the refinement doesn't converge.
    """

    n = len(data)
    data = [np.asarray(d, dtype = float) for d in data]
    for d in data:
        if d.ndim != 2 or d.shape[1] != n:
            raise util.CytoflowOpError(None,
                                       "Each control must have one column "
                                       "per channel")

    transform = lambda x: np.arcsinh(x / cofactor)
    inverse = lambda x: np.sinh(x) * cofactor
    eye = np.eye(n)

    # the initial spillover: regression slopes in the uncompensated data
    initial = _slopes(data, None, None, trim, rlm_max_iter)

    converged = False
    done = False
    it = 0
    last = False
    lam = lambda_coarse
    delta = -1.0
    threshold = threshold_linear
    history = np.full(history_n, -1.0)
    linear = True
    log = []

    spillover = eye.copy()
    update = initial - eye
    slope_error = np.zeros((n, n))

    while not done:
        # update the spillover matrix and compensate the controls
        spillover = spillover + update
        spillover = spillover / np.diag(spillover)[:, None]
        try:
            compensation = np.linalg.inv(spillover)
        except np.linalg.LinAlgError as e:
            raise util.CytoflowOpError(None,
                                       "The spillover matrix became singular. "
                                       "Are two of the controls the same, or "
                                       "is one of them unstained?") from e

        comp = [d @ compensation for d in data]

        # the compensation error: the slopes left in the compensated data
        if linear:
            slope_error = _slopes(comp, None, None, trim, rlm_max_iter) - eye
        else:
            slope_error = _slopes(comp, transform, inverse, trim, rlm_max_iter) - eye

        delta_prev = delta
        delta = float(np.std(slope_error, ddof = 1))
        delta_max = float(np.max(np.abs(slope_error)))

        history[it % history_n] = (delta - delta_prev) if delta_prev >= 0 else -1.0
        delta_change = float(np.mean(history))

        log.append(dict(iteration = it,
                        scale = "linear" if linear else "arcsinh",
                        lambda_ = lam,
                        delta = delta,
                        delta_max = delta_max,
                        delta_change = delta_change))

        if linear and delta_max < threshold:
            # switch to the log-like scale
            linear = False
            threshold = threshold_transformed
            lam = lambda_coarse
            delta = -1.0
            history[:] = -1.0
            delta_change = -1.0

        if delta_change > -threshold_change and lam == lambda_coarse:
            # no longer improving: take smaller steps
            lam = lambda_fine
            delta = -1.0
            history[:] = -1.0
            delta_change = -1.0

        converged = (not linear) and (delta_max < threshold
                                      or delta_change > -threshold_change)

        done = ((converged and last)
                or (delta_change > -threshold_change and linear)
                or (not converged and it == max_iter)
                or it > max_iter)

        last = converged
        it += 1

        update = lam * (slope_error @ spillover)

    convergence = pd.DataFrame(log).rename(columns = {'lambda_' : 'lambda'})

    if not converged:
        raise util.CytoflowOpError(None,
                                   "AutoSpill didn't converge after {} "
                                   "iterations (largest remaining slope "
                                   "error: {:.2g}). Check that each control "
                                   "is stained with one color only, and that "
                                   "it has both dim and bright events."
                                   .format(it, delta_max))

    return spillover, convergence, slope_error


### the operation

@provides(IOperation)
class AutoSpillOp(CompensationOp):
    """
    Compute the spillover matrix from single-color controls with the
    *AutoSpill* algorithm (Roca et al., Nat Commun 2021), then compensate.

    Unlike the traditional approach, the controls don't need to be gated
    into positive and negative populations: AutoSpill fits robust linear
    regressions to all of each control's events, then iteratively refines
    the matrix until the compensated controls show no remaining slope.

    This is a kind of `CompensationOp`: after `estimate`, the matrix is in
    `spillover` (where you can inspect or edit it) and `apply` compensates
    the data in the same way.

    .. note::
        Don't use this together with `CompensationOp` or
        `BleedthroughLinearOp` -- that would compensate the data twice.

    Attributes
    ----------
    controls : Dict(Str, File)
        The channel names to compensate, mapped to the FCS files of their
        single-color controls.  Every channel in the matrix needs a control.

    cofactor : Float (default = 150)
        The second stage of the refinement works on ``arcsinh(x / cofactor)``
        transformed data.

    trim : Float (default = 0.01)
        The fraction of events with the most extreme values that is left out
        of each regression.

    max_events : Int (default = 100000)
        Use at most this many events from each control (a random sample.)
        Set to 0 to use them all.

    max_iterations : Int (default = 100)
        Give up if the refinement hasn't converged after this many rounds.

    Notes
    -----
    The operations already applied to the experiment (eg, a scatter gate)
    are applied to the controls too, and ``subset`` in `estimate` then
    selects the events to use -- so gate the cells (or beads) first, then
    pass eg ``subset = "Cells == True"``.

    Examples
    --------

    >>> as_op = flow.AutoSpillOp()
    >>> as_op.controls = {"FITC-A" : "controls/fitc.fcs",
    ...                   "PE-A" : "controls/pe.fcs",
    ...                   "APC-A" : "controls/apc.fcs"}
    >>> as_op.estimate(ex)
    >>> as_op.default_view().plot(ex)
    >>> ex2 = as_op.apply(ex)
    """

    id = Constant('cytoflow.operations.autospill')
    friendly_id = Constant("AutoSpill Compensation")
    name = Constant("AutoSpill")

    controls = Dict(Str, File)
    cofactor = Float(150.0)
    trim = Float(0.01)
    max_events = Int(100000)
    max_iterations = Int(100)

    _convergence = Instance(pd.DataFrame, transient = True)
    _slope_error = Any(transient = True)

    def estimate(self, experiment, subset = None):
        """
        Compute the spillover matrix from the single-color controls in
        `controls`.

        Parameters
        ----------
        experiment : Experiment
            The experiment the controls belong to.  Its channel names and
            previous operations are used to read and pre-process the
            controls.

        subset : Str
            Which events of the controls to use (after the experiment's
            previous operations have been applied to them.)
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        channels = list(self.controls.keys())

        if len(channels) < 2:
            raise util.CytoflowOpError('controls',
                                       "Need at least two channels (and "
                                       "their controls) to compute spillover.")

        if not self.cofactor > 0:
            raise util.CytoflowOpError('cofactor', "cofactor must be > 0")

        if not 0 <= self.trim < 0.5:
            raise util.CytoflowOpError('trim', "trim must be between 0 and 0.5")

        for channel in channels:
            if channel not in experiment.channels:
                raise util.CytoflowOpError('controls',
                                           "Channel {} isn't in the experiment"
                                           .format(channel))
            if not self.controls[channel] or not os.path.isfile(self.controls[channel]):
                raise util.CytoflowOpError('controls',
                                           "Can't find file {0} for channel {1}."
                                           .format(self.controls[channel], channel))

        if len(set(os.path.abspath(f) for f in self.controls.values())) != len(channels):
            raise util.CytoflowOpError('controls',
                                       "Each channel needs its own "
                                       "single-color control file.")

        self.spillover = {}
        self._convergence = None
        self._slope_error = None

        data = []
        for channel in channels:
            check_tube(self.controls[channel], experiment)

            tube_exp = ImportOp(tubes = [Tube(file = self.controls[channel])],
                                channels = {experiment.metadata[c]["fcs_name"] : c
                                            for c in experiment.channels},
                                name_metadata = experiment.metadata['name_metadata']).apply()

            # apply previous operations (except the first, which is import)
            for op in experiment.history[1:]:
                if hasattr(op, 'by'):
                    for by in op.by:
                        if 'experiment' in experiment.metadata[by]:
                            raise util.CytoflowOpError('experiment',
                                                       "Prior to applying this operation, "
                                                       "you must not apply any operation with 'by' "
                                                       "set to an experimental condition.")
                tube_exp = op.apply(tube_exp)

            if subset:
                try:
                    tube_exp = tube_exp.query(subset)
                except Exception as exc:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{0}' isn't valid"
                                               .format(subset)) from exc

                if len(tube_exp.data) == 0:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{0}' returned no events"
                                               .format(subset))

            tube_data = tube_exp.data[channels].dropna()
            if self.max_events and len(tube_data) > self.max_events:
                tube_data = tube_data.sample(n = self.max_events, random_state = 0)

            if len(tube_data) < 100:
                raise util.CytoflowOpError('controls',
                                           "The control for {} has fewer than "
                                           "100 events".format(channel))

            data.append(tube_data.values)

        # is each control brightest in its own channel?  if not, the files
        # are probably mixed up.
        for i, channel in enumerate(channels):
            spread = np.percentile(data[i], 99, axis = 0) - np.percentile(data[i], 1, axis = 0)
            brightest = channels[int(np.argmax(spread))]
            if brightest != channel:
                warnings.warn("The control for {} varies more in channel {} "
                              "than in {}. Is it the right file?"
                              .format(channel, brightest, channel),
                              util.CytoflowOpWarning)

        try:
            spillover, convergence, error = autospill(data,
                                                      cofactor = self.cofactor,
                                                      trim = self.trim,
                                                      max_iter = self.max_iterations)
        except util.CytoflowOpError as e:
            raise util.CytoflowOpError('controls', e.args[-1]) from e

        self._convergence = convergence
        self._slope_error = pd.DataFrame(error, index = channels, columns = channels)

        # set atomically to support the GUI
        self.spillover = pd.DataFrame(spillover, index = channels, columns = channels)

    def default_view(self, **kwargs):
        """
        Returns a diagnostic view.  By default, a heat map of the spillover
        matrix; pass ``convergence = True`` for the convergence plot or
        ``qc = True`` for the before/after QC heatmaps.
        """
        if kwargs.pop('convergence', False):
            v = AutoSpillConvergenceView(op = self)
            v.trait_set(**kwargs)
            return v
        if kwargs.pop('qc', False):
            v = AutoSpillQCView(op = self)
            v.trait_set(**kwargs)
            return v
        return super().default_view(**kwargs)


@provides(cytoflow.views.IView)
class AutoSpillConvergenceView(HasStrictTraits):
    """
    Plots the convergence of the AutoSpill refinement: the largest
    remaining slope error (and its standard deviation) at each iteration.

    Attributes
    ----------
    op : Instance(AutoSpillOp)
        The operation whose convergence to plot.
    """

    id = Constant('cytoflow.views.autospill_convergence')
    friendly_id = Constant("AutoSpill Convergence")

    op = Instance(AutoSpillOp)

    def plot(self, experiment = None, **kwargs):
        """
        Plot the convergence.  ``experiment`` isn't used.
        """
        if self.op is None or self.op._convergence is None:
            raise util.CytoflowViewError('op',
                                         "Convergence isn't available. Did "
                                         "you forget to call estimate()?")

        conv = self.op._convergence
        plt.figure()
        ax = plt.gca()
        for scale, marker in (("linear", "^"), ("arcsinh", "o")):
            c = conv[conv["scale"] == scale]
            if len(c):
                ax.plot(c["iteration"], c["delta_max"], marker = marker,
                        color = "C1", label = "largest error ({})".format(scale))
                ax.plot(c["iteration"], c["delta"], marker = marker,
                        color = "C0", label = "sd of errors ({})".format(scale))
        ax.set_yscale("log")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Remaining slope error")
        ax.legend()
        if kwargs.get('title'):
            ax.set_title(kwargs['title'])


@provides(cytoflow.views.IView)
class AutoSpillQCView(HasStrictTraits):
    """
    Shows before-and-after QC heatmaps for the AutoSpill compensation.

    The left panel shows the computed spillover matrix (how much of channel
    *i*'s signal bleeds into detector *j*); the right panel shows the residual
    slope error after compensation (ideally all zeros).  Any cell above 1%
    spillover is annotated with its value.

    Attributes
    ----------
    op : Instance(AutoSpillOp)
        The operation whose QC to display.
    """

    id = Constant('cytoflow.views.autospill_qc')
    friendly_id = Constant("AutoSpill QC")

    op = Instance(AutoSpillOp)

    def plot(self, experiment=None, **kwargs):
        """
        Plot the spillover matrix and slope-error heatmaps side-by-side.
        ``experiment`` isn't used.
        """
        import seaborn as sns

        if self.op is None:
            raise util.CytoflowViewError('op', "No operation specified")

        if not self.op.spillover:
            raise util.CytoflowViewError('op',
                                         "Spillover matrix isn't available. "
                                         "Did you forget to call estimate()?")

        channels = list(self.op.spillover.keys()) \
            if isinstance(self.op.spillover, dict) else list(self.op.spillover.index)

        import pandas as pd
        if isinstance(self.op.spillover, dict):
            spill_df = pd.DataFrame(self.op.spillover)
        else:
            spill_df = pd.DataFrame(self.op.spillover)

        # off-diagonal only (diagonal is always 1 by construction)
        spill_pct = (spill_df * 100.0).round(2)
        np.fill_diagonal(spill_pct.values, np.nan)

        fig, axes = plt.subplots(1, 2, figsize=(max(8, len(channels) * 1.2),
                                                 max(4, len(channels) * 0.9)))

        # Spillover matrix
        ax = axes[0]
        mask = np.isnan(spill_pct.values)
        sns.heatmap(spill_pct.fillna(0),
                    ax=ax, annot=True, fmt='.1f',
                    cmap='YlOrRd', mask=mask,
                    linewidths=0.5, linecolor='white',
                    vmin=0, vmax=max(10, spill_pct.fillna(0).values.max()))
        ax.set_title("Spillover matrix (%)")
        ax.set_xlabel("Detector (leaks into)")
        ax.set_ylabel("Channel (control)")

        # Slope-error heatmap
        ax2 = axes[1]
        if self.op._slope_error is not None:
            err_pct = (self.op._slope_error * 100.0).round(2)
            np.fill_diagonal(err_pct.values, np.nan)
            sns.heatmap(err_pct.fillna(0),
                        ax=ax2, annot=True, fmt='.2f',
                        cmap='coolwarm', center=0,
                        linewidths=0.5, linecolor='white')
            ax2.set_title("Residual slope error after compensation (%)")
        else:
            ax2.set_title("Residual slope error (not available)")
        ax2.set_xlabel("Detector")
        ax2.set_ylabel("Channel")

        plt.tight_layout()

        if kwargs.get('title'):
            fig.suptitle(kwargs['title'], y=1.02)


util.expand_class_attributes(AutoSpillOp)
util.expand_method_parameters(AutoSpillOp, AutoSpillOp.estimate)
