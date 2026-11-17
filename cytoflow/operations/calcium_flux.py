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
cytoflow.operations.calcium_flux
---------------------------------

Calcium flux kinetics analysis using fluorescent indicators.

`CalciumFluxOp` -- normalises fluorescence to a pre-stimulation baseline
(F/F₀), bins events by time, and classifies responding cells.

`CalciumFluxView` -- plots mean F/F₀ vs time with optional SEM shading.
"""

import numpy as np
import pandas as pd

from traits.api import (HasStrictTraits, Str, Float, Bool, Constant,
                        Instance, provides)

import cytoflow.utility as util

from .i_operation import IOperation


@provides(IOperation)
class CalciumFluxOp(HasStrictTraits):
    """
    Normalise a calcium-flux fluorescence channel to its pre-stimulation
    baseline (F/F₀) and add derived conditions and statistics to the
    experiment.

    For each event the operation computes::

        F/F₀ = signal / baseline_mean

    where ``baseline_mean`` is the mean signal across all events whose time
    value is ≤ ``baseline_end``.

    The operation adds:

    * A channel ``{name}_ratio`` containing F/F₀ for every event.
    * A condition ``{name}_time_bin`` (float) giving the left edge of each
      event's time bin (bin width = ``bin_width`` seconds).
    * A condition ``{name}_responding`` (category) equal to ``"Responding"``
      when ``f/f0 >= response_threshold``, else ``"Not responding"``.
    * Statistics keyed ``(name, 'mean_ratio')``, ``(name, 'sem_ratio')``
      and ``(name, 'pct_responding')`` — mean F/F₀, SEM, and % responding
      cells — indexed by time bin.

    Attributes
    ----------
    name : Str
        Output prefix used for new channels, conditions and statistics.
    time_channel : Str
        Channel recording time (e.g. ``"Time"``).
    signal_channel : Str
        Fluorescence channel (e.g. ``"Fluo4_A"``).
    ratio_denominator : Str
        For ratiometric dyes (e.g. Indo-1) set this to the denominator
        channel.  When non-empty the signal is computed as
        ``signal_channel / ratio_denominator`` before normalisation.
        Leave empty (default) for single-channel dyes.
    baseline_end : Float (default 60.0)
        Events with time ≤ this value (seconds) form the baseline population.
    response_threshold : Float (default 1.5)
        F/F₀ value above which a cell is classified as "Responding".
    bin_width : Float (default 2.0)
        Width of each time bin in seconds.

    Examples
    --------
    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import numpy as np, pandas as pd
        >>> ex = flow.ExperimentBuilder()  # illustrative; see test suite for real usage

    """

    id = Constant('cytoflow.operations.calcium_flux')
    friendly_id = Constant('Calcium Flux')

    name = Str
    time_channel = Str
    signal_channel = Str
    ratio_denominator = Str("")
    baseline_end = Float(60.0)
    response_threshold = Float(1.5)
    bin_width = Float(2.0)

    def apply(self, experiment):
        """
        Apply the calcium-flux analysis to *experiment*.

        Parameters
        ----------
        experiment : `Experiment`
            The experiment to process.

        Returns
        -------
        `Experiment`
            A new experiment with the ``{name}_ratio`` channel,
            ``{name}_time_bin`` condition, ``{name}_responding`` condition,
            and three statistics entries added.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.name:
            raise util.CytoflowOpError('name', "Must specify a name")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, "
                                       "numbers and underscores.")

        if not self.time_channel:
            raise util.CytoflowOpError('time_channel',
                                       "Must specify a time channel")

        if self.time_channel not in experiment.data.columns:
            raise util.CytoflowOpError('time_channel',
                                       "time_channel '{}' not found in the "
                                       "experiment".format(self.time_channel))

        if not self.signal_channel:
            raise util.CytoflowOpError('signal_channel',
                                       "Must specify a signal channel")

        if self.signal_channel not in experiment.data.columns:
            raise util.CytoflowOpError('signal_channel',
                                       "signal_channel '{}' not found in the "
                                       "experiment".format(self.signal_channel))

        if self.ratio_denominator and \
                self.ratio_denominator not in experiment.data.columns:
            raise util.CytoflowOpError('ratio_denominator',
                                       "ratio_denominator '{}' not found in "
                                       "the experiment"
                                       .format(self.ratio_denominator))

        ratio_channel = "{}_ratio".format(self.name)
        time_bin_condition = "{}_time_bin".format(self.name)
        responding_condition = "{}_responding".format(self.name)

        for col in (ratio_channel, time_bin_condition, responding_condition):
            if col in experiment.data.columns:
                raise util.CytoflowOpError('name',
                                           "'{}' is already in the experiment"
                                           .format(col))

        # ---- compute signal ----
        if self.ratio_denominator:
            signal = (experiment.data[self.signal_channel] /
                      experiment.data[self.ratio_denominator])
        else:
            signal = experiment.data[self.signal_channel]

        # ---- convert time to seconds ----
        time_step = experiment.metadata.get('time_step', 1.0)
        time = experiment.data[self.time_channel] * time_step

        # ---- baseline normalisation ----
        baseline_mask = time <= self.baseline_end
        if baseline_mask.sum() == 0:
            raise util.CytoflowOpError('baseline_end',
                                       "No events found before baseline_end "
                                       "= {}".format(self.baseline_end))

        baseline_mean = signal[baseline_mask].mean()
        if baseline_mean == 0 or np.isnan(baseline_mean):
            raise util.CytoflowOpError('signal_channel',
                                       "Baseline mean is zero or NaN; check "
                                       "signal_channel and baseline_end")

        f_over_f0 = signal / baseline_mean

        # ---- time bins ----
        time_bin = np.floor(time / self.bin_width) * self.bin_width

        # ---- responding classification ----
        responding = np.where(f_over_f0 >= self.response_threshold,
                              "Responding", "Not responding")

        # ---- build new experiment ----
        new_experiment = experiment.clone(deep=False)

        new_experiment.add_channel(ratio_channel, f_over_f0)
        new_experiment.add_condition(time_bin_condition, "float64", time_bin)
        new_experiment.add_condition(responding_condition, "category",
                                     pd.Categorical(responding,
                                                    categories=["Responding",
                                                                "Not responding"]))

        # ---- statistics ----
        agg_data = new_experiment.data[[time_bin_condition, ratio_channel,
                                        responding_condition]].copy()
        agg_data['_responding_int'] = (
            agg_data[responding_condition] == "Responding").astype(float)

        groupby = agg_data.groupby(time_bin_condition, observed=True)
        agg = groupby[ratio_channel].agg(['mean', 'sem']).rename(
            columns={'mean': 'mean_ratio', 'sem': 'sem_ratio'})
        agg['sem_ratio'] = agg['sem_ratio'].fillna(0.0)
        agg['pct_responding'] = groupby['_responding_int'].mean() * 100.0
        agg.index.name = time_bin_condition

        stat_mean = agg[['mean_ratio']]
        stat_sem = agg[['sem_ratio']]
        stat_pct = agg[['pct_responding']]

        new_experiment.statistics["{}_mean_ratio".format(self.name)] = stat_mean
        new_experiment.statistics["{}_sem_ratio".format(self.name)] = stat_sem
        new_experiment.statistics["{}_pct_responding".format(self.name)] = stat_pct

        new_experiment.history.append(
            self.clone_traits(transient=lambda _: True))

        return new_experiment

    def default_view(self, **kwargs):
        """
        Return a `CalciumFluxView` configured for this operation.

        Returns
        -------
        `CalciumFluxView`
        """
        v = CalciumFluxView(op=self)
        v.trait_set(**kwargs)
        return v


from cytoflow.views.i_view import IView as _IView


@provides(_IView)
class CalciumFluxView(HasStrictTraits):
    """
    Plot mean F/F₀ vs time for a calcium-flux experiment.

    Attributes
    ----------
    op : `CalciumFluxOp`
        The operation whose results to display.
    huefacet : Str
        A condition to colour separate lines by.
    rowfacet : Str
        A condition to arrange plots in rows.
    colfacet : Str
        A condition to arrange plots in columns.
    show_error_band : Bool (default True)
        Shade ±SEM around the mean line.
    show_stimulation_line : Bool (default True)
        Draw a vertical dashed line at ``op.baseline_end``.
    """

    id = Constant('cytoflow.views.calcium_flux')
    friendly_id = Constant('Calcium Flux')

    op = Instance(CalciumFluxOp)

    huefacet = Str("")
    rowfacet = Str("")
    colfacet = Str("")
    show_error_band = Bool(True)
    show_stimulation_line = Bool(True)

    def plot(self, experiment, **kwargs):
        """
        Plot mean F/F₀ vs time.

        Parameters
        ----------
        experiment : `Experiment`
            The experiment (after `CalciumFluxOp.apply`) to plot.
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")

        if self.op is None:
            raise util.CytoflowViewError('op', "No operation specified")

        name = self.op.name
        time_bin_condition = "{}_time_bin".format(name)
        ratio_channel = "{}_ratio".format(name)

        mean_key = "{}_mean_ratio".format(name)
        sem_key = "{}_sem_ratio".format(name)

        if mean_key not in experiment.statistics:
            raise util.CytoflowViewError('op',
                                         "Statistics not found; did you apply "
                                         "the operation first?")

        stat_mean = experiment.statistics[mean_key]
        stat_sem = experiment.statistics[sem_key]

        # determine facet structure
        rowfacet = self.rowfacet or None
        colfacet = self.colfacet or None
        huefacet = self.huefacet or None

        row_values = ([None] if rowfacet is None
                      else sorted(experiment.data[rowfacet].unique()))
        col_values = ([None] if colfacet is None
                      else sorted(experiment.data[colfacet].unique()))

        nrows = len(row_values)
        ncols = len(col_values)

        fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                                 figsize=(5 * ncols, 4 * nrows))

        hue_values = ([None] if huefacet is None
                      else sorted(experiment.data[huefacet].unique()))
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

        for ri, rval in enumerate(row_values):
            for ci, cval in enumerate(col_values):
                ax = axes[ri, ci]

                # subset data for this facet
                mask = pd.Series(True, index=experiment.data.index)
                if rval is not None:
                    mask &= experiment.data[rowfacet] == rval
                if cval is not None:
                    mask &= experiment.data[colfacet] == cval
                sub_data = experiment.data[mask]

                for hi, hval in enumerate(hue_values):
                    color = color_cycle[hi % len(color_cycle)]

                    if hval is not None:
                        hue_mask = sub_data[huefacet] == hval
                        grp_data = sub_data[hue_mask]
                    else:
                        grp_data = sub_data

                    # aggregate per time bin in one pass
                    grouped = grp_data.groupby(time_bin_condition, observed=True)
                    agg_v = grouped[ratio_channel].agg(['mean', 'sem'])
                    agg_v['sem'] = agg_v['sem'].fillna(0.0)
                    times = np.array(sorted(agg_v.index))
                    means = agg_v.loc[times, 'mean'].values
                    sems = agg_v.loc[times, 'sem'].values

                    label = str(hval) if hval is not None else None
                    ax.plot(times, means, color=color, label=label)

                    if self.show_error_band:
                        ax.fill_between(times,
                                        means - sems,
                                        means + sems,
                                        color=color, alpha=0.2)

                if self.show_stimulation_line:
                    ax.axvline(self.op.baseline_end, color='black',
                               linestyle='--', linewidth=1,
                               label='Stimulation')

                ax.set_xlabel("Time (s)")
                ax.set_ylabel("F/F₀")
                ax.axhline(1.0, color='grey', linestyle=':', linewidth=0.8)

                title_parts = []
                if rval is not None:
                    title_parts.append("{} = {}".format(rowfacet, rval))
                if cval is not None:
                    title_parts.append("{} = {}".format(colfacet, cval))
                if title_parts:
                    ax.set_title(" | ".join(title_parts))

                if huefacet:
                    ax.legend(title=huefacet)

        plt.tight_layout()
