#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.operations.gate_freq_test
----------------------------------

Paired statistical test on gate frequencies (% parent) across two conditions.

`GateFrequencyTestOp` -- the operation
"""

import warnings

import numpy as np
import pandas as pd
import scipy.stats

from traits.api import HasStrictTraits, Str, List, Enum, Float, Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation


@provides(IOperation)
class GateFrequencyTestOp(HasStrictTraits):
    """
    Paired statistical test comparing gate frequencies (% of parent) between
    two levels of a condition.

    For each `by` group (e.g. each donor), the operation computes the
    % parent of `gate` for `condition_a` and `condition_b`, then runs a
    paired test across donors.

    This is the standard immunology analysis: "Is the frequency of CD4+ T
    cells higher in condition A vs condition B, across N donors?"

    Attributes
    ----------
    name : Str
        Output statistic name.
    gate : Str
        Boolean gate condition whose frequency to test.
    parent : List[Str]
        Boolean parent gates (denominator for % parent). Empty = all events.
    by : List[Str]
        Grouping variables (e.g. ['Donor']) — the paired units.
    condition : Str
        The condition being compared (e.g. 'Treatment').
    condition_a : Str
        First level of `condition` (numerator / "treated").
    condition_b : Str
        Second level of `condition` (reference / "control").
    test : Enum('Wilcoxon', 'Paired t-test')
        Statistical test. Wilcoxon signed-rank (default) is non-parametric
        and robust with small donor numbers.
    alpha : Float
        Significance threshold (default 0.05).
    """

    id = Constant('cytoflow.operations.gate_freq_test')
    friendly_id = Constant('Paired Gate Test')

    name       = Str
    gate       = Str
    parent     = List(Str)
    by         = List(Str)
    condition  = Str
    condition_a = Str
    condition_b = Str
    test       = Enum('Wilcoxon', 'Paired t-test')
    alpha      = Float(0.05)

    def apply(self, experiment):
        if not self.name:
            raise util.CytoflowOpError('name', "Must specify a name")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                "Name can only contain letters, numbers and underscores.")

        if not self.gate:
            raise util.CytoflowOpError('gate', "Must specify a gate condition")

        if self.gate not in experiment.conditions:
            raise util.CytoflowOpError('gate',
                f"Gate '{self.gate}' not found in the experiment.")

        if experiment.conditions[self.gate].dtype.kind != 'b':
            raise util.CytoflowOpError('gate',
                f"Gate '{self.gate}' must be a boolean condition.")

        for p in self.parent:
            if p not in experiment.conditions:
                raise util.CytoflowOpError('parent',
                    f"Parent gate '{p}' not found in the experiment.")
            if experiment.conditions[p].dtype.kind != 'b':
                raise util.CytoflowOpError('parent',
                    f"Parent gate '{p}' must be a boolean condition.")

        if not self.by:
            raise util.CytoflowOpError('by',
                "Must specify grouping variable(s) in 'by' (e.g. Donor).")

        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                    f"Condition '{b}' not found in the experiment.")

        if not self.condition:
            raise util.CytoflowOpError('condition',
                "Must specify the condition to compare.")

        if self.condition not in experiment.conditions:
            raise util.CytoflowOpError('condition',
                f"Condition '{self.condition}' not found.")

        if not self.condition_a or not self.condition_b:
            raise util.CytoflowOpError('condition_a',
                "Must specify both condition_a and condition_b.")

        if self.name in experiment.statistics:
            raise util.CytoflowOpError('name',
                f"'{self.name}' already exists in experiment statistics.")

        data = experiment.data

        # Build parent mask
        parent_mask = pd.Series(True, index=data.index)
        for p in self.parent:
            parent_mask &= data[p].astype(bool)

        gate_mask = parent_mask & data[self.gate].astype(bool)

        groupby_cols = self.by + [self.condition]
        groupby = data.groupby(groupby_cols, observed=True)

        # Compute % parent per by×condition cell
        rows = []
        for key, grp in groupby:
            if not isinstance(key, tuple):
                key = (key,)
            by_key = key[:len(self.by)]
            cond_val = key[len(self.by)]

            pc  = float(parent_mask.loc[grp.index].sum())
            gc  = float(gate_mask.loc[grp.index].sum())
            pct = 100.0 * gc / pc if pc > 0 else np.nan

            row = dict(zip(self.by, by_key))
            row[self.condition] = cond_val
            row['pct_parent'] = pct
            rows.append(row)

        df = pd.DataFrame(rows)

        a_vals = df[df[self.condition] == self.condition_a].set_index(self.by)['pct_parent']
        b_vals = df[df[self.condition] == self.condition_b].set_index(self.by)['pct_parent']

        common = a_vals.index.intersection(b_vals.index)
        if len(common) == 0:
            raise util.CytoflowOpError('by',
                f"No paired observations found. Check that each {self.by} "
                f"group appears in both '{self.condition_a}' and '{self.condition_b}'.")

        if len(common) < len(a_vals) or len(common) < len(b_vals):
            warnings.warn(
                f"Some {self.by} groups are not present in both conditions; "
                f"using {len(common)} paired observations.",
                util.CytoflowOpWarning)

        a_paired = a_vals.loc[common].to_numpy(dtype=float)
        b_paired = b_vals.loc[common].to_numpy(dtype=float)

        valid = ~(np.isnan(a_paired) | np.isnan(b_paired))
        if valid.sum() < 2:
            raise util.CytoflowOpError('by',
                "Need at least 2 valid paired observations to run the test.")

        a_paired = a_paired[valid]
        b_paired = b_paired[valid]

        if self.test == 'Wilcoxon':
            try:
                res = scipy.stats.wilcoxon(a_paired, b_paired, alternative='two-sided')
                stat, pval = float(res.statistic), float(res.pvalue)
            except Exception as e:
                raise util.CytoflowOpError('test',
                    f"Wilcoxon test failed: {e}") from e
        else:  # Paired t-test
            res = scipy.stats.ttest_rel(a_paired, b_paired)
            stat, pval = float(res.statistic), float(res.pvalue)

        result = pd.DataFrame(
            {'statistic': [stat],
             'p_value':   [pval],
             'significant': [pval < self.alpha],
             'n_pairs':   [int(valid.sum())],
             'mean_a':    [float(np.mean(a_paired))],
             'mean_b':    [float(np.mean(b_paired))],
             'mean_diff': [float(np.mean(a_paired - b_paired))]},
            index=pd.Index(
                [f"{self.condition}: {self.condition_a} vs {self.condition_b}"],
                name='Comparison'))

        new_experiment = experiment.clone(deep=False)
        new_experiment.history.append(self.clone_traits(transient=lambda _: True))
        new_experiment.statistics[self.name] = result

        # Also store per-pair values for plotting
        by_cols = ({b: common.get_level_values(b) for b in self.by}
                   if len(self.by) > 1
                   else {self.by[0]: list(common)})
        pairs_df = pd.DataFrame({
            **by_cols,
            self.condition_a: a_vals.loc[common].to_numpy(dtype=float),
            self.condition_b: b_vals.loc[common].to_numpy(dtype=float),
        })
        new_experiment.statistics[self.name + '_Pairs'] = pairs_df

        return new_experiment
