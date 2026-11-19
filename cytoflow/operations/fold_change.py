#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.operations.fold_change
--------------------------------

Compute fold-change (or log2 fold-change) of a statistic between two
condition levels.

`FoldChangeOp` -- the operation
"""

import numpy as np
import pandas as pd

from traits.api import HasStrictTraits, Str, List, Bool, Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation


@provides(IOperation)
class FoldChangeOp(HasStrictTraits):
    """
    Compute the fold-change (or log2 fold-change) of a statistic between two
    levels of an experimental condition, per grouping variable.

    For each `by` group (e.g. each donor), divides the value in `condition_a`
    by the value in `condition_b`.

    Standard usage: compare gate frequencies or channel means between a
    treated and an untreated condition, expressed as fold-change per donor.

    Attributes
    ----------
    name : Str
        Output statistic name.
    statistic : Str
        Input statistic (must exist in the experiment's statistics).
    column : Str
        Column of the statistic to use (e.g. '% Parent', 'Mean').
    by : List[Str]
        Grouping variable(s) that define paired observations (e.g. Donor).
    condition : Str
        The condition being compared (e.g. Treatment).
    condition_a : Str
        Numerator (e.g. 'Stimulated').
    condition_b : Str
        Denominator (e.g. 'Unstimulated').
    log2 : Bool
        Report as log2(fold-change) — standard in immunology (default True).
    """

    id = Constant('cytoflow.operations.fold_change')
    friendly_id = Constant('Fold Change')

    name        = Str
    statistic   = Str
    column      = Str
    by          = List(Str)
    condition   = Str
    condition_a = Str
    condition_b = Str
    log2        = Bool(True)

    def apply(self, experiment):
        if not self.name:
            raise util.CytoflowOpError('name', "Must specify a name")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                "Name can only contain letters, numbers and underscores.")

        if not self.statistic:
            raise util.CytoflowOpError('statistic', "Must specify a statistic")

        if self.statistic not in experiment.statistics:
            raise util.CytoflowOpError('statistic',
                f"Statistic '{self.statistic}' not found. "
                f"Available: {sorted(experiment.statistics.keys())}")

        stat = experiment.statistics[self.statistic]
        if isinstance(stat, pd.Series):
            stat = stat.to_frame(name=stat.name or 'value')

        if not self.column:
            raise util.CytoflowOpError('column', "Must specify a column")

        if self.column not in stat.columns:
            raise util.CytoflowOpError('column',
                f"Column '{self.column}' not found in statistic '{self.statistic}'. "
                f"Available: {list(stat.columns)}")

        if not self.condition:
            raise util.CytoflowOpError('condition', "Must specify a condition")

        if self.condition not in stat.index.names:
            raise util.CytoflowOpError('condition',
                f"Condition '{self.condition}' is not an index level of "
                f"statistic '{self.statistic}'. "
                f"Index levels: {list(stat.index.names)}")

        if not self.condition_a or not self.condition_b:
            raise util.CytoflowOpError('condition_a',
                "Must specify both condition_a and condition_b.")

        if not self.by:
            raise util.CytoflowOpError('by',
                "Must specify grouping variable(s) in 'by'.")

        for b in self.by:
            if b not in stat.index.names:
                raise util.CytoflowOpError('by',
                    f"'{b}' is not an index level of statistic '{self.statistic}'. "
                    f"Index levels: {list(stat.index.names)}")

        if self.name in experiment.statistics:
            raise util.CytoflowOpError('name',
                f"'{self.name}' already exists in experiment statistics.")

        vals = stat[self.column].to_frame()
        vals.index = stat.index

        try:
            a_df = vals.xs(self.condition_a, level=self.condition)
            b_df = vals.xs(self.condition_b, level=self.condition)
        except KeyError as e:
            raise util.CytoflowOpError('condition_a',
                f"Could not select condition level: {e}") from e

        # Align on `by` index levels
        common = a_df.index.intersection(b_df.index)
        if len(common) == 0:
            raise util.CytoflowOpError('by',
                "No matching pairs found after filtering by condition. "
                "Check that each 'by' group appears in both conditions.")

        a_vals = a_df.loc[common, self.column].to_numpy(dtype=float)
        b_vals = b_df.loc[common, self.column].to_numpy(dtype=float)

        with np.errstate(divide='ignore', invalid='ignore'):
            fc = np.where(b_vals != 0, a_vals / b_vals, np.nan)
            if self.log2:
                fc_col = np.where(fc > 0, np.log2(fc), np.nan)
                col_name = 'log2FC'
            else:
                fc_col = fc
                col_name = 'FoldChange'

        result = pd.DataFrame(
            {col_name:    fc_col,
             'value_a':   a_vals,
             'value_b':   b_vals},
            index=common)

        new_experiment = experiment.clone(deep=False)
        new_experiment.history.append(self.clone_traits(transient=lambda _: True))
        new_experiment.statistics[self.name] = result

        return new_experiment
