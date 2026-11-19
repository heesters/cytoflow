#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.operations.gate_stat
-----------------------------

Compute gate frequencies: % of parent, % of total, and event count,
for a boolean gate condition, optionally grouped by experimental conditions.

`GateStatisticOp` -- computes gate frequencies
"""

import pandas as pd
import numpy as np

from traits.api import HasStrictTraits, Str, List, Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation


@provides(IOperation)
class GateStatisticOp(HasStrictTraits):
    """
    Compute gate frequencies (% parent, % total, count) for a boolean gate
    condition, grouped by experimental conditions.

    For each combination of `by` groups, computes:

    * **Count** — number of events where the gate condition is True (and all
      parent conditions are True).
    * **% Parent** — Count divided by the number of events where all *parent*
      conditions are True, expressed as a percentage.  If `parent` is empty,
      all events in the group are the denominator (equivalent to % Total).
    * **% Total** — Count divided by the total events in the group.

    Attributes
    ----------
    name : Str
        The name of the new statistic added to the experiment.
    gate : Str
        The name of the boolean gate condition to compute frequencies for.
    parent : List[Str]
        Boolean conditions whose conjunction defines the parent population
        (denominator for % Parent).  Typically the gates applied before this
        one in the analysis hierarchy.
    by : List[Str]
        Conditions to group by (e.g. ``['Donor', 'Treatment']``).  One row in
        the resulting statistic per unique combination.

    Examples
    --------
    >>> gs = flow.GateStatisticOp(name='LymphFreq',
    ...                           gate='Lymphocytes',
    ...                           parent=['Singlets', 'Live'],
    ...                           by=['Donor', 'Treatment'])
    >>> ex2 = gs.apply(ex)
    >>> print(ex2.statistics['LymphFreq'])
    """

    id = Constant('cytoflow.operations.gate_stat')
    friendly_id = Constant('Gate Frequencies')

    name   = Str
    gate   = Str
    parent = List(Str)
    by     = List(Str)

    def apply(self, experiment):
        """
        Apply the operation to an `Experiment`.

        Parameters
        ----------
        experiment : Experiment
            The experiment to apply this operation to.

        Returns
        -------
        Experiment
            A new experiment with a new statistics entry keyed by `name`.
            The statistic is a DataFrame with columns ``Count``,
            ``% Parent``, and ``% Total``, indexed by the `by` conditions.
        """

        if not self.name:
            raise util.CytoflowOpError('name', "Must specify a name")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError(
                'name',
                "Name can only contain letters, numbers and underscores.")

        if not self.gate:
            raise util.CytoflowOpError('gate', "Must specify a gate condition")

        if self.gate not in experiment.conditions:
            raise util.CytoflowOpError(
                'gate',
                f"Gate '{self.gate}' not found in the experiment. "
                f"Available conditions: {sorted(experiment.conditions.keys())}")

        if experiment.conditions[self.gate].dtype.kind != 'b':
            raise util.CytoflowOpError(
                'gate',
                f"Gate '{self.gate}' must be a boolean condition (True/False), "
                f"but it has dtype '{experiment.conditions[self.gate].dtype}'.")

        for p in self.parent:
            if p not in experiment.conditions:
                raise util.CytoflowOpError(
                    'parent',
                    f"Parent gate '{p}' not found in the experiment.")
            if experiment.conditions[p].dtype.kind != 'b':
                raise util.CytoflowOpError(
                    'parent',
                    f"Parent gate '{p}' must be a boolean condition, "
                    f"but it has dtype '{experiment.conditions[p].dtype}'.")

        if not self.by:
            raise util.CytoflowOpError(
                'by',
                "Must specify at least one grouping condition in 'by'. "
                "Use your sample labels (Donor, Treatment, Timepoint, etc.).")

        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError(
                    'by',
                    f"Grouping condition '{b}' not found. "
                    f"Must be one of {sorted(experiment.conditions.keys())}")

        if self.name in experiment.statistics:
            raise util.CytoflowOpError(
                'name',
                f"'{self.name}' is already in the experiment's statistics.")

        data = experiment.data

        # Build parent mask: events where ALL parent conditions are True
        parent_mask = pd.Series(True, index=data.index)
        for p in self.parent:
            parent_mask &= data[p].astype(bool)

        # Gate mask: gate is True AND within parent
        gate_mask = parent_mask & data[self.gate].astype(bool)

        groupby = data.groupby(self.by, observed=True)
        keys = [x if isinstance(x, tuple) else (x,)
                for x in groupby.groups.keys()]
        idx = pd.MultiIndex.from_tuples(keys, names=self.by)

        stat = pd.DataFrame(
            data=np.full((len(idx), 3), np.nan),
            index=idx,
            columns=['Count', '% Parent', '% Total'],
            dtype='float'
        ).sort_index()

        for group_key, group_data in groupby:
            group_idx  = group_data.index
            pc         = float(parent_mask.loc[group_idx].sum())
            gc         = float(gate_mask.loc[group_idx].sum())
            total      = float(len(group_data))
            key        = group_key if isinstance(group_key, tuple) else (group_key,)

            stat.loc[key, 'Count']    = gc
            stat.loc[key, '% Parent'] = 100.0 * gc / pc    if pc    > 0 else np.nan
            stat.loc[key, '% Total']  = 100.0 * gc / total if total > 0 else np.nan

        new_experiment = experiment.clone(deep=False)
        new_experiment.history.append(self.clone_traits(transient=lambda _: True))
        new_experiment.statistics[self.name] = stat

        return new_experiment
