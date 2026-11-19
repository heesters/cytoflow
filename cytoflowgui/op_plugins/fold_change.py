#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Fold Change
-----------

Compute the fold-change (or log2 fold-change) of a statistic between two
levels of an experimental condition, per grouping variable.

For each `by` group (e.g. each donor), divides the value in **Condition A**
by **Condition B**.  The result can be log2-transformed, which is standard
in immunology and shows symmetric up/down regulation.

.. object:: Name
    Output statistic name.

.. object:: Statistic
    Input statistic to compare.

.. object:: Column
    Which column of the statistic to use (e.g. ``% Parent``, ``Mean``).

.. object:: Group By (pairing)
    Variable(s) defining paired observations (e.g. Donor).

.. object:: Condition
    The condition being compared (e.g. Treatment).

.. object:: Condition A / B
    Numerator (A) and denominator (B).  The fold-change is A/B.

.. object:: Log2
    Report as log2(fold-change)? Standard immunology convention: 0 = no
    change, +1 = 2-fold up, -1 = 2-fold down.
"""

import numpy as np
import pandas as pd

from traits.api import provides, List, Property
from traitsui.api import (View, Item, TextEditor, EnumEditor, CheckListEditor,
                          BooleanEditor)
from pyface.api import ImageResource  # @UnresolvedImport
from envisage.api import Plugin

from ..workflow.operations import FoldChangeWorkflowOp

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, PluginHelpMixin, shared_op_traits_view


class FoldChangeHandler(OpHandler):

    stat_columns = Property(
        observe='context.previous_wi.statistics, model.statistic')

    # values of the currently selected condition in the statistic index
    condition_values = Property(
        observe='context.previous_wi.statistics, model.statistic, model.condition')

    stat_conditions = Property(
        observe='context.previous_wi.statistics, model.statistic')

    def _get_stat_columns(self):
        if not (self.context and self.context.previous_wi
                and self.context.previous_wi.statistics
                and self.model and self.model.statistic
                and self.model.statistic in self.context.previous_wi.statistics):
            return []
        stat = self.context.previous_wi.statistics[self.model.statistic]
        if isinstance(stat, pd.Series):
            return [stat.name] if stat.name else []
        return list(stat.columns)

    def _get_stat_conditions(self):
        if not (self.context and self.context.previous_wi
                and self.context.previous_wi.statistics
                and self.model and self.model.statistic
                and self.model.statistic in self.context.previous_wi.statistics):
            return []
        stat = self.context.previous_wi.statistics[self.model.statistic]
        if isinstance(stat, pd.Series):
            stat = stat.to_frame()
        return list(stat.index.names)

    def _get_condition_values(self):
        if not (self.context and self.context.previous_wi
                and self.context.previous_wi.statistics
                and self.model and self.model.statistic and self.model.condition
                and self.model.statistic in self.context.previous_wi.statistics):
            return []
        stat = self.context.previous_wi.statistics[self.model.statistic]
        if isinstance(stat, pd.Series):
            stat = stat.to_frame()
        if self.model.condition not in stat.index.names:
            return []
        vals = stat.index.get_level_values(self.model.condition).unique()
        return sorted(str(v) for v in vals)

    operation_traits_view = \
        View(Item('name',
                  editor=TextEditor(auto_set=False, placeholder="None"),
                  tooltip="Name for the fold-change statistic"),
             Item('statistic',
                  editor=EnumEditor(name='context_handler.previous_statistics_names'),
                  label="Statistic"),
             Item('column',
                  editor=EnumEditor(name='handler.stat_columns'),
                  label="Column",
                  tooltip="Which column of the statistic to use"),
             Item('by',
                  editor=CheckListEditor(cols=2,
                                         name='handler.stat_conditions'),
                  label='Group By\n(pairing)',
                  style='custom',
                  tooltip="Variable(s) defining paired observations, e.g. Donor"),
             Item('condition',
                  editor=EnumEditor(name='handler.stat_conditions'),
                  label="Condition",
                  tooltip="The condition to compare across"),
             Item('condition_a',
                  editor=EnumEditor(name='handler.condition_values'),
                  label="Cond. A\n(numerator)"),
             Item('condition_b',
                  editor=EnumEditor(name='handler.condition_values'),
                  label="Cond. B\n(denominator)"),
             Item('log2',
                  editor=BooleanEditor(),
                  label="log2(FC)",
                  tooltip="Report as log2 fold-change (standard: 0=no change, "
                          "+1=2-fold up, -1=2-fold down)"),
             shared_op_traits_view)


@provides(IOperationPlugin)
class FoldChangePlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.fold_change'
    operation_id = 'cytoflow.operations.fold_change'
    view_id = None

    name = "Fold Change"
    short_name = "Fold Change"
    menu_group = "Statistics"

    def get_operation(self):
        return FoldChangeWorkflowOp()

    def get_handler(self, model, context):
        return FoldChangeHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('channel_stat')

    plugin = List(contributes_to=OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
