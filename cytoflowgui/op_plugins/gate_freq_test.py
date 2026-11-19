#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Paired Gate Frequency Test
--------------------------

Paired statistical test (Wilcoxon or paired t-test) comparing gate frequencies
(% of parent) between two levels of an experimental condition.

This is the standard immunology analysis: "Is the frequency of this population
higher in condition A vs B, across N donors?"

The result is added as two statistics:

* **<Name>** — overall test: statistic, p_value, significant, n_pairs, mean_a,
  mean_b, mean_diff.
* **<Name>_Pairs** — per-donor values in both conditions (useful for plotting
  paired lines).

.. object:: Name
    Output statistic name.

.. object:: Gate
    The boolean gate condition whose frequency to compare.

.. object:: Parent gates
    Gates that define the parent population (denominator for % parent).
    Leave empty to use all events.

.. object:: Group By (pairing)
    The variable(s) that define paired observations (e.g. Donor).
    Each value must appear in both conditions.

.. object:: Condition
    The condition being compared (e.g. Treatment).

.. object:: Condition A / B
    The two levels to compare.  A is the "treated" group, B the reference.

.. object:: Test
    * **Wilcoxon** — non-parametric, recommended for small donor numbers.
    * **Paired t-test** — parametric; assumes normally distributed differences.

.. object:: alpha
    Significance threshold (default 0.05).
"""

import numpy as np
import pandas as pd

from traits.api import provides, List, Property
from traitsui.api import (View, Item, TextEditor, EnumEditor, CheckListEditor)
from pyface.api import ImageResource  # @UnresolvedImport
from envisage.api import Plugin

from ..workflow.operations import GateFrequencyTestWorkflowOp

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, PluginHelpMixin, shared_op_traits_view


class GateFrequencyTestHandler(OpHandler):

    # bool conditions only — for gate / parent dropdowns
    bool_conditions = Property(
        observe='context.previous_wi.conditions')

    # non-bool conditions — for condition / by dropdowns
    nonbool_conditions = Property(
        observe='context.previous_wi.conditions')

    # values of the currently selected condition
    condition_values = Property(
        observe='context.previous_wi.result, model.condition')

    def _get_bool_conditions(self):
        if (self.context and self.context.previous_wi
                and self.context.previous_wi.conditions):
            return sorted(k for k, v in self.context.previous_wi.conditions.items()
                          if v == 'bool')
        return []

    def _get_nonbool_conditions(self):
        if (self.context and self.context.previous_wi
                and self.context.previous_wi.conditions):
            return sorted(k for k, v in self.context.previous_wi.conditions.items()
                          if v != 'bool')
        return []

    def _get_condition_values(self):
        if not (self.context and self.context.previous_wi
                and self.context.previous_wi.result
                and self.model and self.model.condition
                and self.model.condition in self.context.previous_wi.result.conditions):
            return []
        vals = self.context.previous_wi.result.data[self.model.condition].unique()
        return sorted(str(v) for v in vals)

    operation_traits_view = \
        View(Item('name',
                  editor=TextEditor(auto_set=False, placeholder="None"),
                  tooltip="Name for this test result"),
             Item('gate',
                  editor=EnumEditor(name='handler.bool_conditions'),
                  label="Gate",
                  tooltip="The boolean gate condition to test"),
             Item('parent',
                  editor=CheckListEditor(cols=2, name='handler.bool_conditions'),
                  label='Parent\ngates',
                  style='custom',
                  tooltip="Gates defining the parent population (denominator)"),
             Item('by',
                  editor=CheckListEditor(cols=2,
                                         name='context_handler.previous_conditions_names'),
                  label='Group By\n(pairing)',
                  style='custom',
                  tooltip="Variable(s) defining paired observations, e.g. Donor"),
             Item('condition',
                  editor=EnumEditor(name='handler.nonbool_conditions'),
                  label="Condition",
                  tooltip="The condition being compared (e.g. Treatment)"),
             Item('condition_a',
                  editor=EnumEditor(name='handler.condition_values'),
                  label="Cond. A\n(treated)",
                  tooltip="First level (numerator / treated)"),
             Item('condition_b',
                  editor=EnumEditor(name='handler.condition_values'),
                  label="Cond. B\n(control)",
                  tooltip="Second level (reference / control)"),
             Item('test',
                  editor=EnumEditor(values={'Wilcoxon': '1:Wilcoxon (non-parametric)',
                                            'Paired t-test': '2:Paired t-test'}),
                  label="Test"),
             Item('alpha',
                  editor=TextEditor(auto_set=False, evaluate=float,
                                    format_func=lambda x: "" if x is None else str(x)),
                  label="alpha"),
             shared_op_traits_view)


@provides(IOperationPlugin)
class GateFrequencyTestPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.gate_freq_test'
    operation_id = 'cytoflow.operations.gate_freq_test'
    view_id = None

    name = "Paired Gate Test"
    short_name = "Gate Test"
    menu_group = "Statistics"

    def get_operation(self):
        return GateFrequencyTestWorkflowOp()

    def get_handler(self, model, context):
        return GateFrequencyTestHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('stat_test')

    plugin = List(contributes_to=OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
