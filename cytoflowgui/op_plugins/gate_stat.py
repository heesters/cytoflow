#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Gate Frequencies
----------------

Compute gate frequencies (% of parent, % of total, and event count) for
a boolean gate condition, optionally grouped by experimental conditions.

This is the standard flow-cytometry readout — equivalent to FlowJo's
"% of Parent" statistic.

.. object:: Name

    The name of the new statistic.

.. object:: Gate

    The boolean gate condition to compute frequencies for.

.. object:: Parent gates

    The gate(s) whose conjunction defines the parent population
    (denominator for **% Parent**).  Typically the gates applied before
    this one in the analysis hierarchy.  Leave empty to use all events.

.. object:: Group By

    Conditions to stratify results by (e.g. Donor, Treatment).  One row
    is produced per unique combination.

"""

from traits.api import provides, List
from traitsui.api import (View, Item, TextEditor, EnumEditor, CheckListEditor)
from pyface.api import ImageResource  # @UnresolvedImport
from envisage.api import Plugin

from ..workflow.operations import GateStatisticWorkflowOp

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, PluginHelpMixin, shared_op_traits_view


class GateStatisticHandler(OpHandler):
    operation_traits_view = \
        View(Item('name',
                  editor=TextEditor(auto_set=False,
                                    placeholder="None"),
                  tooltip="Name for this statistic (used as column header in exports)"),
             Item('gate',
                  editor=EnumEditor(name='context_handler.previous_bool_conditions_names'),
                  label="Gate",
                  tooltip="The boolean gate condition whose frequency you want to measure"),
             Item('parent',
                  editor=CheckListEditor(cols=2,
                                         name='context_handler.previous_bool_conditions_names'),
                  label='Parent\ngates',
                  style='custom',
                  tooltip="Gates that define the parent population (denominator). "
                          "Tick all gates applied before this one in the hierarchy."),
             Item('by',
                  editor=CheckListEditor(cols=2,
                                         name='context_handler.previous_conditions_names'),
                  label='Group\nBy',
                  style='custom',
                  tooltip="Experimental conditions to stratify by (e.g. Donor, Treatment)"),
             shared_op_traits_view)


@provides(IOperationPlugin)
class GateStatisticPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.gate_stat'
    operation_id = 'cytoflow.operations.gate_stat'
    view_id = None

    name = "Gate Frequencies"
    short_name = "Gate Freq"
    menu_group = "Statistics"

    def get_operation(self):
        return GateStatisticWorkflowOp()

    def get_handler(self, model, context):
        return GateStatisticHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('channel_stat')

    plugin = List(contributes_to=OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
