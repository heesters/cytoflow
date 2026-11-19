#!/usr/bin/env python3.11
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2025
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

'''
Boolean Gate
------------

Combine existing gates into a new gate using ``and``, ``or`` and ``not``
(like FlowJo's "Boolean gates" or FCS Express's "Combination gates").

The new gate is defined by a logical expression of the names of other
gates. For example, if you have a threshold gate named ``CD4_pos`` and a
polygon gate named ``Dead``, then

* ``CD4_pos and not Dead`` selects live CD4+ events;
* ``CD4_pos or CD8_pos`` selects events in either gate;
* ``not (Polygon1 or Range2)`` selects events in neither gate;
* ``Quad == 'Quad_2'`` selects the upper-right quadrant of a quad gate
  named ``Quad``. Categorical gates (quad gates, gate categories,
  clustering) are used by comparing them to one of their values with
  ``==`` or ``!=``;
* ``Y2_high == False`` is the same as ``not Y2_high``;
* ``CD4_pos and Dox > 1`` uses a numeric experimental condition.

The expression may only contain gate (condition) names, ``and``, ``or``,
``not``, parentheses, and comparisons of a condition to a constant
(``==``, ``!=``, and for numbers ``<``, ``<=``, ``>``, ``>=``). ``&``,
``|`` and ``~`` also work, but put comparisons in parentheses when you use
them: ``Gate1 & (Quad == 'Quad_1')``.

The gates you can use are listed below the expression.

The result is a new gate (a boolean condition) that is ``True`` for events
that satisfy the expression. Use it like any other gate -- eg, as the
``Color facet`` of a scatter plot, or in a subset.

.. object:: Name

    The name of the new gate.

.. object:: Expression

    The logical expression that defines the new gate.

.. object:: Available gates

    The gates (and other conditions) that you can use in the expression.
'''

from traits.api import provides, List, Property, Str
from traitsui.api import View, Item, TextEditor
from pyface.api import ImageResource  # @UnresolvedImport
from envisage.api import Plugin

from natsort import natsorted

from ..workflow.operations import BooleanGateWorkflowOp

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, PluginHelpMixin, shared_op_traits_view


def available_conditions_text(conditions):
    """
    Format the conditions in a `WorkflowItem` (a dict of name ->
    `pandas.Series` of unique values) as a list of names that can be used in
    a boolean gate expression.
    """
    if not conditions:
        return "(none -- make some gates first)"

    lines = []
    others = []
    for name in natsorted(conditions.keys()):
        values = conditions[name]
        if values.dtype.kind == 'b':
            lines.append(name)
        elif values.dtype.kind in 'iuf':
            others.append("{} (number)".format(name))
        else:
            vals = natsorted([str(v) for v in values])
            shown = ", ".join("'{}'".format(v) for v in vals[:6])
            if len(vals) > 6:
                shown += ", ..."
            others.append("{} == {}".format(name, shown))

    if not lines:
        lines.append("(no gates yet)")

    return "\n".join(lines + others)


class BooleanGateHandler(OpHandler):

    available_names = Property(Str, observe = "context.previous_wi.conditions")

    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False, placeholder = "None")),
             Item('expression',
                  editor = TextEditor(auto_set = False,
                                      placeholder = "eg Gate1 and not Gate2")),
             Item('handler.available_names',
                  label = "Available\ngates",
                  style = 'readonly',
                  resizable = True),
             shared_op_traits_view)

    def _get_available_names(self):
        if self.context and self.context.previous_wi:
            return available_conditions_text(self.context.previous_wi.conditions)
        return ""


@provides(IOperationPlugin)
class BooleanGatePlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.boolean_gate'
    operation_id = 'cytoflow.operations.boolean_gate'
    view_id = None

    short_name = "Boolean"
    name = "Boolean Gate"
    menu_group = "Gates"

    def get_operation(self):
        return BooleanGateWorkflowOp()

    def get_handler(self, model, context):
        return BooleanGateHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('boolean_gate')

    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
