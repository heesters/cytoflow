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
Derived Parameter
-----------------

Adds a new "channel" to the workflow, computed from other channels with an
arithmetic expression (like FlowJo's "Derive Parameters" or FCS Express's
"Parameter Math"). For example:

* ``FITC_A / PE_A`` -- the ratio of two channels;
* ``log10(B1_A + 1)`` -- a log-transformed channel;
* ``(V2_A - V1_A) * 2``;
* ``arcsinh(Y2_A / 150)`` -- an arcsinh transform with a cofactor of 150;
* ``max(FSC_A, 0)``.

Refer to a channel by its name with every character that isn't a letter,
number or underscore replaced by ``_`` -- so ``FITC-A`` is ``FITC_A`` --
or put the name in backticks: ```FITC-A```. The channels you can use are
listed below the expression.

The expression may contain channel names, numbers, ``+``, ``-``, ``*``,
``/``, ``**`` (power), parentheses, and these functions:

* ``log`` (natural log; also ``ln``), ``log10``, ``log2``
* ``exp``, ``sqrt``, ``abs``
* ``arcsinh`` (also ``asinh``)
* ``min(a, b)`` and ``max(a, b)`` (the smaller or larger of two values,
  event by event)

Events for which the result isn't a finite number -- eg, division by zero,
or the log of a zero or negative value -- are **removed** from the
experiment, and a warning says how many were removed. To keep them, avoid
the problem in the expression: eg, ``log10(max(B1_A, 1))`` instead of
``log10(B1_A)``.

.. object:: Name

    The name of the new channel.

.. object:: Expression

    The arithmetic expression that computes the new channel.

.. object:: Available channels

    The channels that you can use in the expression.
'''

from traits.api import provides, List, Property, Str
from traitsui.api import View, Item, TextEditor
from pyface.api import ImageResource  # @UnresolvedImport
from envisage.api import Plugin

from natsort import natsorted

import cytoflow.utility as util
from cytoflow.operations.derived import FUNCTIONS

from ..workflow.operations import DerivedParameterWorkflowOp

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, PluginHelpMixin, shared_op_traits_view


def available_channels_text(channels):
    """
    Format a list of channels as the names to use in a derived parameter
    expression.
    """
    if not channels:
        return "(none)"

    lines = []
    for channel in natsorted(channels):
        ident = util.sanitize_identifier(channel)
        if ident == channel:
            lines.append(channel)
        else:
            lines.append("{}  (or `{}`)".format(ident, channel))

    # a few functions per line, so the panel doesn't get too wide
    functions = list(FUNCTIONS.keys())
    lines.append("")
    lines.append("Functions:")
    for i in range(0, len(functions), 4):
        lines.append(", ".join(functions[i:i + 4]))

    return "\n".join(lines)


class DerivedParameterHandler(OpHandler):

    available_names = Property(Str, observe = "context.previous_wi.channels")

    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False, placeholder = "None")),
             Item('expression',
                  editor = TextEditor(auto_set = False,
                                      placeholder = "eg FITC_A / PE_A")),
             Item('handler.available_names',
                  label = "Available\nchannels",
                  style = 'readonly',
                  resizable = True),
             shared_op_traits_view)

    def _get_available_names(self):
        if self.context and self.context.previous_wi:
            return available_channels_text(self.context.previous_wi.channels)
        return ""


@provides(IOperationPlugin)
class DerivedParameterPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.derived'
    operation_id = 'cytoflow.operations.derived'
    view_id = None

    short_name = "Derived"
    name = "Derived Parameter"
    menu_group = "Preprocessing"

    def get_operation(self):
        return DerivedParameterWorkflowOp()

    def get_handler(self, model, context):
        return DerivedParameterHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('derived')

    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
