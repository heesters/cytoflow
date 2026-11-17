#!/usr/bin/env python3.8
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

'''
Sample Labels
-------------

Assign human-readable labels to values of an existing condition (e.g. ``File``).

.. object:: Condition to label

    The condition whose values you want to relabel, e.g. ``File``.

.. object:: Label mapping

    A table of original value → new label pairs.  Any value you do **not**
    remap keeps its original value.

.. object:: New condition name (optional)

    If set, the relabelled values are stored in a **new** condition with this
    name and the original condition is kept unchanged.  Leave blank to replace
    the original condition in-place.
'''

from traits.api import (provides, List, Str, observe, Instance, HasStrictTraits)
from traitsui.api import (View, Item, VGroup, EnumEditor, TextEditor,
                          TableEditor, ObjectColumn)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..editors import ColorTextEditor
from ..workflow.operations import SampleLabelsWorkflowOp
from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class _LabelRow(HasStrictTraits):
    """One row in the label mapping table."""
    original = Str
    label = Str


class SampleLabelsHandler(OpHandler):
    """Handler for SampleLabelsWorkflowOp."""

    # Synced list of rows for the TableEditor
    labels_list = List(Instance(_LabelRow))

    @observe('context.previous_wi.result', post_init=True)
    def _refresh_labels_list(self, _=None):
        """Rebuild rows when the upstream result changes."""
        op = self.model
        if (self.context is not None
                and self.context.previous_wi is not None
                and self.context.previous_wi.result is not None
                and op.condition
                and op.condition in self.context.previous_wi.result.conditions):
            result = self.context.previous_wi.result
            vals = sorted(result.data[op.condition].unique().astype(str))
            # Preserve any existing labels
            new_rows = []
            for v in vals:
                existing = op.labels.get(v, '')
                new_rows.append(_LabelRow(original=v, label=existing))
            self.labels_list = new_rows
        else:
            self.labels_list = []

    @observe('model.condition', post_init=True)
    def _condition_changed(self, _=None):
        self._refresh_labels_list()

    @observe('labels_list:items.label')
    def _labels_list_changed(self, _=None):
        """Sync the UI list back to the op's labels dict."""
        new_labels = {}
        for row in self.labels_list:
            if row.label and row.label != row.original:
                new_labels[row.original] = row.label
        # only update if changed to avoid infinite loops
        if new_labels != dict(self.model.labels):
            self.model.labels = new_labels

    operation_traits_view = \
        View(
            VGroup(
                Item('condition',
                     editor=EnumEditor(name='context_handler.previous_conditions_names'),
                     label='Condition'),
                Item('handler.labels_list',
                     editor=TableEditor(
                         columns=[
                             ObjectColumn(name='original',
                                          label='Current value',
                                          editable=False,
                                          width=0.5),
                             ObjectColumn(name='label',
                                          label='New label',
                                          editable=True,
                                          width=0.5,
                                          editor=TextEditor(auto_set=False)),
                         ],
                         rows=8,
                         editable=True,
                         auto_size=True,
                         configurable=False,
                     ),
                     show_label=False),
                Item('new_condition',
                     editor=TextEditor(auto_set=False, placeholder='Leave blank to replace in-place'),
                     label='New condition name'),
                show_border=False,
            ),
            shared_op_traits_view)


@provides(IOperationPlugin)
class SampleLabelsPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.sample_labels'
    operation_id = 'cytoflow.operations.sample_labels'

    short_name = name = "Sample Labels"
    menu_group = "TOP"

    def get_operation(self):
        return SampleLabelsWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, SampleLabelsWorkflowOp):
            return SampleLabelsHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('sample_labels')

    plugin = List(contributes_to=OP_PLUGIN_EXT)

    def _plugin_default(self):
        return [self]
