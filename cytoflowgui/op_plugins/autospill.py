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

'''
AutoSpill Compensation
----------------------

Compute the spillover (compensation) matrix from single-color controls with
the *AutoSpill* algorithm, then compensate the data.

Unlike traditional compensation, you **don't have to gate positive and 
negative populations** in each control.  AutoSpill fits a robust straight 
line through *all* of each control's events to measure how much of each 
color spills into every other channel, then repeatedly compensates the 
controls and corrects the matrix until no spillover is left.  This is 
usually more accurate, and it copes with controls that have no clear 
negative population.

How to use it:

1. Before this operation, gate out debris (eg, a polygon gate on FSC-A vs 
   SSC-A).  The gate is applied to the controls, too.
2. Click **Add a channel** once for every fluorescence channel, choose the
   channel, and choose the FCS file of its single-color control.  Every 
   channel needs a control.
3. Under **Subset**, choose the gate from step 1.
4. Click **Estimate!**  The spillover matrix appears in the table, where 
   you can inspect (and, if you must, edit) it, and the plot shows it as a
   heat map.

.. note::

    Use only one compensation operation (this one, *Compensation (from 
    FCS)*, or *Linear Bleedthrough*) in a workflow -- otherwise the data 
    are compensated twice.

.. object:: Channels

    The channels to compensate, and the FCS file with each one's 
    single-color control.

.. object:: Spillover Matrix

    The computed matrix: the fraction of each channel's ("From") signal that
    is detected in each other channel ("To").

.. object:: Subset

    Which events of the controls to use.

.. object:: Arcsinh cofactor

    The final refinement is done on ``arcsinh(x / cofactor)`` transformed
    data, so that errors in dim events are corrected too.  The default 
    (150) suits conventional cytometers.

.. object:: Trim

    The fraction of the most extreme events left out of each fit (default
    0.01.)

.. object:: Max events

    Use at most this many events from each control (0 = all of them.)

AutoSpill is described in Roca *et al.*, "AutoSpill is a principled 
framework that simplifies the analysis of multichromatic flow cytometry 
data", *Nature Communications* 12:2890 (2021).  This is a Python port of 
the authors' R package.
'''

from natsort import natsorted

from traits.api import provides, Event, Property, List, Str
from traitsui.api import (View, Item, VGroup, ButtonEditor, FileEditor, HGroup, 
                          EnumEditor, Controller, TableEditor, ObjectColumn, 
                          TextEditor)

from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..editors import InstanceHandlerEditor, VerticalListEditor, SubsetListEditor
from ..workflow.operations import (AutoSpillWorkflowOp, CompensationWorkflowView,
                                   BleedthroughChannel, AutoSpillQCWorkflowView)
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin
from .compensation import CompensationViewHandler
from ..view_plugins import ViewHandler


class AutoSpillChannelHandler(Controller):
    channel_view = View(HGroup(Item('channel',
                                    editor = EnumEditor(name = 'context_handler.channels')),
                               Item('file',
                                    editor = FileEditor(dialog_style = 'open'),
                                    show_label = False)))
    
    
class AutoSpillHandler(OpHandler):
    add_channel = Event
    remove_channel = Event
        
    channels = Property(List(Str), observe = 'context.previous_wi.channels')
     
    operation_traits_view = \
        View(
            VGroup(
                VGroup(
                        Item('channels_list',
                             editor = VerticalListEditor(editor = InstanceHandlerEditor(view = 'channel_view',
                                                                                        handler_factory = AutoSpillChannelHandler),
                                                         style = 'custom',
                                                         mutable = False)),
                        Item('handler.add_channel',
                             editor = ButtonEditor(value = True,
                                                   label = "Add a channel")),
                        Item('handler.remove_channel',
                             editor = ButtonEditor(value = True,
                                                   label = "Remove a channel")),
                        label = "Channels",
                        show_labels = False),
                VGroup(
                    Item('spillover_list',
                         editor = TableEditor(columns = [ObjectColumn(name = 'from_channel',
                                                                      label = "From",
                                                                      editable = False),
                                                         ObjectColumn(name = 'to_channel',
                                                                      label = "To",
                                                                      editable = False),
                                                         ObjectColumn(name = 'spillover',
                                                                      label = "Spillover",
                                                                      editor = TextEditor(evaluate = float,
                                                                                          auto_set = False))],
                                              rows = 8,
                                              editable = True,
                                              sortable = False,
                                              auto_size = True,
                                              configurable = False),
                         show_label = False),
                    label = "Spillover Matrix",
                    show_border = False),
                VGroup(Item('subset_list',
                            show_label = False,
                            editor = SubsetListEditor(conditions = "context_handler.previous_conditions",
                                                      editor = InstanceHandlerEditor(view = 'subset_view',
                                                                                     handler_factory = subset_handler_factory))),
                            label = "Subset",
                            show_border = False,
                            show_labels = False),
                Item('cofactor',
                     editor = TextEditor(auto_set = False, evaluate = float),
                     label = "Arcsinh\ncofactor"),
                Item('trim',
                     editor = TextEditor(auto_set = False, evaluate = float),
                     label = "Trim"),
                Item('max_events',
                     editor = TextEditor(auto_set = False, evaluate = int),
                     label = "Max events"),
                Item('do_estimate',
                     editor = ButtonEditor(value = True,
                                           label = "Estimate!"),
                     show_label = False),
                shared_op_traits_view))
        
    # MAGIC: called when add_channel is set
    def _add_channel_fired(self):
        self.model.channels_list.append(BleedthroughChannel())
        
    # MAGIC: called when remove_channel is set
    def _remove_channel_fired(self):
        if self.model.channels_list:
            self.model.channels_list.pop()

    # MAGIC: returns the value of the 'channels' property
    def _get_channels(self):
        try:
            return natsorted(self.context.previous_wi.channels)
        except AttributeError:
            return []


class AutoSpillQCViewHandler(ViewHandler):
    view_traits_view = View(
        Item('context.view_warning',
             resizable = True,
             visible_when = 'context.view_warning',
             editor = ColorTextEditor(foreground_color = "#000000",
                                     background_color = "#ffff99")),
        Item('context.view_error',
             resizable = True,
             visible_when = 'context.view_error',
             editor = ColorTextEditor(foreground_color = "#000000",
                                      background_color = "#ff9191")))

    view_params_view = View()


@provides(IOperationPlugin)
class AutoSpillPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.autospill'
    operation_id = 'cytoflow.operations.autospill'
    view_id = 'cytoflow.view.compensationdiagnostic'

    name = "AutoSpill Compensation"
    short_name = "AutoSpill"
    menu_group = "Preprocessing"

    def get_operation(self):
        return AutoSpillWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, AutoSpillWorkflowOp):
            return AutoSpillHandler(model = model, context = context)
        elif isinstance(model, CompensationWorkflowView):
            return CompensationViewHandler(model = model, context = context)
        elif isinstance(model, AutoSpillQCWorkflowView):
            return AutoSpillQCViewHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('autospill')

    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
