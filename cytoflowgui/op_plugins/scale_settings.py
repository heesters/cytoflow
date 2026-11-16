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
Axis Transforms
---------------

Set the parameters of the axis transforms ("scales") for individual 
channels -- like "customize axis" in FlowJo or FCS Express.  Every plot and 
gate after this operation that uses the **logicle** or **arcsinh** scale on
that channel uses these settings.  The data themselves aren't changed.

Click **Add a channel**, choose the channel, and fill in the parameters you
want to change; leave the others empty to use the defaults.

.. object:: Logicle W

    The width of the linear region around 0, in decades.  By default it is
    estimated from the data.  Increase it if the negative population is 
    split into two "lumps"; decrease it if populations near 0 are squashed
    together.  Typical values: 0.5 - 1.5.
    
.. object:: Logicle M

    The width of the whole axis, in decades (default 4.5).
    
.. object:: Logicle A

    Additional decades of negative values to show (default 0).
    
.. object:: Arcsinh cofactor

    The arcsinh scale shows ``arcsinh(x / cofactor)``: values much smaller
    than the cofactor are shown linearly, much larger ones logarithmically.
    Typical values: 5 for mass cytometry, 150 for conventional flow, and 
    several hundred to several thousand for some spectral cytometers.
    (Default 150.)
'''

from natsort import natsorted

from traits.api import provides, List, Str, Event, Property
from traitsui.api import (View, Item, TextEditor, EnumEditor, HGroup, VGroup,
                          ButtonEditor, Controller)
from pyface.api import ImageResource  # @UnresolvedImport
from envisage.api import Plugin

from ..workflow.operations import ScaleSettingsWorkflowOp
from ..workflow.operations.scale_settings import ChannelScaleSettings
from ..editors import VerticalListEditor, InstanceHandlerEditor
from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT 
from .op_plugin_base import OpHandler, PluginHelpMixin, shared_op_traits_view


def _float_editor():
    return TextEditor(auto_set = False,
                      placeholder = "default",
                      format_func = lambda x: "" if x is None else str(x))


class ChannelScaleHandler(Controller):
    channel_view = View(VGroup(Item('channel',
                                    editor = EnumEditor(name = 'context_handler.channels')),
                               HGroup(Item('logicle_W', label = "Logicle W", 
                                           editor = _float_editor()),
                                      Item('logicle_M', label = "M",
                                           editor = _float_editor()),
                                      Item('logicle_A', label = "A",
                                           editor = _float_editor())),
                               Item('arcsinh_cofactor', label = "Arcsinh\ncofactor",
                                    editor = _float_editor()),
                               show_border = True))


class ScaleSettingsHandler(OpHandler):
    
    add_channel = Event
    remove_channel = Event
    channels = Property(List(Str), observe = 'context.previous_wi.channels')
    
    operation_traits_view = \
        View(VGroup(Item('settings_list',
                         editor = VerticalListEditor(editor = InstanceHandlerEditor(view = 'channel_view',
                                                                                    handler_factory = ChannelScaleHandler),
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
             shared_op_traits_view) 
    
    # MAGIC: called when add_channel is set
    def _add_channel_fired(self):
        self.model.settings_list.append(ChannelScaleSettings())
        
    def _remove_channel_fired(self):
        if self.model.settings_list:
            self.model.settings_list.pop()
            
    def _get_channels(self):
        try:
            return natsorted(self.context.previous_wi.channels)
        except AttributeError:
            return []


@provides(IOperationPlugin)
class ScaleSettingsPlugin(Plugin, PluginHelpMixin):
    
    id = 'cytoflowgui.op_plugins.scale_settings'
    operation_id = 'cytoflow.operations.scale_settings'
    view_id = None
    
    name = "Axis Transforms"
    short_name = "Axis\nTransforms"
    menu_group = "Preprocessing"
    
    def get_operation(self):
        return ScaleSettingsWorkflowOp()
    
    def get_handler(self, model, context):
        return ScaleSettingsHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('scale_settings')

    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
