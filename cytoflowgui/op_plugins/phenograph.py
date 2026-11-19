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
Graph Clustering (PhenoGraph)
-----------------------------

Finds populations automatically with *graph clustering*, the approach of
PhenoGraph: every event is connected to the events most similar to it, and
the resulting network is split into densely connected communities.  Each 
community is a cluster.

Unlike k-means, you don't have to say how many clusters there are, clusters
can have any shape, and you can use as many channels as you like -- so it 
is well suited to finding cell populations in a multi-color panel.

Add the channels to cluster on (with a scale for each -- usually 
*logicle*), then click **Estimate!**  The operation adds a new condition 
with this operation's name whose values are the cluster names (the largest
cluster is number 1), and a statistic with the clusters' centers and 
proportions.  The plot shows the first two channels, colored by cluster.

.. object:: Name

    The operation name; also the name of the new condition.

.. object:: Channels

    The channels to cluster on, and the scale for each.

.. object:: Neighbors

    How many nearest neighbors each event is connected to (default 30).
    Fewer neighbors find smaller, finer clusters.

.. object:: Resolution

    Larger than 1 gives more (smaller) clusters; smaller than 1 gives 
    fewer (larger) ones.

.. object:: Method

    How communities are found: *leiden* (the default) or *louvain* (the 
    method of the original PhenoGraph.)

.. object:: Sample

    The fraction of events used to build the network (default: all of 
    them.)  Building it is slow for more than ~100,000 events; with a 
    smaller sample, the remaining events are each given the cluster most 
    common among their nearest sampled neighbors.

.. object:: Group Estimates By

    Cluster separately for each value of these conditions.

.. object:: Subset

    Only use these events to find the clusters.

The method is from Levine *et al.*, "Data-driven phenotypic dissection of 
AML reveals progenitor-like cells that correlate with prognosis", *Cell* 
162:184 (2015).
'''

from natsort import natsorted

from traits.api import provides, Event, Property, List, Str
from traitsui.api import (View, Item, EnumEditor, HGroup, VGroup, TextEditor, 
                          CheckListEditor, ButtonEditor, Controller)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..view_plugins.scatterplot import ScatterplotParamsHandler
from ..editors import (SubsetListEditor, ColorTextEditor, InstanceHandlerEditor,
                       VerticalListEditor)
from ..workflow.operations import PhenographWorkflowOp, PhenographWorkflowView, tSNEChannel
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class ChannelHandler(Controller):
    channel_view = View(HGroup(Item('channel',
                                    editor = EnumEditor(name = 'context_handler.channels')),
                               Item('scale')))


class PhenographHandler(OpHandler):
    
    add_channel = Event
    remove_channel = Event
    channels = Property(List(Str), observe = 'context.previous_wi.channels')
    
    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False, placeholder = "None")),
             VGroup(Item('channels_list',
                         editor = VerticalListEditor(editor = InstanceHandlerEditor(view = 'channel_view',
                                                                                    handler_factory = ChannelHandler),
                                                     style = 'custom',
                                                     mutable = False)),
             Item('handler.add_channel',
                  editor = ButtonEditor(value = True,
                                        label = "Add a channel"),
                  show_label = False),
             Item('handler.remove_channel',
                  editor = ButtonEditor(value = True,
                                        label = "Remove a channel")),
             label = "Channels",
             show_labels = False),
             VGroup(Item('num_neighbors',
                         editor = TextEditor(auto_set = False),
                         label = "Neighbors"),
                    Item('resolution',
                         editor = TextEditor(auto_set = False),
                         label = "Resolution"),
                    Item('method',
                         label = "Method"),
                    Item('sample',
                         editor = TextEditor(auto_set = False, evaluate = float),
                         label = "Sample"),
                    Item('by',
                         editor = CheckListEditor(cols = 2,
                                                  name = 'context_handler.previous_conditions_names'),
                         label = 'Group\nEstimates\nBy',
                         style = 'custom'),
                    label = "Estimate parameters"),
             VGroup(Item('subset_list',
                         show_label = False,
                         editor = SubsetListEditor(conditions = "context_handler.previous_conditions",
                                                   editor = InstanceHandlerEditor(view = 'subset_view',
                                                                                  handler_factory = subset_handler_factory))),
                    label = "Subset",
                    show_border = False,
                    show_labels = False),
             Item('do_estimate',
                  editor = ButtonEditor(value = True,
                                        label = "Estimate"),
                  show_label = False),
             shared_op_traits_view)

    # MAGIC: called when add_channel is set
    def _add_channel_fired(self):
        self.model.channels_list.append(tSNEChannel())
        
    def _remove_channel_fired(self):
        if self.model.channels_list:
            self.model.channels_list.pop()   
            
    def _get_channels(self):
        try:
            return natsorted(self.context.previous_wi.channels)
        except AttributeError:
            return []
        
        
class PhenographViewHandler(ViewHandler):
    view_traits_view = \
        View(Item('context.view_warning',
                  resizable = True,
                  visible_when = 'context.view_warning',
                  editor = ColorTextEditor(foreground_color = "#000000",
                                          background_color = "#ffff99")),
             Item('context.view_error',
                  resizable = True,
                  visible_when = 'context.view_error',
                  editor = ColorTextEditor(foreground_color = "#000000",
                                           background_color = "#ff9191")))

    view_params_view = \
        View(Item('plot_params',
                  editor = InstanceHandlerEditor(view = 'view_params_view',
                                                 handler_factory = ScatterplotParamsHandler),
                  style = 'custom',
                  show_label = False))
        

@provides(IOperationPlugin)
class PhenographPlugin(Plugin, PluginHelpMixin):
    
    id = 'cytoflowgui.op_plugins.phenograph'
    operation_id = 'cytoflow.operations.phenograph'
    view_id = 'cytoflowgui.op_plugins.phenograph'

    name = "Graph Clustering (PhenoGraph)"
    short_name = "Graph\nClustering"
    menu_group = "Clustering"
    
    def get_operation(self):
        return PhenographWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, PhenographWorkflowOp):
            return PhenographHandler(model = model, context = context)
        elif isinstance(model, PhenographWorkflowView):
            return PhenographViewHandler(model = model, context = context)
    
    def get_icon(self):
        return ImageResource('phenograph')
    
    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
