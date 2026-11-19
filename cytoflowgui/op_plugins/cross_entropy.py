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
Cross Entropy Test
------------------

Tests whether groups of events -- for example samples, or treated vs 
control -- look *different* in a t-SNE plot.  Deciding that by eye is 
subjective; this test puts a number and a p-value on it.

For every event, the test computes the *cross entropy* between the event's
neighborhood in the original data and in the t-SNE plot (this is the 
quantity that t-SNE minimizes.)  Then it compares the distributions of 
cross entropy between every pair of groups with a Kolmogorov-Smirnov test.

How to use it:

1. Make a t-SNE of **all the groups together** (the *tSNE* operation, with 
   nothing ticked under *Group Estimates By*.)
2. Add this operation.  Choose the **same channels and scales** and the 
   **same perplexity** as the t-SNE, and choose the two t-SNE channels as 
   the embedding.
3. Choose the condition whose values are the groups to compare, and click 
   **Estimate!**

The plot shows the cumulative distribution of cross entropy for each group
(or a dendrogram of the groups).  The results are saved in two statistics,
which you can see with a *Table* view:

- the statistic with this operation's name has one row per pair of groups:
  **D**, the largest gap between the two groups' curves (0 = identical, 
  1 = completely different; use it as a distance), the p-value **p**, 
  **p_adj** (corrected for the number of pairs), and **reject** (whether
  p_adj is below alpha);
- the statistic ending in ``_Groups`` has the number of events and the 
  median and mean cross entropy of each group.

.. note::

    With thousands of events per group, even tiny differences are 
    "significant".  Pay attention to **D** as well, and consider a strict
    threshold (the authors suggest p < 0.001 for large data sets.)
    
.. note::

    The test compares how each group's events sit *within the structure* of
    the data.  Groups that differ only in the proportions of populations
    with a similar shape can look alike to it.

.. object:: Name

    The operation name; used to name the new statistics.

.. object:: Channels

    The channels (and their scales) that the t-SNE was made from.

.. object:: Embedding X, Embedding Y

    The two channels made by the t-SNE operation.

.. object:: Compare

    The condition whose values are the groups to compare.

.. object:: Perplexity

    The perplexity the t-SNE was made with.

.. object:: Max events

    The calculation looks at all pairs of events, so at most this many 
    events (a random sample) are used.

.. object:: Correction

    The correction for testing many pairs of groups: Holm (the default), 
    Bonferroni, or none.

.. object:: Alpha

    The significance threshold.

The method is from Roca *et al.*, "A cross entropy test allows quantitative
statistical comparison of t-SNE and UMAP representations", *Cell Reports
Methods* 3:100390 (2023).  This is an independent implementation of the 
t-SNE form of the test as published; it is not the authors' code, and its
numbers may differ slightly from theirs.
'''

from natsort import natsorted

from traits.api import provides, Event, Property, List, Str
from traitsui.api import (View, Item, EnumEditor, HGroup, VGroup, TextEditor, 
                          ButtonEditor, Controller)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..editors import (SubsetListEditor, InstanceHandlerEditor, VerticalListEditor,
                       ColorTextEditor, ExtendableEnumEditor)
from ..workflow.operations import (CrossEntropyTestWorkflowOp, 
                                   CrossEntropyTestWorkflowView, tSNEChannel)
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class ChannelHandler(Controller):
    channel_view = View(HGroup(Item('channel',
                                    editor = EnumEditor(name = 'context_handler.channels')),
                               Item('scale')))


class CrossEntropyTestHandler(OpHandler):
    
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
             VGroup(Item('xembedding',
                         editor = EnumEditor(name = 'handler.channels'),
                         label = "Embedding X"),
                    Item('yembedding',
                         editor = EnumEditor(name = 'handler.channels'),
                         label = "Embedding Y"),
                    Item('by',
                         editor = ExtendableEnumEditor(name = 'context_handler.previous_conditions_names',
                                                       extra_items = {"None" : ""}),
                         label = "Compare"),
                    Item('perplexity',
                         editor = TextEditor(auto_set = False),
                         label = "Perplexity"),
                    Item('max_events',
                         editor = TextEditor(auto_set = False, evaluate = int),
                         label = "Max events"),
                    Item('correction',
                         label = "Correction"),
                    Item('alpha',
                         editor = TextEditor(auto_set = False, evaluate = float),
                         label = "Alpha"),
                    label = "Test parameters"),
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
        
        
class CrossEntropyTestViewHandler(ViewHandler):
    view_traits_view = \
        View(Item('dendrogram',
                  label = "Dendrogram"),
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
class CrossEntropyTestPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.cross_entropy'
    operation_id = 'cytoflow.operations.cross_entropy'
    view_id = 'cytoflow.view.crossentropytestview'

    name = "Cross Entropy Test"
    short_name = "Cross Ent"
    menu_group = "Computational Flow"
    
    def get_operation(self):
        return CrossEntropyTestWorkflowOp()
    
    def get_handler(self, model, context):
        if isinstance(model, CrossEntropyTestWorkflowOp):
            return CrossEntropyTestHandler(model = model, context = context)
        elif isinstance(model, CrossEntropyTestWorkflowView):
            return CrossEntropyTestViewHandler(model = model, context = context)
    
    def get_icon(self):
        return ImageResource('cross_entropy')
    
    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
