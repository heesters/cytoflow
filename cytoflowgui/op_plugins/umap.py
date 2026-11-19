#!/usr/bin/env python3.8
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
Uniform Manifold Approximation and Projection
---------------------------------------------

Use Uniform Manifold Approximation and Projection (UMAP) to embed
multivariate data in a small number of dimensions while maintaining its
underlying structure. Like tSNE, UMAP is a non-linear method; it is
generally faster than tSNE and tends to preserve more of the global
structure of the data.

Creates new "channels" named ``{name}_1``, ``{name}_2``, ... ``{name}_n``,
where ``name`` is the **Name** attribute and ``n`` is the number of
**Components**.

The model is fit to a random **Sample** of the events; then *all* the
events are projected into the embedding.

The same embedding may not be appropriate for different subsets of the data set.
If this is the case, you can use the **By** attribute to specify
metadata by which to aggregate the data before estimating (and applying) a
model.  The UMAP parameters, such as the number of neighbors, are the same
across each subset, though.

.. object:: Name

    The operation name; determines the name of the new columns.

.. object:: Channels

    The channels to apply the embedding to.

.. object:: Scale

    Re-scale the data in the specified channels before fitting.

.. object:: Components

    How many dimensions to embed the data into? Defaults to 2.

.. object:: Neighbors

    The size of the local neighborhood UMAP looks at when learning the
    structure of the data. Smaller values emphasize local structure; larger
    values emphasize global structure. Defaults to 15.

.. object:: Min. distance

    How tightly UMAP is allowed to pack points together. Smaller values give
    tighter clusters; larger values spread the points out more evenly.
    Defaults to 0.1.

.. object:: Metric

    How should we measure "distance"? Euclidean distance makes sense if the
    number of dimensions (ie, channels) is small, but as the number of
    dimensions increases, maybe try ``cosine`` or one of the others.

.. object:: Random seed

    If set, makes the embedding reproducible (but slower, because UMAP then
    runs on a single processor.) Leave it blank for a different random
    embedding each time.

.. object:: Sample

    What proportion of the data set to use for training? Defaults to 1%
    of the dataset to help with runtime.

.. object:: By

    A list of metadata attributes to aggregate the data before estimating
    the model.  For example, if the experiment has two pieces of metadata,
    ``Time`` and ``Dox``, setting **By** to ``["Time", "Dox"]`` will
    fit the model separately to each subset of the data with a unique
    combination of ``Time`` and ``Dox``.

.. note::

    UMAP requires the ``umap-learn`` package.

.. plot::
   :include-source: False

    import cytoflow as flow
    import_op = flow.ImportOp()
    import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
                                 conditions = {'Dox' : 10.0}),
                       flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
                                 conditions = {'Dox' : 1.0})]
    import_op.conditions = {'Dox' : 'float'}
    ex = import_op.apply()


    umap_op = flow.UMAPOp(name = 'UMAP',
                          channels = ['V2-A', 'Y2-A'],
                          scale = {'V2-A' : 'log',
                                   'Y2-A' : 'log'},
                          random_state = 0)
    umap_op.estimate(ex)
    ex2 = umap_op.apply(ex)
    flow.ScatterplotView(xchannel = "UMAP_1",
                         ychannel = "UMAP_2",
                         huefacet = "Dox").plot(ex2)


'''
from natsort import natsorted

from traits.api import provides, Event, Property, List, Str
from traitsui.api import (View, Item, EnumEditor, HGroup, VGroup, TextEditor,
                          CheckListEditor, ButtonEditor, Controller)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..editors import SubsetListEditor, InstanceHandlerEditor, VerticalListEditor
from ..workflow.operations import UMAPWorkflowOp, UMAPChannel
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class ChannelHandler(Controller):
    channel_view = View(HGroup(Item('channel',
                                    editor = EnumEditor(name = 'context_handler.channels')),
                               Item('scale')))

class UMAPHandler(OpHandler):

    add_channel = Event
    remove_channel = Event
    channels = Property(List(Str), observe = 'context.channels')

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
             show_labels = False),
             VGroup(Item('n_components',
                         editor = TextEditor(auto_set = False),
                         label = "Components"),
                    Item('n_neighbors',
                         editor = TextEditor(auto_set = False),
                         label = "Neighbors"),
                    Item('min_dist',
                         editor = TextEditor(auto_set = False),
                         label = "Min. distance"),
                    Item('metric'),
                    Item('random_state',
                         editor = TextEditor(auto_set = False,
                                             format_func = lambda x: "" if x is None else str(x),
                                             placeholder = "None"),
                         label = "Random seed"),
                    Item('sample',
                         editor = TextEditor(auto_set = False,
                                             evaluate = float,
                                             format_func = lambda x: "" if x is None else str(x)),
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
             Item('status',
                  width = 1.0,
                  style = "readonly"),
             shared_op_traits_view)

    # MAGIC: called when add_channel is set
    def _add_channel_fired(self):
        self.model.channels_list.append(UMAPChannel())

    def _remove_channel_fired(self):
        if self.model.channels_list:
            self.model.channels_list.pop()

    def _get_channels(self):
        if self.context and self.context.channels:
            return natsorted(self.context.channels)
        else:
            return []


@provides(IOperationPlugin)
class UMAPPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.umap'
    operation_id = 'cytoflow.operations.umap'
    view_id = None

    name = "Uniform Manifold Approximation and Projection"
    short_name = "UMAP"
    menu_group = "Computational Flow"

    def get_operation(self):
        return UMAPWorkflowOp()

    def get_handler(self, model, context):
        return UMAPHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('umap')

    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
