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
Proliferation
-------------

Dye-dilution proliferation analysis, for cells labelled with CFSE, 
CellTrace Violet or a similar dye.

Every time a cell divides, each daughter gets half of the dye.  On a 
logarithmic axis, the generations are therefore evenly spaced peaks.  This 
operation fits a curve to each generation and computes the usual 
proliferation statistics.

How to use it:

1. Gate on single, live cells first, and choose that gate under **Subset**.
2. Choose the dye's channel.
3. Tell the operation where the *undivided* cells are.  The most reliable 
   way is an unstimulated control: choose it under **Unstimulated 
   control**.  Otherwise, type the peak's position into **Undivided 
   peak**, or leave both empty to use the brightest peak in each sample -- 
   but if nearly all the cells in a sample have divided, that peak is *not* 
   generation 0.
4. Tick the condition that identifies your samples under **Group Estimates 
   By**, and click **Estimate!**

The plot shows each sample's histogram with the fitted generations 
(numbered: 0 is undivided.)  The results are in a statistic with this 
operation's name (see it with a **Table** view), and a new condition with 
each cell's generation, which you can use to color or subset other plots.

The statistic's columns:

- **Generation_0**, **Generation_1**, ...: the percentage of cells in each 
  generation.
- **Percent_divided**: the percentage of the *starting* cells that divided
  at least once.  (Each cell that divided ``g`` times is one of 2 :sup:`g`
  descendants, so the number of starting cells is worked out by dividing 
  each generation by 2 :sup:`g`.)
- **Division_index**: the average number of divisions of *all* the starting
  cells, including the ones that never divided.
- **Proliferation_index**: the average number of divisions of the cells 
  that *did* divide.
- **Expansion_index**: how many-fold the whole culture expanded.
- **Replication_index**: how many-fold the responding cells expanded.
- **Undivided_peak**: the position of generation 0.

These are the same definitions that FlowJo uses.

.. object:: Name

    The operation name; names the new statistic and condition.

.. object:: Channel

    The channel with the proliferation dye.

.. object:: Generations

    The largest number of divisions to look for (default 8.)

.. object:: Undivided peak

    The position of the undivided cells, in the channel's units.  Leave it
    empty to find it from the control or the data.

.. object:: Fit spacing

    The generations are exactly a factor of 2 apart.  Tick this to fit the
    spacing instead (within 15%.)

.. object:: Group Estimates By

    Fit separately for each value of these conditions (eg, each sample.)

.. object:: Unstimulated control

    The events of a sample whose cells have not divided.  Its peak sets the 
    position of generation 0 for all the samples.

.. object:: Subset

    Only use these events.
'''

from traits.api import provides, List
from traitsui.api import (View, Item, EnumEditor, VGroup, TextEditor, 
                          CheckListEditor, ButtonEditor)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..editors import SubsetListEditor, ColorTextEditor, InstanceHandlerEditor
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin
from ..workflow.operations import ProliferationWorkflowOp, ProliferationWorkflowView


class ProliferationHandler(OpHandler):
    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False, placeholder = "None")),
             Item('channel',
                  editor = EnumEditor(name = 'context_handler.previous_channels'),
                  label = "Channel"),
             Item('num_generations',
                  editor = TextEditor(auto_set = False),
                  label = "Generations"),
             Item('undivided_peak',
                  editor = TextEditor(auto_set = False,
                                      placeholder = "From control / data",
                                      format_func = lambda x: "" if x is None else str(x)),
                  label = "Undivided\npeak"),
             Item('fit_spacing',
                  label = "Fit\nspacing"),
             Item('by',
                  editor = CheckListEditor(cols = 2,
                                           name = 'context_handler.previous_conditions_names'),
                  label = 'Group\nEstimates\nBy',
                  style = 'custom'),
             VGroup(Item('control_subset_list',
                         show_label = False,
                         editor = SubsetListEditor(conditions = "context_handler.previous_conditions",
                                                   editor = InstanceHandlerEditor(view = 'subset_view',
                                                                                  handler_factory = subset_handler_factory))),
                    label = "Unstimulated control",
                    show_border = False,
                    show_labels = False),
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


class ProliferationViewHandler(ViewHandler):
    view_traits_view = \
        View(VGroup(Item('subset_list',
                         show_label = False,
                         editor = SubsetListEditor(conditions = "context_handler.conditions",
                                                   editor = InstanceHandlerEditor(view = 'subset_view',
                                                                                  handler_factory = subset_handler_factory))),
                    label = "Subset",
                    show_border = False,
                    show_labels = False),
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
class ProliferationPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.proliferation'
    operation_id = 'cytoflow.operations.proliferation'
    view_id = 'cytoflow.view.proliferationview'

    name = "Proliferation"
    short_name = "Proliferation"
    menu_group = "Clustering"
    
    def get_operation(self):
        return ProliferationWorkflowOp()
    
    def get_handler(self, model, context):
        if isinstance(model, ProliferationWorkflowOp):
            return ProliferationHandler(model = model, context = context)
        elif isinstance(model, ProliferationWorkflowView):
            return ProliferationViewHandler(model = model, context = context)
    
    def get_icon(self):
        return ImageResource('proliferation')
    
    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
