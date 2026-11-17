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
Cell Cycle
----------

Cell-cycle analysis from a DNA stain (propidium iodide, DAPI, Hoechst, ...): 
the percentage of cells in G0/G1, S and G2/M.

Cells in G2 or mitosis have twice the DNA of cells in G0/G1, and cells in S
phase are in between.  This operation fits a model to the DNA histogram: a 
peak for G0/G1, a peak at about twice its position for G2/M, and a broad 
S phase between them, all blurred by the same measurement error (a 
"Dean-Jett-Fox" type model.)

How to use it:

1. Use the DNA stain's **area** channel, on a **linear** scale.
2. **Gate out doublets first** -- for example with a polygon gate on the 
   DNA stain's area vs its height or width -- and choose that gate under 
   **Subset**.  Two G1 cells stuck together look exactly like one G2 cell.
3. Tick the condition that identifies your samples under **Group Estimates 
   By**, and click **Estimate!**

The plot shows each sample's histogram with the three fitted phases.  The 
results are in a statistic with this operation's name (see it with a 
**Table** view):

- **G1**, **S**, **G2M**: the percentage of the cycling cells in each phase 
  (they add up to 100.)
- **SubG1**: the percentage of *all* events with less DNA than G0/G1 -- 
  debris and apoptotic cells.
- **AboveG2M**: the percentage of all events with more DNA than G2/M -- 
  aggregates and polyploid cells.
- **G1_mean**, **G1_CV** (in percent; a good DNA stain gives less than 
  about 5) and **G2_G1_ratio** (should be close to 2.)
- **RMSE**: how well the model fits (smaller is better; above about 0.05, 
  look at the plot.)

The operation also adds a condition with each cell's most likely phase 
(G1, S, G2M, SubG1 or AboveG2M), for coloring or subsetting other plots.
It is approximate where the phases overlap -- use the statistic for the 
percentages.

.. note::

    How many cells are counted as S phase rather than G1 or G2/M depends 
    on the model, especially for cells that have just entered or are about
    to leave S phase.  Compare samples analyzed the same way, rather than 
    relying on the absolute numbers.

.. object:: Name

    The operation name; names the new statistic and condition.

.. object:: Channel

    The DNA stain's channel.

.. object:: Group Estimates By

    Fit separately for each value of these conditions (eg, each sample.)

.. object:: Subset

    Only use these events (the single cells.)
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
from ..workflow.operations import CellCycleWorkflowOp, CellCycleWorkflowView


class CellCycleHandler(OpHandler):
    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False, placeholder = "None")),
             Item('channel',
                  editor = EnumEditor(name = 'context_handler.previous_channels'),
                  label = "Channel"),
             Item('by',
                  editor = CheckListEditor(cols = 2,
                                           name = 'context_handler.previous_conditions_names'),
                  label = 'Group\nEstimates\nBy',
                  style = 'custom'),
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
                                        label = "Estimate!"),
                  show_label = False),
             shared_op_traits_view)


class CellCycleViewHandler(ViewHandler):
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
class CellCyclePlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.cell_cycle'
    operation_id = 'cytoflow.operations.cell_cycle'
    view_id = 'cytoflow.view.cellcycleview'

    name = "Cell Cycle"
    short_name = "Cell\nCycle"
    menu_group = "Clustering"
    
    def get_operation(self):
        return CellCycleWorkflowOp()
    
    def get_handler(self, model, context):
        if isinstance(model, CellCycleWorkflowOp):
            return CellCycleHandler(model = model, context = context)
        elif isinstance(model, CellCycleWorkflowView):
            return CellCycleViewHandler(model = model, context = context)
    
    def get_icon(self):
        return ImageResource('cell_cycle')
    
    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
