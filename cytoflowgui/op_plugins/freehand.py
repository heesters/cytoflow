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
Freehand Gate
-------------

Draw a freehand (lasso) gate: press the mouse button (or the trackpad) and
**drag** around the population you want.  The gate is closed when you let
go -- you don't need to finish where you started.

- To change the gate, just draw it again.
- A stray click (or a very short drag) doesn't erase the gate.

The gate is stored as a polygon (simplified to at most 100 vertices), so it
works exactly like a polygon gate.

.. object:: Name

    The operation name.  Used to name the new metadata field that's created by 
    this module.
    
.. object:: X Channel

    The name of the channel on the gate's X axis.

.. object:: Y Channel

    The name of the channel on the gate's Y axis.
    
.. object:: Adjust Per

    Makes this a "magnetic" gate.  Tick one or more conditions (eg, the one 
    that identifies each tube or well), and the gate is shifted separately 
    for each of their values, so it follows the population when it sits in
    a slightly different place in each sample.  The plot still shows the 
    gate as you drew it; the shift applied to each sample is saved as a
    statistic with the gate's name, which you can look at with a table view.
    
.. object:: Max Adjustment

    The furthest the gate may move, as a fraction of its own width and 
    height (default 0.25).  Keep it small enough that the gate can't jump to
    a neighboring population.
    
.. object:: X Scale

    The scale of the X axis for the interactive plot.
    
.. object:: Y Scale

    The scale of the Y axis for the interactive plot
    
.. object:: Hue facet

    Show different experimental conditions in different colors.
    
.. object:: Subset

    Show only a subset of the data.
   
.. plot::
   :include-source: False

    import numpy as np
    import cytoflow as flow
    import_op = flow.ImportOp()
    import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
                                 conditions = {'Dox' : 10.0}),
                       flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
                                 conditions = {'Dox' : 1.0})]
    import_op.conditions = {'Dox' : 'float'}
    ex = import_op.apply()
    
    # a wobbly loop, like one drawn by hand
    t = np.linspace(0, 2 * np.pi, 60, endpoint = False)
    r = 1 + 0.15 * np.sin(3 * t) + 0.05 * np.cos(7 * t)
    x = 10 ** (2.0 + 1.0 * r * np.cos(t))
    y = 10 ** (3.8 + 0.5 * r * np.sin(t))

    f = flow.FreehandOp(name = "Freehand",
                        xchannel = "V2-A",
                        ychannel = "Y2-A",
                        xscale = "log",
                        yscale = "log",
                        vertices = list(zip(x, y)))

    f.default_view(huefacet = "Dox").plot(ex)

'''

from traits.api import provides, List
from traitsui.api import View, Item, EnumEditor, VGroup, TextEditor, Controller, CheckListEditor
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..view_plugins.view_plugin_base import Data2DPlotParamsView
from ..editors import SubsetListEditor, ColorTextEditor, ExtendableEnumEditor, InstanceHandlerEditor
from ..workflow.operations import FreehandWorkflowOp, FreehandSelectionView
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin

class FreehandViewParamsHandler(Controller):
    view_params_view = \
        View(Item('density'),
             # density-specific
             Item('gridsize',
                  editor = TextEditor(auto_set = False),
                  visible_when = 'density',
                  label = "Grid size"),
             Item('smoothed',
                  visible_when = 'density',
                  label = "Smooth"),
             Item('smoothed_sigma',
                  editor = TextEditor(auto_set = False),
                  label = "Smooth\nsigma",
                  visible_when = "density and smoothed"),
             
             # scatterplot-specific
             Item('alpha',
                  editor = TextEditor(auto_set = False),
                  visible_when = 'not density'),
             Item('s',
                  editor = TextEditor(auto_set = False),
                  visible_when = 'not density',
                  label = "Size"),
             Item('marker',
                  visible_when = 'not density'),
             Data2DPlotParamsView.content)

class FreehandHandler(OpHandler):
    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False)),
             Item('xchannel',
                  editor=EnumEditor(name='context_handler.previous_channels'),
                  label = "X Channel"),
             Item('ychannel',
                  editor=EnumEditor(name='context_handler.previous_channels'),
                  label = "Y Channel"),
             Item('adjust_by',
                  editor = CheckListEditor(cols = 2,
                                           name = 'context_handler.previous_conditions_names'),
                  label = 'Adjust\nPer',
                  style = 'custom',
                  tooltip = "Magnetic gate: shift the gate separately for "
                            "each value of these conditions, to follow the "
                            "population"),
             Item('max_adjust',
                  editor = TextEditor(auto_set = False),
                  label = "Max\nAdjustment",
                  visible_when = 'adjust_by'),
             shared_op_traits_view) 
        
        
class FreehandViewHandler(ViewHandler):
    view_traits_view = \
        View(VGroup(
             VGroup(Item('xchannel', 
                         label = "X Channel", 
                         style = 'readonly'),
                    Item('xscale',
                         label = "X Scale"),
                    Item('ychannel',
                         label = "Y Channel",
                         style = 'readonly'),
                    Item('yscale',
                         label = "Y Scale"),
                    Item('huefacet',
                         editor=ExtendableEnumEditor(name='context_handler.previous_conditions_names',
                                                     extra_items = {"None" : ""}),
                         visible_when = 'not plot_params.density',
                         label="Color\nFacet"),
                    Item('huescale',
                         visible_when = 'plot_params.density',
                         label = "Color\nScale"),
                    label = "Plot settings",
                    show_border = False),
             VGroup(Item('subset_list',
                         show_label = False,
                         editor = SubsetListEditor(conditions = "context_handler.previous_conditions",
                                                   editor = InstanceHandlerEditor(view = 'subset_view',
                                                                                  handler_factory = subset_handler_factory),
                                                   mutable = False)),
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
                                           background_color = "#ff9191"))))
        
    view_params_view = \
        View(Item('plot_params',
                  editor = InstanceHandlerEditor(view = 'view_params_view',
                                                 handler_factory = FreehandViewParamsHandler),
                  style = 'custom',
                  show_label = False))
        



@provides(IOperationPlugin)
class FreehandPlugin(Plugin, PluginHelpMixin):
    
    id = 'cytoflowgui.op_plugins.freehand'
    operation_id = 'cytoflow.operations.freehand'
    view_id = 'cytoflowgui.workflow.operations.freehandview'

    name = "Freehand Gate"
    short_name = "Freehand"
    menu_group = "Gates"
    
    def get_operation(self):
        return FreehandWorkflowOp()
    
    def get_handler(self, model, context):
        if isinstance(model, FreehandWorkflowOp):
            return FreehandHandler(model = model, context = context)
        elif isinstance(model, FreehandSelectionView):
            return FreehandViewHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('freehand')
    
    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
