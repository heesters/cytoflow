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
Population Tree
---------------

Shows your gates as a **tree of nested populations**, with the number of 
events in each, as a percentage of its parent and of all events -- like 
the workspace tree in FlowJo or the Gate Navigator in FCS Express.

In *Cytoflow*, a gate applies to all the events.  A **population** is a 
gate applied *inside its parent population*: for example "CD4+" inside 
"T cells" inside "Singlets".  This operation doesn't change your data; it 
adds a statistic with the numbers (see it with a **Table** view, where you 
can also export it), and its plot draws the tree.

The quickest way to set it up is to click **Build from workflow**.  Every 
gate before this operation becomes a population.  Its parent is the gate 
that it was drawn inside: either the parent recorded when the gates were 
imported from FlowJo, or the gate that you had selected under *Subset* when 
you drew it.  Check the result -- you can edit any row, and add or remove 
populations.

.. object:: Name

    The operation name; becomes the name of the new statistic.

.. object:: Populations

    For each population: its **name**; its **gate** -- the name of a gate, 
    or an expression such as ``Quad == 'Quad_2'`` or ``CD4_pos and not 
    CD8_pos``; and its **parent** population (leave empty for a population 
    at the top of the tree.)

.. object:: Build from workflow

    Fill in the populations from the gates in the workflow.  This replaces
    the current list.

.. object:: Group By

    Compute the numbers separately for each value of these conditions -- 
    for example, for each sample.  (The plot pools all the events; use its
    *Subset* to look at one sample.)

The statistic has the columns **Count**, **Percent_of_parent**, 
**Percent_of_total** and **Depth** (how deep the population is in the 
tree.)
'''

from traits.api import provides, Event, List
from traitsui.api import (View, Item, VGroup, HGroup, TextEditor, ButtonEditor,
                          CheckListEditor, Controller)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..editors import (SubsetListEditor, ColorTextEditor, InstanceHandlerEditor,
                       VerticalListEditor)
from ..workflow.operations import (PopulationTreeWorkflowOp, PopulationTreeWorkflowView,
                                   PopulationTreePopulation, populations_from_workflow,
                                   BooleanGateWorkflowOp)
from ..workflow.workflow_item import WorkflowItem
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class PopulationHandler(Controller):
    population_view = View(VGroup(HGroup(Item('name',
                                              editor = TextEditor(auto_set = False,
                                                                  placeholder = "Name"),
                                              show_label = False),
                                         Item('parent',
                                              editor = TextEditor(auto_set = False,
                                                                  placeholder = "Parent (or empty)"),
                                              show_label = False)),
                                  Item('gate',
                                       editor = TextEditor(auto_set = False,
                                                           placeholder = "Gate, eg  Singlets"),
                                       show_label = False),
                                  show_border = True))


class PopulationTreeHandler(OpHandler):

    add_population = Event
    remove_population = Event
    build_from_workflow = Event
    create_gates = Event

    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False, placeholder = "None")),
             VGroup(Item('handler.build_from_workflow',
                         editor = ButtonEditor(value = True,
                                               label = "Build from workflow")),
                    Item('populations_list',
                         editor = VerticalListEditor(editor = InstanceHandlerEditor(view = 'population_view',
                                                                                    handler_factory = PopulationHandler),
                                                     style = 'custom',
                                                     mutable = False)),
                    Item('handler.add_population',
                         editor = ButtonEditor(value = True,
                                               label = "Add a population")),
                    Item('handler.remove_population',
                         editor = ButtonEditor(value = True,
                                               label = "Remove a population")),
                    Item('handler.create_gates',
                         editor = ButtonEditor(value = True,
                                               label = "Create Boolean Gates from populations"),
                         tooltip = "Insert a BooleanGate op for every population in the tree"),
                    label = "Populations",
                    show_labels = False),
             Item('by',
                  editor = CheckListEditor(cols = 2,
                                           name = 'context_handler.previous_conditions_names'),
                  label = 'Group\nBy',
                  style = 'custom'),
             shared_op_traits_view)

    # MAGIC: called when add_population is set
    def _add_population_fired(self):
        self.model.populations_list.append(PopulationTreePopulation())
        
    def _remove_population_fired(self):
        if self.model.populations_list:
            self.model.populations_list.pop()
            
    def _build_from_workflow_fired(self):
        # walk back through the workflow, collecting each operation and
        # the subsets of the view its gate was drawn on
        items = []
        wi = self.context.previous_wi if self.context else None
        while wi is not None:
            view = wi.default_view
            subsets = list(view.subset_list) if view is not None and hasattr(view, 'subset_list') else []
            items.insert(0, (wi.operation, subsets))
            wi = wi.previous_wi

        self.model.populations_list = [PopulationTreePopulation(name = n, gate = g, parent = p)
                                       for n, g, p in populations_from_workflow(items)]

    def _create_gates_fired(self):
        """Insert a BooleanGateWorkflowOp for every population in the tree."""
        if not self.model.populations_list:
            return
        workflow = self.context
        # find the WorkflowItem that owns self.model so we can insert after it
        my_wi = next((wi for wi in workflow.workflow if wi.operation is self.model), None)
        if my_wi is None:
            return
        insert_idx = workflow.workflow.index(my_wi) + 1
        for pop in self.model.populations_list:
            if not pop.gate:
                continue
            gate_name = (pop.name + '_gate') if pop.name else None
            if not gate_name:
                continue
            op = BooleanGateWorkflowOp(name=gate_name, expression=pop.gate)
            new_wi = WorkflowItem(operation=op, workflow=workflow)
            workflow.workflow.insert(insert_idx, new_wi)
            insert_idx += 1
        # select the first inserted item (if any)
        if insert_idx > workflow.workflow.index(my_wi) + 1:
            workflow.selected = workflow.workflow[workflow.workflow.index(my_wi) + 1]
        

class PopulationTreeViewHandler(ViewHandler):
    view_traits_view = \
        View(VGroup(Item('subset_list',
                         show_label = False,
                         editor = SubsetListEditor(conditions = "context_handler.conditions",
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
                                           background_color = "#ff9191")))

    view_params_view = View()


@provides(IOperationPlugin)
class PopulationTreePlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.population_tree'
    operation_id = 'cytoflow.operations.population_tree'
    view_id = 'cytoflow.view.populationtreeview'

    name = "Population Tree"
    short_name = "Population\nTree"
    menu_group = "Gates"
    
    def get_operation(self):
        return PopulationTreeWorkflowOp()
    
    def get_handler(self, model, context):
        if isinstance(model, PopulationTreeWorkflowOp):
            return PopulationTreeHandler(model = model, context = context)
        elif isinstance(model, PopulationTreeWorkflowView):
            return PopulationTreeViewHandler(model = model, context = context)
    
    def get_icon(self):
        return ImageResource('population_tree')
    
    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
