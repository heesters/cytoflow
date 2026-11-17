#!/usr/bin/env python3.8

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


"""
cytoflowgui.workflow.operations.population_tree
-----------------------------------------------

"""

from traits.api import (HasTraits, provides, Str, Property, observe, Instance,
                        List, Tuple)

from cytoflow.operations.population_tree import (PopulationTreeOp, PopulationTreeView,
                                                 infer_populations, BOOLEAN_GATE_IDS)
import cytoflow.utility as util

from .. import Changed
from ..views import IWorkflowView, WorkflowView
from ..serialization import camel_registry, cytoflow_class_repr, traits_repr, traits_str, dedent

from .operation_base import IWorkflowOperation, WorkflowOperation

PopulationTreeOp.__repr__ = cytoflow_class_repr


class Population(HasTraits):
    """One population, for the GUI's list"""
    name = Str
    gate = Str
    parent = Str
    
    def __repr__(self):
        return traits_repr(self)
    
    
def populations_from_workflow(items):
    """
    Work out the population tree from the operations in a workflow.
    
    Parameters
    ----------
    items : list of (operation, subsets)
        The workflow's operations in order, each with the subsets of the 
        view that its gate was drawn on (a list of `ISubset`, possibly 
        empty.)
        
    Returns
    -------
    list of (name, gate, parent)
    
    Notes
    -----
    Every gate that makes a True/False condition becomes a population.  Its
    parent is, in order of preference: the parent recorded by a gate 
    imported from FlowJo or Gating-ML (``X_gated = X and Parent``); or the 
    gate that its plot was subset to, if the plot showed exactly one gate's
    "True" events.
    """
    ops = [op for op, _ in items]
    gate_names = [op.name for op in ops if getattr(op, 'id', '') in BOOLEAN_GATE_IDS]
    
    parents = {}
    for op, subsets in items:
        if getattr(op, 'id', '') not in BOOLEAN_GATE_IDS:
            continue
        selected = []
        for s in subsets or []:
            if not getattr(s, 'str', ''):
                continue
            # only "this gate == True" counts as a parent
            if getattr(s, 'selected_t', False) and not getattr(s, 'selected_f', False) \
               and s.name in gate_names and s.name != op.name:
                selected.append(s.name)
            else:
                selected.append(None)
        if len(selected) == 1 and selected[0] is not None:
            parents[op.name] = selected[0]
            
    return infer_populations(ops, parents = parents)


@provides(IWorkflowOperation)    
class PopulationTreeWorkflowOp(WorkflowOperation, PopulationTreeOp):
    name = Str(apply = True)
    by = List(Str, apply = True)
    
    populations_list = List(Population, apply = True)
    populations = Property(List(Tuple(Str, Str, Str)),
                           observe = '[populations_list.items,populations_list.items.name,populations_list.items.gate,populations_list.items.parent]')
    
    @observe('[populations_list:items,populations_list:items.name,populations_list:items.gate,populations_list:items.parent]')
    def _populations_updated(self, _):
        self.changed = 'populations_list'
        
    def _get_populations(self):
        return [(p.name, p.gate, p.parent) for p in self.populations_list]
    
    def default_view(self, **kwargs):
        return PopulationTreeWorkflowView(op = self, **kwargs)
    
    def clear_estimate(self):
        # no-op
        return
    
    def get_notebook_code(self, idx):
        op = PopulationTreeOp()
        op.name = self.name
        op.by = self.by
        op.populations = self.populations

        return dedent("""
        op_{idx} = {repr}
        
        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1))
        
        
@provides(IWorkflowView)
class PopulationTreeWorkflowView(WorkflowView, PopulationTreeView):
    op = Instance(IWorkflowOperation, fixed = True)
    plot_params = Instance(HasTraits, ())
    
    def should_plot(self, changed, payload):
        return changed in (Changed.VIEW, Changed.RESULT, Changed.PREV_RESULT)
    
    def get_notebook_code(self, idx):
        view = PopulationTreeView()
        view.copy_traits(self, view.copyable_trait_names())

        return dedent("""
        op_{idx}.default_view({traits}).plot(ex_{idx})
        """
        .format(traits = traits_str(view),
                idx = idx))
        

@camel_registry.dumper(PopulationTreeWorkflowOp, 'population-tree', version = 1)
def _dump(op):
    return dict(name = op.name,
                populations_list = op.populations_list,
                by = op.by)
    
@camel_registry.loader('population-tree', version = 1)
def _load(data, version):
    return PopulationTreeWorkflowOp(**data)

@camel_registry.dumper(Population, 'population-tree-population', version = 1)
def _dump_population(p):
    return dict(name = p.name,
                gate = p.gate,
                parent = p.parent)
    
@camel_registry.loader('population-tree-population', version = 1)
def _load_population(data, version):
    return Population(**data)

@camel_registry.dumper(PopulationTreeWorkflowView, 'population-tree-view', version = 1)
def _dump_view(view):
    return dict(op = view.op,
                subset_list = view.subset_list)

@camel_registry.loader('population-tree-view', version = 1)
def _load_view(data, version):
    return PopulationTreeWorkflowView(**data)
