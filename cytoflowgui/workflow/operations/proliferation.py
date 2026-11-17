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
cytoflowgui.workflow.operations.proliferation
---------------------------------------------

"""

from traits.api import (HasTraits, provides, Str, Property, observe, Instance,
                        List, Dict, Any, Bool, Int)

from cytoflow.operations.proliferation import ProliferationOp, ProliferationView
import cytoflow.utility as util

from .. import Changed
from ..views import IWorkflowView, WorkflowView
from ..serialization import camel_registry, cytoflow_class_repr, traits_str, dedent
from ..subset import ISubset

from .operation_base import IWorkflowOperation, WorkflowOperation

ProliferationOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)    
class ProliferationWorkflowOp(WorkflowOperation, ProliferationOp):
    # add 'estimate', 'apply' metadata
    name = Str(apply = True)
    channel = Str(estimate = True)
    by = List(Str, estimate = True)
    num_generations = util.PositiveCInt(8, allow_zero = False, estimate = True)
    undivided_peak = util.CFloatOrNone(None, estimate = True)
    fit_spacing = Bool(False, estimate = True)
    
    # the unstimulated control: "control_subset" is generated from this list
    control_subset = Property(Str, observe = "control_subset_list.items.str")
    control_subset_list = List(ISubset, estimate = True)

    # add the 'estimate_result' metadata
    _fits = Dict(Any, Any, estimate_result = True, transient = True)
    
    # "subset" is generated from subset_list
    subset = Property(Str, observe = "subset_list.items.str")
    subset_list = List(ISubset, estimate = True)
    
    @observe('subset_list:items.str')
    def _on_subset_changed(self, _):
        self.changed = 'subset_list'
        
    # MAGIC - returns the value of the "subset" Property, above
    def _get_subset(self):
        return " and ".join([subset.str for subset in self.subset_list if subset.str])
    
    @observe('control_subset_list:items.str')
    def _on_control_subset_changed(self, _):
        self.changed = 'control_subset_list'
        
    # MAGIC - returns the value of the "control_subset" Property, above
    def _get_control_subset(self):
        return " and ".join([subset.str for subset in self.control_subset_list if subset.str])

    def default_view(self, **kwargs):
        return ProliferationWorkflowView(op = self, **kwargs)
    
    def estimate(self, experiment):
        super().estimate(experiment, subset = self.subset)
        
    def apply(self, experiment):
        if not self._fits:
            raise util.CytoflowOpError(None, 'Click "Estimate"!')
        
        return super().apply(experiment)
        
    def clear_estimate(self):
        self._fits = {}
    
    def get_notebook_code(self, idx):
        op = ProliferationOp()
        op.copy_traits(self, op.copyable_trait_names())
        op.control_subset = self.control_subset

        return dedent("""
        op_{idx} = {repr}
        
        op_{idx}.estimate(ex_{prev_idx}{subset})
        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1,
                subset = ", subset = " + repr(self.subset) if self.subset else ""))
        
        
@provides(IWorkflowView)
class ProliferationWorkflowView(WorkflowView, ProliferationView):
    op = Instance(IWorkflowOperation, fixed = True)
    plot_params = Instance(HasTraits, ())
    
    def should_plot(self, changed, payload):
        return changed in (Changed.VIEW, Changed.ESTIMATE_RESULT, Changed.PREV_RESULT)
    
    def get_notebook_code(self, idx):
        view = ProliferationView()
        view.copy_traits(self, view.copyable_trait_names())

        return dedent("""
        op_{idx}.default_view({traits}).plot(ex_{idx})
        """
        .format(traits = traits_str(view),
                idx = idx))
        

@camel_registry.dumper(ProliferationWorkflowOp, 'proliferation', version = 1)
def _dump(op):
    return dict(name = op.name,
                channel = op.channel,
                by = op.by,
                num_generations = op.num_generations,
                undivided_peak = op.undivided_peak,
                fit_spacing = op.fit_spacing,
                control_subset_list = op.control_subset_list,
                subset_list = op.subset_list)
    
@camel_registry.loader('proliferation', version = 1)
def _load(data, version):
    return ProliferationWorkflowOp(**data)

@camel_registry.dumper(ProliferationWorkflowView, 'proliferation-view', version = 1)
def _dump_view(view):
    return dict(op = view.op,
                subset_list = view.subset_list)

@camel_registry.loader('proliferation-view', version = 1)
def _load_view(data, version):
    return ProliferationWorkflowView(**data)
