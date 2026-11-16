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
cytoflowgui.workflow.operations.autospill
-----------------------------------------

"""

import warnings

from traits.api import provides, Str, observe, List, Property, Float, Int

from cytoflow.operations.autospill import AutoSpillOp
import cytoflow.utility as util

from .. import Changed
from ..serialization import camel_registry, cytoflow_class_repr, dedent
from ..subset import ISubset

from .operation_base import IWorkflowOperation
from .bleedthrough_linear import Channel
from .compensation import CompensationWorkflowOp, CompensationWorkflowView
from cytoflow.operations.compensation import CompensationOp

AutoSpillOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class AutoSpillWorkflowOp(CompensationWorkflowOp, AutoSpillOp):
    """
    The GUI's AutoSpill operation.  It is a `CompensationWorkflowOp` (so 
    the estimated matrix shows up in the same editable table, and is 
    applied the same way) whose `estimate` runs AutoSpill on the controls.
    """
    
    # the controls: a list of (channel, file)
    channels_list = List(Channel, estimate = True)
    
    cofactor = Float(150.0, estimate = True)
    trim = Float(0.01, estimate = True)
    max_events = Int(100000, estimate = True)
    max_iterations = Int(100, estimate = True)
    
    @observe('channels_list:items, channels_list:items.+type', post_init = True)
    def _on_channels_items_changed(self, _):
        self.changed = 'channels_list'
        
    # "subset" is generated from subset_list
    subset = Property(Str, observe = "subset_list.items.str")
    subset_list = List(ISubset, estimate = True)
        
    @observe('subset_list:items.str')
    def _on_subset_changed(self, _):
        self.changed = 'subset_list'
        
    # MAGIC - returns the value of the "subset" Property, above
    def _get_subset(self):
        return " and ".join([subset.str for subset in self.subset_list if subset.str])
    
    def estimate(self, experiment):
        channels = [c.channel for c in self.channels_list]
        for c in channels:
            if not c:
                raise util.CytoflowOpError('controls', "Choose a channel for every control")
            if channels.count(c) > 1:
                raise util.CytoflowOpError('controls',
                                           "Channel {0} is included more than once"
                                           .format(c))
                
        self.controls = {c.channel : c.file for c in self.channels_list}
        
        if not self.subset:
            warnings.warn("Are you sure you don't want to specify a subset "
                          "used to estimate the model?",
                          util.CytoflowOpWarning)
            
        return AutoSpillOp.estimate(self, experiment, subset = self.subset)
    
    def apply(self, experiment):
        if not self.spillover:
            raise util.CytoflowOpError('spillover',
                                       "No spillover matrix. Add the "
                                       "single-color controls, then click "
                                       "\"Estimate!\"")
            
        return CompensationOp.apply(self, experiment)
    
    def should_clear_estimate(self, changed, payload):
        if changed == Changed.ESTIMATE:
            return True
        
        return False
    
    def get_notebook_code(self, idx):
        op = AutoSpillOp()
        op.controls = {c.channel : c.file for c in self.channels_list}
        op.cofactor = self.cofactor
        op.trim = self.trim
        op.max_events = self.max_events
        op.max_iterations = self.max_iterations
        
        # the matrix may have been edited after it was estimated, so set it
        # rather than estimating it again
        return dedent("""
        op_{idx} = {repr}
        
        # estimated with op_{idx}.estimate(ex_{prev_idx}{subset}) 
        op_{idx}.spillover = {spillover}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1,
                subset = ", subset = " + repr(self.subset) if self.subset else "",
                spillover = repr({k : dict(v) for k, v in self.spillover.items()})))
        

@camel_registry.dumper(AutoSpillWorkflowOp, 'autospill', version = 1)
def _dump(op):
    return dict(channels_list = op.channels_list,
                subset_list = op.subset_list,
                cofactor = op.cofactor,
                trim = op.trim,
                max_events = op.max_events,
                max_iterations = op.max_iterations,
                spillover = op.spillover)

@camel_registry.loader('autospill', version = 1)
def _load(data, version):
    return AutoSpillWorkflowOp(**data)
