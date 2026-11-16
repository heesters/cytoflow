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
cytoflowgui.workflow.operations.cross_entropy
---------------------------------------------

"""

from traits.api import (HasTraits, provides, Str, Property, observe, Instance,
                        List, Dict, Any, Enum, Bool, Float, Int)

from cytoflow.operations.cross_entropy import CrossEntropyTestOp, CrossEntropyTestView
import cytoflow.utility as util

from .. import Changed
from ..views import IWorkflowView, WorkflowView
from ..serialization import camel_registry, cytoflow_class_repr, traits_str, dedent
from ..subset import ISubset

from .operation_base import IWorkflowOperation, WorkflowOperation
from .tsne import Channel

CrossEntropyTestOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)    
class CrossEntropyTestWorkflowOp(WorkflowOperation, CrossEntropyTestOp):
    # a list of (channel, scale) instead of separate channels / scale
    channels_list = List(Channel, estimate = True)
    channels = Property(List(Str), 
                        observe = '[channels_list.items,channels_list.items.channel,channels_list.items.scale]')
    scale = Property(Dict(Str, util.ScaleEnum),
                     observe = '[channels_list.items,channels_list.items.channel,channels_list.items.scale]')
    
    # the embedding's two channels
    xembedding = Str(estimate = True)
    yembedding = Str(estimate = True)
    embedding = Property(List(Str), observe = '[xembedding,yembedding]')
    
    # add 'estimate', 'apply' metadata
    name = Str(apply = True)
    by = Str(estimate = True)
    perplexity = util.PositiveFloat(30, estimate = True)
    max_events = Int(20000, estimate = True)
    correction = Enum("holm", "bonferroni", "none", estimate = True)
    alpha = Float(0.05, estimate = True)
    
    # add the 'estimate_result' metadata
    _cross_entropy = Dict(Any, Any, estimate_result = True, transient = True)
    
    # "subset" is generated from subset_list
    subset = Property(Str, observe = "subset_list.items.str")
    subset_list = List(ISubset, estimate = True)
    
    @observe('[channels_list:items,channels_list:items.channel,channels_list:items.scale]')
    def _channels_updated(self, event):
        self.changed = 'channels_list'
        
    def _get_channels(self):
        return [c.channel for c in self.channels_list]
    
    def _get_scale(self):
        return {c.channel : c.scale for c in self.channels_list}
    
    def _get_embedding(self):
        return [c for c in (self.xembedding, self.yembedding) if c]
    
    @observe('subset_list:items.str')
    def _on_subset_changed(self, _):
        self.changed = 'subset_list'
        
    # MAGIC - returns the value of the "subset" Property, above
    def _get_subset(self):
        return " and ".join([subset.str for subset in self.subset_list if subset.str])
    
    def default_view(self, **kwargs):
        return CrossEntropyTestWorkflowView(op = self, **kwargs)
    
    def estimate(self, experiment):
        channels = [c.channel for c in self.channels_list]
        for c in channels:
            if channels.count(c) > 1:
                raise util.CytoflowOpError('channels',
                                           "Channel {0} is included more than once"
                                           .format(c))
        if self.xembedding and self.xembedding == self.yembedding:
            raise util.CytoflowOpError('embedding',
                                       "The two embedding channels must be different")
            
        super().estimate(experiment, subset = self.subset)
        
    def apply(self, experiment):
        if not self._cross_entropy:
            raise util.CytoflowOpError(None, 'Click "Estimate"!')
        
        return super().apply(experiment)
        
    def clear_estimate(self):
        self._pairs = None
        self._groups = None
        self._cross_entropy = {}
    
    def get_notebook_code(self, idx):
        op = CrossEntropyTestOp()
        op.copy_traits(self, op.copyable_trait_names())
        
        op.channels = [c.channel for c in self.channels_list]
        op.scale = {c.channel : c.scale for c in self.channels_list}
        op.embedding = self.embedding

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
class CrossEntropyTestWorkflowView(WorkflowView, CrossEntropyTestView):
    op = Instance(IWorkflowOperation, fixed = True)
    dendrogram = Bool(False)
    plot_params = Instance(HasTraits, ())
    
    def should_plot(self, changed, payload):
        return changed in (Changed.VIEW, Changed.ESTIMATE_RESULT)
    
    def get_notebook_code(self, idx):
        view = CrossEntropyTestView()
        view.copy_traits(self, view.copyable_trait_names())

        return dedent("""
        op_{idx}.default_view({traits}).plot(ex_{idx})
        """
        .format(traits = traits_str(view),
                idx = idx))
        

@camel_registry.dumper(CrossEntropyTestWorkflowOp, 'cross-entropy-test', version = 1)
def _dump(op):
    return dict(name = op.name,
                channels_list = op.channels_list,
                xembedding = op.xembedding,
                yembedding = op.yembedding,
                by = op.by,
                perplexity = op.perplexity,
                max_events = op.max_events,
                correction = op.correction,
                alpha = op.alpha,
                subset_list = op.subset_list)
    
@camel_registry.loader('cross-entropy-test', version = 1)
def _load(data, version):
    return CrossEntropyTestWorkflowOp(**data)

@camel_registry.dumper(CrossEntropyTestWorkflowView, 'cross-entropy-test-view', version = 1)
def _dump_view(view):
    return dict(op = view.op,
                dendrogram = view.dendrogram)

@camel_registry.loader('cross-entropy-test-view', version = 1)
def _load_view(data, version):
    return CrossEntropyTestWorkflowView(**data)
