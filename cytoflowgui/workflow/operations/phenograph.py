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
cytoflowgui.workflow.operations.phenograph
------------------------------------------

"""

from traits.api import (provides, Instance, Str, Property, observe, Constant,
                        List, Dict, Any, Enum, DelegatesTo)

import numpy as np

from cytoflow.operations.phenograph import PhenographOp, Phenograph2DView, Phenograph1DView
import cytoflow.utility as util

from ..views import IWorkflowView, WorkflowByView, ScatterplotPlotParams
from ..serialization import camel_registry, traits_str, cytoflow_class_repr, dedent
from ..subset import ISubset

from .operation_base import IWorkflowOperation, WorkflowOperation
from .tsne import Channel

PhenographOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class PhenographWorkflowOp(WorkflowOperation, PhenographOp):
    # a list of (channel, scale) instead of separate channels / scale
    channels_list = List(Channel, estimate = True)
    channels = Property(List(Str), 
                        observe = '[channels_list.items,channels_list.items.channel,channels_list.items.scale]')
    scale = Property(Dict(Str, util.ScaleEnum),
                     observe = '[channels_list.items,channels_list.items.channel,channels_list.items.scale]')
    
    # add "apply", "estimate" metadata
    name = Str(apply = True)
    num_neighbors = util.PositiveCInt(30, allow_zero = False, estimate = True)
    resolution = util.PositiveCFloat(1.0, allow_zero = False, estimate = True)
    method = Enum("leiden", "louvain", estimate = True)
    sample = util.UnitFloat(1.0, estimate = True)
    by = List(Str, estimate = True)
    
    # add the 'estimate_result' metadata.  (this is the last thing that
    # PhenographOp.estimate() sets.)
    _centers = Dict(Any, Instance(np.ndarray), transient = True, estimate_result = True)
    
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
        
    @observe('subset_list:items.str')
    def _on_subset_changed(self, _):
        self.changed = 'subset_list'
        
    # MAGIC - returns the value of the "subset" Property, above
    def _get_subset(self):
        return " and ".join([subset.str for subset in self.subset_list if subset.str])
            
    def default_view(self, **kwargs):
        return PhenographWorkflowView(op = self, **kwargs)
    
    def estimate(self, experiment):
        channels = [c.channel for c in self.channels_list]
        for c in channels:
            if channels.count(c) > 1:
                raise util.CytoflowOpError('channels',
                                           "Channel {0} is included more than once"
                                           .format(c))
        
        super().estimate(experiment, subset = self.subset)
        
    def apply(self, experiment):
        if not self._centers:
            raise util.CytoflowOpError(None, 'Click "Estimate"!')
        return PhenographOp.apply(self, experiment)
    
    def clear_estimate(self):
        self._nn = {}
        self._labels = {}
        self._scale = {}
        self._centers = {}      
    
    def get_notebook_code(self, idx):
        op = PhenographOp()
        op.copy_traits(self, op.copyable_trait_names())     
        op.channels = [c.channel for c in self.channels_list]
        op.scale = {c.channel : c.scale for c in self.channels_list}

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
class PhenographWorkflowView(WorkflowByView):
    """
    The diagnostic plot: the first two channels of the clustering (or a 
    histogram, if there's only one), colored by cluster.
    """
    op = Instance(IWorkflowOperation, fixed = True)
    subset = DelegatesTo('op')
    by = DelegatesTo('op')
    
    plot_params = Instance(ScatterplotPlotParams, ())
    
    id = "cytoflowgui.op_plugins.phenograph"
    friendly_id = "Graph Clustering" 
    
    name = Constant("Graph Clustering")
    
    def _view(self):
        channels = [c for c in self.op.channels if c]
        scale = self.op.scale
        if len(channels) == 0:
            raise util.CytoflowViewError(None, "Add at least one channel")
        elif len(channels) == 1:
            return Phenograph1DView(op = self.op,
                                    channel = channels[0],
                                    scale = scale[channels[0]])
        else:
            return Phenograph2DView(op = self.op,
                                    xchannel = channels[0],
                                    ychannel = channels[1],
                                    xscale = scale[channels[0]],
                                    yscale = scale[channels[1]])
    
    def plot(self, experiment, **kwargs):
        view = self._view()
        if isinstance(view, Phenograph1DView):
            # the scatterplot's parameters don't apply to a histogram
            for k in ('alpha', 's', 'marker'):
                kwargs.pop(k, None)
        view.plot(experiment, **kwargs)
                     
    def enum_plots(self, experiment):
        return self._view().enum_plots(experiment)
            
    def get_notebook_code(self, idx):
        plot_params_str = traits_str(self.plot_params)
        channels = [c for c in self.op.channels if c][:2]
        
        return dedent("""
        op_{idx}.default_view(channels = {channels}).plot(ex_{idx}{plot}{plot_params})
        """
        .format(channels = repr(channels),
                idx = idx,
                plot = ", plot_name = " + repr(self.current_plot) if self.current_plot else "",
                plot_params = ", " + plot_params_str if plot_params_str else ""))
        
        
@camel_registry.dumper(PhenographWorkflowOp, 'phenograph', version = 1)
def _dump(op):
    return dict(name = op.name,
                channels_list = op.channels_list,
                num_neighbors = op.num_neighbors,
                resolution = op.resolution,
                method = op.method,
                sample = op.sample,
                by = op.by,
                subset_list = op.subset_list)
    
@camel_registry.loader('phenograph', version = 1)
def _load(data, version):
    return PhenographWorkflowOp(**data)

@camel_registry.dumper(PhenographWorkflowView, 'phenograph-view', version = 1)
def _dump_view(view):
    return dict(op = view.op,
                plot_params = view.plot_params,
                current_plot = view.current_plot)

@camel_registry.loader('phenograph-view', version = 1)
def _load_view(data, ver):
    return PhenographWorkflowView(**data)
