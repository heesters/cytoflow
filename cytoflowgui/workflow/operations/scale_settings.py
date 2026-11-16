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
cytoflowgui.workflow.operations.scale_settings
----------------------------------------------

"""

from textwrap import dedent
from traits.api import HasTraits, Str, List, provides, observe

import cytoflow.utility as util
from cytoflow import ScaleSettingsOp
                       
from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr, traits_repr
from .operation_base import IWorkflowOperation, WorkflowOperation

ScaleSettingsOp.__repr__ = cytoflow_class_repr


class ChannelScaleSettings(HasTraits):
    """The settings for one channel.  Empty values aren't set."""
    
    channel = Str
    logicle_W = util.CFloatOrNone(None)
    logicle_M = util.CFloatOrNone(None)
    logicle_A = util.CFloatOrNone(None)
    arcsinh_cofactor = util.CFloatOrNone(None)
    
    def __repr__(self):
        return traits_repr(self)
    
    def to_dict(self):
        ret = {}
        for scale, param, value in (('logicle', 'W', self.logicle_W),
                                    ('logicle', 'M', self.logicle_M),
                                    ('logicle', 'A', self.logicle_A),
                                    ('arcsinh', 'cofactor', self.arcsinh_cofactor)):
            if value is not None:
                ret.setdefault(scale, {})[param] = float(value)
        return ret


@provides(IWorkflowOperation)
class ScaleSettingsWorkflowOp(WorkflowOperation, ScaleSettingsOp):
    
    settings_list = List(ChannelScaleSettings, apply = True)
    
    @observe('settings_list:items,settings_list:items.+type', post_init = True)
    def _on_settings_changed(self, _):
        self.changed = 'settings_list'
        
    def _make_settings(self):
        settings = {}
        for s in self.settings_list:
            if s.channel:
                d = s.to_dict()
                if d:
                    settings.setdefault(s.channel, {}).update(d)
        return settings
        
    def apply(self, experiment):
        self.settings = self._make_settings()
        return super().apply(experiment)
    
    def clear_estimate(self):
        # no-op
        return
    
    def get_notebook_code(self, idx):
        op = ScaleSettingsOp()
        op.settings = self._make_settings()

        return dedent("""
        op_{idx} = {repr}
                
        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1))
        
    
@camel_registry.dumper(ScaleSettingsWorkflowOp, 'scale-settings', version = 1)
def _dump(op):
    return dict(settings_list = op.settings_list)
    
@camel_registry.loader('scale-settings', version = 1)
def _load(data, version):
    return ScaleSettingsWorkflowOp(**data)

@camel_registry.dumper(ChannelScaleSettings, 'scale-settings-channel', version = 1)
def _dump_channel(s):
    return dict(channel = s.channel,
                logicle_W = s.logicle_W,
                logicle_M = s.logicle_M,
                logicle_A = s.logicle_A,
                arcsinh_cofactor = s.arcsinh_cofactor)
    
@camel_registry.loader('scale-settings-channel', version = 1)
def _load_channel(data, version):
    return ChannelScaleSettings(**data)
