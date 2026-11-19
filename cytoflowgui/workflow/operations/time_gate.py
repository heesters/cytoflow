#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.workflow.operations.time_gate
-----------------------------------------
"""

from traits.api import provides, Str, Tuple, Property, observe

from cytoflow.operations.time_gate import TimeGateOp, TimeGateSelection
import cytoflow.utility as util

from ..views import IWorkflowView, WorkflowView, HistogramPlotParams
from ..serialization import camel_registry, traits_str, cytoflow_class_repr, dedent

from .operation_base import IWorkflowOperation, WorkflowOperation
from .range import RangeSelectionView

TimeGateOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowView)
class TimeGateSelectionView(RangeSelectionView):
    """Workflow view for the time gate — reuses the range selection view."""
    pass


@provides(IWorkflowOperation)
class TimeGateWorkflowOp(WorkflowOperation, TimeGateOp):
    name = Str(apply=True)
    # channel is set automatically; allow override
    channel = Str(apply=True)

    _range = Tuple(util.FloatOrNone(None), util.FloatOrNone(None), apply=True)
    low  = Property(util.FloatOrNone(None), observe='_range')
    high = Property(util.FloatOrNone(None), observe='_range')

    def _get_low(self):
        return self._range[0]

    def _set_low(self, val):
        self._range = (val, self._range[1])

    def _get_high(self):
        return self._range[1]

    def _set_high(self, val):
        self._range = (self._range[0], val)

    def default_view(self, **kwargs):
        return TimeGateSelectionView(op=self, **kwargs)

    def get_notebook_code(self, idx):
        op = TimeGateOp()
        op.copy_traits(self, op.copyable_trait_names())
        return dedent("""
        op_{idx} = {repr}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """.format(repr=repr(op), idx=idx, prev_idx=idx - 1))


### Serialization

@camel_registry.dumper(TimeGateWorkflowOp, 'time-gate', version=1)
def _dump(op):
    return dict(name=op.name,
                channel=op.channel,
                low=op.low,
                high=op.high)


@camel_registry.loader('time-gate', version=1)
def _load(data, version):
    return TimeGateWorkflowOp(**data)


@camel_registry.dumper(TimeGateSelectionView, 'time-gate-view', version=1)
def _dump_view(view):
    return dict(op=view.op,
                scale=view.scale,
                huefacet=view.huefacet,
                subset_list=view.subset_list,
                plot_params=view.plot_params,
                current_plot=view.current_plot)


@camel_registry.loader('time-gate-view', version=1)
def _load_view(data, version):
    return TimeGateSelectionView(**data)
