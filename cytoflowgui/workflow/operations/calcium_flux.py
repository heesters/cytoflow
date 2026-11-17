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

"""
cytoflowgui.workflow.operations.calcium_flux
--------------------------------------------

"""

from traits.api import provides, Str, Float, Bool, Instance

from cytoflow.operations.calcium_flux import CalciumFluxOp, CalciumFluxView

from ..views import IWorkflowView, WorkflowView
from ..serialization import camel_registry, traits_str, cytoflow_class_repr, dedent

from .operation_base import IWorkflowOperation, WorkflowOperation

CalciumFluxOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class CalciumFluxWorkflowOp(WorkflowOperation, CalciumFluxOp):

    name = Str(apply=True)
    time_channel = Str(apply=True)
    signal_channel = Str(apply=True)
    ratio_denominator = Str("", apply=True)
    baseline_end = Float(60.0, apply=True)
    response_threshold = Float(1.5, apply=True)
    bin_width = Float(2.0, apply=True)

    def default_view(self, **kwargs):
        return CalciumFluxWorkflowView(op=self, **kwargs)

    def clear_estimate(self):
        # CalciumFluxOp has no estimate step
        return

    def should_apply(self, changed, payload):
        from .. import Changed
        if changed in (Changed.APPLY, Changed.PREV_RESULT):
            return True
        return False

    def should_clear_estimate(self, changed, payload):
        return False

    def get_notebook_code(self, idx):
        op = CalciumFluxOp()
        op.copy_traits(self, op.copyable_trait_names())

        return dedent("""
        op_{idx} = {repr}
        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr=repr(op),
                idx=idx,
                prev_idx=idx - 1))


@provides(IWorkflowView)
class CalciumFluxWorkflowView(WorkflowView, CalciumFluxView):
    plot_params = Instance(object, ())

    def should_plot(self, changed, payload):
        from .. import Changed
        if changed in (Changed.APPLY, Changed.PREV_RESULT):
            return True
        return False

    def get_notebook_code(self, idx):
        view = CalciumFluxView()
        view.copy_traits(self, view.copyable_trait_names())

        return dedent("""
        op_{idx}.default_view({traits}).plot(ex_{idx})
        """
        .format(traits=traits_str(view),
                idx=idx))


### Serialization
@camel_registry.dumper(CalciumFluxWorkflowOp, 'calcium-flux', version=1)
def _dump(op):
    return dict(name=op.name,
                time_channel=op.time_channel,
                signal_channel=op.signal_channel,
                ratio_denominator=op.ratio_denominator,
                baseline_end=op.baseline_end,
                response_threshold=op.response_threshold,
                bin_width=op.bin_width)


@camel_registry.loader('calcium-flux', version=1)
def _load(data, version):
    return CalciumFluxWorkflowOp(**data)


@camel_registry.dumper(CalciumFluxWorkflowView, 'calcium-flux-view', version=1)
def _dump_view(view):
    return dict(op=view.op,
                huefacet=view.huefacet,
                rowfacet=view.rowfacet,
                colfacet=view.colfacet,
                show_error_band=view.show_error_band,
                show_stimulation_line=view.show_stimulation_line)


@camel_registry.loader('calcium-flux-view', version=1)
def _load_view(data, version):
    return CalciumFluxWorkflowView(**data)
