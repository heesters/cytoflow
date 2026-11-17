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

"""
cytoflowgui.workflow.operations.sample_labels
----------------------------------------------
"""

from traits.api import provides, Str, Dict

from cytoflow.operations.sample_labels import SampleLabelsOp

from ..serialization import camel_registry, cytoflow_class_repr, dedent

from .operation_base import IWorkflowOperation, WorkflowOperation

SampleLabelsOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class SampleLabelsWorkflowOp(WorkflowOperation, SampleLabelsOp):
    condition = Str(apply=True)
    labels = Dict(Str, Str, apply=True)
    new_condition = Str("", apply=True)

    def clear_estimate(self):
        pass

    def get_notebook_code(self, idx):
        op = SampleLabelsOp()
        op.copy_traits(self, op.copyable_trait_names())

        return dedent("""
        op_{idx} = {repr}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr=repr(op),
                idx=idx,
                prev_idx=idx - 1))


### Serialization

@camel_registry.dumper(SampleLabelsWorkflowOp, 'sample-labels', version=1)
def _dump(op):
    return dict(condition=op.condition,
                labels=dict(op.labels),
                new_condition=op.new_condition)


@camel_registry.loader('sample-labels', version=1)
def _load(data, version):
    return SampleLabelsWorkflowOp(**data)
