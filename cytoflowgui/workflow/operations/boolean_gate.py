#!/usr/bin/env python3.11

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2025
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
cytoflowgui.workflow.operations.boolean_gate
--------------------------------------------

"""

from textwrap import dedent
from traits.api import Str, provides

from cytoflow import BooleanGateOp

from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr
from .operation_base import IWorkflowOperation, WorkflowOperation

BooleanGateOp.__repr__ = cytoflow_class_repr

@provides(IWorkflowOperation)
class BooleanGateWorkflowOp(WorkflowOperation, BooleanGateOp):

    name = Str(apply = True)
    expression = Str(apply = True)

    def clear_estimate(self):
        # no-op
        return

    def get_notebook_code(self, idx):
        op = BooleanGateOp()
        op.copy_traits(self, op.copyable_trait_names())

        return dedent("""
        op_{idx} = {repr}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1))


### Serialization
@camel_registry.dumper(BooleanGateWorkflowOp, 'boolean-gate', version = 1)
def _dump(op):
    return dict(name = op.name,
                expression = op.expression)

@camel_registry.loader('boolean-gate', version = 1)
def _load(data, version):
    return BooleanGateWorkflowOp(**data)
