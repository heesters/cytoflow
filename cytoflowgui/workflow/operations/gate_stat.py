#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.workflow.operations.gate_stat
-----------------------------------------
"""

from traits.api import Str, List, provides

from cytoflow.operations.gate_stat import GateStatisticOp

from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr
from .operation_base import IWorkflowOperation, WorkflowOperation

GateStatisticOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class GateStatisticWorkflowOp(WorkflowOperation, GateStatisticOp):
    # all three traits trigger re-computation when changed
    name   = Str(apply=True)
    gate   = Str(apply=True)
    parent = List(Str, apply=True)
    by     = List(Str, apply=True)

    def default_view(self, **kwargs):
        return None


# ---------- serialization ----------

@camel_registry.dumper(GateStatisticWorkflowOp, 'GateStatisticWorkflowOp', version=1)
def _dump(op):
    return dict(name=op.name, gate=op.gate, parent=op.parent, by=op.by)


@camel_registry.loader('GateStatisticWorkflowOp', version=1)
def _load(data, version):
    return GateStatisticWorkflowOp(
        name=data['name'],
        gate=data.get('gate', ''),
        parent=data.get('parent', []),
        by=data.get('by', []))
