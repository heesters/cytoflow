#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.workflow.operations.fold_change
-------------------------------------------
"""

from traits.api import Str, List, Bool, provides

from cytoflow.operations.fold_change import FoldChangeOp

from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr
from .operation_base import IWorkflowOperation, WorkflowOperation

FoldChangeOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class FoldChangeWorkflowOp(WorkflowOperation, FoldChangeOp):
    name        = Str(apply=True)
    statistic   = Str(apply=True)
    column      = Str(apply=True)
    by          = List(Str, apply=True)
    condition   = Str(apply=True)
    condition_a = Str(apply=True)
    condition_b = Str(apply=True)
    log2        = Bool(True, apply=True)

    def default_view(self, **kwargs):
        return None

    def clear_estimate(self):
        return


@camel_registry.dumper(FoldChangeWorkflowOp, 'FoldChangeWorkflowOp', version=1)
def _dump(op):
    return dict(name=op.name, statistic=op.statistic, column=op.column,
                by=op.by, condition=op.condition, condition_a=op.condition_a,
                condition_b=op.condition_b, log2=op.log2)


@camel_registry.loader('FoldChangeWorkflowOp', version=1)
def _load(data, version):
    return FoldChangeWorkflowOp(
        name=data.get('name', ''),
        statistic=data.get('statistic', ''),
        column=data.get('column', ''),
        by=data.get('by', []),
        condition=data.get('condition', ''),
        condition_a=data.get('condition_a', ''),
        condition_b=data.get('condition_b', ''),
        log2=data.get('log2', True))
