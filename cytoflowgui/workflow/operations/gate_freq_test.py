#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.workflow.operations.gate_freq_test
----------------------------------------------
"""

from traits.api import Str, List, Enum, Float, provides

from cytoflow.operations.gate_freq_test import GateFrequencyTestOp

from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr
from .operation_base import IWorkflowOperation, WorkflowOperation

GateFrequencyTestOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class GateFrequencyTestWorkflowOp(WorkflowOperation, GateFrequencyTestOp):
    name        = Str(apply=True)
    gate        = Str(apply=True)
    parent      = List(Str, apply=True)
    by          = List(Str, apply=True)
    condition   = Str(apply=True)
    condition_a = Str(apply=True)
    condition_b = Str(apply=True)
    test        = Enum('Wilcoxon', 'Paired t-test', apply=True)
    alpha       = Float(0.05, apply=True)

    def default_view(self, **kwargs):
        return None

    def clear_estimate(self):
        return


@camel_registry.dumper(GateFrequencyTestWorkflowOp, 'GateFrequencyTestWorkflowOp', version=1)
def _dump(op):
    return dict(name=op.name, gate=op.gate, parent=op.parent, by=op.by,
                condition=op.condition, condition_a=op.condition_a,
                condition_b=op.condition_b, test=op.test, alpha=op.alpha)


@camel_registry.loader('GateFrequencyTestWorkflowOp', version=1)
def _load(data, version):
    return GateFrequencyTestWorkflowOp(
        name=data.get('name', ''),
        gate=data.get('gate', ''),
        parent=data.get('parent', []),
        by=data.get('by', []),
        condition=data.get('condition', ''),
        condition_a=data.get('condition_a', ''),
        condition_b=data.get('condition_b', ''),
        test=data.get('test', 'Wilcoxon'),
        alpha=data.get('alpha', 0.05))
