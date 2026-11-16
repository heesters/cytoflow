#!/usr/bin/env python3.11
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
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
cytoflowgui.gate_interchange
----------------------------

The glue between the GUI's workflow and `cytoflow.operations.gatingml`:
importing gates from Gating-ML files and FlowJo workspaces, and exporting
the workflow's gates to Gating-ML.

The conversions need the `Experiment` (to look up channel names and the
parameters of the scales), and that only exists in the remote process.  So
`import_gates` and `export_gates`, which are called in the local (GUI)
process, ask the remote process to run `remote_import` and `remote_export`.
"""

import warnings

import cytoflow.utility as util
from cytoflow.operations.gatingml import (export_gatingml, import_gatingml,
                                          import_flowjo_workspace, is_gate)

from .workflow.workflow_item import WorkflowItem
from .workflow.operations import (ThresholdWorkflowOp, RangeWorkflowOp,
                                  Range2DWorkflowOp, QuadWorkflowOp,
                                  PolygonWorkflowOp, EllipseWorkflowOp,
                                  BooleanGateWorkflowOp)

# operation id --> the workflow operation, and the traits to copy (in the
# order to set them)
GATE_OPS = {
    'cytoflow.operations.threshold' :
        (ThresholdWorkflowOp, ['name', 'channel', 'threshold']),
    'cytoflow.operations.range' :
        (RangeWorkflowOp, ['name', 'channel', 'low', 'high']),
    'cytoflow.operations.range2d' :
        (Range2DWorkflowOp, ['name', 'xchannel', 'ychannel', 'xlow', 'xhigh', 'ylow', 'yhigh']),
    'cytoflow.operations.quad' :
        (QuadWorkflowOp, ['name', 'xchannel', 'ychannel', 'xthreshold', 'ythreshold']),
    'cytoflow.operations.polygon' :
        (PolygonWorkflowOp, ['name', 'xchannel', 'ychannel', 'xscale', 'yscale', 'vertices']),
    'cytoflow.operations.ellipse' :
        (EllipseWorkflowOp, ['name', 'xchannel', 'ychannel', 'xscale', 'yscale',
                             'center', 'width', 'height', 'angle']),
    'cytoflow.operations.boolean_gate' :
        (BooleanGateWorkflowOp, ['name', 'expression'])}


def is_flowjo_workspace(path):
    return str(path).lower().endswith('.wsp')


def _plain(value):
    """Plain python values, to send between the processes"""
    if isinstance(value, (list, tuple)):
        return type(value)([_plain(v) for v in value]) if isinstance(value, tuple) \
               else [_plain(v) for v in value]
    elif isinstance(value, str):
        return str(value)
    return float(value)


def _collect(w):
    return [str(x.message) for x in w if issubclass(x.category, util.CytoflowWarning)]


def remote_import(workflow, path, sample = None, group = None):
    """
    Runs in the remote process.  Returns a dict with ``ops`` (a list of
    ``(operation id, {trait : value})``), ``warnings`` (a list of strings) and
    ``error`` (a string, or `None`)
    """
    try:
        experiment = workflow.workflow[-1].result if workflow.workflow else None
        if experiment is None:
            return dict(ops = [], warnings = [],
                        error = "The last operation in the workflow doesn't have "
                                "a result. Import some data, and fix any errors, "
                                "before importing gates.")

        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter('always')
            if is_flowjo_workspace(path):
                ops = import_flowjo_workspace(path, experiment, sample = sample, group = group)
            else:
                ops = import_gatingml(path, experiment)

        return dict(ops = [(op.id, {t : _plain(getattr(op, t)) for t in GATE_OPS[op.id][1]})
                           for op in ops],
                    warnings = _collect(w),
                    error = None)
    except Exception as e:
        return dict(ops = [], warnings = [], error = str(e))


def remote_export(workflow, path):
    """
    Runs in the remote process.  Returns a dict with ``exported`` (the names
    of the gates), ``warnings`` (a list of strings) and ``error``
    """
    try:
        ops = []
        experiments = []
        not_ready = []
        for wi in workflow.workflow:
            if not is_gate(wi.operation):
                continue
            # the experiment that the gate is applied to
            experiment = wi.previous_wi.result if wi.previous_wi else None
            if experiment is None:
                not_ready.append(wi.operation.name or wi.operation.friendly_id)
            else:
                ops.append(wi.operation)
                experiments.append(experiment)

        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter('always')
            exported = export_gatingml(ops, experiments, path)

        warns = _collect(w)
        if not_ready:
            warns.append("These gates were not exported, because the operations "
                         "before them have errors: {}".format(", ".join(not_ready)))
        return dict(exported = exported, warnings = warns, error = None)
    except Exception as e:
        return dict(exported = [], warnings = [], error = str(e))


def _remote(workflow, func, *args):
    result = workflow.remote_eval(
        "__import__('cytoflowgui.gate_interchange', fromlist = ['{0}']).{0}(self, {1})"
        .format(func, ", ".join([repr(a) for a in args])))
    if result['error']:
        raise util.CytoflowError(result['error'])
    return result


def import_gates(workflow, path, sample = None, group = None):
    """
    Import the gates in a Gating-ML file or a FlowJo workspace, and add
    them to the end of the workflow.

    Parameters
    ----------
    workflow : LocalWorkflow

    path : str
        A Gating-ML file, or (if it ends in ``.wsp``) a FlowJo workspace

    sample, group : str
        For a FlowJo workspace: which sample's, or group's, gates to import.

    Returns
    -------
    (list of WorkflowItem, list of str)
        The new workflow items, and any warnings from the import.
    """
    result = _remote(workflow, 'remote_import', str(path), sample, group)

    new_items = []
    for op_id, traits in result['ops']:
        op_class, names = GATE_OPS[op_id]
        op = op_class()
        for name in names:
            setattr(op, name, traits[name])

        wi = WorkflowItem(operation = op, workflow = workflow)
        workflow.workflow.append(wi)
        new_items.append(wi)

    return new_items, result['warnings']


def export_gates(workflow, path):
    """
    Write all of the gates in the workflow to a Gating-ML file.

    Returns
    -------
    (list of str, list of str)
        The names of the gates that were written, and any warnings.
    """
    result = _remote(workflow, 'remote_export', str(path))
    return result['exported'], result['warnings']
