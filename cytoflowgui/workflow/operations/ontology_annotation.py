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
cytoflowgui.workflow.operations.ontology_annotation
-----------------------------------------------------
"""

from traits.api import provides, Str, Dict, Tuple

from cytoflow.operations.ontology_annotation import OntologyAnnotationOp

from ..serialization import camel_registry, cytoflow_class_repr, dedent

from .operation_base import IWorkflowOperation, WorkflowOperation

OntologyAnnotationOp.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class OntologyAnnotationWorkflowOp(WorkflowOperation, OntologyAnnotationOp):
    annotations = Dict(Str, Tuple(Str, Str), apply=True)

    def should_apply(self, changed, payload):
        return True

    def should_clear_estimate(self, changed, payload):
        return False

    def clear_estimate(self):
        pass

    def get_notebook_code(self, idx):
        op = OntologyAnnotationOp()
        op.copy_traits(self, op.copyable_trait_names())

        return dedent("""
        op_{idx} = {repr}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr=repr(op),
                idx=idx,
                prev_idx=idx - 1))


### Serialization

@camel_registry.dumper(OntologyAnnotationWorkflowOp, 'ontology-annotation', version=1)
def _dump(op):
    # Tuple values need to be serialized as lists for YAML compatibility
    return dict(annotations={k: list(v) for k, v in op.annotations.items()})


@camel_registry.loader('ontology-annotation', version=1)
def _load(data, version):
    annotations = {k: tuple(v) for k, v in data.get('annotations', {}).items()}
    return OntologyAnnotationWorkflowOp(annotations=annotations)
