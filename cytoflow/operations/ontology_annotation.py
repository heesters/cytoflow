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
cytoflow.operations.ontology_annotation
-----------------------------------------

Tag gates/populations with standardized Cell Ontology (CL) identifiers.

`OntologyAnnotationOp` -- the `IOperation` class.
`search_cl` -- fuzzy search the bundled CL phenotype lookup table.
`get_cl_phenotypes` -- return the full lookup table.
"""

import json
import os
import warnings

from traits.api import HasStrictTraits, Str, Dict, Tuple, Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation

# ---------------------------------------------------------------------------
# Module-level lookup table (loaded once)
# ---------------------------------------------------------------------------

_CL_PHENOTYPES = None


def _load_cl_phenotypes():
    global _CL_PHENOTYPES
    if _CL_PHENOTYPES is None:
        resources_dir = os.path.join(os.path.dirname(__file__), '..', 'resources')
        json_path = os.path.join(resources_dir, 'cl_phenotypes.json')
        with open(os.path.normpath(json_path), 'r', encoding='utf-8') as fh:
            _CL_PHENOTYPES = json.load(fh)
    return _CL_PHENOTYPES


def get_cl_phenotypes():
    """Return the full bundled CL phenotype lookup table.

    Returns
    -------
    list of dict
        Each dict has keys: ``cl_id``, ``label``, ``synonyms``, ``markers``.
    """
    return _load_cl_phenotypes()


def search_cl(query, max_results=10):
    """Search the bundled CL phenotypes by label, synonym or marker string.

    Uses simple case-insensitive substring matching; no external dependencies.

    Parameters
    ----------
    query : str
        Search string.
    max_results : int, optional
        Maximum number of results to return (default 10).

    Returns
    -------
    list of dict
        Matching entries from the lookup table.
    """
    query_lower = query.lower()
    phenotypes = _load_cl_phenotypes()
    results = []
    seen_ids = set()
    for entry in phenotypes:
        cl_id = entry.get('cl_id', '')
        if cl_id in seen_ids:
            continue
        label = entry.get('label', '')
        synonyms = entry.get('synonyms', [])
        markers = entry.get('markers', '')
        hit = (query_lower in label.lower()
               or any(query_lower in s.lower() for s in synonyms)
               or query_lower in markers.lower()
               or query_lower in cl_id.lower())
        if hit:
            results.append(entry)
            seen_ids.add(cl_id)
        if len(results) >= max_results:
            break
    return results


# ---------------------------------------------------------------------------
# Operation
# ---------------------------------------------------------------------------

@provides(IOperation)
class OntologyAnnotationOp(HasStrictTraits):
    """
    Tag gates or populations in an experiment with standardized Cell Ontology
    (CL) identifiers, making the analysis FAIR-compliant.

    This operation does **not** modify ``experiment.data``; it only stores
    annotations in ``experiment.metadata['ontology_annotations']``.

    Attributes
    ----------
    annotations : Dict(Str, Tuple(Str, Str))
        Mapping from a gate/condition/channel name in the experiment to a
        ``(CL_ID, CL_label)`` tuple.  For example::

            {"Live_CD4": ("CL:0000492", "CD4-positive helper T cell")}

    Examples
    --------

    >>> import cytoflow as flow
    >>> op = flow.OntologyAnnotationOp(
    ...     annotations={"Live_CD4": ("CL:0000492", "CD4-positive helper T cell")})
    >>> ex2 = op.apply(ex)
    >>> ex2.metadata['ontology_annotations']
    {'Live_CD4': ('CL:0000492', 'CD4-positive helper T cell')}
    """

    id = Constant('cytoflow.operations.ontology_annotation')
    friendly_id = Constant('Cell Ontology Annotation')
    name = Constant('Cell Ontology Annotation')

    # Dict mapping condition/channel name -> (CL_ID, CL_label)
    annotations = Dict(Str, Tuple(Str, Str))

    def apply(self, experiment):
        """
        Store ontology annotations in experiment metadata.

        Parameters
        ----------
        experiment : Experiment
            The experiment to annotate.

        Returns
        -------
        Experiment
            A clone of *experiment* with ``metadata['ontology_annotations']``
            populated.

        Raises
        ------
        CytoflowOpError
            If *experiment* is ``None``.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        # Warn about unknown names but do not raise an error
        all_names = (set(experiment.conditions.keys())
                     | set(experiment.data.columns))
        for name in self.annotations:
            if name not in all_names:
                warnings.warn(
                    "OntologyAnnotationOp: '{}' is not a known condition or "
                    "channel in the experiment.".format(name),
                    util.CytoflowOpWarning)

        new_experiment = experiment.clone(deep=False)
        new_experiment.metadata['ontology_annotations'] = dict(self.annotations)
        new_experiment.history.append(
            self.clone_traits(transient=lambda _: True))
        return new_experiment


util.expand_class_attributes(OntologyAnnotationOp)
util.expand_method_parameters(OntologyAnnotationOp, OntologyAnnotationOp.apply)
