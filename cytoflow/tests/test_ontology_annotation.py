#!/usr/bin/env python3.8
# coding: latin-1

"""
Tests for OntologyAnnotationOp, search_cl and get_cl_phenotypes.
"""

import os
import tempfile
import unittest

import cytoflow as flow
from cytoflow.operations.ontology_annotation import (
    OntologyAnnotationOp, search_cl, get_cl_phenotypes
)


class TestGetClPhenotypes(unittest.TestCase):
    """Test the get_cl_phenotypes helper."""

    def test_returns_list(self):
        result = get_cl_phenotypes()
        self.assertIsInstance(result, list)

    def test_at_least_80_entries(self):
        result = get_cl_phenotypes()
        self.assertGreaterEqual(len(result), 80,
                                "Expected at least 80 CL phenotype entries, got {}".format(len(result)))

    def test_entry_structure(self):
        """Each entry must have cl_id, label, synonyms and markers keys."""
        for entry in get_cl_phenotypes():
            self.assertIn('cl_id', entry)
            self.assertIn('label', entry)
            self.assertIn('synonyms', entry)
            self.assertIn('markers', entry)


class TestSearchCl(unittest.TestCase):
    """Test the search_cl helper."""

    def test_cd4_returns_cd4_t_cell(self):
        results = search_cl('CD4')
        self.assertTrue(len(results) > 0, "search_cl('CD4') returned no results")
        labels = [r['label'] for r in results]
        synonyms_flat = [s for r in results for s in r.get('synonyms', [])]
        combined = ' '.join(labels + synonyms_flat).lower()
        self.assertIn('cd4', combined)

    def test_neutrophil_returns_neutrophil(self):
        results = search_cl('neutrophil')
        self.assertTrue(len(results) > 0, "search_cl('neutrophil') returned no results")
        labels = [r['label'].lower() for r in results]
        self.assertTrue(any('neutrophil' in l for l in labels),
                        "No neutrophil entry found in results: {}".format(labels))

    def test_max_results_respected(self):
        results = search_cl('cell', max_results=5)
        self.assertLessEqual(len(results), 5)

    def test_empty_query_returns_results(self):
        results = search_cl('', max_results=10)
        self.assertIsInstance(results, list)

    def test_treg_found(self):
        results = search_cl('Treg')
        self.assertTrue(len(results) > 0, "search_cl('Treg') returned no results")
        all_text = ' '.join(
            r['label'] + ' ' + ' '.join(r.get('synonyms', []))
            for r in results
        ).lower()
        self.assertIn('treg', all_text)

    def test_nk_cell_found(self):
        results = search_cl('NK cell')
        self.assertTrue(len(results) > 0)

    def test_case_insensitive(self):
        results_upper = search_cl('CD8')
        results_lower = search_cl('cd8')
        self.assertEqual(len(results_upper), len(results_lower))


class _ExperimentMixin:
    """Mixin to set up a minimal experiment."""

    def setUp(self):
        cwd = os.path.dirname(os.path.abspath(__file__))
        self.ex = flow.ImportOp(
            conditions={'Dox': 'float'},
            tubes=[
                flow.Tube(file=os.path.join(cwd, 'data/Plate01/RFP_Well_A3.fcs'),
                          conditions={'Dox': 10.0}),
                flow.Tube(file=os.path.join(cwd, 'data/Plate01/RFP_Well_A6.fcs'),
                          conditions={'Dox': 1.0}),
            ]).apply()


class TestOntologyAnnotationOp(_ExperimentMixin, unittest.TestCase):
    """Test OntologyAnnotationOp.apply()."""

    def setUp(self):
        super().setUp()

    def test_apply_stores_annotations_in_metadata(self):
        op = OntologyAnnotationOp(
            annotations={'Dox': ('CL:0000492', 'CD4-positive helper T cell')})
        ex2 = op.apply(self.ex)
        self.assertIn('ontology_annotations', ex2.metadata)
        self.assertEqual(
            ex2.metadata['ontology_annotations'],
            {'Dox': ('CL:0000492', 'CD4-positive helper T cell')})

    def test_apply_does_not_modify_data(self):
        op = OntologyAnnotationOp(
            annotations={'Dox': ('CL:0000492', 'CD4-positive helper T cell')})
        ex2 = op.apply(self.ex)
        self.assertTrue(self.ex.data.equals(ex2.data))

    def test_apply_multiple_annotations(self):
        op = OntologyAnnotationOp(
            annotations={
                'Dox': ('CL:0000492', 'CD4-positive helper T cell'),
                'B1-A': ('CL:0000084', 'T cell'),
            })
        ex2 = op.apply(self.ex)
        ann = ex2.metadata['ontology_annotations']
        self.assertIn('Dox', ann)
        self.assertIn('B1-A', ann)

    def test_apply_no_experiment_raises(self):
        import cytoflow.utility as util
        op = OntologyAnnotationOp(
            annotations={'Dox': ('CL:0000492', 'CD4-positive helper T cell')})
        with self.assertRaises(util.CytoflowOpError):
            op.apply(None)

    def test_unknown_name_warns_but_does_not_raise(self):
        import warnings
        import cytoflow.utility as util
        op = OntologyAnnotationOp(
            annotations={'NonExistentGate': ('CL:0000492', 'CD4-positive helper T cell')})
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            ex2 = op.apply(self.ex)
        self.assertIsNotNone(ex2)
        # Should have issued a warning
        self.assertTrue(any(issubclass(x.category, util.CytoflowOpWarning) for x in w),
                        "Expected a CytoflowOpWarning for unknown gate name")

    def test_empty_annotations(self):
        op = OntologyAnnotationOp(annotations={})
        ex2 = op.apply(self.ex)
        self.assertEqual(ex2.metadata.get('ontology_annotations', {}), {})

    def test_history_updated(self):
        op = OntologyAnnotationOp(
            annotations={'Dox': ('CL:0000492', 'CD4-positive helper T cell')})
        ex2 = op.apply(self.ex)
        self.assertTrue(len(ex2.history) > 0)


class TestOntologyAnnotationSerialization(_ExperimentMixin, unittest.TestCase):
    """Test camel round-trip for OntologyAnnotationWorkflowOp."""

    def setUp(self):
        super().setUp()

    def test_camel_roundtrip(self):
        """Workflow op round-trips through camel save/load correctly."""
        try:
            from cytoflowgui.workflow.operations import OntologyAnnotationWorkflowOp
            from cytoflowgui.workflow.serialization import load_yaml, save_yaml
        except (ImportError, Exception) as e:
            self.skipTest("GUI not available: {}".format(e))

        op = OntologyAnnotationWorkflowOp(
            annotations={
                'Dox': ('CL:0000492', 'CD4-positive helper T cell'),
                'B1-A': ('CL:0000084', 'T cell'),
            })

        with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as f:
            fname = f.name

        try:
            save_yaml(op, fname)
            loaded = load_yaml(fname)
        finally:
            os.unlink(fname)

        self.assertIsInstance(loaded, OntologyAnnotationWorkflowOp)
        self.assertEqual(dict(loaded.annotations), dict(op.annotations))


if __name__ == '__main__':
    unittest.main()
