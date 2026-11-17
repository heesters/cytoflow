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

'''
Cell Ontology Annotation
--------------------------

Tag gates or populations with standardized Cell Ontology (CL) identifiers,
making the analysis FAIR (Findable, Accessible, Interoperable, Reusable).

.. object:: Edit Annotations...

    Opens a dialog where you can search the bundled Cell Ontology lookup table
    and assign a CL term to any gate, condition or channel in the current
    experiment.

    The dialog shows a search field on the right to filter CL terms by name or
    marker expression.  Select a population on the left and a CL term on the
    right then click **OK** to record the mapping.

.. note::

    This operation does not modify your data; it stores the annotations in
    experiment metadata for downstream export and reporting.
'''

from traits.api import provides, List
from traitsui.api import View, Item, Label, VGroup, Handler
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

try:
    from pyface.qt import QtWidgets, QtCore  # @UnresolvedImport
except ImportError:
    QtWidgets = None
    QtCore = None

from cytoflow.operations.ontology_annotation import search_cl
from ..workflow.operations import OntologyAnnotationWorkflowOp
from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class OntologyAnnotationHandler(OpHandler):
    """Handler for OntologyAnnotationWorkflowOp."""

    def on_edit_annotations(self, info=None):
        """Open the annotation editor dialog."""
        if QtWidgets is None:
            return

        op = self.model

        # Gather population/condition/channel names from the upstream experiment
        pop_names = []
        if (self.context is not None
                and self.context.previous_wi is not None
                and self.context.previous_wi.result is not None):
            result = self.context.previous_wi.result
            pop_names = sorted(set(result.conditions.keys())
                               | set(result.data.columns))

        dlg = _AnnotationDialog(op, pop_names)
        if dlg.exec_() == QtWidgets.QDialog.Accepted:
            op.annotations = dlg.get_annotations()

    operation_traits_view = \
        View(
            VGroup(
                Label('Use the Edit button to annotate populations with Cell Ontology terms.'),
                Item('handler.on_edit_annotations',
                     show_label=False,
                     editor=__import__('traitsui.api', fromlist=['ButtonEditor']).ButtonEditor(
                         label='Edit Annotations...')),
                show_border=False,
            ),
            shared_op_traits_view)


class _AnnotationDialog(object):
    """Qt dialog for editing ontology annotations."""

    def __init__(self, op, pop_names):
        if QtWidgets is None:
            self._dlg = None
            self._annotations = dict(op.annotations)
            return

        self._op = op
        self._annotations = dict(op.annotations)

        self._dlg = QtWidgets.QDialog()
        self._dlg.setWindowTitle('Cell Ontology Annotation')
        self._dlg.resize(800, 500)

        layout = QtWidgets.QHBoxLayout(self._dlg)

        # ----- Left panel: existing annotations table -----
        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        left_layout.addWidget(QtWidgets.QLabel('<b>Population Annotations</b>'))

        self._table = QtWidgets.QTableWidget(0, 3)
        self._table.setHorizontalHeaderLabels(['Population', 'CL ID', 'CL Label'])
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self._table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        left_layout.addWidget(self._table)

        btn_row = QtWidgets.QHBoxLayout()
        self._pop_combo = QtWidgets.QComboBox()
        for name in pop_names:
            self._pop_combo.addItem(name)
        btn_row.addWidget(self._pop_combo)

        add_btn = QtWidgets.QPushButton('Add / Update')
        add_btn.clicked.connect(self._add_annotation)
        btn_row.addWidget(add_btn)

        del_btn = QtWidgets.QPushButton('Remove')
        del_btn.clicked.connect(self._remove_annotation)
        btn_row.addWidget(del_btn)
        left_layout.addLayout(btn_row)

        layout.addWidget(left, stretch=1)

        # ----- Right panel: CL search -----
        right = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(right)
        right_layout.addWidget(QtWidgets.QLabel('<b>Search Cell Ontology</b>'))

        search_row = QtWidgets.QHBoxLayout()
        self._search_edit = QtWidgets.QLineEdit()
        self._search_edit.setPlaceholderText('Type to search (e.g. CD4, neutrophil)...')
        self._search_edit.textChanged.connect(self._do_search)
        search_row.addWidget(self._search_edit)
        right_layout.addLayout(search_row)

        self._cl_table = QtWidgets.QTableWidget(0, 3)
        self._cl_table.setHorizontalHeaderLabels(['CL ID', 'Label', 'Markers'])
        self._cl_table.horizontalHeader().setStretchLastSection(True)
        self._cl_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self._cl_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        right_layout.addWidget(self._cl_table)

        layout.addWidget(right, stretch=1)

        # ----- Dialog buttons -----
        btn_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        btn_box.accepted.connect(self._dlg.accept)
        btn_box.rejected.connect(self._dlg.reject)

        outer = QtWidgets.QVBoxLayout()
        outer.addLayout(layout)
        outer.addWidget(btn_box)
        self._dlg.setLayout(outer)

        # Populate table from existing annotations
        self._refresh_table()
        # Show all entries initially
        self._do_search('')

    def _refresh_table(self):
        self._table.setRowCount(0)
        for pop, (cl_id, cl_label) in sorted(self._annotations.items()):
            row = self._table.rowCount()
            self._table.insertRow(row)
            self._table.setItem(row, 0, QtWidgets.QTableWidgetItem(pop))
            self._table.setItem(row, 1, QtWidgets.QTableWidgetItem(cl_id))
            self._table.setItem(row, 2, QtWidgets.QTableWidgetItem(cl_label))

    def _do_search(self, text):
        results = search_cl(text, max_results=50) if text.strip() else search_cl('', max_results=50)
        self._cl_table.setRowCount(0)
        for entry in results:
            row = self._cl_table.rowCount()
            self._cl_table.insertRow(row)
            self._cl_table.setItem(row, 0, QtWidgets.QTableWidgetItem(entry['cl_id']))
            self._cl_table.setItem(row, 1, QtWidgets.QTableWidgetItem(entry['label']))
            self._cl_table.setItem(row, 2, QtWidgets.QTableWidgetItem(entry.get('markers', '')))

    def _add_annotation(self):
        pop = self._pop_combo.currentText()
        if not pop:
            return
        rows = self._cl_table.selectedItems()
        if not rows:
            return
        selected_row = self._cl_table.currentRow()
        cl_id = self._cl_table.item(selected_row, 0).text()
        cl_label = self._cl_table.item(selected_row, 1).text()
        self._annotations[pop] = (cl_id, cl_label)
        self._refresh_table()

    def _remove_annotation(self):
        selected_row = self._table.currentRow()
        if selected_row < 0:
            return
        pop = self._table.item(selected_row, 0).text()
        self._annotations.pop(pop, None)
        self._refresh_table()

    def exec_(self):
        if self._dlg is None:
            return 0
        return self._dlg.exec_()

    def get_annotations(self):
        return dict(self._annotations)


@provides(IOperationPlugin)
class OntologyAnnotationPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.ontology_annotation'
    operation_id = 'cytoflow.operations.ontology_annotation'

    short_name = name = "Cell Ontology"
    menu_group = "TOOLS"

    def get_operation(self):
        return OntologyAnnotationWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, OntologyAnnotationWorkflowOp):
            return OntologyAnnotationHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('ontology_annotation')

    plugin = List(contributes_to=OP_PLUGIN_EXT)

    def _plugin_default(self):
        return [self]
