#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.stats_pane
-----------------------

A dock pane that shows live gate frequencies (% parent, % total, count)
for every boolean gate in the workflow, automatically updated as the
analysis runs.  This mirrors the FlowJo statistics table.
"""

import math
import logging

import numpy as np

from traits.api import Instance, observe
from pyface.tasks.api import TraitsDockPane  # @UnresolvedImport
from pyface.qt import QtGui, QtCore

from .workflow import LocalWorkflow
from .workflow_controller import WorkflowController

logger = logging.getLogger(__name__)


class StatsDockPane(TraitsDockPane):
    """Live gate-statistics panel.

    Shows % Parent, % Total and Count for every boolean gate in the
    currently-selected WorkflowItem, broken down by sample-group
    conditions (Donor, Treatment, Timepoint, …).
    """

    id   = 'cytoflowgui.stats_pane'
    name = 'Gate Statistics'

    model   = Instance(LocalWorkflow)
    handler = Instance(WorkflowController)

    closable  = True
    dock_area = 'right'
    floatable = True
    movable   = True
    visible   = True

    # ------------------------------------------------------------------ #
    #  Qt widget creation
    # ------------------------------------------------------------------ #

    def create_contents(self, parent):
        container = QtGui.QWidget(parent)                # @UndefinedVariable
        layout    = QtGui.QVBoxLayout(container)         # @UndefinedVariable
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        # ---- header row ------------------------------------------------
        header = QtGui.QHBoxLayout()                     # @UndefinedVariable

        self._stat_label = QtGui.QLabel("Show:")         # @UndefinedVariable
        header.addWidget(self._stat_label)

        self._stat_combo = QtGui.QComboBox()             # @UndefinedVariable
        self._stat_combo.addItems(['% Parent', '% Total', 'Count'])
        self._stat_combo.currentIndexChanged.connect(self._refresh_table)
        header.addWidget(self._stat_combo)
        header.addStretch(1)

        self._source_label = QtGui.QLabel("Source:")     # @UndefinedVariable
        header.addWidget(self._source_label)

        self._source_combo = QtGui.QComboBox()           # @UndefinedVariable
        self._source_combo.addItems(['Selected step', 'Last step'])
        self._source_combo.currentIndexChanged.connect(self._refresh_table)
        header.addWidget(self._source_combo)

        layout.addLayout(header)

        # ---- table -----------------------------------------------------
        self._table = QtGui.QTableWidget(0, 0)           # @UndefinedVariable
        self._table.setAlternatingRowColors(True)
        self._table.setEditTriggers(
            QtGui.QTableWidget.NoEditTriggers)           # @UndefinedVariable
        self._table.horizontalHeader().setStretchLastSection(False)
        self._table.horizontalHeader().setSectionResizeMode(
            QtGui.QHeaderView.ResizeToContents)          # @UndefinedVariable
        self._table.verticalHeader().setVisible(False)
        layout.addWidget(self._table)

        # ---- placeholder -----------------------------------------------
        self._placeholder = QtGui.QLabel(                # @UndefinedVariable
            "No gates found.\n\nAdd a gate to see frequencies here.")
        self._placeholder.setAlignment(
            QtCore.Qt.AlignCenter)                       # @UndefinedVariable
        self._placeholder.setStyleSheet("color: #888; font-style: italic;")
        layout.addWidget(self._placeholder)
        self._placeholder.hide()

        self.control = container

        # ---- observe model changes ------------------------------------
        if self.model:
            self.model.observe(self._on_selected_changed,
                               'selected', dispatch='ui')
            self.model.observe(self._on_stats_changed,
                               'workflow:items:auto_gate_stats',
                               dispatch='ui')

        self._refresh_table()
        return container

    # ------------------------------------------------------------------ #
    #  Trait observers
    # ------------------------------------------------------------------ #

    @observe('model')
    def _model_changed(self, _):
        if self.model:
            self.model.observe(self._on_selected_changed,
                               'selected', dispatch='ui')
            self.model.observe(self._on_stats_changed,
                               'workflow:items:auto_gate_stats',
                               dispatch='ui')
        self._refresh_table()

    def _on_selected_changed(self, _):
        self._refresh_table()

    def _on_stats_changed(self, _):
        self._refresh_table()

    # ------------------------------------------------------------------ #
    #  Table population
    # ------------------------------------------------------------------ #

    def _get_df(self):
        """Return the auto_gate_stats DataFrame to display, or None."""
        if self.model is None:
            return None

        use_last = getattr(self, '_source_combo', None) and \
                   self._source_combo.currentIndex() == 1

        if use_last and self.model.workflow:
            wi = self.model.workflow[-1]
        elif self.model.selected:
            wi = self.model.selected
        else:
            return None

        return getattr(wi, 'auto_gate_stats', None)

    def _refresh_table(self, *_):
        if not hasattr(self, '_table'):
            return

        df = self._get_df()

        if df is None or df.empty:
            self._table.hide()
            self._placeholder.show()
            return

        stat_col = self._stat_combo.currentText() \
                   if hasattr(self, '_stat_combo') else '% Parent'

        # determine pivot structure
        stat_cols  = ['Count', '% Total', '% Parent']
        group_cols = [c for c in df.columns
                      if c not in stat_cols and c != 'Gate']

        self._placeholder.hide()
        self._table.show()

        if group_cols:
            try:
                pivot = df.pivot_table(index='Gate',
                                       columns=group_cols,
                                       values=stat_col,
                                       aggfunc='first')
            except Exception:
                pivot = None
        else:
            pivot = None

        if pivot is not None:
            gates   = list(pivot.index)
            samples = [str(c) for c in pivot.columns]

            self._table.setRowCount(len(gates))
            self._table.setColumnCount(len(samples) + 1)
            self._table.setHorizontalHeaderLabels(['Gate'] + samples)

            for r, gate in enumerate(gates):
                gate_item = QtGui.QTableWidgetItem(str(gate))  # @UndefinedVariable
                gate_item.setFont(_bold_font())
                self._table.setItem(r, 0, gate_item)
                for c, sample in enumerate(pivot.columns):
                    val = pivot.loc[gate, sample]
                    self._table.setItem(r, c + 1, _fmt_item(val, stat_col))

        else:
            # no group conditions: simple two-column table (Gate, value)
            gates = list(df['Gate'].unique())
            self._table.setRowCount(len(gates))
            self._table.setColumnCount(2)
            self._table.setHorizontalHeaderLabels(['Gate', stat_col])

            for r, gate in enumerate(gates):
                row = df[df['Gate'] == gate].iloc[0]
                gate_item = QtGui.QTableWidgetItem(str(gate))  # @UndefinedVariable
                gate_item.setFont(_bold_font())
                self._table.setItem(r, 0, gate_item)
                self._table.setItem(r, 1, _fmt_item(row[stat_col], stat_col))

        self._table.resizeColumnsToContents()


# ------------------------------------------------------------------ #
#  Helpers
# ------------------------------------------------------------------ #

def _bold_font():
    f = QtGui.QFont()                                    # @UndefinedVariable
    f.setBold(True)
    return f


def _fmt_item(val, stat_col):
    if val is None or (isinstance(val, float) and math.isnan(val)):
        txt = '—'
    elif stat_col == 'Count':
        txt = str(int(val))
    else:
        txt = f'{val:.1f} %'
    item = QtGui.QTableWidgetItem(txt)                   # @UndefinedVariable
    item.setTextAlignment(QtCore.Qt.AlignRight |          # @UndefinedVariable
                          QtCore.Qt.AlignVCenter)        # @UndefinedVariable
    return item
