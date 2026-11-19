#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.population_tree_pane
---------------------------------

A dock pane showing all boolean gate populations as a hierarchy tree with
event counts, % parent and % total.  Auto-updates when the workflow changes.
"""

import logging

from traits.api import Instance, observe
from pyface.tasks.api import TraitsDockPane  # @UnresolvedImport
from pyface.qt import QtGui, QtCore

from .workflow import LocalWorkflow
from .workflow_controller import WorkflowController

logger = logging.getLogger(__name__)


class PopulationTreeDockPane(TraitsDockPane):
    """Population hierarchy tree panel.

    Shows all boolean gate conditions in the currently-selected WorkflowItem
    as a tree.  Root nodes are top-level gates (no inferred parent); child
    nodes are gates whose immediate predecessor boolean gate is the parent.
    Columns: Gate name | Count | % Parent | % Total.
    """

    id   = 'cytoflowgui.population_tree_pane'
    name = 'Population Tree'

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
        container = QtGui.QWidget(parent)
        layout    = QtGui.QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        # ---- tree widget -----------------------------------------------
        self._tree = QtGui.QTreeWidget()
        self._tree.setColumnCount(4)
        self._tree.setHeaderLabels(['Gate', 'Count', '% Parent', '% Total'])
        self._tree.setAlternatingRowColors(True)
        self._tree.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        header = self._tree.header()
        header.setSectionResizeMode(0, QtGui.QHeaderView.Stretch)
        for col in (1, 2, 3):
            header.setSectionResizeMode(col, QtGui.QHeaderView.ResizeToContents)
        layout.addWidget(self._tree)

        # ---- placeholder -----------------------------------------------
        self._placeholder = QtGui.QLabel(
            "No gates found.\n\nAdd a gate operation to see\nthe population hierarchy here.")
        self._placeholder.setAlignment(QtCore.Qt.AlignCenter)
        self._placeholder.setStyleSheet("color: #888; font-style: italic;")
        layout.addWidget(self._placeholder)
        self._placeholder.hide()

        if self.model:
            self.model.observe(self._on_selected_changed,
                               'selected', dispatch='ui')
            self.model.observe(self._on_stats_changed,
                               'workflow:items:auto_gate_stats',
                               dispatch='ui')

        self._refresh_tree()
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
        self._refresh_tree()

    def _on_selected_changed(self, _):
        self._refresh_tree()

    def _on_stats_changed(self, _):
        self._refresh_tree()

    # ------------------------------------------------------------------ #
    #  Tree population
    # ------------------------------------------------------------------ #

    def _get_df(self):
        if self.model is None:
            return None
        if self.model.selected:
            wi = self.model.selected
        elif self.model.workflow:
            wi = self.model.workflow[-1]
        else:
            return None
        return getattr(wi, 'auto_gate_stats', None)

    def _refresh_tree(self, *_):
        if not hasattr(self, '_tree'):
            return

        df = self._get_df()

        if df is None or df.empty:
            self._tree.hide()
            self._placeholder.show()
            return

        self._placeholder.hide()
        self._tree.show()
        self._tree.clear()

        # determine stat and group columns
        stat_cols  = {'Count', '% Total', '% Parent'}
        group_cols = [c for c in df.columns if c not in stat_cols and c != 'Gate']

        # aggregate: for each gate sum/mean across sample groups (or use first
        # group if no grouping) so we show one tree row per gate
        agg = {}
        for gate, gdf in df.groupby('Gate', sort=False):
            count   = gdf['Count'].sum()   if 'Count'    in gdf.columns else float('nan')
            pct_p   = gdf['% Parent'].mean() if '% Parent' in gdf.columns else float('nan')
            pct_t   = gdf['% Total'].mean()  if '% Total'  in gdf.columns else float('nan')
            agg[gate] = (count, pct_p, pct_t)

        gates = list(agg.keys())

        # build parent map from df order (simple sequential hierarchy:
        # each gate's parent is the most recently seen gate before it)
        parent_map = {}
        seen = []
        for gate in gates:
            parent_map[gate] = seen[-1] if seen else None
            seen.append(gate)

        # build QTreeWidgetItems
        items = {}
        for gate in gates:
            count, pct_p, pct_t = agg[gate]
            item = QtGui.QTreeWidgetItem([
                gate,
                '{:.0f}'.format(count)       if not _isnan(count) else '—',
                '{:.1f}%'.format(pct_p)      if not _isnan(pct_p) else '—',
                '{:.1f}%'.format(pct_t)      if not _isnan(pct_t) else '—',
            ])
            # right-align numeric columns
            for col in (1, 2, 3):
                item.setTextAlignment(col, QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            items[gate] = item

        # attach items to their parents (or to tree root)
        for gate, item in items.items():
            parent = parent_map.get(gate)
            if parent and parent in items:
                items[parent].addChild(item)
            else:
                self._tree.addTopLevelItem(item)

        self._tree.expandAll()


def _isnan(v):
    try:
        import math
        return math.isnan(v)
    except Exception:
        return False
