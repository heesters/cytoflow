#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.editors.spillover_matrix_editor
--------------------------------------------

A TraitsUI editor that renders a List(Spillover) as a proper square
compensation matrix.  The diagonal is grayed out and shows "1.0";
off-diagonal cells are editable inline as percentages.
"""

from natsort import natsorted

from pyface.qt import QtGui, QtCore

from traits.api import Bool
from traitsui.editor_factory import EditorFactory
from traitsui.qt.editor import Editor


# ---------------------------------------------------------------------------
# Actual Qt editor
# ---------------------------------------------------------------------------

class _SpilloverMatrixEditor(Editor):
    """Qt editor that displays spillover_list as a square grid."""

    _updating = Bool(False)

    def init(self, parent):
        table = QtGui.QTableWidget(0, 0)                    # @UndefinedVariable
        table.setAlternatingRowColors(False)
        table.horizontalHeader().setSectionResizeMode(
            QtGui.QHeaderView.ResizeToContents)             # @UndefinedVariable
        table.verticalHeader().setSectionResizeMode(
            QtGui.QHeaderView.ResizeToContents)             # @UndefinedVariable
        table.setSelectionMode(
            QtGui.QAbstractItemView.SingleSelection)        # @UndefinedVariable
        table.setTabKeyNavigation(True)
        table.itemChanged.connect(self._on_item_changed)

        self.control = table
        self._rebuild_table()

    def update_editor(self):
        self._rebuild_table()

    # ------------------------------------------------------------------ #

    def _channels(self):
        """Return sorted channel names from current spillover_list."""
        channels = set()
        for s in (self.value or []):
            channels.add(s.from_channel)
            channels.add(s.to_channel)
        return natsorted(channels)

    def _rebuild_table(self):
        table = self.control
        if table is None:
            return

        self._updating = True
        try:
            channels = self._channels()
            n = len(channels)

            table.blockSignals(True)
            table.setRowCount(n)
            table.setColumnCount(n)
            table.setHorizontalHeaderLabels(channels)
            table.setVerticalHeaderLabels(channels)

            # build lookup: (from, to) → Spillover object
            lookup = {}
            for s in (self.value or []):
                lookup[(s.from_channel, s.to_channel)] = s

            gray = QtGui.QColor(220, 220, 220)              # @UndefinedVariable
            diagonal_flags = (QtCore.Qt.ItemIsEnabled |     # @UndefinedVariable
                              QtCore.Qt.ItemIsSelectable)   # @UndefinedVariable
            editable_flags = (QtCore.Qt.ItemIsEnabled |     # @UndefinedVariable
                              QtCore.Qt.ItemIsSelectable |  # @UndefinedVariable
                              QtCore.Qt.ItemIsEditable)     # @UndefinedVariable

            for r, from_ch in enumerate(channels):
                for c, to_ch in enumerate(channels):
                    if from_ch == to_ch:
                        item = QtGui.QTableWidgetItem("1.0")  # @UndefinedVariable
                        item.setFlags(diagonal_flags)
                        item.setBackground(gray)
                        item.setForeground(
                            QtGui.QColor(140, 140, 140))    # @UndefinedVariable
                        item.setData(QtCore.Qt.UserRole, None)  # @UndefinedVariable
                    else:
                        spill_obj = lookup.get((from_ch, to_ch))
                        val = spill_obj.spillover if spill_obj else 0.0
                        item = QtGui.QTableWidgetItem(f"{val:.4g}")  # @UndefinedVariable
                        item.setFlags(editable_flags)
                        item.setTextAlignment(
                            QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)  # @UndefinedVariable
                        item.setData(QtCore.Qt.UserRole,    # @UndefinedVariable
                                     (from_ch, to_ch))

                    table.setItem(r, c, item)

            table.blockSignals(False)
        finally:
            self._updating = False

    def _on_item_changed(self, item):
        if self._updating:
            return

        key = item.data(QtCore.Qt.UserRole)                 # @UndefinedVariable
        if key is None:
            return  # diagonal — ignore

        from_ch, to_ch = key
        try:
            new_val = float(item.text())
        except (ValueError, TypeError):
            # Revert to old value
            self._updating = True
            for s in (self.value or []):
                if s.from_channel == from_ch and s.to_channel == to_ch:
                    item.setText(f"{s.spillover:.4g}")
                    break
            self._updating = False
            return

        for s in (self.value or []):
            if s.from_channel == from_ch and s.to_channel == to_ch:
                if s.spillover != new_val:
                    s.spillover = new_val
                break

    def dispose(self):
        if self.control is not None:
            try:
                self.control.itemChanged.disconnect(self._on_item_changed)
            except Exception:
                pass
        super().dispose()


# ---------------------------------------------------------------------------
# Editor factory
# ---------------------------------------------------------------------------

class SpilloverMatrixEditor(EditorFactory):
    """Factory for the spillover matrix editor."""

    def simple_editor(self, ui, object, name, description, parent):
        return _SpilloverMatrixEditor(
            parent, factory=self, ui=ui, object=object,
            name=name, description=description)

    def custom_editor(self, ui, object, name, description, parent):
        return self.simple_editor(ui, object, name, description, parent)

    def text_editor(self, ui, object, name, description, parent):
        return self.simple_editor(ui, object, name, description, parent)

    def readonly_editor(self, ui, object, name, description, parent):
        return self.simple_editor(ui, object, name, description, parent)
