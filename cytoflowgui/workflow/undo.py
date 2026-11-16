#!/usr/bin/env python3
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
cytoflowgui.workflow.undo
-------------------------

Undo / redo for the GUI's workflow.

`UndoManager` watches a `LocalWorkflow` and records:

- changes to operation parameters (any trait with ``apply`` or ``estimate``
  metadata -- including gates drawn on a plot, which end up setting the
  operation's traits in the GUI process),
- changes to view parameters (non-transient, non-status traits), and
- operations being added to or removed from the workflow.

Undoing a step simply sets the old values again (or re-inserts / removes the
workflow item), and the normal machinery sends the change to the remote
process and re-runs whatever needs re-running.

Changes that arrive within `MERGE_SECONDS` of each other are merged into
one undo step, so that (for example) dragging a range gate -- which sets
both ``low`` and ``high`` -- is undone in one go.
"""

import copy
import logging
import threading
import time

from traits.api import HasStrictTraits, Instance, Bool, List, Any, Property, observe

logger = logging.getLogger(__name__)

MERGE_SECONDS = 0.5
MAX_STEPS = 100


class _Change(object):
    """One trait change: set ``obj.name`` from ``old`` to ``new``"""

    def __init__(self, obj, name, old, new):
        self.obj = obj
        self.name = name
        self.old = old
        self.new = new

    def undo(self, manager):
        setattr(self.obj, self.name, copy.deepcopy(self.old))

    def redo(self, manager):
        setattr(self.obj, self.name, copy.deepcopy(self.new))

    def describe(self):
        return self.name.replace('_', ' ').strip()


class _ItemChange(object):
    """A workflow item was added (``added = True``) or removed"""

    def __init__(self, index, wi, added):
        self.index = index
        self.wi = wi
        self.added = added

    def _add(self, manager):
        workflow = manager.model.workflow
        workflow.insert(min(self.index, len(workflow)), self.wi)

    def _remove(self, manager):
        workflow = manager.model.workflow
        if self.wi in workflow:
            workflow.remove(self.wi)

    def undo(self, manager):
        self._remove(manager) if self.added else self._add(manager)

    def redo(self, manager):
        self._add(manager) if self.added else self._remove(manager)

    def describe(self):
        name = getattr(self.wi.operation, 'friendly_id', 'operation')
        return ("add " if self.added else "delete ") + name


class _Step(object):
    def __init__(self):
        self.changes = []
        self.time = time.monotonic()

    def describe(self):
        if len(self.changes) == 1:
            return self.changes[0].describe()
        names = []
        for c in self.changes:
            d = c.describe()
            if d not in names:
                names.append(d)
        return ", ".join(names[:3]) + ("..." if len(names) > 3 else "")


def _snapshot(value):
    try:
        return copy.deepcopy(value)
    except Exception:
        return value


class UndoManager(HasStrictTraits):
    """
    Undo / redo history for a `LocalWorkflow`.  Create it with
    ``UndoManager(model = local_workflow)``, then call `undo` and `redo`.
    """

    model = Instance('cytoflowgui.workflow.workflow.LocalWorkflow')

    can_undo = Property(Bool, observe = 'undo_stack.items')
    can_redo = Property(Bool, observe = 'redo_stack.items')
    undo_name = Property(observe = 'undo_stack.items')
    redo_name = Property(observe = 'redo_stack.items')

    undo_stack = List(Any)
    redo_stack = List(Any)

    # True while we're undoing or redoing: don't record those changes
    _applying = Bool(False)

    # last-known parameter values, for traits whose change notifications
    # don't include the old value (list / dict editing, the 'changed' event)
    _shadow = Any
    _lock = Any

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._shadow = {}
        self._lock = threading.RLock()
        self.clear()

    def _get_can_undo(self):
        return len(self.undo_stack) > 0

    def _get_can_redo(self):
        return len(self.redo_stack) > 0

    def _get_undo_name(self):
        return self.undo_stack[-1].describe() if self.undo_stack else ""

    def _get_redo_name(self):
        return self.redo_stack[-1].describe() if self.redo_stack else ""

    def clear(self):
        with self._lock:
            self.undo_stack = []
            self.redo_stack = []
            self._shadow = {}
            if self.model:
                for wi in self.model.workflow:
                    self._shadow_object(wi.operation)
                    for v in wi.views:
                        self._shadow_object(v)

    ### recording

    @staticmethod
    def _recordable_op_traits(op):
        return [n for n in op.trait_names()
                if (op.trait(n).apply or op.trait(n).estimate)
                and not op.trait(n).transient
                and not op.trait(n).status]

    @staticmethod
    def _recordable_view_traits(view):
        return [n for n in view.trait_names()
                if not view.trait(n).transient
                and not view.trait(n).status
                and view.trait(n).type not in ('delegate', 'event', 'property')
                and not n.startswith('_')
                and n not in ('trait_added', 'trait_modified')]

    def _shadow_object(self, obj):
        if obj is None:
            return
        try:
            names = self._recordable_op_traits(obj) \
                    if hasattr(obj, 'should_apply') \
                    else self._recordable_view_traits(obj)
            self._shadow[id(obj)] = {n : _snapshot(getattr(obj, n)) for n in names}
        except Exception:
            logger.debug("UndoManager: couldn't snapshot {}".format(obj))

    def _record(self, change):
        if self._applying:
            return

        with self._lock:
            now = time.monotonic()
            if (self.undo_stack
                and now - self.undo_stack[-1].time < MERGE_SECONDS
                and not self.redo_stack):
                step = self.undo_stack[-1]
            else:
                step = _Step()
                self.undo_stack.append(step)
                if len(self.undo_stack) > MAX_STEPS:
                    self.undo_stack.pop(0)

            # merge repeated changes to the same trait in one step
            if isinstance(change, _Change):
                prev = next((c for c in step.changes
                             if isinstance(c, _Change)
                             and c.obj is change.obj
                             and c.name == change.name), None)
                if prev is not None:
                    prev.new = change.new
                    step.time = now
                    self.redo_stack = []
                    return

            step.changes.append(change)
            step.time = now
            self.redo_stack = []

    def _record_trait_change(self, obj, name, old, new, has_old):
        shadow = self._shadow.setdefault(id(obj), {})
        if not has_old:
            if name not in shadow:
                # we don't know what it was before; start tracking it now
                shadow[name] = _snapshot(new)
                return
            old = shadow[name]
        new = _snapshot(new)
        shadow[name] = new
        if old is new or (type(old) is type(new) and _safe_eq(old, new)):
            return
        self._record(_Change(obj, name, _snapshot(old), new))

    @observe('model:workflow', post_init = True)
    def _on_new_workflow(self, _):
        self.clear()

    @observe('model:workflow:items')
    def _on_items(self, event):
        for wi in event.added:
            self._shadow_object(wi.operation)
            for v in wi.views:
                self._shadow_object(v)
        if self._applying:
            return
        for wi in event.removed:
            self._record(_ItemChange(event.index, wi, added = False))
        for wi in event.added:
            self._record(_ItemChange(event.index, wi, added = True))

    @observe('model:workflow:items:views:items')
    def _on_views_added(self, event):
        for v in event.added:
            self._shadow_object(v)

    @observe('model:workflow:items:operation:[+apply,+estimate]')
    def _on_operation_changed(self, event):
        obj = event.object
        name = event.name
        has_old = True
        if name == 'changed':
            name = event.new
            has_old = False
        trait = obj.trait(name)
        if trait is None or trait.transient or trait.status:
            return
        self._record_trait_change(obj, name, event.old, getattr(obj, name), has_old)

    @observe('model:workflow:items:views:items:+type')
    def _on_view_changed(self, event):
        obj = event.object
        name = event.name
        has_old = True
        if name == 'changed':
            name = event.new
            has_old = False
        if name.startswith('_') or name in ('trait_added', 'trait_modified'):
            return
        trait = obj.trait(name)
        if (trait is None or trait.transient or trait.status
            or trait.type in ('delegate', 'event', 'property')):
            return
        self._record_trait_change(obj, name, event.old, getattr(obj, name), has_old)

    ### undo / redo

    def undo(self):
        with self._lock:
            if not self.undo_stack:
                return
            step = self.undo_stack.pop()
            self._apply(step, undo = True)
            self.redo_stack.append(step)

    def redo(self):
        with self._lock:
            if not self.redo_stack:
                return
            step = self.redo_stack.pop()
            self._apply(step, undo = False)
            step.time = 0     # never merge new changes into a redone step
            self.undo_stack.append(step)

    def _apply(self, step, undo):
        self._applying = True
        try:
            changes = reversed(step.changes) if undo else step.changes
            for c in changes:
                try:
                    c.undo(self) if undo else c.redo(self)
                    if isinstance(c, _Change):
                        self._shadow.setdefault(id(c.obj), {})[c.name] = \
                            _snapshot(c.old if undo else c.new)
                except Exception:
                    logger.exception("Couldn't {} {}"
                                     .format("undo" if undo else "redo", c.describe()))
        finally:
            self._applying = False


def _safe_eq(a, b):
    try:
        return bool(a == b)
    except Exception:
        return False
