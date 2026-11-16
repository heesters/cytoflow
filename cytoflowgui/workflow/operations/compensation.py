#!/usr/bin/env python3.8
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2022
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
cytoflowgui.workflow.operations.compensation
--------------------------------------------

"""

from natsort import natsorted

from traits.api import (HasTraits, provides, Str, observe, Instance,
                        List, Float, File)

from cytoflow.operations.compensation import (CompensationOp,
                                              CompensationDiagnostic,
                                              SpilloverMatrix)
import cytoflow.utility as util

from .. import Changed
from ..views import IWorkflowView, WorkflowView
from ..serialization import camel_registry, traits_str, traits_repr, cytoflow_class_repr, dedent

from .operation_base import IWorkflowOperation, WorkflowOperation

CompensationOp.__repr__ = cytoflow_class_repr


class Spillover(HasTraits):
    """One (off-diagonal) entry of the spillover matrix, for the GUI's table"""
    from_channel = Str
    to_channel = Str
    spillover = Float

    def __repr__(self):
        return traits_repr(self)


def _spillover_to_list(spillover):
    channels = natsorted(spillover.keys())
    return [Spillover(from_channel = f,
                      to_channel = t,
                      spillover = spillover[f].get(t, 0.0))
            for f in channels for t in channels if f != t]


@provides(IWorkflowOperation)
class CompensationWorkflowOp(WorkflowOperation, CompensationOp):
    # add some metadata
    spillover_file = File(estimate = True)
    spillover = SpilloverMatrix(estimate_result = True, status = True)

    # the off-diagonal entries of the spillover matrix, as a list (for the
    # GUI's table editor.)  edited in the local process; sent to the remote
    # process by setting `changed` to 'spillover_list'
    spillover_list = List(Spillover)

    def _list_to_spillover(self):
        channels = list(self.spillover.keys())
        for s in self.spillover_list:
            for c in (s.from_channel, s.to_channel):
                if c not in channels:
                    channels.append(c)

        # keep the diagonal (usually 1.0) from the current matrix
        spillover = {c : {c : self.spillover.get(c, {}).get(c, 1.0)}
                     for c in channels}
        for s in self.spillover_list:
            spillover[s.from_channel][s.to_channel] = s.spillover
        return spillover

    # when the spillover matrix changes (ie, from estimate()), update the list
    @observe('spillover')
    def _on_spillover_changed(self, _):
        new_list = _spillover_to_list(self.spillover)
        if [(s.from_channel, s.to_channel, s.spillover) for s in new_list] != \
           [(s.from_channel, s.to_channel, s.spillover) for s in self.spillover_list]:
            self.spillover_list = new_list

    # when the spillover list is replaced (eg, when it's sent from the local
    # process to the remote process), update the spillover matrix
    @observe('spillover_list', post_init = True)
    def _on_spillover_list_changed(self, _):
        if not self.spillover_list and not self.spillover:
            return

        spillover = self._list_to_spillover()
        if spillover != self.spillover:
            self.spillover = spillover

    # when a value in the table is edited (in the local process), send the
    # list to the remote process
    @observe('spillover_list:items:spillover', post_init = True)
    def _on_spillover_value_changed(self, _):
        self.changed = 'spillover_list'

    def default_view(self, **kwargs):
        return CompensationWorkflowView(op = self, **kwargs)

    def estimate(self, experiment):
        return super().estimate(experiment)

    def apply(self, experiment):
        if not self.spillover:
            raise util.CytoflowOpError('spillover',
                                       "No spillover matrix. Click "
                                       "\"Read spillover matrix\" to read it "
                                       "from the FCS file (or the CSV file).")

        return super().apply(experiment)

    def should_clear_estimate(self, changed, payload):
        if changed == Changed.ESTIMATE:
            return True

        return False

    def clear_estimate(self):
        self.spillover = {}

    def get_notebook_code(self, idx):
        # don't re-read the matrix in the notebook -- the user may have
        # edited it.
        op = CompensationOp()
        op.spillover = self.spillover

        return dedent("""
        op_{idx} = {repr}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1))


@provides(IWorkflowView)
class CompensationWorkflowView(WorkflowView, CompensationDiagnostic):
    plot_params = Instance(HasTraits, ())

    def should_plot(self, changed, payload):
        if changed == Changed.ESTIMATE_RESULT:
            return True

        return False

    def get_notebook_code(self, idx):
        view = CompensationDiagnostic()
        view.copy_traits(self, view.copyable_trait_names())

        return dedent("""
        op_{idx}.default_view({traits}).plot(ex_{prev_idx})
        """
        .format(traits = traits_str(view),
                idx = idx,
                prev_idx = idx - 1))


### Serialization
@camel_registry.dumper(CompensationWorkflowOp, 'compensation', version = 1)
def _dump(op):
    return dict(spillover_file = op.spillover_file,
                spillover = op.spillover)

@camel_registry.loader('compensation', version = 1)
def _load(data, version):
    return CompensationWorkflowOp(**data)

@camel_registry.dumper(CompensationWorkflowView, 'compensation-view', version = 1)
def _dump_view(view):
    return dict(op = view.op)

@camel_registry.loader('compensation-view', version = 1)
def _load_view(data, version):
    return CompensationWorkflowView(**data)
