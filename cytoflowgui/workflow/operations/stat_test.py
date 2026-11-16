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
cytoflowgui.workflow.operations.stat_test
-----------------------------------------

"""

from traits.api import (HasTraits, Str, Enum, Float, Bool, Property, List,
                        Instance, provides, observe)

from cytoflow.operations.stat_test import StatisticalTestOp, StatisticalTestView

from .. import Changed
from ..views import IWorkflowView, WorkflowView
from ..serialization import camel_registry, cytoflow_class_repr, traits_str, dedent
from ..subset import ISubset

from .operation_base import IWorkflowOperation, WorkflowOperation

StatisticalTestOp.__repr__ = cytoflow_class_repr
StatisticalTestView.__repr__ = cytoflow_class_repr


@provides(IWorkflowOperation)
class StatisticalTestWorkflowOp(WorkflowOperation, StatisticalTestOp):
    name = Str(apply = True)
    observations = Enum("statistic", "events", apply = True)
    statistic = Str(apply = True)
    feature = Str(apply = True)
    channel = Str(apply = True)
    by = Str(apply = True)
    test = Enum("anova", "kruskal", "welch", "mannwhitney", apply = True)
    posthoc = Enum("none", "tukey", apply = True)
    alpha = Float(0.05, apply = True)

    # override the base class's "subset" with one that is dynamically generated /
    # updated from subset_list
    subset = Property(Str, observe = "subset_list.items.str")
    subset_list = List(ISubset, apply = True)

    # bits to support the subset editor
    @observe('subset_list:items.str')
    def _on_subset_changed(self, _):
        self.changed = 'subset_list'

    # MAGIC - returns the value of the "subset" Property, above
    def _get_subset(self):
        return " and ".join([subset.str for subset in self.subset_list if subset.str])

    def default_view(self, **kwargs):
        return StatisticalTestWorkflowView(op = self, **kwargs)

    def clear_estimate(self):
        # no-op
        return

    def get_notebook_code(self, idx):
        op = StatisticalTestOp()
        op.copy_traits(self, op.copyable_trait_names())

        return dedent("""
        op_{idx} = {repr}

        ex_{idx} = op_{idx}.apply(ex_{prev_idx})
        """
        .format(repr = repr(op),
                idx = idx,
                prev_idx = idx - 1))


@provides(IWorkflowView)
class StatisticalTestWorkflowView(WorkflowView, StatisticalTestView):
    op = Instance(IWorkflowOperation, fixed = True)
    show_ns = Bool(False)
    plot_params = Instance(HasTraits, ())

    def should_plot(self, changed, payload):
        # this view plots the operation's result
        return changed in (Changed.VIEW, Changed.RESULT)

    def get_notebook_code(self, idx):
        view = StatisticalTestView()
        view.copy_traits(self, view.copyable_trait_names())

        return dedent("""
        op_{idx}.default_view({traits}).plot(ex_{idx})
        """
        .format(traits = traits_str(view),
                idx = idx))


### Serialization
@camel_registry.dumper(StatisticalTestWorkflowOp, 'statistical-test', version = 1)
def _dump(op):
    return dict(name = op.name,
                observations = op.observations,
                statistic = op.statistic,
                feature = op.feature,
                channel = op.channel,
                by = op.by,
                test = op.test,
                posthoc = op.posthoc,
                alpha = op.alpha,
                subset_list = op.subset_list)

@camel_registry.loader('statistical-test', version = 1)
def _load(data, version):
    return StatisticalTestWorkflowOp(**data)

@camel_registry.dumper(StatisticalTestWorkflowView, 'statistical-test-view', version = 1)
def _dump_view(view):
    return dict(op = view.op,
                show_ns = view.show_ns)

@camel_registry.loader('statistical-test-view', version = 1)
def _load_view(data, version):
    return StatisticalTestWorkflowView(**data)
