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
cytoflow.operations.sample_labels
----------------------------------

Assign human-readable labels to sample/condition values.

`SampleLabelsOp` -- the `IOperation` class.
"""

from traits.api import HasStrictTraits, Str, Dict, Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation


@provides(IOperation)
class SampleLabelsOp(HasStrictTraits):
    """
    Assign human-readable labels to values of an existing condition (such as
    ``File``).  This is useful when file names are long or uninformative: map
    each filename to a short, descriptive label.

    If ``new_condition`` is set the relabelled values are stored in a **new**
    condition with that name; the original condition is left unchanged.
    If ``new_condition`` is empty (the default) the existing condition is
    replaced in-place.

    Any value that is **not** in ``labels`` keeps its original value unchanged.

    Attributes
    ----------
    condition : Str
        The condition to relabel (e.g. ``"File"``).

    labels : Dict(Str, Str)
        Mapping of original value → new label, e.g.
        ``{"ctrl.fcs": "Control", "treated.fcs": "Day 3"}``.

    new_condition : Str
        Name of the new condition to create. Leave empty to replace the
        original condition in-place.

    Examples
    --------

    >>> op = flow.SampleLabelsOp(
    ...     condition = "File",
    ...     labels = {"ctrl.fcs": "Control", "treated.fcs": "Treated"},
    ...     new_condition = "Group")
    >>> ex2 = op.apply(ex)
    """

    id = Constant('cytoflow.operations.sample_labels')
    friendly_id = Constant("Sample Labels")
    name = Constant("Sample Labels")

    condition = Str
    labels = Dict(Str, Str)
    new_condition = Str("")

    def apply(self, experiment):
        """
        Apply the label mapping to ``experiment``.

        Parameters
        ----------
        experiment : Experiment
            The old `Experiment` to which this operation is applied.

        Returns
        -------
        Experiment
            A new `Experiment` with the relabelled condition.

        Raises
        ------
        CytoflowOpError
            If the condition is missing or the new condition name is invalid.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.condition:
            raise util.CytoflowOpError('condition',
                                       "Must specify a condition to relabel")

        if self.condition not in experiment.conditions:
            raise util.CytoflowOpError('condition',
                                       "Condition '{}' is not in the experiment"
                                       .format(self.condition))

        if not self.labels:
            raise util.CytoflowOpError('labels',
                                       "Must specify at least one label mapping")

        dest = self.new_condition if self.new_condition else self.condition

        if self.new_condition:
            if dest in experiment.data.columns and dest != self.condition:
                raise util.CytoflowOpError('new_condition',
                                           "Condition '{}' already exists"
                                           .format(dest))
            if dest != util.sanitize_identifier(dest):
                raise util.CytoflowOpError('new_condition',
                                           "Name '{}' can only contain letters, "
                                           "numbers and underscores.".format(dest))

        new_experiment = experiment.clone(deep=False)

        col = new_experiment.data[self.condition].astype(str)
        new_col = col.map(lambda v: self.labels.get(v, v))

        # always store as categorical strings
        from pandas.api.types import CategoricalDtype
        new_cats = sorted(set(new_col.unique()))
        new_col = new_col.astype(CategoricalDtype(categories=new_cats, ordered=False))

        if dest == self.condition:
            new_experiment.data[dest] = new_col
        else:
            # Add as a brand-new condition (also registers metadata)
            new_experiment.add_condition(dest, 'category', new_col)

        new_experiment.history.append(
            self.clone_traits(transient=lambda _: True))
        return new_experiment

    def default_view(self, **kwargs):
        from cytoflow.views import HistogramView
        return HistogramView(**kwargs)


util.expand_class_attributes(SampleLabelsOp)
util.expand_method_parameters(SampleLabelsOp, SampleLabelsOp.apply)
