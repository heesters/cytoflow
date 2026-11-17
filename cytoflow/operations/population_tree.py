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
cytoflow.operations.population_tree
-----------------------------------

A "population tree": nested populations (gates applied inside their parent
gates), with the number of events in each, the percentage of the parent and
the percentage of all events -- like the workspace tree in FlowJo or the
gate navigator in FCS Express.  `population_tree` has:

`PopulationTreeOp` -- the `IOperation` that computes the tree's statistics.

`infer_populations` -- builds the list of populations from a list of 
operations (eg, the gates imported from a FlowJo workspace.)

`PopulationTreeView` -- draws the tree.
"""

import ast

from traits.api import (HasStrictTraits, Str, List, Tuple, Constant, Instance,
                        provides)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cytoflow.views
import cytoflow.utility as util

from .i_operation import IOperation
from .boolean_gate import evaluate_boolean_expression


def _order_and_check(populations):
    """
    Validates a list of ``(name, gate, parent)`` and returns it ordered so
    that parents come before their children (otherwise keeping the order),
    with each population's depth: a list of ``(name, gate, parent, depth)``.
    """
    names = [p[0] for p in populations]
    for name, gate, parent in populations:
        if not name:
            raise util.CytoflowOpError('populations', "Every population needs a name")
        if names.count(name) > 1:
            raise util.CytoflowOpError('populations',
                                       "Population '{}' is defined more than once".format(name))
        if not gate:
            raise util.CytoflowOpError('populations',
                                       "Population '{}' needs a gate".format(name))
        if parent and parent not in names:
            raise util.CytoflowOpError('populations',
                                       "The parent of '{}', '{}', isn't a population"
                                       .format(name, parent))
        if parent == name:
            raise util.CytoflowOpError('populations',
                                       "Population '{}' can't be its own parent".format(name))
            
    parents = {p[0] : p[2] for p in populations}
    depth = {}
    for name in names:
        seen = [name]
        p = parents[name]
        while p:
            if p in seen:
                raise util.CytoflowOpError('populations',
                                           "The populations' parents form a loop: {}"
                                           .format(" -> ".join(seen + [p])))
            seen.append(p)
            p = parents[p]
        depth[name] = len(seen) - 1
        
    # depth-first, keeping the given order among siblings
    children = {}
    for name, _, parent in populations:
        children.setdefault(parent, []).append(name)
    gates = {p[0] : p[1] for p in populations}
    
    ordered = []
    def visit(parent):
        for name in children.get(parent, []):
            ordered.append((name, gates[name], parents[name], depth[name]))
            visit(name)
    visit("")
    return ordered


@provides(IOperation)
class PopulationTreeOp(HasStrictTraits):
    """
    Compute a *population tree*: for a set of nested populations, the number
    of events in each one, as a percentage of its parent population and of 
    all the events.
    
    In `cytoflow`, a gate applies to all the events: it makes a new condition
    that says whether each event is inside the gate.  A **population** is a 
    gate applied *inside its parent population* -- eg, "CD4+" inside 
    "T cells" inside "Singlets".  This operation doesn't change the data;
    it adds a statistic with the numbers, which you can view as a table or 
    export, and its default view draws the tree.
    
    Attributes
    ----------
    name : Str
        The operation name; becomes the name of the new statistic.
        
    populations : List((Str, Str, Str))
        The populations, as ``(name, gate, parent)``:
        
        - ``name``: what to call the population.
        - ``gate``: which events are in the gate.  Usually the name of a 
          True/False condition made by a gate (``"Singlets"``), but it can 
          be any expression that `BooleanGateOp` understands: 
          ``"Quad == 'Quad_2'"``, ``"CD4_pos and not CD8_pos"``.
        - ``parent``: the name of the parent population, or ``""`` if the 
          population is at the top of the tree (its parent is all events.)
          
        An event is in a population if it is in the population's gate and 
        in the gates of all its ancestors.
        
    by : List(Str)
        Compute the tree separately for each combination of these 
        conditions' values -- eg, for each sample or each treatment.
        
    Notes
    -----
    The statistic is indexed by the `by` conditions and ``Population``,
    with the populations in tree order, and has the columns:
    
    - ``Count``: the number of events in the population.
    - ``Percent_of_parent``: ``Count`` as a percentage of the parent 
      population's count (not-a-number if the parent is empty.)
    - ``Percent_of_total``: ``Count`` as a percentage of all the events (in
      that group.)
    - ``Depth``: 0 for populations at the top of the tree, 1 for their
      children, and so on.
      
    Examples
    --------
    
    >>> ex2 = flow.ThresholdOp(name = "Big", channel = "FSC-A", threshold = 50000).apply(ex)
    >>> ex3 = flow.ThresholdOp(name = "Bright", channel = "Y2-A", threshold = 1000).apply(ex2)
    >>> tree_op = flow.PopulationTreeOp(name = "Tree", by = ["Dox"],
    ...                                 populations = [("Cells", "Big", ""),
    ...                                                ("Positive", "Bright", "Cells"),
    ...                                                ("Negative", "not Bright", "Cells")])
    >>> ex4 = tree_op.apply(ex3)
    >>> ex4.statistics["Tree"]
    >>> tree_op.default_view().plot(ex4)
    """
    
    id = Constant('cytoflow.operations.population_tree')
    friendly_id = Constant("Population Tree")
    
    name = Str
    populations = List(Tuple(Str, Str, Str))
    by = List(Str)
    
    def compute(self, experiment, by = None):
        """
        Compute the tree's numbers.
        
        Parameters
        ----------
        experiment : Experiment
        
        by : List(Str)
            Overrides `by`; pass ``[]`` to pool all the events.
            
        Returns
        -------
        pandas.DataFrame
            The statistic described in the class documentation.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")
        
        if not self.populations:
            raise util.CytoflowOpError('populations', "Must specify at least one population")
        
        by = list(self.by if by is None else by)
        for b in by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "Condition '{}' isn't in the experiment".format(b))
                
        ordered = _order_and_check(list(self.populations))
        
        # which events are in each population?
        member = {}
        n = len(experiment)
        for name, gate, parent, _ in ordered:
            try:
                mask = np.asarray(evaluate_boolean_expression(gate, experiment), dtype = bool)
            except util.CytoflowOpError as e:
                raise util.CytoflowOpError('populations',
                                           "Population '{}': {}".format(name, e.args[-1])) from e
            if mask.shape != (n,):
                raise util.CytoflowOpError('populations',
                                           "Population '{}': the gate must be True "
                                           "or False for each event".format(name))
            member[name] = mask & member[parent] if parent else mask
            
        names = [o[0] for o in ordered]
        
        if by:
            groups = experiment.data.groupby(by, observed = True).indices
            keys = sorted(groups.keys())
        else:
            groups = {() : np.arange(n)}
            keys = [()]
            
        rows = []
        index = []
        for key in keys:
            idx = groups[key]
            total = len(idx)
            counts = {name : int(member[name][idx].sum()) for name in names}
            for name, _, parent, depth in ordered:
                parent_count = counts[parent] if parent else total
                rows.append((counts[name],
                             100.0 * counts[name] / parent_count if parent_count else np.nan,
                             100.0 * counts[name] / total if total else np.nan,
                             depth))
                k = key if isinstance(key, tuple) else (key,)
                index.append(tuple(k) + (name,))
                
        stat = pd.DataFrame(rows, 
                            columns = ["Count", "Percent_of_parent", "Percent_of_total", "Depth"],
                            index = pd.MultiIndex.from_tuples(index, names = by + ["Population"]))
        
        if not by:
            # a plain index is easier to use than a one-level MultiIndex
            stat.index = pd.Index([i[0] for i in index], name = "Population")
        
        return stat
    
    def apply(self, experiment):
        """
        Add the tree's statistic to the experiment.
        
        Returns
        -------
        Experiment
            A new experiment (sharing ``experiment``'s data) with the 
            statistic `name`.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")
        
        if not self.name:
            raise util.CytoflowOpError('name', "Must specify a name")
        
        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores.")
            
        if self.name in experiment.statistics:
            raise util.CytoflowOpError('name',
                                       "Experiment already has a statistic named {}"
                                       .format(self.name))
            
        if "Population" in self.by:
            raise util.CytoflowOpError('by', "Can't group by a condition named 'Population'")
            
        stat = self.compute(experiment)
        
        new_experiment = experiment.clone(deep = False)
        new_experiment.statistics[self.name] = stat
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment
    
    def default_view(self, **kwargs):
        """
        Returns a `PopulationTreeView` that draws the tree.
        """
        v = PopulationTreeView(op = self)
        v.trait_set(**kwargs)
        return v
    
    
def infer_populations(ops, parents = None):
    """
    Build a list of populations from a list of operations.
    
    Gates imported from FlowJo or Gating-ML (`import_flowjo_workspace`, 
    `import_gatingml`) keep their hierarchy as `BooleanGateOp` s: the gate 
    ``X`` is followed by ``X_gated`` with the expression ``X and P`` (or 
    ``X and P_gated``.)  This function reads the tree back from them.  
    
    Parameters
    ----------
    ops : list of `IOperation`
        The operations, in order.
        
    parents : dict
        Optional: ``{gate name : parent gate name}`` for gates whose parent 
        isn't recorded in a `BooleanGateOp` (the GUI works these out from 
        the subset that each gate was drawn on.)
        
    Returns
    -------
    list of (name, gate, parent)
        Suitable for `PopulationTreeOp.populations`.  Gates with a True / 
        False condition become populations; other operations are ignored.
    """
    parents = dict(parents or {})
    
    # from "X_gated = X and P[_gated]"
    gated_of = {}        # "X_gated" -> "X"
    hidden = set()
    for op in ops:
        if getattr(op, 'id', None) != 'cytoflow.operations.boolean_gate':
            continue
        try:
            tree = ast.parse(op.expression.strip(), mode = 'eval').body
        except SyntaxError:
            continue
        if (isinstance(tree, ast.BoolOp) and isinstance(tree.op, ast.And)
                and len(tree.values) == 2
                and all(isinstance(v, ast.Name) for v in tree.values)
                and op.name.endswith("_gated")):
            child, parent = tree.values[0].id, tree.values[1].id
            gated_of[op.name] = child
            parents.setdefault(child, parent)
            hidden.add(op.name)
            
    # parents may be given as "P_gated"; point them at "P"
    def resolve(p):
        seen = set()
        while p in gated_of and p not in seen:
            seen.add(p)
            p = gated_of[p]
        return p
    
    populations = []
    for op in ops:
        op_id = getattr(op, 'id', '')
        name = getattr(op, 'name', '')
        if not name or name in hidden:
            continue
        if op_id in BOOLEAN_GATE_IDS:
            populations.append((name, name, resolve(parents.get(name, ""))))
            
    names = [p[0] for p in populations]
    return [(n, g, p if p in names and p != n else "") for n, g, p in populations]


# operations that add one True / False condition named after the operation
BOOLEAN_GATE_IDS = ('cytoflow.operations.threshold',
                    'cytoflow.operations.range',
                    'cytoflow.operations.range2d',
                    'cytoflow.operations.polygon',
                    'cytoflow.operations.freehand',
                    'cytoflow.operations.ellipse',
                    'cytoflow.operations.boolean_gate')


@provides(cytoflow.views.IView)
class PopulationTreeView(HasStrictTraits):
    """
    Draws a population tree: one row for each population, indented under 
    its parent, with a bar showing its share of all the events and its 
    count and percentage of its parent.
    
    Attributes
    ----------
    op : Instance(PopulationTreeOp)
        The operation whose tree to draw.
        
    subset : Str
        Only count these events (eg, ``"Dox == 10.0"`` to see one sample.)
        By default, all the events are pooled.
    """
    
    id = Constant('cytoflow.view.populationtreeview')
    friendly_id = Constant("Population Tree")
    
    op = Instance(PopulationTreeOp)
    subset = Str
    
    def plot(self, experiment, **kwargs):
        """
        Draw the tree.
        
        Parameters
        ----------
        title : str
            The plot title
            
        color : a matplotlib color
            The color of the bars.  Default: ``"C0"``
        """
        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")
        
        if self.op is None:
            raise util.CytoflowViewError('op', "No operation specified")
        
        if self.subset:
            try:
                experiment = experiment.query(self.subset)
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' isn't valid"
                                             .format(self.subset)) from e
            if len(experiment) == 0:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' returned no events"
                                             .format(self.subset))
                
        try:
            stat = self.op.compute(experiment, by = [])
        except util.CytoflowOpError as e:
            raise util.CytoflowViewError(e.args[0], e.args[-1]) from e
        
        names = list(stat.index.get_level_values("Population"))
        n = len(names)
        total = len(experiment)
        color = kwargs.get('color', 'C0')
        
        fig = plt.figure(figsize = (7, max(1.5, 0.42 * (n + 1) + 0.6)))
        ax = fig.add_subplot(111)
        
        max_depth = int(stat["Depth"].max())
        indent = 0.035
        label_x = 0.02
        bar_x0 = min(0.62, 0.36 + indent * max_depth)
        bar_w = 0.97 - bar_x0
        
        rows = [("All events", total, np.nan, 100.0, -1)] + \
               [(name, int(r.Count), r.Percent_of_parent, r.Percent_of_total, int(r.Depth))
                for name, r in zip(names, stat.itertuples())]
               
        for i, (name, count, pct_parent, pct_total, depth) in enumerate(rows):
            y = len(rows) - 1 - i
            x = label_x + indent * (depth + 1)
            ax.text(x, y, name, va = 'center', ha = 'left',
                    fontweight = 'bold' if depth < 0 else 'normal',
                    transform = ax.get_yaxis_transform())
            # the connector to the parent
            if depth >= 0:
                ax.plot([x - indent * 0.75, x - indent * 0.15], [y, y], color = '0.6', lw = 1,
                        transform = ax.get_yaxis_transform(), clip_on = False)
            frac = 0.0 if not np.isfinite(pct_total) else pct_total / 100.0
            ax.barh(y, bar_w, left = bar_x0, height = 0.62, color = '0.93', 
                    transform = ax.get_yaxis_transform(), zorder = 1)
            ax.barh(y, bar_w * frac, left = bar_x0, height = 0.62, 
                    color = '0.55' if depth < 0 else color,
                    transform = ax.get_yaxis_transform(), zorder = 2)
            if depth < 0:
                label = "{:,}".format(count)
            elif np.isfinite(pct_parent):
                label = "{:,}  ({:.1f}% of parent, {:.1f}% of all)".format(count, pct_parent, pct_total)
            else:
                label = "{:,}  (parent is empty)".format(count)
            ax.text(bar_x0 + 0.01, y, label, va = 'center', ha = 'left', fontsize = 'small',
                    transform = ax.get_yaxis_transform(), zorder = 3)
            
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xlim(0, 1)
        ax.set_axis_off()
        if kwargs.get('title'):
            ax.set_title(kwargs['title'])
            

util.expand_class_attributes(PopulationTreeOp)
util.expand_method_parameters(PopulationTreeOp, PopulationTreeOp.apply)
