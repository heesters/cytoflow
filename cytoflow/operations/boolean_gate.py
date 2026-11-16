#!/usr/bin/env python3.11
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2025
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
cytoflow.operations.boolean_gate
--------------------------------

Combine existing gates with AND, OR and NOT into a new gate (like FlowJo's
"Boolean gates" or FCS Express's "Combination gates").  `boolean_gate` has
one class:

`BooleanGateOp` -- applies the operation.
"""

import ast
import numbers

import numpy as np
import pandas as pd

from traits.api import HasStrictTraits, Str, Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation

# comparison operators allowed between a condition and a constant
_COMPARISONS = {ast.Eq : lambda a, b: a == b,
                ast.NotEq : lambda a, b: a != b,
                ast.Lt : lambda a, b: a < b,
                ast.LtE : lambda a, b: a <= b,
                ast.Gt : lambda a, b: a > b,
                ast.GtE : lambda a, b: a >= b}

# ordering comparisons only make sense against numbers
_ORDERING = (ast.Lt, ast.LtE, ast.Gt, ast.GtE)


@provides(IOperation)
class BooleanGateOp(HasStrictTraits):
    """
    Combine existing gates (boolean conditions) into a new gate using
    ``and``, ``or`` and ``not``.

    The new gate is defined by a logical `expression`, for example
    ``"CD4_pos and not Dead"``. Every name in the expression must be an
    existing condition in the `Experiment` -- usually the name of another
    gate, such as a `ThresholdOp`, `RangeOp`, `Range2DOp` or `PolygonOp`.
    The expression may contain:

    - the names of **boolean** conditions (eg, the output of a threshold,
      range or polygon gate)
    - ``and``, ``or``, ``not`` (or their equivalents ``&``, ``|``, ``~``)
      and parentheses
    - comparisons of a condition to a constant, such as
      ``Quad == 'Quad_2'`` (useful for categorical gates like `QuadOp`),
      ``Y2_high == False`` or ``Dox > 1.0``. ``==`` and ``!=`` work with any
      condition; ``<``, ``<=``, ``>`` and ``>=`` need a number.

    Nothing else is allowed. The expression is parsed and checked, never
    passed to Python's `eval`.

    The result is a new boolean condition named `name`, which is ``True``
    for events that satisfy the expression.

    Attributes
    ----------
    name : Str
        The operation name.  Used to name the new boolean condition in the
        experiment that's created by `apply`.

    expression : Str
        The logical expression defining the new gate.

    Examples
    --------
    Make a little data set.

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/RFP_Well_A6.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Make two gates.

    .. plot::
        :context: close-figs

        >>> ex2 = flow.ThresholdOp(name = 'Y2_high',
        ...                        channel = 'Y2-A',
        ...                        threshold = 300).apply(ex)
        >>> ex3 = flow.RangeOp(name = 'B1_mid',
        ...                    channel = 'B1-A',
        ...                    low = 100,
        ...                    high = 1000).apply(ex2)

    Combine them into a new gate.

    .. plot::
        :context: close-figs

        >>> bool_op = flow.BooleanGateOp(name = "Y2_high_not_B1_mid",
        ...                              expression = "Y2_high and not B1_mid")
        >>> ex4 = bool_op.apply(ex3)

    Plot the new gate.

    .. plot::
        :context: close-figs

        >>> bool_op.default_view(xchannel = "B1-A",
        ...                      xscale = "log",
        ...                      ychannel = "Y2-A",
        ...                      yscale = "log").plot(ex4)
    """

    # traits
    id = Constant('cytoflow.operations.boolean_gate')
    friendly_id = Constant("Boolean Gate")

    name = Str
    expression = Str

    def apply(self, experiment):
        """
        Applies the boolean gate to an experiment.

        Parameters
        ----------
        experiment : `Experiment`
            the `Experiment` to which this operation is applied

        Returns
        -------
        Experiment
            a new `Experiment`, the same as ``experiment`` but with a new
            boolean condition named `name`.

            The new condition also has the following new metadata:

            - **expression** : Str
                The expression that defined it.

        Raises
        ------
        CytoflowOpError
            if the name is missing or already used, or if the expression
            is empty, can't be parsed, refers to a condition that isn't in
            the experiment, or uses something other than the allowed
            operators.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if not self.name:
            raise util.CytoflowOpError('name',
                                       "You have to set the gate's name "
                                       "before applying it!")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name {} can only contain letters, numbers and underscores."
                                       .format(self.name))

        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "Experiment already has a column named {0}"
                                       .format(self.name))

        gate = evaluate_boolean_expression(self.expression, experiment)

        new_experiment = experiment.clone(deep = False)
        new_experiment.add_condition(self.name, "bool", gate)
        new_experiment.metadata[self.name]['expression'] = self.expression
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))

        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns a view showing the events in the new gate.

        If ``ychannel`` is set, returns a `ScatterplotView`; otherwise, a
        `HistogramView`. In both cases, ``huefacet`` is set to the new gate,
        so events in the gate are a different color than events outside it.
        Set ``channel`` (or ``xchannel`` and ``ychannel``) and anything else
        with the keyword arguments.

        Returns
        -------
            IView : a view on the new gate
        """
        from cytoflow.views import HistogramView, ScatterplotView

        kwargs = dict(kwargs)
        kwargs.setdefault('huefacet', self.name)

        if 'ychannel' in kwargs:
            return ScatterplotView(**kwargs)

        if 'xchannel' in kwargs:
            kwargs['channel'] = kwargs.pop('xchannel')
        if 'xscale' in kwargs:
            kwargs['scale'] = kwargs.pop('xscale')

        return HistogramView(**kwargs)


def boolean_condition_names(experiment):
    """
    Returns the names of the boolean conditions in an `Experiment` -- the
    names that can be used on their own in a `BooleanGateOp` expression.
    """
    return [c for c in experiment.conditions
            if experiment.data[c].dtype.kind == 'b']


def evaluate_boolean_expression(expression, experiment):
    """
    Safely evaluates a `BooleanGateOp` expression against an `Experiment`.

    The expression is parsed with `ast`; only a small set of node types is
    allowed, and the tree is walked by hand (nothing is passed to `eval`).

    Returns
    -------
    pandas.Series
        A boolean series, indexed like ``experiment.data``.

    Raises
    ------
    CytoflowOpError
        if the expression is empty, isn't valid, refers to an unknown
        condition or uses a forbidden construct.
    """

    if not expression or not expression.strip():
        raise util.CytoflowOpError('expression',
                                   "Must specify an expression, such as "
                                   "'Gate1 and not Gate2'")

    try:
        tree = ast.parse(expression.strip(), mode = 'eval')
    except SyntaxError as e:
        raise util.CytoflowOpError('expression',
                                   "Couldn't parse expression '{}': {}"
                                   .format(expression, e.msg)) from e
    except (RecursionError, MemoryError, ValueError) as e:
        raise util.CytoflowOpError('expression',
                                   "Couldn't parse expression '{}'"
                                   .format(expression)) from e

    try:
        return _BooleanEvaluator(experiment).evaluate(tree.body)
    except RecursionError as e:
        raise util.CytoflowOpError('expression',
                                   "Expression '{}' is too deeply nested"
                                   .format(expression)) from e


class _BooleanEvaluator(object):

    def __init__(self, experiment):
        self.experiment = experiment

    def _error(self, msg):
        return util.CytoflowOpError('expression', msg)

    def _condition(self, node):
        if not isinstance(node.ctx, ast.Load):
            raise self._error("Can't assign to '{}'".format(node.id))

        if node.id not in self.experiment.conditions:
            bools = boolean_condition_names(self.experiment)
            raise self._error("'{}' isn't a condition (gate) in the experiment. "
                              "Available gates: {}"
                              .format(node.id, ", ".join(bools) if bools else "(none)"))
        return self.experiment.data[node.id]

    def evaluate(self, node):
        if isinstance(node, ast.BoolOp):
            values = [self.evaluate(v) for v in node.values]
            result = values[0]
            for v in values[1:]:
                if isinstance(node.op, ast.And):
                    result = result & v
                else:
                    result = result | v
            return result

        elif isinstance(node, ast.BinOp) and isinstance(node.op, (ast.BitAnd, ast.BitOr)):
            left = self.evaluate(node.left)
            right = self.evaluate(node.right)
            return (left & right) if isinstance(node.op, ast.BitAnd) else (left | right)

        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.Not, ast.Invert)):
            return ~self.evaluate(node.operand)

        elif isinstance(node, ast.Name):
            col = self._condition(node)
            if col.dtype.kind != 'b':
                raise self._error("'{0}' isn't a boolean gate, so it can't be "
                                  "used on its own. Compare it to a value "
                                  "instead, eg {0} == {1!r}"
                                  .format(node.id,
                                          self._example_value(node.id)))
            return col.astype(bool)

        elif isinstance(node, ast.Compare):
            return self._compare(node)

        elif isinstance(node, ast.Constant) and isinstance(node.value, bool):
            return pd.Series(node.value, index = self.experiment.data.index)

        else:
            raise self._error("'{}' isn't allowed in a boolean gate expression. "
                              "Use gate names, 'and', 'or', 'not', parentheses "
                              "and comparisons like Quad == 'Quad_1'."
                              .format(self._source(node)))

    def _compare(self, node):
        if len(node.ops) != 1 or len(node.comparators) != 1:
            raise self._error("Chained comparisons like '{}' aren't allowed"
                              .format(self._source(node)))

        op = node.ops[0]
        left = node.left
        right = node.comparators[0]

        if type(op) not in _COMPARISONS:
            raise self._error("Comparison '{}' isn't allowed: use ==, !=, "
                              "<, <=, > or >=".format(self._source(node)))

        # put the condition on the left
        if isinstance(left, ast.Constant) and isinstance(right, ast.Name):
            left, right = right, left
            flipped = {ast.Lt : ast.Gt, ast.Gt : ast.Lt,
                       ast.LtE : ast.GtE, ast.GtE : ast.LtE}
            op = flipped.get(type(op), type(op))()

        if not (isinstance(left, ast.Name) and isinstance(right, ast.Constant)):
            raise self._error("Comparisons must be between a condition and a "
                              "constant, eg Quad == 'Quad_1'. '{}' isn't."
                              .format(self._source(node)))

        col = self._condition(left)
        value = right.value

        if not isinstance(value, (str, numbers.Number)):
            raise self._error("Can't compare to '{}'".format(self._source(right)))

        if isinstance(op, _ORDERING):
            if isinstance(value, (str, bool)) or col.dtype.kind not in 'iuf':
                raise self._error("'{}': <, <=, > and >= only work on numeric "
                                  "conditions and numbers"
                                  .format(self._source(node)))
        else:
            # catch typos in categorical values, which otherwise would
            # silently select nothing
            if col.dtype.kind == 'b':
                if not isinstance(value, (bool, np.bool_)) and value not in (0, 1):
                    raise self._error("'{}' is a boolean gate; compare it to "
                                      "True or False".format(left.id))
            else:
                values = self._values(left.id)
                if value not in values:
                    raise self._error("Condition '{}' has no value {!r}. Its "
                                      "values are: {}"
                                      .format(left.id, value,
                                              ", ".join([repr(v) for v in values])))

        result = _COMPARISONS[type(op)](col, value)
        return pd.Series(np.asarray(result, dtype = bool), index = col.index)

    def _values(self, name):
        meta = self.experiment.metadata.get(name, {})
        if 'values' in meta:
            values = list(meta['values'])
        else:
            values = list(self.experiment.data[name].unique())
        return values

    def _example_value(self, name):
        values = self._values(name)
        return values[0] if values else 'value'

    @staticmethod
    def _source(node):
        try:
            return ast.unparse(node)
        except Exception:
            return type(node).__name__
