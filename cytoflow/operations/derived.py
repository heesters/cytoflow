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
cytoflow.operations.derived
---------------------------

Creates a new "channel" from an arithmetic expression of other channels
(like FlowJo's "Derive Parameters" or FCS Express's "Parameter Math").
`derived` has one class:

`DerivedParameterOp` -- applies the operation.
"""

import ast
import re
import numbers
from warnings import warn

import numpy as np
import pandas as pd

from traits.api import HasStrictTraits, Str, Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation

_BINARY_OPS = {ast.Add : np.add,
               ast.Sub : np.subtract,
               ast.Mult : np.multiply,
               ast.Div : np.true_divide,
               ast.Pow : np.power}

_UNARY_OPS = {ast.USub : np.negative,
              ast.UAdd : np.positive}

# name -> (function, number of arguments)
FUNCTIONS = {'log' : (np.log, 1),
             'ln' : (np.log, 1),
             'log10' : (np.log10, 1),
             'log2' : (np.log2, 1),
             'exp' : (np.exp, 1),
             'sqrt' : (np.sqrt, 1),
             'abs' : (np.abs, 1),
             'arcsinh' : (np.arcsinh, 1),
             'asinh' : (np.arcsinh, 1),
             'min' : (np.minimum, 2),
             'max' : (np.maximum, 2)}
"""The functions that can be used in a `DerivedParameterOp` expression."""

_BACKTICK = re.compile(r'`([^`]*)`')
_PLACEHOLDER = "__cytoflow_channel_{}__"


@provides(IOperation)
class DerivedParameterOp(HasStrictTraits):
    """
    Create a new "channel" by computing an arithmetic expression of other
    channels, for example ``"FITC_A / PE_A"``, ``"log10(B1_A + 1)"`` or
    ``"(V2_A - V1_A) * 2"``.

    Channels are referred to by name. Because a name like ``FITC-A`` isn't
    a valid identifier (the ``-`` would be a minus sign), you can write the
    name with every character that isn't a letter, number or underscore
    replaced by ``_`` (eg ``FITC_A``), or quote the name in backticks
    (eg ```FITC-A```).

    The expression may contain:

    - channel names
    - numbers
    - ``+``, ``-``, ``*``, ``/``, ``**`` (power), unary ``-`` and parentheses
    - these functions: ``log`` (natural log; also ``ln``), ``log10``,
      ``log2``, ``exp``, ``sqrt``, ``abs``, ``arcsinh`` (also ``asinh``),
      ``min(a, b)`` and ``max(a, b)`` (element-wise).

    Nothing else is allowed. The expression is parsed and checked, never
    passed to Python's `eval`.

    Events for which the expression isn't a finite number -- eg, division
    by zero, or the log of a negative number -- are **dropped** from the
    experiment (the same as `RatioOp`), and a warning says how many were
    dropped.

    Attributes
    ----------
    name : Str
        The operation name.  Also becomes the name of the new channel.

    expression : Str
        The arithmetic expression that computes the new channel.

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

    Compute a new channel.

    .. plot::
        :context: close-figs

        >>> derived_op = flow.DerivedParameterOp(name = "Y2_over_B1",
        ...                                      expression = "Y2_A / B1_A")
        >>> ex2 = derived_op.apply(ex)

    Plot it.

    .. plot::
        :context: close-figs

        >>> derived_op.default_view(scale = "log").plot(ex2)
    """

    # traits
    id = Constant('cytoflow.operations.derived')
    friendly_id = Constant("Derived Parameter")

    name = Str
    expression = Str

    def apply(self, experiment):
        """
        Applies the operation to an experiment.

        Parameters
        ----------
        experiment : `Experiment`
            the `Experiment` to which this operation is applied

        Returns
        -------
        Experiment
            a new experiment with the new channel. Events where the new
            channel isn't a finite number are dropped.

            The new channel also has the following new metadata:

            - **expression** : Str
                The expression that computed the new channel.

            - **range** : Float
                The largest absolute (finite) value of the new channel.

        Raises
        ------
        CytoflowOpError
            if the name is missing or already used, or if the expression
            is empty, can't be parsed, refers to a channel that isn't in
            the experiment, or uses something that isn't allowed.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if not self.name:
            raise util.CytoflowOpError('name',
                                       "Must specify a name for the derived channel")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores.")

        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "New channel {0} is already in the experiment"
                                       .format(self.name))

        values = evaluate_channel_expression(self.expression, experiment)

        nonfinite = ~np.isfinite(values)
        n_bad = int(nonfinite.sum())

        if n_bad == len(values):
            raise util.CytoflowOpError('expression',
                                       "'{}' isn't a finite number for any "
                                       "event (eg division by zero, or the "
                                       "log of a negative number)"
                                       .format(self.expression))

        new_experiment = experiment.clone(deep = False)
        new_experiment.add_channel(self.name, values)

        if n_bad > 0:
            warn("'{}' was not a finite number (eg division by zero, or the "
                 "log of a negative number) for {} of {} events ({:.2%}); "
                 "those events were dropped."
                 .format(self.expression, n_bad, len(values), n_bad / len(values)),
                 util.CytoflowOpWarning)
            new_experiment.data = \
                new_experiment.data[~nonfinite.to_numpy()].reset_index(drop = True)

        finite = new_experiment.data[self.name]
        new_experiment.metadata[self.name]['range'] = \
            float(finite.abs().max()) if len(finite) > 0 else 1.0
        new_experiment.metadata[self.name]['expression'] = self.expression
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))

        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns a `HistogramView` of the new channel.

        Returns
        -------
            IView : a histogram of the new channel
        """
        from cytoflow.views import HistogramView
        kwargs = dict(kwargs)
        kwargs.setdefault('channel', self.name)
        return HistogramView(**kwargs)


def resolve_channel_name(name, channels):
    """
    Find the channel that ``name`` refers to: either an exact match, or a
    channel whose name is ``name`` once non-identifier characters are replaced
    by ``_`` (so ``FITC_A`` refers to ``FITC-A``). Returns ``None`` if there's
    no match.

    Raises
    ------
    CytoflowOpError
        if ``name`` matches more than one channel.
    """
    if name in channels:
        return name

    matches = [c for c in channels if util.sanitize_identifier(c) == name]
    if len(matches) > 1:
        raise util.CytoflowOpError('expression',
                                   "'{}' could mean any of channels {}; quote "
                                   "the one you want in backticks, eg `{}`"
                                   .format(name, ", ".join(matches), matches[0]))
    return matches[0] if matches else None


def evaluate_channel_expression(expression, experiment):
    """
    Safely evaluates a `DerivedParameterOp` expression against an `Experiment`.

    The expression is parsed with `ast`; only a small set of node types is
    allowed, and the tree is walked by hand (nothing is passed to `eval`).
    Floating-point errors (division by zero, log of a negative number,
    overflow) don't raise; they produce ``inf`` or ``NaN``.

    Returns
    -------
    pandas.Series
        A float series, indexed like ``experiment.data``.

    Raises
    ------
    CytoflowOpError
        if the expression is empty, isn't valid, refers to an unknown
        channel or uses a forbidden construct.
    """

    if not expression or not expression.strip():
        raise util.CytoflowOpError('expression',
                                   "Must specify an expression, such as "
                                   "'FITC_A / PE_A'")

    # replace `quoted names` with placeholders
    quoted = {}
    def _quote(m):
        key = _PLACEHOLDER.format(len(quoted))
        quoted[key] = m.group(1)
        return key

    source = _BACKTICK.sub(_quote, expression.strip())

    if '`' in source:
        raise util.CytoflowOpError('expression',
                                   "Unmatched backtick in '{}'".format(expression))

    try:
        tree = ast.parse(source, mode = 'eval')
    except SyntaxError as e:
        raise util.CytoflowOpError('expression',
                                   "Couldn't parse expression '{}': {}"
                                   .format(expression, e.msg)) from e
    except (RecursionError, MemoryError, ValueError) as e:
        raise util.CytoflowOpError('expression',
                                   "Couldn't parse expression '{}'"
                                   .format(expression)) from e

    evaluator = _ArithmeticEvaluator(experiment, quoted)
    try:
        with np.errstate(all = 'ignore'):
            result = evaluator.evaluate(tree.body)
    except (RecursionError, OverflowError) as e:
        raise util.CytoflowOpError('expression',
                                   "Couldn't evaluate expression '{}': it is "
                                   "too deeply nested or a number is too large"
                                   .format(expression)) from e

    if not evaluator.used_channels:
        raise util.CytoflowOpError('expression',
                                   "The expression '{}' doesn't use any channels"
                                   .format(expression))

    return pd.Series(np.asarray(result, dtype = 'float64'),
                     index = experiment.data.index)


class _ArithmeticEvaluator(object):

    def __init__(self, experiment, quoted):
        self.experiment = experiment
        self.channels = experiment.channels
        self.quoted = quoted
        self.used_channels = []

    def _error(self, msg):
        return util.CytoflowOpError('expression', msg)

    def _source(self, node):
        try:
            s = ast.unparse(node)
        except Exception:
            return type(node).__name__
        for k, v in self.quoted.items():
            s = s.replace(k, '`' + v + '`')
        return s

    def evaluate(self, node):
        if isinstance(node, ast.BinOp) and type(node.op) in _BINARY_OPS:
            return _BINARY_OPS[type(node.op)](self.evaluate(node.left),
                                              self.evaluate(node.right))

        elif isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY_OPS:
            return _UNARY_OPS[type(node.op)](self.evaluate(node.operand))

        elif isinstance(node, ast.Constant):
            if isinstance(node.value, bool) or not isinstance(node.value, numbers.Real):
                raise self._error("'{}' isn't allowed: only numbers can be "
                                  "constants".format(self._source(node)))
            return float(node.value)

        elif isinstance(node, ast.Name):
            if not isinstance(node.ctx, ast.Load):
                raise self._error("Can't assign to '{}'".format(node.id))

            if node.id in self.quoted:
                name = self.quoted[node.id]
                channel = name if name in self.channels else None
            else:
                name = node.id
                channel = resolve_channel_name(name, self.channels)

            if channel is None:
                if name in FUNCTIONS:
                    raise self._error("'{}' is a function; call it like {}(...)"
                                      .format(name, name))
                raise self._error("'{}' isn't a channel in the experiment. "
                                  "Available channels: {}"
                                  .format(name, ", ".join(self.channels)))

            if channel not in self.used_channels:
                self.used_channels.append(channel)
            return self.experiment.data[channel].to_numpy(dtype = 'float64')

        elif isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in FUNCTIONS:
                raise self._error("'{}' isn't an allowed function. Allowed "
                                  "functions are: {}"
                                  .format(self._source(node.func),
                                          ", ".join(FUNCTIONS.keys())))

            fn, nargs = FUNCTIONS[node.func.id]
            if node.keywords or len(node.args) != nargs \
                or any(isinstance(a, ast.Starred) for a in node.args):
                raise self._error("{}() takes exactly {} argument{}"
                                  .format(node.func.id, nargs,
                                          "" if nargs == 1 else "s"))

            return fn(*[self.evaluate(a) for a in node.args])

        else:
            raise self._error("'{}' isn't allowed in a derived parameter "
                              "expression. Use channel names, numbers, "
                              "+ - * / **, parentheses and the functions {}."
                              .format(self._source(node),
                                      ", ".join(FUNCTIONS.keys())))
