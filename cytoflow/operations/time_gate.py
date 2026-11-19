#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflow.operations.time_gate
-----------------------------

Gate on the Time channel to exclude debris acquired at the start of a run
or events acquired after a clog.

`TimeGateOp` -- extends RangeOp with auto-detection of the time channel.
"""

from traits.api import Constant, provides

import cytoflow.utility as util
from .i_operation import IOperation
from .range import RangeOp, RangeSelection


# FCS files use various names for the time channel
_TIME_CHANNEL_NAMES = ('Time', 'TIME', 'time', 'TIME_MS', 'time_ms',
                       'TimeMSec', 'TIMEMS')


def find_time_channel(experiment):
    """Return the time channel name from an experiment, or None."""
    for name in _TIME_CHANNEL_NAMES:
        if name in experiment.channels:
            return name
    for ch in experiment.channels:
        if 'time' in ch.lower():
            return ch
    return None


@provides(IOperation)
class TimeGateOp(RangeOp):
    """
    Gate on the instrument Time channel to exclude unstable acquisition
    periods (debris at the start of a run, clogs, pressure fluctuations).

    Extends `RangeOp` and auto-detects the time channel if `channel` is not
    set.  The gate keeps events where ``low <= Time <= high``.

    Attributes
    ----------
    name : Str
        Name for the boolean gate condition added to the experiment.
    channel : Str
        The time channel.  Auto-detected from the experiment if not set.
    low : Float
        Lower bound (in the channel's native units, usually ms or seconds).
    high : Float
        Upper bound.

    Examples
    --------
    >>> tg = flow.TimeGateOp(name='TimeGate', low=5000, high=55000)
    >>> ex2 = tg.apply(ex)
    """

    id = Constant('cytoflow.operations.time_gate')
    friendly_id = Constant('Time Gate')

    def apply(self, experiment):
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        # Auto-detect time channel if not explicitly set
        if not self.channel:
            ch = find_time_channel(experiment)
            if ch is None:
                raise util.CytoflowOpError(
                    'channel',
                    "Could not find a Time channel in the experiment. "
                    "Available channels: {}. "
                    "Set channel explicitly.".format(sorted(experiment.channels)))
            self.channel = ch

        return super().apply(experiment)

    def default_view(self, **kwargs):
        self._selection_view = TimeGateSelection(op=self)
        self._selection_view.trait_set(**kwargs)
        return self._selection_view


class TimeGateSelection(RangeSelection):
    """Range selection view for the time channel — always uses linear scale."""

    id = Constant('cytoflow.views.time_gate')
    friendly_id = Constant('Time Gate Selection')

    def plot(self, experiment, **kwargs):
        # Ensure linear scale for time
        kwargs.setdefault('scale', 'linear')
        # Auto-detect channel
        if not self.op.channel and experiment is not None:
            ch = find_time_channel(experiment)
            if ch:
                self.op.channel = ch
        super().plot(experiment, **kwargs)
