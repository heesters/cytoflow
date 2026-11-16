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
cytoflow.operations.scale_settings
----------------------------------

Set the parameters of the axis transforms (logicle, arcsinh, log, hlog)
for individual channels.  `scale_settings` has one class:

`ScaleSettingsOp` -- the `IOperation` that stores the settings.
"""

from traits.api import (HasStrictTraits, Str, Dict, Float, Constant, provides)

import cytoflow.utility as util
from .i_operation import IOperation

# the parameters that can be set for each scale, and a description
SCALE_PARAMETERS = {
    'logicle' : {'W' : "width of the linear region, in decades",
                 'M' : "width of the whole display, in decades",
                 'A' : "additional decades of negative values",
                 'r' : "quantile of the negative values used to estimate W"},
    'arcsinh' : {'cofactor' : "values much smaller than this are linear, much "
                              "larger are logarithmic"},
    'hlog' : {'b' : "location of the transition from linear to log"},
    'log' : {'threshold' : "smallest value to show"}
}


@provides(IOperation)
class ScaleSettingsOp(HasStrictTraits):
    """
    Set the parameters of the axis transforms ("scales") for individual
    channels -- eg, the logicle ``W`` (width of the linear region) or the
    arcsinh ``cofactor``.  Every plot and gate downstream of this operation
    that uses that scale on that channel uses these parameters.  (This is
    the equivalent of "customize axis" / "transform settings" in FlowJo or
    FCS Express.)

    The data aren't changed: the settings are stored in the experiment's
    channel metadata, under ``scale_params``.

    Attributes
    ----------
    settings : Dict(Str, Dict(Str, Dict(Str, Float)))
        ``{channel : {scale : {parameter : value}}}``.  The parameters for
        each scale are:

        - ``logicle``: ``W``, ``M``, ``A``, ``r``
        - ``arcsinh``: ``cofactor``
        - ``hlog``: ``b``
        - ``log``: ``threshold``

    Examples
    --------

    >>> op = flow.ScaleSettingsOp(settings = {"B1-A" : {"logicle" : {"W" : 1.0}},
    ...                                       "Y2-A" : {"arcsinh" : {"cofactor" : 500}}})
    >>> ex2 = op.apply(ex)
    >>> flow.HistogramView(channel = "Y2-A", scale = "arcsinh").plot(ex2)
    """

    id = Constant('cytoflow.operations.scale_settings')
    friendly_id = Constant("Axis Transforms")
    name = Constant("Axis Transforms")

    settings = Dict(Str, Dict(Str, Dict(Str, Float)))

    def apply(self, experiment):
        """
        Store the scale settings in a copy of ``experiment``'s metadata.

        Returns
        -------
        Experiment
            A new experiment (sharing ``experiment``'s data), with each
            configured channel's ``scale_params`` metadata updated.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.settings:
            raise util.CytoflowOpError('settings', "No settings specified")

        for channel, scales in self.settings.items():
            if channel not in experiment.channels:
                raise util.CytoflowOpError('settings',
                                           "Channel '{}' isn't in the experiment"
                                           .format(channel))
            for scale, params in scales.items():
                if scale not in SCALE_PARAMETERS:
                    raise util.CytoflowOpError('settings',
                                               "Can't set parameters for scale "
                                               "'{}'. Choose from {}"
                                               .format(scale, ", ".join(SCALE_PARAMETERS)))
                for param, value in params.items():
                    if param not in SCALE_PARAMETERS[scale]:
                        raise util.CytoflowOpError('settings',
                                                   "Scale '{}' has no parameter '{}'. "
                                                   "Choose from {}"
                                                   .format(scale, param,
                                                           ", ".join(SCALE_PARAMETERS[scale])))
                    _check_value(scale, param, value)

        new_experiment = experiment.clone(deep = False)

        for channel, scales in self.settings.items():
            md = new_experiment.metadata[channel]
            current = {s : dict(p) for s, p in md.get('scale_params', {}).items()}
            for scale, params in scales.items():
                current.setdefault(scale, {}).update(params)
            md['scale_params'] = current

        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment

    def default_view(self, **kwargs):
        """
        A histogram of one channel, to preview the settings.  Pass
        ``channel`` and ``scale``.
        """
        from cytoflow.views import HistogramView
        return HistogramView(**kwargs)


def _check_value(scale, param, value):
    ok = True
    if scale == 'logicle':
        if param in ('W', 'A'):
            ok = value >= 0
        elif param == 'M':
            ok = value > 0
        elif param == 'r':
            ok = 0 < value < 1
    elif scale in ('arcsinh', 'hlog', 'log'):
        ok = value > 0

    if not ok:
        raise util.CytoflowOpError('settings',
                                   "Bad value {} for {} parameter '{}'"
                                   .format(value, scale, param))

util.expand_class_attributes(ScaleSettingsOp)
util.expand_method_parameters(ScaleSettingsOp, ScaleSettingsOp.apply)
