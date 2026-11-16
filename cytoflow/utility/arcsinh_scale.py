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
cytoflow.utility.arcsinh_scale
------------------------------

A scale that transforms data with the inverse hyperbolic sine, 
``arcsinh(x / cofactor)``.  This is the standard transform for mass 
cytometry (CyTOF, cofactor 5) and is widely used for spectral flow data 
(cofactor ~ 150 - 6000, depending on the instrument).  Like logicle, it is
linear near 0 (so it handles negative values) and logarithmic for large 
values; the ``cofactor`` sets where the transition happens.

`ArcsinhScale` -- implements `IScale`, the `cytoflow` interface for the scale.

`MatplotlibArcsinhScale` -- the matplotlib scale (built on matplotlib's own
``asinh`` scale, which provides the axis ticks.)
"""

import numpy as np
import pandas as pd

import matplotlib.scale
import matplotlib.colors

from traits.api import (Instance, Str, provides, Constant, Float, List, Array)

from .scale import IScale, ScaleMixin, register_scale
from .cytoflow_errors import CytoflowError

@provides(IScale)
class ArcsinhScale(ScaleMixin):
    """
    A scale that transforms the data with ``arcsinh(x / cofactor)``.
    
    Attributes
    ----------
    cofactor : Float (default = 150)
        Values much smaller than ``cofactor`` are displayed (nearly) 
        linearly; values much larger, (nearly) logarithmically.  Typical 
        values: 5 for mass cytometry; 150 for conventional flow data; larger
        (eg, 1000 - 6000) for some spectral cytometers.
    """
    
    id = Constant("cytoflow.utility.arcsinh_scale")
    name = "arcsinh"
    
    experiment = Instance("cytoflow.Experiment")
    
    # must set one of these.  they're considered in order.
    channel = Str
    condition = Str
    statistic = Str
    features = List(Str)
    data = Array
    
    cofactor = Float(150.0)
    
    def _check(self):
        if not self.cofactor > 0:
            raise CytoflowError("arcsinh cofactor must be > 0")
    
    def __call__(self, data):
        """Transforms ``data`` using this scale."""
        self._check()
        if isinstance(data, (int, float, np.number)):
            return float(np.arcsinh(data / self.cofactor))
        elif isinstance(data, (list, tuple)):
            ret = [self.__call__(x) for x in data]
            return tuple(ret) if isinstance(data, tuple) else ret
        elif isinstance(data, pd.Series):
            return np.arcsinh(data / self.cofactor)
        elif isinstance(data, np.ndarray):
            return np.arcsinh(data / self.cofactor)
        else:
            raise CytoflowError("Unknown type {} passed to ArcsinhScale.__call__"
                                .format(type(data)))
            
    def inverse(self, data):
        """Transforms ``data`` using the inverse of this scale."""
        self._check()
        if isinstance(data, (int, float, np.number)):
            return float(np.sinh(data) * self.cofactor)
        elif isinstance(data, (list, tuple)):
            ret = [self.inverse(x) for x in data]
            return tuple(ret) if isinstance(data, tuple) else ret
        elif isinstance(data, (np.ndarray, pd.Series)):
            return np.sinh(data) * self.cofactor
        else:
            raise CytoflowError("Unknown type {} passed to ArcsinhScale.inverse"
                                .format(type(data)))
            
    def clip(self, data):
        """arcsinh is defined everywhere, so there's nothing to clip."""
        return data
    
    def norm(self, vmin = None, vmax = None):
        """
        A factory function that returns a `matplotlib.colors.Normalize` 
        instance, which normalizes values for a `matplotlib` color palette.
        """
        if vmin is not None and vmax is not None:
            pass
        elif self.channel:
            vmin = self.experiment[self.channel].min()
            vmax = self.experiment[self.channel].max()
        elif self.condition:
            vmin = self.experiment[self.condition].min()
            vmax = self.experiment[self.condition].max()
        elif self.statistic in self.experiment.statistics:
            stat = self.experiment.statistics[self.statistic][self.features]
            vmin = stat.min().min()
            vmax = stat.max().max()
        elif self.data.size > 0:
            vmin = self.data.min()
            vmax = self.data.max()
        else:
            raise CytoflowError("Must set one of 'channel', 'condition' "
                                "or 'statistic'.")
            
        return matplotlib.colors.AsinhNorm(linear_width = self.cofactor,
                                           vmin = vmin, vmax = vmax)
    
    def get_mpl_params(self, ax):
        """
        Returns a dict with the parameters for `MatplotlibArcsinhScale`
        """
        return {"cofactor" : self.cofactor}
    
register_scale(ArcsinhScale)


class MatplotlibArcsinhScale(matplotlib.scale.AsinhScale):
    """
    matplotlib's ``asinh`` scale, which plots ``a * arcsinh(x / a)``, with
    ``a`` = the cofactor.  That's ``arcsinh(x / cofactor)`` times a constant,
    so shapes (eg, gates) look the same as in the `ArcsinhScale` space.
    """
    
    name = "arcsinh"
    
    def __init__(self, axis, *, cofactor = 150.0, **kwargs):
        super().__init__(axis, linear_width = cofactor, **kwargs)
        
matplotlib.scale.register_scale(MatplotlibArcsinhScale)
