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
cytoflow.views.highlight
------------------------

"Back-gating": draw a population (eg, the events in a gate) on top of all 
the events, in a contrasting color.  Used by `ScatterplotView` and 
`DensityView` (the ``highlight`` attribute).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cytoflow.utility as util

HIGHLIGHT_COLOR = '#e4572e'
BACKGROUND_COLOR = '#b0b0b0'


def highlight_mask(experiment, highlight):
    """
    Return a boolean `pandas.Series` (indexed like ``experiment.data``) 
    that's True for the events to highlight.
    
    ``highlight`` is either the name of a boolean condition (eg, the name of
    a gate) or an expression (eg, ``"Quad == 'Quad_2'"`` or 
    ``"CD4_pos and not Dead"``).
    """
    
    if highlight in experiment.conditions:
        col = experiment.data[highlight]
        if not pd.api.types.is_bool_dtype(col):
            raise util.CytoflowViewError('highlight',
                                         "Condition '{0}' isn't True/False. Use an "
                                         "expression to choose which events to "
                                         "highlight, eg \"{0} == {1!r}\""
                                         .format(highlight, col.iloc[0]))
        return col
    
    try:
        mask = experiment.data.eval(highlight)
    except Exception as e:
        raise util.CytoflowViewError('highlight',
                                     "Couldn't evaluate highlight expression '{}': {}"
                                     .format(highlight, e)) from e
        
    if not isinstance(mask, pd.Series) or not pd.api.types.is_bool_dtype(mask):
        raise util.CytoflowViewError('highlight',
                                     "Highlight expression '{}' must be True or "
                                     "False for each event".format(highlight))
    return mask


def draw_highlight(x, y, mask = None, highlight_kwargs = None, **kwargs):
    """
    Draw the highlighted events in one facet.  Called through
    ``FacetGrid.map``, so ``x`` and ``y`` are the facet's data and ``kwargs``
    are seaborn's (color, label), which are ignored.
    """
    m = mask.reindex(x.index).fillna(False).to_numpy(dtype = bool)
    if not m.any():
        return
    plt.scatter(x.to_numpy()[m], y.to_numpy()[m], **(highlight_kwargs or {}))


def pop_highlight_kwargs(kwargs):
    """
    Pull the highlight plot parameters out of ``kwargs``:
    ``highlight_color``, ``highlight_size`` and ``highlight_alpha``.
    """
    return dict(color = kwargs.pop('highlight_color', HIGHLIGHT_COLOR),
                s = kwargs.pop('highlight_size', 4),
                alpha = kwargs.pop('highlight_alpha', 0.8),
                marker = 'o',
                linewidths = 0,
                zorder = 3,
                rasterized = None)
