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
cytoflow.operations._magnetic
-----------------------------

"Magnetic" (per-sample adjusted) gates, like FlowJo's magnetic gates: a
gate is drawn once (the template), then for each group of events (eg, each
tube) it is shifted so that it follows that group's population.  Used by 
`PolygonOp` / `FreehandOp` and `EllipseOp` (their ``adjust_by`` attribute.)

The shift is found with a mean-shift-like iteration in scaled space: in
the pooled data, the centroid of the events inside the template gate sits
at some position ``c0``.  For each group, the gate is moved by ``offset``,
the centroid ``c`` of the group's events inside the moved gate is computed,
and the offset becomes ``c - c0`` -- until it stops changing.  So the gate
ends up in the same place *relative to the population* as in the pooled
data.  The shift is limited to ``max_adjust`` times the gate's size on 
each axis, so a gate can't wander off to a different population.
"""

import numpy as np
import pandas as pd

import cytoflow.utility as util


def magnetic_membership(experiment, xy, contains, size, adjust_by, max_adjust,
                        max_iter = 20, tol = 1e-3):
    """
    Compute gate membership with a per-group shift.
    
    Parameters
    ----------
    experiment : Experiment
        Used for the grouping conditions.
        
    xy : numpy.ndarray, shape (n_events, 2)
        The (scaled) data, in the same order as ``experiment.data``.
        
    contains : callable
        ``contains(points, offset)`` -> boolean array: which ``points`` are
        in the gate when it's moved by ``offset`` (a length-2 array).
        
    size : (float, float)
        The gate's size (eg, its bounding box) on each axis, in scaled units.
        
    adjust_by : list of str
        The conditions defining the groups.
        
    max_adjust : float
        The largest shift allowed, as a fraction of ``size``.
        
    Returns
    -------
    (numpy.ndarray of bool, pandas.DataFrame)
        The membership of each event, and the shift applied to each group
        (columns ``x_shift``, ``y_shift``, in scaled units, indexed by the
        ``adjust_by`` conditions.)
    """
    
    for c in adjust_by:
        if c not in experiment.conditions:
            raise util.CytoflowOpError('adjust_by',
                                       "Can't find condition '{}' to adjust the "
                                       "gate by".format(c))
    if max_adjust < 0:
        raise util.CytoflowOpError('max_adjust', "max_adjust must be >= 0")
    
    finite = np.all(np.isfinite(xy), axis = 1)
    zero = np.zeros(2)
    
    in_template = contains(xy, zero) & finite
    if in_template.sum() == 0:
        raise util.CytoflowOpError('adjust_by',
                                   "There aren't any events in the gate, so it "
                                   "can't be adjusted.")
    c0 = xy[in_template].mean(axis = 0)
    limit = np.abs(np.asarray(size, dtype = float)) * max_adjust
    
    membership = np.zeros(len(xy), dtype = bool)
    shifts = {}
    
    groups = experiment.data.groupby(list(adjust_by), observed = True).indices
    for group, idx in groups.items():
        idx = np.asarray(idx)
        g_xy = xy[idx]
        g_finite = finite[idx]
        offset = zero.copy()
        for _ in range(max_iter):
            m = contains(g_xy, offset) & g_finite
            if m.sum() < 3:
                break
            new_offset = np.clip(g_xy[m].mean(axis = 0) - c0, -limit, limit)
            if np.all(np.abs(new_offset - offset) <= tol * np.maximum(limit, 1e-12)):
                offset = new_offset
                break
            offset = new_offset
            
        membership[idx] = contains(g_xy, offset) & g_finite
        key = group if isinstance(group, tuple) else (group,)
        shifts[key] = offset
        
    index = pd.MultiIndex.from_tuples(list(shifts.keys()), names = list(adjust_by)) \
            if len(adjust_by) > 1 \
            else pd.Index([k[0] for k in shifts.keys()], name = adjust_by[0])
    stat = pd.DataFrame([s for s in shifts.values()], 
                        index = index, 
                        columns = ['x_shift', 'y_shift']).sort_index()
    
    return membership, stat
