#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

'''
Singlets Gate (Doublet Exclusion)
----------------------------------

Gate out doublets (two cells stuck together) using the FSC-H vs FSC-A
plot.  Single cells fall on a diagonal line; doublets fall below it and
to the right.

This is a Polygon gate with FSC-A on the X axis and FSC-H on the Y axis.
The channels are auto-selected if the experiment contains channels whose
names include ``FSC-A`` and ``FSC-H`` (case-insensitive match for
``fsc-a``/``fsc-h`` or ``FSC_A``/``FSC_H``).

Draw the gate the same way as any polygon gate: single-click to add
vertices, double-click (or click the first vertex) to close.

.. object:: Name

    The operation name.

.. object:: X Channel / Y Channel

    The scatter channels for the doublet plot.  Defaults to FSC-A (X) and
    FSC-H (Y).

.. object:: Adjust Per / Max Adjustment / Scale / Hue facet / Subset

    Same as for the standard Polygon Gate.
'''

from traits.api import provides, List
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..workflow.operations import PolygonWorkflowOp, PolygonSelectionView
from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import PluginHelpMixin
from .polygon import PolygonHandler, PolygonViewHandler

# Fraction of axis range used for the default diagonal gate vertices
# (0 = min, 1 = max of axis range)
_DIAG_FRACS = [0.02, 0.02, 0.02, 0.98, 0.98, 0.98, 0.98, 0.02]


def _fsc_channel(channels, suffix):
    """Return the channel name containing the given FSC suffix."""
    suffix_lower = suffix.lower()
    for ch in channels:
        if ch.lower().replace('_', '-') == 'fsc-' + suffix_lower.lstrip('-'):
            return ch
        if suffix_lower in ch.lower():
            return ch
    return ''


@provides(IOperationPlugin)
class DoubletGatePlugin(Plugin, PluginHelpMixin):
    """
    Doublet-exclusion (singlets) gate: a Polygon gate pre-configured for
    FSC-H vs FSC-A.
    """

    id = 'cytoflowgui.op_plugins.doublet_gate'
    operation_id = 'cytoflow.operations.polygon'   # reuses PolygonOp
    view_id = 'cytoflowgui.workflow.operations.polygonview'

    name = 'Singlets Gate'
    short_name = 'Singlets'
    menu_group = 'Gates'

    def get_operation(self):
        op = PolygonWorkflowOp()
        # Pre-select FSC-A / FSC-H — will be overridden once an experiment
        # is available; these are just reasonable defaults.
        op.xchannel = 'FSC-A'
        op.ychannel = 'FSC-H'
        return op

    def get_handler(self, model, context):
        if isinstance(model, PolygonWorkflowOp):
            return PolygonHandler(model=model, context=context)
        elif isinstance(model, PolygonSelectionView):
            return PolygonViewHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('polygon')

    plugin = List(contributes_to=OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
